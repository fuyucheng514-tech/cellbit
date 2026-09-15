#include "bubbles.hpp"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace subass {
namespace {
constexpr int SOLID_LEN=10, SIMPLE_LEN=4, MAX_BUBBLE=500, MAX_BRANCHES=50;
constexpr int MAX_COVERAGE=1000, MIN_POLISH_ALN=500;
constexpr double SOLID_MISMATCH=.3, SOLID_INDEL=.3;
struct Aln { std::string qid, qseq, tseq; int qs=0,qe=0,ts=0,te=0,tl=0,mapq=0; bool secondary=false; double err=0; };
using InsertionMap = std::unordered_map<std::string,std::string>;
// Most reference positions have no insertion at all.  Keeping a full empty
// unordered_map at every base dominates memory on long contigs, so allocate it
// only at positions where an insertion is actually observed.
struct Profile { std::unique_ptr<InsertionMap> ins; int prop=0,del=0,mis=0,cov=0; };
struct Bubble { std::string id,cons; int pos=0,sub=0; std::vector<std::string> branches; };

std::map<std::string,std::string> fasta(const std::filesystem::path& p) {
  std::ifstream in(p); if(!in) throw std::runtime_error("cannot open FASTA "+p.string());
  std::map<std::string,std::string> r; std::string l,id;
  while(std::getline(in,l)) { if(!l.empty()&&l.back()=='\r') l.pop_back(); if(l.empty()) continue;
    if(l[0]=='>') { id=l.substr(1); auto s=id.find_first_of(" \t"); if(s!=std::string::npos) id.resize(s); if(!r.emplace(id,"").second) throw std::runtime_error("duplicate FASTA id "+id); }
    else { if(id.empty()) throw std::runtime_error("sequence before header"); for(char c:l) if(!std::isspace((unsigned char)c)) r[id]+=std::toupper((unsigned char)c); }
  } return r;
}
char bam_base(const uint8_t* seq,int i){ static const char* n="=ACMGRSVTWYHKDBN"; return n[bam_seqi(seq,i)]; }

Aln decode(const bam1_t* b,const bam_hdr_t* h,const std::map<std::string,std::string>& refs){
  Aln a; a.qid=bam_get_qname(b); const std::string target=h->target_name[b->core.tid]; a.tl=h->target_len[b->core.tid]; a.ts=b->core.pos; a.mapq=b->core.qual;
  a.secondary=b->core.flag&BAM_FSECONDARY;
  // BAM stores SEQ in the alignment orientation; match Flye's SAM parser and
  // do not reverse-complement it here.  The reverse flag is metadata only.
  std::string read; read.reserve(b->core.l_qseq); auto* sq=bam_get_seq(b); for(int i=0;i<b->core.l_qseq;i++) read+=bam_base(sq,i);
  auto it=refs.find(target); if(it==refs.end()) throw std::runtime_error("BAM target absent in FASTA: "+target); const std::string& ref=it->second;
  a.qseq.reserve(read.size()); a.tseq.reserve(read.size());
  int rp=a.ts,qp=0,hardL=0,softL=0,softR=0; bool leftH=true,leftS=true; auto* cg=bam_get_cigar(b);
  for(uint32_t i=0;i<b->core.n_cigar;i++){ int op=bam_cigar_op(cg[i]), n=bam_cigar_oplen(cg[i]);
    if(op==BAM_CHARD_CLIP){if(leftH)hardL+=n;}
    else if(op==BAM_CSOFT_CLIP){qp+=n;if(leftS)softL+=n;else softR+=n;}
    else if(op==BAM_CMATCH||op==BAM_CEQUAL||op==BAM_CDIFF){a.qseq.append(read,qp,n);a.tseq.append(ref,rp,n);qp+=n;rp+=n;}
    else if(op==BAM_CINS){a.qseq.append(read,qp,n);a.tseq.append(n,'-');qp+=n;}
    else if(op==BAM_CDEL){a.qseq.append(n,'-');a.tseq.append(ref,rp,n);rp+=n;}
    else if(op==BAM_CREF_SKIP){rp+=n;} else if(op!=BAM_CPAD) throw std::runtime_error("unsupported CIGAR op");
    leftH=false; if(op!=BAM_CHARD_CLIP)leftS=false;
  }
  a.te=rp; a.qs=hardL+softL; a.qe=qp+hardL-softR; int match=0; for(size_t i=0;i<a.tseq.size();i++) match+=a.tseq[i]==a.qseq[i]; a.err=a.tseq.empty()?1.0:1.0-double(match)/a.tseq.size(); return a;
}

class BamReader {
public:
  BamReader(const std::filesystem::path& bam,const std::map<std::string,std::string>& refs)
      : path_(bam.string()), refs_(refs) {
    file_=sam_open(path_.c_str(),"r");
    if(!file_)throw std::runtime_error("open BAM failed");
    header_=sam_hdr_read(file_);
    if(!header_){sam_close(file_);file_=nullptr;throw std::runtime_error("BAM header required");}
    index_=sam_index_load(file_,path_.c_str());
    if(!index_){bam_hdr_destroy(header_);header_=nullptr;sam_close(file_);file_=nullptr;throw std::runtime_error("indexed BAM required");}
    record_=bam_init1();
    if(!record_){hts_idx_destroy(index_);index_=nullptr;bam_hdr_destroy(header_);header_=nullptr;sam_close(file_);file_=nullptr;throw std::runtime_error("cannot allocate BAM record");}
  }
  ~BamReader(){bam_destroy1(record_);hts_idx_destroy(index_);bam_hdr_destroy(header_);if(file_)sam_close(file_);}
  BamReader(const BamReader&)=delete;
  BamReader& operator=(const BamReader&)=delete;

  std::vector<Aln> alignments(const std::string& id){
    int tid=bam_name2id(header_,id.c_str());
    if(tid<0)throw std::runtime_error("FASTA target absent in BAM: "+id);
    using Iterator=std::unique_ptr<hts_itr_t,decltype(&hts_itr_destroy)>;
    Iterator itr(sam_itr_queryi(index_,tid,0,header_->target_len[tid]),&hts_itr_destroy);
    if(!itr)throw std::runtime_error("cannot query indexed BAM target: "+id);
    std::vector<Aln> v;
    while(sam_itr_next(file_,itr.get(),record_)>=0){
      if((record_->core.flag&BAM_FUNMAP)||record_->core.l_qseq==0)continue;
      v.push_back(decode(record_,header_,refs_));
    }
    std::mt19937 gen(42); std::shuffle(v.begin(),v.end(),gen);
    std::vector<Aln> capped; long long bases=0; int len=refs_.at(id).size();
    for(auto& a:v){capped.push_back(std::move(a));bases+=capped.back().qe-capped.back().qs;if(bases/std::max(1,len)>MAX_COVERAGE)break;}
    std::sort(capped.begin(),capped.end(),[](const auto&x,const auto&y){return x.qid==y.qid?(x.qe-x.qs)>(y.qe-y.qs):x.qid<y.qid;});
    return capped;
  }

private:
  std::string path_;
  const std::map<std::string,std::string>& refs_;
  samFile* file_=nullptr;
  bam_hdr_t* header_=nullptr;
  hts_idx_t* index_=nullptr;
  bam1_t* record_=nullptr;
};
double median(std::vector<int> x){if(x.empty())return 0;std::sort(x.begin(),x.end());return x[x.size()/2];}
std::pair<std::vector<Aln>,int> uniform(std::vector<Aln> a){ if(a.empty())return{{},0}; int n=a[0].tl/100+1; std::vector<int> cov(n); auto reliable=[](const Aln&x){return !x.secondary&&x.mapq>=20;};
  for(auto&x:a)if(reliable(x))for(int i=x.ts/100;i<=x.te/100&&i<n;i++)cov[i]++;
  int thr=std::max(20,(int)median(cov));std::vector<Aln> out;std::unordered_map<std::string,Aln> sec;
  for(auto&x:a){if(reliable(x))out.push_back(std::move(x));else{std::string key=x.qid;sec.insert_or_assign(std::move(key),std::move(x));}}
  auto score=[&](const Aln&x){int g=0,b=0;for(int i=x.ts/100;i<=x.te/100&&i<n;i++)(cov[i]<thr?g:b)++;return std::pair<int,int>{g-2*b,x.te-x.ts};};std::vector<Aln> s;s.reserve(sec.size());for(auto&kv:sec)s.push_back(std::move(kv.second));sec.clear();sec.rehash(0);std::sort(s.begin(),s.end(),[&](const auto&x,const auto&y){return score(x)>score(y);});
  for(auto&x:s){int g=0,b=0;for(int i=x.ts/100;i<=x.te/100&&i<n;i++)(cov[i]<thr?g:b)++;if(g+b&&double(g)/(g+b)>.66){out.push_back(std::move(x));for(int i=out.back().ts/100;i<=out.back().te/100&&i<n;i++)cov[i]++;}}
  return{std::move(out),(int)median(cov)};
}
std::string shift(std::string trg,std::string qry){trg.reserve(trg.size()+2);trg.insert(trg.begin(),'$');trg.push_back('$');qry.reserve(qry.size()+2);qry.insert(qry.begin(),'$');qry.push_back('$');bool gap=false;size_t st=0;for(size_t i=0;i<trg.size();i++){if(gap&&qry[i]!='-'){gap=false;int l=(int)st-1,r=(int)i-1;while(l>0&&r>=(int)st&&qry[l]==trg[r])std::swap(qry[l--],qry[r--]);}if(!gap&&qry[i]=='-'){gap=true;st=i;}}qry.erase(qry.begin());qry.pop_back();return qry;}
std::vector<Profile> profile(const std::vector<Aln>& as,const std::string& ref,double& error_sum,size_t& error_count){std::vector<Profile> p(ref.size());int minlen=std::min<int>(MIN_POLISH_ALN,ref.size()/2);
 for(auto&a:as){if((int)a.qseq.size()<minlen)continue;error_sum+=a.err;++error_count;auto q=shift(a.tseq,a.qseq),t=shift(q,a.tseq);int pos=a.ts;for(size_t i=0;i<t.size();i++){if(t[i]=='-')pos--;if(pos<0||pos>=(int)p.size()){pos++;continue;}auto&z=p[pos];if(t[i]=='-'){if(!z.ins)z.ins=std::make_unique<InsertionMap>();(*z.ins)[a.qid]+=q[i];}else{z.cov++;if(q[i]=='-')z.del++;else if(t[i]!=q[i])z.mis++;}pos++;}}
 for(int i=0;i<(int)p.size();i++)if(p[i].ins)for(auto&kv:*p[i].ins){p[i].prop++;int n=kv.second.size();for(int j=std::max(0,i-n);j<i;j++)p[j].prop++;for(int j=i+1;j<std::min<int>(p.size(),i+n+1);j++)p[j].prop++;}
 return p;}
bool solid(const std::vector<Profile>&p,int x){for(int i=x;i<x+SOLID_LEN;i++){if(!p[i].cov)return false;if(double(p[i].mis+p[i].del)/p[i].cov>SOLID_MISMATCH||double(p[i].prop)/p[i].cov>SOLID_INDEL)return false;}return true;}
bool simple(const std::string&ref,int x){int ext=SIMPLE_LEN*2;std::string s;for(int i=x-ext/2;i<x+ext/2;i++)s+=ref[i];for(int i=ext/2-SIMPLE_LEN/2;i<ext/2+SIMPLE_LEN/2-1;i++)if(s[i]==s[i+1])return false;for(int sh:{0,1})for(int i=0;i<SIMPLE_LEN-sh-1;i++){int z=ext/2-SIMPLE_LEN+sh+i*2;if(s.substr(z,2)==s.substr(z+2,2))return false;}return true;}
std::vector<int> partition(const std::vector<Profile>&p,const std::string&ref){std::vector<char> f(p.size());for(int x=0;x<(int)p.size()-SOLID_LEN;){if(solid(p,x)){for(int i=x;i<x+SOLID_LEN;i++)f[i]=1;x+=SOLID_LEN;}else x++;}std::vector<int> r;int prev=SOLID_LEN;for(int x=SOLID_LEN;x<(int)p.size()-SOLID_LEN;){int cur=x+SIMPLE_LEN/2;bool mark=std::all_of(f.begin()+x,f.begin()+x+SIMPLE_LEN,[](char c){return c;})&&simple(ref,cur);if(mark||x-prev>MAX_BUBBLE){r.push_back(cur);prev=cur;x+=SOLID_LEN;}else x++;}return r;}
std::vector<Bubble> bubbles(const std::vector<Aln>&as,const std::string&ref,const std::vector<int>&part,const std::string&id){if(part.empty()||as.empty())return{};std::vector<int> ep{0};ep.insert(ep.end(),part.begin(),part.end());ep.push_back(ref.size());std::vector<Bubble>b;b.reserve(ep.size()-1);for(size_t i=0;i+1<ep.size();i++){Bubble z;z.id=id;z.pos=ep[i];z.cons.assign(ref,ep[i],ep[i+1]-ep[i]);b.push_back(std::move(z));}
 for(auto&a:as){int bi=std::upper_bound(ep.begin(),ep.end(),a.ts)-ep.begin()-1,next=ep[bi+1],tp=a.ts,bs=0;bool incomplete=a.ts>ep[bi];for(int i=0;i<(int)a.tseq.size();i++){if(a.tseq[i]=='-')continue;if(tp>=next){if(!incomplete){std::string q;q.reserve(i-bs);for(int j=bs;j<i;j++)if(a.qseq[j]!='-')q+=a.qseq[j];b[bi].branches.push_back(std::move(q));}incomplete=false;bi=std::upper_bound(ep.begin(),ep.end(),tp)-ep.begin()-1;next=ep[bi+1];bs=i;}tp++;}if(a.te>=ep.back()){std::string q;q.reserve(a.qseq.size()-bs);for(int j=bs;j<(int)a.qseq.size();j++)if(a.qseq[j]!='-')q+=a.qseq[j];b.back().branches.push_back(std::move(q));}}
std::vector<Bubble> out;for(auto&z:b){if(z.branches.empty())continue;std::vector<const std::string*> lens;lens.reserve(z.branches.size());for(const auto&s:z.branches)lens.push_back(&s);std::sort(lens.begin(),lens.end(),[](const auto*a,const auto*b){return a->size()<b->size();});const std::string* med=lens[lens.size()/2];const size_t med_size=med->size();if(!med_size)continue;if(std::abs((long)med_size-(long)z.cons.size())>(long)med_size/2)z.cons=*med;std::vector<std::string> good;for(auto&s:z.branches)if(!s.empty()&&std::abs((long)s.size()-(long)med_size)/double(med_size)<.5)good.push_back(std::move(s));if(good.empty())continue;if(good.size()>MAX_BRANCHES){std::sort(good.begin(),good.end(),[](const auto&a,const auto&b){return a.size()<b.size();});int l=good.size()/2-MAX_BRANCHES/2;std::vector<std::string> mid;mid.reserve(MAX_BRANCHES);for(int i=l;i<l+MAX_BRANCHES;i++)mid.push_back(std::move(good[i]));good.swap(mid);}z.branches=std::move(good);int chunks=med_size/MAX_BUBBLE;if(chunks>1){for(int c=0;c<chunks;c++){Bubble y;y.id=z.id;y.pos=z.pos;y.sub=c;for(auto&s:z.branches){int n=s.size()/chunks,st=c*n,en=c==chunks-1?s.size():(c+1)*n;y.branches.push_back(s.substr(st,en-st));}y.cons=y.branches[0];out.push_back(std::move(y));}}else out.push_back(std::move(z));}return out;}
}
CoverageResult make_bubbles(const std::filesystem::path& bam,const std::filesystem::path& reference,const std::filesystem::path& output,int){auto refs=fasta(reference);std::ofstream out(output);if(!out)throw std::runtime_error("cannot create bubbles");CoverageResult r;double error_sum=0;size_t error_count=0;BamReader reader(bam,refs);
 for(auto&kv:refs){std::vector<Bubble> bs;{auto raw=reader.alignments(kv.first);if(raw.empty())continue;auto u=uniform(std::move(raw));r.mean_coverage[kv.first]=u.second;auto p=profile(u.first,kv.second,error_sum,error_count);auto part=partition(p,kv.second);bs=bubbles(u.first,kv.second,part,kv.first);}for(auto&b:bs){out<<'>'<<b.id<<' '<<b.pos<<' '<<b.branches.size()<<' '<<b.sub<<'\n'<<b.cons<<'\n';for(size_t i=0;i<b.branches.size();i++)out<<'>'<<i<<'\n'<<b.branches[i]<<'\n';}}
 r.mean_alignment_error=error_count?error_sum/(error_count+1):0;return r;}
}
