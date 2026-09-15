#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

struct Config {
  int k=15, window=500, min_window=100, min_seed_hits=3, max_occ=64, align_band=128;
  int sketch_scale=64, query_sketch_size=1000, max_posting=1000;
  int min_candidate_hits=2, max_candidates=20, top=5, threads=1;
  int triangle_min_shared=1, triangle_max_candidates=128;
  int triangle_sketch_size=0, triangle_max_posting=1000;
  // The fused block implementation is retained as an explicitly selectable
  // experiment, but the exact production default stays on the byte-exact
  // reference-at-a-time path until a canary demonstrates a wall-time win.
  int triangle_reference_block_size=1;
  // Performance-only directory width for the exact global seed join.  It
  // changes lookup granularity, never which seed occurrences are retained.
  int triangle_global_prefix_bits=20;
  double min_identity=70.0, min_af=0.0;
  bool k_explicit=false, sketch_scale_explicit=false, query_sketch_explicit=false, max_posting_explicit=false;
  bool deep_verify=false;
  bool compact_index=true;
  bool triangle_fast=false;
  bool triangle_candidates_only=false;
  bool triangle_global_seed_join=false;
  fs::path query, reference, list, manifest, taxonomy, radii, release_metadata;
  fs::path out_dir, index, queries, out, triangle_candidates_out, stats_out;
};
struct Contig { std::string name, seq; };
struct Genome {
  fs::path path; std::vector<Contig> contigs; std::uint64_t callable=0;
  std::uint64_t source_bytes=0; std::int64_t source_mtime=0;
  std::string content_sha256;
};
struct Estimate { double ani=0, af_ref=0, af_query=0; std::uint64_t aligned=0, matches=0; };
struct Occ { std::uint32_t contig=0, pos=0; };
struct Candidate {
  std::uint32_t qc=0, rc=0, q0=0, q1=0; int diag=0; bool reverse=false;
  int seed_votes=0, compared=0, matches=0; double identity=0;
  std::vector<std::pair<std::uint32_t,std::uint32_t>> aligned;
};
struct RefRow {
  std::uint32_t id=0; fs::path path; std::uint64_t bytes=0, callable=0;
  std::int64_t mtime=0; std::string content_sha256;
  std::uint32_t sketch_count=0, posting_count=0;
};

// Exact triangle accelerator.  The original estimator rebuilds a retained
// seed table for reference i and then decodes every later query once for that
// reference.  For an N-node triangle that repeatedly decodes each query O(N)
// times.  This compact global posting table contains precisely the same
// per-reference retained forward-seed occurrences.  A query can therefore be
// decoded once and its hits dispatched to every earlier reference without
// changing the downstream vote keys, gapped alignments, ANI, or AF formulas.
struct GlobalSeedOcc {
  std::uint32_t seed=0;
  std::uint16_t reference=0, contig=0;
  std::uint32_t pos=0;
};
static_assert(sizeof(GlobalSeedOcc)==12,"unexpected GlobalSeedOcc padding");
struct GlobalSeedIndex {
  int key_bits=0,prefix_bits=0,prefix_shift=0;
  // The prefix directory narrows an exact lower_bound to one contiguous seed
  // bucket.  Its width is tunable only for performance canaries; lookup and
  // retained occurrence semantics are invariant.
  std::vector<std::uint32_t> offsets;
  std::vector<GlobalSeedOcc> occurrences;
  std::pair<const GlobalSeedOcc*,const GlobalSeedOcc*> equal_range(std::uint64_t raw_seed)const{
    const auto seed=std::uint32_t(raw_seed);
    const auto prefix=prefix_shift?seed>>prefix_shift:seed;
    const auto begin_index=offsets[prefix],end_index=offsets[std::size_t(prefix)+1];
    const auto*begin=occurrences.data()+begin_index;const auto*end=occurrences.data()+end_index;
    auto first=std::lower_bound(begin,end,seed,[](const GlobalSeedOcc&x,std::uint32_t value){return x.seed<value;});
    auto last=std::upper_bound(first,end,seed,[](std::uint32_t value,const GlobalSeedOcc&x){return value<x.seed;});
    return {first,last};
  }
};

static std::string trim(std::string s){while(!s.empty()&&std::isspace(static_cast<unsigned char>(s.front())))s.erase(s.begin());while(!s.empty()&&std::isspace(static_cast<unsigned char>(s.back())))s.pop_back();return s;}
static std::vector<std::string> split(const std::string&s,char d){std::vector<std::string>v;std::stringstream x(s);std::string z;while(std::getline(x,z,d))v.push_back(z);if(!s.empty()&&s.back()==d)v.emplace_back();return v;}
static std::string json_escape(const std::string&s){std::string r;for(unsigned char c:s){if(c=='\\')r+="\\\\";else if(c=='"')r+="\\\"";else if(c=='\n')r+="\\n";else if(c=='\r')r+="\\r";else if(c=='\t')r+="\\t";else r+=char(c);}return r;}
static int base_code(char c){switch(c){case'A':return 0;case'C':return 1;case'G':return 2;case'T':return 3;default:return -1;}}
static bool callable(char c){return base_code(c)>=0;}
static char complement(char c){switch(c){case'A':return'T';case'C':return'G';case'G':return'C';case'T':return'A';default:return'N';}}
static std::string reverse_complement(const std::string&s){std::string r(s.size(),'N');for(size_t i=0;i<s.size();++i)r[s.size()-1-i]=complement(s[i]);return r;}
static std::uint64_t mix64(std::uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}

class Sha256 {
  std::uint32_t st[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};std::uint8_t b[64]={};size_t used=0;std::uint64_t bits=0;
  static std::uint32_t rr(std::uint32_t x,int n){return(x>>n)|(x<<(32-n));}
  void block(){static const std::uint32_t k[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};std::uint32_t w[64];for(int i=0;i<16;++i)w[i]=(std::uint32_t(b[4*i])<<24)|(std::uint32_t(b[4*i+1])<<16)|(std::uint32_t(b[4*i+2])<<8)|b[4*i+3];for(int i=16;i<64;++i){auto s0=rr(w[i-15],7)^rr(w[i-15],18)^(w[i-15]>>3),s1=rr(w[i-2],17)^rr(w[i-2],19)^(w[i-2]>>10);w[i]=w[i-16]+s0+w[i-7]+s1;}auto a=st[0],c1=st[1],c2=st[2],d=st[3],e=st[4],f=st[5],g=st[6],h=st[7];for(int i=0;i<64;++i){auto s1=rr(e,6)^rr(e,11)^rr(e,25),ch=(e&f)^((~e)&g),t1=h+s1+ch+k[i]+w[i],s0=rr(a,2)^rr(a,13)^rr(a,22),maj=(a&c1)^(a&c2)^(c1&c2),t2=s0+maj;h=g;g=f;f=e;e=d+t1;d=c2;c2=c1;c1=a;a=t1+t2;}st[0]+=a;st[1]+=c1;st[2]+=c2;st[3]+=d;st[4]+=e;st[5]+=f;st[6]+=g;st[7]+=h;}
public:void update(const void*p,size_t n){auto*x=static_cast<const std::uint8_t*>(p);bits+=std::uint64_t(n)*8;while(n){size_t take=std::min(n,64-used);std::copy(x,x+take,b+used);used+=take;x+=take;n-=take;if(used==64){block();used=0;}}}void update(const std::string&s){update(s.data(),s.size());}std::string finish(){auto original=bits;std::uint8_t one=0x80,zero=0;update(&one,1);while(used!=56)update(&zero,1);std::uint8_t len[8];for(int i=0;i<8;++i)len[7-i]=std::uint8_t(original>>(8*i));update(len,8);std::ostringstream o;o<<std::hex<<std::setfill('0');for(auto x:st)o<<std::setw(8)<<x;return o.str();}};
static std::string sha_text(const std::string&s){Sha256 h;h.update(s);return h.finish();}
static std::string sha_file(const fs::path&p){std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot hash: "+p.string());Sha256 h;char b[1<<16];while(f){f.read(b,sizeof(b));if(f.gcount())h.update(b,size_t(f.gcount()));}return h.finish();}
static std::uint64_t exact_pair_count(std::size_t n){
  if(n<2)return 0;
  std::uint64_t a=static_cast<std::uint64_t>(n),b=a-1;
  if(a/2*2==a)a/=2;else b/=2;
  if(b&&a>std::numeric_limits<std::uint64_t>::max()/b)throw std::runtime_error("triangle pair count overflows uint64");
  return a*b;
}
static std::string triangle_input_content_set_sha256(const std::vector<fs::path>&paths){
  Sha256 h;h.update("gtdb-ani-af-triangle-input-content-set-v1\n");
  for(const auto&p:paths){auto canonical=fs::weakly_canonical(p);h.update(canonical.string());h.update("\t");h.update(sha_file(canonical));h.update("\n");}
  return h.finish();
}

static std::string read_source(const fs::path&p){
  auto ext=p.extension().string();for(char&c:ext)c=char(std::tolower(static_cast<unsigned char>(c)));if(ext==".gz"){gzFile z=gzopen(p.string().c_str(),"rb");if(!z)throw std::runtime_error("cannot open gzip FASTA: "+p.string());std::string out;char b[1<<16];int n=0;while((n=gzread(z,b,sizeof(b)))>0)out.append(b,n);int err=Z_OK;const char*msg=gzerror(z,&err);if(err!=Z_OK&&err!=Z_STREAM_END){std::string m=msg?msg:"gzip read error";gzclose(z);throw std::runtime_error(m);}if(gzclose(z)!=Z_OK)throw std::runtime_error("gzip close/checksum failure: "+p.string());return out;}
  std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot read FASTA: "+p.string());std::ostringstream o;o<<f.rdbuf();return o.str();
}
static Genome parse_fasta_text(const fs::path&p,const std::string&text){Genome g;g.path=fs::absolute(p).lexically_normal();std::istringstream f(text);std::string line;Contig*cur=nullptr;std::set<std::string>names;while(std::getline(f,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();if(!line.empty()&&line[0]=='>'){auto name=trim(line.substr(1));if(name.empty()||!names.insert(name).second)throw std::runtime_error("empty/duplicate FASTA header: "+p.string());g.contigs.push_back({name,{}});cur=&g.contigs.back();continue;}if(!cur&&!trim(line).empty())throw std::runtime_error("sequence before FASTA header: "+p.string());if(!cur)continue;for(unsigned char x:line)if(!std::isspace(x)){char c=char(std::toupper(x));if(c=='U')c='T';if(base_code(c)<0)c='N';else++g.callable;cur->seq+=c;}if(cur->seq.size()>size_t(std::numeric_limits<int>::max()))throw std::runtime_error("FASTA contig exceeds supported coordinate range: "+p.string());}g.contigs.erase(std::remove_if(g.contigs.begin(),g.contigs.end(),[](const Contig&c){return c.seq.empty();}),g.contigs.end());if(g.contigs.empty()||!g.callable||g.contigs.size()>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("FASTA has no callable bases or exceeds coordinate range: "+p.string());g.content_sha256=sha_text(text);return g;}
static Genome read_fasta(const fs::path&p){if(!fs::is_regular_file(p))throw std::runtime_error("FASTA missing: "+p.string());auto text=read_source(p);auto g=parse_fasta_text(fs::weakly_canonical(p),text);g.source_bytes=fs::file_size(p);g.source_mtime=fs::last_write_time(p).time_since_epoch().count();return g;}
static std::string triangle_loaded_content_set_sha256(const std::vector<Genome>&genomes){
  Sha256 h;h.update("gtdb-ani-af-triangle-input-content-set-v1\n");
  for(const auto&g:genomes){h.update(g.path.string());h.update("\t");h.update(g.content_sha256);h.update("\n");}
  return h.finish();
}

template<class F>static void each_kmer(const std::string&s,int k,F fn){std::uint64_t mask=(1ULL<<(2*k))-1,fw=0,rv=0;int valid=0;for(size_t i=0;i<s.size();++i){int x=base_code(s[i]);if(x<0){fw=rv=0;valid=0;continue;}fw=((fw<<2)|std::uint64_t(x))&mask;rv=(rv>>2)|(std::uint64_t(3-x)<<(2*(k-1)));if(++valid>=k)fn(std::uint32_t(i-k+1),fw,rv);}}
static std::vector<std::uint64_t> make_scaled_sketch(const Genome&g,int k,int scale,int cap){
  const auto threshold=std::numeric_limits<std::uint64_t>::max()/std::uint64_t(scale);
  // A flat vector has much lower allocator overhead than one unordered_set
  // node per retained k-mer. Sorting and uniquing produces exactly the same
  // ordered sketch as the former set implementation.
  std::vector<std::uint64_t>x;
  const auto expected=std::uint64_t(scale)>0?g.callable/std::uint64_t(scale):g.callable;
  if(expected<std::uint64_t(std::numeric_limits<size_t>::max()))x.reserve(size_t(expected));
  for(const auto&c:g.contigs)each_kmer(c.seq,k,[&](std::uint32_t,std::uint64_t fw,std::uint64_t rv){auto h=mix64(std::min(fw,rv));if(h<=threshold)x.push_back(h);});
  std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());if(cap>0&&x.size()>size_t(cap))x.resize(size_t(cap));return x;
}

using SeedIndex=std::unordered_map<std::uint64_t,std::vector<Occ>>;
static SeedIndex reference_seeds(const Genome&r,const Config&c){SeedIndex idx;if(r.callable<std::uint64_t(std::numeric_limits<size_t>::max()))idx.reserve(size_t(r.callable));for(size_t ci=0;ci<r.contigs.size();++ci)each_kmer(r.contigs[ci].seq,c.k,[&](std::uint32_t p,std::uint64_t fw,std::uint64_t){auto&v=idx[fw];if(int(v.size())<=c.max_occ)v.push_back({std::uint32_t(ci),p});});for(auto i=idx.begin();i!=idx.end();)if(int(i->second.size())>c.max_occ)i=idx.erase(i);else++i;return idx;}

static GlobalSeedIndex build_global_seed_index(const std::vector<Genome>&genomes,const Config&c,std::size_t threads){
  if(c.k>15)throw std::runtime_error("global exact seed join currently requires k<=15");
  if(genomes.size()>std::size_t(std::numeric_limits<std::uint16_t>::max())+1)throw std::runtime_error("global exact seed join supports at most 65536 references");
  std::vector<std::vector<GlobalSeedOcc>> local(genomes.size());std::atomic<std::size_t>next{0};std::exception_ptr error;std::mutex error_mutex;
  auto build_one=[&](){for(;;){const auto reference=next.fetch_add(1);if(reference>=genomes.size())return;try{
    const auto&g=genomes[reference];std::vector<GlobalSeedOcc>rows;if(g.callable<std::uint64_t(std::numeric_limits<std::size_t>::max()))rows.reserve(std::size_t(g.callable));
    if(g.contigs.size()>std::size_t(std::numeric_limits<std::uint16_t>::max())+1)throw std::runtime_error("global exact seed join supports at most 65536 contigs per reference");
    for(std::size_t contig=0;contig<g.contigs.size();++contig)each_kmer(g.contigs[contig].seq,c.k,[&](std::uint32_t pos,std::uint64_t fw,std::uint64_t){rows.push_back({std::uint32_t(fw),std::uint16_t(reference),std::uint16_t(contig),pos});});
    std::sort(rows.begin(),rows.end(),[](const GlobalSeedOcc&a,const GlobalSeedOcc&b){return std::tie(a.seed,a.contig,a.pos)<std::tie(b.seed,b.contig,b.pos);});
    std::size_t write=0;for(std::size_t begin=0;begin<rows.size();){std::size_t end=begin+1;while(end<rows.size()&&rows[end].seed==rows[begin].seed)++end;if(end-begin<=std::size_t(c.max_occ)){for(auto i=begin;i<end;++i)rows[write++]=rows[i];}begin=end;}rows.resize(write);local[reference]=std::move(rows);
  }catch(...){std::lock_guard<std::mutex>guard(error_mutex);if(!error)error=std::current_exception();return;}}};
  std::vector<std::thread>pool;const auto worker_count=std::min<std::size_t>(std::max<std::size_t>(1,threads),genomes.size());for(std::size_t i=0;i<worker_count;++i)pool.emplace_back(build_one);for(auto&t:pool)t.join();if(error)std::rethrow_exception(error);

  GlobalSeedIndex index;index.key_bits=2*c.k;index.prefix_bits=std::min(c.triangle_global_prefix_bits,index.key_bits);index.prefix_shift=index.key_bits-index.prefix_bits;const auto prefix_count=std::size_t(1)<<index.prefix_bits;index.offsets.assign(prefix_count+1,0);
  std::uint64_t total=0;for(const auto&rows:local){if(rows.size()>std::numeric_limits<std::uint64_t>::max()-total)throw std::runtime_error("global seed occurrence count overflow");total+=rows.size();for(const auto&row:rows){const auto prefix=index.prefix_shift?row.seed>>index.prefix_shift:row.seed;++index.offsets[std::size_t(prefix)+1];}}
  if(total>std::uint64_t(std::numeric_limits<std::uint32_t>::max()))throw std::runtime_error("global seed index requires more than 2^32-1 retained occurrences");
  std::uint64_t running=0;
  for(std::size_t i=1;i<index.offsets.size();++i){running+=index.offsets[i];index.offsets[i]=std::uint32_t(running);}
  if(running!=total)throw std::runtime_error("global seed index prefix count mismatch");
  index.occurrences.resize(std::size_t(total));
  for(auto&rows:local){for(const auto&row:rows){const auto prefix=index.prefix_shift?row.seed>>index.prefix_shift:row.seed;const auto destination=index.offsets[prefix]++;index.occurrences[destination]=row;}std::vector<GlobalSeedOcc>().swap(rows);}std::vector<std::vector<GlobalSeedOcc>>().swap(local);
  for(std::size_t i=prefix_count;i>0;--i)index.offsets[i]=index.offsets[i-1];
  index.offsets[0]=0;
  // Sorting independent high-prefix buckets yields exact
  // (seed,reference,contig,pos) order without one monolithic sort.
  std::atomic<std::size_t>bucket_next{0};error=nullptr;pool.clear();constexpr std::size_t bucket_chunk=4096;
  auto sort_buckets=[&](){for(;;){const auto first=bucket_next.fetch_add(bucket_chunk);if(first>=prefix_count)return;const auto last=std::min(prefix_count,first+bucket_chunk);try{for(auto prefix=first;prefix<last;++prefix){auto begin=index.occurrences.begin()+std::ptrdiff_t(index.offsets[prefix]),end=index.occurrences.begin()+std::ptrdiff_t(index.offsets[prefix+1]);if(end-begin>1)std::sort(begin,end,[](const GlobalSeedOcc&a,const GlobalSeedOcc&b){return std::tie(a.seed,a.reference,a.contig,a.pos)<std::tie(b.seed,b.reference,b.contig,b.pos);});}}catch(...){std::lock_guard<std::mutex>guard(error_mutex);if(!error)error=std::current_exception();return;}}};
  for(std::size_t i=0;i<worker_count;++i)pool.emplace_back(sort_buckets);
  for(auto&t:pool)t.join();
  if(error)std::rethrow_exception(error);
  return index;
}
using ReverseContigs=std::vector<std::string>;
static ReverseContigs make_reverse_contigs(const Genome&g){ReverseContigs out;out.reserve(g.contigs.size());for(const auto&c:g.contigs)out.push_back(reverse_complement(c.seq));return out;}
struct VoteSpan { int count=0; std::uint32_t first=std::numeric_limits<std::uint32_t>::max(),last=0; };
struct VoteKey {
  std::uint32_t qi=0,ri=0;
  int diag=0;
  bool reverse=false;
  bool operator==(const VoteKey&o)const{return qi==o.qi&&ri==o.ri&&diag==o.diag&&reverse==o.reverse;}
};
struct VoteKeyHash {
  std::size_t operator()(const VoteKey&x)const noexcept {
    std::uint64_t h=(std::uint64_t(x.qi)<<32)|x.ri;
    h^=std::uint64_t(std::uint32_t(x.diag))*0x9e3779b97f4a7c15ULL;
    h^=std::uint64_t(x.reverse)*0xbf58476d1ce4e5b9ULL;
    return std::size_t(mix64(h));
  }
};
struct VoteEntry {VoteKey first;VoteSpan second;};
class VoteMap {
  static constexpr std::uint32_t empty_slot=std::numeric_limits<std::uint32_t>::max();
  std::vector<VoteEntry>entries_;
  std::vector<std::uint32_t>buckets_;
  void rehash(std::size_t requested){std::size_t capacity=16;while(capacity<requested)capacity*=2;std::vector<std::uint32_t>next(capacity,empty_slot);const auto mask=capacity-1;VoteKeyHash hash;for(std::uint32_t i=0;i<entries_.size();++i){auto slot=hash(entries_[i].first)&mask;while(next[slot]!=empty_slot)slot=(slot+1)&mask;next[slot]=i;}buckets_.swap(next);}
public:
  VoteMap()=default;
  explicit VoteMap(std::size_t expected_entries){
    if(!expected_entries)return;
    entries_.reserve(expected_entries);
    rehash((expected_entries*10+6)/7);
  }
  VoteSpan&operator[](const VoteKey&key){if(buckets_.empty())rehash(16);if((entries_.size()+1)*10>buckets_.size()*7)rehash(buckets_.size()*2);const auto mask=buckets_.size()-1;VoteKeyHash hash;auto slot=hash(key)&mask;for(;;){const auto index=buckets_[slot];if(index==empty_slot){if(entries_.size()>=std::size_t(empty_slot))throw std::runtime_error("vote table exceeds uint32 index space");entries_.push_back({key,{}});buckets_[slot]=std::uint32_t(entries_.size()-1);return entries_.back().second;}if(entries_[index].first==key)return entries_[index].second;slot=(slot+1)&mask;}}
  auto begin()const{return entries_.begin();}
  auto end()const{return entries_.end();}
};

static Candidate gapped_candidate(const std::string&q,const std::string&r,std::uint32_t q0,std::uint32_t q1,int expected_diag,int band){
  struct Workspace {std::vector<int>prev,cur;std::vector<unsigned char>trace;};
  static thread_local Workspace ws;
  Candidate out;out.q0=q0;out.q1=q1;out.diag=expected_diag;if(q1>q.size()||q0>q1||r.size()>size_t(std::numeric_limits<int>::max()))throw std::runtime_error("alignment coordinates out of range");auto rb64=std::max<std::int64_t>(0,std::int64_t(q0)+expected_diag-std::int64_t(band)),re64=std::min<std::int64_t>(std::int64_t(r.size()),std::int64_t(q1)+expected_diag+std::int64_t(band));if(re64<=rb64||q1<=q0)return out;size_t rb=size_t(rb64),n=q1-q0,m=size_t(re64-rb64);if(n>size_t(std::numeric_limits<int>::max()/4)||m==std::numeric_limits<size_t>::max()||n==std::numeric_limits<size_t>::max()||n+1>std::numeric_limits<size_t>::max()/(m+1))throw std::runtime_error("alignment workspace size overflow");size_t stride=m+1,cells=(n+1)*stride;if(cells>std::numeric_limits<size_t>::max()-3)throw std::runtime_error("alignment traceback size overflow");constexpr int gap=-4;
  ws.prev.assign(stride,0);ws.cur.resize(stride);ws.trace.resize(cells);
  auto get_trace=[&](size_t cell){return ws.trace[cell];};
  auto*prev=ws.prev.data();auto*cur=ws.cur.data();auto*trace=ws.trace.data();const auto*rbase=r.data()+rb;
  for(size_t i=1;i<=n;++i){const auto row=i*stride;const char qb=q[q0+i-1];cur[0]=int(i)*gap;trace[row]=1;for(size_t j=1;j<=m;++j){const char rbv=rbase[j-1];int diag=prev[j-1]+((qb==rbv&&qb!='N')?2:-3),up=prev[j]+gap,left=cur[j-1]+gap;int best=diag;unsigned char move=0;if(up>best){best=up;move=1;}if(left>best){best=left;move=2;}cur[j]=best;trace[row+j]=move;}std::swap(prev,cur);}
  out.aligned.reserve(n);
  size_t i=n,j=0;for(size_t x=1;x<=m;++x)if(prev[x]>prev[j])j=x;while(i){auto move=get_trace(i*stride+j);if(move==0){if(!j)throw std::runtime_error("gapped traceback underflow");auto qp=std::uint32_t(q0+i-1),rp=std::uint32_t(rb+j-1);out.aligned.push_back({qp,rp});--i;--j;}else if(move==1){--i;}else{if(!j)throw std::runtime_error("gapped traceback underflow");--j;}}std::reverse(out.aligned.begin(),out.aligned.end());for(auto[qp,rp]:out.aligned)if(callable(q[qp])&&callable(r[rp])){++out.compared;if(q[qp]==r[rp])++out.matches;}out.identity=out.compared?100.0*out.matches/out.compared:0;return out;
}
template<class VoteIterator>
static Estimate estimate_from_vote_range(const Genome&q,const Genome&r,const ReverseContigs&q_reverse,const Config&c,VoteIterator vote_begin,VoteIterator vote_end){
  if(q_reverse.size()!=q.contigs.size())throw std::runtime_error("reverse-contig cache does not match query genome");
  if(vote_begin==vote_end)return {};
  struct Potential {VoteKey key;VoteSpan span;int lo=0,hi=0;};
  std::vector<Potential>potential;potential.reserve(std::size_t(vote_end-vote_begin));
  std::vector<std::vector<std::pair<std::uint32_t,std::uint32_t>>>q_intervals(q.contigs.size()),r_intervals(r.contigs.size());
  for(auto vote=vote_begin;vote!=vote_end;++vote){const auto&key=vote->first;const auto&v=vote->second;
    if(v.count<c.min_seed_hits)continue;
    const auto qi=key.qi,ri=key.ri;const bool rev=key.reverse;const int diag=key.diag;
    const auto&qseq=rev?q_reverse[qi]:q.contigs[qi].seq;const auto&rseq=r.contigs[ri].seq;
    const auto overlap_lo64=std::max<std::int64_t>(0,-std::int64_t(diag));
    const auto overlap_hi64=std::min<std::int64_t>(std::int64_t(qseq.size()),std::int64_t(rseq.size())-diag);
    const int lo=int(std::max(overlap_lo64,std::int64_t(v.first)-c.align_band));
    const int hi=int(std::min(overlap_hi64,std::int64_t(v.last)+c.k+c.align_band));
    if(lo>=hi)continue;
    potential.push_back({key,v,lo,hi});
    if(rev)q_intervals[qi].push_back({std::uint32_t(qseq.size()-std::size_t(hi)),std::uint32_t(qseq.size()-std::size_t(lo))});
    else q_intervals[qi].push_back({std::uint32_t(lo),std::uint32_t(hi)});
    const auto rb=std::uint32_t(std::max<std::int64_t>(0,std::int64_t(lo)+diag-c.align_band));
    const auto re=std::uint32_t(std::min<std::int64_t>(std::int64_t(rseq.size()),std::int64_t(hi)+diag+c.align_band));
    if(rb<re)r_intervals[ri].push_back({rb,re});
  }
  auto union_bases=[](auto&groups){std::uint64_t total=0;for(auto&intervals:groups){if(intervals.empty())continue;std::sort(intervals.begin(),intervals.end());auto begin=intervals.front().first,end=intervals.front().second;for(std::size_t i=1;i<intervals.size();++i){if(intervals[i].first>end){total+=end-begin;begin=intervals[i].first;end=intervals[i].second;}else end=std::max(end,intervals[i].second);}total+=end-begin;}return total;};
  if(c.min_af>0){const auto possible_aligned=std::min(union_bases(q_intervals),union_bases(r_intervals));const auto denominator=std::min(q.callable,r.callable);const auto possible_af=denominator?100.0L*static_cast<long double>(possible_aligned)/static_cast<long double>(denominator):0.0L;if(possible_af+1e-12L<static_cast<long double>(c.min_af))return {};}
  std::vector<Candidate>cand;for(const auto&p:potential){const auto&key=p.key;const auto&v=p.span;auto qi=key.qi,ri=key.ri;const bool rev=key.reverse;const int diag=key.diag;const auto&qseq=rev?q_reverse[qi]:q.contigs[qi].seq;const auto&rseq=r.contigs[ri].seq;for(int s=p.lo;s<p.hi;){int end=int(std::min<std::int64_t>(p.hi,std::int64_t(s)+c.window));auto x=gapped_candidate(qseq,rseq,std::uint32_t(s),std::uint32_t(end),diag,c.align_band);x.qc=qi;x.rc=ri;x.reverse=rev;x.seed_votes=v.count;if(x.compared>=c.min_window&&x.identity>=c.min_identity)cand.push_back(std::move(x));s=end;}}
  std::sort(cand.begin(),cand.end(),[](const Candidate&a,const Candidate&b){if(a.identity!=b.identity)return a.identity>b.identity;if(a.compared!=b.compared)return a.compared>b.compared;if(a.seed_votes!=b.seed_votes)return a.seed_votes>b.seed_votes;return std::tie(a.qc,a.rc,a.reverse,a.q0,a.diag)<std::tie(b.qc,b.rc,b.reverse,b.q0,b.diag);});
  if(cand.empty())return {};
  std::vector<std::vector<unsigned char>>qcov(q.contigs.size()),rcov(r.contigs.size());std::uint64_t aligned=0,matches=0;
  for(const auto&x:cand){const auto&qseq=x.reverse?q_reverse[x.qc]:q.contigs[x.qc].seq;const auto&rseq=r.contigs[x.rc].seq;auto&qmask=qcov[x.qc];auto&rmask=rcov[x.rc];if(qmask.empty())qmask.assign(q.contigs[x.qc].seq.size(),0);if(rmask.empty())rmask.assign(rseq.size(),0);std::size_t fresh=0;int fresh_match=0;for(auto[qp,rp]:x.aligned){if(!callable(qseq[qp])||!callable(rseq[rp]))continue;size_t oq=x.reverse?qseq.size()-1-qp:qp;if(!qmask[oq]&&!rmask[rp]){++fresh;if(qseq[qp]==rseq[rp])++fresh_match;}}if(int(fresh)<c.min_window/2||int(fresh)*2<x.compared)continue;for(auto[qp,rp]:x.aligned){if(!callable(qseq[qp])||!callable(rseq[rp]))continue;size_t oq=x.reverse?qseq.size()-1-qp:qp;if(!qmask[oq]&&!rmask[rp]){qmask[oq]=1;rmask[rp]=1;}}aligned+=fresh;matches+=fresh_match;}
  Estimate e;e.aligned=aligned;e.matches=matches;e.ani=aligned?100.0*matches/aligned:0;e.af_query=100.0*aligned/q.callable;e.af_ref=100.0*aligned/r.callable;return e;
}
static Estimate estimate_from_votes(const Genome&q,const Genome&r,const ReverseContigs&q_reverse,const Config&c,const VoteMap&votes){return estimate_from_vote_range(q,r,q_reverse,c,votes.begin(),votes.end());}
static Estimate estimate_with_index(const Genome&q,const Genome&r,const SeedIndex&idx,const ReverseContigs&q_reverse,const Config&c){
  if(q_reverse.size()!=q.contigs.size())throw std::runtime_error("reverse-contig cache does not match query genome");
  VoteMap votes;
  for(size_t qi=0;qi<q.contigs.size();++qi){
    const auto&seq=q.contigs[qi].seq;
    each_kmer(seq,c.k,[&](std::uint32_t qp,std::uint64_t fw,std::uint64_t rv){
      auto add_votes=[&](std::uint64_t seed,std::uint32_t query_pos,bool reverse){
        auto it=idx.find(seed);if(it==idx.end())return;
        for(const auto&o:it->second){auto&v=votes[{std::uint32_t(qi),o.contig,int(o.pos)-int(query_pos),reverse}];++v.count;v.first=std::min(v.first,query_pos);v.last=std::max(v.last,query_pos);}
      };
      add_votes(fw,qp,false);
      const auto reverse_pos=std::uint32_t(seq.size()-std::size_t(c.k)-std::size_t(qp));
      add_votes(rv,reverse_pos,true);
    });
  }
  return estimate_from_votes(q,r,q_reverse,c,votes);
}
static Estimate estimate(const Genome&q,const Genome&r,const Config&c){auto idx=reference_seeds(r,c);auto q_reverse=make_reverse_contigs(q);return estimate_with_index(q,r,idx,q_reverse,c);}

static void put_u32(std::ostream&o,std::uint32_t x){char b[4];for(int i=0;i<4;++i)b[i]=char(x>>(8*i));o.write(b,4);}
static void put_u64(std::ostream&o,std::uint64_t x){char b[8];for(int i=0;i<8;++i)b[i]=char(x>>(8*i));o.write(b,8);}
static std::uint32_t get_u32(std::istream&i){unsigned char b[4];if(!i.read(reinterpret_cast<char*>(b),4))throw std::runtime_error("truncated binary index");return std::uint32_t(b[0])|(std::uint32_t(b[1])<<8)|(std::uint32_t(b[2])<<16)|(std::uint32_t(b[3])<<24);}
static std::uint64_t get_u64(std::istream&i){unsigned char b[8];if(!i.read(reinterpret_cast<char*>(b),8))throw std::runtime_error("truncated binary index");std::uint64_t x=0;for(int n=0;n<8;++n)x|=std::uint64_t(b[n])<<(8*n);return x;}
static void magic(std::ostream&o,const char*s){o.write(s,8);}

static std::uint32_t read_le32(const unsigned char*p){return std::uint32_t(p[0])|(std::uint32_t(p[1])<<8)|(std::uint32_t(p[2])<<16)|(std::uint32_t(p[3])<<24);}
static std::uint64_t read_le64(const unsigned char*p){std::uint64_t x=0;for(int n=0;n<8;++n)x|=std::uint64_t(p[n])<<(8*n);return x;}

class MappedFile {
  const unsigned char*data_=nullptr;size_t size_=0;
#if defined(_WIN32)
  HANDLE file_=INVALID_HANDLE_VALUE,map_=nullptr;
  BY_HANDLE_FILE_INFORMATION identity_{};
#else
  int fd_=-1;
  struct stat identity_{};
#endif
public:
  explicit MappedFile(const fs::path&p){
#if defined(_WIN32)
    file_=CreateFileW(p.c_str(),GENERIC_READ,FILE_SHARE_READ|FILE_SHARE_DELETE,nullptr,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);if(file_==INVALID_HANDLE_VALUE)throw std::runtime_error("cannot open mapped index: "+p.string());LARGE_INTEGER n;if(!GetFileSizeEx(file_,&n)||!GetFileInformationByHandle(file_,&identity_)||n.QuadPart<=0||std::uint64_t(n.QuadPart)>std::uint64_t(std::numeric_limits<size_t>::max())){CloseHandle(file_);file_=INVALID_HANDLE_VALUE;throw std::runtime_error("invalid/too-large mapped index: "+p.string());}size_=size_t(n.QuadPart);map_=CreateFileMappingW(file_,nullptr,PAGE_READONLY,0,0,nullptr);if(!map_){CloseHandle(file_);file_=INVALID_HANDLE_VALUE;throw std::runtime_error("cannot create index mapping: "+p.string());}data_=static_cast<const unsigned char*>(MapViewOfFile(map_,FILE_MAP_READ,0,0,0));if(!data_){CloseHandle(map_);CloseHandle(file_);map_=nullptr;file_=INVALID_HANDLE_VALUE;throw std::runtime_error("cannot map index: "+p.string());}
#else
    fd_=::open(p.c_str(),O_RDONLY);if(fd_<0)throw std::runtime_error("cannot open mapped index: "+p.string());if(::flock(fd_,LOCK_SH|LOCK_NB)!=0||::fstat(fd_,&identity_)!=0||identity_.st_size<=0||std::uint64_t(identity_.st_size)>std::uint64_t(std::numeric_limits<size_t>::max())){::close(fd_);fd_=-1;throw std::runtime_error("cannot lock or invalid/too-large mapped index: "+p.string());}size_=size_t(identity_.st_size);void*m=::mmap(nullptr,size_,PROT_READ,MAP_PRIVATE,fd_,0);if(m==MAP_FAILED){::flock(fd_,LOCK_UN);::close(fd_);fd_=-1;throw std::runtime_error("cannot map index: "+p.string());}data_=static_cast<const unsigned char*>(m);
#if defined(MADV_SEQUENTIAL)
    (void)::madvise(const_cast<unsigned char*>(data_),size_,MADV_SEQUENTIAL);
#endif
#endif
  }
  ~MappedFile(){
#if defined(_WIN32)
    if(data_)UnmapViewOfFile(data_);
    if(map_)CloseHandle(map_);
    if(file_!=INVALID_HANDLE_VALUE)CloseHandle(file_);
#else
    if(data_)::munmap(const_cast<unsigned char*>(data_),size_);
    if(fd_>=0){(void)::flock(fd_,LOCK_UN);::close(fd_);}
#endif
  }
  MappedFile(const MappedFile&)=delete;MappedFile&operator=(const MappedFile&)=delete;
  const unsigned char*data()const{return data_;}size_t size()const{return size_;}
  void verify_identity()const{
#if defined(_WIN32)
    BY_HANDLE_FILE_INFORMATION now{};if(file_==INVALID_HANDLE_VALUE||!GetFileInformationByHandle(file_,&now)||now.dwVolumeSerialNumber!=identity_.dwVolumeSerialNumber||now.nFileIndexHigh!=identity_.nFileIndexHigh||now.nFileIndexLow!=identity_.nFileIndexLow||now.nFileSizeHigh!=identity_.nFileSizeHigh||now.nFileSizeLow!=identity_.nFileSizeLow||CompareFileTime(&now.ftLastWriteTime,&identity_.ftLastWriteTime)!=0)throw std::runtime_error("mapped index identity changed");
#else
    struct stat now{};if(fd_<0||::fstat(fd_,&now)!=0||now.st_dev!=identity_.st_dev||now.st_ino!=identity_.st_ino||now.st_size!=identity_.st_size
#if defined(__APPLE__)
      ||now.st_mtimespec.tv_sec!=identity_.st_mtimespec.tv_sec||now.st_mtimespec.tv_nsec!=identity_.st_mtimespec.tv_nsec||now.st_ctimespec.tv_sec!=identity_.st_ctimespec.tv_sec||now.st_ctimespec.tv_nsec!=identity_.st_ctimespec.tv_nsec
#else
      ||now.st_mtim.tv_sec!=identity_.st_mtim.tv_sec||now.st_mtim.tv_nsec!=identity_.st_mtim.tv_nsec||now.st_ctim.tv_sec!=identity_.st_ctim.tv_sec||now.st_ctim.tv_nsec!=identity_.st_ctim.tv_nsec
#endif
    )throw std::runtime_error("mapped index identity changed");
#endif
  }
  void release_pages()const{
#if !defined(_WIN32) && defined(MADV_DONTNEED)
    if(data_)(void)::madvise(const_cast<unsigned char*>(data_),size_,MADV_DONTNEED);
#endif
  }
  void random_access()const{
#if !defined(_WIN32) && defined(MADV_RANDOM)
    if(data_)(void)::madvise(const_cast<unsigned char*>(data_),size_,MADV_RANDOM);
#endif
  }
};

static std::string sha_bytes(const unsigned char*p,size_t n){Sha256 h;constexpr size_t block=1u<<20;while(n){size_t take=std::min(n,block);h.update(p,take);p+=take;n-=take;}return h.finish();}
static void require_bytes(size_t size,std::uint64_t off,std::uint64_t need,const std::string&what){if(off>size||need>std::uint64_t(size)-off)throw std::runtime_error("truncated "+what);}

struct PostingSpan {
  const unsigned char*refs=nullptr;std::uint32_t count=0;
  std::uint32_t at(std::uint32_t i)const{if(i>=count)throw std::runtime_error("posting span access out of range");return read_le32(refs+std::uint64_t(i)*4);}
};
struct PostingIndex {
  enum class Format {LegacyV3,CompactV4};Format format=Format::LegacyV3;
  std::shared_ptr<MappedFile>postings,lookup;
  std::vector<std::uint64_t>legacy_hashes,legacy_count_offsets;
  std::uint64_t groups=0,total_refs=0;
  void verify_snapshot()const{if(!postings)throw std::runtime_error("missing postings snapshot");postings->verify_identity();if(lookup)lookup->verify_identity();}
  PostingSpan find(std::uint64_t h)const{
    if(format==Format::LegacyV3){auto it=std::lower_bound(legacy_hashes.begin(),legacy_hashes.end(),h);if(it==legacy_hashes.end()||*it!=h)return{};size_t i=size_t(it-legacy_hashes.begin());auto count_off=legacy_count_offsets[i],begin=count_off+4,end=i+1<legacy_count_offsets.size()?legacy_count_offsets[i+1]-8:postings->size();if(end<begin||(end-begin)%4||(end-begin)/4>std::uint64_t(std::numeric_limits<std::uint32_t>::max()))throw std::runtime_error("invalid legacy posting span");return{postings->data()+begin,std::uint32_t((end-begin)/4)};}
    std::uint64_t lo=0,hi=groups;const auto*p=lookup->data()+24;while(lo<hi){auto mid=lo+(hi-lo)/2,found=read_le64(p+mid*16);if(found<h)lo=mid+1;else hi=mid;}if(lo==groups||read_le64(p+lo*16)!=h)return{};auto begin=read_le64(p+lo*16+8),end=lo+1<groups?read_le64(p+(lo+1)*16+8):total_refs;if(begin>total_refs||end>total_refs||end<begin||end-begin>std::uint64_t(std::numeric_limits<std::uint32_t>::max()))throw std::runtime_error("invalid compact posting span");return{postings->data()+16+begin*4,std::uint32_t(end-begin)};
  }
};

static std::vector<fs::path> read_paths(const fs::path&p,bool manifest){std::ifstream f(p);if(!f)throw std::runtime_error("cannot read path list: "+p.string());std::vector<fs::path>v;std::set<std::string>seen;std::string line;bool first=true;while(std::getline(f,line)){line=trim(line);if(line.empty()||line[0]=='#')continue;auto x=split(line,'\t');if(first&&manifest&&(x[0]=="Ref_file"||x[0]=="path"||x[0]=="reference")){first=false;continue;}first=false;fs::path q=x[0];if(q.is_relative())q=p.parent_path()/q;if(!fs::is_regular_file(q))throw std::runtime_error("listed FASTA missing: "+q.string());q=fs::weakly_canonical(q);if(!seen.insert(q.string()).second)throw std::runtime_error("duplicate reference/query path: "+q.string());v.push_back(q);}if(v.empty())throw std::runtime_error("empty path list: "+p.string());return v;}
struct Built { RefRow row; std::vector<std::uint64_t> sketch; };
static Built build_one(const fs::path&p,const Config&c){auto g=read_fasta(p);auto s=make_scaled_sketch(g,c.k,c.sketch_scale,0);if(!std::is_sorted(s.begin(),s.end())||std::adjacent_find(s.begin(),s.end())!=s.end())throw std::runtime_error("internal sketch ordering failure");if(s.size()>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("reference sketch exceeds index format");RefRow r;r.path=g.path;r.bytes=g.source_bytes;r.callable=g.callable;r.mtime=g.source_mtime;r.content_sha256=std::move(g.content_sha256);r.sketch_count=std::uint32_t(s.size());return{std::move(r),std::move(s)};}

static void flush_pairs(std::vector<std::pair<std::uint64_t,std::uint32_t>>&pairs,const fs::path&dir,std::vector<fs::path>&chunks){if(pairs.empty())return;std::sort(pairs.begin(),pairs.end());auto p=dir/("postings_chunk_"+std::to_string(chunks.size())+".bin");std::array<char,1u<<20>buffer{};std::ofstream o;o.rdbuf()->pubsetbuf(buffer.data(),std::streamsize(buffer.size()));o.open(p,std::ios::binary);if(!o)throw std::runtime_error("cannot create postings chunk: "+p.string());for(auto[h,r]:pairs){put_u64(o,h);put_u32(o,r);}o.close();if(!o||fs::file_size(p)!=pairs.size()*12ULL)throw std::runtime_error("failed writing postings chunk: "+p.string());chunks.push_back(p);pairs.clear();}
struct ChunkReader {
  std::array<char,1u<<16>buffer{};std::ifstream in;std::uint64_t h=0;std::uint32_t r=0;bool valid=false;
  explicit ChunkReader(const fs::path&p){in.rdbuf()->pubsetbuf(buffer.data(),std::streamsize(buffer.size()));in.open(p,std::ios::binary);if(!in)throw std::runtime_error("cannot open postings chunk: "+p.string());next();}
  void next(){if(in.peek()==EOF){if(!in.eof())throw std::runtime_error("failed reading postings chunk");valid=false;return;}h=get_u64(in);r=get_u32(in);valid=true;}
};
using ChunkQueueItem=std::tuple<std::uint64_t,std::uint32_t,size_t>;
static void merge_pair_chunk_group(const std::vector<fs::path>&inputs,const fs::path&out){
  std::vector<ChunkReader>readers;readers.reserve(inputs.size());for(const auto&p:inputs)readers.emplace_back(p);
  std::priority_queue<ChunkQueueItem,std::vector<ChunkQueueItem>,std::greater<ChunkQueueItem>>q;
  for(size_t i=0;i<readers.size();++i)if(readers[i].valid)q.push({readers[i].h,readers[i].r,i});
  std::array<char,1u<<20>buffer{};std::ofstream o;o.rdbuf()->pubsetbuf(buffer.data(),std::streamsize(buffer.size()));o.open(out,std::ios::binary);if(!o)throw std::runtime_error("cannot create merged postings chunk: "+out.string());
  bool have=false;std::uint64_t last_h=0,records=0;std::uint32_t last_r=0;
  while(!q.empty()){auto[h,r,i]=q.top();q.pop();if(!have||h!=last_h||r!=last_r){put_u64(o,h);put_u32(o,r);last_h=h;last_r=r;have=true;++records;}readers[i].next();if(readers[i].valid)q.push({readers[i].h,readers[i].r,i});}
  o.close();if(!o||fs::file_size(out)!=records*12ULL)throw std::runtime_error("failed writing merged postings chunk: "+out.string());
  readers.clear();for(const auto&p:inputs)if(!fs::remove(p))throw std::runtime_error("failed removing verified temporary chunk: "+p.string());
}
struct MergeResult {std::vector<std::uint32_t>counts;std::uint64_t groups=0,total_refs=0;};
static MergeResult merge_chunks(std::vector<fs::path>chunks,const fs::path&out,const fs::path&lookup,size_t nrefs,int max_posting,bool compact){
  // GTDB contains enough references to create thousands of sorted chunks.
  // Merge them in bounded fan-in rounds so this never depends on a raised
  // process file-descriptor limit.
  constexpr size_t fan_in=128;size_t round=0;
  while(chunks.size()>fan_in){std::vector<fs::path>next;for(size_t base=0;base<chunks.size();base+=fan_in){size_t end=std::min(chunks.size(),base+fan_in);std::vector<fs::path>group(chunks.begin()+base,chunks.begin()+end);auto p=out.parent_path()/("postings_merge_r"+std::to_string(round)+"_"+std::to_string(next.size())+".bin");merge_pair_chunk_group(group,p);next.push_back(p);}chunks=std::move(next);++round;}
  std::vector<ChunkReader>readers;readers.reserve(chunks.size());for(const auto&p:chunks)readers.emplace_back(p);std::priority_queue<ChunkQueueItem,std::vector<ChunkQueueItem>,std::greater<ChunkQueueItem>>q;for(size_t i=0;i<readers.size();++i)if(readers[i].valid)q.push({readers[i].h,readers[i].r,i});
  std::array<char,1u<<20>out_buffer{},lookup_buffer{};std::ofstream o;o.rdbuf()->pubsetbuf(out_buffer.data(),std::streamsize(out_buffer.size()));o.open(out,std::ios::binary);if(!o)throw std::runtime_error("cannot create postings index: "+out.string());std::ofstream lu;if(compact){lu.rdbuf()->pubsetbuf(lookup_buffer.data(),std::streamsize(lookup_buffer.size()));lu.open(lookup,std::ios::binary);if(!lu)throw std::runtime_error("cannot create postings lookup: "+lookup.string());magic(o,"GAAFRF4\0");put_u64(o,0);magic(lu,"GAAFLU4\0");put_u64(lu,0);put_u64(lu,0);}else{magic(o,"GAAFPO3\0");put_u64(o,0);}
  MergeResult result;result.counts.assign(nrefs,0);while(!q.empty()){auto h=std::get<0>(q.top());std::vector<std::uint32_t>refs;while(!q.empty()&&std::get<0>(q.top())==h){auto[hh,r,i]=q.top();q.pop();if(refs.empty()||refs.back()!=r)refs.push_back(r);readers[i].next();if(readers[i].valid)q.push({readers[i].h,readers[i].r,i});}if(refs.size()>size_t(max_posting))continue;if(compact){put_u64(lu,h);put_u64(lu,result.total_refs);}else{put_u64(o,h);put_u32(o,std::uint32_t(refs.size()));}for(auto r:refs){put_u32(o,r);++result.counts[r];++result.total_refs;}++result.groups;}
  o.seekp(8);put_u64(o,compact?result.total_refs:result.groups);if(compact){lu.seekp(8);put_u64(lu,result.groups);put_u64(lu,result.total_refs);}o.close();if(compact)lu.close();if(!o||(compact&&!lu))throw std::runtime_error("failed writing postings index");auto expected_out=compact?16ULL+result.total_refs*4ULL:16ULL+result.groups*12ULL+result.total_refs*4ULL;if(fs::file_size(out)!=expected_out)throw std::runtime_error("postings size closure failure");if(compact&&fs::file_size(lookup)!=24ULL+result.groups*16ULL)throw std::runtime_error("postings lookup size closure failure");
  readers.clear();for(const auto&p:chunks)if(!fs::remove(p))throw std::runtime_error("failed removing verified temporary chunk: "+p.string());return result;
}

static fs::path temporary_sibling(const fs::path&out){auto stamp=std::chrono::high_resolution_clock::now().time_since_epoch().count();return out.parent_path()/(out.filename().string()+".tmp."+std::to_string(stamp));}
static std::array<unsigned char,32> raw_digest(const std::string&hex){std::array<unsigned char,32>out{};auto nibble=[](char c)->unsigned char{return static_cast<unsigned char>(c<='9'?c-'0':c-'a'+10);};for(size_t i=0;i<out.size();++i)out[i]=static_cast<unsigned char>((nibble(hex[2*i])<<4)|nibble(hex[2*i+1]));return out;}
static const std::string&sketch_digest_domain(){static const std::string domain="gtdb-ani-af-sketch-v4-le64";return domain;}
static std::array<unsigned char,32> sketch_digest(const std::vector<std::uint64_t>&sketch){Sha256 h;h.update(sketch_digest_domain());for(auto x:sketch){unsigned char b[8];for(int i=0;i<8;++i)b[i]=static_cast<unsigned char>(x>>(8*i));h.update(b,8);}return raw_digest(h.finish());}
static void build_index(const Config&c){
  if(fs::exists(c.out_dir)){throw std::runtime_error("index destination exists: "+c.out_dir.string());}for(const auto&p:{c.manifest,c.taxonomy,c.radii,c.release_metadata}){if(!fs::is_regular_file(p))throw std::runtime_error("required provenance file missing: "+p.string());}auto refs=read_paths(c.manifest,true);if(refs.size()>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("too many references for index format");auto tmp=temporary_sibling(c.out_dir);if(!fs::create_directories(tmp))throw std::runtime_error("temporary index directory already exists: "+tmp.string());try{
    std::array<char,1u<<20>sketch_buffer{};std::ofstream sketches;sketches.rdbuf()->pubsetbuf(sketch_buffer.data(),std::streamsize(sketch_buffer.size()));sketches.open(tmp/"SKETCHES.bin",std::ios::binary);if(!sketches)throw std::runtime_error("cannot create SKETCHES.bin");magic(sketches,c.compact_index?"GAAFSK4\0":"GAAFSK3\0");put_u32(sketches,std::uint32_t(refs.size()));std::vector<RefRow>rows(refs.size());std::vector<std::pair<std::uint64_t,std::uint32_t>>pairs;pairs.reserve(2000000);std::vector<fs::path>chunks;size_t batch=size_t(std::max(1,c.threads));
    for(size_t base=0;base<refs.size();base+=batch){size_t n=std::min(batch,refs.size()-base);std::vector<Built>built(n);std::atomic<size_t>next{0};std::exception_ptr error;std::mutex em;auto worker=[&](){for(;;){size_t i=next.fetch_add(1);if(i>=n)return;try{built[i]=build_one(refs[base+i],c);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();return;}}};std::vector<std::thread>ts;for(int t=0;t<std::min<int>(c.threads,n);++t)ts.emplace_back(worker);for(auto&t:ts)t.join();if(error)std::rethrow_exception(error);
      for(size_t i=0;i<n;++i){auto id=std::uint32_t(base+i);auto&b=built[i];b.row.id=id;if(c.compact_index){put_u32(sketches,b.row.sketch_count);auto digest=sketch_digest(b.sketch);sketches.write(reinterpret_cast<const char*>(digest.data()),std::streamsize(digest.size()));}else{put_u32(sketches,b.row.sketch_count);for(auto h:b.sketch)put_u64(sketches,h);}for(auto h:b.sketch)pairs.push_back({h,id});rows[id]=std::move(b.row);std::vector<std::uint64_t>().swap(b.sketch);if(pairs.size()>=2000000)flush_pairs(pairs,tmp,chunks);}}
    flush_pairs(pairs,tmp,chunks);sketches.close();if(!sketches)throw std::runtime_error("failed writing SKETCHES.bin");auto merged=merge_chunks(chunks,tmp/"POSTINGS.bin",tmp/"POSTINGS_LOOKUP.bin",refs.size(),c.max_posting,c.compact_index);std::ofstream meta(tmp/"REFS.tsv");meta<<"ref_id\tRef_file\tFile_bytes\tMtime_ticks\tContent_SHA256\tCallable_bases\tSketch_count\tPosting_count\n";std::uint64_t kept_postings=0;for(size_t i=0;i<rows.size();++i){auto&r=rows[i];r.posting_count=merged.counts[i];kept_postings+=r.posting_count;meta<<r.id<<'\t'<<r.path.string()<<'\t'<<r.bytes<<'\t'<<r.mtime<<'\t'<<r.content_sha256<<'\t'<<r.callable<<'\t'<<r.sketch_count<<'\t'<<r.posting_count<<'\n';}meta.close();if(!meta||kept_postings!=merged.total_refs)throw std::runtime_error("failed writing REFS.tsv/posting closure");
    std::ostringstream complete;complete<<"{\n  \"schema\": \""<<(c.compact_index?"gtdb-ani-af-index-v4":"gtdb-ani-af-index-v3")<<"\",\n  \"status\": \"PASS\",\n  \"k\": "<<c.k<<",\n  \"sketch_scale\": "<<c.sketch_scale<<",\n  \"query_sketch_size\": "<<c.query_sketch_size<<",\n  \"max_posting\": "<<c.max_posting<<",\n  \"reference_count\": "<<refs.size()<<",\n  \"kept_posting_count\": "<<kept_postings<<",\n";if(c.compact_index)complete<<"  \"posting_group_count\": "<<merged.groups<<",\n";complete<<"  \"reference_manifest_sha256\": \""<<sha_file(c.manifest)<<"\",\n  \"taxonomy_sha256\": \""<<sha_file(c.taxonomy)<<"\",\n  \"radii_sha256\": \""<<sha_file(c.radii)<<"\",\n  \"release_metadata_sha256\": \""<<sha_file(c.release_metadata)<<"\",\n  \"refs_sha256\": \""<<sha_file(tmp/"REFS.tsv")<<"\",\n  \"sketches_sha256\": \""<<sha_file(tmp/"SKETCHES.bin")<<"\",\n  \"postings_sha256\": \""<<sha_file(tmp/"POSTINGS.bin")<<"\",\n";if(c.compact_index)complete<<"  \"postings_lookup_sha256\": \""<<sha_file(tmp/"POSTINGS_LOOKUP.bin")<<"\",\n";complete<<"  \"manifest_path\": \""<<json_escape(fs::weakly_canonical(c.manifest).string())<<"\",\n  \"taxonomy_path\": \""<<json_escape(fs::weakly_canonical(c.taxonomy).string())<<"\",\n  \"radii_path\": \""<<json_escape(fs::weakly_canonical(c.radii).string())<<"\",\n  \"release_metadata_path\": \""<<json_escape(fs::weakly_canonical(c.release_metadata).string())<<"\"\n}\n";std::ofstream co(tmp/"COMPLETE.json");co<<complete.str();co.close();if(!co)throw std::runtime_error("failed writing COMPLETE.json");fs::rename(tmp,c.out_dir);
  }catch(...){std::error_code ec;fs::remove_all(tmp,ec);throw;}}

static std::string json_string(const std::string&s,const std::string&key){auto n="\""+key+"\": \"";auto p=s.find(n);if(p==std::string::npos)throw std::runtime_error("COMPLETE missing "+key);p+=n.size();std::string out;for(;p<s.size();++p){char ch=s[p];if(ch=='"')return out;if(ch!='\\'){if(static_cast<unsigned char>(ch)<0x20)throw std::runtime_error("malformed COMPLETE string");out+=ch;continue;}if(++p>=s.size())throw std::runtime_error("malformed COMPLETE escape");switch(s[p]){case'"':out+='"';break;case'\\':out+='\\';break;case'/':out+='/';break;case'b':out+='\b';break;case'f':out+='\f';break;case'n':out+='\n';break;case'r':out+='\r';break;case't':out+='\t';break;default:throw std::runtime_error("unsupported COMPLETE escape");}}throw std::runtime_error("malformed COMPLETE");}
static std::uint64_t json_uint(const std::string&s,const std::string&key){auto n="\""+key+"\": ";auto p=s.find(n);if(p==std::string::npos)throw std::runtime_error("COMPLETE missing "+key);p+=n.size();size_t e=p;while(e<s.size()&&std::isdigit(static_cast<unsigned char>(s[e])))++e;if(e==p)throw std::runtime_error("malformed COMPLETE integer");size_t tail=e;while(tail<s.size()&&std::isspace(static_cast<unsigned char>(s[tail])))++tail;if(tail>=s.size()||(s[tail]!=','&&s[tail]!='}'))throw std::runtime_error("malformed COMPLETE integer suffix");try{return std::stoull(s.substr(p,e-p));}catch(const std::out_of_range&){throw std::runtime_error("COMPLETE integer overflow: "+key);}}
static int json_int_range(const std::string&s,const std::string&key,std::uint64_t lo,std::uint64_t hi){auto x=json_uint(s,key);if(x<lo||x>hi)throw std::runtime_error("COMPLETE integer out of range: "+key);return int(x);}
static std::uint64_t parse_u64_exact(const std::string&s,const std::string&what){if(s.empty()||!std::all_of(s.begin(),s.end(),[](unsigned char c){return std::isdigit(c);}))throw std::runtime_error("malformed unsigned integer in "+what);try{size_t used=0;auto x=std::stoull(s,&used);if(used!=s.size())throw std::runtime_error("malformed unsigned integer in "+what);return x;}catch(const std::invalid_argument&){throw std::runtime_error("malformed unsigned integer in "+what);}catch(const std::out_of_range&){throw std::runtime_error("unsigned integer overflow in "+what);}}
static std::uint32_t parse_u32_exact(const std::string&s,const std::string&what){auto x=parse_u64_exact(s,what);if(x>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("32-bit integer overflow in "+what);return std::uint32_t(x);}
static std::int64_t parse_i64_exact(const std::string&s,const std::string&what){size_t first=!s.empty()&&s[0]=='-'?1:0;if(first==s.size()||!std::all_of(s.begin()+std::string::difference_type(first),s.end(),[](unsigned char c){return std::isdigit(c);}))throw std::runtime_error("malformed signed integer in "+what);try{size_t used=0;auto x=std::stoll(s,&used);if(used!=s.size())throw std::runtime_error("malformed signed integer in "+what);return std::int64_t(x);}catch(const std::invalid_argument&){throw std::runtime_error("malformed signed integer in "+what);}catch(const std::out_of_range&){throw std::runtime_error("signed integer overflow in "+what);}}
static bool valid_sha256(const std::string&s){return s.size()==64&&std::all_of(s.begin(),s.end(),[](unsigned char c){return std::isdigit(c)||(c>='a'&&c<='f');});}

// Taxonomy is provenance for this estimator; it is deliberately absent from
// the numerical search/triangle kernels.  A taxonomy rebind must nevertheless
// prove that it describes exactly the accessions frozen in REFS.tsv.  Keep the
// parser here (rather than trusting line counts) so a duplicated, missing, or
// release-mismatched accession fails closed before a new index is published.
static std::string lower_ascii(std::string s){for(char&c:s)c=char(std::tolower(static_cast<unsigned char>(c)));return s;}
static std::string strip_utf8_bom(std::string s){if(s.size()>=3&&static_cast<unsigned char>(s[0])==0xef&&static_cast<unsigned char>(s[1])==0xbb&&static_cast<unsigned char>(s[2])==0xbf)s.erase(0,3);return s;}
static std::string taxonomy_accession(std::string s){
  s=trim(strip_utf8_bom(std::move(s)));if(s.size()>=2&&s.front()=='"'&&s.back()=='"')s=s.substr(1,s.size()-2);s=fs::path(s).filename().string();
  if(s.rfind("RS_",0)==0||s.rfind("GB_",0)==0)s.erase(0,3);
  for(const auto&prefix:{std::string("GCF_"),std::string("GCA_")}){auto p=s.find(prefix);if(p==std::string::npos)continue;auto i=p+prefix.size(),digits=i;while(i<s.size()&&std::isdigit(static_cast<unsigned char>(s[i])))++i;if(i==digits||i>=s.size()||s[i]!='.')continue;++i;auto version=i;while(i<s.size()&&std::isdigit(static_cast<unsigned char>(s[i])))++i;if(i>version)return s.substr(p,i-p);}
  for(const auto&suffix:{"_genomic.fna.gz","_genomic.fna",".fna.gz",".fna",".fa.gz",".fa",".gz"}){auto n=std::strlen(suffix);if(s.size()>=n&&s.compare(s.size()-n,n,suffix)==0){s.resize(s.size()-n);break;}}
  return trim(s);
}
static std::vector<std::string> quoted_fields(const std::string&line,char delimiter){
  std::vector<std::string>out;std::string field;bool quoted=false;
  for(size_t i=0;i<line.size();++i){char ch=line[i];if(ch=='"'){if(quoted&&i+1<line.size()&&line[i+1]=='"'){field+='"';++i;}else quoted=!quoted;}else if(ch==delimiter&&!quoted){out.push_back(field);field.clear();}else field+=ch;}
  if(quoted)throw std::runtime_error("unterminated quote in taxonomy row");
  out.push_back(field);return out;
}
struct TaxonomyAccessions {std::set<std::string> values;std::uint64_t rows=0;std::string set_sha256;};
static std::string accession_set_sha256(const std::set<std::string>&values){Sha256 h;for(const auto&x:values){h.update(x);h.update("\n");}return h.finish();}
static TaxonomyAccessions load_taxonomy_accessions(const fs::path&path){
  std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("cannot read taxonomy: "+path.string());std::string first;while(std::getline(in,first)){if(!first.empty()&&first.back()=='\r')first.pop_back();if(!trim(first).empty())break;}if(trim(first).empty())throw std::runtime_error("empty taxonomy: "+path.string());first=strip_utf8_bom(first);
  const char delimiter=first.find('\t')!=std::string::npos?'\t':',';if(delimiter==','&&first.find(',')==std::string::npos)throw std::runtime_error("taxonomy must be comma- or tab-delimited: "+path.string());auto head=quoted_fields(first,delimiter);for(auto&x:head)x=lower_ascii(trim(x));
  static const std::set<std::string> accession_names={"accession","access_id","genome_id","ncbi_genbank_assembly_accession","ncbi_refseq_assembly_accession"};
  static const std::set<std::string> taxonomy_names={"gtdb_taxonomy","taxonomy","classification"};
  bool has_header=!head.empty()&&accession_names.count(head[0]);size_t accession_column=0,taxonomy_column=1;
  if(has_header){bool have_accession=false,have_taxonomy=false;for(size_t i=0;i<head.size();++i){if(!have_accession&&accession_names.count(head[i])){accession_column=i;have_accession=true;}if(!have_taxonomy&&taxonomy_names.count(head[i])){taxonomy_column=i;have_taxonomy=true;}}if(!have_accession||!have_taxonomy)throw std::runtime_error("taxonomy header lacks accession/taxonomy columns: "+path.string());}
  TaxonomyAccessions result;auto add=[&](const std::string&line,std::uint64_t line_number){if(trim(line).empty())return;auto fields=quoted_fields(strip_utf8_bom(line),delimiter);if(fields.size()<=std::max(accession_column,taxonomy_column))throw std::runtime_error(path.string()+":"+std::to_string(line_number)+" lacks accession/taxonomy fields");auto accession=taxonomy_accession(fields[accession_column]);auto tax=trim(fields[taxonomy_column]);if(accession.empty()||tax.empty())throw std::runtime_error(path.string()+":"+std::to_string(line_number)+" has empty accession/taxonomy");if(!result.values.insert(accession).second)throw std::runtime_error("duplicate taxonomy accession: "+accession);++result.rows;};
  std::uint64_t line_number=1;if(!has_header)add(first,line_number);std::string line;while(std::getline(in,line)){++line_number;if(!line.empty()&&line.back()=='\r')line.pop_back();add(line,line_number);}if(!in.eof())throw std::runtime_error("error reading taxonomy: "+path.string());result.set_sha256=accession_set_sha256(result.values);return result;
}
static std::set<std::string> reference_accessions(const std::vector<RefRow>&refs){std::set<std::string>out;for(const auto&r:refs){auto a=taxonomy_accession(r.path.string());if(a.empty())throw std::runtime_error("cannot derive accession from REFS path: "+r.path.string());if(!out.insert(a).second)throw std::runtime_error("duplicate normalized accession in REFS: "+a);}return out;}
static void require_same_accessions(const std::set<std::string>&refs,const std::set<std::string>&taxonomy){
  if(refs==taxonomy)return;
  std::vector<std::string>missing,extra;std::set_difference(refs.begin(),refs.end(),taxonomy.begin(),taxonomy.end(),std::back_inserter(missing));std::set_difference(taxonomy.begin(),taxonomy.end(),refs.begin(),refs.end(),std::back_inserter(extra));std::ostringstream msg;msg<<"taxonomy/REFS accession sets differ: missing="<<missing.size()<<" extra="<<extra.size();for(size_t i=0;i<std::min<size_t>(missing.size(),3);++i)msg<<" missing["<<i<<"]="<<missing[i];for(size_t i=0;i<std::min<size_t>(extra.size(),3);++i)msg<<" extra["<<i<<"]="<<extra[i];throw std::runtime_error(msg.str());
}
static std::string read_file_bytes(const fs::path&path){std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot read: "+path.string());std::ostringstream out;out<<f.rdbuf();if(f.bad())throw std::runtime_error("error reading: "+path.string());return out.str();}
static void validate_taxonomy_rebind_contract(const fs::path&dir,const std::string&complete,const std::vector<RefRow>&refs){
  if(complete.find("\"taxonomy_rebound_from_schema\"")==std::string::npos)return;
  const auto source_schema=json_string(complete,"taxonomy_rebound_from_schema");if(source_schema!="gtdb-ani-af-index-v3"&&source_schema!="gtdb-ani-af-index-v4")throw std::runtime_error("unsupported taxonomy rebind source schema");
  if(json_string(complete,"source_receipt_file")!="SOURCE_COMPLETE.json"||json_string(complete,"taxonomy_bound_file")!="TAXONOMY.tsv")throw std::runtime_error("taxonomy rebind uses an unsafe/noncanonical bound filename");
  if(json_string(complete,"source_index_path").empty()||!valid_sha256(json_string(complete,"source_complete_sha256"))||!valid_sha256(json_string(complete,"source_taxonomy_sha256"))||!valid_sha256(json_string(complete,"taxonomy_accession_set_sha256")))throw std::runtime_error("invalid taxonomy rebind provenance");
  auto source=read_file_bytes(dir/"SOURCE_COMPLETE.json");if(sha_text(source)!=json_string(complete,"source_complete_sha256"))throw std::runtime_error("embedded source COMPLETE SHA256 mismatch");if(json_string(source,"schema")!=source_schema||json_string(source,"status")!="PASS")throw std::runtime_error("embedded source COMPLETE is not the declared PASS index");
  for(const auto&key:{"k","sketch_scale","query_sketch_size","max_posting","reference_count","kept_posting_count"})if(json_uint(source,key)!=json_uint(complete,key))throw std::runtime_error("taxonomy rebind changed source integer field: "+std::string(key));
  for(const auto&key:{"reference_manifest_sha256","radii_sha256","release_metadata_sha256","refs_sha256","sketches_sha256","postings_sha256"})if(json_string(source,key)!=json_string(complete,key))throw std::runtime_error("taxonomy rebind changed source digest field: "+std::string(key));
  if(json_string(complete,"schema")=="gtdb-ani-af-index-v4"&&json_string(source,"postings_lookup_sha256")!=json_string(complete,"postings_lookup_sha256"))throw std::runtime_error("taxonomy rebind changed source postings lookup digest");
  if(json_string(source,"taxonomy_sha256")!=json_string(complete,"source_taxonomy_sha256"))throw std::runtime_error("taxonomy rebind source taxonomy digest mismatch");
  const auto bound=dir/"TAXONOMY.tsv";if(!fs::is_regular_file(bound)||sha_file(bound)!=json_string(complete,"taxonomy_sha256"))throw std::runtime_error("bound taxonomy SHA256 mismatch");auto taxonomy=load_taxonomy_accessions(bound);auto refset=reference_accessions(refs);require_same_accessions(refset,taxonomy.values);if(taxonomy.rows!=json_uint(complete,"taxonomy_entry_count")||taxonomy.rows!=refs.size()||taxonomy.set_sha256!=json_string(complete,"taxonomy_accession_set_sha256"))throw std::runtime_error("taxonomy rebind accession closure mismatch");
  const auto materialization=json_string(complete,"core_materialization");if(materialization!="hardlink"&&materialization!="copy"&&materialization!="mixed")throw std::runtime_error("unknown taxonomy rebind core materialization");
}
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
struct IndexSourceSnapshot {std::string complete_bytes,refs_bytes;};
static std::pair<std::vector<RefRow>,PostingIndex> load_index(const fs::path&dir,Config&c,std::vector<std::array<unsigned char,32>>*legacy_digests=nullptr,IndexSourceSnapshot*source_snapshot=nullptr){
  for(const auto&n:{"COMPLETE.json","REFS.tsv","SKETCHES.bin","POSTINGS.bin"})if(!fs::is_regular_file(dir/n))throw std::runtime_error("incomplete index: missing "+std::string(n));
  auto complete_map=std::make_shared<MappedFile>(dir/"COMPLETE.json");std::string complete(reinterpret_cast<const char*>(complete_map->data()),complete_map->size());complete_map->verify_identity();auto schema=json_string(complete,"schema");bool compact=schema=="gtdb-ani-af-index-v4";if(!compact&&schema!="gtdb-ani-af-index-v3")throw std::runtime_error("index COMPLETE is not supported v3/v4");if(json_string(complete,"status")!="PASS")throw std::runtime_error("index COMPLETE is not PASS");if(compact&&!fs::is_regular_file(dir/"POSTINGS_LOOKUP.bin"))throw std::runtime_error("incomplete compact index: missing POSTINGS_LOOKUP.bin");
  int ik=json_int_range(complete,"k",5,31),iscale=json_int_range(complete,"sketch_scale",1,std::uint64_t(std::numeric_limits<int>::max())),iqs=json_int_range(complete,"query_sketch_size",10,1000),imp=json_int_range(complete,"max_posting",1,std::uint64_t(std::numeric_limits<int>::max()));if(c.k_explicit&&c.k!=ik)throw std::runtime_error("--k disagrees with index");if(c.sketch_scale_explicit&&c.sketch_scale!=iscale)throw std::runtime_error("--sketch-scale disagrees with index");if(c.query_sketch_explicit&&c.query_sketch_size!=iqs)throw std::runtime_error("--query-sketch-size disagrees with index");if(c.max_posting_explicit&&c.max_posting!=imp)throw std::runtime_error("--max-posting disagrees with index");c.k=ik;c.sketch_scale=iscale;c.query_sketch_size=iqs;c.max_posting=imp;
  for(const auto&key:{"reference_manifest_sha256","taxonomy_sha256","radii_sha256","release_metadata_sha256","refs_sha256","sketches_sha256","postings_sha256"})if(!valid_sha256(json_string(complete,key)))throw std::runtime_error("COMPLETE has invalid SHA256: "+std::string(key));if(compact&&!valid_sha256(json_string(complete,"postings_lookup_sha256")))throw std::runtime_error("COMPLETE has invalid compact lookup SHA256");if(complete.find("\"converted_from_schema\"")!=std::string::npos&&(schema!="gtdb-ani-af-index-v4"||json_string(complete,"converted_from_schema")!="gtdb-ani-af-index-v3"||!valid_sha256(json_string(complete,"source_complete_sha256"))||json_string(complete,"source_index_path").empty()))throw std::runtime_error("invalid compact conversion provenance");auto refs_map=std::make_shared<MappedFile>(dir/"REFS.tsv");std::string refs_text(reinterpret_cast<const char*>(refs_map->data()),refs_map->size());refs_map->verify_identity();if(json_string(complete,"refs_sha256")!=sha_text(refs_text))throw std::runtime_error("index REFS SHA256 mismatch");
  std::istringstream rf(refs_text);std::string line;if(!std::getline(rf,line))throw std::runtime_error("missing REFS header");if(!line.empty()&&line.back()=='\r')line.pop_back();if(line!="ref_id\tRef_file\tFile_bytes\tMtime_ticks\tContent_SHA256\tCallable_bases\tSketch_count\tPosting_count")throw std::runtime_error("wrong REFS header");std::vector<RefRow>refs;std::set<std::string>paths;while(std::getline(rf,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();auto x=split(line,'\t');if(x.size()!=8)throw std::runtime_error("malformed REFS row");RefRow r;r.id=parse_u32_exact(x[0],"REFS ref_id");r.path=x[1];r.bytes=parse_u64_exact(x[2],"REFS File_bytes");r.mtime=parse_i64_exact(x[3],"REFS Mtime_ticks");r.content_sha256=x[4];r.callable=parse_u64_exact(x[5],"REFS Callable_bases");r.sketch_count=parse_u32_exact(x[6],"REFS Sketch_count");r.posting_count=parse_u32_exact(x[7],"REFS Posting_count");if(r.id!=refs.size()||r.path.empty()||!r.bytes||!r.callable||!paths.insert(r.path.string()).second||!valid_sha256(r.content_sha256)||r.posting_count>r.sketch_count)throw std::runtime_error("duplicate/out-of-order/invalid REFS row");refs.push_back(std::move(r));}if(refs.size()!=json_uint(complete,"reference_count"))throw std::runtime_error("reference count disagrees with COMPLETE");
  if(compact){auto sk=std::make_shared<MappedFile>(dir/"SKETCHES.bin");if(sk->size()!=12ULL+refs.size()*36ULL||std::memcmp(sk->data(),"GAAFSK4\0",8)!=0||read_le32(sk->data()+8)!=refs.size())throw std::runtime_error("compact SKETCHES structure/count mismatch");if(json_string(complete,"sketches_sha256")!=sha_bytes(sk->data(),sk->size()))throw std::runtime_error("compact SKETCHES SHA256 mismatch");for(size_t i=0;i<refs.size();++i)if(read_le32(sk->data()+12+i*36)!=refs[i].sketch_count)throw std::runtime_error("compact SKETCHES per-reference count mismatch");sk->verify_identity();sk->release_pages();}
  else{auto sk=std::make_shared<MappedFile>(dir/"SKETCHES.bin");auto*p=sk->data();auto size=sk->size();require_bytes(size,0,12,"legacy SKETCHES");if(std::memcmp(p,"GAAFSK3\0",8)!=0||read_le32(p+8)!=refs.size())throw std::runtime_error("legacy SKETCHES structure/count mismatch");if(legacy_digests){legacy_digests->clear();legacy_digests->reserve(refs.size());}Sha256 whole;whole.update(p,12);std::uint64_t pos=12;for(size_t i=0;i<refs.size();++i){require_bytes(size,pos,4,"legacy SKETCHES count");auto n=read_le32(p+pos);if(n!=refs[i].sketch_count)throw std::runtime_error("SKETCHES per-reference count mismatch");require_bytes(size,pos,4ULL+std::uint64_t(n)*8,"legacy SKETCHES hashes");whole.update(p+pos,size_t(4ULL+std::uint64_t(n)*8));pos+=4;if(legacy_digests){Sha256 digest;digest.update(sketch_digest_domain());digest.update(p+pos,size_t(std::uint64_t(n)*8));legacy_digests->push_back(raw_digest(digest.finish()));}std::uint64_t prev=0;for(std::uint32_t j=0;j<n;++j){auto h=read_le64(p+pos+std::uint64_t(j)*8);if(j&&h<=prev)throw std::runtime_error("SKETCHES hashes not strictly ordered");prev=h;}pos+=std::uint64_t(n)*8;}if(pos!=size)throw std::runtime_error("trailing SKETCHES bytes");if(json_string(complete,"sketches_sha256")!=whole.finish())throw std::runtime_error("legacy SKETCHES SHA256 mismatch");sk->verify_identity();sk->release_pages();}
  PostingIndex pi;pi.postings=std::make_shared<MappedFile>(dir/"POSTINGS.bin");if(json_string(complete,"postings_sha256")!=sha_bytes(pi.postings->data(),pi.postings->size()))throw std::runtime_error("POSTINGS SHA256 mismatch");std::vector<std::uint64_t>counts(refs.size(),0);
  if(!compact){pi.format=PostingIndex::Format::LegacyV3;auto*p=pi.postings->data();auto size=pi.postings->size();require_bytes(size,0,16,"legacy POSTINGS");if(std::memcmp(p,"GAAFPO3\0",8)!=0)throw std::runtime_error("wrong legacy POSTINGS magic");pi.groups=read_le64(p+8);if(pi.groups>(size-16)/12)throw std::runtime_error("impossible legacy posting group count");pi.legacy_hashes.reserve(size_t(pi.groups));pi.legacy_count_offsets.reserve(size_t(pi.groups));std::uint64_t pos=16,prevh=0;for(std::uint64_t g=0;g<pi.groups;++g){require_bytes(size,pos,12,"legacy POSTINGS record");auto h=read_le64(p+pos);pos+=8;if(g&&h<=prevh)throw std::runtime_error("POSTINGS hashes not strictly ordered");prevh=h;pi.legacy_hashes.push_back(h);pi.legacy_count_offsets.push_back(pos);auto n=read_le32(p+pos);pos+=4;if(!n||n>std::uint32_t(c.max_posting))throw std::runtime_error("empty/over-frequency posting");require_bytes(size,pos,std::uint64_t(n)*4,"legacy POSTINGS refs");std::uint32_t prevr=0;for(std::uint32_t j=0;j<n;++j){auto r=read_le32(p+pos+std::uint64_t(j)*4);if(r>=refs.size()||(j&&r<=prevr))throw std::runtime_error("POSTINGS refs not strictly ordered/in range");prevr=r;++counts[r];}pos+=std::uint64_t(n)*4;pi.total_refs+=n;}if(pos!=size)throw std::runtime_error("trailing legacy POSTINGS bytes");}
  else{pi.format=PostingIndex::Format::CompactV4;pi.lookup=std::make_shared<MappedFile>(dir/"POSTINGS_LOOKUP.bin");if(json_string(complete,"postings_lookup_sha256")!=sha_bytes(pi.lookup->data(),pi.lookup->size()))throw std::runtime_error("POSTINGS_LOOKUP SHA256 mismatch");auto*d=pi.postings->data();auto*l=pi.lookup->data();require_bytes(pi.postings->size(),0,16,"compact POSTINGS");require_bytes(pi.lookup->size(),0,24,"compact POSTINGS lookup");if(std::memcmp(d,"GAAFRF4\0",8)!=0||std::memcmp(l,"GAAFLU4\0",8)!=0)throw std::runtime_error("wrong compact postings magic");pi.total_refs=read_le64(d+8);pi.groups=read_le64(l+8);if(read_le64(l+16)!=pi.total_refs||pi.total_refs!=json_uint(complete,"kept_posting_count")||pi.groups!=json_uint(complete,"posting_group_count"))throw std::runtime_error("compact postings header/COMPLETE mismatch");if(pi.total_refs>(std::numeric_limits<std::uint64_t>::max()-16)/4||pi.postings->size()!=16ULL+pi.total_refs*4ULL||pi.groups>(std::numeric_limits<std::uint64_t>::max()-24)/16||pi.lookup->size()!=24ULL+pi.groups*16ULL)throw std::runtime_error("compact postings size closure failure");std::uint64_t prevh=0,prevoff=0;for(std::uint64_t g=0;g<pi.groups;++g){auto h=read_le64(l+24+g*16),off=read_le64(l+24+g*16+8),end=g+1<pi.groups?read_le64(l+24+(g+1)*16+8):pi.total_refs;if((g&&h<=prevh)||off!=prevoff||end<=off||end-off>std::uint64_t(c.max_posting)||end>pi.total_refs)throw std::runtime_error("compact postings lookup ordering/range failure");prevh=h;prevoff=end;std::uint32_t prevr=0;for(std::uint64_t j=off;j<end;++j){auto r=read_le32(d+16+j*4);if(r>=refs.size()||(j>off&&r<=prevr))throw std::runtime_error("compact POSTINGS refs not strictly ordered/in range");prevr=r;++counts[r];}}if((pi.groups&&prevoff!=pi.total_refs)||(!pi.groups&&pi.total_refs))throw std::runtime_error("compact postings final offset mismatch");}
  if(pi.total_refs!=json_uint(complete,"kept_posting_count"))throw std::runtime_error("posting total disagrees with COMPLETE");for(size_t i=0;i<refs.size();++i)if(counts[i]!=refs[i].posting_count)throw std::runtime_error("posting per-reference count mismatch");validate_taxonomy_rebind_contract(dir,complete,refs);complete_map->verify_identity();refs_map->verify_identity();pi.verify_snapshot();if(source_snapshot){source_snapshot->complete_bytes=complete;source_snapshot->refs_bytes=std::move(refs_text);}pi.postings->release_pages();pi.postings->random_access();if(pi.lookup){pi.lookup->release_pages();pi.lookup->random_access();}return{std::move(refs),std::move(pi)};
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

static void compact_legacy_index(const Config&request){
  if(request.index.empty()||request.out_dir.empty())throw std::runtime_error("compact-index requires --index and --out-dir");
  if(fs::exists(request.out_dir))throw std::runtime_error("compact-index destination exists: "+request.out_dir.string());
  auto source_index_path=fs::weakly_canonical(request.index).string();Config c=request;std::vector<std::array<unsigned char,32>>digests;IndexSourceSnapshot source;auto loaded=load_index(c.index,c,&digests,&source);auto&refs=loaded.first;auto&postings=loaded.second;if(postings.format!=PostingIndex::Format::LegacyV3)throw std::runtime_error("compact-index source must be a legacy v3 index");if(digests.size()!=refs.size())throw std::runtime_error("legacy sketch digest closure failure");auto&source_complete=source.complete_bytes;
  auto tmp=temporary_sibling(request.out_dir);if(!fs::create_directories(tmp))throw std::runtime_error("temporary compact-index directory already exists: "+tmp.string());try{
    {std::ofstream ro(tmp/"REFS.tsv",std::ios::binary);ro.write(source.refs_bytes.data(),std::streamsize(source.refs_bytes.size()));ro.close();if(!ro||sha_file(tmp/"REFS.tsv")!=json_string(source_complete,"refs_sha256"))throw std::runtime_error("failed writing verified REFS snapshot");}
    std::array<char,1u<<20>sketch_buffer{};std::ofstream sk;sk.rdbuf()->pubsetbuf(sketch_buffer.data(),std::streamsize(sketch_buffer.size()));sk.open(tmp/"SKETCHES.bin",std::ios::binary);magic(sk,"GAAFSK4\0");put_u32(sk,std::uint32_t(refs.size()));for(size_t i=0;i<refs.size();++i){put_u32(sk,refs[i].sketch_count);sk.write(reinterpret_cast<const char*>(digests[i].data()),std::streamsize(digests[i].size()));}sk.close();if(!sk||fs::file_size(tmp/"SKETCHES.bin")!=12ULL+refs.size()*36ULL)throw std::runtime_error("failed writing compact sketch summary");
    std::array<char,1u<<20>post_buffer{},lookup_buffer{};std::ofstream po,lu;po.rdbuf()->pubsetbuf(post_buffer.data(),std::streamsize(post_buffer.size()));lu.rdbuf()->pubsetbuf(lookup_buffer.data(),std::streamsize(lookup_buffer.size()));po.open(tmp/"POSTINGS.bin",std::ios::binary);lu.open(tmp/"POSTINGS_LOOKUP.bin",std::ios::binary);if(!po||!lu)throw std::runtime_error("cannot create compact postings artifacts");magic(po,"GAAFRF4\0");put_u64(po,postings.total_refs);magic(lu,"GAAFLU4\0");put_u64(lu,postings.groups);put_u64(lu,postings.total_refs);std::uint64_t emitted=0;for(auto h:postings.legacy_hashes){auto span=postings.find(h);put_u64(lu,h);put_u64(lu,emitted);for(std::uint32_t j=0;j<span.count;++j)put_u32(po,span.at(j));emitted+=span.count;}postings.verify_snapshot();po.close();lu.close();if(!po||!lu||emitted!=postings.total_refs||fs::file_size(tmp/"POSTINGS.bin")!=16ULL+emitted*4ULL||fs::file_size(tmp/"POSTINGS_LOOKUP.bin")!=24ULL+postings.groups*16ULL)throw std::runtime_error("compact postings write/size closure failure");
    std::ostringstream complete;complete<<"{\n  \"schema\": \"gtdb-ani-af-index-v4\",\n  \"status\": \"PASS\",\n  \"k\": "<<c.k<<",\n  \"sketch_scale\": "<<c.sketch_scale<<",\n  \"query_sketch_size\": "<<c.query_sketch_size<<",\n  \"max_posting\": "<<c.max_posting<<",\n  \"reference_count\": "<<refs.size()<<",\n  \"kept_posting_count\": "<<postings.total_refs<<",\n  \"posting_group_count\": "<<postings.groups<<",\n  \"reference_manifest_sha256\": \""<<json_string(source_complete,"reference_manifest_sha256")<<"\",\n  \"taxonomy_sha256\": \""<<json_string(source_complete,"taxonomy_sha256")<<"\",\n  \"radii_sha256\": \""<<json_string(source_complete,"radii_sha256")<<"\",\n  \"release_metadata_sha256\": \""<<json_string(source_complete,"release_metadata_sha256")<<"\",\n  \"refs_sha256\": \""<<sha_file(tmp/"REFS.tsv")<<"\",\n  \"sketches_sha256\": \""<<sha_file(tmp/"SKETCHES.bin")<<"\",\n  \"postings_sha256\": \""<<sha_file(tmp/"POSTINGS.bin")<<"\",\n  \"postings_lookup_sha256\": \""<<sha_file(tmp/"POSTINGS_LOOKUP.bin")<<"\",\n  \"manifest_path\": \""<<json_escape(json_string(source_complete,"manifest_path"))<<"\",\n  \"taxonomy_path\": \""<<json_escape(json_string(source_complete,"taxonomy_path"))<<"\",\n  \"radii_path\": \""<<json_escape(json_string(source_complete,"radii_path"))<<"\",\n  \"release_metadata_path\": \""<<json_escape(json_string(source_complete,"release_metadata_path"))<<"\",\n  \"converted_from_schema\": \"gtdb-ani-af-index-v3\",\n  \"source_complete_sha256\": \""<<sha_text(source_complete)<<"\",\n  \"source_index_path\": \""<<json_escape(source_index_path)<<"\"\n}\n";std::ofstream co(tmp/"COMPLETE.json",std::ios::binary);co<<complete.str();co.close();if(!co)throw std::runtime_error("failed writing compact COMPLETE.json");{Config verify;auto verified=load_index(tmp,verify);if(verified.first.size()!=refs.size()||verified.second.groups!=postings.groups||verified.second.total_refs!=postings.total_refs)throw std::runtime_error("converted compact index verification closure failure");}fs::rename(tmp,request.out_dir);
  }catch(...){std::error_code ec;fs::remove_all(tmp,ec);throw;}
}

static void write_binary_exact(const fs::path&path,const std::string&bytes){std::ofstream out(path,std::ios::binary);if(!out)throw std::runtime_error("cannot create: "+path.string());out.write(bytes.data(),std::streamsize(bytes.size()));out.close();if(!out)throw std::runtime_error("failed writing: "+path.string());}
static std::string materialize_verified_file(const fs::path&source,const fs::path&destination,const std::string&expected_sha){
  std::error_code link_error;fs::create_hard_link(source,destination,link_error);std::string mode="hardlink";
  if(link_error){std::error_code copy_error;if(!fs::copy_file(source,destination,fs::copy_options::none,copy_error)||copy_error)throw std::runtime_error("cannot hardlink or copy core artifact "+source.string()+": hardlink="+link_error.message()+" copy="+copy_error.message());mode="copy";}
  if(!fs::is_regular_file(destination)||sha_file(destination)!=expected_sha)throw std::runtime_error("materialized core artifact SHA256 mismatch: "+destination.string());
  return mode;
}
static void rebind_taxonomy(const Config&request){
  if(request.index.empty()||request.taxonomy.empty()||request.out_dir.empty())throw std::runtime_error("rebind-taxonomy requires --index --taxonomy and --out-dir");
  if(fs::exists(request.out_dir))throw std::runtime_error("rebind-taxonomy destination exists: "+request.out_dir.string());
  if(!fs::is_regular_file(request.taxonomy))throw std::runtime_error("new taxonomy is not a regular file: "+request.taxonomy.string());
  Config source_config=request;IndexSourceSnapshot source_snapshot;std::vector<RefRow>refs;PostingIndex::Format source_format=PostingIndex::Format::LegacyV3;std::uint64_t source_groups=0,source_total_refs=0;
  // Finish and destroy all mmap snapshots before creating hardlinks: creating
  // a link legitimately changes inode metadata and would otherwise trip the
  // source snapshot's mutation detector even though no bytes changed.
  {auto loaded=load_index(request.index,source_config,nullptr,&source_snapshot);refs=std::move(loaded.first);source_format=loaded.second.format;source_groups=loaded.second.groups;source_total_refs=loaded.second.total_refs;}
  auto refset=reference_accessions(refs);auto taxonomy=load_taxonomy_accessions(request.taxonomy);require_same_accessions(refset,taxonomy.values);if(taxonomy.rows!=refs.size())throw std::runtime_error("new taxonomy row count differs from REFS count");
  const auto&source=source_snapshot.complete_bytes;const auto source_schema=json_string(source,"schema");if(source_schema!="gtdb-ani-af-index-v3"&&source_schema!="gtdb-ani-af-index-v4")throw std::runtime_error("rebind-taxonomy source is not a supported index");const bool compact=source_schema=="gtdb-ani-af-index-v4";if((source_format==PostingIndex::Format::CompactV4)!=compact)throw std::runtime_error("source schema/index format disagreement");
  auto tmp=temporary_sibling(request.out_dir);if(!fs::create_directories(tmp))throw std::runtime_error("temporary rebind directory already exists: "+tmp.string());try{
    std::vector<std::string>modes;auto materialize=[&](const char*name,const char*digest_key){modes.push_back(materialize_verified_file(request.index/name,tmp/name,json_string(source,digest_key)));};materialize("REFS.tsv","refs_sha256");materialize("SKETCHES.bin","sketches_sha256");materialize("POSTINGS.bin","postings_sha256");if(compact)materialize("POSTINGS_LOOKUP.bin","postings_lookup_sha256");
    std::string core_mode=modes.front();if(!std::all_of(modes.begin(),modes.end(),[&](const std::string&x){return x==core_mode;}))core_mode="mixed";
    write_binary_exact(tmp/"SOURCE_COMPLETE.json",source);if(sha_file(tmp/"SOURCE_COMPLETE.json")!=sha_text(source))throw std::runtime_error("failed to preserve source COMPLETE receipt exactly");
    std::error_code copy_error;if(!fs::copy_file(request.taxonomy,tmp/"TAXONOMY.tsv",fs::copy_options::none,copy_error)||copy_error)throw std::runtime_error("cannot copy bound taxonomy: "+copy_error.message());const auto taxonomy_sha=sha_file(request.taxonomy);if(sha_file(tmp/"TAXONOMY.tsv")!=taxonomy_sha)throw std::runtime_error("bound taxonomy copy SHA256 mismatch");
    // Reopen the source after materialization.  This catches any byte change
    // that raced the link/copy operation without depending on inode metadata.
    {Config source_recheck;auto reopened=load_index(request.index,source_recheck);if(reopened.first.size()!=refs.size()||reopened.second.groups!=source_groups||reopened.second.total_refs!=source_total_refs||sha_file(request.index/"COMPLETE.json")!=sha_text(source))throw std::runtime_error("source index changed while taxonomy was rebound");}
    std::ostringstream complete;complete<<"{\n  \"schema\": \""<<source_schema<<"\",\n  \"status\": \"PASS\",\n  \"k\": "<<source_config.k<<",\n  \"sketch_scale\": "<<source_config.sketch_scale<<",\n  \"query_sketch_size\": "<<source_config.query_sketch_size<<",\n  \"max_posting\": "<<source_config.max_posting<<",\n  \"reference_count\": "<<refs.size()<<",\n  \"kept_posting_count\": "<<source_total_refs<<",\n";if(compact)complete<<"  \"posting_group_count\": "<<source_groups<<",\n";complete
      <<"  \"reference_manifest_sha256\": \""<<json_string(source,"reference_manifest_sha256")<<"\",\n"
      <<"  \"taxonomy_sha256\": \""<<taxonomy_sha<<"\",\n"
      <<"  \"radii_sha256\": \""<<json_string(source,"radii_sha256")<<"\",\n"
      <<"  \"release_metadata_sha256\": \""<<json_string(source,"release_metadata_sha256")<<"\",\n"
      <<"  \"refs_sha256\": \""<<sha_file(tmp/"REFS.tsv")<<"\",\n"
      <<"  \"sketches_sha256\": \""<<sha_file(tmp/"SKETCHES.bin")<<"\",\n"
      <<"  \"postings_sha256\": \""<<sha_file(tmp/"POSTINGS.bin")<<"\",\n";if(compact)complete<<"  \"postings_lookup_sha256\": \""<<sha_file(tmp/"POSTINGS_LOOKUP.bin")<<"\",\n";complete
      <<"  \"manifest_path\": \""<<json_escape(json_string(source,"manifest_path"))<<"\",\n"
      <<"  \"taxonomy_path\": \""<<json_escape(fs::weakly_canonical(request.taxonomy).string())<<"\",\n"
      <<"  \"radii_path\": \""<<json_escape(json_string(source,"radii_path"))<<"\",\n"
      <<"  \"release_metadata_path\": \""<<json_escape(json_string(source,"release_metadata_path"))<<"\",\n"
      <<"  \"taxonomy_rebound_from_schema\": \""<<source_schema<<"\",\n"
      <<"  \"source_index_path\": \""<<json_escape(fs::weakly_canonical(request.index).string())<<"\",\n"
      <<"  \"source_receipt_file\": \"SOURCE_COMPLETE.json\",\n"
      <<"  \"source_complete_sha256\": \""<<sha_text(source)<<"\",\n"
      <<"  \"source_taxonomy_sha256\": \""<<json_string(source,"taxonomy_sha256")<<"\",\n"
      <<"  \"taxonomy_bound_file\": \"TAXONOMY.tsv\",\n"
      <<"  \"taxonomy_entry_count\": "<<taxonomy.rows<<",\n"
      <<"  \"taxonomy_accession_set_sha256\": \""<<taxonomy.set_sha256<<"\",\n"
      <<"  \"core_materialization\": \""<<core_mode<<"\"\n}\n";
    write_binary_exact(tmp/"COMPLETE.json",complete.str());{Config verify;auto rebound=load_index(tmp,verify);if(rebound.first.size()!=refs.size()||rebound.second.groups!=source_groups||rebound.second.total_refs!=source_total_refs)throw std::runtime_error("taxonomy rebound index verification closure failure");}fs::rename(tmp,request.out_dir);
  }catch(...){std::error_code ec;fs::remove_all(tmp,ec);throw;}
}

static void header(std::ostream&o){o<<"Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n";}
static void row(std::ostream&o,const fs::path&r,const fs::path&q,const Estimate&e){o<<r.string()<<'\t'<<q.string()<<'\t'<<std::fixed<<std::setprecision(6)<<e.ani<<'\t'<<e.af_ref<<'\t'<<e.af_query<<'\n';}
struct SearchScratch {std::vector<std::uint16_t>counts;std::vector<std::uint32_t>touched;};
static std::string search_one(const fs::path&query,const std::vector<RefRow>&refs,const PostingIndex&postings,const Config&c,SearchScratch&scratch){
  postings.verify_snapshot();
  if(scratch.counts.size()!=refs.size())scratch.counts.assign(refs.size(),0);
  if(scratch.touched.capacity()<std::min<size_t>(refs.size(),size_t(c.query_sketch_size)*64))scratch.touched.reserve(std::min<size_t>(refs.size(),size_t(c.query_sketch_size)*64));
  auto qg=read_fasta(query);auto qs=make_scaled_sketch(qg,c.k,c.sketch_scale,c.query_sketch_size);for(auto h:qs){auto p=postings.find(h);for(std::uint32_t j=0;j<p.count;++j){auto id=p.at(j);if(id>=scratch.counts.size())throw std::runtime_error("posting reference ID changed/out of range");if(!scratch.counts[id])scratch.touched.push_back(id);if(scratch.counts[id]==std::numeric_limits<std::uint16_t>::max())throw std::runtime_error("candidate counter overflow");++scratch.counts[id];}}
  struct Rank{std::uint32_t id;int shared;double hit_fraction;};std::vector<Rank>rank;rank.reserve(scratch.touched.size());for(auto id:scratch.touched){auto n=scratch.counts[id];scratch.counts[id]=0;if(n>=c.min_candidate_hits)rank.push_back({id,n,qs.empty()?0.0:double(n)/qs.size()});}scratch.touched.clear();std::sort(rank.begin(),rank.end(),[&](const Rank&a,const Rank&b){if(a.hit_fraction!=b.hit_fraction)return a.hit_fraction>b.hit_fraction;if(a.shared!=b.shared)return a.shared>b.shared;return refs[a.id].path.string()<refs[b.id].path.string();});if(int(rank.size())>c.max_candidates)rank.resize(c.max_candidates);
  struct Result{Estimate e;fs::path ref;};std::vector<Result>results;for(const auto&x:rank){const auto&m=refs[x.id];if(!fs::is_regular_file(m.path)||fs::file_size(m.path)!=m.bytes||fs::last_write_time(m.path).time_since_epoch().count()!=m.mtime)throw std::runtime_error("reference metadata changed: "+m.path.string());auto rg=read_fasta(m.path);if(c.deep_verify&&rg.content_sha256!=m.content_sha256)throw std::runtime_error("reference decompressed-content SHA256 changed: "+m.path.string());results.push_back({estimate(qg,rg,c),rg.path});}postings.verify_snapshot();std::sort(results.begin(),results.end(),[](const Result&a,const Result&b){if(a.e.ani!=b.e.ani)return a.e.ani>b.e.ani;if(a.e.af_query!=b.e.af_query)return a.e.af_query>b.e.af_query;if(a.e.af_ref!=b.e.af_ref)return a.e.af_ref>b.e.af_ref;return a.ref.string()<b.ref.string();});if(int(results.size())>c.top)results.resize(c.top);std::ostringstream o;for(const auto&r:results)row(o,r.ref,qg.path,r.e);return o.str();
}
static void write_once(const fs::path&p,const std::string&s){if(p.empty()){std::cout<<s;return;}if(fs::exists(p))throw std::runtime_error("output exists: "+p.string());std::ofstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot write: "+p.string());f.write(s.data(),static_cast<std::streamsize>(s.size()));f.close();if(!f)throw std::runtime_error("failed writing: "+p.string());}
static Config parse(int argc,char**argv,int start){
  Config c;
  for(int i=start;i<argc;++i){
    std::string a=argv[i];
    auto v=[&](){if(++i>=argc)throw std::runtime_error("missing value for "+a);return std::string(argv[i]);};
    if(a=="--query")c.query=v();
    else if(a=="--reference")c.reference=v();
    else if(a=="--list"||a=="-l")c.list=v();
    else if(a=="--manifest")c.manifest=v();
    else if(a=="--taxonomy")c.taxonomy=v();
    else if(a=="--radii")c.radii=v();
    else if(a=="--release-metadata")c.release_metadata=v();
    else if(a=="--out-dir")c.out_dir=v();
    else if(a=="--index")c.index=v();
    else if(a=="--queries")c.queries=v();
    else if(a=="--out")c.out=v();
    else if(a=="--stats-out")c.stats_out=v();
    else if(a=="--index-format"){auto f=v();if(f=="compact-v4")c.compact_index=true;else if(f=="legacy-v3")c.compact_index=false;else throw std::runtime_error("--index-format must be compact-v4 or legacy-v3");}
    else if(a=="--k"){c.k=std::stoi(v());c.k_explicit=true;}
    else if(a=="--sketch-scale"){c.sketch_scale=std::stoi(v());c.sketch_scale_explicit=true;}
    else if(a=="--query-sketch-size"){c.query_sketch_size=std::stoi(v());c.query_sketch_explicit=true;}
    else if(a=="--max-posting"){c.max_posting=std::stoi(v());c.max_posting_explicit=true;}
    else if(a=="--min-candidate-hits")c.min_candidate_hits=std::stoi(v());
    else if(a=="--window")c.window=std::stoi(v());
    else if(a=="--min-window")c.min_window=std::stoi(v());
    else if(a=="--min-seed-hits")c.min_seed_hits=std::stoi(v());
    else if(a=="--max-occ")c.max_occ=std::stoi(v());
    else if(a=="--align-band")c.align_band=std::stoi(v());
    else if(a=="--min-identity")c.min_identity=std::stod(v());
    else if(a=="--min-af")c.min_af=std::stod(v());
    else if(a=="--profile"){if(v()!="medium")throw std::runtime_error("only --profile medium is defined");}
    else if(a=="--edge-mode"){if(v()!="sparse")throw std::runtime_error("only --edge-mode sparse is defined");}
    else if(a=="--triangle-mode"){auto mode=v();if(mode=="exact")c.triangle_fast=false;else if(mode=="sketch")c.triangle_fast=true;else throw std::runtime_error("--triangle-mode must be exact or sketch");}
    else if(a=="--triangle-min-shared")c.triangle_min_shared=std::stoi(v());
    else if(a=="--triangle-max-candidates")c.triangle_max_candidates=std::stoi(v());
    else if(a=="--triangle-sketch-size")c.triangle_sketch_size=std::stoi(v());
    else if(a=="--triangle-max-posting")c.triangle_max_posting=std::stoi(v());
    else if(a=="--triangle-reference-block-size")c.triangle_reference_block_size=std::stoi(v());
    else if(a=="--triangle-global-prefix-bits")c.triangle_global_prefix_bits=std::stoi(v());
    else if(a=="--triangle-global-seed-join")c.triangle_global_seed_join=true;
    else if(a=="--triangle-candidates-out")c.triangle_candidates_out=v();
    else if(a=="--triangle-candidates-only")c.triangle_candidates_only=true;
    else if(a=="--max-candidates")c.max_candidates=std::stoi(v());
    else if(a=="--top")c.top=std::stoi(v());
    else if(a=="--threads")c.threads=std::stoi(v());
    else if(a=="--deep-verify")c.deep_verify=true;
    else throw std::runtime_error("unknown option: "+a);
  }
  if(c.k<5||c.k>31||c.window<50||c.min_window<20||c.min_window>c.window||c.min_seed_hits<1||c.max_occ<1||c.align_band<1||c.sketch_scale<1||c.query_sketch_size<10||c.query_sketch_size>1000||c.max_posting<1||c.min_candidate_hits<1||c.max_candidates<1||c.top<1||c.threads<1||c.min_identity<0||c.min_identity>100||c.min_af<0||c.min_af>100||c.triangle_min_shared<1||c.triangle_max_candidates<0||c.triangle_sketch_size<0||c.triangle_max_posting<0||c.triangle_reference_block_size<1||c.triangle_reference_block_size>64||c.triangle_global_prefix_bits<12||c.triangle_global_prefix_bits>28)throw std::runtime_error("invalid numeric option");
  if(!c.triangle_fast&&!c.triangle_candidates_out.empty())throw std::runtime_error("--triangle-candidates-out requires --triangle-mode sketch");
  if(c.triangle_fast&&c.triangle_global_seed_join)throw std::runtime_error("--triangle-global-seed-join requires --triangle-mode exact");
  if(c.triangle_global_seed_join&&c.triangle_reference_block_size!=1)throw std::runtime_error("--triangle-global-seed-join cannot be combined with fused reference blocks");
  if(c.triangle_candidates_only&&(!c.triangle_fast||c.triangle_candidates_out.empty()))throw std::runtime_error("--triangle-candidates-only requires --triangle-mode sketch and --triangle-candidates-out");
  return c;
}

static Candidate gapped_candidate_reference_for_test(const std::string&q,const std::string&r,std::uint32_t q0,std::uint32_t q1,int expected_diag,int band){Candidate out;out.q0=q0;out.q1=q1;out.diag=expected_diag;int rb=std::max(0,int(q0)+expected_diag-band),re=std::min<int>(r.size(),int(q1)+expected_diag+band);if(re<=rb||q1<=q0)return out;size_t n=q1-q0,m=size_t(re-rb),stride=m+1;std::vector<int>dp((n+1)*stride,0);std::vector<unsigned char>tr((n+1)*stride,0);constexpr int gap=-4;for(size_t i=1;i<=n;++i){dp[i*stride]=int(i)*gap;tr[i*stride]=1;}for(size_t i=1;i<=n;++i)for(size_t j=1;j<=m;++j){char qb=q[q0+i-1],rbv=r[size_t(rb)+j-1];int diag=dp[(i-1)*stride+j-1]+((callable(qb)&&callable(rbv)&&qb==rbv)?2:-3),up=dp[(i-1)*stride+j]+gap,left=dp[i*stride+j-1]+gap;int best=diag;unsigned char move=0;if(up>best){best=up;move=1;}if(left>best){best=left;move=2;}dp[i*stride+j]=best;tr[i*stride+j]=move;}size_t i=n,j=0;for(size_t x=1;x<=m;++x)if(dp[n*stride+x]>dp[n*stride+j])j=x;while(i){auto move=tr[i*stride+j];if(move==0){auto qp=std::uint32_t(q0+i-1),rp=std::uint32_t(size_t(rb)+j-1);out.aligned.push_back({qp,rp});--i;--j;}else if(move==1){--i;}else{if(!j)throw std::runtime_error("reference traceback underflow");--j;}}std::reverse(out.aligned.begin(),out.aligned.end());for(auto[qp,rp]:out.aligned)if(callable(q[qp])&&callable(r[rp])){++out.compared;if(q[qp]==r[rp])++out.matches;}out.identity=out.compared?100.0*out.matches/out.compared:0;return out;}

static std::string random_dna(size_t n,std::uint64_t seed=0x123456789abcdef0ULL){std::uint64_t x=seed;std::string s;s.reserve(n);for(size_t i=0;i<n;++i){x^=x<<13;x^=x>>7;x^=x<<17;s+="ACGT"[x&3];}return s;}
static std::string mutate_95(std::string s,std::uint64_t seed){std::vector<size_t>p(s.size());std::iota(p.begin(),p.end(),0);std::mt19937_64 rng(seed);std::shuffle(p.begin(),p.end(),rng);for(size_t i=0;i<s.size()/20;++i){auto&b=s[p[i]];b="ACGT"[(base_code(b)+1)%4];}return s;}
static size_t shared_hashes(const std::vector<std::uint64_t>&a,const std::vector<std::uint64_t>&b){size_t n=0,i=0,j=0;while(i<a.size()&&j<b.size()){if(a[i]<b[j])++i;else if(b[j]<a[i])++j;else{++n;++i;++j;}}return n;}

struct TriangleCandidatePair {
  std::uint32_t i=0,j=0,shared=0;
  bool selected_by_i=false,selected_by_j=false;
};
struct TriangleCandidatePlan {
  std::vector<std::vector<std::uint32_t>> by_reference;
  std::vector<TriangleCandidatePair> pairs;
  std::uint64_t sketch_memberships=0,skipped_high_frequency_lookups=0;
};
struct TriangleSketchMember {std::uint64_t hash=0;std::uint32_t genome=0;};
struct TriangleSketchRange {std::uint64_t hash=0;std::size_t begin=0,end=0;};

static TriangleCandidatePlan build_triangle_candidate_plan(const std::vector<Genome>&genomes,const Config&c){
  if(genomes.size()>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("too many triangle genomes for candidate index");
  const std::size_t n=genomes.size();std::vector<std::vector<std::uint64_t>>sketches(n);std::atomic<std::size_t>next{0};std::exception_ptr error;std::mutex em;
  auto sketch_worker=[&](){for(;;){auto i=next.fetch_add(1);if(i>=n)return;try{sketches[i]=make_scaled_sketch(genomes[i],c.k,c.sketch_scale,c.triangle_sketch_size);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();return;}}};
  std::vector<std::thread>workers;const auto nw=std::min<std::size_t>(std::size_t(c.threads),n);for(std::size_t w=0;w<nw;++w)workers.emplace_back(sketch_worker);for(auto&t:workers)t.join();if(error)std::rethrow_exception(error);
  std::size_t total=0;for(const auto&s:sketches){if(s.size()>std::numeric_limits<std::size_t>::max()-total)throw std::runtime_error("triangle sketch membership count overflow");total+=s.size();}
  std::vector<TriangleSketchMember>members;members.reserve(total);for(std::size_t i=0;i<n;++i)for(auto h:sketches[i])members.push_back({h,std::uint32_t(i)});
  std::sort(members.begin(),members.end(),[](const auto&a,const auto&b){return std::tie(a.hash,a.genome)<std::tie(b.hash,b.genome);});
  std::vector<TriangleSketchRange>ranges;ranges.reserve(members.size());for(std::size_t begin=0;begin<members.size();){std::size_t end=begin+1;while(end<members.size()&&members[end].hash==members[begin].hash)++end;ranges.push_back({members[begin].hash,begin,end});begin=end;}
  // A 16-bit prefix directory preserves exact lookup semantics while reducing
  // each sorted-range lookup from the whole index to one small hash bucket.
  std::vector<std::size_t>prefix(65537);std::size_t rp=0;for(std::size_t p=0;p<65536;++p){while(rp<ranges.size()&&(ranges[rp].hash>>48)<p)++rp;prefix[p]=rp;}prefix[65536]=ranges.size();
  auto find_range=[&](std::uint64_t h)->const TriangleSketchRange*{auto p=std::size_t(h>>48);auto first=ranges.begin()+std::ptrdiff_t(prefix[p]),last=ranges.begin()+std::ptrdiff_t(prefix[p+1]);auto it=std::lower_bound(first,last,h,[](const auto&r,std::uint64_t value){return r.hash<value;});return it!=last&&it->hash==h?&*it:nullptr;};
  struct Rank {std::uint32_t id=0,shared=0;};std::vector<std::vector<Rank>>directed(n);std::vector<std::string>path_keys;path_keys.reserve(n);for(const auto&g:genomes)path_keys.push_back(g.path.string());
  next=0;error=nullptr;workers.clear();std::atomic<std::uint64_t>skipped_groups{0};
  auto rank_worker=[&](){
    std::vector<std::uint32_t>counts(n,0),touched;std::uint64_t local_skipped=0;
    try{for(;;){auto i=next.fetch_add(1);if(i>=n)break;touched.clear();for(auto h:sketches[i]){auto*r=find_range(h);if(!r)continue;auto posting=r->end-r->begin;if(c.triangle_max_posting>0&&posting>std::size_t(c.triangle_max_posting)){++local_skipped;continue;}for(auto p=r->begin;p<r->end;++p){auto id=members[p].genome;if(id==i)continue;if(!counts[id])touched.push_back(id);if(counts[id]==std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("triangle shared-sketch counter overflow");++counts[id];}}auto&ranked=directed[i];ranked.reserve(touched.size());for(auto id:touched){auto shared=counts[id];counts[id]=0;if(shared>=std::uint32_t(c.triangle_min_shared))ranked.push_back({id,shared});}std::sort(ranked.begin(),ranked.end(),[&](const Rank&a,const Rank&b){if(a.shared!=b.shared)return a.shared>b.shared;return path_keys[a.id]<path_keys[b.id];});if(c.triangle_max_candidates>0&&ranked.size()>std::size_t(c.triangle_max_candidates))ranked.resize(std::size_t(c.triangle_max_candidates));}}
    catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();}
    skipped_groups.fetch_add(local_skipped,std::memory_order_relaxed);
  };
  for(std::size_t w=0;w<nw;++w)workers.emplace_back(rank_worker);
  for(auto&t:workers)t.join();
  if(error)std::rethrow_exception(error);
  struct Choice {std::uint64_t key=0;std::uint32_t shared=0;std::uint8_t mask=0;};std::vector<Choice>choices;for(std::size_t i=0;i<n;++i){if(directed[i].size()>std::numeric_limits<std::size_t>::max()-choices.size())throw std::runtime_error("triangle candidate choice count overflow");for(const auto&r:directed[i]){auto a=std::min<std::uint32_t>(std::uint32_t(i),r.id),b=std::max<std::uint32_t>(std::uint32_t(i),r.id);choices.push_back({(std::uint64_t(a)<<32)|b,r.shared,std::uint8_t(i==a?1:2)});}}
  std::sort(choices.begin(),choices.end(),[](const Choice&a,const Choice&b){return std::tie(a.key,a.mask)<std::tie(b.key,b.mask);});TriangleCandidatePlan plan;plan.by_reference.resize(n);plan.sketch_memberships=members.size();plan.skipped_high_frequency_lookups=skipped_groups.load();
  for(std::size_t begin=0;begin<choices.size();){std::size_t end=begin+1;auto shared=choices[begin].shared;std::uint8_t mask=choices[begin].mask;while(end<choices.size()&&choices[end].key==choices[begin].key){if(choices[end].shared!=shared)throw std::runtime_error("asymmetric shared-sketch count");mask=std::uint8_t(mask|choices[end].mask);++end;}auto i=std::uint32_t(choices[begin].key>>32),j=std::uint32_t(choices[begin].key);plan.pairs.push_back({i,j,shared,bool(mask&1),bool(mask&2)});plan.by_reference[i].push_back(j);begin=end;}
  return plan;
}

static void write_triangle_candidates(const fs::path&p,const TriangleCandidatePlan&plan,const std::vector<Genome>&genomes,const Config&c){
  if(p.empty())return;
  if(fs::exists(p))throw std::runtime_error("output exists: "+p.string());
  fs::path tmp=p;tmp+=".tmp";
  if(fs::exists(tmp))throw std::runtime_error("stale candidate-output temporary: "+tmp.string());
  std::ofstream out(tmp);if(!out)throw std::runtime_error("cannot write triangle candidates: "+tmp.string());
  out<<"Ref_file\tQuery_file\tshared_sketches\tselected_by_ref_topk\tselected_by_query_topk\tsketch_scale\tsketch_size_cap\tmax_posting\tmin_shared\tmax_candidates_per_genome\n";
  for(const auto&x:plan.pairs)out<<genomes[x.i].path.string()<<'\t'<<genomes[x.j].path.string()<<'\t'<<x.shared<<'\t'<<(x.selected_by_i?1:0)<<'\t'<<(x.selected_by_j?1:0)<<'\t'<<c.sketch_scale<<'\t'<<c.triangle_sketch_size<<'\t'<<c.triangle_max_posting<<'\t'<<c.triangle_min_shared<<'\t'<<c.triangle_max_candidates<<'\n';
  out.close();if(!out)throw std::runtime_error("failed writing triangle candidates: "+tmp.string());fs::rename(tmp,p);
}
static void triangle(const Config&c);
static void self_test(){
  Config c;c.k=13;c.window=300;c.min_window=80;auto refseq=random_dna(10000),q1=refseq.substr(200,1400),q2=reverse_complement(refseq.substr(1900,1300));for(size_t i=37;i<q1.size();i+=101)q1[i]="ACGT"[(base_code(q1[i])+1)%4];q1.replace(500,8,"NNNNNNNN");
  auto packed=gapped_candidate(q1,refseq,0,std::uint32_t(q1.size()),200,128),reference=gapped_candidate_reference_for_test(q1,refseq,0,std::uint32_t(q1.size()),200,128);if(std::tie(packed.compared,packed.matches,packed.identity,packed.aligned)!=std::tie(reference.compared,reference.matches,reference.identity,reference.aligned))throw std::runtime_error("packed traceback differs from full-matrix recurrence/tie rules");
  for(int t=0;t<24;++t){auto rr=random_dna(size_t(420+t),0xe500000000000000ULL+std::uint64_t(t));auto start=size_t(20+t%31),len=size_t(90+(t*17)%180);auto qq=rr.substr(start,len);for(size_t i=size_t(3+t%7);i<qq.size();i+=size_t(19+t%11))qq[i]="ACGT"[(base_code(qq[i])+1)%4];if(t%3==0&&qq.size()>40)qq.replace(31,5,"NNNNN");auto x=gapped_candidate(qq,rr,0,std::uint32_t(qq.size()),int(start),32),y=gapped_candidate_reference_for_test(qq,rr,0,std::uint32_t(qq.size()),int(start),32);if(std::tie(x.compared,x.matches,x.identity,x.aligned)!=std::tie(y.compared,y.matches,y.identity,y.aligned))throw std::runtime_error("randomized packed traceback equivalence failed");}for(const auto&tie_case:std::vector<std::pair<std::string,std::string>>{{"AAAAAAAAAAAAAAAA","AAAAAAAAAAAAAAAAAAAA"},{"ACACACACACACACAC","ACACACACACACACACACAC"},{"NNNNACGTACGTNNNN","TTTTACGTACGTTTTT"}}){auto x=gapped_candidate(tie_case.first,tie_case.second,0,std::uint32_t(tie_case.first.size()),0,8),y=gapped_candidate_reference_for_test(tie_case.first,tie_case.second,0,std::uint32_t(tie_case.first.size()),0,8);if(std::tie(x.compared,x.matches,x.identity,x.aligned)!=std::tie(y.compared,y.matches,y.identity,y.aligned))throw std::runtime_error("tie-case packed traceback equivalence failed");}
  Genome r=parse_fasta_text("ref",">r1\n"+refseq.substr(0,5000)+"\n>r2\n"+refseq.substr(5000)+"\n"),q=parse_fasta_text("query",">forward\n"+q1+"\n>reverse\n"+q2+"\n");auto a=estimate(q,r,c),b=estimate(q,r,c);auto reusable_index=reference_seeds(r,c);auto q_reverse=make_reverse_contigs(q);auto indexed=estimate_with_index(q,r,reusable_index,q_reverse,c);if(a.aligned!=std::uint64_t(a.af_query*q.callable/100.0+0.5)||a.aligned!=std::uint64_t(a.af_ref*r.callable/100.0+0.5))throw std::runtime_error("one-to-one AF accounting failed");if(std::tie(a.ani,a.af_ref,a.af_query,a.aligned,a.matches)!=std::tie(b.ani,b.af_ref,b.af_query,b.aligned,b.matches))throw std::runtime_error("determinism failed");if(std::tie(a.ani,a.af_ref,a.af_query,a.aligned,a.matches)!=std::tie(indexed.ani,indexed.af_ref,indexed.af_query,indexed.aligned,indexed.matches))throw std::runtime_error("reusable reference index path differs from one-shot estimate");std::ostringstream one_shot_row,indexed_row;row(one_shot_row,r.path,q.path,a);row(indexed_row,r.path,q.path,indexed);if(one_shot_row.str()!=indexed_row.str())throw std::runtime_error("reusable reference index row is not byte-exact");Config fc=c;fc.sketch_scale=1;fc.triangle_sketch_size=0;fc.triangle_max_posting=0;fc.triangle_min_shared=1;fc.triangle_max_candidates=0;fc.threads=2;Genome unrelated=parse_fasta_text("unrelated",">u\n"+random_dna(4000,0xf57c9e3779b9ULL)+"\n");auto candidate_plan=build_triangle_candidate_plan(std::vector<Genome>{r,q,unrelated},fc);bool found_candidate=false;for(const auto&x:candidate_plan.pairs)if(x.i==0&&x.j==1)found_candidate=true;if(!found_candidate)throw std::runtime_error("sketch candidate generator missed deterministic related pair");
  auto longseq=random_dna(3100000,0xa17c9e3779b9ULL);Genome longref=parse_fasta_text("longref",">r\n"+longseq+"\n"),d1=parse_fasta_text("d1",">r\n"+random_dna(3100000,0xb27c9e3779b9ULL)+"\n"),d2=parse_fasta_text("d2",">r\n"+random_dna(3100000,0xc37c9e3779b9ULL)+"\n");auto exact1=parse_fasta_text("exact1",">q\n"+longseq.substr(777777,1000)+"\n"),one=parse_fasta_text("one",">q\n"+mutate_95(longseq.substr(777777,1000),11)+"\n"),ten=parse_fasta_text("ten",">q\n"+mutate_95(longseq.substr(1777777,10000),17)+"\n");auto ls=make_scaled_sketch(longref,c.k,c.sketch_scale,0),s1=make_scaled_sketch(one,c.k,c.sketch_scale,c.query_sketch_size),s10=make_scaled_sketch(ten,c.k,c.sketch_scale,c.query_sketch_size),se=make_scaled_sketch(exact1,c.k,c.sketch_scale,c.query_sketch_size),sd1=make_scaled_sketch(d1,c.k,c.sketch_scale,0),sd2=make_scaled_sketch(d2,c.k,c.sketch_scale,0);auto hit1=shared_hashes(s1,ls),hit10=shared_hashes(s10,ls);if(shared_hashes(se,ls)!=se.size()||hit1<size_t(c.min_candidate_hits)||hit10<size_t(c.min_candidate_hits)||hit1<=std::max(shared_hashes(s1,sd1),shared_hashes(s1,sd2))||hit10<=std::max(shared_hashes(s10,sd1),shared_hashes(s10,sd2)))throw std::runtime_error("3.1Mb 1kb/10kb 95%-subset candidate top recall failed");if(sha_text("abc")!="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")throw std::runtime_error("SHA256 failed");
  auto tmp=fs::temp_directory_path()/("gaaf_compact_selftest_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));fs::create_directory(tmp);std::uintmax_t legacy_bytes=0,compact_bytes=0;try{
    std::vector<std::string>seqs={refseq.substr(0,3200),refseq.substr(1600,3200),random_dna(3200,0xd47c9e3779b9ULL)};for(size_t i=0;i<seqs.size();++i){std::ofstream f(tmp/("r"+std::to_string(i)+".fna"),std::ios::binary);f<<">r"<<i<<"\n"<<seqs[i]<<"\n";}
    std::string fa=">x\n"+refseq.substr(0,1200)+"\n";{std::ofstream f(tmp/"x.fna",std::ios::binary);f<<fa;}gzFile z=gzopen((tmp/"x.fna.gz").string().c_str(),"wb");if(!z||gzwrite(z,fa.data(),unsigned(fa.size()))!=int(fa.size())||gzclose(z)!=Z_OK)throw std::runtime_error("gzip write failed");auto plain=read_fasta(tmp/"x.fna"),gz=read_fasta(tmp/"x.fna.gz");if(plain.content_sha256!=gz.content_sha256||plain.callable!=gz.callable||plain.contigs[0].seq!=gz.contigs[0].seq)throw std::runtime_error("gzip/plain parity failed");
    {std::ofstream f(tmp/"manifest.tsv");f<<"Ref_file\n";for(size_t i=0;i<seqs.size();++i)f<<"r"<<i<<".fna\n";}for(const auto&name:{"taxonomy.tsv","radii.tsv","release.txt"}){std::ofstream f(tmp/name);f<<name<<"\n";}
    Config ic;ic.k=13;ic.sketch_scale=4;ic.query_sketch_size=100;ic.max_posting=10;ic.threads=2;ic.manifest=tmp/"manifest.tsv";ic.taxonomy=tmp/"taxonomy.tsv";ic.radii=tmp/"radii.tsv";ic.release_metadata=tmp/"release.txt";ic.out_dir=tmp/"v3";ic.compact_index=false;build_index(ic);auto direct=ic;direct.out_dir=tmp/"v4_direct";direct.compact_index=true;build_index(direct);Config convert;convert.index=tmp/"v3";convert.out_dir=tmp/"v4";compact_legacy_index(convert);for(const auto&name:{"REFS.tsv","SKETCHES.bin","POSTINGS.bin","POSTINGS_LOOKUP.bin"})if(sha_file(tmp/"v4"/name)!=sha_file(tmp/"v4_direct"/name))throw std::runtime_error("direct-v4 and converted-v4 core artifacts differ");{std::ifstream v3_complete(tmp/"v3"/"COMPLETE.json"),v4_complete(tmp/"v4"/"COMPLETE.json");std::ostringstream as,bs;as<<v3_complete.rdbuf();bs<<v4_complete.rdbuf();for(const auto&key:{"reference_manifest_sha256","taxonomy_sha256","radii_sha256","release_metadata_sha256","manifest_path","taxonomy_path","radii_path","release_metadata_path"})if(json_string(as.str(),key)!=json_string(bs.str(),key))throw std::runtime_error("converted v4 provenance differs from v3");if(json_string(bs.str(),"converted_from_schema")!="gtdb-ani-af-index-v3"||json_string(bs.str(),"source_complete_sha256")!=sha_file(tmp/"v3"/"COMPLETE.json")||json_string(bs.str(),"source_index_path")!=fs::weakly_canonical(tmp/"v3").string())throw std::runtime_error("converted v4 source provenance closure failed");}
    Config l3,l4;auto v3=load_index(tmp/"v3",l3),v4=load_index(tmp/"v4",l4);if(v3.first.size()!=v4.first.size()||l3.k!=l4.k||l3.sketch_scale!=l4.sketch_scale||l3.query_sketch_size!=l4.query_sketch_size||l3.max_posting!=l4.max_posting)throw std::runtime_error("v3/v4 metadata equivalence failed");std::set<std::uint64_t>hashes;for(size_t i=0;i<seqs.size();++i){auto g=read_fasta(tmp/("r"+std::to_string(i)+".fna"));auto s=make_scaled_sketch(g,ic.k,ic.sketch_scale,0);hashes.insert(s.begin(),s.end());}for(auto h:hashes){auto x=v3.second.find(h),y=v4.second.find(h);if(x.count!=y.count)throw std::runtime_error("v3/v4 posting count differs");for(std::uint32_t j=0;j<x.count;++j)if(x.at(j)!=y.at(j))throw std::runtime_error("v3/v4 posting IDs differ");}SearchScratch sc3,sc4;std::ostringstream t3,t4;header(t3);header(t4);t3<<search_one(tmp/"x.fna",v3.first,v3.second,l3,sc3);t4<<search_one(tmp/"x.fna",v4.first,v4.second,l4,sc4);if(t3.str()!=t4.str())throw std::runtime_error("v3/v4 full search TSV is not byte-exact");
    for(const auto&name:{"SKETCHES.bin","POSTINGS.bin"})legacy_bytes+=fs::file_size(tmp/"v3"/name);
    for(const auto&name:{"SKETCHES.bin","POSTINGS.bin","POSTINGS_LOOKUP.bin"})compact_bytes+=fs::file_size(tmp/"v4"/name);
    if(compact_bytes>=legacy_bytes)throw std::runtime_error("compact toy index is not smaller");
    fs::copy(tmp/"v4",tmp/"mutable",fs::copy_options::recursive);{Config stable;auto snapshot=load_index(tmp/"mutable",stable);bool protected_or_rejected=false;{std::ofstream f(tmp/"mutable"/"POSTINGS_LOOKUP.bin",std::ios::binary|std::ios::app);if(!f)protected_or_rejected=true;else{f.put('\0');f.close();try{snapshot.second.verify_snapshot();}catch(const std::exception&){protected_or_rejected=true;}}}if(!protected_or_rejected)throw std::runtime_error("mapped index mutation was not blocked/rejected");}
    {std::ofstream f(tmp/"triangle.list",std::ios::binary);for(size_t i=0;i<seqs.size();++i)f<<(tmp/("r"+std::to_string(i)+".fna")).string()<<"\n";}Config tc;tc.list=tmp/"triangle.list";tc.out=tmp/"triangle.tsv";tc.stats_out=tmp/"triangle.stats.json";tc.threads=2;tc.min_af=0;triangle(tc);auto st=read_file_bytes(tc.stats_out);if(json_string(st,"schema")!="gtdb-ani-af-triangle-stats-v1"||json_string(st,"status")!="PASS"||json_string(st,"triangle_mode")!="exact"||json_uint(st,"nodes")!=3||json_uint(st,"pairs_expected")!=3||json_uint(st,"pairs_evaluated")!=3||json_uint(st,"rows_emitted")!=3||json_string(st,"input_list_sha256")!=sha_file(tc.list)||json_string(st,"input_content_set_sha256")!=triangle_input_content_set_sha256(read_paths(tc.list,false))||json_string(st,"output_sha256")!=sha_file(tc.out))throw std::runtime_error("triangle stats/receipt self-test failed");
    fs::copy(tmp/"v4",tmp/"corrupt",fs::copy_options::recursive);{std::fstream f(tmp/"corrupt"/"POSTINGS_LOOKUP.bin",std::ios::binary|std::ios::in|std::ios::out);f.seekg(24);char ch=0;f.get(ch);f.seekp(24);ch^=1;f.put(ch);}bool rejected=false;try{Config bad;auto ignored=load_index(tmp/"corrupt",bad);(void)ignored;}catch(const std::exception&){rejected=true;}if(!rejected)throw std::runtime_error("corrupt compact index was not rejected");fs::remove_all(tmp);
  }catch(...){std::error_code ec;fs::remove_all(tmp,ec);throw;}
  if(a.ani<97||a.af_query<90||a.af_ref<20)throw std::runtime_error("mapping thresholds failed");
  std::cout<<"PASS gtdb-ani-af self-test ANI="<<a.ani<<" AFq="<<a.af_query<<" AFr="<<a.af_ref<<" subset1k_hits="<<hit1<<"/"<<s1.size()<<" subset10k_hits="<<hit10<<"/"<<s10.size()<<" long_reference_bp="<<longref.callable<<" traceback=equivalent reusable_reference_index=equivalent sketch_candidate_estimate=byte-exact v3_v4_postings=equivalent v3_v4_tsv=byte-exact conversion=no-reference-rebuild compact_ratio="<<std::fixed<<std::setprecision(4)<<double(compact_bytes)/double(legacy_bytes)<<" snapshot_mutation=fail-closed corruption=fail-closed\n";
}

static void triangle(const Config&c){
  if(c.list.empty()||(!c.triangle_candidates_only&&(c.out.empty()||c.stats_out.empty())))throw std::runtime_error("triangle requires --list, --out and --stats-out (unless --triangle-candidates-only is used)");
  if(!c.triangle_candidates_only&&fs::absolute(c.out).lexically_normal()==fs::absolute(c.stats_out).lexically_normal())throw std::runtime_error("--out and --stats-out must be distinct write-once paths");
  const auto input_list_sha=sha_file(c.list);auto paths=read_paths(c.list,false);if(sha_file(c.list)!=input_list_sha)throw std::runtime_error("triangle input list changed while being parsed");const auto pairs_expected=exact_pair_count(paths.size());std::vector<Genome>genomes(paths.size());std::vector<ReverseContigs>reverse_genomes(paths.size());std::atomic<size_t>load_next{0};std::exception_ptr error;std::mutex em;
  auto loader=[&](){for(;;){auto i=load_next.fetch_add(1);if(i>=paths.size())return;try{genomes[i]=read_fasta(paths[i]);if(!c.triangle_candidates_only)reverse_genomes[i]=make_reverse_contigs(genomes[i]);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();return;}}};
  std::vector<std::thread>workers;auto nw=std::min<size_t>(size_t(c.threads),paths.size());for(size_t i=0;i<nw;++i)workers.emplace_back(loader);for(auto&t:workers)t.join();if(error)std::rethrow_exception(error);const auto input_content_set_sha=triangle_loaded_content_set_sha256(genomes);
  TriangleCandidatePlan candidate_plan;if(c.triangle_fast){candidate_plan=build_triangle_candidate_plan(genomes,c);write_triangle_candidates(c.triangle_candidates_out,candidate_plan,genomes,c);std::cerr<<"triangle sketch candidates: nodes="<<genomes.size()<<" memberships="<<candidate_plan.sketch_memberships<<" undirected_pairs="<<candidate_plan.pairs.size()<<" skipped_high_frequency_lookups="<<candidate_plan.skipped_high_frequency_lookups<<" min_shared="<<c.triangle_min_shared<<" max_candidates_per_genome="<<c.triangle_max_candidates<<" sketch_size_cap="<<c.triangle_sketch_size<<" max_posting="<<c.triangle_max_posting<<"\n";if(c.triangle_candidates_only)return;}
  std::vector<std::string>lines(paths.size());std::atomic<size_t>pair_next{0};std::atomic<std::uint64_t>pairs_evaluated{0},rows_emitted{0};error=nullptr;workers.clear();
  auto comparer=[&](){for(;;){auto i=pair_next.fetch_add(1);if(i>=paths.size())return;try{const bool has_pairs=c.triangle_fast?!candidate_plan.by_reference[i].empty():i+1<paths.size();if(!has_pairs){lines[i].clear();continue;}auto reference_index=reference_seeds(genomes[i],c);std::ostringstream out;auto compare_one=[&](std::size_t j){auto e=estimate_with_index(genomes[j],genomes[i],reference_index,reverse_genomes[j],c);pairs_evaluated.fetch_add(1,std::memory_order_relaxed);if(std::max(e.af_ref,e.af_query)+1e-12>=c.min_af){row(out,genomes[i].path,genomes[j].path,e);rows_emitted.fetch_add(1,std::memory_order_relaxed);}};if(c.triangle_fast){for(auto j:candidate_plan.by_reference[i])compare_one(j);}else for(size_t j=i+1;j<paths.size();++j)compare_one(j);lines[i]=out.str();}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();return;}}};
  struct BlockOcc {std::uint16_t reference_slot=0;std::uint32_t contig=0,pos=0;};
  using BlockSeedIndex=std::unordered_map<std::uint64_t,std::vector<BlockOcc>>;
  const std::size_t reference_block_size=std::size_t(c.triangle_reference_block_size);
  const std::size_t block_count=(paths.size()+reference_block_size-1)/reference_block_size;
  const bool use_reference_blocks=!c.triangle_fast&&reference_block_size>1;
  if(c.triangle_global_seed_join){
    const auto index_started=std::chrono::steady_clock::now();auto global_index=build_global_seed_index(genomes,c,nw);const auto index_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-index_started).count();std::cerr<<"triangle global exact seed index: occurrences="<<global_index.occurrences.size()<<" prefix_bits="<<global_index.prefix_bits<<" offset_entries="<<global_index.offsets.size()<<" build_seconds="<<std::fixed<<std::setprecision(3)<<index_seconds<<"\n";
    struct EmittedRow {std::uint32_t reference=0;std::string text;};
    std::vector<std::vector<EmittedRow>>query_rows(paths.size());std::vector<std::size_t>query_order;query_order.reserve(paths.size());for(std::size_t j=1;j<paths.size();++j)query_order.push_back(j);
    std::sort(query_order.begin(),query_order.end(),[&](std::size_t a,std::size_t b){if(genomes[a].callable!=genomes[b].callable)return genomes[a].callable>genomes[b].callable;return a<b;});
    std::atomic<std::size_t>query_next{0};workers.clear();
    auto global_query_worker=[&](){for(;;){const auto task=query_next.fetch_add(1);if(task>=query_order.size())return;const auto j=query_order[task];try{
      const auto&q=genomes[j];std::vector<std::unique_ptr<VoteMap>>votes(j);
      for(std::size_t qi=0;qi<q.contigs.size();++qi){const auto&seq=q.contigs[qi].seq;each_kmer(seq,c.k,[&](std::uint32_t qp,std::uint64_t fw,std::uint64_t rv){
        auto add_votes=[&](std::uint64_t seed,std::uint32_t query_pos,bool reverse){const auto range=global_index.equal_range(seed);for(auto occurrence=range.first;occurrence!=range.second;++occurrence){if(occurrence->reference>=j)break;auto&target=votes[occurrence->reference];if(!target){const auto keyspace=std::uint64_t(1)<<(2*c.k);const auto random_hits=2.0L*static_cast<long double>(q.callable)*static_cast<long double>(genomes[occurrence->reference].callable)/static_cast<long double>(keyspace);const auto expected=std::size_t(std::max<long double>(16.0L,std::min<long double>(256.0L,random_hits)));target=std::make_unique<VoteMap>(expected);}auto&v=(*target)[{std::uint32_t(qi),occurrence->contig,int(occurrence->pos)-int(query_pos),reverse}];++v.count;v.first=std::min(v.first,query_pos);v.last=std::max(v.last,query_pos);}};
        add_votes(fw,qp,false);const auto reverse_pos=std::uint32_t(seq.size()-std::size_t(c.k)-std::size_t(qp));add_votes(rv,reverse_pos,true);
      });}
      const VoteMap empty_votes;std::vector<EmittedRow>emitted;for(std::size_t reference=0;reference<j;++reference){const auto&pair_votes=votes[reference]?*votes[reference]:empty_votes;auto e=estimate_from_votes(q,genomes[reference],reverse_genomes[j],c,pair_votes);if(std::max(e.af_ref,e.af_query)+1e-12>=c.min_af){std::ostringstream output;row(output,genomes[reference].path,genomes[j].path,e);emitted.push_back({std::uint32_t(reference),output.str()});rows_emitted.fetch_add(1,std::memory_order_relaxed);}}
      query_rows[j]=std::move(emitted);pairs_evaluated.fetch_add(j,std::memory_order_relaxed);
    }catch(...){std::lock_guard<std::mutex>guard(em);if(!error)error=std::current_exception();return;}}};
    const auto query_started=std::chrono::steady_clock::now();const auto query_workers=std::min<std::size_t>(nw,query_order.size());for(std::size_t i=0;i<query_workers;++i)workers.emplace_back(global_query_worker);for(auto&t:workers)t.join();if(error)std::rethrow_exception(error);const auto query_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-query_started).count();std::cerr<<"triangle global exact query/refine seconds="<<std::fixed<<std::setprecision(3)<<query_seconds<<"\n";
    for(std::size_t j=1;j<query_rows.size();++j)for(auto&emitted:query_rows[j])lines[emitted.reference]+=emitted.text;
  }else if(use_reference_blocks){
    // Keep several immutable reference blocks resident and schedule all of
    // their (block, query) units through one shared queue.  This overlaps the
    // long-tail genomes of different blocks while each task still decodes a
    // query k-mer stream only once for every reference in its block.
    struct BlockWork {std::size_t begin=0,count=0;BlockSeedIndex index;std::vector<std::vector<std::string>>pair_rows;};
    struct BlockQueryTask {std::size_t block=0,query=0;};
    constexpr std::size_t resident_blocks=8;
    for(std::size_t batch_begin=0;batch_begin<block_count;batch_begin+=resident_blocks){
      const auto batch_end=std::min(block_count,batch_begin+resident_blocks);std::vector<std::unique_ptr<BlockWork>>batch;batch.reserve(batch_end-batch_begin);
      for(std::size_t block_id=batch_begin;block_id<batch_end;++block_id){
        auto work=std::make_unique<BlockWork>();work->begin=block_id*reference_block_size;const auto end=std::min(paths.size(),work->begin+reference_block_size);work->count=end-work->begin;
        for(std::size_t slot=0;slot<work->count;++slot){auto index=reference_seeds(genomes[work->begin+slot],c);work->index.reserve(work->index.size()+index.size());for(const auto&[seed,occurrences]:index){auto&target=work->index[seed];target.reserve(target.size()+occurrences.size());for(const auto&o:occurrences)target.push_back({std::uint16_t(slot),o.contig,o.pos});}}
        work->pair_rows.assign(work->count,std::vector<std::string>(paths.size()));batch.push_back(std::move(work));
      }
      std::vector<BlockQueryTask>tasks;for(std::size_t b=0;b<batch.size();++b)for(std::size_t j=batch[b]->begin+1;j<paths.size();++j)tasks.push_back({b,j});
      std::sort(tasks.begin(),tasks.end(),[&](const BlockQueryTask&a,const BlockQueryTask&b){if(genomes[a.query].callable!=genomes[b.query].callable)return genomes[a.query].callable>genomes[b.query].callable;if(batch[a.block]->begin!=batch[b.block]->begin)return batch[a.block]->begin<batch[b.block]->begin;return a.query<b.query;});
      std::atomic<std::size_t>task_next{0};std::exception_ptr batch_error;std::mutex batch_error_mutex;workers.clear();
      auto query_worker=[&](){for(;;){const auto task_index=task_next.fetch_add(1);if(task_index>=tasks.size())return;const auto task=tasks[task_index];auto&work=*batch[task.block];const auto j=task.query;try{
        const auto active=std::min(work.count,j-work.begin);std::vector<VoteMap>votes(active);const auto&q=genomes[j];
        for(std::size_t qi=0;qi<q.contigs.size();++qi){const auto&seq=q.contigs[qi].seq;each_kmer(seq,c.k,[&](std::uint32_t qp,std::uint64_t fw,std::uint64_t rv){
          auto add_votes=[&](std::uint64_t seed,std::uint32_t query_pos,bool reverse){auto found=work.index.find(seed);if(found==work.index.end())return;for(const auto&o:found->second){if(o.reference_slot>=active)continue;auto&v=votes[o.reference_slot][{std::uint32_t(qi),o.contig,int(o.pos)-int(query_pos),reverse}];++v.count;v.first=std::min(v.first,query_pos);v.last=std::max(v.last,query_pos);}};
          add_votes(fw,qp,false);const auto reverse_pos=std::uint32_t(seq.size()-std::size_t(c.k)-std::size_t(qp));add_votes(rv,reverse_pos,true);
        });}
        pairs_evaluated.fetch_add(active,std::memory_order_relaxed);
        for(std::size_t slot=0;slot<active;++slot){auto e=estimate_from_votes(q,genomes[work.begin+slot],reverse_genomes[j],c,votes[slot]);if(std::max(e.af_ref,e.af_query)+1e-12>=c.min_af){std::ostringstream output;row(output,genomes[work.begin+slot].path,genomes[j].path,e);work.pair_rows[slot][j]=output.str();rows_emitted.fetch_add(1,std::memory_order_relaxed);}}
      }catch(...){std::lock_guard<std::mutex>g(batch_error_mutex);if(!batch_error)batch_error=std::current_exception();return;}}};
      const auto batch_workers=std::min<std::size_t>(nw,tasks.size());for(std::size_t w=0;w<batch_workers;++w)workers.emplace_back(query_worker);for(auto&t:workers)t.join();if(batch_error)std::rethrow_exception(batch_error);
      for(const auto&owned:batch){const auto&work=*owned;for(std::size_t slot=0;slot<work.count;++slot){std::ostringstream output;for(std::size_t j=work.begin+slot+1;j<paths.size();++j)output<<work.pair_rows[slot][j];lines[work.begin+slot]=output.str();}}
    }
  }else{
    for(size_t i=0;i<nw;++i)workers.emplace_back(comparer);
    for(auto&t:workers)t.join();
  }
  if(error)std::rethrow_exception(error);
  const auto evaluated=pairs_evaluated.load(),emitted=rows_emitted.load();if(sha_file(c.list)!=input_list_sha)throw std::runtime_error("triangle input list changed during execution");
  if(!c.triangle_fast&&evaluated!=pairs_expected)throw std::runtime_error("exact triangle did not evaluate N*(N-1)/2 pairs");
  if(emitted>evaluated)throw std::runtime_error("triangle emitted more rows than evaluated pairs");
  std::ostringstream out;header(out);for(const auto&s:lines)out<<s;write_once(c.out,out.str());
  const auto output_sha=sha_file(c.out);std::ostringstream stats;stats
    <<"{\n  \"schema\": \"gtdb-ani-af-triangle-stats-v1\",\n  \"status\": \"PASS\",\n"
    <<"  \"triangle_mode\": \""<<(c.triangle_fast?"sketch":"exact")<<"\",\n"
    <<"  \"edge_mode\": \"sparse\",\n  \"nodes\": "<<paths.size()<<",\n"
    <<"  \"pairs_expected\": "<<pairs_expected<<",\n  \"pairs_evaluated\": "<<evaluated<<",\n"
    <<"  \"rows_emitted\": "<<emitted<<",\n  \"min_af\": "<<std::setprecision(17)<<c.min_af<<",\n"
    <<"  \"exact_reference_block_size\": "<<(c.triangle_fast?1:c.triangle_reference_block_size)<<",\n"
    <<"  \"exact_execution\": \""<<(c.triangle_fast?"sketch-candidate":(c.triangle_global_seed_join?"global-seed-join":(use_reference_blocks?"fused-reference-blocks":"reference-at-a-time")))<<"\",\n"
    <<"  \"input_list_path\": \""<<json_escape(fs::weakly_canonical(c.list).string())<<"\",\n"
    <<"  \"input_list_sha256\": \""<<input_list_sha<<"\",\n"
    <<"  \"input_content_set_sha256\": \""<<input_content_set_sha<<"\",\n"
    <<"  \"output_path\": \""<<json_escape(fs::weakly_canonical(c.out).string())<<"\",\n"
    <<"  \"output_sha256\": \""<<output_sha<<"\"\n}\n";
  write_once(c.stats_out,stats.str());
}

static void index_info(const Config&request){if(request.index.empty())throw std::runtime_error("index-info requires --index");auto start=std::chrono::steady_clock::now();Config c=request;auto loaded=load_index(c.index,c);auto&pi=loaded.second;std::ifstream cf(c.index/"COMPLETE.json");std::ostringstream cs;cs<<cf.rdbuf();auto schema=json_string(cs.str(),"schema");std::uint64_t disk=fs::file_size(c.index/"COMPLETE.json")+fs::file_size(c.index/"REFS.tsv")+fs::file_size(c.index/"SKETCHES.bin")+fs::file_size(c.index/"POSTINGS.bin");if(fs::exists(c.index/"POSTINGS_LOOKUP.bin"))disk+=fs::file_size(c.index/"POSTINGS_LOOKUP.bin");auto heap=pi.legacy_hashes.capacity()*sizeof(std::uint64_t)+pi.legacy_count_offsets.capacity()*sizeof(std::uint64_t);auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"schema\treferences\tposting_groups\tkept_postings\tartifact_bytes\tpersistent_posting_heap_bytes\tload_validate_seconds\n"<<schema<<'\t'<<loaded.first.size()<<'\t'<<pi.groups<<'\t'<<pi.total_refs<<'\t'<<disk<<'\t'<<heap<<'\t'<<std::fixed<<std::setprecision(6)<<elapsed<<'\n';}

static void help(){std::cout<<"gtdb-ani-af: experimental independent ANI/AF estimator (not skani; not production)\n"
" pair --query Q --reference R [--out TSV]\n"
" index --manifest REFS --taxonomy TAX --radii RADII --release-metadata RELEASE --out-dir DB [--index-format compact-v4|legacy-v3 --threads N --sketch-scale 64 --query-sketch-size 1000 --max-posting 1000]\n"
" search --index DB --queries LIST --out TSV [--threads N --min-candidate-hits 2 --max-candidates 20 --top 5 --deep-verify]\n"
" compact-index --index LEGACY_V3 --out-dir COMPACT_V4\n"
" rebind-taxonomy --index PASS_DB --taxonomy MATCHING_TAXONOMY --out-dir NEW_DB\n"
" index-info --index DB\n"
 " triangle --list FASTA_LIST --out TSV --stats-out JSON [--threads N --profile medium --edge-mode sparse --min-af 15]\n"
 "   default: --triangle-mode exact (all pairs; unchanged reference semantics)\n"
 "            [--triangle-reference-block-size 1]  # experimental exact fused lookup, 1..64\n"
 "            [--triangle-global-seed-join] [--triangle-global-prefix-bits 20] # exact one-pass query dispatch; k<=15\n"
 "   fast: --triangle-mode sketch [--sketch-scale 64 --triangle-min-shared 1 --triangle-max-candidates 128]\n"
 "         [--triangle-sketch-size 0 --triangle-max-posting 1000 --triangle-candidates-out TSV]\n"
 "         [--triangle-candidates-only]  # build/audit candidates without ANI/AF evaluation\n"
 "   sketch-size/max-candidates/max-posting 0 mean unlimited; per-genome top-K is symmetrized by union\n"
" self-test\n";}

int main(int argc,char**argv)try{
  if(argc<2){help();return 2;}std::string mode=argv[1];if(mode=="help"||mode=="--help"||mode=="-h"){help();return 0;}if(mode=="self-test"||mode=="--self-test"){self_test();return 0;}Config c=parse(argc,argv,2);
  if(mode=="pair"){if(c.query.empty()||c.reference.empty())throw std::runtime_error("pair requires --query/--reference");auto q=read_fasta(c.query),r=read_fasta(c.reference);auto e=estimate(q,r,c);std::ostringstream o;header(o);row(o,r.path,q.path,e);write_once(c.out,o.str());return 0;}
  if(mode=="triangle"){triangle(c);return 0;}
  if(mode=="index"){if(c.manifest.empty()||c.taxonomy.empty()||c.radii.empty()||c.release_metadata.empty()||c.out_dir.empty())throw std::runtime_error("index requires --manifest --taxonomy --radii --release-metadata --out-dir");build_index(c);return 0;}
  if(mode=="compact-index"){compact_legacy_index(c);return 0;}
  if(mode=="rebind-taxonomy"){rebind_taxonomy(c);return 0;}
  if(mode=="index-info"){index_info(c);return 0;}
  if(mode=="search"){
    if(c.index.empty()||c.queries.empty()||c.out.empty())throw std::runtime_error("search requires --index --queries --out");
    auto loaded=load_index(c.index,c);auto&refs=loaded.first;auto&postings=loaded.second;auto queries=read_paths(c.queries,false);std::vector<std::string>outputs(queries.size());std::atomic<size_t>next{0};std::exception_ptr error;std::mutex em;
    auto worker=[&](){SearchScratch scratch;for(;;){size_t qi=next.fetch_add(1);if(qi>=queries.size())return;try{outputs[qi]=search_one(queries[qi],refs,postings,c,scratch);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();return;}}};
    std::vector<std::thread>ts;for(int i=0;i<std::min<int>(c.threads,queries.size());++i)ts.emplace_back(worker);for(auto&t:ts)t.join();if(error)std::rethrow_exception(error);std::ostringstream o;header(o);for(const auto&s:outputs)o<<s;write_once(c.out,o.str());return 0;
  }
  throw std::runtime_error("unknown mode: "+mode);
}catch(const std::exception&e){std::cerr<<"fatal: "<<e.what()<<"\n";return 2;}
