#include <algorithm>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef _WIN32
#include <csignal>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

// Stage 3B uses this repository's experimental gtdb-ani-af executable only
// as the direct SAG-vs-SAG ANI/AF engine.  It deliberately does not search a
// GTDB reference database: every exact Dna2bit-negative SAG that passes the
// two historical quality gates enters the all-pairs triangle.  GTDB-Tk's
// bac120 marker map remains an independent source of conflict evidence.
// BLAST+, leidenalg and cpp-subass-full remain external implementations.

struct Config {
  fs::path manifest, quality, marker_map, out;
  fs::path ani_engine="gtdb-ani-af", makeblastdb="makeblastdb", blastn="blastn";
  fs::path python="python3", leiden_backend, subass="cpp-subass";
  fs::path flye_root="/home/data/fyc/biosoft/miniconda3/envs/assemble";
  int threads=16;
  bool dry=false, resume=false, self_test=false, allow_experimental_ani=false;
  std::vector<fs::path> flye_provenance;
};

struct FastaStats { std::uint64_t total=0, max_contig=0, gc=0, acgt=0; };
struct Sag {
  std::string id, reason;
  fs::path assembly;
  double declared_gc=std::numeric_limits<double>::quiet_NaN();
  double contamination=std::numeric_limits<double>::quiet_NaN();
  std::uint64_t declared_total=0, declared_max=0;
  FastaStats actual;
};
struct MarkerAgg {
  double sum=0.0, min=101.0;
  std::uint32_t count=0;
  void add(double pid){sum+=pid;min=std::min(min,pid);++count;}
};

struct RuntimeMetrics {
  std::uint64_t tsv_rows=0, sha_cache_hits=0, sha_cache_misses=0;
  std::uint64_t sha_bytes_read=0, fasta_bytes_read=0;
  std::size_t marker_peak_jobs=0, subass_peak_jobs=0;
};
static RuntimeMetrics metrics;

static std::string lower(std::string s) {
  for(char& c:s)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
static std::string trim(std::string s) {
  std::size_t first=0,last=s.size();
  while(first<last&&std::isspace(static_cast<unsigned char>(s[first])))++first;
  while(last>first&&std::isspace(static_cast<unsigned char>(s[last-1])))--last;
  if(first==0&&last==s.size())return s;
  return s.substr(first,last-first);
}
static std::string strip_bom(std::string s){if(s.size()>=3&&static_cast<unsigned char>(s[0])==0xef&&static_cast<unsigned char>(s[1])==0xbb&&static_cast<unsigned char>(s[2])==0xbf)s.erase(0,3);return s;}
static std::vector<std::string> split(const std::string&s,char d) {
  std::vector<std::string> v; std::stringstream ss(s); std::string x;
  while(std::getline(ss,x,d))v.push_back(x);
  if(!s.empty()&&s.back()==d)v.emplace_back();
  return v;
}

class TsvReader {
  fs::path path_;
  std::ifstream in_;
  std::vector<std::string> header_;
  std::size_t line_number_=1;
public:
  explicit TsvReader(const fs::path&p):path_(p),in_(p){
    if(!in_)throw std::runtime_error("cannot read TSV: "+p.string());
    std::string line;if(!std::getline(in_,line))throw std::runtime_error("empty TSV: "+p.string());
    line=strip_bom(line);header_=split(line,'\t');for(auto&x:header_)x=trim(std::move(x));
    std::set<std::string> seen;for(const auto&x:header_)if(x.empty()||!seen.insert(lower(x)).second)throw std::runtime_error("empty/duplicate TSV header in "+p.string());
  }
  const std::vector<std::string>& header()const{return header_;}
  std::size_t column(std::initializer_list<const char*>names,bool required=true)const{
    for(const auto*n:names)for(std::size_t i=0;i<header_.size();++i)if(lower(header_[i])==lower(n))return i;
    if(!required)return std::numeric_limits<std::size_t>::max();
    std::string msg="missing required column (one of:";for(const auto*n:names)msg+=" "+std::string(n);throw std::runtime_error(msg+") in "+path_.string());
  }
  bool next(std::vector<std::string>&row){
    std::string line;
    while(std::getline(in_,line)){
      ++line_number_;if(line.empty())continue;
      row=split(line,'\t');if(row.size()!=header_.size())throw std::runtime_error(path_.string()+":"+std::to_string(line_number_)+" wrong field count");
      for(auto&x:row)x=trim(std::move(x));
      ++metrics.tsv_rows;return true;
    }
    return false;
  }
};
static std::string safe_id(const std::string& raw) {
  if(raw.empty()||raw=="."||raw=="..")throw std::runtime_error("empty/unsafe identifier");
  std::string s=raw;
  for(char& c:s)if(!(std::isalnum(static_cast<unsigned char>(c))||c=='_'||c=='-'||c=='.'))c='_';
  if(s!=raw)throw std::runtime_error("identifier contains unsafe characters: "+raw);
  return s;
}
static std::string q(const fs::path&p) {
  std::string r="'"; for(char c:p.string())r+=c=='\''?"'\\''":std::string(1,c); return r+"'";
}
static std::string json(const std::string&s) {
  std::string r; for(unsigned char c:s){
    if(c=='\\')r+="\\\\";else if(c=='\"')r+="\\\"";else if(c=='\n')r+="\\n";
    else if(c=='\r')r+="\\r";else if(c=='\t')r+="\\t";else if(c<32){std::ostringstream x;x<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c);r+=x.str();}else r+=char(c);
  }return r;
}

class Sha256 {
  std::uint32_t state_[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                           0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  std::uint8_t block_[64]={};std::size_t used_=0;std::uint64_t bits_=0;
  static std::uint32_t rr(std::uint32_t x,int n){return (x>>n)|(x<<(32-n));}
  void transform(){
    static const std::uint32_t k[64]={
      0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
      0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
      0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
      0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
      0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
      0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
      0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
      0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::uint32_t w[64];for(int i=0;i<16;++i)w[i]=(std::uint32_t(block_[4*i])<<24)|(std::uint32_t(block_[4*i+1])<<16)|(std::uint32_t(block_[4*i+2])<<8)|block_[4*i+3];
    for(int i=16;i<64;++i){auto s0=rr(w[i-15],7)^rr(w[i-15],18)^(w[i-15]>>3);auto s1=rr(w[i-2],17)^rr(w[i-2],19)^(w[i-2]>>10);w[i]=w[i-16]+s0+w[i-7]+s1;}
    auto a=state_[0],b=state_[1],c=state_[2],d=state_[3],e=state_[4],f=state_[5],g=state_[6],h=state_[7];
    for(int i=0;i<64;++i){auto s1=rr(e,6)^rr(e,11)^rr(e,25);auto ch=(e&f)^((~e)&g);auto t1=h+s1+ch+k[i]+w[i];auto s0=rr(a,2)^rr(a,13)^rr(a,22);auto maj=(a&b)^(a&c)^(b&c);auto t2=s0+maj;h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;}
    state_[0]+=a;state_[1]+=b;state_[2]+=c;state_[3]+=d;state_[4]+=e;state_[5]+=f;state_[6]+=g;state_[7]+=h;
  }
public:
  void update(const void*data,std::size_t n){auto*p=static_cast<const std::uint8_t*>(data);bits_+=std::uint64_t(n)*8;while(n){auto take=std::min(n,64-used_);std::copy(p,p+take,block_+used_);used_+=take;p+=take;n-=take;if(used_==64){transform();used_=0;}}}
  void update(const std::string&s){update(s.data(),s.size());}
  std::string finish(){auto original_bits=bits_;std::uint8_t one=0x80;update(&one,1);std::uint8_t zero=0;while(used_!=56)update(&zero,1);std::uint8_t len[8];for(int i=0;i<8;++i)len[7-i]=std::uint8_t(original_bits>>(8*i));update(len,8);std::ostringstream o;o<<std::hex<<std::setfill('0');for(auto x:state_)o<<std::setw(8)<<x;return o.str();}
};
static std::string sha256_text(const std::string&s){Sha256 h;h.update(s);return h.finish();}

struct FileIdentity {
  std::uintmax_t size=0;
  fs::file_time_type mtime{};
#ifndef _WIN32
  dev_t device=0;ino_t inode=0;std::int64_t mtime_s=0,mtime_ns=0,ctime_s=0,ctime_ns=0;
#endif
};
static FileIdentity file_identity(const fs::path&p){
  if(!fs::is_regular_file(p))throw std::runtime_error("fingerprint input missing: "+p.string());
  FileIdentity x;x.size=fs::file_size(p);x.mtime=fs::last_write_time(p);
#ifndef _WIN32
  struct stat st{};if(::stat(p.c_str(),&st)!=0)throw std::runtime_error("cannot stat file: "+p.string());
  x.device=st.st_dev;x.inode=st.st_ino;x.mtime_s=st.st_mtim.tv_sec;x.mtime_ns=st.st_mtim.tv_nsec;x.ctime_s=st.st_ctim.tv_sec;x.ctime_ns=st.st_ctim.tv_nsec;
#endif
  return x;
}
static bool same_identity(const FileIdentity&a,const FileIdentity&b){
  if(a.size!=b.size||a.mtime!=b.mtime)return false;
#ifndef _WIN32
  return a.device==b.device&&a.inode==b.inode&&a.mtime_s==b.mtime_s&&a.mtime_ns==b.mtime_ns&&a.ctime_s==b.ctime_s&&a.ctime_ns==b.ctime_ns;
#else
  return true;
#endif
}
struct CachedDigest {FileIdentity identity;std::string sha256;};
static std::unordered_map<std::string,CachedDigest> digest_cache;
static std::string canonical_key(const fs::path&p){return fs::weakly_canonical(p).string();}
static void cache_digest(const fs::path&p,const FileIdentity&before,const FileIdentity&after,const std::string&sha){
  if(!same_identity(before,after))throw std::runtime_error("file changed while being read: "+p.string());
  digest_cache[canonical_key(p)]={after,sha};
}
static std::string sha256_file(const fs::path&p){
  auto key=canonical_key(p);auto before=file_identity(p);auto it=digest_cache.find(key);
  if(it!=digest_cache.end()){
    if(!same_identity(before,it->second.identity))throw std::runtime_error("file changed after its validated fingerprint (fail-closed): "+p.string());
    ++metrics.sha_cache_hits;return it->second.sha256;
  }
  ++metrics.sha_cache_misses;std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot hash file: "+p.string());Sha256 h;char b[1<<16];
  while(f){f.read(b,sizeof(b));auto n=f.gcount();if(n>0){h.update(b,static_cast<std::size_t>(n));metrics.sha_bytes_read+=static_cast<std::uint64_t>(n);}}
  auto sha=h.finish();auto after=file_identity(p);cache_digest(p,before,after,sha);return sha;
}
struct FingerprintedFile {fs::path canonical;std::uintmax_t bytes=0;std::string sha256;};
static std::pair<std::string,std::vector<FingerprintedFile>> fingerprint_files(const std::vector<fs::path>&paths,bool require_nonempty=false){
  Sha256 set_hash;std::vector<FingerprintedFile> files;files.reserve(paths.size());
  for(const auto&p:paths){auto id=file_identity(p);if(require_nonempty&&id.size==0)throw std::runtime_error("fingerprint input empty: "+p.string());FingerprintedFile x{fs::weakly_canonical(p),id.size,sha256_file(p)};std::string row=x.canonical.string()+"\t"+std::to_string(x.bytes)+"\t"+x.sha256+"\n";set_hash.update(row);files.push_back(std::move(x));}
  return {set_hash.finish(),std::move(files)};
}
static std::string file_set_fingerprint(const std::vector<fs::path>&paths,bool require_nonempty=false){return fingerprint_files(paths,require_nonempty).first;}
static std::string triangle_input_content_set_sha256(const std::vector<fs::path>&paths){
  Sha256 h;h.update("gtdb-ani-af-triangle-input-content-set-v1\n");
  for(const auto&p:paths){auto canonical=fs::weakly_canonical(p);h.update(canonical.string());h.update("\t");h.update(sha256_file(canonical));h.update("\n");}
  return h.finish();
}
static std::string receipt_field(const std::string&s,const std::string&key){auto needle="\""+key+"\": \"";auto p=s.find(needle);if(p==std::string::npos)throw std::runtime_error("receipt missing field: "+key);p+=needle.size();auto e=s.find('"',p);if(e==std::string::npos)throw std::runtime_error("malformed receipt field: "+key);return s.substr(p,e-p);}
static std::uint64_t receipt_uint(const std::string&s,const std::string&key){auto needle="\""+key+"\": ";auto p=s.find(needle);if(p==std::string::npos)throw std::runtime_error("receipt missing numeric field: "+key);p+=needle.size();auto e=p;while(e<s.size()&&std::isdigit(static_cast<unsigned char>(s[e])))++e;if(e==p)throw std::runtime_error("malformed receipt numeric field: "+key);auto token=s.substr(p,e-p);try{std::size_t used=0;auto value=std::stoull(token,&used);if(used!=token.size())throw std::runtime_error("");return value;}catch(...){throw std::runtime_error("invalid receipt numeric field: "+key);}}
static void verify_receipt_fields(const std::string&rs,const std::string&stage_sha,const std::string&command_sha,const std::string&input_sha,const std::string&config_sha,const std::string&contract_sha,const std::string&outputs_sha){if(receipt_field(rs,"schema")!="sag-stage3b-receipt-v2"||receipt_field(rs,"status")!="PASS"||receipt_field(rs,"stage_sha256")!=stage_sha||receipt_field(rs,"command_sha256")!=command_sha||receipt_field(rs,"input_fingerprint_sha256")!=input_sha||receipt_field(rs,"config_sha256")!=config_sha||receipt_field(rs,"contract_sha256")!=contract_sha||receipt_field(rs,"outputs_fingerprint_sha256")!=outputs_sha)throw std::runtime_error("receipt field mismatch");}
static double number(const std::string&s,const std::string&what) {
  try{size_t n=0;double v=std::stod(s,&n);if(n!=s.size()||!std::isfinite(v))throw std::runtime_error("");return v;}
  catch(...){throw std::runtime_error("invalid numeric "+what+": "+s);}
}
static std::uint64_t integer(const std::string&s,const std::string&what) {
  double v=number(s,what);if(v<0||std::floor(v)!=v)throw std::runtime_error("invalid integer "+what+": "+s);return static_cast<std::uint64_t>(v);
}
static bool dna2bit_negative_reason(const std::string&s){static const std::set<std::string> allowed={"dna2bit_rejected_or_no_hit","dna2bit_negative","dna2bit_no_hit","dna2bit_rejected"};return allowed.count(lower(trim(s)))!=0;}
static bool stage3b_quality_pass(const Sag&s){return s.actual.max_contig>=1000&&s.contamination<5.0;}
static std::uint64_t exact_pair_count(std::size_t n){
  if(n<2)return 0;
  std::uint64_t a=static_cast<std::uint64_t>(n),b=a-1;
  if(a/2*2==a)a/=2;else b/=2;
  if(b&&a>std::numeric_limits<std::uint64_t>::max()/b)throw std::runtime_error("exact triangle pair count overflows uint64");
  return a*b;
}
static FastaStats fasta_stats(const fs::path&p) {
  auto before=file_identity(p);std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot read FASTA: "+p.string());
  FastaStats z;std::string line;bool seen_header=false;std::uint64_t cur=0;Sha256 hash;
  auto consume=[&](const std::string&value){
    if(!value.empty()&&value[0]=='>'){if(seen_header)z.max_contig=std::max(z.max_contig,cur);seen_header=true;cur=0;return;}
    if(!seen_header&&!trim(value).empty())throw std::runtime_error("FASTA sequence before header: "+p.string());
    for(unsigned char c:value)if(!std::isspace(c)){char u=static_cast<char>(std::toupper(c));
      if(!(u=='A'||u=='C'||u=='G'||u=='T'||u=='N'||u=='R'||u=='Y'||u=='S'||u=='W'||u=='K'||u=='M'||u=='B'||u=='D'||u=='H'||u=='V'||u=='-'||u=='.'))throw std::runtime_error("invalid FASTA character in "+p.string());
      ++z.total;++cur;if(u=='A'||u=='C'||u=='G'||u=='T'){++z.acgt;if(u=='G'||u=='C')++z.gc;}
    }
  };
  char block[1<<16];while(f){f.read(block,sizeof(block));auto n=f.gcount();if(n<=0)continue;hash.update(block,static_cast<std::size_t>(n));metrics.fasta_bytes_read+=static_cast<std::uint64_t>(n);for(std::streamsize i=0;i<n;++i){char ch=block[i];if(ch=='\n'){consume(line);line.clear();}else line.push_back(ch);}}
  if(!line.empty())consume(line);
  if(seen_header)z.max_contig=std::max(z.max_contig,cur);
  if(!seen_header||z.total==0)throw std::runtime_error("missing/empty FASTA (fail-closed): "+p.string());
  auto after=file_identity(p);cache_digest(p,before,after,hash.finish());
  return z;
}
static fs::path canonical_existing(const fs::path&p) {
  if(!fs::is_regular_file(p))throw std::runtime_error("required file not found: "+p.string());
  return fs::weakly_canonical(p);
}
static bool command_is_path(const fs::path&p){return p.has_parent_path()||p.is_absolute();}
static fs::path resolve_tool(const fs::path&p,bool dry){
  if(p.empty())throw std::runtime_error("empty tool path");
  if(command_is_path(p)){auto resolved=canonical_existing(p);return resolved;}
  if(dry)return p;
  const char*raw=std::getenv("PATH");if(!raw)throw std::runtime_error("PATH is unset while resolving tool: "+p.string());
#ifdef _WIN32
  constexpr char separator=';';const std::vector<std::string>suffixes={"",".exe",".cmd",".bat"};
#else
  constexpr char separator=':';const std::vector<std::string>suffixes={""};
#endif
  for(const auto&dir:split(raw,separator)){if(dir.empty())continue;for(const auto&suffix:suffixes){auto candidate=fs::path(dir)/(p.string()+suffix);if(fs::is_regular_file(candidate))return fs::canonical(candidate);}}
  throw std::runtime_error("tool not found on PATH: "+p.string());
}
static std::vector<fs::path> collect_flye_provenance(const fs::path&root){
  if(!fs::is_directory(root))throw std::runtime_error("Flye root not found: "+root.string());
  const auto package=root/"lib/python3.9/site-packages/flye";if(!fs::is_directory(package))throw std::runtime_error("Flye package root not found: "+package.string());
  std::vector<fs::path> files={root/"bin/flye",root/"bin/flye-modules",root/"bin/flye-minimap2",root/"bin/flye-samtools",package/"config/bin_cfg/asm_subasm.cfg"};
  for(const auto&entry:fs::recursive_directory_iterator(package)){if(!entry.is_regular_file())continue;auto p=entry.path();if(p.extension()==".pyc"||p.string().find("__pycache__")!=std::string::npos)continue;files.push_back(p);}
  for(auto&p:files)p=canonical_existing(p);
  std::sort(files.begin(),files.end());files.erase(std::unique(files.begin(),files.end()),files.end());return files;
}
static std::string read_text(const fs::path&p){std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("cannot read: "+p.string());std::ostringstream out;out<<f.rdbuf();if(f.bad())throw std::runtime_error("failed reading: "+p.string());return out.str();}

static void write_once(const fs::path&p,const std::string&content) {
  if(fs::exists(p))throw std::runtime_error("refusing to overwrite: "+p.string());
  fs::create_directories(p.parent_path());fs::path tmp=p;tmp+=".tmp";
  if(fs::exists(tmp))throw std::runtime_error("stale temporary file: "+tmp.string());
  {std::ofstream f(tmp,std::ios::binary);if(!f)throw std::runtime_error("cannot write: "+tmp.string());f<<content;if(!f)throw std::runtime_error("write failed: "+tmp.string());}
  fs::rename(tmp,p);
}
static void write_or_verify(const fs::path&p,const std::string&content,const Config&c) {
  if(!fs::exists(p)){write_once(p,content);return;}
  if(!c.resume)throw std::runtime_error("refusing to overwrite: "+p.string());
  if(!fs::is_regular_file(p)||fs::file_size(p)!=content.size()||sha256_file(p)!=sha256_text(content))
    throw std::runtime_error("resume derived artifact differs from deterministic recomputation: "+p.string());
}
template<class Writer>static void write_generated(const fs::path&p,const Config&c,Writer writer){
  if(fs::exists(p)&&!c.resume)throw std::runtime_error("refusing to overwrite: "+p.string());
  fs::create_directories(p.parent_path());fs::path tmp=p;tmp+=".recomputed.tmp";
  if(fs::exists(tmp))throw std::runtime_error("stale generated-file temporary: "+tmp.string());
  {std::ofstream out(tmp,std::ios::binary);if(!out)throw std::runtime_error("cannot write: "+tmp.string());writer(out);if(!out)throw std::runtime_error("write failed: "+tmp.string());}
  if(fs::exists(p)){
    bool same=fs::is_regular_file(p)&&fs::file_size(p)==fs::file_size(tmp)&&sha256_file(p)==sha256_file(tmp);fs::remove(tmp);
    if(!same)throw std::runtime_error("resume derived artifact differs from deterministic streaming recomputation: "+p.string());
  }else fs::rename(tmp,p);
}
static void verify_outputs(const std::vector<fs::path>&outs,bool allow_empty=false) {
  for(const auto&p:outs)if(!fs::is_regular_file(p)||(!allow_empty&&fs::file_size(p)==0))throw std::runtime_error("stage output missing/empty: "+p.string());
}
static int execute_argv(const std::vector<std::string>&argv,const std::string&fallback_command){
  if(argv.empty())return std::system(fallback_command.c_str());
#ifdef _WIN32
  // The production target is Linux. Keep a buildable Windows fallback for
  // self-tests; all untrusted paths remain shell-quoted in fallback_command.
  return std::system(fallback_command.c_str());
#else
  pid_t pid=::fork();if(pid<0)throw std::runtime_error("fork failed");
  if(pid==0){std::vector<char*>args;args.reserve(argv.size()+1);for(const auto&s:argv)args.push_back(const_cast<char*>(s.c_str()));args.push_back(nullptr);::execvp(args[0],args.data());std::cerr<<"execvp failed for "<<argv[0]<<": errno="<<errno<<"\n";::_exit(127);}
  int status=0;while(::waitpid(pid,&status,0)<0){if(errno==EINTR)continue;throw std::runtime_error("waitpid failed");}
  if(WIFEXITED(status))return WEXITSTATUS(status);
  if(WIFSIGNALED(status))return 128+WTERMSIG(status);
  return 255;
#endif
}
static void run_stage(const std::string&name,const std::string&cmd,const std::vector<fs::path>&outs,const fs::path&receipt,const Config&c,const std::vector<fs::path>&inputs={},const std::string&config_contract="",bool allow_empty=false,const std::vector<std::string>&argv={}) {
  std::cerr<<"+ ["<<name<<"] "<<cmd<<"\n";if(c.dry)return;
  const auto input_fingerprint=fingerprint_files(inputs);const auto&input_sha=input_fingerprint.first;const auto&input_files=input_fingerprint.second;
  const auto stage_sha=sha256_text(name),command_sha=sha256_text(cmd),config_sha=sha256_text(config_contract);
  const auto contract_sha=sha256_text("stage="+stage_sha+"\ncommand="+command_sha+"\ninputs="+input_sha+"\nconfig="+config_sha+"\n");
  if(fs::exists(receipt)){
    if(!c.resume)throw std::runtime_error("receipt already exists (use --resume): "+receipt.string());
    std::ifstream rf(receipt,std::ios::binary);std::stringstream rb;rb<<rf.rdbuf();auto rs=rb.str();
    verify_outputs(outs,allow_empty);auto current_outputs=file_set_fingerprint(outs);
    try{verify_receipt_fields(rs,stage_sha,command_sha,input_sha,config_sha,contract_sha,current_outputs);}catch(const std::exception&e){throw std::runtime_error("receipt command/input/config/output size+SHA256 mismatch: "+receipt.string()+" ("+e.what()+")");}
    std::cerr<<"  resume: verified PASS, command/input/config and output size+SHA256\n";return;
  }
  for(const auto&p:outs)if(fs::exists(p))throw std::runtime_error("unreceipted output exists; refusing overwrite: "+p.string());
  fs::create_directories(receipt.parent_path());int rc=execute_argv(argv,cmd);
  if(rc!=0)throw std::runtime_error("stage failed rc="+std::to_string(rc)+": "+name);
  const auto post_execution_input_sha=file_set_fingerprint(inputs);if(post_execution_input_sha!=input_sha)throw std::runtime_error("stage input/tool fingerprint changed during execution: "+name);
  verify_outputs(outs,allow_empty);
  auto fingerprint=fingerprint_files(outs);const auto&outputs_sha=fingerprint.first;const auto&output_files=fingerprint.second;std::ostringstream r;r<<"{\n  \"schema\": \"sag-stage3b-receipt-v2\",\n  \"status\": \"PASS\",\n  \"stage\": \""<<json(name)<<"\",\n  \"stage_sha256\": \""<<stage_sha<<"\",\n  \"command\": \""<<json(cmd)<<"\",\n  \"command_sha256\": \""<<command_sha<<"\",\n  \"input_fingerprint_sha256\": \""<<input_sha<<"\",\n  \"config_sha256\": \""<<config_sha<<"\",\n  \"contract_sha256\": \""<<contract_sha<<"\",\n  \"outputs_fingerprint_sha256\": \""<<outputs_sha<<"\",\n  \"inputs\": [";
  for(size_t i=0;i<input_files.size();++i){const auto&x=input_files[i];if(i)r<<",";r<<"\n    {\"path\": \""<<json(x.canonical.string())<<"\", \"bytes\": "<<x.bytes<<", \"sha256\": \""<<x.sha256<<"\"}";}r<<"\n  ],\n  \"outputs\": [";
  for(size_t i=0;i<output_files.size();++i){const auto&x=output_files[i];if(i)r<<",";r<<"\n    {\"path\": \""<<json(x.canonical.string())<<"\", \"bytes\": "<<x.bytes<<", \"sha256\": \""<<x.sha256<<"\"}";}r<<"\n  ]\n}\n";
  write_once(receipt,r.str());
}

static Config parse(int argc,char**argv) {
  Config c;
  for(int i=1;i<argc;++i){std::string a=argv[i];auto val=[&](){if(++i>=argc)throw std::runtime_error("missing value for "+a);return fs::path(argv[i]);};
    if(a=="--manifest")c.manifest=val();else if(a=="--quality-manifest")c.quality=val();else if(a=="--marker-map")c.marker_map=val();
    else if(a=="--out")c.out=val();
    else if(a=="--ani-engine")c.ani_engine=val();else if(a=="--makeblastdb")c.makeblastdb=val();else if(a=="--blastn")c.blastn=val();
    else if(a=="--python")c.python=val();else if(a=="--leiden-backend")c.leiden_backend=val();else if(a=="--subass")c.subass=val();else if(a=="--flye-root")c.flye_root=val();
    else if(a=="--threads"){auto s=val().string();c.threads=static_cast<int>(integer(s,"threads"));if(c.threads<1)throw std::runtime_error("threads must be >=1");}
    else if(a=="--dry-run")c.dry=true;else if(a=="--resume")c.resume=true;else if(a=="--self-test")c.self_test=true;else if(a=="--allow-experimental-ani-engine")c.allow_experimental_ani=true;
    else if(a=="--help"||a=="-h"){
      std::cout<<"sag-stage3b-tractor --manifest 03B_unclassified_pending.tsv --quality-manifest quality.tsv\n"
        "  --marker-map bac120_marker_nt_map.tsv\n"
        "  --ani-engine gtdb-ani-af --allow-experimental-ani-engine\n"
        "  --leiden-backend stage3b_signed_leiden.py --out DIR [--threads N --resume --dry-run]\n"
        "  Stage3B does not search GTDB references: every Dna2bit-negative SAG passing\n"
        "  max_contig>=1000 and CheckM2 contamination<5 enters an exact all-pairs triangle.\n";
      std::exit(0);
    }else throw std::runtime_error("unknown option: "+a);
  }
  if(c.self_test)return c;
  if(c.manifest.empty()||c.quality.empty()||c.marker_map.empty()||c.out.empty()||c.leiden_backend.empty())
    throw std::runtime_error("required: --manifest --quality-manifest --marker-map --leiden-backend --out");
  if(!c.allow_experimental_ani)throw std::runtime_error("gtdb-ani-af is experimental and not production-validated; explicit --allow-experimental-ani-engine is required");
  if(c.dry&&c.resume)throw std::runtime_error("--dry-run and --resume are mutually exclusive");
  return c;
}

static std::map<std::string,Sag> load_sags(const Config&c) {
  TsvReader manifest(c.manifest);auto ci=manifest.column({"sag_id","SAG_id"});auto ca=manifest.column({"assembly_fasta"});auto cr=manifest.column({"reason"});
  std::map<std::string,Sag>sags;std::vector<std::string> row;
  while(manifest.next(row)){Sag s;s.id=safe_id(row[ci]);s.assembly=canonical_existing(row[ca]);s.reason=trim(row[cr]);auto lr=lower(s.reason);
    if(!dna2bit_negative_reason(lr))throw std::runtime_error("reason must be an exact DNA2bit-negative enum for "+s.id+": "+s.reason);
    auto id=s.id;if(!sags.emplace(id,std::move(s)).second)throw std::runtime_error("duplicate SAG id: "+id);
  }
  if(sags.empty())throw std::runtime_error("Stage 3B manifest has no SAGs");
  TsvReader quality(c.quality);auto qi=quality.column({"sag_id","SAG_id"});auto qgc=quality.column({"gc_pct","GC","gc"});
  auto qmax=quality.column({"max_contig","max_contig_bp"});auto qcon=quality.column({"checkm2_contamination","contamination"});auto qtot=quality.column({"total_bp","total_len","assembly_total_bp"},false);
  std::set<std::string> got;
  while(quality.next(row)){auto id=row[qi];auto it=sags.find(id);if(it==sags.end())continue;if(!got.insert(id).second)throw std::runtime_error("duplicate quality row: "+id);
    it->second.declared_gc=number(row[qgc],id+" gc_pct");it->second.declared_max=integer(row[qmax],id+" max_contig");it->second.contamination=number(row[qcon],id+" CheckM2 contamination");if(qtot!=std::numeric_limits<std::size_t>::max())it->second.declared_total=integer(row[qtot],id+" total_bp");
  }
  for(auto&[id,s]:sags){if(!got.count(id))throw std::runtime_error("quality manifest missing SAG: "+id);s.actual=fasta_stats(s.assembly);
    if(s.declared_total&&s.declared_total!=s.actual.total)throw std::runtime_error("declared total_bp disagrees with FASTA for "+id);
    if(s.declared_max!=s.actual.max_contig)throw std::runtime_error("declared max_contig disagrees with FASTA for "+id);
    if(!(s.declared_gc>=0&&s.declared_gc<=100&&s.contamination>=0))throw std::runtime_error("invalid GC/contamination for "+id);
    // Historical fasta_stats uses all sequence characters (including N and
    // ambiguity symbols) in the denominator, not A/C/G/T-only length.
    double actual_gc=100.0*double(s.actual.gc)/double(s.actual.total);
    if(std::fabs(actual_gc-s.declared_gc)>0.00011)throw std::runtime_error("declared 4-decimal gc_pct disagrees with FASTA for "+id);
    // The upstream total-bp gate is additional; it is not a replacement for
    // the historical Stage 3B max-contig and CheckM2 gates below.
    if(s.actual.total<1000)throw std::runtime_error("Stage 1 total_bp<1000 SAG leaked into 3B: "+id);
  }
  return sags;
}

static std::unordered_map<std::string,std::string> path_to_sag(const std::map<std::string,Sag>&sags,const std::set<std::string>*only=nullptr) {
  std::unordered_map<std::string,std::string> m;
  for(const auto&[id,s]:sags){if(only&&!only->count(id))continue;auto full=fs::weakly_canonical(s.assembly).string();if(m.count(full))throw std::runtime_error("two SAGs share assembly path");m[full]=id;auto base=s.assembly.filename().string();auto key="BASE:"+base;if(m.count(key))m[key]="";else m[key]=id;}return m;
}
static std::string resolve_sag(const std::string&p,const std::unordered_map<std::string,std::string>&m) {
  std::error_code ec;auto full=fs::weakly_canonical(fs::path(p),ec);if(!ec){auto i=m.find(full.string());if(i!=m.end())return i->second;}
  auto i=m.find("BASE:"+fs::path(p).filename().string());if(i!=m.end()&&!i->second.empty())return i->second;
  throw std::runtime_error("cannot unambiguously map external result path to SAG: "+p);
}

static std::map<std::string,std::map<std::string,std::string>> load_markers(const fs::path&p,const std::set<std::string>&keep) {
  TsvReader rows(p);auto cs=rows.column({"SAG_id","sag_id"});auto cm=rows.column({"marker_id"});auto cq=rows.column({"sequence"});
  std::map<std::string,std::map<std::string,std::string>> per;
  std::vector<std::string> row;while(rows.next(row)){auto s=row[cs];if(!keep.count(s))continue;auto seq=std::move(row[cq]);seq.erase(std::remove_if(seq.begin(),seq.end(),[](unsigned char c){return std::isspace(c);}),seq.end());for(char&c:seq)c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if(seq.size()<150)continue;
    for(char c:seq)if(!(c=='A'||c=='C'||c=='G'||c=='T'||c=='N'||c=='R'||c=='Y'||c=='S'||c=='W'||c=='K'||c=='M'||c=='B'||c=='D'||c=='H'||c=='V'))throw std::runtime_error("invalid marker sequence for "+s);
    auto&slot=per[row[cm]][s];if(seq.size()>slot.size())slot=std::move(seq);
  }return per;
}
static void validate_marker_map_header(const fs::path&p) {
  TsvReader rows(p);rows.column({"SAG_id","sag_id"});rows.column({"marker_id"});rows.column({"sequence"});
}
static std::string marker_name(const std::string&mid,size_t idx){std::string s;for(char c:mid)s+=std::isalnum(static_cast<unsigned char>(c))?c:'_';if(s.size()>50)s.resize(50);return std::string("M")+(idx<10?"00":idx<100?"0":"")+std::to_string(idx)+"_"+s;}
struct MarkerJob {std::string marker_id;fs::path dir,fasta,db,blast_out,nhr,nin,nsq;};
static void run_marker_job(const MarkerJob&j,const Config&c){
  std::string mk=q(c.makeblastdb)+" -in "+q(j.fasta)+" -dbtype nucl -out "+q(j.db);
  run_stage("makeblastdb:"+j.marker_id,mk,{j.nhr,j.nin,j.nsq},j.dir/"MAKEBLASTDB.PASS.json",c,{j.fasta,c.makeblastdb},"dbtype=nucl\nmarker="+j.marker_id,false,
    {c.makeblastdb.string(),"-in",j.fasta.string(),"-dbtype","nucl","-out",j.db.string()});
  std::string bl=q(c.blastn)+" -query "+q(j.fasta)+" -db "+q(j.db)+" -evalue 1e-5 -max_target_seqs 10000 -outfmt '6 qseqid sseqid pident length' -num_threads 2 -out "+q(j.blast_out);
  run_stage("blastn:"+j.marker_id,bl,{j.blast_out},j.dir/"BLASTN.PASS.json",c,{j.fasta,j.nhr,j.nin,j.nsq,c.blastn},"evalue=1e-5\nmax_target_seqs=10000\noutfmt=6 qseqid sseqid pident length\nthreads=2",true,
    {c.blastn.string(),"-query",j.fasta.string(),"-db",j.db.string(),"-evalue","1e-5","-max_target_seqs","10000","-outfmt","6 qseqid sseqid pident length","-num_threads","2","-out",j.blast_out.string()});
}
static void run_marker_jobs(const std::vector<MarkerJob>&jobs,const Config&c){
  if(jobs.empty())return;
  const std::size_t limit=std::max<std::size_t>(1,static_cast<std::size_t>(c.threads)/2);
#ifdef _WIN32
  for(const auto&j:jobs)run_marker_job(j,c);metrics.marker_peak_jobs=1;
#else
  std::map<pid_t,std::size_t> running;std::size_t next=0;bool failed=false;std::string failed_marker;
  auto launch=[&](std::size_t index){pid_t pid=::fork();if(pid<0)throw std::runtime_error("fork failed while launching marker BLAST");if(pid==0){::setpgid(0,0);try{run_marker_job(jobs[index],c);::_exit(0);}catch(const std::exception&e){std::cerr<<"marker worker failed ["<<jobs[index].marker_id<<"]: "<<e.what()<<"\n";::_exit(2);}}::setpgid(pid,pid);running.emplace(pid,index);metrics.marker_peak_jobs=std::max(metrics.marker_peak_jobs,running.size());};
  while(next<jobs.size()||!running.empty()){
    while(!failed&&next<jobs.size()&&running.size()<limit)launch(next++);
    int status=0;pid_t pid;do{pid=::waitpid(-1,&status,0);}while(pid<0&&errno==EINTR);if(pid<0)throw std::runtime_error("waitpid failed for marker workers");
    auto it=running.find(pid);if(it==running.end())throw std::runtime_error("reaped unknown marker worker");bool ok=WIFEXITED(status)&&WEXITSTATUS(status)==0;if(!ok&&!failed){failed=true;failed_marker=jobs[it->second].marker_id;for(const auto&[other,_]:running)if(other!=pid)::kill(-other,SIGTERM);}running.erase(it);
  }
  if(failed)throw std::runtime_error("bounded marker BLAST worker failed: "+failed_marker);
#endif
}

struct SubassJob {std::string cluster;fs::path input,run;std::size_t n=0;};
static void run_subass_job(const SubassJob&j,const Config&c,int job_threads,const std::string&phase){
  auto fr=c.flye_root,modules=fr/"bin/flye-modules",minimap=fr/"bin/flye-minimap2",samtools=fr/"bin/flye-samtools",package=fr/"lib/python3.9/site-packages/flye",config=fr/"lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg";
  std::string cmd=q(c.subass)+" --reads "+q(j.input)+" --out-dir "+q(j.run)+" --flye-modules "+q(modules)+" --minimap2 "+q(minimap)+" --samtools "+q(samtools)+" --package-root "+q(package)+" --config "+q(config)+" --threads "+std::to_string(job_threads)+" --assemble-threads 1 --phase "+phase+" --no-overlap-policy passthrough";
  std::vector<fs::path> inputs={j.input,c.subass};inputs.insert(inputs.end(),c.flye_provenance.begin(),c.flye_provenance.end());
  if(phase=="finish"){
    inputs.push_back(j.run/"00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE");
    auto draft=j.run/"00-assembly/draft_assembly.fasta";
    if(fs::is_regular_file(draft)&&fs::file_size(draft)>0)inputs.push_back(draft);
    else {inputs.push_back(j.run/"assembly.fasta");inputs.push_back(j.run/"NO_OVERLAP_PASSTHROUGH.PASS.json");}
  }
  const auto output=phase=="assemble"?j.run/"00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE":j.run/"40-polishing/CPP_FINISH_PHASE_COMPLETE";
  const auto receipt=phase=="assemble"?j.run/"STAGE3B_SUBASS_ASSEMBLE.PASS.json":j.run/"STAGE3B_SUBASS.PASS.json";
  run_stage("cpp-subass-"+phase+":"+j.cluster,cmd,{output},receipt,c,inputs,"flye_root="+fs::absolute(fr).string()+"\nphase="+phase+"\nassemble_threads=1\ndownstream_threads="+std::to_string(job_threads)+"\nno_overlap_policy=passthrough",false,
    {c.subass.string(),"--reads",j.input.string(),"--out-dir",j.run.string(),"--flye-modules",modules.string(),"--minimap2",minimap.string(),"--samtools",samtools.string(),"--package-root",package.string(),"--config",config.string(),"--threads",std::to_string(job_threads),"--assemble-threads","1","--phase",phase,"--no-overlap-policy","passthrough"});
}
static void run_subass_phase(const std::vector<SubassJob>&scheduled,const Config&c,int job_threads,std::size_t limit,const std::string&phase){
#ifdef _WIN32
  for(const auto&j:scheduled)run_subass_job(j,c,job_threads,phase);metrics.subass_peak_jobs=1;
#else
  std::map<pid_t,std::size_t> running;std::size_t next=0;bool failed=false;std::string failed_cluster;
  auto launch=[&](std::size_t index){pid_t pid=::fork();if(pid<0)throw std::runtime_error("fork failed while launching subassembly");if(pid==0){::setpgid(0,0);try{run_subass_job(scheduled[index],c,job_threads,phase);::_exit(0);}catch(const std::exception&e){std::cerr<<"subassembly "<<phase<<" worker failed ["<<scheduled[index].cluster<<"]: "<<e.what()<<"\n";::_exit(2);}}::setpgid(pid,pid);running.emplace(pid,index);metrics.subass_peak_jobs=std::max(metrics.subass_peak_jobs,running.size());};
  while(next<scheduled.size()||!running.empty()){
    while(!failed&&next<scheduled.size()&&running.size()<limit)launch(next++);
    int status=0;pid_t pid;do{pid=::waitpid(-1,&status,0);}while(pid<0&&errno==EINTR);if(pid<0)throw std::runtime_error("waitpid failed for subassembly workers");
    auto it=running.find(pid);if(it==running.end())throw std::runtime_error("reaped unknown subassembly worker");bool ok=WIFEXITED(status)&&WEXITSTATUS(status)==0;if(!ok&&!failed){failed=true;failed_cluster=scheduled[it->second].cluster;for(const auto&[other,_]:running)if(other!=pid)::kill(-other,SIGTERM);}running.erase(it);
  }
  if(failed)throw std::runtime_error("bounded subassembly "+phase+" worker failed: "+failed_cluster);
#endif
}
static void run_subass_jobs(const std::vector<SubassJob>&jobs,const Config&c){
  if(jobs.empty())return;
  // First saturate the budget with deterministic one-thread assembly jobs;
  // then use verified two-thread downstream jobs at half the concurrency.
  // Both passes are LPT ordered and their declared aggregate budget is no
  // greater than the caller's --threads value.
  std::vector<SubassJob>scheduled=jobs;
  std::sort(scheduled.begin(),scheduled.end(),[](const SubassJob&a,const SubassJob&b){auto as=fs::file_size(a.input),bs=fs::file_size(b.input);if(as!=bs)return as>bs;return a.cluster<b.cluster;});
  const int budget=std::max(1,c.threads);
  const auto assemble_limit=std::min<std::size_t>(scheduled.size(),static_cast<std::size_t>(budget));
  run_subass_phase(scheduled,c,1,assemble_limit,"assemble");
  const int finish_threads=budget>=2?2:1;
  const auto finish_limit=std::min<std::size_t>(scheduled.size(),static_cast<std::size_t>(std::max(1,budget/finish_threads)));
  run_subass_phase(scheduled,c,finish_threads,finish_limit,"finish");
}
static void parse_blast(const fs::path&p,std::map<std::pair<std::string,std::string>,MarkerAgg>&agg) {
  std::ifstream f(p);if(!f)throw std::runtime_error("cannot read blastn result: "+p.string());std::map<std::pair<std::string,std::string>,std::pair<double,int>>best;std::string line;size_t ln=0;
  while(std::getline(f,line)){++ln;if(line.empty())continue;auto x=split(line,'\t');if(x.size()!=4)throw std::runtime_error(p.string()+":"+std::to_string(ln)+" expected outfmt 6 qseqid sseqid pident length");if(x[0]==x[1])continue;auto key=std::minmax(x[0],x[1]);double pid=number(x[2],"pident");int len=static_cast<int>(integer(x[3],"alignment length"));auto i=best.find(key);if(i==best.end()||len>i->second.second)best[key]={pid,len};}
  for(const auto&[k,v]:best)agg[k].add(v.first);
}
static void append_fasta(const fs::path&src,std::ostream&dst,const std::string&sag) {
  std::ifstream f(src);if(!f)throw std::runtime_error("cannot read FASTA: "+src.string());std::string line;bool header=false;
  while(std::getline(f,line)){if(!line.empty()&&line[0]=='>'){header=true;auto h=trim(line.substr(1));if(h.empty())throw std::runtime_error("empty FASTA header: "+src.string());dst<<">"<<sag<<"__"<<h<<"\n";}else dst<<line<<"\n";}
  if(!header)throw std::runtime_error("FASTA lacks header: "+src.string());
}
static void materialize_merged_fasta(const fs::path&out,const std::vector<std::string>&members,const std::map<std::string,Sag>&sags,const Config&c){
  fs::create_directories(out.parent_path());fs::path receipt=out.parent_path()/"MERGED_INPUT.PASS.json",tmp=out;tmp+=".recomputed.tmp";
  std::vector<fs::path> sources;sources.reserve(members.size());for(const auto&id:members)sources.push_back(sags.at(id).assembly);const auto source_fp=file_set_fingerprint(sources);
  if(fs::exists(out)&&!c.resume)throw std::runtime_error("refusing overwrite merged FASTA: "+out.string());
  if(!fs::exists(out)&&fs::exists(receipt))throw std::runtime_error("merged-input receipt exists but output is missing: "+receipt.string());
  if(fs::exists(out)&&fs::exists(receipt)){
    std::ifstream f(receipt,std::ios::binary);std::stringstream b;b<<f.rdbuf();auto r=b.str();auto output_fp=file_set_fingerprint({out});
    if(receipt_field(r,"schema")!="sag-stage3b-merged-input-v1"||receipt_field(r,"status")!="PASS"||receipt_field(r,"source_fingerprint_sha256")!=source_fp||receipt_field(r,"output_fingerprint_sha256")!=output_fp||receipt_field(r,"member_count")!=std::to_string(members.size()))throw std::runtime_error("merged-input source/output receipt mismatch: "+receipt.string());
    return;
  }
  if(fs::exists(tmp))throw std::runtime_error("stale merged FASTA temporary: "+tmp.string());
  {std::ofstream f(tmp,std::ios::binary);if(!f)throw std::runtime_error("cannot write merged FASTA temporary");for(const auto&id:members)append_fasta(sags.at(id).assembly,f,id);if(!f)throw std::runtime_error("failed writing merged FASTA");}
  if(fs::exists(out)){if(fs::file_size(out)!=fs::file_size(tmp)||sha256_file(out)!=sha256_file(tmp)){fs::remove(tmp);throw std::runtime_error("resume merged FASTA differs from deterministic recomputation: "+out.string());}fs::remove(tmp);}else fs::rename(tmp,out);
  const auto output_fp=file_set_fingerprint({out});std::ostringstream r;r<<"{\n  \"schema\": \"sag-stage3b-merged-input-v1\",\n  \"status\": \"PASS\",\n  \"member_count\": \""<<members.size()<<"\",\n  \"source_fingerprint_sha256\": \""<<source_fp<<"\",\n  \"output_fingerprint_sha256\": \""<<output_fp<<"\"\n}\n";write_once(receipt,r.str());
}

static void self_test() {
  if(!dna2bit_negative_reason("dna2bit_rejected_or_no_hit")||dna2bit_negative_reason("dna2bit_not_negative")||dna2bit_negative_reason("positive_not_rejected"))throw std::runtime_error("DNA2bit reason enum self-test failed");
  Sag gate;gate.actual.max_contig=1000;gate.contamination=4.999999;if(!stage3b_quality_pass(gate))throw std::runtime_error("quality gate rejected inclusive max-contig/strict contamination pass");gate.contamination=5.0;if(stage3b_quality_pass(gate))throw std::runtime_error("quality gate accepted contamination=5");gate.contamination=0;gate.actual.max_contig=999;if(stage3b_quality_pass(gate))throw std::runtime_error("quality gate accepted max_contig<1000");
  if(exact_pair_count(0)!=0||exact_pair_count(1)!=0||exact_pair_count(10)!=45)throw std::runtime_error("exact N choose 2 pair-count self-test failed");
  if(sha256_text("abc")!="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")throw std::runtime_error("SHA256 self-test failed");
  std::string receipt="{\"schema\": \"sag-stage3b-receipt-v2\", \"status\": \"PASS\", \"stage_sha256\": \"s\", \"command_sha256\": \"c\", \"input_fingerprint_sha256\": \"i\", \"config_sha256\": \"g\", \"contract_sha256\": \"k\", \"outputs_fingerprint_sha256\": \"o\"}";verify_receipt_fields(receipt,"s","c","i","g","k","o");bool rejected=false;try{verify_receipt_fields(receipt,"s","tampered","i","g","k","o");}catch(...){rejected=true;}if(!rejected)throw std::runtime_error("tampered receipt self-test failed");
  double w=(97.5-95.0)/5.0*0.4+0.01;if(std::fabs(w-0.21)>1e-12)throw std::runtime_error("positive weight self-test failed");
  if(std::fabs((97.0-96.0)/97.0-0.010309278350515464)>1e-12)throw std::runtime_error("negative weight self-test failed");
  MarkerAgg marker;marker.add(98.0);marker.add(96.0);if(marker.count!=2||marker.sum!=194.0||marker.min!=96.0)throw std::runtime_error("constant-memory marker aggregate self-test failed");
  auto nonce=std::chrono::high_resolution_clock::now().time_since_epoch().count();auto base=fs::temp_directory_path()/("stage3b_stream_selftest_"+std::to_string(nonce));auto tsv=base;tsv+=".tsv";auto fasta=base;fasta+=".fa";
  {std::ofstream f(tsv,std::ios::binary);f<<"id\tvalue\nA\t 1 \nB\t2\n";}{TsvReader r(tsv);auto id=r.column({"id"}),value=r.column({"value"});std::vector<std::string> row;std::size_t n=0;while(r.next(row)){++n;if(row[id]!=(n==1?"A":"B")||row[value]!=(n==1?"1":"2"))throw std::runtime_error("streaming TSV self-test failed");}if(n!=2)throw std::runtime_error("streaming TSV row-count self-test failed");}
  const std::string fasta_text=">a\nACGN\n>b\nGG\n";{std::ofstream f(fasta,std::ios::binary);f<<fasta_text;}auto stats=fasta_stats(fasta);if(stats.total!=6||stats.max_contig!=4||stats.gc!=4||sha256_file(fasta)!=sha256_text(fasta_text))throw std::runtime_error("single-pass FASTA/stat/hash self-test failed");fs::remove(tsv);fs::remove(fasta);
  std::cout<<"PASS stage3b C++ self-test reference_search=disabled triangle_mode=exact\n";
}

int main(int argc,char**argv)try{
  Config c=parse(argc,argv);if(c.self_test){self_test();return 0;}
  auto sags=load_sags(c);validate_marker_map_header(c.marker_map);
  c.ani_engine=resolve_tool(c.ani_engine,c.dry);c.makeblastdb=resolve_tool(c.makeblastdb,c.dry);c.blastn=resolve_tool(c.blastn,c.dry);c.python=resolve_tool(c.python,c.dry);c.subass=resolve_tool(c.subass,c.dry);c.leiden_backend=canonical_existing(c.leiden_backend);
  if(!c.dry){c.flye_root=fs::canonical(c.flye_root);c.flye_provenance=collect_flye_provenance(c.flye_root);}
  if(!c.dry){if(fs::exists(c.out)&&!c.resume)throw std::runtime_error("output directory exists; write-once policy requires a new path or --resume: "+c.out.string());fs::create_directories(c.out);}

  // Dna2bit negativity is already proven by the exact manifest reason enum.
  // Apply the two remaining historical quality gates.  No GTDB reference
  // search or second taxonomy-positive exclusion exists in this route.
  std::map<std::string,Sag> eligible;fs::path quality_audit=c.out/"quality_audit.tsv";
  if(!c.dry)write_generated(quality_audit,c,[&](std::ostream&qa){qa<<"sag_id\tassembly_fasta\ttotal_bp\tmax_contig\tgc_pct\tcheckm2_contamination\tstage3b_quality_pass\treason\n";
    for(const auto&[id,s]:sags){bool pass=stage3b_quality_pass(s);qa<<id<<"\t"<<s.assembly.string()<<"\t"<<s.actual.total<<"\t"<<s.actual.max_contig<<"\t"<<std::fixed<<std::setprecision(5)<<s.declared_gc<<"\t"<<s.contamination<<"\t"<<(pass?"1":"0")<<"\t"<<(s.actual.max_contig<1000?"max_contig_lt_1000":s.contamination>=5.0?"checkm2_contamination_ge_5":"pass")<<"\n";if(pass)eligible.emplace(id,s);}
  });else for(const auto&[id,s]:sags)if(stage3b_quality_pass(s))eligible.emplace(id,s);
  if(eligible.empty())throw std::runtime_error("no SAG passes max_contig>=1000 and CheckM2 contamination<5.0");
  if(c.dry)std::cerr<<"DRY-RUN would write quality_audit.tsv (eligible="<<eligible.size()<<")\n";

  std::set<std::string> graph_nodes;for(const auto&[id,_]:eligible)graph_nodes.insert(id);
  const auto graph_node_count=graph_nodes.size();
  const std::uint64_t requested_pairs=exact_pair_count(graph_node_count);
  fs::path node_dir=c.out/"01_cellbit_negative_quality_pass",node_tsv=node_dir/"cellbit_negative_quality_pass.tsv",node_list=node_dir/"cellbit_negative_quality_pass.list";
  if(c.dry){
    std::cerr<<"GTDB reference search: disabled; reference_positive_exclusion=false\n";
    std::cerr<<"exact pairwise evaluations: "<<requested_pairs<<" = "<<graph_node_count<<"*("<<graph_node_count<<"-1)/2\n";
    std::cerr<<"+ [sag_pairwise_ani_af_triangle_exact] "<<q(c.ani_engine)<<" triangle --list "<<q(node_list)<<" --out "<<q(c.out/"02_pairwise_ani_af/ani_triangle_sparse.tsv")<<" --stats-out "<<q(c.out/"02_pairwise_ani_af/triangle_stats.json")<<" --threads "<<c.threads<<" --profile medium --edge-mode sparse --min-af 15 --triangle-mode exact --triangle-global-seed-join --triangle-global-prefix-bits 22\n";
    std::cerr<<"DRY-RUN stops before data-dependent BLAST/Leiden/subassembly planning; no files or commands executed\n";return 0;
  }
  write_generated(node_tsv,c,[&](std::ostream&ds){ds<<"sag_id\tassembly_fasta\tgc_pct\n";for(const auto&id:graph_nodes){const auto&s=eligible.at(id);ds<<id<<"\t"<<s.assembly.string()<<"\t"<<s.declared_gc<<"\n";}});
  write_generated(node_list,c,[&](std::ostream&dl){for(const auto&id:graph_nodes)dl<<eligible.at(id).assembly.string()<<"\n";});

  fs::path tri_dir=c.out/"02_pairwise_ani_af",tri=tri_dir/"ani_triangle_sparse.tsv",tri_stats=tri_dir/"triangle_stats.json";std::string tri_cmd=q(c.ani_engine)+" triangle --list "+q(node_list)+" --out "+q(tri)+" --stats-out "+q(tri_stats)+" --threads "+std::to_string(c.threads)+" --profile medium --edge-mode sparse --min-af 15 --triangle-mode exact --triangle-global-seed-join --triangle-global-prefix-bits 22";std::vector<fs::path> graph_assemblies;graph_assemblies.reserve(graph_nodes.size());for(const auto&id:graph_nodes)graph_assemblies.push_back(eligible.at(id).assembly);std::vector<fs::path> tri_inputs={node_list,c.ani_engine};tri_inputs.insert(tri_inputs.end(),graph_assemblies.begin(),graph_assemblies.end());run_stage("sag_pairwise_ani_af_triangle_exact",tri_cmd,{tri,tri_stats},tri_dir/"PASS.json",c,tri_inputs,"threads="+std::to_string(c.threads)+"\nprofile=medium\nedge_mode=sparse\nmin_af=max(AF_ref,AF_query)>=15\ntriangle_mode=exact\nexact_execution=global-seed-join\nglobal_prefix_bits=22\npair_selection=all_N_choose_2\nreference_positive_exclusion=false",false,
    {c.ani_engine.string(),"triangle","--list",node_list.string(),"--out",tri.string(),"--stats-out",tri_stats.string(),"--threads",std::to_string(c.threads),"--profile","medium","--edge-mode","sparse","--min-af","15","--triangle-mode","exact","--triangle-global-seed-join","--triangle-global-prefix-bits","22"});
  // Accept the N choose 2 claim only from the ANI engine's write-once receipt.
  // The sparse table intentionally retains only rows passing min-af=15.
  const auto stats=read_text(tri_stats);if(receipt_field(stats,"schema")!="gtdb-ani-af-triangle-stats-v1"||receipt_field(stats,"status")!="PASS"||receipt_field(stats,"triangle_mode")!="exact"||receipt_field(stats,"edge_mode")!="sparse"||receipt_field(stats,"exact_execution")!="global-seed-join")throw std::runtime_error("invalid exact-triangle stats identity");
  if(receipt_uint(stats,"nodes")!=graph_node_count||receipt_uint(stats,"pairs_expected")!=requested_pairs||receipt_uint(stats,"pairs_evaluated")!=requested_pairs)throw std::runtime_error("exact-triangle stats do not prove N*(N-1)/2 evaluations");
  if(receipt_field(stats,"input_list_path")!=fs::weakly_canonical(node_list).string()||receipt_field(stats,"output_path")!=fs::weakly_canonical(tri).string()||receipt_field(stats,"input_list_sha256")!=sha256_file(node_list)||receipt_field(stats,"input_content_set_sha256")!=triangle_input_content_set_sha256(graph_assemblies)||receipt_field(stats,"output_sha256")!=sha256_file(tri))throw std::runtime_error("exact-triangle stats input/output path or SHA256 mismatch");
  const std::uint64_t evaluated_pairs=receipt_uint(stats,"pairs_evaluated"),stats_rows_emitted=receipt_uint(stats,"rows_emitted");
  TsvReader tr(tri);auto tref=tr.column({"Ref_file","Reference","Ref_name"});auto tq=tr.column({"Query_file","Query"});auto tani=tr.column({"ANI"});auto tafr=tr.column({"Align_fraction_ref","AF_ref"});auto tafq=tr.column({"Align_fraction_query","AF_query"});auto pm=path_to_sag(eligible,&graph_nodes);std::map<std::pair<std::string,std::string>,std::tuple<double,double,double,double,double>> pos;std::set<std::pair<std::string,std::string>> emitted_pairs;std::vector<std::string> tr_row;std::uint64_t triangle_rows=0;
  while(tr.next(tr_row)){++triangle_rows;auto a=resolve_sag(tr_row[tref],pm),b=resolve_sag(tr_row[tq],pm);if(a==b)throw std::runtime_error("exact triangle emitted a self pair: "+a);auto key=std::minmax(a,b);if(!emitted_pairs.insert(key).second)throw std::runtime_error("exact triangle emitted duplicate pair: "+key.first+" / "+key.second);double ani=number(tr_row[tani],"triangle ANI"),afr=number(tr_row[tafr],"triangle AF_ref"),afq=number(tr_row[tafq],"triangle AF_query");if(ani<95)continue;double gd=std::fabs(eligible.at(a).declared_gc-eligible.at(b).declared_gc);if(gd>2.0)continue;double w=(ani-95.0)/5.0*(std::max(afr,afq)/100.0)+0.01;auto i=pos.find(key);if(i==pos.end()||std::get<4>(i->second)<w)pos[key]={ani,afr,afq,gd,w};}
  if(triangle_rows!=stats_rows_emitted||triangle_rows>evaluated_pairs)throw std::runtime_error("sparse triangle row count disagrees with exact-triangle stats receipt");
  // Apply the existing 12-significant-digit contract before the first
  // numeric field.  std::setprecision is persistent; placing it immediately
  // before weight left only the first row's ANI/AF/GC at the default six
  // significant digits even though weight used the unrounded values.
  fs::path pos_tsv=tri_dir/"positive_edges.tsv";write_generated(pos_tsv,c,[&](std::ostream&pe){pe<<"SAG1\tSAG2\tANI\tAF_ref\tAF_query\tgc_diff\tweight\n"<<std::setprecision(12);for(const auto&[k,v]:pos)pe<<k.first<<"\t"<<k.second<<"\t"<<std::get<0>(v)<<"\t"<<std::get<1>(v)<<"\t"<<std::get<2>(v)<<"\t"<<std::get<3>(v)<<"\t"<<std::get<4>(v)<<"\n";});

  // Historical marker targeting is a separate positive-only Leiden pass.
  // Only members of its size>=10 communities enter marker BLAST; the final
  // signed Leiden below still contains every quality-passing Cellbit-negative
  // node, including nodes outside those marker-target communities.
  fs::path pre_dir=c.out/"03_positive_precluster",marker_targets=pre_dir/"marker_target_membership.tsv",pre_report=pre_dir/"precluster_report.json";
  std::string pc=q(c.python)+" "+q(c.leiden_backend)+" --precluster --nodes "+q(node_tsv)+" --positive "+q(pos_tsv)+" --membership "+q(marker_targets)+" --report "+q(pre_report);
  run_stage("positive_only_leiden_marker_targets",pc,{marker_targets,pre_report},pre_dir/"PASS.json",c,{node_tsv,pos_tsv,c.python,c.leiden_backend},"seed=20260811\nresolution=1\nn_iterations=-1",false,
    {c.python.string(),c.leiden_backend.string(),"--precluster","--nodes",node_tsv.string(),"--positive",pos_tsv.string(),"--membership",marker_targets.string(),"--report",pre_report.string()});
  TsvReader ptr(marker_targets);auto pts=ptr.column({"SAG_id","sag_id"});auto ptc=ptr.column({"cluster_size"});std::set<std::string> marker_keep;std::vector<std::string> pre_row;
  while(ptr.next(pre_row)){auto id=pre_row[pts];if(!graph_nodes.count(id))throw std::runtime_error("precluster returned non-eligible SAG: "+id);if(integer(pre_row[ptc],"precluster cluster_size")<10)throw std::runtime_error("precluster emitted cluster below size 10");if(!marker_keep.insert(id).second)throw std::runtime_error("duplicate precluster SAG: "+id);}
  auto markers=load_markers(c.marker_map,marker_keep);marker_keep.clear();fs::path md=c.out/"04_marker_blastn";std::vector<MarkerJob> marker_jobs;size_t mi=0;
  for(auto&[mid,seqs]:markers){if(seqs.size()<2){std::map<std::string,std::string>().swap(seqs);continue;}++mi;std::string tag=marker_name(mid,mi);MarkerJob j;j.marker_id=mid;j.dir=md/tag;j.fasta=j.dir/"marker.fna";j.db=j.dir/"db";j.blast_out=j.dir/"blastn.tsv";j.nhr=fs::path(j.db.string()+".nhr");j.nin=fs::path(j.db.string()+".nin");j.nsq=fs::path(j.db.string()+".nsq");
    write_generated(j.fasta,c,[&](std::ostream&out){for(const auto&[s,seq]:seqs)out<<">"<<s<<"\n"<<seq<<"\n";});marker_jobs.push_back(std::move(j));std::map<std::string,std::string>().swap(seqs);
  }
  std::map<std::string,std::map<std::string,std::string>>().swap(markers);const size_t jobs=marker_jobs.size();run_marker_jobs(marker_jobs,c);
  std::map<std::pair<std::string,std::string>,MarkerAgg> agg;for(const auto&j:marker_jobs)parse_blast(j.blast_out,agg);
  if(jobs==0)std::cerr<<"marker evidence: zero eligible jobs; continuing with header-only evidence and no negative edges\n";
  fs::path marker_pairs=md/"marker_pairwise_blastn.tsv",neg_edges=md/"negative_edges.tsv";size_t neg_n=0;
  write_generated(marker_pairs,c,[&](std::ostream&mp){mp<<"SAG1\tSAG2\tn_markers\tmean_pident\tmin_pident\n";for(const auto&[k,v]:agg){double mean=std::round((v.sum/v.count)*1000.0)/1000.0;mp<<k.first<<"\t"<<k.second<<"\t"<<v.count<<"\t"<<std::fixed<<std::setprecision(3)<<mean<<"\t"<<v.min<<"\n";}});
  write_generated(neg_edges,c,[&](std::ostream&ne){ne<<"SAG1\tSAG2\tn_markers\tmean_pident\tweight\n";for(const auto&[k,v]:agg){double mean=std::round((v.sum/v.count)*1000.0)/1000.0;if(v.count>=3&&mean<97.0){++neg_n;ne<<k.first<<"\t"<<k.second<<"\t"<<v.count<<"\t"<<std::fixed<<std::setprecision(3)<<mean<<"\t"<<std::setprecision(12)<<(97.0-mean)/97.0<<"\n";}}});

  // Every frozen parameter point has an independent optimiser and the backend
  // uses spawned processes to keep both the native RNG and graph state isolated.
  // Running the six points concurrently therefore changes scheduling only, not
  // the seed, parameter order, tie-breaking, or selected membership.
  const int leiden_workers=std::min(6,c.threads);
  const std::string leiden_workers_text=std::to_string(leiden_workers);
  fs::path ld=c.out/"05_signed_leiden",membership=ld/"chosen_membership.tsv",report=ld/"signed_report.json";std::string lc=q(c.python)+" "+q(c.leiden_backend)+" --nodes "+q(node_tsv)+" --positive "+q(pos_tsv)+" --marker-pairs "+q(marker_pairs)+" --negative "+q(neg_edges)+" --membership "+q(membership)+" --report "+q(report)+" --workers "+leiden_workers_text;run_stage("signed_leiden_sweep",lc,{membership,report},ld/"PASS.json",c,{node_tsv,pos_tsv,marker_pairs,neg_edges,c.python,c.leiden_backend},"seed=20260811\nn_iterations=-1\nparameters=1:1,1:3,1:10,2:3,2:10,2:30\nworkers="+leiden_workers_text,false,
    {c.python.string(),c.leiden_backend.string(),"--nodes",node_tsv.string(),"--positive",pos_tsv.string(),"--marker-pairs",marker_pairs.string(),"--negative",neg_edges.string(),"--membership",membership.string(),"--report",report.string(),"--workers",leiden_workers_text});
  TsvReader mr(membership);auto mc=mr.column({"cluster"});auto ms=mr.column({"SAG_id","sag_id"});auto mz=mr.column({"cluster_size"});std::map<std::string,std::vector<std::string>> groups;std::map<std::string,std::uint64_t> declared_sizes;std::set<std::string> assigned;std::vector<std::string> membership_row;
  while(mr.next(membership_row)){auto id=membership_row[ms],cl=safe_id(membership_row[mc]);auto declared=integer(membership_row[mz],"cluster_size");if(!graph_nodes.count(id))throw std::runtime_error("Leiden returned unknown/non-eligible SAG: "+id);if(declared<10)throw std::runtime_error("backend emitted cluster below size 10");if(!assigned.insert(id).second)throw std::runtime_error("SAG occurs in multiple output clusters: "+id);auto [it,inserted]=declared_sizes.emplace(cl,declared);if(!inserted&&it->second!=declared)throw std::runtime_error("inconsistent declared cluster_size for "+cl);groups[cl].push_back(id);}
  for(const auto&[cl,v]:groups)if(declared_sizes.at(cl)!=v.size())throw std::runtime_error("cluster_size mismatch for "+cl);
  std::set<std::string> unaggregated_ids;for(const auto&id:graph_nodes)if(!assigned.count(id))unaggregated_ids.insert(id);
  std::set<std::string> partition=assigned;partition.insert(unaggregated_ids.begin(),unaggregated_ids.end());
  if(partition!=graph_nodes||assigned.size()+unaggregated_ids.size()!=graph_nodes.size())throw std::runtime_error("assigned/unaggregated partition does not exactly close over quality-pass graph nodes");
  fs::path unaggregated=ld/"unaggregated_sags.tsv";write_generated(unaggregated,c,[&](std::ostream&ua){ua<<"sag_id\tassembly_fasta\treason\n";for(const auto&id:unaggregated_ids)ua<<id<<"\t"<<eligible.at(id).assembly.string()<<"\tnot_in_chosen_size_ge10_cluster\n";});
  const auto assigned_count=assigned.size(),unaggregated_count=unaggregated_ids.size(),group_count=groups.size(),positive_edge_count=pos.size();assigned.clear();unaggregated_ids.clear();partition.clear();declared_sizes.clear();pos.clear();agg.clear();
  fs::path cd=c.out/"06_subassemble",cluster_manifest=cd/"clusters.tsv";std::vector<SubassJob> subass_jobs;subass_jobs.reserve(groups.size());
  for(auto&[cl,v]:groups){const auto n=v.size();auto d=cd/cl,input=d/"input_subassemblies.fasta",run=d/"run";materialize_merged_fasta(input,v,eligible,c);subass_jobs.push_back({cl,input,run,n});std::vector<std::string>().swap(v);}
  run_subass_jobs(subass_jobs,c);
  write_generated(cluster_manifest,c,[&](std::ostream&cm){cm<<"cluster\tn_sags\tn_found\tinput_fasta\tbin_fasta\n";for(const auto&j:subass_jobs)cm<<j.cluster<<"\t"<<j.n<<"\t"<<j.n<<"\t"<<j.input.string()<<"\t"<<(j.run/"assembly.fasta").string()<<"\n";});
  std::ostringstream done;done<<"{\n  \"schema\": \"sag-stage3b-complete-v2\",\n  \"status\": \"PASS\",\n  \"pipeline\": \"cellbit-negative-sag-stage3b-tractor-v0.2\",\n  \"entry_rule\": \"Cellbit/Dna2bit negative && max_contig>=1000 && CheckM2 contamination<5\",\n  \"reference_search\": \"disabled\",\n  \"reference_positive_exclusion\": false,\n  \"ani_engine\": \"experimental-gtdb-ani-af\",\n  \"ani_engine_path\": \""<<json(c.ani_engine.string())<<"\",\n  \"ani_engine_sha256\": \""<<sha256_file(c.ani_engine)<<"\",\n  \"makeblastdb_path\": \""<<json(c.makeblastdb.string())<<"\",\n  \"makeblastdb_sha256\": \""<<sha256_file(c.makeblastdb)<<"\",\n  \"blastn_path\": \""<<json(c.blastn.string())<<"\",\n  \"blastn_sha256\": \""<<sha256_file(c.blastn)<<"\",\n  \"python_path\": \""<<json(c.python.string())<<"\",\n  \"python_sha256\": \""<<sha256_file(c.python)<<"\",\n  \"leiden_backend_path\": \""<<json(c.leiden_backend.string())<<"\",\n  \"leiden_backend_sha256\": \""<<sha256_file(c.leiden_backend)<<"\",\n  \"subass_path\": \""<<json(c.subass.string())<<"\",\n  \"subass_sha256\": \""<<sha256_file(c.subass)<<"\",\n  \"flye_provenance_fingerprint_sha256\": \""<<file_set_fingerprint(c.flye_provenance)<<"\",\n  \"triangle_mode\": \"exact_all_pairs\",\n  \"candidate_or_sketch_prefilter\": false,\n  \"triangle_min_af\": 15,\n  \"triangle_stats_sha256\": \""<<sha256_file(tri_stats)<<"\",\n  \"positive_edge_rule\": \"ANI>=95 && abs(delta_GC)<=2; triangle output requires max(AF_ref,AF_query)>=15\",\n  \"cellbit_negative_input\": "<<sags.size()<<",\n  \"quality_failed\": "<<(sags.size()-eligible.size())<<",\n  \"quality_pass_graph_nodes\": "<<graph_node_count<<",\n  \"pairwise_pairs_requested\": "<<requested_pairs<<",\n  \"pairwise_pairs_evaluated\": "<<evaluated_pairs<<",\n  \"triangle_rows_emitted\": "<<triangle_rows<<",\n  \"positive_edges\": "<<positive_edge_count<<",\n  \"marker_jobs\": "<<jobs<<",\n  \"negative_edges\": "<<neg_n<<",\n  \"clusters_ge10\": "<<group_count<<",\n  \"assigned_sags\": "<<assigned_count<<",\n  \"unaggregated_sags\": "<<unaggregated_count<<",\n  \"assigned_plus_unaggregated_equals_quality_pass\": true\n}\n";write_or_verify(c.out/"COMPLETE.json",done.str(),c);
  std::cerr<<"optimization metrics: tsv_rows="<<metrics.tsv_rows<<" sha_cache_hits="<<metrics.sha_cache_hits<<" sha_cache_misses="<<metrics.sha_cache_misses<<" sha_bytes="<<metrics.sha_bytes_read<<" fasta_bytes_single_pass="<<metrics.fasta_bytes_read<<" marker_peak_jobs="<<metrics.marker_peak_jobs<<" subass_peak_jobs="<<metrics.subass_peak_jobs<<"\n";
  std::cout<<"PASS Stage3B: quality_pass_graph_nodes="<<graph_node_count<<" clusters_ge10="<<group_count<<"\n";return 0;
}catch(const std::exception&e){std::cerr<<"fatal: "<<e.what()<<"\n";return 2;}
