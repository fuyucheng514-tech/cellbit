#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <process.h>
#include <intrin.h>
#else
#include <dirent.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr std::size_t kProductionRowBytes = 6912;
constexpr std::size_t kGtdb232ReferenceCount = 199923;
constexpr const char* kManifestSchema = "dna2bit-packed-manifest-v1";
constexpr const char* kPackName = "references.pack";
constexpr const char* kManifestName = "references.tsv";
constexpr const char* kShaName = "references.pack.sha256";
constexpr const char* kManifestShaName = "references.tsv.sha256";
constexpr const char* kTaxonomyName = "genome_taxonomy.csv";
constexpr const char* kTaxonomyShaName = "genome_taxonomy.csv.sha256";
constexpr const char* kCompleteName = "COMPLETE.json";

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::uint64_t process_id() {
#ifdef _WIN32
  return static_cast<std::uint64_t>(_getpid());
#else
  return static_cast<std::uint64_t>(::getpid());
#endif
}

bool has_trailing_separator(const std::string& path) {
  if (path.empty()) return false;
  const char last = path.back();
  return last == '/' || last == '\\';
}

void require_legacy_directory_argument(const std::string& path, const char* flag) {
  if (!has_trailing_separator(path)) {
    fail(std::string(flag) + " must end in '/' exactly as required by teacher Dna2bit: " + path);
  }
}

bool unsafe_manifest_text(const std::string& value) {
  return value.find_first_of("\t\r\n\0") != std::string::npos;
}

std::string json_escape(const std::string& value) {
  std::ostringstream out;
  for (unsigned char c : value) {
    switch (c) {
      case '\\': out << "\\\\"; break;
      case '"': out << "\\\""; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<unsigned>(c) << std::dec;
        } else {
          out << static_cast<char>(c);
        }
    }
  }
  return out.str();
}

// Small self-contained SHA-256 implementation.  The index builder hashes the
// exact byte stream while writing it; search verifies that digest by default.
// This keeps a same-sized corrupted or accidentally replaced pack fail-closed.
class Sha256 {
 public:
  Sha256() { reset(); }

  void update(const void* bytes, std::size_t length) {
    const auto* input = static_cast<const std::uint8_t*>(bytes);
    total_bytes_ += length;
    while (length != 0) {
      const std::size_t take = std::min(length, block_.size() - buffered_);
      std::memcpy(block_.data() + buffered_, input, take);
      buffered_ += take;
      input += take;
      length -= take;
      if (buffered_ == block_.size()) {
        transform(block_.data());
        buffered_ = 0;
      }
    }
  }

  std::string final_hex() {
    const std::uint64_t message_bits = total_bytes_ * 8;
    block_[buffered_++] = 0x80;
    if (buffered_ > 56) {
      std::fill(block_.begin() + static_cast<std::ptrdiff_t>(buffered_), block_.end(), 0);
      transform(block_.data());
      buffered_ = 0;
    }
    std::fill(block_.begin() + static_cast<std::ptrdiff_t>(buffered_), block_.begin() + 56, 0);
    for (unsigned i = 0; i < 8; ++i) {
      block_[63 - i] = static_cast<std::uint8_t>(message_bits >> (8 * i));
    }
    transform(block_.data());

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::uint32_t word : state_) out << std::setw(8) << word;
    const std::string digest = out.str();
    reset();
    return digest;
  }

 private:
  static constexpr std::array<std::uint32_t, 64> k_ = {
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
      0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
      0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
      0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
      0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  static std::uint32_t rotate_right(std::uint32_t value, unsigned count) {
    return (value >> count) | (value << (32 - count));
  }

  void reset() {
    state_ = {0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
              0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U};
    total_bytes_ = 0;
    buffered_ = 0;
    block_.fill(0);
  }

  void transform(const std::uint8_t* block) {
    std::array<std::uint32_t, 64> words{};
    for (unsigned i = 0; i < 16; ++i) {
      words[i] = (static_cast<std::uint32_t>(block[4 * i]) << 24) |
                 (static_cast<std::uint32_t>(block[4 * i + 1]) << 16) |
                 (static_cast<std::uint32_t>(block[4 * i + 2]) << 8) |
                 static_cast<std::uint32_t>(block[4 * i + 3]);
    }
    for (unsigned i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotate_right(words[i - 15], 7) ^
                               rotate_right(words[i - 15], 18) ^ (words[i - 15] >> 3);
      const std::uint32_t s1 = rotate_right(words[i - 2], 17) ^
                               rotate_right(words[i - 2], 19) ^ (words[i - 2] >> 10);
      words[i] = words[i - 16] + s0 + words[i - 7] + s1;
    }

    std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
    std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
    for (unsigned i = 0; i < 64; ++i) {
      const std::uint32_t sum1 = rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
      const std::uint32_t choose = (e & f) ^ (~e & g);
      const std::uint32_t temp1 = h + sum1 + choose + k_[i] + words[i];
      const std::uint32_t sum0 = rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
  }

  std::array<std::uint32_t, 8> state_{};
  std::array<std::uint8_t, 64> block_{};
  std::uint64_t total_bytes_ = 0;
  std::size_t buffered_ = 0;
};

std::string sha256_of_file(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) fail("cannot open file for SHA-256: " + path.string());
  Sha256 hash;
  std::array<char, 1U << 20> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    if (count > 0) hash.update(buffer.data(), static_cast<std::size_t>(count));
  }
  if (!input.eof()) fail("error hashing file: " + path.string());
  return hash.final_hex();
}

struct LegacyPath {
  std::string basename;
  std::string full_path;
};

std::vector<LegacyPath> legacy_readdir(const std::string& directory) {
  require_legacy_directory_argument(directory, "directory");
  std::vector<LegacyPath> result;
#ifdef _WIN32
  // This fallback exists for local smoke builds.  Production equivalence is
  // defined against Linux opendir/readdir, which is what the teacher binary
  // itself calls.  FindFirstFile order is intentionally left unsorted.
  WIN32_FIND_DATAA data{};
  HANDLE handle = FindFirstFileA((directory + "*").c_str(), &data);
  if (handle == INVALID_HANDLE_VALUE) fail("cannot open directory: " + directory);
  do {
    const std::string name = data.cFileName;
    if (name != "." && name != "..") result.push_back({name, directory + name});
  } while (FindNextFileA(handle, &data));
  FindClose(handle);
#else
  DIR* dir = ::opendir(directory.c_str());
  if (!dir) fail("cannot open directory " + directory + ": " + std::strerror(errno));
  errno = 0;
  while (dirent* entry = ::readdir(dir)) {
    const std::string name = entry->d_name;
    if (name != "." && name != "..") result.push_back({name, directory + name});
    errno = 0;
  }
  const int read_errno = errno;
  ::closedir(dir);
  if (read_errno != 0) fail("readdir failed for " + directory + ": " + std::strerror(read_errno));
#endif
  return result;
}

std::string extract_accession(const std::string& path) {
  // Exact pattern embedded in teacher binary at .rodata 0x304c0.
  static const std::regex pattern(R"(.*(GC[FA]_\d+\.\d+)_)");
  std::smatch match;
  if (!std::regex_search(path, match, pattern) || match.size() < 2) {
    fail("reference filename does not match teacher accession regex: " + path);
  }
  return match[1].str();
}

std::uint64_t regular_file_size(const fs::path& path) {
  std::error_code error;
  if (!fs::is_regular_file(path, error) || error) fail("not a regular file: " + path.string());
  const std::uint64_t size = fs::file_size(path, error);
  if (error) fail("cannot stat file: " + path.string());
  return size;
}

struct IndexEntry {
  std::string basename;
  std::string accession;
};

struct IndexManifest {
  std::size_t row_bytes = 0;
  std::vector<IndexEntry> entries;
};

struct BoundTaxonomy {
  bool present = false;
  std::size_t entry_count = 0;
  std::string sha256;
};

std::vector<std::string> split_tabs(const std::string& line) {
  std::vector<std::string> fields;
  std::size_t begin = 0;
  for (;;) {
    const std::size_t end = line.find('\t', begin);
    if (end == std::string::npos) {
      fields.push_back(line.substr(begin));
      return fields;
    }
    fields.push_back(line.substr(begin, end - begin));
    begin = end + 1;
  }
}

std::size_t parse_size(const std::string& value, const std::string& label) {
  std::size_t used = 0;
  unsigned long long parsed = 0;
  try {
    parsed = std::stoull(value, &used);
  } catch (...) {
    fail("invalid " + label + ": " + value);
  }
  if (used != value.size() || parsed > std::numeric_limits<std::size_t>::max()) {
    fail("invalid " + label + ": " + value);
  }
  return static_cast<std::size_t>(parsed);
}

BoundTaxonomy copy_and_validate_taxonomy(const fs::path& source, const fs::path& destination,
                                         const std::vector<IndexEntry>& references) {
  if (!fs::is_regular_file(source)) fail("taxonomy is not a regular file: " + source.string());

  std::unordered_set<std::string> unmatched;
  unmatched.reserve(references.size() * 2);
  for (const IndexEntry& reference : references) {
    if (!unmatched.insert(reference.accession).second) {
      fail("duplicate accession in reference inventory: " + reference.accession);
    }
  }

  // Preserve the exact taxonomy byte stream.  The deployed teacher reader
  // observably retains CR from CRLF records in output column 3.
  std::ifstream input(source, std::ios::binary);
  std::ofstream output(destination, std::ios::binary | std::ios::trunc);
  if (!input) fail("cannot open taxonomy: " + source.string());
  if (!output) fail("cannot create bound taxonomy: " + destination.string());
  Sha256 digest;
  std::array<char, 1U << 20> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    if (count > 0) {
      output.write(buffer.data(), count);
      if (!output) fail("failed writing bound taxonomy: " + destination.string());
      digest.update(buffer.data(), static_cast<std::size_t>(count));
    }
  }
  if (!input.eof()) fail("error reading taxonomy: " + source.string());
  output.flush();
  if (!output) fail("failed flushing bound taxonomy: " + destination.string());
  output.close();
  const std::string taxonomy_sha = digest.final_hex();

  // A production index is valid only when taxonomy and packed references are
  // the exact same accession universe.  This catches the observed failure
  // mode where an older 113,104-row table was paired with 199,923 GTDB232
  // reference bits, while leaving taxonomy out of the distance computation.
  std::ifstream validate(destination, std::ios::binary);
  if (!validate) fail("cannot reopen bound taxonomy: " + destination.string());
  std::string line;
  std::size_t count = 0;
  while (std::getline(validate, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) fail("empty row in taxonomy: " + source.string());
    const std::size_t comma = line.find(',');
    if (comma == std::string::npos || comma == 0 || comma + 1 == line.size()) {
      fail("invalid accession,taxonomy row in: " + source.string());
    }
    const std::string accession = line.substr(0, comma);
    if (unmatched.erase(accession) != 1) {
      fail("taxonomy has duplicate or non-reference accession: " + accession);
    }
    ++count;
  }
  if (!validate.eof()) fail("error validating taxonomy: " + destination.string());
  if (!unmatched.empty()) {
    fail("taxonomy is missing " + std::to_string(unmatched.size()) +
         " packed reference accessions (example: " + *unmatched.begin() + ")");
  }
  if (count != references.size()) {
    fail("taxonomy/reference entry count mismatch");
  }
  return {true, count, taxonomy_sha};
}

void build_index(const std::string& database, const fs::path& output,
                 std::size_t expected_row_bytes, const fs::path& taxonomy) {
  require_legacy_directory_argument(database, "--database");
  if (fs::exists(output)) fail("index destination already exists (write-once): " + output.string());
  const fs::path staging = fs::path(output.string() + ".incomplete-" + std::to_string(process_id()));
  if (fs::exists(staging)) fail("staging destination already exists: " + staging.string());
  fs::create_directories(staging);

  const std::vector<LegacyPath> paths = legacy_readdir(database);
  if (paths.empty()) fail("database directory is empty: " + database);
  if (taxonomy.empty() && paths.size() == kGtdb232ReferenceCount) {
    fail("a 199,923-reference GTDB232 index must bind its matching taxonomy; pass --taxonomy");
  }

  const fs::path pack_path = staging / kPackName;
  const fs::path manifest_path = staging / kManifestName;
  std::ofstream pack(pack_path, std::ios::binary | std::ios::trunc);
  if (!pack) fail("cannot create pack: " + pack_path.string());

  std::vector<IndexEntry> entries;
  entries.reserve(paths.size());
  Sha256 digest;
  std::array<char, 1U << 20> buffer{};
  std::size_t row_bytes = 0;
  for (std::size_t index = 0; index < paths.size(); ++index) {
    const LegacyPath& item = paths[index];
    if (unsafe_manifest_text(item.basename)) fail("unsafe reference filename: " + item.basename);
    const std::uint64_t size64 = regular_file_size(item.full_path);
    if (size64 == 0 || size64 % sizeof(std::uint64_t) != 0 ||
        size64 > std::numeric_limits<std::size_t>::max()) {
      fail("reference bit file must contain a positive whole number of uint64 words: " + item.full_path);
    }
    const std::size_t size = static_cast<std::size_t>(size64);
    if (index == 0) row_bytes = size;
    if (size != row_bytes) fail("reference bit files have unequal byte lengths: " + item.full_path);
    if (expected_row_bytes != 0 && size != expected_row_bytes) {
      fail("reference bit file has " + std::to_string(size) + " bytes, expected " +
           std::to_string(expected_row_bytes) + ": " + item.full_path);
    }

    std::ifstream input(item.full_path, std::ios::binary);
    if (!input) fail("cannot open reference bit file: " + item.full_path);
    std::size_t copied = 0;
    while (copied < size) {
      const std::size_t want = std::min(buffer.size(), size - copied);
      input.read(buffer.data(), static_cast<std::streamsize>(want));
      if (input.gcount() != static_cast<std::streamsize>(want)) {
        fail("short read from reference bit file: " + item.full_path);
      }
      pack.write(buffer.data(), static_cast<std::streamsize>(want));
      if (!pack) fail("failed writing pack: " + pack_path.string());
      digest.update(buffer.data(), want);
      copied += want;
    }
    char extra = 0;
    if (input.read(&extra, 1)) fail("reference grew while packing: " + item.full_path);
    entries.push_back({item.basename, extract_accession(item.full_path)});
  }
  pack.flush();
  if (!pack) fail("failed flushing pack: " + pack_path.string());
  pack.close();
  const std::string pack_sha = digest.final_hex();

  std::ofstream manifest(manifest_path, std::ios::binary | std::ios::trunc);
  if (!manifest) fail("cannot create manifest: " + manifest_path.string());
  manifest << '#' << kManifestSchema << '\n';
  manifest << "#row_bytes\t" << row_bytes << '\n';
  manifest << "#entry_count\t" << entries.size() << '\n';
  manifest << "index\tbasename\taccession\n";
  for (std::size_t i = 0; i < entries.size(); ++i) {
    manifest << i << '\t' << entries[i].basename << '\t' << entries[i].accession << '\n';
  }
  manifest.flush();
  if (!manifest) fail("failed writing manifest: " + manifest_path.string());
  manifest.close();
  const std::string manifest_sha = sha256_of_file(manifest_path);

  BoundTaxonomy bound_taxonomy;
  if (!taxonomy.empty()) {
    bound_taxonomy = copy_and_validate_taxonomy(taxonomy, staging / kTaxonomyName, entries);
  }

  {
    std::ofstream sha_file(staging / kShaName, std::ios::binary | std::ios::trunc);
    if (!sha_file) fail("cannot create pack digest file");
    sha_file << pack_sha << "  " << kPackName << '\n';
  }
  {
    std::ofstream sha_file(staging / kManifestShaName, std::ios::binary | std::ios::trunc);
    if (!sha_file) fail("cannot create manifest digest file");
    sha_file << manifest_sha << "  " << kManifestName << '\n';
  }
  if (bound_taxonomy.present) {
    std::ofstream sha_file(staging / kTaxonomyShaName, std::ios::binary | std::ios::trunc);
    if (!sha_file) fail("cannot create taxonomy digest file");
    sha_file << bound_taxonomy.sha256 << "  " << kTaxonomyName << '\n';
  }
  {
    std::ofstream complete(staging / kCompleteName, std::ios::binary | std::ios::trunc);
    if (!complete) fail("cannot create completion receipt");
    complete << "{\n"
             << "  \"schema\": \"dna2bit-packed-index-complete-v1\",\n"
             << "  \"status\": \"PASS\",\n"
             << "  \"source_directory\": \"" << json_escape(database) << "\",\n"
             << "  \"entry_count\": " << entries.size() << ",\n"
             << "  \"row_bytes\": " << row_bytes << ",\n"
             << "  \"payload_bytes\": " << (entries.size() * row_bytes) << ",\n"
             << "  \"pack_sha256\": \"" << pack_sha << "\",\n"
             << "  \"manifest_sha256\": \"" << manifest_sha << "\",\n"
             << "  \"taxonomy_bound\": " << (bound_taxonomy.present ? "true" : "false") << ",\n"
             << "  \"taxonomy_entry_count\": " << bound_taxonomy.entry_count << ",\n"
             << "  \"taxonomy_sha256\": \"" << bound_taxonomy.sha256 << "\"\n"
             << "}\n";
  }

  fs::rename(staging, output);
  std::cerr << "packed " << entries.size() << " references x " << row_bytes << " bytes into "
            << output << " (sha256 " << pack_sha << ")\n";
}

IndexManifest load_manifest(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) fail("cannot open manifest: " + path.string());
  std::string line;
  if (!std::getline(input, line) || line != std::string("#") + kManifestSchema) {
    fail("unsupported or missing packed manifest schema: " + path.string());
  }
  if (!std::getline(input, line) || line.rfind("#row_bytes\t", 0) != 0) {
    fail("missing row_bytes in manifest: " + path.string());
  }
  IndexManifest manifest;
  manifest.row_bytes = parse_size(line.substr(11), "row_bytes");
  if (manifest.row_bytes == 0 || manifest.row_bytes % sizeof(std::uint64_t) != 0) {
    fail("manifest row_bytes is not a positive multiple of 8");
  }
  if (!std::getline(input, line) || line.rfind("#entry_count\t", 0) != 0) {
    fail("missing entry_count in manifest: " + path.string());
  }
  const std::size_t expected_count = parse_size(line.substr(13), "entry_count");
  if (!std::getline(input, line) || line != "index\tbasename\taccession") {
    fail("invalid manifest table header: " + path.string());
  }
  manifest.entries.reserve(expected_count);
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto fields = split_tabs(line);
    if (fields.size() != 3) fail("invalid manifest row: " + line);
    const std::size_t index = parse_size(fields[0], "manifest index");
    if (index != manifest.entries.size()) fail("manifest indices are not contiguous and ordered");
    if (unsafe_manifest_text(fields[1]) || unsafe_manifest_text(fields[2])) {
      fail("unsafe text in manifest row");
    }
    if (extract_accession(fields[1]) != fields[2]) {
      fail("manifest accession no longer matches teacher regex for " + fields[1]);
    }
    manifest.entries.push_back({fields[1], fields[2]});
  }
  if (manifest.entries.size() != expected_count) {
    fail("manifest entry_count mismatch: declared " + std::to_string(expected_count) +
         ", read " + std::to_string(manifest.entries.size()));
  }
  return manifest;
}

std::string read_expected_sha(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) fail("cannot open pack digest: " + path.string());
  std::string digest;
  input >> digest;
  if (digest.size() != 64 || !std::all_of(digest.begin(), digest.end(), [](unsigned char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
      })) {
    fail("invalid SHA-256 digest file: " + path.string());
  }
  return digest;
}

class MappedFile {
 public:
  explicit MappedFile(const fs::path& path) {
#ifdef _WIN32
    file_ = CreateFileW(path.wstring().c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                        FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file_ == INVALID_HANDLE_VALUE) fail("cannot open pack: " + path.string());
    LARGE_INTEGER size{};
    if (!GetFileSizeEx(file_, &size) || size.QuadPart <= 0) fail("cannot size pack: " + path.string());
    size_ = static_cast<std::size_t>(size.QuadPart);
    mapping_ = CreateFileMappingW(file_, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (!mapping_) fail("cannot create pack mapping: " + path.string());
    data_ = static_cast<const std::uint8_t*>(MapViewOfFile(mapping_, FILE_MAP_READ, 0, 0, 0));
    if (!data_) fail("cannot map pack: " + path.string());
#else
    fd_ = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0) fail("cannot open pack " + path.string() + ": " + std::strerror(errno));
    struct stat status {};
    if (::fstat(fd_, &status) != 0 || status.st_size <= 0) fail("cannot size pack: " + path.string());
    size_ = static_cast<std::size_t>(status.st_size);
    void* mapped = ::mmap(nullptr, size_, PROT_READ, MAP_SHARED, fd_, 0);
    if (mapped == MAP_FAILED) fail("cannot mmap pack " + path.string() + ": " + std::strerror(errno));
    data_ = static_cast<const std::uint8_t*>(mapped);
#ifdef MADV_SEQUENTIAL
    ::madvise(const_cast<std::uint8_t*>(data_), size_, MADV_SEQUENTIAL);
#endif
#endif
  }

  ~MappedFile() {
#ifdef _WIN32
    if (data_) UnmapViewOfFile(data_);
    if (mapping_) CloseHandle(mapping_);
    if (file_ != INVALID_HANDLE_VALUE) CloseHandle(file_);
#else
    if (data_) ::munmap(const_cast<std::uint8_t*>(data_), size_);
    if (fd_ >= 0) ::close(fd_);
#endif
  }

  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;
  const std::uint8_t* data() const { return data_; }
  std::size_t size() const { return size_; }

 private:
  const std::uint8_t* data_ = nullptr;
  std::size_t size_ = 0;
#ifdef _WIN32
  HANDLE file_ = INVALID_HANDLE_VALUE;
  HANDLE mapping_ = nullptr;
#else
  int fd_ = -1;
#endif
};

std::string sha256_of_mapping(const MappedFile& mapping) {
  Sha256 hash;
  constexpr std::size_t chunk = 1U << 24;
  std::size_t offset = 0;
  while (offset < mapping.size()) {
    const std::size_t take = std::min(chunk, mapping.size() - offset);
    hash.update(mapping.data() + offset, take);
    offset += take;
  }
  return hash.final_hex();
}

std::vector<std::uint64_t> load_bit_words(const std::string& path, std::size_t expected_bytes) {
  const std::uint64_t size64 = regular_file_size(path);
  if (size64 != expected_bytes) {
    fail("query/reference bit byte-length mismatch: " + path + " has " + std::to_string(size64) +
         ", expected " + std::to_string(expected_bytes));
  }
  std::vector<std::uint64_t> words(expected_bytes / sizeof(std::uint64_t));
  std::ifstream input(path, std::ios::binary);
  if (!input) fail("cannot open query bit file: " + path);
  input.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(expected_bytes));
  if (input.gcount() != static_cast<std::streamsize>(expected_bytes)) {
    fail("short read from query bit file: " + path);
  }
  return words;
}

std::unordered_map<std::string, std::string> load_teacher_taxonomy(const fs::path& path) {
  // Teacher classify::load uses fgets(key,16), fgetc(comma), then
  // fgets(value,1024) and removes exactly the final byte.  In the deployed
  // CRLF taxonomy this intentionally retains '\r'; save() adds '\n', yielding
  // byte-identical CRLF output.
  std::FILE* file = std::fopen(path.string().c_str(), "rb");
  if (!file) fail("cannot open taxonomy file: " + path.string());
  std::unordered_map<std::string, std::string> result;
  std::array<char, 16> key{};
  std::array<char, 1024> value{};
  while (std::fgets(key.data(), static_cast<int>(key.size()), file)) {
    const int delimiter = std::fgetc(file);
    if (delimiter == EOF || !std::fgets(value.data(), static_cast<int>(value.size()), file)) {
      std::fclose(file);
      fail("truncated taxonomy row in: " + path.string());
    }
    std::string taxonomy(value.data());
    if (!taxonomy.empty()) taxonomy.pop_back();
    result[std::string(key.data())] = std::move(taxonomy);
  }
  if (std::ferror(file)) {
    std::fclose(file);
    fail("error reading taxonomy file: " + path.string());
  }
  std::fclose(file);
  return result;
}

unsigned popcount64(std::uint64_t value) {
#if defined(_MSC_VER) && defined(_M_X64)
  return static_cast<unsigned>(__popcnt64(value));
#elif defined(__GNUC__) || defined(__clang__)
  return static_cast<unsigned>(__builtin_popcountll(value));
#else
  unsigned count = 0;
  while (value != 0) {
    value &= value - 1;
    ++count;
  }
  return count;
#endif
}

struct SearchResult {
  int reference_index = -1;
  int best_distance = INT_MAX;
  int legacy_second = 0;
};

[[maybe_unused]] SearchResult legacy_search_one(const std::vector<std::uint64_t>& query,
                               const std::uint8_t* pack,
                               std::size_t reference_count,
                               std::size_t row_bytes,
                               float min_ratio) {
  const std::size_t words_per_row = row_bytes / sizeof(std::uint64_t);
  int best = INT_MAX;
  int second = 0;
  int winner = -1;
  for (std::size_t reference = 0; reference < reference_count; ++reference) {
    const auto* row = reinterpret_cast<const std::uint64_t*>(pack + reference * row_bytes);
    int distance = 0;
    for (std::size_t word = 0; word < words_per_row; ++word) {
      distance += static_cast<int>(popcount64(query[word] ^ row[word]));
    }

    // Deliberately preserve the teacher binary's observable quirk: second is
    // only the previous record minimum when a STRICT new minimum is found.
    // It is not the globally second-smallest distance.  Strict comparison
    // also makes the earliest raw-readdir reference win a tie.
    if (distance < best) {
      second = best;
      best = distance;
      winner = static_cast<int>(reference);
    }
  }

  if (winner != -1) {
    const float ratio = static_cast<float>(second - best) / static_cast<float>(best);
    // Teacher rejects margin <= threshold; equivalently acceptance is strict
    // `ratio > min_ratio`.  Keep this float (not double) expression intact.
    if (!(ratio > min_ratio)) winner = -1;
  }
  return {winner, best, second};
}

int hamming_distance_bounded(const std::vector<std::uint64_t>& query,
                             const std::uint64_t* row,
                             std::size_t words_per_row,
                             int strict_upper_bound) {
  int distance = 0;
  std::size_t word = 0;
  // Eight independent POPCNT inputs let Zen3 overlap their execution.  The
  // bound is checked every 64 words; once the partial sum is >= current best,
  // the remaining non-negative terms cannot satisfy the teacher's strict
  // `distance < best` update and are therefore irrelevant.
  for (; word + 64 <= words_per_row; word += 64) {
    for (std::size_t offset = 0; offset < 64; offset += 8) {
      distance += static_cast<int>(popcount64(query[word + offset + 0] ^ row[word + offset + 0]));
      distance += static_cast<int>(popcount64(query[word + offset + 1] ^ row[word + offset + 1]));
      distance += static_cast<int>(popcount64(query[word + offset + 2] ^ row[word + offset + 2]));
      distance += static_cast<int>(popcount64(query[word + offset + 3] ^ row[word + offset + 3]));
      distance += static_cast<int>(popcount64(query[word + offset + 4] ^ row[word + offset + 4]));
      distance += static_cast<int>(popcount64(query[word + offset + 5] ^ row[word + offset + 5]));
      distance += static_cast<int>(popcount64(query[word + offset + 6] ^ row[word + offset + 6]));
      distance += static_cast<int>(popcount64(query[word + offset + 7] ^ row[word + offset + 7]));
    }
    if (distance >= strict_upper_bound) return strict_upper_bound;
  }
  for (; word < words_per_row; ++word) {
    distance += static_cast<int>(popcount64(query[word] ^ row[word]));
  }
  return distance;
}

void apply_legacy_margin(SearchResult& result, float min_ratio) {
  if (result.reference_index == -1) return;
  const float ratio = static_cast<float>(result.legacy_second - result.best_distance) /
                      static_cast<float>(result.best_distance);
  if (!(ratio > min_ratio)) result.reference_index = -1;
}

struct SearchOptions {
  fs::path index;
  std::string queries;
  fs::path taxonomy;
  fs::path output = "search_result.csv";
  float min_ratio = 0.01F;
  unsigned threads = 1;
  unsigned query_block_size = 16;
  bool verify_pack_sha = true;
};

void search_index(const SearchOptions& options) {
  require_legacy_directory_argument(options.queries, "--search");
  if (options.threads == 0) fail("--nthreads must be positive");
  if (!fs::is_regular_file(options.index / kCompleteName)) {
    fail("packed index has no PASS completion receipt: " + options.index.string());
  }

  {
    std::ifstream complete(options.index / kCompleteName, std::ios::binary);
    const std::string receipt((std::istreambuf_iterator<char>(complete)),
                              std::istreambuf_iterator<char>());
    if (complete.bad() || receipt.find("\"status\": \"PASS\"") == std::string::npos ||
        receipt.find("\"schema\": \"dna2bit-packed-index-complete-v1\"") == std::string::npos) {
      fail("packed index completion receipt is invalid or not PASS");
    }
  }
  if (options.verify_pack_sha) {
    const std::string expected_manifest = read_expected_sha(options.index / kManifestShaName);
    const std::string observed_manifest = sha256_of_file(options.index / kManifestName);
    if (observed_manifest != expected_manifest) {
      fail("manifest SHA-256 mismatch: expected " + expected_manifest + ", got " + observed_manifest);
    }
  }

  const IndexManifest manifest = load_manifest(options.index / kManifestName);
  if (manifest.entries.empty()) fail("packed index contains no references");
  if (manifest.entries.size() > static_cast<std::size_t>(INT_MAX)) {
    fail("packed index has too many references for teacher-compatible int result indices");
  }
  const MappedFile pack(options.index / kPackName);
  if (manifest.row_bytes > std::numeric_limits<std::size_t>::max() / manifest.entries.size() ||
      pack.size() != manifest.row_bytes * manifest.entries.size()) {
    fail("pack size does not equal row_bytes * entry_count");
  }
  if (options.verify_pack_sha) {
    const std::string expected = read_expected_sha(options.index / kShaName);
    const std::string observed = sha256_of_mapping(pack);
    if (observed != expected) fail("pack SHA-256 mismatch: expected " + expected + ", got " + observed);
  }

  const fs::path bound_taxonomy = options.index / kTaxonomyName;
  const fs::path bound_taxonomy_sha = options.index / kTaxonomyShaName;
  const bool has_bound_taxonomy = fs::is_regular_file(bound_taxonomy);
  if (has_bound_taxonomy != fs::is_regular_file(bound_taxonomy_sha)) {
    fail("packed index has an incomplete bound-taxonomy contract");
  }
  fs::path taxonomy_path = options.taxonomy;
  if (taxonomy_path.empty()) {
    if (!has_bound_taxonomy) {
      fail("search requires --tax because this developer index has no bound taxonomy");
    }
    taxonomy_path = bound_taxonomy;
  }
  if (!fs::is_regular_file(taxonomy_path)) {
    fail("taxonomy is not a regular file: " + taxonomy_path.string());
  }
  if (has_bound_taxonomy) {
    const std::string expected = read_expected_sha(bound_taxonomy_sha);
    const std::string observed = sha256_of_file(taxonomy_path);
    if (observed != expected) {
      fail("taxonomy SHA-256 does not match the table bound to this packed index: expected " +
           expected + ", got " + observed);
    }
  }

  const std::vector<LegacyPath> query_paths = legacy_readdir(options.queries);
  std::vector<std::vector<std::uint64_t>> query_words;
  query_words.reserve(query_paths.size());
  for (const LegacyPath& query : query_paths) {
    query_words.push_back(load_bit_words(query.full_path, manifest.row_bytes));
  }
  const auto taxonomy = load_teacher_taxonomy(taxonomy_path);

  std::vector<SearchResult> results(query_words.size());
  const std::size_t block_count = (query_words.size() + options.query_block_size - 1) /
                                  options.query_block_size;
  const unsigned worker_count = std::min<unsigned>(options.threads,
      static_cast<unsigned>(std::max<std::size_t>(1, block_count)));
  std::atomic<std::size_t> next_query{0};
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (unsigned worker = 0; worker < worker_count; ++worker) {
    workers.emplace_back([&]() {
      for (;;) {
        const std::size_t begin = next_query.fetch_add(options.query_block_size);
        if (begin >= query_words.size()) return;
        const std::size_t end = std::min(query_words.size(),
                                         begin + options.query_block_size);
        for (std::size_t reference = 0; reference < manifest.entries.size(); ++reference) {
          const auto* row = reinterpret_cast<const std::uint64_t*>(
              pack.data() + reference * manifest.row_bytes);
#if defined(__GNUC__) || defined(__clang__)
          if (reference + 1 < manifest.entries.size())
            __builtin_prefetch(pack.data() + (reference + 1) * manifest.row_bytes, 0, 0);
#endif
          for (std::size_t query = begin; query < end; ++query) {
            SearchResult& state = results[query];
            const int distance = hamming_distance_bounded(query_words[query], row,
                                                           manifest.row_bytes / sizeof(std::uint64_t),
                                                           state.best_distance);
            if (distance < state.best_distance) {
              state.legacy_second = state.best_distance;
              state.best_distance = distance;
              state.reference_index = static_cast<int>(reference);
            }
          }
        }
        for (std::size_t query = begin; query < end; ++query)
          apply_legacy_margin(results[query], options.min_ratio);
      }
    });
  }
  for (auto& worker : workers) worker.join();

  std::FILE* output = std::fopen(options.output.string().c_str(), "wb");
  if (!output) fail("cannot open result file: " + options.output.string());
  for (std::size_t query = 0; query < query_paths.size(); ++query) {
    const int winner = results[query].reference_index;
    if (winner == -1) continue;
    const std::string& accession = manifest.entries[static_cast<std::size_t>(winner)].accession;
    const auto found = taxonomy.find(accession);
    const std::string empty;
    const std::string& tax = found == taxonomy.end() ? empty : found->second;
    if (std::fputs(query_paths[query].full_path.c_str(), output) < 0 ||
        std::fputc(',', output) == EOF || std::fputs(accession.c_str(), output) < 0 ||
        std::fputc(',', output) == EOF || std::fputs(tax.c_str(), output) < 0 ||
        std::fputc('\n', output) == EOF) {
      std::fclose(output);
      fail("failed writing result file: " + options.output.string());
    }
  }
  if (std::fclose(output) != 0) fail("failed closing result file: " + options.output.string());
  std::cerr << "searched " << query_paths.size() << " queries against " << manifest.entries.size()
            << " packed references; query_block_size=" << options.query_block_size
            << "; result=" << options.output << "\n";
}

void print_usage(const char* program) {
#if defined(DNA2BIT_PACKED_BUILDER_ONLY)
  std::cerr << "One-time teacher-compatible Dna2bit pack builder\n\n  " << program
            << " --database DB_DIR/ --index INDEX_DIR [--taxonomy TAXONOMY.csv]"
               " [--expected-row-bytes 6912]\n";
#elif defined(DNA2BIT_PACKED_SEARCH_ONLY)
  std::cerr << "Teacher-compatible packed Dna2bit query search\n\n  " << program
            << " --database INDEX_DIR --search QUERY_DIR/ --tax TAXONOMY.csv"
               " [--min_ratio 0.01] [--nthreads N] [--result_file out.csv]"
               " [--query-block-size 16] [--no-verify-pack-sha]\n";
#else
  std::cerr
      << "Teacher-compatible packed Dna2bit search\n\n"
      << "Build (one time; preserves raw opendir/readdir order):\n  " << program
      << " build --database DB_DIR/ --index INDEX_DIR [--taxonomy TAXONOMY.csv]"
         " [--expected-row-bytes 6912]\n\n"
      << "Search:\n  " << program
      << " search --index INDEX_DIR --search QUERY_DIR/ --tax TAXONOMY.csv"
         " [--min-ratio 0.01] [--nthreads N] [--result-file out.csv]"
         " [--query-block-size 16] [--no-verify-pack-sha]\n";
#endif
}

std::string require_value(int& index, int argc, char** argv, const std::string& option) {
  if (++index >= argc) fail("missing value after " + option);
  return argv[index];
}

int main_impl(int argc, char** argv) {
#if defined(DNA2BIT_PACKED_BUILDER_ONLY)
  const std::string command = "build";
  const int argument_begin = 1;
#elif defined(DNA2BIT_PACKED_SEARCH_ONLY)
  const std::string command = "search";
  const int argument_begin = 1;
#else
  if (argc < 2) {
    print_usage(argv[0]);
    return 2;
  }
  const std::string command = argv[1];
  const int argument_begin = 2;
#endif
  if (command == "build") {
    std::string database;
    fs::path index;
    fs::path taxonomy;
    std::size_t expected = kProductionRowBytes;
    for (int i = argument_begin; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument == "--database" || argument == "-d") database = require_value(i, argc, argv, argument);
      else if (argument == "--index" || argument == "-o") index = require_value(i, argc, argv, argument);
      else if (argument == "--taxonomy" || argument == "--tax" || argument == "-t") {
        taxonomy = require_value(i, argc, argv, argument);
      }
      else if (argument == "--expected-row-bytes") {
        expected = parse_size(require_value(i, argc, argv, argument), "expected-row-bytes");
      } else if (argument == "--help" || argument == "-h") {
        print_usage(argv[0]);
        return 0;
      } else {
        fail("unknown build option: " + argument);
      }
    }
    if (database.empty() || index.empty()) fail("build requires --database and --index");
    if (taxonomy.empty()) {
      const fs::path database_path = fs::path(database).lexically_normal();
      const fs::path sibling = database_path.parent_path() / "genome_taxonomy_1.csv";
      if (fs::is_regular_file(sibling)) {
        taxonomy = sibling;
        std::cerr << "auto-binding sibling taxonomy: " << taxonomy << '\n';
      }
    }
    build_index(database, index, expected, taxonomy);
    return 0;
  }

  if (command == "search") {
    SearchOptions options;
    for (int i = argument_begin; i < argc; ++i) {
      const std::string argument = argv[i];
      if (argument == "--index" || argument == "--database" || argument == "-d") {
        options.index = require_value(i, argc, argv, argument);
      }
      else if (argument == "--search" || argument == "-s") options.queries = require_value(i, argc, argv, argument);
      else if (argument == "--tax" || argument == "-t") options.taxonomy = require_value(i, argc, argv, argument);
      else if (argument == "--result-file" || argument == "--result_file" || argument == "-o") {
        options.output = require_value(i, argc, argv, argument);
      }
      else if (argument == "--min-ratio" || argument == "--min_ratio" || argument == "-m") {
        const std::string value = require_value(i, argc, argv, argument);
        std::size_t used = 0;
        try {
          options.min_ratio = std::stof(value, &used);
        } catch (...) {
          fail("invalid min-ratio: " + value);
        }
        if (used != value.size()) fail("invalid min-ratio: " + value);
      } else if (argument == "--nthreads" || argument == "-n") {
        const std::size_t parsed = parse_size(require_value(i, argc, argv, argument), "nthreads");
        if (parsed == 0 || parsed > std::numeric_limits<unsigned>::max()) fail("invalid nthreads");
        options.threads = static_cast<unsigned>(parsed);
      } else if (argument == "--query-block-size") {
        const std::size_t parsed = parse_size(require_value(i, argc, argv, argument), "query-block-size");
        if (parsed == 0 || parsed > 64) fail("query-block-size must be in 1..64");
        options.query_block_size = static_cast<unsigned>(parsed);
      } else if (argument == "--no-verify-pack-sha") {
        options.verify_pack_sha = false;
      } else if (argument == "--help" || argument == "-h") {
        print_usage(argv[0]);
        return 0;
      } else {
        fail("unknown search option: " + argument);
      }
    }
    if (options.index.empty() || options.queries.empty()) {
      fail("search requires --index and --search; --tax is optional only for a bound index");
    }
    search_index(options);
    return 0;
  }

  if (command == "--help" || command == "-h") {
    print_usage(argv[0]);
    return 0;
  }
  fail("unknown command: " + command);
}

}  // namespace

#ifndef DNA2BIT_PACKED_LIBRARY
int main(int argc, char** argv) {
  try {
    return main_impl(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "dna2bit-packed: ERROR: " << error.what() << '\n';
    return 1;
  }
}
#endif  // DNA2BIT_PACKED_LIBRARY
