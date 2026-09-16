#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <zlib.h>
#include "dna2bit_embedded_api.hpp"
#include "dna2bit_packed_api.hpp"

namespace fs = std::filesystem;

enum class InputKind { paired_reads, contigs };
enum class SequenceFormat { fasta, fastq };

static fs::path environment_path(const char* name) {
  const char* value = std::getenv(name);
  return value && *value ? fs::path(value) : fs::path();
}

static fs::path default_flye_root() {
  fs::path configured = environment_path("MICROSAGS_FLYE_ROOT");
  return configured.empty() ? environment_path("CONDA_PREFIX") : configured;
}

struct Sag {
  std::string id;
  InputKind input_kind = InputKind::paired_reads;
  fs::path source1, source2;
  fs::path clean_r1, clean_r2, assembly, bit;
  std::string bit_search_token;
  std::uint64_t assembly_bp = 0;
  std::uint64_t assembly_max_contig = 0;
  std::uint64_t assembly_gc = 0;
  std::uint64_t assembly_acgt = 0;
  bool assembly_bp_known = false;
};

struct Config {
  fs::path manifest, out;
  fs::path fastp = "fastp";
  fs::path spades = "spades.py";
  fs::path dna_tax = environment_path("MICROSAGS_DNA_TAX");
  std::string dna_search_engine = "packed";
  fs::path dna_packed_db = environment_path("MICROSAGS_DNA_PACKED_DB");
  fs::path subass = "cpp-subass";
  fs::path flye_root = default_flye_root();
  int threads = 24;
  int memory_gb = 128;
  bool dry = false;
  bool resume = false;
  bool stop_after_annotation = false;
};

class SequenceLineReader {
 public:
  explicit SequenceLineReader(const fs::path& path) : path_(path) {
    file_ = gzopen(path.string().c_str(), "rb");
    if (!file_) throw std::runtime_error("cannot open sequence file: " + path.string());
    // zlib defaults to a very small buffer.  Sequence inputs are streamed
    // sequentially, so a bounded 1 MiB inflate buffer materially reduces
    // syscall overhead without making memory scale with file size.
    if (gzbuffer(file_, 1U << 20) != 0) {
      gzclose(file_);
      file_ = nullptr;
      throw std::runtime_error("cannot configure gzip buffer: " + path.string());
    }
  }

  ~SequenceLineReader() {
    if (file_) gzclose(file_);
  }

  SequenceLineReader(const SequenceLineReader&) = delete;
  SequenceLineReader& operator=(const SequenceLineReader&) = delete;

  bool getline(std::string& out) {
    out.clear();
    while (true) {
      char* got = gzgets(file_, buffer_.data(), static_cast<int>(buffer_.size()));
      if (!got) {
        int error_number = Z_OK;
        const char* message = gzerror(file_, &error_number);
        if (error_number != Z_OK && error_number != Z_STREAM_END) {
          throw std::runtime_error("error reading " + path_.string() + ": " +
                                   (message ? message : "unknown zlib error"));
        }
        return !out.empty();
      }
      out += buffer_.data();
      if (!out.empty() && out.back() == '\n') {
        out.pop_back();
        if (!out.empty() && out.back() == '\r') out.pop_back();
        return true;
      }
      if (gzeof(file_)) return !out.empty();
    }
  }

 private:
  fs::path path_;
  gzFile file_ = nullptr;
  std::array<char, 65536> buffer_{};
};

static std::string q(const fs::path& p) {
  std::string result = "'";
  for (char c : p.string()) result += c == '\'' ? "'\\''" : std::string(1, c);
  return result + "'";
}

static std::string clean_id(std::string value) {
  for (char& c : value) {
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')) c = '_';
  }
  if (value.empty() || value == "." || value == "..") throw std::runtime_error("invalid empty SAG id");
  return value;
}

static void run(const std::string& command, const fs::path& receipt, bool dry) {
  fs::create_directories(receipt.parent_path());
  std::cerr << "+ " << command << "\n";
  if (dry) return;
  const int rc = std::system(command.c_str());
  if (rc != 0) throw std::runtime_error("command failed rc=" + std::to_string(rc) + ": " + command);
  std::ofstream(receipt) << "PASS\n" << command << "\n";
}

struct CommandTask {
  std::string command;
  fs::path receipt;
  std::uintmax_t work_bytes = 0;
  std::string stable_key;
};

// Execute independent external jobs without allowing the sum of their child
// thread budgets to exceed --threads.  The caller prepares all deterministic
// inputs first, then supplies tasks in canonical order.  Individual receipts
// remain exactly the normal run() receipt and are only written after that job
// exits successfully; a failed child prevents later queued children from
// starting and makes the whole stage fail closed.
static void run_bounded_parallel(const std::vector<CommandTask>& tasks, std::size_t workers, bool dry) {
  if (tasks.empty()) return;
  if (dry || workers <= 1 || tasks.size() == 1) {
    for (const auto& task : tasks) run(task.command, task.receipt, dry);
    return;
  }

  workers = std::min(workers, tasks.size());
  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::mutex error_mutex;
  std::string first_error;
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (std::size_t worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&]() {
      while (!failed.load(std::memory_order_acquire)) {
        const std::size_t task_index = next.fetch_add(1, std::memory_order_relaxed);
        if (task_index >= tasks.size()) return;
        try {
          run(tasks[task_index].command, tasks[task_index].receipt, false);
        } catch (const std::exception& error) {
          {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (first_error.empty()) first_error = error.what();
          }
          failed.store(true, std::memory_order_release);
          return;
        }
      }
    });
  }
  for (auto& worker : pool) worker.join();
  if (failed.load(std::memory_order_acquire)) {
    throw std::runtime_error("parallel external stage failed: " + first_error);
  }
}

// The two expensive child stages have different process startup and memory
// profiles.  These targets keep each child useful while deriving the number
// of simultaneous children solely from the single global --threads budget.
// Multiplying returned workers by returned threads-per-worker never exceeds
// the requested budget.
static std::pair<std::size_t, int> bounded_parallel_shape(std::size_t jobs, int total_threads,
                                                          int target_threads_per_child) {
  if (jobs == 0) return {0, 1};
  const auto wanted_workers = static_cast<std::size_t>(
      std::max(1, total_threads / std::max(1, target_threads_per_child)));
  const auto workers = std::min(jobs, wanted_workers);
  const int threads_per_worker = std::max(1, total_threads / static_cast<int>(workers));
  return {workers, threads_per_worker};
}

static std::vector<std::string> split(const std::string& value, char delimiter) {
  std::vector<std::string> values;
  if (value.empty()) return values;
  std::size_t begin = 0;
  while (true) {
    const std::size_t end = value.find(delimiter, begin);
    if (end == std::string::npos) {
      values.emplace_back(value.substr(begin));
      break;
    }
    values.emplace_back(value.substr(begin, end - begin));
    begin = end + 1;
  }
  return values;
}

static std::string lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

static void require_file(const fs::path& path, const std::string& name) {
  if (!fs::is_regular_file(path)) throw std::runtime_error(name + " not found: " + path.string());
}

static fs::path resolve_program(const fs::path& configured, const fs::path& launcher,
                                const std::string& name) {
  auto accept = [](const fs::path& candidate) -> fs::path {
    std::error_code error;
    if (!fs::is_regular_file(candidate, error) || error) return {};
    return fs::weakly_canonical(candidate, error);
  };
  if (configured.has_parent_path() || configured.is_absolute()) {
    auto found = accept(configured);
    if (!found.empty()) return found;
    throw std::runtime_error(name + " not found: " + configured.string());
  }

  // Installed binaries are siblings.  In an uninstalled CMake build the
  // bundled driver lives one subdirectory lower.
  std::error_code error;
  const auto launcher_dir = fs::absolute(launcher, error).parent_path();
  if (!error) {
    for (const auto& candidate : {launcher_dir / configured, launcher_dir / "subass" / configured}) {
      auto found = accept(candidate);
      if (!found.empty()) return found;
    }
  }

  const char* path_value = std::getenv("PATH");
  if (path_value) {
#ifdef _WIN32
    constexpr char path_separator = ';';
#else
    constexpr char path_separator = ':';
#endif
    for (const auto& directory : split(path_value, path_separator)) {
      auto found = accept(fs::path(directory) / configured);
      if (!found.empty()) return found;
    }
  }
  throw std::runtime_error(name + " not found beside the launcher or on PATH: " + configured.string());
}

static fs::path resolve_manifest_path(const fs::path& manifest, const std::string& value) {
  if (value.empty()) throw std::runtime_error("empty input path in manifest");
  fs::path path(value);
  if (path.is_relative()) path = manifest.parent_path() / path;
  require_file(path, "input sequence file");
  return fs::weakly_canonical(path);
}

static bool valid_nucleotide_line(const std::string& line) {
  bool saw_base = false;
  for (unsigned char c : line) {
    if (std::isspace(c)) continue;
    if (!(std::isalpha(c) || c == '-' || c == '.' || c == '*')) return false;
    saw_base = true;
  }
  return saw_base;
}

static SequenceFormat detect_sequence_format(const fs::path& path) {
  SequenceLineReader reader(path);
  std::string first;
  while (reader.getline(first) && first.empty()) {}
  if (first.size() >= 3 && static_cast<unsigned char>(first[0]) == 0xef &&
      static_cast<unsigned char>(first[1]) == 0xbb && static_cast<unsigned char>(first[2]) == 0xbf) {
    first.erase(0, 3);
  }
  if (first.empty()) throw std::runtime_error("empty sequence file: " + path.string());

  if (first[0] == '>') {
    bool saw_sequence = false;
    std::string line;
    std::size_t inspected = 0;
    while (inspected < 256 && reader.getline(line)) {
      if (line.empty()) continue;
      ++inspected;
      if (line[0] == '>') continue;
      if (!valid_nucleotide_line(line)) {
        throw std::runtime_error("FASTA contains a non-nucleotide sequence line near the start: " + path.string());
      }
      saw_sequence = true;
    }
    if (!saw_sequence) throw std::runtime_error("FASTA has no sequence after its header: " + path.string());
    return SequenceFormat::fasta;
  }

  if (first[0] == '@') {
    std::string header = first;
    for (int record = 0; record < 8; ++record) {
      std::string sequence, plus, quality;
      if (!reader.getline(sequence) || !reader.getline(plus) || !reader.getline(quality)) {
        throw std::runtime_error("truncated FASTQ record: " + path.string());
      }
      if (header.empty() || header[0] != '@' || plus.empty() || plus[0] != '+' || sequence.empty() ||
          sequence.size() != quality.size() || !valid_nucleotide_line(sequence)) {
        throw std::runtime_error("invalid four-line FASTQ structure: " + path.string());
      }
      if (!reader.getline(header)) break;
      if (header.empty() || header[0] != '@') {
        throw std::runtime_error("invalid FASTQ header after record " + std::to_string(record + 1) + ": " +
                                 path.string());
      }
    }
    return SequenceFormat::fastq;
  }

  throw std::runtime_error("cannot identify sequence content as FASTA or FASTQ: " + path.string());
}

static bool is_manifest_header(const std::vector<std::string>& columns) {
  if (columns.empty() || lower(columns[0]) != "sag_id") return false;
  if (columns.size() == 2) {
    const std::string second = lower(columns[1]);
    return second == "assembly_fasta" || second == "contigs" || second == "input";
  }
  if (columns.size() != 3) return false;
  const std::string second = lower(columns[1]);
  const std::string third = lower(columns[2]);
  return (second == "r1" && third == "r2") ||
         ((second == "input1" || second == "source1") && (third == "input2" || third == "source2"));
}

static std::vector<Sag> read_manifest(const fs::path& manifest) {
  std::ifstream input(manifest);
  if (!input) throw std::runtime_error("cannot read manifest: " + manifest.string());
  std::vector<Sag> sags;
  std::set<std::string> ids;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    const auto columns = split(line, '\t');
    if (is_manifest_header(columns)) continue;
    if (columns.size() != 2 && columns.size() != 3) {
      throw std::runtime_error("manifest line " + std::to_string(line_number) +
                               " must have 2 columns (SAG_ID, contigs FASTA) or 3 columns (SAG_ID, R1, R2)");
    }

    Sag sag;
    sag.id = clean_id(columns[0]);
    if (!ids.insert(sag.id).second) throw std::runtime_error("duplicate SAG id after ID normalization: " + sag.id);
    sag.source1 = resolve_manifest_path(manifest, columns[1]);
    const auto format1 = detect_sequence_format(sag.source1);

    if (columns.size() == 2) {
      if (format1 != SequenceFormat::fasta) {
        throw std::runtime_error("manifest line " + std::to_string(line_number) +
                                 " has one FASTQ input; paired reads require both R1 and R2 columns");
      }
      sag.input_kind = InputKind::contigs;
    } else {
      sag.source2 = resolve_manifest_path(manifest, columns[2]);
      const auto format2 = detect_sequence_format(sag.source2);
      if (format1 != SequenceFormat::fastq || format2 != SequenceFormat::fastq) {
        throw std::runtime_error("manifest line " + std::to_string(line_number) +
                                 " has three columns, but both R1 and R2 must be FASTQ by content");
      }
      if (fs::equivalent(sag.source1, sag.source2)) {
        throw std::runtime_error("R1 and R2 resolve to the same file for SAG " + sag.id);
      }
      sag.input_kind = InputKind::paired_reads;
    }
    sags.push_back(std::move(sag));
  }
  if (sags.empty()) throw std::runtime_error("manifest has no SAGs");
  return sags;
}

static Config parse(int argc, char** argv) {
  Config config;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    auto value = [&]() {
      if (++i >= argc) throw std::runtime_error("missing value for " + option);
      return std::string(argv[i]);
    };
    if (option == "--manifest") config.manifest = value();
    else if (option == "--out") config.out = value();
    else if (option == "--threads") config.threads = std::stoi(value());
    else if (option == "--memory-gb") config.memory_gb = std::stoi(value());
    else if (option == "--fastp") config.fastp = value();
    else if (option == "--spades") config.spades = value();
    else if (option == "--dna-tax") config.dna_tax = value();
    else if (option == "--dna-search-engine") config.dna_search_engine = value();
    else if (option == "--dna-packed-db") config.dna_packed_db = value();
    else if (option == "--subass") config.subass = value();
    else if (option == "--flye-root") config.flye_root = value();
    else if (option == "--dry-run") config.dry = true;
    else if (option == "--resume") config.resume = true;
    else if (option == "--stop-after") {
      const auto stage = value();
      if (stage != "annotation") {
        throw std::runtime_error("--stop-after currently accepts only: annotation");
      }
      config.stop_after_annotation = true;
    }
    else if (option == "--help" || option == "-h") {
      std::cout << "dna2bit-sag-pipeline --manifest SAGs.tsv --out DIR [options]\n\n"
                << "Input is detected from file contents, not filename extensions:\n"
                << "  SAG_ID<TAB>R1<TAB>R2       paired FASTQ; run fastp -> SPAdes -> >=1000-bp gate\n"
                << "  SAG_ID<TAB>contigs.fasta   FASTA; skip fastp/SPAdes -> >=1000-bp gate\n"
                << "A matching optional header (sag_id/r1/r2 or sag_id/assembly_fasta) is accepted.\n\n"
                << "Search engine is embedded teacher-compatible packed search:\n"
                << "  --dna-search-engine packed --dna-packed-db INDEX_DIR\n\n"
                << "Workflow endpoint:\n"
                << "  --stop-after annotation   write DNA2bit labels/pending files and skip Stage 3A\n\n"
                << "Portable dependency paths:\n"
                << "  --fastp PATH --spades PATH --subass PATH --flye-root PREFIX\n"
                << "  MICROSAGS_DNA_TAX, MICROSAGS_DNA_PACKED_DB and MICROSAGS_FLYE_ROOT\n"
                << "  may be used instead of repeating data/runtime paths.\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (config.manifest.empty() || config.out.empty()) throw std::runtime_error("required: --manifest and --out");
  if (config.threads < 1 || config.memory_gb < 1) throw std::runtime_error("threads and memory-gb must be positive");
  if (config.dna_search_engine != "packed") {
    throw std::runtime_error("embedded package supports only --dna-search-engine packed");
  }
  if (config.dna_search_engine == "packed" && config.dna_packed_db.empty()) {
    throw std::runtime_error("packed search requires --dna-packed-db or MICROSAGS_DNA_PACKED_DB");
  }
  if (config.dna_tax.empty()) {
    throw std::runtime_error("taxonomy requires --dna-tax or MICROSAGS_DNA_TAX");
  }
  if (!config.stop_after_annotation && config.flye_root.empty()) {
    throw std::runtime_error("Flye runtime requires --flye-root, MICROSAGS_FLYE_ROOT, or an active Conda environment");
  }
  return config;
}

static bool is_gzip(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  unsigned char bytes[2] = {0, 0};
  input.read(reinterpret_cast<char*>(bytes), 2);
  return input.gcount() == 2 && bytes[0] == 0x1f && bytes[1] == 0x8b;
}

static void ensure_symlink(const fs::path& source, const fs::path& link) {
  require_file(source, "symlink source");
  fs::create_directories(link.parent_path());
  if (fs::exists(link)) {
    if (!fs::equivalent(source, link)) {
      throw std::runtime_error("existing link target differs from required source: " + link.string());
    }
    return;
  }
  fs::create_symlink(fs::absolute(source), link);
}

struct FastaStats {
  std::uint64_t total = 0;
  std::uint64_t max_contig = 0;
  std::uint64_t gc = 0;
  std::uint64_t acgt = 0;
};

static void consume_fasta_sequence_line(const std::string& line, FastaStats& stats,
                                        std::uint64_t& current_contig) {
  for (unsigned char c : line) {
    if (std::isspace(c)) continue;
    ++stats.total;
    ++current_contig;
    const unsigned char base = static_cast<unsigned char>(std::toupper(c));
    if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
      ++stats.acgt;
      if (base == 'C' || base == 'G') ++stats.gc;
    }
  }
}

// Returns true when all FASTA statistics were obtained while a gzip input was
// already being materialized.  Plain FASTA stays a zero-copy symlink and is
// counted in the common bounded scan below.
static bool materialize_contigs(const fs::path& source, const fs::path& destination,
                                FastaStats& stats) {
  fs::create_directories(destination.parent_path());
  if (fs::exists(destination)) throw std::runtime_error("refusing to overwrite existing contig view: " + destination.string());
  if (!is_gzip(source)) {
    fs::create_symlink(fs::absolute(source), destination);
    return false;
  }

  const fs::path temporary = destination.string() + ".incomplete";
  if (fs::exists(temporary)) throw std::runtime_error("stale incomplete contig copy: " + temporary.string());
  SequenceLineReader reader(source);
  std::ofstream output(temporary, std::ios::binary);
  if (!output) throw std::runtime_error("cannot create decompressed contig view: " + temporary.string());
  std::string line;
  stats = {};
  std::uint64_t current_contig = 0;
  while (reader.getline(line)) {
    output << line << '\n';
    if (!line.empty() && line[0] == '>') {
      stats.max_contig = std::max(stats.max_contig, current_contig);
      current_contig = 0;
      continue;
    }
    consume_fasta_sequence_line(line, stats, current_contig);
  }
  stats.max_contig = std::max(stats.max_contig, current_contig);
  output.close();
  if (!output) throw std::runtime_error("failed writing decompressed contig view: " + temporary.string());
  fs::rename(temporary, destination);
  return true;
}

static std::string species_key(const std::string& taxonomy) {
  auto ranks = split(taxonomy, ';');
  std::string species = ranks.empty() ? taxonomy : ranks.back();
  if (species.rfind("s__", 0) == 0) species = species.substr(3);
  return clean_id(species);
}

static void append_fasta(const fs::path& source, std::ofstream& destination, const std::string& sag_id) {
  SequenceLineReader input(source);
  std::string line;
  while (input.getline(line)) {
    if (!line.empty() && line[0] == '>') destination << '>' << sag_id << '|' << line.substr(1) << '\n';
    else destination << line << '\n';
  }
}

static FastaStats fasta_stats(const fs::path& source) {
  SequenceLineReader input(source);
  FastaStats stats;
  std::uint64_t current_contig = 0;
  std::string line;
  while (input.getline(line)) {
    if (!line.empty() && line[0] == '>') {
      stats.max_contig = std::max(stats.max_contig, current_contig);
      current_contig = 0;
      continue;
    }
    consume_fasta_sequence_line(line, stats, current_contig);
  }
  stats.max_contig = std::max(stats.max_contig, current_contig);
  return stats;
}

// Length is derived from each assembly FASTA itself, never trusted from a
// manifest field.  Independent SAGs are scanned concurrently, but the fixed
// worker bound avoids turning a 220-thread run into a metadata/seek storm.
// No audit artifact is emitted until every reader has joined successfully;
// the first exception is then rethrown so the stage remains fail-closed.
static void populate_missing_assembly_bases(std::vector<Sag>& sags, std::size_t workers) {
  std::vector<std::size_t> pending;
  pending.reserve(sags.size());
  for (std::size_t i = 0; i < sags.size(); ++i) {
    if (!sags[i].assembly_bp_known) pending.push_back(i);
  }
  if (pending.empty()) return;
  workers = std::max<std::size_t>(1, std::min(workers, pending.size()));

  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::mutex error_mutex;
  std::exception_ptr first_error;
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (std::size_t worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&]() {
      while (!failed.load(std::memory_order_acquire)) {
        const std::size_t task = next.fetch_add(1, std::memory_order_relaxed);
        if (task >= pending.size()) return;
        try {
          Sag& sag = sags[pending[task]];
          const FastaStats stats = fasta_stats(sag.assembly);
          sag.assembly_bp = stats.total;
          sag.assembly_max_contig = stats.max_contig;
          sag.assembly_gc = stats.gc;
          sag.assembly_acgt = stats.acgt;
          sag.assembly_bp_known = true;
        } catch (...) {
          {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!first_error) first_error = std::current_exception();
          }
          failed.store(true, std::memory_order_release);
          return;
        }
      }
    });
  }
  for (auto& thread : pool) thread.join();
  if (first_error) std::rethrow_exception(first_error);
}


struct EmbeddedSketchTask {
  std::vector<fs::path> inputs;
  fs::path output;
  fs::path receipt;
  std::string stable_key;
};

static void run_embedded_sketches(const std::vector<EmbeddedSketchTask>& tasks,
                                  std::size_t workers, bool dry) {
  if (tasks.empty() || dry) return;
  workers = std::max<std::size_t>(1, std::min(workers, tasks.size()));
  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::mutex error_mutex;
  std::string first_error;
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (std::size_t worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&]() {
      while (!failed.load(std::memory_order_acquire)) {
        const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
        if (index >= tasks.size()) return;
        const auto& task = tasks[index];
        try {
          dna2bit_embedded::sketch_files(task.inputs, task.output, 17, 55296, 0);
          std::ofstream receipt(task.receipt);
          receipt << "PASS\nimplementation=embedded_dna2bit_source\nk=17\nbit_len=55296\nhash_type=0\n";
          for (const auto& input : task.inputs) receipt << "input=" << input.string() << "\n";
          receipt << "output=" << task.output.string() << "\n";
        } catch (const std::exception& error) {
          std::lock_guard<std::mutex> lock(error_mutex);
          if (first_error.empty()) first_error = task.stable_key + ": " + error.what();
          failed.store(true, std::memory_order_release);
          return;
        }
      }
    });
  }
  for (auto& worker : pool) worker.join();
  if (failed.load(std::memory_order_acquire)) {
    throw std::runtime_error("embedded DNA2bit sketch stage failed: " + first_error);
  }
}

static const char* input_kind_name(InputKind kind) {
  return kind == InputKind::paired_reads ? "paired_reads" : "contigs";
}

int main(int argc, char** argv) try {
  using SteadyClock = std::chrono::steady_clock;
  const auto total_started = SteadyClock::now();
  auto elapsed_seconds = [](SteadyClock::time_point begin, SteadyClock::time_point end) {
    return std::chrono::duration<double>(end - begin).count();
  };
  Config config = parse(argc, argv);
  require_file(config.dna_tax, "taxonomy table");
  if (config.dna_search_engine != "packed") {
    throw std::runtime_error("embedded package supports only --dna-search-engine packed");
  }
  require_file(config.dna_packed_db / "COMPLETE.json", "packed dna2bit database receipt");
  std::error_code packed_error;
  config.dna_packed_db = fs::weakly_canonical(config.dna_packed_db, packed_error);
  if (packed_error) throw std::runtime_error("cannot canonicalize packed dna2bit database");
  if (!config.stop_after_annotation) {
    config.subass = resolve_program(config.subass, argv[0], "C++ subassemble");
  }
  auto sags = read_manifest(config.manifest);
  fs::create_directories(config.out);

  std::size_t read_inputs = 0;
  std::size_t contig_inputs = 0;
  for (auto& sag : sags) {
    const auto directory = config.out / "01_assembly" / sag.id;
    const auto clean = directory / "clean";
    fs::create_directories(clean);
    sag.assembly = directory / (sag.id + ".fasta");

    if (sag.input_kind == InputKind::paired_reads) {
      ++read_inputs;
      sag.clean_r1 = clean / "R1.fastq.gz";
      sag.clean_r2 = clean / "R2.fastq.gz";
      const auto stage_pass = directory / "STAGE1.PASS";
      const auto spades_assembly = directory / "spades" / "scaffolds.fasta";
      if (!(config.resume && fs::exists(stage_pass))) {
        const std::string qc = q(config.fastp) + " -i " + q(sag.source1) + " -I " + q(sag.source2) +
                               " -o " + q(sag.clean_r1) + " -O " + q(sag.clean_r2) + " --thread " +
                               std::to_string(std::min(config.threads, 16)) + " --json " +
                               q(directory / "fastp.json") + " --html " + q(directory / "fastp.html");
        run(qc, directory / "FASTP.PASS", config.dry);
        const std::string spades = q(config.spades) + " --sc --careful -1 " + q(sag.clean_r1) + " -2 " +
                                   q(sag.clean_r2) + " -o " + q(directory / "spades") + " -t " +
                                   std::to_string(config.threads) + " -m " + std::to_string(config.memory_gb);
        run(spades, stage_pass, config.dry);
      }
      if (!config.dry && !fs::exists(sag.assembly)) {
        require_file(spades_assembly, "SPAdes scaffolds");
        fs::create_symlink(fs::absolute(spades_assembly), sag.assembly);
      }
    } else {
      ++contig_inputs;
      const auto stage_pass = directory / "STAGE1.PASS";
      if (!config.dry && !fs::exists(sag.assembly)) {
        FastaStats stats;
        sag.assembly_bp_known = materialize_contigs(sag.source1, sag.assembly, stats);
        if (sag.assembly_bp_known) {
          sag.assembly_bp = stats.total;
          sag.assembly_max_contig = stats.max_contig;
          sag.assembly_gc = stats.gc;
          sag.assembly_acgt = stats.acgt;
        }
      }
      if (!config.dry) {
        require_file(sag.assembly, "auto-routed contig FASTA");
        if (!fs::exists(stage_pass)) {
          std::ofstream(stage_pass) << "PASS\ninput_type=contigs\nsource=" << sag.source1.string()
                                    << "\nskipped=fastp,SPAdes\nstart_stage=assembly_length_gate\n";
        }
      }
      std::cerr << "+ [auto-route] " << sag.id << ": FASTA contigs; skip fastp and SPAdes\n";
    }
  }

  std::vector<Sag*> eligible;
  const auto excluded_path = config.out / "01_assembly" / "excluded_lt1000bp.tsv";
  // Bounded independent scans remove the old serial 8,785-file bottleneck.
  // A hard cap of 16 avoids turning a high --threads value into an I/O storm.
  const std::size_t length_scan_workers = static_cast<std::size_t>(std::min(config.threads, 16));
  if (!config.dry) {
    populate_missing_assembly_bases(sags, length_scan_workers);
    std::ofstream excluded(excluded_path);
    excluded << "sag_id\tassembly_fasta\ttotal_bp\treason\n";
    for (auto& sag : sags) {
      if (sag.assembly_bp < 1000) {
        excluded << sag.id << '\t' << sag.assembly.string() << '\t' << sag.assembly_bp
                 << "\ttotal_assembly_bp_lt_1000\n";
      } else {
        eligible.push_back(&sag);
      }
    }
    if (eligible.empty()) throw std::runtime_error("all SAGs were excluded by the <1000 bp hard cutoff");

    std::ofstream audit(config.out / "INPUT_AUDIT.tsv");
    audit << "sag_id\tdetected_input_type\tsource_1\tsource_2\tstart_stage\tskipped_stages\tassembly_fasta\tassembly_bp\teligible\n";
    for (const auto& sag : sags) {
      audit << sag.id << '\t' << input_kind_name(sag.input_kind) << '\t' << sag.source1.string() << '\t'
            << (sag.source2.empty() ? "NA" : sag.source2.string()) << '\t'
            << (sag.input_kind == InputKind::paired_reads ? "fastp" : "assembly_length_gate") << '\t'
            << (sag.input_kind == InputKind::paired_reads ? "none" : "fastp,SPAdes") << '\t'
            << sag.assembly.string() << '\t' << sag.assembly_bp << '\t'
            << (sag.assembly_bp >= 1000 ? "yes" : "no") << '\n';
    }

    // This is the one-pass, auditable FASTA-stat sidecar for Stage 3B.  It is
    // deliberately separate from INPUT_AUDIT to preserve that frozen schema.
    // CheckM2 contamination must be joined by the wrapper after a fresh run;
    // it cannot be inferred or safely filled during sequence scanning.
    std::ofstream stats(config.out / "STAGE3B_FASTA_STATS.tsv");
    stats << "sag_id\tassembly_fasta\ttotal_bp\tmax_contig\tgc_pct\tgc_bases\tacgt_bases\n";
    stats << std::fixed << std::setprecision(5);
    for (const auto& sag : sags) {
      const double gc_pct = sag.assembly_bp == 0
                                ? 0.0
                                : 100.0 * static_cast<double>(sag.assembly_gc) /
                                      static_cast<double>(sag.assembly_bp);
      stats << sag.id << '\t' << sag.assembly.string() << '\t' << sag.assembly_bp << '\t'
            << sag.assembly_max_contig << '\t' << gc_pct << '\t' << sag.assembly_gc << '\t'
            << sag.assembly_acgt << '\n';
    }
  } else {
    for (auto& sag : sags) eligible.push_back(&sag);
  }

  const double preflight_length_seconds = elapsed_seconds(total_started, SteadyClock::now());
  const auto sketch_started = SteadyClock::now();

  const auto bit_directory = config.out / "02_dna2bit" / "bits";
  fs::create_directories(bit_directory);
  // dna2bit sketches are independent per SAG.  Four threads is enough to
  // avoid spending most time in process/setup overhead on the Lake inputs,
  // while the product of worker count and -n stays within --threads.
  const auto [sketch_workers, sketch_threads] =
      bounded_parallel_shape(eligible.size(), config.threads, 4);
  std::vector<EmbeddedSketchTask> sketch_tasks;
  sketch_tasks.reserve(eligible.size());
  for (auto* sag : eligible) {
    const auto directory = bit_directory / sag->id;
    const auto pass = directory / "SKETCH.PASS";
    fs::create_directories(directory);
    if (config.resume && fs::exists(pass)) continue;

    if (sag->input_kind == InputKind::paired_reads) {
      ensure_symlink(sag->clean_r1, directory / "R1.fastq.gz");
      ensure_symlink(sag->clean_r2, directory / "R2.fastq.gz");
      const auto list = directory / "inputs.list";
      if (!config.dry) {
        std::ofstream output(list);
        output << "R1.fastq.gz\nR2.fastq.gz\n";
      }
      sketch_tasks.push_back({{sag->clean_r1, sag->clean_r2},
                              directory / "paired_reads.k.17.l.55296.bit", pass, sag->id});
    } else {
      ensure_symlink(sag->assembly, directory / "input.fasta");
      const auto list = directory / "inputs.list";
      if (!config.dry) {
        std::ofstream output(list);
        output << "input.fasta\n";
      }
      sketch_tasks.push_back({{sag->assembly},
                              directory / "input.fasta.k.17.l.55296.bit", pass, sag->id});
    }
  }
  run_embedded_sketches(sketch_tasks, static_cast<std::size_t>(sketch_workers), config.dry);
  if (!config.dry) {
    std::ofstream(config.out / "02_dna2bit" / "SKETCH.PASS")
        << "PASS\nimplementation=embedded_dna2bit_source\nmode=per-SAG_auto_routed\nk=17\nbit_len=55296\nhash_type=0\n";
  }

  const auto bit_list = config.out / "02_dna2bit" / "bits.list";
  if (!config.dry) {
    std::ofstream output(bit_list);
    const auto search_inputs = config.out / "02_dna2bit" / "search_inputs";
    fs::create_directories(search_inputs);
    std::size_t index = 0;
    for (auto* sag : eligible) {
      std::vector<fs::path> hits;
      for (const auto& entry : fs::directory_iterator(bit_directory / sag->id)) {
        if (entry.path().extension() == ".bit") hits.push_back(entry.path());
      }
      if (hits.size() != 1) {
        throw std::runtime_error("dna2bit must produce exactly one bit for " + sag->id + ", observed " +
                                 std::to_string(hits.size()));
      }
      sag->bit = hits[0];
      std::ostringstream name;
      name << 'b';
      name.width(6);
      name.fill('0');
      name << ++index << ".bit";
      const fs::path alias = search_inputs / name.str();
      ensure_symlink(sag->bit, alias);
      sag->bit_search_token = (fs::path("search_inputs") / name.str()).generic_string();
      output << sag->bit_search_token << '\n';
    }
  }

  const double sketch_seconds = elapsed_seconds(sketch_started, SteadyClock::now());
  const auto search_started = SteadyClock::now();

  const auto result = config.out / "02_dna2bit" / "search_result.csv";
  const auto search_pass = config.out / "02_dna2bit" / "SEARCH.PASS";
  if (!(config.resume && fs::exists(search_pass))) {
    if (config.dry) {
      std::cerr << "+ [embedded] packed DNA2bit search\n";
    } else {
      const auto query_directory = (config.out / "02_dna2bit" / "search_inputs").generic_string() + "/";
      dna2bit_packed::search_index_for_pipeline(
          config.dna_packed_db, query_directory, config.dna_tax, result, 0.01F,
          static_cast<unsigned>(config.threads));
      std::ofstream receipt(search_pass);
      receipt << "PASS\nimplementation=embedded_dna2bit_packed_search\nmin_ratio=0.01\n";
    }
  }

  if (config.dry) {
    std::cout << "DRY-RUN complete: auto-detected paired_reads=" << read_inputs << " contigs=" << contig_inputs
              << "; no scientific output written\n";
    return 0;
  }

  std::unordered_map<std::string, Sag*> by_bit;
  by_bit.reserve(eligible.size() * 3);
  for (auto* sag : eligible) {
    by_bit[sag->bit_search_token] = sag;
    by_bit[fs::path(sag->bit_search_token).filename().generic_string()] = sag;
    by_bit[fs::weakly_canonical(sag->bit).string()] = sag;
  }
  std::map<std::string, std::vector<Sag*>> groups;
  std::unordered_set<std::string> labeled;
  labeled.reserve(eligible.size());
  std::ifstream result_file(result);
  if (!result_file) throw std::runtime_error("cannot read dna2bit search result: " + result.string());
  std::string line;
  std::ofstream labels(config.out / "02_dna2bit" / "labels.tsv");
  labels << "sag_id\treference\ttaxonomy\tspecies_group\n";
  while (std::getline(result_file, line)) {
    const auto columns = split(line, ',');
    if (columns.size() < 3) continue;
    auto found = by_bit.find(fs::path(columns[0]).lexically_normal().generic_string());
    if (found == by_bit.end() && fs::path(columns[0]).is_relative()) {
      const auto relative_to_search = config.out / "02_dna2bit" / columns[0];
      if (fs::exists(relative_to_search)) found = by_bit.find(fs::weakly_canonical(relative_to_search).string());
    }
    if (found == by_bit.end() && fs::path(columns[0]).is_absolute()) {
      found = by_bit.find(fs::weakly_canonical(columns[0]).string());
    }
    if (found == by_bit.end()) continue;
    std::string taxonomy = columns[2];
    for (std::size_t i = 3; i < columns.size(); ++i) taxonomy += "," + columns[i];
    const auto key = species_key(taxonomy);
    groups[key].push_back(found->second);
    labeled.insert(found->second->id);
    labels << found->second->id << '\t' << columns[1] << '\t' << taxonomy << '\t' << key << '\n';
  }

  std::ofstream pending(config.out / "03B_unclassified_pending.tsv");
  pending << "sag_id\tassembly_fasta\treason\n";
  for (auto* sag : eligible) {
    if (!labeled.count(sag->id)) {
      pending << sag->id << '\t' << sag->assembly.string() << "\tdna2bit_rejected_or_no_hit\n";
    }
  }

  labels.close();
  pending.close();

  const double search_seconds = elapsed_seconds(search_started, SteadyClock::now());
  if (config.stop_after_annotation) {
    const auto scientific_finished = SteadyClock::now();
    const double total_seconds = elapsed_seconds(total_started, scientific_finished);
    {
      std::ofstream timing(config.out / "TIMING.tsv");
      timing << "phase\tseconds\n" << std::fixed << std::setprecision(6)
             << "preflight_length_gate\t" << preflight_length_seconds << '\n'
             << "sketch\t" << sketch_seconds << '\n'
             << "search\t" << search_seconds << '\n'
             << "total\t" << total_seconds << '\n';
    }
    {
      std::ofstream timing(config.out / "TIMING.json");
      timing << std::fixed << std::setprecision(6)
             << "{\"clock\":\"std::chrono::steady_clock\","
             << "\"workflow_endpoint\":\"annotation\","
             << "\"dna_search_engine\":\"" << config.dna_search_engine << "\","
             << "\"length_scan_workers\":" << length_scan_workers << ','
             << "\"preflight_length_gate_seconds\":" << preflight_length_seconds << ','
             << "\"sketch_seconds\":" << sketch_seconds << ','
             << "\"search_seconds\":" << search_seconds << ','
             << "\"total_seconds\":" << total_seconds << "}\n";
    }
    std::ofstream(config.out / "COMPLETE.json")
        << "{\"status\":\"PASS\",\"pipeline\":\"dna2bit-original-embedded-annotation-v1-auto-input\","
        << "\"workflow_endpoint\":\"annotation\",\"dna_search_engine\":\"" << config.dna_search_engine << "\","
        << "\"input_sags\":" << sags.size() << ",\"paired_read_inputs\":" << read_inputs
        << ",\"contig_inputs\":" << contig_inputs << ",\"eligible_sags\":" << eligible.size()
        << ",\"excluded_lt1000bp\":" << (sags.size() - eligible.size()) << ",\"labeled\":" << labeled.size()
        << ",\"pending\":" << (eligible.size() - labeled.size()) << "}\n";
    std::cout << "PASS annotation: input=" << sags.size() << " paired_reads=" << read_inputs
              << " contigs=" << contig_inputs << " eligible=" << eligible.size()
              << " excluded_lt1000bp=" << (sags.size() - eligible.size())
              << " labeled=" << labeled.size() << " pending=" << (eligible.size() - labeled.size()) << "\n";
    return 0;
  }
  const auto subass_started = SteadyClock::now();

  fs::create_directories(config.out / "03A_subassemble");
  std::ofstream group_manifest(config.out / "03A_subassemble" / "groups.tsv");
  group_manifest << "species_group\tsag_count\tinput_fasta\tbin_fasta\n";
  // Materialize each species input in map (lexicographic) order before any
  // assembly starts.  This preserves byte-for-byte FASTA member ordering even
  // though the independent C++ subassemble invocations below run concurrently.
  // Execute deterministic one-thread Flye assemble for all groups first, then
  // the independently verified two-thread repeat/contigger/polish phase.  The
  // split lets assembly use up to --threads groups concurrently while the
  // finish phase uses at most floor(--threads/2) x 2 threads.  Both phases
  // remain within the caller's global budget, and LPT order minimizes the
  // long tail in each phase.
  const int subass_budget = std::max(1, config.threads);
  const int finish_threads = subass_budget >= 2 ? 2 : 1;
  const std::size_t assemble_workers = std::min<std::size_t>(
      groups.size(), static_cast<std::size_t>(subass_budget));
  const std::size_t finish_workers = std::min<std::size_t>(
      groups.size(), static_cast<std::size_t>(std::max(1, subass_budget / finish_threads)));
  std::vector<CommandTask> assemble_tasks, finish_tasks;
  assemble_tasks.reserve(groups.size());
  finish_tasks.reserve(groups.size());
  for (auto& [key, members] : groups) {
    const auto group_directory = config.out / "03A_subassemble" / key;
    const auto merged_fasta = group_directory / "input_subassemblies.fasta";
    const auto output = group_directory / "run";
    fs::create_directories(group_directory);
    std::ofstream merged(merged_fasta);
    for (auto* sag : members) append_fasta(sag->assembly, merged, sag->id);
    merged.close();

    const fs::path flye = config.flye_root;
    const std::string command_prefix =
        q(config.subass) + " --reads " + q(merged_fasta) + " --out-dir " + q(output) + " --flye-modules " +
        q(flye / "bin/flye-modules") + " --minimap2 " + q(flye / "bin/flye-minimap2") + " --samtools " +
        q(flye / "bin/flye-samtools") + " --package-root " + q(flye / "lib/python3.9/site-packages/flye") +
        " --config " + q(flye / "lib/python3.9/site-packages/flye/config/bin_cfg/asm_subasm.cfg");
    if (!(config.resume && fs::exists(output / "CPP_FULL_PIPELINE_PASS"))) {
      const auto bytes = fs::file_size(merged_fasta);
      if (!(config.resume && fs::exists(output / "00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE"))) {
        assemble_tasks.push_back({command_prefix +
                                      " --threads 1 --assemble-threads 1 --phase assemble"
                                      " --no-overlap-policy passthrough",
                                  output / "ASSEMBLY.PASS", bytes, key});
      }
      finish_tasks.push_back({command_prefix + " --threads " + std::to_string(finish_threads) +
                                  " --assemble-threads 1 --phase finish --no-overlap-policy passthrough",
                              output / "PIPELINE.PASS", bytes, key});
    }
    group_manifest << key << '\t' << members.size() << '\t' << merged_fasta.string() << '\t' << (output / "assembly.fasta").string()
                   << '\n';
  }
  group_manifest.close();
  auto lpt = [](const CommandTask& a, const CommandTask& b) {
    if (a.work_bytes != b.work_bytes) return a.work_bytes > b.work_bytes;
    return a.stable_key < b.stable_key;
  };
  std::sort(assemble_tasks.begin(), assemble_tasks.end(), lpt);
  std::sort(finish_tasks.begin(), finish_tasks.end(), lpt);
  run_bounded_parallel(assemble_tasks, assemble_workers, config.dry);
  const auto subass_assembly_finished = SteadyClock::now();
  run_bounded_parallel(finish_tasks, finish_workers, config.dry);

  const auto scientific_finished = SteadyClock::now();
  const double subass_assemble_seconds = elapsed_seconds(subass_started, subass_assembly_finished);
  const double subass_finish_seconds = elapsed_seconds(subass_assembly_finished, scientific_finished);
  const double subass_seconds = elapsed_seconds(subass_started, scientific_finished);
  const double total_seconds = elapsed_seconds(total_started, scientific_finished);
  {
    std::ofstream timing(config.out / "TIMING.tsv");
    timing << "phase\tseconds\n" << std::fixed << std::setprecision(6)
           << "preflight_length_gate\t" << preflight_length_seconds << '\n'
           << "sketch\t" << sketch_seconds << '\n'
           << "search\t" << search_seconds << '\n'
           << "stage3a_assemble_1t\t" << subass_assemble_seconds << '\n'
           << "stage3a_finish_2t\t" << subass_finish_seconds << '\n'
           << "stage3a_subassemble\t" << subass_seconds << '\n'
           << "total\t" << total_seconds << '\n';
  }
  {
    std::ofstream timing(config.out / "TIMING.json");
    timing << std::fixed << std::setprecision(6)
           << "{\"clock\":\"std::chrono::steady_clock\","
           << "\"dna_search_engine\":\"" << config.dna_search_engine << "\","
           << "\"length_scan_workers\":" << length_scan_workers << ','
           << "\"preflight_length_gate_seconds\":" << preflight_length_seconds << ','
           << "\"sketch_seconds\":" << sketch_seconds << ','
           << "\"search_seconds\":" << search_seconds << ','
           << "\"stage3a_assemble_1t_seconds\":" << subass_assemble_seconds << ','
           << "\"stage3a_finish_2t_seconds\":" << subass_finish_seconds << ','
           << "\"stage3a_subassemble_seconds\":" << subass_seconds << ','
           << "\"total_seconds\":" << total_seconds << "}\n";
  }

  std::ofstream(config.out / "COMPLETE.json")
      << "{\"status\":\"PASS\",\"pipeline\":\"dna2bit-original-embedded-stage1-3A-v1-auto-input\","
      << "\"dna_search_engine\":\"" << config.dna_search_engine << "\","
      << "\"input_sags\":" << sags.size() << ",\"paired_read_inputs\":" << read_inputs
      << ",\"contig_inputs\":" << contig_inputs << ",\"eligible_sags\":" << eligible.size()
      << ",\"excluded_lt1000bp\":" << (sags.size() - eligible.size()) << ",\"labeled\":" << labeled.size()
      << ",\"groups\":" << groups.size() << "}\n";
  std::cout << "PASS: input=" << sags.size() << " paired_reads=" << read_inputs << " contigs=" << contig_inputs
            << " eligible=" << eligible.size() << " excluded_lt1000bp=" << (sags.size() - eligible.size())
            << " labeled=" << labeled.size() << " groups=" << groups.size() << "\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << "fatal: " << error.what() << "\n";
  return 2;
}
