#include "dna2bit_embedded_api.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>
#ifdef DNA2BIT_HAVE_LIBDEFLATE
#include <libdeflate.h>
#endif

#include "dna2bit_embedded_vendor/MurmurHash3.h"
#include "dna2bit_embedded_vendor/kseq.h"
#include "dna2bit_embedded_vendor/rollinghash.h"
#include "dna2bit_embedded_vendor/wyhash.h"

namespace {

// A bounded, optional gzip fast path. Keep the same kseq parser and use zlib
// for plain input, multi-member gzip, large streams, or unsupported containers.
// CRC/size checks are performed by libdeflate before any sequence is exposed.
class SequenceInput {
 public:
  explicit SequenceInput(const std::filesystem::path& path) {
#ifdef DNA2BIT_HAVE_LIBDEFLATE
    if (load_gzip(path)) return;
#endif
    stream_ = gzopen(path.string().c_str(), "rb");
    if (!stream_) throw std::runtime_error("cannot open sequence for embedded DNA2bit: " + path.string());
  }
  ~SequenceInput() { if (stream_) gzclose(stream_); }
  SequenceInput(const SequenceInput&) = delete;
  SequenceInput& operator=(const SequenceInput&) = delete;
  int read(void* output, unsigned count) {
    if (stream_) return gzread(stream_, output, count);
    const std::size_t available = std::min<std::size_t>(count, decoded_.size() - position_);
    if (available) std::memcpy(output, decoded_.data() + position_, available);
    position_ += available;
    return static_cast<int>(available);
  }
 private:
  gzFile stream_ = nullptr;
  std::vector<unsigned char> decoded_;
  std::size_t position_ = 0;
#ifdef DNA2BIT_HAVE_LIBDEFLATE
  bool load_gzip(const std::filesystem::path& path) {
    constexpr std::uintmax_t limit = 64U << 20;
    if (path.extension() != ".gz") return false;
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error) || error) return false;
    const auto size = std::filesystem::file_size(path, error);
    if (error || size < 18 || size > limit) return false;
    std::ifstream input(path, std::ios::binary);
    if (!input) return false;
    std::vector<unsigned char> compressed(static_cast<std::size_t>(size));
    if (!input.read(reinterpret_cast<char*>(compressed.data()), compressed.size())) return false;
    if (compressed[0] != 31 || compressed[1] != 139 || compressed[2] != 8) return false;
    const auto* trailer = compressed.data() + compressed.size() - 4;
    const std::uint32_t expected = static_cast<std::uint32_t>(trailer[0]) |
        (static_cast<std::uint32_t>(trailer[1]) << 8) |
        (static_cast<std::uint32_t>(trailer[2]) << 16) |
        (static_cast<std::uint32_t>(trailer[3]) << 24);
    if (expected > limit) return false;
    using Decoder = std::unique_ptr<libdeflate_decompressor, decltype(&libdeflate_free_decompressor)>;
    thread_local Decoder decoder(libdeflate_alloc_decompressor(), libdeflate_free_decompressor);
    if (!decoder) return false;
    std::vector<unsigned char> decoded(std::max<std::size_t>(1, expected));
    std::size_t consumed = 0;
    const auto result = libdeflate_gzip_decompress_ex(decoder.get(), compressed.data(),
        compressed.size(), decoded.data(), expected, &consumed, nullptr);
    // In particular, do not silently discard additional concatenated members.
    if (result != LIBDEFLATE_SUCCESS || consumed != compressed.size()) return false;
    decoded.resize(expected);
    decoded_ = std::move(decoded);
    return true;
  }
#endif
};

int read_sequence(SequenceInput* input, void* output, unsigned count) {
  return input->read(output, count);
}

}  // namespace

KSEQ_INIT(SequenceInput*, read_sequence)

namespace {

template <class Count>
void read2dis_wy(const char* read, const std::size_t sliding_len,
                 std::vector<Count>& dis) {
  const std::size_t length = std::strlen(read);
  if (length <= sliding_len) return;
  const std::size_t vec_size = dis.size();
  for (std::size_t i = 0; i < length - sliding_len; ++i) {
    const std::uint64_t hash = wyhash(read + i, sliding_len, 0, _wyp);
    dis[(hash & 0x7FFFFFFFFFFFFFFFULL) % vec_size] +=
        2 * static_cast<long>(hash >> 63) - 1;
  }
}

void read2dis_mu(const char* read, const std::size_t sliding_len,
                 std::vector<long>& dis) {
  const std::size_t length = std::strlen(read);
  if (length <= sliding_len) return;
  const std::size_t vec_size = dis.size();
  char data[16];
  std::uint64_t hash = 0;
  for (std::size_t i = 0; i < length - sliding_len; ++i) {
    (void)i;  // Preserve the archived implementation's fixed-window hash call.
    MurmurHash3_x64_128(read, static_cast<int>(sliding_len), 0, data);
    hash = *reinterpret_cast<std::uint64_t*>(data);
    dis[(hash & 0x7FFFFFFFFFFFFFFFULL) % vec_size] +=
        2 * static_cast<long>(hash >> 63) - 1;
  }
}

template <class Count>
void reads2dis_wy(const char* code, const char* read,
                  const std::size_t sliding_len, std::vector<Count>& dis) {
  const std::size_t length = std::strlen(read);
  if (length <= sliding_len) return;
  std::vector<char> complement(length + 1, '\0');
  for (std::size_t i = 0; i < length - 1; ++i) {
    complement[i] = code[static_cast<unsigned char>(read[length - 2 - i])];
  }
  // The archived code writes a newline sentinel, which is not a nucleotide
  // and is intentionally retained for teacher-compatible k-mer traversal.
  complement[length - 1] = '\n';
  read2dis_wy(read, sliding_len, dis);
  read2dis_wy(complement.data(), sliding_len, dis);
}

void reads2dis_mu(const char* code, const char* read,
                  const std::size_t sliding_len, std::vector<long>& dis) {
  const std::size_t length = std::strlen(read);
  if (length <= sliding_len) return;
  std::vector<char> complement(length + 1, '\0');
  for (std::size_t i = 0; i < length - 1; ++i) {
    complement[i] = code[static_cast<unsigned char>(read[length - 2 - i])];
  }
  complement[length - 1] = '\n';
  read2dis_mu(read, sliding_len, dis);
  read2dis_mu(complement.data(), sliding_len, dis);
}

void dis2bit(std::vector<long>& bit) {
  for (long& value : bit) value = (value >> 63) + 1;
}

void handle_fasta_wy(SequenceInput* input, const char* code, const std::size_t kmer,
                     std::vector<long>& bit) {
  kseq_t* sequence = kseq_init(input);
  // Smaller counters keep the random-access accumulator cache-friendly. A
  // conservative window budget proves that no int32 counter can overflow.
  // Flush exact signed counts into the legacy wide accumulator before the
  // budget is exhausted; unusually long records use the original wide path.
  std::vector<std::int32_t> narrow(bit.size(), 0);
#ifdef DNA2BIT_TEST_COUNTER_BUDGET
  constexpr std::size_t max_updates = DNA2BIT_TEST_COUNTER_BUDGET;
  static_assert(max_updates > 0 && max_updates <= std::numeric_limits<std::int32_t>::max(),
                "invalid test counter budget");
#else
  constexpr std::size_t max_updates = std::numeric_limits<std::int32_t>::max();
#endif
  std::size_t available = max_updates;
  auto flush = [&] {
    for (std::size_t i = 0; i < bit.size(); ++i) {
      bit[i] += narrow[i];
      narrow[i] = 0;
    }
    available = max_updates;
  };
  while (kseq_read(sequence) >= 0) {
    const std::size_t length = sequence->seq.l;
    if (length <= kmer) continue;
    if (length - kmer > max_updates / 2) {
      flush();
      reads2dis_wy(code, sequence->seq.s, kmer, bit);
    } else {
      const std::size_t updates = 2 * (length - kmer);
      if (updates > available) flush();
      reads2dis_wy(code, sequence->seq.s, kmer, narrow);
      available -= updates;
    }
  }
  flush();
  kseq_destroy(sequence);
}

void handle_fasta_ro(SequenceInput* input, const std::size_t kmer,
                     std::vector<long>& bit) {
  const std::size_t bit_len = bit.size();
  kseq_t* sequence = kseq_init(input);
  while (kseq_read(sequence) >= 0) {
    const std::string read(sequence->seq.s);
    if (read.size() <= kmer) continue;
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    for (int i = static_cast<int>(kmer) - 1; i >= 0; --i) {
      forward = r33(forward) ^ Tab[static_cast<unsigned char>(read[kmer - 1 - i])];
      reverse = r33(reverse) ^ Tab[static_cast<unsigned char>(read[i]) & doff];
    }
    bit[(forward & 0x7FFFFFFFFFFFFFFFULL) % bit_len] +=
        2 * static_cast<long>(forward >> 63) - 1;
    bit[(reverse & 0x7FFFFFFFFFFFFFFFULL) % bit_len] +=
        2 * static_cast<long>(reverse >> 63) - 1;
    for (std::size_t i = 0; i < read.length() - kmer; ++i) {
      forward = r33(forward) ^ Tab[static_cast<unsigned char>(read[i + kmer])] ^
                opchar(static_cast<unsigned char>(read[i]), kmer);
      reverse = r3263(reverse ^
                      opchar(static_cast<unsigned char>(read[i + kmer]) & doff, kmer) ^
                      Tab[static_cast<unsigned char>(read[i]) & doff]);
      bit[(forward & 0x7FFFFFFFFFFFFFFFULL) % bit_len] +=
          2 * static_cast<long>(forward >> 63) - 1;
      bit[(reverse & 0x7FFFFFFFFFFFFFFFULL) % bit_len] +=
          2 * static_cast<long>(reverse >> 63) - 1;
    }
  }
  kseq_destroy(sequence);
}

void handle_fasta_mu(SequenceInput* input, const char* code, const std::size_t kmer,
                     std::vector<long>& bit) {
  kseq_t* sequence = kseq_init(input);
  while (kseq_read(sequence) >= 0) {
    reads2dis_mu(code, sequence->seq.s, kmer, bit);
  }
  kseq_destroy(sequence);
}

void write_bits(std::FILE* output, const std::vector<long>& bit) {
  const std::size_t count = bit.size() / 64;
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t offset = i * 64;
    std::uint64_t value = 0;
    for (std::size_t j = 0; j < 64; ++j) {
      value = (value << 1) | static_cast<std::uint64_t>(bit[offset + j]);
    }
    if (std::fwrite(&value, sizeof(value), 1, output) != 1) {
      throw std::runtime_error("failed writing embedded DNA2bit sketch");
    }
  }
  if (bit.size() % 64 != 0) {
    std::uint64_t value = 0;
    for (std::size_t i = count * 64; i < bit.size(); ++i) {
      value = (value << 1) | static_cast<std::uint64_t>(bit[i]);
    }
    if (std::fwrite(&value, sizeof(value), 1, output) != 1) {
      throw std::runtime_error("failed writing embedded DNA2bit tail");
    }
  }
}

}  // namespace

namespace dna2bit_embedded {

void sketch_file(const std::filesystem::path& input,
                 const std::filesystem::path& output,
                 const std::size_t kmer_len,
                 const std::size_t bit_len,
                 const std::size_t hash_type) {
  sketch_files({input}, output, kmer_len, bit_len, hash_type);
}

void sketch_files(const std::vector<std::filesystem::path>& inputs,
                  const std::filesystem::path& output,
                  const std::size_t kmer_len,
                  const std::size_t bit_len,
                  const std::size_t hash_type) {
  if (kmer_len == 0 || bit_len == 0) {
    throw std::runtime_error("embedded DNA2bit requires positive kmer and bit lengths");
  }
  if (inputs.empty()) throw std::runtime_error("embedded DNA2bit received no inputs");
  std::filesystem::create_directories(output.parent_path());
  std::vector<long> bit(bit_len, 0);
  char code[256];
  std::fill(std::begin(code), std::end(code), 'N');
  code[static_cast<unsigned char>('A')] = 'T';
  code[static_cast<unsigned char>('T')] = 'A';
  code[static_cast<unsigned char>('G')] = 'C';
  code[static_cast<unsigned char>('C')] = 'G';

  std::FILE* output_file = std::fopen(output.string().c_str(), "wb");
  if (!output_file) {
    throw std::runtime_error("cannot create embedded DNA2bit sketch: " + output.string());
  }
  try {
    for (const auto& input : inputs) {
      SequenceInput file(input);
      switch (hash_type) {
        case 0: handle_fasta_wy(&file, code, kmer_len, bit); break;
        case 1: handle_fasta_ro(&file, kmer_len, bit); break;
        case 2: handle_fasta_mu(&file, code, kmer_len, bit); break;
        default: throw std::runtime_error("unsupported embedded DNA2bit hash type");
      }
    }
    dis2bit(bit);
    write_bits(output_file, bit);
    if (std::fclose(output_file) != 0) {
      output_file = nullptr;
      throw std::runtime_error("failed closing embedded DNA2bit sketch: " + output.string());
    }
  } catch (...) {
    std::fclose(output_file);
    throw;
  }
}

}  // namespace dna2bit_embedded
