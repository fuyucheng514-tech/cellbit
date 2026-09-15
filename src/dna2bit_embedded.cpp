#include "dna2bit_embedded_api.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "dna2bit_embedded_vendor/MurmurHash3.h"
#include "dna2bit_embedded_vendor/kseq.h"
#include "dna2bit_embedded_vendor/rollinghash.h"
#include "dna2bit_embedded_vendor/wyhash.h"

KSEQ_INIT(gzFile, gzread)

namespace {

void read2dis_wy(const char* read, const std::size_t sliding_len,
                 std::vector<long>& dis) {
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

void reads2dis_wy(const char* code, const char* read,
                  const std::size_t sliding_len, std::vector<long>& dis) {
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

void handle_fasta_wy(gzFile input, const char* code, const std::size_t kmer,
                     std::vector<long>& bit) {
  kseq_t* sequence = kseq_init(input);
  while (kseq_read(sequence) >= 0) {
    reads2dis_wy(code, sequence->seq.s, kmer, bit);
  }
  kseq_destroy(sequence);
}

void handle_fasta_ro(gzFile input, const std::size_t kmer,
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

void handle_fasta_mu(gzFile input, const char* code, const std::size_t kmer,
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
      gzFile file = gzopen(input.string().c_str(), "rb");
      if (!file) {
        throw std::runtime_error("cannot open sequence for embedded DNA2bit: " + input.string());
      }
      try {
        switch (hash_type) {
          case 0: handle_fasta_wy(file, code, kmer_len, bit); break;
          case 1: handle_fasta_ro(file, kmer_len, bit); break;
          case 2: handle_fasta_mu(file, code, kmer_len, bit); break;
          default: throw std::runtime_error("unsupported embedded DNA2bit hash type");
        }
      } catch (...) {
        gzclose(file);
        throw;
      }
      gzclose(file);
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
