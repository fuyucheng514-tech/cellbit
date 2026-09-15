#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

// In-process implementation of the original teacher DNA2bit sketch step.
// The bit layout, k=17/bit_len=55296/hash_type=0 defaults and reverse
// complement handling are kept byte-compatible with the archived source.
namespace dna2bit_embedded {

void sketch_file(const std::filesystem::path& input,
                 const std::filesystem::path& output,
                 std::size_t kmer_len = 17,
                 std::size_t bit_len = 55296,
                 std::size_t hash_type = 0);

void sketch_files(const std::vector<std::filesystem::path>& inputs,
                  const std::filesystem::path& output,
                  std::size_t kmer_len = 17,
                  std::size_t bit_len = 55296,
                  std::size_t hash_type = 0);

}  // namespace dna2bit_embedded
