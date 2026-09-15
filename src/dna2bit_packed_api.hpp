#pragma once

#include <filesystem>
#include <string>

namespace dna2bit_packed {

// Calls the frozen teacher-compatible packed search implementation in-process.
// No dna2bit executable or search helper process is launched by the pipeline.
void search_index_for_pipeline(const std::filesystem::path& index,
                               const std::string& query_directory_with_slash,
                               const std::filesystem::path& taxonomy,
                               const std::filesystem::path& output,
                               float min_ratio,
                               unsigned threads);

}  // namespace dna2bit_packed
