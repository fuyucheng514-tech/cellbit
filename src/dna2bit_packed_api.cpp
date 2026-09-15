// Reuse the audited teacher-compatible packed search implementation as a
// library translation unit.  The standalone CLI remains available for
// independent index/search audits, while the main pipeline links this file
// and invokes search_index_for_pipeline directly.
#define DNA2BIT_PACKED_LIBRARY 1
#include "dna2bit_packed.cpp"

namespace dna2bit_packed {

void search_index_for_pipeline(const std::filesystem::path& index,
                               const std::string& query_directory_with_slash,
                               const std::filesystem::path& taxonomy,
                               const std::filesystem::path& output,
                               const float min_ratio,
                               const unsigned threads) {
  SearchOptions options;
  options.index = index;
  options.queries = query_directory_with_slash;
  options.taxonomy = taxonomy;
  options.output = output;
  options.min_ratio = min_ratio;
  options.threads = threads;
  search_index(options);
}

}  // namespace dna2bit_packed
