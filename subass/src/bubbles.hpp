#pragma once
#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace subass {
struct CoverageResult {
  std::unordered_map<std::string, int> mean_coverage;
  double mean_alignment_error = 0.0;
};

CoverageResult make_bubbles(const std::filesystem::path& bam,
                            const std::filesystem::path& reference,
                            const std::filesystem::path& output,
                            int threads);
}
