#include "polish_glue.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

static std::unordered_map<std::string, int> load_coverage(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("open coverage stats: " + path.string());
    std::unordered_map<std::string, int> coverage;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        const auto first = line.find('\t');
        const auto second = first == std::string::npos ? first : line.find('\t', first + 1);
        if (first == std::string::npos || second == std::string::npos)
            throw std::runtime_error("invalid coverage stats line");
        coverage[line.substr(0, first)] = std::stoi(line.substr(second + 1));
    }
    return coverage;
}

int main(int argc, char** argv) try {
    if (argc < 2) throw std::runtime_error("mode required: compose, compose-files, filter, or gfa");
    const std::string mode = argv[1];
    if (mode == "compose") {
        if (argc != 8)
            throw std::runtime_error("compose: consensus coverage.tsv out.fa out.stats out.bed.gz expected.tsv");
        auto coverage = load_coverage(argv[3]);
        const auto result = subass::compose_consensus(argv[2], argv[4], coverage, argv[5], argv[6]);
        std::ofstream summary(argv[7]);
        if (!summary) throw std::runtime_error("create compose summary");
        size_t bases = 0;
        for (const auto& item : result.lengths) bases += item.second;
        summary << "sequences\t" << result.sequences.size() << '\n'
                << "lengths\t" << result.lengths.size() << '\n'
                << "bases\t" << bases << '\n';
    } else if (mode == "compose-files") {
        if (argc != 7)
            throw std::runtime_error("compose-files: consensus coverage.tsv out.fa out.stats out.bed.gz");
        auto coverage = load_coverage(argv[3]);
        subass::compose_consensus_to_files(argv[2], argv[4], coverage, argv[5], argv[6]);
    } else if (mode == "filter") {
        if (argc != 6) throw std::runtime_error("filter: stats.in contigs.in stats.out contigs.out");
        subass::filter_subassembly_coverage(argv[2], argv[3], argv[4], argv[5]);
    } else if (mode == "gfa") {
        if (argc != 7) throw std::runtime_error("gfa: edges.fa graph.gfa polished.fa alignments.bam out.gfa");
        subass::generate_polished_gfa(argv[2], argv[3], argv[4], argv[5], argv[6]);
    } else {
        throw std::runtime_error("unknown mode: " + mode);
    }
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
