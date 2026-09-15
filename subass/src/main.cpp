#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "bubbles.hpp"
#include "no_overlap_policy.hpp"
#include "polish_glue.hpp"

namespace fs = std::filesystem;

enum class Phase { all, assemble, finish };

struct Options {
    fs::path reads, out, flye_modules, minimap2, samtools, package_root, config;
    int threads = 1;
    int assemble_threads = 1;
    int min_overlap = 1000;
    Phase phase = Phase::all;
    subass::NoOverlapPolicy no_overlap_policy = subass::NoOverlapPolicy::strict;
};

static std::string quote(const fs::path& p) {
    std::string s = p.string(), out = "'";
    out.reserve(s.size() + 2);
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    return out + "'";
}

static void run(const std::string& command, const fs::path& receipt) {
    std::cerr << "+ " << command << "\n";
    const int rc = std::system(command.c_str());
    if (rc != 0) throw std::runtime_error("stage failed, rc=" + std::to_string(rc));
    std::ofstream(receipt) << "PASS\n" << command << "\n";
}

static void write_marker(const fs::path& path, const std::string& content) {
    if (fs::exists(path)) throw std::runtime_error("refusing to overwrite marker: " + path.string());
    std::ofstream output(path, std::ios::binary);
    if (!output) throw std::runtime_error("cannot write marker: " + path.string());
    output << content;
    output.close();
    if (!output) throw std::runtime_error("cannot finish marker: " + path.string());
}

static Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("missing value for " + a);
            return argv[i];
        };
        if (a == "--reads") o.reads = value();
        else if (a == "--out-dir") o.out = value();
        else if (a == "--flye-modules") o.flye_modules = value();
        else if (a == "--config") o.config = value();
        else if (a == "--minimap2") o.minimap2 = value();
        else if (a == "--samtools") o.samtools = value();
        else if (a == "--package-root") o.package_root = value();
        else if (a == "--threads") o.threads = std::stoi(value());
        else if (a == "--assemble-threads") o.assemble_threads = std::stoi(value());
        else if (a == "--phase") {
            const auto phase = value();
            if (phase == "all") o.phase = Phase::all;
            else if (phase == "assemble") o.phase = Phase::assemble;
            else if (phase == "finish") o.phase = Phase::finish;
            else throw std::runtime_error("--phase must be all, assemble, or finish");
        }
        else if (a == "--no-overlap-policy") o.no_overlap_policy = subass::parse_no_overlap_policy(value());
        else throw std::runtime_error("unknown argument: " + a);
    }
    if (o.reads.empty() || o.out.empty() || o.flye_modules.empty() || o.minimap2.empty() ||
        o.samtools.empty() || o.package_root.empty() || o.config.empty())
        throw std::runtime_error("required: --reads --out-dir --flye-modules --minimap2 --samtools --package-root --config");
    if (o.threads < 1) throw std::runtime_error("--threads must be positive");
    if (o.assemble_threads != 1) {
        throw std::runtime_error("--assemble-threads must be exactly 1 in the deterministic speed release");
    }
    return o;
}

int main(int argc, char** argv) try {
    const Options o = parse(argc, argv);
    fs::create_directories(o.out / "00-assembly");
    fs::create_directories(o.out / "20-repeat");
    fs::create_directories(o.out / "30-contigger");
    fs::create_directories(o.out / "40-polishing");
    const auto draft = o.out / "00-assembly/draft_assembly.fasta";
    const auto assembly_phase_marker = o.out / "00-assembly/CPP_ASSEMBLY_PHASE_COMPLETE";
    const auto finish_phase_marker = o.out / "40-polishing/CPP_FINISH_PHASE_COMPLETE";
    const auto full_pipeline_marker = o.out / "CPP_FULL_PIPELINE_PASS";

    if (o.phase == Phase::finish) {
        if (!fs::is_regular_file(assembly_phase_marker)) {
            throw std::runtime_error("finish phase requires CPP_ASSEMBLY_PHASE_COMPLETE");
        }
        // An explicit no-overlap assembly already published the byte-exact
        // passthrough and its full scientific receipt.  The finish invocation
        // only closes the scheduler phase; it must not touch that output.
        if (fs::is_regular_file(full_pipeline_marker)) {
            if (!fs::is_regular_file(o.out / "NO_OVERLAP_PASSTHROUGH.PASS.json") ||
                !fs::is_regular_file(o.out / "assembly.fasta")) {
                throw std::runtime_error("unexpected pre-existing full pipeline marker");
            }
            write_marker(finish_phase_marker, "PASS\nmode=no_overlap_passthrough_already_complete\n");
            return 0;
        }
        if (!fs::is_regular_file(draft) || fs::file_size(draft) == 0) {
            throw std::runtime_error("finish phase requires non-empty deterministic draft assembly");
        }
    }

    if (o.phase != Phase::finish) {
    const auto assembly_log = o.no_overlap_policy == subass::NoOverlapPolicy::strict
        ? o.out / "cpp-subass.log"
        : o.out / "00-assembly/no_overlap_evidence.log";
    if (o.no_overlap_policy == subass::NoOverlapPolicy::passthrough && fs::exists(assembly_log)) {
        throw std::runtime_error("refusing stale no-overlap evidence log: " + assembly_log.string());
    }
    const std::string assembly_command = quote(o.flye_modules) + " assemble --reads " + quote(o.reads) +
        " --out-asm " + quote(draft) + " --config " + quote(o.config) +
        " --log " + quote(assembly_log) + " --threads " +
        std::to_string(o.assemble_threads) + " --min-ovlp " + std::to_string(o.min_overlap);
    std::cerr << "+ " << assembly_command << "\n";
    const int assembly_rc = std::system(assembly_command.c_str());
    if (subass::handle_no_overlap({o.no_overlap_policy, assembly_rc, o.reads, draft,
                                  assembly_log, o.out, assembly_command})) {
        write_marker(assembly_phase_marker, "PASS\nassemble_threads=1\nmode=no_overlap_passthrough\n");
        if (o.phase == Phase::all) {
            write_marker(finish_phase_marker, "PASS\nmode=no_overlap_passthrough_already_complete\n");
        }
        std::cerr << "PASS: no-overlap input promoted byte-exactly to assembly.fasta; see "
                  << (o.out / "NO_OVERLAP_PASSTHROUGH.PASS.json") << "\n";
        return 0;
    }
    std::ofstream(o.out / "00-assembly/CPP_STAGE_PASS") << "PASS\n" << assembly_command << "\n";
    write_marker(assembly_phase_marker, "PASS\nassemble_threads=1\nmode=assembled_draft\n");
    if (o.phase == Phase::assemble) return 0;
    }

    run(quote(o.flye_modules) + " repeat --disjointigs " + quote(draft) +
        " --reads " + quote(o.reads) + " --out-dir " + quote(o.out / "20-repeat") +
        " --config " + quote(o.config) + " --log " + quote(o.out / "cpp-subass.log") +
        " --threads " + std::to_string(o.threads) + " --min-ovlp " + std::to_string(o.min_overlap),
        o.out / "20-repeat/CPP_STAGE_PASS");
    run(quote(o.flye_modules) + " contigger --graph-edges " +
        quote(o.out / "20-repeat/repeat_graph_edges.fasta") + " --reads " + quote(o.reads) +
        " --out-dir " + quote(o.out / "30-contigger") + " --config " + quote(o.config) +
        " --repeat-graph " + quote(o.out / "20-repeat/repeat_graph_dump") +
        " --graph-aln " + quote(o.out / "20-repeat/read_alignment_dump") +
        " --log " + quote(o.out / "cpp-subass.log") + " --threads " +
        std::to_string(o.threads) + " --min-ovlp " + std::to_string(o.min_overlap),
        o.out / "30-contigger/CPP_STAGE_PASS");
    const auto polish_dir=o.out/"40-polishing";
    auto align = [&](const fs::path& reference,const fs::path& query,const fs::path& bam,const fs::path& receipt){
      const std::string smem=fs::file_size(reference)>100ull*1024*1024?"4G":"1G";
      const std::string batch=fs::file_size(reference)>100ull*1024*1024?"5G":"1G";
      const int samtools_aux_threads=std::max(0,std::min(4,o.threads-1));
      std::string cmd="set -eo pipefail; "+quote(o.minimap2)+" "+quote(reference)+" "+quote(query)+
        " -x map-pb -t "+std::to_string(o.threads)+" -a -p 0.5 -N 10 --sam-hit-only -L -K "+batch+
        " -z 1000 -Q --secondary-seq -I 64G | "+quote(o.samtools)+" view -T "+quote(reference)+
        " -u - | "+quote(o.samtools)+" sort -T "+quote(polish_dir/"sort_tmp")+
        " -O bam -@ "+std::to_string(samtools_aux_threads)+" -l 1 -m "+smem+" -o "+quote(bam)+
        " && "+quote(o.samtools)+" index -c -@ "+std::to_string(samtools_aux_threads)+" "+quote(bam);
      run("bash -c \""+cmd+"\"",receipt);
    };
    const auto contigs=o.out/"30-contigger/contigs.fasta";
    const auto bam=polish_dir/"minimap_1.bam";
    align(contigs,o.reads,bam,polish_dir/"CPP_ALIGN_PASS");
    const auto bubbles=polish_dir/"bubbles_1.fasta";
    auto cov=subass::make_bubbles(bam,contigs,bubbles,o.threads);
    std::ofstream(polish_dir/"CPP_BUBBLES_PASS")<<"PASS\nmean_error="<<cov.mean_alignment_error<<"\n";
    if(fs::file_size(bubbles)==0)throw std::runtime_error("no polishing bubbles");
    const auto consensus=polish_dir/"consensus_1.fasta";
    run(quote(o.flye_modules)+" polisher --bubbles "+quote(bubbles)+" --subs-mat "+
        quote(o.package_root/"config/bin_cfg/pacbio_chm13_substitutions.mat")+" --hopo-mat "+
        quote(o.package_root/"config/bin_cfg/pacbio_chm13_homopolymers.mat")+" --out "+quote(consensus)+
        " --threads "+std::to_string(o.threads),polish_dir/"CPP_POLISHER_PASS");
    const auto polished_raw=polish_dir/"polished_1.raw.fasta";
    const auto raw_stats=polish_dir/"contigs_stats.raw.txt";
    subass::compose_consensus_to_files(consensus,polished_raw,cov.mean_coverage,raw_stats,polish_dir/"base_coverage.bed.gz");
    cov.mean_coverage.clear();
    cov.mean_coverage.rehash(0);
    const auto polished=polish_dir/"polished_1.fasta";
    const auto stats=polish_dir/"contigs_stats.txt";
    subass::filter_subassembly_coverage(raw_stats,polished_raw,stats,polished);
    std::ofstream(polish_dir/"CPP_COMPOSE_FILTER_PASS")<<"PASS\n";
    const auto edge_bam=polish_dir/"edges_aln.bam";
    align(polished,o.out/"30-contigger/graph_final.fasta",edge_bam,polish_dir/"CPP_EDGE_ALIGN_PASS");
    const auto polished_gfa=polish_dir/"polished_edges.gfa";
    subass::generate_polished_gfa(o.out/"30-contigger/graph_final.fasta",o.out/"30-contigger/graph_final.gfa",polished,edge_bam,polished_gfa);
    subass::finalize_no_scaffold(polished,o.out/"30-contigger/graph_final.gv",polished_gfa,o.out/"30-contigger/contigs_stats.txt",o.out);
    write_marker(finish_phase_marker, "PASS\ndownstream_threads=" + std::to_string(o.threads) + "\n");
    std::ofstream(full_pipeline_marker)
        << "PASS\nassemble_threads=" << o.assemble_threads
        << "\ndownstream_threads=" << o.threads << "\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << "\n";
    return 2;
}
