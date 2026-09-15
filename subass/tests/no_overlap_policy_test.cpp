#include "no_overlap_policy.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

static void write(const fs::path& path, const std::string& content) {
    fs::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    output << content;
    output.close();
    if (!output) throw std::runtime_error("failed writing fixture");
}

static std::string read(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    std::ostringstream output;
    output << input.rdbuf();
    return output.str();
}

static void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
static void require_throws(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

static subass::NoOverlapContext context(const fs::path& root,
                                        subass::NoOverlapPolicy policy,
                                        int return_code = 0) {
    return {policy,
            return_code,
            root / "input.fasta",
            root / "00-assembly/draft_assembly.fasta",
            root / "00-assembly/no_overlap_evidence.log",
            root,
            "fake-flye-modules assemble --fixture"};
}

int main() try {
    const auto nonce = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto temporary = fs::temp_directory_path() / ("cpp-subass-no-overlap-test-" + std::to_string(nonce));
    fs::create_directories(temporary);

    const std::string single = ">SAG01__contig_1 description\nACGTNN\n>SAG01__contig_2\nGGCC\n";
    const std::string multiple = ">SAG01__a\nAAAA\n>SAG02__b\nCCCC\n>SAG03__c\nNNNN\n";

    // Default/strict mode retains the original behavior: an empty draft is a
    // hard error even when the log explicitly says there were no overlaps.
    const auto strict = temporary / "strict";
    write(strict / "input.fasta", single);
    write(strict / "00-assembly/draft_assembly.fasta", "");
    write(strict / "00-assembly/no_overlap_evidence.log", "No overlaps found!\n");
    require_throws([&] { subass::handle_no_overlap(context(strict, subass::NoOverlapPolicy::strict)); },
                   "strict mode accepted an empty draft");
    require(!fs::exists(strict / "assembly.fasta"), "strict mode published an assembly");

    // A single-SAG/no-overlap group is promoted only under the explicit policy.
    const auto one = temporary / "one";
    write(one / "input.fasta", single);
    write(one / "00-assembly/draft_assembly.fasta", "");
    write(one / "00-assembly/no_overlap_evidence.log", "No overlaps found - unable to estimate parameters\n");
    require(subass::handle_no_overlap(context(one, subass::NoOverlapPolicy::passthrough)),
            "single-SAG passthrough did not trigger");
    require(read(one / "assembly.fasta") == single, "single-SAG passthrough was not byte-exact");
    const auto one_receipt = read(one / "NO_OVERLAP_PASSTHROUGH.PASS.json");
    require(one_receipt.find("\"status\": \"PASS\"") != std::string::npos,
            "single-SAG PASS receipt missing");
    require(one_receipt.find("\"records\": 2") != std::string::npos,
            "single-SAG receipt has wrong record count");
    require(one_receipt.find("\"byte_exact\": true") != std::string::npos,
            "single-SAG receipt lacks byte closure");

    // Multiple source SAGs in one species group follow the same auditable rule.
    const auto many = temporary / "many";
    write(many / "input.fasta", multiple);
    write(many / "00-assembly/draft_assembly.fasta", "");
    write(many / "00-assembly/no_overlap_evidence.log", "Assembled 0 disjointigs\n");
    require(subass::handle_no_overlap(context(many, subass::NoOverlapPolicy::passthrough)),
            "multi-SAG passthrough did not trigger");
    require(read(many / "assembly.fasta") == multiple, "multi-SAG passthrough was not byte-exact");
    require(read(many / "NO_OVERLAP_PASSTHROUGH.PASS.json").find("\"records\": 3") !=
                std::string::npos,
            "multi-SAG receipt has wrong record count");

    // A normal non-empty draft must not enter the fallback.  Returning false
    // leaves the original repeat/contigger/polishing code path unchanged.
    const auto normal = temporary / "normal";
    write(normal / "input.fasta", single);
    write(normal / "00-assembly/draft_assembly.fasta", ">disjointig_1\nACGTACGT\n");
    require(!subass::handle_no_overlap(context(normal, subass::NoOverlapPolicy::passthrough)),
            "normal draft was incorrectly intercepted");
    require(!fs::exists(normal / "assembly.fasta"), "normal path unexpectedly published fallback output");

    // Empty output is never enough on its own, and a non-zero child rc always
    // fails even if the three files otherwise resemble a no-overlap run.
    const auto no_evidence = temporary / "no-evidence";
    write(no_evidence / "input.fasta", single);
    write(no_evidence / "00-assembly/draft_assembly.fasta", "");
    write(no_evidence / "00-assembly/no_overlap_evidence.log", "assembly ended quietly\n");
    require_throws(
        [&] { subass::handle_no_overlap(context(no_evidence, subass::NoOverlapPolicy::passthrough)); },
        "empty draft without evidence was accepted");

    const auto failed = temporary / "failed";
    write(failed / "input.fasta", single);
    write(failed / "00-assembly/draft_assembly.fasta", "");
    write(failed / "00-assembly/no_overlap_evidence.log", "No overlaps found!\n");
    require_throws(
        [&] { subass::handle_no_overlap(context(failed, subass::NoOverlapPolicy::passthrough, 256)); },
        "non-zero assembly rc was accepted");

    fs::remove_all(temporary);
    std::cout << "PASS strict rejection, single/multi passthrough closure, normal-path preservation, fail-closed guards\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
}
