#include "no_overlap_policy.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace subass {
namespace {

class Sha256 {
    std::uint32_t state_[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                               0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::uint8_t block_[64] = {};
    std::size_t used_ = 0;
    std::uint64_t bits_ = 0;

    static std::uint32_t rotate_right(std::uint32_t value, int bits) {
        return (value >> bits) | (value << (32 - bits));
    }

    void transform() {
        static const std::uint32_t constants[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
            0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
            0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
            0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
            0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
            0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
            0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
            0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        std::uint32_t words[64];
        for (int i = 0; i < 16; ++i) {
            words[i] = (std::uint32_t(block_[4 * i]) << 24) |
                       (std::uint32_t(block_[4 * i + 1]) << 16) |
                       (std::uint32_t(block_[4 * i + 2]) << 8) | block_[4 * i + 3];
        }
        for (int i = 16; i < 64; ++i) {
            const auto s0 = rotate_right(words[i - 15], 7) ^ rotate_right(words[i - 15], 18) ^
                            (words[i - 15] >> 3);
            const auto s1 = rotate_right(words[i - 2], 17) ^ rotate_right(words[i - 2], 19) ^
                            (words[i - 2] >> 10);
            words[i] = words[i - 16] + s0 + words[i - 7] + s1;
        }
        auto a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        auto e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int i = 0; i < 64; ++i) {
            const auto s1 = rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
            const auto choose = (e & f) ^ ((~e) & g);
            const auto first = h + s1 + choose + constants[i] + words[i];
            const auto s0 = rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
            const auto majority = (a & b) ^ (a & c) ^ (b & c);
            const auto second = s0 + majority;
            h = g;
            g = f;
            f = e;
            e = d + first;
            d = c;
            c = b;
            b = a;
            a = first + second;
        }
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += h;
    }

public:
    void update(const void* data, std::size_t size) {
        auto* bytes = static_cast<const std::uint8_t*>(data);
        bits_ += std::uint64_t(size) * 8;
        while (size) {
            const auto take = std::min(size, 64 - used_);
            std::copy(bytes, bytes + take, block_ + used_);
            used_ += take;
            bytes += take;
            size -= take;
            if (used_ == 64) {
                transform();
                used_ = 0;
            }
        }
    }

    void update(const std::string& text) { update(text.data(), text.size()); }

    std::string finish() {
        const auto original_bits = bits_;
        const std::uint8_t one = 0x80;
        update(&one, 1);
        const std::uint8_t zero = 0;
        while (used_ != 56) update(&zero, 1);
        std::uint8_t length[8];
        for (int i = 0; i < 8; ++i) length[7 - i] = std::uint8_t(original_bits >> (8 * i));
        update(length, 8);
        std::ostringstream output;
        output << std::hex << std::setfill('0');
        for (const auto value : state_) output << std::setw(8) << value;
        return output.str();
    }
};

std::string sha256_text(const std::string& text) {
    Sha256 digest;
    digest.update(text);
    return digest.finish();
}

std::string sha256_file(const fs::path& path) {
    if (!fs::is_regular_file(path)) throw std::runtime_error("cannot hash missing file: " + path.string());
    const auto size_before = fs::file_size(path);
    const auto time_before = fs::last_write_time(path);
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot hash file: " + path.string());
    Sha256 digest;
    char buffer[1 << 16];
    while (input) {
        input.read(buffer, sizeof(buffer));
        const auto size = input.gcount();
        if (size > 0) digest.update(buffer, static_cast<std::size_t>(size));
    }
    if (!input.eof()) throw std::runtime_error("error while hashing file: " + path.string());
    if (size_before != fs::file_size(path) || time_before != fs::last_write_time(path)) {
        throw std::runtime_error("file changed while being hashed: " + path.string());
    }
    return digest.finish();
}

std::string json_escape(const std::string& text) {
    std::string output;
    for (const unsigned char value : text) {
        if (value == '\\') output += "\\\\";
        else if (value == '"') output += "\\\"";
        else if (value == '\n') output += "\\n";
        else if (value == '\r') output += "\\r";
        else if (value == '\t') output += "\\t";
        else if (value < 32) {
            std::ostringstream escaped;
            escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(value);
            output += escaped.str();
        } else {
            output += char(value);
        }
    }
    return output;
}

std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    return text;
}

struct FastaStats {
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
    std::uint64_t gc_bases = 0;
    std::uint64_t min_length = 0;
    std::uint64_t max_length = 0;
    std::string id_set_sha256;
    std::string record_set_sha256;
};

FastaStats fasta_stats(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read passthrough FASTA: " + path.string());
    FastaStats stats;
    std::set<std::string> seen_ids;
    std::vector<std::string> ids;
    std::vector<std::string> record_digests;
    Sha256 record_digest;
    std::uint64_t current_length = 0;
    bool have_record = false;

    auto finish_record = [&]() {
        if (!have_record) return;
        if (current_length == 0) throw std::runtime_error("empty FASTA record in " + path.string());
        record_digests.push_back(record_digest.finish());
        stats.min_length = stats.records == 0 ? current_length : std::min(stats.min_length, current_length);
        stats.max_length = std::max(stats.max_length, current_length);
        ++stats.records;
    };

    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line.front() == '>') {
            finish_record();
            const auto header = line.substr(1);
            const auto begin = header.find_first_not_of(" \t");
            if (begin == std::string::npos) throw std::runtime_error("empty FASTA header in " + path.string());
            const auto end = header.find_first_of(" \t", begin);
            const auto id = header.substr(begin, end == std::string::npos ? end : end - begin);
            if (!seen_ids.insert(id).second) throw std::runtime_error("duplicate FASTA primary ID: " + id);
            ids.push_back(id);
            record_digest = Sha256{};
            record_digest.update(id);
            const char separator = '\0';
            record_digest.update(&separator, 1);
            current_length = 0;
            have_record = true;
            continue;
        }
        if (!have_record) throw std::runtime_error("sequence before FASTA header in " + path.string());
        for (const unsigned char raw : line) {
            if (std::isspace(raw)) continue;
            if (!(std::isalpha(raw) || raw == '-' || raw == '.' || raw == '*')) {
                throw std::runtime_error("invalid FASTA sequence character in " + path.string());
            }
            const char base = static_cast<char>(std::toupper(raw));
            record_digest.update(&base, 1);
            ++current_length;
            ++stats.bases;
            if (base == 'G' || base == 'C') ++stats.gc_bases;
        }
    }
    if (!input.eof()) throw std::runtime_error("error while reading FASTA: " + path.string());
    finish_record();
    if (stats.records == 0) throw std::runtime_error("passthrough FASTA has no records: " + path.string());

    std::sort(ids.begin(), ids.end());
    Sha256 id_set;
    for (const auto& id : ids) {
        id_set.update(std::to_string(id.size()));
        id_set.update(":");
        id_set.update(id);
        id_set.update("\n");
    }
    stats.id_set_sha256 = id_set.finish();
    std::sort(record_digests.begin(), record_digests.end());
    Sha256 record_set;
    for (const auto& digest : record_digests) {
        record_set.update(digest);
        record_set.update("\n");
    }
    stats.record_set_sha256 = record_set.finish();
    return stats;
}

bool same_stats(const FastaStats& left, const FastaStats& right) {
    return left.records == right.records && left.bases == right.bases &&
           left.gc_bases == right.gc_bases && left.min_length == right.min_length &&
           left.max_length == right.max_length && left.id_set_sha256 == right.id_set_sha256 &&
           left.record_set_sha256 == right.record_set_sha256;
}

std::string read_text(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read assembly evidence log: " + path.string());
    std::ostringstream text;
    text << input.rdbuf();
    if (!input && !input.eof()) throw std::runtime_error("error reading assembly evidence log: " + path.string());
    return text.str();
}

std::string no_overlap_evidence(const fs::path& log) {
    const auto content = lower(read_text(log));
    if (content.find("no overlaps found") != std::string::npos) return "flye_no_overlaps_found";
    if (content.find("assembled 0 disjointig") != std::string::npos) return "flye_assembled_zero_disjointigs";
    if (content.find("0 disjointigs were assembled") != std::string::npos)
        return "flye_zero_disjointigs_were_assembled";
    return {};
}

void write_once(const fs::path& path, const std::string& content) {
    if (fs::exists(path)) throw std::runtime_error("refusing to overwrite receipt: " + path.string());
    const fs::path temporary = path.string() + ".tmp";
    if (fs::exists(temporary)) throw std::runtime_error("stale receipt temporary exists: " + temporary.string());
    {
        std::ofstream output(temporary, std::ios::binary);
        if (!output) throw std::runtime_error("cannot create receipt: " + temporary.string());
        output.write(content.data(), static_cast<std::streamsize>(content.size()));
        output.close();
        if (!output) throw std::runtime_error("cannot finish receipt: " + temporary.string());
    }
    fs::rename(temporary, path);
}

std::string canonical_string(const fs::path& path) { return fs::weakly_canonical(path).string(); }

}  // namespace

NoOverlapPolicy parse_no_overlap_policy(const std::string& value) {
    if (value == "strict") return NoOverlapPolicy::strict;
    if (value == "passthrough") return NoOverlapPolicy::passthrough;
    throw std::runtime_error("--no-overlap-policy must be strict or passthrough");
}

const char* no_overlap_policy_name(NoOverlapPolicy policy) {
    return policy == NoOverlapPolicy::strict ? "strict" : "passthrough";
}

bool handle_no_overlap(const NoOverlapContext& context) {
    // Never reinterpret a failed executable as a biological no-overlap case.
    if (context.assembly_return_code != 0) {
        throw std::runtime_error("assembly stage failed, rc=" +
                                 std::to_string(context.assembly_return_code));
    }
    if (!fs::is_regular_file(context.draft_assembly)) {
        throw std::runtime_error("assembly returned success without draft_assembly.fasta");
    }
    if (fs::file_size(context.draft_assembly) != 0) return false;

    if (context.policy == NoOverlapPolicy::strict) {
        throw std::runtime_error("no disjointigs were assembled (strict no-overlap policy)");
    }

    // The log is dedicated to this exact assembly invocation.  Empty draft
    // alone is insufficient: require one of Flye's explicit no-overlap/zero-
    // disjointig messages before promoting the input.
    if (!fs::is_regular_file(context.assembly_log)) {
        throw std::runtime_error("empty draft has no assembly evidence log");
    }
    const auto evidence = no_overlap_evidence(context.assembly_log);
    if (evidence.empty()) {
        throw std::runtime_error("empty draft lacks explicit Flye no-overlap evidence");
    }

    if (!fs::is_regular_file(context.input_subassemblies) ||
        fs::file_size(context.input_subassemblies) == 0) {
        throw std::runtime_error("cannot passthrough missing/empty input subassemblies");
    }
    fs::create_directories(context.output_directory);
    const auto input_stats = fasta_stats(context.input_subassemblies);
    const auto input_sha = sha256_file(context.input_subassemblies);
    const auto input_bytes = fs::file_size(context.input_subassemblies);

    const auto assembly = context.output_directory / "assembly.fasta";
    const fs::path temporary = assembly.string() + ".passthrough.tmp";
    if (fs::exists(assembly) || fs::exists(temporary)) {
        throw std::runtime_error("refusing to overwrite prior passthrough assembly");
    }
    if (!fs::copy_file(context.input_subassemblies, temporary, fs::copy_options::none)) {
        throw std::runtime_error("failed to copy no-overlap passthrough assembly");
    }
    const auto output_stats = fasta_stats(temporary);
    const auto output_sha = sha256_file(temporary);
    const auto output_bytes = fs::file_size(temporary);
    if (input_bytes != output_bytes || input_sha != output_sha || !same_stats(input_stats, output_stats)) {
        throw std::runtime_error("no-overlap passthrough input/output closure failure");
    }
    fs::rename(temporary, assembly);
    if (sha256_file(assembly) != input_sha || !same_stats(fasta_stats(assembly), input_stats)) {
        throw std::runtime_error("published passthrough assembly changed after atomic rename");
    }

    const auto log_sha = sha256_file(context.assembly_log);
    const auto draft_sha = sha256_file(context.draft_assembly);
    std::ostringstream receipt;
    receipt << "{\n"
            << "  \"schema\": \"cpp-subass-no-overlap-passthrough-v1\",\n"
            << "  \"status\": \"PASS\",\n"
            << "  \"mode\": \"byte_exact_input_passthrough\",\n"
            << "  \"policy\": \"passthrough\",\n"
            << "  \"trigger\": \"rc0_empty_draft_explicit_no_overlap_evidence\",\n"
            << "  \"evidence_kind\": \"" << evidence << "\",\n"
            << "  \"assembly_command\": \"" << json_escape(context.assembly_command) << "\",\n"
            << "  \"assembly_command_sha256\": \"" << sha256_text(context.assembly_command) << "\",\n"
            << "  \"assembly_log\": {\"path\": \"" << json_escape(canonical_string(context.assembly_log))
            << "\", \"bytes\": " << fs::file_size(context.assembly_log) << ", \"sha256\": \""
            << log_sha << "\"},\n"
            << "  \"empty_draft\": {\"path\": \"" << json_escape(canonical_string(context.draft_assembly))
            << "\", \"bytes\": 0, \"sha256\": \"" << draft_sha << "\"},\n"
            << "  \"input\": {\"path\": \"" << json_escape(canonical_string(context.input_subassemblies))
            << "\", \"bytes\": " << input_bytes << ", \"sha256\": \"" << input_sha << "\"},\n"
            << "  \"output\": {\"path\": \"" << json_escape(canonical_string(assembly))
            << "\", \"bytes\": " << output_bytes << ", \"sha256\": \"" << output_sha << "\"},\n"
            << "  \"fasta_stats\": {\"records\": " << input_stats.records
            << ", \"bases\": " << input_stats.bases << ", \"gc_bases\": " << input_stats.gc_bases
            << ", \"min_record_bases\": " << input_stats.min_length
            << ", \"max_record_bases\": " << input_stats.max_length
            << ", \"id_set_sha256\": \"" << input_stats.id_set_sha256
            << "\", \"record_set_sha256\": \"" << input_stats.record_set_sha256 << "\"},\n"
            << "  \"closure\": {\"byte_exact\": true, \"id_set_equal\": true, "
               "\"record_set_equal\": true, \"record_count_equal\": true, \"base_count_equal\": true}\n"
            << "}\n";
    const auto receipt_path = context.output_directory / "NO_OVERLAP_PASSTHROUGH.PASS.json";
    write_once(receipt_path, receipt.str());

    std::ostringstream marker;
    marker << "PASS\n"
           << "mode=byte_exact_input_passthrough\n"
           << "policy=passthrough\n"
           << "receipt=NO_OVERLAP_PASSTHROUGH.PASS.json\n"
           << "assembly_sha256=" << output_sha << "\n"
           << "records=" << input_stats.records << "\n"
           << "bases=" << input_stats.bases << "\n";
    write_once(context.output_directory / "CPP_FULL_PIPELINE_PASS", marker.str());
    return true;
}

}  // namespace subass
