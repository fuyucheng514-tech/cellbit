#pragma once

#include <filesystem>
#include <string>

namespace subass {

enum class NoOverlapPolicy { strict, passthrough };

NoOverlapPolicy parse_no_overlap_policy(const std::string& value);
const char* no_overlap_policy_name(NoOverlapPolicy policy);

struct NoOverlapContext {
    NoOverlapPolicy policy = NoOverlapPolicy::strict;
    int assembly_return_code = 0;
    std::filesystem::path input_subassemblies;
    std::filesystem::path draft_assembly;
    std::filesystem::path assembly_log;
    std::filesystem::path output_directory;
    std::string assembly_command;
};

// Return true only when an empty-draft/no-overlap run was safely promoted to
// byte-exact input passthrough.  Return false for a normal non-empty draft so
// the unchanged Flye repeat/contigger/polishing path can continue.  Every
// other abnormal state throws and therefore fails closed.
bool handle_no_overlap(const NoOverlapContext& context);

}  // namespace subass
