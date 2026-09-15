#pragma once
#include <filesystem>
#include <string>
#include <unordered_map>
namespace subass {
struct ComposeResult { std::unordered_map<std::string,int> lengths; std::unordered_map<std::string,std::string> sequences; };
ComposeResult compose_consensus(const std::filesystem::path&, const std::filesystem::path&,
 const std::unordered_map<std::string,int>&, const std::filesystem::path&, const std::filesystem::path&);
void compose_consensus_to_files(const std::filesystem::path&, const std::filesystem::path&,
 const std::unordered_map<std::string,int>&, const std::filesystem::path&, const std::filesystem::path&);
void filter_subassembly_coverage(const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&);
void generate_polished_gfa(const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&);
void finalize_no_scaffold(const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&,const std::filesystem::path&);
}
