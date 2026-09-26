#include "dna2bit_embedded_api.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>

namespace {

void usage(const char* program) {
  std::cerr << "usage: " << program << " INPUT [INPUT ...] OUTPUT\n"
            << "       " << program
            << " --batch-list FILE --output-dir DIR --threads N\n";
}

// Each reference remains an independent sketch. Only the scheduling changes:
// persistent native workers replace one executable launch per reference.
int batch_main(int argc, char** argv) {
  std::filesystem::path list, directory;
  unsigned thread_count = 1;
  std::unordered_set<std::string> options;
  for (int i = 1; i < argc; ++i) {
    const std::string option = argv[i];
    if (option != "--batch-list" && option != "--output-dir" && option != "--threads")
      throw std::runtime_error("unknown batch option: " + option);
    if (!options.insert(option).second)
      throw std::runtime_error("duplicate batch option: " + option);
    if (++i == argc) throw std::runtime_error("missing value for " + option);
    const std::string value = argv[i];
    if (option == "--batch-list") list = value;
    else if (option == "--output-dir") directory = value;
    else {
      if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("--threads must be a positive integer");
      const auto parsed = std::stoull(value);
      if (parsed == 0 || parsed > std::numeric_limits<unsigned>::max())
        throw std::runtime_error("--threads is out of range");
      thread_count = static_cast<unsigned>(parsed);
    }
  }
  if (list.empty() || directory.empty())
    throw std::runtime_error("batch mode requires --batch-list and --output-dir");
  std::ifstream input(list);
  if (!input) throw std::runtime_error("cannot open reference list: " + list.string());
  struct Job { std::filesystem::path input, output; };
  std::vector<Job> jobs;
  std::unordered_set<std::string> names;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    std::filesystem::path reference(line);
    if (reference.is_relative()) reference = list.parent_path() / reference;
    if (!std::filesystem::is_regular_file(reference))
      throw std::runtime_error("not a reference file: " + reference.string());
    const auto name = reference.filename().string() + ".k.17.l.55296.bit";
    if (!names.insert(name).second)
      throw std::runtime_error("duplicate reference sketch filename: " + name);
    const auto output = directory / name;
    if (std::filesystem::exists(output))
      throw std::runtime_error("sketch destination already exists: " + output.string());
    jobs.push_back({reference, output});
  }
  if (!input.eof()) throw std::runtime_error("error reading reference list");
  if (jobs.empty()) throw std::runtime_error("empty reference list");
  std::filesystem::create_directories(directory);
  const auto worker_count = std::min<std::size_t>(thread_count, jobs.size());
  std::atomic<std::size_t> next{0}, completed{0};
  std::atomic<bool> stopped{false};
  std::exception_ptr error;
  std::mutex error_mutex;
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  const auto start = std::chrono::steady_clock::now();
  std::cerr << "batch-sketch start references=" << jobs.size()
            << " workers=" << worker_count << '\n';
  auto worker = [&] {
    while (!stopped.load(std::memory_order_relaxed)) {
      const auto index = next.fetch_add(1, std::memory_order_relaxed);
      if (index >= jobs.size()) break;
      try {
        dna2bit_embedded::sketch_file(jobs[index].input, jobs[index].output);
        completed.fetch_add(1, std::memory_order_relaxed);
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!error) error = std::current_exception();
        stopped.store(true, std::memory_order_relaxed);
      }
    }
  };
  try {
    for (std::size_t i = 0; i < worker_count; ++i) workers.emplace_back(worker);
  } catch (...) {
    stopped.store(true, std::memory_order_relaxed);
    for (auto& thread : workers) thread.join();
    throw;
  }
  for (auto& thread : workers) thread.join();
  if (error) std::rethrow_exception(error);
  if (completed != jobs.size()) throw std::runtime_error("incomplete batch sketch");
  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  std::cerr << "batch-sketch complete references=" << completed
            << " workers=" << worker_count << " seconds=" << elapsed << '\n';
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
    usage(argv[0]);
    return 0;
  }
  if (argc < 3) {
    usage(argv[0]);
    return 2;
  }
  try {
    if (std::string(argv[1]) == "--batch-list") return batch_main(argc, argv);
    std::vector<std::filesystem::path> inputs;
    for (int index = 1; index < argc - 1; ++index) inputs.emplace_back(argv[index]);
    dna2bit_embedded::sketch_files(inputs, argv[argc - 1]);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "embedded-dna2bit-sketch: ERROR: " << error.what() << '\n';
    return 1;
  }
}
