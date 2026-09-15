#include "dna2bit_embedded_api.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " INPUT [INPUT ...] OUTPUT\n";
    return 2;
  }
  try {
    std::vector<std::filesystem::path> inputs;
    for (int index = 1; index < argc - 1; ++index) inputs.emplace_back(argv[index]);
    dna2bit_embedded::sketch_files(inputs, argv[argc - 1]);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "embedded-dna2bit-sketch: ERROR: " << error.what() << '\n';
    return 1;
  }
}
