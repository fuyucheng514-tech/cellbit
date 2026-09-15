cmake_minimum_required(VERSION 3.16)
project(dna2bit_sag_pipeline VERSION 0.1.0 LANGUAGES CXX)

set(HTS_ROOT "" CACHE PATH "Optional htslib/zlib installation prefix for bundled cpp-subass")
if(HTS_ROOT)
  list(PREPEND CMAKE_PREFIX_PATH "${HTS_ROOT}")
  set(ZLIB_ROOT "${HTS_ROOT}")
endif()

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# This repository is a production command-line package.  A bare CMake
# invocation should therefore produce the optimized build; callers can still
# request Debug explicitly.  We deliberately avoid -march=native so that the
# packaged binaries remain portable between x86-64 servers.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(DNA2BIT_ENABLE_IPO "Enable supported interprocedural optimization in Release builds" ON)
if(DNA2BIT_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT DNA2BIT_IPO_SUPPORTED OUTPUT DNA2BIT_IPO_ERROR)
endif()

function(dna2bit_configure_target target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive- /Gy /Gw)
    target_link_options(${target} PRIVATE
      $<$<CONFIG:Release>:/OPT:REF>
      $<$<CONFIG:Release>:/OPT:ICF>
      $<$<CONFIG:MinSizeRel>:/OPT:REF>
      $<$<CONFIG:MinSizeRel>:/OPT:ICF>)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE
      $<$<CONFIG:Release>:-ffunction-sections;-fdata-sections>
      $<$<CONFIG:MinSizeRel>:-ffunction-sections;-fdata-sections>)
  endif()
  if(APPLE)
    target_link_options(${target} PRIVATE
      $<$<CONFIG:Release>:-Wl,-dead_strip>
      $<$<CONFIG:MinSizeRel>:-Wl,-dead_strip>)
  elseif(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_options(${target} PRIVATE
      $<$<CONFIG:Release>:-Wl,--gc-sections;-Wl,--as-needed>
      $<$<CONFIG:MinSizeRel>:-Wl,--gc-sections;-Wl,--as-needed>)
  endif()
  if(DNA2BIT_ENABLE_IPO AND DNA2BIT_IPO_SUPPORTED)
    set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL TRUE)
  endif()
endfunction()

add_executable(dna2bit-sag-pipeline src/main.cpp src/dna2bit_embedded.cpp src/dna2bit_packed_api.cpp)
add_executable(dna2bit-embedded-sketch src/dna2bit_embedded_cli.cpp src/dna2bit_embedded.cpp)
target_link_libraries(dna2bit-embedded-sketch PRIVATE ZLIB::ZLIB Threads::Threads)
dna2bit_configure_target(dna2bit-embedded-sketch)
target_link_libraries(dna2bit-sag-pipeline PRIVATE ZLIB::ZLIB Threads::Threads)
dna2bit_configure_target(dna2bit-sag-pipeline)

add_executable(sag-stage3b-tractor src/stage3b.cpp)
dna2bit_configure_target(sag-stage3b-tractor)

add_executable(gtdb-ani-af src/gtdb_ani_af.cpp)
target_link_libraries(gtdb-ani-af PRIVATE ZLIB::ZLIB Threads::Threads)
dna2bit_configure_target(gtdb-ani-af)

# Drop-in search-side compatibility implementation for the teacher Dna2bit
# binary.  It freezes the legacy opendir/readdir order once and stores all
# fixed-width .bit rows in one contiguous file, avoiding ~200k tiny opens on
# every pipeline invocation.  It is intentionally a separate executable
# until full-dataset byte-for-byte equivalence has been signed off.
add_executable(dna2bit-packed src/dna2bit_packed.cpp)
target_link_libraries(dna2bit-packed PRIVATE Threads::Threads)
dna2bit_configure_target(dna2bit-packed)

# Explicit production-facing entry points.  They share the exact same core as
# the combined developer CLI but do not require a build/search subcommand.
add_executable(dna2bit-pack-builder src/dna2bit_packed.cpp)
target_compile_definitions(dna2bit-pack-builder PRIVATE DNA2BIT_PACKED_BUILDER_ONLY=1)
target_link_libraries(dna2bit-pack-builder PRIVATE Threads::Threads)
dna2bit_configure_target(dna2bit-pack-builder)

add_executable(dna2bit-packed-search src/dna2bit_packed.cpp)
target_compile_definitions(dna2bit-packed-search PRIVATE DNA2BIT_PACKED_SEARCH_ONLY=1)
target_link_libraries(dna2bit-packed-search PRIVATE Threads::Threads)
dna2bit_configure_target(dna2bit-packed-search)

# The deployed teacher Dna2bit binary itself requires x86 POPCNT.  Enable the
# same ISA only for the three packed-search executables so each 64-bit XOR is
# counted by the hardware instruction rather than a libgcc helper, without
# raising the ISA requirement of the rest of this portable package.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
  target_compile_options(dna2bit-packed PRIVATE -mpopcnt)
  target_compile_options(dna2bit-pack-builder PRIVATE -mpopcnt)
  target_compile_options(dna2bit-packed-search PRIVATE -mpopcnt)
endif()

option(DNA2BIT_BUILD_SUBASS "Build the bundled C++ Flye subassembly driver" ON)
if(DNA2BIT_BUILD_SUBASS)
  # Register the subdirectory's CTest entries in the top-level build tree.
  # Without a root CTestTestfile.cmake, invoking `ctest --test-dir build`
  # silently reports "No tests were found" even though subass defines tests.
  include(CTest)
  add_subdirectory(subass)
endif()

if(BUILD_TESTING AND UNIX)
  add_test(NAME gtdb_ani_af_rebind_taxonomy
           COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/gtdb_ani_af_rebind_taxonomy.sh
                        $<TARGET_FILE:gtdb-ani-af>)
  set_tests_properties(gtdb_ani_af_rebind_taxonomy PROPERTIES TIMEOUT 120)
endif()

install(TARGETS dna2bit-sag-pipeline dna2bit-embedded-sketch sag-stage3b-tractor gtdb-ani-af dna2bit-packed
                dna2bit-pack-builder dna2bit-packed-search
        RUNTIME DESTINATION bin)
if(TARGET cpp-subass)
  install(TARGETS cpp-subass RUNTIME DESTINATION bin)
endif()
