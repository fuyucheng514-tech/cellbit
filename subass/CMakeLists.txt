cmake_minimum_required(VERSION 3.16)
project(cpp_subass_full LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Keep command-line builds fast by default without overriding an explicitly
# selected build type or a multi-config generator (Visual Studio, Xcode, ...).
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(HTS_ROOT "" CACHE PATH "Optional htslib installation prefix")
if(HTS_ROOT)
  # Let FindZLIB select the library from the same self-contained environment.
  list(PREPEND CMAKE_PREFIX_PATH "${HTS_ROOT}")
  set(ZLIB_ROOT "${HTS_ROOT}")
endif()

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(HTS QUIET IMPORTED_TARGET htslib)
endif()

if(TARGET PkgConfig::HTS)
  set(SUBASS_HTS_TARGET PkgConfig::HTS)
else()
  find_path(HTS_INCLUDE_DIR
    NAMES htslib/sam.h
    HINTS "${HTS_ROOT}" ENV HTS_ROOT
    PATH_SUFFIXES include
  )
  find_library(HTS_LIBRARY
    NAMES hts libhts
    HINTS "${HTS_ROOT}" ENV HTS_ROOT
    PATH_SUFFIXES lib lib64
  )
  if(NOT HTS_INCLUDE_DIR OR NOT HTS_LIBRARY)
    message(FATAL_ERROR
      "htslib was not found. Install its pkg-config metadata, or configure "
      "with -DHTS_ROOT=/path/to/htslib-prefix")
  endif()
  add_library(subass_hts UNKNOWN IMPORTED)
  set_target_properties(subass_hts PROPERTIES
    IMPORTED_LOCATION "${HTS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTS_INCLUDE_DIR}"
  )
  set(SUBASS_HTS_TARGET subass_hts)
endif()

add_executable(cpp-subass
  src/main.cpp
  src/bubbles.cpp
  src/no_overlap_policy.cpp
  src/polish_glue.cpp
)
target_link_libraries(cpp-subass PRIVATE ${SUBASS_HTS_TARGET} ZLIB::ZLIB Threads::Threads)

if(MSVC)
  target_compile_options(cpp-subass PRIVATE /W4 /permissive- /Gy /Gw)
  target_link_options(cpp-subass PRIVATE $<$<CONFIG:Release>:/OPT:REF> $<$<CONFIG:Release>:/OPT:ICF>)
else()
  target_compile_options(cpp-subass PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-ffunction-sections>
    $<$<CONFIG:Release>:-fdata-sections>
  )
  if(APPLE)
    target_link_options(cpp-subass PRIVATE $<$<CONFIG:Release>:-Wl,-dead_strip>)
  elseif(UNIX)
    target_link_options(cpp-subass PRIVATE $<$<CONFIG:Release>:-Wl,--gc-sections>)
  endif()
endif()

# IPO/LTO changes neither the algorithms nor floating-point rules used here.
# Enable it when the selected compiler/linker pair supports it, and otherwise
# retain a normal portable Release build.
option(SUBASS_ENABLE_IPO "Enable link-time optimization for Release builds" ON)
if(SUBASS_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT SUBASS_IPO_SUPPORTED OUTPUT SUBASS_IPO_ERROR LANGUAGES CXX)
  if(SUBASS_IPO_SUPPORTED)
    set_property(TARGET cpp-subass PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  else()
    message(STATUS "cpp-subass: IPO/LTO unavailable: ${SUBASS_IPO_ERROR}")
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(cpp-subass-bubble-test src/bubble_test.cpp src/bubbles.cpp)
  target_link_libraries(cpp-subass-bubble-test PRIVATE ${SUBASS_HTS_TARGET} Threads::Threads)
  add_executable(cpp-subass-polish-glue-test tests/polish_glue_test.cpp src/polish_glue.cpp)
  target_include_directories(cpp-subass-polish-glue-test PRIVATE src)
  target_link_libraries(cpp-subass-polish-glue-test PRIVATE ${SUBASS_HTS_TARGET} ZLIB::ZLIB)
  add_executable(cpp-subass-no-overlap-policy-test
    tests/no_overlap_policy_test.cpp
    src/no_overlap_policy.cpp
  )
  target_include_directories(cpp-subass-no-overlap-policy-test PRIVATE src)
  if(MSVC)
    target_compile_options(cpp-subass-bubble-test PRIVATE /W4 /permissive-)
    target_compile_options(cpp-subass-polish-glue-test PRIVATE /W4 /permissive-)
    target_compile_options(cpp-subass-no-overlap-policy-test PRIVATE /W4 /permissive-)
  else()
    target_compile_options(cpp-subass-bubble-test PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cpp-subass-polish-glue-test PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cpp-subass-no-overlap-policy-test PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME cpp-subass-bubble-cli-guard COMMAND cpp-subass-bubble-test)
  set_tests_properties(cpp-subass-bubble-cli-guard PROPERTIES WILL_FAIL TRUE)
  add_test(NAME cpp-subass-polish-glue-cli-guard COMMAND cpp-subass-polish-glue-test)
  set_tests_properties(cpp-subass-polish-glue-cli-guard PROPERTIES WILL_FAIL TRUE)
  add_test(NAME cpp-subass-no-overlap-policy COMMAND cpp-subass-no-overlap-policy-test)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME cpp-subass-no-overlap-cli
      COMMAND ${Python3_EXECUTABLE}
              ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_no_overlap_cli.py
              $<TARGET_FILE:cpp-subass>)
  endif()
endif()
