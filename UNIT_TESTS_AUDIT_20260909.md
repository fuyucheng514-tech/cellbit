# Unit-test audit — DNA2bit-SAG Original Integrated v1

Date (UTC): 2026-09-09T14:03:32Z

The CMake test targets were explicitly built before execution. All five registered tests passed:

Internal ctest changing into directory: /home/data/fyc/dna2bit_sag_original_integrated_20260909/build_embedded_original
Test project /home/data/fyc/dna2bit_sag_original_integrated_20260909/build_embedded_original
    Start 1: gtdb_ani_af_rebind_taxonomy
1/5 Test #1: gtdb_ani_af_rebind_taxonomy ........   Passed    0.24 sec
    Start 2: cpp-subass-bubble-cli-guard
2/5 Test #2: cpp-subass-bubble-cli-guard ........   Passed    0.00 sec
    Start 3: cpp-subass-polish-glue-cli-guard
3/5 Test #3: cpp-subass-polish-glue-cli-guard ...   Passed    0.00 sec
    Start 4: cpp-subass-no-overlap-policy
4/5 Test #4: cpp-subass-no-overlap-policy .......   Passed    0.01 sec
    Start 5: cpp-subass-no-overlap-cli
5/5 Test #5: cpp-subass-no-overlap-cli ..........   Passed    0.33 sec

100% tests passed, 0 tests failed out of 5

Total Test time (real) =   0.59 sec

Test binary SHA256:

4c1b26eedfe5e961c01eefb3a9257260c79a110d955691ba032c30419e51aec0  /home/data/fyc/dna2bit_sag_original_integrated_20260909/build_embedded_original/subass/cpp-subass-bubble-test
ac7672f2f9ee02abc0f92816c1c644f12f4a147a5aaabc2cc242fdd04cf067c9  /home/data/fyc/dna2bit_sag_original_integrated_20260909/build_embedded_original/subass/cpp-subass-polish-glue-test
e8f82be5b6a7fd6a02b43cf079ccb3ede8bea1ae9cc3c14c8c259a22bffa839f  /home/data/fyc/dna2bit_sag_original_integrated_20260909/build_embedded_original/subass/cpp-subass-no-overlap-policy-test

The initial ctest invocation reported three tests as Not Run because those optional test executables had not been built by the targeted production build. This was a build/test invocation gap, not a code failure; the explicit test-target build closed it and the complete suite is now 5/5 PASS.
