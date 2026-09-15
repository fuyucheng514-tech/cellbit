# DNA2bit-SAG Original Embedded

This is a new, write-once derivative of the first unmodified speed-v1 workflow.
It keeps the original k=17, bit_len=55296, hash_type=0 sketch and strict margin >0.01 semantics.

The teacher DNA2bit executable is not launched by the pipeline. The archived DNA2bit sketch source is compiled into dna2bit-sag-pipeline; the audited packed-search implementation is linked in-process as well. dna2bit-embedded-sketch is a small standalone self-test/utility using the same embedded source.

The original Stage3A/Stage3B runner and scientific rules remain separate and are not silently changed. External tools that remain external are documented in EXTERNAL_COMPONENTS.tsv.

Package name: dna2bit_sag_original_integrated_20260909
