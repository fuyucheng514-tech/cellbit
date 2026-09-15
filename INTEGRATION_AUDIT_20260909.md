# DNA2bit-SAG Original Integrated v1 — final integration audit

Date (UTC): 2026-09-09T13:42:31Z
Package: /home/data/fyc/dna2bit_sag_original_integrated_20260909
Scope: source integration and small real-contig end-to-end acceptance smoke; no 13,742-SAG production rerun.

## Scientific contract

- Frozen semantic lineage: original v1 pipeline, not Cellbit57 and not a speed-v5 scientific redesign.
- Embedded DNA2bit parameters: k=17, bit_len=55296, hash_type=0, reverse-complement canonical sketch.
- Search semantics: teacher-compatible packed index, strict min_ratio > 0.01, original winner/tie/record-low behavior retained in the packed implementation.
- The existing <1000 bp hard gate, Stage3A species aggregation and C++ subassemble/Flye handoff are unchanged by this integration.

## Integration checks

1. Main pipeline links src/dna2bit_embedded.cpp and src/dna2bit_packed_api.cpp in CMake.
2. The main executable calls dna2bit_embedded::sketch_files and dna2bit_packed::search_index_for_pipeline in-process.
3. No external teacher DNA2bit executable or dna2bit-packed-search process is referenced by the executable path or the formal runner.
4. The standalone packed-search executable remains available only as an audit oracle; it is not required by the main pipeline.
5. Main binary SHA256: 735f66cf351ef80a461012fbbfed5ff306d6075303645d467919b14157fad5fa
6. Embedded sketch CLI SHA256: 253f5304ff0cc8f1ee70db63a36061a11d6a712b1c7ac3d6f1a80f465dcb63cd
7. Standalone packed-search oracle SHA256: 2b511e4a2b27c41aa69d21998af1a6dfd308de1d1604a7af2cf4d0af56358bec

## Functional evidence

- Teacher-vs-embedded FASTA sketch: byte-identical (validated before this audit).
- Teacher-vs-embedded paired FASTQ sketch: byte-identical (validated before this audit).
- In-process packed search vs standalone packed-search oracle: byte-identical; smoke result SHA256 69f6c946d4e46bc431199b99d9c7c09c6cab69d2aa5ca4703b9cd4350fed9cd8.
- Smoke input: 4 real contig SAGs from the authoritative Lake contig manifest.
- Smoke closure: 4 input, 4 eligible, 0 excluded below 1000 bp, 2 DNA2bit labels, 2 pending, 2 Stage3A groups, 2 assembled bin FASTAs, overall COMPLETE status PASS.
- Smoke resource/time evidence: /home/data/fyc/dna2bit_sag_original_integrated_smoke_20260909/attempt_003/MAIN_TIME.txt; full stdout/stderr and all intermediate files are preserved under the write-once smoke root.

## Defects found and corrected during audit

- Defect 1: groups.tsv streamed std::filesystem::path values with implementation-added quotes. This could make downstream path validation fail. Corrected by writing .string() values; rebuilt main SHA is recorded above.
- Defect 2: top-level 02_dna2bit/SKETCH.PASS omitted its implementation identity. Corrected by adding implementation=embedded_dna2bit_source; per-SAG receipts already carried the identity.
- The first two smoke attempts are intentionally preserved as failed audit evidence under /home/data/fyc/dna2bit_sag_original_integrated_smoke_20260909/attempt_001 and attempt_002; no output was overwritten.

## External components (not silently claimed as embedded)

The original v1 wrapper still invokes the external fastp/SPAdes path for paired reads and the configured Flye modules through the C++ subassemble wrapper. CheckM2, GTDB-Tk and BLAST+ remain external evaluators/helpers in the full wrapper. They are listed in EXTERNAL_COMPONENTS.tsv. Only DNA2bit sketching and packed search were integrated in-process in this package.

## Scientific acceptance boundary

This audit establishes that the new executable is a reproducible, teacher-compatible in-process integration and that the contig route completes on real data without the discovered receipt/path bugs. It does not replace a fresh 13,742-SAG full scientific run; the fixed full runner remains write-once and retains its original closure assertions and authority hashes.
