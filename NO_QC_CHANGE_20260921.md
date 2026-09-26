# Microsags v0.5.0 direct-read annotation contract

Microsags annotation does not execute a read-preprocessing stage.

- Paired FASTQ: R1 and R2 are paired by SAG name and the original files are
  passed directly to the embedded DNA2bit sketcher.
- Singleton FASTQ: the original file is passed directly to the embedded
  DNA2bit sketcher.
- Contig FASTA: the existing 1,000 bp total-length eligibility gate is retained.
- Preflight checks only path existence, filename extension, unique SAG IDs, and
  valid R1/R2 pairing. It does not open or decompress sequence files.
- Sequence readability is checked naturally when DNA2bit performs its first and
  only read of the input.
- No transformed FASTQ or read-preprocessing report/pass marker is produced.
- The public command-line interface has no read-preprocessing executable option.

The scientific DNA2bit parameters are unchanged: k=17, bit length=55,296 and
hash type=0. Assembly mode and its Stage 3A/3B algorithms are unchanged.
