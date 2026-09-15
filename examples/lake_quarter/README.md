# Lake quarter-scale integration example

`select_quarter.py` deterministically selects every fourth SAG after sorting the
two-column contig manifest by `sag_id`. For the audited 13,742-SAG lake
manifest, it selects 3,436 SAGs.

`verified_result.json` records the result closure, timings and core output
checksums from the completed public-release integration run.

The biological FASTA files are not redistributed in this repository. See the
[Lake quarter tutorial](../../docs/lake-quarter-tutorial.md) for the tested
manifest checksums, command and output interpretation.
