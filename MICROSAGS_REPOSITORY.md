# Microsags repository

This repository is the source distribution of Microsags (DNA2bit-SAG Tractor original integrated release).

## Contents

- `src/`, `subass/`: C++ source and embedded DNA2bit interfaces.
- `python/`: Stage3B helper scripts.
- `tools/`: preparation, audit, and reproduction utilities.
- `docs/`: workflow and reproducibility documentation.
- `examples/`: minimal manifest examples.
- `config/`: external database path template.
- `website/`: static project documentation site.

Large sequencing inputs, result trees, databases, build directories, and credentials are intentionally not included.

## Build

```bash
conda env create -f environment.yml
conda activate dna2bit-sag-original
JOBS=8 PREFIX="$HOME/.local/microsags" bash install.sh
```

See `README_RELEASE.md` and `config/paths.env.example` for the complete installation and database configuration contract.

The same static site files are mirrored at docs/ so GitHub Pages can publish from the main branch /docs folder.
