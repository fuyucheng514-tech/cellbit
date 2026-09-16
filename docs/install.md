# Install

Microsags supports Linux x86-64. Conda or Mamba is recommended.

```bash
git clone https://github.com/fuyucheng514-tech/cellbit.git Microsags
cd Microsags
conda env create -f environment.yml
conda activate microsags
PREFIX="$CONDA_PREFIX" JOBS=8 bash install.sh
microsags --help
```

## Download the database

The database is distributed as a GitHub Release asset:

```bash
mkdir -p "$HOME/microsags-data"
cd "$HOME/microsags-data"
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
tar -xzf Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
```

Use the extracted directory with `-d`:

```bash
microsags annotate SAGs/ -d database -o annotation_output
```

`assemble` also requires CheckM2 and GTDB-Tk data:

```bash
microsags assemble SAGs/ -d database -o assembly_output \
  --checkm2-database checkm2_database \
  --gtdbtk-data gtdbtk_database
```
