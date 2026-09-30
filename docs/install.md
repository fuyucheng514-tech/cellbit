# Install

Microsags runs on Linux x86-64. Install the program and download its reference
database once. Annotation needs nothing else.

## 1. Install Microsags

Download the prebuilt package; no compiler or Conda installation is needed:

```bash
mkdir -p "$HOME/.local/microsags"
cd "$HOME/.local/microsags"
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.5.2-linux1/Microsags-v0.5.2-linux-x86_64.tar.gz
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.5.2-linux1/Microsags-v0.5.2-linux-x86_64.tar.gz.sha256
sha256sum -c Microsags-v0.5.2-linux-x86_64.tar.gz.sha256
tar -xzf Microsags-v0.5.2-linux-x86_64.tar.gz
export PATH="$HOME/.local/microsags/bin:$PATH"
conda-unpack
microsags --help
```

`conda-unpack` is included in the package; run it once after extraction. To use
`microsags` in future terminal sessions, add the `export PATH=...` line to your
shell startup file. Do not move the extracted directory after this step.

## 2. Download the database

```bash
mkdir -p "$HOME/microsags-data"
cd "$HOME/microsags-data"
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
wget https://github.com/fuyucheng514-tech/cellbit/releases/download/v0.1.0/Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
sha256sum -c Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz.sha256
tar -xzf Microsags-GTDB232-DNA2bit-k17-packed-v1.tar.gz
```

Pass the extracted database directory with `-d` in the
[annotation examples](usage.md). If the server cannot reach GitHub, download
both files on another machine, copy them to the server, and run the checksum
and extraction commands there.

## 3. Assembly dependencies (optional)

Skip this section if you only need annotation. Assembly also calls CheckM2 and
GTDB-Tk, which each need their own reference database. Install those tools in
separate environments:

```bash
mamba create -n microsags-checkm2 -c conda-forge -c bioconda checkm2=1.0.1
mamba create -n microsags-gtdbtk -c conda-forge -c bioconda gtdbtk=2.7.2
```

Create the configuration directory with
`mkdir -p "$HOME/.config/microsags"`, then create
`$HOME/.config/microsags/paths.env` with the actual locations of the two
executables and their downloaded databases:

```text
CHECKM2_BIN=/absolute/path/to/checkm2
CHECKM2DB=/absolute/path/to/checkm2_database.dmnd
GTDBTK_BIN=/absolute/path/to/gtdbtk
GTDBTK_DATA_PATH=/absolute/path/to/gtdbtk_data
```

You then run Assembly without repeating these four paths. To compile
Microsags from source or build a database for a new GTDB release, see the
[GitHub README](https://github.com/fuyucheng514-tech/cellbit#install) and
[database-building section](https://github.com/fuyucheng514-tech/cellbit#build-a-database).
