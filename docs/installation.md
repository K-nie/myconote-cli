# Installation

## Prerequisites

- **macOS** (Apple Silicon or Intel) or **Linux** (x86-64)
- [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html)
- Rust ≥ 1.74 (installed automatically if missing)
- ~20 GB disk space for core tools; ~75 GB if downloading eggNOG databases

---

## Option A: One-shot installer (recommended)

```bash
git clone https://github.com/K-nie/myconote-cli.git
cd myconote-cli
bash install.sh
```

This script handles everything:

1. Installs Miniconda if conda is not found
2. Creates a dedicated `myconote` conda environment
3. Installs all bioinformatics dependencies (Augustus, RepeatMasker, IQ-TREE2, etc.)
4. Compiles the Myconote_CLI binary with `cargo build --release`
5. Adds `myconote` to your PATH

Estimated time: 20–40 minutes on first install (database downloads excluded).

---

## Option B: Manual conda + cargo

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate myconote

# 2. Build binary
cargo build --release

# 3. Add to PATH
export PATH="$PWD/target/release:$PATH"
```

---

## Option C: Pre-built binary

Download the latest release from the [GitHub Releases page](https://github.com/K-nie/myconote-cli/releases), extract, and place the `myconote` binary somewhere on your PATH.

---

## Verify installation

```bash
myconote --version
myconote check
```

`myconote check` reports which external tools are found and which are missing.

---

## Database setup

```bash
myconote db setup --kingdom fungi
```

Downloads and indexes required databases (dbCAN, MEROPS, GO, Pfam). The `--kingdom fungi` flag fetches the minimal set; omit it for the full suite including eggNOG (~50 GB).

!!! warning "eggNOG database"
    The eggNOG database is ~50 GB. Download only if you need functional annotation beyond BLAST + dbCAN.
