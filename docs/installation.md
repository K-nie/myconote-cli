# Installation

## Prerequisites

- **macOS** (Apple Silicon or Intel) or **Linux** (x86-64)
- [Miniconda or Anaconda](https://docs.conda.io/en/latest/miniconda.html)
- ~20 GB disk space for core tools; ~75 GB if downloading eggNOG databases

---

## Option A: Quick install (pre-built binary)

Download a pre-built binary for your platform. No Rust compiler needed.

```bash
curl -fsSL https://raw.githubusercontent.com/K-nie/myconote-cli/main/quick-install.sh | bash
```

Then install external tools and databases:
```bash
myconote-cli install --yes   # installs 30+ bioinformatics tools via conda
myconote-cli setup           # downloads annotation databases (~2.5 GB)
```

---

## Option B: Full installer (recommended for first-time setup)

Handles everything: Miniconda, Rust, binary build, tool installation, and database downloads.

```bash
git clone https://github.com/K-nie/myconote-cli.git
cd myconote-cli
bash install.sh
```

Estimated time: 20-40 minutes on first install.

---

## Option C: Docker

```bash
docker pull ghcr.io/k-nie/myconote-cli:latest
docker run -v $(pwd):/data myconote-cli predict /data/genome.fa --kingdom fungi
```

All tools and databases are pre-installed in the container.

---

## Option D: Singularity / Apptainer (HPC)

```bash
singularity pull myconote-cli.sif docker://ghcr.io/k-nie/myconote-cli:latest
singularity run myconote-cli.sif predict genome.fa --kingdom fungi
```

---

## Option E: Build from source

If you have Rust >= 1.74:

```bash
git clone https://github.com/K-nie/myconote-cli.git
cd myconote-cli
cargo build --release
cp target/release/myconote-cli ~/.local/bin/
myconote-cli install --yes
myconote-cli setup
```

---

## Verify installation

```bash
myconote-cli --version       # shows version and banner
myconote-cli check           # reports status of all 30+ external tools
myconote-cli setup --check   # shows database download status
```

---

## Database setup

```bash
myconote-cli setup                        # download all core databases
myconote-cli setup --db swiss-prot pfam   # download specific databases
myconote-cli setup --check                # verify what's installed
```

Core databases (~2.5 GB):

| Database | Size | Used by |
|----------|------|---------|
| Swiss-Prot (MMseqs2) | ~1 GB | `annotate` (product names) |
| Pfam-A HMMs | ~300 MB | `annotate` (domain search) |
| BUSCO lineages | ~500 MB | `annotate` (completeness) |
| dbCAN (CAZyme) | ~200 MB | `annotate --cazyme` |
| MEROPS | ~100 MB | `annotate --merops` |

!!! warning "eggNOG database"
    The eggNOG database is ~50 GB and is NOT downloaded by default. Download only if you need COG/KEGG annotation:
    ```bash
    download_eggnog_data.py -y --data_dir ~/.eggnog_mapper/data
    ```

---

## Ollama setup (for `explain` AI interpreter)

The `explain` command uses a local LLM via Ollama for privacy-first result interpretation. No data leaves your machine.

```bash
myconote-cli setup ollama        # install Ollama, detect RAM, pull best model
```

This will:

1. Install Ollama (latest version from GitHub releases)
2. Start the Ollama server
3. Auto-detect your system RAM and select the most capable model
4. Pull and verify the model

Available model tiers:

| RAM | Model | Quality |
|-----|-------|---------|
| 48+ GB | llama3.3:70b-instruct-q4_K_M | Best |
| 24+ GB | qwen2.5:32b-instruct-q4_K_M | Excellent |
| 16+ GB | mistral-small:22b | Strong |
| 12+ GB | qwen2.5:14b | Good |
| 8+ GB | llama3.1:8b | Baseline |

Override with `MYCONOTE_CHAT_MODEL=<model>` or `--model <model>` on the `explain` command.

Optionally download Q1 open-access papers for grounded citations:

```bash
myconote-cli setup chat-corpus   # ~50 MB, CC-BY licensed papers only
```

---

## Getting started

After installation:

```bash
myconote-cli learn    # interactive tutorial (8 lessons, ~1 hour)
```

Or jump straight to the [Quick Start](quickstart.md) guide.
