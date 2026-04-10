# myconote-cli Docker image
#
# Single-stage build on mambaforge: install Rust inside the same image
# to avoid glibc version mismatches between builder and runtime.
#
# Build:
#   docker build -t myconote-cli .
#
# Run:
#   docker run -v $(pwd):/data myconote-cli predict /data/genome.fa --kingdom fungi
#   docker run -v $(pwd):/data myconote-cli annotate /data/genes.gff3 --fasta /data/genome.fa
#
# With databases pre-mounted:
#   docker run -v $(pwd):/data -v ~/.myconote/dbs:/root/.myconote/dbs myconote-cli annotate ...

FROM condaforge/mambaforge:latest

LABEL maintainer="Benjamin Narh-Madey <narhmadey@wisc.edu>"
LABEL org.opencontainers.image.description="myconote-cli: Blazing-fast genome annotation pipeline"
LABEL org.opencontainers.image.version="0.1.0"
LABEL org.opencontainers.image.source="https://github.com/K-nie/myconote-cli"
LABEL org.opencontainers.image.licenses="MIT"

# Install build dependencies and compile myconote-cli
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl build-essential pkg-config libfontconfig1-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Build myconote-cli
WORKDIR /build
COPY Cargo.toml Cargo.lock* ./
COPY src/ src/
COPY tests/ tests/

RUN cargo build --release && \
    strip target/release/myconote-cli && \
    cp target/release/myconote-cli /usr/local/bin/myconote-cli && \
    rm -rf /build /root/.cargo/registry /root/.cargo/git && \
    rustup self uninstall -y || true

# Verify binary runs
RUN myconote-cli --version

# ── Install bioinformatics tools in groups ──────────────────────────────────

# Group 1: Core alignment & sequence tools
RUN mamba install -y -c bioconda -c conda-forge \
    minimap2 samtools hmmer diamond blast mafft muscle \
    && mamba clean -afy

# Group 2: Gene predictors
RUN mamba install -y -c bioconda -c conda-forge \
    augustus snap glimmerhmm \
    && mamba clean -afy

# Group 3: Search & specialized annotation
RUN mamba install -y -c bioconda -c conda-forge \
    mmseqs2 miniprot trnascan-se \
    && mamba clean -afy

# Group 4: Repeat masking
RUN mamba install -y -c bioconda -c conda-forge \
    repeatmasker repeatmodeler \
    && mamba clean -afy

# Group 5: Phylogenetics
RUN mamba install -y -c bioconda -c conda-forge \
    iqtree fasttree \
    && mamba clean -afy

# Group 6: Heavy tools (large dependency trees)
RUN mamba install -y -c bioconda -c conda-forge \
    busco && mamba clean -afy

RUN mamba install -y -c bioconda -c conda-forge \
    eggnog-mapper && mamba clean -afy

RUN mamba install -y -c bioconda -c conda-forge \
    trinity transdecoder && mamba clean -afy \
    || echo "WARNING: Trinity/TransDecoder not available for this platform"

RUN mamba install -y -c bioconda -c conda-forge \
    pasa && mamba clean -afy \
    || echo "WARNING: PASA not available for this platform"

# Group 7: pip tools
RUN pip install --no-cache-dir pybiolib

# Create database directory
RUN mkdir -p /root/.myconote/dbs

WORKDIR /data

ENTRYPOINT ["myconote-cli"]
CMD ["--help"]
