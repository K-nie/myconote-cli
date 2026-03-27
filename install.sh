#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  myconote-cli — one-shot installer
#
#  What this script does:
#    1. Detects OS (macOS / Linux)
#    2. Installs Miniconda if conda is not found
#    3. Creates a dedicated 'myconote' conda environment
#    4. Installs Rust (via rustup) if missing or too old (need ≥ 1.74)
#    5. Builds and installs the myconote-cli binary
#    6. Calls `myconote-cli install` to install all bioinformatics tools
#    7. Downloads annotation databases (~2.5 GB)
#    8. Adds a shell alias so you can run 'myconote-cli' from anywhere
#
#  Usage:
#    bash install.sh
#
#  Skip individual steps with environment variables:
#    SKIP_CONDA=1       bash install.sh   # skip Miniconda install
#    SKIP_RUST=1        bash install.sh   # skip Rust install / version check
#    SKIP_BUILD=1       bash install.sh   # skip building the Rust binary
#    SKIP_TOOLS=1       bash install.sh   # skip bioinformatics tool install
#    SKIP_DATABASES=1   bash install.sh   # skip database downloads
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

ok()   { echo -e "${GREEN}  ✓${RESET}  $*"; }
info() { echo -e "${CYAN}  →${RESET}  $*"; }
warn() { echo -e "${YELLOW}  ⚠${RESET}  $*"; }
fail() { echo -e "${RED}  ✗${RESET}  $*"; exit 1; }
step() { echo -e "\n${BOLD}${CYAN}══ $* ══${RESET}"; }

# ── Configuration ─────────────────────────────────────────────────────────────
ENV_NAME="myconote"
MINICONDA_DIR="$HOME/miniconda3"
DB_DIR="$HOME/.myconote/dbs"
INSTALL_DIR="$HOME/.local/bin"
MIN_RUST_VERSION="1.74"   # driven by clap 4.5

# Detect OS + architecture
OS="$(uname -s)"
ARCH="$(uname -m)"
case "$OS" in
    Linux)  PLATFORM="Linux"  ;;
    Darwin) PLATFORM="MacOSX" ;;
    *)      fail "Unsupported OS: $OS  (only Linux and macOS are supported)" ;;
esac

# ── Banner ────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}╔═══════════════════════════════════════════════════════╗${RESET}"
echo -e "${BOLD}║         myconote-cli  —  Full Installer               ║${RESET}"
echo -e "${BOLD}╚═══════════════════════════════════════════════════════╝${RESET}"
echo ""
echo -e "  Platform : ${BOLD}${PLATFORM} (${ARCH})${RESET}"
echo -e "  Env name : ${BOLD}${ENV_NAME}${RESET}"
echo -e "  Databases: ${BOLD}${DB_DIR}${RESET}"
echo ""
echo -e "  This will install:"
echo -e "    • Miniconda3             (if not already present)"
echo -e "    • Rust ≥ ${MIN_RUST_VERSION}              (via rustup — no root required)"
echo -e "    • myconote-cli binary    (compiled from source)"
echo -e "    • All 27 bioinformatics tools  (via conda/mamba + pip)"
echo -e "    • Annotation databases   (~2.5 GB)"
echo ""
read -r -p "  Continue? [y/N] " CONFIRM
[[ "$CONFIRM" =~ ^[Yy]$ ]] || { echo "  Aborted."; exit 0; }

# ─────────────────────────────────────────────────────────────────────────────
# Helper: compare version strings (e.g. "1.74" vs "1.70.0")
# Returns 0 (true) if $1 >= $2
# ─────────────────────────────────────────────────────────────────────────────
version_gte() {
    # Strip any suffix like "-nightly" and compare numerically
    local a b
    IFS='.' read -ra a <<< "$(echo "$1" | tr -cs '0-9.' '.' | sed 's/^\.//' | sed 's/\.$//')"
    IFS='.' read -ra b <<< "$(echo "$2" | tr -cs '0-9.' '.' | sed 's/^\.//' | sed 's/\.$//')"
    local i
    for (( i=0; i<${#b[@]}; i++ )); do
        local av="${a[$i]:-0}"
        local bv="${b[$i]:-0}"
        (( 10#$av > 10#$bv )) && return 0
        (( 10#$av < 10#$bv )) && return 1
    done
    return 0  # equal
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Install Miniconda if needed
# ─────────────────────────────────────────────────────────────────────────────
step "Step 1/6 — Conda"

if [[ "${SKIP_CONDA:-0}" == "1" ]]; then
    warn "SKIP_CONDA=1, skipping Miniconda install"
elif command -v conda &>/dev/null; then
    ok "conda already installed: $(conda --version)"
else
    info "conda not found — installing Miniconda3 to ${MINICONDA_DIR}"

    case "${PLATFORM}-${ARCH}" in
        Linux-x86_64)   MC_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;;
        Linux-aarch64)  MC_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh" ;;
        MacOSX-x86_64)  MC_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh" ;;
        MacOSX-arm64)   MC_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh" ;;
        *) fail "No Miniconda installer for ${PLATFORM}-${ARCH}" ;;
    esac

    MC_INSTALLER="/tmp/miniconda_install.sh"
    info "Downloading: $MC_URL"
    curl -fsSL "$MC_URL" -o "$MC_INSTALLER"
    bash "$MC_INSTALLER" -b -p "$MINICONDA_DIR"
    rm -f "$MC_INSTALLER"

    # shellcheck source=/dev/null
    source "${MINICONDA_DIR}/etc/profile.d/conda.sh"

    for RC in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.zshrc"; do
        if [[ -f "$RC" ]] && ! grep -q "miniconda3/etc/profile.d/conda.sh" "$RC" 2>/dev/null; then
            echo "" >> "$RC"
            echo "# >>> conda initialize >>>" >> "$RC"
            echo "source \"${MINICONDA_DIR}/etc/profile.d/conda.sh\"" >> "$RC"
            echo "# <<< conda initialize <<<" >> "$RC"
            ok "Added conda init to $RC"
        fi
    done

    ok "Miniconda3 installed at ${MINICONDA_DIR}"
fi

CONDA_BASE="$(conda info --base 2>/dev/null || echo "$MINICONDA_DIR")"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh" 2>/dev/null || true

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Create the myconote conda environment
# ─────────────────────────────────────────────────────────────────────────────
step "Step 2/6 — Conda environment '${ENV_NAME}'"

if conda env list | grep -q "^${ENV_NAME} "; then
    warn "Environment '${ENV_NAME}' already exists — reusing it"
    warn "To start fresh: conda env remove -n ${ENV_NAME} && bash install.sh"
else
    info "Creating Python 3.11 environment…"
    conda create -y -n "$ENV_NAME" python=3.11
    ok "Environment '${ENV_NAME}' created"
fi

conda activate "$ENV_NAME" 2>/dev/null || true
CONDA_ENV_BIN="$(conda run -n "$ENV_NAME" bash -c 'echo $CONDA_PREFIX' 2>/dev/null)/bin"
export PATH="${CONDA_ENV_BIN}:$PATH"
ok "Environment active (bin: ${CONDA_ENV_BIN})"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Install Rust (with minimum version enforcement)
# ─────────────────────────────────────────────────────────────────────────────
step "Step 3/6 — Rust toolchain (need ≥ ${MIN_RUST_VERSION})"

if [[ "${SKIP_RUST:-0}" == "1" ]]; then
    warn "SKIP_RUST=1, skipping Rust install / version check"
else
    # Source cargo env if it exists (covers previous installs in the same session)
    # shellcheck source=/dev/null
    [[ -f "$HOME/.cargo/env" ]] && source "$HOME/.cargo/env"
    export PATH="$HOME/.cargo/bin:$PATH"

    NEED_INSTALL=false
    NEED_UPDATE=false

    if command -v rustc &>/dev/null; then
        # Parse version string: "rustc 1.77.2 (25ef9e3d8 2024-04-09)"
        RUSTC_VER="$(rustc --version 2>/dev/null | awk '{print $2}')"
        if version_gte "$RUSTC_VER" "$MIN_RUST_VERSION"; then
            ok "Rust ${RUSTC_VER} — meets requirement (≥ ${MIN_RUST_VERSION})"
        else
            warn "Rust ${RUSTC_VER} is too old — minimum required is ${MIN_RUST_VERSION}"
            if command -v rustup &>/dev/null; then
                NEED_UPDATE=true
            else
                warn "rustup not found; will reinstall Rust from scratch"
                NEED_INSTALL=true
            fi
        fi
    else
        info "Rust not found — will install via rustup"
        NEED_INSTALL=true
    fi

    if [[ "$NEED_UPDATE" == "true" ]]; then
        info "Updating Rust to latest stable via rustup…"
        rustup update stable
        # shellcheck source=/dev/null
        source "$HOME/.cargo/env"
        RUSTC_VER="$(rustc --version 2>/dev/null | awk '{print $2}')"
        ok "Rust updated to ${RUSTC_VER}"
    fi

    if [[ "$NEED_INSTALL" == "true" ]]; then
        info "Downloading and installing Rust via rustup (no root required)…"
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
            | sh -s -- -y --no-modify-path --default-toolchain stable
        # shellcheck source=/dev/null
        source "$HOME/.cargo/env"
        export PATH="$HOME/.cargo/bin:$PATH"
        RUSTC_VER="$(rustc --version 2>/dev/null | awk '{print $2}')"
        ok "Rust ${RUSTC_VER} installed"

        # Persist cargo env in shell rc files
        for RC in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.zshrc"; do
            if [[ -f "$RC" ]] && ! grep -q '\.cargo/env' "$RC" 2>/dev/null; then
                echo "" >> "$RC"
                echo '# Rust / cargo' >> "$RC"
                echo '. "$HOME/.cargo/env"' >> "$RC"
                ok "Added cargo env to $RC"
            fi
        done
    fi

    # Final sanity check
    command -v cargo &>/dev/null || fail "cargo still not found after Rust install. Open a new terminal and re-run."
    ok "cargo $(cargo --version | awk '{print $2}')"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Build and install myconote-cli binary
# ─────────────────────────────────────────────────────────────────────────────
step "Step 4/6 — Building myconote-cli"

if [[ "${SKIP_BUILD:-0}" == "1" ]]; then
    warn "SKIP_BUILD=1, skipping binary build"
else
    export PATH="$HOME/.cargo/bin:$PATH"
    command -v cargo &>/dev/null || fail "cargo not found. Run without SKIP_RUST=1 or install Rust manually."

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    info "Building release binary (this takes ~1–2 min on first run)…"
    cd "$SCRIPT_DIR"

    # Show only high-signal lines; re-run verbosely if it fails
    if ! cargo build --release 2>&1 | grep -E "^(error|warning\[E|Compiling myconote|Finished)" ; then
        info "Retrying with full output to diagnose errors…"
        cargo build --release
    fi

    [[ -f "target/release/myconote-cli" ]] || fail "Binary not produced — build failed."
    ok "Build complete: ${SCRIPT_DIR}/target/release/myconote-cli"

    # Install to ~/.local/bin (no sudo required)
    mkdir -p "$INSTALL_DIR"
    cp "target/release/myconote-cli" "$INSTALL_DIR/myconote-cli"
    chmod +x "$INSTALL_DIR/myconote-cli"
    ok "Binary installed to: ${INSTALL_DIR}/myconote-cli"

    # Add ~/.local/bin to PATH in shell rc files if not already there
    for RC in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.zshrc"; do
        if [[ -f "$RC" ]] && ! grep -q "${INSTALL_DIR}" "$RC" 2>/dev/null; then
            echo "" >> "$RC"
            echo "# myconote-cli" >> "$RC"
            echo "export PATH=\"${INSTALL_DIR}:\$PATH\"" >> "$RC"
            ok "Added ${INSTALL_DIR} to PATH in $RC"
        fi
    done
    export PATH="${INSTALL_DIR}:$PATH"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — Install all 27 bioinformatics tools via myconote-cli install
# ─────────────────────────────────────────────────────────────────────────────
step "Step 5/6 — Installing bioinformatics tools"

if [[ "${SKIP_TOOLS:-0}" == "1" ]]; then
    warn "SKIP_TOOLS=1, skipping bioinformatics tool installation"
else
    MYCLI="${INSTALL_DIR}/myconote-cli"
    [[ -x "$MYCLI" ]] || MYCLI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/target/release/myconote-cli"
    [[ -x "$MYCLI" ]] || fail "myconote-cli binary not found — did Step 4 succeed?"

    # Prefer mamba/micromamba if available; otherwise fall back to conda
    if command -v mamba &>/dev/null || command -v micromamba &>/dev/null; then
        info "Using mamba for faster installs…"
        "$MYCLI" install --yes --mamba
    else
        info "Using conda for installs (mamba not found)…"
        "$MYCLI" install --yes
    fi

    echo ""
    ok "All automatable tools installed."
    echo ""
    warn "The following tools require a MANUAL download (licence restrictions):"
    echo -e "    ${RED}✗${RESET}  signalp6       Free with registration:"
    echo -e "               https://services.healthtech.dtu.dk/services/SignalP-6.0/"
    echo -e "    ${YELLOW}⚠${RESET}  gmes_petap.pl  Academic licence (Perl deps via conda):"
    echo -e "               http://topaz.gatech.edu/GeneMark/"
    echo -e "               After download: conda install -c bioconda perl-hash-merge"
    echo ""
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — Download annotation databases (~2.5 GB)
# ─────────────────────────────────────────────────────────────────────────────
step "Step 6/6 — Downloading annotation databases (~2.5 GB)"

if [[ "${SKIP_DATABASES:-0}" == "1" ]]; then
    warn "SKIP_DATABASES=1, skipping database downloads"
else
    mkdir -p "$DB_DIR"

    # ── Swiss-Prot + MMseqs2 index ───────────────────────────────────────────
    SPROT_FA="${DB_DIR}/uniprot_sprot.fasta.gz"
    SPROT_DB="${DB_DIR}/swissprot/swissprot"

    if [[ -f "${SPROT_DB}.index" ]]; then
        ok "Swiss-Prot MMseqs2 database already exists"
    else
        info "Downloading UniProt Swiss-Prot (~270 MB compressed)…"
        curl -fL \
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
            -o "$SPROT_FA" --progress-bar
        info "Building MMseqs2 database…"
        mkdir -p "${DB_DIR}/swissprot"
        TMP_SPROT=$(mktemp -d)
        mmseqs createdb "$SPROT_FA" "$SPROT_DB" --threads 4
        mmseqs createindex "$SPROT_DB" "$TMP_SPROT" --threads 4
        rm -rf "$TMP_SPROT"
        ok "Swiss-Prot MMseqs2 database ready"
    fi

    # ── Pfam-A ──────────────────────────────────────────────────────────────
    PFAM_HMM="${DB_DIR}/pfam/Pfam-A.hmm"

    if [[ -f "${PFAM_HMM}.h3i" ]]; then
        ok "Pfam-A already pressed"
    else
        info "Downloading Pfam-A (~500 MB compressed)…"
        mkdir -p "${DB_DIR}/pfam"
        curl -fL \
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" \
            -o "${PFAM_HMM}.gz" --progress-bar
        info "Decompressing and pressing Pfam-A…"
        gunzip -f "${PFAM_HMM}.gz"
        hmmpress "$PFAM_HMM"
        ok "Pfam-A ready"
    fi

    # ── BUSCO lineages ───────────────────────────────────────────────────────
    BUSCO_DB_DIR="${DB_DIR}/busco_lineages"
    mkdir -p "${BUSCO_DB_DIR}/lineages"

    download_busco_lineage() {
        local LIN="$1"
        local LINDIR="${BUSCO_DB_DIR}/lineages/${LIN}"
        if [[ -d "$LINDIR" ]]; then
            ok "BUSCO lineage already present: ${LIN}"; return
        fi
        info "Downloading BUSCO lineage: ${LIN}…"
        if busco --download "$LIN" --download_path "$BUSCO_DB_DIR" --quiet 2>/dev/null; then
            ok "BUSCO lineage ready: ${LIN}"
        else
            local URL="https://busco-data.ezlab.org/v5/data/lineages/${LIN}.tar.gz"
            local TGZ="${BUSCO_DB_DIR}/${LIN}.tar.gz"
            if curl -fsSL "$URL" -o "$TGZ" --progress-bar 2>/dev/null; then
                tar -xzf "$TGZ" -C "${BUSCO_DB_DIR}/lineages/"
                rm -f "$TGZ"
                ok "BUSCO lineage ready: ${LIN}"
            else
                warn "Could not pre-download '${LIN}' — BUSCO will fetch it on first use"
            fi
        fi
    }

    download_busco_lineage "fungi_odb10"
    download_busco_lineage "viridiplantae_odb10"
    download_busco_lineage "metazoa_odb10"
    download_busco_lineage "insecta_odb10"

    ok "All databases ready in ${DB_DIR}"

    # ── eggNOG-mapper (optional, large) ──────────────────────────────────────
    echo ""
    warn "eggNOG-mapper databases (~50 GB) are NOT downloaded by default."
    warn "Only needed if you plan to use --eggnog in annotate."
    warn "To download them later, run:"
    echo -e "      download_eggnog_data.py -y --data_dir ~/.eggnog_mapper/data"
    echo ""
fi

# ─────────────────────────────────────────────────────────────────────────────
# Done
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}╔═══════════════════════════════════════════════════════╗${RESET}"
echo -e "${BOLD}${GREEN}║   Installation complete!                              ║${RESET}"
echo -e "${BOLD}${GREEN}╚═══════════════════════════════════════════════════════╝${RESET}"
echo ""
echo -e "  ${BOLD}Verify everything is ready:${RESET}"
echo -e "    myconote-cli check"
echo ""
echo -e "  ${BOLD}Quick-start (full pipeline):${RESET}"
echo ""
echo -e "    # 1. Mask repeats"
echo -e "    myconote-cli mask genome.fa --engine repeatmodeler --threads 8"
echo ""
echo -e "    # 2. Predict genes"
echo -e "    myconote-cli predict genome_masked.fa --kingdom fungi \\"
echo -e "                 --locus-prefix MYFUN --train --threads 8"
echo ""
echo -e "    # 3. Annotate"
echo -e "    myconote-cli annotate predict_out/consensus.gff3 \\"
echo -e "                 --fasta genome_masked.fa --kingdom fungi \\"
echo -e "                 --cazyme --antismash --threads 8"
echo ""
echo -e "  ${BOLD}Open a new terminal${RESET} (or run ${BOLD}source ~/.zshrc${RESET}) for PATH changes to take effect."
echo ""
