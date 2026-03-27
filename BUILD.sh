#!/usr/bin/env bash
# Build helper for myconote-cli
# Run this from the project root: bash BUILD.sh

set -e

MIN_RUST_VERSION="1.74"   # driven by clap 4.5

echo "═══════════════════════════════════════════════════════"
echo "  myconote-cli build script"
echo "═══════════════════════════════════════════════════════"
echo

# ── Rust presence check ───────────────────────────────────────────────────────
if ! command -v cargo &>/dev/null; then
    echo "✗ cargo not found."
    echo ""
    echo "  Install Rust (no root required) with:"
    echo "    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    echo ""
    echo "  Then open a new terminal and re-run this script."
    exit 1
fi

# ── Rust version check ────────────────────────────────────────────────────────
RUSTC_VER="$(rustc --version 2>/dev/null | awk '{print $2}')"

version_gte() {
    local a b
    IFS='.' read -ra a <<< "$(echo "$1" | sed 's/[^0-9.]//g')"
    IFS='.' read -ra b <<< "$(echo "$2" | sed 's/[^0-9.]//g')"
    local i
    for (( i=0; i<${#b[@]}; i++ )); do
        local av="${a[$i]:-0}" bv="${b[$i]:-0}"
        (( 10#$av > 10#$bv )) && return 0
        (( 10#$av < 10#$bv )) && return 1
    done
    return 0
}

if ! version_gte "$RUSTC_VER" "$MIN_RUST_VERSION"; then
    echo "✗ Rust ${RUSTC_VER} is too old. Need ≥ ${MIN_RUST_VERSION}."
    echo ""
    if command -v rustup &>/dev/null; then
        echo "  Update with:  rustup update stable"
    else
        echo "  Install a newer version:"
        echo "    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    fi
    exit 1
fi

echo "✓ $(cargo --version)"
echo "✓ $(rustc --version)"
echo

# ── Development build ─────────────────────────────────────────────────────────
echo "─── cargo build ────────────────────────────────────────"
cargo build 2>&1
echo

echo "─── Build complete ─────────────────────────────────────"
echo "  Binary: ./target/debug/myconote-cli"
echo ""
echo "  Quick test:"
echo "    ./target/debug/myconote-cli --help"
echo "    ./target/debug/myconote-cli check"
echo ""
echo "  Release build (optimised, ~1 min):"
echo "    cargo build --release"
echo "  Binary: ./target/release/myconote-cli"
echo
