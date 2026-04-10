#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  myconote-cli — quick binary installer
#
#  Downloads a pre-built binary for your platform. No Rust required.
#  External bioinformatics tools (Augustus, HMMER, etc.) must be installed
#  separately — run `myconote-cli install` after this script.
#
#  Usage:
#    curl -fsSL https://raw.githubusercontent.com/K-nie/myconote-cli/main/quick-install.sh | bash
#
#  Or:
#    bash quick-install.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

REPO="K-nie/myconote-cli"
INSTALL_DIR="${HOME}/.local/bin"

RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

ok()   { echo -e "${GREEN}  ✓${RESET}  $*"; }
info() { echo -e "${CYAN}  →${RESET}  $*"; }
fail() { echo -e "${RED}  ✗${RESET}  $*"; exit 1; }

echo ""
echo -e "${BOLD}myconote-cli — quick installer${RESET}"
echo ""

# ── Detect platform ──────────────────────────────────────────────────────────
OS="$(uname -s)"
ARCH="$(uname -m)"

case "${OS}-${ARCH}" in
    Linux-x86_64)       ARTIFACT="myconote-cli-linux-amd64" ;;
    Linux-aarch64)      ARTIFACT="myconote-cli-linux-arm64" ;;
    Darwin-arm64)       ARTIFACT="myconote-cli-macos-arm64" ;;
    Darwin-x86_64)      ARTIFACT="myconote-cli-macos-amd64" ;;
    *) fail "Unsupported platform: ${OS}-${ARCH}. Build from source: cargo build --release" ;;
esac

info "Platform: ${OS} ${ARCH} → ${ARTIFACT}"

# ── Get latest release URL ───────────────────────────────────────────────────
info "Fetching latest release from GitHub..."

RELEASE_URL=$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" \
    | grep "browser_download_url.*${ARTIFACT}" \
    | head -1 \
    | cut -d '"' -f 4)

if [[ -z "${RELEASE_URL:-}" ]]; then
    # No release yet — fall back to building from source
    echo ""
    echo -e "${CYAN}  No pre-built release found.${RESET}"
    echo -e "  This is expected for v0.1.0 — the first release hasn't been tagged yet."
    echo ""
    echo -e "  ${BOLD}Install from source instead:${RESET}"
    echo ""
    echo "    git clone https://github.com/${REPO}.git"
    echo "    cd myconote-cli"
    echo "    cargo build --release"
    echo "    cp target/release/myconote-cli ~/.local/bin/"
    echo ""
    echo -e "  Or use the full installer: ${BOLD}bash install.sh${RESET}"
    echo ""
    exit 0
fi

# ── Download and install ─────────────────────────────────────────────────────
info "Downloading ${ARTIFACT}..."
TMP_DIR=$(mktemp -d)
curl -fsSL "${RELEASE_URL}" -o "${TMP_DIR}/${ARTIFACT}.tar.gz"

info "Extracting..."
tar xzf "${TMP_DIR}/${ARTIFACT}.tar.gz" -C "${TMP_DIR}"

mkdir -p "${INSTALL_DIR}"
cp "${TMP_DIR}/myconote-cli" "${INSTALL_DIR}/myconote-cli"
chmod +x "${INSTALL_DIR}/myconote-cli"
rm -rf "${TMP_DIR}"

ok "Installed to ${INSTALL_DIR}/myconote-cli"

# ── Add to PATH if needed ────────────────────────────────────────────────────
if ! echo "$PATH" | grep -q "${INSTALL_DIR}"; then
    SHELL_RC=""
    if [[ -f "$HOME/.zshrc" ]]; then
        SHELL_RC="$HOME/.zshrc"
    elif [[ -f "$HOME/.bashrc" ]]; then
        SHELL_RC="$HOME/.bashrc"
    elif [[ -f "$HOME/.bash_profile" ]]; then
        SHELL_RC="$HOME/.bash_profile"
    fi

    if [[ -n "$SHELL_RC" ]] && ! grep -q "${INSTALL_DIR}" "$SHELL_RC" 2>/dev/null; then
        echo "" >> "$SHELL_RC"
        echo "# myconote-cli" >> "$SHELL_RC"
        echo "export PATH=\"${INSTALL_DIR}:\$PATH\"" >> "$SHELL_RC"
        ok "Added ${INSTALL_DIR} to PATH in ${SHELL_RC}"
    fi
    export PATH="${INSTALL_DIR}:$PATH"
fi

# ── Verify ────────────────────────────────────────────────────────────────────
echo ""
"${INSTALL_DIR}/myconote-cli" --version 2>/dev/null | head -10
echo ""

ok "myconote-cli is ready!"
echo ""
echo -e "  ${BOLD}Next steps:${RESET}"
echo ""
echo "    # Install external bioinformatics tools"
echo "    myconote-cli install --yes"
echo ""
echo "    # Download annotation databases"
echo "    myconote-cli setup"
echo ""
echo "    # Learn how to use it (interactive tutorial)"
echo "    myconote-cli learn"
echo ""
echo "    # Or jump straight in"
echo "    myconote-cli predict genome.fa --kingdom fungi"
echo ""

if ! echo "$PATH" | grep -q "${INSTALL_DIR}"; then
    echo -e "  ${CYAN}Note:${RESET} Open a new terminal for PATH changes to take effect."
    echo ""
fi
