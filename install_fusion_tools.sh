#!/usr/bin/env bash
#
# install_fusion_tools.sh — install the external CLI tools used by the Step-3
# fusion-design stage (fuse.py) on the Linux server.
#
# Tools and how they install:
#   AGGRESCAN3D   Gate 2 (aggregation)  — free, pip-installable → fully automated here
#   NetMHCpan 4.2 Gate 1 (MHC-I)        — DTU academic licence; you download the tarball,
#                                         this script extracts + configures it
#   NetCTLpan 1.1 Gate 1 (DEPRECATED)   — discontinued by DTU (superseded by NetMHCpan);
#                                         configured only if a legacy tarball is in vendor/
#   FoldX 5       structure phase       — academic licence tarball; wired to FOLDX_PATH
#   CamSol/UbPred Gate2 alt / Gate3 opt — web-server only, no CLI (reported, not installed)
#
# Licence-gated tools cannot be auto-downloaded. Register, download, and drop the
# tarballs into the vendor dir, then run this script. Anything missing is skipped
# with a note, so you can run it repeatedly as you obtain each tool.
#
# Usage:
#   conda activate ecr
#   bash install_fusion_tools.sh                 # uses ./vendor and ~/opt/ecr_tools
#   VENDOR_DIR=/path/to/tarballs bash install_fusion_tools.sh
#
# Download pages (academic licence — free):
#   NetMHCpan 4.2  https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/
#                  (NetCTLpan 1.1 is discontinued by DTU; NetMHCpan supersedes it —
#                   it covers MHC-I binding but not proteasomal cleavage/TAP)
#   FoldX 5        https://foldxsuite.crg.eu/
#
set -euo pipefail

# ── Configuration (override via env vars) ───────────────────────────────────
VENDOR_DIR="${VENDOR_DIR:-$(pwd)/vendor}"          # where you drop downloaded tarballs
INSTALL_DIR="${INSTALL_DIR:-$HOME/opt/ecr_tools}"  # where tools get installed
BIN_DIR="${BIN_DIR:-$INSTALL_DIR/bin}"             # wrappers/symlinks go here (add to PATH)
TMP_DIR="${INSTALL_DIR}/tmp"

mkdir -p "$VENDOR_DIR" "$INSTALL_DIR" "$BIN_DIR" "$TMP_DIR"

msg()  { printf '\033[1;34m[install]\033[0m %s\n' "$*"; }
ok()   { printf '\033[1;32m[ ok ]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[skip]\033[0m %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }

# First tarball matching a glob in VENDOR_DIR, or empty string.
find_tarball() {
    local f
    f=$(ls -1 "$VENDOR_DIR"/$1 2>/dev/null | head -n1 || true)
    printf '%s' "$f"
}

# ── Pre-flight checks ───────────────────────────────────────────────────────
msg "vendor dir : $VENDOR_DIR  (drop licence-gated tarballs here)"
msg "install dir: $INSTALL_DIR"
msg "bin dir    : $BIN_DIR  (add this to PATH)"
echo

if [[ "${CONDA_DEFAULT_ENV:-}" != "ecr" ]]; then
    warn "conda env 'ecr' is not active (current: '${CONDA_DEFAULT_ENV:-none}')."
    warn "AGGRESCAN3D will install into whatever Python is on PATH. Run 'conda activate ecr' first if that's not intended."
    echo
fi

# tcsh is required by the DTU wrapper scripts (netMHCpan / netCTLpan are tcsh).
if ! have tcsh; then
    warn "tcsh not found — NetMHCpan/NetCTLpan wrappers need it."
    if have conda; then
        msg "installing tcsh via conda..."
        conda install -y -c conda-forge tcsh >/dev/null && ok "tcsh installed"
    else
        warn "install tcsh with your package manager (e.g. 'sudo apt-get install tcsh')."
    fi
    echo
fi

# ════════════════════════════════════════════════════════════════════════════
# AGGRESCAN3D (Gate 2) — pip, fully automated
# ════════════════════════════════════════════════════════════════════════════
msg "AGGRESCAN3D (Gate 2 aggregation)"
if have aggrescan3d; then
    ok "aggrescan3d already on PATH ($(command -v aggrescan3d))"
else
    # freesasa is the SASA backend A3D needs; numpy/scipy come as deps.
    pip install --quiet aggrescan3d freesasa && ok "aggrescan3d installed via pip" \
        || warn "pip install aggrescan3d failed — check pip/conda env and retry."
fi
echo

# ════════════════════════════════════════════════════════════════════════════
# NetMHCpan 4.2 (Gate 1) — licence-gated tarball
# ════════════════════════════════════════════════════════════════════════════
install_netmhcpan() {
    local tarball pkgdir
    tarball=$(find_tarball 'netMHCpan-*.Linux.tar.gz')
    [[ -z "$tarball" ]] && tarball=$(find_tarball 'netMHCpan-*.tar.gz')
    if [[ -z "$tarball" ]]; then
        warn "NetMHCpan: no tarball in $VENDOR_DIR (expected netMHCpan-4.2c.Linux.tar.gz)."
        warn "  Get it: https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/"
        return
    fi
    msg "NetMHCpan: extracting $(basename "$tarball")"
    tar -xzf "$tarball" -C "$INSTALL_DIR"
    pkgdir=$(ls -d "$INSTALL_DIR"/netMHCpan-* 2>/dev/null | head -n1)

    # Data files may ship separately; fetch only if the package didn't include them
    # (the 4.2c static build often bundles data, in which case this is skipped).
    if [[ ! -d "$pkgdir/data" ]]; then
        msg "NetMHCpan: data/ missing — fetching data files..."
        if have wget; then
            wget -q "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/data.tar.gz" \
                -O "$pkgdir/data.tar.gz" && tar -xzf "$pkgdir/data.tar.gz" -C "$pkgdir" \
                && rm -f "$pkgdir/data.tar.gz" && ok "data installed" \
                || warn "data download failed — fetch data.tar.gz manually into $pkgdir and untar."
        else
            warn "wget missing — download data.tar.gz manually into $pkgdir and untar."
        fi
    fi

    # Patch the tcsh wrapper: NMHOME must point at the real install dir, TMPDIR writable.
    sed -i.bak -E "s|^setenv[[:space:]]+NMHOME.*|setenv NMHOME $pkgdir|" "$pkgdir/netMHCpan"
    sed -i     -E "s|^setenv[[:space:]]+TMPDIR.*|setenv TMPDIR $TMP_DIR|" "$pkgdir/netMHCpan"
    chmod +x "$pkgdir/netMHCpan"
    ln -sf "$pkgdir/netMHCpan" "$BIN_DIR/netMHCpan"
    ok "NetMHCpan wired → $BIN_DIR/netMHCpan"
}
install_netmhcpan
echo

# ════════════════════════════════════════════════════════════════════════════
# NetCTLpan 1.1 (Gate 1) — DEPRECATED: discontinued by DTU, superseded by
# NetMHCpan. We do NOT point new users at it; only configure a legacy tarball if
# one is already present (for pipelines that still depend on it).
# ════════════════════════════════════════════════════════════════════════════
install_netctlpan() {
    local tarball pkgdir
    tarball=$(find_tarball 'netCTLpan-*.Linux.tar.gz')
    [[ -z "$tarball" ]] && tarball=$(find_tarball 'netCTLpan-*.tar.gz')
    if [[ -z "$tarball" ]]; then
        msg "NetCTLpan: discontinued by DTU — skipping (use NetMHCpan). Drop a legacy"
        msg "  netCTLpan-*.tar.gz in $VENDOR_DIR only if an existing pipeline needs it."
        return
    fi
    warn "NetCTLpan is DEPRECATED (discontinued by DTU) — configuring legacy tarball anyway."
    msg "NetCTLpan: extracting $(basename "$tarball")"
    tar -xzf "$tarball" -C "$INSTALL_DIR"
    pkgdir=$(ls -d "$INSTALL_DIR"/netCTLpan-* 2>/dev/null | head -n1)

    # NetCTLpan's wrapper sets NETCTLpan (install root) and TMPDIR.
    sed -i.bak -E "s|^setenv[[:space:]]+NETCTLpan.*|setenv NETCTLpan $pkgdir|" "$pkgdir/netCTLpan"
    sed -i     -E "s|^setenv[[:space:]]+TMPDIR.*|setenv TMPDIR $TMP_DIR|" "$pkgdir/netCTLpan"
    chmod +x "$pkgdir/netCTLpan"
    ln -sf "$pkgdir/netCTLpan" "$BIN_DIR/netCTLpan"
    ok "NetCTLpan wired → $BIN_DIR/netCTLpan"
}
install_netctlpan
echo

# ════════════════════════════════════════════════════════════════════════════
# FoldX 5 (structure phase) — licence-gated tarball
# ════════════════════════════════════════════════════════════════════════════
install_foldx() {
    local tarball pkgdir binpath
    tarball=$(find_tarball 'foldx*5*.tar.gz')
    [[ -z "$tarball" ]] && tarball=$(find_tarball 'foldx*.tar.gz')
    if [[ -z "$tarball" ]]; then
        warn "FoldX: no tarball in $VENDOR_DIR. Get it: https://foldxsuite.crg.eu/"
        return
    fi
    msg "FoldX: extracting $(basename "$tarball")"
    pkgdir="$INSTALL_DIR/foldx5"
    mkdir -p "$pkgdir"
    tar -xzf "$tarball" -C "$pkgdir"
    # The executable is named foldx / foldx_<date>; symlink to a stable name.
    binpath=$(find "$pkgdir" -maxdepth 2 -type f -name 'foldx*' -perm -u+x 2>/dev/null | head -n1)
    if [[ -z "$binpath" ]]; then
        binpath=$(find "$pkgdir" -maxdepth 2 -type f -name 'foldx*' 2>/dev/null | head -n1)
        [[ -n "$binpath" ]] && chmod +x "$binpath"
    fi
    if [[ -z "$binpath" ]]; then
        warn "FoldX: extracted but no 'foldx*' binary found under $pkgdir — set FOLDX_PATH manually."
        return
    fi
    ln -sf "$binpath" "$BIN_DIR/foldx"
    ok "FoldX wired → $BIN_DIR/foldx  (FOLDX_PATH below)"
}
install_foldx
echo

# ════════════════════════════════════════════════════════════════════════════
# Web-only tools — reported, not installable
# ════════════════════════════════════════════════════════════════════════════
warn "CamSol — web server only (no CLI). Use AGGRESCAN3D for Gate 2, or set the"
warn "         camsol 'api' backend in config.yaml if you wrap the service yourself."
warn "UbPred — web server only and OPTIONAL. Gate 3 (N-end rule + degron scan)"
warn "         runs without it; leave fusion.tools.ubpred.backend: disabled."
echo

# ── Summary ─────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════════"
ok   "Done. Next steps:"
echo
echo "1. Put the tool bin dir on PATH (add to ~/.bashrc to persist):"
echo "     export PATH=\"$BIN_DIR:\$PATH\""
if [[ -L "$BIN_DIR/foldx" ]]; then
echo
echo "2. Point the refine/fusion stages at FoldX:"
echo "     export FOLDX_PATH=\"$BIN_DIR/foldx\""
fi
echo
echo "3. Enable the tools you installed in config.yaml (fusion.tools.<tool>.backend: local):"
echo "     - one MHC-I tool   → netmhcpan (NetMHCpan 4.2; NetCTLpan is deprecated)"
echo "     - one aggregation  → aggrescan3d"
echo "   Leave camsol / ubpred as 'disabled'."
echo
echo "4. Verify:"
for t in aggrescan3d netMHCpan netCTLpan foldx; do
    if PATH="$BIN_DIR:$PATH" have "$t"; then
        printf "     \033[1;32m✓\033[0m %s\n" "$t"
    else
        printf "     \033[1;33m–\033[0m %s (not installed)\n" "$t"
    fi
done
echo "════════════════════════════════════════════════════════════════════════"
