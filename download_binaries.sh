#!/bin/bash
# Download pre-built binary files from GitHub Releases -- gh-CLI-FREE.
#
# Per the "No GitHub CLI" policy (CLAUDE.md), this uses
# tools/download_release_asset.py (raw GitHub REST API via urllib), NOT
# `gh release download`.  An optional $GH_TOKEN / $GITHUB_TOKEN lifts the
# anonymous 60 req/hr API limit but is not required for this public repo.
#
# Usage: ./download_binaries.sh
set -e

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

echo "Downloading pre-built binaries from GitHub Releases (gh-free)..."

# .pyd modules (multiple) -- pattern replaces `gh release download -p "*.pyd"`
echo "  -> src/radia/*.pyd"
python tools/download_release_asset.py --repo ksugahar/Radia --tag binaries \
    --pattern "*.pyd" --dest src/radia

# Standalone C ABI libraries (for Simulink/FMI/native consumers).
for pattern in "*.dll" "*.so" "*.dylib"; do
    echo "  -> src/radia/$pattern (optional by platform)"
    python tools/download_release_asset.py --repo ksugahar/Radia --tag binaries \
        --pattern "$pattern" --dest src/radia \
        || echo "  ($pattern not present for this platform)"
done

# fmm3d.lib is legacy (ExaFMM-t was removed); tolerate its absence.
echo "  -> src/ext/fmm3d/lib/fmm3d.lib (optional)"
python tools/download_release_asset.py --repo ksugahar/Radia --tag binaries \
    --name fmm3d.lib --dest src/ext/fmm3d/lib --min-size 1000 \
    || echo "  (fmm3d.lib not in release -- ok, FMM was removed)"

echo "Done. Downloaded files:"
ls -la src/radia/*.pyd src/radia/*.dll src/radia/*.so src/radia/*.dylib \
    src/ext/fmm3d/lib/fmm3d.lib 2>/dev/null || true
