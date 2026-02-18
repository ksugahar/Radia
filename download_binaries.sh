#!/bin/bash
# Download pre-built binary files from GitHub Releases
# Usage: ./download_binaries.sh

set -e

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

echo "Downloading pre-built binaries from GitHub Releases..."

# Download .pyd modules
echo "  -> src/radia/*.pyd"
gh release download binaries -D src/radia/ -p "*.pyd" --clobber

# Download fmm3d library
echo "  -> src/ext/fmm3d/lib/fmm3d.lib"
mkdir -p src/ext/fmm3d/lib
gh release download binaries -D src/ext/fmm3d/lib/ -p "fmm3d.lib" --clobber

echo "Done. Downloaded files:"
ls -la src/radia/*.pyd src/ext/fmm3d/lib/fmm3d.lib 2>/dev/null || true
