"""Synchronize MCP's standalone numerical gate from its sole source of truth.

Run after editing src/radia/_residual_gate.py. --check is read-only and fails
if the independently distributed MCP copy is stale. No runtime fallback is used.
"""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'src/radia/_residual_gate.py'
TARGET = ROOT / 'packages/radia-mcp/src/radia_mcp/radia_ngsolve/_vendor/residual_gate.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    source = SOURCE.read_bytes()
    if args.check:
        if not TARGET.exists() or TARGET.read_bytes() != source:
            parser.exit(1, 'MCP residual gate is stale; run tools/sync_mcp_residual_gate.py\n')
    else:
        TARGET.parent.mkdir(parents=True, exist_ok=True)
        TARGET.write_bytes(source)
    print('MCP residual gate matches its canonical source')


if __name__ == '__main__':
    main()
