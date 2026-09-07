"""Read-only development editable-source verification.

Defaults to the canonical LAB checkout. Use --mcp-source explicitly for an
approved clean MCP runtime. Never infer the expectation from the installed
package: that would accept drift as its own source of truth.
"""
from __future__ import annotations

import argparse
import release_quad


def expected_packages(mcp_source=None):
    packages = release_quad._canonical_lab_editable_packages()
    if mcp_source:
        packages = [(name, mcp_source if name == "radia-mcp" else path)
                    for name, path in packages]
    return packages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcp-source", help="Approved editable MCP project root")
    args = parser.parse_args(argv)
    return 4 if release_quad._verify_lab_editable(
        expected_packages(args.mcp_source)) else 0


if __name__ == "__main__":
    raise SystemExit(main())
