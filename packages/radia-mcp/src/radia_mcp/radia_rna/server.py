"""Startup entry point for radia-mcp.rna."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.rna")


if __name__ == "__main__":
    main()
