"""Startup entry point for radia-mcp-analysis."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp-analysis")


if __name__ == "__main__":
    main()
