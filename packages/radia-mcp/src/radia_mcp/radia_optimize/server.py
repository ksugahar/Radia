"""Startup entry point for radia-mcp.optimize."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.optimize")


if __name__ == "__main__":
    main()
