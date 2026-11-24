"""Startup entry point for radia-mcp.postprocess."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.postprocess")


if __name__ == "__main__":
    main()
