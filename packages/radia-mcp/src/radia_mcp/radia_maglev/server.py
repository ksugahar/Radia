"""Startup entry point for radia-mcp-maglev."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp-maglev")


if __name__ == "__main__":
    main()
