"""Startup entry point for radia-mcp-ih."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp-ih")


if __name__ == "__main__":
    main()
