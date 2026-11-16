"""Startup entry point for radia-mcp-acoustic."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp-acoustic")


if __name__ == "__main__":
    main()
