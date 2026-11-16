"""Startup entry point for radia-mcp-electromagnet."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp-electromagnet")


if __name__ == "__main__":
    main()
