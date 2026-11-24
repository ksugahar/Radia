"""Startup entry point for radia-mcp.learning."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.learning")


if __name__ == "__main__":
    main()
