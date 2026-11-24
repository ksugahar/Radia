"""Startup entry point for radia-mcp.motor."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.motor")


if __name__ == "__main__":
    main()
