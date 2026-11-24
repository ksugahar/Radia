"""Startup entry point for radia-mcp.pcb."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-mcp.pcb")


if __name__ == "__main__":
    main()
