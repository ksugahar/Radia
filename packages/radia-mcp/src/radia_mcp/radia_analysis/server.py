"""Startup entry point for radia-analysis."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-analysis")


if __name__ == "__main__":
    main()
