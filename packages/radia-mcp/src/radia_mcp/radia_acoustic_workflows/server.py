"""Startup entry point for radia-acoustic-workflows."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-acoustic-workflows")


if __name__ == "__main__":
    main()
