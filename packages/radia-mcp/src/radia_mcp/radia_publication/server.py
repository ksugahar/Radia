"""Startup entry point for radia-publication."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-publication")


if __name__ == "__main__":
    main()
