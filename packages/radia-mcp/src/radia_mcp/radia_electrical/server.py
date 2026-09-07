"""Startup entry point for radia-electrical."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-electrical")


if __name__ == "__main__":
    main()
