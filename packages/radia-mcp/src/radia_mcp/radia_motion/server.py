"""Startup entry point for radia-motion."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-motion")


if __name__ == "__main__":
    main()
