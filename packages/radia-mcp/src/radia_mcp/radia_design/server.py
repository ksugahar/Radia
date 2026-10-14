"""Startup entry point for radia-design."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("radia-design")


if __name__ == "__main__":
    main()
