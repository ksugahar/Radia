"""Startup entry point for document-ops."""

from ..capability_packs.server import main as run_pack


def main():
    run_pack("document-ops")


if __name__ == "__main__":
    main()
