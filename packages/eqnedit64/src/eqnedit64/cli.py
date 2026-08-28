"""Console launcher for the bundled standalone Eqnedit64 application."""
from __future__ import annotations

import subprocess
import sys

from .api import backend_path


def main() -> int:
    executable = str(backend_path())
    arguments = sys.argv[1:]
    if not arguments:
        subprocess.Popen([executable])
        return 0
    return subprocess.run([executable, *arguments], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
