"""Console launcher for the bundled standalone Eqnedit64 application."""
from __future__ import annotations

import subprocess
import sys

from .api import backend_path


_HELP = """\
Eqnedit64 - TeX equation editor and converter

Usage:
  eqnedit64
  eqnedit64 equation.tex
  eqnedit64 INPUT OUTPUT

INPUT:
  equation.tex   UTF-8 TeX file
  clipboard      TeX currently on the clipboard

OUTPUT:
  office         editable PowerPoint / Word clipboard data
  slides         Google Slides PNG + HTML clipboard data
  png            PNG clipboard data
  equation.png   PNG image file
  equation.emf   EMF image file

Examples:
  eqnedit64 equation.tex office
  eqnedit64 equation.tex equation.png
  eqnedit64 clipboard png
"""


def main() -> int:
    arguments = sys.argv[1:]
    if arguments in (["--help"], ["-h"], ["/?"]):
        print(_HELP, end="")
        return 0
    executable = str(backend_path())
    if not arguments:
        subprocess.Popen([executable])
        return 0
    return subprocess.run([executable, *arguments], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
