"""Post-install contract for an Eqnedit64 wheel."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile

import eqnedit64


assert eqnedit64.__version__ == "3.0.13"
assert eqnedit64.backend_path().is_file()
assert eqnedit64.web_asset().is_file()
assert eqnedit64.web_asset("equation-editor.fragment.html").is_file()
assert (eqnedit64.web_asset().parents[1] / "licenses" / "GUST-FONT-LICENSE.txt").is_file()
assert eqnedit64.tex_normalize(r"\frac{x}{y}") == r"\frac{x}{y}"
assert "<mfrac>" in eqnedit64.tex_to_mathml(r"\frac{x}{y}")
assert "<svg" in eqnedit64.tex_to_svg(r"E=mc^2")

equation = eqnedit64.Equation()
assert equation.load_latex(r"E=mc^2H^2")
assert equation.backspace()
assert equation.latex() == r"E = mc^{2}H^{}"
assert equation.backspace()
assert equation.latex() == r"E = mc^{2}H"

with tempfile.TemporaryDirectory() as directory:
    rendered = eqnedit64.render_equation(
        r"\frac{x}{y}", Path(directory) / "equation.png"
    )
    assert rendered.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")

help_result = subprocess.run(
    [sys.executable, "-m", "eqnedit64.cli", "--help"],
    capture_output=True,
    check=False,
    text=True,
)
assert help_result.returncode == 0
assert "eqnedit64 INPUT OUTPUT" in help_result.stdout
assert "--copy-tex-file" not in help_result.stdout

print("PASS: installed eqnedit64 wheel exposes core, backend, and Web assets")
