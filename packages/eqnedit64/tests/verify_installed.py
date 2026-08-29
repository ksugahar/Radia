"""Post-install contract for an Eqnedit64 wheel."""
from __future__ import annotations

import eqnedit64


assert eqnedit64.__version__ == "3.0.9"
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

print("PASS: installed eqnedit64 wheel exposes core, backend, and Web assets")
