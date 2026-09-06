"""Every palette cell, exercised: does it draw as itself and survive a save?

The prime cell drew a hat for a whole release because the renderer picked its
mark with an if-chain that fell through to the hat, and no test compared one
decoration against another. Reasoning about which families might share that
shape found four more defects, but reasoning only reaches as far as someone
thought to look. This sweep does not reason: it presses every cell the
palettes offer and asks two questions that do not need a human eye.

  1. Do two different commands produce the same picture? Two cells that draw
     identically are either a duplicate entry or a fall-through, and both are
     defects. This is what would have caught the prime.
  2. Does the TeX the editor emits load back as the same equation? A cell
     whose markup is dropped or rewritten on save loses the user's work
     quietly; the round trip is where that shows.

Distinct rendering is checked WITHIN a palette rather than globally: separate
palettes legitimately offer the same symbol from different categories.

Run:  python tests\\test_palette_sweep.py
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import Equation, palettes, math_font_loaded  # noqa: E402


def build(command: str) -> Equation | None:
    """Apply one palette command to a fresh equation, exactly as the window
    does in eqnedt64_app.cpp: templates, symbols, raw LaTeX and typeface
    changes each have their own entry point, and using the wrong one would
    test something the user never does."""
    equation = Equation()
    if command.startswith("template."):
        if not equation.insert_template(command[len("template."):]):
            return None
    elif command.startswith("symbol."):
        if not equation.insert_symbol(command[len("symbol."):]):
            return None
    elif command.startswith("latex."):
        if not equation.insert_latex(command[len("latex."):]):
            return None
    elif command.startswith("style."):
        equation.insert_text("x")
        equation.select_all()
        if not equation.restyle_selection(command[len("style."):]):
            return None
    else:
        return None
    return equation


def main() -> int:
    if not math_font_loaded():
        print("FAIL  embedded math font unavailable")
        return 1

    failures: list[str] = []
    cells = 0
    unusable: list[str] = []

    for title, _face, _columns, items in palettes():
        drawn: dict[str, str] = {}
        for command, face, _label in items:
            cells += 1
            equation = build(command)
            if equation is None:
                unusable.append(f"{title}/{face}: {command} was rejected")
                continue

            svg = equation.svg()
            if not svg.strip():
                failures.append(f"{title}/{face}: {command} rendered nothing")
                continue
            if svg in drawn and drawn[svg] != command:
                failures.append(
                    f"{title}: {command} draws exactly like {drawn[svg]}")
            else:
                drawn[svg] = command

            # What the editor saves must be what it loads.
            tex = equation.latex()
            again = Equation()
            if not again.load_latex(tex):
                failures.append(
                    f"{title}/{face}: {command} emitted TeX that will not "
                    f"load back: {tex!r}")
            elif again.latex() != tex:
                failures.append(
                    f"{title}/{face}: {command} is not a fixed point: "
                    f"{tex!r} -> {again.latex()!r}")

    if unusable:
        failures.extend(unusable)

    if failures:
        print(f"FAIL  {len(failures)} of {cells} palette cells")
        for failure in failures[:40]:
            print("  " + failure)
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
        return 1
    print(f"ok    {cells} palette cells draw distinctly within their palette "
          "and round-trip through TeX")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
