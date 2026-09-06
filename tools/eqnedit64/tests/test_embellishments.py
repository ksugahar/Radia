"""Every decoration must draw as itself, and none may draw as another.

The renderer chose its accent with an if-chain that fell through to the hat,
and it listed five of nineteen embellishments. So the Decoration palette's
prime, double prime, triple prime, double dot, triple dot, cancel, frown and
smile all drew a hat -- while the LaTeX and the Office MathML they produced
were correct. A wrong picture beside a right paste is the worst shape for a
defect: nothing downstream disagrees, so only a person looking at the screen
can catch it, and one did, after 3.0.15 shipped.

The compiler now refuses an unmapped embellishment (a switch with no default
under /W4 /WX /w14062), which stops the class at its source. These checks
cover what a compiler cannot: that the marks are actually DIFFERENT, and that
a prime lands beside the letter rather than on top of it.

Run:  python tests\\test_embellishments.py
"""
from __future__ import annotations

import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import Equation, math_font_loaded  # noqa: E402

# The Decoration palette, in the order palettes.cpp declares it.
DECORATIONS = (
    "hat", "tilde", "bar", "vec", "dot", "ddot", "dddot",
    "prime", "dprime", "tprime", "strike", "frown", "smile",
)


def render(kind: str | None) -> str:
    """Put one decoration around the same base and return its SVG."""
    equation = Equation()
    if kind is not None:
        equation.insert_template(kind)
    equation.insert_text("a")
    return equation.svg()


def svg_width(svg: str) -> float:
    match = re.search(r'width="([0-9.]+)', svg)
    return float(match.group(1)) if match else -1.0


def main() -> int:
    failures: list[str] = []
    if not math_font_loaded():
        print("FAIL  embedded math font unavailable")
        return 1

    drawn = {kind: render(kind) for kind in DECORATIONS}
    for kind, svg in drawn.items():
        if not svg.strip():
            failures.append(f"{kind} rendered nothing")

    # The fall-through made eight of these identical to the hat.
    seen: dict[str, str] = {}
    for kind, svg in drawn.items():
        if svg in seen:
            failures.append(f"{kind} renders exactly like {seen[svg]}")
        else:
            seen[svg] = kind

    # A prime is a superscript suffix: it advances the line where a centred
    # accent does not. Comparing against the bare base keeps this independent
    # of the font's own advance widths.
    plain = svg_width(render(None))
    for kind in ("prime", "dprime", "tprime"):
        got = svg_width(drawn[kind])
        if got <= plain:
            failures.append(
                f"{kind} did not advance past the base ({got} <= {plain}); "
                "it was drawn as an accent instead of a suffix")

    # More primes must take more room than fewer.
    if not (svg_width(drawn["prime"]) < svg_width(drawn["dprime"])
            < svg_width(drawn["tprime"])):
        failures.append(
            "prime widths do not increase with the number of primes: "
            f"{svg_width(drawn['prime'])}, {svg_width(drawn['dprime'])}, "
            f"{svg_width(drawn['tprime'])}")

    # A template the palette offers must be one the model accepts.
    known = set(Equation.templates())
    missing = [kind for kind in DECORATIONS if kind not in known]
    if missing:
        failures.append(f"palette offers templates the model rejects: {missing}")

    if failures:
        print(f"FAIL  {len(failures)}")
        for failure in failures:
            print("  " + failure)
        return 1
    print(f"ok    {len(DECORATIONS)} decorations draw distinctly; "
          "primes are suffixes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
