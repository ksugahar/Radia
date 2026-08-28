"""Deterministic parser/normalizer/renderer fuzzing for raw TeX input.

The structural-operation fuzzer starts from a valid Equation tree.  This test
attacks the other boundary: arbitrary TeX pasted from papers, browsers, and
half-edited source.  Every generated input must parse without a crash, reach a
one-pass normalization fixed point, produce finite geometry, and emit valid
SVG XML.  Seeds are fixed so a failure is exactly reproducible.
"""
from __future__ import annotations

import math
import os
import random
import re
import sys
import xml.etree.ElementTree as ET

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import (  # noqa: E402
    Equation,
    tex_normalize,
    tex_to_mathml,
    tex_to_svg,
)


ATOMS = (
    "a", "x", "0", "12", "+", "-", "=", "(", ")", "[", "]", ",",
    "π", "θ", "速度", "Δ", "∂", "∞",
    r"\alpha", r"\Gamma", r"\partial", r"\nabla", r"\infty",
    r"\leq", r"\neq", r"\rightarrow", r"\cdot", r"\times",
    r"\operatorname{curl}", r"\text{speed}", r"\mathbf{x}",
    # A named glyph whose command must survive verbatim, an operator name
    # that must stay upright, and a variant Greek letter that used to share
    # its code point with the plain form.
    r"\ldots", r"\rangle", r"\parallel", r"\div", r"\lim", r"\varphi",
)

ENVIRONMENTS = (
    r"\begin{aligned}a&=b\\c&=d\end{aligned}",
    r"\begin{matrix}a&b\\c&d\end{matrix}",
    r"\begin{bmatrix}1&0\\0&1\end{bmatrix}",
    r"\begin{cases}x&x>0\\-x&x\leq0\end{cases}",
    # Shapes whose parse rules are special: a bra-ket, a fence whose sides
    # differ, and one with no visible delimiter at all.  Leaving them out is
    # how the \rangle round-trip defect reached a release with a green suite.
    r"\left\langle a\middle|b\right\rangle",
    r"\left(a,b\right]",
    r"\left.\frac{x}{y}\right.",
    r"\left\{\begin{matrix}a&b\end{matrix}\right.",
)

_VIEWBOX = re.compile(r'viewBox="0 0 ([^ ]+) ([^"]+)"')


def expression(rng: random.Random, depth: int = 0) -> str:
    if depth >= 5 or rng.randrange(100) < 34:
        return rng.choice(ATOMS)
    choice = rng.randrange(12)
    if choice == 0:
        return rf"\frac{{{expression(rng, depth + 1)}}}{{{expression(rng, depth + 1)}}}"
    if choice == 1:
        return rf"\sqrt{{{expression(rng, depth + 1)}}}"
    if choice == 2:
        return rf"\sqrt[{expression(rng, depth + 1)}]{{{expression(rng, depth + 1)}}}"
    if choice == 3:
        return rf"{{{expression(rng, depth + 1)}}}_{{{expression(rng, depth + 1)}}}"
    if choice == 4:
        return rf"{{{expression(rng, depth + 1)}}}^{{{expression(rng, depth + 1)}}}"
    if choice == 5:
        return rf"\sum_{{{expression(rng, depth + 1)}}}^{{{expression(rng, depth + 1)}}}{expression(rng, depth + 1)}"
    if choice == 6:
        return rf"\int_{{{expression(rng, depth + 1)}}}^{{{expression(rng, depth + 1)}}}{expression(rng, depth + 1)}"
    if choice == 7:
        return rf"\left({expression(rng, depth + 1)}\right)"
    if choice == 8:
        return rf"\hat{{{expression(rng, depth + 1)}}}"
    if choice == 9:
        return rng.choice(ENVIRONMENTS)
    if choice == 10:
        return expression(rng, depth + 1) + expression(rng, depth + 1)
    return "{" + expression(rng, depth + 1) + "}"


def damage(tex: str, rng: random.Random) -> str:
    """Make realistic incomplete/incorrect paste text without invalid UTF-8."""
    if not tex:
        return rng.choice(("{", "\\", r"\begin{matrix}"))
    choice = rng.randrange(8)
    position = rng.randrange(len(tex) + 1)
    if choice == 0:
        return tex[:position] + tex[position + 1:]
    if choice == 1:
        return tex[:position] + rng.choice(("{", "}", "_", "^", "&", "\\")) + tex[position:]
    if choice == 2:
        return tex[:position]
    if choice == 3:
        return tex + rng.choice(("{", "}", "_", "^", "\\"))
    if choice == 4:
        return tex.replace("\\", "\\unknowncommand ", 1)
    if choice == 5:
        return r"\begin{not-an-environment}" + tex + r"\end{different}"
    if choice == 6:
        return rng.choice(("$", "$$", r"\(", r"\[")) + tex
    return tex[:position] + "日本語🙂" + tex[position:]


def check(tex: str, seed: int, case: int, kind: str) -> None:
    normalized = tex_normalize(tex)
    normalized_twice = tex_normalize(normalized)
    if normalized != normalized_twice:
        raise AssertionError(
            f"normalization not fixed: seed={seed} case={case} kind={kind}\n"
            f"input={tex!r}\nonce={normalized!r}\ntwice={normalized_twice!r}"
        )

    equation = Equation()
    equation.load_latex(tex)
    canonical = equation.latex()
    if tex_normalize(canonical) != canonical:
        raise AssertionError(
            f"Equation output not canonical: seed={seed} case={case} kind={kind}\n"
            f"input={tex!r}\noutput={canonical!r}"
        )

    width, height, baseline = equation.metrics()
    if not all(math.isfinite(value) for value in (width, height, baseline)):
        raise AssertionError(
            f"non-finite metrics: seed={seed} case={case} kind={kind} tex={tex!r}"
        )
    if width < 0 or height <= 0 or baseline < 0 or baseline > height:
        raise AssertionError(
            f"invalid metrics {(width, height, baseline)!r}: "
            f"seed={seed} case={case} kind={kind} tex={tex!r}"
        )

    svg = tex_to_svg(tex)
    ET.fromstring(svg)
    viewbox = _VIEWBOX.search(svg)
    if not viewbox:
        raise AssertionError(f"missing SVG viewBox: seed={seed} case={case} kind={kind}")
    dimensions = tuple(float(value) for value in viewbox.groups())
    if not all(math.isfinite(value) and 0 < value < 1_000_000 for value in dimensions):
        raise AssertionError(
            f"invalid SVG dimensions {dimensions!r}: seed={seed} case={case} kind={kind}"
        )

    mathml = tex_to_mathml(tex)
    math_root = ET.fromstring(mathml)
    if math_root.tag != "{http://www.w3.org/1998/Math/MathML}math" or \
            math_root.attrib.get("mathsize") != "24pt":
        raise AssertionError(
            f"invalid 24 pt MathML root: seed={seed} case={case} kind={kind}"
        )


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")
    seeds = 64
    cases_per_seed = 24
    checks = 0
    for seed in range(seeds):
        rng = random.Random(0xE6410000 + seed)
        for case in range(cases_per_seed):
            valid = expression(rng)
            try:
                check(valid, seed, case, "generated")
                check(damage(valid, rng), seed, case, "damaged")
            except Exception as exc:
                print(f"FAIL seed={seed} case={case}: {exc}")
                return 1
            checks += 2

    # Pathologically deep nesting must not take the process down.  Parse,
    # layout, and emit are all recursive over the tree, so a deep enough paste
    # overflowed the stack -- measured fine to ~1200 levels of \sqrt and dead
    # by ~1400.  The parser now caps nesting (Equation Editor 3.0 does the
    # same, its string 16044), so every depth here must return, not crash.
    # A crash shows up as this test's process dying, not as an assertion.
    for depth in (500, 2000, 20000):
        deep = "x"
        for _ in range(depth):
            deep = r"\sqrt{" + deep + "}"
        try:
            check(deep, 0, 0, "deep-nest")
        except Exception as exc:
            print(f"FAIL deep nesting {depth}: {exc}")
            return 1
        checks += 1

    print(
        f"ok    {checks} raw TeX parse/normalize/metrics/SVG/MathML fuzz checks"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
