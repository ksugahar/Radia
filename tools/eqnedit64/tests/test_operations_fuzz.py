"""Deterministic structural-operation stress test for Eqnedit64.

This is intentionally headless: it drives the same Equation model used by the
native canvas without moving the user's mouse or sending keys to Windows.
"""
from __future__ import annotations

import math
import os
import random
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import Equation, MAX_NESTING_DEPTH  # noqa: E402


TEXT = tuple("abcxyz0123456789+-=(),") + ("π", "θ")
SYMBOLS = (
    r"\alpha", r"\beta", r"\gamma", r"\delta", r"\theta",
    r"\lambda", r"\mu", r"\pi", r"\sigma", r"\phi", r"\omega",
    r"\partial", r"\nabla", r"\infty", r"\times", r"\leq", r"\geq",
)
MOVES = (
    "caret.left", "caret.right", "caret.next_slot", "caret.prev_slot",
    "caret.out", "caret.home", "caret.end", "caret.up", "caret.down",
)
EDITS = (
    "edit.backspace", "edit.delete", "edit.undo", "edit.redo",
    "edit.new_line", "edit.alignment",
)
PASTES = (
    r"\frac{a}{b}", r"\sqrt{x+1}", r"\sum_{i=0}^{n}x_i",
    r"\left( a+b \right)",
    r"\begin{aligned}F&=ma\\E&=mc^{2}\end{aligned}",
)


def check_state(eq: "Equation", seed: int, step: int, trail: list[str]) -> None:
    tex = eq.latex()
    if not tex:
        return
    first = Equation()
    first.load_latex(tex)
    normalized = first.latex()
    second = Equation()
    second.load_latex(normalized)
    normalized_twice = second.latex()
    if normalized.strip() != normalized_twice.strip():
        raise AssertionError(
            f"normalization did not reach a fixed point at seed={seed} step={step}\n"
            f"trail={trail[-12:]}\nfirst={normalized!r}\nsecond={normalized_twice!r}"
        )
    width, height, baseline = eq.metrics()
    if not all(math.isfinite(v) for v in (width, height, baseline)):
        raise AssertionError(f"non-finite metrics at seed={seed} step={step}: {tex!r}")
    if width < 0 or height <= 0 or baseline < 0:
        raise AssertionError(
            f"invalid metrics at seed={seed} step={step}: "
            f"{(width, height, baseline)!r} tex={tex!r}"
        )
    if ":" not in eq.caret():
        raise AssertionError(f"invalid caret at seed={seed} step={step}: {eq.caret()!r}")


def run_seed(seed: int, steps: int) -> int:
    rng = random.Random(seed)
    eq = Equation()
    trail: list[str] = []
    templates = Equation.templates()

    for step in range(steps):
        # Keep the stress case bounded while still allowing deeply nested
        # templates.  Resetting is itself a realistic New-document boundary.
        if (len(eq.latex()) > 5000 or
                eq.caret().count("/") >= MAX_NESTING_DEPTH):
            eq = Equation()
            eq.insert_text("x")
            trail.append("reset")

        choice = rng.randrange(100)
        if choice < 30:
            value = rng.choice(TEXT)
            eq.insert_text(value)
            trail.append(f"text:{value}")
        elif choice < 45:
            kind = rng.choice(templates)
            eq.insert_template(kind)
            trail.append(f"template:{kind}")
        elif choice < 55:
            symbol = rng.choice(SYMBOLS)
            if not eq.insert_symbol(symbol):
                raise AssertionError(f"known symbol rejected: {symbol}")
            trail.append(f"symbol:{symbol}")
        elif choice < 76:
            command = rng.choice(MOVES)
            eq.command(command)
            trail.append(command)
        elif choice < 91:
            command = rng.choice(EDITS)
            eq.command(command)
            trail.append(command)
        elif choice < 96:
            if rng.randrange(2):
                eq.select_all()
                trail.append("select_all")
            else:
                eq.clear_selection()
                trail.append("clear_selection")
        else:
            pasted = rng.choice(PASTES)
            eq.insert_latex(pasted)
            trail.append(f"paste:{pasted}")

        if step % 25 == 0:
            check_state(eq, seed, step, trail)

    check_state(eq, seed, steps, trail)
    return steps


def main() -> int:
    seeds = 100
    steps = 250
    operations = 0
    for seed in range(seeds):
        try:
            operations += run_seed(seed, steps)
        except Exception as exc:
            print(f"FAIL seed={seed}: {exc}")
            return 1
    print(f"ok    {operations} deterministic structural operations across {seeds} seeds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
