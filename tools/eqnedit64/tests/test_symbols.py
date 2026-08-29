"""Every named symbol, swept for the two ways a saved .tex silently rots.

The raw-TeX fuzzer already checks that normalization reaches a fixed point,
but its corpus mentions twelve commands out of the symbol table's hundred and
fifty.  Both defects this file exists for lived in the gap:

  * `\\rangle` was the one multi-letter row in the emitter's code-point table
    written without a trailing space, so `\\langle a \\rangle x` serialized to
    `\\langle a\\ranglex` and came back from the file as the word "ranglex".
    Undo restores a snapshot by re-parsing its own output, so the same hole
    corrupted Undo, not just save and reopen;
  * forty-five commands had no reverse entry at all and were written as bare
    Unicode glyphs -- `⟹`, `⌊`, `†`, `…` -- which are not valid pdfLaTeX
    math.  The equation still looked right on the canvas and the file no
    longer compiled.

So the invariants are: a symbol survives its own round trip unchanged, one
pass of normalization is a fixed point, and the serialized form is plain
ASCII TeX rather than a glyph that only a Unicode-aware engine can set.

Run:  python tests\\test_symbols.py
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import symbol_commands, tex_normalize  # noqa: E402


def sweep_symbols(failures: list) -> int:
    commands = symbol_commands()
    if len(commands) < 100:
        failures.append(f"symbol table looks truncated: {len(commands)} entries")
    for command in commands:
        # 1. one pass of normalization is a fixed point, in a context where a
        #    missing separator would run the command into the next atom.
        once = tex_normalize(command + " x")
        twice = tex_normalize(once)
        if once != twice:
            failures.append(
                f"{command}: not a fixed point, {once!r} -> {twice!r}")

        # 2. the saved form is TeX, not a raw glyph pdfLaTeX cannot set.
        alone = tex_normalize(command)
        raw = [c for c in alone if ord(c) > 126]
        if raw:
            failures.append(
                f"{command}: emitted as raw glyph {''.join(raw)!r} ({alone!r})")

        # 3. the author's own command comes back, not a synonym that happens
        #    to share the code point.
        elif alone.strip() != command:
            failures.append(f"{command}: round-tripped to {alone.strip()!r}")
    return len(commands)


def check_named_glyphs_beat_operator_names(failures: list) -> None:
    """\\div is the division sign; \\Re and \\Im are Fraktur glyphs."""
    for command in (r"\div", r"\Re", r"\Im"):
        out = tex_normalize(command)
        if out.strip() != command:
            failures.append(
                f"{command} should stay a glyph command, got {out!r}")
    # Genuine operator names must still be set upright.
    for command in (r"\sin", r"\log", r"\lim"):
        out = tex_normalize(command + " x")
        if command not in out:
            failures.append(f"{command} lost its operator form: {out!r}")


def check_structures(failures: list) -> None:
    cases = [
        # cases keeps both columns and its environment.
        (r"\begin{cases} x & x > 0 \\ -x & x \leq 0 \end{cases}",
         [r"\begin{cases}", "&"]),
        # a mismatched closing delimiter is meaningful, not a typo.
        (r"\left( a , b \right]", [r"\left(", r"\right]"]),
        (r"\left[ a , b \right)", [r"\left[", r"\right)"]),
        # \left. x \right. is a sizing group; it must not grow parentheses.
        (r"\left. x \right.", []),
        # matched pairs are unchanged.
        (r"\left( x \right)", [r"\left(", r"\right)"]),
    ]
    for tex, required in cases:
        out = tex_normalize(tex)
        if out != tex_normalize(out):
            failures.append(f"{tex!r}: not a fixed point, {out!r}")
        for token in required:
            if token not in out:
                failures.append(f"{tex!r}: lost {token!r}, got {out!r}")

    invisible = tex_normalize(r"\left. x \right.")
    if "(" in invisible or ")" in invisible:
        failures.append(
            rf"\left. x \right. invented parentheses: {invisible!r}")

    columns = tex_normalize(
        r"\begin{cases} x & x > 0 \\ -x & x \leq 0 \end{cases}")
    if columns.count("&") != 2:
        failures.append(f"cases lost a column separator: {columns!r}")

    # A limits operator must not be wrapped in braces: {\lim}_{x} is an
    # ordinary symbol, so the subscript moves from under the operator to
    # beside it.  Compound bases still need their braces.
    limits = tex_normalize(r"\lim_{x \to 0} f(x)")
    if limits.startswith("{"):
        failures.append(f"\\lim was demoted to an ordinary atom: {limits!r}")
    for tex in (r"{a+b}^{2}", r"{x^{2}}_{3}"):
        out = tex_normalize(tex)
        if not out.startswith("{"):
            failures.append(f"{tex}: compound base lost its braces: {out!r}")

    # \middle used to fall through to the unknown-command path and print the
    # letters of its own name.  An angle-fenced bra-ket is now a Dirac node,
    # whose bar stretches; any other \middle keeps the delimiter inline.
    braket = tex_normalize(r"\left\langle a \middle| b \right\rangle")
    if r"\text{middle}" in braket or "middle}" in braket:
        failures.append(f"\\middle leaked its command name: {braket!r}")
    if r"\middle|" not in braket:
        failures.append(f"bra-ket lost its stretching bar: {braket!r}")
    for tex in (r"\left\langle a \middle| b \right\rangle",
                r"\left( a \middle| b \right)",
                r"\left\langle a \middle| b \middle| c \right\rangle"):
        out = tex_normalize(tex)
        if "text{middle" in out:
            failures.append(f"{tex!r}: \\middle leaked as text: {out!r}")
        if "|" not in out:
            failures.append(f"{tex!r}: lost its separator: {out!r}")
        if out != tex_normalize(out):
            failures.append(f"{tex!r}: not a fixed point, {out!r}")

    # LaTeX's \phi is the straight symbol and \varphi the ordinary letter;
    # six var- commands used to share a code point with their plain form, so
    # the canvas and Office both showed the wrong glyph.
    variants = {r"\phi": r"\varphi", r"\epsilon": r"\varepsilon",
                r"\theta": r"\vartheta", r"\pi": r"\varpi",
                r"\kappa": r"\varkappa", r"\rho": r"\varrho"}
    for plain, variant in variants.items():
        both = tex_normalize(plain + " " + variant)
        if plain not in both or variant not in both:
            failures.append(f"{plain}/{variant} did not both survive: {both!r}")
        if tex_normalize(plain) == tex_normalize(variant):
            failures.append(f"{plain} and {variant} still share one glyph")

    # A damaged escaped literal can pass through the legacy ASCII-symbol map
    # before becoming a named command.  Outer padding must be normalized in
    # that first pass, rather than disappearing only after save/reopen.
    damaged = tex_normalize(r"\^leq")
    if damaged != tex_normalize(damaged):
        failures.append(
            f"escaped symbol was not a one-pass fixed point: {damaged!r}")
    if damaged[:1].isspace():
        failures.append(f"escaped symbol kept outer whitespace: {damaged!r}")


def check_serialisation_traps(failures: list) -> None:
    """Places where the emitter can write TeX it cannot read back."""
    # TeX ends an optional argument at the first unbraced ']', so a root index
    # containing one has to be wrapped.  \sqrt[{]}]{x} used to come back as
    # \sqrt{]}x -- index gone, wrong radicand, x outside the root.
    for tex in (r"\sqrt[{]}]{x}", r"\sqrt[{]}^{2}]{x}", r"\sqrt[{[}]{x}",
                r"\sqrt[n]{x}", r"\sqrt{x}"):
        out = tex_normalize(tex)
        if out != tex_normalize(out):
            failures.append(
                f"{tex!r}: root index not a fixed point, {out!r}")

    # One character node can serialise to several tokens.  A text or vector
    # character under a script came out unbraced, while the same content read
    # back from a file came out braced, so the file changed every other time
    # it was opened.
    for tex in (r"\text{a}^{2}", r"\mathrm{abc}^{2}",
                r"\mathit{abc}^{2}", r"\mathrm{\frac{x}{y}}",
                r"\mathit{\sqrt{x}}", r"\mathsf{abc}", r"\mathtt{abc}",
                r"\mathcal{ABC}", r"\mathbb{R}", r"\mathfrak{F}",
                r"\boldsymbol{\alpha}", r"\bm{\alpha}",
                r"\mathnormal{x}", r"\mathbf{x}^{2}", r"\text{速度}_{1}",
                r"\operatorname{Re}^{2}", r"{a+b}^{2}", r"x^{2}",
                r"\alpha^{2}", r"\lim_{x}"):
        out = tex_normalize(tex)
        if out != tex_normalize(out):
            failures.append(f"{tex!r}: script base not a fixed point, {out!r}")

    # cases owns its brace and left-aligns; a plain matrix in a brace is a
    # different thing and must not be rewritten into one.
    plain = tex_normalize(r"\left\{ \begin{matrix} a & b \end{matrix} \right.")
    if "cases" in plain:
        failures.append(f"a braced matrix was rewritten as cases: {plain!r}")
    if plain != tex_normalize(plain):
        failures.append(f"a braced matrix is not a fixed point: {plain!r}")


def main() -> int:
    failures: list = []
    count = sweep_symbols(failures)
    check_serialisation_traps(failures)
    check_named_glyphs_beat_operator_names(failures)
    check_structures(failures)
    if failures:
        print(f"FAIL  {len(failures)}")
        for failure in failures[:40]:
            print("  " + failure)
        if len(failures) > 40:
            print(f"  ...and {len(failures) - 40} more")
        return 1
    print(f"ok    {count} symbols: round trip, fixed point, ASCII TeX output; "
          "cases columns, mixed and invisible fences")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
