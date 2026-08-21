r"""An integral inside a line still gets the limits written after that line.

Equation Editor writes an operator's limits as a separate block AFTER the thing
they belong to: a switch to subscript size, the limit lines, a switch to symbol
size, the operator's own glyph.  The reader reunites them.

That worked while the operator and its block were siblings and failed the moment
the operator sat one level deeper -- inside a LINE, with the block left in the
parent's list.  The nested path asked for an operator with no limits "yet" by
testing hasLower/hasUpper, but those record which SLOTS the template wrote, and
an integral written with variation 2 -- the common case -- puts its INTEGRAND in
the slot this reader calls `upper`.  Every such integral therefore looked like
it already had limits and could never be given the block that followed its line.

**This lives here and not in tests/ for a reason.**  The shape only occurs in
documents EQNEDT32 itself wrote; our own writer emits the tidy form, where the
limits are real slots.  A round-trip fixture built with `tex_to_mtef` passes
against the BROKEN reader too -- it was written, it passed, and it was deleted
for claiming a coverage it did not have.  Only a real .eqn exercises this.

The corpus is private, so nothing is committed here but the file names and the
properties their LaTeX must have.  Point RADIA_EQN_CORPUS at the directory, or
let it default to the lab path; the test skips when it is not there.
"""

from __future__ import annotations

import os
import pathlib
import sys

import pytest

for _p in pathlib.Path(__file__).resolve().parents:
    if (_p / "src" / "radia").exists():
        sys.path.insert(0, str(_p / "src"))
        break

equation = pytest.importorskip("radia.equation")

DEFAULT_CORPUS = pathlib.Path(
    r"W:\00_CAE\数式エディタ\python\eqnedt32\db\bidirectional_pass")


def corpus_dir():
    env = os.environ.get("RADIA_EQN_CORPUS")
    return pathlib.Path(env) if env else DEFAULT_CORPUS


# The documents these fixes recovered.  Equation Editor draws all of them
# correctly -- checked with validation_test/equation/render_in_ee3.ps1, which is
# what turned them from "malformed input" back into bugs worth fixing.
RECOVERED = [
    "harrington_ch2_integral_equation_conducting_plate.eqn",
    "harrington_ch2_lmn_matrix_element.eqn",
    "harrington_ch2_operator_l_plate.eqn",
]


@pytest.mark.parametrize("name", RECOVERED)
def test_the_limit_block_reaches_an_integral_inside_a_line(name):
    path = corpus_dir() / name
    if not path.exists():
        pytest.skip(f"corpus not present: {path}")

    latex = equation.mtef_to_latex(path.read_bytes())

    assert r"\scriptstyle" not in latex, (
        "a stray size marker survived, which is what a limit block that found "
        f"no owner leaves behind:\n  {latex}")
    assert r"\limits" in latex, (
        f"the limits went missing entirely:\n  {latex}")


def test_a_nudged_record_does_not_desynchronise_the_stream():
    """The content of a nudged equation is read, not skipped past.

    Nudge is TWO bytes; the reader skipped four, ate the start of whatever came
    next, and carried on emitting plausible records from misaligned bytes.  The
    equation below lost BOTH inner products of its fraction that way, and
    nothing said so -- the output was well-formed LaTeX for a different
    equation.
    """
    path = corpus_dir() / "perturbation_alpha_native.eqn"
    if not path.exists():
        pytest.skip(f"corpus not present: {path}")

    latex = equation.mtef_to_latex(path.read_bytes())

    # What Equation Editor draws: alpha = 1 - <f0,g> / <f0,M f0>
    assert "f_{0}" in latex, latex
    assert latex.count(r"\langle") >= 2, (
        f"an inner product went missing -- the stream desynchronised:\n  {latex}")
    assert "Mf_{0}" in latex or "M f_{0}" in latex, latex
    assert r"\scriptscriptstyle" not in latex, (
        f"size markers from misaligned bytes:\n  {latex}")

    # And the denominator is INSIDE the fraction.  Equation Editor writes
    # this one as numerator-in-template, denominator-immediately-after, with
    # no size markers between; the reader needed both halves of that before
    # the equation matched what EE3 draws: alpha = 1 - <f0,g> / <f0,M f0>
    assert "}{}" not in latex, (
        f"the denominator escaped the fraction:\n  {latex}")
    whole = (r"\dfrac{\left\langle  f_{0},g \right\rangle }"
             r"{\left\langle  f_{0},Mf_{0} \right\rangle }")
    assert whole in latex, f"expected the whole fraction:\n  got {latex}"


def test_the_corpus_carries_no_unowned_limit_block():
    """Nothing anywhere in the corpus leaves a limit block unclaimed.

    This is the sweep: it is what catches the same defect arriving in a
    document nobody has looked at.
    """
    root = corpus_dir()
    if not root.exists():
        pytest.skip(f"corpus not present: {root}")

    stray = [p.name for p in sorted(root.glob("*.eqn"))
             if r"\scriptstyle" in equation.mtef_to_latex(p.read_bytes())]
    assert not stray, f"documents with an unclaimed limit block: {stray}"
