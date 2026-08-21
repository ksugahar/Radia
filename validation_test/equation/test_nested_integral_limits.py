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


# The documents the fix recovered.  Equation Editor draws both correctly --
# checked with validation_test/equation/render_in_ee3.ps1, which is what turned
# these from "malformed input" back into a bug worth fixing.
RECOVERED = [
    "harrington_ch2_integral_equation_conducting_plate.eqn",
    "harrington_ch2_lmn_matrix_element.eqn",
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


def test_the_corpus_carries_no_unowned_limit_block():
    """Nothing anywhere in the corpus leaves a limit block unclaimed.

    The per-document checks above name two files; this one is the sweep, and it
    is the check that would catch the same defect arriving somewhere else.
    """
    root = corpus_dir()
    if not root.exists():
        pytest.skip(f"corpus not present: {root}")

    stray = []
    for path in sorted(root.glob("*.eqn")):
        latex = equation.mtef_to_latex(path.read_bytes())
        if r"\scriptstyle" in latex:
            stray.append(path.name)

    # Three remain, each its own unfinished investigation; they are listed so
    # that a NEW one fails this test instead of hiding in the count.
    known = {
        "harrington_ch2_operator_l_plate.eqn",
        "harrington_ch2_polarizability_xx_matrix.eqn",
        "perturbation_alpha_native.eqn",
    }
    unexpected = sorted(set(stray) - known)
    assert not unexpected, f"new documents with an unclaimed limit block: {unexpected}"
    fixed = sorted(known - set(stray))
    assert not fixed, (
        f"these are clean now -- remove them from `known`: {fixed}")
