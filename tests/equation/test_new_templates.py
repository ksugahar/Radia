r"""The constructs added on 2026-08-19, from the editor's side.

The appearance half of each is pinned against TeX in test_tex_metrics.py.
This file is the other half: that the template exists, that the caret can
reach its holes, and that what it builds reads back as the LaTeX an author
would have typed.

Several of these were in the model already and drew nothing at all --
OversetNode had no layout case, so \overset was silently discarded, and the
three arrow selectors reached the rule-drawing path and produced an
UNDERLINE.  A missing feature announces itself; a wrong mark does not.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)


def fresh(latex=""):
    e = equation.Equation()
    if latex:
        e.load_latex(latex)
    return e


# ---- everything on the palette can actually be inserted --------------------

def test_every_palette_button_maps_to_a_template():
    known = set(equation.Equation.templates())
    for group in equation.template_palettes():
        for item in group.items:
            if item.is_template:
                assert item.command in known, (group.name, item.command)


@pytest.mark.parametrize("kind", [
    "prescript", "overset", "underset", "braket",
    "xrightarrow", "xleftarrow", "xleftrightarrow", "xrightleftharpoons",
    "dbloverline", "dblunderline",
    "paren_l", "paren_r", "bracket_l", "bracket_r",
    "brace_l", "brace_r", "abs_l", "abs_r", "angle_l", "angle_r",
])
def test_the_new_templates_insert_and_draw(kind):
    e = fresh()
    assert e.insert_template(kind), kind
    out = e.latex()
    assert out.strip(), kind
    # and it survives being written out and read back
    again = fresh(out)
    assert again.latex().strip() == out.strip(), (kind, out)


# ---- the notation rules ----------------------------------------------------

def test_a_vector_is_bold_italic():
    r"""\vec\bm is the base spelling of a vector here, so the letter under the
    arrow carries the face as well."""
    e = fresh("B")
    e.select_all()
    assert e.insert_template("vec")
    assert e.latex() == B + "vec{" + B + "bm{B}}"


def test_bm_and_mathbf_are_different_things():
    r"""\bm is bold italic and is what a vector takes; \mathbf is upright bold
    and is what a matrix name takes.  They were one face until the vector rule
    was written down, so both drew the same letter."""
    a = fresh(B + "bm{B}")
    b = fresh(B + "mathbf{B}")
    assert a.latex() == B + "bm{B}"
    assert b.latex() == B + "mathbf{B}"
    assert equation.tex_to_svg(B + "bm{B}") != equation.tex_to_svg(B + "mathbf{B}")


def test_the_default_fraction_is_dfrac():
    r"""The outermost fraction is DRAWN at display size, so it goes out as
    \dfrac: pasted into running text, a bare \frac would be smaller than the
    picture the author accepted."""
    assert fresh(B + "frac{a}{b}").latex() == B + "dfrac{a}{b}"


def test_but_a_fraction_inside_one_stays_plain():
    r"""LaTeX steps a nested \frac down by itself, so the two rules agree at
    every level and nothing has to be said twice."""
    out = fresh(B + "frac{" + B + "frac{a}{b}}{c}").latex()
    assert out == B + "dfrac{" + B + "frac{a}{b}}{c}"


def test_tfrac_keeps_its_own_name():
    assert fresh(B + "tfrac{a}{b}").latex() == B + "tfrac{a}{b}"


# ---- the ones that used to draw the wrong mark, or nothing -----------------

def test_overrightarrow_draws_an_arrow_and_not_a_rule():
    over = equation.tex_to_svg(B + "overrightarrow{AB}")
    line = equation.tex_to_svg(B + "overline{AB}")
    assert over != line
    # a rule is a <rect>; an arrow is a glyph
    assert line.count("<rect") > over.count("<rect")


def test_overset_survives_at_all():
    out = fresh(B + "overset{a}{=}").latex()
    assert "overset" in out
    assert B + "text" not in out          # the unknown-command path


def test_a_labelled_arrow_reads_back_as_one():
    assert fresh(B + "xrightarrow{f}").latex() == B + "xrightarrow{f}"
    assert fresh(B + "xrightarrow[g]{f}").latex() == B + "xrightarrow[g]{f}"


def test_a_prescript_attaches_to_what_follows():
    r"""{}^{14}_{6}C is one atom, not two scripts floating before a C."""
    e = fresh("{}^{14}_{6}C")
    assert e.latex() == "{}^{14}_{6}C"


def test_a_middle_bar_is_a_delimiter_and_grows():
    small = equation.tex_metrics(
        B + "left" + B + "langle a" + B + "middle| b" + B + "right" + B + "rangle",
        _style())
    tall = equation.tex_metrics(
        B + "left" + B + "langle " + B + "frac{a}{b}" + B + "middle| c"
        + B + "right" + B + "rangle",
        _style())
    # the taller content makes a taller box, which a typed | could not do
    assert tall[1] + tall[2] > small[1] + small[2] + 2.0


def test_one_bar_and_two_bars_are_not_the_same_delimiter():
    r"""``\middle|`` is a single bar and ``\middle\|`` is a double one.

    They were swapped, and silently: the C++ compared against "\|", which is
    not an escape, so the compiler read it as "|" -- the plain bar took the
    double-bar branch and the double bar fell through to the plain one.  The
    LaTeX round-tripped either way, because the emitter had the same broken
    literal, so only the drawing was wrong.
    """
    single = fresh(B + "left( a " + B + "middle| b " + B + "right)")
    double = fresh(B + "left( a " + B + "middle" + B + "| b " + B + "right)")
    assert single.latex() != double.latex()
    assert B + "middle| " in single.latex() or B + "middle|" in single.latex()
    assert B + "middle" + B + "|" in double.latex()

    # and what is drawn differs too -- U+2016 only for the double bar
    assert "‖" not in equation.tex_to_svg(single.latex())
    assert "‖" in equation.tex_to_svg(double.latex())


def test_Vert_is_the_double_bar_it_names():
    e = fresh(B + "left" + B + "langle a " + B + "middle" + B + "Vert b "
              + B + "right" + B + "rangle")
    assert "‖" in equation.tex_to_svg(e.latex())


def _style():
    s = equation.SvgStyle()
    s.padding = 0.0
    return s


def test_a_one_sided_fence_keeps_the_space_of_the_absent_one():
    r"""\left. f \right| is 1.2 pt wider than the bar and the f, because the
    delimiter that is not there still takes \nulldelimiterspace."""
    both = equation.tex_metrics(B + "left| f " + B + "right|", _style())
    one = equation.tex_metrics(B + "left. f " + B + "right|", _style())
    assert one[0] < both[0]
    assert one[0] > 0
