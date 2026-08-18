"""The typesetting parameters read out of the math font itself.

The font is Latin Modern Math -- Computer Modern in OpenType.  Its MATH table
is not merely close to TeX's parameters, it IS them: asked directly, TeX at 12
point wants an axis at 8.124/12 = 0.677 em for a display numerator and a rule
at 0.040, and the table says 677 and 40 per thousand.  Setting in it does not
approximate the way TeX sets mathematics; it uses the same numbers.

The layout used to guess these -- fifty-nine hand-tuned multiples of the point
size, none of them from the font -- and it shows most on the radical, which was
drawn as the plain character with a rule tacked on and so could not meet a tall
radicand.  No constant reaches that: the fix is to ask the font for a bigger
radical, and the font has been offering a family of them all along.

Every expected number here was read independently with fontTools, so this is a
check of the reader against the font, not against itself.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation, "math_constants"):
    pytest.skip("built without the MATH reader", allow_module_level=True)

UPEM = 1000.0          # Latin Modern Math's design units per em
RADICAL = 0x221A


def units(em):
    return em * UPEM


# ---- the constants ----------------------------------------------------------

@pytest.mark.parametrize("key,want", [
    ("axis_height", 250),                     # TeX's axis_height, 3.0/12 em
    ("radical_vertical_gap", 50),
    ("radical_display_vertical_gap", 148),
    ("radical_rule_thickness", 40),           # TeX's default_rule_thickness
    ("radical_extra_ascender", 40),
    ("radical_kern_before_degree", 278),
    ("radical_kern_after_degree", -556),
    ("fraction_rule_thickness", 40),
    ("fraction_num_display_shift_up", 677),   # TeX's num1, 8.12389/12 em
    ("fraction_denom_display_shift_down", 686),   # TeX's denom1, 8.23198/12
])
def test_a_constant_matches_the_font(key, want):
    got = units(equation.math_constants()[key])
    assert abs(got - want) < 0.5, "%s: read %.1f, font says %d" % (key, got, want)


def test_the_script_sizes_are_the_fonts_own():
    """70% and 50%, which is what TeX has always used.  Cambria Math asks for
    73 and 60 instead -- the reason to read the table rather than a textbook is
    that the answer depends on the font."""
    c = equation.math_constants()
    assert abs(c["script_percent"] - 0.70) < 1e-9
    assert abs(c["script_script_percent"] - 0.50) < 1e-9


def test_constants_are_in_em_not_design_units():
    """A caller multiplies by the type size, so the font's units must not leak
    out -- an axis height of 585 rather than 0.29 would put the fraction bar
    six hundred times too high."""
    assert 0.1 < equation.math_constants()["axis_height"] < 0.5


# ---- the radical is a family, not a character -------------------------------

def test_the_radical_resolves():
    assert equation.math_glyph(RADICAL) == 3077


def test_an_absent_character_is_zero_not_a_guess():
    assert equation.math_glyph(0x0E0100) == 0


def test_the_font_ships_five_radicals():
    v = equation.math_stretch(RADICAL)["variants"]
    assert [round(units(h)) for _g, h in v] == [1001, 1201, 1801, 2401, 3001]


def test_the_radicals_get_taller():
    heights = [h for _g, h in equation.math_stretch(RADICAL)["variants"]]
    assert heights == sorted(heights)
    assert len(set(g for g, _h in equation.math_stretch(RADICAL)["variants"])) == 5


def test_there_are_parts_to_build_a_taller_one():
    """Past the largest ready-made size the font hands over three pieces -- a
    top, a repeating middle, a tail -- so a radical can be any height at all."""
    a = equation.math_stretch(RADICAL)["assembly"]
    assert len(a) == 3
    assert [ext for _g, _adv, ext in a] == [False, True, False]


def test_parts_overlap_by_the_amount_the_font_asks():
    """Butting the pieces end to end leaves a seam; the font states how much
    they must overlap."""
    assert units(equation.math_stretch(RADICAL)["min_overlap"]) > 0


# ---- what the layout will actually ask for ----------------------------------

@pytest.mark.parametrize("want_em", [0.5, 0.9, 1.4, 2.5, 3.0])
def test_asking_by_height_gives_something_at_least_that_tall(want_em):
    _g, got = equation.math_variant_for_height(RADICAL, want_em)
    assert got >= want_em, "asked %.1f em, got %.2f" % (want_em, got)


def test_a_taller_radicand_gets_a_different_drawing():
    """The failure this replaces: one radical glyph for every radicand, so a
    fraction inside a root left the bar joined to nothing."""
    small, _ = equation.math_variant_for_height(RADICAL, 0.5)
    large, _ = equation.math_variant_for_height(RADICAL, 3.0)
    assert small != large


def test_beyond_the_largest_it_reports_falling_short():
    """Then the caller must assemble; silently returning the biggest as if it
    fit is how a too-short radical would come back."""
    _g, got = equation.math_variant_for_height(RADICAL, 20.0)
    assert got < 20.0


# ---- other things that stretch ----------------------------------------------

@pytest.mark.parametrize("cp,name", [
    (0x0028, "parenthesis"), (0x005B, "bracket"), (0x007B, "brace"),
    (0x2223, "vertical bar"), (0x27E8, "angle bracket"),
    (0x222B, "integral"), (0x2211, "summation"),
])
def test_the_delimiters_stretch_too(cp, name):
    """One mechanism, every growing thing -- which is why this replaces the
    per-symbol tuning rather than adding to it."""
    s = equation.math_stretch(cp)
    assert s["variants"] or s["assembly"], name


def test_the_mathematical_angle_bracket_stretches():
    """Latin Modern maps the deprecated U+2329 to the same drawing and lets it
    stretch too; Cambria Math leaves that one flat.  U+27E8 is the one to reach
    for either way."""
    assert equation.math_stretch(0x27E8)["variants"]
