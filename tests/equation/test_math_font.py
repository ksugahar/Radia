"""The typesetting parameters read out of the math font itself.

The layout used to guess these -- fifty-nine hand-tuned multiples of the point
size, none of them from the font -- and it shows most on the radical, which was
drawn as the plain character with a rule tacked on and so could not meet a tall
radicand.  No constant reaches that: the fix is to ask the font for a bigger
radical, and Cambria Math has been offering six of them all along.

Every expected number here was read independently with fontTools, so this is a
check of the reader against the font, not against itself.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation, "math_constants"):
    pytest.skip("built without the MATH reader", allow_module_level=True)

UPEM = 2048.0          # Cambria Math's design units per em
RADICAL = 0x221A


def units(em):
    return em * UPEM


# ---- the constants ----------------------------------------------------------

@pytest.mark.parametrize("key,want", [
    ("axis_height", 585),
    ("radical_vertical_gap", 166),
    ("radical_display_vertical_gap", 345),
    ("radical_rule_thickness", 133),
    ("radical_extra_ascender", 133),
    ("radical_kern_before_degree", 133),
    ("radical_kern_after_degree", -640),
    ("fraction_rule_thickness", 133),
    ("fraction_num_display_shift_up", 1550),
    ("fraction_denom_display_shift_down", 1370),
])
def test_a_constant_matches_the_font(key, want):
    got = units(equation.math_constants()[key])
    assert abs(got - want) < 0.5, "%s: read %.1f, font says %d" % (key, got, want)


def test_the_script_sizes_are_the_fonts_own():
    """71% and 50% are TeX's numbers; this font asks for 73% and 60%, and the
    whole point of reading the table is to use the designer's figure."""
    c = equation.math_constants()
    assert abs(c["script_percent"] - 0.73) < 1e-9
    assert abs(c["script_script_percent"] - 0.60) < 1e-9


def test_constants_are_in_em_not_design_units():
    """A caller multiplies by the type size, so the font's units must not leak
    out -- an axis height of 585 rather than 0.29 would put the fraction bar
    six hundred times too high."""
    assert 0.1 < equation.math_constants()["axis_height"] < 0.5


# ---- the radical is a family, not a character -------------------------------

def test_the_radical_resolves():
    assert equation.math_glyph(RADICAL) == 958


def test_an_absent_character_is_zero_not_a_guess():
    assert equation.math_glyph(0x0E0100) == 0


def test_the_font_ships_six_radicals():
    v = equation.math_stretch(RADICAL)["variants"]
    assert [round(units(h)) for _g, h in v] == [1972, 2544, 4569, 6829, 9129, 11429]


def test_the_radicals_get_taller():
    heights = [h for _g, h in equation.math_stretch(RADICAL)["variants"]]
    assert heights == sorted(heights)
    assert len(set(g for g, _h in equation.math_stretch(RADICAL)["variants"])) == 6


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

@pytest.mark.parametrize("want_em", [0.5, 0.9, 1.4, 2.5, 4.0])
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


def test_the_angle_bracket_is_the_mathematical_one():
    """U+2329 is the deprecated compatibility character and the font leaves it
    flat; the stretchy one is U+27E8.  Reaching for the obvious codepoint gets
    a bracket that will not grow."""
    assert not equation.math_stretch(0x2329)["variants"]
    assert equation.math_stretch(0x27E8)["variants"]
