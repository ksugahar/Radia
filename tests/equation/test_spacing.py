"""The space that goes around an operator.

TeX sorts every atom into a class -- ordinary, binary, relation, opening,
closing, punctuation -- and the space between two of them is a table lookup in
eighteenths of an em.  Equation Editor 3.1 does the same thing.  Measured on its
own screen at four times size, it leaves 0.219 em around a relation, 0.156
around a binary operator, nothing before a comma and 0.109 after one, and
nothing at a parenthesis: the same classes, the same asymmetry at the comma,
about three quarters of TeX's amounts.

The table is checked directly rather than through a rendered width.  A missing
entry is invisible at that distance -- the glyph's own advance drowns out the
space around it -- and that is exactly how an ASCII asterisk came to be classed
as an ordinary letter, with 0.03 em around it where Equation Editor gives 0.16,
because only U+2217 was listed and not '*' as anyone would type it.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation, "atom_kind"):
    pytest.skip("built without the atom-class accessor", allow_module_level=True)

K = equation.AtomKind
kind = equation.atom_kind
space = equation.atom_space_mu

THIN, MED, THICK = 3, 4, 5


# ---- what a character is ----------------------------------------------------

@pytest.mark.parametrize("ch,want", [
    ("+", K.Bin), ("-", K.Bin), ("*", K.Bin),
    ("=", K.Rel), ("<", K.Rel), (">", K.Rel), ("|", K.Rel),
    ("(", K.Open), ("[", K.Open), ("{", K.Open),
    (")", K.Close), ("]", K.Close), ("}", K.Close),
    (",", K.Punct), (";", K.Punct), (":", K.Punct),
    ("a", K.Ord), ("1", K.Ord), ("/", K.Ord),
])
def test_a_character_is_classed_as_it_reads(ch, want):
    assert kind(ord(ch)) == want, ch


@pytest.mark.parametrize("ascii_form,typographic", [
    ("*", 0x2217),      # asterisk operator
    ("-", 0x2212),      # minus sign
    ("|", 0x2225),      # parallel to
])
def test_the_keyboard_form_is_classed_like_the_typographic_one(ascii_form,
                                                               typographic):
    """A person types these; the palette inserts the others.  They are the same
    operator and must space the same, which is what was wrong: '*' was ordinary
    while U+2217 was binary."""
    assert kind(ord(ascii_form)) == kind(typographic), ascii_form


# ---- how much room goes between two of them ---------------------------------

def test_a_relation_gets_more_room_than_a_binary_operator():
    """Five eighteenths against four.  Equation Editor measures 0.219 against
    0.156 -- the same ordering, about three quarters the size."""
    assert space(K.Ord, K.Rel) > space(K.Ord, K.Bin) > 0


def test_a_binary_operator_gets_room_a_letter_does_not():
    assert space(K.Ord, K.Bin) == MED
    assert space(K.Ord, K.Ord) == 0


def test_a_relation_gets_the_thick_space():
    assert space(K.Ord, K.Rel) == THICK
    assert space(K.Rel, K.Ord) == THICK


def test_nothing_at_a_parenthesis():
    """Equation Editor leaves 0.016 em there, which is nothing."""
    assert space(K.Ord, K.Open) == 0
    assert space(K.Open, K.Ord) == 0
    assert space(K.Ord, K.Close) == 0


def test_punctuation_is_asymmetric():
    """Nothing before a comma, a thin space after it.  Equation Editor measures
    0.016 before and 0.109 after -- the same shape."""
    assert space(K.Ord, K.Punct) == 0
    assert space(K.Punct, K.Ord) == THIN


def test_two_relations_do_not_push_each_other_apart():
    """In "a = b = c" the space belongs between a relation and its operand, not
    between two relations that happen to meet."""
    assert space(K.Rel, K.Rel) == 0


def test_a_binary_operator_does_not_push_against_a_closing_bracket():
    assert space(K.Bin, K.Close) == 0


def test_an_operator_takes_a_thin_space():
    """A large operator -- an integral, a summation -- against its operand."""
    assert space(K.Op, K.Ord) == THIN
    assert space(K.Ord, K.Op) == THIN


# ---- and it reaches the page ------------------------------------------------

def x_positions(latex):
    import re
    return [float(m) for m in
            re.findall(r'\sx="(-?[0-9.]+)"',
                       equation.tex_to_svg(latex, equation.SvgStyle()))]


@pytest.mark.parametrize("op", ["+", "-", "*", "=", "<"])
def test_the_table_reaches_the_rendered_equation(op):
    """The classes are of no use if the layout does not consult them."""
    xs = x_positions("a%sb" % op)
    plain = x_positions("axb")
    assert len(xs) >= 3 and len(plain) >= 3
    assert (xs[2] - xs[0]) > (plain[2] - plain[0]), op
