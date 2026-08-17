"""The palette bar, following Equation Editor 3.0's two rows.

The palette is BUILT from the shared command table, never listed by hand, and
these tests exist for the two ways that can still go wrong.

It can hide part of the editor: classifying only the ranges someone thought of
left \\ldots, \\cdots, \\dagger and \\bullet insertable from the keyboard but
absent from the bar.  So the classification has to be total.

And it can mislabel itself: the group names are a parallel array indexed by the
classifier's enum, so inserting a group into one and not the other silently
relabels every palette after it -- which is how "Greek" came to hold arrows.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def group(groups, name):
    return next(g for g in groups if g.name == name)


# ---- nothing is hidden -----------------------------------------------------

REACHABLE = {i.command for g in equation.symbol_palettes() for i in g.items}

# A spread of commands the editor accepts; each must be on the bar somewhere.
INSERTABLE = [
    r"\alpha", r"\Omega", r"\infty", r"\partial", r"\nabla",
    r"\perp", r"\parallel", r"\cdots", r"\ldots", r"\vdots", r"\ddots",
    r"\dagger", r"\bullet", r"\star", r"\times", r"\oplus",
    r"\forall", r"\exists", r"\emptyset", r"\subset", r"\leq", r"\neq",
    r"\rightarrow", r"\Leftrightarrow", r"\approx", r"\equiv",
    r"\hookleftarrow", r"\hookrightarrow",
]


@pytest.mark.parametrize("command", INSERTABLE)
def test_what_the_editor_can_insert_is_on_the_bar(command):
    e = Equation()
    assert e.insert_symbol(command), f"{command} is not insertable at all"
    assert command in REACHABLE, f"{command} is insertable but the bar hides it"


def test_the_bar_never_offers_what_cannot_be_inserted():
    """This also guards the command table's sort order.

    The lookup is a binary search, so a single pair out of order makes an
    entry that is IN the table unreachable.  \\hookleftarrow sat behind
    \\hookrightarrow and had never been insertable.  Because the palette
    enumerates the table and the classification is total, walking the bar and
    inserting each item covers every entry -- so any future inversion fails
    here rather than being discovered by a user who wanted an arrow.
    """
    for g in equation.symbol_palettes():
        for item in g.items:
            if item.is_template:
                assert item.command in Equation.templates()
            else:
                e = Equation()
                assert e.insert_symbol(item.command), item.command


def test_every_template_on_the_bar_actually_exists():
    known = set(Equation.templates())
    for g in equation.template_palettes():
        for item in g.items:
            assert item.is_template
            assert item.command in known, item.command


# ---- the labels tell the truth ---------------------------------------------

LABELLED = [
    ("Greek",          "\\alpha"),
    ("Greek capitals", "\\Omega"),
    ("Arrows",         "\\rightarrow"),
    ("Relations",      "\\leq"),
    ("Set theory",     "\\subset"),
    ("Logic",          "\\forall"),
    ("Ellipses",       "\\cdots"),
    ("Operators",      "\\pm"),
]


@pytest.mark.parametrize("name,command", LABELLED)
def test_a_group_contains_what_its_name_says(name, command):
    """The parallel-array bug shows up here and nowhere else."""
    g = group(equation.symbol_palettes(), name)
    assert command in {i.command for i in g.items}


def test_greek_holds_greek_and_nothing_else():
    g = group(equation.symbol_palettes(), "Greek")
    for item in g.items:
        assert 0x03B1 <= item.code <= 0x03C9, f"{item.command} is not Greek"


def test_arrows_hold_arrows():
    g = group(equation.symbol_palettes(), "Arrows")
    for item in g.items:
        assert (0x2190 <= item.code <= 0x21FF or
                0x27F8 <= item.code <= 0x27FA), item.command


# ---- shape ------------------------------------------------------------------

def test_no_group_is_empty():
    """An empty dropdown is worse than an absent one."""
    for g in equation.symbol_palettes() + equation.template_palettes():
        assert g.items, g.name


def test_no_symbol_appears_on_two_buttons():
    seen = set()
    for g in equation.symbol_palettes():
        for item in g.items:
            if item.is_template:
                continue
            assert item.command not in seen, item.command
            seen.add(item.command)


def test_the_two_rows_are_the_shape_of_the_old_editor():
    """Two rows of dropdowns, symbols above templates."""
    assert 6 <= len(equation.symbol_palettes()) <= 12
    assert 5 <= len(equation.template_palettes()) <= 10


def test_a_palette_cell_knows_what_to_draw():
    """Cells are drawn from the code point with our own routine, not from a
    bitmap -- so they stay sharp at any DPI and show exactly what gets
    inserted."""
    for g in equation.symbol_palettes():
        for item in g.items:
            if not item.is_template:
                assert item.code > 0, item.command
