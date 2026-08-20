"""The icon each palette button wears.

A button draws a REAL member of its group, inserted by the real code, so it
cannot start advertising something the template no longer is.  Which member is
named by the group, because "first in the list" is Equation Editor's ordering
and not always the member that says what the group IS -- the matrix palette
opens on a 1x2, which reads as two boxes rather than as a matrix.

A name that is not in its group would silently fall back to the first item and
the button would quietly go back to being wrong.  That is not hypothetical: the
first version of the table reached the compiler as "\neq" with ONE backslash,
which C++ read as a newline followed by "eq", so every symbol icon silently
stayed on the member it was supposed to replace.  Nothing failed; the buttons
just did not change.  Hence this file.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation, "symbol_palettes"):
    pytest.skip("built before the palettes were exposed", allow_module_level=True)


def all_groups():
    return list(equation.symbol_palettes()) + list(equation.template_palettes())


def insert(item):
    e = equation.Equation()
    if item.is_template:
        e.insert_template(item.command)
    else:
        e.insert_symbol(item.command)
    return e


# The icons the groups are expected to wear.  This mirrors the C++ table in
# eq_edit.cpp; a group not listed here wears its first member.
EXPECTED = {
    "Relations": r"\neq",
    "Arrows": r"\rightarrow",
    "Logic": r"\forall",
    "Set theory": r"\subset",
    "Miscellaneous": r"\infty",
    "Embellishments": "vec",
    "Matrices": "matrix2x2",
}


@pytest.mark.parametrize("group_name,command", sorted(EXPECTED.items()))
def test_the_button_wears_the_named_icon(group_name, command):
    """What the WINDOW resolves, not merely what the group contains.

    Checking membership alone would pass whatever the table said, since a
    mis-spelled name falls back to the first item.  `icon` is the same call the
    window makes to decide what to draw.
    """
    groups = {g.name: g for g in all_groups()}
    assert group_name in groups, f"no palette group named {group_name!r}"
    g = groups[group_name]
    assert g.icon.command == command, (
        f"{group_name} draws {g.icon.command!r}, expected {command!r}; "
        f"members are {[it.command for it in g.items]}"
    )


def test_an_unnamed_group_wears_its_first_member():
    """The default, so a new group needs no table entry to look like itself."""
    for g in all_groups():
        if g.name in EXPECTED:
            continue
        assert g.icon.command == g.items[0].command, g.name


def test_every_group_has_something_to_draw():
    """An empty group would give a button with nothing on it."""
    for g in all_groups():
        assert g.items, f"{g.name} has no members"


def test_the_icon_renders_as_more_than_a_bare_point():
    """What the button shows has to have some extent.

    A template's slots are drawn as dotted boxes by the window, so a template
    whose LaTeX is all empty slots still shows something; what this rules out
    is an icon with no width AND no height, which paints an empty button.
    """
    for g in all_groups():
        item = g.items[0]
        for name, want in EXPECTED.items():
            if g.name == name:
                item = next(it for it in g.items if it.command == want)
        w, asc, desc = insert(item).extents()
        assert w > 0 or (asc + desc) > 0, f"{g.name}: {item.command!r} draws nothing"


def test_the_matrix_icon_is_square():
    """A 1x2 reads as two boxes; the word "matrices" means the 2x2."""
    groups = {g.name: g for g in all_groups()}
    tex = insert(next(it for it in groups["Matrices"].items
                      if it.command == "matrix2x2")).latex()
    rows = [r for r in tex.split(r"\\") if r.strip()]
    assert len(rows) == 2, tex
    assert rows[0].count("&") == 1, tex
