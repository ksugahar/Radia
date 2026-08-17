"""Japanese in an equation.

Two separate defects lived here, and both were invisible while only ASCII was
tried.

The parser advanced one BYTE at a time, so a three-byte character became three
nodes and every renderer drew it as three Latin-1 glyphs -- textbook mojibake.
Fixing that revealed the second: Times New Roman and Cambria Math do not
contain a single kana, so correctly decoded Japanese was drawn as correctly
spaced blank paper.

Both are locked here at the tree level, where a wrong answer is a wrong count
rather than a picture someone has to look at.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation

WORDS = ["磁", "磁束密度", "鉄心", "ひらがな", "カタカナ", "全角ＡＢ", "、。"]


# ---- one character, one node -----------------------------------------------

@pytest.mark.parametrize("word", WORDS)
def test_a_japanese_character_is_one_node_not_three(word):
    """A byte-at-a-time parse made this len(utf8) instead of len(text)."""
    tree = equation.tex_dump_tree(word)
    assert tree.count("CHAR") == len(word), tree


@pytest.mark.parametrize("word", WORDS)
def test_japanese_survives_the_round_trip(word):
    e = Equation()
    e.load_latex(word)
    assert e.latex() == word


def test_japanese_mixed_with_maths_keeps_both():
    e = Equation()
    e.load_latex(r"B_{鉄心} = \frac{\Phi}{S}")
    out = e.latex()
    assert "鉄心" in out
    assert r"\Phi" in out


# ---- it reaches every output ------------------------------------------------

@pytest.mark.parametrize("word", ["磁束密度", "鉄心"])
def test_mathml_carries_the_characters_themselves(word):
    """Each character is its own element, so the word is not one substring --
    what matters is that every character arrives, not as an escape."""
    xml = equation.tex_to_mathml(word)
    for ch in word:
        assert ch in xml


@pytest.mark.parametrize("word", ["磁束密度", "鉄心"])
def test_omml_carries_the_characters_themselves(word):
    xml = equation.tex_to_omml(word)
    for ch in word:
        assert ch in xml


def test_office_sets_japanese_upright():
    """A kanji in italic is simply wrong; Office is told plain explicitly."""
    assert 'm:val="p"' in equation.tex_to_omml("鉄")


def test_the_svg_asks_for_a_japanese_face():
    """Neither the text face nor the maths face has a kana in it, so a viewer
    given the default families would render nothing."""
    svg = equation.tex_to_svg("磁")
    assert "Mincho" in svg
    plain = equation.tex_to_svg("B")
    assert "Mincho" not in plain


def test_a_picture_of_japanese_is_not_the_picture_of_nothing():
    """The blank-paper failure produced a valid but nearly empty PNG."""
    st = equation.SvgStyle()
    jp = equation.tex_to_png("磁束密度", st, 4.0)
    blank = equation.tex_to_png(r"\ ", st, 4.0)
    assert len(jp) > len(blank) * 2


# ---- the editing model ------------------------------------------------------

def test_backspace_removes_a_whole_character_not_a_byte():
    e = Equation()
    e.load_latex("磁束")
    e.move_end()
    assert e.backspace()
    assert e.latex() == "磁"


def test_typing_japanese_inserts_whole_characters():
    """What the IME hands the window is a composed character, not bytes."""
    e = Equation()
    e.insert_text("鉄")
    e.insert_text("心")
    assert e.latex() == "鉄心"


def test_the_caret_can_be_placed_inside_japanese():
    e = Equation()
    e.load_latex("磁束密度")
    e.move_home()
    for _ in range(2):
        assert e.move_right()
    e.insert_text("X")
    assert e.latex() == "磁束X密度"
