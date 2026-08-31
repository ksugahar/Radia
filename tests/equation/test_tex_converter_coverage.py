"""Regression coverage for TeX constructs used by presentation sources.

These are converter tests, not source-style recommendations.  A supported TeX
spelling must reach OMML without becoming visible prose, and an unsupported
control sequence must fail before a slide can be generated.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")


@pytest.mark.parametrize(
    "latex",
    [
        r"\bigl(x\bigr)",
        r"\big(x\Big)",
        r"\Bigl[x\Bigr]",
        r"\bigg\{x\biggr\}",
        r"\Biggl\langle x\Biggr\rangle",
    ],
)
def test_delimiter_size_commands_keep_the_delimiter(latex):
    normalized = equation.tex_normalize(latex)
    assert "\\big" not in normalized
    assert "x" in normalized
    assert "<m:oMath" in equation.tex_to_omml(latex)


@pytest.mark.parametrize(
    "latex,forbidden",
    [
        (r"K_{\rm SIBC}", r"\rm"),
        (r"{\bf A}", r"\bf"),
        (r"{\it A}", r"\it"),
        (r"{\sf A}", r"\sf"),
        (r"{\tt A}", r"\tt"),
    ],
)
def test_old_font_declarations_are_semantic_not_literal(latex, forbidden):
    normalized = equation.tex_normalize(latex)
    assert forbidden not in normalized
    if forbidden == r"\rm":
        assert r"\mathrm" in normalized
    omml = equation.tex_to_omml(latex)
    assert forbidden not in omml
    assert "<m:oMath" in omml


@pytest.mark.parametrize(
    "latex,normalized_command,omml_character",
    [
        (r"\mathsf{T}", r"\mathsf", "\U0001d5b3"),
        (r"\mathtt{x}", r"\mathtt", "\U0001d6a1"),
        (r"\mathfrak{g}", r"\mathfrak", "\U0001d524"),
    ],
)
def test_extended_math_alphabets_reach_omml(
    latex, normalized_command, omml_character
):
    assert normalized_command in equation.tex_normalize(latex)
    assert omml_character in equation.tex_to_omml(latex)


def test_dots_alias_is_normalized_to_a_supported_ellipsis():
    normalized = equation.tex_normalize(r"a,\dots,z")
    assert r"\dots" not in normalized
    assert "…" in normalized
    assert "…" in equation.tex_to_omml(r"a,\dots,z")


@pytest.mark.parametrize("command", [r"\text{in }", r"\mbox{in }"])
def test_text_boxes_preserve_their_trailing_space(command):
    normalized = equation.tex_normalize(command)
    assert "in " in normalized
    assert 'xml:space="preserve"' in equation.tex_to_omml(command)


def _width(latex: str) -> float:
    style = equation.SvgStyle()
    style.padding = 0.0
    return equation.tex_metrics(latex, style)[0]


@pytest.mark.parametrize(
    "latex,minimum_extra_points",
    [
        (r"a\qquad b", 23.9),
        (r"a\ b", 3.9),
        (r"a\quad b", 11.9),
        (r"a\hspace{1em}b", 11.9),
    ],
)
def test_explicit_positive_spaces_survive_normalization_and_layout(
    latex, minimum_extra_points
):
    normalized = equation.tex_normalize(latex)
    assert normalized != "ab"
    assert _width(latex) - _width("ab") >= minimum_extra_points
    omml = equation.tex_to_omml(latex)
    assert 'xml:space="preserve"' in omml


def test_negative_thin_space_remains_an_explicitly_accepted_no_op():
    assert equation.tex_normalize(r"a\!b") == "ab"


def test_tilde_is_a_nonbreaking_space_not_the_sim_relation():
    normalized = equation.tex_normalize(r"\text{in}~\Omega")
    assert "~" in normalized
    assert r"\sim" not in normalized
    omml = equation.tex_to_omml(r"\text{in}~\Omega")
    assert "\u00a0" in omml
    assert 'xml:space="preserve"' in omml


def test_unknown_control_sequence_fails_loudly():
    with pytest.raises(ValueError, match=r"unsupported TeX control sequence: \\notacommand"):
        equation.tex_to_omml(r"a+\notacommand+b")
