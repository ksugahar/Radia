"""A figure nobody can read must not reach a file.

Each case here is something that happened on one deck in one afternoon. None
of them is a rendering fault, so the existing paste-scale and dpi checks call
the same figures clean: the picture is correct and the information is not
recoverable from it.

The gate raises rather than warns, because a figure that is merely warned
about gets shipped.
"""

from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from radia_mcp.figure.tools import (  # noqa: E402
    enforce_readable,
    figure_readability_problems,
    lab_savefig,
)


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def test_a_healthy_figure_has_nothing_to_report():
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot([1, 2, 3], [1, 4, 9])
    ax.set_xlabel("N")
    assert figure_readability_problems(fig) == []


def test_a_panel_flattened_into_a_strip_is_refused():
    """Slide 9 was 20.8 cm wide and 6.5 cm tall in 11.4 cm of free height, so
    the curves were squeezed into a band and the separation between solvers
    could not be read."""
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
    problems = figure_readability_problems(fig)
    assert any("flattened into a strip" in p for p in problems), problems
    with pytest.raises(ValueError, match="unreadable"):
        enforce_readable(fig)


def test_text_covering_a_plotted_point_is_refused():
    """Slide 7's inline label sat on the s=0.75 marker, hiding it."""
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot([1, 2, 3], [1, 4, 9])
    ax.text(2, 4, "covers the marker", ha="center", va="center")
    problems = figure_readability_problems(fig)
    assert any("covers" in p and "plotted point" in p for p in problems), problems


def test_text_running_off_the_figure_is_refused():
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot([1, 2, 3], [1, 4, 9])
    ax.annotate("way off", xy=(3, 9), xytext=(400, -400),
                textcoords="offset points")
    problems = figure_readability_problems(fig)
    assert any("runs outside the figure" in p for p in problems), problems


def test_a_label_in_white_space_is_fine():
    """The fix for slide 7: the label moved into the gap the data leaves."""
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot([1, 2, 3], [1, 4, 9])
    ax.text(0.03, 0.92, "comp. charge", transform=ax.transAxes)
    assert figure_readability_problems(fig) == []


def test_lab_savefig_refuses_to_write_it(tmp_path):
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot([1, 2, 3], [1, 4, 9])
    out = tmp_path / "strip.png"
    with pytest.raises(ValueError, match="unreadable"):
        lab_savefig(fig, str(out), medium="presentation")
    assert not out.exists(), "the file must not be written"


def test_the_escape_hatch_has_to_be_asked_for(tmp_path):
    """It exists, and it has to be typed at the call site so the decision is
    visible where the figure is made.

    Note "paper": allow_unreadable waives THIS gate only. The on-page font
    floor is a separate contract and still applies, which is why a
    presentation-medium figure at default font sizes is refused either way.
    """
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot([1, 2, 3], [1, 4, 9])
    out = tmp_path / "strip.png"
    lab_savefig(fig, str(out), medium="paper", allow_unreadable=True)
    assert out.exists()


def test_the_escape_hatch_does_not_waive_the_font_floor(tmp_path):
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot([1, 2, 3], [1, 4, 9])
    with pytest.raises(ValueError):
        lab_savefig(fig, str(tmp_path / "x.png"), medium="presentation",
                    allow_unreadable=True)


def test_the_problems_are_named_not_just_counted():
    """A gate that says 'unreadable' without saying why gets switched off."""
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot([1, 2, 3], [1, 4, 9])
    try:
        enforce_readable(fig)
    except ValueError as exc:
        assert "axes 1" in str(exc)
        assert "floor" in str(exc)
    else:
        pytest.fail("expected a refusal")
