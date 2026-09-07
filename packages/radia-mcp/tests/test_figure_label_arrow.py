"""place_label_arrow keeps its arrow on the side of the label nearest the
target and stores the tail in axes fraction; place_legend_clear says when no
corner is clear instead of hiding the legend outside the canvas."""
from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from radia_mcp.figure.tools import place_label_arrow, place_legend_clear  # noqa: E402


def test_label_arrow_tail_sits_on_the_label_side_nearest_the_target():
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    x = np.linspace(0, 1, 50)
    ax.plot(x, x ** 2)
    t, ann = place_label_arrow(ax, "the knee", (0.7, 0.49), color="r",
                               region=(0.05, 0.5, 0.55, 0.95))
    fig.canvas.draw()
    bb = t.get_window_extent()
    sx, sy = ax.transAxes.transform(ann.xyann)
    assert bb.x0 - 6 <= sx <= bb.x1 + 6
    assert bb.y0 - 6 <= sy <= bb.y1 + 6
    # nearest side: the target is below-right of the label, so the tail is
    # on the bottom or right edge, never on the top or left one
    on_bottom = abs(sy - bb.y0) <= 6
    on_right = abs(sx - bb.x1) <= 6
    assert on_bottom or on_right
    plt.close(fig)


def test_label_arrow_tail_is_stored_in_axes_fraction():
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    ax.plot([0, 1], [0, 1])
    _, ann = place_label_arrow(ax, "here", (0.5, 0.5), color="b")
    tail = tuple(ann.xyann)
    assert ann.anncoords == "axes fraction"
    assert 0.0 <= tail[0] <= 1.0 and 0.0 <= tail[1] <= 1.0
    fig.dpi = 300
    fig.canvas.draw()
    assert tuple(ann.xyann) == tail
    plt.close(fig)


def test_legend_clear_names_the_outside_fallback():
    fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
    x = np.linspace(0, 1, 400)
    for k in range(8):                         # curves everywhere
        ax.plot(x, np.sin(11 * x + k), label=f"c{k}")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        loc = place_legend_clear(ax)
    assert loc == "outside"
    assert any("cut it off" in str(w.message) for w in caught)
    plt.close(fig)


def test_legend_clear_can_refuse_instead():
    fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
    x = np.linspace(0, 1, 400)
    for k in range(8):
        ax.plot(x, np.sin(11 * x + k), label=f"c{k}")
    assert place_legend_clear(ax, outside=False) is None
    assert ax.get_legend() is None
    plt.close(fig)


def test_legend_clear_finds_a_clear_corner_when_there_is_one():
    fig, ax = plt.subplots(figsize=(6, 4), dpi=100)
    ax.plot([0, 1], [1, 1], label="flat, at the top")
    ax.set_ylim(0, 1.2)
    assert place_legend_clear(ax) in ("lower left", "lower right")
    plt.close(fig)
