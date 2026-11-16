"""Flux-line FEM post-processing helpers stay in radia-mcp.

Guards: ``trace_flux_line`` integrates dx/ds = B along a field line and the
trajectory and arclength plot recipes draw on caller-supplied matplotlib axes
without any document or figure-style package.
Assumes: scipy and matplotlib are installed (radia-mcp analysis extras).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from radia_mcp.radia_ngsolve.flux_line_postprocess import (  # noqa: E402
    plot_field_vs_arclength,
    plot_flux_line_trajectory,
    trace_flux_line,
)


def test_uniform_field_line_is_straight():
    s, xy = trace_flux_line(lambda x, y: (1.0, 0.0), x0=(0.0, 0.5),
                            s_span=(0.0, 2.0), n_eval=21)
    xy = np.asarray(xy)
    assert np.allclose(xy[:, 0], np.asarray(s), atol=1e-6), (xy[:, 0], s)
    assert np.allclose(xy[:, 1], 0.5, atol=1e-9), xy[:, 1]


def test_plot_recipes_draw_on_plain_axes():
    import matplotlib.pyplot as plt

    s, xy = trace_flux_line(lambda x, y: (-y, x), x0=(1.0, 0.0),
                            s_span=(0.0, np.pi), n_eval=50)
    xy = np.asarray(xy)
    fig, axes = plt.subplots(1, 2)
    plot_flux_line_trajectory(axes[0], xy[:, 0], xy[:, 1], source_xy=(1.0, 0.0))
    plot_field_vs_arclength(axes[1], s, -xy[:, 1], xy[:, 0])
    assert axes[0].lines and axes[1].lines
    plt.close(fig)
