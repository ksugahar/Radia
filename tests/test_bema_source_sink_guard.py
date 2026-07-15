"""Unit test: the BEM-A source/sink current-path (shortcut) guard.

Background (2026-07-15): the ih_fem_kelvin demo coil's two end-caps are each
SPLIT in z.  The test helper's "two smallest faces by mass" auto-detect
labelled the z<0 and z>0 HALVES OF THE SAME cap as source/sink, so the
impedance-EFIE -- which minimises energy -- drove the terminal current
through the 2.6 mm conductor THICKNESS instead of around the 66 mm ring and
returned **L = 0.28 nH instead of ~90 nH** (PEEC on the same coil: 104.9 nH):
a 373x error that passed SILENTLY into a method comparison.

``check_source_sink_current_path`` catches it BEFORE the solve.  Note that
straight-line distance cannot: on that coil the broken (through-thickness)
source/sink separation is 2.4 mm and the correct (across-the-gap) one is
2.6 mm.  The discriminator is the geodesic ALONG the conductor surface --
the real current path (measured: broken 1.73 mm = 2.6% of extent, correct
202.9 mm = 3.1x extent).

These tests use a synthetic triangulated strip (no meshing / no NGSolve), so
they are fast and pin the guard's decision directly.
"""
from __future__ import annotations

import numpy as np
import pytest

from radia.bem.coil_inductance_ngsolve import check_source_sink_current_path


def _strip(n=20, length=1.0, width=0.05):
    """Triangulated flat strip along x: 2*n triangles, n+1 rungs.

    Returns (verts, tris).  Triangle 2i / 2i+1 live in the i-th quad, so
    the first quad is at x~0 and the last at x~length.
    """
    xs = np.linspace(0.0, length, n + 1)
    verts = []
    for x in xs:
        verts.append([x, 0.0, 0.0])
        verts.append([x, width, 0.0])
    verts = np.asarray(verts, dtype=float)
    tris = []
    for i in range(n):
        a, b = 2 * i, 2 * i + 1          # rung i   (bottom, top)
        c, d = 2 * i + 2, 2 * i + 3      # rung i+1 (bottom, top)
        tris.append([a, b, c])
        tris.append([b, d, c])
    return verts, np.asarray(tris, dtype=np.int64)


def test_guard_passes_when_current_traverses_the_conductor():
    """source at one end, sink at the other -> geodesic ~ extent -> OK."""
    verts, tris = _strip()
    src = np.zeros(len(tris), dtype=bool)
    snk = np.zeros(len(tris), dtype=bool)
    src[0] = src[1] = True        # first quad (x ~ 0)
    snk[-1] = snk[-2] = True      # last quad  (x ~ length)

    info = check_source_sink_current_path(verts, tris, src, snk)
    assert info is not None
    # The path must run the length of the strip, i.e. ~ the bbox extent.
    assert info["ratio"] > 0.5
    assert info["geodesic_m"] == pytest.approx(1.0, rel=0.2)


def test_guard_raises_on_shortcut_source_sink():
    """source and sink on the SAME station -> geodesic ~ 0 -> RAISE.

    This is the ih_fem_kelvin failure mode in miniature: both labels sit on
    one cross-section, so the current shortcuts instead of traversing.
    """
    verts, tris = _strip()
    src = np.zeros(len(tris), dtype=bool)
    snk = np.zeros(len(tris), dtype=bool)
    src[0] = True                 # both in the FIRST quad -- adjacent
    snk[1] = True

    with pytest.raises(ValueError) as exc:
        check_source_sink_current_path(verts, tris, src, snk)
    msg = str(exc.value)
    assert "degenerate" in msg
    assert "SHORTCUT" in msg
    # The diagnostic must name the actionable fix, not just complain.
    assert "FIX:" in msg and "opposite" in msg.lower()


def test_guard_raises_when_source_sink_disconnected():
    """Two disconnected components -> no current can flow at all."""
    v1, t1 = _strip(n=4, length=1.0)
    v2, t2 = _strip(n=4, length=1.0)
    v2 = v2 + np.array([0.0, 10.0, 0.0])          # far away, unconnected
    verts = np.vstack([v1, v2])
    tris = np.vstack([t1, t2 + len(v1)])
    src = np.zeros(len(tris), dtype=bool)
    snk = np.zeros(len(tris), dtype=bool)
    src[0] = True
    snk[len(t1)] = True                            # on the other component

    with pytest.raises(ValueError, match="no surface path connects"):
        check_source_sink_current_path(verts, tris, src, snk)


def test_guard_is_a_noop_without_labels():
    """Empty masks are the caller's own fail-fast check, not this guard's."""
    verts, tris = _strip()
    empty = np.zeros(len(tris), dtype=bool)
    assert check_source_sink_current_path(verts, tris, empty, empty) is None
