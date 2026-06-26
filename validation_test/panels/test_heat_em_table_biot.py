"""Unit tests for the ``--ht-source biot`` core helpers in
calc_heat_with_em_table (_biot_Ht_on_surface /
_polyline_to_filament_segments).

Pure NumPy + analytical Biot-Savart -- no NGSolve mesh, no Cubit, no
STEP file.  These lock the NEW physics the biot ht-source adds:

  * incident-field tangential projection |H_t| = |H - (H.n)n|
  * the explicit --biot-image-factor scaling
  * the centerline polyline -> filament-segment bridge (incl. the
    closed-loop wrap-around)

validated against the closed-form circular-loop ON-AXIS field
``H_z(z) = I a^2 / (2 (a^2 + z^2)^{3/2})``.  The end-to-end heat solve
with a real coil STEP + workpiece .vol is a separate smoke (deferred).
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PANELS = os.path.join(os.path.dirname(os.path.dirname(HERE)),
                      "src", "radia", "panels")
if PANELS not in sys.path:
    sys.path.insert(0, PANELS)

import calc_heat_with_em_table as ht  # noqa: E402


def _loop_segments(a, n_seg=600, z0=0.0):
    """Full circular loop of radius ``a`` in plane z=z0, centered at
    the origin, as an (n_seg, 2, 3) filament-segment array."""
    th = np.linspace(0.0, 2 * math.pi, n_seg + 1)
    p = np.zeros((n_seg + 1, 3))
    p[:, 0] = a * np.cos(th)
    p[:, 1] = a * np.sin(th)
    p[:, 2] = z0
    return np.stack([p[:-1], p[1:]], axis=1)


def _loop_Hz_analytic(a, z, current):
    return current * a * a / (2.0 * (a * a + z * z) ** 1.5)


def test_polyline_to_filament_segments_open_and_closed():
    poly = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]], dtype=float)

    seg_open = ht._polyline_to_filament_segments(poly, closed=False)
    assert seg_open.shape == (2, 2, 3)
    assert np.allclose(seg_open[0], [[0, 0, 0], [1, 0, 0]])
    assert np.allclose(seg_open[1], [[1, 0, 0], [1, 1, 0]])

    seg_closed = ht._polyline_to_filament_segments(poly, closed=True)
    assert seg_closed.shape == (3, 2, 3)
    # wrap-around closes the loop: last seg = (poly[-1] -> poly[0]).
    assert np.allclose(seg_closed[-1], [[1, 1, 0], [0, 0, 0]])

    with pytest.raises(ValueError):
        ht._polyline_to_filament_segments(np.zeros((1, 3)), closed=False)


def test_biot_Ht_axial_field_matches_loop_formula():
    a, I, z = 0.05, 1000.0, 0.03
    segs = _loop_segments(a, n_seg=600)
    obs = np.array([[0.0, 0.0, z]])
    Hz = _loop_Hz_analytic(a, z, I)

    # Normal along the axis (z): on-axis H is purely axial, so the
    # tangential component is ~0.
    n_axial = np.array([[0.0, 0.0, 1.0]])
    Ht_axial = ht._biot_Ht_on_surface(segs, obs, n_axial, current=I)
    assert Ht_axial[0] < 1e-3 * Hz

    # Radial normal (x): the tangent plane contains the z axis, so the
    # full axial field is tangential -> |H_t| ~ |H_z|.
    n_radial = np.array([[1.0, 0.0, 0.0]])
    Ht_radial = ht._biot_Ht_on_surface(segs, obs, n_radial, current=I)
    assert Ht_radial[0] == pytest.approx(Hz, rel=5e-3)


def test_biot_Ht_image_factor_scales_linearly():
    a, I, z = 0.05, 500.0, 0.02
    segs = _loop_segments(a, n_seg=400)
    obs = np.array([[0.0, 0.0, z]])
    n_radial = np.array([[1.0, 0.0, 0.0]])
    Ht1 = ht._biot_Ht_on_surface(segs, obs, n_radial, current=I,
                                 image_factor=1.0)
    Ht2 = ht._biot_Ht_on_surface(segs, obs, n_radial, current=I,
                                 image_factor=2.0)
    assert Ht2[0] == pytest.approx(2.0 * Ht1[0], rel=1e-12)


def test_biot_Ht_zero_normal_returns_full_magnitude():
    a, I, z = 0.05, 800.0, 0.04
    segs = _loop_segments(a, n_seg=400)
    obs = np.array([[0.0, 0.0, z]])
    # A zero normal means no tangent plane is known -> the helper
    # returns the full incident magnitude |H| (here purely axial).
    Ht_zero = ht._biot_Ht_on_surface(segs, obs, np.zeros((1, 3)), current=I)
    Hz = _loop_Hz_analytic(a, z, I)
    assert Ht_zero[0] == pytest.approx(Hz, rel=5e-3)


def test_biot_Ht_current_scales_linearly():
    a, z = 0.05, 0.03
    segs = _loop_segments(a, n_seg=400)
    obs = np.array([[0.0, 0.0, z]])
    n_radial = np.array([[1.0, 0.0, 0.0]])
    Ht_1A = ht._biot_Ht_on_surface(segs, obs, n_radial, current=1.0)
    Ht_kA = ht._biot_Ht_on_surface(segs, obs, n_radial, current=1000.0)
    assert Ht_kA[0] == pytest.approx(1000.0 * Ht_1A[0], rel=1e-9)
