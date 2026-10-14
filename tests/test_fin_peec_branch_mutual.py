"""The fin PEEC's branch mutual is filamentary, not cross-section averaged.

`MutualInductanceRectBar` averages the Neumann kernel over both cross-sections
and exists precisely because, in the kernel's own words, filamentary Neumann
gives a "spurious circulating current artifact ... for close parallel bars".
It is reached only when both segments are sub-filaments of one parent, which
`add_connected_segment` never produces, so every fin branch pair takes the
filamentary formula instead.

These checks pin that behaviour and its size, so gate 2 of the fin PEEC ledger
cannot be closed by accident.  If the fin path is ever routed through the
averaged kernel, the first test here changes and this file is the reminder to
re-measure the ledger's numbers.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation_test" / "induction_heating"))

#: The 48 mm fixture at the accepted discretisation: adjacent perimeter bands
#: touch, and their separation is smaller than their own depth.
BAND_LENGTH_M = 0.0015      # 48 mm over 32 axial segments
BAND_WIDTH_M = 0.0894e-3    # perimeter spacing at n_peri = 256
SHEET_DEPTH_M = 0.1706e-3   # skin depth, copper at 150 kHz
DELIVERY_TOLERANCE = 0.03   # the gate's current/loss fraction limit


def _pair_mutual():
    peec_matrices = pytest.importorskip("radia.peec_matrices")

    builder = peec_matrices.PEECBuilder()
    for shift in (0.0, BAND_WIDTH_M):
        n0 = builder.add_node_at(0.0, shift, 0.0)
        n1 = builder.add_node_at(0.0, shift, BAND_LENGTH_M)
        builder.add_connected_segment(n0, n1, BAND_WIDTH_M, SHEET_DEPTH_M,
                                      sigma=5.8e7)
    L, _R, _P, _M = builder.build()
    return float(L[0][1])


def test_branch_mutual_is_the_bare_filament_formula():
    """Machine-precision equality is the proof that the averaged kernel is
    not reached: an averaged value would differ in the third digit."""
    from fin_partial_element_geometry import filament_mutual

    kernel = _pair_mutual()
    grover = filament_mutual(BAND_LENGTH_M, BAND_WIDTH_M)
    assert kernel == pytest.approx(grover, rel=1e-12)


@pytest.mark.parametrize("gap_factor", [1.0, 2.0, 5.0])
def test_filament_formula_holds_at_every_separation(gap_factor):
    from fin_partial_element_geometry import filament_mutual

    peec_matrices = pytest.importorskip("radia.peec_matrices")

    gap = BAND_WIDTH_M * gap_factor
    builder = peec_matrices.PEECBuilder()
    for shift in (0.0, gap):
        n0 = builder.add_node_at(0.0, shift, 0.0)
        n1 = builder.add_node_at(0.0, shift, BAND_LENGTH_M)
        builder.add_connected_segment(n0, n1, BAND_WIDTH_M, SHEET_DEPTH_M,
                                      sigma=5.8e7)
    L, _R, _P, _M = builder.build()
    assert float(L[0][1]) == pytest.approx(
        filament_mutual(BAND_LENGTH_M, gap), rel=1e-12)


def test_the_filament_error_exceeds_the_delivery_tolerance():
    """Size the gap between what the kernel computes and what the geometry
    asks for, on the bands the accepted gate actually used."""
    import numpy as np

    from fin_partial_element_geometry import bar_mutual, filament_mutual

    axis = np.array([0.0, 0.0, 1.0])
    # width along the perimeter tangent, depth along the outward normal
    frame = (np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    centre_i = np.array([0.0, 0.0, BAND_LENGTH_M / 2.0])
    centre_j = np.array([0.0, BAND_WIDTH_M, BAND_LENGTH_M / 2.0])

    averaged = bar_mutual(centre_i, centre_j, axis, BAND_LENGTH_M,
                          BAND_WIDTH_M, SHEET_DEPTH_M,
                          BAND_WIDTH_M, SHEET_DEPTH_M, frame, frame)
    filament = filament_mutual(BAND_LENGTH_M, BAND_WIDTH_M)
    relative = filament / averaged - 1.0

    assert relative > DELIVERY_TOLERANCE, (
        f"the filamentary mutual is {relative:.2%} above the "
        "cross-section-averaged value on the accepted band geometry; if this "
        "ever drops below the delivery tolerance the gate-2 argument changes")
    assert relative < 0.10, (
        "the error is systematic and bounded; a jump here means the band "
        "geometry constants no longer describe the accepted configuration")


def test_separation_is_smaller_than_the_band_depth():
    """Why the filament formula is the wrong tool here, as one assertion."""
    assert BAND_WIDTH_M < SHEET_DEPTH_M
    assert BAND_LENGTH_M / BAND_WIDTH_M > 10.0, (
        "the bands are long and thin, so the error is in the cross-section "
        "treatment rather than in the end effects")
    assert math.isclose(BAND_WIDTH_M / SHEET_DEPTH_M, 0.524, rel_tol=0.01)
