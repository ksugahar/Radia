"""End-to-end test: STEP -> centerline -> filaments -> PEEC -> inductance.

Verifies that the full pipeline (coil_from_step.extract_centerline ->
coil_from_cad.build_peec_from_path -> PEECBuilder -> L) produces correct
inductance for a circular torus coil with known analytical self-inductance.

Analytical reference (Grover, circular loop, wire radius a, loop radius R):
    L = mu_0 * R * [ln(8R/a) - 2 + Y_i]
where Y_i = 1/4 for DC (uniform current distribution).

Test cases:
  1. Circular cross-section torus: R=30mm, a=3mm -> L_ana = 79.73 nH
  2. Square cross-section torus: R=30mm, w=h=6mm -> L_ana from GMD
"""

import math
import os
import sys
import tempfile

import numpy as np
import pytest

import radia  # noqa: F401 — must be first (sets up MKL DLL paths)

MU_0 = 4e-7 * math.pi


def _analytical_L_circular(R, a):
    """Grover self-inductance of a circular loop with circular cross-section.

    L = mu_0 * R * [ln(8R/a) - 2 + 1/4]
    The 1/4 = DC internal inductance.
    """
    return MU_0 * R * (math.log(8 * R / a) - 2.0 + 0.25)


def _analytical_L_square(R, w):
    """Grover self-inductance of a circular loop with square cross-section.

    Uses the GMD (geometric mean distance) of a square: GMD = 0.44705*w
    L = mu_0 * R * [ln(8R/GMD) - 2]
    (internal inductance is already embedded in GMD for square).

    Reference: Rosa NBS 169 (1908), Grover Dover (1946).
    For a square conductor with uniform current:
      ln(GMD) = ln(w) - 0.5*ln(2) - pi/12 + 0.5*ln(pi)
      GMD/w = exp(-0.5*ln(2) - pi/12 + 0.5*ln(pi)) = 0.447049...
    But we also need to add internal inductance Y_i = 0.25 (DC).
    So: L = mu_0 * R * [ln(8R/a_eq) - 2 + 0.25]
    where a_eq = GMD_square = w * exp(-ln(2)/2 - pi/12 + ln(pi)/2)
    """
    # For square w x w cross-section, equivalent radius via GMD:
    # Rosa NBS 169: for a square bar, GMD from itself = 0.44705 * w
    gmd = 0.44705 * w
    # External + internal (1/4)
    return MU_0 * R * (math.log(8 * R / gmd) - 2.0 + 0.25)


def _build_torus_step(path, R, a):
    """Build circular-section torus STEP via OCC."""
    from netgen.occ import WorkPlane, Axes, Pnt, Dir, Axis
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circ = wp.Circle(a).Face()
    torus = circ.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    torus.WriteStep(path)
    return path


def _build_square_torus_step(path, R, w):
    """Build square-section torus STEP via OCC."""
    from netgen.occ import WorkPlane, Axes, Pnt, Dir, Axis
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    rect = wp.MoveTo(-w / 2, -w / 2).Rectangle(w, w).Face()
    torus = rect.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    torus.WriteStep(path)
    return path


def _step_to_peec_inductance(step_path, R, wire_w, wire_h,
                             sigma=5.8e7, nwinc=1, step_size=None):
    """Full pipeline: STEP -> centerline -> filaments -> PEEC -> L."""
    from radia.coil_from_step import extract_centerline
    from radia.coil_from_cad import build_peec_from_path
    from radia.peec_topology import PEECCircuitSolver

    res = extract_centerline(
        step_path,
        start_hint=(np.array([R, 0, 0]), np.array([0, 1, 0])),
        step_size=step_size or (R * 2 * math.pi / 60),
        verbose=False,
    )

    n_stations = len(res.polyline)
    n_seg = n_stations - (0 if not res.closed else 1)

    # Build path (close the loop by appending first point)
    if res.closed:
        path = np.vstack([res.polyline, res.polyline[0:1]])
    else:
        path = res.polyline

    widths = np.full(path.shape[0] - 1, wire_w)
    heights = np.full(path.shape[0] - 1, wire_h)

    topo = build_peec_from_path(path, widths, heights,
                                sigma=sigma, nwinc=nwinc, nhinc=nwinc)

    # DC inductance (freq=0 -> only L, no R contribution)
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(0.0)
    # At DC, Z = R + j*0, but L = Im(Z)/omega is 0/0.
    # Instead extract L directly from the L matrix.
    L_matrix = topo["L"]
    n = L_matrix.shape[0]

    # Port impedance at low freq to get L
    f_low = 100.0  # 100 Hz
    Z_low = solver.compute_port_impedance(f_low)
    L = Z_low.imag / (2 * math.pi * f_low)

    return L, n_seg, topo


class TestStepToPeecInductance:
    """STEP -> PEEC inductance end-to-end verification."""

    def test_circular_torus(self, tmp_path):
        """Circular cross-section torus: R=30mm, a=3mm."""
        R, a = 0.030, 0.003
        step = str(tmp_path / "torus_circ.step")
        _build_torus_step(step, R, a)

        L_ana = _analytical_L_circular(R, a)
        L_peec, n_seg, _ = _step_to_peec_inductance(
            step, R, wire_w=2 * a, wire_h=2 * a)

        err = (L_peec - L_ana) / L_ana * 100
        print(f"\nCircular torus: R={R*1e3:.0f}mm, a={a*1e3:.1f}mm")
        print(f"  L_analytical = {L_ana*1e9:.2f} nH")
        print(f"  L_PEEC       = {L_peec*1e9:.2f} nH  ({n_seg} segments)")
        print(f"  Error        = {err:+.2f}%")
        assert abs(err) < 5, f"PEEC L error {err:.2f}% exceeds 5% threshold"

    def test_square_torus(self, tmp_path):
        """Square cross-section torus: R=30mm, w=6mm."""
        R, w = 0.030, 0.006
        step = str(tmp_path / "torus_sq.step")
        _build_square_torus_step(step, R, w)

        L_ana = _analytical_L_square(R, w)
        L_peec, n_seg, _ = _step_to_peec_inductance(
            step, R, wire_w=w, wire_h=w)

        err = (L_peec - L_ana) / L_ana * 100
        print(f"\nSquare torus: R={R*1e3:.0f}mm, w={w*1e3:.1f}mm")
        print(f"  L_analytical = {L_ana*1e9:.2f} nH")
        print(f"  L_PEEC       = {L_peec*1e9:.2f} nH  ({n_seg} segments)")
        print(f"  Error        = {err:+.2f}%")
        assert abs(err) < 5, f"PEEC L error {err:.2f}% exceeds 5% threshold"


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        from pathlib import Path
        t = TestStepToPeecInductance()
        t.test_circular_torus(Path(tmp))
        t.test_square_torus(Path(tmp))
        print("\nAll tests passed!")
