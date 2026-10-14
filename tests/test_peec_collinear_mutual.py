"""Mutual inductance must stay finite for touching collinear filaments.

A tapered sweep puts consecutive lane segments on one straight line, and
the parallel-filament kernel evaluated

    F(x, d) = x * log((x + sqrt(x**2 + d**2)) / d) - sqrt(x**2 + d**2)

which is algebraically x * arsinh(x/d) - sqrt(x**2 + d**2) but cancels
catastrophically for x < 0.  Once the perpendicular separation drops
below |x| * sqrt(2 * eps) the sum rounds to exactly 0, the log returns
-inf and the mutual inductance comes back +inf.  The kernel's own guard
tested d < 1e-15 m, five orders of magnitude under the onset for a 10 mm
filament, so it never fired.

Measured 2026-09-20 on a tapered beak fin: 2 of 704**2 entries of L were
non-finite, and nothing downstream looked -- the run finished and put a
NaN inductance in the result.  Both halves are pinned here: the kernel
must return a finite mutual, and the solver must refuse a matrix that is
not finite.
"""

from __future__ import annotations

import numpy as np
import pytest

# Two touching, collinear 10 mm filaments lifted straight out of the
# tapered fixture's lane graph.  Full precision matters: rounding these
# coordinates perturbs the geometry enough to hide the singularity.
SEG_A0 = (-0.0006303055161212521, 0.002, 0.019999999966666668)
SEG_A1 = (-0.0009281921642270711, 0.002, 0.03)
SEG_B1 = (-0.0012260784353539985, 0.002, 0.04000000003333334)
WIDTH_A = 0.00044624231366741324
WIDTH_B = 0.00041994324112576935
HEIGHT = 0.00017063200947021692
PERPENDICULAR_SEPARATION_M = 1.8840587215840723e-10


def _two_segment_topology():
    from radia.peec_matrices import PEECBuilder

    builder = PEECBuilder()
    n0 = builder.add_node_at(*SEG_A0)
    n1 = builder.add_node_at(*SEG_A1)
    n2 = builder.add_node_at(*SEG_B1)
    builder.add_connected_segment(n0, n1, WIDTH_A, HEIGHT, sigma=5.8e7)
    builder.add_connected_segment(n1, n2, WIDTH_B, HEIGHT, sigma=5.8e7)
    builder.add_port(n0, n2)
    return builder.build_topology()


def test_geometry_still_sits_below_the_cancellation_onset():
    """Guard the guard: the fixture is only a test if it is degenerate."""
    a0, a1, b1 = (np.array(p) for p in (SEG_A0, SEG_A1, SEG_B1))
    di, dj = a1 - a0, b1 - a1
    length = float(np.linalg.norm(di))
    axis = di / length
    separation = float(np.linalg.norm((a1 + b1) / 2 - (a0 + a1) / 2
                                      - axis * (((a1 + b1) / 2
                                                 - (a0 + a1) / 2) @ axis)))
    onset = length * np.sqrt(2.0 * np.finfo(float).eps)
    assert separation == pytest.approx(PERPENDICULAR_SEPARATION_M, rel=1e-6)
    assert separation < onset, (
        f"separation {separation:.3e} m is no longer under the {onset:.3e} m "
        f"cancellation onset, so this fixture no longer tests anything")
    assert float(di @ dj / (length * np.linalg.norm(dj))) > 1 - 1e-14


def test_collinear_touching_filaments_have_a_finite_mutual():
    topology = _two_segment_topology()
    inductance = np.asarray(topology["L"])

    assert np.all(np.isfinite(inductance)), (
        f"mutual inductance is not finite: {inductance.tolist()}")
    # Neighbouring, slightly less degenerate configurations of the same
    # pair land at 1.375e-09 .. 1.387e-09 H.
    assert inductance[0, 1] == pytest.approx(1.3869e-09, rel=5e-3)
    assert inductance[0, 1] > 0.0


def test_solver_refuses_a_non_finite_partial_element_matrix():
    from radia.peec_topology import PEECCircuitSolver

    topology = _two_segment_topology()
    poisoned = np.array(topology["L"], dtype=float)
    poisoned[0, 1] = np.inf
    poisoned[1, 0] = np.inf
    topology = dict(topology)
    topology["L"] = poisoned

    with pytest.raises(ValueError, match="non-finite"):
        PEECCircuitSolver(topology)
