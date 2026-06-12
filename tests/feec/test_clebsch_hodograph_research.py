"""Golden lock for the Clebsch hodograph RESEARCH examples
(examples/clebsch_hodograph/): the A-method (vector-potential primary) dual
and the Kelvin (exact open boundary) variant.  Imports each script's
solve() and asserts the headline accuracy + bidirectional consistency bands.

These are research demonstrations (not panel modes); the verified FORWARD
panel mode is locked separately by tests/panels/test_clebsch_golden.py.

Verified 2026-06-12:
  A-method 2-D   field_error ~3e-5  consistency ~4e-4
  Kelvin axisym  field_error ~1e-7  consistency ~2e-5  (exact open boundary)
"""
import sys
from pathlib import Path

EXDIR = Path(__file__).resolve().parents[2] / "examples" / "clebsch_hodograph"
sys.path.insert(0, str(EXDIR))


def test_a_method_clebsch_2d():
    import a_method_clebsch_2d as am
    r = am.solve(mu_r=1000.0, order=3, maxh=0.10)   # coarser for CI speed
    # Vector-potential A-method recovers the 2-D interior field 2 mu_r/(mu_r+1).
    assert r["field_error"] < 2e-3, r
    # A_z / V accuracy vs the exact conjugate pair.
    assert r["Az_error"] < 2e-3 and r["V_error"] < 2e-3, r
    # Hodograph self-consistency B(from A_z) vs B(from V).
    assert 0.0 < r["consistency"] < 2e-3, r


def test_hodograph_kelvin_axisym():
    import hodograph_kelvin_axisym as hk
    r = hk.solve(mu_r=100.0, order=3, maxh=0.05)    # coarser for CI speed
    # Kelvin open boundary is EXACT -> field_error must be tiny (<< the
    # ~3e-3 of the far-truncated panel-mode sphere).
    assert r["field_error"] < 1e-3, r
    assert 0.0 < r["consistency"] < 5e-3, r


def test_cohomology_currentlink():
    """The hodograph's scalar coordinate is a 1st-cohomology class iff a
    current threads a hole (period != 0).  Locks the radia.cohomology
    generator + the grad-fit obstruction contrast."""
    import cohomology_hodograph_currentlink as ch
    r = ch.solve(maxh=0.025)
    assert r["b1_solid"] == 0 and r["b1_washer"] == 1, r
    assert r["curl_rel"] < 1e-6, r                       # generator is curl-free
    assert abs(abs(r["oint_hole"]) - 1.0) < 0.05, r      # unit circulation
    assert abs(r["oint_contractible"]) < 1e-2, r
    # current-linking field is NOT a gradient -> cohomology required.
    assert r["residual_cohomology_field"] > 0.5, r
    # zero-period field IS a gradient -> single-valued scalar coordinate.
    assert r["residual_gradient_field"] < 1e-3, r
