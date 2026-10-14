r"""Demagnetizing factor & uniformly-magnetized-body internal field (#71) -- test.

N tuples sum to 1 (ellipsoid rule); B_in = Br(1-N) gives 2Br/3 (sphere) / Br/2 (transverse cylinder)
/ 0 (thin slab); H_in = -N Br/mu0. The transverse-cylinder B_in=Br/2 is confirmed by an NGSolve
solve within the open-domain truncation. Tool-independent."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (demagnetizing_factor, magnetized_body_internal_field,
                                           demagnetizing_field, MU0)


def test_factors_sum_to_one_and_values():
    for shape in ("sphere", "cylinder_transverse", "cylinder_axial", "thin_slab"):
        N = demagnetizing_factor(shape)
        assert math.isclose(sum(N), 1.0, rel_tol=1e-12)
    assert demagnetizing_factor("sphere") == (1 / 3, 1 / 3, 1 / 3)
    assert demagnetizing_factor("cylinder_transverse") == (0.5, 0.5, 0.0)
    assert demagnetizing_factor("thin_slab") == (1.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        demagnetizing_factor("banana")


def test_internal_field_and_demag_field():
    Br = 1.2
    # B_in = Br(1-N): sphere 2/3, transverse cyl 1/2, slab 0
    assert math.isclose(magnetized_body_internal_field(Br, 1 / 3), 2 * Br / 3, rel_tol=1e-12)
    assert math.isclose(magnetized_body_internal_field(Br, 0.5), Br / 2, rel_tol=1e-12)
    assert math.isclose(magnetized_body_internal_field(Br, 1.0), 0.0, abs_tol=1e-15)
    # axial cylinder (N=0) holds the full remanence
    assert math.isclose(magnetized_body_internal_field(Br, 0.0), Br, rel_tol=1e-12)
    # demag field H_in = -N Br/mu0 (negative, opposes M)
    assert math.isclose(demagnetizing_field(Br, 0.5), -0.5 * Br / MU0, rel_tol=1e-12)
    assert demagnetizing_field(Br, 1 / 3) < 0.0
