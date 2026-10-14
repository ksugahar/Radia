r"""Cylindrical conductor Joule heating (electro-thermal) -- regression test (#41).

A round conductor (radius a) with a uniform Joule source q=J^2/sigma has the parabolic
radial temperature rise dT(r)=q(a^2-r^2)/(4k), centre dT_peak=q a^2/(4k). radia
`solve_heat_steady` on a disk reproduces it; the closed form is the tool-independent gate
(the cylindrical sibling of the slab electro-thermal #19). A reference commercial FE code
matches the same closed form to 0.000% live (recorded internally, not here)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import joule_heated_cylinder_peak_dT

A, K = 0.02, 2.0
Q = 5.0e5          # uniform Joule density [W/m^3] (= J^2/sigma)


def test_peak_dt_formula():
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, K), Q * A * A / (4 * K), rel_tol=1e-12)
    # dT_peak = 25 K for these values; scales with q and a^2, inverse with k
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, K), 25.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(2 * Q, A, K), 50.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, 2 * A, K), 100.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, 4 * K), 6.25, rel_tol=1e-12)
