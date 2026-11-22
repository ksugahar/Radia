"""Device surrogate of the linear volumetric F1 slab against its stored baseline.

Guards: breakage of the device route from mesh arrays to Joule power, at p=1 with complex AMS
and p=2 with BDDC over AMS, using the same model, integration and power evaluation as
``test_ih_ams_f1_surrogate``. Device matrices, the Dirichlet lift, the COCR iteration and the true
residual gate all sit on this path. The ledger must name only permitted payload families.
Assumes: surrogate of adoption record induction_heating/20261024_ams_linear_slab (level m=2). The
baseline is the CPU route's fixture and tolerance, unchanged: user approval
APPROVAL-20261024-IH-AMS-F1-SURROGATE-TOL-01, rtol=1e-8, atol=0. Optional GPU lane: skips visibly
unless RADIA_FEM_NATIVE_TEST_DIR is set. No accuracy gate belongs here.
"""
import cmath
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")


@pytest.fixture(scope="module")
def route():
    if not os.environ.get("RADIA_FEM_NATIVE_TEST_DIR"):
        pytest.skip("Optional GPU lane requires RADIA_FEM_NATIVE_TEST_DIR")
    import radia.fem_gpu as module
    return module


def _compute(route, p):
    from ngsolve.meshes import MakeStructured3DMesh
    mu, sigma, delta, H0 = 4 * math.pi * 1e-7, 5e6, 1e-3, 1000.0
    omega = 2 / (mu * sigma * delta ** 2)
    a, width, k = 2 * delta, delta / 4, (1 + 1j) / delta
    Ay = mu * H0 * (ng.exp(k * ng.x) - ng.exp(-k * ng.x)) / (2 * k * cmath.cosh(k * a))
    boundary_A = ng.CoefficientFunction((0, Ay, 0))
    mesh = MakeStructured3DMesh(hexes=False, nx=32, ny=2, nz=2, secondorder=False,
                                mapping=lambda x, y, z: (4 * delta * x - a, width * (y - .5), width * (z - .5)))
    space = (ng.HCurl(mesh, order=1, nograds=True, complex=True, dirichlet=".*") if p == 1
             else ng.HCurl(mesh, order=2, complex=True, dirichlet=".*"))
    coefficients = route.Coefficients(stiff=1 / mu, mass_im=omega * sigma, stiff_aux=1 / mu, mass_aux=omega * sigma)
    g = ng.GridFunction(space)
    with ng.TaskManager():
        g.Set(boundary_A, definedon=mesh.Boundaries(".*"))
    lift = np.array(g.vec.FV().NumPy(), copy=True)
    result = route.solve_hcurl(space, coefficients, np.zeros(space.ndof, dtype=complex), Backend="gpu", lift=lift,
                               bonus=(2 * p + 2, 2 * p + 2, 0), tolerance=1e-11, cycles=4)
    free = np.array(list(space.FreeDofs()), dtype=bool)
    assert np.array_equal(result["solution"][~free], lift[~free]), "the Dirichlet lift changed"
    g.vec.FV().NumPy()[:] = result["solution"]
    J = -1j * omega * sigma * g
    with ng.TaskManager():
        power = float(ng.Integrate(ng.InnerProduct(J, J).real, mesh, order=2 * p + 12)) / (2 * sigma)
    return power, result


@pytest.mark.parametrize("p", [1, 2])
def test_device_f1_joule_power(route, p):
    power, result = _compute(route, p)
    fixture = json.loads((Path(__file__).parent / "fixtures" / f"ih_f1_m2_p{p}_5_3_3.json").read_text(encoding="utf-8"))
    quantity = fixture["quantities"]["power_W"]
    assert quantity["rtol"] == 1e-8 and quantity["atol"] == 0
    bound = quantity["rtol"] * abs(quantity["expected"])
    assert abs(power - quantity["expected"]) <= bound, (
        f"power {power} W, expected {quantity['expected']} W, bound {bound} W, iterations {result['iterations']}")
    assert result["true_residual"] <= 1e-6, f"true residual {result['true_residual']}"
    assert result["backend"] == "gpu" and result["preconditioner"] == ("ams" if p == 1 else "bddc")
