"""Stage-1 adjoint gate at the research-mesh configuration (golden bands).

Unit ball, maxh 0.35 (~270 tets), log-uniform s in [1e-2, 1]: the adjoint
gradient must match central finite differences to better than 1e-7 relative
along a random direction (measured 8.1e-10 at h=1e-5 on 2026-07-28), and the
dipole-reciprocity load must reproduce the independent C++ analytic charge
evaluator to better than 1e-7 relative (measured 1.1e-10 at bonus_intorder=10).
"""
from math import pi
from types import SimpleNamespace

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("radia")

from netgen.occ import OCCGeometry, Pnt, Sphere  # noqa: E402
from ngsolve import HDiv, InnerProduct, Mesh, SetNumThreads, TaskManager  # noqa: E402

from radia.isochronous_topopt import (  # noqa: E402
    MU0, DensityAdjointVIM, demag_field_from_solution, field_functional_load,
    gradient_pair_points, orbit_arc_points, uniform_field_load,
)


@pytest.fixture(scope="module")
def gate():
    SetNumThreads(4)
    with TaskManager():
        mesh = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.35))
        fes = HDiv(mesh, order=1)
        prob = DensityAdjointVIM(fes, eps=1e-7)
        rng = np.random.default_rng(1234)
        s0 = np.exp(rng.uniform(np.log(1e-2), 0.0, prob.n_el))
        orbit, radial = orbit_arc_points(1.6, 0.3, 6)
        w_orbit = np.array([1.0, -0.7, 1.3, 0.5, -1.1, 0.9])
        pts, wts = gradient_pair_points(orbit, w_orbit, delta=0.1, axis=0)
        f_state = uniform_field_load(fes, (0.0, 0.0, 1.0))
        f_adj = field_functional_load(fes, pts, wts, axis=2, scale=MU0,
                                      bonus_intorder=10)
        base = prob.objective_and_gradient(s0, f_state, f_adj, tol=1e-12)
    assert prob.n_el > 200, "research-mesh class expected"
    return SimpleNamespace(mesh=mesh, fes=fes, prob=prob, s0=s0,
                           dipole_points=pts, dipole_weights=wts,
                           f_state=f_state, f_adj=f_adj, base=base)


def test_directional_adjoint_matches_fd(gate):
    d = np.random.default_rng(99).standard_normal(gate.prob.n_el)
    ds = gate.s0 * d
    predicted = float(gate.base.gradient @ ds)
    h = 1e-5
    with TaskManager():
        gp, _ = gate.prob.solve(gate.s0 + h * ds, gate.f_state, tol=1e-12)
        gm, _ = gate.prob.solve(gate.s0 - h * ds, gate.f_state, tol=1e-12)
    fd = (float(InnerProduct(gate.f_adj.vec, gp.vec))
          - float(InnerProduct(gate.f_adj.vec, gm.vec))) / (2.0 * h)
    rel = abs(fd - predicted) / abs(predicted)
    assert rel < 1e-7, (fd, predicted, rel)


def test_reciprocity_against_cpp_evaluator(gate):
    with TaskManager():
        H_d = demag_field_from_solution(gate.prob.demag, gate.base.gfM,
                                        gate.dipole_points)
    independent = float(MU0 * np.dot(gate.dipole_weights, H_d[:, 2]))
    rel = abs(gate.base.objective - independent) / abs(independent)
    assert rel < 1e-7, (gate.base.objective, independent, rel)
