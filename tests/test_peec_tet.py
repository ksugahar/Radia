"""Native tetrahedral volume PEEC kernel (rad_peec_tet via radia.peec_tet).

Guards: the closed-form tetrahedron potentials of ``peec_matrices.tet_potential`` against a
high-order collapsed Gauss rule and Poisson's equation; the partial inductance and dense
circuit of ``tet_inductance`` / ``tet_kkt_solve`` on a coarse revolved copper ring against the
exact DC resistance of a rectangular ring and the axisymmetric ring kernel; symmetry and
positive definiteness of L; the HACApK operator ``TetPEECHMatrix`` and the preconditioned GMRES
circuit (``operator="hmatrix"``) against that dense route on a ring with low-rank blocks; loud
refusal of invalid input.
Assumes: Netgen meshes the OCC ring; the axisymmetric kernel is the independently validated route.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pm = pytest.importorskip("radia.peec_matrices")
ng = pytest.importorskip("ngsolve")
if not hasattr(pm, "tet_inductance"):
    pytest.fail("radia.peec_matrices lacks tet_inductance: rebuild the native module")

import radia._radia_pybind as rp  # noqa: E402
from radia.peec_axisym import solve_axisym_ring_peec  # noqa: E402
from radia.peec_tet import _face_basis, solve_tet_peec  # noqa: E402

if not hasattr(rp, "TetPEECHMatrix"):
    pytest.fail("radia._radia_pybind lacks TetPEECHMatrix: rebuild the native module")

TET = np.array([[0., 0, 0], [1.3, 0.1, 0], [0.2, 1.1, 0.1], [0.3, 0.4, 0.9]])
VOL = abs(np.linalg.det(TET[1:] - TET[0]))/6
R0, A, SIGMA = 10e-3, 2e-3, 5.8e7


def gauss_potential(x, n=40):
    g, w = np.polynomial.legendre.leggauss(n)
    t, w = (g + 1)/2, w/2
    U, V, S = np.meshgrid(t, t, t, indexing="ij")
    WU, WV, WS = np.meshgrid(w, w, w, indexing="ij")
    lam = np.stack([U, V*(1 - U), S*(1 - U)*(1 - V)], -1).reshape(-1, 3)
    ww = (WU*WV*WS*(1 - U)**2*(1 - V)).ravel()*6*VOL
    y = TET[0] + lam @ (TET[1:] - TET[0])
    r = np.linalg.norm(y - x, axis=1)
    return (ww/r).sum(), ((ww/r)[:, None]*y).sum(0)


@pytest.mark.parametrize("x", [[3.0, 2.0, 1.0], [1.5, 0.5, 0.4], [0.6, 0.5, -0.3]])
def test_potential_matches_gauss_quadrature_off_the_tetrahedron(x):
    x = np.array(x)
    p0, p1 = pm.tet_potential(x, TET)
    q0, q1 = gauss_potential(x)
    assert abs(p0/q0 - 1) < 1e-12, (p0, q0)
    assert np.abs(np.array(p1) - q1).max() < 1e-12*np.abs(q1).max(), (p1, q1)


@pytest.mark.parametrize("x, expected", [(TET.mean(0), -4*math.pi), (np.array([1.5, 0.5, 0.4]), 0.0)])
def test_potential_satisfies_poisson(x, expected):
    h = 1e-3
    f = lambda y: pm.tet_potential(y, TET)[0]
    lap = sum(f(x + h*e) + f(x - h*e) - 2*f(x) for e in np.eye(3))/h**2
    assert abs(lap - expected) < 1e-5, (lap, expected)


def test_potential_is_finite_on_vertices_and_edges():
    for x in (TET[2], 0.5*(TET[0] + TET[3]), TET[1:].mean(0)):
        p0, p1 = pm.tet_potential(x, TET)
        assert np.isfinite(p0) and np.all(np.isfinite(p1)), (x, p0, p1)


def ring_mesh(maxh):
    from netgen.occ import X, Y, Z, Axes, Axis, Glue, OCCGeometry, WorkPlane
    face = WorkPlane(Axes((R0 - A/2, 0, -A/2), n=-Y, h=X)).Rectangle(A, A).Face()
    shape = Glue([face.Revolve(Axis((0, 0, 0), Z), 180), face.Revolve(Axis((0, 0, 0), Z), -180)])
    for f in shape.faces:
        c = f.center
        f.name = ("cut" if c.x > 0 else "joint") if abs(c.y) < 1e-9 else "surface"
    return ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))


@pytest.fixture(scope="module")
def coarse_ring():
    mesh = ring_mesh(2e-3)
    with ng.TaskManager():
        out = solve_tet_peec(mesh, SIGMA, 100.0, "cut")
    return mesh, out


def test_coarse_ring_dc_resistance_is_the_exact_ring_resistance(coarse_ring):
    _, out = coarse_ring
    exact = 2*math.pi/(SIGMA*A*math.log((R0 + A/2)/(R0 - A/2)))   # azimuthal current, rectangular section
    assert abs(out["R_ohm"]/exact - 1) < 3e-3, (out["R_ohm"], exact)
    assert out["power_balance"] < 1e-10 and out["true_relative_residual"] < 1e-10, out


def test_coarse_ring_inductance_matches_the_axisymmetric_kernel(coarse_ring):
    _, out = coarse_ring
    n = 32
    xs, zs = np.linspace(R0 - A/2, R0 + A/2, n + 1), np.linspace(-A/2, A/2, n + 1)
    pts = np.array([(x, z) for z in zs for x in xs])
    idx = lambda i, j: j*(n + 1) + i
    tris = np.array([t for j in range(n) for i in range(n) for t in
                     ((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)), (idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)))],
                    dtype=np.int32)
    ref = solve_axisym_ring_peec(pts, tris, 100.0, SIGMA)
    assert abs(out["L_H"]/ref["L_H"] - 1) < 5e-3, (out["L_H"], ref["L_H"])


def interior_basis(mesh):
    fes = ng.HDiv(mesh, order=0, RT=True)
    pts, tets, order, coef = _face_basis(mesh, fes)
    counts = np.bincount(order.ravel(), minlength=fes.ndof)
    interior = np.flatnonzero(counts == 2)
    compact = -np.ones(fes.ndof, dtype=np.int64)
    compact[interior] = np.arange(len(interior))
    return pts, tets, coef, compact[order].astype(np.int32), len(interior)


def test_inductance_is_symmetric_positive_definite(coarse_ring):
    mesh, _ = coarse_ring
    pts, tets, coef, dof, n = interior_basis(mesh)
    L = pm.tet_inductance(pts, tets, coef, dof, n)
    assert np.array_equal(L, L.T)
    assert np.linalg.eigvalsh(L).min() > 0, np.linalg.eigvalsh(L).min()


@pytest.fixture(scope="module")
def fine_ring():
    return ring_mesh(1e-3)


def test_hmatrix_inductance_matches_the_dense_route(fine_ring):
    pts, tets, coef, dof, n = interior_basis(fine_ring)
    with ng.TaskManager():
        L = pm.tet_inductance(pts, tets, coef, dof, n)
        H = rp.TetPEECHMatrix(pts, tets, coef, dof, n)
    stats = H.stats()
    assert stats["n_lowrank"] > 0, stats                      # the ACA path is exercised
    rng = np.random.default_rng(3)
    x = rng.normal(size=n) + 1j*rng.normal(size=n)
    err = np.linalg.norm(np.asarray(H.matvec(x)) - L @ x)/np.linalg.norm(L @ x)
    assert err < 1e-7, (err, stats)
    for i, j in rng.integers(0, n, size=(50, 2)):
        assert abs(H.entry(int(i), int(j)) - L[i, j]) <= 1e-13*L[i, i], (i, j, H.entry(int(i), int(j)), L[i, j])
    assert np.allclose(np.asarray(H.diagonal()), np.diag(L), rtol=1e-13, atol=0.0)


def test_hmatrix_circuit_matches_the_dense_route(fine_ring):
    with ng.TaskManager():
        dense = solve_tet_peec(fine_ring, SIGMA, [100.0, 17.5e3], "cut")
        comp = solve_tet_peec(fine_ring, SIGMA, [100.0, 17.5e3], "cut", operator="hmatrix")
    for d, c in zip(dense, comp):
        assert c["frequency_hz"] == d["frequency_hz"] and c["assembly_s"] == comp[0]["assembly_s"], (c, d)
        assert c["iterations"] > 0 and c["true_relative_residual"] <= 1e-6, c
        dz = abs(c["impedance_ohm"]/d["impedance_ohm"] - 1)
        assert dz < 1e-7, (c["frequency_hz"], c["impedance_ohm"], d["impedance_ohm"])
        assert c["power_balance"] < 1e-8, c
    assert dense[1]["R_ohm"] > 1.02*dense[0]["R_ohm"], (dense[0]["R_ohm"], dense[1]["R_ohm"])   # skin effect


def test_invalid_input_is_refused():
    tets = np.array([[0, 1, 2, 3]], dtype=np.int32)
    coef, dof = np.ones((1, 4)), np.array([[0, 1, 2, 3]], dtype=np.int32)
    flat = TET.copy()
    flat[3] = flat[:3].mean(0)
    with pytest.raises(ValueError, match="degenerate"):
        pm.tet_inductance(flat, tets, coef, dof, 4)
    with pytest.raises(ValueError, match="references point"):
        pm.tet_inductance(TET, np.array([[0, 1, 2, 4]], dtype=np.int32), coef, dof, 4)
    with pytest.raises(ValueError, match="out of range"):
        pm.tet_inductance(TET, tets, coef, np.array([[0, 1, 2, 9]], dtype=np.int32), 4)
    with pytest.raises(ValueError, match="cut surface"):
        solve_tet_peec(ring_mesh(4e-3), SIGMA, 100.0, "missing")
    with pytest.raises(ValueError, match="operator must be"):
        solve_tet_peec(ring_mesh(4e-3), SIGMA, 100.0, "cut", operator="sparse")
