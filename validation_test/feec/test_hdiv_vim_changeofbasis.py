"""Golden: the HDiv-VIM charge-map change-of-basis lands the monomial charge coefficients in the GRAM's
cell_verts geometry frame (NOT NGSolve's CalcShape frame).

This locks the 2026-06-13 fix for a latent frame-mismatch bug: NGSolve's L2/SurfaceL2 `CalcShape` uses its own
reference-element vertex ordering (ref(0,0,0)->the LAST mesh vertex), but the C++ charge-Gram interprets the
monomials via `cell_verts` in MESH-VERTEX order (ref(0,0,0)->V0).  The old change-of-basis evaluated the
monomial and CalcShape at the SAME reference point, so it produced coefficients in NGSolve's frame -- a fixed
vertex permutation off from the Gram's frame.  This is INVISIBLE to every uniform-M / demag-factor test
(uniform M has div M = 0 => zero volume charge, and a constant charge is frame-invariant), but it scrambles the
VOLUME charge of any non-uniform (high-order) magnetisation, corrupting high-order demag solves.

The fix evaluates the monomial at the cell_verts-frame coordinate g(pt) that maps (via GetTrafo) to the same
physical point as the NGSolve-ref point pt.  Correctness <=> the change-of-basis exactly REPRODUCES the L2 /
SurfaceL2 function in the cell_verts frame.  By the divergence-theorem identity this also makes the discrete
charge carry the exact dipole (sum c_k INT z*m_k == INT M_z), the moment the old code dropped.

NGSolve + Netgen required.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Sphere, OCCGeometry, Pnt  # noqa: E402

from radia.vim._vim import _change_of_basis, _tet_ref, _tri_ref, _monos_vol, _monos_surf  # noqa: E402


def _mesh(maxh=1.5):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


def test_volume_change_of_basis_reproduces_L2_in_gram_frame():
    """Sv reconstructs a random L2(1) field at cell_verts-frame points to machine precision (was ~O(1) wrong in
    the old NGSolve-frame code for any field with non-trivial spatial variation)."""
    mesh = _mesh()
    order = 2
    mons = _monos_vol(max(order - 1, 0))
    with ng.TaskManager():
        L2v = ng.L2(mesh, order=max(order - 1, 0))
        gf = ng.GridFunction(L2v)
        rng = np.random.default_rng(3)
        gf.vec.FV().NumPy()[:] = rng.standard_normal(L2v.ndof)
        e0 = ng.ElementId(ng.VOL, 0)
        V = np.array([mesh[v].point for v in mesh[e0].vertices])          # cell_verts (mesh order)
        Sv = _change_of_basis(L2v.GetFE(e0), mons, *_tet_ref(max(3 * order, 4)), dim=3,
                              trafo=mesh.GetTrafo(e0), Vmesh=V)
        mc = Sv @ np.array(gf.vec)[list(L2v.GetDofNrs(e0))]
        Jm = np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]).T
        worst = 0.0
        for g in [(0.2, 0.2, 0.2), (0.5, 0.1, 0.1), (0.1, 0.5, 0.2), (0.1, 0.1, 0.6)]:
            P = V[0] + Jm @ np.array(g)                                   # physical = cell_verts map of g
            recon = sum(mc[a] * g[0] ** i * g[1] ** j * g[2] ** k for a, (i, j, k) in enumerate(mons))
            worst = max(worst, abs(recon - gf(mesh(*P))))
    assert worst < 1e-9, f"volume change-of-basis not in the Gram (cell_verts) frame: worst recon err {worst:.2e}"


def test_surface_change_of_basis_reproduces_SurfaceL2_in_gram_frame():
    """Ss reproduces a random SurfaceL2(2) field's moments against physical test CFs (1,x,y,z,x^2,xy,...) to
    machine precision.  (SurfaceL2 GridFunction point-eval returns 0, so we check moments, not point values.)"""
    mesh = _mesh()
    order = 2
    mons = _monos_surf(order)
    with ng.TaskManager():
        L2b = ng.SurfaceL2(mesh, order=order)
        gf = ng.GridFunction(L2b)
        rng = np.random.default_rng(5)
        gf.vec.FV().NumPy()[:] = rng.standard_normal(L2b.ndof)
        tests = {"1": ng.CF(1), "x": ng.x, "y": ng.y, "z": ng.z,
                 "x2": ng.x ** 2, "z2": ng.z ** 2, "xy": ng.x * ng.y, "yz": ng.y * ng.z}
        gfm = {k: float(ng.Integrate(gf * t, mesh.Boundaries(".*"), element_wise=True).NumPy()[0]) for k, t in tests.items()}
        f0 = ng.ElementId(ng.BND, 0)
        V = np.array([mesh[v].point for v in mesh[f0].vertices])
        Ss = _change_of_basis(L2b.GetFE(f0), mons, *_tri_ref(max(3 * order, 4)), dim=2,
                              trafo=mesh.GetTrafo(f0), Vmesh=V)
        mc = Ss @ np.array(gf.vec)[list(L2b.GetDofNrs(f0))]
    J2 = np.array([V[1] - V[0], V[2] - V[0]]).T
    twoA = np.linalg.norm(np.cross(J2[:, 0], J2[:, 1]))
    Ptri, Wtri = _tri_ref(10)
    xyz = V[0][None, :] + (J2 @ Ptri.T).T

    def tval(name):
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        return {"1": np.ones_like(x), "x": x, "y": y, "z": z, "x2": x * x, "z2": z * z, "xy": x * y, "yz": y * z}[name]

    worst = 0.0
    for name in tests:
        recon = sum(mc[a] * twoA * np.sum(Wtri * (Ptri[:, 0] ** i * Ptri[:, 1] ** j) * tval(name))
                    for a, (i, j) in enumerate(mons))
        worst = max(worst, abs(recon - gfm[name]))
    assert worst < 1e-9, f"surface change-of-basis not in the Gram (cell_verts) frame: worst moment err {worst:.2e}"


def test_demag_factor_unaffected_by_frame_fix():
    """The frame fix must NOT change the (frame-invariant) demag factor for uniform M -- still ~1/3."""
    from radia.vim import DemagOperator
    mesh = _mesh()
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        N = DemagOperator(fes)
        D = N.DemagFactor(ng.CF((0, 0, 1)))
    assert 0.31 < D < 0.345, f"demag factor {D:.5f} not ~1/3 after the change-of-basis frame fix"
