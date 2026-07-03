"""Golden: the 2D PLANAR HDiv-VIM charge Gram (motor cross-section layer).

`radia.vim._vim.build_charge_gram(HDiv(mesh2d, order=1))` must route to the C++ dim2 log-kernel Gram
(-ln(r)/2pi; charges = -div M on tri/quad cells + M.n on boundary edges, Piola-exact REF measures) and
produce a demag operator N = B^T G B whose generalized spectrum vs the HDiv mass respects the 2D physical
bound [0, 1], with the CLOSED-FORM anchors:
  * DISK demag factor = 1/2 EXACT (uniform in-plane M);
  * ELLIPSE 2:1 -> (1/3, 2/3);
  * linear solve -> 2D Clausius-Mossotti  M/H0 = chi/(1 + chi/2).

This also locks the two quadrature lessons that produced eig > 1 during bring-up (2026-07-03):
  * the OUTER must be the PRODUCT-GAUSS rule -- Dunavant-7's coherent misintegration of the potential's
    C1 kinks, amplified by the div-scaled charge map, leaked the QUAD-mesh spectrum to 1.072 while every
    individual entry agreed to ~3e-5;
  * the EDGE inner must split-grade at the kernel-peak PARAMETER (the self entry's analytic value
    (3/2 - ln L)/2pi is reproduced to 6 digits) -- endpoint grading overestimated edge self-energies.
Self-contained (unit_square / MakeStructured2DMesh / OCC disk); ~20 s.
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from ngsolve.meshes import MakeStructured2DMesh  # noqa: E402
from netgen.occ import WorkPlane, OCCGeometry  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from radia.vim._vim import build_charge_gram  # noqa: E402


def _materialize_N(B, G):
    n = B.shape[1]
    N = np.zeros((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        N[:, j] = B.T @ np.asarray(G.matvec_sym((B @ e).tolist()), float)
        e[j] = 0.0
    return 0.5 * (N + N.T)


def _gate(mesh):
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M_mass = build_charge_gram(fes)
        N = _materialize_N(B, G)
        Md = M_mass.toarray()
        w = sla.eigh(N, Md, eigvals_only=True)
        gf = ng.GridFunction(fes)
        gf.Set(ng.CoefficientFunction((1.0, 0.0)))
        mux = gf.vec.FV().NumPy().copy()
        gf.Set(ng.CoefficientFunction((0.0, 1.0)))
        muy = gf.vec.FV().NumPy().copy()
    chk = G.hex_state_check()
    assert chk["ctor"] == chk["now"], f"2D Gram state canary mismatch {chk}"
    Dx = float(mux @ (N @ mux)) / float(mux @ (Md @ mux))
    Dy = float(muy @ (N @ muy)) / float(muy @ (Md @ muy))
    return w, Dx, Dy, (N, Md, mux)


@pytest.mark.parametrize(
    "name,make",
    [
        ("tri-square", lambda: ng.Mesh(ng.unit_square.GenerateMesh(maxh=0.35))),
        ("quad-square", lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3)),
        ("quad-distorted", lambda: MakeStructured2DMesh(
            quads=True, nx=3, ny=3,
            mapping=lambda x, y: (x + 0.15*x*y*(1 - x), y + 0.12*x*y*(1 - y)))),
    ],
)
def test_2d_spectrum_in_bound_and_square_demag(name, make):
    w, Dx, Dy, _ = _gate(make())
    assert w.min() > -1e-8, f"{name}: N not PSD (min eig {w.min():.2e})"
    assert w.max() < 1.005, f"{name}: demag spectrum escaped [0,1] (max eig {w.max():.6f})"
    assert w.max() > 0.99, f"{name}: max eig {w.max():.6f} unexpectedly low (quadrature regression?)"
    # unit square: Dx + Dy = 1 by symmetry -> each 1/2
    assert abs(Dx - 0.5) < 5e-3 and abs(Dy - 0.5) < 5e-3, f"{name}: square demag ({Dx:.4f},{Dy:.4f}) off 1/2"


def test_2d_disk_demag_half_and_clausius_mossotti():
    wp = WorkPlane().Circle(0, 0, 1.0).Face()
    mesh = ng.Mesh(OCCGeometry(wp, dim=2).GenerateMesh(maxh=0.3))
    w, Dx, Dy, (N, Md, mux) = _gate(mesh)
    assert w.max() < 1.005 and w.min() > -1e-8
    # 2D transverse demag factor of a disk is EXACTLY 1/2
    assert abs(Dx - 0.5) < 2e-3 and abs(Dy - 0.5) < 2e-3, f"disk demag ({Dx:.5f},{Dy:.5f}) off 1/2"
    # linear demag solve vs 2D Clausius-Mossotti M/H0 = chi/(1 + chi/2)
    for chi in (10.0, 999.0):
        A = Md/chi + N
        m = np.linalg.solve(A, Md @ mux)
        Mavg = float(mux @ (Md @ m)) / float(mux @ (Md @ mux))
        ref = chi/(1.0 + chi/2.0)
        assert abs(Mavg - ref)/ref < 3e-3, f"chi={chi}: M/H0={Mavg:.6f} vs CM2D {ref:.6f}"


def test_2d_ellipse_demag_thirds():
    wp = WorkPlane().Ellipse(2.0, 1.0).Face()
    mesh = ng.Mesh(OCCGeometry(wp, dim=2).GenerateMesh(maxh=0.35))
    w, Dx, Dy, _ = _gate(mesh)
    assert w.max() < 1.005 and w.min() > -1e-8
    # ellipse a:b -> N_x = b/(a+b) = 1/3, N_y = 2/3
    assert abs(Dx - 1/3) < 5e-3, f"ellipse Dx={Dx:.5f} off 1/3"
    assert abs(Dy - 2/3) < 5e-3, f"ellipse Dy={Dy:.5f} off 2/3"
