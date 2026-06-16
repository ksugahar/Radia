"""M0-linear parity gate (productionization, docs/hdiv_vim/PRODUCTIONIZATION.md): the HDiv-VIM
production entry `radia.vim.hdiv_demag_solve` vs the shipped yano-type MSC (`rad.Solve`) on the SAME
soft-iron body -- a cube under a uniform applied field.

A cube's volume-averaged demag factor is exactly 1/3 (cubic symmetry, Dx=Dy=Dz, sum=1), but its M is
NON-uniform, so the ellipsoid formula chi/(1+chi/3)H is NOT exact for the volume-average M; the right
parity metric is the two solvers measured head-to-head.

Measured convergence (2026-06-15, H0=1000 A/m, mu_r=100):
    yano (rad.Solve):  n=4 3348.2, n=6 3410.3, n=8 3437.7, n=10 3452.9   (hexes; rising -> ~3465)
    hdiv (VIM):        L/4 3450.5, L/5 3462.7, L/6 3470.4, L/8 3479.2     (tets;  rising -> ~3490)
    finest yano vs finest hdiv = 0.76%.  HDiv +N GMRES iters mesh-robust: 24, 28, 26, 37.

The gate locks (1) yano == HDiv within 3% at a moderate resolution, (2) HDiv demag factor ~ 1/3, and
(3) HDiv mesh-robustness (the M_mass^-1 preconditioner keeps the +N iteration count bounded as the
mesh refines -- a Jacobi diagonal blew past 5000 iters at 7224 faces, the bug this gate guards).

This is a transitional gate: it compares against yano-type while it still ships; it retires when
yano-type is sealed (M5).
"""
import math
import warnings

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402
from radia.netgen_mesh_import import create_hex_mesh_grid  # noqa: E402

MU0 = 4.0e-7 * math.pi
H0 = 1000.0     # uniform applied field, +z (A/m)
L = 0.02        # cube edge (m)
MU_R = 100.0


def _yano_cube_Mz(n):
    """Volume-average M_z of the soft-iron cube via the yano-type MSC rad.Solve (n x n x n hexes,
    equal volume -> mean of per-element M_z).  Applied field via ObjBckg B = mu0 H0 (the free-space
    source field whose H is H0)."""
    rad.UtiDelAll()
    core = create_hex_mesh_grid(center=[0, 0, 0], size=[L, L, L], divisions=[n, n, n],
                                mu_r=MU_R, verbose=False)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    cont = rad.ObjCnt([core, bkg])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)   # yano backend is deprecated; expected here
        rad.Solve(cont, 1e-6, 3000, 0)                        # LU
    res = rad.ObjM(core)
    return float(np.mean([m[2] for (_c, m) in res]))


def _hdiv_cube(maxh):
    box = Box((-L / 2, -L / 2, -L / 2), (L / 2, L / 2, L / 2))
    geo = OCCGeometry(box)
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
        return hdiv_demag_solve(mesh, MU_R, ng.CoefficientFunction((0, 0, H0)))


def test_yano_hdiv_cube_linear_parity():
    """yano-type MSC and HDiv-VIM agree within 3% on the soft-iron cube (volume-average M_z)."""
    mz_yano = _yano_cube_Mz(8)                       # 512 hexes
    res = _hdiv_cube(L / 6)                           # 860 tets
    mz_hdiv = res["M_avg"][2]
    rel = abs(mz_yano - mz_hdiv) / abs(mz_hdiv)
    assert rel < 0.03, f"yano {mz_yano:.1f} vs HDiv {mz_hdiv:.1f} parity rel {rel:.2e}"
    assert abs(res["demag"] - 1.0 / 3.0) < 5e-3, f"HDiv demag factor off: {res['demag']}"
    # transverse average ~ 0 (uniform +z drive)
    assert abs(res["M_avg"][0]) < 1e-2 * mz_hdiv and abs(res["M_avg"][1]) < 1e-2 * mz_hdiv


def test_hdiv_plus_N_mesh_robust():
    """The M_mass^-1-preconditioned +N GMRES iteration count stays BOUNDED as the cube mesh refines
    (the bug this guards: a Jacobi-diagonal preconditioner blew past 5000 iters at 7224 faces)."""
    iters = {}
    for maxh in (L / 4, L / 8):
        res = _hdiv_cube(maxh)
        iters[maxh] = res["iters"]
    assert max(iters.values()) < 100, f"HDiv +N GMRES not mesh-robust: {iters}"
    # refining must not blow up the count (mesh-robustness, not just a one-off bound)
    assert iters[L / 8] <= 3 * iters[L / 4] + 20, f"iters grow too fast with refinement: {iters}"
