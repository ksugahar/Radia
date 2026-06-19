"""IMA (Image Method of Analysis / mirror symmetry) for the FEEC HDiv-VIM demag.

The demag operator is an INTEGRAL operator (charge-Gram), so mirror symmetry enters as IMAGE CHARGES
folded into the Gram (not an FES symmetry BC):

    G_IMA(a,b) = G(a,b) + sum_{S != empty} sign(S) * G(a, R_S(b))

where R_S reflects charge b's geometry on the mirror axes in subset S, and sign(S) = product of the
per-plane image-string signs (parse_image_string / image_group).  The SAME per-image sign applies to
BOTH cell (rho = -div M) and face (sigma = M.n) charges -- pinned empirically below.  This lets the
REDUCED (1/2, 1/4, 1/8) mesh reproduce the FULL model.

Locks (high mu_r -> non-uniform M -> exercises rho AND sigma signs):
  (1) half (1 plane), quarter (2 planes), octant (3 planes) of a cube + IMA reproduce the FULL cube
      M_avg and demag to ~1e-3 (operator to machine precision; M via the iterative solve);
  (2) a NONLINEAR quarter + IMA matches the full nonlinear cube (rho path under saturation);
  (3) parse_image_string / image_group give the documented (axis, sign) + product-sign group;
  (4) hdiv_demag_solve(image=..., scalable=True) is a fail-loud conflict (dense-only for now).
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
import radia.vim as vim  # noqa: E402
from radia.vim import _core as core  # noqa: E402

A = 0.01
MU_R = 200.0
H0 = 1000.0
MU0 = 4e-7 * math.pi


def _box(x0, x1, y0, y1, z0, z1, nx, ny, nz):
    mp = lambda x, y, z: (x0 + (x1 - x0) * x, y0 + (y1 - y0) * y, z0 + (z1 - z0) * z)  # noqa: E731
    with ng.TaskManager():
        return MakeStructured3DMesh(hexes=True, nx=nx, ny=ny, nz=nz, mapping=mp)


def _Hz():
    return ng.CoefficientFunction((0, 0, H0))


def test_image_string_and_group():
    assert core.parse_image_string("+x-z") == [(0, 1), (2, -1)]
    assert core.parse_image_string("-x+y-z") == [(0, -1), (1, 1), (2, -1)]
    # group of 2 planes (+x,-z): {x}:+1, {z}:-1, {x,z}:-1
    g = core.image_group([(0, 1), (2, -1)])
    assert ((0,), 1) in g and ((2,), -1) in g and ((0, 2), -1) in g and len(g) == 3
    with pytest.raises(ValueError):
        core.parse_image_string("xx")          # malformed
    with pytest.raises(ValueError):
        core.parse_image_string("+x+x")        # repeated axis


@pytest.fixture(scope="module")
def _full_linear():
    with ng.TaskManager():
        r = vim.hdiv_demag_solve(_box(-A, A, -A, A, -A, A, 4, 4, 4), mu_r=MU_R, H_ext=_Hz(),
                                 scalable=False)
    return r["M_avg"][2], r["demag"]


@pytest.mark.parametrize("name,boxargs,image,nx,ny,nz", [
    ("half",    (-A, A, -A, A, 0, A), "-z",     4, 4, 2),
    ("quarter", (0, A, -A, A, 0, A),  "+x-z",   2, 4, 2),
    ("octant",  (0, A, 0, A, 0, A),   "+x+y-z", 2, 2, 2),
])
def test_ima_reduced_matches_full_linear(_full_linear, name, boxargs, image, nx, ny, nz):
    """A reduced (1/2,1/4,1/8) cube + IMA reproduces the full cube M_avg + demag at high mu_r."""
    Mf, Df = _full_linear
    with ng.TaskManager():
        r = vim.hdiv_demag_solve(_box(*boxargs, nx, ny, nz), mu_r=MU_R, H_ext=_Hz(),
                                 image=image, scalable=False)
    relM = abs(r["M_avg"][2] - Mf) / abs(Mf)
    relD = abs(r["demag"] - Df) / abs(Df)
    assert relM < 3e-3, f"{name} IMA M_avg {r['M_avg'][2]:.2f} vs full {Mf:.2f} (rel {relM:.2e})"
    assert relD < 3e-3, f"{name} IMA demag {r['demag']:.5f} vs full {Df:.5f} (rel {relD:.2e})"


def test_ima_quarter_matches_full_nonlinear():
    """NONLINEAR quarter + IMA matches the full nonlinear cube (the rho/div-M path under saturation)."""
    Hs = [0, 50, 100, 300, 1000, 3000, 1e4, 3e4, 1e5, 1e6]
    Bs = [0, 0.6, 1.0, 1.5, 1.75, 1.9, 2.0, 2.05, 2.1, 2.1 + MU0 * (1e6 - 1e5)]
    BH = [[h, b] for h, b in zip(Hs, Bs)]
    Hdrive = ng.CoefficientFunction((0, 0, 5000.0))
    with ng.TaskManager():
        rf = vim.hdiv_demag_solve(_box(-A, A, -A, A, -A, A, 3, 3, 3), bh_table=BH, H_ext=Hdrive,
                                  scalable=False)
        rq = vim.hdiv_demag_solve(_box(0, A, -A, A, 0, A, 2, 3, 2), bh_table=BH, H_ext=Hdrive,
                                  image="+x-z", scalable=False)
    rel = abs(rq["M_avg"][2] - rf["M_avg"][2]) / abs(rf["M_avg"][2])
    assert rq["nonlinear"] and rf["nonlinear"]
    assert rel < 2e-2, f"nonlinear quarter+IMA M_avg {rq['M_avg'][2]:.1f} vs full {rf['M_avg'][2]:.1f} (rel {rel:.2e})"


def test_ima_scalable_true_conflict():
    """image=... with an EXPLICIT scalable=True is fail-loud (IMA is dense-only for now)."""
    with ng.TaskManager():
        with pytest.raises(NotImplementedError):
            vim.hdiv_demag_solve(_box(0, A, -A, A, 0, A, 2, 4, 2), mu_r=MU_R, H_ext=_Hz(),
                                 image="+x-z", scalable=True)
