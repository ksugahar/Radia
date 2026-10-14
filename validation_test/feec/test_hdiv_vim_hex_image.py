"""HDiv-VIM hex + wedge IMA (image method): a reduced (1/2, 1/4) symmetric model + image= reproduces the
FULL model.

The FEEC HDiv-VIM folds the mirror-image charge interactions into the charge Gram so a reduced-symmetry
model (cut along a symmetry plane) gives the same magnetization as the full model.  For TET this is the C++
QuadDotRefl->PhiInner fold (2026-07-04); for HEX / WEDGE it is the reflected-block QuadBlockHex/Wedge(mask)
(2026-07-04) -- the source host is reflected on the 3-bit axis mask (reflected sub-centroid drives the
near/far grading + Duffy corner; the outer eval point is reflected; the pair is never treated as SELF).

Geometry: a cube iron [-a,a]^3 in a uniform applied field Hz.  Hz is PARALLEL to the x=0 and y=0 mirror
planes, so a half model (x in [0,a]) uses image='+x' and a quarter model (x,y in [0,a]) image='+x+y'
(symmetric '+').  Metric: the volume-averaged magnetization M_avg (physical, method-independent) and the
Rayleigh demag factor -- both must match the full model.

VALIDATED 2026-07-04 (mask=0 stays byte-identical to the direct block; the image fold reproduces the full
model): hex 1/2 = 6.0e-5, hex 1/4 = 1.4e-4, wedge 1/2 = 7.4e-4 (all << the 5e-3 gate).
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

pytestmark = pytest.mark.compute_host

A = 0.01
MU_R = 1000.0
H0 = 1.0e4
HCF = ng.CoefficientFunction((0.0, 0.0, H0))
TOL = 5e-3          # M_avg rel gate (validated values are 1-2 orders below)


def _box(x0, x1, y0, y1, z0, z1, nx, ny, nz, hexes=True, prism=False):
    return MakeStructured3DMesh(
        hexes=hexes, prism=prism, nx=nx, ny=ny, nz=nz,
        mapping=lambda X, Y, Z: (x0 + (x1 - x0) * X, y0 + (y1 - y0) * Y, z0 + (z1 - z0) * Z))


def _rel(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / np.linalg.norm(np.asarray(b)))


def test_hex_image_half_reproduces_full():
    """hex half (x in [0,A]) + image='+x' == full cube."""
    rad.UtiDelAll()
    with ng.TaskManager():
        full = vim.Solve(_box(-A, A, -A, A, -A, A, 4, 4, 4), mu_r=MU_R, H_ext=HCF)
        half = vim.Solve(_box(0.0, A, -A, A, -A, A, 2, 4, 4), mu_r=MU_R, H_ext=HCF, image='+x')
    rel = _rel(half["M_avg"], full["M_avg"])
    assert rel < TOL, f"hex 1/2 image='+x' M_avg rel {rel:.2e} vs full (Mz {half['M_avg'][2]:.1f} != {full['M_avg'][2]:.1f})"
    assert abs(half["demag"] - full["demag"]) < 1e-3, f"hex 1/2 demag {half['demag']:.5f} != full {full['demag']:.5f}"
    rad.UtiDelAll()


def test_hex_image_quarter_reproduces_full():
    """hex quarter (x,y in [0,A]) + image='+x+y' == full cube (multi-image fold: 3 masks)."""
    rad.UtiDelAll()
    with ng.TaskManager():
        full = vim.Solve(_box(-A, A, -A, A, -A, A, 4, 4, 4), mu_r=MU_R, H_ext=HCF)
        quarter = vim.Solve(_box(0.0, A, 0.0, A, -A, A, 2, 2, 4), mu_r=MU_R, H_ext=HCF, image='+x+y')
    rel = _rel(quarter["M_avg"], full["M_avg"])
    assert rel < TOL, f"hex 1/4 image='+x+y' M_avg rel {rel:.2e} vs full"
    assert abs(quarter["demag"] - full["demag"]) < 1e-3, f"hex 1/4 demag {quarter['demag']:.5f} != full {full['demag']:.5f}"
    rad.UtiDelAll()


def test_wedge_image_half_reproduces_full():
    """wedge (prism) half (x in [0,A]) + image='+x' == full cube."""
    rad.UtiDelAll()
    with ng.TaskManager():
        full = vim.Solve(_box(-A, A, -A, A, -A, A, 3, 3, 3, hexes=False, prism=True),
                                    mu_r=MU_R, H_ext=HCF)
        half = vim.Solve(_box(0.0, A, -A, A, -A, A, 2, 3, 3, hexes=False, prism=True),
                                    mu_r=MU_R, H_ext=HCF, image='+x')
    rel = _rel(half["M_avg"], full["M_avg"])
    assert rel < TOL, f"wedge 1/2 image='+x' M_avg rel {rel:.2e} vs full (Mz {half['M_avg'][2]:.1f} != {full['M_avg'][2]:.1f})"
    assert abs(half["demag"] - full["demag"]) < 1e-3, f"wedge 1/2 demag {half['demag']:.5f} != full {full['demag']:.5f}"
    rad.UtiDelAll()


def test_radsolve_auto_hex_uses_hdiv_with_image():
    """rad.Solve auto-routes a mesh-backed HEX soft iron to the HDiv-VIM (CLAUDE.md DIRECTION 2026-07-04),
    including with image= (the reflected-block IMA)."""
    from radia.vim import _radsolve
    MU0 = 4.0e-7 * math.pi
    rad.set_demag_backend("auto")
    rad.UtiDelAll()
    with ng.TaskManager():
        iron = vim.MeshSoftIron(_box(0.0, A, -A, A, -A, A, 2, 4, 4), mu_r=MU_R)
        assert _radsolve.is_hdiv_eligible(iron)          # mesh-backed hex -> HDiv-eligible (the flip)
        bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
        res = rad.Solve(rad.ObjCnt([iron, bkg]), 1e-6, 2000, 0, image='+x')
    assert isinstance(res, dict) and "demag" in res, "auto hex+image did not route to the HDiv-VIM dispatch"
    assert 0.2 < res["demag"] < 0.45, f"hex demag {res['demag']:.4f} out of physical band"
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
