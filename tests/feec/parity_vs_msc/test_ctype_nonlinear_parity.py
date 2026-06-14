"""M0 nonlinear parity gate (productionization, docs/hdiv_vim/PRODUCTIONIZATION.md): the HDiv-VIM
production entry `radia.vim.hdiv_demag_solve` (NONLINEAR, COIL source) vs the shipped yano-type MSC
(`rad.Solve`) on a coil-driven C-type electromagnet.

C-yoke iron (thin in z, with a gap) driven by a rectangular current loop ENCIRCLING the back leg (loop
in the x-z plane at y=0, threading the window) -> circulating flux around the yoke, across the gap.
Both solvers use the SAME iron mesh + SAME coil + SAME MatSatIsoTab BH table.  HDiv's applied field is
the coil's Biot-Savart H (rad.RadiaField(coil,'h')); yano sees the same coil object in its container.
KELVIN-LESS (iron-only mesh).

Measured (2026-06-15, NI=8000 A, 270 tets, Msat~1.27e6 A/m): HDiv (Anderson-Hantila) and yano agree to
~0.6% on volume-average M (|M| ~ 507 kA/m, ~40% of saturation -- genuinely into the BH knee).  The
Anderson-Hantila fixed point converges in ~120 outer iters (more than Newton's ~20, but each step is a
single cheap M_mass^-1-GMRES with no tangent assembly / line search -> within ~2x of yano wall-clock).

The gate locks: HDiv solves the coil-driven NONLINEAR C-type and agrees with yano within 2% (vector),
with a bounded outer-iteration count.  Transitional gate -- retires when yano-type is sealed (M5).
"""
import math
import warnings

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt as OPnt  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402
from radia import netgen_mesh_import as nmi  # noqa: E402

MU0 = 4e-7 * math.pi
CHI0 = 2000.0
MSAT = 1.6 / MU0                        # M saturation ~ 1.6 T / mu0 ~ 1.27e6 A/m
_Hs = np.concatenate([[0.0], np.logspace(-1, 7, 80)])
_Ms = CHI0 * _Hs / (1.0 + CHI0 * _Hs / MSAT)
_Bs = MU0 * (_Hs + _Ms)
BH = [[float(h), float(b)] for h, b in zip(_Hs, _Bs)]   # MatSatIsoTab data both solvers consume
NI = 8000.0                              # coil ampere-turns


def _cyoke():
    return (Box(OPnt(-0.06, -0.06, -0.02), OPnt(0.06, 0.06, 0.02))
            - Box(OPnt(-0.035, -0.035, -0.03), OPnt(0.035, 0.035, 0.03))
            - Box(OPnt(0.018, -0.07, -0.03), OPnt(0.07, 0.07, 0.03)))


def _coil():
    """Rectangular loop in the x-z plane at y=0 encircling the back leg (threads the window)."""
    pts = [[-0.075, 0.0, -0.035], [-0.025, 0.0, -0.035], [-0.025, 0.0, 0.035],
           [-0.075, 0.0, 0.035], [-0.075, 0.0, -0.035]]
    return rad.ObjFlmCur(pts, NI)


def _yano_Mavg(mesh):
    rad.UtiDelAll()
    coil = _coil()
    cont = nmi.netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m', verbose=False)
    rad.MatApl(cont, rad.MatSatIsoTab(BH))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)    # yano backend deprecated; expected here
        rad.Solve(rad.ObjCnt([cont, coil]), 1e-6, 3000, 0)
    cents, vols = [], []
    for el in mesh.Elements(ng.VOL):
        vs = np.array([mesh[v].point for v in el.vertices])
        cents.append(vs.mean(0)); vols.append(abs(np.linalg.det(vs[1:] - vs[0])) / 6.0)
    cents, vols = np.array(cents), np.array(vols)
    M = np.array(rad.Fld(cont, 'm', cents.tolist())).reshape(-1, 3)
    Mavg = (M * vols[:, None]).sum(0) / vols.sum()
    rad.UtiDelAll()
    return Mavg


def test_ctype_coil_nonlinear_parity():
    """HDiv-VIM and yano-type MSC agree within 2% (vector volume-average M) on the coil-driven
    nonlinear C-type electromagnet; HDiv converges in bounded damped-Newton iters."""
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(_cyoke()).GenerateMesh(maxh=0.016))

    rad.UtiDelAll()
    coil = _coil()
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, bh_table=BH, H_ext=rad.RadiaField(coil, 'h'))
    M_hdiv = res["M_avg"]
    rad.UtiDelAll()

    M_yano = _yano_Mavg(mesh)
    rel = float(np.linalg.norm(M_hdiv - M_yano) / (np.linalg.norm(M_yano) + 1e-30))
    assert res["nonlinear"] is True
    assert rel < 0.02, f"C-type HDiv {M_hdiv} vs yano {M_yano} parity rel {rel:.2e}"
    assert res["iters"] < 200, f"HDiv nonlinear Anderson-Hantila outer iters not bounded: {res['iters']}"
    # the drive is genuinely into the nonlinear knee (|M| a large fraction of saturation)
    assert np.linalg.norm(M_yano) > 0.3 * MSAT
