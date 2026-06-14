"""M0 nonlinear parity gate (productionization, docs/hdiv_vim/PRODUCTIONIZATION.md): the HDiv-VIM
production entry `radia.vim.hdiv_demag_solve` (NONLINEAR, COIL source) vs the shipped yano-type MSC/MMM
(`rad.Solve`) on a coil-driven C-type electromagnet, evaluated by the ENGINEERING quantity -- the
GAP-CENTRE flux density B.

C-yoke iron (thin in z, with a gap/bore) driven by a rectangular current loop ENCIRCLING the back leg.
Both solvers use the SAME iron mesh + SAME coil + SAME MatSatIsoTab BH table, and the gap B is then
evaluated with the SAME field kernel (radia.Fld) -- only the magnetization M differs:
  yano:  rad.Solve(iron tets [MMM] + coil) -> solved M ; B_gap = rad.Fld(iron+coil,'b',gap)
  HDiv:  hdiv_demag_solve -> per-element M ; Radia iron rebuilt from that M (no solve) ;
         B_gap = rad.Fld(iron_hdivM + coil,'b',gap)
HDiv applied field = the coil's Biot-Savart H (rad.RadiaField(coil,'h')); KELVIN-LESS (iron-only mesh).

Measured (2026-06-15, NI=8000 A, 2715 tets, Msat~1.27e6 A/m): gap-centre |B| HDiv 0.0982 T vs yano
0.0979 T agree to 0.32% (bore centre), 0.18% at the gap opening.  The HDiv Anderson-Hantila fixed
point converges in ~120 outer iters.

The gate locks: HDiv and yano agree within 1.5% on the gap-centre flux density (the engineering
quantity), with a bounded outer-iteration count.  Transitional gate -- retires when yano-type is
sealed (M5).
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
GAP_PTS = [[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]]   # bore centre + gap-opening point


def _cyoke():
    return (Box(OPnt(-0.06, -0.06, -0.02), OPnt(0.06, 0.06, 0.02))
            - Box(OPnt(-0.035, -0.035, -0.03), OPnt(0.035, 0.035, 0.03))
            - Box(OPnt(0.018, -0.07, -0.03), OPnt(0.07, 0.07, 0.03)))


def _coil():
    """Rectangular loop in the x-z plane at y=0 encircling the back leg (threads the bore)."""
    pts = [[-0.075, 0.0, -0.035], [-0.025, 0.0, -0.035], [-0.025, 0.0, 0.035],
           [-0.075, 0.0, 0.035], [-0.075, 0.0, -0.035]]
    return rad.ObjFlmCur(pts, NI)


def _yano_gapB(mesh):
    """Gap-centre B (Tesla) via yano-type MMM rad.Solve (iron tets + coil), same field kernel."""
    rad.UtiDelAll()
    coil = _coil()
    cont = nmi.netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m', verbose=False)
    rad.MatApl(cont, rad.MatSatIsoTab(BH))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)    # yano backend deprecated; expected here
        rad.Solve(rad.ObjCnt([cont, coil]), 1e-6, 3000, 0)     # LU
    B = np.array(rad.Fld(rad.ObjCnt([cont, coil]), 'b', GAP_PTS)).reshape(-1, 3)
    rad.UtiDelAll()
    return B


def test_ctype_coil_nonlinear_gap_field():
    """HDiv-VIM and yano-type agree within 1.5% on the GAP-CENTRE flux density B (the engineering
    quantity) of the coil-driven nonlinear C-type; HDiv converges in a bounded outer-iter count."""
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(_cyoke()).GenerateMesh(maxh=0.008))

    # HDiv-VIM solve -> per-element M -> rebuild Radia iron from that M -> gap B (same field kernel)
    rad.UtiDelAll()
    coil = _coil()
    with ng.TaskManager():
        res = hdiv_demag_solve(mesh, bh_table=BH, H_ext=rad.RadiaField(coil, 'h'))
    M_el = res["M"]
    coil2 = _coil()
    iron_hdiv = nmi.netgen_mesh_to_radia(
        mesh, material=lambda i: {'magnetization': M_el[i].tolist()}, units='m', verbose=False)
    B_hdiv = np.array(rad.Fld(rad.ObjCnt([iron_hdiv, coil2]), 'b', GAP_PTS)).reshape(-1, 3)
    rad.UtiDelAll()

    B_yano = _yano_gapB(mesh)

    assert res["nonlinear"] is True
    for k, p in enumerate(GAP_PTS):
        rel = float(np.linalg.norm(B_hdiv[k] - B_yano[k]) / (np.linalg.norm(B_yano[k]) + 1e-30))
        assert rel < 0.015, f"gap B at {p}: HDiv {B_hdiv[k]} vs yano {B_yano[k]} rel {rel:.2e}"
    # bore field is a meaningful fraction of a Tesla (genuinely driven, into the BH knee)
    assert np.linalg.norm(B_yano[0]) > 0.05
    assert res["iters"] < 200, f"HDiv nonlinear Anderson-Hantila outer iters not bounded: {res['iters']}"
