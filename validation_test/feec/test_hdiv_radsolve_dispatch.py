"""Gate-2: rad.Solve dispatches the FEEC HDiv-VIM backend (demag_backend='hdiv').

Previously rad.Solve only ran the legacy six-face surface-charge MSC; selecting 'hdiv' raised
NotImplementedError.  This locks the wiring (radia.vim.soft_iron_from_mesh +
radia.vim._radsolve.dispatch): a soft-iron container built from an NGSolve mesh, solved
via rad.Solve(demag_backend='hdiv'), reproduces a direct radia.vim.hdiv_demag_solve to
machine precision and writes the per-element M back so rad.ObjM/rad.Fld reflect it.

With hdiv wired in, rad.Solve keeps a working soft-iron demag path: a TET iron solves via the HDiv-VIM
(RT1), and a HEX / WEDGE iron is routed by the 'auto' split to the collocation MMMM backend (HDiv-VIM is
tet/RT1-only, 2026-06-29 -- an explicit demag_backend='hdiv' on a non-tet iron fails loud).
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry  # noqa: E402

MU0 = 4.0e-7 * math.pi
H0 = 1000.0
L = 0.02
MU_R = 100.0


def _tet_cube_mesh(maxh=L / 6):
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Box((-L / 2, -L / 2, -L / 2), (L / 2, L / 2, L / 2))).GenerateMesh(maxh=maxh))


def test_radsolve_hdiv_equals_direct():
    """rad.Solve(demag_backend='hdiv') on a soft_iron_from_mesh container == direct hdiv_demag_solve
    (same mesh, same applied field) to machine precision, and writes M back onto the Radia elements."""
    rad.UtiDelAll()
    mesh = _tet_cube_mesh()
    with ng.TaskManager():
        direct = vim.hdiv_demag_solve(mesh, MU_R, ng.CoefficientFunction((0, 0, H0)))

    iron = vim.soft_iron_from_mesh(mesh, mu_r=MU_R)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])      # free-space B whose H is H0
    cont = rad.ObjCnt([iron, bkg])
    # No set_demag_backend: a mesh-backed TET soft iron (soft_iron_from_mesh) auto-routes to HDiv-VIM RT1.
    with ng.TaskManager():
        res = rad.Solve(cont, 1e-6, 1000, 0)

    rel = abs(res["M_avg"][2] - direct["M_avg"][2]) / abs(direct["M_avg"][2])
    assert rel < 1e-9, f"rad.Solve(hdiv) M_avg {res['M_avg'][2]:.3f} != direct {direct['M_avg'][2]:.3f} (rel {rel:.2e})"

    # write-back: rad.ObjM reflects the HDiv per-element M (element-mean ~ volume-avg within tet-size spread)
    objm = rad.ObjM(iron)
    mz_objm = float(np.mean([m[2] for (_c, m) in objm]))
    assert abs(mz_objm - res["M_avg"][2]) / abs(res["M_avg"][2]) < 0.05, \
        f"ObjM write-back mean Mz {mz_objm:.1f} not consistent with solved M_avg {res['M_avg'][2]:.1f}"
    rad.UtiDelAll()


# (test_radsolve_hdiv_image_passed_through removed 2026-06-29: IMA image symmetry is retired from HDiv-VIM
#  (and the mesh was hex) -- reduced symmetric models use collocation MMMM.  Retirement locked by
#  test_hdiv_vim_rt1_contract.py; the hex->collocation MMMM auto-route by its hex test.)


def test_meshless_hex_soft_iron_not_registered_uses_collocation_mmmm():
    """A hex soft iron built the mesh-less way (ObjHexahedron + MatLin, NOT via soft_iron_from_mesh) has
    no mesh association, so rad.Solve does NOT route it to the HDiv-VIM -- it falls through to the C++
    solver, which now solves it with the six-face surface-charge MSC (the Error203 guard was removed 2026-06-19).
    The full surface-charge MSC physics (cube demag ~1/3) is locked by tests/test_demag_backend.py."""
    from radia.vim import _radsolve
    rad.set_demag_backend("auto")
    rad.UtiDelAll()
    hexv = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
    obj = rad.ObjHexahedron([[c * 0.005 for c in v] for v in hexv], [0, 0, 0])
    rad.MatApl(obj, rad.MatLin(1000.0))
    cont = rad.ObjCnt([obj])
    assert not _radsolve.is_registered(cont)        # mesh-less -> not HDiv-registered
    rad.Solve(cont, 1e-6, 100, 0)                   # no Error203; surface-charge MSC solves (M stays 0, no source)
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
