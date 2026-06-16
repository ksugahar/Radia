"""Gate-2: rad.Solve dispatches the FEEC HDiv-VIM backend (demag_backend='hdiv').

Previously rad.Solve only ran the legacy yano-type MSC; selecting 'hdiv' raised
NotImplementedError.  This locks the wiring (radia.vim.soft_iron_from_mesh +
radia.vim._radsolve.dispatch): a soft-iron container built from an NGSolve mesh, solved
via rad.Solve(demag_backend='hdiv'), reproduces a direct radia.vim.hdiv_demag_solve to
machine precision and writes the per-element M back so rad.ObjM/rad.Fld reflect it.

This is the prerequisite for sealing (deleting) the yano-type C++: with hdiv wired in,
rad.Solve keeps a working hex/wedge soft-iron demag path through the HDiv-VIM.
(tet-first: build_demag's volume self-energy is tet-only, so the dispatch requires a tet mesh.)
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
    prev = rad.set_demag_backend("hdiv")
    try:
        with ng.TaskManager():
            res = rad.Solve(cont, 1e-6, 1000, 0)
    finally:
        rad.set_demag_backend(prev)

    rel = abs(res["M_avg"][2] - direct["M_avg"][2]) / abs(direct["M_avg"][2])
    assert rel < 1e-9, f"rad.Solve(hdiv) M_avg {res['M_avg'][2]:.3f} != direct {direct['M_avg'][2]:.3f} (rel {rel:.2e})"

    # write-back: rad.ObjM reflects the HDiv per-element M (element-mean ~ volume-avg within tet-size spread)
    objm = rad.ObjM(iron)
    mz_objm = float(np.mean([m[2] for (_c, m) in objm]))
    assert abs(mz_objm - res["M_avg"][2]) / abs(res["M_avg"][2]) < 0.05, \
        f"ObjM write-back mean Mz {mz_objm:.1f} not consistent with solved M_avg {res['M_avg'][2]:.1f}"
    rad.UtiDelAll()


def test_radsolve_hdiv_unregistered_raises():
    """rad.Solve(demag_backend='hdiv') on a container with no HDiv-registered iron raises
    NotImplementedError redirecting to soft_iron_from_mesh (fail-loud, No-Fallback)."""
    rad.UtiDelAll()
    # a plain hex magnet container -- NOT built via soft_iron_from_mesh, so no mesh association
    hexv = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
    obj = rad.ObjHexahedron([[c * 0.005 for c in v] for v in hexv], [0, 0, 0])
    cont = rad.ObjCnt([obj])
    prev = rad.set_demag_backend("hdiv")
    try:
        with pytest.raises(NotImplementedError):
            rad.Solve(cont, 1e-6, 100, 0)
    finally:
        rad.set_demag_backend(prev)
    rad.UtiDelAll()
