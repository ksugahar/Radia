"""Gate-2 hex/wedge: the FEEC HDiv-VIM demag solve on HEX and WEDGE meshes, via the POLYTOPE charge Gram
(build_demag's analytic Gram generalized from tet/triangle to any flat-faced convex cell + quad face: cell
Newtonian potential = divergence-theorem sum over the convex-hull faces of the exact Wilton triangle
potential; a quad face = two flat triangles; centroid-fan / Dunavant outer quadrature).  Both the dense
Python path (analytic_charge_gram) AND the scalable C++ charge-Gram H-matrix (the polytope triangle-soup
_ChargeGramHMatrix ctor) handle hex/wedge -- scalable=None (default) now uses the C++ H-matrix for ANY
element type.

Locks:
  (1) hex- and wedge-meshed cubes have demag_z ~ 1/3 (the physics gate -- isotropic uniform-M body);
  (2) a hex cube's M_avg matches a tet mesh of the same cube (element-type-independent physics);
  (3) hdiv_demag_solve(scalable=True) on a hex mesh now WORKS (the C++ polytope Gram) and matches the
      dense polytope Gram (scalable=False) -- the Item-1 scalable hex/wedge increment;
  (4) rad.Solve(demag_backend='hdiv') dispatches the HEX solve == direct hdiv_demag_solve to machine
      precision and writes per-element M back (ObjM) onto the Radia hex elements;
  (5) soft_iron_from_mesh on a WEDGE mesh raises a clean NotImplementedError (the container write-back
      has no ObjWedge import yet) -- while direct hdiv_demag_solve(wedge) DOES solve it (locked by (1)).
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

MU0 = 4.0e-7 * math.pi
H0 = 1000.0
L = 0.02
MU_R = 100.0


def _cube(hexes, n):
    """Structured [0,L]^3 cube of `hexes` (else prisms/wedges), n^3 cells."""
    mp = lambda x, y, z: (L * x, L * y, L * z)  # noqa: E731
    with ng.TaskManager():
        return MakeStructured3DMesh(hexes=hexes, prism=(not hexes), nx=n, ny=n, nz=n, mapping=mp)


def _Hz():
    return ng.CoefficientFunction((0, 0, H0))


def test_hex_cube_demag_third():
    """A hex-meshed cube: the polytope-Gram demag factor (Rayleigh quotient) is ~1/3 (isotropic)."""
    with ng.TaskManager():
        res = vim.hdiv_demag_solve(_cube(True, 4), mu_r=MU_R, H_ext=_Hz())
    assert 0.30 < res["demag"] < 0.36, f"hex cube demag {res['demag']:.4f} not ~1/3"


def test_wedge_cube_demag_third():
    """A wedge(prism)-meshed cube through the SAME polytope Gram: demag ~1/3 (mixed tri+quad faces)."""
    with ng.TaskManager():
        res = vim.hdiv_demag_solve(_cube(False, 4), mu_r=MU_R, H_ext=_Hz())
    assert 0.30 < res["demag"] < 0.36, f"wedge cube demag {res['demag']:.4f} not ~1/3"


def test_hex_matches_tet_cube():
    """The solved M_avg is element-type independent: a hex cube agrees with a tet mesh of the same cube
    (the physics is the geometry's, not the mesh's)."""
    from netgen.occ import Box, OCCGeometry
    with ng.TaskManager():
        r_hex = vim.hdiv_demag_solve(_cube(True, 4), mu_r=MU_R, H_ext=_Hz())
        tetm = ng.Mesh(OCCGeometry(Box((0, 0, 0), (L, L, L))).GenerateMesh(maxh=L / 4))
        r_tet = vim.hdiv_demag_solve(tetm, mu_r=MU_R, H_ext=_Hz())
    rel = abs(r_hex["M_avg"][2] - r_tet["M_avg"][2]) / abs(r_tet["M_avg"][2])
    assert rel < 2e-2, f"hex M_avg {r_hex['M_avg'][2]:.1f} vs tet {r_tet['M_avg'][2]:.1f} (rel {rel:.2e})"


def test_hex_nonlinear_matches_tet_cube():
    """NONLINEAR (bh_table) soft iron on a hex cube: the polytope Gram feeds the same element-agnostic
    Newton as tet, so a hex cube's saturated M_avg matches a tet mesh of the same cube."""
    from netgen.occ import Box, OCCGeometry
    Hs = [0, 50, 100, 300, 1000, 3000, 1e4, 3e4, 1e5, 1e6]
    Bs = [0, 0.6, 1.0, 1.5, 1.75, 1.9, 2.0, 2.05, 2.1, 2.1 + MU0 * (1e6 - 1e5)]
    BH = [[h, b] for h, b in zip(Hs, Bs)]
    Hdrive = ng.CoefficientFunction((0, 0, 5000.0))     # strong enough to push into the knee
    with ng.TaskManager():
        r_hex = vim.hdiv_demag_solve(_cube(True, 3), bh_table=BH, H_ext=Hdrive)
        tetm = ng.Mesh(OCCGeometry(Box((0, 0, 0), (L, L, L))).GenerateMesh(maxh=L / 3))
        r_tet = vim.hdiv_demag_solve(tetm, bh_table=BH, H_ext=Hdrive)
    assert r_hex["nonlinear"] and r_tet["nonlinear"]
    rel = abs(r_hex["M_avg"][2] - r_tet["M_avg"][2]) / abs(r_tet["M_avg"][2])
    assert rel < 3e-2, f"nonlinear hex M_avg {r_hex['M_avg'][2]:.1f} vs tet {r_tet['M_avg'][2]:.1f} (rel {rel:.2e})"


def test_scalable_hex_matches_dense():
    """Item-1 scalable hex/wedge Gram: scalable=True (C++ polytope charge-Gram H-matrix) reproduces
    scalable=False (the dense Python polytope Gram) on a hex cube to a tight tolerance -- the C++ entry ==
    the dense analytic_charge_gram polytope entry (same hull triangulation, same built-in quad rules)."""
    with ng.TaskManager():
        r_scal = vim.hdiv_demag_solve(_cube(True, 3), mu_r=MU_R, H_ext=_Hz(), scalable=True)
        r_dense = vim.hdiv_demag_solve(_cube(True, 3), mu_r=MU_R, H_ext=_Hz(), scalable=False)
    assert abs(r_scal["demag"] - r_dense["demag"]) < 1e-3 * abs(r_dense["demag"]) + 1e-4, \
        f"hex scalable demag {r_scal['demag']:.5f} vs dense {r_dense['demag']:.5f}"
    rel = abs(r_scal["M_avg"][2] - r_dense["M_avg"][2]) / abs(r_dense["M_avg"][2])
    assert rel < 1e-3, f"hex scalable M_avg {r_scal['M_avg'][2]:.2f} vs dense {r_dense['M_avg'][2]:.2f} (rel {rel:.2e})"


def test_scalable_wedge_matches_dense():
    """Same Item-1 check on a WEDGE cube (mixed tri+quad faces): the C++ polytope Gram == the dense one."""
    with ng.TaskManager():
        r_scal = vim.hdiv_demag_solve(_cube(False, 3), mu_r=MU_R, H_ext=_Hz(), scalable=True)
        r_dense = vim.hdiv_demag_solve(_cube(False, 3), mu_r=MU_R, H_ext=_Hz(), scalable=False)
    rel = abs(r_scal["M_avg"][2] - r_dense["M_avg"][2]) / abs(r_dense["M_avg"][2])
    assert rel < 1e-3, f"wedge scalable M_avg {r_scal['M_avg'][2]:.2f} vs dense {r_dense['M_avg'][2]:.2f} (rel {rel:.2e})"


def test_radsolve_hdiv_hex_dispatch():
    """rad.Solve(demag_backend='hdiv') on a HEX soft_iron_from_mesh container == direct hdiv_demag_solve
    to machine precision, and writes the per-element M back onto the Radia hex elements (ObjM)."""
    rad.UtiDelAll()
    mesh = _cube(True, 4)
    with ng.TaskManager():
        direct = vim.hdiv_demag_solve(mesh, mu_r=MU_R, H_ext=_Hz())
    iron = vim.soft_iron_from_mesh(mesh, mu_r=MU_R)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])      # free-space B whose H is H0
    cont = rad.ObjCnt([iron, bkg])
    prev = rad.get_demag_backend()        # set_demag_backend returns the NEW value, not prev
    rad.set_demag_backend("hdiv")
    try:
        with ng.TaskManager():
            res = rad.Solve(cont, 1e-6, 1000, 0)
    finally:
        rad.set_demag_backend(prev)
    rel = abs(res["M_avg"][2] - direct["M_avg"][2]) / abs(direct["M_avg"][2])
    assert rel < 1e-9, f"rad.Solve(hdiv) hex M_avg {res['M_avg'][2]:.1f} != direct {direct['M_avg'][2]:.1f} (rel {rel:.2e})"
    objm = rad.ObjM(iron)
    assert len(objm) == direct["n_el"], f"{len(objm)} Radia hex handles vs {direct['n_el']} HDiv elements"
    mz = float(np.mean([m[2] for (_c, m) in objm]))
    assert abs(mz - res["M_avg"][2]) / abs(res["M_avg"][2]) < 0.05, \
        f"ObjM write-back mean Mz {mz:.1f} not consistent with solved M_avg {res['M_avg'][2]:.1f}"
    rad.UtiDelAll()


def test_radsolve_hdiv_wedge_dispatch():
    """rad.Solve(demag_backend='hdiv') on a WEDGE soft_iron_from_mesh container (ObjWedge per element)
    == direct hdiv_demag_solve to machine precision, with per-element M written back (ObjM)."""
    rad.UtiDelAll()
    mesh = _cube(False, 3)
    with ng.TaskManager():
        direct = vim.hdiv_demag_solve(mesh, mu_r=MU_R, H_ext=_Hz())
    iron = vim.soft_iron_from_mesh(mesh, mu_r=MU_R)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    cont = rad.ObjCnt([iron, bkg])
    prev = rad.get_demag_backend()        # set_demag_backend returns the NEW value, not prev
    rad.set_demag_backend("hdiv")
    try:
        with ng.TaskManager():
            res = rad.Solve(cont, 1e-6, 1000, 0)
    finally:
        rad.set_demag_backend(prev)
    rel = abs(res["M_avg"][2] - direct["M_avg"][2]) / abs(direct["M_avg"][2])
    assert rel < 1e-9, f"rad.Solve(hdiv) wedge M_avg {res['M_avg'][2]:.1f} != direct {direct['M_avg'][2]:.1f} (rel {rel:.2e})"
    objm = rad.ObjM(iron)
    assert len(objm) == direct["n_el"], f"{len(objm)} Radia wedge handles vs {direct['n_el']} HDiv elements"
    mz = float(np.mean([m[2] for (_c, m) in objm]))
    assert abs(mz - res["M_avg"][2]) / abs(res["M_avg"][2]) < 0.05, \
        f"ObjM write-back mean Mz {mz:.1f} not consistent with solved M_avg {res['M_avg'][2]:.1f}"
    rad.UtiDelAll()
