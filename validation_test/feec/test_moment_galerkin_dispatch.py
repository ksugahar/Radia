"""Gate: rad.Solve dispatches the SYMMETRIC moment-Galerkin MMMM backend (demag_backend='moment_galerkin').

A soft-iron container built from an NGSolve mesh (radia.vim.soft_iron_from_mesh) now supports THREE
co-valid demag backends on the SAME container: 'yano' (multipole-moment MMM surface-charge MSC),
'hdiv' (FEEC HDiv-VIM), and 'moment_galerkin' (the symmetric N = B^T G B moment-Galerkin MMMM,
radia.moment_galerkin).  soft_iron_from_mesh stores the per-element vertex lists (same mesh.Elements(VOL)
order as the Radia handles) so the dispatch (radia.vim._radsolve.dispatch_moment_galerkin) can solve the
container and write the per-element M back via ObjSetM.

Locks: the dispatch reproduces a DIRECT moment_galerkin_demag_solve on the same elements + applied field
to machine precision (the wiring -- element order, chi, H_ext, write-back), the "mg"/"galerkin" aliases
resolve to the same result, the mesh-less error path raises, and the stray field is the same physical
ballpark as the 'yano' backend (the two differ only by the dipole-vs-6DOF discretisation).

moment-Galerkin is LINEAR-only; nonlinear (bh_table) and IMA image symmetry route through 'hdiv'/'yano'.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.moment_galerkin as mg  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
from radia.netgen_mesh_import import extract_elements  # noqa: E402
from radia.vim import _radsolve  # noqa: E402

MU0 = 4.0e-7 * math.pi
H0 = 1000.0
L = 0.02
MU_R = 1000.0


def _hex_cube_mesh(n=2):
    with ng.TaskManager():
        return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n,
                                    mapping=lambda x, y, z: (L * x, L * y, L * z))


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def test_dispatch_equals_direct_moment_galerkin():
    """rad.Solve(demag_backend='moment_galerkin') on a soft_iron_from_mesh container == a DIRECT
    moment_galerkin_demag_solve on the same elements + uniform applied field (machine precision),
    and writes M back so rad.ObjM reflects it.  (Uniform background -> centroid H = H0 everywhere,
    so the direct call uses H_ext=(0,0,H0).)"""
    mesh = _hex_cube_mesh(2)
    iron = _radsolve.soft_iron_from_mesh(mesh, mu_r=MU_R)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])       # free-space B whose H is H0
    cont = rad.ObjCnt([iron, bkg])
    with ng.TaskManager():
        res = rad.Solve(cont, 1e-9, 3000, 0, demag_backend="moment_galerkin")

    extracted, _ = extract_elements(mesh, allow_hex=True, allow_wedge=True)
    elems = [np.asarray(e["vertices"], float) for e in extracted]
    direct = mg.moment_galerkin_demag_solve(elems, mu_r=MU_R, H_ext=(0.0, 0.0, H0))

    scale = np.linalg.norm(direct["M"])
    rel = np.linalg.norm(np.asarray(res["M"]) - direct["M"]) / scale
    assert rel < 1e-9, f"dispatch M != direct moment_galerkin M (rel {rel:.2e})"

    # write-back: rad.ObjM reflects the dispatched per-element M
    objm = rad.ObjM(iron)
    M_obj = np.array([m for (_c, m) in objm])
    assert np.allclose(M_obj, res["M"], atol=1e-6 * scale), "ObjM write-back != dispatched M"

    # physics: the 2x2x2 cube magnetises along +z, M_z ~ chi/(1+d chi) H0 with d ~ cube demag
    Mz = float(np.mean(res["M"][:, 2]))
    pred = (MU_R - 1.0) / (1.0 + (1.0 / 3.0) * (MU_R - 1.0)) * H0
    assert 0.8 * pred < Mz < 1.2 * pred, f"mean Mz {Mz:.1f} not near cube self-consistency {pred:.1f}"


@pytest.mark.parametrize("alias", ["mg", "galerkin"])
def test_aliases_resolve_to_moment_galerkin(alias):
    """The 'mg' / 'galerkin' aliases dispatch identically to 'moment_galerkin' (set_demag_backend +
    per-call), giving the SAME field."""
    mesh = _hex_cube_mesh(2)
    probe = [0.0, 0.0, 0.1]

    def solve(backend):
        rad.UtiDelAll()
        iron = _radsolve.soft_iron_from_mesh(mesh, mu_r=MU_R)
        cont = rad.ObjCnt([iron, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
        with ng.TaskManager():
            rad.Solve(cont, 1e-9, 3000, 0, demag_backend=backend)
        return np.array(rad.Fld(iron, "b", probe))

    B_canon = solve("moment_galerkin")
    B_alias = solve(alias)
    assert np.allclose(B_alias, B_canon, rtol=0, atol=1e-18), f"alias {alias} field != canonical"
    # global set_demag_backend alias too
    assert rad.set_demag_backend(alias) == "moment_galerkin"


def test_meshless_moment_galerkin_raises():
    """A mesh-LESS soft iron (ObjHexahedron + MatLin, not soft_iron_from_mesh) has no registered element
    vertex lists, so demag_backend='moment_galerkin' raises (fail-loud, with guidance) -- it cannot
    recover vertices from a Radia handle."""
    v = [[0, 0, 0], [L, 0, 0], [L, L, 0], [0, L, 0], [0, 0, L], [L, 0, L], [L, L, L], [0, L, L]]
    h = rad.ObjHexahedron(v, [0, 0, 0])
    rad.MatApl(h, rad.MatLin(MU_R))
    cont = rad.ObjCnt([h, rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    with pytest.raises(ValueError, match="mesh-backed soft iron"):
        rad.Solve(cont, 1e-9, 1000, 0, demag_backend="moment_galerkin")


def test_nonlinear_moment_galerkin_raises():
    """moment-Galerkin is linear-only: a bh_table (nonlinear) soft iron raises with guidance to use
    'hdiv'/'yano' (No-Fallbacks -- it must not silently solve a linear approximation)."""
    mesh = _hex_cube_mesh(2)
    bh = [[0.0, 0.0], [100.0, 0.5], [1000.0, 1.5], [50000.0, 2.0]]
    iron = _radsolve.soft_iron_from_mesh(mesh, bh_table=bh)
    cont = rad.ObjCnt([iron, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    with pytest.raises(NotImplementedError, match="linear-only"):
        with ng.TaskManager():
            rad.Solve(cont, 1e-9, 3000, 0, demag_backend="moment_galerkin")


def test_moment_galerkin_same_ballpark_as_yano():
    """The moment-Galerkin dispatch stray field is the same physical ballpark as the 'yano' backend on
    the SAME container (they differ only by the dipole-3DOF vs 6DOF-MSC discretisation -- a few %, not
    a wrong answer)."""
    mesh = _hex_cube_mesh(2)
    probe = [0.0, 0.0, 0.1]

    def stray(backend):
        rad.UtiDelAll()
        iron = _radsolve.soft_iron_from_mesh(mesh, mu_r=MU_R)
        cont = rad.ObjCnt([iron, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
        with ng.TaskManager():
            rad.Solve(cont, 1e-9, 3000, 0, demag_backend=backend)
        return np.array(rad.Fld(iron, "b", probe))

    B_mg = stray("moment_galerkin")
    B_yano = stray("yano")
    rel = np.linalg.norm(B_mg - B_yano) / np.linalg.norm(B_yano)
    assert rel < 0.12, f"moment-Galerkin vs yano stray field rel {rel:.2e} too large (expected ~discretisation gap)"
    assert B_mg[2] > 0 and B_yano[2] > 0, "both backends magnetise +z"
