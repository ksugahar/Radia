"""The demag-backend API: Radia soft iron uses FEEC HDiv-VIM.

Default is "auto": mesh-BACKED pure TET / HEX / WEDGE soft iron
(radia.vim.MeshSoftIron) is eligible for FEEC HDiv-VIM RT1. Permanent-magnet
field-only objects are unaffected. The
mesh-backed HDiv routing is locked by validation_test/feec/test_hdiv_radsolve_dispatch.py."""
import math

import pytest
import radia as rad

MU0 = 4.0e-7 * math.pi


def test_backend_default_is_auto():
    rad.set_demag_backend("auto")
    assert rad.get_demag_backend() == "auto"


def test_set_hdiv_accepted():
    assert rad.set_demag_backend("hdiv") == "hdiv"
    assert rad.get_demag_backend() == "hdiv"
    assert rad.set_demag_backend("auto") == "auto"
    rad.set_demag_backend("auto")


def test_invalid_backend_raises():
    with pytest.raises(ValueError):
        rad.set_demag_backend("bogus")
    rad.set_demag_backend("auto")


def test_invalid_per_call_backend_raises_before_cpp_dispatch():
    with pytest.raises(ValueError):
        rad.Solve(0, demag_backend="bogus")
    rad.set_demag_backend("auto")


def test_solverconfig_hdiv_and_auto_ok():
    rad.SolverConfig(demag_backend="hdiv")
    rad.SolverConfig(demag_backend="auto")
    rad.set_demag_backend("auto")


def test_meshless_hex_soft_iron_rejected():
    """ObjHexahedron + MatLin soft iron must use a mesh-backed HDiv route."""
    rad.set_demag_backend("auto")
    rad.UtiDelAll()
    L = 0.01
    hexv = [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L],
            [-L, -L, L], [L, -L, L], [L, L, L], [-L, L, L]]
    h = rad.ObjHexahedron(hexv, [0, 0, 0])
    rad.MatApl(h, rad.MatLin(1000.0))
    H0 = 1000.0
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])   # uniform Bz = mu0*H0
    cont = rad.ObjCnt([h, bkg])
    with pytest.raises(RuntimeError, match="[Mm]esh-less soft iron"):
        rad.Solve(cont, 1e-6, 1000, 0)
    rad.UtiDelAll()
    rad.set_demag_backend("auto")


def _mesh_backed_top(mesh, ng):
    from radia.vim import MeshSoftIron

    rad.UtiDelAll()
    with ng.TaskManager():
        iron = MeshSoftIron(mesh, mu_r=20.0)
    return rad.ObjCnt([iron])


def test_mesh_backed_hex_is_hdiv_auto_eligible():
    """Fast routing contract: a mesh-backed pure HEX iron is now an HDiv-VIM auto candidate.

    This stays in tests/ rather than validation_test/ because it does not solve; it only guards the cheap
    dispatch predicate that decides whether ``rad.Solve(..., demag_backend=None)`` enters the HDiv bridge.
    The numerical solve parity lives in validation_test/feec.
    """
    ng = pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim import _radsolve

    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    with ng.TaskManager():
        hex_mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1, mapping=mp)
    assert _radsolve.is_hdiv_eligible(_mesh_backed_top(hex_mesh, ng))

    rad.UtiDelAll()
    rad.set_demag_backend("auto")


def test_mesh_backed_wedge_is_hdiv_auto_eligible():
    """Fast routing contract: a mesh-backed pure WEDGE iron is now an HDiv-VIM auto candidate."""
    ng = pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim import _radsolve

    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    try:
        with ng.TaskManager():
            wedge_mesh = MakeStructured3DMesh(prism=True, nx=1, ny=1, nz=1, mapping=mp)
    except TypeError:
        pytest.skip("this NGSolve MakeStructured3DMesh has no prism= kwarg")
    if {len(el.vertices) for el in wedge_mesh.Elements(ng.VOL)} != {6}:
        pytest.skip("prism mesh did not produce pure wedges")
    assert _radsolve.is_hdiv_eligible(_mesh_backed_top(wedge_mesh, ng))

    rad.UtiDelAll()
    rad.set_demag_backend("auto")


def test_multiple_mesh_backed_irons_are_hdiv_auto_eligible():
    ng = pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim import MeshSoftIron, _radsolve

    rad.UtiDelAll()
    meshes = [
        MakeStructured3DMesh(
            hexes=True, nx=1, ny=1, nz=1,
            mapping=lambda x, y, z, shift=shift: (
                0.01*(x-0.5)+shift, 0.01*(y-0.5), 0.01*(z-0.5)))
        for shift in (-0.02, 0.02)
    ]
    with ng.TaskManager():
        irons = [MeshSoftIron(mesh, mu_r=20.0) for mesh in meshes]
    top = rad.ObjCnt(irons)

    assert _radsolve.registered_iron_count(top) == 2
    assert _radsolve.is_hdiv_eligible(top)

    rad.UtiDelAll()
    rad.set_demag_backend("auto")
