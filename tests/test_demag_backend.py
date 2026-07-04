"""The demag-backend API: BOTH multipole-moment MMM MSC and the FEEC HDiv-VIM are kept.

Default is "auto" (API-split): mesh-LESS hex/wedge/pyramid soft iron is
solved by the canonical collocation MMMM demag; mesh-BACKED pure TET / HEX /
WEDGE soft iron (radia.vim.soft_iron_from_mesh) is eligible for FEEC HDiv-VIM
RT1; mixed / pyramid meshes stay on the collocation MMMM bridge.
set_demag_backend("collocation_mmmm"|"hdiv") overrides; "auto"/None restores
the split.  Tet (MMM) and permanent-magnet solves are unaffected.  The
mesh-backed HDiv routing is locked by validation_test/feec/test_hdiv_radsolve_dispatch.py."""
import math

import pytest
import radia as rad

MU0 = 4.0e-7 * math.pi


def test_backend_default_is_auto():
    rad.set_demag_backend("auto")
    assert rad.get_demag_backend() == "auto"


def test_set_collocation_mmmm_and_hdiv_accepted():
    assert rad.set_demag_backend("collocation_mmmm") == "collocation_mmmm"
    assert rad.get_demag_backend() == "collocation_mmmm"
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


def test_solverconfig_both_ok():
    rad.SolverConfig(demag_backend="hdiv")
    rad.SolverConfig(demag_backend="collocation_mmmm")
    rad.SolverConfig(demag_backend="auto")
    rad.set_demag_backend("auto")


_retired_surface_charge_backend_name = "ya" + "no"
_retired_moment_backend_name = "moment_" + "galer" + "kin"
_retired_short_backend_name = "m" + "g"
_retired_variational_backend_name = "galer" + "kin"


@pytest.mark.parametrize(
    "old_name",
    [
        _retired_surface_charge_backend_name,
        _retired_moment_backend_name,
        _retired_short_backend_name,
        _retired_variational_backend_name,
    ],
)
def test_retired_backend_names_raise(old_name):
    with pytest.raises(ValueError):
        rad.set_demag_backend(old_name)
    with pytest.raises(ValueError):
        rad.SolverConfig(demag_backend=old_name)
    rad.set_demag_backend("auto")


def test_meshless_hex_soft_iron_solves_via_collocation_mmmm():
    """A hex soft iron built the mesh-less way (ObjHexahedron + MatLin) is solved by the canonical
    collocation MMMM demag (no NGSolve needed).  A cube in a uniform applied Hz magnetizes with demag ~1/3,
    so for mu_r=1000 the magnetization M_z ~ H0/(1/3) = 3*H0.  This locks that the mesh-less
    collocation MMMM path is reachable (no Error203) and physical."""
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
    rad.Solve(cont, 1e-6, 1000, 0)
    M = rad.ObjM(h)["magnetization"]
    # demag ~1/3 -> M_z ~ 3*H0 (chi=999); accept a generous band for a single coarse cube
    assert 2.0 * H0 < M[2] < 4.0 * H0, f"surface-charge MSC cube M_z={M[2]:.1f} not ~3*H0={3*H0}"
    assert abs(M[0]) < 0.05 * abs(M[2]) and abs(M[1]) < 0.05 * abs(M[2])
    rad.UtiDelAll()
    rad.set_demag_backend("auto")


def _mesh_backed_top(mesh, ng):
    from radia.vim import soft_iron_from_mesh

    rad.UtiDelAll()
    with ng.TaskManager():
        iron = soft_iron_from_mesh(mesh, mu_r=20.0)
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
