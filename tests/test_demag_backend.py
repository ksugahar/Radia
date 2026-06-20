"""The demag-backend API: BOTH yano-type MSC and the FEEC HDiv-VIM are kept (decision 2026-06-19).

Default is "auto" (API-split): mesh-LESS hex/wedge soft iron (ObjHexahedron/ObjWedge + MatLin) is
solved by the yano-type collocation MSC demag; mesh-BACKED soft iron (radia.vim.soft_iron_from_mesh)
is solved by the FEEC HDiv-VIM.  set_demag_backend("yano"|"hdiv") overrides; "auto"/None restores the
split.  Tet (MMM) and permanent-magnet solves are unaffected.  The mesh-backed HDiv routing is locked
by tests/feec/test_hdiv_radsolve_dispatch.py."""
import math

import pytest
import radia as rad

MU0 = 4.0e-7 * math.pi


def test_backend_default_is_auto():
    rad.set_demag_backend("auto")
    assert rad.get_demag_backend() == "auto"


def test_set_yano_and_hdiv_accepted():
    assert rad.set_demag_backend("yano") == "yano"
    assert rad.get_demag_backend() == "yano"
    assert rad.set_demag_backend("hdiv") == "hdiv"
    assert rad.get_demag_backend() == "hdiv"
    assert rad.set_demag_backend("auto") == "auto"
    rad.set_demag_backend("auto")


def test_invalid_backend_raises():
    with pytest.raises(ValueError):
        rad.set_demag_backend("bogus")
    rad.set_demag_backend("auto")


def test_solverconfig_both_ok():
    rad.SolverConfig(demag_backend="hdiv")
    rad.SolverConfig(demag_backend="yano")
    rad.SolverConfig(demag_backend="auto")
    rad.set_demag_backend("auto")


def test_meshless_hex_soft_iron_solves_via_yano():
    """A hex soft iron built the mesh-less way (ObjHexahedron + MatLin) is solved by the yano-type MSC
    demag (no NGSolve needed).  A cube in a uniform applied Hz magnetizes with demag ~1/3, so for
    mu_r=1000 the magnetization M_z ~ H0/(1/3) = 3*H0.  This locks that the yano path is REACHABLE
    (no Error203) and physical."""
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
    assert 2.0 * H0 < M[2] < 4.0 * H0, f"yano-MSC cube M_z={M[2]:.1f} not ~3*H0={3*H0}"
    assert abs(M[0]) < 0.05 * abs(M[2]) and abs(M[1]) < 0.05 * abs(M[2])
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
