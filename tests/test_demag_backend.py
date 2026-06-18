"""The demag-backend API after the yano-type MSC removal.

The FEEC HDiv-VIM (radia.vim) is the ONLY soft-iron demag backend: get_demag_backend() is always
'hdiv', set_demag_backend('yano') is rejected, and rad.Solve on a mesh-less hexahedral / wedge soft
iron (the old yano path) raises Radia::Error203 directing to radia.vim.soft_iron_from_mesh.
Tetrahedron (MMM) and permanent-magnet solves are unaffected.  The auto-routing of a mesh-backed soft
iron to the HDiv-VIM is locked by tests/feec/test_hdiv_radsolve_dispatch.py."""
import pytest
import radia as rad


def test_backend_is_hdiv():
    assert rad.get_demag_backend() == "hdiv"


def test_set_yano_rejected():
    with pytest.raises(NotImplementedError):
        rad.set_demag_backend("yano")


def test_set_hdiv_is_noop():
    assert rad.set_demag_backend("hdiv") == "hdiv"
    assert rad.get_demag_backend() == "hdiv"


def test_invalid_backend_raises():
    with pytest.raises(ValueError):
        rad.set_demag_backend("bogus")


def test_solverconfig_hdiv_ok_yano_rejected():
    rad.SolverConfig(demag_backend="hdiv")          # no-op, must not raise
    with pytest.raises(NotImplementedError):
        rad.SolverConfig(demag_backend="yano")


def test_meshless_hex_soft_iron_raises():
    """A hex soft iron built the mesh-less way (ObjHexahedron + MatLin) -> rad.Solve refuses it
    (Radia::Error203: the yano-type MSC was removed) directing to soft_iron_from_mesh.  No NGSolve
    needed -- the guard is in the C++ solve entry."""
    rad.UtiDelAll()
    hexv = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
    h = rad.ObjHexahedron([[c * 0.005 for c in v] for v in hexv], [0, 0, 0])
    rad.MatApl(h, rad.MatLin(1000.0))
    cont = rad.ObjCnt([h])
    with pytest.raises(RuntimeError, match="yano-type MSC"):
        rad.Solve(cont, 1e-5, 100, 0)
    rad.UtiDelAll()
