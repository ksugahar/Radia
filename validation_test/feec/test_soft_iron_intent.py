"""radia.SoftIron -- the unified intent-based soft-iron object (2-layer API, 2026-06-19).

Locks that ONE object built from a .vol drives HDiv-VIM with no ObjHexahedron in the user's hands:
build -> solve() -> field()."""
import math

import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

MU0 = 4.0e-7 * math.pi
L = 0.02
MU_R = 100.0
H0 = 1000.0


def _vol(path):
    """Write the iron cube as a ``.vol``: this test's subject is SoftIron's FILE path.

    The mesh comes from OCC rather than ``MakeStructured3DMesh`` because saving a structured mesh
    crashes NGSolve 6.2.2606 (access violation inside ``ngmesh.Save``, reproducible with no radia
    imported at all, for HEX and TET alike, while the OCC and CSG routes save and reload fine).
    """
    with ng.TaskManager():
        geometry = OCCGeometry(Box(Pnt(0.0, 0.0, 0.0), Pnt(L, L, L)))
        geometry.GenerateMesh(maxh=L / 3.0).Save(str(path))


def test_soft_iron_hdiv_from_vol(tmp_path):
    vol = tmp_path / "cube_hex.vol"
    _vol(vol)

    def run(backend="hdiv"):
        rad.UtiDelAll()
        with ng.TaskManager():
            iron = rad.SoftIron(vol, mu_r=MU_R)                  # intent object, no ObjHexahedron
            src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])    # applied H = H0
            iron.solve(source=src, backend=backend)
            mz = sum(m[2] for (_c, m) in iron.magnetization()) / len(iron.magnetization())
            # field() returns total (iron+source); check it is finite at the centre
            bz = iron.field("b", [[0.0, 0.0, 0.0]])[0][2]
        rad.UtiDelAll()
        return mz, bz

    mz_hdiv, bz = run()
    assert 2.5e3 < mz_hdiv < 4.0e3, f"SoftIron(hdiv) unphysical M_avg_z={mz_hdiv:.1f}"
    assert abs(bz) > 0.0


def test_soft_iron_repr_and_auto(tmp_path):
    vol = tmp_path / "c.vol"
    _vol(vol)
    rad.UtiDelAll()
    with ng.TaskManager():
        iron = rad.SoftIron(vol, mu_r=MU_R)
        assert "SoftIron" in repr(iron) and "mu_r=100" in repr(iron)
        # auto on a mesh-backed iron -> HDiv-VIM (returns the solve dict)
        src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
        res = iron.solve(source=src, backend="auto")
    assert isinstance(res, dict) and "M_avg" in res, "auto on a .vol iron should dispatch HDiv (dict result)"
    rad.UtiDelAll()
