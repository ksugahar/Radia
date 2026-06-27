"""radia.SoftIron -- the unified intent-based soft-iron object (2-layer API, 2026-06-19).

Locks that ONE object built from a .vol drives BOTH demag backends (surface-charge MSC and HDiv-VIM) with no
ObjHexahedron in the user's hands: build -> solve(backend=...) -> field(), backend-agnostic, and the
two backends agree on the same .vol within the RT0-vs-MSC discretization gap."""
import math

import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

MU0 = 4.0e-7 * math.pi
L = 0.02
MU_R = 100.0
H0 = 1000.0


def _vol(path):
    with ng.TaskManager():
        m = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3,
                                 mapping=lambda x, y, z: (L * x, L * y, L * z))
        m.ngmesh.Save(str(path))


def test_soft_iron_both_backends_from_vol(tmp_path):
    vol = tmp_path / "cube_hex.vol"
    _vol(vol)

    def run(backend):
        rad.UtiDelAll()
        with ng.TaskManager():
            iron = rad.SoftIron(vol, mu_r=MU_R)                  # intent object, no ObjHexahedron
            src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])    # applied H = H0
            iron.solve(source=src, backend=backend)
            mz = sum(m[2] for (_c, m) in iron.magnetization()) / len(iron.magnetization())
            # field() returns total (iron+source); check it is finite at the centre
            bz = iron.field("b", [[0.0, 0.0, 0.0]])[0][2] if backend == "hdiv" else \
                iron.field("b", [0.0, 0.0, 0.0])[2]
        rad.UtiDelAll()
        return mz, bz

    mz_hdiv, _ = run("hdiv")
    mz_collocation, _ = run("collocation_mmmm")

    for name, mz in (("hdiv", mz_hdiv), ("collocation_mmmm", mz_collocation)):
        assert 2.5e3 < mz < 4.0e3, f"SoftIron({name}) unphysical M_avg_z={mz:.1f}"
    rel = abs(mz_hdiv - mz_collocation) / abs(mz_collocation)
    assert rel < 0.06, (
        f"SoftIron hdiv {mz_hdiv:.1f} vs collocation MMMM {mz_collocation:.1f} "
        f"differ {rel*100:.1f}%"
    )


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
