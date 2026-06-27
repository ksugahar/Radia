"""The unified .vol geometry path (decision 2026-06-19): a soft iron loaded from a netgen .vol FILE
solves with BOTH demag backends -- six-face surface-charge MSC and the FEEC HDiv-VIM -- selected by
set_demag_backend.  .vol is the SOLE Cubit<->NGSolve mesh interchange, so netgen owns the mesh
orientation (no hand-built-mesh boundary-winding pitfalls).  This locks:
  (1) radia.vim.soft_iron_from_vol(path) round-trips a .vol into a both-backend soft-iron container;
  (2) the HDiv path (default) and the collocation MMMM path (set_demag_backend('collocation_mmmm')) BOTH solve on it and
      agree to within the RT0-vs-MSC discretization gap on a structured hex cube (~few %)."""
import math

import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

MU0 = 4.0e-7 * math.pi
L = 0.02
MU_R = 100.0
H0 = 1000.0


def _make_vol(path):
    """A structured hex cube [0,L]^3 saved to a netgen .vol (netgen owns orientation)."""
    with ng.TaskManager():
        m = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3,
                                 mapping=lambda x, y, z: (L * x, L * y, L * z))
        m.ngmesh.Save(str(path))


def _solve_from_vol(path, backend):
    rad.UtiDelAll()
    rad.set_demag_backend(backend)
    try:
        with ng.TaskManager():
            iron = vim.soft_iron_from_vol(path, mu_r=MU_R)        # <- .vol -> both-backend iron
            bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])      # uniform Bz = mu0*H0  (H = H0)
            res = rad.Solve(rad.ObjCnt([iron, bkg]), 1e-6, 2000, 0)
        if isinstance(res, dict) and "M_avg" in res:               # HDiv path returns the solve dict
            return res["M_avg"][2]
        objm = rad.ObjM(iron)                                      # collocation MMMM path -> read back via ObjM
        return sum(m[2] for (_c, m) in objm) / len(objm)
    finally:
        rad.UtiDelAll()
        rad.set_demag_backend("auto")


def test_vol_solves_with_both_backends(tmp_path):
    vol = tmp_path / "cube_hex.vol"
    _make_vol(vol)

    mz_hdiv = _solve_from_vol(vol, "hdiv")
    mz_collocation = _solve_from_vol(vol, "collocation_mmmm")

    # chi = mu_r-1 = 99, cube demag ~1/3 -> M_avg_z ~ 99*H0/(1+99/3) ~ 2912; the discretized non-uniform
    # M_avg runs a bit higher.  Lock a generous physical band for BOTH backends.
    for name, mz in (("hdiv", mz_hdiv), ("collocation_mmmm", mz_collocation)):
        assert 2.5e3 < mz < 4.0e3, f"{name} from .vol gave unphysical M_avg_z={mz:.1f}"

    rel = abs(mz_hdiv - mz_collocation) / abs(mz_collocation)
    assert rel < 0.06, (
        f"HDiv {mz_hdiv:.1f} vs collocation MMMM {mz_collocation:.1f} "
        f"from same .vol differ {rel*100:.1f}% (> RT0-vs-MSC)"
    )
