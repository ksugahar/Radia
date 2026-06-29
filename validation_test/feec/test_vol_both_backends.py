"""The unified .vol geometry path (decision 2026-06-19): a soft iron loaded from a netgen .vol FILE
solves with BOTH demag backends -- six-face surface-charge MSC and the FEEC HDiv-VIM -- selected by
set_demag_backend.  .vol is the SOLE Cubit<->NGSolve mesh interchange, so netgen owns the mesh
orientation (no hand-built-mesh boundary-winding pitfalls).  This locks:
  (1) radia.vim.soft_iron_from_vol(path) round-trips a .vol into a soft-iron container;
  (2) a HEX .vol solves with the collocation MMMM backend (the surface-charge MSC, M_avg ~ the cube fixed
      point), AND the HDiv-VIM backend REFUSES the hex .vol (fail loud) -- HDiv-VIM is tet/RT1-only
      (2026-06-29), so the rad.Solve 'auto' split routes a non-tet mesh-backed iron to collocation MMMM.
  (A tet .vol would solve with BOTH backends, but the cross-backend agreement on a tet body is gated by the
  collocation-MMMM-on-tet accuracy -- a separate concern -- so this test fixes the hex routing instead.)"""
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
    """A structured HEX cube [0,L]^3 saved to a netgen .vol (netgen owns orientation)."""
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


def test_vol_hex_collocation_solves_hdiv_fails_loud(tmp_path):
    vol = tmp_path / "cube_hex.vol"
    _make_vol(vol)

    # collocation MMMM (surface-charge MSC) solves the hex .vol: chi = mu_r-1 = 99, cube demag ~1/3 ->
    # M_avg_z ~ 99*H0/(1+99/3) ~ 2912; the discretized non-uniform M_avg runs higher.  Generous physical band.
    mz_collocation = _solve_from_vol(vol, "collocation_mmmm")
    assert 2.5e3 < mz_collocation < 6.0e3, \
        f"collocation MMMM from hex .vol gave unphysical M_avg_z={mz_collocation:.1f}"

    # HDiv-VIM is tet/RT1-only -> it REFUSES the hex .vol (No-Fallbacks: a clear error naming collocation MMMM).
    with pytest.raises(ValueError, match="TET"):
        _solve_from_vol(vol, "hdiv")
