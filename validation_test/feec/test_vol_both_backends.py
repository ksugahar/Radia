"""The unified .vol geometry path (decision 2026-06-19): a soft iron loaded from a netgen .vol FILE
solves with BOTH demag backends -- six-face surface-charge MSC and the FEEC HDiv-VIM -- selected by
set_demag_backend.  .vol is the SOLE Cubit<->NGSolve mesh interchange, so netgen owns the mesh
orientation (no hand-built-mesh boundary-winding pitfalls).  This locks:
  (1) radia.vim.VolSoftIron(path) round-trips a .vol into a soft-iron container;
  (2) a HEX .vol solves with BOTH demag backends -- the collocation MMMM surface-charge MSC AND the FEEC
      HDiv-VIM (hex unlocked 2026-07-04: the wired hex RT1 charge Gram + the shipped mass-Riesz CG) -- and
      the two AGREE on M_avg_z (cross-method).  The production rad.Solve 'auto' default now routes pure HEX
      mesh-backed irons to HDiv-VIM; this file still selects each backend explicitly to keep the comparison
      unambiguous.
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
            iron = vim.VolSoftIron(path, mu_r=MU_R)        # <- .vol -> both-backend iron
            bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])      # uniform Bz = mu0*H0  (H = H0)
            res = rad.Solve(rad.ObjCnt([iron, bkg]), 1e-6, 2000, 0)
        if isinstance(res, dict) and "M_avg" in res:               # HDiv path returns the solve dict
            return res["M_avg"][2]
        objm = rad.ObjM(iron)                                      # collocation MMMM path -> read back via ObjM
        return sum(m[2] for (_c, m) in objm) / len(objm)
    finally:
        rad.UtiDelAll()
        rad.set_demag_backend("auto")


def test_vol_hex_solves_both_backends(tmp_path):
    vol = tmp_path / "cube_hex.vol"
    _make_vol(vol)

    # collocation MMMM (surface-charge MSC) solves the hex .vol: chi = mu_r-1 = 99, cube demag ~1/3 ->
    # M_avg_z ~ 99*H0/(1+99/3) ~ 2912; the discretized non-uniform M_avg runs higher.  Generous physical band.
    mz_collocation = _solve_from_vol(vol, "collocation_mmmm")
    assert 2.5e3 < mz_collocation < 6.0e3, \
        f"collocation MMMM from hex .vol gave unphysical M_avg_z={mz_collocation:.1f}"

    # HDiv-VIM ALSO solves the hex .vol now (hex unlocked 2026-07-04): set_demag_backend('hdiv') drives the
    # wired hex RT1 charge Gram + the shipped mass-Riesz CG.  Retry the KNOWN bursty GetTrafo first-touch flake.
    last = None
    for _ in range(5):
        try:
            mz_hdiv = _solve_from_vol(vol, "hdiv")
            break
        except RuntimeError as e:
            if "GetTrafo lattice evaluation unstable" in str(e):
                last = e
                continue
            raise
    else:
        raise last
    assert 2.5e3 < mz_hdiv < 6.0e3, f"HDiv from hex .vol gave unphysical M_avg_z={mz_hdiv:.1f}"
    # cross-method agreement on the SAME hex .vol (generous coarse-mesh 3x3x3 band)
    rel = abs(mz_hdiv - mz_collocation) / abs(mz_collocation)
    assert rel < 0.05, f"hex .vol: HDiv {mz_hdiv:.1f} vs collocation MMMM {mz_collocation:.1f} rel {rel:.2e}"
