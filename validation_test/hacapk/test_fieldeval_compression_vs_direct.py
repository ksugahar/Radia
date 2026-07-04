"""Validation: the embed-in-square H-matrix field evaluation (`_FieldEvalHMatrix`) reproduces the direct
`rad.Fld('b', ...)` and its ACA compression accuracy is controlled by the tolerance eps.

`_FieldEvalHMatrix` (src/core/rad_hacapk_field.*) accelerates rad.Fld-style field evaluation by embedding
the rectangular obs x src field operator in a SQUARE SYMMETRIC HACApK H-matrix over the combined [obs; src]
point set:  A = [[0, K],[K^T, 0]],  K[3o+a, 3s+b] = a-component of B at obs_o from src_s per unit M_b.
With x = [0; M],  y = A x = [K M; 0],  so y[:3*n_obs] is the B-field at the observation points.

The kernel is bit-consistent with rad.Fld (same radTg3d::B_genComp), so:
  * `entry(i,j)` is the EXACT, uncompressed, eps-independent embed entry (Tesla), and
  * `matvec(x)`  is the COMPRESSED HACApK apply (ACA at tolerance eps).
Ground truth is the direct `rad.Fld`.  This locks three guarantees:

  1. KERNEL BIT-CONSISTENCY -- the exact entry-matrix applied to M reproduces direct rad.Fld to ~machine.
  2. ACCURACY TRACKS eps    -- the compressed-matvec error vs direct decreases with eps.
  3. COMPRESSION IS REAL    -- at a loose eps the error is well above machine zero (an uncompressed dense
     apply would be exact regardless of eps), and the geometry is genuinely MULTI-CLUSTER (n_elem >> leaf,
     the O(N log N) regime).  This is the case that a naive one-sided embed + un-normalised (O(mu0*H) ~
     1e-13) kernel silently mis-fills; the fix is the symmetric embed + O(1) entry normalisation
     (memory radfld-hmatrix-derisk).

FLAT global-coordinate container only (no ancestor-group transforms) and OUTSIDE observation points (the
field superposition of per-element responses equals rad.Fld only outside the sources).
Run explicitly (validation lane): `pytest validation_test/hacapk/test_fieldeval_compression_vs_direct.py`.
"""
import numpy as np
import pytest

import radia as rad
import radia._radia_pybind as _rp

_EPS_SWEEP = (1e-2, 1e-4, 1e-6, 1e-8)
_DIM = 0.14
_LEAF = 12          # << n_src, n_obs  -> genuinely multi-cluster (the fine-tree / O(N log N) regime)


def _geometry(n_src=3, n_obs=4, obs_lo=5.0, seed=0):
    """Flat ObjCnt of n_src^3 unit-ish ObjRecMag with known random M; n_obs^3 FAR observation points."""
    rng = np.random.default_rng(seed)
    xs = np.linspace(0.5 / n_src, 1 - 0.5 / n_src, n_src)
    objs = [rad.ObjRecMag([float(x), float(y), float(z)], [_DIM, _DIM, _DIM],
                          (rng.standard_normal(3) * 3.0e5).tolist())
            for x in xs for y in xs for z in xs]
    cont = rad.ObjCnt(objs)
    oxs = np.linspace(obs_lo, obs_lo + 1.0, n_obs)      # well-separated from the [0,1]^3 sources
    obs = np.array([[x, y, z] for x in oxs for y in oxs for z in oxs], float)
    return cont, obs


@pytest.fixture(scope="module", params=["b", "a"], ids=["B_flux_density", "A_vector_potential"])
def result(request):
    """Runs for both field_type='b' (flux density, Tesla) and 'a' (vector potential, T*m -- the A-form
    eddy-current FEM coupling source term)."""
    ft = request.param
    rad.UtiDelAll()
    cont, obs = _geometry()
    n_obs = len(obs)
    Fdir = np.asarray(rad.Fld(cont, ft, obs.tolist()), float).reshape(n_obs, 3)   # ground truth FIRST

    obs_flat = obs.reshape(-1).tolist()
    # exact entry-matrix (uncompressed) at a tight eps -> kernel bit-consistency
    Gx = _rp._FieldEvalHMatrix(cont, obs_flat, eps=1e-10, leaf=_LEAF, eta=2.0, field_type=ft)
    ndof, no, ns = Gx.ndof(), Gx.n_obs(), Gx.n_src()
    x = np.array([0.0] * (3 * no) + list(Gx.src_magnetization()))
    D = np.array([[Gx.entry(i, j) for j in range(ndof)] for i in range(ndof)])
    F_entry = (D @ x)[:3 * no].reshape(no, 3)
    err_entry = float(np.linalg.norm(F_entry - Fdir) / np.linalg.norm(Fdir))

    # compressed matvec error vs direct rad.Fld, per eps
    errs = {}
    for eps in _EPS_SWEEP:
        G = _rp._FieldEvalHMatrix(cont, obs_flat, eps=eps, leaf=_LEAF, eta=2.0, field_type=ft)
        y = np.asarray(G.matvec(x.tolist()), float)[:3 * no].reshape(no, 3)
        errs[eps] = float(np.linalg.norm(y - Fdir) / np.linalg.norm(Fdir))

    info = dict(field_type=ft, is_a=Gx.is_a_field(), ndof=ndof, n_obs=no, n_src=ns, n_elem=no + ns,
                err_entry=err_entry, errs=errs)
    rad.UtiDelAll()
    return info


def test_multi_cluster_regime(result):
    """The geometry must be genuinely multi-cluster (n_elem >> leaf) -- otherwise the H-matrix is a single
    dense block and 'compression' is vacuously exact, hiding the fine-tree fill / normalisation contract."""
    assert result["ndof"] == 3 * result["n_elem"], "ndof must be 3*(n_obs + n_src) (uniform 3-DOF embed)"
    assert result["n_elem"] > 4 * _LEAF, f"n_elem={result['n_elem']} not >> leaf={_LEAF} (need multi-cluster)"
    assert result["is_a"] == (result["field_type"] == "a"), "is_a_field() disagrees with field_type"


def test_kernel_bit_consistency(result):
    """The EXACT entry-matrix applied to M reproduces direct rad.Fld(field_type) to ~machine (the kernel is
    bit-consistent with B_genComp; the per-field-type unit scaling -- mu0 for B, unity for the already-
    physical A -- and per-element superposition are correct)."""
    assert result["err_entry"] < 1e-10, \
        f"entry-matrix vs direct rad.Fld('{result['field_type']}') = {result['err_entry']:.3e}"


def test_compression_accuracy_tracks_eps(result):
    """The compressed-matvec error vs direct decreases with eps (down to the ~1e-8 ACA/kernel floor)."""
    errs = [result["errs"][e] for e in _EPS_SWEEP]         # eps decreasing 1e-2 -> 1e-8
    for a, b, ea, eb in zip(errs, errs[1:], _EPS_SWEEP, _EPS_SWEEP[1:]):
        assert b <= a * 3.0, f"error did not decrease from eps={ea:.0e} ({a:.3e}) to eps={eb:.0e} ({b:.3e})"
    for eps, err in result["errs"].items():
        assert err < 200.0 * eps + 1e-7, f"eps={eps:.0e}: matvec error {err:.3e} not tracking eps"


def test_compression_is_real(result):
    """Loose eps error is well above machine zero (a dense apply would be exact regardless of eps -> proves
    ACA low-rank compression is actually happening on this multi-cluster H-matrix); tight eps is accurate."""
    assert result["errs"][1e-2] > 1e-5, f"loose-eps error {result['errs'][1e-2]:.3e} ~ machine (no compression?)"
    assert result["errs"][1e-8] < 1e-6, f"tight-eps error {result['errs'][1e-8]:.3e} should be near-exact"
