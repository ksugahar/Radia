"""Golden test: the symmetric HDiv-type VIM demag operator N = B^T G B built as a HACApK H-matrix
(rad_hacapk_hdiv, RadHACApKHDivManager).  This is the PRODUCTION step that makes the symmetric
operator SCALE -- the dense O(n^2) N (and its O(n_charge^2 nsub^6) assembly) is replaced by an
H-matrix whose ACA entry function evaluates N[i][j] on demand as the charge-cluster Coulomb sum.

The C++ probe (_hdiv_vim_hmatrix_probe) builds the H-matrix AND the dense reference
(rad_hdiv::AssembleN, same nsub) and returns the comparison:
  - matvec_relerr  : max over deterministic probe vectors of ||H x - N_dense x|| / ||N_dense x||
  - symmetry_relerr: |x^T (N y) - y^T (N x)| from the H-matvec
plus the H-matrix stats (compression, n_lowrank, ...).

Three locks:
  (1) when the grid is small enough that ALL leaves are dense (no ACA compression), the H-matvec
      must equal dense N to MACHINE precision -- this proves the entry function + matvec plumbing
      are EXACT, independent of any ACA tolerance;
  (2) on larger grids (with compressed far blocks) the H-matvec matches dense N within a few times
      the ACA tolerance, on regular AND distorted grids, at nsub=0 (fast Gram) and nsub>=1 (accurate);
  (3) the operator is symmetric (the entries are symmetric; ACA preserves it to its tolerance).
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp


def _probe(nx, ny, nz, nsub=0, distort=0.0, eps=1e-4, leaf=32, eta=2.0):
    d = _rp._hdiv_vim_hmatrix_probe(nx, ny, nz, nsub, distort, eps, leaf, eta)
    assert d["ok"], "H-matrix build failed"
    return d


def test_hmatrix_exact_when_all_dense():
    """2x2x2 (ndof=36) with leaf=64 -> a single dense leaf, no ACA compression -> the H-matvec
    must reproduce dense N to MACHINE precision.  Isolates the entry function + matvec plumbing
    from any ACA approximation."""
    d = _probe(2, 2, 2, nsub=0, leaf=64)
    assert d["n_lowrank"] == 0, f"expected all-dense (n_lowrank=0), got {d['n_lowrank']}"
    assert d["matvec_relerr"] < 1e-10, f"all-dense H-matvec not exact: relerr={d['matvec_relerr']:.3e}"
    assert d["symmetry_relerr"] < 1e-10, f"all-dense not symmetric: {d['symmetry_relerr']:.3e}"


@pytest.mark.parametrize("nx,ny,nz,nsub,distort", [
    (3, 3, 3, 0, 0.0),     # the prototype golden size, fast Gram
    (3, 3, 3, 2, 0.0),     # accurate sub-point Gram (entry function reproduces it on demand)
    (4, 4, 4, 0, 0.0),
    (5, 5, 5, 0, 0.0),
    (4, 4, 4, 0, 0.18),    # distorted grid (trilinear sub-points via the shared Gram helpers)
    (4, 4, 4, 2, 0.18),
])
def test_hmatrix_matches_dense_within_aca_tol(nx, ny, nz, nsub, distort):
    """The H-matvec matches dense N = B^T G B within a few times the ACA tolerance (eps=1e-5),
    on regular AND distorted grids, at nsub=0 (centroid-monopole) and nsub>=1 (accurate)."""
    eps = 1e-5
    d = _probe(nx, ny, nz, nsub=nsub, distort=distort, eps=eps, leaf=32)
    assert d["matvec_relerr"] < 100 * eps, (
        f"H-matvec relerr {d['matvec_relerr']:.3e} >> ACA eps {eps:.0e} "
        f"(nx={nx} nsub={nsub} distort={distort})")
    assert d["symmetry_relerr"] < 100 * eps, (
        f"H-matvec not symmetric: {d['symmetry_relerr']:.3e} (nx={nx} nsub={nsub} distort={distort})")


def test_hmatrix_compression_recorded():
    """Record the H-matrix compression on a scaling sweep (nsub=0, fast Gram).  The structured cube
    is a COMPACT geometry (near-field heavy; CLAUDE.md notes most blocks stay dense for such), so we
    assert only the sanity bound (compression <= 1) and a growing DOF count, and print the actual
    numbers -- the hard correctness is locked by the two tests above."""
    rows = []
    for n in (4, 6, 8):
        d = _probe(n, n, n, nsub=0, eps=1e-4, leaf=32)
        rows.append((d["n_dof"], d["n_leaves"], d["n_lowrank"], d["compression"], d["build_time"]))
        print(f"  {n}^3: ndof={d['n_dof']:5d} leaves={d['n_leaves']:4d} "
              f"lowrank={d['n_lowrank']:4d} compression={d['compression']:.3f} "
              f"build={d['build_time']:.2f}s")
        assert d["compression"] <= 1.0 + 1e-9, f"compression > 1: {d['compression']:.4f}"
        assert d["matvec_relerr"] < 1e-2, f"relerr too large at {n}^3: {d['matvec_relerr']:.3e}"
    assert rows[-1][0] > rows[0][0], "DOF count should grow with grid size"
