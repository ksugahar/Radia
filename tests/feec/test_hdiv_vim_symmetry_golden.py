"""Golden test for the symmetric HDiv-type VIM demag operator (Radia C++ core, rad_hdiv_vim).

Locks the structural foundation of the HDiv-type VIM (the symmetric alternative to the collocation
MSC kernel): on a structured nx*ny*nz hex grid (RT0 faces), the demag operator N = B^T G B is
  (1) SYMMETRIC (Galerkin energy form) -- ||N - N^T||/||N|| ~ machine eps, and
  (2) loops are FIELD-NULL BY CONSTRUCTION (loops = ker B; B.loop = 0 => N.loop = 0).
Golden values come from the NGSolve prototype (examples/feec_vim/hdiv_demag_quad_self.json):
regular 3x3x3 -> ndof=108, n_loop=28, asym~1e-16, loop_res~1e-16.  The C++ hand-enumerated topology
(rad_hdiv_vim, no NGSolve) reproduces these exactly (verified standalone before integration).
"""
import numpy as np
import pytest

import radia._radia_pybind as _rp


def _assemble(nx, ny, nz):
    d = _rp._hdiv_vim_assemble(nx, ny, nz)
    nf, n_charge = d["nf"], d["n_charge"]
    N = np.asarray(d["N"], float).reshape(nf, nf)
    B = np.asarray(d["B"], float).reshape(n_charge, nf)
    return d, N, B


@pytest.mark.parametrize("nx,ny,nz,golden_ndof,golden_nloop", [
    (1, 1, 1, 6, 0),
    (2, 2, 2, 36, 5),
    (3, 3, 3, 108, 28),   # the NGSolve prototype golden (hdiv_demag_quad_self.json)
    (4, 4, 4, 240, 81),   # scaling table (hdiv_demag_scaling.json)
    (5, 5, 5, 450, 176),
])
def test_hdiv_vim_symmetry_and_loop_nullity(nx, ny, nz, golden_ndof, golden_nloop):
    d, N, B = _assemble(nx, ny, nz)
    nf = d["nf"]
    assert nf == golden_ndof, f"ndof {nf} != golden {golden_ndof}"

    # (1) symmetry (Galerkin energy form)
    nN = np.linalg.norm(N, 2)
    asym = np.linalg.norm(N - N.T) / nN
    assert asym < 1e-12, f"N not symmetric: ||N-N^T||/||N|| = {asym:.2e}"

    # (2) loops = ker(B); count + field-nullity
    s = np.linalg.svd(B, compute_uv=False)
    rankB = int(np.sum(s > 1e-9 * s.max()))
    n_loop = nf - rankB
    assert n_loop == golden_nloop, f"n_loop {n_loop} != golden {golden_nloop}"

    if n_loop:
        _, _, Vt = np.linalg.svd(B)
        loops = Vt[rankB:, :]                      # ker(B) basis
        loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / nN
        assert loop_res < 1e-10, f"loops not field-null: max||N.loop||/||N|| = {loop_res:.2e}"


def test_hdiv_vim_charge_conservation_gauss():
    """rank(B) = n_charge - 1 (Gauss: the total-charge mode is the single redundancy)."""
    d, N, B = _assemble(3, 3, 3)
    s = np.linalg.svd(B, compute_uv=False)
    rankB = int(np.sum(s > 1e-9 * s.max()))
    assert rankB == d["n_charge"] - 1, f"rank(B)={rankB} != n_charge-1={d['n_charge']-1}"
