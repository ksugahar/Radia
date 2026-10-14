"""Test isolation for the HDiv/collocation demag-backend GLOBAL.

``rad.set_demag_backend(...)`` mutates a module-global (``radia._demag_backend``).
A test that forces a backend must restore it, or it LEAKS into later tests -- which
is exactly how ``validation_test/feec/`` went chronically CI-red: a retired
hex/wedge HDiv test forced ``"hdiv"`` and tried to restore it by storing the
return value of ``set_demag_backend("hdiv")``.  That function returns the NEW
value (not the previous), so the restore re-applied ``"hdiv"`` and later tests
ran with the wrong global backend.

This autouse fixture pins the global to the ``"auto"`` default around EVERY feec
test, so a forced backend can never leak across tests (defence in depth).
"""
import pytest


def hdiv_vim_dense_N_and_loops(mesh, eps=1e-9, leaf=64):
    """Build the DENSE HDiv-VIM demag operator N = B^T G B + the loop basis (ker B), from the production
    C++ charge-Gram H-matrix (the dense Python Gram path was removed).  For SMALL meshes only -- it forms
    N column-by-column via the C++ H-matvec (N v = B^T (H.matvec(B v))) and the ker(B) basis via a dense
    SVD of the sparse charge map B.  Returns (N dense, loops (n_loop x ndof), n_loop).

    This is the structural-test helper: symmetry (||N - N^T||) and loop-nullity (||N @ loop||) of the
    C++-Gram operator on a small mesh.  Callers must open `with ng.TaskManager():`."""
    import numpy as np
    import ngsolve as ng
    from radia.vim import ChargeGram

    fes = ng.HDiv(mesh, order=1)
    B, H, _ = ChargeGram(
        fes, eps=eps, leafsize=leaf, eta=2.0,
        ho_far_factor=float("inf"),
    )

    def N_apply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    ndof = int(fes.ndof)
    N = np.column_stack([N_apply(np.eye(ndof)[:, k]) for k in range(ndof)])
    # loops = right null space of the (small) dense charge map B
    Bd = B.toarray()
    sv = np.linalg.svd(Bd, compute_uv=False)
    rankB = int(np.sum(sv > 1e-9 * sv.max()))
    n_loop = ndof - rankB
    _, _, Vt = np.linalg.svd(Bd)
    loops = Vt[rankB:, :]
    return N, loops, n_loop


@pytest.fixture(autouse=True)
def _isolate_demag_backend():
    try:
        import radia as rad
    except Exception:                      # radia unavailable -> feec tests skip anyway
        yield
        return
    rad.set_demag_backend("auto")          # clean start, regardless of any prior leak
    try:
        yield
    finally:
        rad.set_demag_backend("auto")      # clean end, regardless of test outcome
