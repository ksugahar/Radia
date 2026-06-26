"""Test isolation for the HDiv/yano demag-backend GLOBAL.

``rad.set_demag_backend(...)`` mutates a module-global (``radia._demag_backend``).
A test that forces a backend must restore it, or it LEAKS into later tests -- which
is exactly how ``validation_test/feec/`` went chronically CI-red:

  ``test_hdiv_vim_hex_wedge`` forced ``"hdiv"`` and "restored" it with
  ``prev = rad.set_demag_backend("hdiv"); finally: rad.set_demag_backend(prev)`` --
  but ``set_demag_backend`` returns the NEW value (not the previous), so ``prev``
  was ``"hdiv"`` and the restore RE-APPLIED ``"hdiv"``.  pytest runs ``hex_wedge``
  before ``newton_vs_radia`` + ``pm_mixing`` (alphabetical), which then ran
  ``rad.Solve`` (default -> leaked ``"hdiv"``) on a raw-tet iron and hit
  ``demag_backend='hdiv' needs a mesh-backed soft iron``.

This autouse fixture pins the global to the ``"auto"`` default around EVERY feec
test, so a forced backend can never leak across tests (defence in depth -- the
hex_wedge save/restore is also corrected to capture the previous value directly).
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
    import radia._radia_pybind as _rp
    from radia.vim import _core as _tet

    d = _tet.build_demag(mesh)
    B = d["B_csr"]
    n_el = int(d["n_el"])
    all_tet = len(d["cell_verts"]) > 0 and len(d["face_verts"]) > 0 and d["poly"] is None
    if all_tet:
        H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                   n_el=n_el, eps=eps, leaf=leaf, eta=2.0)
    else:
        p = d["poly"]
        H = _rp._ChargeGramHMatrix(
            cell_tris=list(p["cell_tris"]), cell_troff=list(p["cell_troff"]),
            cell_cent=list(p["cell_cent"]), cell_meas=list(p["cell_meas"]),
            face_tris=list(p["face_tris"]), face_troff=list(p["face_troff"]),
            face_cent=list(p["face_cent"]), face_meas=list(p["face_meas"]),
            n_el=n_el, eps=eps, leaf=leaf, eta=2.0)

    def N_apply(v):
        v = np.asarray(v, float)
        return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)

    ndof = int(d["ndof"])
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
