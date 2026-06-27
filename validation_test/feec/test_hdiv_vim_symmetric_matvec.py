"""Golden: the SYMMETRIC charge-Gram H-matvec (symmetric HACApK).

The HACApK H-matrix stores both the (I,J) and (J,I) admissible leaf blocks but ACA-truncates them
INDEPENDENTLY, so the GENERAL matvec (`matvec`) is only approximately symmetric.  `matvec_sym` applies the
UPPER-triangular leaves only -- each upper leaf supplies its own block AND the mirror as its exact transpose
-- so the operator is EXACTLY symmetric (||G - G^T|| == 0) regardless of the per-block truncation.  This is
what makes the +N CG robust by construction (the symmetric-CG default), replacing the earlier GMRES retreat.

Locks:
  (1) matvec_transpose is the exact transpose of matvec;
  (2) matvec_sym is EXACTLY symmetric (G_sym == G_sym^T to machine precision) and == 0.5*(G + G^T)
      (which proves the upper-triangular leaf partition is symmetric -- the construction's premise);
  (3) the default solve uses the symmetric CG and matches the GMRES cross-check.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

import radia._radia_pybind as rp  # noqa: E402
from radia.vim._core import build_demag  # noqa: E402
from radia.vim import hdiv_demag_solve  # noqa: E402


def _gram(maxh=0.45, near_factor=2.0, far_quad=4):
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))).GenerateMesh(maxh=maxh))
        d = build_demag(mesh)
    H = rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                              n_el=int(d["n_el"]), eps=1e-12, leaf=32, eta=2.0,
                              near_factor=near_factor, image_masks=[], image_signs=[], far_quad=int(far_quad))
    return H


def _materialize(apply_fn, n):
    G = np.zeros((n, n)); e = np.zeros(n)
    for j in range(n):
        e[:] = 0.0; e[j] = 1.0
        G[:, j] = np.asarray(apply_fn(e.tolist()), float)
    return G


def test_transpose_matvec_is_exact_transpose():
    """matvec_transpose(x) == G^T x for the general H-matrix (mirror leaf apply)."""
    H = _gram()
    n = H.ndof()
    G = _materialize(H.matvec, n)
    GT = _materialize(H.matvec_transpose, n)
    err = np.linalg.norm(GT - G.T) / np.linalg.norm(G)
    assert err < 1e-12, f"transpose-matvec err {err:.2e}"


def test_matvec_sym_exactly_symmetric_and_equals_average():
    """matvec_sym is EXACTLY symmetric, and equals 0.5*(G + G^T) -- the latter PROVES the upper-triangular
    leaf partition is symmetric (else upper-only leaves would not reproduce the symmetrized general Gram)."""
    H = _gram()
    n = H.ndof()
    G = _materialize(H.matvec, n)
    GS = _materialize(H.matvec_sym, n)
    asym_general = np.linalg.norm(G - G.T) / np.linalg.norm(G)
    asym_sym = np.linalg.norm(GS - GS.T) / np.linalg.norm(GS)
    diff_avg = np.linalg.norm(GS - 0.5 * (G + G.T)) / np.linalg.norm(G)
    assert asym_sym < 1e-12, f"matvec_sym not symmetric: {asym_sym:.2e}"
    assert diff_avg < 1e-12, f"matvec_sym != 0.5(G+G^T): {diff_avg:.2e} (leaf partition not symmetric?)"
    # the general matvec is at most mildly asymmetric (ACA), but matvec_sym is machine-exact regardless
    assert asym_sym <= asym_general + 1e-15


def test_matvec_sym_symmetry_bilinear_probe():
    """Bilinear-form symmetry probe (cheap, scales to large N where materialization is infeasible):
    x^T (G_sym y) == y^T (G_sym x) to machine precision.  (The general matvec's asymmetry is config-
    dependent -- already machine-level on well-resolved small meshes, ~1e-9 with the monopole far at large
    N -- so we lock only the durable contract: matvec_sym is machine-symmetric for ANY config.)"""
    H = _gram(maxh=0.3)
    n = H.ndof()
    x = np.cos(np.arange(n) * 0.7); y = np.sin(np.arange(n) * 1.3)
    GSx = np.asarray(H.matvec_sym(x.tolist()), float); GSy = np.asarray(H.matvec_sym(y.tolist()), float)
    asym_sym = abs(x @ GSy - y @ GSx) / (abs(x @ GSy) + 1e-300)
    assert asym_sym < 1e-12, f"symmetric bilinear probe not machine-symmetric: {asym_sym:.2e}"


def test_default_solve_is_symmetric_cg():
    """The default hdiv_demag_solve uses the symmetric mass-Riesz CG (the symmetric-HACApK matvec) and
    matches the GMRES cross-check on the uniform-field cube."""
    H_ext = ng.CoefficientFunction((0, 0, 1000.0))
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))).GenerateMesh(maxh=0.4))
        auto = hdiv_demag_solve(mesh, 200.0, H_ext)
        gm = hdiv_demag_solve(mesh, 200.0, H_ext, linear_solver="gmres")
    assert auto["linear_solver"] == "mass-riesz-cg"
    assert gm["linear_solver"] == "mass-riesz-gmres"
    rel = abs(auto["M_avg"][2] - gm["M_avg"][2]) / abs(gm["M_avg"][2])
    assert rel < 1e-6, f"symmetric CG vs GMRES disagree: {rel:.2e}"
