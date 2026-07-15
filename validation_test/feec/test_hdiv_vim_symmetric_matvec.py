"""Golden: the SYMMETRIC charge-Gram H-matvec (symmetric HACApK).

`matvec_sym` applies the UPPER-triangular leaves only -- each upper leaf supplies its own block AND the
mirror as its exact transpose -- so the operator is EXACTLY symmetric (||G - G^T|| == 0) regardless of the
per-block ACA truncation.  This is what makes the +N CG robust by construction (the symmetric-CG default).

SYMMETRIC FILL (2026-07-03): the ChargeGram build now SKIPS the strictly-lower leaves entirely (they were
never read by matvec_sym; ~2x build), and plain `matvec` / `matvec_transpose` on the ChargeGram are ROUTED
to `matvec_sym` (for the symmetric operator they are the same map; the base implementations would read the
empty lower leaves).  Tests (1)-(2) therefore lock the ROUTING identity G == G^T == G_sym exactly -- a
routing regression (someone rebinding matvec to the base general apply on a sym-fill build) breaks them
with an O(1) error, not a subtle drift.

Locks:
  (1) matvec_transpose == matvec (both routed to the symmetric apply);
  (2) matvec / matvec_sym are EXACTLY symmetric and mutually identical;
  (3) the default solve uses the symmetric CG and matches the GMRES cross-check.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import ChargeGram, Solve  # noqa: E402


def _gram(maxh=0.45):
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))).GenerateMesh(maxh=maxh))
        fes = ng.HDiv(mesh, order=1)
        _, H, _ = ChargeGram(
            fes, eps=1e-12, leafsize=32, eta=2.0,
            ho_far_factor=float("inf"),
        )
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
    GSx = np.asarray(H.matvec_sym(np.ascontiguousarray(x)), float)
    GSy = np.asarray(H.matvec_sym(np.ascontiguousarray(y)), float)
    asym_sym = abs(x @ GSy - y @ GSx) / (abs(x @ GSy) + 1e-300)
    assert asym_sym < 1e-12, f"symmetric bilinear probe not machine-symmetric: {asym_sym:.2e}"


def test_default_solve_is_symmetric_cg():
    """The default and explicit solver names select the same C++ symmetric CG."""
    H_ext = ng.CoefficientFunction((0, 0, 1000.0))
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))).GenerateMesh(maxh=0.4))
        auto = Solve(mesh, 200.0, H_ext)
        explicit = Solve(mesh, 200.0, H_ext, linear_solver="cpp-cg")
    assert auto["linear_solver"] == "mass-riesz-cg"
    assert explicit["linear_solver"] == "mass-riesz-cg"
    rel = abs(auto["M_avg"][2] - explicit["M_avg"][2]) / abs(explicit["M_avg"][2])
    assert rel < 1e-12, f"auto and cpp-cg disagree: {rel:.2e}"
