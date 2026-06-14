"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r, H_ext, ...)` is the first shippable HDiv-type VIM Radia API: the LINEAR
soft-iron applied-field demag solve, the candidate replacement for the yano-type MSC hex/wedge
soft-iron demag in `rad.Solve`.  It wraps `build_demag` + the C++ `_ChargeGramHMatrix.solve_linear_material`
(the PHYSICAL +N material system) and returns per-element M.

## Formulation (verified-first, 2026-06-15)
The physical linear-material system is the **+N** system

    ((1/chi) M_mass + N) m = M_mass h_ext ,   N = B^T G B ,   chi = mu_r - 1 ,

with `h_ext` the applied field H_ext L2-projected onto RT0 (HDiv order 0).  Verified against the
analytic uniform-sphere magnetization `M = chi/(1 + chi D) H` (D = demag factor ~ 1/3) to <= 1.4e-4
across mu_r 10..1000.  The **-N** system (`solve_material_minres`) is the loop-field-null mu_r-
independent operator but is NON-physical for the magnetization (wrong sign: -4497 vs +2249 on the
sphere) -- so this entry uses the +N `solve_linear_material` Jacobi-PCG (its iteration count grows
mildly with mu_r: ~50 -> ~162 for mu_r 10 -> 1000, bounded).

KELVIN-LESS: the 1/r charge Gram IS the open boundary (a volume integral method like MMM/MSC); only
the iron is meshed -- no air box / Kelvin needed.

## Scope (M1, this increment)
- LINEAR, single isotropic soft-iron region (scalar `mu_r`).
- Arbitrary applied field `H_ext` (any NGSolve CoefficientFunction: `rad.RadiaField(coil,'h')`,
  analytic, ...).
- Returns per-element M (n_el x 3), the volume-average M, iters, the demag factor.

NONLINEAR (BH table) + per-region materials + the M0 parity gate are the remaining productionization
steps (docs/hdiv_vim/PRODUCTIONIZATION.md).  Until they land, the yano-type MSC stays the `rad.Solve`
default demag backend (`radia.set_demag_backend`); this entry does not touch it.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- this library helper does NOT
open a TaskManager; the caller (calc_*.py / example / test) wraps the call in `with ng.TaskManager():`.
"""
import numpy as np
import scipy.sparse as sp
import ngsolve as ng

import radia._radia_pybind as _rp
from . import _core as _tet


def _approx_diag_N(B_csr, diagG):
    """Approximate Jacobi diagonal of N = B^T G B from the diagonal of G only:
    diag(N)_f ~ sum_a B[a,f]^2 G_aa.  O(nnz(B)); a Jacobi preconditioner tolerates the dropped
    off-diagonal cross terms (the +N solve still converges to the exact solution, just a few more
    iters).  Used on the scalable path where the dense N (hence its exact diagonal) is never formed."""
    B2 = B_csr.multiply(B_csr)                       # elementwise square, (n_charge x n_face)
    return np.asarray(B2.T @ np.asarray(diagG, float)).ravel()


def hdiv_demag_solve(mesh, mu_r, H_ext, *, scalable=False, gram_eps=1e-8, leaf=32, eta=2.0,
                     tol=1e-10, maxit=5000):
    """Linear soft-iron demag solve via the HDiv-type VIM (the +N physical material system).

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Tet mesh of the soft-iron body ONLY (KELVIN-LESS: the 1/r Gram is the open boundary).
    mu_r : float
        Relative permeability (> 1); chi = mu_r - 1, isotropic linear.
    H_ext : ngsolve.CoefficientFunction
        Applied (incident) field H_ext over the mesh (A/m).
    scalable : bool
        False (default): dense analytic charge Gram -> EXACT Jacobi diagonal (small problems,
        golden-exact).  True: skip the dense O(N^2) Gram; demag apply via the C++ charge-Gram
        H-matvec + the approximate Jacobi diagonal (the production-scale path).
    gram_eps, leaf, eta : float, int, float
        H-matrix (ACA) parameters for `_ChargeGramHMatrix`.
    tol, maxit : float, int
        Jacobi-PCG tolerance + iteration cap (RAISES if not converged -- CLAUDE.md No-Fallbacks).

    Returns
    -------
    dict
        M        : (n_el, 3) per-element magnetization (A/m).
        M_avg    : (3,) volume-average magnetization (A/m).
        iters    : Jacobi-PCG iteration count.
        demag    : demag factor D (Rayleigh quotient; ~1/3 for a sphere).
        ndof, n_el, n_charge : sizes.
    """
    chi = float(mu_r) - 1.0
    if chi <= 0.0:
        raise ValueError("hdiv_demag_solve: mu_r must be > 1 for a soft-iron demag solve (got %r)"
                         % (mu_r,))
    inv_chi = 1.0 / chi

    d = _tet.build_demag(mesh, skip_dense_gram=True) if scalable \
        else _tet.build_demag(mesh, analytic_gram=True)
    n_face = int(d["ndof"]); n_el = int(d["n_el"])
    Mm = d["M_mass"]                                  # sparse CSR (scalable) or dense (reference)
    B = d["B_csr"]

    H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                               n_el=n_el, eps=gram_eps, leaf=leaf, eta=eta)

    # source projection: H_ext L2-projected onto RT0 -> h_ext DOFs; physical RHS = M_mass h_ext
    fes = ng.HDiv(mesh, order=0)
    gfH = ng.GridFunction(fes); gfH.Set(H_ext)
    h_ext = gfH.vec.FV().NumPy().copy()
    rhs = np.asarray(Mm @ h_ext).ravel()

    # Jacobi diagonal of the +N system ((1/chi) M_mass + N)
    Mm_diag = Mm.diagonal() if sp.issparse(Mm) else np.diag(Mm)
    diagN = _approx_diag_N(B, d["self_energy"]) if scalable else np.diag(d["N"])
    prec = inv_chi * Mm_diag + diagN

    mco = sp.coo_matrix(Mm)
    res = H.solve_linear_material(
        list(map(int, B.indptr)), list(map(int, B.indices)), list(B.data.astype(float)),
        n_face, list(map(int, mco.row)), list(map(int, mco.col)), list(mco.data.astype(float)),
        inv_chi, list(prec.astype(float)), list(rhs.astype(float)), float(tol), int(maxit))
    iters = int(res["iters"])
    if iters >= maxit:
        raise RuntimeError("hdiv_demag_solve: +N Jacobi-PCG did not converge in %d iters "
                           "(mu_r=%g, n_face=%d)" % (maxit, mu_r, n_face))
    m = np.asarray(res["m"], float)

    # per-element M: L2-project the RT0 magnetization onto VectorL2(0).  VectorL2 order-0 DOFs are
    # COMPONENT-major -> reshape(3, n_el).T gives (n_el, 3).
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol

    # demag factor (Rayleigh quotient of the unit-M projection)
    mu = d["m_unit"]
    if scalable:
        Nmu = B.T @ np.asarray(H.matvec((B @ mu).tolist()), float)
        D = float((mu @ Nmu) / float(mu @ (Mm @ mu)))
    else:
        D = _tet.demag_factor(d)

    return dict(M=M_el, M_avg=M_avg, iters=iters, demag=D,
                ndof=n_face, n_el=n_el, n_charge=int(d["n_charge"]))
