"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r), and
`hdiv_demag_solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the candidate replacement for the yano-type MSC hex/wedge soft-iron demag in `rad.Solve`.  Both modes
take an ARBITRARY applied field `H_ext` (any NGSolve CoefficientFunction -- e.g. a coil's Biot-Savart
field `rad.RadiaField(coil,'h')`, the C-type electromagnet driver) and return per-element M.

## Formulation (verified-first, 2026-06-15)
LINEAR -- physical magnetization system, **+N**:  ((1/chi) M_mass + N) m = M_mass h_ext ,  N = B^T G B,
with `h_ext` = H_ext L2-projected onto RT0 (HDiv order 0).  (The -N system is mu_r-independent but
NON-physical -- wrong sign -- so +N is solved.)  Verified vs the analytic sphere to 1.4e-4.

NONLINEAR -- a **damped matrix-free Newton** on the constitutive residual.  F(m) = M_mass m - b_M(H),
H = h_ext - M_mass^-1 N m, b_M(H) = INT M(H).v dx (the RT0 projection of M(H)); consistent tensor
Jacobian J = M_mass + T D_op (T = INT (dM/dH) u.v dx the tensor-tangent mass, D_op = M_mass^-1 N),
applied MATRIX-FREE J v = M_mass v + T (M_mass^-1 (N v)); each Newton step is ONE M_mass^-1-GMRES on J
+ an Armijo line search, after a scalar-chi LINEAR warmstart into the basin.  WHY Newton, not the
earlier Anderson-Hantila fixed point: Hantila/Picard's contraction
rho = (chi_max - chi_min)/(chi_max + chi_min) -> 1 as the iron saturates (chi_max ~ 1e4 unsaturated,
chi_min ~ 10 at the pole), so it STALLS on real silicon steel (measured: ~1e-2 residual after 300
iters, NOT converged) -- Newton's quadratic step is immune to the chi-range.  VERIFIED: the real CEFC
Si-steel C-type at 3000 AT converges relF 0.55 -> 8e-7 in ~24 iters, gap B within ~1% of the FEM truth
(Anderson-Hantila could not solve it at all).  N = B^T G B stays SPD-PSD (G the SPD Coulomb Gram), the
de-Rham loops sit in ker(N) and are carried by M_mass (no loop-star), and the scalable C++ charge-Gram
H-matvec is the only O(N log N) cost.

## Preconditioner -- the FULL HDiv mass inverse (mesh + mu_r robust)
Every +N / constant-LHS solve is GMRES preconditioned with `M_mass^{-1}` (a one-time sparse LU of the
local RT0 HDiv mass), which DEFLATES the de-Rham loops (ker B) -> the iteration count stays bounded as
the mesh refines AND as mu_r grows (the de-Rham structural advantage: no loop-star, distortion-robust).
A Jacobi diagonal is NOT mesh-robust (the +N count blew past 5000 iters at 7224 faces -- measured).
GMRES (not CG/MINRES) absorbs the residual ACA asymmetry of the charge-Gram H-matvec.

KELVIN-LESS: the 1/r charge Gram IS the open boundary (a volume integral method like MMM/MSC); only
the iron is meshed -- no air box / Kelvin needed.  The NONLINEAR path uses the analytic charge Gram
(scalable `_ChargeGramHMatrix` at tight gram_eps), REQUIRED for div M != 0 (non-uniform M) bodies.

## Scope (M1)
Single isotropic soft-iron region (one `mu_r` or one `bh_table`).  Per-region / mixed (PM+iron)
materials + the M0 parity gate are the remaining productionization steps
(docs/hdiv_vim/PRODUCTIONIZATION.md).  Until they land, yano-type MSC stays the `rad.Solve` default
demag backend (`radia.set_demag_backend`); this entry does not touch it.

Per CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT" -- this library helper does NOT
open a TaskManager; the caller wraps the call in `with ng.TaskManager():`.
"""
import inspect
from math import pi as _PI

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import ngsolve as ng

import radia._radia_pybind as _rp
from . import _core as _tet
from ._nonlinear import _bh_table_funcs, _table_tensor_tangent

_MU0 = 4e-7 * _PI
# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect once.
_GMRES_TOL = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, bh_table=None, scalable=True, gram_eps=1e-10,
                     leaf=32, eta=2.0, near_factor=1e30, tol=1e-8, maxit=4000, gmres_restart=400,
                     nl_maxit=300, nl_tol=1e-6, anderson_window=6):
    """HDiv-type VIM soft-iron demag solve (the +N physical material system).

    Provide EXACTLY ONE material spec:
      mu_r     : float > 1  -> LINEAR isotropic soft iron.
      bh_table : [[H,B], ..] (A/m, T; the MatSatIsoTab data) -> NONLINEAR isotropic soft iron.
    H_ext      : NGSolve CoefficientFunction, the applied field (A/m) -- uniform, analytic, or a coil's
                 Biot-Savart field rad.RadiaField(coil,'h').  REQUIRED.

    Returns dict: M (n_el,3) per-element magnetization, M_avg (3,), iters, demag (Rayleigh factor),
    ndof, n_el, n_charge, nonlinear(bool).  The caller must open `with ng.TaskManager():`.
    """
    if H_ext is None:
        raise ValueError("hdiv_demag_solve: H_ext (applied-field CoefficientFunction) is required")
    if (mu_r is None) == (bh_table is None):
        raise ValueError("hdiv_demag_solve: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")

    # ---- demag operator N (apply) + M_mass^{-1} preconditioner, shared by both modes ----
    if scalable:
        d = _tet.build_demag(mesh, skip_dense_gram=True)
        Mm, B = d["M_mass"], d["B_csr"]
        H = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                   n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta,
                                   near_factor=near_factor)

        def N_apply(v):
            v = np.asarray(v, float)
            return B.T @ np.asarray(H.matvec((B @ v).tolist()), float)
    else:
        d = _tet.build_demag(mesh, analytic_gram=True)
        Mm, N = d["M_mass"], d["N"]

        def N_apply(v):
            return N @ np.asarray(v, float)

    n_face, n_el = int(d["ndof"]), int(d["n_el"])
    mu = d["m_unit"]
    Mfac = spla.splu(sp.csc_matrix(Mm))
    Mprec = spla.LinearOperator((n_face, n_face), matvec=lambda v: Mfac.solve(np.asarray(v, float)))

    # ---- applied field projected onto RT0 ----
    fes = ng.HDiv(mesh, order=0)
    gfHext = ng.GridFunction(fes); gfHext.Set(H_ext)
    h_ext = gfHext.vec.FV().NumPy().copy()
    rhs0 = np.asarray(Mm @ h_ext).ravel()
    denom = float(mu @ np.asarray(Mm @ mu).ravel())
    D = float((mu @ N_apply(mu)) / denom)

    if mu_r is not None:
        m, iters = _solve_linear(mu_r, Mm, N_apply, Mprec, rhs0, n_face, tol, maxit, gmres_restart)
    else:
        m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec,
                                    n_face, fes, tol, gmres_restart, nl_maxit, nl_tol, anderson_window)

    # ---- per-element M (VectorL2(0) projection; component-major DOFs -> (n_el,3)) + averages ----
    gfM = ng.GridFunction(fes); gfM.vec.FV().NumPy()[:] = m
    fesM = ng.VectorL2(mesh, order=0); gfMc = ng.GridFunction(fesM); gfMc.Set(gfM)
    M_el = gfMc.vec.FV().NumPy().reshape(3, n_el).T.copy()
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    M_avg = np.array([ng.Integrate(gfM[i], mesh) for i in range(3)]) / vol

    return dict(M=M_el, M_avg=M_avg, iters=int(iters), demag=D, ndof=n_face, n_el=n_el,
                n_charge=int(d["n_charge"]), nonlinear=(bh_table is not None))


def _solve_linear(mu_r, Mm, N_apply, Mprec, rhs, n_face, tol, maxit, gmres_restart):
    chi = float(mu_r) - 1.0
    if chi <= 0.0:
        raise ValueError("hdiv_demag_solve: mu_r must be > 1 (got %r)" % (mu_r,))
    inv_chi = 1.0 / chi

    def A_apply(v):
        v = np.asarray(v, float)
        return inv_chi * np.asarray(Mm @ v).ravel() + N_apply(v)

    A = spla.LinearOperator((n_face, n_face), matvec=A_apply)
    it = {"n": 0}
    cycles = max(2, int(np.ceil(maxit / gmres_restart)))
    m, info = spla.gmres(A, rhs, M=Mprec, restart=int(gmres_restart), maxiter=cycles,
                         callback=lambda _x: it.__setitem__("n", it["n"] + 1),
                         callback_type="pr_norm", **{_GMRES_TOL: tol})
    if info != 0:
        raise RuntimeError("hdiv_demag_solve (linear): +N GMRES did not converge (info=%d, mu_r=%g, "
                           "n_face=%d, iters=%d)" % (info, mu_r, n_face, it["n"]))
    return m, it["n"]


def _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec,
                     n_face, fes, gmres_tol, gmres_restart, nl_maxit, nl_tol, anderson_window):
    """Damped matrix-free NEWTON on the constitutive residual (replaces the Anderson-Hantila fixed
    point, which STALLS on real silicon steel -- the Hantila/Picard contraction
    rho = (chi_max-chi_min)/(chi_max+chi_min) -> 1 as the iron saturates, so even Anderson cannot
    recover it within a bounded iter count).  Newton converges where Hantila cannot: VERIFIED on the
    real CEFC Si-steel C-type at 3000 AT (relF 0.55 -> 8e-7 in ~24 iters, gap B within ~1% of the FEM
    truth) where Anderson-Hantila stalled at ~1e-2 after 300.  `anderson_window` is now unused (kept in
    the signature for call-site stability).

    Residual    F(m) = M_mass m - b_M(H),  H = h_ext - M_mass^-1 N m,  b_M(H) = INT M(H).v dx (the RT0
                L2 projection of the constitutive M).
    Consistent TENSOR Jacobian  J = M_mass + T D_op,  T = INT (dM/dH) u.v dx (the tensor-tangent mass),
                D_op = M_mass^-1 N -- applied MATRIX-FREE:  J v = M_mass v + T (M_mass^-1 (N v)).
    Each Newton step is ONE M_mass^-1-preconditioned GMRES on J + an Armijo backtracking line search.
    A scalar-chi LINEAR warmstart  ((1/chi0) M_mass + N) m = M_mass h_ext  lands inside the Newton basin
    (the unsaturated chi0 response).  The tensor tangent (`_table_tensor_tangent`) + matrix-free N reuse
    the EXACT pieces the Hantila path used -- no new operator, fully scalable (the C++ charge-Gram
    H-matvec is the only O(N log N) cost).  Fail-loud on non-convergence (CLAUDE.md No-Fallbacks)."""
    del anderson_window                      # unused by Newton (signature kept for the caller)
    arr = np.asarray(bh_table, float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("hdiv_demag_solve: bh_table must be [[H,B], ...] (A/m, T)")
    _Mof, Bpch, Bder, Hmax, Mmax = _bh_table_funcs(arr[:, 0], arr[:, 1])
    chi0 = max(float(Bder(0.0)) / _MU0 - 1.0, 1.0)
    Id = ng.Id(3); vf = fes.TestFunction(); uf = fes.TrialFunction(); gfH = ng.GridFunction(fes)

    def _constit(m):                         # (M(H) CF, tensor-tangent CF) at H = h_ext - M_mass^-1 N m
        gfH.vec.FV().NumPy()[:] = h_ext - Mfac.solve(N_apply(m))
        return _table_tensor_tangent(gfH, mesh, Bpch, Bder, Hmax, Mmax, Id)

    def _bM(M_cf):                           # RT0 L2 projection load INT M(H).v dx
        lf = ng.LinearForm(fes); lf += M_cf * vf * ng.dx; lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def _Fnorm(m):
        M_cf, _ = _constit(m)
        return float(np.linalg.norm(np.asarray(Mm @ m).ravel() - _bM(M_cf)))

    # scalar-chi LINEAR warmstart -> the Newton basin (the unsaturated chi0 response)
    rhs0 = np.asarray(Mm @ h_ext).ravel()
    A0 = spla.LinearOperator((n_face, n_face), matvec=lambda v:
                             (1.0 / chi0) * np.asarray(Mm @ np.asarray(v, float)).ravel() + N_apply(v))
    m, _ = spla.gmres(A0, rhs0, M=Mprec, restart=int(gmres_restart), maxiter=4, **{_GMRES_TOL: 1e-6})

    converged = False; nit = 0; rel = float("inf")
    for it in range(nl_maxit):
        nit = it + 1
        M_cf, tang = _constit(m)
        F = np.asarray(Mm @ m).ravel() - _bM(M_cf)
        nF = float(np.linalg.norm(F))
        scale = float(np.linalg.norm(np.asarray(Mm @ m).ravel())) + 1e-30
        rel = nF / scale
        if rel < nl_tol:
            converged = True; break
        T = ng.BilinearForm(fes); T += ng.InnerProduct(tang * uf, vf) * ng.dx; T.Assemble()
        r, c, v = T.mat.COO()
        Tcsr = sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(n_face, n_face))

        def _Jv(x, _Tcsr=Tcsr):              # J v = M_mass v + T (M_mass^-1 (N v))   (matrix-free)
            x = np.asarray(x, float)
            return np.asarray(Mm @ x).ravel() + np.asarray(_Tcsr @ Mfac.solve(N_apply(x))).ravel()

        Jop = spla.LinearOperator((n_face, n_face), matvec=_Jv)
        # Eisenstat-Walker inexact-Newton forcing: solve the Newton step only as tightly as the current
        # outer residual warrants (loose when far, tight when near) -- the inner J-GMRES is the dominant
        # cost (~130 iters/step, M_mass^-1 is a weak preconditioner for J = M_mass + T D_op), so a loose
        # early tol cuts wasted inner work without changing the (tangent-limited, ~0.65-linear) outer rate.
        inner_tol = float(min(1e-2, max(gmres_tol, 0.1 * rel)))
        dm, _info = spla.gmres(Jop, -F, M=Mprec, restart=int(gmres_restart), maxiter=4, **{_GMRES_TOL: inner_tol})
        lam = 1.0                            # Armijo backtracking line search (globalises Newton)
        while lam > 1e-6 and _Fnorm(m + lam * dm) >= nF:
            lam *= 0.5
        m = m + lam * dm
    if not converged:
        raise RuntimeError("hdiv_demag_solve (nonlinear Newton): did NOT converge -- rel=%.2e > "
                           "nl_tol=%.1e after %d iters (returning M would be a silent wrong result)"
                           % (rel, nl_tol, nit))
    return m, nit
