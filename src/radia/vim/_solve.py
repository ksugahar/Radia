"""radia.vim._solve -- consolidated production demag-solve entry (productionization M1).

`hdiv_demag_solve(mesh, mu_r=.., H_ext=..)`           -- LINEAR soft iron (scalar mu_r), and
`hdiv_demag_solve(mesh, bh_table=.., H_ext=..)`       -- NONLINEAR soft iron (real BH table),

the candidate replacement for the yano-type MSC hex/wedge soft-iron demag in `rad.Solve`.  Both modes
take an ARBITRARY applied field `H_ext` (any NGSolve CoefficientFunction -- e.g. a coil's Biot-Savart
field `rad.RadiaField(coil,'h')`, the C-type electromagnet driver) and return per-element M.

## Formulation (verified-first, 2026-06-15)
Physical magnetization system, **+N**:  ((1/chi) M_mass + N) m = M_mass h_ext ,  N = B^T G B,
with `h_ext` = H_ext L2-projected onto RT0 (HDiv order 0).  (The -N system is mu_r-independent but
NON-physical -- wrong sign -- so +N is solved.)  LINEAR verified vs the analytic sphere to 1.4e-4;
NONLINEAR is the damped tensor-tangent Newton (the same machinery as `_nonlinear.solve_nonlinear_newton`,
generalized to an arbitrary source + a real BH table + per-element M output).

## Preconditioner -- the FULL HDiv mass inverse (mesh + mu_r robust)
The +N solves (linear system, nonlinear Newton step, scalar-chi Picard warmstart) are GMRES
preconditioned with `M_mass^{-1}` (a one-time sparse LU of the local RT0 HDiv mass), which DEFLATES the
de-Rham loops (ker B) -> the iteration count stays bounded as the mesh refines AND as mu_r grows.  A
Jacobi diagonal is NOT mesh-robust (the +N count blew past 5000 iters at 7224 faces -- measured).
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
from ._nonlinear import _bh_table_funcs, _table_tensor_tangent, _bf_to_csr

_MU0 = 4e-7 * _PI
# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect once.
_GMRES_TOL = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def hdiv_demag_solve(mesh, mu_r=None, H_ext=None, *, bh_table=None, scalable=True, gram_eps=1e-10,
                     leaf=32, eta=2.0, tol=1e-8, maxit=4000, gmres_restart=400,
                     newton_maxit=80, newton_tol=1e-7, picard_warmstart=8):
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
                                   n_el=int(d["n_el"]), eps=gram_eps, leaf=leaf, eta=eta)

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
        m, iters = _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec, mu, D, denom,
                                    n_face, fes, gram_eps, tol, gmres_restart,
                                    newton_maxit, newton_tol, picard_warmstart)

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


def _solve_nonlinear(mesh, bh_table, h_ext, Mm, N_apply, Mfac, Mprec, mu, Dscal, denom,
                     n_face, fes, gram_eps, gmres_tol, gmres_restart,
                     newton_maxit, newton_tol, picard_warmstart):
    """Damped tensor-tangent Newton on the constitutive residual F(m)=M_mass m - b_M(H),
    H = h_ext - M_mass^{-1} N m, with a real BH table + an arbitrary applied field.  Scalable: N via
    the charge-Gram H-matvec, every +N solve M_mass^{-1}-preconditioned GMRES, scalar-chi Picard
    warmstart, Armijo line search, fail-loud on non-convergence (CLAUDE.md No-Fallbacks)."""
    arr = np.asarray(bh_table, float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("hdiv_demag_solve: bh_table must be [[H,B], ...] (A/m, T)")
    Harr, Barr = arr[:, 0], arr[:, 1]
    Mof, Bpch, Bder, Hmax, Mmax = _bh_table_funcs(Harr, Barr)
    chi_init = max(float(Bder(0.0)) / _MU0 - 1.0, 1.0)

    def Dop_apply(v):
        return Mfac.solve(N_apply(v))

    Id = ng.Id(3)
    uf, vf = fes.TnT()
    gfH = ng.GridFunction(fes)
    vol = ng.Integrate(ng.CoefficientFunction(1.0), mesh)
    # representative scalar applied-field magnitude for the scalar-chi Picard warmstart
    gfHe = ng.GridFunction(fes); gfHe.vec.FV().NumPy()[:] = h_ext
    H0rep = float(ng.Integrate(ng.sqrt(ng.InnerProduct(gfHe, gfHe) + 1e-30), mesh) / vol)
    rhs0 = np.asarray(Mm @ h_ext).ravel()

    def set_field(m):
        gfH.vec.FV().NumPy()[:] = h_ext - Dop_apply(m)

    def constit():
        return _table_tensor_tangent(gfH, mesh, Bpch, Bder, Hmax, Mmax, Id)

    def bM():
        Mcf, _ = constit()
        lf = ng.LinearForm(fes); lf += Mcf * vf * ng.dx; lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def Mavg(m):
        return float((mu @ np.asarray(Mm @ m).ravel()) / denom)

    def Fnorm(m):
        set_field(m)
        return float(np.linalg.norm(np.asarray(Mm @ m).ravel() - bM()))

    # scalar-chi Picard warmstart (GMRES on the symmetric +N system) -> lands in the Newton basin
    chi = chi_init
    m = np.zeros(n_face)
    for _ in range(picard_warmstart):
        Aop = spla.LinearOperator((n_face, n_face),
                                  matvec=lambda v, c=chi: (1.0 / c) * np.asarray(Mm @ np.asarray(v, float)).ravel()
                                  + N_apply(v))
        m, _ = spla.gmres(Aop, rhs0, M=Mprec, x0=m, maxiter=20, restart=int(gmres_restart),
                          **{_GMRES_TOL: gmres_tol})
        Hi = H0rep - Dscal * Mavg(m)
        chi = 0.5 * chi + 0.5 * (Mof(Hi) / Hi if abs(Hi) > 1e-30 else chi)

    converged = False
    relF = float("inf")
    nit = 0
    for it in range(newton_maxit):
        nit = it + 1
        set_field(m)
        Mcf, tang = constit()
        lf = ng.LinearForm(fes); lf += Mcf * vf * ng.dx; lf.Assemble()
        F = np.asarray(Mm @ m).ravel() - lf.vec.FV().NumPy()
        nF = float(np.linalg.norm(F))
        relF = nF / (np.linalg.norm(np.asarray(Mm @ m).ravel()) + 1e-30)
        if relF < newton_tol:
            converged = True
            break
        T = ng.BilinearForm(fes); T += ng.InnerProduct(tang * uf, vf) * ng.dx; T.Assemble()
        Tcsr = _bf_to_csr(T)

        def Japply(v):
            v = np.asarray(v, float)
            return np.asarray(Mm @ v).ravel() + Tcsr @ Dop_apply(v)

        Jop = spla.LinearOperator((n_face, n_face), matvec=Japply)
        dm, _ = spla.gmres(Jop, -F, M=Mprec, maxiter=20, restart=int(gmres_restart),
                           **{_GMRES_TOL: gmres_tol})
        lam = 1.0
        while lam > 1e-7 and Fnorm(m + lam * dm) >= nF:
            lam *= 0.5
        m = m + lam * dm
    if not converged:
        raise RuntimeError("hdiv_demag_solve (nonlinear): Newton did NOT converge -- relF=%.2e > "
                           "newton_tol=%.1e after %d iters (returning M would be a silent wrong "
                           "result)" % (relF, newton_tol, nit))
    return m, nit
