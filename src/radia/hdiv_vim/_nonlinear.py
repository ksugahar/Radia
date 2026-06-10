"""hdiv_demag_tet_nonlinear.py -- NONLINEAR HDiv-type VIM demag (applied-field solve + BH-curve Picard).

The linear HDiv-VIM (hdiv_demag_tet.py) computes demag FACTORS via eig(N, M_mass).  A NONLINEAR solve
needs an APPLIED-FIELD solve (M for a given H_ext), then a per-element constitutive iteration.

## Applied-field formulation (verify-first result, 2026-06-07)
The eigenvalue framing A_eig = (1/chi) M_mass - N is NOT the applied-field system.  The physical
applied-field weak form is

    A+ m = M_mass h_ext ,   A+ = (1/chi) M_mass + N         (PLUS N)

  because  M = chi (H_ext + H_demag),  H_demag,weak = -N m  =>  (1/chi) M_mass m + N m = M_mass h_ext.
For a sphere this reproduces the analytic  M/H = chi/(1 + chi D)  (D = demag factor) -- VERIFIED here
to <=2.5% for mu_r<=100 (A+); the minus-sign system gives nonsense (negative / divergent).

## Nonlinear (Picard, secant susceptibility)
With a BH curve M(H), iterate:  solve A+(chi^k) m -> M_avg, internal field H_int = H0 - D M_avg,
chi^{k+1} = M(H_int)/H_int (secant), re-solve.  Converges + saturates (M -> M_sat).  This prototype
uses the SCALAR (uniform-M sphere) update; the per-element (non-uniform body) generalization needs the
per-cell field reconstruction (RT0 -> cell H) and an MMM/MSC cross-check -- the next increment.

## Honest accuracy
The nonlinear M value carries the SAME operator accuracy as the linear demag: the centroid-monopole N
UNDER-estimates the sphere demag (D~0.31 vs analytic 1/3), so the high-chi (deep-saturation-onset)
response is off by ~the same amount.  The #3 near-field correction (build_near_correction) raises D
toward 1/3 and tightens it -- demonstrated below.
"""
import os
from math import pi

import numpy as np
import scipy.sparse as sp

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

from radia.hdiv_vim import _core as tet   # M1: core promoted to radia.hdiv_vim (was examples/feec_vim)


def _bf_to_dense(bf):
    """NGSolve BilinearForm -> dense numpy (prototype-scale matrices only)."""
    return _bf_to_csr(bf).toarray()


def _bf_to_csr(bf):
    """NGSolve BilinearForm -> scipy CSR (sparse; for the scalable matrix-free Jacobian apply)."""
    m = bf.mat
    r, c, val = m.COO()
    return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(m.height, m.width))


# scipy renamed the Krylov tolerance kwarg 'tol' -> 'rtol' (scipy >= 1.12); detect it once.
import inspect as _inspect  # noqa: E402
import scipy.sparse.linalg as _spla_probe  # noqa: E402
_GMRES_TOL = "rtol" if "rtol" in _inspect.signature(_spla_probe.gmres).parameters else "tol"
_MINRES_TOL = "rtol" if "rtol" in _inspect.signature(_spla_probe.minres).parameters else "tol"


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def _bh_curve(chi0, Msat):
    """Saturating M(H) = chi0 H / (1 + chi0|H|/Msat): slope chi0 at H=0, asymptote +-Msat."""
    return lambda H: chi0 * H / (1.0 + chi0 * abs(H) / Msat)


def _scalar_fixed_point(Mof, D, H0):
    """Correct analytic uniform-sphere root: solve M = M(H0 - D M) by bisection on f(M)=M-M(H0-DM)."""
    lo, hi = -1.0, 1.0
    f = lambda M: M - Mof(H0 - D * M)
    # widen until sign change (cap 1e12 so a REAL-table Msat ~ 1.5e6 is reachable, not just the
    # analytic curve's Msat=1)
    while f(lo) * f(hi) > 0 and hi < 1e12:
        lo *= 2; hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def solve_nonlinear(mesh, Mof, H0, near_correction=False, nsub=4, maxit=300, tol=1e-9):
    """Nonlinear applied-field HDiv-VIM solve (scalar-chi Picard on the A+ system).  Returns
    (M_avg, n_iter, D_used)."""
    d = tet.build_demag(mesh, nsub)
    N = d["N"].copy()
    M = d["M_mass"]
    mu = d["m_unit"]
    denom = mu @ M @ mu
    if near_correction:
        # demag operator with the #3 near-field correction folded in: N_eff = B^T (G + corr) B
        corr = tet.build_near_correction(mesh, d, nsub=nsub, near_factor=2.0)
        B = d["B_csr"]
        N = N + np.asarray((B.T @ corr @ B).todense())
    D = float((mu @ N @ mu) / denom)
    b0 = M @ mu
    chi, Mavg = 1000.0, 0.0
    for it in range(maxit):
        A = (1.0 / chi) * M + N
        m = np.linalg.solve(A, H0 * b0)
        Mavg = float((mu @ M @ m) / denom)
        Hint = H0 - D * Mavg
        chi_new = Mof(Hint) / Hint if abs(Hint) > 1e-30 else 1000.0
        if abs(chi_new - chi) < tol * chi:
            break
        chi = 0.5 * chi + 0.5 * chi_new
    return Mavg, it + 1, D


def _tensor_tangent_cfs(gfH, chi0, Msat, Id):
    """Constitutive M(H) field + the CONSISTENT TENSOR tangent dM/dH for the saturating curve.

    For the isotropic vector law M(H) = chi_sec(|H|) H, the consistent tangent is the rank-1 + scalar
    split (slope along H, secant perpendicular):

        dM/dH = chi_diff Hhat(x)Hhat + chi_sec (I - Hhat(x)Hhat),
        chi_sec  = chi0 / (1 + chi0|H|/Msat)        (M = chi_sec H),
        chi_diff = chi0 / (1 + chi0|H|/Msat)^2       (d|M|/d|H|).

    The naive scalar tangent chi_diff*I FAILS at moderate drive (where chi_sec != chi_diff): e.g. at
    |H|~4e-4 here chi_sec=725 but chi_diff=510, so omitting the perpendicular secant term gives a
    badly wrong Jacobian and Newton crawls / stalls.
    """
    H2 = ng.InnerProduct(gfH, gfH) + 1e-30
    Hm = ng.sqrt(H2)
    sec = chi0 / (1.0 + chi0 * Hm / Msat)
    dif = chi0 / (1.0 + chi0 * Hm / Msat) ** 2
    par = ng.OuterProduct(gfH, gfH) / H2                 # Hhat (x) Hhat
    return sec * gfH, dif * par + sec * (Id - par)


_MU0 = 4e-7 * pi


def _bh_table_funcs(Harr, Barr):
    """Build M(H) + the PCHIP B(H)/B'(H) from a REAL [[H,B]] table (the same data Radia's MatSatIsoTab
    consumes): monotone C1 interpolation -> smooth dM/dH for Newton.  M(H) = B(|H|)/mu0 - |H|.

    BEYOND the table (|H| > H_max) the curve SATURATES: B extends with slope mu0 (mu_r -> 1), so
    M -> M(H_max) = const and dM/dH -> 0 -- the physical saturation, matching Radia's MatSatIsoTab
    linear-B extension.  (Without this, PCHIP polynomial extrapolation blows up: M >> M_sat.)
    Returns (Mof, Bpch, Bder, Hmax, Mmax)."""
    from scipy.interpolate import PchipInterpolator
    Harr = np.asarray(Harr, float)
    Barr = np.asarray(Barr, float)
    Bpch = PchipInterpolator(Harr, Barr)
    Bder = Bpch.derivative()
    Hmax = float(Harr[-1])
    Mmax = float(Bpch(Hmax) / _MU0 - Hmax)

    def Mof(H):
        a = abs(H)
        m = (float(Bpch(a)) / _MU0 - a) if a <= Hmax else Mmax
        return float(np.sign(H) * m)

    return Mof, Bpch, Bder, Hmax, Mmax


def _table_tensor_tangent(gfH, mesh, Bpch, Bder, Hmax, Mmax, Id):
    """Table version of _tensor_tangent_cfs: per-element chi_sec / chi_diff from the PCHIP BH table,
    set as L2(0) element-constant fields (chi_sec = M(|H|)/|H|, chi_diff = dM/d|H| = B'(|H|)/mu0 - 1).
    Beyond the table (|H| > Hmax) M saturates to Mmax (chi_diff = 0)."""
    l2 = ng.L2(mesh, order=0)
    gfm = ng.GridFunction(l2)
    gfm.Set(ng.sqrt(ng.InnerProduct(gfH, gfH) + 1e-30))
    Hmag = np.maximum(gfm.vec.FV().NumPy(), 1e-30)         # per-element |H|
    sat = Hmag > Hmax
    a_cl = np.minimum(Hmag, Hmax)
    m_e = np.where(sat, Mmax, Bpch(a_cl) / _MU0 - a_cl)    # M(|H|), saturated beyond table
    sec_e = m_e / Hmag                                     # chi_sec = M/|H| (decreases beyond sat)
    dif_e = np.where(sat, 0.0, Bder(a_cl) / _MU0 - 1.0)    # chi_diff = dM/d|H| (0 beyond sat)
    gf_sec = ng.GridFunction(l2)
    gf_sec.vec.FV().NumPy()[:] = sec_e
    gf_dif = ng.GridFunction(l2)
    gf_dif.vec.FV().NumPy()[:] = dif_e
    H2 = ng.InnerProduct(gfH, gfH) + 1e-30
    par = ng.OuterProduct(gfH, gfH) / H2
    return gf_sec * gfH, gf_dif * par + gf_sec * (Id - par)


def solve_nonlinear_newton(mesh, chi0, Msat, H0, near_correction=True, nsub=4,
                           picard_warmstart=8, maxit=200, tol=1e-10, wilton_surface=False,
                           bh_table=None, analytic_gram=False, require_convergence=True):
    """Robust NONLINEAR HDiv-VIM solve via damped (line-search) Newton-Raphson.

    Newton on the constitutive residual in RT0 coefficient form

        F(m) = M_mass m - b_M(H),   H = h_ext - D_op m,   D_op = M_mass^{-1} N,
        b_M(H) = INT M(H).v dx      (L2 projection of the constitutive M onto RT0),

    with the CONSISTENT TENSOR Jacobian (derived, not the scalar approximation)

        J = M_mass + T D_op,   T = INT (dM/dH) u.v dx     (tensor tangent mass).

    Three ingredients make it robust AND accurate (verify-first, 2026-06-07):
      (a) the tensor tangent dM/dH (see _tensor_tangent_cfs) -- the scalar chi_diff*I stalls at
          moderate drive;
      (b) the #3 NEAR-FIELD CORRECTION on N -- the centroid-monopole N has poor PER-ELEMENT (local)
          field accuracy, so per-element Newton converges to a WRONG root (e.g. M=0.09 instead of 0.30
          on the sphere at H0=0.1); the near correction restores the local fields -> correct root;
      (c) line-search DAMPING (Armijo backtracking) for global robustness + a scalar-chi PICARD
          warmstart to enter the basin at the BH knee (chi_sec/chi_diff ~ 8x), where pure Newton
          crawls.

    This is the Newton-Raphson counterpart of Radia's MMM/MSC newton_damping=True path; the outer
    Newton machinery (tangent susceptibility + damping) is SHARED -- only the demag operator N differs
    (here N = B^T G B, the HDiv-VIM charge form).

    Returns (M_avg, n_newton_iter, D_used).  The caller must open `with ng.TaskManager():`.
    """
    # constitutive: analytic saturating curve (chi0, Msat) OR a REAL tabulated BH curve (bh_table =
    # (Harr, Barr), the same data Radia's MatSatIsoTab consumes).
    if bh_table is not None:
        Mof, _Bpch, _Bder, _Hmax, _Mmax = _bh_table_funcs(bh_table[0], bh_table[1])
        chi_init = max(float(_Bder(0.0)) / _MU0 - 1.0, 1.0)   # table's zero-field susceptibility
    else:
        Mof = _bh_curve(chi0, Msat)
        chi_init = chi0
    d = tet.build_demag(mesh, nsub, wilton_surface=wilton_surface, analytic_gram=analytic_gram)
    M_mass = d["M_mass"]
    mu = d["m_unit"]
    denom = mu @ M_mass @ mu
    N = d["N"].copy()
    if near_correction and not analytic_gram:
        # near-field correction.  With wilton_surface the surface-surface block is already exact, so
        # correct only the VOLUME-involving (cell-cell, cell-face) near pairs (skip_surface_surface)
        # -- the per-element nonlinear Newton needs the volume near-field, and double-correcting the
        # surface would over-count.  Without wilton_surface, correct all near pairs (the old path).
        # analytic_gram makes the WHOLE Gram exact -> no near-correction needed (skip it).
        corr = tet.build_near_correction(mesh, d, nsub=nsub, near_factor=2.0,
                                         skip_surface_surface=wilton_surface)
        B = d["B_csr"]
        N = N + np.asarray((B.T @ corr @ B).todense())
    Dop = np.linalg.solve(M_mass, N)
    Dscal = float((mu @ N @ mu) / denom)
    b0 = M_mass @ mu

    def Mavg(m):
        return float((mu @ M_mass @ m) / denom)

    fes = ng.HDiv(mesh, order=0)
    u, v = fes.TnT()
    gfH = ng.GridFunction(fes)
    Id = ng.Id(3)

    def constit():
        """(M(H) field, consistent tensor tangent) at the current gfH -- analytic or table."""
        if bh_table is not None:
            return _table_tensor_tangent(gfH, mesh, _Bpch, _Bder, _Hmax, _Mmax, Id)
        return _tensor_tangent_cfs(gfH, chi0, Msat, Id)

    def set_field(m):
        gfH.vec.FV().NumPy()[:] = H0 * mu - Dop @ m

    def bM():
        Mcf, _ = constit()
        lf = ng.LinearForm(fes)
        lf += Mcf * v * ng.dx
        lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def Fnorm(m):
        set_field(m)
        return np.linalg.norm(M_mass @ m - bM())

    # (c) scalar-chi Picard warmstart: cheap, robust, lands inside the Newton basin (matters at the knee).
    chi = chi_init
    m = np.zeros(M_mass.shape[0])
    for _ in range(picard_warmstart):
        m = np.linalg.solve((1.0 / chi) * M_mass + N, H0 * b0)
        Hi = H0 - Dscal * Mavg(m)
        chi = 0.5 * chi + 0.5 * (Mof(Hi) / Hi if abs(Hi) > 1e-30 else chi)

    # damped operator Newton on the constitutive residual
    nit = 0
    M_prev = Mavg(m)
    converged = False
    for it in range(maxit):
        nit = it + 1
        set_field(m)
        Mcf, tang = constit()
        lf = ng.LinearForm(fes)
        lf += Mcf * v * ng.dx
        lf.Assemble()
        F = M_mass @ m - lf.vec.FV().NumPy()
        nF = np.linalg.norm(F)
        if nF < tol * (np.linalg.norm(M_mass @ m) + 1e-30):
            converged = True; break
        T = ng.BilinearForm(fes)
        T += ng.InnerProduct(tang * u, v) * ng.dx
        T.Assemble()
        dm = np.linalg.solve(M_mass + _bf_to_dense(T) @ Dop, -F)
        lam = 1.0                                   # (c) Armijo backtracking line search
        while lam > 1e-7 and Fnorm(m + lam * dm) >= nF:
            lam *= 0.5
        m = m + lam * dm
        # Physically meaningful convergence: the quantity of interest (M_avg) has stopped moving.
        # Near the knee the line search takes tiny steps, so the tight residual tol is never met even
        # though M_avg is already correct -- key on the observable, not the raw residual norm.
        M_now = Mavg(m)
        if abs(M_now - M_prev) < 1e-8 * (abs(M_now) + 1e-30):
            converged = True; break
        M_prev = M_now
    if require_convergence and not converged:
        # fail loud (no silent non-converged result): the most common cause is the WRONG Gram for a
        # NON-uniform-M nonlinear problem -- the surface-only wilton_surface leaves the volume (cell)
        # blocks crude so Newton stalls; div M != 0 bodies (cube, C-yoke, ...) need analytic_gram=True.
        raise RuntimeError(
            f"solve_nonlinear_newton did NOT converge in maxit={maxit} iters (last M_avg={Mavg(m):.1f}). "
            f"For NON-uniform-M nonlinear problems (div M != 0: cube, C-yoke, any non-ellipsoid) pass "
            f"analytic_gram=True -- the surface-only wilton_surface Gram leaves the volume blocks crude "
            f"and Newton cannot converge. Pass require_convergence=False to inspect the non-converged iterate.")
    return Mavg(m), nit, Dscal


def solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0, nsub=4, gram_eps=1e-10,
                                    picard_warmstart=8, maxit=200, gmres_tol=1e-8, newton_tol=1e-6,
                                    near_factor=1e30, gmres_restart=400, return_timing=False, verbose=False):
    """SCALABLE damped Newton (production #2): the dense O(N^3)/O(N^2) demag is replaced by the C++
    HACApK charge-Gram H-matrix (O(N log N) apply), and the Newton system is solved ITERATIVELY
    (GMRES) -- no dense factorization anywhere.

    M2b: the C++ Gram is the ANALYTIC charge Gram (PhiTet/TriPotential, exact near AND far), so the
    demag apply is simply  N v = B^T ( H.matvec(B v) )  with H = the analytic _ChargeGramHMatrix --
    no separate sparse near-correction.  This matches the dense solve_nonlinear_newton(analytic_gram=
    True) operator entry-by-entry, so the scalable solver reproduces the dense ANALYTIC Newton (the
    correct reference for NON-uniform M: cube, C-yoke, any div M != 0 body), closing the under-
    resolution the old monopole-far + sparse-near split left.  The Jacobian
    J v = M_mass v + T M_mass^{-1} N v is applied matrix-free (M_mass factored once, sparse); GMRES
    (M_mass-preconditioned) solves each Newton step; Armijo line search + Picard warmstart as in the
    dense solver.  Returns (M_avg, n_newton_iter, D_used)."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    import radia._radia_pybind as _rp
    import time

    _t0 = time.perf_counter()
    Mof = _bh_curve(chi0, Msat)
    # SCALABLE build: skip the dense O(N^2) G/N and the loop SVD; M_mass + B come back SPARSE.  The demag
    # apply is the analytic charge-Gram H-matvec (Hg) below, so the dense N is never needed here.
    d = tet.build_demag(mesh, nsub, skip_dense_gram=True)
    M_mass = d["M_mass"]                         # sparse CSR
    mu = d["m_unit"]
    denom = float(mu @ (M_mass @ mu))            # sparse-safe Rayleigh denominator
    B = d["B_csr"]
    # M2b: the C++ ANALYTIC charge Gram H-matrix is exact near AND far (matches dense analytic_gram
    # entry-by-entry), so the demag apply is just B^T (G_analytic-Hmatvec (B v)) -- NO separate Python
    # near-correction.  This closes the non-uniform-M (div M != 0: cube, C-yoke) scalable gap that the
    # old monopole-far + sparse-near split left, while staying O(N log N).
    # near_factor=1e30 (default) = all-analytic (the tight scalable golden); a finite near_factor (e.g. 2.0)
    # enables the near/far split (analytic near, monopole far) -> the BUILD speedup, accurate to ~3% Gram.
    Hg = _rp._ChargeGramHMatrix(cell_verts=list(d["cell_verts"]), face_verts=list(d["face_verts"]),
                                n_el=int(d["n_el"]), eps=gram_eps, leaf=32, eta=2.0, near_factor=near_factor)
    Mcsc = sp.csc_matrix(M_mass)
    Mfac = spla.splu(Mcsc)                       # factor the HDiv mass ONCE (sparse, well-conditioned)

    def N_apply(v):                              # scalable demag: the analytic charge-Gram H-matvec
        q = B @ v
        return B.T @ np.asarray(Hg.matvec(q.tolist()), float)

    def Dop_apply(v):                            # M_mass^{-1} N v  (the weak demag field)
        return Mfac.solve(N_apply(v))

    ndof = d["ndof"]
    Dscal = float((mu @ N_apply(mu)) / denom)    # demag factor via the H-matrix apply
    b0 = M_mass @ mu

    def Mavg(m):
        return float((mu @ (M_mass @ m)) / denom)

    fes = ng.HDiv(mesh, order=0)
    uf, vf = fes.TnT()
    gfH = ng.GridFunction(fes)
    Id = ng.Id(3)

    def set_field(m):
        gfH.vec.FV().NumPy()[:] = H0 * mu - Dop_apply(m)

    def constit(m):
        set_field(m)
        return _tensor_tangent_cfs(gfH, chi0, Msat, Id)

    def bM(Mcf):
        lf = ng.LinearForm(fes)
        lf += Mcf * vf * ng.dx
        lf.Assemble()
        return lf.vec.FV().NumPy().copy()

    def Fnorm(m):
        Mcf, _ = constit(m)
        return np.linalg.norm(M_mass @ m - bM(Mcf))

    Mprec = spla.LinearOperator((ndof, ndof), matvec=lambda v: Mfac.solve(np.asarray(v, float)))

    _t_build = time.perf_counter()        # end of one-time setup (build_demag + Hg + splu); solve phase follows
    # scalar-chi Picard warmstart via the scalable apply (CG on the symmetric A+ = (1/chi)M_mass + N)
    chi = chi0
    m = np.zeros(ndof)
    Aplus_diag = np.maximum(np.abs((1.0 / chi) * M_mass.diagonal()) + 1e-30, 1e-30)
    for _pw in range(picard_warmstart):
        Aop = spla.LinearOperator((ndof, ndof), matvec=lambda v, c=chi: (1.0 / c) * (M_mass @ np.asarray(v, float)) + N_apply(np.asarray(v, float)))
        # GMRES, NOT MINRES/CG.  N is applied via the ACA H-matrix, which is symmetric only to ~ACA tol
        # -- MEASURED ~2-5% asymmetric at scale (each off-diagonal block + its transpose get independent
        # ACA pivots; the asymmetry grows with N).  Symmetric Krylov (MINRES/CG) STALL or DIVERGE on it
        # (the warmstart used to hit its 2000-iter cap, info=2000); GMRES is asymmetry-robust and
        # converges in ~100-150 iters (the Newton step below already uses GMRES for the same reason).
        # Warm-restart from the previous Picard iterate (x0=m).
        # restart MUST exceed the inner-iteration count or GMRES STAGNATES (restarting throws away the
        # Krylov subspace).  The +N inner solve grows mildly with N (~31 @ ndof 8573 -> ~115 @ 38383, the
        # star-space demag spectrum -- NOT the loops, which M_mass^{-1} already deflates).  restart=50 was
        # the entire "44k scaling wall": at ndof 38383 it needs 115 iters and restart=50 never converged
        # (info=20); restart>=200 converges in 115.  Default 400 covers ~100k ndof; very large N benefits
        # from a star-space preconditioner (to BOUND the iters) rather than an ever-larger restart.
        m, _info = spla.gmres(Aop, H0 * b0, M=Mprec, x0=m, maxiter=20, restart=gmres_restart, **{_GMRES_TOL: gmres_tol})
        Hi = H0 - Dscal * Mavg(m)
        chi = 0.5 * chi + 0.5 * (Mof(Hi) / Hi if abs(Hi) > 1e-30 else chi)
        if verbose:
            print(f"    [warmstart {_pw}] gmres_info={_info} chi={chi:.3e} Mavg={Mavg(m):.1f}", flush=True)

    nit = 0
    converged = False
    relF = float("inf")
    for it in range(maxit):
        nit = it + 1
        Mcf, tang = constit(m)
        F = M_mass @ m - bM(Mcf)
        nF = np.linalg.norm(F)
        # SOUND convergence: the actual nonlinear residual relF = ||F|| / ||M_mass m|| is small.  Do NOT
        # break on Mavg stagnation (|M_now - M_prev| small) -- during the slow globalization phase Mavg
        # plateaus while relF is still O(0.1), so an Mavg-stagnation break silently returns an
        # under-converged (wrong) M.  That was the "drift to 509k" bug; the operator/tangent are correct.
        relF = nF / (np.linalg.norm(M_mass @ m) + 1e-30)
        if relF < newton_tol:
            converged = True
            if verbose:
                print(f"    [newton {it:2d}] relF={relF:.2e} CONVERGED", flush=True)
            break
        T = ng.BilinearForm(fes)
        T += ng.InnerProduct(tang * uf, vf) * ng.dx
        T.Assemble()
        Tcsr = _bf_to_csr(T)

        def Japply(v):                           # J v = M_mass v + T M_mass^{-1} N v (matrix-free)
            v = np.asarray(v, float)
            return M_mass @ v + Tcsr @ Dop_apply(v)

        Jop = spla.LinearOperator((ndof, ndof), matvec=Japply)
        dm, ginfo = spla.gmres(Jop, -F, M=Mprec, maxiter=20, restart=gmres_restart, **{_GMRES_TOL: gmres_tol})
        lam = 1.0
        while lam > 1e-7 and Fnorm(m + lam * dm) >= nF:
            lam *= 0.5
        m = m + lam * dm
        if verbose:
            print(f"    [newton {it:2d}] relF={relF:.2e} gmres_info={ginfo} lam={lam:.3e} Mavg={Mavg(m):.1f}",
                  flush=True)
    if not converged:
        # FAIL LOUD (CLAUDE.md "No Fallbacks"): never silently return an under-converged M.  Slow
        # convergence here is the ill-conditioned +N warmstart (loop near-null modes); raise maxit or
        # improve the warmstart / preconditioner (the -N mu_r-independent material formulation).
        raise RuntimeError(
            f"solve_nonlinear_newton_scalable did NOT converge: relF={relF:.2e} > newton_tol={newton_tol:.1e} "
            f"after {nit} Newton iters. Returning M now would be a silent wrong result.")
    if return_timing:
        _now = time.perf_counter()
        return Mavg(m), nit, Dscal, {"t_build_s": _t_build - _t0, "t_solve_s": _now - _t_build,
                                     "t_total_s": _now - _t0, "ndof": int(d["ndof"]), "n_charge": int(d["n_charge"])}
    return Mavg(m), nit, Dscal


def main():
    mesh = _sphere(0.35)
    chi0, Msat = 1000.0, 1.0
    Mof = _bh_curve(chi0, Msat)
    # self-test runs serial: per the caller-wraps policy this HELPER module contains NO `with TaskManager()`.
    print(f"Nonlinear HDiv-VIM sphere demag (chi0={chi0}, Msat={Msat})")
    print(f"{'H0':>8} {'Picard M':>10} {'+nearcorr':>10} {'analytic(1/3)':>14} {'M/Msat':>8}")
    for H0 in (1e-4, 1e-3, 1e-2, 1e-1):
        Mmono, it1, Dm = solve_nonlinear(mesh, Mof, H0, near_correction=False)
        Mcorr, it2, Dc = solve_nonlinear(mesh, Mof, H0, near_correction=True)
        Mana = _scalar_fixed_point(Mof, 1.0 / 3.0, H0)   # analytic uses the EXACT sphere D=1/3
        print(f"{H0:8.0e} {Mmono:10.5f} {Mcorr:10.5f} {Mana:14.5f} {Mcorr/Msat:8.3f}")
    print(f"  (monopole D={Dm:.4f}, near-corrected D={Dc:.4f}, analytic D=0.3333)")
    print("  => Picard converges + saturates; near-correction moves D toward 1/3 -> closer to analytic.\n")

    # ---- damped Newton-Raphson (robust at ALL drive incl. deep saturation) ----
    print("Damped Newton-Raphson (tensor tangent + near corr + line search + Picard warmstart):")
    print(f"{'H0':>8} {'Newton M':>10} {'newton it':>10} {'analytic(1/3)':>14} {'rel':>9}")
    for H0 in (1e-2, 1e-1, 3e-1, 1.0, 5.0):
        Mn, nit, Dn = solve_nonlinear_newton(mesh, chi0, Msat, H0)
        Mana = _scalar_fixed_point(Mof, 1.0 / 3.0, H0)
        print(f"{H0:8.0e} {Mn:10.5f} {nit:10d} {Mana:14.5f} {abs(Mn - Mana) / Mana:9.1e}")
    print("  => Newton matches the analytic at every drive (where the scalar-Picard saturates but"
          " the simple per-element Picard/Hantila of earlier sessions failed).")


if __name__ == "__main__":
    main()
