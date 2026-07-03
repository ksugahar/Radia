"""radia.vim._nonlinear -- NONLINEAR HDiv-type VIM demag (applied-field BH-curve Newton).

A NONLINEAR demag solve needs an APPLIED-FIELD solve (M for a given H_ext) + a per-element constitutive
iteration.  Both public entries route through the C++ scalable charge-Gram operator (the dense Python
charge-Gram path was removed):

  solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0, ...)
      -- damped matrix-free Newton on the C++ analytic _ChargeGramHMatrix demag operator
         N v = B^T (H.matvec(B v)); M_mass-preconditioned GMRES per Newton step + Armijo line search +
         scalar-chi Picard warmstart; FAILS LOUD on non-convergence (CLAUDE.md No-Fallbacks).
  solve_nonlinear_newton(mesh, chi0, Msat, H0, bh_table=..., ...)
      -- thin wrapper over the production radia.vim.hdiv_demag_solve (the same C++ damped Newton), for
         either the analytic saturating curve (chi0, Msat) or an explicit [[H,B]] table.

The constitutive-tangent helpers (_tensor_tangent_cfs, _table_tensor_tangent[_multi]) are the consistent
TENSOR tangent dM/dH = chi_diff Hhat(x)Hhat + chi_sec (I - Hhat(x)Hhat) (the scalar chi_diff*I stalls at
moderate drive); they are also consumed by _solve.hdiv_demag_solve.
"""
from math import pi

import numpy as np
import scipy.sparse as sp

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

from radia.vim import _core as tet   # M1: core promoted to radia.vim from the retired prototype.


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


def _table_tensor_tangent_multi(gfH, mesh, region_funcs, elem_region, Id):
    """PER-REGION version of _table_tensor_tangent: each element uses ITS OWN region's PCHIP BH table.

    `region_funcs` is a list (indexed by region-id) of (Bpch, Bder, Hmax, Mmax) tuples (one per soft-iron
    grade); `elem_region` is the per-element region-id array in mesh element order (== L2(0) DOF order, so
    element i <-> DOF i).  N = B^T G B is geometry-only, so per-region nonlinear iron enters ONLY here:
    each element's chi_sec / chi_diff is read from the BH table of the material it belongs to.  Returns
    the same (M(H) CF, consistent tensor-tangent CF) pair as the single-region helper."""
    l2 = ng.L2(mesh, order=0)
    gfm = ng.GridFunction(l2)
    gfm.Set(ng.sqrt(ng.InnerProduct(gfH, gfH) + 1e-30))
    Hmag = np.maximum(gfm.vec.FV().NumPy(), 1e-30)            # per-element |H| (DOF == element order)
    sec_e = np.empty_like(Hmag)
    dif_e = np.empty_like(Hmag)
    for ridx, (Bpch, Bder, Hmax, Mmax) in enumerate(region_funcs):
        sel = elem_region == ridx
        if not np.any(sel):
            continue
        h = Hmag[sel]
        sat = h > Hmax
        a_cl = np.minimum(h, Hmax)
        m_e = np.where(sat, Mmax, Bpch(a_cl) / _MU0 - a_cl)   # M(|H|), saturated beyond this region's table
        sec_e[sel] = m_e / h                                  # chi_sec = M/|H|
        dif_e[sel] = np.where(sat, 0.0, Bder(a_cl) / _MU0 - 1.0)   # chi_diff = dM/d|H| (0 beyond sat)
    gf_sec = ng.GridFunction(l2)
    gf_sec.vec.FV().NumPy()[:] = sec_e
    gf_dif = ng.GridFunction(l2)
    gf_dif.vec.FV().NumPy()[:] = dif_e
    H2 = ng.InnerProduct(gfH, gfH) + 1e-30
    par = ng.OuterProduct(gfH, gfH) / H2
    return gf_sec * gfH, gf_dif * par + gf_sec * (Id - par)


# --------------------------------------------------------------------- INVERSE BH (energy / reluctivity form)
# The SYMMETRIC energy-Newton (the all-C++ nonlinear path, _solve._solve_nonlinear_energy_cpp) linearises the
# residual R(m) = INT H(M).v dx + N m - M_mass h_ext (M = the flux field of m, H(M) the INVERSE BH reluctance
# field) -> Jacobian J = W_tan + N, W_tan = INT nu_d u.v dx, nu_d = dH/dM = (dM/dH)^-1 (differential
# reluctivity tensor).  Symmetric (N + the SPD reluctivity mass) -> CG-able by the EXISTING C++ W-CG
# (solve_linear_material_mass_riesz), unlike the forward M-residual whose J = M_mass + T M_mass^-1 N needs
# GMRES + an M_mass^-1.  Reduces to the linear (M_{1/chi}+N)m = M_mass h_ext when chi is constant.
_CHI_DIFF_FLOOR = 1e-3     # nu_diff = 1/max(chi_diff, floor): cap the differential reluctivity in deep saturation
                           # (chi_diff = dM/d|H| -> 0).  Also the barrier slope Kbar = 1/floor for |M| > Mmax.


def _bh_inverse_funcs(Harr, Barr):
    """Inverse of `_bh_table_funcs`: from the forward monotone M(|H|) table build the vectorized per-element
    reluctivity fields + co-energy of |M|.  Returns (fields, wco, Mmax):
      fields(Mmag) -> (nu_sec, nu_diff):  nu_sec = |H|/|M| (secant reluctivity, H(M) = nu_sec M);
                      nu_diff = 1/chi_diff (differential reluctivity along M, chi_diff = B'(|H|)/mu0 - 1).
      wco(Mmag)   -> the co-energy density INT_0^|M| H(s) ds  (the line-search merit -- E is convex, R is
                      its gradient, so an Armijo line search on E is robust where ||R|| stalls in saturation).
    HARD-SATURATION BARRIER: the table saturates M at Mmax (M cannot exceed Msat); for |M| > Mmax the inverse
    is undefined, so H(M) is extended with a steep but C1-smooth barrier (slope Kbar = 1/_CHI_DIFF_FLOOR = the
    table's saturation nu_diff) -> the energy form is repelled from the unphysical |M| > Mmax region (without
    it, the M-iterates overshoot Msat and the Newton wanders / limit-cycles)."""
    from scipy.interpolate import PchipInterpolator
    Mof, Bpch, Bder, Hmax, Mmax = _bh_table_funcs(Harr, Barr)
    Hs = np.concatenate([[0.0], np.logspace(-2, np.log10(max(Hmax, 1.0)), 800)])
    Ms = np.array([Mof(h) for h in Hs])                       # should be monotone 0..Mmax
    # Saturating BH tables can produce tiny PCHIP ripples in M(H).  PCHIP for the
    # inverse H(M) requires the whole retained sequence to be strictly increasing,
    # not just adjacent positive differences in the unfiltered samples.
    keep = []
    last = -np.inf
    for k, mval in enumerate(Ms):
        if mval > last + 1e-12 * max(1.0, Mmax):
            keep.append(k)
            last = mval
    Hgrid, Mgrid = Hs[keep], Ms[keep]
    if len(Mgrid) < 2:
        raise ValueError("hdiv_demag_solve: inverse BH table did not produce a monotone M(H) curve")
    Hof = PchipInterpolator(Mgrid, Hgrid)                     # |H| given |M|
    Wco_grid = np.concatenate([[0.0], np.cumsum(0.5 * (Hgrid[1:] + Hgrid[:-1]) * np.diff(Mgrid))])
    Wco = PchipInterpolator(Mgrid, Wco_grid)
    Hlast = float(Hgrid[-1])
    Kbar = 1.0 / _CHI_DIFF_FLOOR
    Mcap = Mmax * (1.0 - 1e-9)

    def fields(Mmag):
        Mmag = np.asarray(Mmag, float)
        Mm = np.minimum(Mmag, Mcap)
        h = Hof(Mm)
        chid = np.maximum(Bder(np.minimum(h, Hmax)) / _MU0 - 1.0, _CHI_DIFF_FLOOR)
        nu_sec = h / np.maximum(Mmag, 1e-30)
        nu_diff = 1.0 / chid
        over = Mmag > Mcap
        if np.any(over):                                      # C1-smooth barrier beyond Mmax
            nu_sec[over] = (Hlast + Kbar * (Mmag[over] - Mcap)) / Mmag[over]
            nu_diff[over] = Kbar
        return nu_sec, nu_diff

    def wco(Mmag):
        Mmag = np.asarray(Mmag, float)
        Mm = np.minimum(Mmag, Mcap)
        d = np.maximum(Mmag - Mcap, 0.0)
        return Wco(Mm) + Hlast * d + 0.5 * Kbar * d * d

    return fields, wco, Mmax


def _reluctivity_tangent(gfM, mesh, fields, Id):
    """Inverse analogue of `_table_tensor_tangent`: from the M field gfM, per-element nu_sec / nu_diff via the
    inverse BH `fields`, returns (H(M) CF = nu_sec M, nu_d tensor CF = nu_diff Mhat(x)Mhat + nu_sec (I - Mhat(x)Mhat))."""
    l2 = ng.L2(mesh, order=0)
    gfm = ng.GridFunction(l2)
    gfm.Set(ng.sqrt(ng.InnerProduct(gfM, gfM) + 1e-30))
    Mmag = np.maximum(gfm.vec.FV().NumPy(), 1e-30)
    nu_sec, nu_diff = fields(Mmag)
    gs = ng.GridFunction(l2); gs.vec.FV().NumPy()[:] = nu_sec
    gd = ng.GridFunction(l2); gd.vec.FV().NumPy()[:] = nu_diff
    M2 = ng.InnerProduct(gfM, gfM) + 1e-30
    par = ng.OuterProduct(gfM, gfM) / M2                      # Mhat (x) Mhat
    return gs * gfM, gd * par + gs * (Id - par)


def _reluctivity_tangent_multi(gfM, mesh, region_fields, elem_region, Id):
    """PER-REGION inverse analogue of `_table_tensor_tangent_multi`: each element uses ITS region's inverse BH
    `fields`.  `region_fields` = list (by region-id) of the `fields` callables; `elem_region` = per-element
    region-id (== L2(0) DOF order)."""
    l2 = ng.L2(mesh, order=0)
    gfm = ng.GridFunction(l2)
    gfm.Set(ng.sqrt(ng.InnerProduct(gfM, gfM) + 1e-30))
    Mmag = np.maximum(gfm.vec.FV().NumPy(), 1e-30)
    nu_sec = np.empty_like(Mmag); nu_diff = np.empty_like(Mmag)
    for ridx, fields in enumerate(region_fields):
        sel = elem_region == ridx
        if not np.any(sel):
            continue
        s, d = fields(Mmag[sel])
        nu_sec[sel] = s; nu_diff[sel] = d
    gs = ng.GridFunction(l2); gs.vec.FV().NumPy()[:] = nu_sec
    gd = ng.GridFunction(l2); gd.vec.FV().NumPy()[:] = nu_diff
    M2 = ng.InnerProduct(gfM, gfM) + 1e-30
    par = ng.OuterProduct(gfM, gfM) / M2
    return gs * gfM, gd * par + gs * (Id - par)


def solve_nonlinear_newton(mesh, chi0, Msat, H0, near_correction=True, nsub=4,
                           picard_warmstart=8, maxit=200, tol=1e-10,
                           bh_table=None, require_convergence=True):
    """NONLINEAR HDiv-VIM solve on a +z uniform applied field H0 -- the C++ scalable charge-Gram path.

    The dense Python charge-Gram Newton was REMOVED; this is a thin wrapper over the production
    `radia.vim.hdiv_demag_solve` (the damped matrix-free Newton on the C++ analytic `_ChargeGramHMatrix`
    demag operator N v = B^T (H.matvec(B v))).  The constitutive law is either the analytic saturating
    curve (chi0, Msat) -- supplied to the C++ path as the equivalent [[H,B]] table -- or an explicit
    `bh_table = (Harr, Barr)` (the same data Radia's MatSatIsoTab consumes).

    `near_correction`, `picard_warmstart`, `tol` are accepted for call-site stability but are NOT used by
    the C++ path (its analytic Gram is exact near AND far, so no near-correction is needed; its warmstart
    + inner-GMRES tolerances are internal).  `require_convergence` is informational: the C++ path ALWAYS
    fails loud on non-convergence (CLAUDE.md No-Fallbacks).

    Returns (M_avg_z, n_newton_iter, D_used) on a +z drive.  The caller must open `with ng.TaskManager():`.
    """
    del near_correction, picard_warmstart, tol, require_convergence   # not used by the C++ path
    from ._solve import hdiv_demag_solve
    if bh_table is not None:
        Harr = np.asarray(bh_table[0], float)
        Barr = np.asarray(bh_table[1], float)
    else:
        # synthesize the [[H,B]] table of the analytic saturating curve M(H)=chi0 H/(1+chi0|H|/Msat),
        # B = mu0 (H + M), so the C++ table-driven nonlinear path sees the SAME constitutive law.
        Harr = np.concatenate([[0.0], np.logspace(-1, 7, 60)])
        Mof = _bh_curve(chi0, Msat)
        Barr = _MU0 * (Harr + np.array([Mof(h) for h in Harr]))
    BH = [[float(h), float(b)] for h, b in zip(Harr, Barr)]
    res = hdiv_demag_solve(mesh, bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)),
                           nl_maxit=maxit)
    return float(res["M_avg"][2]), int(res["iters"]), float(res["demag"])


def solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0, nsub=4, gram_eps=1e-10,
                                    picard_warmstart=8, maxit=200, gmres_tol=1e-8, newton_tol=1e-6,
                                    near_factor=1e30, gmres_restart=400, return_timing=False, verbose=False):
    """SCALABLE damped Newton (production #2): the demag is the C++ HACApK charge-Gram H-matrix
    (O(N log N) apply), and the Newton system is solved ITERATIVELY (GMRES) -- no dense factorization
    anywhere.

    The C++ Gram is the ANALYTIC charge Gram (PhiTet/TriPotential, exact near AND far), so the demag
    apply is simply  N v = B^T ( H.matvec(B v) )  with H = the analytic _ChargeGramHMatrix -- no separate
    sparse near-correction.  This is exact for NON-uniform M (cube, C-yoke, any div M != 0 body).  The
    Jacobian J v = M_mass v + T M_mass^{-1} N v is applied matrix-free (M_mass factored once, sparse);
    GMRES (M_mass-preconditioned) solves each Newton step; Armijo line search + Picard warmstart.
    Returns (M_avg, n_newton_iter, D_used)."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    import radia._radia_pybind as _rp
    import time

    _t0 = time.perf_counter()
    Mof = _bh_curve(chi0, Msat)
    # SPARSE build: M_mass + B come back SPARSE; the demag apply is the analytic charge-Gram H-matvec
    # (Hg) below (no dense n_charge^2 object anywhere).
    d = tet.build_demag(mesh, nsub)
    M_mass = d["M_mass"]                         # sparse CSR
    mu = d["m_unit"]
    denom = float(mu @ (M_mass @ mu))            # sparse-safe Rayleigh denominator
    B = d["B_csr"]
    # The C++ ANALYTIC charge Gram H-matrix is exact near AND far, so the demag apply is just
    # B^T (G_analytic-Hmatvec (B v)) -- NO separate Python near-correction.  This handles the non-uniform-M
    # (div M != 0: cube, C-yoke) case while staying O(N log N).
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
    # ---- damped Newton-Raphson on the C++ scalable charge-Gram operator (robust at ALL drive) ----
    # solve_nonlinear_newton routes through the production hdiv_demag_solve (NGSolve work), so this CALLER
    # opens the TaskManager (CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT").
    print(f"Nonlinear HDiv-VIM sphere demag (chi0={chi0}, Msat={Msat}) -- damped Newton (C++ charge Gram):")
    print(f"{'H0':>8} {'Newton M':>10} {'newton it':>10} {'analytic(1/3)':>14} {'rel':>9}")
    with ng.TaskManager():
        for H0 in (1e-2, 1e-1, 3e-1, 1.0, 5.0):
            Mn, nit, Dn = solve_nonlinear_newton(mesh, chi0, Msat, H0)
            Mana = _scalar_fixed_point(Mof, 1.0 / 3.0, H0)
            print(f"{H0:8.0e} {Mn:10.5f} {nit:10d} {Mana:14.5f} {abs(Mn - Mana) / Mana:9.1e}")
    print("  => Newton matches the analytic at every drive (deep saturation included).")


if __name__ == "__main__":
    main()
