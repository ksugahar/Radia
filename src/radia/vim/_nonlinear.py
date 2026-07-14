"""Constitutive-law helpers for the production HDiv-VIM nonlinear solver.

The public nonlinear entry point is :func:`radia.vim.Solve`.  This module
contains the constitutive curves and consistent tangent helpers consumed by
``radia.vim._solve``; it contains constitutive helpers only.
"""
from math import pi

import numpy as np
import ngsolve as ng


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
        raise ValueError("vim.Solve: inverse BH table did not produce a monotone M(H) curve")
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
