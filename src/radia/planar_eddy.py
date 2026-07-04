"""radia.planar_eddy -- SHARED 2D staggered eddy-current coupling (maglev / induction machine).

Weak (sequential) coupling of an ANALYTIC soft-iron demagnetisation solver -- either the collocation
MMMM (``radia.mmmm2d``) or the HDiv-VIM (``radia.vim._vim2d``), both of which produce a per-element
magnetisation M -- with an NGSolve REDUCED-POTENTIAL complex A_z eddy FEM in a SEPARATE conductor.
This is the 2D instance of the lab's maglev method (CLAUDE.md "Maglev Analysis: Radia + NGSolve"):
Radia/MMMM supplies the open-boundary iron field ANALYTICALLY (no air mesh, magnet may move without
re-meshing); NGSolve solves ONLY the eddy reaction in the conductor.  References: Chadebec 2006 (IEM
open boundary), Biro 2000 (reduced potential).

The coupling is METHOD-AGNOSTIC: ``iron_solve`` is a callback ``H_ext_complex (nEl,2) -> M (nEl,2)``,
so the SAME orchestration serves MMMM and the HDiv-VIM.  The eddy system matrix (grad.grad +
j w mu0 sigma mass) is angle/iteration INDEPENDENT -- only the RHS (the applied + iron source A0)
changes -- so it is assembled + factored ONCE and back-substituted every staggered iteration.

De-risked (C:\\temp\\mmmm2d_eddy) against a MONOLITHIC FEM (iron+conductor+air meshed together): the
staggered result reproduces it to 6e-4 (iron M_avg) / 1.6e-4 (conductor <Bx>) in ~4 iterations, and
the standalone eddy FEM matches the analytic conducting-cylinder Bessel gate to ~1e-4.

TaskManager: this is a HELPER module (caller-wraps policy) -- it never opens ``with TaskManager()``;
call it inside the caller's TaskManager region (the eddy FEM Assemble/Inverse run parallel there).

    import ngsolve as ng
    import radia.planar_eddy as pe
    with ng.TaskManager():
        res = pe.couple_mmmm(iron_mesh, mu_r=100.0, fem_mesh=cond_air_mesh, sigma=3.7e7, freq=150.0)
    res["M"]          # complex per-element iron magnetisation (nEl, 2)
    res["gfu"]        # NGSolve complex GridFunction: the eddy reaction potential A_r
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, TaskManager,  # noqa: F401 (re-export)
                     grad, dx, x, y, atan2, CF)

import radia._radia_pybind as _rp
from radia.mmmm2d import (_extract_geometry, _element_materials, _sub_geometry, _pm_hard_M,
                          _region_chi_for)
from radia.planar_charges import mn_edge_cloud, exterior_field

MU0 = 4e-7 * np.pi


# ---- iron field as an NGSolve source (shared: the same M.n log-charge cloud planar_charges uses) --

def iron_field_cf(iron_mesh, M_elem, ngauss=2):
    """The iron's out-of-plane vector potential A_z as an NGSolve CoefficientFunction:
        A_z(x,y) = mu0/(2 pi) sum_q Q_q atan2(y - y_q, x - x_q),
    Q_q the shared M.n edge-charge cloud (``planar_charges.mn_edge_cloud``; Q complex when M is a
    phasor).  Exact (no interpolation), so it can be used directly in the eddy FEM RHS.

    BRANCH-CUT: atan2 cuts along -x from each charge; place the iron so the conductor lies to one side
    (the net cloud charge is zero, so a global 2 pi jump cancels, but a cut through the conductor would
    corrupt the RHS integral).  The total field B = curl A_z is cut-free regardless.
    """
    Xq, Q = mn_edge_cloud(iron_mesh, M_elem, ngauss)      # Q complex if M complex
    k = MU0 / (2.0 * np.pi)
    acc = CF(0j)
    for i in range(len(Q)):
        acc = acc + complex(Q[i]) * k * atan2(y - float(Xq[i, 1]), x - float(Xq[i, 0]))
    return acc


# ---- reduced-potential complex A_z eddy FEM (factored once) --------------------------------------

def eddy_operator(fem_mesh, sigma, freq, order=4, conductor="conductor", dirichlet="outer"):
    """Assemble + factor the reduced-A_z eddy system ONCE.  Returns (fes, inv, w).

    System (iteration-independent): int grad(u).grad(v) + j w mu0 sigma int_cond u v.  The RHS
    (-j w mu0 sigma int_cond A0 v) is rebuilt each staggered step; see ``eddy_reaction``."""
    w = 2.0 * np.pi * freq
    fem_mesh.Curve(order)
    fes = H1(fem_mesh, order=order, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    af = BilinearForm(fes, symmetric=True)
    af += grad(u) * grad(v) * dx
    af += 1j * w * MU0 * sigma * u * v * dx(conductor)
    af.Assemble()
    inv = af.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
    return fes, inv, w


def eddy_reaction(fes, inv, w, sigma, applied_Az_cf, conductor="conductor"):
    """One reduced eddy solve for a given source potential A0 = applied_Az_cf (applied + iron).
    Returns the complex GridFunction A_r (the eddy reaction)."""
    v = fes.TestFunction()
    lf = LinearForm(fes)
    lf += -1j * w * MU0 * sigma * applied_Az_cf * v * dx(conductor)
    lf.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = inv * lf.vec
    return gfu


# ---- staggered coupling (method-agnostic) -------------------------------------------------------

def couple(iron_mesh, iron_solve, fem_mesh, sigma, freq, *, applied_Az_cf, H_app,
           order=4, tol=1e-6, maxit=40, ngauss=2, conductor="conductor"):
    """Staggered MMMM/VIM-iron <-> reduced-eddy-FEM coupling.

    iron_solve : callable(H_ext_complex (nEl,2)) -> M (nEl,2) complex.  The soft-iron demag with a
                 COMPLEX applied field; MMMM and the HDiv-VIM both fit (see couple_mmmm).
    applied_Az_cf : the APPLIED-field vector potential A_z as an NGSolve CF (e.g. B0*y for a uniform
                 B0 x-hat), harmonic in the FEM domain.
    H_app      : applied H (nEl,2, complex) at the iron centroids, consistent with applied_Az_cf.

    Returns dict: M, M_avg, iters, hist (per-step relative M change), gfu (eddy reaction A_r), fes.
    """
    _, _, centroids, areas = _extract_geometry(iron_mesh)
    nEl = len(areas)
    fes, inv, w = eddy_operator(fem_mesh, sigma, freq, order=order, conductor=conductor)
    M = np.zeros((nEl, 2), complex)
    hist = []
    for it in range(maxit):
        A0 = applied_Az_cf + iron_field_cf(iron_mesh, M, ngauss=ngauss)
        gfu = eddy_reaction(fes, inv, w, sigma, A0, conductor=conductor)
        gB = grad(gfu)
        H_eddy = np.empty((nEl, 2), complex)
        for i, (px, py) in enumerate(centroids):
            g = gB(fem_mesh(px, py))                       # eddy reaction B_r = curl A_r, at the iron
            H_eddy[i] = (g[1] / MU0, -g[0] / MU0)          # H = B/mu0 (Bx=dA_r/dy, By=-dA_r/dx)
        M_new = iron_solve(np.asarray(H_app, complex) + H_eddy)
        d = np.linalg.norm(M_new - M) / max(np.linalg.norm(M_new), 1e-300)
        hist.append(float(d))
        M = M_new
        if d < tol:
            it += 1
            break
    else:
        raise RuntimeError("planar_eddy.couple: staggered iteration NOT converged "
                           "(rel dM=%.2e after %d iters -- returning M would be a silent wrong result)"
                           % (hist[-1] if hist else np.inf, maxit))
    wa = areas / areas.sum()
    M_avg = np.array([wa @ M[:, 0], wa @ M[:, 1]])
    return {"M": M, "M_avg": M_avg, "iters": it, "hist": hist, "gfu": gfu, "fes": fes, "n_el": nEl}


def _mmmm_iron_solve(iron_mesh, mu_r=None, bh_table=None, pm=None,
                     nl_tol=1e-6, nl_maxit=200, nl_damp=0.6, chi_floor=1e-12):
    """Build a MMMM complex iron_solve callback for the staggered eddy coupling.

    * LINEAR (``mu_r`` scalar or {region: mu_r}): chi real -> M = solve(Re H) + j solve(Im H).
    * NONLINEAR (``bh_table`` [[H,B],...]): amplitude-based EFFECTIVE-chi Picard -- 1st-harmonic AC
      approximation: chi_eff = M(|H|)/|H| with |H| the PHASOR MAGNITUDE sqrt(|Hx|^2+|Hy|^2) (which
      reduces to |H| for a real/DC field, so the sigma->0 limit recovers the DC nonlinear demag).
      The chi state warm-starts across staggered calls.  (Captures amplitude-dependent saturation,
      NOT harmonic generation -- that needs time stepping.)

    With ``pm`` = {region: [Mx,My]} the mesh is a PM-motor / ECB rotor: those regions are RIGID
    permanent magnets (fixed REAL M, an in-phase phasor source at omega).  Their field magnetises the
    soft iron (added to the soft drive) AND, because the returned M carries the fixed PM elements,
    drives the conductor eddy through iron_field_cf.  Only the soft subsystem is solved."""
    if (mu_r is None) == (bh_table is None):
        raise ValueError("planar_eddy: provide EXACTLY ONE of mu_r (linear) or bh_table (nonlinear)")
    verts, offsets, centroids, areas = _extract_geometry(iron_mesh)
    nEl = len(areas)
    mats = _element_materials(iron_mesh)
    if pm:
        hard = set(pm)
        soft_ids = np.array([i for i, m in enumerate(mats) if m not in hard], int)
        hard_ids = np.array([i for i, m in enumerate(mats) if m in hard], int)
        if soft_ids.size == 0:
            raise ValueError("planar_eddy: pm covers every region -- no soft iron to solve")
        law_regions = set(mu_r) if isinstance(mu_r, dict) else set()
        if hard & law_regions:
            raise ValueError("planar_eddy: region(s) %s are BOTH pm and soft mu_r"
                             % sorted(hard & law_regions))
        verts_s, offsets_s = _sub_geometry(verts, offsets, soft_ids)
        M_hard = _pm_hard_M(iron_mesh, pm, mats, nEl)
        H_pm = exterior_field(iron_mesh, M_hard, centroids[soft_ids])     # PM field at soft (real)
        mats_s = [mats[i] for i in soft_ids]
    else:
        soft_ids = np.arange(nEl)
        hard_ids = np.empty(0, int)
        verts_s, offsets_s, M_hard, H_pm, mats_s = verts, offsets, None, 0.0, mats

    def _lin(Hs, chi):
        Mr = _rp.Moment2DSolveLinear(verts_s, offsets_s, chi, np.ascontiguousarray(Hs.real, float))
        Mi = _rp.Moment2DSolveLinear(verts_s, offsets_s, chi, np.ascontiguousarray(Hs.imag, float))
        return Mr + 1j * Mi

    def _assemble(Mc):
        M = np.zeros((nEl, 2), complex)
        M[soft_ids] = Mc
        if pm:
            M[hard_ids] = M_hard[hard_ids]
        return M

    if mu_r is not None:                                                  # LINEAR
        if isinstance(mu_r, dict):
            chi = _region_chi_for(mats_s, mu_r)
        else:
            if not mu_r > 1.0:
                raise ValueError("planar_eddy: mu_r must be > 1 (got %r)" % (mu_r,))
            chi = np.full(len(soft_ids), mu_r - 1.0)

        def solve(H):
            return _assemble(_lin(np.ascontiguousarray(H[soft_ids]) + H_pm, chi))
        return solve

    # NONLINEAR: scalar bh_table, effective-chi Picard (warm-started via closure state)
    if isinstance(bh_table, dict):
        raise NotImplementedError("planar_eddy: per-region bh_table not yet wired; use a scalar table")
    from radia.mmmm2d import _law_from_table
    M_of_h, _chi_sec, chi0 = _law_from_table(bh_table)
    state = {"chi": np.full(len(soft_ids), chi0)}

    def solve(H):
        Hs = np.ascontiguousarray(H[soft_ids]) + H_pm                     # applied + eddy + PM (complex)
        chi = state["chi"]
        for _ in range(nl_maxit):
            Mc = _lin(Hs, np.maximum(chi, chi_floor))
            Mmag = np.sqrt(np.abs(Mc[:, 0]) ** 2 + np.abs(Mc[:, 1]) ** 2)  # |M phasor| (== |M| at DC)
            nH = np.maximum(Mmag / np.maximum(chi, chi_floor), 1e-300)     # |H_local| = |M|/chi
            chi_star = np.maximum(M_of_h(nH) / nH, chi_floor)
            r = chi_star - chi
            if np.linalg.norm(r) / max(np.linalg.norm(chi_star), 1e-300) < nl_tol:
                chi = chi_star
                break
            chi = np.maximum(chi + nl_damp * r, chi_floor)
        state["chi"] = chi
        return _assemble(Mc)
    return solve


def couple_mmmm(iron_mesh, fem_mesh, sigma, freq, *, mu_r=None, bh_table=None, pm=None, B0=1.0,
                order=4, tol=1e-6, maxit=40, ngauss=2, conductor="conductor"):
    """Convenience: couple a MMMM soft-iron body to the eddy FEM under a uniform applied field
    B0 x-hat (A0 = B0 y).  Soft-iron law: EXACTLY ONE of ``mu_r`` (scalar or {region: mu_r} dict,
    LINEAR) or ``bh_table`` ([[H,B],...], NONLINEAR effective-chi).  ``pm={region: [Mx,My]}`` adds
    embedded rigid permanent magnets (PM-motor / eddy-current-brake rotor).  Returns couple()."""
    _, _, centroids, _ = _extract_geometry(iron_mesh)
    H_app = np.tile([B0 / MU0, 0.0], (len(centroids), 1)).astype(complex)
    solve = _mmmm_iron_solve(iron_mesh, mu_r=mu_r, bh_table=bh_table, pm=pm)
    return couple(iron_mesh, solve, fem_mesh, sigma, freq, applied_Az_cf=B0 * y, H_app=H_app,
                  order=order, tol=tol, maxit=maxit, ngauss=ngauss, conductor=conductor)
