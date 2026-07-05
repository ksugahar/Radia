"""radia.planar_eddy -- 2D staggered eddy-current coupling (maglev / induction machine).

Weak (sequential) coupling of a planar soft-iron demagnetisation callback, which produces a per-element
magnetisation M, with an NGSolve REDUCED-POTENTIAL complex A_z eddy FEM in a SEPARATE conductor.
This is the 2D instance of the lab's maglev method (CLAUDE.md "Maglev Analysis: Radia + NGSolve"):
Radia supplies the open-boundary iron field analytically (no air mesh, magnet may move without
re-meshing); NGSolve solves ONLY the eddy reaction in the conductor.  References: Chadebec 2006
(IEM open boundary), Biro 2000 (reduced potential).

The coupling is METHOD-AGNOSTIC: ``iron_solve`` is a callback ``H_ext_complex (nEl,2) -> M (nEl,2)``,
so the SAME orchestration serves the HDiv-VIM and any future planar reduced model.  The eddy system
matrix (grad.grad + j w mu0 sigma mass) is angle/iteration INDEPENDENT -- only the RHS (the applied +
iron source A0) changes -- so it is assembled + factored ONCE and back-substituted every staggered
iteration.

TaskManager: this is a HELPER module (caller-wraps policy) -- it never opens ``with TaskManager()``;
call it inside the caller's TaskManager region (the eddy FEM Assemble/Inverse run parallel there).

    import ngsolve as ng
    import radia.planar_eddy as pe
    with ng.TaskManager():
        res = pe.couple(iron_mesh, iron_solve, cond_air_mesh, sigma=3.7e7, freq=150.0,
                        applied_Az_cf=B0*y, H_app=H_app)
    res["M"]          # complex per-element iron magnetisation (nEl, 2)
    res["gfu"]        # NGSolve complex GridFunction: the eddy reaction potential A_r
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, TaskManager,  # noqa: F401 (re-export)
                     grad, dx, x, y, atan2, CF)

from radia.planar_geometry import _extract_geometry
from radia.planar_charges import mn_edge_cloud

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
    """Staggered planar-iron <-> reduced-eddy-FEM coupling.

    iron_solve : callable(H_ext_complex (nEl,2)) -> M (nEl,2) complex.  The soft-iron demag callback
                 sees the applied, eddy, and source field at iron element centroids.
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
