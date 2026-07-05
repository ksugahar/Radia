"""radia.planar_hysteresis -- SHARED 2D play-hysteresis quasi-static demag solver.

A soft-magnetic body with a Prandtl-Ishlinskii play law (``planar_materials.PlayHysteresis``) is stepped
through a sequence of applied fields; at each step the demagnetising fixed point M = play(H0 + N M) is
solved by NEWTON on the dense demag operator N (``planar_aniso.demag_operator``, assembled on the SHARED
planar_charges kernel).  The play's INCREMENTAL susceptibility (>= 0 even on the descending branch) keeps
the Newton Jacobian I - diag(chi_inc) N well-conditioned -- a secant-chi Picard would see negative chi
near coercivity and break.  Because N comes from the shared kernel, this stays aligned with the
HDiv-VIM planar layer.

UNIAXIAL: the hysteresis axis is x (field applied along +x, M along x, M_y = 0 by symmetry for a body
symmetric about the x-axis) -- the standard first 2D hysteresis model.  Verified: the anhysteretic limit
(eta -> 0) recovers the linear demag chi/(1+chi/2), and the full N-Newton disk loop matches the scalar
reference M = play(H_ext - D M) (D = 1/2) to ~1e-4.

    import radia.planar_hysteresis as ph
    from radia.planar_materials import PlayHysteresis
    play = PlayHysteresis(eta=[0.1,0.3,0.6,1.0], w=[3,2,1.2,0.6])
    H = 3.0 * np.sin(np.linspace(0, 4*np.pi, 97))          # 2 cycles
    r = ph.solve_hysteresis_demag(mesh, play, H)
    r["M_avg"]   # the (H_ext, M) hysteresis loop (volume-averaged)
"""
from __future__ import annotations

import numpy as np

from radia.planar_geometry import _extract_geometry
from radia.planar_aniso import demag_operator

MU0 = 4e-7 * np.pi


def _newton(play, Nxx, H0, p, m_init, tol, maxit):
    """Solve M = play.M(H0 + Nxx M, p) by Newton (committed state p fixed); returns M (no commit)."""
    n = len(H0)
    M = m_init.copy()
    for _ in range(maxit):
        H = H0 + Nxx @ M
        F = M - play.M(H, p)
        if np.linalg.norm(F) < tol * max(np.linalg.norm(M), 1.0):
            break
        J = np.eye(n) - (play.chi_inc(H, p)[:, None] * Nxx)
        M = M - np.linalg.solve(J, F)
    else:
        raise RuntimeError("planar_hysteresis: Newton did NOT converge in %d iters (residual %.2e) -- "
                           "returning M would be a silent wrong result" % (maxit, np.linalg.norm(F)))
    return M


def solve_hysteresis_demag(mesh, play, H_ext_seq, ngauss=6, newton_tol=1e-10, newton_maxit=50):
    """Quasi-static uniaxial (x-axis) play-hysteresis demag over the applied-field sequence
    ``H_ext_seq`` (A/m along +x).  ``play`` is a planar_materials.PlayHysteresis.

    Returns dict: H_ext (nStep,), M_avg (nStep,) the volume-averaged M loop, M (n,) the final
    per-element magnetisation, n_el, ndof(=n).  The caller wraps TaskManager (mesh construction)."""
    _, _, centroids, areas = _extract_geometry(mesh)
    n = len(areas)
    Nxx = demag_operator(mesh, centroids, ngauss)[0::2, 0::2]     # x-field from x-magnetisation block
    w = areas / areas.sum()
    p = play.fresh_state(n)
    M = np.zeros(n)
    M_avg = np.empty(len(H_ext_seq))
    for i, Hext in enumerate(H_ext_seq):
        H0 = np.full(n, float(Hext))
        M = _newton(play, Nxx, H0, p, M, newton_tol, newton_maxit)
        p = play.advance(H0 + Nxx @ M, p)                        # commit the play state at this step
        M_avg[i] = w @ M
    return {"H_ext": np.asarray(H_ext_seq, float), "M_avg": M_avg, "M": M, "n_el": n, "ndof": n}
