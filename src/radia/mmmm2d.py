"""radia.mmmm2d -- 2D planar collocation MMMM (Multipole Magnetic Moment Method).

The 2D twin of the 3D MMMM (rad.Solve moment path) and the sibling of the 2D HDiv-VIM
(radia.vim._vim2d): a per-unit-length motor-cross-section soft-iron demag solver on a mesh of
triangles / quadrilaterals with the 2D Laplace kernel G = -ln(r)/(2 pi).  Each element carries one
uniform line-charge DOF per EDGE; the constitutive law M = chi H is imposed on the field MOMENTS
about the element centroid (1 monopole + 2 dipole + (nEdge-3) quadrupole rows).

The NUMERICAL CORE is C++ (rad_moment2d.cpp, exposed as radia._radia_pybind.Moment2DSolveLinear):
the segment log-kernel field/gradient closed forms, the moment-system assembly, and the dense LU
solve.  This module is the thin NGSolve-mesh ingestion + linear / nonlinear (scalar-chi Picard)
driver + demag-factor postprocessing.

API mirrors radia.vim._vim2d.solve_planar_demag.

TaskManager: per the caller-wraps policy this module never opens a TaskManager (the C++ solve does
not use it); wrap the mesh construction in the caller's `with TaskManager():` as usual.
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng

import radia._radia_pybind as _rp

MU0 = 4e-7 * np.pi


def _extract_geometry(mesh):
    """NGSolve 2D mesh -> (verts (nVert,2) f64, offsets (nElem+1,) i32, centroids (nElem,2),
    areas (nElem,)).  Vertices are passed in NGSolve order; the C++ re-orients each element CCW.
    The centroid/area (used for the applied field + averaging) are orientation-independent."""
    if mesh.dim != 2:
        raise ValueError("mmmm2d: mesh.dim must be 2 (got %d)" % mesh.dim)
    pts = np.array([list(mesh[v].point)[:2] for v in mesh.vertices])
    verts, offsets, centroids, areas = [], [0], [], []
    for el in mesh.Elements(ng.VOL):
        V = pts[[v.nr for v in el.vertices]]
        x, y = V[:, 0], V[:, 1]
        A2 = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)     # signed 2*area
        # area centroid (orientation-independent: numerator & A2 both flip under reversal)
        cx = np.sum((x + np.roll(x, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3 * A2)
        cy = np.sum((y + np.roll(y, -1)) * (x * np.roll(y, -1) - np.roll(x, -1) * y)) / (3 * A2)
        verts.append(V)
        offsets.append(offsets[-1] + len(V))
        centroids.append((cx, cy))
        areas.append(0.5 * abs(A2))
    return (np.vstack(verts).astype(np.float64), np.asarray(offsets, dtype=np.int32),
            np.asarray(centroids, float), np.asarray(areas, float))


def _eval_Hext(H_ext, centroids, mesh):
    """Applied field at the element centroids -> (nElem, 2).  H_ext is a 2-tuple/array (uniform)
    or a 2-component NGSolve CoefficientFunction (evaluated at the centroids)."""
    arr = np.asarray(H_ext, dtype=object) if not isinstance(H_ext, ng.CoefficientFunction) else None
    if arr is not None:
        a = np.asarray(H_ext, float)
        if a.shape == (2,):
            return np.tile(a, (len(centroids), 1))
        if a.shape == (len(centroids), 2):
            return a
        raise ValueError("mmmm2d: H_ext array must be (2,) uniform or (nElem, 2)")
    # CoefficientFunction: evaluate per centroid
    out = np.zeros((len(centroids), 2))
    for i, (px, py) in enumerate(centroids):
        val = H_ext(mesh(px, py))
        out[i] = (val[0], val[1])
    return out


def _M_avg(M, areas):
    w = areas / areas.sum()
    return np.array([float(w @ M[:, 0]), float(w @ M[:, 1])])


def _law_from_table(bh_table):
    """[[H,B],...] (A/m, T) -> (M_of_h, chi_sec, chi0) with saturation clamp beyond Hmax."""
    tab = np.asarray(bh_table, float)
    if tab.ndim != 2 or tab.shape[1] != 2 or tab.shape[0] < 3:
        raise ValueError("bh_table must be [[H, B], ...] (A/m, T) with >= 3 rows")
    H, Bt = tab[:, 0].copy(), tab[:, 1].copy()
    if H[0] != 0.0:
        H = np.concatenate([[0.0], H]); Bt = np.concatenate([[0.0], Bt])
    if np.any(np.diff(H) <= 0):
        raise ValueError("bh_table H column must be strictly increasing")
    M = Bt / MU0 - H
    if np.any(M < -1e-9):
        raise ValueError("bh_table implies negative magnetization (B < mu0 H) -- not a soft iron")
    chi0 = M[1] / H[1]

    def M_of_h(h):
        return np.interp(h, H, M)

    def chi_sec(h):
        h = np.asarray(h, float)
        out = np.full_like(h, chi0)
        big = h >= H[1]
        out[big] = np.interp(h[big], H, M) / h[big]
        return out
    return M_of_h, chi_sec, chi0


def solve_planar_demag(mesh, mu_r=None, H_ext=None, bh_table=None, *,
                       nl_tol=1e-6, nl_maxit=300, nl_damp=0.6, chi_floor=1e-12):
    """Single-region planar soft-iron demag solve (the C++ 2D moment core + this driver).

    EXACTLY ONE of ``mu_r`` (linear) or ``bh_table`` (nonlinear [[H,B],...]) must be given.
    ``H_ext`` is a 2-tuple (uniform) or a 2-component NGSolve CoefficientFunction.

    Returns dict: M (nElem,2), M_avg (2,), demag_factors (Dx,Dy), iters, residual, ndof (edge DOF),
    n_el, nonlinear (bool), linear_solver='dense-2d-cpp'.
    """
    if H_ext is None:
        raise ValueError("solve_planar_demag: H_ext is required")
    if (mu_r is None) == (bh_table is None):
        raise ValueError("solve_planar_demag: provide EXACTLY ONE of mu_r (linear) or bh_table")
    if isinstance(mu_r, dict) or isinstance(bh_table, dict):
        raise NotImplementedError("solve_planar_demag: per-region (dict) materials not wired yet "
                                  "(the first increment is a single soft-iron region)")
    verts, offsets, centroids, areas = _extract_geometry(mesh)
    nElem = len(areas)
    ndof = int(offsets[-1])
    Hc = _eval_Hext(H_ext, centroids, mesh)

    if mu_r is not None:
        if not mu_r > 1.0:
            raise ValueError("solve_planar_demag: mu_r must be > 1 (got %r)" % (mu_r,))
        chi = np.full(nElem, mu_r - 1.0)
        M = _rp.Moment2DSolveLinear(verts, offsets, chi, Hc)
        iters, res, nonlinear = 1, 0.0, False
    else:
        M_of_h, chi_sec, chi0 = _law_from_table(bh_table)
        # scalar-chi Picard: linear solve with per-element chi, then H(c)=M/chi, update chi (secant)
        chi = np.full(nElem, chi0)
        prev = None
        res = np.inf
        M = None
        for it in range(nl_maxit):
            M = _rp.Moment2DSolveLinear(verts, offsets, np.maximum(chi, chi_floor), Hc)
            He = M / np.maximum(chi, chi_floor)[:, None]          # H(c) = M/chi (dipole row is exact)
            nH = np.maximum(np.linalg.norm(He, axis=1), 1e-300)
            chi_star = np.maximum(M_of_h(nH) / nH, chi_floor)
            r = chi_star - chi
            res = np.linalg.norm(r) / max(np.linalg.norm(chi_star), 1e-300)
            if res < nl_tol:
                iters = it + 1
                break
            chi_next = chi + nl_damp * r
            if prev is not None:                                 # safeguarded Anderson(1)
                chi_p, r_p = prev
                dr = r - r_p
                den = float(dr @ dr)
                if den > 1e-300:
                    th = float(r @ dr) / den
                    cand = (1 - th) * (chi + r) + th * (chi_p + r_p)
                    if np.all(cand > 0):
                        chi_next = np.maximum(cand, chi_floor)
            prev = (chi.copy(), r.copy())
            chi = np.maximum(chi_next, chi_floor)
        else:
            raise RuntimeError("solve_planar_demag: Picard NOT converged (res=%.2e after %d)"
                               % (res, nl_maxit))
        nonlinear = True

    Mavg = _M_avg(M, areas)
    return {
        "M": M, "M_avg": Mavg, "demag_factors": _demag_factors(verts, offsets, areas, mu_r),
        "iters": iters, "residual": float(res), "ndof": ndof, "n_el": nElem,
        "nonlinear": nonlinear, "linear_solver": "dense-2d-cpp",
    }


def _demag_factors(verts, offsets, areas, mu_r):
    """(Dx, Dy) from two unit-field LINEAR solves.  Uses mu_r if given, else a default chi=1
    probe (the demag factor is chi-dependent for non-ellipsoidal bodies, so it is reported at the
    solve's own chi when linear)."""
    chi_val = (mu_r - 1.0) if (mu_r is not None and mu_r > 1.0) else 1.0
    nElem = len(areas)
    chi = np.full(nElem, chi_val)
    out = []
    for axis in range(2):
        H = np.zeros((nElem, 2)); H[:, axis] = 1.0
        M = _rp.Moment2DSolveLinear(verts, offsets, chi, H)
        Ma = _M_avg(M, areas)[axis]
        out.append(1.0 / Ma - 1.0 / chi_val)
    return (out[0], out[1])


def demag_factors(mesh, chi=1.0):
    """Convenience: (Dx, Dy) demag factors of a planar body at susceptibility ``chi`` (two unit
    LINEAR solves).  For an ellipse a(x):b(y) the exact values are b/(a+b), a/(a+b)."""
    verts, offsets, _, areas = _extract_geometry(mesh)
    return _demag_factors(verts, offsets, areas, chi + 1.0)
