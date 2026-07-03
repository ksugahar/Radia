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
    H, M, chi0 = _hm_arrays(bh_table)

    def M_of_h(h):
        return np.interp(h, H, M)

    def chi_sec(h):
        h = np.asarray(h, float)
        out = np.full_like(h, chi0)
        big = h >= H[1]
        out[big] = np.interp(h[big], H, M) / h[big]
        return out
    return M_of_h, chi_sec, chi0


def _hm_arrays(bh_table):
    """[[H,B],...] (A/m,T) -> (H, M, chi0), 0-anchored + validated (shared by single + per-region)."""
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
    return H, M, M[1] / H[1]


def _element_materials(mesh):
    """Per-element material name (same VOL element order as _extract_geometry)."""
    return [el.mat for el in mesh.Elements(ng.VOL)]


def _region_ids(mats):
    d = {}
    for i, m in enumerate(mats):
        d.setdefault(m, []).append(i)
    return {k: np.asarray(v, int) for k, v in d.items()}


def _check_regions(mats, provided, what):
    missing = set(mats) - set(provided)
    if missing:
        raise ValueError("solve_planar_demag: mesh regions %s have no %s; provided: %s"
                         % (sorted(missing), what, sorted(provided)))


def _per_region_chi(mesh, mu_r_dict):
    """Per-element chi from a {region_name: mu_r} dict (multi-grade soft iron)."""
    mats = _element_materials(mesh)
    _check_regions(mats, mu_r_dict, "mu_r")
    chi = np.empty(len(mats))
    for name, ids in _region_ids(mats).items():
        mr = mu_r_dict[name]
        if not mr > 1.0:
            raise ValueError("solve_planar_demag: mu_r[%r] must be > 1 (got %r)" % (name, mr))
        chi[ids] = mr - 1.0
    return chi


def _per_region_law(mesh, bh_dict):
    """Per-element (M_of_h, chi_sec, chi0) from a {region_name: [[H,B],...]} dict."""
    mats = _element_materials(mesh)
    _check_regions(mats, bh_dict, "bh_table")
    rid = _region_ids(mats)
    law = {name: _hm_arrays(bh_dict[name]) for name in rid}
    chi0_e = np.empty(len(mats))
    for name, ids in rid.items():
        chi0_e[ids] = law[name][2]

    def M_of_h(h):
        out = np.empty(len(mats))
        for name, ids in rid.items():
            H, M, _ = law[name]
            out[ids] = np.interp(h[ids], H, M)
        return out

    def chi_sec(h):
        out = np.empty(len(mats))
        for name, ids in rid.items():
            H, M, c0 = law[name]
            hi = h[ids]; ci = np.full_like(hi, c0)
            big = hi >= H[1]
            ci[big] = np.interp(hi[big], H, M) / hi[big]
            out[ids] = ci
        return out
    return M_of_h, chi_sec, chi0_e


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
    per_region = isinstance(mu_r, dict) or isinstance(bh_table, dict)
    verts, offsets, centroids, areas = _extract_geometry(mesh)
    nElem = len(areas)
    ndof = int(offsets[-1])
    Hc = _eval_Hext(H_ext, centroids, mesh)

    if mu_r is not None:
        # LINEAR: uniform (scalar mu_r) or per-region ({region: mu_r} dict, multi-grade soft iron)
        if isinstance(mu_r, dict):
            chi = _per_region_chi(mesh, mu_r)
        else:
            if not mu_r > 1.0:
                raise ValueError("solve_planar_demag: mu_r must be > 1 (got %r)" % (mu_r,))
            chi = np.full(nElem, mu_r - 1.0)
        M = _rp.Moment2DSolveLinear(verts, offsets, chi, Hc)
        iters, res, nonlinear = 1, 0.0, False
    else:
        # NONLINEAR: single table or per-region ({region: [[H,B],...]} dict)
        if isinstance(bh_table, dict):
            M_of_h, chi_sec, chi0 = _per_region_law(mesh, bh_table)
            chi = np.array(chi0, float)                          # per-element initial chi
        else:
            M_of_h, chi_sec, chi0 = _law_from_table(bh_table)
            chi = np.full(nElem, chi0)
        # scalar-chi Picard: linear solve with per-element chi, then H(c)=M/chi, update chi (secant)
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
    # demag factors are a single-body concept -> None for a heterogeneous (per-region) body
    dfac = None if per_region else _demag_factors(verts, offsets, areas, mu_r)
    return {
        "M": M, "M_avg": Mavg, "demag_factors": dfac,
        "iters": iters, "residual": float(res), "ndof": ndof, "n_el": nElem,
        "nonlinear": nonlinear, "linear_solver": "dense-2d-cpp", "per_region": per_region,
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


# ---- exterior field / Maxwell torque -- delegate to the SHARED planar_charges layer (the same
#      routine the HDiv-VIM uses): given the solved per-element M, the exterior field is the field
#      of that magnetization (M.n equivalent charges, C++ log kernel). --------------------------

def exterior_field(mesh, M_elem, P, ngauss=4):
    """Exterior H (in air) at points P (n,2) of the solved magnetization M_elem on ``mesh``."""
    from radia.planar_charges import exterior_field as _ef
    return _ef(mesh, M_elem, P, ngauss=ngauss)


def maxwell_torque(mesh, M_elem, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440, ngauss=4):
    """Maxwell-stress torque per unit length on a circle Rc (in air, enclosing the body) for the
    solved magnetization M_elem in a UNIFORM applied field H_ext.  Reluctance torque = mu0 A (M x H0)."""
    from radia.planar_charges import maxwell_torque as _mt
    return _mt(mesh, M_elem, Rc, H_ext=H_ext, center=center, n=n, ngauss=ngauss)


def torque_angle_sweep(mesh, H0, angles_rad, Rc, *, mu_r=None, bh_table=None, center=(0.0, 0.0),
                       n=1440, ngauss=4, **solve_kw):
    """Reluctance-torque-vs-angle sweep -- the planar-motor headline.  For each applied-field angle
    ``theta`` the body is re-solved with H_ext = H0 (cos theta, sin theta) and the Maxwell-stress
    torque about ``center`` on the circle Rc is returned.

    Returns dict: angles (rad), torque (nAng,), M_avg (nAng,2).  (LINEAR: mu_r; NONLINEAR: bh_table.)
    For a LINEAR material the moment matrix is angle-INDEPENDENT (only the RHS rotates), so it is
    assembled + LU-factored ONCE and back-substituted for all angles (C++ Moment2DSolveMulti).
    NONLINEAR (bh_table) re-solves per angle (chi varies -> the matrix changes)."""
    angles = np.asarray(angles_rad, float)
    T = np.zeros(len(angles))
    Mavg = np.zeros((len(angles), 2))
    if mu_r is not None and bh_table is None:
        # LINEAR: factor ONCE, solve all angles
        verts, offsets, _, areas = _extract_geometry(mesh)
        nEl = len(areas)
        if isinstance(mu_r, dict):
            chi = _per_region_chi(mesh, mu_r)
        else:
            if not mu_r > 1.0:
                raise ValueError("torque_angle_sweep: mu_r must be > 1 (got %r)" % (mu_r,))
            chi = np.full(nEl, mu_r - 1.0)
        Hmulti = np.zeros((len(angles), nEl, 2))
        Hmulti[:, :, 0] = (H0 * np.cos(angles))[:, None]
        Hmulti[:, :, 1] = (H0 * np.sin(angles))[:, None]
        Mmulti = _rp.Moment2DSolveMulti(verts, offsets, chi, Hmulti)      # (nAng, nEl, 2)
        w = areas / areas.sum()
        for i, th in enumerate(angles):
            H_ext = (H0 * np.cos(th), H0 * np.sin(th))
            T[i] = maxwell_torque(mesh, Mmulti[i], Rc, H_ext=H_ext, center=center, n=n, ngauss=ngauss)
            Mavg[i] = (float(w @ Mmulti[i, :, 0]), float(w @ Mmulti[i, :, 1]))
        return {"angles": angles, "torque": T, "M_avg": Mavg,
                "factored_once": True, "per_region": isinstance(mu_r, dict)}
    for i, th in enumerate(angles):                                       # NONLINEAR: per angle
        H_ext = (H0 * np.cos(th), H0 * np.sin(th))
        r = solve_planar_demag(mesh, mu_r=mu_r, bh_table=bh_table, H_ext=H_ext, **solve_kw)
        T[i] = maxwell_torque(mesh, r["M"], Rc, H_ext=H_ext, center=center, n=n, ngauss=ngauss)
        Mavg[i] = r["M_avg"]
    return {"angles": angles, "torque": T, "M_avg": Mavg,
            "factored_once": False, "per_region": isinstance(bh_table, dict)}
