"""radia.planar_aniso -- SHARED 2D anisotropic-susceptibility linear demag solver.

Grain-oriented (GO) silicon steel and other uniaxially-anisotropic laminations obey M = X.H with X a
2x2 susceptibility tensor (``planar_materials.chi_tensor``: chi_par along the easy axis, chi_perp
across it).  Neither the scalar collocation MMMM (``Moment2DSolveLinear``) nor the HDiv-VIM handles a
tensor chi, and a matrix-free Picard on M is ILL-CONDITIONED for a real chi (the demag operator X.N
has spectral radius ~ chi_max/2, so the Picard diverges for chi>>2).  So this assembles the demag
operator N DENSELY on the SHARED planar_charges kernel (N[i,j] = field at centroid i from a unit
magnetisation on element j, the Gauss-sampled M.n log-charge cloud) and solves the well-conditioned
dense system directly:

    (I - X_blockdiag N) M = X_blockdiag H0 .

Because N is built from ``planar_charges`` (the same kernel MMMM and the HDiv-VIM already use for field
evaluation), this is METHOD-AGNOSTIC -- one anisotropic solver for both.  Verified: the isotropic
special case (X = chi I) reproduces the exact ``Moment2DSolveLinear`` to ~1e-4 (Gauss-N self-term), and
the anisotropic disk matches the analytic  M = (I + D X)^-1 X H0  (D = 1/2) to ~3e-4.

2D moment methods are FEW-element (motor cross-section ~1e2-1e3 DOF), so the dense O(n^3) solve is
cheap; for very large n use the scalar solvers (this is the tensor-chi niche).

    import radia.planar_aniso as pa
    r = pa.solve_anisotropic_demag(mesh, chi_par=5000, chi_perp=200, easy_deg=0.0, H0=(H0, 0))
    r["M"]       # per-element magnetisation (n,2)     r["M_avg"]   # volume average (2,)
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng

from radia.planar_charges import charge_field, _gauss01
from radia.planar_materials import chi_tensor, region_ids, check_regions
from radia.mmmm2d import _extract_geometry, _element_materials

MU0 = 4e-7 * np.pi


def _element_clouds(mesh, ngauss):
    """Per element: (Xq (m,2), NWL (m,2)) -- the edge Gauss points and, at each, the OUTWARD normal
    times the edge length * Gauss weight, so the M.n charge of a unit e_d magnetisation is NWL[:,d]."""
    if mesh.dim != 2:
        raise ValueError("planar_aniso: mesh.dim must be 2 (got %d)" % mesh.dim)
    pts = np.array([list(mesh[v].point)[:2] for v in mesh.vertices])
    tg, wg = _gauss01(ngauss)
    clouds = []
    for el in mesh.Elements(ng.VOL):
        V = pts[[v.nr for v in el.vertices]]
        x, y = V[:, 0], V[:, 1]
        A2 = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        if A2 < 0:
            V = V[::-1]
        c = V.mean(axis=0)
        nv = len(V)
        Xs, NW = [], []
        for i in range(nv):
            P0 = V[i]; P1 = V[(i + 1) % nv]
            t = P1 - P0; L = np.hypot(t[0], t[1]); th = t / L
            nout = np.array([th[1], -th[0]])
            mid = 0.5 * (P0 + P1)
            if nout @ (mid - c) < 0:
                nout = -nout
            X = P0[None, :] + tg[:, None] * (P1 - P0)[None, :]            # (ngauss, 2)
            Xs.append(X)
            NW.append(np.tile(nout, (ngauss, 1)) * (L * wg)[:, None])     # (ngauss, 2): n_out * L * w
        clouds.append((np.vstack(Xs), np.vstack(NW)))
    return clouds


def demag_operator(mesh, centroids=None, ngauss=6):
    """Dense demag operator N (2n x 2n): (N @ M_flat)[i] = field at centroid i from magnetisation M
    (M_flat = [Mx0,My0,Mx1,My1,...]).  Built column-by-column via the shared charge_field kernel."""
    if centroids is None:
        _, _, centroids, _ = _extract_geometry(mesh)
    clouds = _element_clouds(mesh, ngauss)
    n = len(clouds)
    N = np.zeros((2 * n, 2 * n))
    for j, (Xq, NWL) in enumerate(clouds):
        for d in range(2):
            H = charge_field(Xq, np.ascontiguousarray(NWL[:, d]), centroids)   # (n,2) field at centroids
            N[0::2, 2 * j + d] = H[:, 0]
            N[1::2, 2 * j + d] = H[:, 1]
    return N


def _X_per_element(mesh, chi_par, chi_perp, easy_deg):
    """(n,2,2) per-element susceptibility tensor.  Scalars -> uniform; dicts -> per material region."""
    mats = _element_materials(mesh)
    n = len(mats)
    if isinstance(chi_par, dict) or isinstance(chi_perp, dict) or isinstance(easy_deg, dict):
        cp = chi_par if isinstance(chi_par, dict) else {m: chi_par for m in set(mats)}
        cq = chi_perp if isinstance(chi_perp, dict) else {m: chi_perp for m in set(mats)}
        ea = easy_deg if isinstance(easy_deg, dict) else {m: easy_deg for m in set(mats)}
        check_regions(mats, cp, "chi_par"); check_regions(mats, cq, "chi_perp")
        check_regions(mats, ea, "easy_deg")
        X = np.empty((n, 2, 2))
        for name, ids in region_ids(mats).items():
            X[ids] = chi_tensor(cp[name], cq[name], ea[name])
        return X
    return np.tile(chi_tensor(chi_par, chi_perp, easy_deg), (n, 1, 1))


def solve_anisotropic_demag(mesh, chi_par, chi_perp, easy_deg=0.0, H0=(0.0, 0.0), ngauss=6):
    """Anisotropic linear soft-iron demag (M = X.H), X = uniaxial tensor (chi_par easy / chi_perp
    across, easy axis ``easy_deg`` from +x).  Scalars = uniform; {region: value} dicts = per grade.
    ``H0`` is the uniform applied field (2,).  Returns dict: M (n,2), M_avg (2,), n_el, ndof(=2n)."""
    _, _, centroids, areas = _extract_geometry(mesh)
    n = len(areas)
    X = _X_per_element(mesh, chi_par, chi_perp, easy_deg)
    N = demag_operator(mesh, centroids, ngauss)
    Xbd = np.zeros((2 * n, 2 * n))
    for k in range(n):
        Xbd[2 * k:2 * k + 2, 2 * k:2 * k + 2] = X[k]
    H0f = np.tile(np.asarray(H0, float), n)
    M = np.linalg.solve(np.eye(2 * n) - Xbd @ N, Xbd @ H0f).reshape(n, 2)
    w = areas / areas.sum()
    return {"M": M, "M_avg": np.array([w @ M[:, 0], w @ M[:, 1]]), "n_el": n, "ndof": 2 * n,
            "linear_solver": "dense-aniso-2d"}
