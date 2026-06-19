"""hdiv_demag_tet.py -- HDiv-type VIM demag on UNSTRUCTURED TET meshes (production step #1a).

The structured-hex prototype (hdiv_demag_quad_self.py) is generalized to real tet RT0 meshes from
NGSolve.  The NGSolve HDiv(order=0) extraction is element-AGNOSTIC:
  B = [ -div(u) tested on L2 (volume charge rho) ;  u.Trace().n tested on SurfaceL2 (surface charge sigma) ]
  M_mass = the HDiv(0) mass (the physical demag factors are eig(N, M_mass), basis-invariant).
Only the Gram self-energy geometry is element-specific -> here TET (barycentric sub-points) and
TRIANGLE faces, with the self constants c_tet / c_tri (the tet/tri analogs of c_cube / c_sq).

Verification (the physics gate): a tet-meshed SPHERE has demag factor EXACTLY 1/3 (isotropic);
a CUBE has the three axis demag factors equal to 1/3 each (sum to 1 by symmetry).  We compute the
demag factor of a UNIFORM M_z as the Rayleigh quotient (m^T N m)/(m^T M_mass m) with m = the RT0
L2-projection of (0,0,1).  Sphere -> ~1/3 confirms the unstructured tet operator is physical.

This is the Python reference; production step #1b moves the dense Gram G to the C++ HACApK charge
H-matrix (scalable) and checks it reproduces this.
"""
import json
import os
from math import pi, sqrt, log, atan2

import numpy as np
import scipy.sparse as sp
from scipy.spatial import ConvexHull

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt, OrthoBrick

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 1e-6

# ---- self constants c = INT INT_{unit-measure shape} 1/|x-y| (sub-cell self correction) ----
# cube/square via exact quadrature (as the hex prototype); tet/tri extrapolated (compute_self_constants.py)
from scipy import integrate  # noqa: E402
_c_sq, _ = integrate.dblquad(lambda v, u: (1 - u) * (1 - v) / sqrt(u * u + v * v + 1e-300), 0, 1, 0, 1)
C_SQ = 4.0 * _c_sq
C_TET = 1.77   # INT INT 1/|x-y| over a unit-volume regular tet (Richardson, ~0.3% uncertainty)
C_TRI = 2.89   # INT INT 1/|x-y| over a unit-area equilateral triangle


def _bary_tet(nsub):
    """Barycentric sub-point lattice of a reference tet: (i,j,k,l) lattice nodes with i+j+k+l=nsub-1,
    point (i+1/4, j+1/4, k+1/4, l+1/4)/nsub.  The 4 coords sum to EXACTLY 1 (inside the tet) and
    nsub=1 -> the centroid (1/4,1/4,1/4,1/4).  Equal sub-weights.  (The earlier [.+0.5,...]/nsub form
    summed to 1+1/nsub -> sub-points OUTSIDE the tet -- a quadrature bug; fixed 2026-06-07.)"""
    lam = []
    for i in range(nsub):
        for j in range(nsub - i):
            for k in range(nsub - i - j):
                l = nsub - 1 - i - j - k
                lam.append([(i + 0.25) / nsub, (j + 0.25) / nsub, (k + 0.25) / nsub, (l + 0.25) / nsub])
    return np.array(lam)   # (m, 4) barycentric, each row sums to 1


def _bary_tri(nsub):
    """Barycentric sub-point lattice of a reference triangle: (i,j,k) with i+j+k=nsub-1, point
    (i+1/3, j+1/3, k+1/3)/nsub (sums to EXACTLY 1; nsub=1 -> centroid).  Equal sub-weights."""
    lam = []
    for i in range(nsub):
        for j in range(nsub - i):
            k = nsub - 1 - i - j
            lam.append([(i + 1.0 / 3) / nsub, (j + 1.0 / 3) / nsub, (k + 1.0 / 3) / nsub])
    return np.array(lam)   # (m, 3), each row sums to 1


def _gauss_duffy_tet(o):
    """Gauss-Duffy collapsed-cube tet rule: o Gauss-Legendre pts/dim (o^3 nodes), ref-tet coords
    (lam1,lam2,lam3) + weights summing to 1/6.  For the SMOOTH outer integral of an EXACT analytic inner
    (phi_tet/tri_potential) this converges to ~machine precision -- it REPLACES the crude equal-weight
    _bary_tet outer rule, which under-integrated the volume self-energy by ~6.5% (golden-invisible because
    every uniform-M demag has div M = 0).  Same rule as radia.vim._vim._tet_ref(o)."""
    x, w = np.polynomial.legendre.leggauss(o)
    s, ws = 0.5 * (x + 1.0), 0.5 * w
    P, W = [], []
    for a, wa in zip(s, ws):
        for b, wb in zip(s, ws):
            for c, wc in zip(s, ws):
                P.append((a, b * (1 - a), c * (1 - a) * (1 - b)))
                W.append(wa * wb * wc * (1 - a) ** 2 * (1 - b))
    return np.array(P), np.array(W)   # (m,3) ref pts (lam1,lam2,lam3), weights sum 1/6


def tet_self_energy(V, vol, nsub):
    """G_aa for a tet (4 verts V) of volume `vol`: cross sub-point sum + sub-cell self (c_tet).

    NOTE (2026-06-12): this empirical cross-subpoint + C_TET*w^(5/3) self is ~10-15% LOW -- the same
    volume-self quadrature defect that was fixed in the C++ analytic ctor (the production order-0 Gram) and
    in analytic_charge_gram (the dense reference).  It is NOT fixed here because it feeds only the dense-path
    PRECONDITIONER diagonal + the monopole-path Gram diagonal, never the production scalable operator (which
    uses the C++ ctor's own QuadDot(a,a)).  Correcting it makes the wilton_surface=True non-uniform solve
    CONVERGE (to ~0.6% of analytic_gram) instead of diverging, which retires the surface-only fail-loud
    guard (test_hdiv_vim_cyoke_nonlinear) -- a behaviour change deferred for a separate decision."""
    lam = _bary_tet(nsub)
    C = lam @ V                                   # (m,3) sub-point positions
    w = np.full(len(C), vol / len(C))             # equal sub-weights, sum = vol
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2)
    np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(w, w) / (4 * pi * D))
    selfsub = np.sum(C_TET * w ** (5.0 / 3.0) / (4 * pi))
    return cross + selfsub


def tri_self_energy(V, area, nsub):
    lam = _bary_tri(nsub)
    C = lam @ V
    w = np.full(len(C), area / len(C))
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2)
    np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(w, w) / (4 * pi * D))
    selfsub = np.sum(C_TRI * w ** 1.5 / (4 * pi))
    return cross + selfsub


# --- Wilton analytic surface Gram (the #1 production accuracy fix, 2026-06-07) ---
# For a UNIFORM M the demag is ENTIRELY surface charge (div M = 0 -> zero volume charge), so the demag
# factor is governed by the surface-surface (boundary-triangle) Gram block.  The centroid-MONOPOLE
# off-diagonal under-resolves adjacent/near boundary triangles -> the demag factor comes out ~5-6% low
# (cube 0.311, sphere 0.314 vs the exact 1/3), and the sub-point near-correction does NOT fix the cube.
# Wilton's exact analytic potential of a uniformly-charged flat triangle closes this: cube/sphere demag
# -> 1/3 to <0.15%.  Reference: Wilton et al., IEEE TAP 32(3):276 (1984); Graglia, IEEE TAP 41(10):1448
# (1993).  Off-diagonal only -- the validated tri_self_energy keeps the diagonal.

# Dunavant degree-5 symmetric triangle rule (7 points, barycentric; weights sum to 1) for the OUTER
# integral over the observation triangle (the inner integral is the exact Wilton potential).
_DUN5 = [(1.0 / 3, 1.0 / 3, 1.0 / 3, 0.225)]
for _a, _w in [(0.0597158717, 0.1323941527), (0.7974269853, 0.1259391805)]:
    _b = (1.0 - _a) / 2.0
    _DUN5 += [(_a, _b, _b, _w), (_b, _a, _b, _w), (_b, _b, _a, _w)]
_DUN5 = np.array(_DUN5)


def tri_potential(V, r):
    """Exact INT_T 1/|r - r'| dA' over a flat triangle T (vertices V, 3x3) at observation point(s) r
    (constant unit density) -- the Wilton/Graglia analytic potential integral.  Verified vs fine
    numerical quadrature to the reference's own accuracy (~1e-4, limited by the reference near the
    plane).  Handles r on either side of / on the triangle plane.

    Vectorized: r may be a single point (3,) -> float, or a batch (M,3) -> (M,) -- the batch form is
    what makes the O(n_bf^2) surface Gram practical."""
    r = np.asarray(r, float)
    single = (r.ndim == 1)
    R = r.reshape(1, 3) if single else r            # (M,3)
    v0, v1, v2 = V[0], V[1], V[2]
    n = np.cross(v1 - v0, v2 - v0)
    n = n / np.linalg.norm(n)
    d = (R - v0) @ n                                # (M,) signed heights above the plane
    p = R - d[:, None] * n                          # (M,3) projections onto the plane
    ad = np.abs(d)
    I = np.zeros(len(R))
    vs = (v0, v1, v2)
    for i in range(3):
        a = vs[i]
        b = vs[(i + 1) % 3]
        lh = (b - a) / np.linalg.norm(b - a)
        uh = np.cross(lh, n)                         # in-plane unit normal to the edge
        P0 = (a - p) @ uh                            # (M,) signed perp distance projection -> edge
        sm = (a - p) @ lh                            # (M,) endpoint projections along the edge
        sp = (b - p) @ lh
        Rm = np.linalg.norm(R - a, axis=1)           # (M,) = |r - a|
        Rp = np.linalg.norm(R - b, axis=1)
        R0sq = P0 * P0 + d * d
        dm, dp = Rm + sm, Rp + sp
        safe = (dp > 1e-300) & (dm > 1e-300)
        ratio = np.where(safe, dp, 1.0) / np.where(safe, dm, 1.0)
        f = np.where(safe, np.log(ratio), 0.0)
        beta = np.arctan2(P0 * sp, R0sq + ad * Rp) - np.arctan2(P0 * sm, R0sq + ad * Rm)
        I += P0 * f - ad * beta
    return float(I[0]) if single else I


def wilton_surface_block(bf_V):
    """Surface-surface Gram block (n_bf x n_bf) by the Wilton analytic inner integral + a Dunavant
    outer rule: G[a][b] = (1/4pi) INT_{tri_a} (INT_{tri_b} 1/|r-r'| dA') dA.  Symmetric; the diagonal
    here is the Dunavant-outer self (shape-exact for any triangle).  Vectorized over the observation
    points: for each SOURCE triangle b, evaluate its Wilton potential at ALL outer quad points at once
    -> the Python loop is O(n_bf), not O(n_bf^2)."""
    nb = len(bf_V)
    area = np.array([0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0])) for V in bf_V])
    QP = np.array([_DUN5[:, :3] @ V for V in bf_V])  # (nb, 7, 3) outer quad points per triangle
    allP = QP.reshape(-1, 3)                          # (nb*7, 3)
    wq = _DUN5[:, 3]                                  # (7,)
    G = np.zeros((nb, nb))
    for b in range(nb):
        phi = tri_potential(bf_V[b], allP).reshape(nb, 7)   # potential of source b at every outer point
        G[:, b] = (phi * wq[None, :]).sum(axis=1) * area / (4 * pi)
    return 0.5 * (G + G.T)


def phi_tet(V, P):
    """Newtonian potential INT_tet 1/|P - r'| dV' of a uniform tetrahedron (vertices V, 4x3) at P --
    the VOLUME analog of tri_potential, via the divergence theorem (nabla'^2 R = 2/R):
        INT_V 1/R dV = (1/2) sum_{4 faces} d_face * INT_face 1/R dA'
    with d_face = (r' - P).n_hat (constant on a flat face) and the inner face integral = tri_potential.
    Reuses the exact Wilton triangle potential -> the volume Gram needs NO new singular quadrature.
    P: (3,) -> float, or (M,3) -> (M,).  Verified vs fine volume quadrature to the reference's accuracy."""
    P = np.asarray(P, float)
    single = (P.ndim == 1)
    R = P.reshape(1, 3) if single else P
    cen = V.mean(0)
    tot = np.zeros(len(R))
    for f in ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)):
        Fv = V[list(f)]
        n = np.cross(Fv[1] - Fv[0], Fv[2] - Fv[0])
        n = n / np.linalg.norm(n)
        if np.dot(Fv.mean(0) - cen, n) < 0:
            n = -n                                          # outward face normal
        d = (Fv[0] - R) @ n                                 # (M,) signed distance P -> face plane
        tot += d * tri_potential(Fv, R)
    tot *= 0.5
    return float(tot[0]) if single else tot


# --- polytope (hex / wedge cell, quad face) charge-Gram geometry ---------------------------------
# The tet/triangle analytic charge Gram (phi_tet / tri_potential + Gauss-Duffy / Dunavant) generalizes
# to ANY flat-faced convex cell + quad face with NO new singular quadrature: the cell Newtonian
# potential is the divergence-theorem sum over the cell's (convex-hull) triangular faces of the SAME
# exact Wilton triangle potential (tri_potential), and a quad face is two flat triangles.  This is the
# hex/wedge path for build_demag(analytic_gram=True); tets + triangles keep the EXACT phi_tet /
# tri_potential calls (bit-identical -- the C++ analytic charge-Gram golden depends on it).


def _cell_hull_tris(V):
    """(triangle (3,3), outward unit normal) faces of a convex cell (hex/wedge/...) via ConvexHull.
    Feeds the polytope volume-charge Newtonian potential; tets use phi_tet directly (this is only
    reached for cells with != 4 vertices)."""
    V = np.asarray(V, float)
    hull = ConvexHull(V)
    cen = V.mean(0)
    tris = []
    for simplex in hull.simplices:
        P = V[simplex]
        nrm = np.cross(P[1] - P[0], P[2] - P[0])
        ln = np.linalg.norm(nrm)
        if ln < 1e-300:
            continue
        nrm = nrm / ln
        if np.dot(nrm, P.mean(0) - cen) < 0:
            nrm = -nrm                                       # outward
        tris.append((P, nrm))
    return tris


def _polytope_potential(tris, R):
    """VECTORIZED INT_cell 1/|R-r'| dV' over a flat-faced convex cell (hull faces `tris`) at points R,
    via the divergence theorem on the (vectorized) Wilton triangle potential -- the hex/wedge analog of
    phi_tet (the two agree on a tet to the reference's accuracy, verified).  R: (3,)->float, (M,3)->(M,)."""
    R = np.asarray(R, float)
    single = (R.ndim == 1)
    RR = R.reshape(1, 3) if single else R
    tot = np.zeros(len(RR))
    for (P, nrm) in tris:
        d = (P[0] - RR) @ nrm                                # (M,) signed distance R -> face plane
        tot += d * tri_potential(P, RR)
    tot *= 0.5
    return float(tot[0]) if single else tot


def _face_subtris(V):
    """Sub-triangles of a boundary face: triangle (3 verts) -> [V]; quad (4 verts) -> 2-triangle fan."""
    V = np.asarray(V, float)
    if len(V) == 3:
        return [V]
    if len(V) == 4:
        return [V[[0, 1, 2]], V[[0, 2, 3]]]
    raise ValueError("HDiv-VIM: boundary face with %d vertices is not supported (triangle or quad only)"
                     % len(V))


def _cell_outer_quad(V, vol, lam_t, w_t):
    """Interior quadrature (points (m,3), weights (m,) summing to `vol`) over a convex cell: centroid-fan
    tets (centroid + each hull face), each filled by the Gauss-Duffy tet rule (lam_t, w_t)."""
    V = np.asarray(V, float)
    cen = V.mean(0)
    Ps, Ws = [], []
    for (T, _n) in _cell_hull_tris(V):
        tet = np.vstack([cen, T])
        tvol = abs(np.linalg.det(tet[1:] - tet[0])) / 6.0
        Ps.append(tet[0] + lam_t @ (tet[1:] - tet[0]))
        Ws.append(w_t * 6.0 * tvol)
    return np.vstack(Ps), np.concatenate(Ws)


def _face_outer_quad(V, area):
    """Outer quadrature (points, weights summing to `area`) over a boundary face: Dunavant-5 per sub-tri."""
    Ps, Ws = [], []
    for T in _face_subtris(V):
        a = 0.5 * np.linalg.norm(np.cross(T[1] - T[0], T[2] - T[0]))
        Ps.append(_DUN5[:, :3] @ T)
        Ws.append(_DUN5[:, 3] * a)
    return np.vstack(Ps), np.concatenate(Ws)


def _poly_cell_self_energy(V, vol):
    """Analytic volume self-energy G_aa = (1/4pi) INT_cell INT_cell 1/r for a convex polytope cell (its
    own outer quadrature against its own Newtonian potential) -- the hex/wedge analog of tet_self_energy,
    but EXACT-quadrature rather than the C_TET extrapolation.  Used only for the monopole-path
    preconditioner diagonal (the analytic Gram computes its own diagonal)."""
    lam_t, w_t = _gauss_duffy_tet(4)
    P, W = _cell_outer_quad(V, vol, lam_t, w_t)
    return float(np.sum(W * _polytope_potential(_cell_hull_tris(V), P)) / (4.0 * pi))


def _poly_face_self_energy(V, area):
    """Analytic surface self-energy G_aa for a quad/triangle boundary face (Dunavant outer against its
    sub-triangle Wilton potentials).  Monopole-path preconditioner diagonal only."""
    P, W = _face_outer_quad(V, area)
    pot = sum(tri_potential(T, P) for T in _face_subtris(V))
    return float(np.sum(W * pot) / (4.0 * pi))


def analytic_charge_gram(el_V, bf_V, el_vol, bf_area, nsub_tet=3):
    """FULL analytic charge Gram (n_charge x n_charge, cells then boundary faces): every block via the
    EXACT analytic potential -- cell sources via phi_tet (tet) / _polytope_potential (hex, wedge, ...),
    face sources via tri_potential (triangle) / sum-of-tri (quad) -- integrated against an outer
    quadrature on the target cell (tet sub-points / centroid-fan) / face (Dunavant).  Replaces the
    centroid-monopole G everywhere, so the vol-vol and vol-surf (cell-cell, cell-face) blocks -- which
    matter for NON-uniform M (nonlinear, div M != 0) -- are analytic, not monopole+sub-point.  This is
    the #B 'volume Gram' that closes the residual sharp-body non-uniform gap (cube ~8.7%, C-yoke ~4%).
    Symmetric; O(n_charge) Python loop (vectorized over target points).

    TET cells + TRIANGLE faces take the EXACT phi_tet / tri_potential + Gauss-Duffy / Dunavant path, so
    an all-tet/all-triangle mesh is BIT-IDENTICAL to the original (the C++ analytic charge-Gram golden,
    test_hdiv_vim_chargegram_analytic, locks C++ == this dense G).  HEX / WEDGE cells + QUAD faces take
    the convex-hull polytope path -- verified demag_z -> 1/3 on hex- and wedge-meshed cubes."""
    n_el, n_bf = len(el_V), len(bf_V)
    n = n_el + n_bf
    lam_t, w_t = _gauss_duffy_tet(4)                         # Gauss-Duffy tet outer rule (nsub_tet now unused:
    QP, QW = [], []                                          # the old _bary_tet outer under-integrated vol-self ~6.5%)
    src = []                                                 # per-source potential evaluator descriptor
    for k, V in enumerate(el_V):
        if len(V) == 4:
            QP.append(V[0] + lam_t @ (V[1:] - V[0]))         # physical outer pts (V0 + l1 e1 + l2 e2 + l3 e3)
            QW.append(w_t * 6.0 * el_vol[k])                 # phys weight = w_ref * |J|, |J| = 6*vol
            src.append(("tet", V))
        else:
            P, W = _cell_outer_quad(V, el_vol[k], lam_t, w_t)
            QP.append(P); QW.append(W); src.append(("poly", _cell_hull_tris(V)))
    for k, V in enumerate(bf_V):
        if len(V) == 3:
            QP.append(_DUN5[:, :3] @ V); QW.append(_DUN5[:, 3] * bf_area[k])
            src.append(("tri", V))
        else:
            P, W = _face_outer_quad(V, bf_area[k])
            QP.append(P); QW.append(W); src.append(("face", _face_subtris(V)))
    allP = np.vstack(QP)
    allW = np.concatenate(QW)
    tidx = np.concatenate([np.full(len(QP[i]), i) for i in range(n)])
    inv4pi = 1.0 / (4.0 * pi)
    G = np.zeros((n, n))
    for s in range(n):
        kind, data = src[s]
        if kind == "tet":
            pot = phi_tet(data, allP)
        elif kind == "tri":
            pot = tri_potential(data, allP)
        elif kind == "poly":
            pot = _polytope_potential(data, allP)
        else:                                                # "face": quad -> sum of its sub-triangle potentials
            pot = sum(tri_potential(T, allP) for T in data)
        G[:, s] = np.bincount(tidx, weights=allW * pot * inv4pi, minlength=n)
    return 0.5 * (G + G.T)


def parse_image_string(image):
    """Parse a Radia IMA image string ('+x-z', '-x+y-z', ...) into a list of (axis, sign): axis in
    {0,1,2} for x/y/z, sign +1 (symmetric across that plane) / -1 (antisymmetric).  The string is a
    concatenation of (+|-)(x|y|z) tokens.  Matches the rad.Solve(image=...) convention (CLAUDE.md
    IMA Sign Selection: field PARALLEL to the mirror -> '+', PERPENDICULAR -> '-')."""
    s = image.strip().lower().replace(" ", "")
    planes, i = [], 0
    axis_of = {"x": 0, "y": 1, "z": 2}
    while i < len(s):
        if s[i] not in "+-" or i + 1 >= len(s) or s[i + 1] not in axis_of:
            raise ValueError("bad IMA image string %r (expected tokens like '+x','-z')" % image)
        planes.append((axis_of[s[i + 1]], 1 if s[i] == "+" else -1))
        i += 2
    axes = [a for a, _ in planes]
    if len(set(axes)) != len(axes):
        raise ValueError("IMA image string %r repeats an axis" % image)
    return planes


def image_group(planes):
    """Image-method reflection group from parsed planes [(axis,sign),...]: every NON-EMPTY subset of
    the mirror planes gives one image, whose reflection negates the coords on those axes and whose
    scalar sign is the PRODUCT of the per-plane signs (standard image method; the paper eq 23
    G^IMA = G + sum_m s_m G_{i,m(j)}).  Returns [(axes_tuple, sign), ...] (2^P - 1 images)."""
    out = []
    P = len(planes)
    for mask in range(1, 1 << P):
        axes, sign = [], 1
        for k in range(P):
            if mask & (1 << k):
                axes.append(planes[k][0]); sign *= planes[k][1]
        out.append((tuple(sorted(axes)), sign))
    return out


def _reflect(V, axes):
    """Reflect vertex coords V (n,3) by negating the columns in `axes`."""
    V = np.asarray(V, float).copy()
    for a in axes:
        V[:, a] = -V[:, a]
    return V


def image_charge_gram(el_V, bf_V, el_vol, bf_area, axes):
    """Cross charge Gram G_img[a][b] = (1/4pi) INT_{target a} Phi_{R_axes(b)} -- the analytic Coulomb
    interaction of each real charge a with charge b's geometry REFLECTED on `axes` (negate those
    coords).  Same outer quadrature + source potentials as analytic_charge_gram, but the SOURCE
    geometry is reflected.  Symmetric (G(a,R(b))=G(b,R(a)) by kernel symmetry + R isometry), so the
    image term keeps N_IMA symmetric.  Used to fold mirror symmetry (IMA) into the demag operator."""
    n_el, n_bf = len(el_V), len(bf_V)
    n = n_el + n_bf
    lam_t, w_t = _gauss_duffy_tet(4)
    QP, QW = [], []
    for k, V in enumerate(el_V):
        if len(V) == 4:
            QP.append(V[0] + lam_t @ (V[1:] - V[0])); QW.append(w_t * 6.0 * el_vol[k])
        else:
            P, W = _cell_outer_quad(V, el_vol[k], lam_t, w_t); QP.append(P); QW.append(W)
    for k, V in enumerate(bf_V):
        if len(V) == 3:
            QP.append(_DUN5[:, :3] @ V); QW.append(_DUN5[:, 3] * bf_area[k])
        else:
            P, W = _face_outer_quad(V, bf_area[k]); QP.append(P); QW.append(W)
    allP = np.vstack(QP); allW = np.concatenate(QW)
    tidx = np.concatenate([np.full(len(QP[i]), i) for i in range(n)])
    inv4pi = 1.0 / (4.0 * pi)

    def refl_src_pot(s):                                  # potential of source s reflected on `axes`
        if s < n_el:
            V = _reflect(el_V[s], axes)
            return phi_tet(V, allP) if len(V) == 4 else _polytope_potential(_cell_hull_tris(V), allP)
        V = _reflect(bf_V[s - n_el], axes)
        return sum(tri_potential(T, allP) for T in _face_subtris(V))

    G = np.zeros((n, n))
    for s in range(n):
        G[:, s] = np.bincount(tidx, weights=allW * refl_src_pot(s) * inv4pi, minlength=n)
    return 0.5 * (G + G.T)


def polytope_flat_geom(el_V, bf_V, el_vol, bf_area):
    """Flatten the per-cell convex-hull triangles + per-face sub-triangles (+ vertex-mean centroids and
    measures) into the flat triangle-soup arrays the C++ POLYTOPE _ChargeGramHMatrix constructor consumes
    (hex/wedge scalable path).  Cells -> _cell_hull_tris, boundary faces -> _face_subtris; the C++ rebuilds
    the centroid-fan (cells) / Dunavant (faces) outer quadrature and the divergence-theorem source potential
    from these, so it matches the dense analytic_charge_gram polytope path entry-by-entry (SAME triangles,
    SAME built-in quad rules).  Returns flat float arrays + int CSR offsets (triangle soup = 9 doubles/tri)."""
    cell_tris, cell_troff, cell_cent = [], [0], []
    for V in el_V:
        V = np.asarray(V, float)
        cell_cent.append(V.mean(0))
        tris = _cell_hull_tris(V)
        for (P, _n) in tris:
            cell_tris.append(np.asarray(P, float).ravel())
        cell_troff.append(cell_troff[-1] + len(tris))
    face_tris, face_troff, face_cent = [], [0], []
    for V in bf_V:
        V = np.asarray(V, float)
        face_cent.append(V.mean(0))
        subs = _face_subtris(V)
        for T in subs:
            face_tris.append(np.asarray(T, float).ravel())
        face_troff.append(face_troff[-1] + len(subs))
    return dict(
        cell_tris=(np.concatenate(cell_tris) if cell_tris else np.zeros(0)),
        cell_troff=np.asarray(cell_troff, np.int32),
        cell_cent=(np.asarray(cell_cent, float).ravel() if cell_cent else np.zeros(0)),
        cell_meas=np.asarray(el_vol, float),
        face_tris=(np.concatenate(face_tris) if face_tris else np.zeros(0)),
        face_troff=np.asarray(face_troff, np.int32),
        face_cent=(np.asarray(face_cent, float).ravel() if face_cent else np.zeros(0)),
        face_meas=np.asarray(bf_area, float),
    )


def _csr_sp(bf):
    """NGSolve BilinearForm.mat -> scipy CSR (SPARSE). The scalable (skip_dense_gram) build keeps the
    charge maps + HDiv mass sparse; only the dense reference path densifies via _csr()."""
    m = bf.mat
    r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width))


def _csr(bf):
    return _csr_sp(bf).toarray()


def build_demag(mesh, nsub=4, wilton_surface=False, analytic_gram=False, skip_dense_gram=False,
                image=None):
    """Assemble the HDiv-type VIM demag operator N = B^T G B on a tet mesh + the HDiv mass M_mass.
    Returns N, M_mass, the charge map B, the loop basis, and diagnostics.

    image (IMA mirror symmetry, dense analytic path only): a Radia image string ('+x-z', ...).  The
    mesh is the reduced (1/2, 1/4, 1/8) model; the demag operator folds in the mirror-image charge
    interactions G_IMA = G + sum_S sign(S) image_charge_gram(.., S) so the reduced solve reproduces the
    full model (validated quarter+IMA == full to machine precision).  Requires analytic_gram=True.

    wilton_surface=True replaces the surface-surface (boundary-triangle) OFF-diagonal Gram block with
    the exact Wilton analytic integral (the diagonal keeps the validated tri_self_energy).  This makes
    the demag factor exact to <0.15% on cube AND sphere (vs the ~5-6% monopole error that the sub-point
    near-correction cannot fix on the cube) -- the #1 production accuracy fix.  O(n_bf^2) dense; the
    scalable C++ HACApK charge-Gram path is the next step."""
    # NGSolve assembly (FES, charge-map BilinearForms, HDiv mass).  Per the caller-wraps policy this HELPER
    # does NOT open a TaskManager; the CALLER wraps build_demag in `with ng.TaskManager():` (CLAUDE.md
    # "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT").
    fes = ng.HDiv(mesh, order=0)
    ndof = fes.ndof
    nn = ng.specialcf.normal(mesh.dim)
    L2v, L2b = ng.L2(mesh, order=0), ng.SurfaceL2(mesh, order=0)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v)
    bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx
    bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b)
    bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds
    bb.Assemble()
    Bv_sp, Bb_sp = _csr_sp(bv), _csr_sp(bb)   # charge maps as sparse CSR (densified only on the ref path)
    massv = ng.BilinearForm(L2v); massv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; massv.Assemble()
    massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; massb.Assemble()
    el_vol = _csr_sp(massv).diagonal(); bf_area = _csr_sp(massb).diagonal()   # L2-0 mass is diagonal == measures
    # HDiv mass (the physical demag-factor metric) -- kept SPARSE; the skip_dense_gram path never densifies it
    vh = fes.TestFunction()
    mh = ng.BilinearForm(fes); mh += u * vh * ng.dx; mh.Assemble()
    M_mass_sp = _csr_sp(mh)
    el_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.VOL)]
    bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
    # uniform M_z, L2-projected onto RT0 (for the demag Rayleigh quotient)
    gfu = ng.GridFunction(fes); gfu.Set(ng.CoefficientFunction((0, 0, 1))); m_unit = np.array(gfu.vec)
    el_c = np.array([V.mean(0) for V in el_V]); bf_c = np.array([V.mean(0) for V in bf_V])
    n_el, n_bf = len(el_c), len(bf_c)
    cent = np.vstack([el_c, bf_c]); meas = np.concatenate([el_vol, bf_area])
    # scaled charge map B (each row / its measure) as sparse CSR -- the single source for both B_csr and,
    # on the dense reference path, the dense B.  (sparse: each charge cell touches only its own faces.)
    _Brows = [sp.diags(1.0 / el_vol) @ Bv_sp]
    if n_bf:
        _Brows.append(sp.diags(1.0 / bf_area) @ Bb_sp)
    B_sp = sp.vstack(_Brows).tocsr()

    # diagonal charge self-energies (O(N); kept even on the scalable path, for the preconditioner).
    # tet/triangle use the validated tet_self_energy/tri_self_energy (unchanged); hex/wedge cells +
    # quad faces use the analytic polytope self-energy (EXACT-quadrature) so a non-tet mesh does not
    # crash here and gets a physical monopole-path diagonal.
    diagG = np.empty(n_el + n_bf)
    for k, V in enumerate(el_V):
        diagG[k] = tet_self_energy(V, el_vol[k], nsub) if len(V) == 4 else _poly_cell_self_energy(V, el_vol[k])
    for k, V in enumerate(bf_V):
        diagG[n_el + k] = (tri_self_energy(V, bf_area[k], nsub) if len(V) == 3
                           else _poly_face_self_energy(V, bf_area[k]))
    if skip_dense_gram:
        # SCALABLE path: SKIP the O(N^2) dense G/N and the O(N^2) loop SVD (the biggest dense costs;
        # G/N are n_charge^2).  The C++ analytic _ChargeGramHMatrix + SolveMaterialMINRES use
        # cent / meas / cell_verts / face_verts / diagG + the sparse B_csr only -- never the dense G.
        # M_mass + B are kept SPARSE here (no n_charge^2 dense), so this path has NO dense N^2 object at
        # all -- the analytic-Gram BUILD needs only the verts (straight from the mesh).
        G = N = loops = None
        n_loop = None
        M_mass = M_mass_sp        # sparse CSR (RT0 HDiv mass is local -> sparse)
        B = B_sp                  # sparse CSR scaled charge map
    else:
        M_mass = M_mass_sp.toarray()                       # dense reference path: small N
        B = B_sp.toarray()
        Dd = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2); np.fill_diagonal(Dd, np.inf)
        G = (meas[:, None] * meas[None, :]) / (4 * pi * Dd)
        np.fill_diagonal(G, diagG)
        if analytic_gram:
            # FULL analytic charge Gram (#B volume Gram): every block (vol-vol, vol-surf, surf-surf) via
            # the exact phi_tet / tri_potential potentials -- closes the NON-uniform-M (nonlinear) gap.
            G = analytic_charge_gram(el_V, bf_V, el_vol, bf_area)
            if image:
                # IMA: fold in the mirror-image charge interactions (G_IMA = G + sum_S sign(S) G_img_S),
                # so the reduced (1/2,1/4,1/8) mesh reproduces the full model.  Same per-image sign for
                # cell rho + face sigma (pinned empirically: quarter+IMA == full to 1e-8).
                for axes, sign in image_group(parse_image_string(image)):
                    G = G + sign * image_charge_gram(el_V, bf_V, el_vol, bf_area, axes)
        elif image:
            raise ValueError("build_demag: image= (IMA) requires analytic_gram=True")
        elif wilton_surface and n_bf > 0:
            # exact Wilton surface-surface block (shape-exact for ANY triangle; demag 1/3 to <0.15%).
            G[n_el:, n_el:] = wilton_surface_block(bf_V)
        N = B.T @ G @ B
        Q = np.vstack([Bv_sp.toarray(), Bb_sp.toarray()]) if n_bf else Bv_sp.toarray()
        sv = np.linalg.svd(Q, compute_uv=False)
        rankQ = int(np.sum(sv > 1e-9 * sv.max()))
        n_loop = ndof - rankQ
        _, _, Vt = np.linalg.svd(Q)
        loops = Vt[rankQ:, :]
    # charge geometry (for the C++ HACApK charge-Gram H-matrix path, #1b): centroids, measures,
    # diagonal self-energies, the dense Gram G, and the sparse charge map B as scipy CSR.
    # cell_verts [n_el*12] (tets, 4 verts) + face_verts [n_bf*9] (tris, 3 verts) feed the ANALYTIC
    # C++ charge Gram (M2b: _ChargeGramHMatrix analytic mode == this dense analytic_gram path).
    # cell_verts/face_verts feed ONLY the C++ scalable charge-Gram H-matrix (tet cells + triangle
    # faces); ravel only when every cell/face has the same vertex count (uniform tet/triangle mesh).
    # A hex/wedge or mixed mesh leaves these empty -- it takes the dense analytic polytope path instead.
    cell_verts = (np.asarray(el_V, float).ravel()
                  if n_el and len({len(V) for V in el_V}) == 1 else np.zeros(0))
    face_verts = (np.asarray(bf_V, float).ravel()
                  if n_bf and len({len(V) for V in bf_V}) == 1 else np.zeros(0))
    # POLYTOPE flat geometry for the C++ scalable HEX/WEDGE charge-Gram H-matrix (the polytope
    # _ChargeGramHMatrix ctor): computed only when the mesh has non-tet cells / non-triangle faces (a pure
    # tet mesh takes the bit-identical cell_verts/face_verts path, so the per-cell ConvexHull is skipped).
    need_poly = any(len(V) != 4 for V in el_V) or any(len(V) != 3 for V in bf_V)
    poly = (polytope_flat_geom(el_V, bf_V, el_vol, bf_area) if need_poly else None)
    return dict(N=N, M_mass=M_mass, B=B, ndof=ndof, n_loop=n_loop, loops=loops, m_unit=m_unit,
                cent=cent, meas=meas, self_energy=diagG, G=G,
                B_csr=B_sp, n_charge=n_el + n_bf,
                cell_verts=cell_verts, face_verts=face_verts, n_el=n_el, poly=poly)


def demag_factor(d):
    m = d["m_unit"]
    return float((m @ d["N"] @ m) / (m @ d["M_mass"] @ m))


def _charge_subpoints(mesh, nsub):
    """Sub-points + equal weights per charge (cells then boundary faces), element-type aware:
    tet cells use the (fixed) barycentric tet lattice, boundary tri faces the triangle lattice."""
    el_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.VOL)]
    bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
    lam_t, lam_f = _bary_tet(nsub), _bary_tri(nsub)
    SP, SW, kind = [], [], []
    for V in el_V:
        vol = abs(np.linalg.det(V[1:] - V[0])) / 6.0
        SP.append(lam_t @ V); SW.append(np.full(len(lam_t), vol / len(lam_t))); kind.append(0)
    for V in bf_V:
        ar = 0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0]))
        SP.append(lam_f @ V); SW.append(np.full(len(lam_f), ar / len(lam_f))); kind.append(1)
    return SP, SW, kind


def build_near_correction(mesh, d, nsub=4, near_factor=2.0, skip_surface_surface=False):
    """Sparse near-field Gram correction (exact sub-point MINUS centroid-monopole) for NEAR charge
    pairs -- the standard H-matrix near-field correction.  The scalable Gram is then the compressed
    monopole H-matrix (far, cheap) PLUS this sparse local correction (near, exact): it lifts the
    sphere demag from the monopole ~0.31 to the Gram-exact ~0.33 (-> analytic 1/3) at O(N) extra cost.
    Cells contribute ~0 for uniform M (div M = 0), but the correction is computed for all near pairs
    for generality.  Returns a scipy CSR (n_charge x n_charge).

    skip_surface_surface=True omits surface-surface pairs -- use this together with the Wilton surface
    Gram (build_demag(wilton_surface=True)), which already makes the surface-surface block exact; the
    near-correction then only fixes the VOLUME-involving (cell-cell, cell-face) near pairs that the
    per-element NONLINEAR Newton needs (without it the volume near-field is un-corrected and Newton
    finds a wrong root)."""
    import scipy.sparse as sp
    inv4pi = 1.0 / (4.0 * pi)
    cent, meas = d["cent"], d["meas"]
    n_cell = d["n_charge"] - sum(1 for _ in mesh.Elements(ng.BND))
    SP, SW, _ = _charge_subpoints(mesh, nsub)
    n = len(cent)
    size = np.array([meas[a] ** (1.0 / 3.0) if a < n_cell else meas[a] ** 0.5 for a in range(n)])
    rows, cols, vals = [], [], []
    for a in range(n):
        ca, sa, wa = cent[a], SP[a], SW[a]
        for b in range(a + 1, n):
            if skip_surface_surface and a >= n_cell and b >= n_cell:
                continue                         # surface-surface handled exactly by the Wilton Gram
            dx = ca - cent[b]
            r = float(np.sqrt(dx @ dx))
            if r < near_factor * (size[a] + size[b]):
                D = np.linalg.norm(sa[:, None, :] - SP[b][None, :, :], axis=2)
                exact = float(np.sum(np.outer(wa, SW[b]) * inv4pi / D))
                mono = meas[a] * meas[b] * inv4pi / r
                delta = exact - mono
                rows += [a, b]; cols += [b, a]; vals += [delta, delta]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def report(mesh, tag, nsub=4):
    d = build_demag(mesh, nsub)
    Nn = np.linalg.norm(d["N"], 2)
    asym = np.linalg.norm(d["N"] - d["N"].T) / Nn
    loop_res = (max(np.linalg.norm(d["N"] @ d["loops"][k]) for k in range(d["n_loop"])) / Nn
                if d["n_loop"] else 0.0)
    Dz = demag_factor(d)
    print(f"[{tag}] ndof={d['ndof']} n_loop={d['n_loop']} asym={asym:.1e} loop_res={loop_res:.1e} "
          f"demag_z={Dz:.4f}")
    return {"tag": tag, "ndof": d["ndof"], "n_loop": d["n_loop"], "asym": asym,
            "loop_res": loop_res, "demag_z": Dz}


if __name__ == "__main__":
    res = {}
    # cube (tet-meshed): each axis demag factor ~ 1/3
    geo = CSGeometry(); geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    res["cube"] = report(ng.Mesh(geo.GenerateMesh(maxh=0.25)), "cube tet", nsub=4)
    # sphere (tet-meshed): demag factor EXACTLY 1/3
    for h in (0.5, 0.35, 0.25):
        geo = CSGeometry(); geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
        res[f"sphere_h{h}"] = report(ng.Mesh(geo.GenerateMesh(maxh=h)), f"sphere tet h={h}", nsub=4)
    with open(os.path.join(HERE, "hdiv_demag_tet.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_demag_tet.json"))
