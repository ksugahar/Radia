"""Loop-free collocation surface-charge (yano-type MSC) prototype: pyramid-cloud kernel + star
projection.  Python validation gate BEFORE the C++ port into the yano-MSC kernel.

Motivation (from `yano_hdiv_loop_bridge.py`): the yano-type collocation matrix carries the cell-graph
cycle space as a latent near-null subspace, ill-conditioned at high mu_r.  Two levers fix it -- and
this prototype implements BOTH and validates the conditioning win on a distorted multi-hex patch:

  1. **Star projection** -- restrict the unknowns from the 6-sigma/hex local-face space to the
     internal-flux subspace q_internal = B^T a (B = cell/internal-face incidence).  The loop
     coefficients ker(B) are then NOT unknowns, so an iterative solve cannot grow them.  Sparse,
     O(n_face).  This is exactly the HDiv-VIM RT0 DOF space, reached from the collocation side.

  2. **Pyramid-cloud kernel** -- replace the single-point compensation (one point charge at the
     element center) with an ELEMENT-COMMON cloud distributed on the 12 (centroid, edge) partition
     triangles, area-weighted, evaluated at the pyramid centroid (0.75*face_area_centroid +
     0.25*trilinear_volume_centroid).  Element-common (shared by all 6 face DOFs) -> the loop SOURCE
     field cancels at machine precision; charge-neutral per DOF -> NO spurious monopole (physics
     preserved); and it conditions the self-block far better than the single point.

Validated here (all on the SAME distorted patch, single-point vs pyramid-cloud):
  (a) source clouds are charge-neutral per DOF (sum of weights ~ 0)  -> physics-preserving;
  (b) the full local-face matrix near-null dim == the cell-graph cycle count (the latent loops);
  (c) star projection drives near-null -> 0 and the condition number from ~1e16 to O(10-100);
  (d) the pyramid-cloud star-projected condition number is markedly LOWER than the single-point one
      (the lab reference on a 100-element distorted patch: single-point 107 vs pyramid-cloud 20.5),
      with correspondingly fewer BiCGStab iterations.

The C++ port then puts this kernel + star projection behind the yano-MSC `Use6DOF_MSC` solve, gated by
the existing demag-factor goldens (sphere/cube -> 1/3) for the physics.

KERNEL CORRECTNESS (cross-validated): on a canonical 100-element distorted patch this exact kernel +
star projection gives star-projected condition number **20.49** (single-point compensation gives
107.4) with BiCGStab converging in **26-31** iterations (single-point 56-60) -- the published
reference figures.  The numbers PRINTED BELOW use this example's OWN (different) distortion, so they
show the same TREND (pyramid-cloud markedly beats single-point, both loop-free after star projection)
at this example's geometry rather than the reference 20.49.


Reference: H. Yano and K. Sugahara, "Magnetic Moment Method with the Idea of Magnetic Surface Charge
Method", J. Magn. Soc. Jpn. 47(2), 52-56, 2023.  The loop / de Rham reading is the FEEC RT0 structure
(see `yano_hdiv_loop_bridge.py`, README).
"""
import os
import json
import math

import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab

HERE = os.path.dirname(os.path.abspath(__file__))

# hex local topology (VTK-style): 6 quad faces + 12 edges, consistent winding
HEX_FACE_NODES = np.array([[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                           [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]], dtype=int)
HEX_EDGES = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4],
                      [0, 4], [1, 5], [2, 6], [3, 7]], dtype=int)
_BARY3 = np.array([[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])


# ----- hex geometry -----
def _tri_area(a, b, c):
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def _face_quadrature(V, face):
    """6-point quadrature over a quad face (2 triangles x 3 bary points)."""
    q = V[HEX_FACE_NODES[face]]
    pts, wts = [], []
    for tri in ((0, 1, 2), (0, 2, 3)):
        a, b, c = (q[i] for i in tri)
        area = _tri_area(a, b, c)
        for w in _BARY3:
            pts.append(w[0] * a + w[1] * b + w[2] * c)
            wts.append(area / 3.0)
    return np.array(pts), np.array(wts)


def _face_area(V, face):
    return float(_face_quadrature(V, face)[1].sum())


def _face_area_centroid(V, face):
    pts, w = _face_quadrature(V, face)
    return (w[:, None] * pts).sum(0) / w.sum()


def _face_normal(V, face, center):
    p = V[HEX_FACE_NODES[face]]
    n = np.cross(p[1] - p[0], p[2] - p[0])
    nn = np.linalg.norm(n)
    if nn <= 1e-300:
        return np.zeros(3)
    n = n / nn
    return -n if np.dot(n, p.mean(0) - center) < 0 else n


def _trilinear_centroid(V):
    """2x2x2 Gauss volume centroid (the trilinear hex centroid)."""
    g = 1.0 / math.sqrt(3.0)
    s = np.array([-g, g])
    num, den = np.zeros(3), 0.0
    for a in s:
        for b in s:
            for c in s:
                N = 0.125 * np.array([(1 - a) * (1 - b) * (1 - c), (1 + a) * (1 - b) * (1 - c),
                                      (1 + a) * (1 + b) * (1 - c), (1 - a) * (1 + b) * (1 - c),
                                      (1 - a) * (1 - b) * (1 + c), (1 + a) * (1 - b) * (1 + c),
                                      (1 + a) * (1 + b) * (1 + c), (1 - a) * (1 + b) * (1 + c)])
                dN = 0.125 * np.array([
                    [-(1 - b) * (1 - c), -(1 - a) * (1 - c), -(1 - a) * (1 - b)],
                    [(1 - b) * (1 - c), -(1 + a) * (1 - c), -(1 + a) * (1 - b)],
                    [(1 + b) * (1 - c), (1 + a) * (1 - c), -(1 + a) * (1 + b)],
                    [-(1 + b) * (1 - c), (1 - a) * (1 - c), -(1 - a) * (1 + b)],
                    [-(1 - b) * (1 + c), -(1 - a) * (1 + c), (1 - a) * (1 - b)],
                    [(1 - b) * (1 + c), -(1 + a) * (1 + c), (1 + a) * (1 - b)],
                    [(1 + b) * (1 + c), (1 + a) * (1 + c), (1 + a) * (1 + b)],
                    [-(1 + b) * (1 + c), (1 - a) * (1 + c), (1 - a) * (1 + b)]])
                J = V.T @ dN
                detJ = abs(np.linalg.det(J))
                num += detJ * (N @ V)
                den += detJ
    return num / den if den > 1e-300 else V.mean(0)


def _surface_cloud_from_triangles(tris):
    """3-pt quadrature cloud over triangles, normalized to unit total charge (element-common)."""
    pts, wts = [], []
    for a, b, c in tris:
        area = _tri_area(a, b, c)
        if area <= 1e-300:
            continue
        for w in _BARY3:
            pts.append(w[0] * a + w[1] * b + w[2] * c)
            wts.append(area / 3.0)
    w = np.array(wts)
    return np.array(pts), w / w.sum()


def _compensation_cloud(V, center, layout):
    """Element-COMMON compensation cloud (source-face-independent -> loop-source-null)."""
    if layout == "single":
        return center.reshape(1, 3), np.array([1.0])
    if layout == "pyr_faces12":               # 12 (centroid, edge) partition triangles, area-weighted
        tris = [(center, V[a], V[b]) for a, b in HEX_EDGES]
        return _surface_cloud_from_triangles(tris)
    raise ValueError(layout)


def _source_cloud(V, face, center, layout):
    """One unit-density face DOF: physical face charge (+, total=area) minus area*common_cloud (-).
    Net charge == 0 (monopole-free) -> physics-preserving."""
    fp, fw = _face_quadrature(V, face)             # physical face charge, total = area
    area = fw.sum()
    cp, cw = _compensation_cloud(V, center, layout)
    pts = np.vstack([fp, cp])
    wts = np.concatenate([fw, -area * cw])         # compensation total = -area  -> neutral
    return pts, wts


# ----- collocation matrix (monopole field; the 1/4pi cancels in the conditioning) -----
def collocation_matrix(hexes, layout):
    n = len(hexes)
    centers = [_trilinear_centroid(V) for V in hexes]
    eval_pts = np.zeros((6 * n, 3))
    normals = np.zeros((6 * n, 3))
    for e, V in enumerate(hexes):
        for f in range(6):
            r = 6 * e + f
            eval_pts[r] = 0.75 * _face_area_centroid(V, f) + 0.25 * centers[e]   # pyramid centroid
            normals[r] = _face_normal(V, f, centers[e])
    clouds = [_source_cloud(hexes[e], f, centers[e], layout) for e in range(n) for f in range(6)]
    A = np.zeros((6 * n, 6 * n))
    neutral = max(abs(w.sum()) for _, w in clouds)         # physics gate: per-DOF net charge ~ 0
    for col, (sp, sw) in enumerate(clouds):
        diff = eval_pts[:, None, :] - sp[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        m = dist > 1e-14
        inv = np.zeros_like(dist)
        inv[m] = 1.0 / dist[m] ** 3
        fld = np.sum(sw[None, :, None] * diff * inv[:, :, None], axis=1)
        A[:, col] = -np.einsum("ij,ij->i", normals, fld)
    return A, neutral


# ----- star projection T (boundary identity + internal sign*row_space_basis/area) -----
def face_topology(elements):
    adj = {}
    for ei, el in enumerate(elements):
        for fi, ln in enumerate(HEX_FACE_NODES):
            key = tuple(sorted(int(el[n]) for n in ln))
            adj.setdefault(key, []).append((ei, fi))
    internal = [k for k, a in adj.items() if len(a) == 2]
    boundary = [k for k, a in adj.items() if len(a) == 1]
    return internal, boundary, adj


def star_transform(hexes, elements, internal, boundary, adj):
    B = np.zeros((len(elements), len(internal)))
    for col, key in enumerate(internal):
        for sgn, (ei, _f) in zip((1.0, -1.0), adj[key]):
            B[ei, col] = sgn
    _u, sv, vh = np.linalg.svd(B, full_matrices=True)
    rank = int(np.sum(sv > max(B.shape) * np.finfo(float).eps * (sv[0] if sv.size else 1.0)))
    star_basis = vh[:rank, :]
    T = np.zeros((6 * len(elements), len(boundary) + rank))
    for col, key in enumerate(boundary):
        ei, fi = adj[key][0]
        T[6 * ei + fi, col] = 1.0
    for ic, key in enumerate(internal):
        for sgn, (ei, fi) in zip((1.0, -1.0), adj[key]):
            area = _face_area(hexes[ei], fi)
            if area <= 1e-300:
                continue
            for sc in range(rank):
                T[6 * ei + fi, len(boundary) + sc] += sgn * star_basis[sc, ic] / area
    return T


# ----- distorted multi-hex patch (shared global nodes; per-node twist+jitter -- pure Python, no
#       planar-face constraint, so genuine hex distortion is allowed here) -----
def distorted_patch(nx, ny, nz, strength, seed=7):
    rng = np.random.default_rng(seed)
    nid = {}
    coords = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid[(i, j, k)] = len(coords)
                coords.append([i - 0.5 * nx, j - 0.5 * ny, k - 0.5 * nz])
    P = np.array(coords, float)
    # twist about z + shear + interior jitter
    for n, (x, y, z) in enumerate(P):
        ang = strength * 0.10 * z
        c, s = math.cos(ang), math.sin(ang)
        P[n, 0] = c * x - s * y + strength * 0.05 * y
        P[n, 1] = s * x + c * y - strength * 0.04 * x
        P[n, 2] = z * (1 + strength * 0.02 * x) + strength * 0.04 * math.sin(math.pi * x / max(nx, 1))
    for (i, j, k), n in nid.items():
        if 0 < i < nx and 0 < j < ny and 0 < k < nz:
            P[n] += strength * 0.10 * (rng.random(3) - 0.5)
    P -= P.mean(0)
    elements = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elements.append([nid[(i, j, k)], nid[(i + 1, j, k)], nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
                                 nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)], nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]])
    hexes = [P[np.array(el)] for el in elements]
    return hexes, elements


def _bicgstab_iters(A, rng):
    cnt = {"n": 0}
    op = LinearOperator(A.shape, matvec=lambda v: A @ v)
    res = []
    for kind, b in (("random", rng.random(A.shape[0]) - 0.5),
                    ("uniform", np.ones(A.shape[0])),
                    ("quad", np.linspace(-1, 1, A.shape[0]))):
        cnt["n"] = 0
        _x, info = bicgstab(op, b, rtol=1e-10, maxiter=20 * A.shape[0],
                            callback=lambda xk: cnt.__setitem__("n", cnt["n"] + 1))
        res.append((kind, cnt["n"], info))
    return res


def evaluate(nx, ny, nz, strength):
    hexes, elements = distorted_patch(nx, ny, nz, strength)
    internal, boundary, adj = face_topology(elements)
    cycle = len(internal) - (len(elements) - 1)
    print(f"\n=== {nx}x{ny}x{nz} strength={strength} : n_elem={len(elements)} "
          f"internal={len(internal)} boundary={len(boundary)} cycle(loops)={cycle} ===")
    rng = np.random.default_rng(0)
    rows = {}
    for layout in ("single", "pyr_faces12"):
        A, neutral = collocation_matrix(hexes, layout)
        sv = np.linalg.svd(A, compute_uv=False)
        near_null = int(np.sum(sv < 1e-9 * sv[0]))
        cond_full = sv[0] / sv[-1] if sv[-1] > 0 else float("inf")
        T = star_transform(hexes, elements, internal, boundary, adj)
        Ar = T.T @ A @ T
        svr = np.linalg.svd(Ar, compute_uv=False)
        cond_star = svr[0] / svr[-1] if svr[-1] > 0 else float("inf")
        nn_star = int(np.sum(svr < 1e-9 * svr[0]))
        iters = _bicgstab_iters(Ar, rng)
        it_str = " ".join(f"{k}={n}({i})" for k, n, i in iters)
        print(f"  [{layout:>11}] neutral(max|sum w|)={neutral:.1e}  full: cond={cond_full:.2e} "
              f"near_null={near_null}  | star({Ar.shape[0]}): cond={cond_star:.2f} near_null={nn_star}  "
              f"bicgstab[{it_str}]")
        rows[layout] = dict(neutral=neutral, cond_full=cond_full, near_null=near_null,
                            star_dim=Ar.shape[0], cond_star=cond_star, near_null_star=nn_star,
                            bicgstab={k: n for k, n, _ in iters})
    s, p = rows["single"]["cond_star"], rows["pyr_faces12"]["cond_star"]
    print(f"  [verdict   ] loop bridge near_null==cycle: "
          f"{rows['single']['near_null']==cycle==rows['pyr_faces12']['near_null']};  "
          f"star near_null->0: {rows['pyr_faces12']['near_null_star']==0};  "
          f"pyr cond {p:.1f} < single cond {s:.1f}: {p < s} ({s/p:.2f}x)")
    return dict(grid=f"{nx}x{ny}x{nz}", strength=strength, n_elem=len(elements), cycle=cycle, **rows)


def main():
    out = [evaluate(3, 3, 2, 1.0), evaluate(5, 5, 4, 4.0)]
    ok = all(r["pyr_faces12"]["near_null_star"] == 0
             and r["pyr_faces12"]["near_null"] == r["cycle"]
             and abs(r["single"]["neutral"]) < 1e-9 and abs(r["pyr_faces12"]["neutral"]) < 1e-9
             and r["pyr_faces12"]["cond_star"] < r["single"]["cond_star"] for r in out)
    with open(os.path.join(HERE, "yano_pyr_faces12_star.json"), "w") as f:
        json.dump({"all_ok": ok, "cases": out}, f, indent=2, default=float)
    print(f"\nall_ok (neutral + loop-free + pyr beats single): {ok}")
    return ok


if __name__ == "__main__":
    main()
