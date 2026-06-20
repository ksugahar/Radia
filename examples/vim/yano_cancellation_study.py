"""yano_cancellation_study.py -- RESEARCH (user idea): put the cancellation charge on the hex EDGES
instead of the centroid.  The previous study (yano_loopfree_kernel_study.py) found the eval point and
the cancellation are SEPARABLE knobs: eval_alpha=0.5 (EIEM2) preserves accuracy; the cancellation
cloud shape sets the star-block conditioning (single cond_star~84, pyr_faces12~14.6).  Lower cond_star
-> fewer Krylov iters.  This sweeps cancellation-charge LAYOUTS at the accuracy-preserving EIEM2 eval
point (alpha=0.5), all element-common + charge-neutral, to find the best-conditioned cancellation:

  single        : 1 point at the trilinear centroid (EIEM2 default)
  pyr_faces12   : cloud on the 12 (centroid,edge) triangles, area-weighted (the prior best)
  edges_mid     : 12 edge MIDPOINTS, equal weight                       <- user's idea
  edges_len     : 12 edge midpoints, edge-LENGTH weighted               <- user's idea (len-aware)
  edges_q       : 3-pt quadrature ALONG each of the 12 edges            <- user's idea (distributed)
  faces6        : 6 face-area-centroids, area-weighted (for contrast)

Reports cond_star + BiCGStab iters (lower=better), and confirms (a) charge-neutral per DOF and
(b) loops stay topological (near_null==cycle; cancellation never lifts them -- prior finding).
Monopole field => CONDITIONING only (accuracy of the winner needs the C++ exact-field demag golden).
"""
import os, json, math
import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from yano_pyr_faces12_star import (
    HEX_FACE_NODES, HEX_EDGES, _tri_area, _face_quadrature, _face_area, _face_area_centroid,
    _trilinear_centroid, _face_normal, _surface_cloud_from_triangles, face_topology,
    star_transform, distorted_patch, _BARY3,
)


def comp_cloud(V, center, layout):
    """element-common cancellation cloud, normalized to unit total charge (sum w = 1)."""
    if layout == "single":
        return center.reshape(1, 3), np.array([1.0])
    if layout == "pyr_faces12":
        tris = [(center, V[a], V[b]) for a, b in HEX_EDGES]
        return _surface_cloud_from_triangles(tris)
    if layout == "faces6":
        tris = []
        for f in range(6):
            q = V[HEX_FACE_NODES[f]]
            tris += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
        return _surface_cloud_from_triangles(tris)
    if layout == "edges_mid":
        pts = np.array([0.5 * (V[a] + V[b]) for a, b in HEX_EDGES])
        return pts, np.full(len(pts), 1.0 / len(pts))
    if layout == "edges_len":
        pts, w = [], []
        for a, b in HEX_EDGES:
            pts.append(0.5 * (V[a] + V[b])); w.append(np.linalg.norm(V[b] - V[a]))
        w = np.array(w); return np.array(pts), w / w.sum()
    if layout == "edges_q":   # 3-pt Gauss along each edge
        g = np.array([0.5 - 0.5 * math.sqrt(0.6), 0.5, 0.5 + 0.5 * math.sqrt(0.6)])
        gw = np.array([5.0, 8.0, 5.0]) / 18.0
        pts, w = [], []
        for a, b in HEX_EDGES:
            L = np.linalg.norm(V[b] - V[a])
            for t, ww in zip(g, gw):
                pts.append((1 - t) * V[a] + t * V[b]); w.append(ww * L)
        w = np.array(w); return np.array(pts), w / w.sum()
    raise ValueError(layout)


def source_cloud(V, face, center, layout):
    fp, fw = _face_quadrature(V, face)
    area = fw.sum()
    cp, cw = comp_cloud(V, center, layout)
    return np.vstack([fp, cp]), np.concatenate([fw, -area * cw])


def collocation_matrix(hexes, layout, eval_alpha=0.5):
    n = len(hexes)
    centers = [_trilinear_centroid(V) for V in hexes]
    eval_pts = np.zeros((6 * n, 3)); normals = np.zeros((6 * n, 3))
    for e, V in enumerate(hexes):
        for f in range(6):
            r = 6 * e + f
            eval_pts[r] = eval_alpha * _face_area_centroid(V, f) + (1 - eval_alpha) * centers[e]
            normals[r] = _face_normal(V, f, centers[e])
    clouds = [source_cloud(hexes[e], f, centers[e], layout) for e in range(n) for f in range(6)]
    neutral = max(abs(w.sum()) for _, w in clouds)
    A = np.zeros((6 * n, 6 * n))
    for col, (sp, sw) in enumerate(clouds):
        diff = eval_pts[:, None, :] - sp[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        m = dist > 1e-14
        inv = np.zeros_like(dist); inv[m] = 1.0 / dist[m] ** 3
        fld = np.sum(sw[None, :, None] * diff * inv[:, :, None], axis=1)
        A[:, col] = -np.einsum("ij,ij->i", normals, fld)
    return A, neutral


def _iters(Ar, rng):
    op = LinearOperator(Ar.shape, matvec=lambda v: Ar @ v)
    out = []
    for b in (rng.random(Ar.shape[0]) - 0.5, np.ones(Ar.shape[0])):
        cnt = {"n": 0}
        bicgstab(op, b, rtol=1e-10, maxiter=20 * Ar.shape[0],
                 callback=lambda xk: cnt.__setitem__("n", cnt["n"] + 1))
        out.append(cnt["n"])
    return out


def study(nx, ny, nz, strength, tag):
    hexes, elements = distorted_patch(nx, ny, nz, strength)
    internal, boundary, adj = face_topology(elements)
    cycle = len(internal) - (len(elements) - 1)
    T = star_transform(hexes, elements, internal, boundary, adj)
    rng = np.random.default_rng(0)
    print(f"\n=== {tag}: {nx}x{ny}x{nz} strength={strength}  n_elem={len(elements)} cycle={cycle} ===")
    print(f"  {'cancellation':>13} | {'neutral':>8} {'near_null':>9} | {'cond_star':>9} | {'bicgstab(rand,unif)':>20}")
    rows = []
    for layout in ("single", "pyr_faces12", "faces6", "edges_mid", "edges_len", "edges_q"):
        A, neutral = collocation_matrix(hexes, layout)
        sv = np.linalg.svd(A, compute_uv=False)
        near_null = int(np.sum(sv < 1e-9 * sv[0]))
        Ar = T.T @ A @ T
        svr = np.linalg.svd(Ar, compute_uv=False)
        cond_star = svr[0] / svr[-1] if svr[-1] > 0 else float("inf")
        it = _iters(Ar, rng)
        print(f"  {layout:>13} | {neutral:>8.1e} {near_null:>9} | {cond_star:>9.1f} | {str(it):>20}")
        rows.append(dict(layout=layout, neutral=float(neutral), near_null=near_null,
                         cond_star=float(cond_star), bicgstab=it))
    best = min(rows, key=lambda r: r["cond_star"])
    print(f"  -> best cond_star: {best['layout']} ({best['cond_star']:.1f})")
    return dict(tag=tag, cycle=cycle, rows=rows, best=best["layout"])


def main():
    out = [study(3, 3, 2, 0.0, "REGULAR"), study(3, 3, 2, 2.0, "distorted"),
           study(5, 5, 4, 4.0, "distorted_big")]
    with open(os.path.join(HERE, "yano_cancellation_study.json"), "w") as f:
        json.dump({"cases": out}, f, indent=2, default=float)
    print("\nNote: monopole field -> CONDITIONING only. The lowest-cond_star cancellation is the")
    print("candidate; its ACCURACY (does it preserve EIEM2 demag at alpha=0.5?) needs the C++ exact field.")
    print("saved", os.path.join(HERE, "yano_cancellation_study.json"))


if __name__ == "__main__":
    main()
