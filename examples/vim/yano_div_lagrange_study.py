"""yano_div_lagrange_study.py -- RESEARCH (user idea): impose the PHYSICAL div(B)=0 constraint
(sum of an element's 6 face charges = 0, i.e. magnetic charge neutrality / no monopole) via LAGRANGE
multipliers, instead of via the CANCELLATION CHARGE that the yano-MSC kernel currently uses.

Context: yano_cancellation_study.py showed the cancellation cloud (single/pyr12/edges) sets the
conditioning; pyr12 is best.  But the cancellation is an AD-HOC per-DOF neutralizer.  The user asks:
the physical statement is div(B)=0 -> sum_f (sigma_f * area_f) = 0 per element.  Impose THAT directly
via Lagrange on the RAW (un-cancelled) collocation, [[A_raw, C^T],[C,0]][sigma;lambda]=[b;0] with
C[e, dof(e,f)] = area_f.  Does it (1) match pyr12's conditioning, (2) remove the monopole modes,
(3) relate to the loops (cell-graph cycles)?  div=0 (n_el constraints) != loops (cycle count) -- they
are different subspaces (physics-Gauss vs field-null), so this measures whether div=0-Lagrange is a
viable PHYSICAL replacement for the cancellation charge.

Monopole field -> CONDITIONING/structure only (the prototype's valid regime).
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from yano_pyr_faces12_star import (
    _face_quadrature, _face_area, _face_area_centroid, _trilinear_centroid, _face_normal,
    _source_cloud, face_topology, star_transform, distorted_patch,
)


def collocation(hexes, layout, eval_alpha=0.5):
    """yano-type collocation at EIEM2 eval point; layout='raw' (face charge only, NO cancellation),
    'single', or 'pyr_faces12'."""
    n = len(hexes)
    centers = [_trilinear_centroid(V) for V in hexes]
    eval_pts = np.zeros((6 * n, 3)); normals = np.zeros((6 * n, 3))
    for e, V in enumerate(hexes):
        for f in range(6):
            r = 6 * e + f
            eval_pts[r] = eval_alpha * _face_area_centroid(V, f) + (1 - eval_alpha) * centers[e]
            normals[r] = _face_normal(V, f, centers[e])
    if layout == "raw":
        clouds = [(_face_quadrature(hexes[e], f)) for e in range(n) for f in range(6)]  # (pts, w) face only
    else:
        clouds = [_source_cloud(hexes[e], f, centers[e], layout) for e in range(n) for f in range(6)]
    A = np.zeros((6 * n, 6 * n))
    for col, (sp, sw) in enumerate(clouds):
        diff = eval_pts[:, None, :] - sp[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        m = dist > 1e-14
        inv = np.zeros_like(dist); inv[m] = 1.0 / dist[m] ** 3
        fld = np.sum(sw[None, :, None] * diff * inv[:, :, None], axis=1)
        A[:, col] = -np.einsum("ij,ij->i", normals, fld)
    return A


def div_constraint(hexes):
    """C[e, 6e+f] = area_f : sum of the element's 6 face TOTAL charges (div B = 0 / neutrality)."""
    n = len(hexes)
    C = np.zeros((n, 6 * n))
    for e, V in enumerate(hexes):
        for f in range(6):
            C[e, 6 * e + f] = _face_area(V, f)
    return C


def _near_null(M):
    sv = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(sv < 1e-9 * sv[0]))


def _ker(C):
    """orthonormal basis of ker(C) (the div-free subspace)."""
    _u, s, vt = np.linalg.svd(C, full_matrices=True)
    rank = int(np.sum(s > 1e-9 * s[0]))
    return vt[rank:, :].T, rank


def _star_cond(N, mu_r):
    """system cond after removing the loops (near-null of N): (1/chi)I - N on the non-loop subspace."""
    chi = mu_r - 1.0
    sv = np.linalg.svd(N, compute_uv=False)
    nl = int(np.sum(sv < 1e-9 * sv[0]))
    _u, _s, vt = np.linalg.svd(N)
    P = vt[:N.shape[0] - nl, :].T            # complement of the near-null (loops)
    A = (1.0 / chi) * np.eye(N.shape[0]) - N
    Ar = P.T @ A @ P
    return float(np.linalg.cond(Ar)), nl


def main():
    out = []
    for (nx, ny, nz, strength) in [(3, 3, 2, 0.0), (3, 3, 2, 2.0)]:
        hexes, elements = distorted_patch(nx, ny, nz, strength)
        internal, boundary, adj = face_topology(elements)
        cycle = len(internal) - (len(elements) - 1)
        C = div_constraint(hexes)
        Q, rankC = _ker(C)                    # div-free subspace basis (ndof x (ndof-rankC))
        n_el = len(elements)
        tag = "REGULAR" if strength == 0 else "distorted"
        print(f"\n=== {tag} {nx}x{ny}x{nz}: n_el={n_el}, ndof={6*n_el}, rank(div C)={rankC}, "
              f"loops(cycle)={cycle} ===")

        N_raw = collocation(hexes, "raw")
        N_pyr = collocation(hexes, "pyr_faces12")
        nn_raw = _near_null(N_raw)
        nn_pyr = _near_null(N_pyr)
        # div=0: restrict raw N to the div-free subspace -> does it drop the monopole modes?
        N_raw_df = Q.T @ N_raw @ Q
        nn_raw_df = _near_null(N_raw_df)
        print(f"  near_null:  raw={nn_raw} (monopole+loops)   raw|div-free={nn_raw_df}   "
              f"pyr12={nn_pyr} (loops)")
        # system conditioning at mu_r=1e4 after loop removal (the usable metric)
        c_pyr, l_pyr = _star_cond(N_pyr, 1e4)
        # raw|div-free, then remove loops too: build the system in the div-free subspace
        chi = 1e4 - 1.0
        svdf = np.linalg.svd(N_raw_df, compute_uv=False)
        nldf = int(np.sum(svdf < 1e-9 * svdf[0]))
        _u, _s, vtdf = np.linalg.svd(N_raw_df)
        Pdf = vtdf[:N_raw_df.shape[0] - nldf, :].T
        Adf = (1.0 / chi) * np.eye(N_raw_df.shape[0]) - N_raw_df
        c_rawdf = float(np.linalg.cond(Pdf.T @ Adf @ Pdf))
        print(f"  system cond @mu_r=1e4 (loops removed):  pyr12={c_pyr:.1f} (loops={l_pyr})   "
              f"raw|div-free|loop-free={c_rawdf:.1f} (loops_df={nldf})")
        verdict = ("div=0 removes the monopole modes (raw|div-free near_null == pyr12) -> physical "
                   "Gauss constraint REPLACES the cancellation charge"
                   if nn_raw_df == nn_pyr else
                   "div=0 does NOT reduce to the cancellation near_null -- different effect")
        print(f"  -> {verdict}")
        out.append(dict(tag=tag, n_el=n_el, ndof=6 * n_el, cycle=cycle, rankC=rankC,
                        nn_raw=nn_raw, nn_raw_divfree=nn_raw_df, nn_pyr=nn_pyr,
                        cond_pyr12=c_pyr, cond_raw_divfree_loopfree=c_rawdf))
    with open(os.path.join(HERE, "yano_div_lagrange_study.json"), "w") as f:
        json.dump({"cases": out}, f, indent=2, default=float)
    print("\nReadout: div=0 (Gauss, n_el constraints) removes the MONOPOLE near-null modes (the")
    print("cancellation charge's job); the LOOPS (cell-graph cycles) are a SEPARATE near-null that still")
    print("needs the star/loop-Lagrange. If raw|div-free|loop-free system cond ~ pyr12, the physical")
    print("div=0-Lagrange + loop-Lagrange REPLACES the ad-hoc cancellation cloud.")
    print("saved", os.path.join(HERE, "yano_div_lagrange_study.json"))


if __name__ == "__main__":
    main()
