"""yano_loopfree_kernel_study.py -- RESEARCH: can the (evaluation point, central cancellation charge)
make the yano-type collocation matrix INTRINSICALLY loop-free (loops lifted off 0 WITHOUT the star
projection), while reducing to EIEM2 on undistorted elements (accuracy preserved)?

Context: yano_pyr_faces12_star.py showed the loop-free-ness there comes from the STAR PROJECTION
(post-hoc); the kernel (single vs pyramid cloud) only changes the star-block conditioning.  The user
asks the deeper question: design the KERNEL (eval point + cancellation) so the OPERATOR itself has no
spread loop near-null space -- then plain Krylov is loop-free, no star projection needed.

DECISIVE FIRST QUESTION (this script): are the loop modes TOPOLOGICAL near-null (the cell-graph cycles,
intrinsically ~0 field for ANY collocation choice -> kernel tweaks cannot lift them -> star projection
mandatory), OR can some (eval_alpha, cancellation) push the `cycle` smallest singular values UP off 0?

Measured on a REGULAR grid (strength=0; eval_alpha=0.5 + single == EIEM2, the accuracy-preserving
reference) and a distorted grid:
  - near_null(A) at thresholds 1e-9/1e-6/1e-3  (does it stay == cycle for all kernels?)
  - the loop/non-loop singular-value GAP: sv[-(cycle)] (largest loop sv) vs sv[-(cycle+1)] (smallest
    star sv) -- if the loop svs stay ~0 (orders below the star svs) for every (alpha, cancellation),
    the loops are topological and the kernel CANNOT make it intrinsically loop-free.

NOTE: monopole field (this prototype is valid for CONDITIONING/loop structure, NOT accuracy -- exact
analytic field / C++ is needed for the demag factor).  This measures the LOOP STRUCTURE only.
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from yano_pyr_faces12_star import (
    HEX_FACE_NODES, _face_area_centroid, _trilinear_centroid, _face_normal,
    _source_cloud, face_topology, star_transform, distorted_patch, _face_area,
)


def collocation_matrix_alpha(hexes, eval_alpha, layout):
    """yano-type collocation matrix with a PARAMETERIZED eval point:
        eval = eval_alpha * face_area_centroid + (1 - eval_alpha) * trilinear_centroid
    (EIEM2 = 0.5; pyramid = 0.75).  layout = central cancellation cloud ('single' / 'pyr_faces12')."""
    n = len(hexes)
    centers = [_trilinear_centroid(V) for V in hexes]
    eval_pts = np.zeros((6 * n, 3)); normals = np.zeros((6 * n, 3))
    for e, V in enumerate(hexes):
        for f in range(6):
            r = 6 * e + f
            eval_pts[r] = eval_alpha * _face_area_centroid(V, f) + (1 - eval_alpha) * centers[e]
            normals[r] = _face_normal(V, f, centers[e])
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


def study(nx, ny, nz, strength, tag):
    hexes, elements = distorted_patch(nx, ny, nz, strength)
    internal, boundary, adj = face_topology(elements)
    cycle = len(internal) - (len(elements) - 1)
    T = star_transform(hexes, elements, internal, boundary, adj)
    print(f"\n=== {tag}: {nx}x{ny}x{nz} strength={strength}  n_elem={len(elements)} cycle(loops)={cycle} ===")
    print(f"  {'eval_alpha':>10} {'cancel':>12} | {'nn@1e-9':>7} {'nn@1e-6':>7} {'nn@1e-3':>7} | "
          f"{'loop_sv_max':>11} {'star_sv_min':>11} {'gap(star/loop)':>14} | {'cond_star':>9}")
    rows = []
    for layout in ("single", "pyr_faces12"):
        for alpha in (0.5, 0.6, 0.75, 0.9, 1.0):
            A = collocation_matrix_alpha(hexes, alpha, layout)
            sv = np.linalg.svd(A, compute_uv=False)
            s0 = sv[0]
            nn = {t: int(np.sum(sv < t * s0)) for t in (1e-9, 1e-6, 1e-3)}
            loop_sv_max = sv[-cycle] if cycle > 0 else 0.0          # largest of the `cycle` smallest svs
            star_sv_min = sv[-(cycle + 1)] if cycle > 0 and len(sv) > cycle else sv[-1]
            gap = star_sv_min / loop_sv_max if loop_sv_max > 0 else float("inf")
            Ar = T.T @ A @ T
            svr = np.linalg.svd(Ar, compute_uv=False)
            cond_star = svr[0] / svr[-1] if svr[-1] > 0 else float("inf")
            flag = " <-EIEM2" if (abs(alpha - 0.5) < 1e-9 and layout == "single") else ""
            print(f"  {alpha:>10.2f} {layout:>12} | {nn[1e-9]:>7} {nn[1e-6]:>7} {nn[1e-3]:>7} | "
                  f"{loop_sv_max/s0:>11.2e} {star_sv_min/s0:>11.2e} {gap:>14.1f} | {cond_star:>9.1f}{flag}")
            rows.append(dict(alpha=alpha, layout=layout, nn=nn, loop_sv_max=float(loop_sv_max/s0),
                             star_sv_min=float(star_sv_min/s0), gap=float(gap), cond_star=float(cond_star)))
    return dict(tag=tag, cycle=cycle, rows=rows)


def main():
    out = [study(3, 3, 2, 0.0, "REGULAR"), study(3, 3, 2, 2.0, "distorted")]
    # decisive read: does ANY (alpha, cancel) drop nn@1e-3 below cycle (loops lifted off 0)?
    print("\nDECISIVE READ:")
    for o in out:
        lifted = [r for r in o["rows"] if r["nn"][1e-3] < o["cycle"] and r["alpha"] < 1.0]
        if lifted:
            best = min(lifted, key=lambda r: r["nn"][1e-3])
            print(f"  [{o['tag']}] cycle={o['cycle']}: SOME kernel LIFTS loops -> nn@1e-3={best['nn'][1e-3]} "
                  f"(alpha={best['alpha']}, {best['layout']}, loop_sv={best['loop_sv_max']:.1e}) "
                  f"=> intrinsic loop-freeness POSSIBLE")
        else:
            mn = min(o["rows"], key=lambda r: r["loop_sv_max"])
            mx = max(o["rows"], key=lambda r: r["loop_sv_max"])
            print(f"  [{o['tag']}] cycle={o['cycle']}: NO kernel lifts loops (nn@1e-3 == cycle for all); "
                  f"loop_sv stays {mn['loop_sv_max']:.1e}..{mx['loop_sv_max']:.1e} below star -> loops are "
                  f"TOPOLOGICAL => star projection (or HDiv) MANDATORY; kernel tweak alone cannot.")
    with open(os.path.join(HERE, "yano_loopfree_kernel_study.json"), "w") as f:
        json.dump({"cases": out}, f, indent=2, default=float)
    print("\nsaved", os.path.join(HERE, "yano_loopfree_kernel_study.json"))


if __name__ == "__main__":
    main()
