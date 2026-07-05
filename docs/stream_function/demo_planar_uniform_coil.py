"""SF -> planar uniform-Bz coil (the easy end of single-stroke complexity).

Source surface: square z=0 plane of half-extent `plane_half`.
Target: uniform Bz = B0 at z=target_z over a square of half-extent `target_half`.

Pipeline (same shape as demo_sf_to_peec_gx.py, planar geometry):

  [1] stream function method   target Bz=B0 (constant) over (x,y) at z=h
                                -> psi(x,y) via (ACA+)+TSVD pseudo-inverse
  [2] equal-current contours    psi level set -> nested closed curves
                                (for uniform target -> CONCENTRIC, no
                                fingerprint / no Maxwell-pair complexity)
  [3] single-stroke SPIRAL      open each contour at its rightmost point,
                                connect via short straight blend at the cut
                                line.  Topologically trivial (axisymmetric
                                contour family + Kuijpers Method-1 rung) so
                                the chain has NO criss-cross diagonals.
  [4] CAD                       sweep a round wire along the spiral -> STEP
  [5] PEEC                      chain polyline -> PEECBuilder -> L, R at f
  [6] field                     Bz over the target plane + on-axis profile
  [7] verify                    target uniformity vs manufactured-coil Bz

Path A (compensated iteration) is also available here -- the axisymmetric
contour family DOES vary smoothly with psi (no level-set topology jumps
like the Gx fingerprint), so the iteration is expected to actually
converge here.  This validates the principle "simple coil -> easy
single-stroke -> Path A converges" from the
``radia-mcp aca_tsvd(single_stroke)`` discussion.

Run:  python demo_planar_uniform_coil.py [--n 21] [--nlevels 12]
                                          [--compensated-iter N] [--with-peec]
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.stream_function import aca_tsvd, pseudo_inverse_solve
from radia.biot_savart import h_segments_batch, MU0


# --------------------------------------------------------------------------
# Same exact Biot-Savart segment Hz as the Gx demo (consistent kernel for
# matrix entries and verification field).
# --------------------------------------------------------------------------
def _seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1
    R21 = float(d @ d)
    o1 = O - P1
    o2 = O - P2
    OR1 = np.sqrt(o1 @ o1)
    OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d)
    O221 = float(o2 @ d)
    Rc12 = O121 / OR1
    Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f


def _loop_Hz(obs, corners):
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, corners[c], corners[n])
    return hz


def planar_loop_corners(x, y, dx, dy, z=0.0):
    """4 corners of a quad loop on the z-plane, CCW as seen from +z.

    Order (-,-), (+,-), (+,+), (-,+).  CCW from +z gives a loop whose
    magnetic moment points in +z when carrying positive current.
    """
    return np.array([
        [x - dx / 2, y - dy / 2, z],
        [x + dx / 2, y - dy / 2, z],
        [x + dx / 2, y + dy / 2, z],
        [x - dx / 2, y + dy / 2, z],
    ])


# --------------------------------------------------------------------------
# Contour extraction in (x, y).  No phi-periodicity wrap (unlike the
# cylinder case); matplotlib emits each closed contour as one polyline.
# --------------------------------------------------------------------------
def contour_polylines_xy(psi_xy, x_grid, y_grid, n_levels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, hi = float(psi_xy.min()), float(psi_xy.max())
    if hi <= lo:
        return [], 0.0
    dI = (hi - lo) / n_levels
    levels = lo + (np.arange(n_levels) + 0.5) * dI

    X, Y = np.meshgrid(x_grid, y_grid, indexing="ij")
    fig = plt.figure()
    cs = plt.contour(X, Y, psi_xy, levels=levels)
    plt.close(fig)

    polylines = []
    allsegs = getattr(cs, "allsegs", None)
    if allsegs is None:
        return [], dI
    for segs_at_level in allsegs:
        for poly in segs_at_level:
            if len(poly) < 3:
                continue
            polylines.append(np.asarray(poly, dtype=float))
    return polylines, dI


# --------------------------------------------------------------------------
# Single-stroke SPIRAL: open every contour at its rightmost point (max x),
# sort contours by enclosed area (outer first), connect via straight (x, y)
# rung at the cut.  This is Kuijpers Method-1 specialised to an
# axisymmetric contour family on a planar surface -- the cut line is the
# positive x-axis ray.
# --------------------------------------------------------------------------
def _dedupe_path(path, tol=1.0e-9):
    diffs = np.linalg.norm(np.diff(path, axis=0), axis=1)
    keep = np.concatenate([[True], diffs > tol])
    return path[keep]


def single_stroke_spiral_xy(polylines, n_blend=1):
    """Spiral chain: outer-to-inner, opened at rightmost (max-x) point."""
    n = len(polylines)
    if n == 0:
        return np.zeros((0, 2))

    def cut_index(poly):
        return int(np.argmax(poly[:, 0]))

    def cut_radius_sq(poly):
        i = cut_index(poly)
        return float(poly[i, 0] ** 2 + poly[i, 1] ** 2)

    indices = sorted(range(n), key=lambda k: -cut_radius_sq(polylines[k]))

    chain_parts = []
    for k in indices:
        poly = polylines[k]
        i = cut_index(poly)
        closed = np.linalg.norm(poly[0] - poly[-1]) < 1e-9
        if closed:
            body = poly[:-1]
            opened = np.vstack([body[i:], body[:i]])
        else:
            opened = poly if i == 0 else np.vstack([poly[i:], poly[:i]])

        if chain_parts:
            prev_end = chain_parts[-1][-1]
            t = np.linspace(0.0, 1.0, n_blend + 2)[1:-1]
            if len(t) > 0:
                blend = np.column_stack([
                    prev_end[0] + t * (opened[0, 0] - prev_end[0]),
                    prev_end[1] + t * (opened[0, 1] - prev_end[1]),
                ])
                chain_parts.append(blend)
        chain_parts.append(opened)

    return np.vstack(chain_parts) if chain_parts else np.zeros((0, 2))


def chain_xy_to_3d(chain_xy, z=0.0, dedupe_tol=1.0e-9):
    path = np.column_stack([
        chain_xy[:, 0],
        chain_xy[:, 1],
        z * np.ones(len(chain_xy)),
    ])
    return _dedupe_path(path, tol=dedupe_tol)


def bz_at(path_3d, current, obs):
    segs = np.stack([path_3d[:-1], path_3d[1:]], axis=1)
    H = h_segments_batch(segs, obs, current=current)
    return MU0 * H[:, 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=21, help="basis grid n x n")
    ap.add_argument("--plane-half", type=float, default=0.25,
                    help="source plane half-extent [m] (default 250 mm)")
    ap.add_argument("--target-half", type=float, default=0.05,
                    help="target region half-extent [m] (default 50 mm)")
    ap.add_argument("--target-z", type=float, default=0.10,
                    help="target plane z above source [m] (default 100 mm)")
    ap.add_argument("--n-target", type=int, default=5,
                    help="target grid n x n (default 5 -> M=25)")
    ap.add_argument("--B0", type=float, default=1.0e-3,
                    help="target uniform Bz [T] (default 1 mT)")
    ap.add_argument("--nlevels", type=int, default=12,
                    help="equal-current contour levels")
    ap.add_argument("--k", type=int, default=0,
                    help="TSVD modes (0 = auto)")
    ap.add_argument("--compensated-iter", type=int, default=0,
                    help="Path-A compensated iteration count (0 = disable). "
                         "On this axisymmetric topology the iteration is "
                         "expected to actually converge.")
    ap.add_argument("--compensated-step", type=float, default=1.0,
                    help="phi update step size (0 < step <= 1)")
    ap.add_argument("--with-peec", action="store_true",
                    help="also run STEP export + PEEC L, R")
    args = ap.parse_args()

    print("=== SF -> planar uniform-Bz coil ===")

    # ------- basis grid: n x n quad loops on z = 0 plane ----------------
    g = np.linspace(-args.plane_half, args.plane_half, args.n)
    dx = dy = (2 * args.plane_half) / (args.n - 1)
    corners_list = []
    for xi in g:
        for yi in g:
            corners_list.append(planar_loop_corners(xi, yi, dx, dy))
    N = len(corners_list)

    # ------- target 3D points at z = target_z ---------------------------
    gt = np.linspace(-args.target_half, args.target_half, args.n_target)
    Xt, Yt = np.meshgrid(gt, gt, indexing="ij")
    obs = np.column_stack([Xt.ravel(), Yt.ravel(),
                            args.target_z * np.ones(Xt.size)])
    M = obs.shape[0]
    B_target = args.B0 * np.ones(M)

    def entry(i, j):
        return _loop_Hz(obs[i], corners_list[j])

    # ------- [1] SF solve via (ACA+)+TSVD --------------------------------
    modes = min(M, N) if args.k <= 0 else min(args.k, M, N)
    res = aca_tsvd(M, N, entry, modes=modes, kmax=min(M, N),
                   aca_eps=1.0e-8)
    k_use = res.modes if args.k <= 0 else min(args.k, res.modes)
    psi = pseudo_inverse_solve(res, B_target, k_mode=k_use)

    print(f"[1] SF design: plane {args.n}x{args.n} -> N={N} basis loops,"
          f" M={M} target obs at z={args.target_z*1e3:.0f}mm")
    print(f"    ACA+ k_aca = {res.k_aca} (of min(M,N)={min(M, N)}),"
          f" TSVD modes used = {k_use}")

    # Continuous residual = how well the SF design hits B_target
    A_psi = np.array([sum(entry(i, j) * psi[j] for j in range(N))
                      for i in range(M)])
    cont_res = float(np.linalg.norm(A_psi - B_target) /
                     (np.linalg.norm(B_target) + 1e-30))
    print(f"    continuous SF: ||A psi - B0||/||B0|| = {cont_res:.3e}")

    psi_xy = psi.reshape(args.n, args.n)

    # ------- [2] equal-current contours ---------------------------------
    polylines, dI = contour_polylines_xy(psi_xy, g, g, args.nlevels)
    n_wire = len(polylines)
    if n_wire == 0:
        print("[2] no usable contours (psi is constant). aborting.")
        return
    print(f"[2] contours: {n_wire} closed curves, dI = {dI:.4g}")

    # ------- [3a] Path-A compensated iteration (optional) ---------------
    if args.compensated_iter > 0:
        print(f"[3a] Path-A compensated iteration "
              f"({args.compensated_iter} iters, step={args.compensated_step})")
        best_psi = psi.copy()
        best_polylines = polylines
        best_res_norm = float("inf")
        for it in range(args.compensated_iter):
            chain_it = single_stroke_spiral_xy(polylines)
            path_it = chain_xy_to_3d(chain_it)
            Bz_unit_it = bz_at(path_it, 1.0, obs)
            denom_it = float(np.dot(Bz_unit_it, Bz_unit_it))
            I_w_it = float(np.dot(Bz_unit_it, B_target) / denom_it) \
                if denom_it > 0 else 0.0
            residual = B_target - I_w_it * Bz_unit_it
            res_norm = float(np.linalg.norm(residual) /
                             (np.linalg.norm(B_target) + 1e-30))
            tag = ""
            if res_norm < best_res_norm:
                best_res_norm = res_norm
                best_psi = psi.copy()
                best_polylines = polylines
                tag = " <-- best"
            print(f"     iter {it+1}: I_w = {I_w_it:.3e},"
                  f" |B_target - chain|/|B_target| = {res_norm:.3e}{tag}")
            delta_phi = pseudo_inverse_solve(res, residual, k_mode=k_use)
            psi = psi + args.compensated_step * delta_phi
            psi_xy = psi.reshape(args.n, args.n)
            polylines, dI = contour_polylines_xy(
                psi_xy, g, g, args.nlevels)
            if not polylines:
                print("     iter aborted: psi became flat")
                break
        psi = best_psi
        polylines = best_polylines
        psi_xy = psi.reshape(args.n, args.n)
        print(f"     final: best residual = {best_res_norm:.3e}")

    # ------- [3] single-stroke spiral chain ------------------------------
    chain_xy = single_stroke_spiral_xy(polylines)
    path = chain_xy_to_3d(chain_xy)
    seg_lens = np.linalg.norm(np.diff(path, axis=0), axis=1)
    length = float(seg_lens.sum())
    print(f"[3] spiral chain: {len(path)} points,"
          f" {len(path) - 1} segments, wire length = {length:.3f} m"
          f" (one continuous conductor)")

    # ------- [6][7] field of single-stroke coil over the target ---------
    Bz_unit = bz_at(path, 1.0, obs)
    denom = float(np.dot(Bz_unit, Bz_unit))
    I_w = float(np.dot(Bz_unit, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit
    rms = float(np.linalg.norm(Bz_dsv - B_target) /
                (np.linalg.norm(B_target) + 1e-30))
    Bz_mean = float(np.mean(Bz_dsv))
    Bz_max = float(np.max(Bz_dsv))
    Bz_min = float(np.min(Bz_dsv))
    nonuniformity = (Bz_max - Bz_min) / (abs(Bz_mean) + 1e-30)
    print(f"[6] best-fit single current I_w = {I_w:.4g} A,"
          f" target-plane RMS = {rms:.3e}")
    print(f"    Bz over target: mean = {Bz_mean:.3e} T,"
          f" range = [{Bz_min:.3e}, {Bz_max:.3e}] T,"
          f" peak-to-peak / mean = {nonuniformity:.3e}")
    print(f"    (target B0 = {args.B0:.3e} T)")

    # On-axis Bz vs z to visualise coil's region of uniformity
    zv = np.linspace(0.5 * args.target_z, 1.5 * args.target_z, 41)
    obs_axis = np.column_stack([np.zeros_like(zv), np.zeros_like(zv), zv])
    Bz_axis = bz_at(path, I_w, obs_axis)
    print(f"[7] on-axis Bz at z={args.target_z*1e3:.0f}mm: {Bz_axis[20]:.3e} T")

    if args.with_peec:
        run_peec_chain(path, I_w)

    # ------- plot -------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig = plt.figure(figsize=(13, 4))

        ax1 = fig.add_subplot(1, 3, 1)
        cf = ax1.contourf(g * 1e3, g * 1e3, psi_xy.T, levels=20, cmap="RdBu_r")
        lo, hi = psi_xy.min(), psi_xy.max()
        if hi > lo:
            lv = lo + (np.arange(args.nlevels) + 0.5) * (hi - lo) / args.nlevels
            ax1.contour(g * 1e3, g * 1e3, psi_xy.T, levels=lv,
                        colors="k", linewidths=0.5)
        ax1.set_xlabel("$x$ (mm)")
        ax1.set_ylabel("$y$ (mm)")
        ax1.set_aspect("equal")
        ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=10)
        fig.colorbar(cf, ax=ax1, fraction=0.045)

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, "b-", lw=0.5)
        ax2.set_xlabel("$x$ (mm)"); ax2.set_ylabel("$y$ (mm)")
        ax2.set_aspect("equal")
        ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
        ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.plot(zv * 1e3, Bz_axis * 1e3, "o-", ms=3, label="single-stroke coil")
        ax3.axhline(args.B0 * 1e3, color="k", linestyle="--",
                    label=f"target {args.B0*1e3:.2f} mT")
        ax3.axvline(args.target_z * 1e3, color="g", linestyle=":",
                    alpha=0.5, label=f"target z={args.target_z*1e3:.0f}mm")
        ax3.set_xlabel("$z$ (mm)"); ax3.set_ylabel("$B_z$ (mT)")
        ax3.text(0.02, 0.96, "(c)", transform=ax3.transAxes, fontsize=10)
        ax3.legend(loc="best", fontsize=8); ax3.grid(alpha=0.3)

        out = HERE / "demo_planar_uniform_coil.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


def run_peec_chain(path, I_w):
    """[4] CAD STEP via build123d Spline + Frenet sweep + [5] PEEC L, R."""
    stride = max(1, len(path) // 500)
    cpath = path[::stride]
    if not np.allclose(cpath[-1], path[-1]):
        cpath = np.vstack([cpath, path[-1]])
    diffs = np.linalg.norm(np.diff(cpath, axis=0), axis=1)
    cpath = cpath[np.concatenate([[True], diffs > 1.0e-4])]

    step_path = str(HERE / "planar_uniform_coil.step")
    try:
        import build123d as b3d
        from radia.coil_from_cad import export_step_with_labels
        mm = 1000.0
        pmm = cpath * mm
        seg_lens = np.linalg.norm(np.diff(pmm, axis=0), axis=1)
        r_wire = float(max(0.2, min(1.0, 0.25 * np.median(seg_lens))))
        t0 = time.time()
        try:
            wire = b3d.Spline(*[b3d.Vector(float(x), float(y), float(z))
                                for x, y, z in pmm])
            # tangent at first segment
            tan = pmm[1] - pmm[0]
            tan = tan / np.linalg.norm(tan)
            pl = b3d.Plane(origin=b3d.Vector(*pmm[0]), z_dir=b3d.Vector(*tan))
            sec = pl.location * b3d.Circle(radius=r_wire).face()
            solid = b3d.sweep(sec, path=wire, is_frenet=True)
            try:
                solid.label = "coil"
            except Exception:
                pass
            export_step_with_labels([solid], step_path)
            print(f"[4] CAD: round wire r={r_wire:.2f}mm swept along spiral -> "
                  f"{os.path.basename(step_path)} ({os.path.getsize(step_path)//1024} KB,"
                  f" vol={solid.volume:.0f} mm^3, {time.time()-t0:.1f}s)")
        except Exception as se:
            poly = b3d.Polyline(*[b3d.Vector(float(x), float(y), float(z))
                                  for x, y, z in pmm])
            b3d.export_step(poly, step_path)
            print(f"[4] CAD: solid sweep failed ({type(se).__name__}); "
                  f"exported Polyline WIRE -> {os.path.basename(step_path)}"
                  f" ({os.path.getsize(step_path)//1024} KB)")
    except Exception as e:
        print(f"[4] CAD STEP export FAILED: {type(e).__name__}: {str(e)[:120]}")

    try:
        from radia.peec_matrices import PEECBuilder
        from radia.peec_topology import PEECCircuitSolver
        b = PEECBuilder()
        nodes = [b.add_node_at(float(p[0]), float(p[1]), float(p[2]))
                 for p in cpath]
        seg_lens = np.linalg.norm(np.diff(cpath, axis=0), axis=1)
        w = hh = float(max(2.0e-4, min(2.0e-3, 0.2 * np.median(seg_lens))))
        for i in range(len(nodes) - 1):
            b.add_connected_segment(nodes[i], nodes[i + 1], w, hh, sigma=5.8e7)
        b.add_port(nodes[0], nodes[-1])
        topo = b.build_topology()
        solver = PEECCircuitSolver(topo)
        freq = 1.0e3
        t0 = time.time()
        Z = solver.compute_port_impedance(freq)
        Zc = complex(np.asarray(Z).ravel()[0]) if np.ndim(Z) else complex(Z)
        L = Zc.imag / (2 * np.pi * freq)
        R = Zc.real
        length = float(np.sum(np.linalg.norm(np.diff(cpath, axis=0), axis=1)))
        print(f"[5] PEEC ({len(nodes) - 1} filament segs, {length:.2f} m @ {freq:.0f} Hz,"
              f" {time.time()-t0:.1f}s): L = {L * 1e6:.3f} uH, R = {R * 1e3:.3f} mOhm")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[5] PEEC FAILED: {type(e).__name__}: {str(e)[:120]}")


if __name__ == "__main__":
    main()
