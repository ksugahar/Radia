#!/usr/bin/env python3
"""mini-REGCOIL / NESCOIL: a fusion coil-design (stellarator Stage-2) demo.

The stream-function (current-potential) machinery that designs an MRI gradient
coil or an induction-heating work coil is, mathematically, the SAME object the
fusion community calls a **winding-surface current potential** -- the unknown of
NESCOIL (Merkel 1987), REGCOIL (Landreman 2017) and the surface part of FOCUS
(Zhu 2018).  The Stage-2 coil problem is:

    given a target NORMAL field  B.n  on the PLASMA boundary,
    find a current potential psi on a WINDING SURFACE around it
    whose Biot-Savart field reproduces that B.n  (its iso-contours = the coils).

This demo runs that forward problem with Radia's existing surface-FE stream
function.  ``A_n psi = B_n`` where the rows of ``A_n`` are the plasma-boundary
normal component ``n.B`` of the winding-surface Biot-Savart kernel (we already
assemble the 3-component ``A``; we just dot it with the plasma normal).  Two
demonstrations:

  1. FORWARD WORKS -- two producible targets reproduced to MACHINE PRECISION:
       * uniform vertical field B.n   (a PF / equilibrium / vertical-field coil)
       * non-axisymmetric  sin(theta) cos(2 phi)  (a stellarator-like shape)
     Both give ||A_n psi - B_n|| / ||B_n|| ~ 1e-8 -- the SF designer does the
     stellarator Stage-2 forward problem exactly.

  2. THE REGCOIL TRADE-OFF -- on a GENUINELY HARD target the winding surface
     cannot reproduce cheaply (a high (theta, phi) mode that decays across the
     plasma-coil gap), sweeping the regularisation weight alpha traces the
     classic REGCOIL L-curve: large alpha -> smooth coil, high B.n residual;
     small alpha -> low residual, high peak current density |grad psi| (the
     coil gets sharp / hard to wind).  (field error, coil complexity) is the
     fundamental Stage-2 trade-off REGCOIL is built to expose.

HONEST SCOPE (what this is and is NOT):
  * This is the FORWARD problem with a SINGLE-VALUED current potential psi --
    correct for PF / RMP / shaping / shim fields (no net poloidal current
    through the winding torus hole).  A net-current (TF-type) coil needs the
    MULTI-VALUED secular term  G*zeta + I*theta ; that is the cohomology
    generator the lab already has in ``cohomology_cut.py`` (GMSH cohomology),
    not wired into this demo.
  * The producible targets reproduce to ~1e-8 precisely BECAUSE they are smooth
    and within the winding surface's reach; that is honest (the forward map is
    exact), not a hidden tolerance.  The REGCOIL residual trade-off only appears
    for a target the surface canNOT exactly produce (panel 2).
  * A production stellarator run additionally needs: B.n from a free-boundary
    VMEC equilibrium (here it is an analytic model field), the winding-surface
    geometry optimisation (FOCUS), and coil force / stress.  Those are the named
    next steps, not part of this PoC.

Run (standalone -- no Cubit, no panel UI):
    python demo_regcoil_fusion.py
    python demo_regcoil_fusion.py --no-plot          # JSON only
    python demo_regcoil_fusion.py --n-alpha 12       # finer REGCOIL L-curve

Writes ``demo_regcoil_fusion.json`` (Data Persistence Policy) and, if matplotlib
+ radia-mcp are present, ``demo_regcoil_fusion.{png,pdf}`` next to this script.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from datetime import datetime

import numpy as np

# major radius shared by the winding torus and the plasma torus (a simple,
# reproducible stand-in for a VMEC boundary -- swap pts/nrm for a real one).
R_MAJOR = 0.30
A_PLASMA = 0.06       # plasma-boundary minor radius
A_WIND = 0.12         # winding-surface minor radius (the coil lives here)


def _src_paths():
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "..", "..", "src", "radia"))
    sys.path.insert(0, os.path.join(here, "..", "..", "src", "radia", "panels"))


def _torus_surface_vol(a, maxh, path):
    """Save a torus SURFACE mesh (revolved wire = dim-2 surface, what the
    surface stream function lives on -- revolving a Face would give a solid)."""
    from netgen.occ import WorkPlane, Axes, Pnt, Y, Z, Axis, OCCGeometry
    wire = WorkPlane(Axes(Pnt(R_MAJOR, 0, 0), n=Y)).Circle(a).Wire()
    surf = wire.Revolve(Axis(Pnt(0, 0, 0), Z), 360)
    OCCGeometry(surf).GenerateMesh(maxh=maxh).Save(path)
    return path


def _plasma_points_normals(plasma, eval_max):
    """Plasma-boundary sample points + analytic outward normals (torus)."""
    pts = np.array([list(v.point) for v in plasma.vertices])
    if len(pts) > eval_max:
        pts = pts[np.linspace(0, len(pts) - 1, eval_max).astype(int)]
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    nrm = np.column_stack([pts[:, 0] - R_MAJOR * np.cos(phi),
                           pts[:, 1] - R_MAJOR * np.sin(phi), pts[:, 2]])
    nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-30)
    # poloidal angle theta about the tube centre (for shaped targets)
    rho = np.hypot(pts[:, 0], pts[:, 1]) - R_MAJOR
    theta = np.arctan2(pts[:, 2], rho)
    return pts, nrm, theta, phi


def _design(C, fes, coil, A_n, Bn, regularize, alpha_rel):
    """Solve min psi^T S psi s.t. (regularised) A_n psi = Bn; return psi + diag.

    Returns (psi, B.n relative residual, peak |grad psi| = coil complexity)."""
    from ngsolve import GridFunction
    S = C._seminorm(fes, regularize).toarray()
    AtA = A_n.T @ A_n
    alpha = alpha_rel * np.trace(AtA) / fes.ndof
    psi = np.linalg.solve(AtA + alpha * S, A_n.T @ Bn)
    res = float(np.linalg.norm(A_n @ psi - Bn) / (np.linalg.norm(Bn) + 1e-30))
    gfu = GridFunction(fes)
    gfu.vec.FV().NumPy()[:] = psi
    return psi, res, float(C._peak_current_density(fes, coil, gfu))


def _contours(C, coil, psi, n_levels):
    verts = np.array([list(v.point) for v in coil.vertices])
    tris = C._surface_triangles(coil)
    tol = max(1e-12, 1e-6 * C._char_len(verts, tris))
    loops = []
    for lev in C._contour_levels(psi[:coil.nv], n_levels):
        loops.extend(C._segs_to_polylines(
            C._contour_segments(verts, psi[:coil.nv], tris, lev), tol))
    return loops


def main():
    ap = argparse.ArgumentParser(
        description="mini-REGCOIL / NESCOIL fusion coil-design demo")
    ap.add_argument("--out-dir",
                    default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--work-dir",
                    default=os.path.join(os.environ.get("TEMP", "/tmp"),
                                         "sf_regcoil"))
    ap.add_argument("--eval-max", type=int, default=220,
                    help="plasma-boundary B.n sample points")
    ap.add_argument("--wind-maxh", type=float, default=0.045)
    ap.add_argument("--plasma-maxh", type=float, default=0.035)
    ap.add_argument("--n-levels", type=int, default=14,
                    help="current-potential contours (= coil turns) to draw")
    ap.add_argument("--n-alpha", type=int, default=9,
                    help="REGCOIL L-curve points (alpha sweep)")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    _src_paths()
    import calc_streamfunction as C
    from ngsolve import Mesh, H1, specialcf, TaskManager

    os.makedirs(args.work_dir, exist_ok=True)

    with TaskManager():
        coil = Mesh(_torus_surface_vol(
            A_WIND, args.wind_maxh, os.path.join(args.work_dir, "wind.vol")))
        plasma = Mesh(_torus_surface_vol(
            A_PLASMA, args.plasma_maxh, os.path.join(args.work_dir, "plasma.vol")))
        fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        n_cf = specialcf.normal(3)

        pts, nrm, theta, phi = _plasma_points_normals(plasma, args.eval_max)
        # A (3-component) at the plasma points -> dot with plasma normal = A_n.
        A3 = C._assemble_biot_savart(fes, n_cf, pts, [0, 1, 2]).reshape(
            len(pts), 3, fes.ndof)
        A_n = np.einsum("mc,mcj->mj", nrm, A3)

        # ---- 1. forward works: two producible targets -> machine precision ----
        targets = {
            "uniform_vertical_PF": nrm[:, 2],
            "stellarator_sinTheta_cos2phi": np.sin(theta) * np.cos(2 * phi),
        }
        designs = {}
        ref_loops = None
        for name, Bn in targets.items():
            psi, res, peakJ = _design(C, fes, coil, A_n, Bn, "h1", 1e-7)
            loops = _contours(C, coil, psi, args.n_levels)
            designs[name] = {"bn_residual_rel": res,
                             "peak_grad_psi": peakJ,
                             "n_contours": len(loops)}
            if name == "uniform_vertical_PF":
                ref_loops = loops
            print(f"[forward] {name:32s}  B.n resid={res:.2e}  "
                  f"peak|grad psi|={peakJ:.2e}  contours={len(loops)}")

        # ---- 2. REGCOIL L-curve on a GENUINELY HARD target ----
        # a high (theta, phi) mode decays across the plasma-coil gap, so it is
        # NOT cheaply producible: alpha now genuinely trades B.n residual against
        # peak |grad psi| (coil complexity).  alpha large -> smooth coil, high
        # residual; alpha small -> low residual, current density saturates at the
        # surface's representation limit (the vertical leg of the L-curve).
        Bn_hard = np.sin(3 * theta) * np.cos(5 * phi)
        alpha_rel_grid = np.logspace(2.0, -8.0, args.n_alpha)
        lcurve = []
        for a_rel in alpha_rel_grid:
            _psi, res, peakJ = _design(C, fes, coil, A_n, Bn_hard, "h1",
                                       float(a_rel))
            lcurve.append({"alpha_rel": float(a_rel),
                           "bn_residual_rel": res, "peak_grad_psi": peakJ})
            print(f"[L-curve] alpha_rel={a_rel:.1e}  B.n resid={res:.3e}  "
                  f"peak|grad psi|={peakJ:.3e}")

    data = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "benchmark": "regcoil_fusion_forward",
        "problem": {
            "winding_torus": {"R": R_MAJOR, "a": A_WIND, "maxh": args.wind_maxh},
            "plasma_torus": {"R": R_MAJOR, "a": A_PLASMA, "maxh": args.plasma_maxh},
            "winding_ndof": int(fes.ndof),
            "plasma_eval_points": int(len(pts)),
            "regularize": "h1",
            "note": "single-valued psi (no net-current secular term); "
                    "B.n from analytic model fields (swap for free-boundary "
                    "VMEC for a production run).",
        },
        "producible_targets": designs,
        "regcoil_lcurve": {
            "hard_target": "sin(3 theta) cos(5 phi)",
            "points": lcurve,
        },
    }
    jpath = os.path.join(args.out_dir, "demo_regcoil_fusion.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"\nResults -> {jpath}")

    if not args.no_plot:
        _plot(args, ref_loops, lcurve)


def _plot(args, loops, lcurve):
    """(a) 3D current-potential coil around the plasma boundary;
       (b) the REGCOIL L-curve (B.n residual vs peak current density)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Lab figure convention (radia_mcp.figure): authored AT the 16 cm embed
    # width, 10 pt, Times New Roman, TrueType; NO in-figure title (the README /
    # caption carries it).  apply_lab_style sets the fonts; a 3D panel is not
    # served by the 2D-grid lab_figure() helper, so we build the mixed
    # 3D + 2D figure ourselves and let save_lab_figure() gate + write it.
    try:
        from radia_mcp.figure import apply_lab_style, save_lab_figure
    except Exception as e:
        print(f"(plot skipped: radia-mcp figure unavailable: {e})")
        return
    plt.rcParams["pdf.fonttype"] = 42
    w_in, h_in = apply_lab_style(embed_width_cm=16.0, aspect=0.46)
    fig = plt.figure(figsize=(w_in, h_in))

    # (a) the designed coil = current-potential contours on the winding torus
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, 2 * np.pi, 40)
    U, V = np.meshgrid(u, v)
    X = (R_MAJOR + A_PLASMA * np.cos(V)) * np.cos(U)
    Yy = (R_MAJOR + A_PLASMA * np.cos(V)) * np.sin(U)
    Zz = A_PLASMA * np.sin(V)
    ax.plot_surface(X, Yy, Zz, color="orange", alpha=0.18, linewidth=0)
    for p in loops:
        ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.8, color="C0")
    ax.text2D(0.02, 0.97, "current-potential coil", transform=ax.transAxes)
    ax.text2D(0.02, 0.90, "plasma boundary", transform=ax.transAxes,
              color="darkorange")
    for s in ("x", "y", "z"):
        getattr(ax, f"set_{s}ticks")([])
    ax.set_box_aspect((1, 1, 0.5))

    # (b) the REGCOIL L-curve: residual vs coil complexity as alpha sweeps
    ax2 = fig.add_subplot(1, 2, 2)
    res = [pt["bn_residual_rel"] for pt in lcurve]
    pk = [pt["peak_grad_psi"] for pt in lcurve]
    ax2.loglog(pk, res, "o-", color="C3")
    ax2.set_xlabel(r"peak $|\nabla\psi|$  (coil complexity)")
    ax2.set_ylabel(r"$B{\cdot}n$ residual  (rel.)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.annotate(r"large $\alpha$", xy=(pk[0], res[0]),
                 xytext=(0.30, 0.85), textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=0.7))
    ax2.annotate(r"small $\alpha$", xy=(pk[-1], res[-1]),
                 xytext=(0.55, 0.18), textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=0.7))

    info = save_lab_figure(fig, os.path.join(args.out_dir,
                                             "demo_regcoil_fusion"),
                           embed_width_cm=16.0)
    print(f"Plot    -> {', '.join(info['wrote'])}")


if __name__ == "__main__":
    main()
