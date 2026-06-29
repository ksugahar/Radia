#!/usr/bin/env python3
"""Shim-coil harmonic PURITY demo (spherical-harmonic target + decomposition).

Designs SF surface coils on ONE cylindrical former for a sequence of named
spherical-harmonic targets -- a gradient (Gz, Gx), three 2nd-order shims
(Z2, X2-Y2, ZX) and a 4th-order shim (Z4) -- and reports, for each, how PURE
the achieved field is in its target harmonic and which harmonic is the largest
residual CONTAMINANT.  These are the canonical MRI gradient/shim quality
metrics.

It runs the SHIPPED panel calc end-to-end, exactly as the GUI does:

    former (cylinder lateral surface, OCC) + DSV (sphere volume)
      ->  for each target T in {Z, X, Z2, C2, ZX, Z4}:
            python calc_streamfunction.py --target-harmonic T --harmonic-lmax 4
      ->  JSON ``harmonics``: purity (target-harmonic field fraction over the
          DSV) + the named max contaminant + the per-(l,m) spectrum.

Why it is a good demonstration
------------------------------
The SF designer fits every target to near-unit purity (l<=3: residual
<0.1 %), and the SAME solid-harmonic table both BUILDS the target and ANALYSES
the result, so the round-trip is exact.  The difficulty grows with the
harmonic order: over a fixed small DSV the field of an l-th harmonic scales as
r^l, so the high-order Z4 shim is the hardest -- yet it still reaches ~0.9998
purity, with the decomposition naming its ~1 % Z3 contaminant.  This is the
l=4 extension of the harmonic basis in action.

Mesh generation runs in a SUBPROCESS (re-entrant ``--make-mesh``) so this
orchestrating process never imports NGSolve; each design also runs in its own
calc subprocess.  Standalone:

    python demo_shim_coil_purity.py [--order 2] [--lmax 4]

Outputs (committed next to this script, per the Data Persistence Policy):
    demo_shim_coil_purity.json   aggregated purity/contamination per target
    demo_shim_coil_purity.png    impurity-vs-order + hardest-target spectrum
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CALC = REPO / "src" / "radia" / "panels" / "calc_streamfunction.py"

# Showcase targets: (name, harmonic order l, human label).  One gradient per
# transverse/axial direction, the three independent 2nd-order shims, and a
# 4th-order shim to exercise the l=4 basis extension.
TARGETS = [
    ("Z",  1, "Gz gradient"),
    ("X",  1, "Gx gradient"),
    ("Z2", 2, "Z2 shim"),
    ("C2", 2, "X2-Y2 shim"),
    ("ZX", 2, "ZX shim"),
    ("Z4", 4, "Z4 shim"),
]
# colour each bar by its harmonic order l (l=1 / l=2 / l=4)
_L_COLOR = {1: "#2c7fb8", 2: "#41ab5d", 3: "#d95f0e", 4: "#993404"}


def _make_meshes(outdir):
    """Re-entrant mesh-gen subprocess target: cylinder lateral former + DSV
    sphere volume (OCC, Cubit-free).  Isolated in a subprocess so the
    orchestrating ``main()`` -- and any import of this module -- stays
    NGSolve-free."""
    from netgen.occ import Cylinder, Sphere, Pnt, Z, OCCGeometry
    from ngsolve import TaskManager
    os.makedirs(outdir, exist_ok=True)
    with TaskManager():            # caller wraps mesh-gen (TaskManager-Only)
        a, L = 0.15, 0.50
        cyl = Cylinder(Pnt(0, 0, -L / 2), Z, r=a, h=L)
        lateral = max(cyl.faces, key=lambda f: f.mass)   # 2*pi*a*L > cap area
        OCCGeometry(lateral).GenerateMesh(maxh=0.04).Save(
            os.path.join(outdir, "former.vol"))
        OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.025).Save(
            os.path.join(outdir, "dsv.vol"))


def _run_target(coil, evalv, target, lmax, order, workdir):
    """Design one harmonic target via the shipped calc; return its JSON dict."""
    out = os.path.join(workdir, f"design_{target}.json")
    cmd = [sys.executable, str(CALC), "--coil-vol", coil, "--eval-vol", evalv,
           "--order", str(order), "--confine", "abe",
           "--target-harmonic", target, "--harmonic-lmax", str(lmax),
           "--output", out]
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env,
                       timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f"calc_streamfunction.py failed for {target} "
                           f"(rc={r.returncode}):\n{r.stderr[-1200:]}")
    with open(out) as f:
        return json.load(f)


def _design_all(coil, evalv, lmax, order, workdir):
    rows = []
    for target, l, label in TARGETS:
        print(f"  designing {target:4s} ({label}) ...", end=" ", flush=True)
        d = _run_target(coil, evalv, target, lmax, order, workdir)
        if "error" in d:
            raise RuntimeError(f"{target}: {d['error']}")
        h = d["harmonics"]
        mc = h.get("max_contaminant") or {}
        row = {"target": target, "l": l, "label": label,
               "dominant": h["dominant"]["name"],
               "purity": h["purity"],
               "residual_fraction": h["residual_fraction"],
               "homogeneity_rms": d["rms"],
               "peak_J": d["peak_J"],
               "contaminant": mc.get("name"),
               "contaminant_fraction": mc.get("fraction", 0.0),
               "spectrum": h["spectrum"][:8]}
        rows.append(row)
        print(f"dominant {row['dominant']:5s} purity {row['purity']:.5f} "
              f"contam {row['contaminant']} {row['contaminant_fraction']:.2e}")
    return rows


def _print_table(rows, lmax):
    print(f"\n=== Shim-coil harmonic purity "
          f"(achieved-Bz decomposition over the DSV, l<={lmax}) ===")
    print(f"{'target':9s}{'dominant':10s}{'purity':>9s}"
          f"{'impurity':>11s}{'contaminant':>16s}{'peak_J':>11s}")
    for r in rows:
        contam = f"{r['contaminant'] or '-'} {r['contaminant_fraction']:.1e}"
        print(f"{r['target']:9s}{r['dominant']:10s}{r['purity']:>9.5f}"
              f"{r['residual_fraction']:>11.2e}{contam:>16s}"
              f"{r['peak_J']:>11.2e}")


def _plot(rows, png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
    except ImportError:
        print("(matplotlib not installed; skipped plot)")
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (a) impurity (1 - purity, via the LSQ residual fraction) per target, log
    #     y, coloured by harmonic order l; annotate the named contaminant.
    names = [r["target"] for r in rows]
    imp = [max(r["residual_fraction"], 1e-6) for r in rows]
    cols = [_L_COLOR.get(r["l"], "#777") for r in rows]
    xpos = np.arange(len(names))
    ax1.bar(xpos, imp, color=cols, width=0.62)
    ax1.set_yscale("log")
    ax1.set_xticks(xpos)
    ax1.set_xticklabels(names)
    ax1.set_xlabel("target spherical harmonic")
    ax1.set_ylabel("field impurity  (LSQ residual fraction)")
    ax1.grid(axis="y", alpha=0.3, which="both")
    for xi, r in zip(xpos, rows):
        ax1.annotate(r["contaminant"] or "",
                     (xi, max(r["residual_fraction"], 1e-6)),
                     textcoords="offset points", xytext=(0, 4),
                     ha="center", fontsize=8, color="#333")
    # legend = harmonic order colour key
    from matplotlib.patches import Patch
    seen = sorted({r["l"] for r in rows})
    ax1.legend(handles=[Patch(facecolor=_L_COLOR[l], label=f"l = {l}")
                        for l in seen], loc="upper center", fontsize=9)
    ax1.text(0.02, 0.96, "(a)", transform=ax1.transAxes, fontsize=11,
             va="top")

    # (b) harmonic spectrum (field fraction per named harmonic) of the HARDEST
    #     target (largest residual) -- the dominant target + small contaminants.
    hard = max(rows, key=lambda r: r["residual_fraction"])
    spec = [s for s in hard["spectrum"] if s["fraction"] > 1e-6][:8]
    labels = [s["name"] for s in spec][::-1]
    fracs = [max(s["fraction"], 1e-6) for s in spec][::-1]
    bcol = ["#993404" if s["name"] == hard["dominant"] else "#bbbbbb"
            for s in spec][::-1]
    ypos = np.arange(len(labels))
    ax2.barh(ypos, fracs, color=bcol)
    ax2.set_xscale("log")
    ax2.set_yticks(ypos)
    ax2.set_yticklabels(labels)
    ax2.set_xlabel("achieved-field fraction")
    ax2.set_ylabel(f"harmonic (target = {hard['target']})")
    ax2.grid(axis="x", alpha=0.3, which="both")
    ax2.text(0.97, 0.06, "(b)", transform=ax2.transAxes, fontsize=11,
             ha="right", va="bottom")

    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(f"Saved plot: {os.path.basename(png)}")
    return png


def main():
    ap = argparse.ArgumentParser(
        description="Shim-coil harmonic purity demo (SF spherical-harmonic "
                    "target + decomposition).")
    ap.add_argument("--order", type=int, default=2, help="psi FE order")
    ap.add_argument("--lmax", type=int, default=4,
                    help="max degree l for the achieved-field decomposition")
    ap.add_argument("--out-dir", default=str(HERE),
                    help="where the .json + .png land (default: this dir)")
    ap.add_argument("--work-dir", default=None,
                    help="mesh/JSON scratch (default: a fresh system temp dir)")
    ap.add_argument("--make-mesh", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.make_mesh:                      # re-entrant mesh-gen subprocess
        _make_meshes(args.make_mesh)
        return

    work = args.work_dir or tempfile.mkdtemp(prefix="sf_shim_")
    os.makedirs(work, exist_ok=True)
    coil = os.path.join(work, "former.vol")
    evalv = os.path.join(work, "dsv.vol")
    if not (os.path.exists(coil) and os.path.exists(evalv)):
        print("generating former + DSV meshes (subprocess) ...")
        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"
        r = subprocess.run([sys.executable, __file__, "--make-mesh", work],
                           capture_output=True, text=True, env=env,
                           timeout=900)
        if r.returncode != 0:
            sys.stderr.write(r.stderr[-1500:] + "\n")
            sys.exit("mesh generation failed (NGSolve/OCC unavailable?)")

    print(f"former a=150 mm L=500 mm | DSV r=50 mm | order {args.order} | "
          f"l<={args.lmax} | confine abe")
    rows = _design_all(coil, evalv, args.lmax, args.order, work)
    _print_table(rows, args.lmax)

    jpath = os.path.join(args.out_dir, "demo_shim_coil_purity.json")
    with open(jpath, "w") as f:
        json.dump({"order": args.order, "lmax": args.lmax,
                   "former": {"radius_m": 0.15, "length_m": 0.50},
                   "dsv_radius_m": 0.05, "rows": rows}, f, indent=2)
    print(f"\nSaved data: {os.path.basename(jpath)}")

    _plot(rows, os.path.join(args.out_dir, "demo_shim_coil_purity.png"))


if __name__ == "__main__":
    main()
