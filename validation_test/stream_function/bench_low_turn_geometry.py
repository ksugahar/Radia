"""bench_low_turn_geometry.py -- is the Z2 (l=2 shim) SINGLE-WIRE limit fixable
by GEOMETRY?

bench_low_turn_budget.py found Z2 fails as a single series wire on the stubby
golden former (a=0.15, L=0.5, aspect L/2a=1.7) -- only a driven array works.
A Z2 (zonal/axial) shim physically needs AXIAL extent, so this sweeps the
cylinder LENGTH L (the aspect ratio) and the DSV radius and asks: does the best
achievable SINGLE-WIRE Z2 residual drop into the few-% range with a better-
matched former?  If yes, Z2 single-wire "works" and the earlier failure was
geometry; if it stays high, the single-wire l=2 limit needs more research.

Per geometry/target it reports:
  cont_homo    continuous regularised design homogeneity (the design ceiling)
  best_wire    best single-series-current residual over N (optimize-levels)
  wire_deliv   delivered single-stroke wire residual at that N
  indep_floor  N independent currents (driven-array floor) -- for contrast

--confine abe (closes contours).  Metric = relative field rms.  Targets: Z2
(the hard l=2 case) + Gx (l=1 sanity contrast).

CONCLUSION (measured): geometry does NOT change the Z2 single-wire residual --
it is flat across aspect L/2a = 1.7 -> 6.7 and for smaller DSV.  The lever was
NEVER geometry: the Z2 CONTINUOUS design is ~0.9% everywhere (order 1 and 2,
purity 0.9999) and INDEPENDENT currents on the same contour loops reach ~5e-4
(the loop basis is perfect).  The single-wire bottleneck was the equal-current
contour ORIENTATION, now FIXED in production (commit 630ece06): loops are
oriented by the consistent grad-psi winding K = n_hat x grad_s(psi), not the
sign(f.B) flip that flattened the saddle's mixed current senses.  POST-FIX this
bench's best_wire for Z2 is a few-% (was ~70%); see verify_gradpsi_orientation.py
for the controlled before/after (Z2 88% -> 5.6% loop-set).  Gx (l=1, monotone psi)
is unchanged (~5e-3) and improves with length.

NOTE on wire_deliv (the DELIVERED single-stroke column): best_wire is the
connector-FREE loop set; wire_deliv adds the field-aware single-stroke connectors.
For a multi-region (saddle) psi the connectors are the bottleneck and were
chain-resolution-limited -- the field-aware cut-opt now AUTO-scales its resolution
with the loop count (calc_streamfunction --chain-ncut/--chain-passes), which cut a
clustered Z2 delivered wire ~0.20 -> ~0.067 on a longer former.  This bench's
wire_deliv reflects that auto-resolution; it stays geometry-limited (the stubby
aspect 1.7 former holds a ~0.16 connector floor that a longer former relaxes).

Outputs (committed): bench_low_turn_geometry.json + .png next to this script.
Run:  python bench_low_turn_geometry.py
"""
from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CALC = os.path.join(REPO, "src", "radia", "panels", "calc_streamfunction.py")

# parameterised cylinder-former + DSV-sphere generator (run as a subprocess so
# the NGSolve/Netgen import stays out of this process)
_GEN = r"""
import os, sys
from netgen.occ import Cylinder, Sphere, Pnt, Z, OCCGeometry
from ngsolve import TaskManager
outdir, a, L, dsv = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
os.makedirs(outdir, exist_ok=True)
with TaskManager():
    cyl = Cylinder(Pnt(0, 0, -L/2), Z, r=a, h=L)
    lateral = max(cyl.faces, key=lambda f: f.mass)
    OCCGeometry(lateral).GenerateMesh(maxh=0.04).Save(os.path.join(outdir, "coil.vol"))
    OCCGeometry(Sphere(Pnt(0, 0, 0), dsv)).GenerateMesh(maxh=dsv*0.5).Save(
        os.path.join(outdir, "eval.vol"))
print("ok")
"""

TARGETS = {"Z2": "z*z-(x*x+y*y)/2", "Gx": "x"}
N_SET = [6, 10, 16]            # optimize-levels turn counts to scan for the best
EVAL_MAX = 40                  # >=40: fewer points under-sample Z2 (a 24-pt design
                               # reads a false ~10% Z2 ceiling; 40 pts -> ~0.9%)

# (label, a, L, dsv, aspect=L/2a)
GEOMS = [
    ("L=0.5",          0.15, 0.5, 0.05),
    ("L=1.0",          0.15, 1.0, 0.05),
    ("L=1.5",          0.15, 1.5, 0.05),
    ("L=2.0",          0.15, 2.0, 0.05),
    ("L=1.0 DSV0.035", 0.15, 1.0, 0.035),
    ("L=1.0 DSV0.025", 0.15, 1.0, 0.025),
]
# Gx is a cheap sanity contrast -- only the extremes
GX_GEOMS = {"L=0.5", "L=2.0"}


def _env():
    e = os.environ.copy()
    e["PYTHONIOENCODING"] = "utf-8"
    return e


def gen(outdir, a, L, dsv):
    r = subprocess.run([sys.executable, "-c", _GEN, outdir, str(a), str(L),
                        str(dsv)], capture_output=True, text=True, env=_env(),
                       timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"mesh gen failed (a={a},L={L},dsv={dsv}):\n"
                           f"{r.stderr[-800:]}")
    return os.path.join(outdir, "coil.vol"), os.path.join(outdir, "eval.vol")


def calc(coil, evalv, cf, extra):
    cmd = [sys.executable, CALC, "--coil-vol", coil, "--eval-vol", evalv,
           "--order", "1", "--target-cf", cf, "--eval-max", str(EVAL_MAX),
           "--confine", "abe"] + extra
    r = subprocess.run(cmd, capture_output=True, text=True, env=_env(),
                       timeout=1200)
    if r.returncode != 0:
        raise RuntimeError(f"calc failed {extra}:\n{r.stderr[-1200:]}")
    return json.loads([ln for ln in r.stdout.splitlines() if ln.strip()][-1])


def measure(coil, evalv, cf):
    cont = calc(coil, evalv, cf, ["--method", "design"])
    best = None
    for N in N_SET:
        o = calc(coil, evalv, cf,
                 ["--method", "manufacture", "--nlevels", str(N),
                  "--optimize-levels"])
        rec = (o.get("equal_current_rms_optimized"),
               o.get("wire_homogeneity_rms"),
               o.get("loops_homogeneity_rms"), N)
        if best is None or (rec[0] is not None and rec[0] < best[0]):
            best = rec
    return {
        "cont_homo": cont.get("homogeneity_rms"),
        "peak_J": cont.get("peak_J"),
        "best_wire": best[0], "wire_deliv": best[1],
        "indep_floor": best[2], "best_N": best[3],
    }


def main():
    rows = []
    with tempfile.TemporaryDirectory() as td:
        for i, (label, a, L, dsv) in enumerate(GEOMS):
            gd = os.path.join(td, f"g{i}")
            coil, evalv = gen(gd, a, L, dsv)
            aspect = L / (2 * a)
            for tname, cf in TARGETS.items():
                if tname == "Gx" and label not in GX_GEOMS:
                    continue
                m = measure(coil, evalv, cf)
                m.update({"geom": label, "target": tname, "a": a, "L": L,
                          "dsv": dsv, "aspect": aspect})
                rows.append(m)
                print(f"  {label:16s} {tname}: cont={m['cont_homo']:.4f} "
                      f"best_wire={m['best_wire']:.4f} (N={m['best_N']}) "
                      f"deliv={m['wire_deliv']:.4f} indep={m['indep_floor']:.4f}")

    out = {
        "benchmark": "low_turn_geometry",
        "hostname": platform.node(),
        "problem": {"a": 0.15, "L_sweep": [0.5, 1.0, 1.5, 2.0],
                    "dsv_sweep": [0.05, 0.035, 0.025], "confine": "abe",
                    "order": 1, "eval_max": EVAL_MAX, "n_set": N_SET},
        "metric": "relative field rms; best_wire = best single-series-current "
                  "(optimize-levels) over N; indep_floor = N independent currents",
        "rows": rows,
    }
    json_path = os.path.join(HERE, "bench_low_turn_geometry.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    # plot: Z2 single-wire vs aspect ratio (L-sweep, DSV=0.05) + Gx contrast
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        Lsweep = [r for r in rows if r["dsv"] == 0.05]
        z2 = sorted([r for r in Lsweep if r["target"] == "Z2"],
                    key=lambda r: r["aspect"])
        fig, ax = plt.subplots(figsize=(5.4, 4.0))
        ax.plot([r["aspect"] for r in z2], [r["cont_homo"] for r in z2],
                "o-", color="C0", label="Z2 continuous design")
        ax.plot([r["aspect"] for r in z2], [r["best_wire"] for r in z2],
                "s-", color="C1", label="Z2 best single wire (opt-levels)")
        ax.plot([r["aspect"] for r in z2], [r["wire_deliv"] for r in z2],
                "^-", color="C3", label="Z2 delivered single-stroke wire")
        ax.plot([r["aspect"] for r in z2], [r["indep_floor"] for r in z2],
                ":", color="0.5", label="Z2 indep-current floor (array)")
        gx = sorted([r for r in Lsweep if r["target"] == "Gx"],
                    key=lambda r: r["aspect"])
        if gx:
            ax.plot([r["aspect"] for r in gx], [r["best_wire"] for r in gx],
                    "d--", color="C2", label="Gx best single wire (contrast)")
        ax.axhline(0.05, color="0.85", lw=0.8)
        ax.set_yscale("log")
        ax.set_xlabel("former aspect ratio  L / 2a")
        ax.set_ylabel("relative field residual (rms)")
        ax.set_title("does a longer former fix Z2 single-wire?")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7)
        fig.tight_layout()
        png = os.path.join(HERE, "bench_low_turn_geometry.png")
        fig.savefig(png, dpi=150)
        plt.close(fig)
        print(f"Saved {png}")
    except Exception as e:                                  # noqa: BLE001
        print(f"plot skipped: {e}")

    print("\n========== Z2 single-wire vs geometry ==========")
    print("  %-16s %8s %10s %10s %10s" % ("geom", "aspect", "cont", "best_wire",
                                          "indep"))
    for r in rows:
        if r["target"] == "Z2":
            print("  %-16s %8.2f %10.4f %10.4f %10.4f" %
                  (r["geom"], r["aspect"], r["cont_homo"], r["best_wire"],
                   r["indep_floor"]))
    print(f"\nSaved {json_path}")


if __name__ == "__main__":
    main()
