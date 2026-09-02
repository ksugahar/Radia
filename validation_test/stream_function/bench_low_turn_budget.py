"""bench_low_turn_budget.py -- does the stream-function method succeed at LOW
turn count?  Quantify "fewest turns to meet a field spec" per strategy.

Sweeps the turn count across the low-turn manufacture strategies on a
representative SF coil design (the golden cylinder former, MRI-gradient-scale
DSV; Gx gradient + Z2 shim targets, --confine abe so the contours close) and
reports, for a set of relative-field-rms specs, the FEWEST turns each strategy
needs.  Answers the question with hard numbers.

METRIC (apples-to-apples): the EQUAL-CURRENT single-series field residual of the
N-turn loop set (||I*sum_k s_k f_k - B|| / ||B||, one common current I) --
exactly what a single wound wire delivers per turn count:
  * equal_dI_contour   N equal-deltaI psi iso-contours (the classic SF discretisation)
  * optimize_levels    the N contour LEVELS optimised (--optimize-levels)
  * greedy_contour/pin/bubble   greedy constructive (--greedy-turns, monotone trace)
Context curves (JSON + reported, not the headline): loops_homogeneity_rms (N
INDEPENDENT currents -- the driven-array floor) and wire_homogeneity_rms (the
delivered single-stroke wire AFTER chaining, optimised levels).

HONEST scope: ONE representative geometry (the golden cylinder fixture), TWO
targets.  Few turns = coarse quantisation of a continuous current density, so
the achievable field accuracy is bounded; this bench measures that bound.

POST-FIX NOTE (grad-psi orientation, commit 630ece06): the production single-wire
orientation now uses K = n_hat x grad_s(psi), not sign(f.B), so Z2 (l=2) single-
wire NOW works -- optimize_levels reaches a few-% equal-current residual instead
of ~70%.  The earlier "Z2 single-wire fails / needs a driven array" reading was an
ORIENTATION BUG, not a physics limit (see verify_gradpsi_orientation.py).  Gx
(l=1) is unchanged.  These numbers are POST-fix.

Outputs (committed, per the Data Persistence Policy): bench_low_turn_budget.json
+ bench_low_turn_budget.png next to this script.  The .vol meshes are gitignored
and regenerated into a temp dir by validation_test/panels/fixtures/make_streamfunction_vol.py.
Run:  python bench_low_turn_budget.py
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
FIXTURE = os.path.join(REPO, "tests", "panels", "fixtures",
                       "make_streamfunction_vol.py")

TARGETS = [("Gx", "x"), ("Z2", "z*z-(x*x+y*y)/2")]
N_SET = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16]
GREEDY_MAX = 16
SPECS = [0.30, 0.20, 0.10, 0.05, 0.02]
CONFINE = "abe"          # the recommended BC: closes the contours on a finite former
EVAL_MAX = 40            # >=40 so Z2 (l=2) is not under-sampled (24 reads false)


def _env():
    e = os.environ.copy()
    e["PYTHONIOENCODING"] = "utf-8"
    return e


def run_calc(coil, evalv, target_cf, extra):
    cmd = [sys.executable, CALC, "--coil-vol", coil, "--eval-vol", evalv,
           "--order", "1", "--target-cf", target_cf, "--method", "manufacture",
           "--eval-max", str(EVAL_MAX), "--confine", CONFINE] + extra
    r = subprocess.run(cmd, capture_output=True, text=True, env=_env(),
                       timeout=1200)
    if r.returncode != 0:
        raise RuntimeError(f"calc failed (rc={r.returncode}) for {extra}:\n"
                           f"{r.stderr[-1500:]}")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return json.loads(lines[-1])


def first_turns_meeting(curve, spec):
    """curve = list of (n_turns, rms); return the FEWEST n with rms <= spec
    (None if never).  For a non-monotone curve this is the first crossing."""
    for n, rms in curve:
        if rms is not None and rms <= spec:
            return n
    return None


def gen_fixture(outdir):
    r = subprocess.run([sys.executable, FIXTURE, outdir],
                       capture_output=True, text=True, env=_env(), timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"fixture gen failed:\n{r.stderr[-800:]}")
    coil = os.path.join(outdir, "coil_cyl_surf.vol")
    evalv = os.path.join(outdir, "eval_dsv.vol")
    if not (os.path.exists(coil) and os.path.exists(evalv)):
        raise RuntimeError("fixture .vol not produced")
    return coil, evalv


def sweep_one_target(coil, evalv, name, cf):
    print(f"\n=== target {name} ({cf}) ===")
    # contour + optimize-levels: one --optimize-levels call per N gives BOTH the
    # uniform (equal-deltaI) and the optimised equal-current residual, plus the
    # independent-current floor and the delivered single-stroke wire.
    equal_dI, opt_lv, loops_floor, wire_deliv = [], [], [], []
    for N in N_SET:
        o = run_calc(coil, evalv, cf, ["--nlevels", str(N), "--optimize-levels"])
        equal_dI.append((N, o.get("equal_current_rms_uniform")))
        opt_lv.append((N, o.get("equal_current_rms_optimized")))
        loops_floor.append((N, o.get("loops_homogeneity_rms")))
        wire_deliv.append((N, o.get("wire_homogeneity_rms")))
        print(f"  N={N:2d}  equal-dI={o.get('equal_current_rms_uniform'):.4f}  "
              f"opt-levels={o.get('equal_current_rms_optimized'):.4f}  "
              f"loops(indep)={o.get('loops_homogeneity_rms'):.4f}  "
              f"wire(deliv)={o.get('wire_homogeneity_rms'):.4f}")
    # greedy: one call per dict gives the full monotone trace n=1..GREEDY_MAX
    greedy = {}
    for dct in ("contour", "pin", "bubble"):
        extra = ["--greedy-turns", str(GREEDY_MAX), "--greedy-dict", dct]
        if dct == "bubble":
            extra += ["--pin-tiling-pins", "60"]
        o = run_calc(coil, evalv, cf, extra)
        tr = [(t["n_turns"], t["rms"]) for t in o.get("greedy_trace", [])]
        greedy[dct] = tr
        floor = tr[-1][1] if tr else None
        print(f"  greedy-{dct:7s} {len(tr)} turns, floor rms="
              f"{floor:.4f}" if floor is not None else f"  greedy-{dct}: empty")
    return {
        "target": name, "target_cf": cf,
        "equal_dI_contour": equal_dI,
        "optimize_levels": opt_lv,
        "loops_indep_floor": loops_floor,
        "wire_delivered": wire_deliv,
        "greedy": greedy,
    }


def spec_table(res):
    """Fewest turns to meet each spec, per strategy (the headline answer)."""
    rows = {}
    curves = {
        "equal_dI_contour": res["equal_dI_contour"],
        "optimize_levels": res["optimize_levels"],
        "greedy_contour": res["greedy"].get("contour", []),
        "greedy_pin": res["greedy"].get("pin", []),
        "greedy_bubble": res["greedy"].get("bubble", []),
    }
    for strat, curve in curves.items():
        rows[strat] = {f"{int(s*100)}%": first_turns_meeting(curve, s)
                       for s in SPECS}
    return rows


def plot(results, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(results), figsize=(5.2 * len(results), 4.0),
                             squeeze=False)
    styles = [("equal_dI_contour", "equal-dI contour", "o-", "C0"),
              ("optimize_levels", "optimize-levels", "s-", "C1"),
              ("greedy_contour", "greedy contour", "^--", "C2"),
              ("greedy_pin", "greedy pin", "v--", "C3"),
              ("greedy_bubble", "greedy bubble", "d--", "C4")]
    for ax, res in zip(axes[0], results):
        curves = {
            "equal_dI_contour": res["equal_dI_contour"],
            "optimize_levels": res["optimize_levels"],
            "greedy_contour": res["greedy"].get("contour", []),
            "greedy_pin": res["greedy"].get("pin", []),
            "greedy_bubble": res["greedy"].get("bubble", []),
        }
        for key, lbl, sty, col in styles:
            c = [(n, r) for (n, r) in curves[key] if r is not None]
            if c:
                ax.plot([n for n, _ in c], [r for _, r in c], sty, color=col,
                        label=lbl, ms=4, lw=1.2)
        floor = [(n, r) for (n, r) in res["loops_indep_floor"] if r is not None]
        if floor:
            ax.plot([n for n, _ in floor], [r for _, r in floor], ":",
                    color="0.5", label="indep-current floor", lw=1.0)
        for s in SPECS:
            ax.axhline(s, color="0.85", lw=0.8, zorder=0)
        ax.set_yscale("log")
        ax.set_xlabel("turns (N)")
        ax.set_ylabel("equal-current field residual (relative rms)")
        ax.set_title(f"target {res['target']}")    # subplot label, NOT in-figure title
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    with tempfile.TemporaryDirectory() as td:
        coil, evalv = gen_fixture(td)
        results = [sweep_one_target(coil, evalv, name, cf)
                   for name, cf in TARGETS]

    out = {
        "benchmark": "low_turn_budget",
        "hostname": platform.node(),
        "problem": {
            "former": "cylinder r=0.15 m L=0.50 m (golden fixture)",
            "dsv": "sphere r=0.05 m",
            "confine": CONFINE, "order": 1, "eval_max": EVAL_MAX,
            "n_set": N_SET, "greedy_max": GREEDY_MAX, "specs": SPECS,
        },
        "metric": "equal-current single-series field residual (relative rms) of "
                  "the N-turn loop set; one common current",
        "results": results,
        "spec_table": {r["target"]: spec_table(r) for r in results},
    }
    json_path = os.path.join(HERE, "bench_low_turn_budget.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    png_path = os.path.join(HERE, "bench_low_turn_budget.png")
    try:
        plot(results, png_path)
    except Exception as e:                                  # noqa: BLE001
        print(f"plot skipped: {e}")
        png_path = None

    print("\n================ FEWEST TURNS TO MEET SPEC ================")
    for tgt, tbl in out["spec_table"].items():
        print(f"\n  target {tgt}  (relative field rms spec -> fewest turns)")
        header = "    %-18s " % "strategy" + "".join(
            "%7s" % f"{int(s*100)}%" for s in SPECS)
        print(header)
        for strat, row in tbl.items():
            cells = "".join("%7s" % ("-" if row[f"{int(s*100)}%"] is None
                                     else row[f"{int(s*100)}%"]) for s in SPECS)
            print("    %-18s %s" % (strat, cells))
    print(f"\nSaved {json_path}")
    if png_path:
        print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
