"""ESIM Karl benchmark — 3 production paths × frequency sweep.

Runs all three radia v4.46+ ESIM-enabled calc scripts on the same .vol
+ coil + BH curve at a sweep of frequencies, captures key outputs +
timing into JSON, and prints a side-by-side summary table.

The three paths:

1. calc_inductance.py --coil-solver peec  (PEEC + scalar BEM-SIBC)
2. calc_fem_kelvin.py  --impedance esim   (PEEC + FEM-HCurl Robin + Kelvin)
3. calc_fem_coilmesh.py --impedance-model esim (FEM A-V + Robin + Kelvin)

All three use the ih_fem_kelvin_demo bundled sample (coil + steel
workpiece in air + Kelvin sphere).  calc_inductance extracts the
workpiece BND from the same .vol via --wp-label sibc.

Output: results.json (sweep × path matrix of metrics).

Usage:
    python benchmark.py [--frequencies "1e4,5e4,1e5,5e5"]
                        [--output results.json]
                        [--paths inductance,fem_kelvin,fem_coilmesh]
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SAMPLES = REPO / "src" / "radia" / "panels" / "samples"
PANELS = REPO / "src" / "radia" / "panels"

VOL = SAMPLES / "ih_fem_kelvin_demo.vol"
COIL_STEP = SAMPLES / "ih_fem_kelvin_demo_coil.step"
BH = SAMPLES / "em_sample_bh.txt"

SIGMA = 2.0e6
MU_R = 100.0  # initial linear-mu guess (BH governs ESIM)
HALF_THICKNESS = 5.0e-3
CURRENT = 1.0
COIL_SIGMA = 5.8e7

DEFAULT_FREQS = [1.0e4, 5.0e4, 1.0e5, 5.0e5]  # Hz
DEFAULT_PATHS = ["inductance", "fem_kelvin", "fem_coilmesh"]


def _normalize(r, path):
    """Map per-script field names onto a canonical schema.

    `calc_fem_kelvin` predates the unified Karl JSON schema and emits
    different field names than `calc_inductance` / `calc_fem_coilmesh`.
    Translate so the benchmark, plot, and summary table see one shape.
    """
    if "error" in r:
        return r
    if path == "fem_kelvin":
        # Z_s is serialised as a Python complex-repr string, e.g.
        # "(0.01734+0.03294j)" — parse it.
        z_str = r.get("Z_s")
        if isinstance(z_str, str):
            try:
                z = complex(z_str.strip("()"))
                r["Z_s_wp_real"] = float(z.real)
                r["Z_s_wp_imag"] = float(z.imag)
            except ValueError:
                pass
        if "L" in r and "L_total_nH" not in r:
            r["L_total_nH"] = float(r["L"]) * 1e9
        if "L_peec" in r and "L_peec_nH" not in r:
            r["L_peec_nH"] = float(r["L_peec"]) * 1e9
        if "P_total" in r and "P_wp_W" not in r:
            r["P_wp_W"] = float(r["P_total"])
        if "iterations" in r and "esim_iterations" not in r:
            r["esim_iterations"] = int(r["iterations"])
    return r


def _run(cmd, label):
    print(f"  [{label}] starting", flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    wall_s = time.perf_counter() - t0
    if proc.returncode != 0:
        print(f"  [{label}] FAILED (rc={proc.returncode}): "
              f"{proc.stderr[-500:]}", flush=True)
        return {"error": proc.stderr[-500:], "wall_s": wall_s,
                "returncode": proc.returncode}
    # The calc scripts emit a JSON line at the end of stdout.
    json_lines = [ln for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("{")]
    if not json_lines:
        return {"error": "no JSON in stdout", "wall_s": wall_s,
                "stdout_tail": proc.stdout[-500:]}
    result = json.loads(json_lines[-1])
    # path label is "<path> @ <freq>"; extract path word.
    path_token = label.split(" ", 1)[0]
    result = _normalize(result, path_token)
    result["_wall_s"] = wall_s
    zr = result.get("Z_s_wp_real", 0.0)
    zi = result.get("Z_s_wp_imag", 0.0)
    Zmag = abs(complex(zr, zi))
    iters = result.get("esim_iterations")
    print(f"  [{label}] done in {wall_s:.1f}s "
          f"(esim_iter={iters}, |Z_s|={Zmag:.3e})", flush=True)
    return result


def run_inductance(freq, output_json, per_panel=False):
    cmd = [sys.executable, str(PANELS / "calc_inductance.py"),
           "--coil-solver", "peec",
           "--coil-step", str(COIL_STEP),
           "--peec-n-peri", "16",
           "--vol", str(VOL),
           "--wp-label", "sibc",
           "--sigma", str(SIGMA),
           "--mu-r", str(MU_R),
           "--impedance-model", "esim",
           "--bh-file", str(BH),
           "--half-thickness", str(HALF_THICKNESS),
           "--esim-max-iter", "15",
           "--esim-tol", "1e-3",
           "--frequency", f"{freq:.6e}",
           "--current", str(CURRENT),
           "--output", str(output_json)]
    if per_panel:
        cmd += ["--esim-per-panel", "--wp-bem-backend", "intree-dense"]
        label = f"inductance_per_panel @ {freq:.2e} Hz"
    else:
        label = f"inductance @ {freq:.2e} Hz"
    return _run(cmd, label)


def run_inductance_per_panel(freq, output_json):
    return run_inductance(freq, output_json, per_panel=True)


def run_fem_kelvin(freq, output_json):
    cmd = [sys.executable, str(PANELS / "calc_fem_kelvin.py"),
           "--vol", str(VOL),
           "--peec-step", str(COIL_STEP),
           "--frequency", f"{freq:.6e}",
           "--sigma", str(SIGMA),
           "--mu-r", str(MU_R),
           "--impedance", "esim",
           "--bh-file", str(BH),
           "--half-thickness", str(HALF_THICKNESS),
           "--solver", "pardiso",
           "--output", str(output_json)]
    return _run(cmd, f"fem_kelvin @ {freq:.2e} Hz")


def run_fem_coilmesh(freq, output_json):
    cmd = [sys.executable, str(PANELS / "calc_fem_coilmesh.py"),
           "--vol", str(VOL),
           "--frequency", f"{freq:.6e}",
           "--current", str(CURRENT),
           "--coil-sigma", str(COIL_SIGMA),
           "--sigma", str(SIGMA),
           "--mu-r", str(MU_R),
           "--half-thickness", str(HALF_THICKNESS),
           "--solver", "pardiso",
           "--impedance-model", "esim",
           "--bh-file", str(BH),
           "--esim-max-iter", "15",
           "--esim-tol", "1e-3",
           "--output", str(output_json)]
    return _run(cmd, f"fem_coilmesh @ {freq:.2e} Hz")


PATH_RUNNERS = {
    "inductance": run_inductance,
    "inductance_per_panel": run_inductance_per_panel,
    "fem_kelvin": run_fem_kelvin,
    "fem_coilmesh": run_fem_coilmesh,
}


def _post_normalize_results(results):
    """Apply _normalize() to every (path, freq) entry in a results dict.

    Safe to call on an already-loaded results.json from a previous run
    that pre-dates the schema unification.
    """
    for row in results.get("sweep", []):
        for path, r in row.get("results", {}).items():
            _normalize(r, path)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequencies", default=None,
                        help="Comma-separated Hz list, e.g. 1e4,1e5,1e6")
    parser.add_argument("--paths", default=",".join(DEFAULT_PATHS),
                        help="Comma-separated subset of "
                             "{inductance,fem_kelvin,fem_coilmesh}")
    parser.add_argument("--output", default=str(Path(__file__).parent / "results.json"))
    parser.add_argument("--tmpdir", default=None,
                        help="Directory for per-run output JSONs (kept for "
                             "inspection); default: temp/")
    args = parser.parse_args()

    if args.frequencies:
        freqs = [float(s) for s in args.frequencies.split(",")]
    else:
        freqs = DEFAULT_FREQS
    paths = [p.strip() for p in args.paths.split(",")]
    for p in paths:
        if p not in PATH_RUNNERS:
            parser.error(f"unknown path {p!r}; supported: "
                         f"{list(PATH_RUNNERS)}")

    for f in (VOL, COIL_STEP, BH):
        if not f.is_file():
            parser.error(f"missing sample: {f}")

    tmpdir = Path(args.tmpdir) if args.tmpdir else Path(__file__).parent / "tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)

    results = {
        "schema": "radia.validation.ih-esim-cross-path-benchmark.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "radia_version": None,
        "problem": {
            "vol": str(VOL.relative_to(REPO)).replace("\\", "/"),
            "coil_step": str(COIL_STEP.relative_to(REPO)).replace("\\", "/"),
            "bh_file": str(BH.relative_to(REPO)).replace("\\", "/"),
            "sigma_wp_Spm": SIGMA,
            "mu_r_wp_initial": MU_R,
            "half_thickness_m": HALF_THICKNESS,
            "current_A": CURRENT,
            "coil_sigma_Spm": COIL_SIGMA,
        },
        "frequencies_Hz": freqs,
        "paths": paths,
        "sweep": [],
    }
    try:
        import radia
        results["radia_version"] = radia.__version__
    except Exception:
        pass

    for f in freqs:
        row = {"frequency_Hz": f, "results": {}}
        for p in paths:
            out_json = tmpdir / f"{p}_f{f:.2e}.json"
            row["results"][p] = PATH_RUNNERS[p](f, out_json)
        results["sweep"].append(row)
        # incremental flush so partial sweeps survive interrupts
        with open(args.output, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)

    print()
    print(f"Wrote {args.output} ({len(freqs)} freqs × {len(paths)} paths)")
    _print_summary_table(results)


def _print_summary_table(results):
    print()
    print("=" * 90)
    print(f"Summary (radia {results['radia_version']}, {results['hostname']})")
    print("=" * 90)
    paths = results["paths"]
    header = f"{'freq[Hz]':>10s}"
    for p in paths:
        header += f"  {p}: t[s]  P_wp[mW]  L[nH]  iter"
    print(header)
    print("-" * len(header))
    for row in results["sweep"]:
        f = row["frequency_Hz"]
        line = f"{f:>10.2e}"
        for p in paths:
            r = row["results"][p]
            if "error" in r:
                line += f"  {p}: FAILED"
                continue
            t = r.get("_wall_s", 0.0)
            P = r.get("P_wp_W") or r.get("P_total_W") or r.get("P_total")
            if P is None:
                P = 0.0
            L = r.get("L_total_nH") or (r.get("L_total", 0) * 1e9)
            ni = r.get("esim_iterations", "?")
            line += f"  t={t:5.1f}  P={P*1e3:7.3e}  L={L:6.2f}  i={ni}"
        print(line)
    print("=" * 90)


if __name__ == "__main__":
    main()
