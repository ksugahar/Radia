"""IH per-element ESIM (f, I_port) sweep — 4 freqs x 4 currents = 16 cases.

For each combination of frequency and port current, run BOTH the
scalar-Karl and per-element-Karl variants and collect P_wp, |H_t|
range, |Z_s| contrast, Karl convergence info into a single results
JSON.

The output is a 2D heatmap of (per-element-vs-scalar P_wp gap) vs
(f, I) showing the operating regime where per-element ESIM matters
most for IH design.  Figure target: IGTE journal paper.

Run time: ~5 min per ESIM case (cell-solver loop dominates),
~5 min per scalar case -> 16 x 2 x 5 min = ~2.5 hours total.

Usage:
    python examples/ih_esim_benchmark/sweep_f_I.py [OUT_DIR]

Output:
    <OUT_DIR>/sweep_results.json    -- aggregated results
    <OUT_DIR>/I*_f*_{scalar,per_panel}.json -- per-case JSONs
"""
import sys
import json
import os
import time
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CALC = REPO / "src" / "radia" / "panels" / "calc_inductance.py"
SAMPLES = REPO / "src" / "radia" / "panels" / "samples"

FREQS_HZ = [10000, 50000, 100000, 500000]
CURRENTS_A = [1.0, 10.0, 100.0, 300.0]

CMD_BASE = [
    sys.executable, str(CALC),
    "--coil-step", str(SAMPLES / "ih_fem_kelvin_demo_coil.step"),
    "--coil-solver", "peec",
    "--vol", str(SAMPLES / "ih_bem_sample_p1.vol"),
    "--wp-label", "sibc",
    "--sigma", "2e6", "--mu-r", "100", "--half-thickness", "0.005",
    "--coil-sigma", "5.8e7",
    "--impedance-model", "esim",
    "--bh-file", str(SAMPLES / "em_sample_bh.txt"),
    "--esim-max-iter", "30", "--esim-tol", "1e-3",
    "--esim-relax", "0.5",
    "--esim-anderson-m", "5",     # safeguarded Anderson
    "--h1-order", "1",
    "--wp-bem-backend", "intree-dense",
]


def run_one(freq_hz: float, current_A: float, per_panel: bool,
            out_path: Path) -> dict:
    cmd = CMD_BASE + [
        "--frequency", str(freq_hz),
        "--current", str(current_A),
        "--output", str(out_path),
    ]
    if per_panel:
        cmd.append("--esim-per-panel")
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        print(f"  ERROR ({dt:.0f}s): {r.stderr[-500:]}")
        return {"error": r.stderr[-500:], "elapsed_s": dt}
    with open(out_path) as f:
        d = json.load(f)
    return {
        "P_wp_W": d.get("P_wp_W"),
        "L_total_uH": d.get("L_total_uH"),
        "H_t_rms": d.get("H_t_rms_A_per_m"),
        "esim_iterations": d.get("esim_iterations"),
        "esim_converged": d.get("esim_converged"),
        "esim_anderson_restarts": d.get("esim_anderson_restarts"),
        "esim_anderson_clips": d.get("esim_anderson_clips"),
        "elapsed_s": dt,
    }


def main():
    # Data Persistence Policy: canonical run writes JSON into the repo
    # (committed alongside the figure), NOT into transient C:/temp.
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "sweep_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Sweep output dir: {out_dir}")
    print(f"Cases: {len(FREQS_HZ) * len(CURRENTS_A) * 2}")
    print()

    results = []
    t_total = time.time()
    for I_A in CURRENTS_A:
        for f_Hz in FREQS_HZ:
            for per_panel in [False, True]:
                tag = f"I{int(I_A)}_f{int(f_Hz/1000)}k_{'per_panel' if per_panel else 'scalar'}"
                out_path = out_dir / f"{tag}.json"
                t0 = time.time()
                print(f"[{len(results)+1:2d}/32] {tag} ... ", end="", flush=True)
                if out_path.exists():
                    with open(out_path) as f:
                        d = json.load(f)
                    summary = {
                        "P_wp_W": d.get("P_wp_W"),
                        "L_total_uH": d.get("L_total_uH"),
                        "H_t_rms": d.get("H_t_rms_A_per_m"),
                        "esim_iterations": d.get("esim_iterations"),
                        "esim_converged": d.get("esim_converged"),
                        "esim_anderson_restarts": d.get("esim_anderson_restarts"),
                        "esim_anderson_clips": d.get("esim_anderson_clips"),
                        "elapsed_s": -1,
                        "from_cache": True,
                    }
                    print(f"cached, P_wp={summary['P_wp_W']:.3f} W")
                else:
                    summary = run_one(f_Hz, I_A, per_panel, out_path)
                    if "error" not in summary:
                        print(f"P_wp={summary['P_wp_W']:.3f} W "
                              f"({summary['elapsed_s']:.0f}s)")
                results.append({
                    "frequency_Hz": f_Hz,
                    "current_A": I_A,
                    "per_panel": per_panel,
                    **summary,
                })
                with open(out_dir / "sweep_results.json", "w") as f:
                    json.dump({"runs": results}, f, indent=2)

    elapsed = time.time() - t_total
    print()
    print(f"Total time: {elapsed/60:.1f} min ({elapsed:.0f}s)")
    print(f"Aggregated results: {out_dir / 'sweep_results.json'}")


if __name__ == "__main__":
    main()
