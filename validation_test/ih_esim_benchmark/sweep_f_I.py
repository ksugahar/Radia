"""IH per-element ESIM (f, I_port) dense sweep.

For each combination of frequency and port current, run BOTH the
scalar-Karl and per-element-Karl variants and collect P_wp, |H_t|
range, |Z_s| contrast, Karl convergence info into a single results
JSON.

The output is a 2D heatmap of (per-element-vs-scalar P_wp gap) vs
(f, I) showing the operating regime where per-element ESIM matters
most for IH design.  Figure target: IGTE 2026 digest.

Grid: 9 currents x 6 frequencies x 2 modes = 108 cases.
Run time on LAB is typically 3--4 hours, depending on high-current
stall behaviour.

Usage:
    python validation_test/ih_esim_benchmark/sweep_f_I.py [OUT_DIR]

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

FREQS_HZ = [10000, 20000, 50000, 100000, 200000, 500000]
CURRENTS_A = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]

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
    # calc_inductance.py may exit 0 but write {"error": "..."} into the
    # per-case JSON when an internal exception is caught (e.g. a transient
    # NAS-share import flicker on S:).  Treat that as an error so the
    # sweep records the case and moves on rather than crashing on
    # P_wp_W = None.  Also delete the placeholder JSON so a future restart
    # re-runs the case rather than re-caching the error.
    with open(out_path) as f:
        d = json.load(f)
    if "error" in d and d.get("P_wp_W") is None:
        try:
            out_path.unlink()
        except OSError:
            pass
        return {"error": str(d["error"])[-500:], "elapsed_s": dt}
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
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "sweep_data_dense"
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
                from_cache = False
                if out_path.exists():
                    with open(out_path) as f:
                        d = json.load(f)
                    if "error" in d and d.get("P_wp_W") is None:
                        # Stale error-placeholder from a previous
                        # transient failure -- delete so this case
                        # re-runs instead of caching the error.
                        try:
                            out_path.unlink()
                        except OSError:
                            pass
                        summary = run_one(f_Hz, I_A, per_panel, out_path)
                    else:
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
                        from_cache = True
                else:
                    summary = run_one(f_Hz, I_A, per_panel, out_path)
                # Single print: cached vs fresh, with P_wp_W=None guard.
                if "error" not in summary:
                    pw = summary.get("P_wp_W")
                    pw_s = f"{pw:.3f} W" if pw is not None else "None"
                    if from_cache:
                        print(f"cached, P_wp={pw_s}")
                    else:
                        print(f"P_wp={pw_s} "
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
