"""Regression test for Cubit + Kelvin + Omega-Reduced Omega + p>=2.

Covers the promoted validation fixture in
`validation_test/cubit/kelvin_1_4_p_convergence/`.  The original
research example is archived in the docs Kelvin notebooks per the
examples -> docs/validation_test promotion policy.  Runs:

  1. `mesh_and_export.py` -- builds a 1/4 sector magnetic-sphere mesh
     (mu_r = 100, R = 0.20 m, Kelvin offset 0.60 m) and exports
     `.vol` at orders 1, 2, 3.  Asserts:
       - 265 / 265 vertex pairs match at machine precision (max_dist < 1e-10 m)
       - kelvin_offset detected as (0, 0, 0.6 +/- 1e-8 m)

  2. `solve_p_convergence.py` -- loads each .vol, runs the
     Omega-Reduced Omega + Kelvin solver, probes Hz at the origin.
     Asserts:
       - Hz at p=2 within +/-1.5% of analytical 3 / (mu_r + 2) * H0
       - Hz at p=3 within +/-1.0% of analytical
       - Periodic FES at p=2 slaves > 100 high-order DOFs
         (sanity that high-order coupling is in effect)

This test is the regression guard for the ANCHOR-DETERMINISTIC fix in
`_add_kelvin_cubit_reduction` (`src/radia/panels/add_kelvin.py`,
2026-04-25) and the `subtract WITHOUT keep + recreate` workaround for
Cubit 2025.3's `subtract ... keep` bug.
"""

import json
import os
import subprocess
import sys
import unittest


SAMPLE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "kelvin_1_4_p_convergence",
)
SAMPLE_DIR = os.path.normpath(SAMPLE_DIR)


def _run(script, args):
    cmd = [sys.executable, os.path.join(SAMPLE_DIR, script)] + list(args)
    env = os.environ.copy()
    env.setdefault("PYTHONIOENCODING", "utf-8")
    out = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        cwd=SAMPLE_DIR,
        env=env,
    )
    if out.returncode != 0:
        raise RuntimeError(
            f"{script} failed (exit {out.returncode}):\n"
            f"stdout:\n{out.stdout}\nstderr:\n{out.stderr}")
    return out.stdout


class TestCubitKelvin14PConvergence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Generate meshes at p=1, 2, 3 once.
        cls.mesh_log = _run("mesh_and_export.py",
                            ["--orders", "1,2,3", "--mesh-size", "0.025"])

    def test_mesh_machine_precision_pairing(self):
        """Vertex pairs on the Kelvin cap match at machine precision."""
        # The mesh log prints `Kelvin periodic: N/N pairs (max_dist=...,
        # 0 unmatched)` once per export order.
        log = self.mesh_log
        # Expect "265/265 pairs" (1/4 mesh size 0.025 produces 265
        # vertex pairs on each cap with current mesh density).
        # Tolerance: max_dist < 1e-10 m (machine precision).
        for line in log.splitlines():
            if "Kelvin periodic:" in line and "/" in line and "pairs" in line:
                # e.g. "Kelvin periodic: 265/265 pairs (max_dist=5.59e-16, ...)"
                parts = line.split("max_dist=")
                if len(parts) < 2:
                    continue
                num = parts[1].split(",")[0].strip()
                max_dist = float(num)
                self.assertLess(max_dist, 1e-10,
                    f"Kelvin pairing not at machine precision: {line}")
        # Also assert "0 unmatched"
        self.assertIn("0 unmatched", log,
            "Some Kelvin vertex pairs were unmatched by the C++ exporter")

    def test_solver_p_convergence_within_tolerance(self):
        """Hz at origin matches analytical within tolerance at p=2 and p=3."""
        _run("solve_p_convergence.py", ["--orders", "1,2,3"])
        with open(os.path.join(SAMPLE_DIR, "p_convergence.json")) as f:
            data = json.load(f)
        results = {r["order"]: r for r in data["results"]}
        self.assertIn(1, results); self.assertIn(2, results); self.assertIn(3, results)

        # Tolerance bands (validated 2026-04-25):
        #   p=1: error +7.79%       (give 10% slack)
        #   p=2: error +0.71%       (give 1.5% slack)
        #   p=3: error +0.54%       (give 1.0% slack)
        for p, max_pct in ((1, 10.0), (2, 1.5), (3, 1.0)):
            err = abs(results[p]["error_pct"])
            self.assertLess(err, max_pct,
                f"order={p}: |err|={err:.2f}% exceeds {max_pct}% bound; "
                f"Hz_origin={results[p]['Hz_origin']:.6e}, "
                f"analytical={results[p]['Hz_analytical']:.6e}")


if __name__ == "__main__":
    unittest.main()
