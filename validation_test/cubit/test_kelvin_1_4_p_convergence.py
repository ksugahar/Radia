"""Regression test for Cubit + Kelvin + Omega-Reduced Omega + p>=2.

Covers the promoted validation fixture in
`validation_test/cubit/kelvin_1_4_p_convergence/`.  The original
research example copy has been retired; this fixture is the maintained
examples -> validation_test promotion target.  Runs:

  1. `mesh_and_export.py` -- builds a 1/4 sector magnetic-sphere mesh
     (mu_r = 100, R = 0.20 m, Kelvin offset 0.60 m) and exports
     `.vol` at orders 1, 2, 3 into a temporary directory, with the
     radia-mcp Cubit plugin (RADIA_CUBIT_PLUGIN_DIR; skipped visibly
     when unset).

  2. `solve_p_convergence.py` -- loads each .vol, runs the
     Omega-Reduced Omega + Kelvin solver, probes Hz at the origin.
     Asserts:
       - each exported .vol carries the Kelvin pairs, equal to an independent
         post-hoc pairing of its bare copy (unmatched = 0, max_dist < 1e-10 m)
       - Hz at p=1 within +/-10% of analytical
       - Hz at p=2 within +/-1.5% of analytical 3 / (mu_r + 2) * H0
       - Hz at p=3 within +/-1.0% of analytical

This test is the regression guard for the ANCHOR-DETERMINISTIC fix in
`_add_kelvin_cubit_reduction` (`radia_mcp.cubit.kelvin`,
2026-04-25) and the `subtract WITHOUT keep + recreate` workaround for
Cubit 2025.3's `subtract ... keep` bug.
"""

import json
import os
import subprocess
import sys
import tempfile
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


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kelvin_pair_oracle  # noqa: E402


class TestCubitKelvin14PConvergence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        plugin_dir = os.environ.get("RADIA_CUBIT_PLUGIN_DIR")
        if not plugin_dir:
            raise unittest.SkipTest(
                "RADIA_CUBIT_PLUGIN_DIR (directory holding radia_cubit_export.dll) "
                "is not set; this validation needs licensed Cubit and the plugin")
        cls.tmp = tempfile.TemporaryDirectory()
        out = cls.tmp.name
        cls.mesh_log = _run("mesh_and_export.py",
                            ["--orders", "1,2,3", "--mesh-size", "0.025",
                             "--plugin-dir", plugin_dir, "--out-dir", out])
        cls.json_path = os.path.join(out, "p_convergence.json")
        cls.solve_log = _run("solve_p_convergence.py",
                             ["--orders", "1,2,3", "--in-dir", out,
                              "--out-json", cls.json_path])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_mesh_machine_precision_pairing(self):
        """Every order's .vol carries the Kelvin pairs the post-hoc oracle finds."""
        for order in (1, 2, 3):
            vol = os.path.join(self.tmp.name, f"sphere_1_4_p{order}.vol")
            info = kelvin_pair_oracle.compare(vol, self.tmp.name)
            self.assertEqual(info["groups"], 1, info)
            self.assertTrue(info["covers_kelvin"], info)
            self.assertEqual(info["bare_groups_before_pairing"], 0, info)
            self.assertEqual(info["posthoc_unmatched"], 0, info)
            self.assertLess(info["posthoc_max_dist"], 1e-10, info)
            self.assertGreater(info["exported_pairs"], 0, info)
            self.assertEqual(info["exported_pairs"], info["posthoc_pairs"], info)
            self.assertEqual(info["differing_pairs"], 0, info)
            self.assertLess(max(abs(a - b) for a, b in zip(info["posthoc_offset"], (0.0, 0.0, 0.6))),
                            1e-9, info)

    def test_solver_p_convergence_within_tolerance(self):
        """Hz at origin matches analytical within tolerance at p=1, 2 and 3."""
        with open(self.json_path) as f:
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
