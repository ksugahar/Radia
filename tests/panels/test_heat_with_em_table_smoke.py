"""Smoke tests for calc_heat_with_em_table.py.

These cover the CLI surface, helper functions, and error paths that
do NOT require a real workpiece .vol or a calc_fem_kelvin reference
output.  The full 60s heat-up validation against a Phase 7 FEM ground
truth is the next-cycle gate (user-selected scope 2026-05-25).
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

_panels_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia", "panels"))
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

import calc_heat_with_em_table as mod  # noqa: E402


class TestCurrentTrajectory:
    def test_constant_value(self):
        I = mod._load_current_trajectory(123.0, "")
        assert I(0.0) == 123.0
        assert I(99.0) == 123.0
        assert I.kind == "constant"

    def test_constant_requires_one(self):
        with pytest.raises(ValueError, match="--coil-current"):
            mod._load_current_trajectory(None, "")

    def test_csv_path_interp_and_clamp(self, tmp_path):
        csv = tmp_path / "I_traj.csv"
        csv.write_text("0.0, 100.0\n30.0, 200.0\n60.0, 50.0\n",
                       encoding="utf-8")
        I = mod._load_current_trajectory(None, str(csv))
        assert I.kind == "csv"
        assert I(0.0) == 100.0
        assert I(30.0) == 200.0
        assert I(60.0) == 50.0
        # Linear interp midpoint.
        assert abs(I(15.0) - 150.0) < 1e-9
        # Clamp below and above.
        assert I(-10.0) == 100.0
        assert I(99999.0) == 50.0

    def test_csv_bad_shape_raises(self, tmp_path):
        csv = tmp_path / "bad.csv"
        csv.write_text("1.0\n2.0\n3.0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="2-col"):
            mod._load_current_trajectory(None, str(csv))


class TestCLI:
    def _run(self, *extra_args):
        cmd = [sys.executable,
               os.path.join(_panels_dir, "calc_heat_with_em_table.py"),
               *extra_args]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=30)

    def test_help_exits_zero(self):
        proc = self._run("--help")
        assert proc.returncode == 0
        assert "--em-table" in proc.stdout
        assert "--ht-source" in proc.stdout
        assert "--I-ref" in proc.stdout

    def test_missing_required_args(self):
        proc = self._run()
        assert proc.returncode != 0
        # argparse writes "the following arguments are required" on stderr.
        assert "required" in proc.stderr.lower()

    def test_biot_path_returns_clear_error(self, tmp_path):
        # We need files to satisfy argparse + the early existence checks
        # before the ht-source dispatch raises "not implemented".
        # Easiest path: build a 3x3 table, point --wp-vol at it (the
        # script returns an error JSON early when --wp-vol doesn't load
        # as a mesh -- that's the error we expect, NOT a Python crash).
        npz = tmp_path / "tab.npz"
        subprocess.run([
            sys.executable,
            os.path.join(_panels_dir, "calc_em_table.py"),
            "--material", "steel", "--frequency", "50e3",
            "--H-grid-min", "1e3", "--H-grid-max", "1e5",
            "--n-H", "3", "--T-grid-min", "20", "--T-grid-max", "300",
            "--n-T", "3", "--em-table-output", str(npz),
        ], capture_output=True, text=True, timeout=60, check=True)

        # Make a 1-line "fake vol" that exists but does not parse as a
        # mesh -- the script must error CLEANLY, not crash with a stack
        # trace at parse time.
        fake_vol = tmp_path / "fake.vol"
        fake_vol.write_text("not_a_real_vol\n", encoding="utf-8")

        proc = self._run(
            "--wp-vol", str(fake_vol),
            "--em-table", str(npz),
            "--ht-source", "biot",
            "--coil-step", "/nonexistent.step",
            "--coil-current", "100",
            "--I-ref", "100")
        # calc_main writes JSON to stdout AND exits non-zero on error
        # (radia >= 4.92.0, commit 5b88b67f, "No Fallbacks -- Fail Fast"
        # policy; soft-fail with returncode 0 was an upstream bug that
        # caused the IGTE 2026 sweep to miscount results).
        assert proc.returncode != 0, \
            "calc_main must exit non-zero on internal error"
        # Last line of stdout is the JSON dump.
        last = proc.stdout.strip().splitlines()[-1]
        import json
        result = json.loads(last)
        assert "error" in result
        # The error message should mention either the missing wp_vol
        # (caught first) OR the biot-not-implemented (if NGSolve loads
        # the bad file and the dispatch is reached).  Both are acceptable
        # -- the contract is "no Python traceback, JSON 'error' key".
        assert isinstance(result["error"], str) and result["error"]
