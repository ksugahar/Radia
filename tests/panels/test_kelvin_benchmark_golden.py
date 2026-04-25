"""Golden test for the EM panel's Kelvin Benchmark mode.

End-to-end:
  1. Build EMPanel, switch to "Kelvin Benchmark" formulation
  2. Set knobs to the canonical values (mu_r=100, H0=1.0, axis=z,
     fes_order=2, R_kelvin=0.20)
  3. Build the CLI command with the bundled sample .vol
  4. Run the CLI as subprocess
  5. Assert: error_pct within +/-1.5% (matches Cubit_1_4_p_convergence
     verification)

This is the **panel-level** companion to
`tests/cubit/test_kelvin_1_4_p_convergence.py` (which exercises the
mesh+solver scripts directly).  This test exercises the FULL panel
chain: widget values -> build_command -> subprocess -> JSON output.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[2]
SAMPLE_VOL = (REPO / "src" / "radia" / "panels" / "samples"
              / "kelvin_benchmark_sphere_1_4.vol")


@pytest.fixture(scope="module")
def panel_module():
    sys.path.insert(0, str(REPO / "src" / "radia"))
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    yield app


@pytest.fixture(scope="module")
def kelvin_panel(panel_module):
    sys.path.insert(0, str(REPO / "src" / "radia"))
    from radia_em import EMPanel, FORM_KELVIN_BENCH
    panel = EMPanel()
    panel._widgets["formulation"].setCurrentText(FORM_KELVIN_BENCH)
    panel._on_formulation_changed(FORM_KELVIN_BENCH)
    panel._widgets["mu_r"].setText("100")
    panel._widgets["H0"].setText("1.0")
    panel._widgets["field_axis"].setCurrentText("z")
    panel._widgets["fes_order"].setValue(2)
    panel._widgets["R_kelvin"].setText("0.20")
    yield panel


def test_sample_vol_present():
    """The packaged Kelvin Benchmark sample .vol must ship with the wheel."""
    if not SAMPLE_VOL.exists():
        pytest.skip(
            f"Sample {SAMPLE_VOL.name} not built locally.  Run "
            f"`python src/radia/panels/samples/kelvin_benchmark_sphere_1_4_build.py "
            f"--out-dir src/radia/panels/samples --orders 2` to generate "
            f"(requires Cubit).  CI / wheel build pulls this file from "
            f"the GitHub Releases binaries tag (per CLAUDE.md Binary File "
            f"Policy).")


def test_panel_runnable_with_sample(kelvin_panel):
    """Panel must report runnable=True when Kelvin Benchmark mode is active."""
    assert kelvin_panel.is_runnable(), \
        "Kelvin Benchmark needs no coil_script -- panel should be runnable"


def test_panel_build_command_minimal(kelvin_panel):
    """Verify all expected CLI flags are emitted with sane values."""
    cmd = kelvin_panel.build_command(str(SAMPLE_VOL))
    cmd_str = " ".join(cmd)
    for required_flag in (
            "calc_kelvin_benchmark.py",
            "--vol", str(SAMPLE_VOL),
            "--fes-order", "2",
            "--mu-r", "100",
            "--H0", "1.0",
            "--field-axis", "z",
            "--R-kelvin", "0.20"):
        assert required_flag in cmd_str, \
            f"build_command missing {required_flag!r}: {cmd_str}"


def test_panel_end_to_end_golden(kelvin_panel, tmp_path):
    """Run the panel-built CLI on the bundled sample; lock the numerics.

    Tolerance band (validated 2026-04-25 for fes_order=2):
        analytical Hz_inside = 3 / (mu_r + 2) * H0 = 0.029412
        |error_pct| <= 1.5%   (Cubit_1_4_p_convergence reports 0.71%)
    """
    if not SAMPLE_VOL.exists():
        pytest.skip(f"Sample {SAMPLE_VOL.name} not built locally.")
    cmd = kelvin_panel.build_command(str(SAMPLE_VOL))
    # Redirect outputs into tmp_path so we don't pollute the repo
    out_json = tmp_path / "kelvin_bench_out.json"
    out_msh = tmp_path / "kelvin_bench_out.msh"
    new_cmd = []
    skip = 0
    for tok in cmd:
        if skip > 0:
            skip -= 1
            continue
        if tok == "--output":
            new_cmd += ["--output", str(out_json)]; skip = 1
            continue
        if tok == "--msh-output":
            new_cmd += ["--msh-output", str(out_msh)]; skip = 1
            continue
        new_cmd.append(tok)
    res = subprocess.run(new_cmd, capture_output=True, text=True, cwd=str(REPO))
    assert res.returncode == 0, (
        f"calc CLI failed (exit {res.returncode})\n"
        f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}")
    assert out_json.exists(), \
        f"Expected --output JSON at {out_json}; stdout:\n{res.stdout}"
    data = json.loads(out_json.read_text())
    assert data.get("converged") is True, f"solver did not converge: {data}"
    err_pct = abs(data["error_pct"])
    assert err_pct < 1.5, (
        f"Kelvin Benchmark error {err_pct:.2f}% exceeds 1.5% bound; "
        f"Hi_origin={data['Hi_origin']:.6e}, "
        f"Hi_analytical={data['Hi_analytical']:.6e}")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
