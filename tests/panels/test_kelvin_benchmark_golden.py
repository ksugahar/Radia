"""Golden test for the EM panel's Kelvin Benchmark mode.

End-to-end:
  1. Build EMPanel, switch to "Kelvin Benchmark" formulation
  2. Set knobs to the canonical values (mu_r=100, H0=1.0, axis=z,
     fes_order=2, R_kelvin=0.20)
  3. Build the CLI command with the bundled sample .vol
  4. Run the CLI as subprocess
  5. Assert: |error_pct| within the per-fraction golden band
     (verified 2026-04-26 -- see VARIANTS table below)

This is the **panel-level** companion to
`tests/cubit/test_kelvin_1_4_p_convergence.py` (which exercises the
mesh+solver scripts directly).  This test exercises the FULL panel
chain: widget values -> build_command -> subprocess -> JSON output.

Variants:
  1_2 (half model, y >= 0):     +1.07% at p=2 -> band 1.5%
  1_4 (quarter model, x,y >= 0): +0.71% at p=2 -> band 1.5%

The 1/8 variant is intentionally NOT covered: the kelvin_far
Dirichlet plane passes through r'=0 (Kelvin centre) where the
reluctivity Mu = mu0*(R/r')^2 is singular, and the linear solve
gives Hz ~ 0 instead of analytical.  Tracked in
`memory/feedback_kelvin_1_8_blocker.md`.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[2]
SAMPLES_DIR = REPO / "src" / "radia" / "panels" / "samples"


# (variant_id, sample_file, golden_band_pct)
VARIANTS = [
    ("1_2", "kelvin_benchmark_sphere_1_2.vol", 1.5),
    ("1_4", "kelvin_benchmark_sphere_1_4.vol", 1.5),
]


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
    panel._widgets["method"].setCurrentText(FORM_KELVIN_BENCH)
    panel._on_method_changed(FORM_KELVIN_BENCH)
    sub = panel._sub_panels[FORM_KELVIN_BENCH]
    sub._widgets["mu_r"].setText("100")
    sub._widgets["H0"].setText("1.0")
    sub._widgets["field_axis"].setCurrentText("z")
    sub._widgets["fes_order"].setValue(2)
    sub._widgets["R_kelvin"].setText("0.20")
    yield panel


@pytest.mark.parametrize("variant,sample_name,band_pct", VARIANTS,
                         ids=[v[0] for v in VARIANTS])
def test_sample_vol_present(variant, sample_name, band_pct):
    """The packaged Kelvin Benchmark sample .vol must ship with the wheel."""
    sample = SAMPLES_DIR / sample_name
    if not sample.exists():
        pytest.skip(
            f"Sample {sample_name} not built locally.  Run "
            f"`python src/radia/panels/samples/_build_all_kelvin_benchmarks.py "
            f"--fracs {variant}` to generate (requires Cubit).  CI / wheel "
            f"build pulls the file from the GitHub Releases binaries tag "
            f"(per CLAUDE.md Binary File Policy).")


def test_panel_runnable_with_sample(kelvin_panel):
    """Panel must report runnable=True when Kelvin Benchmark mode is active."""
    assert kelvin_panel.is_runnable(), \
        "Kelvin Benchmark needs no coil_script -- panel should be runnable"


@pytest.mark.parametrize("variant,sample_name,band_pct", VARIANTS,
                         ids=[v[0] for v in VARIANTS])
def test_panel_build_command_minimal(kelvin_panel, variant, sample_name, band_pct):
    """Verify all expected CLI flags are emitted with sane values."""
    sample = SAMPLES_DIR / sample_name
    cmd = kelvin_panel.build_command(str(sample))
    cmd_str = " ".join(cmd)
    for required_flag in (
            "calc_kelvin_benchmark.py",
            "--vol", str(sample),
            "--fes-order", "2",
            "--mu-r", "100",
            "--H0", "1.0",
            "--field-axis", "z",
            "--R-kelvin", "0.20"):
        assert required_flag in cmd_str, \
            f"build_command missing {required_flag!r}: {cmd_str}"


@pytest.mark.parametrize("variant,sample_name,band_pct", VARIANTS,
                         ids=[v[0] for v in VARIANTS])
def test_panel_end_to_end_golden(kelvin_panel, tmp_path,
                                  variant, sample_name, band_pct):
    """Run the panel-built CLI on the bundled sample; lock the numerics."""
    sample = SAMPLES_DIR / sample_name
    if not sample.exists():
        pytest.skip(f"Sample {sample_name} not built locally.")
    cmd = kelvin_panel.build_command(str(sample))
    out_json = tmp_path / f"kelvin_bench_{variant}_out.json"
    out_msh = tmp_path / f"kelvin_bench_{variant}_out.msh"
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
    assert err_pct < band_pct, (
        f"[{variant}] Kelvin Benchmark error {err_pct:.2f}% exceeds "
        f"{band_pct}% bound; "
        f"Hi_origin={data['Hi_origin']:.6e}, "
        f"Hi_analytical={data['Hi_analytical']:.6e}")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
