"""PEEC scalar strong-coupling CLI and port-power accounting.

The coupled route redistributes filament currents through complete body
reaction. Body heat and reaction are checked independently. The change
in port loss includes body reaction AND the change in coil self loss.
Current scaling and finite-mesh reaction checks use self-authored cases
in run_strong_reaction_ring.py; this file also checks legacy CLI wiring.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

import calc_inductance as ci

_REPO = Path(__file__).resolve().parents[2]
_SAMPLES = _REPO / "src" / "radia" / "panels" / "samples"
_CALC = _SAMPLES.parent / "calc_inductance.py"
_DEMO_COIL_STEP = _SAMPLES / "ih_fem_kelvin_demo_coil.step"
_DEMO_VOL = _SAMPLES / "ih_fem_kelvin_demo.vol"


# ----------------------------------------------------------------------
# 1. argparse + guard surface
# ----------------------------------------------------------------------
def test_argparse_accepts_peec_strong():
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--frequency", "7000",
        "--coil-step", "c.step", "--vol", "w.vol", "--sigma", "5.8e6",
        "--coupling-mode", "strong",
    ])
    assert ns.coil_solver == "peec"
    assert ns.coupling_mode == "strong"


def test_strong_peec_requires_workpiece_vol():
    """strong needs a workpiece --vol even for the peec coil."""
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--frequency", "7000",
        "--coil-step", "c.step", "--sigma", "5.8e6",
        "--coupling-mode", "strong", "--coil-only",
    ])
    out = ci.run_inductance(ns)
    assert out.get("status") == "error"
    assert "workpiece" in out["error"] and "--vol" in out["error"]


def test_strong_peec_rejects_proximity_model_mismatch():
    """Strong PEEC must not mix proximity and isolated-wire baselines."""
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--frequency", "7000",
        "--coil-step", "c.step", "--vol", "w.vol", "--sigma", "5.8e6",
        "--coupling-mode", "strong",
    ])
    out = ci.run_inductance(ns)
    assert out.get("status") == "error"
    assert "--no-peec-proximity" in out["error"]


# ----------------------------------------------------------------------
# 2. end-to-end self-consistency on the committed demo
# ----------------------------------------------------------------------
_SKIP = not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file())


def _run_cli(mode, tmp_path, extra=None, current=1.0):
    out_json = tmp_path / f"peec_{mode}_I{current:g}.json"
    cmd = [
        sys.executable, str(_CALC),
        "--coil-solver", "peec", "--coupling-mode", mode,
        "--coil-step", str(_DEMO_COIL_STEP),
        "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
        "--frequency", "7000", "--current", str(current),
        "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
        "--no-peec-proximity",
        "--output", str(out_json),
    ] + (extra or [])
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=1200, env=env)
    assert proc.returncode == 0, \
        f"[{mode}]\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
    assert out_json.is_file(), f"[{mode}] no json:\n{proc.stdout[-2000:]}"
    return json.loads(out_json.read_text(encoding="utf-8"))


@pytest.mark.skipif(_SKIP, reason="demo coil STEP / workpiece .vol not present")
def test_peec_strong_end_to_end_self_consistent(tmp_path):
    d = _run_cli("strong", tmp_path,
                 extra=["--coupling-max-iter", "4", "--coupling-tol", "5e-3"])
    assert d.get("status") == "ok", d
    assert d["coupling_mode"] == "strong"
    assert d["method"] == "peec-bem-strong"
    assert d["coil_bem_backend"] == "peec-loop-bundle"
    assert int(d["coupling_iterations"]) >= 1
    assert d["coupling_converged"] is True
    assert d["coupling_residual"] <= 5e-3

    # Output-assembly identities (hold by construction, like the BEM-A path).
    assert math.isclose(d["L_total_nH"], d["L_coil_nH"] + d["delta_L_nH"],
                        rel_tol=0, abs_tol=1e-9)
    assert math.isclose(d["R_total_mOhm"], d["R_coil_mOhm"] + d["delta_R_mOhm"],
                        rel_tol=0, abs_tol=1e-9)
    # A changed filament-current distribution changes coil self loss too.
    assert math.isclose(d["delta_R_mOhm"],
                        2.0*(d["body_reaction_power_W"]+d["coil_loss_change_W"])*1e3,
                        rel_tol=1e-6, abs_tol=1e-12)
    assert math.isclose(d["port_power_W"],
                        d["body_reaction_power_W"]+d["coil_loss_W"], rel_tol=1e-6)
    assert d["body_power_balance_relative_error"] <= .1
    assert d["body_residual"] <= 1e-6
    assert d["P_wp_W"] >= 0.0
    assert math.isfinite(d["delta_L_nH"])

    # EXPERIMENTAL flag + diagnostic reaction-R must be surfaced.
    assert d.get("experimental") is True
    assert "experimental_note" in d
    assert "coupled_delta_R_reaction_mOhm" in d
    assert d["coupled_n_filaments"] >= 1

    for key in ("coupled_L_air_nH", "coupled_L_total_nH", "wp_ndof",
                "H_t_rms_A_per_m", "t_coupled_solve_s"):
        assert key in d, f"missing coupled key {key!r}"


@pytest.mark.skipif(_SKIP, reason="demo coil STEP / workpiece .vol not present")
def test_peec_strong_reduces_to_weak_on_weak_coupling(tmp_path):
    """Legacy fixture compatibility check; no general strong/weak error bound."""
    weak = _run_cli("weak", tmp_path)
    strong = _run_cli("strong", tmp_path,
                      extra=["--coupling-max-iter", "4", "--coupling-tol", "5e-3"])
    assert weak["P_wp_W"] > 0 and strong["P_wp_W"] > 0
    rel = abs(strong["P_wp_W"] - weak["P_wp_W"]) / weak["P_wp_W"]
    assert rel < 0.05, (
        f"strong P_wp={strong['P_wp_W']:.4e} vs weak P_wp={weak['P_wp_W']:.4e} "
        f"differ by {rel:.1%} on a weakly-coupled demo (expected <5%)")


@pytest.mark.skipif(_SKIP, reason="demo coil STEP / workpiece .vol not present")
def test_peec_strong_output_scales_with_terminal_current(tmp_path):
    """P_wp ~ I^2, H_t ~ I; dL and dR current-independent (PEEC driver).

    Companion of test_coupling_strong.py::
    test_strong_output_scales_with_terminal_current (the BEM-A driver,
    which forgot the unit-current rescale -- Kubota 2026-07-16, P_wp
    4.5e7x low at 6700 A).  The PEEC driver passes I_port to its solver,
    so it was already real-current-scaled; this locks the shared
    "strong drivers return REAL-current-scaled P_total / H_t" contract
    on the PEEC side too.
    """
    knobs = ["--coupling-max-iter", "4", "--coupling-tol", "5e-3"]
    d1 = _run_cli("strong", tmp_path, extra=knobs, current=1.0)
    d5 = _run_cli("strong", tmp_path, extra=knobs, current=5.0)

    assert d1["P_wp_W"] > 0
    assert math.isclose(d5["P_wp_W"], 25.0 * d1["P_wp_W"], rel_tol=1e-6)
    assert math.isclose(d5["H_t_rms_A_per_m"], 5.0 * d1["H_t_rms_A_per_m"],
                        rel_tol=1e-6)
    assert math.isclose(d5["delta_L_nH"], d1["delta_L_nH"],
                        rel_tol=1e-6, abs_tol=1e-12)
    assert math.isclose(d5["delta_R_mOhm"], d1["delta_R_mOhm"],
                        rel_tol=1e-6, abs_tol=1e-12)
