"""Contract + self-consistency smoke for PEEC-coil strong coupling.

The weak Telegen path evaluates workpiece heating from the INCIDENT
(bare-coil) field, so it can over-estimate a strongly loaded magnetic
workpiece.  The BEM-A coil has a strong (self-consistent) path via
``CoupledBEMSolver``; the PEEC coil shares the SAME incident-field weak
path, so its P_wp is over by the same factor.  ``CoupledPEECBEMSolver``
gives the PEEC coil the same strong coupling (filament back-EMF from the
workpiece reaction redistributes the loop-bundle currents).  It is exposed
as ``calc_inductance.py --coil-solver peec --coupling-mode strong``.

These tests lock the WIRING + SELF-CONSISTENCY only:
  1. argparse accepts ``--coil-solver peec --coupling-mode strong``.
  2. the strong guard now allows peec (and still needs a workpiece --vol).
  3. the CLI runs end-to-end and returns self-consistent output
     (L_total = L_coil + dL, R_total = R_coil + dR, dR = 2 P_wp / I^2),
     the peec-loop-bundle backend label, and the EXPERIMENTAL flag.
  4. on the (weakly-coupled) demo, strong P_wp ~= weak P_wp -- locks that
     the coupled solve reduces to the weak forward when the workpiece
     barely loads the coil (and does not diverge).

The strong-loading P_wp response is NOT validated here:
the committed demo is weakly coupled (BEM-A weak == strong to ~2%), so it
cannot exercise the coil-current redistribution.  That validation needs
a durable strongly-coupled reference case on a compute host -- see the
``radia.peec_coupled_bem_solver`` module docstring (VALIDATION STATUS).
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


def _run_cli(mode, tmp_path, extra=None):
    out_json = tmp_path / f"peec_{mode}.json"
    cmd = [
        sys.executable, str(_CALC),
        "--coil-solver", "peec", "--coupling-mode", mode,
        "--coil-step", str(_DEMO_COIL_STEP),
        "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
        "--frequency", "7000", "--current", "1",
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
    # R is taken from energy: dR = 2 P_wp / I^2.
    assert math.isclose(d["delta_R_mOhm"], 2.0 * d["P_wp_W"] / 1.0 * 1e3,
                        rel_tol=1e-6, abs_tol=1e-12)
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
    """On the weakly-coupled demo, strong P_wp ~= weak P_wp.

    The demo workpiece barely loads the coil (BEM-A weak == strong to
    ~2%), so the coil-current redistribution is negligible and the strong
    forward must reduce to the weak one.  This locks (a) the strong
    forward is the validated weak forward and (b) the coupled solve does
    not diverge / blow up on a weakly-coupled input.  It does NOT validate
    the strong REDUCTION (needs a strongly-coupled case; see module doc).
    """
    weak = _run_cli("weak", tmp_path)
    strong = _run_cli("strong", tmp_path,
                      extra=["--coupling-max-iter", "4", "--coupling-tol", "5e-3"])
    assert weak["P_wp_W"] > 0 and strong["P_wp_W"] > 0
    rel = abs(strong["P_wp_W"] - weak["P_wp_W"]) / weak["P_wp_W"]
    assert rel < 0.05, (
        f"strong P_wp={strong['P_wp_W']:.4e} vs weak P_wp={weak['P_wp_W']:.4e} "
        f"differ by {rel:.1%} on a weakly-coupled demo (expected <5%)")
