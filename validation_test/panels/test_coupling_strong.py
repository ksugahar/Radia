"""Contract + end-to-end smoke for the re-exposed BEM-A strong coupling.

Background (2026-07-03): ``bem_coupled_solver.CoupledBEMSolver`` -- the
validated iterative per-DOF back-reaction coil<->workpiece solver
(cross-checked vs FEM-Kelvin SIBC: copper +0.3%, steel mu_r=100 +1.7% on
L, see validation_test/induction_heating/bem_reference) -- had been
orphaned to the validation lane: the production panel/CLI exposed only
the weak one-way Telegen path.  It is now re-exposed as
``calc_inductance.py --coupling-mode strong`` and the notebook-workbench
method ``METHOD_BEMA_BEM_STRONG``.

These tests lock:
  1. argparse accepts ``--coupling-mode strong`` + the coupling knobs.
  2. the fail-fast guards (strong needs bem-a + a workpiece --vol).
  3. the notebook DesignSpec (IHDesignSpec) build_command for the strong
     method emits only flags the calc argparse accepts.
  4. the CLI path runs end-to-end and returns a self-consistent result
     (L_total = L_coil + dL, R_total = R_coil + dR, dR = 2 P_wp / I^2).

The heavy numeric golden (vs FEM) lives in the bem_reference lane; here
the end-to-end is a fast self-consistency smoke on the committed demo.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

# conftest.py puts src/radia, src/radia/panels and validation_test/panels
# on sys.path.
import calc_inductance as ci
from radia import ih_design as ihd

_REPO = Path(__file__).resolve().parents[2]
_SAMPLES = _REPO / "src" / "radia" / "panels" / "samples"
_CALC = _REPO / "src" / "radia" / "panels" / "calc_inductance.py"
_DEMO_COIL_STEP = _SAMPLES / "ih_fem_kelvin_demo_coil.step"
_DEMO_VOL = _SAMPLES / "ih_fem_kelvin_demo.vol"


# ----------------------------------------------------------------------
# 1. argparse surface
# ----------------------------------------------------------------------
def test_argparse_accepts_strong_and_knobs():
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "bem-a", "--frequency", "7000",
        "--coil-vol", "c.vol", "--vol", "w.vol", "--sigma", "5.8e6",
        "--coupling-mode", "strong",
    ])
    assert ns.coupling_mode == "strong"
    assert ns.coupling_max_iter == 10
    assert ns.coupling_tol == 1e-3
    assert ns.coupling_relax == 0.5
    # default stays weak
    ns_w = p.parse_args(["--coil-solver", "bem-a", "--frequency", "7000"])
    assert ns_w.coupling_mode == "weak"
    # strong is a real choice
    choices = [a.choices for a in p._actions if a.dest == "coupling_mode"][0]
    assert set(choices) == {"weak", "strong"}


# ----------------------------------------------------------------------
# 2. fail-fast guards (return before any heavy solve)
# ----------------------------------------------------------------------
def test_strong_requires_bem_a_coil():
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--frequency", "7000", "--sigma", "5.8e6",
        "--coupling-mode", "strong", "--coil-step", "c.step", "--vol", "w.vol",
    ])
    r = ci.run_inductance(ns)
    assert r.get("status") == "error"
    assert "requires --coil-solver bem-a" in r["error"]


def test_strong_requires_workpiece_vol():
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "bem-a", "--frequency", "7000",
        "--coupling-mode", "strong", "--coil-vol", "c.vol",
    ])
    r = ci.run_inductance(ns)
    assert r.get("status") == "error"
    assert "requires a workpiece" in r["error"]


# ----------------------------------------------------------------------
# 3. notebook DesignSpec build_command <-> calc argparse
# ----------------------------------------------------------------------
def test_designspec_strong_build_command_parses():
    p = ci.build_argparser()
    # calc_main adds the shared --output flag at runtime.
    known = {a.option_strings[0] for a in p._actions if a.option_strings}
    known.add("--output")

    assert ihd.METHOD_BEMA_BEM_STRONG in ihd.IH_METHODS
    assert ihd.METHOD_BEMA_BEM_STRONG in ihd.WORKPIECE_METHODS
    assert ihd.METHOD_BEMA_BEM_STRONG in ihd.BEMA_COIL_VOL_METHODS

    spec = ihd.IHDesignSpec(
        method=ihd.METHOD_BEMA_BEM_STRONG,
        coil_vol="coil.vol", wp_vol="wp.vol",
        frequency="7000", current="6700", coil_sigma="5.8e7",
        wp_sigma="5.8e6", mu_r="100", half_thickness="0.005")
    assert spec.coil_solver_cli() == "bem-a"

    cmd = [str(c) for c in spec.build_command(python="python", panels_dir=None)]
    assert cmd[cmd.index("--coupling-mode") + 1] == "strong"
    assert cmd[cmd.index("--coil-solver") + 1] == "bem-a"
    emitted = [c for c in cmd if c.startswith("--")]
    unknown = [f for f in emitted if f not in known]
    assert not unknown, f"strong build_command emits flags argparse rejects: {unknown}"

    # strong is linear-SIBC only: no ESIM / fes-order / impedance-model.
    vf = spec.visible_fields()
    assert "impedance_model" not in vf
    assert "fes_order" not in vf
    assert {"coil_vol", "wp_vol", "wp_sigma", "mu_r", "half_thickness"} <= vf

    missing = ihd.IHDesignSpec(
        method=ihd.METHOD_BEMA_BEM_STRONG).missing_required_inputs()
    assert "Coil .vol" in missing and "Workpiece .vol" in missing


# ----------------------------------------------------------------------
# 4. end-to-end self-consistency smoke on the committed demo
# ----------------------------------------------------------------------
@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo coil STEP / workpiece .vol fixtures not present")
def test_strong_end_to_end_self_consistent(tmp_path):
    """Run --coupling-mode strong on the demo (coil STEP -> coil.vol,
    workpiece = sibc hole) and assert the coupled solve is self-consistent.

    Values are tiny (the demo geometry is small/coarse); the physics
    magnitude golden is in the bem_reference lane.  Here we lock the
    wiring: status ok, coupling_mode strong, converged, and the output
    assembly identities L_total = L_coil + dL, R_total = R_coil + dR,
    dR = 2 P_wp / I^2.
    """
    from _bema_coil_vol_helper import coil_vol_for

    coil_vol = coil_vol_for(str(_DEMO_COIL_STEP), cache_dir=str(tmp_path))
    out_json = tmp_path / "strong.json"
    current = 1.0
    cmd = [
        sys.executable, str(_CALC),
        "--coil-solver", "bem-a", "--coupling-mode", "strong",
        "--coil-vol", coil_vol,
        "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
        "--frequency", "7000", "--current", str(current),
        "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
        "--coupling-max-iter", "4", "--coupling-tol", "5e-3",
        "--output", str(out_json),
    ]
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=1200, env=env)
    assert proc.returncode == 0, f"calc failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
    assert out_json.is_file(), f"no output json; stdout tail:\n{proc.stdout[-2000:]}"

    d = json.loads(out_json.read_text(encoding="utf-8"))
    assert d.get("status") == "ok", d
    assert d["coupling_mode"] == "strong"
    assert d["method"] == "bem-a-bem-strong"
    assert int(d["coupling_iterations"]) >= 1

    # Output-assembly identities.
    assert math.isclose(d["L_total_nH"], d["L_coil_nH"] + d["delta_L_nH"],
                        rel_tol=0, abs_tol=1e-9)
    assert math.isclose(d["R_total_mOhm"], d["R_coil_mOhm"] + d["delta_R_mOhm"],
                        rel_tol=0, abs_tol=1e-9)
    # dR = 2 P_wp / I^2  (P_wp = 1/2 I^2 dR).
    assert math.isclose(d["delta_R_mOhm"],
                        2.0 * d["P_wp_W"] / (current * current) * 1e3,
                        rel_tol=1e-6, abs_tol=1e-12)
    assert d["P_wp_W"] >= 0.0
    assert math.isfinite(d["delta_L_nH"])

    # The re-exposed path must surface the coupled-solver context keys.
    for key in ("coupled_L_air_nH", "coupled_L_total_nH", "wp_ndof",
                "coupled_n_J_coil", "H_t_rms_A_per_m", "t_coupled_solve_s"):
        assert key in d, f"missing coupled key {key!r}"
