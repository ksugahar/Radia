"""Required fixture-free strong-coupling CLI contracts (native numerical tier).

Demo-dependent validation is explicitly optional and lives in
optional_coupling_strong.py and optional_coupling_strong_peec.py.
"""
from __future__ import annotations
import pytest
import calc_inductance as ci
from radia import ih_design as ihd


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



def test_strong_accepts_peec_coil():
    """strong now supports the peec coil (CoupledPEECBEMSolver) as well as
    bem-a; the coil-solver guard no longer rejects peec (2026-07-15).  With
    dummy inputs the run still errors -- but DOWNSTREAM (STEP read), not at
    the old ``requires --coil-solver bem-a`` guard.
    """
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "peec", "--frequency", "7000", "--sigma", "5.8e6",
        "--coupling-mode", "strong", "--coil-step", "c.step", "--vol", "w.vol",
        "--no-peec-proximity",
    ])
    # peec passes the coil-solver guard and proceeds to the coil build,
    # which raises on the dummy STEP -- proving it got PAST the guard
    # (the old contract returned an error dict rejecting peec here).
    with pytest.raises(Exception) as exc:
        ci.run_inductance(ns)
    msg = str(exc.value)
    assert "requires --coil-solver bem-a" not in msg
    assert "c.step" in msg or "STEP" in msg



def test_strong_requires_workpiece_vol():
    p = ci.build_argparser()
    ns = p.parse_args([
        "--coil-solver", "bem-a", "--frequency", "7000",
        "--coupling-mode", "strong", "--coil-vol", "c.vol",
    ])
    r = ci.run_inductance(ns)
    assert r.get("status") == "error"
    assert "requires a workpiece" in r["error"]



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
    # wp BIE backend is surfaced (HACApK by default -> scalable workpiece).
    assert "--wp-bem-backend" in cmd
    # coil saddle backend is surfaced so the "HACApK (large)" preset can
    # select loop-COCR for the CoupledBEMSolver coil EFIE.
    assert cmd[cmd.index("--coil-saddle-solver") + 1] == "auto"
    emitted = [c for c in cmd if c.startswith("--")]
    unknown = [f for f in emitted if f not in known]
    assert not unknown, f"strong build_command emits flags argparse rejects: {unknown}"

    large = ihd.IHDesignSpec(
        method=ihd.METHOD_BEMA_BEM_STRONG,
        solver="HACApK (large)",
        coil_vol="coil.vol", wp_vol="wp.vol",
        frequency="7000", current="6700", coil_sigma="5.8e7",
        wp_sigma="5.8e6", mu_r="100", half_thickness="0.005")
    large_cmd = [str(c) for c in large.build_command(
        python="python", panels_dir=None)]
    assert large_cmd[large_cmd.index("--coil-saddle-solver") + 1] == "hacapk_cocr"

    # strong is linear-SIBC only: no ESIM / fes-order / impedance-model.
    vf = spec.visible_fields()
    assert "impedance_model" not in vf
    assert "fes_order" not in vf
    assert {"coil_vol", "wp_vol", "wp_sigma", "mu_r", "half_thickness"} <= vf

    missing = ihd.IHDesignSpec(
        method=ihd.METHOD_BEMA_BEM_STRONG).missing_required_inputs()
    assert "Coil .vol" in missing and "Workpiece .vol" in missing



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
