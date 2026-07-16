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


@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo coil STEP / workpiece .vol fixtures not present")
def test_strong_output_scales_with_terminal_current(tmp_path):
    """P_wp ~ I^2, H_t ~ I; dL and dR current-independent.

    Regression for the Kubota 2026-07-16 incident: CoupledBEMSolver drives
    the coil EFIE at UNIT terminal current (its solve() takes no current),
    and the strong driver forgot to rescale -- so at 6700 A the panel
    reported P_wp = 8.3e-4 W / H_t = 9.84 A/m (1 A-drive values, I^2 =
    4.5e7x low) while Delta_L came out sane and masked the bug.  The
    problem is linear, so run the demo at I = 1 A and I = 5 A and require
    exact linear-scaling ratios.
    """
    from _bema_coil_vol_helper import coil_vol_for

    coil_vol = coil_vol_for(str(_DEMO_COIL_STEP), cache_dir=str(tmp_path))
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")

    def _run(current):
        out_json = tmp_path / f"strong_I{current:g}.json"
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
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=1200, env=env)
        assert proc.returncode == 0, \
            f"[I={current}]\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
        return json.loads(out_json.read_text(encoding="utf-8"))

    d1 = _run(1.0)
    d5 = _run(5.0)

    assert d1["P_wp_W"] > 0
    # Dissipation scales as I^2, fields as I.
    assert math.isclose(d5["P_wp_W"], 25.0 * d1["P_wp_W"], rel_tol=1e-6)
    assert math.isclose(d5["H_t_rms_A_per_m"], 5.0 * d1["H_t_rms_A_per_m"],
                        rel_tol=1e-6)
    # Inductance / resistance deltas are current-independent.
    assert math.isclose(d5["delta_L_nH"], d1["delta_L_nH"],
                        rel_tol=1e-6, abs_tol=1e-12)
    assert math.isclose(d5["delta_R_mOhm"], d1["delta_R_mOhm"],
                        rel_tol=1e-6, abs_tol=1e-12)
    # And the energy identity holds at BOTH currents.
    for d, I in ((d1, 1.0), (d5, 5.0)):
        assert math.isclose(d["delta_R_mOhm"],
                            2.0 * d["P_wp_W"] / (I * I) * 1e3,
                            rel_tol=1e-6, abs_tol=1e-12)


@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo coil STEP / workpiece .vol fixtures not present")
def test_strong_wp_hacapk_matches_dense(tmp_path):
    """Q1: the workpiece HACApK backend reproduces the dense result.

    CoupledBEMSolver's workpiece BIE gained an O(N log N) intree-HACApK
    backend (--wp-bem-backend hacapk, the default) so the coupled solve
    scales past the ~12k-tri dense-assembly wall (e.g. the 20k-tri
    Takahashi workpiece).  On the demo it must agree with the dense path
    (--wp-bem-backend intree-dense) to within the ACA compression accuracy.
    """
    from _bema_coil_vol_helper import coil_vol_for

    coil_vol = coil_vol_for(str(_DEMO_COIL_STEP), cache_dir=str(tmp_path))

    def _run(backend):
        out_json = tmp_path / f"strong_{backend}.json"
        cmd = [
            sys.executable, str(_CALC),
            "--coil-solver", "bem-a", "--coupling-mode", "strong",
            "--coil-vol", coil_vol,
            "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
            "--frequency", "7000", "--current", "1",
            "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
            "--coupling-max-iter", "4", "--coupling-tol", "5e-3",
            "--wp-bem-backend", backend,
            "--output", str(out_json),
        ]
        env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=1200, env=env)
        assert proc.returncode == 0, \
            f"[{backend}]\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
        return json.loads(out_json.read_text(encoding="utf-8"))

    dense = _run("intree-dense")
    hac = _run("hacapk")
    assert dense["wp_bem_backend"] == "intree-dense"
    assert hac["wp_bem_backend"] == "hacapk"

    # Same physics to well within the HACApK ACA compression accuracy
    # (measured worst rel.diff ~7e-4 on the demo's tiny P_wp).
    for k in ("L_total_nH", "coupled_L_total_nH", "R_total_mOhm",
              "P_wp_W", "H_t_rms_A_per_m"):
        a, b = float(dense[k]), float(hac[k])
        rel = abs(a - b) / max(abs(a), 1e-30)
        assert rel < 5e-3, f"{k}: dense={a:.6e} hacapk={b:.6e} rel={rel:.2e}"


@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo coil STEP / workpiece .vol fixtures not present")
def test_strong_coil_hacapk_matches_dense(tmp_path):
    """Q1: the coil HACApK (loop-COCR) saddle backend reproduces dense LU.

    CoupledBEMSolver's coil EFIE gained an O(N log N) H-matrix + loop-COCR
    saddle solve (--coil-saddle-solver hacapk_cocr; O(N r) storage, built
    once and re-solved each Picard iteration).  Holding the workpiece
    backend fixed (intree-dense), the coil H-matrix path must reproduce the
    dense-LU coil solve to loop-COCR tolerance (~1e-8).
    """
    from _bema_coil_vol_helper import coil_vol_for

    coil_vol = coil_vol_for(str(_DEMO_COIL_STEP), cache_dir=str(tmp_path))

    def _run(coil_saddle):
        out_json = tmp_path / f"strong_coil_{coil_saddle}.json"
        cmd = [
            sys.executable, str(_CALC),
            "--coil-solver", "bem-a", "--coupling-mode", "strong",
            "--coil-vol", coil_vol,
            "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
            "--frequency", "7000", "--current", "1",
            "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
            "--coupling-max-iter", "4", "--coupling-tol", "5e-3",
            "--wp-bem-backend", "intree-dense",   # hold wp fixed -> isolate coil
            "--coil-saddle-solver", coil_saddle,
            "--output", str(out_json),
        ]
        env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=1200, env=env)
        assert proc.returncode == 0, \
            f"[{coil_saddle}]\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
        return json.loads(out_json.read_text(encoding="utf-8"))

    dense = _run("auto")            # auto -> dense LU on this small coil
    hac = _run("hacapk_cocr")
    assert dense["coil_bem_backend"] == "dense-lu"
    assert hac["coil_bem_backend"] == "hacapk_cocr"

    # The coil H-matrix matvec (aca_eps 1e-10) is tighter than the loop-COCR
    # tolerance, so the coil solve matches dense LU to ~1e-8 (measured ~4e-11).
    for k in ("L_coil_nH", "L_total_nH", "coupled_L_total_nH", "R_total_mOhm",
              "P_wp_W", "H_t_rms_A_per_m"):
        a, b = float(dense[k]), float(hac[k])
        rel = abs(a - b) / max(abs(a), 1e-30)
        assert rel < 1e-6, f"{k}: dense={a:.6e} hacapk={b:.6e} rel={rel:.2e}"


@pytest.mark.skipif(
    not (_DEMO_COIL_STEP.is_file() and _DEMO_VOL.is_file()),
    reason="demo coil STEP / workpiece .vol fixtures not present")
def test_strong_accepts_volume_coil_vol(tmp_path):
    """A VOLUME coil .vol (tets present, the common Cubit export) must run.

    Regression for the Takahashi coil_only.vol failure (2026-07-16, same
    class as keiko gapped_torus 2026-05-12): the strong driver loaded
    ``Mesh(args.coil_vol)`` raw, so on a volume .vol HDivSurface picked up
    thousands of extra null modes (n_J 17,799 instead of 6,204 on
    Takahashi), the coupled saddle LU went singular, and the solve died
    with "array must not contain infs or NaNs".  The driver now routes
    through ``_build_bema_coil_mesh`` (label validation + volume->surface
    extraction), the same loader as the vacuum BEM-A path.
    """
    from _bema_coil_vol_helper import step_to_coil_vol

    coil_vol = str(tmp_path / "coil_volume.vol")
    step_to_coil_vol(str(_DEMO_COIL_STEP), coil_vol, volume_mesh=True)

    # Confirm the fixture actually has volume elements (else this test
    # silently degenerates to the surface case).
    from ngsolve import Mesh
    assert Mesh(coil_vol).ne > 0, "fixture must contain volume tets"

    out_json = tmp_path / "strong_volvol.json"
    cmd = [
        sys.executable, str(_CALC),
        "--coil-solver", "bem-a", "--coupling-mode", "strong",
        "--coil-vol", coil_vol,
        "--vol", str(_DEMO_VOL), "--wp-label", "sibc",
        "--frequency", "7000", "--current", "1",
        "--sigma", "5.8e6", "--mu-r", "100", "--half-thickness", "0.005",
        "--coupling-max-iter", "4", "--coupling-tol", "5e-3",
        "--output", str(out_json),
    ]
    env = dict(os.environ, MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=1200, env=env)
    assert proc.returncode == 0, \
        f"volume coil .vol strong run failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
    d = json.loads(out_json.read_text(encoding="utf-8"))
    assert d.get("status") == "ok", d
    assert math.isfinite(d["P_wp_W"]) and d["P_wp_W"] >= 0.0
    assert math.isfinite(d["delta_L_nH"])
