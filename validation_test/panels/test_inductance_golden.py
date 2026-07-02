"""Golden regression tests for calc_inductance.py.

Consolidated test for the unified inductance CLI which replaced three
older scoped CLIs in 2026-05-02 (Refactor A1+B+C):
  * tests for calc_peec_inductance.py        -> tests below for --coil-solver peec
  * tests for calc_peec_bem.py               -> tests below for --coil-solver peec --vol
  * tests for calc_coil_bem_a_workpiece.py   -> tests below for --coil-solver bem-a --vol

The unified CLI dispatches on:
  * ``--coil-solver peec | bem-a``   (REQUIRED)
  * presence/absence of ``--vol``     (vacuum vs weak-coupled)

giving 4 modes; this file has 4 corresponding @pytest.mark.slow tests.

Coupling terminology (corrected 2026-05-02): the old "1-way forward"
terminology was misleading because Telegen Δ L IS computed (the
back-reaction is captured at the port level via reciprocity, even
though coil J is fixed).  This is **weak coupling**, not strict
one-way; the latter would skip Δ L entirely.  See calc_inductance.py
docstring.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CALC = os.path.join(ROOT, "src", "radia", "panels", "calc_inductance.py")


def _load_golden(fname: str) -> dict:
    with open(os.path.join(HERE, "golden", fname), "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_sample(path: str) -> str:
    if os.path.isabs(path):
        return path
    return os.path.join(ROOT, path)


def _run_calc(args_list, timeout_s: int = 900) -> dict:
    """Run calc_inductance.py with the supplied arg list, return JSON."""
    proc = subprocess.run([sys.executable, CALC] + args_list,
                          capture_output=True, text=True,
                          timeout=timeout_s)
    assert proc.returncode == 0, (
        f"calc_inductance returned {proc.returncode}\n"
        f"CMD: {[CALC] + args_list}\n"
        f"STDOUT (tail):\n{proc.stdout[-4000:]}\n"
        f"STDERR (tail):\n{proc.stderr[-2000:]}")
    json_lines = [ln for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("{")]
    assert json_lines, f"no JSON in stdout:\n{proc.stdout}"
    return json.loads(json_lines[-1])


def _assert_in_range(value, lo, hi, label, hint=""):
    assert lo <= value <= hi, (
        f"{label}: {value:.4g} outside regression guard [{lo}, {hi}].  "
        f"{hint}")


def _assert_close(value, expected, tol_pct, label):
    err_pct = abs(value - expected) / max(abs(expected), 1e-30) * 100.0
    assert err_pct <= tol_pct, (
        f"{label}: got {value:.4g}, expected {expected:.4g} "
        f"(err {err_pct:.2f}%, tol {tol_pct}%)")


# ============================================================================
# Mode 1: --coil-solver peec, vacuum (--vol omitted)
# ============================================================================

@pytest.mark.slow
def test_inductance_peec_vacuum_gapped_torus():
    """355deg gapped torus, PEEC vacuum.  Path: TORUS analytical sweep."""
    g = _load_golden("peec_inductance_torus_50kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        pytest.skip(
            f"Sample STEP missing: {step}.  Regenerate via Cubit "
            f"playback of {g['sample']['generator_jou']}.")

    phys = g["physics"]
    result = _run_calc([
        "--coil-step", step,
        "--coil-solver", "peec",
        "--peec-n-peri", str(phys["n_peri"]),
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", str(phys["sigma_S_per_m"]),
    ])
    assert result["status"] == "ok"
    assert result["coil_solver_used"] == "peec"

    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                     guard["max_L_nH"], "L_coil_nH (hard band)",
                     hint=guard["comment"])
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")
    _assert_close(result["R_coil_mOhm"], exp["R_coil_mOhm"],
                  exp["tolerance_R_pct"], "R_coil_mOhm (golden)")


@pytest.mark.slow
def test_inductance_peec_vacuum_3turn_loft():
    """3-turn pancake loft, PEEC vacuum.  Path: cross-section centroid chain."""
    g = _load_golden("peec_inductance_3turn_150kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        if g["sample"].get("skip_if_missing", False):
            pytest.skip(f"Sample STEP not available: {step}")
        pytest.fail(f"Sample STEP missing: {step}")

    phys = g["physics"]
    # Copy STEP to an ASCII tmp path with NO sibling .jou (sibling auto-
    # preference must not interfere; we exercise the STEP-only path).
    with tempfile.TemporaryDirectory() as td:
        isolated_step = os.path.join(td, "coil.stp")
        shutil.copy2(step, isolated_step)
        result = _run_calc([
            "--coil-step", isolated_step,
            "--coil-solver", "peec",
            "--peec-n-peri", str(phys["n_peri"]),
            "--frequency", str(phys["frequency_Hz"]),
            "--current", str(phys["current_A"]),
            "--coil-sigma", str(phys["sigma_S_per_m"]),
        ])

    assert result["status"] == "ok"
    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                     guard["max_L_nH"], "L_coil_nH (hard band)",
                     hint=guard["comment"])
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")
    _assert_close(result["R_coil_mOhm"], exp["R_coil_mOhm"],
                  exp["tolerance_R_pct"], "R_coil_mOhm (golden)")


@pytest.mark.slow
def test_inductance_peec_vacuum_rect_united():
    """Phase C-heavy regression: united multi-loft + rect cross-section."""
    g = _load_golden("peec_inductance_rect_united_50kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        pytest.fail(f"Sample STEP missing: {step}")

    phys = g["physics"]
    result = _run_calc([
        "--coil-step", step,
        "--coil-solver", "peec",
        "--peec-n-peri", str(phys["n_peri"]),
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", str(phys["sigma_S_per_m"]),
    ])
    assert result["status"] == "ok"
    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                     guard["max_L_nH"], "L_coil_nH (hard band)",
                     hint=guard["comment"])
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")


# ============================================================================
# Mode 2: --coil-solver bem-a, vacuum (--vol omitted)
# ============================================================================

@pytest.mark.slow
def test_inductance_bem_a_vacuum_rect_united():
    """Cross-method: same rect_torus geometry via BEM-A.

    BEM-A converged value 153.5 nH (vs PEEC asymptote ~149 nH at high
    n_peri).  ~3.5% spread is real modeling difference (BEM-A surface
    RWG captures rect-corner crowding; PEEC perimeter filaments don't).
    """
    g = _load_golden("peec_inductance_rect_united_50kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        pytest.fail(f"Sample STEP missing: {step}")

    phys = g["physics"]
    # BEM-A reads pre-meshed .vol (4.33.0 redesign).  Pre-mesh the STEP
    # via the test helper -- this replaces the retired on-the-fly
    # _build_bema_coil_mesh() inside calc_inductance.py.
    from _bema_coil_vol_helper import coil_vol_for
    coil_vol = coil_vol_for(step, maxh=0.012)
    result = _run_calc([
        "--coil-vol", coil_vol,
        "--coil-solver", "bem-a",
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", str(phys["sigma_S_per_m"]),
    ])
    assert result["status"] == "ok"
    assert result["coil_solver_used"] == "bem-a"
    # BEM-A converged at maxh=0.012: 153.5 nH (cross-validation against
    # PEEC golden 145.30 with a 5.7% modeling-difference gap).  Hard
    # band 145-160 covers convergence drift across mesh densities.
    _assert_in_range(result["L_coil_nH"], 145.0, 160.0,
                     "L_coil_nH (BEM-A vacuum band)",
                     hint="BEM-A converged ~153.5 nH; ~3.5% above PEEC "
                          "due to rect-corner current crowding modelled "
                          "by surface RWG.")
    # Panel-level R lock for the unified impedance-EFIE (the ONLY
    # place the CLI-path BEM-A R is asserted since the PEC removal;
    # captured 0.6785 mOhm at 50 kHz / maxh=0.012 on LAB 2026-07-02,
    # +-20% band for per-machine re-mesh variation).
    _assert_in_range(result["R_coil_mOhm"], 0.54, 0.82,
                     "R_coil_mOhm (BEM-A impedance-EFIE band)",
                     hint="Captured 0.6785 mOhm (impedance-EFIE, Zs in "
                          "the saddle).  A big jump suggests the PEC "
                          "post-hoc over-concentration came back; a "
                          "collapse suggests Zs fell out of the "
                          "(1,1)=jw*mu0*SL+Zs*M block.")
    assert result["bem_a_residual"] < 1e-12


# ============================================================================
# Mode 3 + 4: --coil-solver {peec, bem-a}, weak-coupled (--vol present)
# ============================================================================

_WEAK_VOL = ("src", "radia", "panels", "samples", "ih_peec_bem_coarse.vol")
_WEAK_COIL = ("src", "radia", "panels", "samples",
              "ih_peec_inductance_coil.step")


def _have_weak_fixtures():
    vol = os.path.join(ROOT, *_WEAK_VOL)
    step = os.path.join(ROOT, *_WEAK_COIL)
    return os.path.isfile(vol) and os.path.isfile(step)


@pytest.mark.slow
@pytest.mark.skipif(not _have_weak_fixtures(),
                     reason="weak-coupled fixtures (.vol + coil .step) missing")
def test_inductance_peec_weak_gapped_torus():
    """Weak-coupled: PEEC coil + scalar BEM-SIBC workpiece.

    Cross-method check: BEM-A weak (next test) should agree on P_wp
    within ~5% modeling difference.
    """
    vol = os.path.join(ROOT, *_WEAK_VOL)
    step = os.path.join(ROOT, *_WEAK_COIL)
    result = _run_calc([
        "--coil-step", step,
        "--coil-solver", "peec",
        "--peec-n-peri", "16",
        "--frequency", "7000",
        "--current", "1.0",
        "--coil-sigma", "5.8e7",
        "--vol", vol,
        "--wp-label", "sibc",
        "--sigma", "5.8e7",
        "--mu-r", "1.0",
        "--half-thickness", "0.0125",
        "--impedance-model", "sibc",
        "--h1-order", "1",
        "--wp-bem-backend", "hacapk",
    ])
    assert result["status"] == "ok"
    assert result["coil_solver_used"] == "peec"
    assert result["coupling_mode"] == "weak"
    # Hard band: gapped torus R=30 a=3 in air gives PEEC L_air ~85 nH
    # (n_peri=16); workpiece gives P_wp ~6e-5 W at 7 kHz Cu.
    _assert_in_range(result["L_coil_nH"], 80.0, 95.0,
                     "L_coil_nH (peec weak)",
                     hint="PEEC vacuum L for gapped torus R=30 a=3")
    _assert_in_range(result["P_wp_W"], 4e-5, 9e-5,
                     "P_wp_W (peec weak)",
                     hint="2D SIBC ref ~6.6e-5 W at 7 kHz Cu R=25mm")
    # delta_L should be Lenz-correct (negative for non-magnetic Cu wp)
    assert result["delta_L_nH"] < 0


@pytest.mark.slow
@pytest.mark.skipif(not _have_weak_fixtures(),
                     reason="weak-coupled fixtures (.vol + coil .step) missing")
def test_inductance_bem_a_weak_gapped_torus():
    """Weak-coupled: BEM-A coil + scalar BEM-SIBC workpiece.

    BEM-A version of the test above.  L_coil agrees with PEEC weak to
    ~3% (mesh limit on circular cross-section); P_wp agrees to ~5%.
    """
    g = _load_golden("coil_bem_a_workpiece_7kHz_Cu.json")
    coil_step = _resolve_sample(g["sample"]["coil_step"])
    vol = _resolve_sample(g["sample"]["vol"])
    if not (os.path.isfile(coil_step) and os.path.isfile(vol)):
        pytest.skip(f"missing fixture: {coil_step} or {vol}")

    phys = g["physics"]
    # BEM-A 4.33.0 redesign: pre-mesh coil STEP -> .vol (with
    # source/sink labels) before calling calc_inductance.py.
    from _bema_coil_vol_helper import coil_vol_for
    coil_vol = coil_vol_for(coil_step, maxh=phys["coil_maxh_m"])
    result = _run_calc([
        "--coil-vol", coil_vol,
        "--coil-solver", "bem-a",
        "--coil-rwg-quad-degree", str(phys["coil_rwg_quad_degree"]),
        "--coil-rwg-singular-nq", str(phys["coil_rwg_singular_nq"]),
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", "5.8e7",
        "--vol", vol,
        "--wp-label", "sibc",
        "--sigma", "5.8e7",
        "--mu-r", "1.0",
        "--half-thickness", "0.0125",
        "--impedance-model", "sibc",
        "--h1-order", str(phys["h1_order"]),
        "--wp-bem-backend", phys["wp_bem_backend"],
        "--wp-aca-eps", str(phys["wp_aca_eps"]),
    ])
    assert result["status"] == "ok"
    assert result["coil_solver_used"] == "bem-a"
    assert result["coupling_mode"] == "weak"

    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"],
                     guard["min_L_coil_nH"], guard["max_L_coil_nH"],
                     "L_coil_nH (bem-a weak hard band)",
                     hint=guard["comment"])
    _assert_in_range(result["P_wp_W"],
                     guard["min_P_wp_W"], guard["max_P_wp_W"],
                     "P_wp_W (bem-a weak hard band)")
    _assert_in_range(result["H_t_rms_A_per_m"],
                     guard["min_H_t_rms"], guard["max_H_t_rms"],
                     "H_t_rms (bem-a weak hard band)")

    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                   exp["tolerance_L_coil_pct"], "L_coil_nH (golden)")
    _assert_close(result["P_wp_W"], exp["P_wp_W"],
                   exp["tolerance_P_pct"], "P_wp_W (golden)")
    _assert_close(result["H_t_rms_A_per_m"], exp["H_t_rms_A_per_m"],
                   exp["tolerance_Ht_rms_pct"], "H_t_rms (golden)")
