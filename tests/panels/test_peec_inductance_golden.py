"""Golden regression tests for calc_peec_inductance.py.

Locks in the two end-to-end paths of the PEEC-inductance (coil only,
STEP) panel mode:

  1. Gapped torus (single-loop) STEP -> TORUS analytical sweep ->
     n_peri=16 perimeter PEEC -> L ~ 85 nH
  2. Multi-turn loft STEP        -> cross-section-centroid chain ->
     n_peri=16 perimeter PEEC -> L ~ 431 nH

Both tests run the real calc_peec_inductance.py subprocess (no direct
Python API calls) so the JSON protocol + argparse contract are part
of what's locked.

Run::

    pytest tests/panels/test_peec_inductance_golden.py -v

These tests are marked ``slow`` (2 - 20 s each); the deploy skill's
Panel Mode Matrix calls them on Stage 2 (100号機) and Stage 3 (mdx).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CALC = os.path.join(ROOT, "src", "radia", "panels",
                     "calc_peec_inductance.py")


def _load_golden(fname: str) -> dict:
    with open(os.path.join(HERE, "golden", fname), "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_sample(path: str) -> str:
    """Resolve ``path`` as either an absolute file or repo-relative."""
    if os.path.isabs(path):
        return path
    return os.path.join(ROOT, path)


def _run_peec(step_or_jou: str, n_peri: int, freq_hz: float,
              current_A: float, sigma: float, timeout_s: int = 600) -> dict:
    """Run calc_peec_inductance.py as a subprocess, return the JSON result."""
    cmd = [sys.executable, CALC,
           "--peec-step", step_or_jou,
           "--peec-n-peri", str(n_peri),
           "--frequency", str(freq_hz),
           "--current", str(current_A),
           "--coil-sigma", str(sigma)]
    proc = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout_s)
    assert proc.returncode == 0, (
        f"calc_peec_inductance returned {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout[-2000:]}\n"
        f"STDERR:\n{proc.stderr[-2000:]}")
    # Last JSON line of stdout (progress messages go to stderr)
    json_lines = [ln for ln in proc.stdout.splitlines()
                  if ln.strip().startswith("{")]
    assert json_lines, f"no JSON in stdout:\n{proc.stdout}"
    return json.loads(json_lines[-1])


def _assert_in_range(value: float, lo: float, hi: float, label: str,
                      hint: str = "") -> None:
    assert lo <= value <= hi, (
        f"{label}: {value:.3f} outside regression guard [{lo}, {hi}].  "
        f"{hint}")


def _assert_close(value: float, expected: float, tol_pct: float,
                   label: str) -> None:
    err_pct = abs(value - expected) / max(abs(expected), 1e-30) * 100.0
    assert err_pct <= tol_pct, (
        f"{label}: got {value:.4f}, expected {expected:.4f} "
        f"(err {err_pct:.2f}%, tol {tol_pct}%)")


# ---------------------------------------------------------------------------
# Test 1: gapped torus single-loop (TORUS analytical sweep)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_peec_inductance_gapped_torus():
    """355deg gapped torus via STEP-only: TORUS analytical path."""
    g = _load_golden("peec_inductance_torus_50kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        pytest.skip(
            f"Sample STEP missing: {step}.  Regenerate via Cubit "
            f"playback of {g['sample']['generator_jou']}.")

    result = _run_peec(step,
                       n_peri=g["physics"]["n_peri"],
                       freq_hz=g["physics"]["frequency_Hz"],
                       current_A=g["physics"]["current_A"],
                       sigma=g["physics"]["sigma_S_per_m"])

    assert result["status"] == "ok", f"solver status != ok: {result}"
    assert result["input_kind"] == g["expected"]["input_kind"], (
        f"wrong input_kind: got {result['input_kind']}, "
        f"expected {g['expected']['input_kind']} "
        f"(sibling .jou auto-preference may have misfired)")

    # Regression hard-band first (catches order-of-magnitude breaks)
    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                      guard["max_L_nH"], "L_coil_nH (hard band)",
                      hint=guard["comment"])

    # Golden tolerance (catches finer drift)
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")
    _assert_close(result["R_coil_mOhm"], exp["R_coil_mOhm"],
                  exp["tolerance_R_pct"], "R_coil_mOhm (golden)")


# ---------------------------------------------------------------------------
# Test 2: 3-turn pancake loft (cross-section centroid chain)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_peec_inductance_3turn_loft():
    """3-turn pancake loft STEP-only: cross-section-centroid chain path."""
    g = _load_golden("peec_inductance_3turn_150kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        if g["sample"].get("skip_if_missing", False):
            pytest.skip(
                f"Sample STEP not available on this machine: {step}. "
                f"This test is designed to run on LAB / 100号機 where "
                f"Kubota's W: share is mounted.")
        else:
            pytest.fail(f"Sample STEP missing: {step}")

    # Copy STEP to an ASCII path with NO sibling .jou so sibling-
    # preference does NOT kick in and we exercise the STEP-only
    # centerline extraction.
    import shutil
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        isolated_step = os.path.join(td, "coil.stp")
        shutil.copy2(step, isolated_step)
        result = _run_peec(isolated_step,
                           n_peri=g["physics"]["n_peri"],
                           freq_hz=g["physics"]["frequency_Hz"],
                           current_A=g["physics"]["current_A"],
                           sigma=g["physics"]["sigma_S_per_m"])

    assert result["status"] == "ok", f"solver status != ok: {result}"
    assert result["input_kind"] == "STEP", (
        f"expected STEP path; got {result['input_kind']} "
        f"(sibling .jou might have been auto-copied)")

    # Hard band first
    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                      guard["max_L_nH"], "L_coil_nH (hard band)",
                      hint=guard["comment"])

    # Golden tolerance
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")
    _assert_close(result["R_coil_mOhm"], exp["R_coil_mOhm"],
                  exp["tolerance_R_pct"], "R_coil_mOhm (golden)")


# ---------------------------------------------------------------------------
# Test 3: united multi-loft RECT cross-section torus -- Phase C-heavy
# (Path 2c, _filaments_from_section_planes).  Locks in v4.24.0+
# behaviour: previously NaN, now 138 nH.
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_peec_inductance_rect_united():
    """Phase C-heavy regression: united multi-loft + rect cross-section.

    The fixture is a synthetic 355deg torus with a 8x6mm rect cross-
    section, built as 12 ruled lofts then UNITED into one solid.  Pre-
    v4.24.0 this returned NaN (Tier 3 longest-edge picked a lateral
    rect edge as spine).  v4.24.0 added Path 2c (section-planes)
    which recovers a proper cross-section at each of n_stations sample
    points along a rotation-axis-derived spine.
    """
    g = _load_golden("peec_inductance_rect_united_50kHz_Cu.json")
    step = _resolve_sample(g["sample"]["step"])
    if not os.path.isfile(step):
        pytest.fail(f"Sample STEP missing: {step}")

    result = _run_peec(step,
                       n_peri=g["physics"]["n_peri"],
                       freq_hz=g["physics"]["frequency_Hz"],
                       current_A=g["physics"]["current_A"],
                       sigma=g["physics"]["sigma_S_per_m"])

    assert result["status"] == "ok", f"solver status != ok: {result}"
    assert result["input_kind"] == g["expected"]["input_kind"]

    # Hard band first
    guard = g["regression_guard"]
    _assert_in_range(result["L_coil_nH"], guard["min_L_nH"],
                      guard["max_L_nH"], "L_coil_nH (hard band)",
                      hint=guard["comment"])

    # Golden tolerance
    exp = g["expected"]
    _assert_close(result["L_coil_nH"], exp["L_coil_nH"],
                  exp["tolerance_L_pct"], "L_coil_nH (golden)")
