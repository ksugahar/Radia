"""Golden regression test for calc_peec_bem.py.

Purpose: lock in the known-good P output of the PEEC+BEM 1-way panel
path on a fixed coarse sample, so future edits to the BND extractor /
winding-orientation / BEM assembly that regress the result fail here,
not during user-facing validation 3 months later.

Background: during the 2026-04-19 session the panel integration of
calc_peec_bem was developed while the demo script
peec_bundle_bem_sibc_demo.py was known to agree with the 2D axisym
reference within 5-11%.  Porting the demo's OCC-based wp extraction
to a .vol-based sibc-sideset extraction introduced a winding-orientation
bug (FD-per-FD inconsistency from the webcut-split air volumes) that
silently gave P 20x too low or 4x too high depending on mesh.  The
user pointed out the core problem: knowledge proven in the demo was
not captured in a way that made it reproducible in the panel path.
This test closes that gap.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(HERE, "golden", "peec_bem_coarse_7kHz_Cu.json")
CALC = os.path.join(ROOT, "src", "radia", "panels", "calc_peec_bem.py")


def _load_golden():
    with open(GOLDEN, "r", encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.slow
def test_peec_bem_matches_2d_ref_coarse():
    g = _load_golden()
    vol = os.path.join(ROOT, g["sample"]["vol"])
    step = os.path.join(ROOT, g["sample"]["step"])
    if not (os.path.isfile(vol) and os.path.isfile(step)):
        pytest.skip(
            f"Sample files missing -- regenerate via Cubit batch using "
            f"{g['sample']['generator_jou']}:\n"
            f"  coreform_cubit -batch -nojournal -nographics "
            f"-input {g['sample']['generator_jou']}\n"
            f"Missing: {vol if not os.path.isfile(vol) else step}")

    phys = g["physics"]
    cmd = [
        sys.executable, CALC,
        "--peec-step", step,
        "--peec-nwinc", str(phys["peec_nwinc"]),
        "--peec-nhinc", str(phys["peec_nhinc"]),
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", "5.8e7",
        "--vol", vol,
        "--wp-label", "sibc",
        "--sigma", "5.8e7",
        "--half-thickness", str(phys["half_thickness_m"]),
        "--mu-r", "1.0",
        "--impedance-model", "sibc",
        "--h1-order", str(phys["h1_order"]),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, (
        f"calc_peec_bem failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")

    # calc_main emits JSON to stdout.
    result = json.loads(proc.stdout)

    exp = g["expected"]

    # PEEC coil impedance (geometry-only; tight tolerance).
    assert abs(result["L_coil_nH"] - exp["L_coil_peec_nH"]) \
        / exp["L_coil_peec_nH"] * 100 < exp["tolerance_L_coil_pct"], (
            f"L_coil drift: got {result['L_coil_nH']:.3f}, "
            f"expected {exp['L_coil_peec_nH']:.3f} "
            f"(tol {exp['tolerance_L_coil_pct']}%)")

    # Workpiece dissipation vs 2D SIBC reference (mesh-sensitive; 15% band).
    p_ref = g["reference_2d_axisym_sibc"]["P_W"]
    p_got = result["P_wp_W"]
    err_pct = abs(p_got - p_ref) / p_ref * 100
    assert err_pct < exp["tolerance_P_pct"], (
        f"P_wp drift: got {p_got:.4e} W, 2D SIBC ref {p_ref:.4e} W, "
        f"err {err_pct:.2f}% (tol {exp['tolerance_P_pct']}%).  "
        f"If this fails, the normal-winding logic in "
        f"calc_peec_bem._extract_bnd_only likely regressed.")

    # Surface heat-flux distribution (Phase A optimization filter).
    # q_surf is the same physical quantity as P_wp/A_wp at the mean,
    # but the max/p95 values lock the SPATIAL distribution: a regress
    # that makes heat land in the wrong place (e.g. wrong tangential
    # projection) shows up here even when the integrated P matches.
    qmax_got = result["q_surf_max"]
    qmean_got = result["q_surf_mean"]
    qmax_ref = exp["q_surf_max_Wm2"]
    qmean_ref = exp["q_surf_mean_Wm2"]
    err_qmax = abs(qmax_got - qmax_ref) / qmax_ref * 100
    err_qmean = abs(qmean_got - qmean_ref) / qmean_ref * 100
    assert err_qmax < exp["tolerance_qsurf_max_pct"], (
        f"q_surf_max drift: got {qmax_got:.3e}, ref {qmax_ref:.3e}, "
        f"err {err_qmax:.2f}% (tol {exp['tolerance_qsurf_max_pct']}%).")
    assert err_qmean < exp["tolerance_qsurf_mean_pct"], (
        f"q_surf_mean drift: got {qmean_got:.3e}, ref {qmean_ref:.3e}, "
        f"err {err_qmean:.2f}% (tol {exp['tolerance_qsurf_mean_pct']}%).")
    # Cross-check: surface integral of q_surf must equal P_wp.
    p_check = result["P_total_check"]
    err_p_int = abs(p_check - p_got) / p_got * 100
    assert err_p_int < 0.1, (
        f"P_total_check / P_wp inconsistency: int q_surf dA = {p_check:.4e} W, "
        f"P_wp = {p_got:.4e} W, err {err_p_int:.4f}% (tol 0.1%). "
        f"If this fails, the q_surf CF construction has a different "
        f"normalization than the bem.solve() P_density.")
