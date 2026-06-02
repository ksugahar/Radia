"""Standalone golden for calc_heat_axisym.py (Cubit-free, uniform-q).

Locks the 2D-axisymmetric thermal solver and its 2*pi*r surface
integral WITHOUT the EM chain.  test_heat_chain_golden.py's axisym
leg exercises the cross-mesh spatial projection, but it SKIPS when
the coil STEP sample (ih_peec_bem_coarse_coil.step) is absent -- so
the axisym solver had no CI-runnable golden until this one.

The headline assertion is the analytic cross-check
``q_surf_int_W == q_uniform * surface_area``: a uniform surface flux
integrated over the 2*pi*r-weighted lateral surface must equal the
flux times the true lateral area, which catches any regression in
the axisymmetric weighting of the boundary integral.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(HERE, "golden", "heat_axisym_uniform.json")
FIXTURE_DIR = os.path.join(HERE, "fixtures")
WP_VOL_AXI = os.path.join(FIXTURE_DIR,
                          "heat_workpiece_cylinder_R25_H25_axisym.vol")
GEN_AXI = os.path.join(FIXTURE_DIR, "generate_heat_cylinder_axisym.py")
CALC = os.path.join(ROOT, "src", "radia", "panels", "calc_heat_axisym.py")


def _load_golden():
    with open(GOLDEN, "r", encoding="utf-8") as f:
        return json.load(f)


def _ensure_fixture():
    if os.path.isfile(WP_VOL_AXI):
        return
    proc = subprocess.run([sys.executable, GEN_AXI],
                          capture_output=True, text=True, timeout=120)
    if proc.returncode != 0 or not os.path.isfile(WP_VOL_AXI):
        pytest.skip(f"could not generate axisym fixture:\n"
                    f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def test_heat_axisym_uniform_golden():
    g = _load_golden()
    _ensure_fixture()
    c = g["cli"]
    exp = g["expected"]
    cmd = [
        sys.executable, CALC,
        "--wp-vol", WP_VOL_AXI,
        "--surface-label", c["surface_label"],
        "--material", c["material"],
        "--q-uniform", str(c["q_uniform_Wm2"]),
        "--h-conv", str(c["h_conv_Wm2K"]),
        "--t-ext", str(c["t_ext_C"]),
        "--t-initial", str(c["t_initial_C"]),
        "--dt", str(c["dt_s"]),
        "--t-end", str(c["t_end_s"]),
        "--probe-point", f"{c['probe_point_rz'][0]},{c['probe_point_rz'][1]}",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    assert proc.returncode == 0, (
        f"calc_heat_axisym failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")
    r = json.loads(proc.stdout.strip().splitlines()[-1])

    assert r["mesh_type"] == "axisymmetric"

    # Mesh fingerprint -- catches a silent fixture-regen drift.
    assert r["ne"] == exp["ne"], f"ne {r['ne']} != {exp['ne']}"
    assert r["ndof"] == exp["ndof"], f"ndof {r['ndof']} != {exp['ndof']}"

    # Headline: analytic axisym surface-integral cross-check.
    q_int = r["q_surf_int_W"]
    analytic = c["q_uniform_Wm2"] * r["surface_area_m2"]
    assert abs(q_int - analytic) / analytic < 1e-6, (
        f"q_surf_int {q_int} != q_uniform*area {analytic} -- the axisym "
        f"2*pi*r surface integral regressed.")

    tol = exp["tolerance_rel_pct"] / 100.0
    for key in ("q_surf_int_W", "surface_area_m2"):
        assert abs(r[key] - exp[key]) / exp[key] < tol, (
            f"{key} drift: {r[key]} vs {exp[key]}")

    tK = exp["tolerance_temp_K"]
    assert abs(r["T_max_C"] - exp["T_max_C"]) < tK, (
        f"T_max drift: {r['T_max_C']} vs {exp['T_max_C']}")
    assert abs(r["T_min_C"] - exp["T_min_C"]) < tK, (
        f"T_min drift: {r['T_min_C']} vs {exp['T_min_C']}")

    # Monotone heating (uniform source, backward Euler, no cooling
    # below the initial temperature).
    hist = r["T_probe_history_C"]
    for i in range(1, len(hist)):
        assert hist[i] >= hist[i - 1] - 1e-9, (
            f"non-monotone T_probe at step {i}: "
            f"{hist[i-1]} -> {hist[i]}")
