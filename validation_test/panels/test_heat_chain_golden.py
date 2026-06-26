"""Golden regression test for the IH thermal pipeline (Phase B).

Chain under test::

    calc_fem_kelvin.py  --vol <em_vol> --peec-step <step> ... \
                        --msh-output <em_msh>           # writes <em_msh stem>_qsurf.sol
                                                        # and <em_msh stem>_fem.vol
        |
        v
    calc_heat.py        --wp-vol <wp_vol>                # workpiece thermal mesh
                        --qsurf-sol <qsurf.sol>          # cross-mesh source
                        --em-vol <em fem.vol>            # geometry context
                        ...                              # thermal params

Why this test exists
====================

The cross-mesh projection in ``calc_heat._build_qsurf_cf`` walks every
workpiece-mesh surface vertex and queries the EM mesh's GridFunction
at that point.  A regression in any of:

  - calc_fem_kelvin's qsurf.sol layout (FES order / dim / mesh
    correspondence),
  - the .sol companion path convention
    (``<msh stem>_qsurf.sol`` and ``<msh stem>_fem.vol``),
  - calc_heat's surface-vertex enumeration (when filtering by
    ``--surface-label``),

would silently leave most of the workpiece surface at zero flux and
the thermal solve would look "right but cold".  The integral cross-
check ``int q_surf dA`` vs ``P_total`` from the EM run is the tight
detector.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(HERE, "golden", "heat_chain_coarse_7kHz.json")
FIXTURE_DIR = os.path.join(HERE, "fixtures")
WP_VOL = os.path.join(FIXTURE_DIR, "heat_workpiece_cylinder_R25_H25.vol")
GEN_SCRIPT = os.path.join(FIXTURE_DIR, "generate_heat_cylinder.py")
WP_VOL_AXI = os.path.join(FIXTURE_DIR,
                          "heat_workpiece_cylinder_R25_H25_axisym.vol")
GEN_SCRIPT_AXI = os.path.join(FIXTURE_DIR,
                              "generate_heat_cylinder_axisym.py")
CALC_FEM = os.path.join(ROOT, "src", "radia", "panels", "calc_fem_kelvin.py")
CALC_HEAT = os.path.join(ROOT, "src", "radia", "panels", "calc_heat.py")
CALC_HEAT_AXI = os.path.join(ROOT, "src", "radia", "panels",
                              "calc_heat_axisym.py")


def _load_golden():
    with open(GOLDEN, "r", encoding="utf-8") as f:
        return json.load(f)


def _ensure_wp_fixture(tmp_path):
    """Generate the workpiece-volume cylinder .vol fixture on demand
    via the OCC-only helper.  No-op when the fixture already exists.
    Returns the path."""
    if os.path.isfile(WP_VOL):
        return WP_VOL
    proc = subprocess.run(
        [sys.executable, GEN_SCRIPT],
        capture_output=True, text=True, timeout=120)
    if proc.returncode != 0 or not os.path.isfile(WP_VOL):
        pytest.skip(
            f"Could not generate workpiece fixture {WP_VOL}.\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return WP_VOL


def _ensure_wp_axisym_fixture(tmp_path):
    """Generate the 2D axisym workpiece .vol fixture on demand."""
    if os.path.isfile(WP_VOL_AXI):
        return WP_VOL_AXI
    proc = subprocess.run(
        [sys.executable, GEN_SCRIPT_AXI],
        capture_output=True, text=True, timeout=120)
    if proc.returncode != 0 or not os.path.isfile(WP_VOL_AXI):
        pytest.skip(
            f"Could not generate axisym fixture {WP_VOL_AXI}.\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return WP_VOL_AXI


@pytest.mark.slow
def test_fem_kelvin_to_heat_chain(tmp_path):
    g = _load_golden()
    em_vol = os.path.join(ROOT, g["sample"]["em_vol"])
    step = os.path.join(ROOT, g["sample"]["em_step"])
    if not (os.path.isfile(em_vol) and os.path.isfile(step)):
        pytest.skip(
            f"IH sample missing -- regenerate via Cubit batch using "
            f"src/radia/panels/samples/ih_peec_bem_coarse.jou:\n"
            f"Missing: {em_vol if not os.path.isfile(em_vol) else step}")

    wp_vol = _ensure_wp_fixture(tmp_path)

    phys = g["physics"]
    em_msh = str(tmp_path / "em.msh")
    em_qsurf = str(tmp_path / "em_qsurf.sol")
    em_fem_vol = str(tmp_path / "em_fem.vol")

    # Step 1: calc_fem_kelvin -> qsurf.sol + em.vol companion.
    em_cmd = [
        sys.executable, CALC_FEM,
        "--vol", em_vol,
        "--fes-order", str(phys["fes_order_em"]),
        "--frequency", str(phys["frequency_Hz"]),
        "--material", "custom",
        "--sigma", str(phys["wp_sigma_Sm"]),
        "--mu-r", str(phys["wp_mu_r"]),
        "--impedance", "sibc",
        "--formulation", "total",
        "--current", str(phys["current_A"]),
        "--half-thickness", str(phys["half_thickness_m"]),
        "--solver", "pardiso",
        "--peec-step", step,
        "--peec-sigma", str(phys["wp_sigma_Sm"]),
        "--peec-nwinc", "3",
        "--peec-nhinc", "3",
        "--msh-output", em_msh,
    ]
    proc = subprocess.run(em_cmd, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, (
        f"calc_fem_kelvin failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")
    em_result = json.loads(proc.stdout)

    # The .sol + companion .vol naming is part of the contract; the
    # rename helper in calc_fem_kelvin uses ``<msh-stem>_qsurf.sol``
    # and ``<msh-stem>_fem.vol`` so the chain can find them by
    # convention.  Hard-fail if either is missing.
    msh_stem = os.path.splitext(em_msh)[0]
    assert os.path.isfile(msh_stem + "_qsurf.sol"), (
        f"calc_fem_kelvin did not write _qsurf.sol next to "
        f"{em_msh}.  Phase A's .sol output regressed?")
    assert os.path.isfile(msh_stem + "_fem.vol"), (
        f"calc_fem_kelvin did not write _fem.vol next to "
        f"{em_msh}.  vol2msh save_vol_sol_pair regressed?")
    em_qsurf = msh_stem + "_qsurf.sol"
    em_fem_vol = msh_stem + "_fem.vol"

    # EM-side cross-checks (loose -- this test exists to lock the
    # CHAIN, not to re-litigate calc_fem_kelvin's own physics).
    exp_em = g["expected_em"]
    tol_em = exp_em["tolerance_em_pct"]
    assert abs(em_result["P_total"] - exp_em["P_total_W"]) \
        / exp_em["P_total_W"] * 100 < tol_em, (
        f"P_total drift: {em_result['P_total']:.4e} W, "
        f"ref {exp_em['P_total_W']:.4e} W")
    assert em_result["q_surf_max"] is not None
    assert em_result["q_surf_mean"] is not None

    # Step 2: calc_heat consumes qsurf.sol and the workpiece mesh.
    heat_msh = str(tmp_path / "heat_T.msh")
    heat_csv = str(tmp_path / "heat_probe.csv")
    heat_cmd = [
        sys.executable, CALC_HEAT,
        "--wp-vol", wp_vol,
        "--surface-label", "outer",
        "--material", phys["thermal_material"],
        "--qsurf-sol", em_qsurf,
        "--em-vol", em_fem_vol,
        "--qsurf-order", str(phys["fes_order_em"]),
        "--h-conv", str(phys["h_conv_Wm2K"]),
        "--t-ext", str(phys["t_ext_C"]),
        "--t-initial", str(phys["t_initial_C"]),
        "--dt", str(phys["dt_s"]),
        "--t-end", str(phys["t_end_s"]),
        "--probe-point",
        f"{phys['probe_point_m'][0]},{phys['probe_point_m'][1]},"
        f"{phys['probe_point_m'][2]}",
        "--csv-output", heat_csv,
        "--msh-output", heat_msh,
    ]
    proc = subprocess.run(heat_cmd, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, (
        f"calc_heat failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")
    heat_result = json.loads(proc.stdout)

    exp_h = g["expected_heat"]

    # Cross-mesh transfer integrity: q_surf integrated on the WP
    # thermal surface should match the EM-side P_total within mesh
    # interpolation noise.  This is the headline check: a broken
    # surface-label filter or a broken vertex enumeration leaves
    # most vertices at zero and the integral collapses.
    qsurf_int = heat_result["q_surf_int_W"]
    err_pct = abs(qsurf_int - exp_h["qsurf_int_W"]) \
        / exp_h["qsurf_int_W"] * 100
    assert err_pct < exp_h["tolerance_qsurf_int_pct"], (
        f"qsurf_int drift: got {qsurf_int:.4e} W, "
        f"ref {exp_h['qsurf_int_W']:.4e} W, "
        f"err {err_pct:.2f}% (tol {exp_h['tolerance_qsurf_int_pct']}%).  "
        f"Most likely cause: cross-mesh transfer regressed -- "
        f"check the surface-vertex enumeration in "
        f"calc_heat._build_qsurf_cf.")

    # Cross-mesh sanity vs the EM-side P_total.  Looser bracket
    # because the two surfaces are meshed independently.
    p_em = em_result["P_total"]
    rel = abs(qsurf_int - p_em) / max(p_em, 1e-30) * 100
    assert rel < 5.0, (
        f"q_surf integrated on the THERMAL surface ({qsurf_int:.4e} W) "
        f"differs from the EM-side P_total ({p_em:.4e} W) by "
        f"{rel:.2f}% -- expected <5% (mesh interpolation noise).")

    # Final-state temperature checks.  Absolute K bracket because
    # the test current is 1 A and the resulting temperature rise
    # is in the ~1e-5 K range (linear in q_surf, so kA-current
    # production runs land in hundreds of K).
    T_max = heat_result["T_max_C"]
    T_probe = heat_result["T_probe_history_C"][-1]
    tol_K = exp_h["tolerance_temp_K"]
    assert abs(T_max - exp_h["T_max_C"]) < tol_K, (
        f"T_max drift: got {T_max:.10f} C, "
        f"ref {exp_h['T_max_C']:.10f} C, "
        f"abs err {abs(T_max - exp_h['T_max_C']):.2e} K (tol {tol_K} K).")
    assert abs(T_probe - exp_h["T_probe_final_C"]) < tol_K, (
        f"T_probe drift: got {T_probe:.10f} C, "
        f"ref {exp_h['T_probe_final_C']:.10f} C, "
        f"abs err {abs(T_probe - exp_h['T_probe_final_C']):.2e} K "
        f"(tol {tol_K} K).")

    # Sanity: monotone heating from initial temperature.
    history = heat_result["T_probe_history_C"]
    for i in range(1, len(history)):
        assert history[i] >= history[i - 1] - 1e-12, (
            f"T_probe history is non-monotone at step {i}: "
            f"{history[i-1]:.10f} -> {history[i]:.10f}")

    # ----------------------------------------------------------------
    # Axisymmetric leg: reuse the same EM-side qsurf.sol on a 2D
    # axisym mesh and verify the v1.5 path agrees with the 3D leg
    # within mesh-interpolation noise.  This is the regression catch
    # for calc_heat_axisym.py specifically (phi-averaging, axisym
    # weight 2*pi*r in the bilinear forms, 2D mesh handling).
    # ----------------------------------------------------------------
    exp_a = g["expected_heat_axisym"]
    wp_vol_axi = _ensure_wp_axisym_fixture(tmp_path)
    heat_axi_csv = str(tmp_path / "heat_probe_axi.csv")
    heat_axi_cmd = [
        sys.executable, CALC_HEAT_AXI,
        "--wp-vol", wp_vol_axi,
        "--surface-label", "outer",
        "--material", phys["thermal_material"],
        "--qsurf-sol", em_qsurf,
        "--em-vol", em_fem_vol,
        "--qsurf-order", str(phys["fes_order_em"]),
        "--n-phi-samples", str(exp_a["n_phi_samples"]),
        "--h-conv", str(phys["h_conv_Wm2K"]),
        "--t-ext", str(phys["t_ext_C"]),
        "--t-initial", str(phys["t_initial_C"]),
        "--dt", str(phys["dt_s"]),
        "--t-end", str(phys["t_end_s"]),
        "--probe-point",
        f"{exp_a['probe_point_rz'][0]},{exp_a['probe_point_rz'][1]}",
        "--csv-output", heat_axi_csv,
    ]
    proc = subprocess.run(heat_axi_cmd, capture_output=True,
                           text=True, timeout=300)
    assert proc.returncode == 0, (
        f"calc_heat_axisym failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")
    heat_axi_result = json.loads(proc.stdout)

    # Schema marker.
    assert heat_axi_result.get("mesh_type") == "axisymmetric"

    # qsurf integral: lock the captured value AND verify it agrees
    # with the EM-side P_total.  Two checks because the 3D and
    # axisym surfaces are meshed independently and the axisym side
    # carries an extra phi-average step.
    qsurf_int_axi = heat_axi_result["q_surf_int_W"]
    err_pct = abs(qsurf_int_axi - exp_a["qsurf_int_W"]) \
        / exp_a["qsurf_int_W"] * 100
    assert err_pct < exp_a["tolerance_qsurf_int_pct"], (
        f"axisym qsurf_int drift: got {qsurf_int_axi:.4e} W, "
        f"ref {exp_a['qsurf_int_W']:.4e} W, "
        f"err {err_pct:.2f}% (tol {exp_a['tolerance_qsurf_int_pct']}%).")
    rel_em = abs(qsurf_int_axi - p_em) / max(p_em, 1e-30) * 100
    assert rel_em < 5.0, (
        f"axisym qsurf vs EM P_total: got {qsurf_int_axi:.4e} W, "
        f"P_total {p_em:.4e} W, err {rel_em:.2f}% (>5%).")

    # Final-state temperature.
    T_max_axi = heat_axi_result["T_max_C"]
    T_probe_axi = heat_axi_result["T_probe_history_C"][-1]
    assert abs(T_max_axi - exp_a["T_max_C"]) < exp_a["tolerance_temp_K"]
    assert abs(T_probe_axi - exp_a["T_probe_final_C"]) < exp_a["tolerance_temp_K"]

    # Cross-leg sanity: the 3D and axisym legs solve the SAME
    # physical problem (cylindrical workpiece + cylindrical EM
    # mesh).  Their final-state probe temperatures must agree
    # within mesh-interpolation noise.
    assert abs(T_probe_axi - T_probe) < 5.0e-6, (
        f"3D vs axisym disagree on final probe temperature: "
        f"3D={T_probe:.10f} C, axisym={T_probe_axi:.10f} C, "
        f"diff {abs(T_probe_axi - T_probe):.2e} K (tol 5e-6 K).")
