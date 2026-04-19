"""Golden regression test for calc_fem_coilmesh.py (A-V formulation).

Locks in the A-V coil (volumetric sigma + phi Dirichlet source/sink)
+ wp SIBC + Kelvin formulation against the 2D axisym reference.

The non-obvious parts captured here (to prevent regressions):
  1. A-V compound FES with phi defined only on coil material, Dirichlet
     1 on 'source' / 0 on 'sink'.  Solve via Dirichlet lift.
  2. Current extraction is a VOLUME integral
     I_out = int_coil J . grad(psi_n) dV, with psi_n a scalar H1 test
     function with psi_n=1 on source, 0 on sink.  This is
     FEM-consistent (per MCP ih_knowledge), unlike cut-plane surface
     flux.
  3. P_wp from wp SIBC Robin surface formula (matches 2D ref to <2%).
  4. P_coil is volumetric |J|^2/sigma, mesh-sensitive — the Dowell
     thin-wire analytical formula was demoted from panel to
     examples/ 2026-04-19.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(HERE, "golden",
                       "fem_coilmesh_gapped_fine_7kHz_Cu.json")
CALC = os.path.join(ROOT, "src", "radia", "panels",
                     "calc_fem_coilmesh.py")


def _load_golden():
    with open(GOLDEN, "r", encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.slow
def test_fem_coilmesh_matches_2d_ref_gapped():
    g = _load_golden()
    vol = os.path.join(ROOT, g["sample"]["vol"])
    if not os.path.isfile(vol):
        pytest.skip(
            f"Sample .vol missing.  Regenerate via Cubit headless:\n"
            f"  coreform_cubit -batch -nojournal -nographics -input "
            f"{g['sample']['generator_jou']}\n"
            f"Missing: {vol}")

    phys = g["physics"]
    cmd = [
        sys.executable, CALC,
        "--vol", vol,
        "--frequency", str(phys["frequency_Hz"]),
        "--current", str(phys["current_A"]),
        "--coil-sigma", "5.8e7",
        "--sigma", "5.8e7",
        "--mu-r", "1.0",
        "--half-thickness", str(phys["half_thickness_m"]),
        "--fes-order", str(phys["fes_order"]),
        "--solver", "pardiso",
        "--sibc-bnd", "sibc",
        "--source-bnd", "source",
        "--sink-bnd", "sink",
        "--coil-mat", "coil",
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, (
        f"calc_fem_coilmesh failed:\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}")

    result = json.loads(proc.stdout)
    exp = g["expected"]

    # PRIMARY metric for IH engineering: P_wp (workpiece heating)
    # against 2D axisym SIBC reference.  Both PEEC+BEM and FEM A-V
    # agree with 2D to <2% on this — that's the value of having two
    # independent 3D methods.
    P_wp_got = result["P_wp_W"]
    P_wp_ref = g["reference_2d_axisym_sibc"]["P_wp_W"]
    P_wp_err = abs(P_wp_got - P_wp_ref) / P_wp_ref * 100
    assert P_wp_err < exp["tolerance_Pwp_pct"], (
        f"P_wp drift: got {P_wp_got:.4e}, 2D SIBC ref {P_wp_ref:.4e}, "
        f"err {P_wp_err:.2f}% (tol {exp['tolerance_Pwp_pct']}%).  "
        f"If this fails, suspect the wp SIBC Robin or Kelvin.")

    # L against 2D axisym full-eddy reference.  Wider tolerance since
    # FEM L depends on coil mesh resolution and is inherently slightly
    # low on under-resolved coil.
    L_got = result["L_total_nH"]
    L_ref = g["reference_2d_axisym_full_eddy"]["L_nH"]
    L_err = abs(L_got - L_ref) / L_ref * 100
    assert L_err < 10.0, (
        f"L drift: got {L_got:.2f} nH, 2D ref {L_ref:.2f} nH, "
        f"err {L_err:.2f}% (tol 10%, mesh-sensitive).")

    # P_coil is mesh-sensitive; assert against captured value, not 2D.
    P_coil_got = result["P_coil_W"]
    P_coil_exp = exp["P_coil_W"]
    P_coil_err = abs(P_coil_got - P_coil_exp) / P_coil_exp * 100
    assert P_coil_err < exp["tolerance_Pcoil_pct"], (
        f"P_coil regressed vs captured: got {P_coil_got:.4e}, "
        f"captured {P_coil_exp:.4e} (tol {exp['tolerance_Pcoil_pct']}%)."
    )
