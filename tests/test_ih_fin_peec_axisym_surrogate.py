"""Surrogate test of the adopted axisymmetric volume PEEC on the beak fin.

Guards: the native ring partial inductance, the loop-voltage circuit solve and
the regional loss split of ``radia.peec_axisym.solve_axisym_ring_peec`` on the
coarsest level of the adopted fin series (R, partial L, root and tip heat).
Assumes: surrogate test of adoption record induction_heating/20261116_fin_peec_axisym_p2 (convergence level 1).
The 0.06 mm meridian mesh is a stored fixture written by the validation
script's ``section_mesh``; meshing is not exercised here.  150 kHz copper,
radius 200 mm, 1 A terminal current.  The baseline checks historical output of
the production entry, not accuracy against the P2 reference; changed kernels
need reviewed regeneration, not tolerance relaxation.  Runtime is metadata.
APPROVAL-20261116-FIN-PEEC-SURROGATE-TOL-01 sets rtol=1e-6 and atol=0; relative
perturbations of 1e-5 in the inductance or conductivity fail it, 1e-6 survive.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

pm = pytest.importorskip("radia.peec_matrices")
if not hasattr(pm, "axisym_ring_inductance"):
    pytest.fail("radia.peec_matrices lacks axisym_ring_inductance: rebuild the native module")

from radia.peec_axisym import solve_axisym_ring_peec  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures" / "ih_fin"


def _lock():
    spec = importlib.util.spec_from_file_location("fin_regression_lock", ROOT / "tools" / "regression_lock.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_coarsest_fin_level_matches_adopted_history():
    mesh = json.loads((FIXTURES / "beak_section_r200_h060_mesh.json").read_text(encoding="utf-8"))
    points = np.array(mesh["points"], dtype=float)
    tris = np.array(mesh["triangles"], dtype=np.int32)
    region = np.array(mesh["region"])
    assert (len(points), len(tris)) == (1555, 2726), (len(points), len(tris))
    assert set(region) == {"body", "root", "tip"}, sorted(set(region))
    lock = _lock()
    baseline = lock.load_baseline(FIXTURES / "fin_peec_axisym_h060_76fb6b3a2.json")
    conditions = dict(baseline["conditions"])
    out = solve_axisym_ring_peec(points, tris, conditions["frequency_Hz"], conditions["sigma_S_per_m"], region=region)
    assert out["true_relative_residual"] <= 1e-6, out["true_relative_residual"]
    assert abs(out["power_balance"]) <= 1e-10, out["power_balance"]
    loss = out["region_loss_W_per_A2"]
    lock.assert_regression({"R_ohm": float(out["R_ohm"]), "L_H": float(out["L_H"]),
                            "root_heat_W_per_A2": float(loss["root"]), "tip_heat_W_per_A2": float(loss["tip"])},
                           baseline, conditions=conditions)
