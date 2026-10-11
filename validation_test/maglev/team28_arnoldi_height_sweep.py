"""TEAM 28 Arnoldi-Galerkin reduction across the 25-height sweep.

At every height of the full-FEM sweep the coil-driven eddy-current disk is
reduced by the block Krylov / congruence projection of
docs/maglev/demos/team28/team28_arnoldi_force.py (no circuit synthesis), and
its levitation force is compared with

* the direct solve of the same discrete system (the projection error),
* the lab reference force (an independent FEM, lab .mat Fz1), and
* the repository full-FEM sweep team28_axisym_fem_height_sweep.json.

The reduced order is chosen by the measured operating-band error: the
smallest order whose maximum force error against the direct solve is below
ORDER_SELECTION_N at all 25 heights.  The gates are the ones the retired
reduced sweep carried: reduced force within 1 mN of the lab reference at every
height, and the reduced equilibrium height within 0.6 mm of the published
stationary height (official Model A, 11.3 mm).

Run:  python validation_test/maglev/team28_arnoldi_height_sweep.py
"""

from __future__ import annotations

import hashlib
import json
import platform
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = REPO_ROOT / "docs" / "maglev" / "demos" / "team28" / "team28_arnoldi_force.py"
SWEEP = HERE / "team28_axisym_fem_height_sweep.json"
REFERENCE = HERE / "team28_reference.json"
DEFAULT_RESULT = HERE / "team28_arnoldi_height_sweep.json"
MAX_ORDER = 10
ORDER_SELECTION_N = 1.0e-4
LAB_GATE_N = 1.0e-3
PUBLISHED_HEIGHT_GATE_MM = 0.6

sys.path.insert(0, str(HERE))
from team28_axisym_fem_height_sweep import (  # noqa: E402
    DISK_BOTTOM_DZ0_MM, DISK_WEIGHT_N, DZ_MM, LAB_LEGACY_N, _equilibrium_mm)


def _sha(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def run():
    sys.path.insert(0, str(MODEL.parent))
    import ngsolve
    import scipy
    import team28_arnoldi_force as model

    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    published_mm = float(reference["stationary_height_mm"])
    sweep = json.loads(SWEEP.read_text(encoding="utf-8"))
    assert sweep["dZ_mm"] == list(DZ_MM)

    direct, by_order, seconds = [], [], []
    for dz in DZ_MM:
        start = time.perf_counter()
        fz_full, order_forces = model.arnoldi_forces(
            model.aluminium_z + dz * 1.0e-3, max_order=MAX_ORDER)
        seconds.append(time.perf_counter() - start)
        direct.append(fz_full)
        by_order.append(order_forces)
    direct = np.asarray(direct)
    orders = min(len(row) for row in by_order)
    reduced = np.asarray([row[:orders] for row in by_order])      # (height, order)
    order_error = np.abs(reduced - direct[:, None]).max(axis=0)   # per order
    eligible = [m for m in range(orders) if order_error[m] < ORDER_SELECTION_N]
    if not eligible:
        raise RuntimeError("no Arnoldi order reaches the selection error at every height")
    selected = eligible[0]                                         # 0-based
    legacy = reduced[:, selected]
    lab = np.asarray(LAB_LEGACY_N)
    full_fem = np.asarray(sweep["legacy_N"])
    upward = -0.5 * legacy
    equilibrium = _equilibrium_mm(DZ_MM, upward)
    lab_equilibrium = _equilibrium_mm(DZ_MM, -0.5 * lab)
    height_mm = DISK_BOTTOM_DZ0_MM + equilibrium
    checks = {
        "selected_order_below_selection_error_at_every_height":
            bool(order_error[selected] < ORDER_SELECTION_N),
        "reduced_vs_lab_max_abs_below_1mN": bool(np.abs(legacy - lab).max() < LAB_GATE_N),
        "reduced_vs_full_fem_sweep_max_abs_below_1mN":
            bool(np.abs(legacy - full_fem).max() < LAB_GATE_N),
        "reduced_equilibrium_within_0p6_mm_of_published":
            bool(abs(height_mm - published_mm) < PUBLISHED_HEIGHT_GATE_MM),
        "reduced_equilibrium_vs_lab_below_0p05_mm":
            bool(abs(equilibrium - lab_equilibrium) < 0.05),
        "lift_monotonically_nonincreasing": bool(np.all(np.diff(upward) <= 1.0e-9)),
    }
    return {
        "schema": "radia.team28-arnoldi-height-sweep.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform_class": platform.system(),
        "versions": {"python": platform.python_version(), "ngsolve": ngsolve.__version__,
                     "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): _sha(path)
            for path in (MODEL, Path(__file__).resolve(), REFERENCE, SWEEP)
        },
        "method": "axisymmetric mixed phi-B, p=2; block Krylov K_m(K^-1 N, K^-1 F) "
                  "orthonormalised by modified Gram-Schmidt and congruence projection; "
                  "one fresh mesh per height",
        "force_note": "legacy_N is the TEAM 28 integral Re[B_r J_t] = 2x the physical "
                      "time-averaged force; the lift is -legacy_N/2",
        "dZ_mm": list(DZ_MM),
        "direct_legacy_N": direct.tolist(),
        "reduced_legacy_N_by_order": reduced.T.tolist(),
        "max_abs_reduced_minus_direct_N_by_order": order_error.tolist(),
        "order_selection_N": ORDER_SELECTION_N,
        "selected_order": selected + 1,
        "selected_legacy_N": legacy.tolist(),
        "lab_legacy_N": list(LAB_LEGACY_N),
        "max_abs_reduced_minus_lab_N": float(np.abs(legacy - lab).max()),
        "max_abs_reduced_minus_full_fem_sweep_N": float(np.abs(legacy - full_fem).max()),
        "disk_weight_N": DISK_WEIGHT_N,
        "equilibrium_dZ_mm": equilibrium,
        "equilibrium_abs_height_mm": height_mm,
        "lab_equilibrium_dZ_mm": lab_equilibrium,
        "published_stationary_height_mm": published_mm,
        "solve_seconds": seconds,
        "checks": checks,
        "pass": bool(all(checks.values())),
    }


def main():
    started = time.perf_counter()
    result = run()
    result["elapsed_s"] = time.perf_counter() - started
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_RESULT
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "selected_order", "max_abs_reduced_minus_lab_N",
        "max_abs_reduced_minus_full_fem_sweep_N", "equilibrium_abs_height_mm",
        "pass")}, indent=2))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
