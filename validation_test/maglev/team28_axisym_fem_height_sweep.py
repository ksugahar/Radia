"""TEAM 28 full-FEM force-height sweep (the reduced-model-free reference curve).

Every height is a fresh axisymmetric mesh and a full mixed phi-B solve of
docs/maglev/demos/team28/team28_axisym_fem.py; no reduced model is involved.
The stored lab regression (legacy integral, lab .mat Fz1) is an independent
second FEM route on the same heights.

Run:  python validation_test/maglev/team28_axisym_fem_height_sweep.py
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
MODEL = REPO_ROOT / "docs" / "maglev" / "demos" / "team28" / "team28_axisym_fem.py"
DEFAULT_RESULT = HERE / "team28_axisym_fem_height_sweep.json"
DZ_MM = tuple(range(-7, 18))
DISK_WEIGHT_N = 1.055
DISK_BOTTOM_DZ0_MM = 10.8
# Lab .mat Fz1 (legacy integral Re[B_r J_t], N) at DZ_MM; independent lab FEM.
LAB_LEGACY_N = (
    -6.5794, -5.6713, -4.877, -4.1833, -3.578, -3.0505, -2.5916, -2.1928,
    -1.8469, -1.5475, -1.2887, -1.0655, -0.8736, -0.709, -0.5683, -0.4483,
    -0.3466, -0.2606, -0.1883, -0.1279, -0.0779, -0.0367, -0.0033, 0.0236,
    0.0447,
)


def _equilibrium_mm(dz_mm, upward_N, weight_N=DISK_WEIGHT_N):
    residual = np.asarray(upward_N) - weight_N
    for i in range(len(dz_mm) - 1):
        if residual[i] * residual[i + 1] <= 0.0:
            frac = -residual[i] / (residual[i + 1] - residual[i])
            return float(dz_mm[i] + frac * (dz_mm[i + 1] - dz_mm[i]))
    raise RuntimeError("the sweep does not bracket lift equals weight")


def run():
    sys.path.insert(0, str(MODEL.parent))
    import ngsolve
    import team28_axisym_fem as model

    legacy, physical, seconds = [], [], []
    for dz in DZ_MM:
        start = time.perf_counter()
        f_legacy, f_phys = model.solve_force_pair(dz * 1.0e-3)
        seconds.append(time.perf_counter() - start)
        legacy.append(f_legacy)
        physical.append(f_phys)
    legacy = np.asarray(legacy)
    # the physical time-averaged force is positive upward; the legacy integral
    # uses the opposite sign (negative = lift)
    upward = np.asarray(physical)
    lab = np.asarray(LAB_LEGACY_N)
    twice_residual = np.abs(legacy + 2.0 * upward)
    lab_abs = np.abs(legacy - lab)
    equilibrium = _equilibrium_mm(DZ_MM, upward)
    lab_equilibrium = _equilibrium_mm(DZ_MM, -0.5 * lab)
    checks = {
        "legacy_vs_lab_max_abs_below_1mN": bool(lab_abs.max() < 1.0e-3),
        # the ratio is ill-defined near the force zero crossing (dZ~15 mm), so
        # bound the residual of legacy = -2 physical against the peak force
        "legacy_is_minus_twice_physical_within_5e-4_of_peak": bool(
            twice_residual.max() < 5.0e-4 * np.abs(legacy).max()),
        "lift_monotonically_nonincreasing": bool(np.all(np.diff(upward) <= 1.0e-9)),
        "lift_equals_weight_bracketed": bool(upward[0] > DISK_WEIGHT_N > upward[-1]),
        "equilibrium_vs_lab_below_0p05_mm": bool(abs(equilibrium - lab_equilibrium) < 0.05),
    }
    return {
        "schema": "radia.team28-axisym-fem-height-sweep.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "versions": {
            "python": platform.python_version(),
            "ngsolve": ngsolve.__version__,
            "numpy": np.__version__,
        },
        "source_sha256": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"):
                hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
            for path in (MODEL, Path(__file__).resolve())
        },
        "method": "axisymmetric mixed phi-B full FEM, p=2, anisotropic-nu "
                  "infinite-element shell, one fresh mesh per height",
        "height_datum": "dZ offset from the 10.8 mm disk-bottom position",
        "force_note": "legacy_N is the TEAM 28 integral Re[B_r J_t] = 2x the physical "
                      "time-averaged force; upward_physical_N is the physical lift",
        "frequency_Hz": float(model.FREQ),
        "coil_current_peak_A": float(model.I1),
        "disk_weight_N": DISK_WEIGHT_N,
        "dZ_mm": list(DZ_MM),
        "legacy_N": legacy.tolist(),
        "upward_physical_N": upward.tolist(),
        "lab_legacy_N": list(LAB_LEGACY_N),
        "max_abs_legacy_minus_lab_N": float(lab_abs.max()),
        "max_abs_legacy_plus_twice_physical_N": float(twice_residual.max()),
        "equilibrium_dZ_mm": equilibrium,
        "equilibrium_abs_height_mm": DISK_BOTTOM_DZ0_MM + equilibrium,
        "lab_equilibrium_dZ_mm": lab_equilibrium,
        "solve_seconds": seconds,
        "checks": checks,
        "pass": bool(all(checks.values())),
    }


def main():
    started = time.perf_counter()
    result = run()
    result["elapsed_s"] = time.perf_counter() - started
    DEFAULT_RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("checks", "pass", "equilibrium_dZ_mm",
                                             "max_abs_legacy_minus_lab_N", "elapsed_s")},
                     indent=2))
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
