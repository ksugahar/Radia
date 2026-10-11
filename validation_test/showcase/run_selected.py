"""Recompute selected application data; docs notebooks only read these records."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import ngsolve as ng
import radia
from recording import runtime_metadata

ROOT = Path(__file__).resolve().parents[2]


def run(case):
    ng.SetNumThreads(2)
    if case == "motor":
        helper = ROOT / "docs/electric_machine/cogging_skew_demo.py"
        sys.path.insert(0, str(helper.parent))
        from cogging_skew_demo import run_demo
        data = run_demo(verbose=True)
        stack_length_m, length_unit_m, symmetry_factor = 0.050, 1e-3, 1.0
        data.update(stack_length_m=stack_length_m, length_unit_m=length_unit_m, symmetry_factor=symmetry_factor,
                    torque_Nm=(np.array(data["torque_2d_per_depth"])*length_unit_m**2*stack_length_m*symmetry_factor).tolist(),
                    scope="Two-dimensional torque ripple and skew ratios for a 50 mm stack; no drive controller.")
        if not all(data["checks"].values()):
            raise AssertionError(data["checks"])
    elif case == "complex_coil":
        helper = ROOT / "docs/complex_coil_geometry/coil_model.py"
        sys.path.insert(0, str(helper.parent))
        from coil_model import create_beam_steering_coil
        coil, parameters = create_beam_steering_coil()
        xs, ys = np.linspace(-0.40, 0.80, 55), np.linspace(-0.80, 0.80, 70)
        field = np.array([[radia.Fld(coil, "b", [float(x), float(y), 0.0])
                           for x in xs] for y in ys])
        probes = np.array([[-0.3, 0.1, 0.0], [0.0, 0.0, 0.0], [0.7, -0.1, 0.0]])
        base = np.array([radia.Fld(coil, "b", p.tolist()) for p in probes])
        radia.ObjScaleCur(coil, -1.0)
        reverse = np.array([radia.Fld(coil, "b", p.tolist()) for p in probes])
        error = float(np.linalg.norm(reverse + base) / np.linalg.norm(base))
        checks = {"finite_field": bool(np.isfinite(field).all()),
                  "nonzero_field": bool(np.linalg.norm(field) > 0),
                  "current_reversal": error < 1e-10}
        if not all(checks.values()):
            raise AssertionError(checks)
        data = {"parameters": parameters, "parameter_units": {"current": "A"}, "x_m": xs.tolist(), "y_m": ys.tolist(),
                "B_T": field.tolist(), "current_reversal_relative_error": error,
                "checks": checks, "scope": "Prescribed-current solid coil; reversal is a consistency check, not an independent accuracy estimate."}
    else:
        raise ValueError(case)
    return {"schema": "radia.selected_showcase.v1", "case": case,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "platform_class": platform.system(),
            **runtime_metadata(__file__, [helper], threads=2),
            "helper": helper.relative_to(ROOT).as_posix(),
            "helper_sha256": hashlib.sha256(helper.read_bytes()).hexdigest(),
            "result": data}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=["complex_coil", "motor"])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = run(args.case)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(args.output)
