"""Independent impressed-source/Joule power gate, solid and bored conductors.

Explicit in-memory OCC/Netgen validation geometry (not a production mesh export).
This is an FEM self-consistency gate, NOT a BEM agreement certificate.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import ngsolve as ng
from reference_2d_axisym import solve_2d_axisym, EMMaterial


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coil-radius", type=float, default=0.003,
                        help="Circular wire cross-section radius, metres")
    args = parser.parse_args()
    ng.SetNumThreads(4)
    rows = []
    # Modest skin ratio first: refinement establishes the independent oracle
    # without a costly production-frequency mesh or external user data.
    with ng.TaskManager():
        for bore in (0.0, 0.0075):
            for order in (2, 3):
                result = solve_2d_axisym(
                    mat=EMMaterial.from_name("copper"), frequency=1000,
                    r_wp_inner=bore, order=order, skin_ratio=3,
                    maxh_air=0.015, sigma_coil=0.0, a_coil=args.coil_radius)
                result["order"] = order
                expected_volume = math.pi * (0.025**2 - bore**2) * 0.025
                result["volume_relative_error"] = abs(result["workpiece_volume"] / expected_volume - 1)
                rows.append(result)
                print(json.dumps(result), flush=True)
    errors = [row["power_balance_relative_error"] for row in rows]
    passed = all(row["P_joule_workpiece"] > 0 and row["volume_relative_error"] < 1e-10
                 for row in rows) and max(errors) < 1e-4
    sources = [Path(__file__), Path(__file__).with_name("reference_2d_axisym.py")]
    report = {
        "scope": "axisymmetric FEM impressed-source versus volume Joule power",
        "bem_validated": False,
        "phasor": "peak, exp(+i omega t)",
        "coil_cross_section_radius": args.coil_radius,
        "host": platform.node(), "python": sys.version,
        "ngsolve": ng.__version__,
        "sources_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "rows": rows, "power_balance_tolerance": 1e-4, "passed": passed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not passed:
        raise RuntimeError("FEM independent power balance failed; inspect JSON")


if __name__ == "__main__":
    main()
