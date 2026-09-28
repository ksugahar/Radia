"""Benchmark curved HCurl-VIM geometry against FEM and analytic sphere measures.

The sphere is remapped at several geometric orders while preserving one mesh
topology.  For every order this script records FEM area/volume, VIM sampled
surface area, surface-current tangentiality, and the conductor-cycle bridge
dual-volume measure.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SRC = REPO_ROOT / "src"
DEFAULT_OUTPUT = HERE / "curved_sphere_geometry_benchmark.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ngsolve as ng  # noqa: E402

from evrs_current_field_compare import _geometry_diagnostics  # noqa: E402
from evrs_pn_convergence import _parse_ints  # noqa: E402
from evrs_sibc_mixed_schur import _make_skin_mesh  # noqa: E402


def run(args: argparse.Namespace) -> dict[str, object]:
    started = time.perf_counter()
    cases = []
    for curve_order in args.curve_orders:
        mesh = _make_skin_mesh(args.maxh, "sphere", curve_order)
        case = _geometry_diagnostics(mesh, "sphere", args.intorder)
        case["requested_curve_order"] = curve_order
        cases.append(case)

    bridge_reference = cases[-1]["bridge_dual_volume_sum"]
    for case in cases:
        case["bridge_dual_volume_relative_to_highest_order"] = float(
            abs(case["bridge_dual_volume_sum"] - bridge_reference)
            / abs(bridge_reference)
        )

    return {
        "schema": "radia.validation.curved_sphere_geometry.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "configuration": {
            "curve_orders": args.curve_orders,
            "maxh": args.maxh,
            "intorder": args.intorder,
        },
        "total_wall_seconds": float(time.perf_counter() - started),
        "cases": cases,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-orders", type=_parse_ints, default=[1, 2, 3, 4])
    parser.add_argument("--maxh", type=float, default=2.0)
    parser.add_argument("--intorder", type=int, default=10)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if any(order < 1 for order in args.curve_orders):
        parser.error("--curve-orders must contain positive integers")
    if args.maxh <= 0.0:
        parser.error("--maxh must be positive")
    if args.intorder < 0:
        parser.error("--intorder must be non-negative")

    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("Curved sphere FEM/VIM geometry benchmark")
    print(f"  output: {args.output}")
    print("  q  elements  area error  volume error  tangent defect  bridge delta")
    for case in result["cases"]:
        print(
            f"  {case['requested_curve_order']:>1} "
            f"{case['elements']:>9} "
            f"{case['surface_area_relative_error_to_exact']:.3e} "
            f"{case['volume_relative_error_to_exact']:.3e} "
            f"{case['surface_tangential_defect']:.3e} "
            f"{case['bridge_dual_volume_relative_to_highest_order']:.3e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
