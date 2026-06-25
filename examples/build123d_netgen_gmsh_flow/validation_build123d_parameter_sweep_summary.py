"""Validation-class build123d parameter sweep measurement summary.

Run:

    python examples/build123d_netgen_gmsh_flow/validation_build123d_parameter_sweep_summary.py

This example creates a small 2 x 3 x h box sweep and turns the measurement rows
into a design table before meshing or solver setup.  For this geometry,
``volume = 6 h`` and ``area = 12 + 10 h``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from build123d import Box


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.build123d.modeling import (  # noqa: E402
    shape_measurement_row,
    shape_parameter_sweep_summary,
)


OUT_JSON = HERE / "validation_build123d_parameter_sweep_summary.json"
HEIGHTS = [1.0, 2.0, 3.0, 4.0]


def _measurement_rows() -> list[dict[str, object]]:
    rows = []
    for height in HEIGHTS:
        shape = Box(2.0, 3.0, height).solid()
        shape.label = f"box_h_{height:g}"
        row = shape_measurement_row(shape)
        row["height"] = height
        row["expected_volume"] = 6.0 * height
        row["expected_area"] = 12.0 + 10.0 * height
        rows.append(row)
    return rows


def _assert_close(actual: float, expected: float, atol: float = 1.0e-12) -> float:
    error = abs(actual - expected)
    if error > atol:
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    rows = _measurement_rows()
    sweep = shape_parameter_sweep_summary(
        rows,
        "height",
        metric_keys=("volume", "area"),
        limits_by_metric={"volume": {"min": 12.0, "max": 24.0}},
    )
    clean_sweep = shape_parameter_sweep_summary(rows, "height", metric_keys=("volume", "area"))
    metrics = {row["metric"]: row for row in sweep["metric_rows"]}

    errors = {}
    for row in rows:
        height = row["height"]
        errors[f"volume_h_{height:g}"] = _assert_close(row["volume"], row["expected_volume"])
        errors[f"area_h_{height:g}"] = _assert_close(row["area"], row["expected_area"])

    checks = {
        "n_cases": sweep["n_cases"],
        "parameter_values": sweep["parameter_values"],
        "volume_min": metrics["volume"]["min"],
        "volume_max": metrics["volume"]["max"],
        "area_first": metrics["area"]["first"],
        "area_last": metrics["area"]["last"],
        "volume_monotonic": metrics["volume"]["monotonic_non_decreasing"],
        "area_monotonic": metrics["area"]["monotonic_non_decreasing"],
        "constrained_status": sweep["status"],
        "constraint_violation_count": sweep["constraint_violation_count"],
        "clean_status": clean_sweep["status"],
        "max_abs_error": max(errors.values()),
    }

    assert checks["n_cases"] == 4
    assert checks["parameter_values"] == HEIGHTS
    _assert_close(checks["volume_min"], 6.0)
    _assert_close(checks["volume_max"], 24.0)
    _assert_close(checks["area_first"], 22.0)
    _assert_close(checks["area_last"], 52.0)
    assert checks["volume_monotonic"] is True
    assert checks["area_monotonic"] is True
    assert checks["constrained_status"] == "needs_attention"
    assert checks["constraint_violation_count"] == 1
    assert checks["clean_status"] == "ok"

    return {
        "kind": "build123d_parameter_sweep_summary_validation",
        "validation_class": True,
        "learning_theme": (
            "CAD parameter sweeps should be reduced to measurement design tables "
            "before meshing or optimization"
        ),
        "checks": checks,
        "errors": errors,
        "rows": rows,
        "constrained_sweep": sweep,
        "clean_sweep": clean_sweep,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[build123d parameter sweep]")
    print(f"  cases={checks['n_cases']} heights={checks['parameter_values']}")
    print(f"  volume: min={checks['volume_min']:.12g} max={checks['volume_max']:.12g}")
    print(f"  area: first={checks['area_first']:.12g} last={checks['area_last']:.12g}")
    print(
        f"  constrained_status={checks['constrained_status']} "
        f"violations={checks['constraint_violation_count']}"
    )
    print(f"  clean_status={checks['clean_status']} max_abs_error={checks['max_abs_error']:.3e}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
