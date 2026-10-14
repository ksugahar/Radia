"""Validation-class sphere-over-ground-plane capacitance example.

This example exercises the image-series capacitance helper on a sweep that is
large enough to catch truncation mistakes:

* near the plane, the image series converges slowly and C rises well above the
  isolated-sphere value;
* far from the plane, C approaches 4*pi*eps0*a with the first correction a/(2h);
* the electrostatic attraction inferred from -dC/dh is positive.

Run:

    python validation_test/electrostatics/validation_sphere_ground_plane_capacitance.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.electrostatics import (  # noqa: E402
    EPS0,
    sphere_above_plane_capacitance,
)
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_sphere_ground_plane_capacitance_summary.json"
RADIUS_A = 0.01
EPS_R = 1.0
VOLTAGE = 100.0
HEIGHT_RATIOS = (1.05, 1.10, 1.25, 1.50, 2.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0)
TERM_COUNTS = (5, 10, 20, 40, 80)
REFERENCE_TERMS = 1500


def _image_term(alpha: float, n: int) -> float:
    x = n * alpha
    numerator = math.sinh(alpha)
    if x < 40.0:
        return numerator / math.sinh(x)
    return 2.0 * numerator * math.exp(-x)


def reference_capacitance(radius_a: float, height_h: float, eps_r: float = 1.0,
                          n_terms: int = REFERENCE_TERMS) -> float:
    """Overflow-safe high-term image-series reference."""
    if not height_h > radius_a > 0.0:
        raise ValueError("require height_h > radius_a > 0")
    alpha = math.acosh(height_h / radius_a)
    series = sum(_image_term(alpha, n) for n in range(1, n_terms + 1))
    return 4.0 * math.pi * EPS0 * eps_r * radius_a * series


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


def _finite_difference_force(radius_a: float, height_h: float, eps_r: float, voltage: float) -> dict:
    # At fixed voltage, W = 1/2 C V^2.  The attractive force magnitude is
    # -dW/dh = -1/2 V^2 dC/dh, positive because C decreases as h grows.
    step = max(1.0e-5 * height_h, 1.0e-6 * radius_a)
    c_minus = reference_capacitance(radius_a, height_h - step, eps_r)
    c_plus = reference_capacitance(radius_a, height_h + step, eps_r)
    dcdh = (c_plus - c_minus) / (2.0 * step)
    return {
        "height_m": height_h,
        "step_m": step,
        "dCdh_F_per_m": dcdh,
        "attractive_force_N": -0.5 * voltage * voltage * dcdh,
    }


def build_rows() -> list[dict]:
    c_isolated = 4.0 * math.pi * EPS0 * EPS_R * RADIUS_A
    rows = []
    for ratio in HEIGHT_RATIOS:
        height = ratio * RADIUS_A
        cref = reference_capacitance(RADIUS_A, height, EPS_R)
        by_terms = {}
        rel_errors = {}
        for n_terms in TERM_COUNTS:
            value = sphere_above_plane_capacitance(
                RADIUS_A,
                height,
                eps_r=EPS_R,
                n_terms=n_terms,
            )["C"]
            by_terms[str(n_terms)] = value
            rel_errors[str(n_terms)] = _relative_error(value, cref)
        rows.append({
            "height_over_radius": ratio,
            "height_m": height,
            "reference_C_F": cref,
            "C_over_isolated": cref / c_isolated,
            "excess_over_isolated": cref / c_isolated - 1.0,
            "capacitance_by_terms_F": by_terms,
            "relative_errors_by_terms": rel_errors,
        })
    return rows


def _far_field_checks(rows: list[dict]) -> list[dict]:
    checks = []
    for row in rows:
        ratio = row["height_over_radius"]
        if ratio < 25.0:
            continue
        asymptotic = 1.0 / (2.0 * ratio)
        actual = row["excess_over_isolated"]
        checks.append({
            "height_over_radius": ratio,
            "actual_excess": actual,
            "first_image_asymptotic": asymptotic,
            "relative_mismatch": abs(actual - asymptotic) / asymptotic,
        })
    return checks


def validate(rows: list[dict]) -> dict:
    normalized = [row["C_over_isolated"] for row in rows]
    force = _finite_difference_force(RADIUS_A, 2.0 * RADIUS_A, EPS_R, VOLTAGE)
    far_field = _far_field_checks(rows)
    checks = {
        "near_plane_C_over_isolated": normalized[0],
        "mid_sweep_C_over_isolated": next(row["C_over_isolated"] for row in rows
                                          if row["height_over_radius"] == 2.0),
        "far_C_over_isolated": normalized[-1],
        "monotone_decrease_with_height": all(a > b for a, b in zip(normalized, normalized[1:])),
        "max_rel_error_20_terms": max(row["relative_errors_by_terms"]["20"] for row in rows),
        "max_rel_error_40_terms": max(row["relative_errors_by_terms"]["40"] for row in rows),
        "max_rel_error_80_terms": max(row["relative_errors_by_terms"]["80"] for row in rows),
        "far_field": far_field,
        "force_at_height_over_radius_2": force,
    }
    assert checks["near_plane_C_over_isolated"] > 2.4
    assert checks["mid_sweep_C_over_isolated"] > 1.3
    assert abs(checks["far_C_over_isolated"] - 1.0) < 6.0e-3
    assert checks["monotone_decrease_with_height"]
    assert checks["max_rel_error_20_terms"] < 1.5e-3
    assert checks["max_rel_error_40_terms"] < 3.0e-6
    assert checks["max_rel_error_80_terms"] < 1.0e-10
    assert max(item["relative_mismatch"] for item in far_field) < 3.0e-2
    assert force["dCdh_F_per_m"] < 0.0
    assert force["attractive_force_N"] > 0.0
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "sphere_ground_plane_capacitance_validation",
        "validation_class": True,
        "radius_m": RADIUS_A,
        "eps_r": EPS_R,
        "reference_terms": REFERENCE_TERMS,
        "term_counts": list(TERM_COUNTS),
        "rows": rows,
        "checks": checks,
    }
    summary = add_result_metadata(summary, __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[sphere above grounded plane capacitance]")
    for row in rows:
        print(
            f"  h/a={row['height_over_radius']:6.2f}  "
            f"C/C_iso={row['C_over_isolated']:.9f}  "
            f"err20={row['relative_errors_by_terms']['20']:.3e}  "
            f"err40={row['relative_errors_by_terms']['40']:.3e}"
        )
    force = checks["force_at_height_over_radius_2"]
    print(
        "[checks] "
        f"near C/C_iso={checks['near_plane_C_over_isolated']:.6f}, "
        f"far C/C_iso={checks['far_C_over_isolated']:.6f}, "
        f"F(h/a=2,V={VOLTAGE:g})={force['attractive_force_N']:.6e} N"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
