"""Validation-class low-frequency Helmholtz Green-kernel example.

This is an example/validation run, not a pytest test.  It records the analytic
split used by low-frequency acoustic BEM:

    exp(-i k r)/(4 pi r)
      = 1/(4 pi r) + regular smooth corrections.

The static Laplace singularity is isolated, while the frequency corrections are
regular panel integrals.  This is the small kernel-level gate before building a
readable FEM/BEM acoustic coupling example.

Run:

    python examples/acoustic_bem/validation_low_frequency_helmholtz_kernel.py
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import (  # noqa: E402
    helmholtz_green_3d,
    helmholtz_green_low_frequency_series,
)


OUT_JSON = Path(__file__).with_name("validation_low_frequency_helmholtz_kernel_summary.json")
DISTANCE = 0.75
KR_VALUES = [1.0e-6, 1.0e-4, 1.0e-2, 5.0e-2, 1.0e-1, 2.0e-1, 5.0e-1]
ORDERS = [0, 1, 2, 4, 6, 8]


def _complex_record(value):
    z = complex(value)
    return {
        "real": z.real,
        "imag": z.imag,
        "abs": abs(z),
        "phase_rad": cmath.phase(z),
    }


def _case(kr: float) -> dict:
    k = kr / DISTANCE
    exact = helmholtz_green_3d(DISTANCE, k)
    by_order = []
    previous_error = None
    monotone_from_order2 = True
    for order in ORDERS:
        out = helmholtz_green_low_frequency_series(DISTANCE, k, order=order)
        err = out["abs_error"]
        if previous_error is not None and order >= 2 and err > previous_error:
            monotone_from_order2 = False
        previous_error = err
        by_order.append({
            "order": order,
            "approx": _complex_record(out["approx"]),
            "abs_error": err,
            "rel_error": err / max(abs(exact), 1.0e-300),
        })

    terms = helmholtz_green_low_frequency_series(DISTANCE, k, order=3)["terms"]
    return {
        "kr": kr,
        "distance_m": DISTANCE,
        "wavenumber_per_m": k,
        "exact": _complex_record(exact),
        "laplace_term": _complex_record(terms[0]),
        "first_regular_term": _complex_record(terms[1]),
        "second_regular_term": _complex_record(terms[2]),
        "third_regular_term": _complex_record(terms[3]),
        "by_order": by_order,
        "monotone_from_order2": monotone_from_order2,
    }


def build_summary():
    rows = [_case(kr) for kr in KR_VALUES]
    max_order6_error_small = max(
        row["by_order"][ORDERS.index(6)]["abs_error"]
        for row in rows
        if row["kr"] <= 0.2
    )
    max_order8_rel_error = max(row["by_order"][ORDERS.index(8)]["rel_error"] for row in rows)
    low = rows[0]
    k_low = low["wavenumber_per_m"]
    checks = {
        "distance_m": DISTANCE,
        "laplace_term_real": low["laplace_term"]["real"],
        "laplace_term_expected": 1.0 / (4.0 * math.pi * DISTANCE),
        "first_regular_imag_over_minus_k_over_4pi": (
            low["first_regular_term"]["imag"] / (-k_low / (4.0 * math.pi))
        ),
        "second_regular_real_over_minus_k2r_over_8pi": (
            low["second_regular_term"]["real"] / (-(k_low * k_low) * DISTANCE / (8.0 * math.pi))
        ),
        "max_order6_abs_error_for_kr_le_0p2": max_order6_error_small,
        "max_order8_rel_error_for_sweep": max_order8_rel_error,
        "all_monotone_from_order2": all(row["monotone_from_order2"] for row in rows),
    }

    assert abs(checks["laplace_term_real"] - checks["laplace_term_expected"]) < 1.0e-15
    assert abs(checks["first_regular_imag_over_minus_k_over_4pi"] - 1.0) < 1.0e-15
    assert abs(checks["second_regular_real_over_minus_k2r_over_8pi"] - 1.0) < 1.0e-15
    assert checks["max_order6_abs_error_for_kr_le_0p2"] < 1.0e-9
    assert checks["max_order8_rel_error_for_sweep"] < 1.0e-8
    assert checks["all_monotone_from_order2"]

    return {
        "kind": "low_frequency_helmholtz_green_kernel_validation",
        "validation_class": True,
        "time_convention": "peak phasors with exp(+i omega t); outgoing exp(-i k r)",
        "checks": checks,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[low-frequency Helmholtz Green kernel]")
    print(
        f"  r={summary['checks']['distance_m']:.6f} m, "
        f"Laplace={summary['checks']['laplace_term_real']:.15e}"
    )
    print(
        f"  first regular ratio={summary['checks']['first_regular_imag_over_minus_k_over_4pi']:.12f}, "
        f"second regular ratio={summary['checks']['second_regular_real_over_minus_k2r_over_8pi']:.12f}"
    )
    print(
        f"  max order-6 abs error (kr<=0.2)="
        f"{summary['checks']['max_order6_abs_error_for_kr_le_0p2']:.3e}, "
        f"max order-8 rel error={summary['checks']['max_order8_rel_error_for_sweep']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
