"""Validation-class planar acoustic DtN-symbol example.

This example records the exact half-space Helmholtz DtN symbol for planar
trace modes.  It is the flat-boundary companion to the spherical DtN example:
the FEM side supplies pressure on the boundary, and the exterior BEM/radiation
side supplies ``partial_n p = lambda p``.

Run:

    python examples/acoustic_bem/validation_planar_dtn_symbol.py
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
    planar_helmholtz_dtn_symbol,
    planar_mode_radiation_impedance,
)


OUT_JSON = Path(__file__).with_name("validation_planar_dtn_symbol_summary.json")
FREQUENCY = 1000.0
RHO = 1.2041
C = 343.0
INCIDENCE_DEGREES = [0.0, 15.0, 30.0, 60.0, 80.0]
EVANESCENT_KT_OVER_K = [1.05, 1.5, 3.0]


def _complex_record(value):
    z = complex(value)
    return {
        "real": z.real,
        "imag": z.imag,
        "abs": abs(z),
        "phase_rad": cmath.phase(z),
    }


def _propagating_row(degrees: float) -> dict:
    theta = math.radians(degrees)
    out = planar_mode_radiation_impedance(
        FREQUENCY,
        incidence_angle_rad=theta,
        rho=RHO,
        c=C,
    )
    expected = 1.0 / math.cos(theta)
    return {
        "incidence_degrees": degrees,
        "kt_over_k": math.sin(theta),
        "normal_wavenumber": _complex_record(out["normal_wavenumber"]),
        "dtn_eigenvalue": _complex_record(out["dtn_eigenvalue"]),
        "normalized_impedance": _complex_record(out["normalized_impedance"]),
        "expected_normalized_impedance": expected,
        "normalized_impedance_abs_error": abs(out["normalized_impedance"] - expected),
    }


def _evanescent_row(kt_over_k: float) -> dict:
    k = 2.0 * math.pi * FREQUENCY / C
    out = planar_mode_radiation_impedance(
        FREQUENCY,
        tangential_wavenumber=kt_over_k * k,
        rho=RHO,
        c=C,
    )
    expected_reactance = 1.0 / math.sqrt(kt_over_k * kt_over_k - 1.0)
    return {
        "kt_over_k": kt_over_k,
        "normal_wavenumber": _complex_record(out["normal_wavenumber"]),
        "dtn_eigenvalue": _complex_record(out["dtn_eigenvalue"]),
        "normalized_impedance": _complex_record(out["normalized_impedance"]),
        "expected_reactance_ratio": expected_reactance,
        "reactance_abs_error": abs(out["reactance_ratio"] - expected_reactance),
        "radiation_efficiency": out["radiation_efficiency"],
    }


def build_summary():
    k = 2.0 * math.pi * FREQUENCY / C
    propagating = [_propagating_row(deg) for deg in INCIDENCE_DEGREES]
    evanescent = [_evanescent_row(ratio) for ratio in EVANESCENT_KT_OVER_K]

    symbol_rows = []
    for kt_over_k in [0.0, 0.25, 0.75, 1.05, 2.0, 4.0]:
        symbol = planar_helmholtz_dtn_symbol(k, kt_over_k * k)
        residual = symbol["symbol_identity_residual"]
        symbol_rows.append({
            "kt_over_k": kt_over_k,
            "regime": symbol["regime"],
            "dtn_eigenvalue": _complex_record(symbol["dtn_eigenvalue"]),
            "identity_residual_abs": abs(residual),
        })

    checks = {
        "frequency_hz": FREQUENCY,
        "wavenumber_per_m": k,
        "normal_incidence_impedance_error": propagating[0]["normalized_impedance_abs_error"],
        "max_propagating_impedance_error": max(row["normalized_impedance_abs_error"] for row in propagating),
        "max_evanescent_reactance_error": max(row["reactance_abs_error"] for row in evanescent),
        "max_symbol_identity_residual": max(row["identity_residual_abs"] for row in symbol_rows),
        "all_evanescent_zero_active_power": all(
            abs(row["radiation_efficiency"]) < 1.0e-14 for row in evanescent
        ),
        "all_evanescent_dtn_real_negative": all(
            row["dtn_eigenvalue"]["real"] < 0.0 and abs(row["dtn_eigenvalue"]["imag"]) < 1.0e-14
            for row in evanescent
        ),
    }

    assert checks["normal_incidence_impedance_error"] < 1.0e-14
    assert checks["max_propagating_impedance_error"] < 1.0e-13
    assert checks["max_evanescent_reactance_error"] < 1.0e-13
    assert checks["max_symbol_identity_residual"] < 1.0e-12
    assert checks["all_evanescent_zero_active_power"]
    assert checks["all_evanescent_dtn_real_negative"]

    return {
        "kind": "planar_helmholtz_dtn_symbol_validation",
        "validation_class": True,
        "time_convention": "peak phasors with exp(+i omega t); outgoing exp(-i q n)",
        "checks": checks,
        "propagating": propagating,
        "evanescent": evanescent,
        "symbol_rows": symbol_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[planar Helmholtz DtN symbol]")
    print(
        f"  f={summary['checks']['frequency_hz']:.1f} Hz, "
        f"k={summary['checks']['wavenumber_per_m']:.9f} 1/m"
    )
    print(
        f"  max propagating impedance error="
        f"{summary['checks']['max_propagating_impedance_error']:.3e}, "
        f"max evanescent reactance error="
        f"{summary['checks']['max_evanescent_reactance_error']:.3e}"
    )
    print(
        f"  max DtN identity residual="
        f"{summary['checks']['max_symbol_identity_residual']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
