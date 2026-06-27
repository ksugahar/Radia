"""Validation-class spherical acoustic DtN-mode example.

This is an example/validation run, not a pytest test.  On a spherical
artificial boundary, each pressure trace mode ``Y_l^m`` has an exact outgoing
exterior Dirichlet-to-Neumann value:

    partial_n p / p = k h_l^(2)'(k a) / h_l^(2)(k a).

That is the smallest readable scalar acoustic FEM/BEM coupling gate: the FEM
side provides pressure on the boundary, and the BEM/radiation side provides the
normal derivative.  The degree-zero impedance also matches the uniformly
pulsating sphere formula.

Run:

    python validation_test/acoustic_bem/validation_spherical_dtn_modes.py
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
    pulsating_sphere_radiation,
    spherical_hankel2,
    spherical_helmholtz_dtn_eigenvalue,
    spherical_mode_radiation_impedance,
)


OUT_JSON = Path(__file__).with_name("validation_spherical_dtn_modes_summary.json")
RADIUS = 0.18
RHO = 1.2041
C = 343.0
DEGREES = [0, 1, 2, 3, 4]
KA_VALUES = [0.05, 0.1, 0.3, 1.0, 3.0, 20.0]


def _complex_record(value):
    z = complex(value)
    return {
        "real": z.real,
        "imag": z.imag,
        "abs": abs(z),
        "phase_rad": cmath.phase(z),
    }


def _frequency_for_ka(ka: float) -> float:
    return ka * C / (2.0 * math.pi * RADIUS)


def _finite_difference_dtn(degree: int, ka: float) -> complex:
    k = ka / RADIUS
    delta = 1.0e-6 * RADIUS
    h_boundary = spherical_hankel2(degree, ka)
    plus = spherical_hankel2(degree, k * (RADIUS + delta)) / h_boundary
    minus = spherical_hankel2(degree, k * (RADIUS - delta)) / h_boundary
    return (plus - minus) / (2.0 * delta)


def _row(degree: int, ka: float) -> dict:
    f = _frequency_for_ka(ka)
    mode = spherical_mode_radiation_impedance(RADIUS, f, degree, rho=RHO, c=C)
    dtn = spherical_helmholtz_dtn_eigenvalue(RADIUS, ka / RADIUS, degree)
    fd = _finite_difference_dtn(degree, ka)
    fd_abs_error = abs(dtn - fd)
    return {
        "degree": degree,
        "ka": ka,
        "frequency_hz": f,
        "dtn_eigenvalue": _complex_record(dtn),
        "specific_impedance": _complex_record(mode["specific_impedance"]),
        "radiation_efficiency": mode["radiation_efficiency"],
        "reactance_ratio": mode["reactance_ratio"],
        "finite_difference_dtn": _complex_record(fd),
        "finite_difference_abs_error": fd_abs_error,
        "finite_difference_rel_error": fd_abs_error / max(abs(dtn), 1.0e-300),
    }


def build_summary() -> dict:
    rows = [_row(degree, ka) for degree in DEGREES for ka in KA_VALUES]
    by_degree_ka = {(row["degree"], row["ka"]): row for row in rows}

    mono = by_degree_ka[(0, 1.0)]
    sphere = pulsating_sphere_radiation(RADIUS, _frequency_for_ka(1.0), 1.0, rho=RHO, c=C)
    mono_impedance_error = abs(
        complex(mono["specific_impedance"]["real"], mono["specific_impedance"]["imag"])
        - sphere["specific_impedance"]
    )

    low_frequency_ratios = []
    for degree in DEGREES:
        low = by_degree_ka[(degree, 0.05)]["radiation_efficiency"]
        high = by_degree_ka[(degree, 0.1)]["radiation_efficiency"]
        low_frequency_ratios.append({
            "degree": degree,
            "radiation_efficiency_ka_0p05": low,
            "radiation_efficiency_ka_0p1": high,
            "doubling_ratio": high / low,
        })

    high_ka = 80.0
    high_ka_errors = []
    for degree in DEGREES:
        k = high_ka / RADIUS
        dtn = spherical_helmholtz_dtn_eigenvalue(RADIUS, k, degree)
        sommerfeld_curvature = -1.0 / RADIUS - 1j * k
        high_ka_errors.append({
            "degree": degree,
            "relative_error_vs_minus_ik_minus_1_over_a": abs(dtn - sommerfeld_curvature)
            / abs(sommerfeld_curvature),
        })

    max_fd_rel_error = max(row["finite_difference_rel_error"] for row in rows if row["degree"] > 0)
    checks = {
        "monopole_dtn_real_abs_error": abs(
            by_degree_ka[(0, 1.0)]["dtn_eigenvalue"]["real"] - (-1.0 / RADIUS)
        ),
        "monopole_dtn_imag_abs_error": abs(
            by_degree_ka[(0, 1.0)]["dtn_eigenvalue"]["imag"] - (-1.0 / RADIUS)
        ),
        "monopole_impedance_abs_error_vs_pulsating_sphere": mono_impedance_error,
        "max_finite_difference_dtn_rel_error": max_fd_rel_error,
        "low_frequency_ratios_increase_with_degree": all(
            low_frequency_ratios[i + 1]["doubling_ratio"] > low_frequency_ratios[i]["doubling_ratio"]
            for i in range(len(low_frequency_ratios) - 1)
        ),
        "max_high_ka_sommerfeld_curvature_rel_error": max(
            row["relative_error_vs_minus_ik_minus_1_over_a"] for row in high_ka_errors
        ),
    }

    assert checks["monopole_dtn_real_abs_error"] < 1.0e-14
    assert checks["monopole_dtn_imag_abs_error"] < 1.0e-14
    assert checks["monopole_impedance_abs_error_vs_pulsating_sphere"] < 1.0e-12
    assert checks["max_finite_difference_dtn_rel_error"] < 1.0e-8
    assert checks["low_frequency_ratios_increase_with_degree"]
    assert checks["max_high_ka_sommerfeld_curvature_rel_error"] < 1.0e-2

    return {
        "kind": "spherical_helmholtz_dtn_mode_validation",
        "validation_class": True,
        "time_convention": "peak phasors with exp(+i omega t); outgoing h_l^(2)(k r)",
        "radius_m": RADIUS,
        "rho_kg_m3": RHO,
        "sound_speed_m_s": C,
        "degrees": DEGREES,
        "ka_values": KA_VALUES,
        "checks": checks,
        "low_frequency_ratios": low_frequency_ratios,
        "high_ka_sommerfeld_curvature_errors": high_ka_errors,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[spherical Helmholtz DtN modes]")
    print(
        f"  radius={summary['radius_m']:.3f} m degrees={summary['degrees']} "
        f"ka={summary['ka_values']}"
    )
    print(
        f"  monopole dtn errors real/imag="
        f"{summary['checks']['monopole_dtn_real_abs_error']:.3e}/"
        f"{summary['checks']['monopole_dtn_imag_abs_error']:.3e}"
    )
    print(
        f"  max FD rel error={summary['checks']['max_finite_difference_dtn_rel_error']:.3e}, "
        f"max high-ka Sommerfeld rel error="
        f"{summary['checks']['max_high_ka_sommerfeld_curvature_rel_error']:.3e}"
    )
    print(
        "  low-frequency resistance doubling ratios="
        + ", ".join(
            f"l={row['degree']}:{row['doubling_ratio']:.3f}"
            for row in summary["low_frequency_ratios"]
        )
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
