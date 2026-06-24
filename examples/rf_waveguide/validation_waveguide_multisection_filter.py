"""Validation-class rectangular-waveguide multi-section filter example.

This is an example/validation run, not a pytest test.  It exercises the TE10
ABCD cascade helper on a small but realistic microwave structure:

* one quarter-wave dielectric slab must match the closed-form slab helper
* a 3-period high/low permittivity quarter-wave stack behaves as a Bragg reflector
* lossless sections preserve unitarity and ABCD determinant one

Run:

    python examples/rf_waveguide/validation_waveguide_multisection_filter.py
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

from radia_mcp.radia_ngsolve.waveguide import (  # noqa: E402
    C0,
    reflection_metrics,
    waveguide_cascade_sparams,
    waveguide_dielectric_slab_sparams,
)


OUT_JSON = HERE / "validation_waveguide_multisection_filter_summary.json"
WIDTH_A = 0.02286
CENTER_FREQUENCY = 10.0e9
FREQUENCIES_HZ = [8.0e9 + 0.25e9 * i for i in range(17)]


def _beta_te10(frequency: float, width_a: float, eps_r: float) -> float:
    k0 = 2.0 * math.pi * frequency / C0
    kc = math.pi / width_a
    arg = eps_r * k0 * k0 - kc * kc
    if arg <= 0.0:
        raise ValueError("section is below TE10 cutoff for this validation")
    return math.sqrt(arg)


def quarter_wave_length(width_a: float, frequency: float, eps_r: float) -> float:
    return math.pi / (2.0 * _beta_te10(frequency, width_a, eps_r))


def _quarter_slab_record() -> dict:
    eps_r = 2.2
    length = quarter_wave_length(WIDTH_A, CENTER_FREQUENCY, eps_r)
    slab = waveguide_dielectric_slab_sparams(CENTER_FREQUENCY, WIDTH_A, eps_r, length)
    cascade = waveguide_cascade_sparams(CENTER_FREQUENCY, WIDTH_A, [(length, eps_r)])
    gamma = slab["gamma"]
    expected_peak = 2.0 * abs(gamma) / (1.0 + gamma * gamma)
    return {
        "eps_r": eps_r,
        "length_m": length,
        "theta_rad": slab["theta"],
        "interface_gamma": gamma,
        "slab_S11_mag": slab["S11_mag"],
        "expected_quarter_wave_S11_mag": expected_peak,
        "slab_vs_cascade_S11_abs_error": abs(slab["S11"] - cascade["S11"]),
        "slab_vs_cascade_S21_abs_error": abs(slab["S21"] - cascade["S21"]),
        "unitarity": slab["unitarity"],
    }


def bragg_stack_sections(repeats: int = 3, eps_hi: float = 4.0, eps_lo: float = 1.2) -> list[tuple[float, float]]:
    hi = quarter_wave_length(WIDTH_A, CENTER_FREQUENCY, eps_hi)
    lo = quarter_wave_length(WIDTH_A, CENTER_FREQUENCY, eps_lo)
    sections: list[tuple[float, float]] = []
    for _ in range(repeats):
        sections.append((hi, eps_hi))
        sections.append((lo, eps_lo))
    return sections


def _complex_pair(z: complex) -> dict:
    return {"real": z.real, "imag": z.imag, "abs": abs(z)}


def _bragg_records(sections: list[tuple[float, float]]) -> list[dict]:
    rows = []
    for f in FREQUENCIES_HZ:
        s = waveguide_cascade_sparams(f, WIDTH_A, sections)
        metrics = reflection_metrics(s["S11"])
        rows.append({
            "frequency_hz": f,
            "S11": _complex_pair(s["S11"]),
            "S21": _complex_pair(s["S21"]),
            "S11_mag": s["S11_mag"],
            "S21_mag": s["S21_mag"],
            "reflected_power": s["S11_mag"] ** 2,
            "transmitted_power": s["S21_mag"] ** 2,
            "return_loss_db": metrics["return_loss_db"],
            "delivered_power_fraction": metrics["delivered_power_fraction"],
            "unitarity": s["unitarity"],
            "abcd_det": _complex_pair(s["abcd_det"]),
        })
    return rows


def _validate(quarter: dict, rows: list[dict]) -> dict:
    center = min(rows, key=lambda r: abs(r["frequency_hz"] - CENTER_FREQUENCY))
    low_edge = rows[0]
    high_edge = rows[-1]
    max_unitarity_error = max(abs(row["unitarity"] - 1.0) for row in rows)
    max_det_error = max(abs(complex(row["abcd_det"]["real"], row["abcd_det"]["imag"]) - 1.0) for row in rows)
    checks = {
        "quarter_theta_error": abs(quarter["theta_rad"] - math.pi / 2.0),
        "quarter_formula_abs_error": abs(quarter["slab_S11_mag"] - quarter["expected_quarter_wave_S11_mag"]),
        "quarter_slab_vs_cascade_max_abs_error": max(
            quarter["slab_vs_cascade_S11_abs_error"],
            quarter["slab_vs_cascade_S21_abs_error"],
        ),
        "bragg_center_S11_mag": center["S11_mag"],
        "bragg_center_S21_mag": center["S21_mag"],
        "bragg_center_reflected_power": center["reflected_power"],
        "bragg_low_edge_S11_mag": low_edge["S11_mag"],
        "bragg_high_edge_S11_mag": high_edge["S11_mag"],
        "max_unitarity_error": max_unitarity_error,
        "max_abcd_det_error": max_det_error,
    }
    assert checks["quarter_theta_error"] < 1.0e-12
    assert checks["quarter_formula_abs_error"] < 1.0e-12
    assert checks["quarter_slab_vs_cascade_max_abs_error"] < 1.0e-12
    assert checks["bragg_center_S11_mag"] > 0.97
    assert checks["bragg_center_S21_mag"] < 0.21
    assert checks["bragg_center_reflected_power"] > 0.95
    assert checks["bragg_low_edge_S11_mag"] < checks["bragg_center_S11_mag"]
    assert checks["bragg_high_edge_S11_mag"] < checks["bragg_center_S11_mag"]
    assert max_unitarity_error < 1.0e-12
    assert max_det_error < 1.0e-12
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    quarter = _quarter_slab_record()
    sections = bragg_stack_sections()
    rows = _bragg_records(sections)
    checks = _validate(quarter, rows)
    summary = {
        "kind": "waveguide_multisection_filter_validation",
        "validation_class": True,
        "width_a_m": WIDTH_A,
        "center_frequency_hz": CENTER_FREQUENCY,
        "quarter_slab": quarter,
        "bragg_sections": [{"length_m": length, "eps_r": eps_r} for length, eps_r in sections],
        "frequency_rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[quarter-wave slab]")
    print(
        f"  eps_r={quarter['eps_r']:.3g}, length={quarter['length_m']*1e3:.6f} mm, "
        f"|S11|={quarter['slab_S11_mag']:.6f}, "
        f"closed-form={quarter['expected_quarter_wave_S11_mag']:.6f}"
    )
    print("[3-period Bragg stack]")
    for row in rows:
        print(
            f"  f={row['frequency_hz']/1e9:5.2f} GHz  "
            f"|S11|={row['S11_mag']:.6f}  |S21|={row['S21_mag']:.6f}  "
            f"unitarity={row['unitarity']:.12f}"
        )
    print(
        "[checks] "
        f"center |S11|={checks['bragg_center_S11_mag']:.6f}, "
        f"center |S21|={checks['bragg_center_S21_mag']:.6f}, "
        f"max det err={checks['max_abcd_det_error']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
