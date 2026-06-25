"""Validation-class two-port scattering sweep momentum force.

Run:

    python examples/rf_waveguide/validation_two_port_scattering_sweep_force.py

This example audits a small frequency sweep of power-normalized two-port
S-parameters.  For straight-through propagation, the axial momentum force is

    F = P_inc (1 + |S11|^2 - |S21|^2) / c.

The validation checks force extrema, mean force, absorbed power fraction, and a
separate passivity probe with ``|S11|^2 + |S21|^2 > 1``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    two_port_scattering_sweep_momentum_force_summary,
)


OUT_JSON = HERE / "validation_two_port_scattering_sweep_force_summary.json"


def _assert_close(actual: float, expected: float, atol: float = 1.0e-21) -> float:
    error = abs(actual - expected)
    if error > atol:
        raise AssertionError(f"{actual!r} != {expected!r}; error={error!r}")
    return error


def build_summary() -> dict[str, object]:
    power_W = 3.0
    frequencies = [3.0e9, 4.0e9, 5.0e9, 6.0e9]
    s11 = [0.05 + 0.0j, 0.20 + 0.0j, 0.45 + 0.0j, 0.10 + 0.0j]
    s21 = [0.95 + 0.0j, 0.80 + 0.0j, 0.55 + 0.0j, 0.90 + 0.0j]
    sweep = two_port_scattering_sweep_momentum_force_summary(
        frequencies,
        s11,
        s21,
        power_incident_W=power_W,
        incident_direction=(1.0, 0.0, 0.0),
        transmitted_direction=(1.0, 0.0, 0.0),
    )
    passivity_probe = two_port_scattering_sweep_momentum_force_summary(
        [7.0e9],
        [0.8 + 0.0j],
        [0.7 + 0.0j],
        power_incident_W=power_W,
        passivity_tolerance=1.0e-6,
    )

    expected_factors = [
        1.0 + abs(gamma) ** 2 - abs(tau) ** 2
        for gamma, tau in zip(s11, s21)
    ]
    expected_forces = [power_W * factor / C0 for factor in expected_factors]
    max_expected = max(expected_forces)
    min_expected = min(expected_forces)
    mean_expected = sum(expected_forces) / len(expected_forces)
    force_errors = [
        abs(row["force_magnitude_N"] - expected)
        for row, expected in zip(sweep["rows"], expected_forces)
    ]

    checks = {
        "n_points": sweep["n_points"],
        "frequency_monotonic_increasing": sweep["frequency_monotonic_increasing"],
        "passive_status": sweep["status"],
        "max_force_frequency_Hz": sweep["max_force_frequency_Hz"],
        "min_force_frequency_Hz": sweep["min_force_frequency_Hz"],
        "max_force_magnitude_N": sweep["max_force_magnitude_N"],
        "min_force_magnitude_N": sweep["min_force_magnitude_N"],
        "mean_force_magnitude_N": sweep["mean_force_magnitude_N"],
        "expected_max_force_N": max_expected,
        "expected_min_force_N": min_expected,
        "expected_mean_force_N": mean_expected,
        "max_abs_force_error_N": max(force_errors),
        "passivity_probe_status": passivity_probe["status"],
        "passivity_probe_excess": passivity_probe["max_passivity_excess_power_fraction"],
    }

    assert checks["n_points"] == 4
    assert checks["frequency_monotonic_increasing"] is True
    assert checks["passive_status"] == "ok"
    assert checks["max_force_frequency_Hz"] == 5.0e9
    assert checks["min_force_frequency_Hz"] == 3.0e9
    _assert_close(checks["max_force_magnitude_N"], max_expected)
    _assert_close(checks["min_force_magnitude_N"], min_expected)
    _assert_close(checks["mean_force_magnitude_N"], mean_expected)
    assert checks["max_abs_force_error_N"] < 1.0e-21
    assert checks["passivity_probe_status"] == "needs_attention"
    assert checks["passivity_probe_excess"] > 0.0

    return {
        "kind": "two_port_scattering_sweep_force_validation",
        "validation_class": True,
        "rf_learning": (
            "two-port S-parameter sweeps can be audited as momentum-flux tables "
            "with passivity diagnostics"
        ),
        "checks": checks,
        "sweep": sweep,
        "passivity_probe": passivity_probe,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[two-port scattering sweep force]")
    print(f"  n_points={checks['n_points']} status={checks['passive_status']}")
    print(
        f"  min_force={checks['min_force_magnitude_N']:.12g} N "
        f"at {checks['min_force_frequency_Hz']:.12g} Hz"
    )
    print(
        f"  max_force={checks['max_force_magnitude_N']:.12g} N "
        f"at {checks['max_force_frequency_Hz']:.12g} Hz"
    )
    print(
        "  passivity_probe="
        f"{checks['passivity_probe_status']} excess={checks['passivity_probe_excess']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
