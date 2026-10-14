"""Validation-class acoustic impedance reflection sweep.

Run:

    python validation_test/acoustic_bem/validation_acoustic_impedance_sweep.py

For a local acoustic impedance load under normal incidence,

    Gamma = (Z - Z0) / (Z + Z0),
    A = 1 - |Gamma|^2,
    p_mom = (1 + |Gamma|^2) I_inc / c.

The example records absorption and momentum-pressure extrema over a small
impedance sweep, plus a separate non-passive probe.
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

from radia_mcp.radia_ngsolve.acoustics import acoustic_impedance_reflection_sweep_summary  # noqa: E402


OUT_JSON = HERE / "validation_acoustic_impedance_sweep_summary.json"


def build_summary() -> dict[str, object]:
    rho = 1.2041
    c = 343.0
    z0 = rho * c
    area = 0.3
    incident_pressure = 2.0
    frequencies = [125.0, 250.0, 500.0, 1000.0]
    impedances = [z0, 2.0 * z0, z0 * (1.0 + 0.5j), 1.0j * z0]
    sweep = acoustic_impedance_reflection_sweep_summary(
        frequencies,
        impedances,
        area=area,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
    )
    passivity_probe = acoustic_impedance_reflection_sweep_summary(
        [2000.0],
        [-2.0 * z0],
        rho=rho,
        c=c,
        passivity_tolerance=1.0e-6,
    )

    checks = {
        "n_points": sweep["n_points"],
        "frequency_monotonic_increasing": sweep["frequency_monotonic_increasing"],
        "status": sweep["status"],
        "max_absorption_frequency_Hz": sweep["max_absorption_frequency_Hz"],
        "min_absorption_frequency_Hz": sweep["min_absorption_frequency_Hz"],
        "max_force_frequency_Hz": sweep["max_force_frequency_Hz"],
        "min_force_frequency_Hz": sweep["min_force_frequency_Hz"],
        "matched_absorption": sweep["rows"][0]["absorption_coefficient"],
        "double_resistance_reflectance": sweep["rows"][1]["power_reflection_coefficient"],
        "reactive_reflectance": sweep["rows"][3]["power_reflection_coefficient"],
        "passivity_probe_status": passivity_probe["status"],
        "passivity_probe_excess_absorption": passivity_probe["max_passivity_excess_absorption"],
    }

    assert checks["n_points"] == 4
    assert checks["frequency_monotonic_increasing"] is True
    assert checks["status"] == "ok"
    assert checks["max_absorption_frequency_Hz"] == 125.0
    assert checks["min_absorption_frequency_Hz"] == 1000.0
    assert checks["max_force_frequency_Hz"] == 1000.0
    assert checks["min_force_frequency_Hz"] == 125.0
    assert checks["matched_absorption"] == 1.0
    assert checks["double_resistance_reflectance"] == 1.0 / 9.0
    assert checks["reactive_reflectance"] == 1.0
    assert checks["passivity_probe_status"] == "needs_attention"
    assert checks["passivity_probe_excess_absorption"] == 8.0

    return {
        "kind": "acoustic_impedance_sweep_validation",
        "validation_class": True,
        "learning_theme": (
            "acoustic impedance sweeps should report reflection, absorption, "
            "momentum pressure, and passive-load diagnostics together"
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
    print("[acoustic impedance sweep]")
    print(f"  n_points={checks['n_points']} status={checks['status']}")
    print(
        f"  max_absorption_frequency={checks['max_absorption_frequency_Hz']:.12g} Hz "
        f"max_force_frequency={checks['max_force_frequency_Hz']:.12g} Hz"
    )
    print(
        "  passivity_probe="
        f"{checks['passivity_probe_status']} "
        f"excess_absorption={checks['passivity_probe_excess_absorption']:.3e}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
