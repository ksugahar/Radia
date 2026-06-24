"""Validation-class two-port scattering momentum force.

Run:

    python examples/rf_waveguide/validation_two_port_scattering_momentum_force.py

The helper is a closed-form RF momentum-balance gate.  It converts incident
power, reflected/transmitted power fractions, and port propagation directions
into the force on the scattering object:

    F = P/c * ((1 + R) k_inc - T k_out)

This keeps the common straight-through, short, absorber, and 90-degree bend
checks readable before a full time-harmonic Maxwell-stress integration is used.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    C0,
    two_port_scattering_momentum_force_summary,
)

OUT_JSON = Path(__file__).with_name("validation_two_port_scattering_momentum_force_summary.json")


def _assert_close(actual: float, expected: float, name: str, tol: float = 1.0e-18) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{name}: got {actual:.16e}, expected {expected:.16e}")


def _assert_vector_close(actual, expected, name: str, tol: float = 1.0e-18) -> None:
    for i, (a, e) in enumerate(zip(actual, expected)):
        _assert_close(float(a), float(e), f"{name}[{i}]", tol=tol)


def main() -> None:
    power = 3.0
    kin = (1.0, 0.0, 0.0)
    ky = (0.0, 1.0, 0.0)
    p_over_c = power / C0

    cases = {
        "transparent_straight": two_port_scattering_momentum_force_summary(
            power,
            kin,
            reflectance=0.0,
            transmittance=1.0,
            transmitted_direction=kin,
        ),
        "perfect_short": two_port_scattering_momentum_force_summary(
            power,
            kin,
            reflectance=1.0,
            transmittance=0.0,
        ),
        "matched_absorber": two_port_scattering_momentum_force_summary(
            power,
            kin,
            reflectance=0.0,
            transmittance=0.0,
        ),
        "lossless_90deg_bend": two_port_scattering_momentum_force_summary(
            power,
            kin,
            reflectance=0.0,
            transmittance=1.0,
            transmitted_direction=ky,
        ),
        "lossy_straight": two_port_scattering_momentum_force_summary(
            power,
            kin,
            reflectance=0.1,
            transmittance=0.6,
            transmitted_direction=kin,
        ),
    }

    _assert_vector_close(cases["transparent_straight"]["force_N"], (0.0, 0.0, 0.0), "straight")
    _assert_vector_close(cases["perfect_short"]["force_N"], (2.0 * p_over_c, 0.0, 0.0), "short")
    _assert_vector_close(cases["matched_absorber"]["force_N"], (p_over_c, 0.0, 0.0), "absorber")
    _assert_vector_close(cases["lossless_90deg_bend"]["force_N"], (p_over_c, -p_over_c, 0.0), "bend")
    _assert_close(
        cases["lossless_90deg_bend"]["force_magnitude_N"],
        math.sqrt(2.0) * p_over_c,
        "bend magnitude",
    )
    _assert_vector_close(cases["lossy_straight"]["force_N"], (0.5 * p_over_c, 0.0, 0.0), "lossy")

    summary = {
        "kind": "two_port_scattering_momentum_force_validation",
        "power_incident_W": power,
        "p_over_c_N": p_over_c,
        "cases": cases,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] two-port scattering momentum force validated; wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
