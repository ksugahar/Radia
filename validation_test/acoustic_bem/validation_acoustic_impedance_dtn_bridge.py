"""Validation-class acoustic impedance <-> DtN/Robin bridge.

This example makes the scalar acoustic FEM/BEM sign convention explicit.  A
boundary impedance or admittance row can be written as a Helmholtz DtN/Robin
coefficient, and vice versa, using only Euler's equation.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.acoustics import (  # noqa: E402
    acoustic_dtn_from_impedance,
    acoustic_impedance_from_dtn,
    baffled_circular_piston_radiation,
    planar_mode_radiation_impedance,
    spherical_mode_radiation_impedance,
)


OUT_JSON = HERE / "validation_acoustic_impedance_dtn_bridge_summary.json"


def _complex_record(value):
    z = complex(value)
    return {"real": z.real, "imag": z.imag, "abs": abs(z)}


def _row(label, frequency, impedance, reference_dtn, rho):
    dtn = acoustic_dtn_from_impedance(frequency, specific_impedance=impedance, rho=rho)
    inv = acoustic_impedance_from_dtn(frequency, dtn["dtn_eigenvalue"], rho=rho)
    return {
        "label": label,
        "frequency": frequency,
        "specific_impedance": _complex_record(impedance),
        "dtn_eigenvalue": _complex_record(dtn["dtn_eigenvalue"]),
        "reference_dtn_eigenvalue": _complex_record(reference_dtn),
        "dtn_abs_error": abs(dtn["dtn_eigenvalue"] - reference_dtn),
        "roundtrip_impedance_abs_error": abs(inv["specific_impedance"] - impedance),
    }


def main() -> int:
    rho = 1.2041
    c = 343.0
    frequency = 1000.0
    k = 2.0 * math.pi * frequency / c

    normal_impedance = rho * c
    normal = _row("planar_normal", frequency, normal_impedance, -1j * k, rho)

    theta = math.radians(55.0)
    oblique_mode = planar_mode_radiation_impedance(
        frequency, incidence_angle_rad=theta, rho=rho, c=c
    )
    oblique = _row(
        "planar_oblique_55deg",
        frequency,
        oblique_mode["specific_impedance"],
        oblique_mode["dtn_eigenvalue"],
        rho,
    )

    sphere_mode = spherical_mode_radiation_impedance(0.18, frequency, 1, rho=rho, c=c)
    sphere = _row(
        "spherical_degree_1",
        frequency,
        sphere_mode["specific_impedance"],
        sphere_mode["dtn_eigenvalue"],
        rho,
    )

    piston = baffled_circular_piston_radiation(0.075, frequency, 0.01, rho=rho, c=c)
    piston_dtn = acoustic_dtn_from_impedance(
        frequency, specific_impedance=piston["specific_impedance"], rho=rho
    )
    piston_roundtrip = acoustic_impedance_from_dtn(
        frequency, piston_dtn["dtn_eigenvalue"], rho=rho
    )
    piston_row = {
        "label": "baffled_piston_average_impedance",
        "frequency": frequency,
        "ka": piston["ka"],
        "specific_impedance": _complex_record(piston["specific_impedance"]),
        "dtn_eigenvalue": _complex_record(piston_dtn["dtn_eigenvalue"]),
        "roundtrip_impedance_abs_error": abs(
            piston_roundtrip["specific_impedance"] - piston["specific_impedance"]
        ),
        "radiated_active_power_W": piston["radiated_power"],
    }

    rows = [normal, oblique, sphere]
    checks = {
        "normal_dtn_abs_error": normal["dtn_abs_error"],
        "max_reference_dtn_abs_error": max(row["dtn_abs_error"] for row in rows),
        "max_roundtrip_impedance_abs_error": max(
            [row["roundtrip_impedance_abs_error"] for row in rows]
            + [piston_row["roundtrip_impedance_abs_error"]]
        ),
        "piston_radiated_active_power_W": piston_row["radiated_active_power_W"],
    }
    assert checks["normal_dtn_abs_error"] < 1.0e-14
    assert checks["max_reference_dtn_abs_error"] < 1.0e-12
    assert checks["max_roundtrip_impedance_abs_error"] < 1.0e-12
    assert checks["piston_radiated_active_power_W"] > 0.0

    summary = {
        "kind": "acoustic_impedance_dtn_bridge_validation",
        "validation_class": True,
        "rho": rho,
        "c": c,
        "frequency": frequency,
        "rows": rows,
        "piston_row": piston_row,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[acoustic impedance <-> DtN bridge]")
    for row in rows:
        print(
            f"  {row['label']}: dtn={row['dtn_eigenvalue']['real']:+.6e}"
            f"{row['dtn_eigenvalue']['imag']:+.6e}j "
            f"err={row['dtn_abs_error']:.3e}"
        )
    print(
        f"  piston: ka={piston_row['ka']:.6f}, "
        f"roundtrip={piston_row['roundtrip_impedance_abs_error']:.3e}, "
        f"P={piston_row['radiated_active_power_W']:.6e} W"
    )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
