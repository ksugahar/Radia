"""Validation-class TEM transmission-line geometry sweep.

This example is a quasi-static RF preflight for canonical transmission lines:

* coax and two-wire lines satisfy the geometry-independent TEM identity L*C;
* coax/two-wire characteristic impedance grows with conductor spacing;
* wire-over-ground capacitance matches the image-symmetry two-wire result;
* microstrip effective permittivity stays between air and substrate limits.

Run:

    python examples/rf_waveguide/validation_tem_line_geometry_sweep.py
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

from radia_mcp.radia_ngsolve.transmission_line import (  # noqa: E402
    C0,
    coaxial_line_parameters,
    microstrip_line_parameters,
    tem_lc_identity_summary,
    two_wire_line_parameters,
    wire_capacitance_per_length,
)


OUT_JSON = HERE / "validation_tem_line_geometry_sweep_summary.json"

COAX_INNER_RADIUS_M = 0.5e-3
COAX_OUTER_RADII_M = (1.5e-3, 2.5e-3, 5.0e-3, 10.0e-3)
COAX_EPS_R = 2.1
TWO_WIRE_RADIUS_M = 0.5e-3
TWO_WIRE_SEPARATIONS_M = (3.0e-3, 6.0e-3, 12.0e-3, 24.0e-3)
MICROSTRIP_H_M = 0.8e-3
MICROSTRIP_EPS_R = 4.2
MICROSTRIP_WIDTHS_M = (0.24e-3, 0.8e-3, 1.6e-3, 4.0e-3, 8.0e-3)


def _identity_for(row: dict, eps_r: float, mu_r: float = 1.0) -> dict:
    return tem_lc_identity_summary(row["C_per_m"], row["L_per_m"], eps_r=eps_r, mu_r=mu_r)


def build_coax_rows() -> list[dict]:
    rows = []
    for outer in COAX_OUTER_RADII_M:
        params = coaxial_line_parameters(COAX_INNER_RADIUS_M, outer, eps_r=COAX_EPS_R)
        ident = _identity_for(params, eps_r=COAX_EPS_R)
        rows.append({
            "inner_radius_m": COAX_INNER_RADIUS_M,
            "outer_radius_m": outer,
            "radius_ratio": outer / COAX_INNER_RADIUS_M,
            **params,
            "LC_relative_error": ident["LC_relative_error"],
            "vp_relative_error": ident["vp_relative_error"],
        })
    return rows


def build_two_wire_rows() -> list[dict]:
    rows = []
    for separation in TWO_WIRE_SEPARATIONS_M:
        params = two_wire_line_parameters(TWO_WIRE_RADIUS_M, separation)
        ident = _identity_for(params, eps_r=1.0)
        rows.append({
            "wire_radius_m": TWO_WIRE_RADIUS_M,
            "separation_m": separation,
            "separation_over_diameter": separation / (2.0 * TWO_WIRE_RADIUS_M),
            **params,
            "LC_relative_error": ident["LC_relative_error"],
            "vp_relative_error": ident["vp_relative_error"],
        })
    return rows


def build_microstrip_rows() -> list[dict]:
    rows = []
    for width in MICROSTRIP_WIDTHS_M:
        params = microstrip_line_parameters(width, MICROSTRIP_H_M, MICROSTRIP_EPS_R)
        rows.append({
            "width_m": width,
            "substrate_h_m": MICROSTRIP_H_M,
            "eps_r": MICROSTRIP_EPS_R,
            **params,
            "vp_over_c0": params["vp"] / C0,
        })
    return rows


def build_image_symmetry_rows() -> list[dict]:
    rows = []
    for height in (2.0e-3, 5.0e-3, 10.0e-3):
        plane = wire_capacitance_per_length("wire_plane", TWO_WIRE_RADIUS_M, height)["C_per_m"]
        two = wire_capacitance_per_length("two_wire", TWO_WIRE_RADIUS_M, 2.0 * height)["C_per_m"]
        rows.append({
            "wire_radius_m": TWO_WIRE_RADIUS_M,
            "height_m": height,
            "wire_plane_C_per_m": plane,
            "two_wire_D_2h_C_per_m": two,
            "image_symmetry_abs_error": abs(plane - 2.0 * two),
            "image_symmetry_rel_error": abs(plane - 2.0 * two) / plane,
        })
    return rows


def _strictly_increasing(values: list[float]) -> bool:
    return all(a < b for a, b in zip(values, values[1:]))


def _strictly_decreasing(values: list[float]) -> bool:
    return all(a > b for a, b in zip(values, values[1:]))


def validate(coax_rows: list[dict], two_wire_rows: list[dict],
             microstrip_rows: list[dict], image_rows: list[dict]) -> dict:
    checks = {
        "coax_Z0_monotone_increasing": _strictly_increasing([row["Z0"] for row in coax_rows]),
        "coax_C_monotone_decreasing": _strictly_decreasing([row["C_per_m"] for row in coax_rows]),
        "coax_L_monotone_increasing": _strictly_increasing([row["L_per_m"] for row in coax_rows]),
        "coax_max_abs_LC_relative_error": max(abs(row["LC_relative_error"]) for row in coax_rows),
        "coax_vp_constant": max(row["vp"] for row in coax_rows) - min(row["vp"] for row in coax_rows),
        "two_wire_Z0_monotone_increasing": _strictly_increasing([row["Z0"] for row in two_wire_rows]),
        "two_wire_C_monotone_decreasing": _strictly_decreasing([row["C_per_m"] for row in two_wire_rows]),
        "two_wire_max_abs_LC_relative_error": max(abs(row["LC_relative_error"]) for row in two_wire_rows),
        "microstrip_Z0_monotone_decreasing": _strictly_decreasing([row["Z0"] for row in microstrip_rows]),
        "microstrip_eps_eff_monotone_increasing": _strictly_increasing([row["eps_eff"] for row in microstrip_rows]),
        "microstrip_eps_eff_in_bounds": all(
            1.0 < row["eps_eff"] < MICROSTRIP_EPS_R for row in microstrip_rows
        ),
        "image_symmetry_max_rel_error": max(row["image_symmetry_rel_error"] for row in image_rows),
        "coax_first_Z0_ohm": coax_rows[0]["Z0"],
        "coax_last_Z0_ohm": coax_rows[-1]["Z0"],
        "two_wire_first_Z0_ohm": two_wire_rows[0]["Z0"],
        "two_wire_last_Z0_ohm": two_wire_rows[-1]["Z0"],
        "microstrip_narrow_Z0_ohm": microstrip_rows[0]["Z0"],
        "microstrip_wide_Z0_ohm": microstrip_rows[-1]["Z0"],
    }

    assert checks["coax_Z0_monotone_increasing"]
    assert checks["coax_C_monotone_decreasing"]
    assert checks["coax_L_monotone_increasing"]
    assert checks["coax_max_abs_LC_relative_error"] < 1.0e-12
    assert checks["coax_vp_constant"] < 1.0e-6
    assert checks["two_wire_Z0_monotone_increasing"]
    assert checks["two_wire_C_monotone_decreasing"]
    assert checks["two_wire_max_abs_LC_relative_error"] < 1.0e-12
    assert checks["microstrip_Z0_monotone_decreasing"]
    assert checks["microstrip_eps_eff_monotone_increasing"]
    assert checks["microstrip_eps_eff_in_bounds"]
    assert checks["image_symmetry_max_rel_error"] < 1.0e-12
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    coax_rows = build_coax_rows()
    two_wire_rows = build_two_wire_rows()
    microstrip_rows = build_microstrip_rows()
    image_rows = build_image_symmetry_rows()
    checks = validate(coax_rows, two_wire_rows, microstrip_rows, image_rows)

    summary = {
        "kind": "tem_transmission_line_geometry_sweep_validation",
        "validation_class": True,
        "coax_eps_r": COAX_EPS_R,
        "coax_rows": coax_rows,
        "two_wire_rows": two_wire_rows,
        "microstrip_rows": microstrip_rows,
        "image_symmetry_rows": image_rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[coax TEM geometry sweep]")
    for row in coax_rows:
        print(
            f"  b/a={row['radius_ratio']:5.1f}  "
            f"Z0={row['Z0']:8.3f} ohm  "
            f"C={row['C_per_m'] * 1e12:8.3f} pF/m  "
            f"L={row['L_per_m'] * 1e9:8.3f} nH/m"
        )
    print("[two-wire TEM geometry sweep]")
    for row in two_wire_rows:
        print(
            f"  D/2a={row['separation_over_diameter']:5.1f}  "
            f"Z0={row['Z0']:8.3f} ohm  "
            f"C={row['C_per_m'] * 1e12:8.3f} pF/m"
        )
    print("[microstrip sweep]")
    for row in microstrip_rows:
        print(
            f"  w/h={row['u']:5.2f}  "
            f"eps_eff={row['eps_eff']:7.4f}  "
            f"Z0={row['Z0']:8.3f} ohm  "
            f"vp/c0={row['vp_over_c0']:7.4f}"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
