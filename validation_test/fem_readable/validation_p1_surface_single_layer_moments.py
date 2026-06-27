"""Validation-class P1 surface single-layer moment gate.

Run:

    python validation_test/fem_readable/validation_p1_surface_single_layer_moments.py

For a P1 surface density on one triangle, the total source and first moment
are enough to check the far-field of a Laplace single-layer operator:

    Q = int_T sigma dS
    m = int_T sigma x dS
    phi_far(x) = Q/(4*pi*r) + (xhat.m)/(4*pi*r^2) + O(r^-3)

This is a tiny readable FEM/BEM coupling block.  The same surface mass matrix
formula can be mirrored directly in MATLAB or Gypsilab-style teaching code.
"""

from __future__ import annotations

import json
from result_metadata import add_result_metadata
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    laplace_single_layer_far_potential,
    p1_surface_triangle_density_moments,
    p1_surface_triangle_geometry,
)

OUT_JSON = Path(__file__).with_name("validation_p1_surface_single_layer_moments_summary.json")

TRIANGLE = [(1.0, 0.0, 0.0), (1.0, 2.0, 0.0), (1.0, 0.0, 3.0)]
DENSITY = [2.0, 2.0, 2.0]
OBSERVATION = (10.0, 0.0, 0.0)


def _max_abs(values) -> float:
    return max(abs(float(value)) for value in values)


def main() -> None:
    geom = p1_surface_triangle_geometry(TRIANGLE)
    moments = p1_surface_triangle_density_moments(TRIANGLE, DENSITY)
    far = laplace_single_layer_far_potential(
        OBSERVATION,
        moments["total_source"],
        moments["first_moment"],
    )

    area = geom["area"]
    centroid = tuple(sum(point[i] for point in TRIANGLE) / 3.0 for i in range(3))
    expected_q = 2.0 * area
    expected_m = tuple(expected_q * value for value in centroid)
    expected_monopole = expected_q / (4.0 * math.pi * far["radius"])
    expected_dipole = expected_m[0] / (4.0 * math.pi * far["radius"] ** 2)

    q_error = abs(moments["total_source"] - expected_q)
    m_error = _max_abs(moments["first_moment"][i] - expected_m[i] for i in range(3))
    far_error = abs(far["far_potential"] - (expected_monopole + expected_dipole))
    if q_error > 1.0e-15 or m_error > 1.0e-15 or far_error > 1.0e-15:
        raise AssertionError("P1 surface single-layer moment gate failed")

    summary = {
        "kind": "p1_surface_single_layer_moments",
        "triangle_vertices_m": TRIANGLE,
        "nodal_density": DENSITY,
        "observation_point_m": OBSERVATION,
        "area_m2": area,
        "centroid_m": centroid,
        "moments": moments,
        "far_potential": far,
        "checks": {
            "expected_total_source": expected_q,
            "expected_first_moment": expected_m,
            "expected_monopole_potential": expected_monopole,
            "expected_dipole_potential": expected_dipole,
            "total_source_abs_error": q_error,
            "first_moment_max_abs_error": m_error,
            "far_potential_abs_error": far_error,
        },
    }
    OUT_JSON.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2) + "\n", encoding="utf-8")
    print("[P1 surface single-layer moments]")
    print(f"  area={area:.12g} Q={moments['total_source']:.12g}")
    print(f"  first_moment={moments['first_moment']}")
    print(f"  far_potential={far['far_potential']:.12e}")
    print(f"  wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
