"""Planar BDM/RT HDiv-MMM response-reduction smoke on an L-shaped body.

The reduced model is compared with independent parent-space solves in magnetic
energy, element-average magnetization, and the re-entrant-corner neighborhood.
This is a desktop correctness smoke, not a timing benchmark.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SRC = REPO_ROOT / "src"
DEFAULT_OUTPUT = HERE / "planar_hdiv_mmm_response_smoke.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import radia.vim as vim


def _make_l_mesh(maxh: float):
    import ngsolve as ng
    from netgen.geom2d import SplineGeometry

    geometry = SplineGeometry()
    coordinates = ((0, 0), (1, 0), (1, 0.4), (0.4, 0.4), (0.4, 1), (0, 1))
    points = [geometry.AppendPoint(*point) for point in coordinates]
    for index in range(len(points)):
        geometry.Append(
            ["line", points[index], points[(index + 1) % len(points)]],
            leftdomain=1,
            rightdomain=0,
            bc="surface",
        )
    geometry.SetMaterial(1, "body")
    return ng.Mesh(geometry.GenerateMesh(maxh=maxh))


def _element_centroids(mesh) -> np.ndarray:
    import ngsolve as ng

    centroids = np.empty((mesh.ne, 2), dtype=float)
    for element in mesh.Elements(ng.VOL):
        points = np.asarray(
            [tuple(mesh[vertex].point)[:2] for vertex in element.vertices],
            dtype=float,
        )
        centroids[element.nr] = points.mean(axis=0)
    return centroids


def _weighted_relative_l2(
    difference: np.ndarray,
    reference: np.ndarray,
    weights: np.ndarray,
) -> float:
    numerator = np.sum(weights * np.sum(np.abs(difference) ** 2, axis=1))
    denominator = np.sum(weights * np.sum(np.abs(reference) ** 2, axis=1))
    return float(np.sqrt(numerator / max(denominator, np.finfo(float).tiny)))


def run(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    started = time.perf_counter()
    mesh = _make_l_mesh(args.maxh)
    fields = (
        ng.CoefficientFunction((1.0, 0.0)),
        ng.CoefficientFunction((0.0, 1.0)),
    )
    ports = vim.NgsolvePlanarHarmonicPorts(
        mesh,
        max_degree=args.harmonic_degree,
    )
    with ng.TaskManager():
        model = vim.NgsolvePlanarHDivMMMResponseReduction(
            mesh,
            mu_r=args.mu_r,
            order=args.order,
            rt=args.family == "rt",
            external_fields=fields,
            external_names=("H_x", "H_y"),
            training_fields=ports,
            pod_rtol=args.pod_rtol,
            cg_tol=args.cg_tol,
        )
        reduced = model.solve()
        projected = [model.body.project(field) for field in fields]
        parent = np.column_stack(
            [model.body.solve_linear(args.mu_r - 1.0, values) for values in projected]
        )

    inv_chi = 1.0 / (args.mu_r - 1.0)

    def apply_operator(values: np.ndarray) -> np.ndarray:
        return np.asarray(model.body.apply_demag(values)) + inv_chi * np.asarray(
            model.body.Mm @ values
        ).reshape(-1)

    centroids = _element_centroids(mesh)
    corner_distance = np.linalg.norm(centroids - np.array((0.4, 0.4)), axis=1)
    corner_mask = corner_distance <= args.corner_radius
    if not np.any(corner_mask):
        corner_mask[np.argmin(corner_distance)] = True

    rows = []
    for index, name in enumerate(("H_x", "H_y")):
        parent_values = parent[:, index]
        reconstructed = reduced.parent_coefficients[:, index]
        difference = reconstructed - parent_values
        energy_error = float(
            np.sqrt(
                max(float(difference @ apply_operator(difference)), 0.0)
                / max(
                    float(parent_values @ apply_operator(parent_values)),
                    np.finfo(float).tiny,
                )
            )
        )
        parent_m = model.body.M_elem(parent_values)
        reduced_m = model.body.M_elem(reconstructed)
        rows.append(
            {
                "field": name,
                "parent_relative_energy_error": energy_error,
                "element_magnetization_relative_l2_error": _weighted_relative_l2(
                    reduced_m - parent_m,
                    parent_m,
                    model.body.areas,
                ),
                "corner_magnetization_relative_l2_error": _weighted_relative_l2(
                    (reduced_m - parent_m)[corner_mask],
                    parent_m[corner_mask],
                    model.body.areas[corner_mask],
                ),
                "parent_average_magnetization": list(model.body.M_avg(parent_values)),
                "reduced_average_magnetization": list(
                    model.body.M_avg(reconstructed)
                ),
            }
        )

    generation = model.basis_generation
    checks = {
        "rt_order_supported": args.order in (1, 2),
        "harmonic_order_admissible": args.order >= args.harmonic_degree - 1,
        "physical_responses_protected": generation["protected_physical_modes"] == 2,
        "snapshot_solves_converged": (
            generation["max_snapshot_relative_residual"] < 1.0e-8
        ),
        "parent_energy_errors_below_1e-8": all(
            row["parent_relative_energy_error"] < 1.0e-8 for row in rows
        ),
        "corner_errors_below_1e-8": all(
            row["corner_magnetization_relative_l2_error"] < 1.0e-8
            for row in rows
        ),
        "reduced_residual_below_1e-12": reduced.residual_relative_norm < 1.0e-12,
    }
    checks["passed"] = all(checks.values())
    return {
        "schema": "radia.validation.planar_hdiv_mmm_response.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": "Desktop correctness smoke; elapsed time is not benchmark data.",
        "configuration": {
            "geometry": "L-shaped re-entrant corner",
            "hdiv_order": args.order,
            "hdiv_family": args.family,
            "harmonic_degree": args.harmonic_degree,
            "maxh": args.maxh,
            "mu_r": args.mu_r,
            "pod_rtol": args.pod_rtol,
            "cg_tol": args.cg_tol,
            "corner_radius": args.corner_radius,
        },
        "mesh": {
            "elements": int(mesh.ne),
            "corner_elements": int(np.count_nonzero(corner_mask)),
        },
        "reduction": model.diagnostics(),
        "solution": reduced.diagnostics(),
        "field_rows": rows,
        "checks": checks,
        "elapsed_seconds_desktop_smoke": float(time.perf_counter() - started),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--family", choices=("bdm", "rt"), default="bdm")
    parser.add_argument("--harmonic-degree", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--maxh", type=float, default=0.35)
    parser.add_argument("--mu-r", type=float, default=1001.0)
    parser.add_argument("--pod-rtol", type=float, default=1.0e-10)
    parser.add_argument("--cg-tol", type=float, default=1.0e-10)
    parser.add_argument("--corner-radius", type=float, default=0.35)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if args.order < args.harmonic_degree - 1:
        parser.error("--order must be at least --harmonic-degree - 1")
    if args.mu_r <= 1.0 or args.pod_rtol <= 0.0 or args.cg_tol <= 0.0:
        parser.error("material and solver tolerances must be positive")
    if args.maxh <= 0.0 or args.corner_radius <= 0.0:
        parser.error("mesh size and corner radius must be positive")

    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    reduction = result["reduction"]
    print("Planar HDiv-MMM response-reduction smoke")
    print(f"  output: {args.output}")
    print(
        f"  {args.family.upper()}{args.order}: {reduction['parent_ndof']} DoF -> "
        f"{reduction['reduced_modes']} response modes"
    )
    for row in result["field_rows"]:
        print(
            f"  {row['field']}: energy={row['parent_relative_energy_error']:.3e}  "
            f"corner={row['corner_magnetization_relative_l2_error']:.3e}"
        )
    print(f"  passed: {result['checks']['passed']}")
    return 0 if result["checks"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
