"""Export corner-refined HCurl-FEM/EVRS fields and local basis rank to Gmsh."""

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
DEFAULT_OUTPUT = HERE / "hcurl_corner_fields.msh"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import ngsolve as ng  # noqa: E402
import radia.vim as vim  # noqa: E402
from radia.gmsh_post_export import GmshPostExport  # noqa: E402

from evrs_current_field_compare import (  # noqa: E402
    _assemble_parent_on_mesh,
    _field_metrics,
    _free_mask,
    _project_reduced,
    _refine_reentrant_corner,
    _response_basis,
    _sample_current,
    _solve_full,
    _solve_reduced,
)
from evrs_sibc_mixed_schur import _make_skin_mesh  # noqa: E402


def _grid_function(fes, coefficients: np.ndarray):
    values = np.asarray(coefficients)
    if np.max(np.abs(values.imag)) > 1.0e-11 * max(np.max(np.abs(values.real)), 1.0):
        raise ValueError("visualization currently expects real HCurl coefficients")
    gf = ng.GridFunction(fes)
    gf.vec.FV().NumPy()[:] = values.real
    return gf


def _local_basis_diagnostics(mesh, fes, response, intorder: int, rtol: float):
    mode_curls = []
    for column in np.asarray(response.vectors).T:
        mode_curls.append(ng.curl(_grid_function(fes, column)))

    local_rank = []
    energy_density = []
    element_size = []
    corner_distance = []
    for element in mesh.Elements(ng.VOL):
        trafo = mesh.GetTrafo(element)
        rows = []
        volume = 0.0
        centroid_moment = np.zeros(3)
        for ip in ng.IntegrationRule(element.type, intorder):
            mip = trafo(ip)
            weight = float(ip.weight * mip.measure)
            volume += weight
            centroid_moment += weight * np.asarray(mip.point, dtype=float)
            rows.extend(
                [
                    np.sqrt(weight) * np.asarray(mode(mip), dtype=float)
                    for mode in mode_curls
                ]
            )
        matrix = np.asarray(rows).reshape(-1, len(mode_curls), 3).transpose(0, 2, 1)
        matrix = matrix.reshape(-1, len(mode_curls))
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        threshold = rtol * singular_values[0] if singular_values.size else 0.0
        local_rank.append(int(np.count_nonzero(singular_values > threshold)))
        energy_density.append(float(np.sum(singular_values**2) / volume))
        element_size.append(float(volume ** (1.0 / 3.0)))
        centroid = centroid_moment / volume
        corner_distance.append(float(np.linalg.norm(centroid[:2] - (0.45, 0.45))))
    return {
        "local_effective_rank": np.asarray(local_rank, dtype=float),
        "basis_energy_density": np.asarray(energy_density),
        "element_size": np.asarray(element_size),
        "corner_distance": np.asarray(corner_distance),
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    started = time.perf_counter()
    mesh = _make_skin_mesh(
        args.maxh,
        "notched-box",
        corner_edge_maxh=args.corner_edge_maxh,
    )
    refinement = _refine_reentrant_corner(
        mesh,
        "notched-box",
        args.corner_refinements,
        args.refinement_radius,
    )
    fes, stiffness, mass, ports = _assemble_parent_on_mesh(mesh, args.order, False)
    free_dofs = fes.FreeDofs(False)
    full_coefficients, full_solve_seconds = _solve_full(
        stiffness,
        mass,
        ports,
        free_dofs,
        args.shift,
        args.inverse,
    )
    response = _response_basis(
        fes,
        stiffness,
        mass,
        ports,
        steps=args.steps,
        condense=False,
        inverse=args.inverse,
        rtol=args.rtol,
    )
    kr, mr, br = _project_reduced(stiffness, mass, ports, response)
    reduced_coefficients, _ = _solve_reduced(
        response,
        kr,
        mr,
        br,
        args.shift,
        1,
    )

    port = args.port
    full_gf = _grid_function(fes, full_coefficients[:, port])
    reduced_gf = _grid_function(fes, reduced_coefficients[:, port])
    current_fem = ng.curl(full_gf)
    current_evrs = ng.curl(reduced_gf)
    current_error = current_evrs - current_fem
    magnitude_fem = ng.sqrt(ng.InnerProduct(current_fem, current_fem))
    magnitude_evrs = ng.sqrt(ng.InnerProduct(current_evrs, current_evrs))
    magnitude_error = ng.sqrt(ng.InnerProduct(current_error, current_error))

    local = _local_basis_diagnostics(
        mesh,
        fes,
        response,
        args.activity_intorder,
        args.local_rank_rtol,
    )
    post = GmshPostExport(mesh)
    post.add_vector_field("J_FEM", current_fem, cell_data=True)
    post.add_scalar_field("J_FEM_magnitude", magnitude_fem, cell_data=True)
    post.add_vector_field("J_EVRS", current_evrs, cell_data=True)
    post.add_scalar_field("J_EVRS_magnitude", magnitude_evrs, cell_data=True)
    post.add_scalar_field("J_error_magnitude", magnitude_error, cell_data=True)
    post.add_scalar_field(
        "EVRS_local_effective_rank",
        local["local_effective_rank"],
        cell_data=True,
    )
    post.add_scalar_field(
        "EVRS_basis_energy_density",
        local["basis_energy_density"],
        cell_data=True,
    )
    post.add_scalar_field("element_size", local["element_size"], cell_data=True)
    post.add_scalar_field("corner_distance", local["corner_distance"], cell_data=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output = Path(post.write(str(args.output)))

    reference_current = _sample_current(
        mesh,
        fes,
        full_coefficients,
        intorder=args.metric_intorder,
        prefix="fem",
    )
    reduced_current = _sample_current(
        mesh,
        fes,
        reduced_coefficients,
        intorder=args.metric_intorder,
        prefix="evrs",
    )
    metrics = _field_metrics(
        reduced_current,
        reference_current,
        sigma=args.sigma,
        geometry="notched-box",
        corner_radius=args.corner_radius,
    )
    info = response.diagnostics()
    corner_mask = local["corner_distance"] <= args.corner_radius
    far_mask = ~corner_mask
    corner_energy = float(np.mean(local["basis_energy_density"][corner_mask]))
    far_energy = float(np.mean(local["basis_energy_density"][far_mask]))
    result = {
        "schema": "radia.validation.hcurl_corner_gmsh.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "configuration": vars(args) | {"output": str(output)},
        "mesh": {
            "elements": int(mesh.ne),
            "refinement_history": refinement,
        },
        "basis": info,
        "local_basis": {
            "minimum_effective_rank": int(np.min(local["local_effective_rank"])),
            "maximum_effective_rank": int(np.max(local["local_effective_rank"])),
            "corner_mean_effective_rank": float(
                np.mean(
                    local["local_effective_rank"][
                        local["corner_distance"] <= args.corner_radius
                    ]
                )
            ),
            "far_mean_effective_rank": float(
                np.mean(
                    local["local_effective_rank"][far_mask]
                )
            ),
            "corner_mean_basis_energy_density": corner_energy,
            "far_mean_basis_energy_density": far_energy,
            "corner_to_far_basis_energy_ratio": float(corner_energy / far_energy),
        },
        "full_fem_solve_seconds": float(full_solve_seconds),
        "field_metrics": metrics,
        "total_wall_seconds": float(time.perf_counter() - started),
    }
    summary = output.with_suffix(".json")
    summary.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--steps", type=int, default=22)
    parser.add_argument("--shift", type=float, default=1.0)
    parser.add_argument("--port", type=int, choices=(0, 1), default=0)
    parser.add_argument("--maxh", type=float, default=2.0)
    parser.add_argument("--corner-refinements", type=int, default=0)
    parser.add_argument("--refinement-radius", type=float, default=0.18)
    parser.add_argument("--corner-edge-maxh", type=float, default=0.2)
    parser.add_argument("--corner-radius", type=float, default=0.35)
    parser.add_argument("--metric-intorder", type=int, default=10)
    parser.add_argument("--activity-intorder", type=int, default=4)
    parser.add_argument("--local-rank-rtol", type=float, default=1.0e-3)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--inverse", default="sparsecholesky")
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if args.steps < 1 or args.order < 1:
        parser.error("--order and --steps must be positive")
    if args.shift <= 0.0:
        parser.error("--shift must be positive")
    if args.corner_edge_maxh <= 0.0:
        parser.error("--corner-edge-maxh must be positive")
    result = run(args)
    print("HCurl corner Gmsh visualization")
    print(f"  output: {args.output}")
    print(f"  parent/reduced: {result['basis']['active_dofs']}/{result['basis']['rank']}")
    print(
        "  local rank corner/far: "
        f"{result['local_basis']['corner_mean_effective_rank']:.2f}/"
        f"{result['local_basis']['far_mean_effective_rank']:.2f}"
    )
    print(
        "  basis energy corner/far ratio: "
        f"{result['local_basis']['corner_to_far_basis_energy_ratio']:.3f}"
    )
    print(
        "  current error global/corner: "
        f"{result['field_metrics']['max_relative_current_l2']:.3e}/"
        f"{result['field_metrics']['max_corner_relative_current_l2']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
