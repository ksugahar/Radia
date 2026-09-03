"""TEAM 28 acceptance summary for HCurl Eddy Bubble + CLN.

This validation locks four independent facts:

1. The 3 mm aluminium disk at 50 Hz is volumetric (skin depth about 12.2 mm),
   even though every exterior face is adjacent to air.
2. A real p=6 HCurl disk mesh can be reduced by EVRS plus the conductor-graph
   cycle basis while retaining the loop bridge class.
3. The existing full-FEM/6-stage-CLN TEAM 28 force curve is a stable numerical
   acceptance target for the new route.
4. The epsilon-free affine-tetrahedron HCurl-VIM interaction reproduces the
   physical force at the reference position on three meshes, with an
   independent outer-quadrature check.

The fixed-position 3-D force and the 25-position CLN curve are separate gates.
An end-to-end moving 3-D sweep and curved-tetrahedron moments remain explicit
production gates.
"""

from __future__ import annotations

import argparse
import json
import platform
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SRC = REPO_ROOT / "src"
TEAM28_JSON = HERE / "demos" / "team28" / "team28_cln_sweep_results.json"
TEAM28_HCURL_FORCE_JSON = HERE / "team28_hcurl_vim_force_summary.json"
DEFAULT_OUTPUT = HERE / "team28_hcurl_eddy_bubble_summary.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import radia
import radia.vim as vim
from radia.maglev import PositionForceCurve


def _team28_curves(data):
    position = np.asarray(data["dZ_mm"], dtype=float) * 1.0e-3
    return {
        "full_fem": PositionForceCurve(position, data["fz_full_N"], "full-fem"),
        "cln": PositionForceCurve(position, data["fz_cln_N"], "cln-6-stage"),
        "reference": PositionForceCurve(position, data["fz_lab_N"], "reference"),
    }


def _p6_disk_topology(maxh_m: float, evrs_rank: int):
    import ngsolve as ng
    import netgen.occ as occ

    disk = occ.Cylinder(occ.Pnt(0.0, 0.0, 0.0), occ.Z, 0.065, 0.003)
    disk.mat("Al")
    for face in disk.faces:
        face.name = "disk_air"
    mesh = ng.Mesh(occ.OCCGeometry(disk).GenerateMesh(maxh=maxh_m))
    fes = ng.HCurl(mesh, order=6, nograds=True)
    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="Al")
    policy = vim.NgsolveEddyDofPolicy(mesh, fes, topology)
    graph = topology.conductor_graph()
    plan = policy.reduction_plan(
        evrs_rank=evrs_rank,
        surface_modes=0,
        loop_bridge_modes=graph.cycle_rank,
        bridge_strategy="cycle-basis",
    )
    return {
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "parent_family": "HCurl",
        "parent_order": 6,
        "parent_ndof": int(fes.ndof),
        "topology": topology.diagnostics(),
        "conductor_graph": graph.diagnostics(),
        "dof_policy": policy.diagnostics(),
        "reduction_plan": plan.diagnostics(),
    }


def run(maxh_m: float = 0.025, evrs_rank: int = 6) -> dict[str, object]:
    started = time.perf_counter()
    with TEAM28_JSON.open("r", encoding="utf-8") as stream:
        source = json.load(stream)
    with TEAM28_HCURL_FORCE_JSON.open("r", encoding="utf-8") as stream:
        hcurl_force = json.load(stream)
    curves = _team28_curves(source)
    cln_vs_full = curves["cln"].compare(curves["full_fem"])
    cln_vs_reference = curves["cln"].compare(curves["reference"])

    physical_force = PositionForceCurve(
        curves["cln"].positions_m,
        0.5 * curves["cln"].force_N,
        "cln-physical-time-average",
    )
    equilibrium = physical_force.crossings(-float(source["disk_weight_N"]))
    equilibrium_dz_m = float(equilibrium[0]) if equilibrium.size else None
    equilibrium_abs_height_m = (
        10.8e-3 + equilibrium_dz_m
        if equilibrium_dz_m is not None
        else None
    )

    sibc = vim.EddySIBCApplicability(
        frequency_hz=50.0,
        sigma=3.4e7,
        characteristic_thickness_m=3.0e-3,
    )
    p6 = _p6_disk_topology(maxh_m=maxh_m, evrs_rank=evrs_rank)
    candidate_faces = int(p6["topology"]["sibc_face_count"])
    selected_sibc_faces = candidate_faces if sibc.sibc_applicable else 0
    reduced_modes = int(p6["reduction_plan"]["estimated_reduced_modes"])
    hcurl_force_passed = bool(hcurl_force["hcurl_vim_force_acceptance_complete"])

    checks = {
        "team28_routes_to_volumetric_hcurl": sibc.selected_model == "volumetric",
        "air_adjacency_creates_sibc_candidates": candidate_faces > 0,
        "skin_gate_disables_all_sibc_faces": selected_sibc_faces == 0,
        "disk_conductor_is_connected": p6["conductor_graph"]["component_count"] == 1,
        "cycle_bridge_basis_is_nonempty": p6["conductor_graph"]["cycle_rank"] > 0,
        "p6_reduction_is_below_one_percent": (
            p6["reduction_plan"]["estimated_reduction_ratio"] < 0.01
        ),
        "cln_force_curve_matches_full_fem_below_5uN": (
            cln_vs_full["max_abs_error_N"] < 5.0e-6
        ),
        "cln_force_curve_matches_reference_below_0p5mN": (
            cln_vs_reference["max_abs_error_N"] < 5.0e-4
        ),
        "physical_equilibrium_found": equilibrium_abs_height_m is not None,
        "physical_equilibrium_within_10_percent_of_11p5mm": (
            equilibrium_abs_height_m is not None
            and abs(equilibrium_abs_height_m - 11.5e-3) / 11.5e-3 < 0.10
        ),
        "fixed_position_3d_hcurl_vim_force_passed": hcurl_force_passed,
    }
    return {
        "schema": "radia.team28.hcurl-eddy-bubble-acceptance.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "hostname": socket.gethostname(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "radia_version": getattr(radia, "__version__", "unknown"),
            "ngsolve_version": getattr(__import__("ngsolve"), "__version__", "unknown"),
            "elapsed_seconds_lab_smoke": time.perf_counter() - started,
        },
        "problem": {
            "name": "TEAM Workshop Problem 28",
            "frequency_hz": 50.0,
            "disk_radius_m": 0.065,
            "disk_thickness_m": 0.003,
            "conductivity_S_per_m": 3.4e7,
            "parent_hcurl_order": 6,
            "evrs_rank": int(evrs_rank),
            "force_convention": "stored force is 2x physical time-average",
        },
        "surface_model": {
            **sibc.diagnostics(),
            "adjacency_candidate_face_count": candidate_faces,
            "selected_sibc_face_count": selected_sibc_faces,
            "volumetric_boundary_face_count": candidate_faces - selected_sibc_faces,
        },
        "p6_spatial_reduction": p6,
        "cln_reference_acceptance": {
            "cln_stages": int(source["cln_stages"]),
            "cln_vs_full_fem": cln_vs_full,
            "cln_vs_reference": cln_vs_reference,
            "equilibrium_dz_m": equilibrium_dz_m,
            "equilibrium_abs_height_m": equilibrium_abs_height_m,
            "published_steady_height_m": 11.5e-3,
        },
        "reduced_mode_summary": {
            "parent_ndof": int(p6["parent_ndof"]),
            "evrs_modes": int(evrs_rank),
            "cycle_bridge_modes": int(p6["conductor_graph"]["cycle_rank"]),
            "sibc_modes": 0,
            "estimated_total_modes": reduced_modes,
            "estimated_reduction_ratio": float(
                p6["reduction_plan"]["estimated_reduction_ratio"]
            ),
        },
        "hcurl_vim_force_acceptance": {
            "result_file": TEAM28_HCURL_FORCE_JSON.name,
            "generated_at_utc": hcurl_force["generated_at_utc"],
            "validation_host": hcurl_force["runtime"]["hostname"],
            "mesh_case_count": len(hcurl_force["cases"]),
            "maximum_force_relative_error": hcurl_force[
                "maximum_force_relative_error"
            ],
            "outer_quadrature_relative_force_change": hcurl_force[
                "outer_quadrature_relative_force_change"
            ],
            "checks": hcurl_force["checks"],
        },
        "checks": checks,
        "structural_and_reference_acceptance_passed": all(checks.values()),
        "hcurl_vim_force_acceptance_complete": hcurl_force_passed,
        "next_numeric_gate": (
            "Run the same HCurl-VIM basis through the 25-position moving CLN "
            "sweep, then keep the reduced HCurl operator in HACApK form instead "
            "of materializing its final dense matrix"
        ),
    }


def main() -> None:
    import ngsolve as ng

    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, default=0.025)
    parser.add_argument("--evrs-rank", type=int, default=6)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    with ng.TaskManager():
        result = run(maxh_m=args.maxh, evrs_rank=args.evrs_rank)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result["structural_and_reference_acceptance_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
