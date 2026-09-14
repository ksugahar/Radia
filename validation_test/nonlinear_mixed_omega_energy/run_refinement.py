"""Run a three-level nonlinear mixed-Omega field/energy refinement study."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import ngsolve as ng
import numpy as np
from netgen.occ import Box, Glue, OCCGeometry, Pnt, X

from radia.kelvin_solver import (
    project_source_interface_potential,
    solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin,
)
from radia_mcp.radia_ngsolve.field_profile_gate import (
    nonlinear_magnetic_refinement_energy_gate,
)


MU0 = 4.0e-7 * math.pi
MESH_SIZES_M = (0.50, 0.32, 0.25)
BH_TABLE = (
    (0.0, 0.0),
    (0.5, 0.5 * MU0 * 2000.0),
    (2.0, 1.4e-3),
    (8.0, 2.0e-3),
    (40.0, 2.4e-3),
)
REFINEMENT_PARENT_IDENTITY = {
    "schema": "radia.nonlinear-magnetic-refinement-parent.v1",
    "geometry": "two-unit-boxes-sharing-x0-interface",
    "materials": {
        "reduced": "vacuum",
        "total": "single-valued-isotropic-soft-magnetic",
    },
    "excitation": "fixed-analytic-source-field-and-interface-trace",
    "coordinate_system": "right-handed Cartesian",
    "unit_system": "SI",
}


def build_case(maxh: float):
    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "outer"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "outer"
    total.faces.Min(X).name = "source_total_interface"
    mesh = ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=maxh))

    r_x = ng.x - 2.5
    r2 = r_x * r_x + ng.y * ng.y + ng.z * ng.z
    exact = ng.CoefficientFunction(
        (r_x / r2**1.5, ng.y / r2**1.5, ng.z / r2**1.5)
    )
    a = (1.0 - ng.x**2) ** 2
    b = (1.0 - ng.y**2) ** 2
    c = (1.0 - ng.z**2) ** 2
    da_dx = -4.0 * ng.x * (1.0 - ng.x**2)
    db_dy = -4.0 * ng.y * (1.0 - ng.y**2)
    source_h = exact + ng.CoefficientFunction(
        (a * db_dy * c, -da_dx * b * c, 0.0)
    )
    trace = project_source_interface_potential(
        mesh,
        source_h,
        "source_total_interface",
        order=2,
        relative_tolerance=0.04,
    )
    return mesh, source_h, trace


def solve_level(maxh: float) -> dict:
    started = time.perf_counter()
    with ng.TaskManager():
        mesh, source_h, trace = build_case(maxh)
        solved = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh,
            source_h,
            trace["potential"],
            1.0,
            (3.0, 0.0, 0.0),
            bh_table=BH_TABLE,
            nonlinear_materials=("total",),
            reduced_materials=("reduced",),
            total_materials=("total",),
            interface_boundary="source_total_interface",
            order=2,
            material_update_order=1,
            dirichlet_bbbnd="outer",
            tolerance=1.0e-6,
            max_iterations=60,
            relaxation=0.3,
            anderson_depth=0,
            refinement_parent_identity=REFINEMENT_PARENT_IDENTITY,
        )
        selector = mesh.Materials("total")
        volume = float(ng.Integrate(1.0, mesh, definedon=selector, order=8))
        average = [
            float(ng.Integrate(solved["B_cf"][index], mesh, definedon=selector, order=8))
            / volume
            for index in range(3)
        ]
        rms = math.sqrt(
            float(
                ng.Integrate(
                    ng.InnerProduct(solved["B_cf"], solved["B_cf"]),
                    mesh,
                    definedon=selector,
                    order=8,
                )
            )
            / volume
        )
    stats = solved["nonlinear_stats"]
    energy = solved["energy_observables"]
    return {
        "mesh_size_m": maxh,
        "mesh_vertices": int(mesh.nv),
        "mesh_volume_elements": int(sum(1 for _ in mesh.Elements(ng.VOL))),
        "solver_converged": bool(stats["converged"]),
        "iterations": int(stats["iterations"]),
        "response_order": int(stats["response_order"]),
        "material_update_order": int(stats["material_update_order"]),
        "physical_relative_permeability_bounds": stats[
            "physical_permeability_bounds"
        ],
        "volume_m3": volume,
        "average_field_T": average,
        "rms_magnitude_T": rms,
        "magnetic_energy_J": float(energy["energy_J"]),
        "magnetic_coenergy_J": float(energy["coenergy_J"]),
        "h_dot_b_integral_J": float(energy["h_dot_b_integral_J"]),
        "legendre_residual_relative": float(
            energy["legendre_residual_relative"]
        ),
        "field_identity": solved["field_observable_identity"],
        "energy_identity": energy["identity"],
        "elapsed_seconds": time.perf_counter() - started,
        "source_trace_relative_residual": float(
            trace["relative_tangential_residual"]
        ),
    }


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_hdf5(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as artifact:
        artifact.attrs.update(
            {
                "schema": payload["schema"],
                "artifact_id": payload["artifact_id"],
                "created_at_utc": payload["created_at_utc"],
                "producer": "validation_test/nonlinear_mixed_omega_energy",
                "producer_version": payload["source_commit"],
                "solver": "radia-ngsolve mixed total/reduced Omega",
                "coordinate_system": "right-handed Cartesian",
                "unit_system": "SI",
                "source_commit": payload["source_commit"],
                "total_compute_seconds": payload["total_compute_seconds"],
                "timing_breakdown_s": json.dumps(
                    {
                        f"level_{index}": level["elapsed_seconds"]
                        for index, level in enumerate(payload["levels"])
                    },
                    sort_keys=True,
                ),
                "gate_status": payload["gate"]["status"],
                "dimension_order": "C",
            }
        )
        levels = artifact.create_group("levels")
        levels.attrs["axes"] = "refinement_level"
        for key, unit in (
            ("mesh_size_m", "m"),
            ("volume_m3", "m^3"),
            ("average_field_T", "T"),
            ("rms_magnitude_T", "T"),
            ("magnetic_energy_J", "J"),
            ("magnetic_coenergy_J", "J"),
            ("h_dot_b_integral_J", "J"),
            ("legendre_residual_relative", "1"),
            ("elapsed_seconds", "s"),
        ):
            dataset = levels.create_dataset(
                key, data=np.asarray([level[key] for level in payload["levels"]])
            )
            dataset.attrs.update(
                {
                    "unit": unit,
                    "location": "global",
                    "axes": "refinement_level"
                    + (",component" if key == "average_field_T" else ""),
                    "dimension_order": "C",
                }
            )
        identity_dataset = levels.create_dataset(
            "identity_json",
            data=np.asarray(
                [
                    json.dumps(level["energy_identity"], sort_keys=True)
                    for level in payload["levels"]
                ],
                dtype=h5py.string_dtype("utf-8"),
            ),
        )
        identity_dataset.attrs.update(
            {
                "unit": "1",
                "location": "global",
                "axes": "refinement_level",
                "dimension_order": "C",
            }
        )
        gate_dataset = artifact.create_dataset(
            "gate_json",
            data=json.dumps(payload["gate"], sort_keys=True),
            dtype=h5py.string_dtype("utf-8"),
        )
        gate_dataset.attrs.update(
            {
                "unit": "1",
                "location": "global",
                "axes": "scalar",
                "dimension_order": "C",
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    levels = [solve_level(maxh) for maxh in MESH_SIZES_M]
    gate = nonlinear_magnetic_refinement_energy_gate({"levels": levels})
    repo = Path(__file__).resolve().parents[2]
    payload = {
        "schema": "cae-ai-lab.nonlinear-mixed-omega-energy-refinement.v2",
        "artifact_id": "nonlinear-mixed-omega-energy-three-level",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True
        ).strip(),
        "levels": levels,
        "gate": gate,
        "total_compute_seconds": time.perf_counter() - started,
        "claim_boundary": (
            "Self-convergence and field/energy identity only; this artifact is not "
            "cross-solver parity evidence."
        ),
    }
    write_hdf5(args.output, payload)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "sha256": file_sha256(args.output),
                "gate_status": gate["status"],
                "issues": gate["issues"],
                "total_compute_seconds": payload["total_compute_seconds"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
