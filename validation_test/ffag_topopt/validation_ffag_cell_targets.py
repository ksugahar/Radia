"""Bell--Abell soft-edge FFAG target-family validation.

This fast lane checks the reduced target contract before a costly HDiv-MMM
material solve: every momentum must advance by exactly one cell, meet the
periodic position/tangent planes, and reproduce its first-order map through
the analytic response transform.  It does not claim to reproduce the paper's
unpublished PTC placement or to validate a realized 3-D magnet.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from radia.accelerator_magnet_topopt import static_magnet_symplectic_residual
from radia.ffag_topopt import (
    FFAGSoftEdgeCellSpec,
    build_ffag_cell_target_family,
)


DEFAULT_RESULT = Path(__file__).with_name("results_ffag_cell_targets.json")


def _machine_tune(trace, cell_count):
    trace = float(trace)
    if abs(trace) >= 2.0:
        return None
    return float(cell_count * np.arccos(0.5 * trace) / (2.0 * np.pi))


def run(output=DEFAULT_RESULT):
    spec = FFAGSoftEdgeCellSpec.bell_abell(full_gap_m=0.10)
    energies = np.asarray([31.0, 50.0, 80.0, 120.0,
                           160.0, 205.0, 250.0])
    family = build_ffag_cell_target_family(
        energies, spec=spec, n_segments=256,
        transfer_matrix_band=2.0e-3, bend_field_band=2.0e-3)
    raw = np.concatenate([
        reference.field_response for reference in family.references])
    transform_residual = float(np.max(np.abs(
        family.objective.transform(raw) - family.objective.response_target)))
    records = []
    for reference in family.references:
        radius = np.linalg.norm(reference.orbit.positions[:, :2], axis=1)
        radial = reference.orbit.positions[0, :2]
        radial /= np.linalg.norm(radial)
        azimuthal = np.array([-radial[1], radial[0]])
        tangent = reference.orbit.tangents[0, :2]
        incidence = float(np.arctan2(
            tangent @ radial, tangent @ azimuthal))
        records.append({
            "kinetic_energy_mev": reference.kinetic_energy_mev,
            "magnetic_rigidity_tm": reference.magnetic_rigidity_tm,
            "reduced_transverse_offset_m": reference.transverse_offset_m,
            "minimum_orbit_radius_m": float(np.min(radius)),
            "maximum_orbit_radius_m": float(np.max(radius)),
            "entrance_incidence_angle_deg": float(np.degrees(incidence)),
            "bend_angle_error_rad": float(
                reference.bend_angle_rad - spec.cell_bend_angle_rad),
            "periodic_position_residual_m": (
                reference.periodic_position_residual_m),
            "periodic_tangent_residual": (
                reference.periodic_tangent_residual),
            "radial_trace": reference.transfer.optics.radial_trace,
            "vertical_trace": reference.transfer.optics.vertical_trace,
            "radial_machine_tune": _machine_tune(
                reference.transfer.optics.radial_trace, spec.cell_count),
            "vertical_machine_tune": _machine_tune(
                reference.transfer.optics.vertical_trace, spec.cell_count),
            "symplectic_residual": static_magnet_symplectic_residual(
                reference.transfer.matrix),
        })
    maximum_bend_error = max(abs(item["bend_angle_error_rad"])
                             for item in records)
    maximum_position_residual = max(
        item["periodic_position_residual_m"] for item in records)
    maximum_tangent_residual = max(
        item["periodic_tangent_residual"] for item in records)
    maximum_symplectic_residual = max(
        item["symplectic_residual"] for item in records)
    gates = {
        "cell_bend": maximum_bend_error < 1.0e-12,
        "periodic_position": maximum_position_residual < 1.0e-12,
        "periodic_tangent": maximum_tangent_residual < 1.0e-12,
        "analytic_map_transform": transform_residual < 1.0e-10,
        "symplectic_map": maximum_symplectic_residual < 1.0e-10,
    }
    result = {
        "schema": "radia.ffag-soft-edge-cell-targets/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "Reduced soft-edge periodic target family; the realized 3-D "
            "HDiv-MMM closed-orbit solve is a separate acceptance lane."),
        "source_model": "Bell--Abell arXiv:1202.0805 Table 1",
        "cell_count": spec.cell_count,
        "cell_length_m": spec.cell_length_m,
        "cell_bend_angle_rad": spec.cell_bend_angle_rad,
        "fringe_epsilon_m": spec.fringe_epsilon_m,
        "assumed_full_gap_m": spec.full_gap_m,
        "enge_i1": family.fringe_integrals.i1,
        "enge_i2": family.fringe_integrals.i2,
        "enge_equal_integral_residual": (
            family.fringe_integrals.equal_integral_residual),
        "raw_field_row_count": family.objective.raw_field_response_size,
        "design_response_count": family.objective.response_target.size,
        "maximum_bend_angle_error_rad": maximum_bend_error,
        "maximum_periodic_position_residual_m": maximum_position_residual,
        "maximum_periodic_tangent_residual": maximum_tangent_residual,
        "maximum_analytic_map_transform_residual": transform_residual,
        "maximum_symplectic_residual": maximum_symplectic_residual,
        "gates": gates,
        "operating_points": records,
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
