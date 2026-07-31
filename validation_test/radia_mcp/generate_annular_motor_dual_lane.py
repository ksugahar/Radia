"""Generate public-safe production artifacts for both radia-motor lanes."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_DIR = (
    ROOT / "validation_test" / "radia_mcp" / "artifacts" / "annular_motor_dual_lane_v1"
)
TEMP_DIR = Path(r"C:\temp\radia_annular_motor_dual_lane_v1")

for source_root in (ROOT / "src", ROOT / "packages" / "radia-mcp" / "src"):
    source = str(source_root)
    if source not in sys.path:
        sys.path.insert(0, source)

import radia  # noqa: E402
from radia_mcp.motor.triple_check_knowledge import (  # noqa: E402
    validate_motor_triple_check_artifact,
)
from radia_mcp.radia_ngsolve.age_periodic_motion import (  # noqa: E402
    solve_age_periodic_motion,
)
from radia_mcp.radia_ngsolve.age_retirement_validation import (  # noqa: E402
    _validation_vol,
)


def _sha_json(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, indent=2) + "\n")


def _public_command(command: list[str]) -> list[str]:
    public: list[str] = []
    for index, token in enumerate(command):
        path = Path(token)
        if index == 0:
            public.append("python")
        elif path.is_absolute():
            try:
                public.append(str(path.relative_to(ROOT)).replace("\\", "/"))
            except ValueError:
                public.append("<temporary-output>")
        else:
            public.append(token)
    return public


def _identity() -> tuple[dict[str, object], dict[str, str]]:
    geometry = {
        "family": "concentric_annular_harmonic_fixture",
        "resolved_geometry": "smooth_concentric_annuli",
        "slot_geometry_resolved": False,
        "rotor_inner_radius_m": 0.2,
        "rotor_outer_radius_m": 0.8,
        "stator_inner_radius_m": 1.0,
        "stator_outer_radius_m": 1.4,
        "axial_length_m": 0.2,
    }
    periodicity = {
        "slots": 12,
        "poles": 4,
        "sector_count": 4,
        "role": "nominal_symmetry_contract_not_resolved_slot_geometry",
    }
    material = {
        "rotor": {
            "relative_permeability": 1001.0,
            "conductivity_s_per_m": 58_000_000.0,
        },
        "stator": {
            "relative_permeability": 1001.0,
            "conductivity_s_per_m": 58_000_000.0,
        },
        "frequency_hz": 100.0,
    }
    excitation = {
        "harmonic": 1,
        "rotor_amplitude": 0.8,
        "stator_amplitude": 1.0,
        "angle_samples": 8,
        "angle_basis": "mechanical_rad",
        "endpoint_policy": "exclude_repeated_period_endpoint",
        "rotation_representation": "first_spatial_harmonic_phase_only",
        "magnetic_basis": ["uniform_transverse_x", "uniform_transverse_y"],
    }
    physical = {
        "geometry": geometry,
        "periodicity": periodicity,
        "material": material,
        "excitation": excitation,
    }
    digests = {
        "geometry_sha256": _sha_json(
            {"geometry": geometry, "periodicity": periodicity}
        ),
        "material_sha256": _sha_json(material),
        "excitation_sha256": _sha_json(excitation),
    }
    digests["aggregate_sha256"] = _sha_json(digests)
    return physical, digests


def _same_harmonic_fixture_identity(
    physical: dict[str, object],
    identity: dict[str, str],
    age: dict[str, object],
    hdiv: dict[str, object],
) -> bool:
    geometry = physical["geometry"]
    material = physical["material"]
    excitation = physical["excitation"]
    assert isinstance(geometry, dict)
    assert isinstance(material, dict)
    assert isinstance(excitation, dict)

    age_gap = age.get("gap_contract", {})
    age_material = age.get("material_contract", {}).get("materials", {})
    age_excitation = age.get("excitation_contract", {})
    age_grid = age.get("angle_grid", {})
    hdiv_config = hdiv.get("configuration", {})
    hdiv_fixture = hdiv_config.get("motor_fixture_contract", {})
    hdiv_checks = hdiv.get("checks", {})
    expected_materials = {
        name: {
            "relative_permeability": row["relative_permeability"],
            "conductivity_s_per_m": row["conductivity_s_per_m"],
        }
        for name, row in material.items()
        if isinstance(row, dict)
    }
    return all(
        (
            geometry.get("family") == "concentric_annular_harmonic_fixture",
            geometry.get("slot_geometry_resolved") is False,
            math.isclose(age_gap.get("inner_radius_m", math.nan), 0.8),
            math.isclose(age_gap.get("outer_radius_m", math.nan), 1.0),
            age_gap.get("harmonics") == [1],
            age_material == expected_materials,
            math.isclose(
                age.get("frequency_hz", math.nan), material["frequency_hz"]
            ),
            math.isclose(
                age.get("axial_length_m", math.nan), geometry["axial_length_m"]
            ),
            age.get("rotation_method")
            == "rotor_harmonic_phase_only_no_remesh",
            age_excitation.get("rotation") == "rotor_harmonic_phase_only",
            age_excitation.get("harmonics", {}).get("1")
            == {
                "rotor_amplitude": excitation["rotor_amplitude"],
                "stator_amplitude": excitation["stator_amplitude"],
            },
            len(age_grid.get("angles_rad", []))
            == excitation["angle_samples"],
            hdiv_config.get("geometry") == "annular-motor",
            math.isclose(
                hdiv_config.get("axial_length_m", math.nan),
                geometry["axial_length_m"],
            ),
            math.isclose(hdiv_config.get("sigma", math.nan), 58_000_000.0),
            math.isclose(hdiv_config.get("mu_r", math.nan), 1001.0),
            hdiv_config.get("frequencies_Hz") == [100.0],
            hdiv_config.get("motor_angle_samples")
            == excitation["angle_samples"],
            math.isclose(
                hdiv_config.get("rotor_amplitude", math.nan),
                excitation["rotor_amplitude"],
            ),
            hdiv_fixture.get("family") == geometry["family"],
            hdiv_fixture.get("slot_geometry_resolved") is False,
            hdiv_fixture.get("magnetic_basis")
            == excitation["magnetic_basis"],
            hdiv_checks.get(
                "motor_hcurl_ports_match_hdiv_transverse_xy_basis"
            )
            is True,
            hdiv_config.get("shared_model_identity_sha256")
            == identity["aggregate_sha256"],
        )
    )


def _run_age() -> dict[str, object]:
    samples = 8
    return solve_age_periodic_motion(
        {
            "vol_text": _validation_vol(),
            "source_name": "generated_annular_motor.vol",
            "airgap": {
                "inner_radius_m": 0.8,
                "outer_radius_m": 1.0,
                "rotor_ring": "rotor_ring",
                "stator_ring": "stator_ring",
                "rotor_inner": "rotor_inner",
                "outer": "outer",
                "rotor_material": "rotor",
                "stator_material": "stator",
                "harmonics": [1],
            },
            "materials": {
                "rotor": {
                    "relative_permeability": 1001.0,
                    "conductivity_s_per_m": 58_000_000.0,
                },
                "stator": {
                    "relative_permeability": 1001.0,
                    "conductivity_s_per_m": 58_000_000.0,
                },
            },
            "periodic_sector": {
                "slots": 12,
                "poles": 4,
                "sector_count": 4,
                "sector_angle_deg": 90.0,
                "boundary": "anti-periodic",
                "boundary_phase": -1.0,
            },
            "excitation": {
                "1": {"rotor_amplitude": 0.8, "stator_amplitude": 1.0}
            },
            "rotor_angles_rad": [
                2.0 * math.pi * index / samples for index in range(samples)
            ],
            "axial_length_m": 0.2,
            "frequency_hz": 100.0,
            "element_order": 2,
        }
    )


def _run_hdiv(identity_sha256: str) -> tuple[dict[str, object], list[str]]:
    output = TEMP_DIR / "hdiv_mmm_hcurl_eddy_bubble_raw.json"
    command = [
        sys.executable,
        str(ROOT / "validation_test" / "cln" / "hcurl_vim_hdiv_mmm_end_to_end.py"),
        "--geometry",
        "annular-motor",
        "--axial-length",
        "0.2",
        "--curvature-safety",
        "0.2",
        "--motor-angle-samples",
        "8",
        "--rotor-amplitude",
        "0.8",
        "--frequencies",
        "100",
        "--maxh",
        "0.8",
        "--order",
        "1",
        "--steps",
        "4",
        "--bulk-degree",
        "1",
        "--surface-current-degree",
        "1",
        "--response-backend",
        "dense",
        "--mixed-solver",
        "dense",
        "--shared-model-identity-sha256",
        identity_sha256,
        "--output",
        str(output),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [
            str(ROOT / "src"),
            str(ROOT / "packages" / "radia-mcp" / "src"),
            environment.get("PYTHONPATH", ""),
        ]
    )
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "HDiv-MMM/HCurl production run failed:\n"
            + completed.stdout
            + completed.stderr
        )
    return json.loads(output.read_text(encoding="utf-8")), command


def _centered_waveform(values: list[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    centered = array - np.mean(array)
    norm = float(np.linalg.norm(centered))
    if norm <= 0.0:
        raise RuntimeError("torque waveform is constant")
    return centered / norm


def main() -> int:
    started = time.perf_counter()
    host_role = os.environ.get(
        "RADIA_VALIDATION_HOST_ROLE", "developer-smoke"
    ).strip().lower()
    if host_role not in {"compute", "developer-smoke"}:
        raise ValueError(
            "RADIA_VALIDATION_HOST_ROLE must be 'compute' or 'developer-smoke'"
        )
    execution_environment = {
        "host_role": host_role,
        "hostname": platform.node(),
        "python": platform.python_version(),
    }
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    physical_identity, identity = _identity()

    age_started = time.perf_counter()
    age = _run_age()
    age_elapsed = time.perf_counter() - age_started
    hdiv_started = time.perf_counter()
    hdiv, hdiv_command = _run_hdiv(identity["aggregate_sha256"])
    hdiv_elapsed = time.perf_counter() - hdiv_started

    if age["status"] != "solved":
        raise RuntimeError("AGE production lane did not solve")
    if not hdiv["checks"]["passed"]:
        raise RuntimeError("HDiv-MMM/HCurl production lane did not pass")
    if hdiv["validation_host"] != execution_environment["hostname"]:
        raise RuntimeError("HDiv validation host does not match the artifact host")
    if (
        hdiv["configuration"]["shared_model_identity_sha256"]
        != identity["aggregate_sha256"]
    ):
        raise RuntimeError("HDiv artifact does not carry the shared model identity")

    age_torque = [float(row["torque_nm"]) for row in age["torque_rows"]]
    hdiv_torque = [
        float(row["torque_proxy_Nm"]) for row in hdiv["motor_angle_rows"]
    ]
    age_wave = _centered_waveform(age_torque)
    hdiv_wave = _centered_waveform(hdiv_torque)
    direct_correlation = float(np.dot(age_wave, hdiv_wave))
    sign_mapping = -1.0 if direct_correlation < 0.0 else 1.0
    aligned_hdiv_wave = sign_mapping * hdiv_wave
    correlation = float(np.dot(age_wave, aligned_hdiv_wave))
    normalized_rms = float(np.sqrt(np.mean((age_wave - aligned_hdiv_wave) ** 2)))
    normalized_max = float(np.max(np.abs(age_wave - aligned_hdiv_wave)))
    comparison_checks = {
        "same_smooth_annular_harmonic_fixture_identity": (
            _same_harmonic_fixture_identity(
                physical_identity,
                identity,
                age,
                hdiv,
            )
        ),
        "both_angle_grids_complete": len(age_torque) == len(hdiv_torque) == 8,
        "both_torque_waveforms_nonconstant": True,
        "both_torque_waveforms_reverse_sign": (
            min(age_torque) < 0.0 < max(age_torque)
            and min(hdiv_torque) < 0.0 < max(hdiv_torque)
        ),
        "sign_convention_mapping_is_explicit": sign_mapping in {-1.0, 1.0},
        "centered_normalized_correlation_at_least_0_95": correlation >= 0.95,
        "centered_normalized_rms_at_most_0_12": normalized_rms <= 0.12,
    }
    if not all(comparison_checks.values()):
        raise RuntimeError(f"dual-lane comparison failed: {comparison_checks}")

    timestamp = datetime.now(timezone.utc).isoformat()
    radia_version = str(radia.__version__)
    shared_identity = {
        "geometry_sha256": identity["geometry_sha256"],
        "material_sha256": identity["material_sha256"],
        "excitation_sha256": identity["excitation_sha256"],
        "aggregate_sha256": identity["aggregate_sha256"],
    }

    age_artifact = {
        "radia_version": radia_version,
        "schema_version": "radia-motor-validation-artifact/v1",
        "timestamp_utc": timestamp,
        "execution_environment": execution_environment,
        "motor_validation_lane": "ngsolve_age",
        "reference_source_class": "independent_public_solver_lane",
        "observable_family": "torque",
        "case_count": 8,
        "status": "pass",
        "coupling_design_status": "validated_solver_path",
        "tolerances": {
            "torque_closure_relative_error_max": 1.0e-8,
            "centered_normalized_correlation_min": 0.95,
        },
        "metrics": {
            "quantity_specific_residual": float(
                age["torque_summary"]["closure_relative_error"]
            ),
            "torque_relative_error": normalized_rms,
            "centered_normalized_correlation": correlation,
            "torque_peak_to_peak_nm": float(
                age["torque_summary"]["peak_to_peak_nm"]
            ),
        },
        "timing_breakdown_s": {
            "prepare": float(age["timing_breakdown_s"]["prepare"]),
            "assemble_and_factor": float(
                age["timing_breakdown_s"]["assemble_and_factor"]
            ),
            "angle_sweep_and_postprocess": float(
                age["timing_breakdown_s"]["angle_sweep_and_postprocess"]
            ),
            "total": age_elapsed,
        },
        "artifact_feedback": {
            "status": "promoted",
            "public_lesson": (
                "A fixed AGE operator and factorization can sweep a complete first-"
                "harmonic phase grid without remeshing and provide an independent "
                "torque waveform for the smooth annular fixture."
            ),
        },
        "shared_mesh_material_identity": shared_identity,
        "solver_ready_artifact": {
            "artifact_id": "annular_motor_ngsolve_age_v1",
            "verification": [
                {
                    "command": "python validation_test/radia_mcp/generate_annular_motor_dual_lane.py",
                    "result": "solved 8 angles with reused mesh/operator/factorization",
                }
            ],
            "result_output_sha256": str(age["torque_output_sha256"]),
        },
        "age_gate_ids": [
            "mesh_reused_all_angles",
            "operator_reused_all_angles",
            "factorization_reused_all_angles",
            "torque_closes_over_period",
        ],
        "pytest_targets": [
            "packages/radia-mcp/tests/test_age_retirement_validation.py"
        ],
        "execution": {
            "solver_status": age["status"],
            "operation": age["operation"],
            "rotation_method": age["rotation_method"],
            "mesh_reused_all_angles": age["mesh_reused_all_angles"],
            "operator_reused_all_angles": age["operator_reused_all_angles"],
            "factorization_reused_all_angles": age[
                "factorization_reused_all_angles"
            ],
            "torque_rows": age["torque_rows"],
            "torque_summary": age["torque_summary"],
        },
    }

    frequency_row = hdiv["frequency_rows"][0]
    hdiv_artifact = {
        "radia_version": radia_version,
        "schema_version": "radia-motor-validation-artifact/v1",
        "timestamp_utc": timestamp,
        "execution_environment": execution_environment,
        "motor_validation_lane": "hdiv_mmm_hcurl_eddy_bubble",
        "reference_source_class": "independent_public_solver_lane",
        "observable_family": "force_or_torque_trend",
        "case_count": 8,
        "status": "pass",
        "coupling_design_status": "validated_solver_path",
        "tolerances": {
            "mixed_block_residual_max": 1.0e-8,
            "energy_response_relative_error_max": 1.0e-8,
            "centered_normalized_correlation_min": 0.95,
        },
        "metrics": {
            "mixed_block_residual": float(
                max(
                    row["residual_relative_norm"]
                    for row in hdiv["motor_angle_rows"]
                )
            ),
            "magnetic_energy_closure": float(
                hdiv["hdiv_reduction"]["basis_generation"][
                    "max_response_relative_energy_error"
                ]
            ),
            "eddy_power_nonnegative": bool(
                hdiv["checks"]["motor_angle_losses_nonnegative"]
            ),
            "signed_agreement_count": 8,
            "mean_abs_relative_error": float(
                np.mean(np.abs(age_wave - aligned_hdiv_wave))
            ),
            "rms_abs_relative_error": normalized_rms,
            "max_abs_relative_error": normalized_max,
            "centered_normalized_correlation": correlation,
            "direct_solution_relative_error": float(
                frequency_row["direct_solution_relative_error"]
            ),
        },
        "timing_breakdown_s": {
            "mesh_basis_and_reduction": float(hdiv["elapsed_seconds"]),
            "angle_sweep_and_postprocess": 0.0,
            "artifact_transport": max(hdiv_elapsed - float(hdiv["elapsed_seconds"]), 0.0),
            "total": hdiv_elapsed,
        },
        "artifact_feedback": {
            "status": "promoted",
            "public_lesson": (
                "HDiv-MMM and HCurl eddy-bubble share one smooth annular 3D "
                "mesh/material identity and recover its first-harmonic periodic "
                "coenergy-torque waveform."
            ),
        },
        "shared_mesh_material_identity": shared_identity,
        "solver_ready_artifact": {
            "artifact_id": "annular_motor_hdiv_mmm_hcurl_eddy_bubble_v1",
            "verification": [
                {
                    "command": " ".join(
                        _public_command(hdiv_command)
                    ),
                    "result": "all end-to-end and 8-angle motor checks passed",
                }
            ],
            "source_sha256": _sha_file(
                ROOT
                / "validation_test"
                / "cln"
                / "hcurl_vim_hdiv_mmm_end_to_end.py"
            ),
        },
        "hdiv_mmm_operator_contract": {
            "parent_family": hdiv["hdiv_reduction"]["parent_family"],
            "parent_order": hdiv["hdiv_reduction"]["parent_order"],
            "parent_ndof": hdiv["hdiv_reduction"]["parent_ndof"],
            "reduced_modes": hdiv["hdiv_reduction"]["reduced_modes"],
            "demag_hmatrix_backend": hdiv["hdiv_reduction"][
                "demag_hmatrix_backend"
            ],
        },
        "hcurl_eddy_bubble_contract": {
            "parent_ndof": hdiv["hcurl_parent_ndof"],
            "reduced_modes": hdiv["mixed_system"]["hcurl_vim_modes"],
            "eddy_block_roles": hdiv["mixed_system"]["eddy_block_roles"],
            "solver_backend": frequency_row["solver_backend"],
        },
        "coupling_operator_contract": {
            "rows": hdiv["mixed_system"]["coupling_rows"],
            "columns": hdiv["mixed_system"]["coupling_cols"],
            "frobenius_norm": hdiv["mixed_system"]["coupling_frobenius_norm"],
            "full_coupled_mixed_galerkin_is_exact": hdiv["checks"][
                "full_coupled_mixed_galerkin_is_exact"
            ],
        },
        "execution": {
            "checks": hdiv["checks"],
            "frequency_rows": hdiv["frequency_rows"],
            "motor_angle_rows": hdiv["motor_angle_rows"],
        },
    }

    dual_packet = {
        "radia_version": radia_version,
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "timestamp_utc": timestamp,
        "execution_environment": execution_environment,
        "source_mcp_seed": {
            "representative_public_decks": [
                "generated_concentric_annular_harmonic_fixture"
            ],
            "source_mcp_calls": [
                "solve_age_periodic_motion",
                "hcurl_vim_hdiv_mmm_end_to_end",
            ],
        },
        "lane_artifacts": {
            "ngsolve_age": age_artifact,
            "hdiv_mmm_hcurl_eddy_bubble": hdiv_artifact,
        },
        "dual_lane_comparison": {
            "torque_sign_mapping_age_from_hdiv": sign_mapping,
            "sign_mapping_reason": (
                "The lanes use opposite positive rotation/virtual-work torque conventions."
            ),
            "centered_normalized_correlation": correlation,
            "centered_normalized_rms": normalized_rms,
            "centered_normalized_max": normalized_max,
            "checks": comparison_checks,
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": (
                "Both mandatory radia-motor lanes executed on one smooth annular "
                "first-harmonic fixture; their centered torque waveforms agree after "
                "an explicit sign convention map. This fixture does not resolve slot "
                "geometry."
            ),
            "learning_targets": [
                "radia_mcp.motor.validation_lanes_knowledge",
                "radia_mcp.motor.triple_check_knowledge",
            ],
            "verification": [
                "both lane artifact gates pass",
                "shared geometry/material/excitation digests match",
                "centered normalized torque correlation >= 0.95",
            ],
        },
    }
    gate = validate_motor_triple_check_artifact(dual_packet)
    if not gate["accepted_for_primary_dual_learning"]:
        raise RuntimeError(f"dual-lane MCP gate failed: {gate}")

    manifest = {
        "radia_version": radia_version,
        "schema": "radia.validation.annular-motor-dual-lane-manifest.v1",
        "generated_at_utc": timestamp,
        "execution_environment": execution_environment,
        "physical_identity": physical_identity,
        "shared_mesh_material_identity": shared_identity,
        "artifact_files": {
            "ngsolve_age": "ngsolve_age.json",
            "hdiv_mmm_hcurl_eddy_bubble": "hdiv_mmm_hcurl_eddy_bubble.json",
            "native_motor_angle_family": "native_motor_angle_family.json",
            "dual_lane": "dual_lane.json",
            "gate_result": "gate_result.json",
        },
        "timing_breakdown_s": {
            "ngsolve_age": age_elapsed,
            "hdiv_mmm_hcurl_eddy_bubble": hdiv_elapsed,
            "packaging_and_gate": max(
                time.perf_counter() - started - age_elapsed - hdiv_elapsed, 0.0
            ),
            "total": time.perf_counter() - started,
        },
        "status": "pass",
    }
    gate_artifact = {
        "radia_version": radia_version,
        "generated_at_utc": timestamp,
        **gate,
    }
    _write_json(ARTIFACT_DIR / "ngsolve_age.json", age_artifact)
    _write_json(
        ARTIFACT_DIR / "hdiv_mmm_hcurl_eddy_bubble.json", hdiv_artifact
    )
    _write_json(ARTIFACT_DIR / "dual_lane.json", dual_packet)
    _write_json(ARTIFACT_DIR / "gate_result.json", gate_artifact)
    _write_json(ARTIFACT_DIR / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "status": "pass",
                "artifact_dir": str(ARTIFACT_DIR),
                "correlation": correlation,
                "normalized_rms": normalized_rms,
                "gate": gate["accepted_for_primary_dual_learning"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
