"""Jointed assembly STEP closure and source-replay validation gates."""

from __future__ import annotations

import math
from datetime import datetime
from typing import Mapping


def _relative_error(measured: float, reference: float) -> float:
    if not math.isfinite(reference) or reference <= 0.0 or not math.isfinite(measured):
        return math.inf
    return abs(measured - reference) / reference


def _timestamp(value: object) -> datetime | None:
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None


def jointed_assembly_step_closure_gate(
    summary: Mapping[str, object],
    *,
    tessellated_volume_rtol: float = 2.0e-4,
    self_roundtrip_volume_rtol: float = 1.0e-6,
    external_volume_rtol: float = 2.0e-5,
) -> dict[str, object]:
    """Diagnose component-level solid loss in a jointed assembly STEP handoff.

    ``status=ok`` means the diagnosis is supported by a portable component,
    a rejected component, and independent volume evidence.  It does not mean
    the assembly is solver-ready; that is reported separately.
    """

    tolerances = (
        float(tessellated_volume_rtol),
        float(self_roundtrip_volume_rtol),
        float(external_volume_rtol),
    )
    if any(value < 0.0 or not math.isfinite(value) for value in tolerances):
        raise ValueError("volume tolerances must be finite and nonnegative")

    raw_components = summary.get("components") or []
    if not isinstance(raw_components, list) or len(raw_components) < 2:
        raise ValueError("components must contain at least two rows")
    external_rows = summary.get("external_rows") or []
    if not isinstance(external_rows, list):
        raise ValueError("external_rows must be a list")
    external_by_name = {str(row.get("name", "")): row for row in external_rows}

    components = []
    portable = []
    rejected = []
    for raw in raw_components:
        name = str(raw.get("name", "")).strip()
        native_volume = float(raw.get("native_volume_mm3", math.nan))
        tessellated_volume = float(raw.get("tessellated_volume_mm3", math.nan))
        self_roundtrip = raw.get("self_roundtrip") or {}
        self_volume = float(self_roundtrip.get("total_volume_mm3", math.nan))
        self_solid_count = int(self_roundtrip.get("solid_count", -1))
        external = external_by_name.get(f"{name}.step", {})
        external_volume = float(external.get("total_volume_mm3", math.nan))
        external_volume_count = int(external.get("volume_count", -1))
        tess_error = _relative_error(tessellated_volume, native_volume)
        self_error = _relative_error(self_volume, native_volume)
        external_error = _relative_error(external_volume, native_volume)
        native_supported = (
            bool(name)
            and raw.get("native_valid") is True
            and int(raw.get("native_solid_count", 0)) == 1
            and native_volume > 0.0
            and tess_error <= tolerances[0]
        )
        is_portable = (
            native_supported
            and self_solid_count == 1
            and self_error <= tolerances[1]
            and external_volume_count >= 1
            and external_volume > 0.0
            and external_error <= tolerances[2]
        )
        is_closure_loss = (
            native_supported
            and self_solid_count == 0
            and self_volume == 0.0
            and external_volume_count >= 1
            and external_volume == 0.0
        )
        expected = str(raw.get("expected_disposition", ""))
        disposition = (
            "portable_control"
            if is_portable
            else "reject_solid_closure_loss"
            if is_closure_loss
            else "unresolved"
        )
        row = {
            "name": name,
            "expected_disposition": expected,
            "observed_disposition": disposition,
            "native_volume_mm3": native_volume,
            "tessellated_volume_relative_error": tess_error,
            "self_roundtrip_volume_relative_error": self_error,
            "external_volume_relative_error": external_error,
            "self_roundtrip_solid_count": self_solid_count,
            "external_volume_count": external_volume_count,
            "external_volume_mm3": external_volume,
            "disposition_matches": disposition == expected,
        }
        components.append(row)
        if is_portable:
            portable.append(row)
        if is_closure_loss:
            rejected.append(row)

    assembly = summary.get("assembly") or {}
    assembly_external = external_by_name.get(str(assembly.get("step_name", "")), {})
    native_total = float(assembly.get("native_total_volume_mm3", math.nan))
    self_total = float((assembly.get("self_roundtrip") or {}).get("total_volume_mm3", math.nan))
    external_total = float(assembly_external.get("total_volume_mm3", math.nan))
    rejected_native_total = sum(float(row["native_volume_mm3"]) for row in rejected)
    self_lost = native_total - self_total
    external_lost = native_total - external_total
    names = [str(row["name"]) for row in components]
    checks = {
        "component_names_are_nonempty_and_unique": all(names) and len(set(names)) == len(names),
        "external_rows_cover_every_component_and_assembly": all(
            f"{name}.step" in external_by_name for name in names
        )
        and str(assembly.get("step_name", "")) in external_by_name,
        "native_brep_supported_by_tessellated_volume": all(
            float(row["tessellated_volume_relative_error"]) <= tolerances[0]
            for row in components
        ),
        "portable_component_control_present": bool(portable),
        "solid_closure_loss_component_present": bool(rejected),
        "component_dispositions_match_expectations": all(
            bool(row["disposition_matches"]) for row in components
        ),
        "assembly_self_loss_matches_rejected_components": rejected_native_total > 0.0
        and _relative_error(self_lost, rejected_native_total) <= tolerances[1],
        "assembly_external_loss_matches_rejected_components": rejected_native_total > 0.0
        and _relative_error(external_lost, rejected_native_total) <= tolerances[2],
    }
    issues = [name for name, ok in checks.items() if not ok]
    diagnosis = "component_solid_closure_loss" if not issues else "incomplete_component_evidence"
    return {
        "policy": "build123d_jointed_assembly_step_closure_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "diagnosis": diagnosis,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "components": components,
        "assembly": {
            "native_total_volume_mm3": native_total,
            "self_roundtrip_total_volume_mm3": self_total,
            "external_total_volume_mm3": external_total,
            "rejected_component_native_volume_mm3": rejected_native_total,
            "self_roundtrip_lost_volume_mm3": self_lost,
            "external_lost_volume_mm3": external_lost,
        },
        "notes": [
            "Assembly total volume alone cannot identify which source component lost solid closure.",
            "Require a portable component control and an independently tessellated source volume before blaming STEP translation.",
            "Do not promote an assembly when any positive-volume source component returns as a zero-volume external entity.",
        ],
    }


def jointed_assembly_heal_invariance_gate(
    summary: Mapping[str, object],
    *,
    volume_rtol: float = 2.0e-5,
) -> dict[str, object]:
    """Check whether external STEP healing changes a component closure diagnosis.

    A CAD import may create a nominal ``volume`` entity whose measured volume is
    zero.  The gate therefore uses positive measured volume, not entity count,
    and requires both ``noheal`` and ``heal`` evidence for every component.
    """

    tolerance = float(volume_rtol)
    if tolerance < 0.0 or not math.isfinite(tolerance):
        raise ValueError("volume_rtol must be finite and nonnegative")
    components = summary.get("components") or []
    rows = summary.get("external_rows") or []
    if not isinstance(components, list) or len(components) < 2:
        raise ValueError("components must contain at least two rows")
    if not isinstance(rows, list):
        raise ValueError("external_rows must be a list")

    indexed = {
        (str(row.get("name", "")), str(row.get("import_mode", ""))): row
        for row in rows
    }
    observed = []
    portable_total = 0.0
    rejected_total = 0.0
    for component in components:
        name = str(component.get("name", "")).strip()
        native_volume = float(component.get("native_volume_mm3", math.nan))
        expected = str(component.get("expected_disposition", ""))
        mode_rows = [indexed.get((f"{name}.step", mode), {}) for mode in ("noheal", "heal")]
        measured = [float(row.get("total_volume_mm3", math.nan)) for row in mode_rows]
        positive_counts = [
            int(row.get("positive_volume_count", -1)) for row in mode_rows
        ]
        errors = [_relative_error(value, native_volume) for value in measured]
        portable = (
            native_volume > 0.0
            and all(count >= 1 for count in positive_counts)
            and all(error <= tolerance for error in errors)
        )
        persistent_loss = (
            native_volume > 0.0
            and positive_counts == [0, 0]
            and measured == [0.0, 0.0]
        )
        disposition = (
            "portable_control"
            if portable
            else "reject_persistent_solid_closure_loss"
            if persistent_loss
            else "heal_sensitive_or_unresolved"
        )
        if portable:
            portable_total += native_volume
        if persistent_loss:
            rejected_total += native_volume
        observed.append(
            {
                "name": name,
                "expected_disposition": expected,
                "observed_disposition": disposition,
                "native_volume_mm3": native_volume,
                "noheal_volume_mm3": measured[0],
                "heal_volume_mm3": measured[1],
                "noheal_relative_error": errors[0],
                "heal_relative_error": errors[1],
                "disposition_matches": disposition == expected,
            }
        )

    assembly = summary.get("assembly") or {}
    assembly_name = str(assembly.get("step_name", ""))
    assembly_native = float(assembly.get("native_total_volume_mm3", math.nan))
    assembly_rows = [indexed.get((assembly_name, mode), {}) for mode in ("noheal", "heal")]
    assembly_measured = [
        float(row.get("total_volume_mm3", math.nan)) for row in assembly_rows
    ]
    names = [str(row["name"]) for row in observed]
    checks = {
        "component_names_are_nonempty_and_unique": all(names)
        and len(set(names)) == len(names),
        "both_import_modes_cover_every_artifact": all(
            (f"{name}.step", mode) in indexed
            for name in names
            for mode in ("noheal", "heal")
        )
        and all((assembly_name, mode) in indexed for mode in ("noheal", "heal")),
        "portable_control_present": portable_total > 0.0,
        "persistent_closure_loss_present": rejected_total > 0.0,
        "component_dispositions_match_expectations": all(
            bool(row["disposition_matches"]) for row in observed
        ),
        "native_assembly_is_component_sum": _relative_error(
            portable_total + rejected_total, assembly_native
        )
        <= tolerance,
        "assembly_loss_persists_after_heal": all(
            _relative_error(value, portable_total) <= tolerance
            for value in assembly_measured
        ),
        "assembly_lost_volume_matches_rejected_components": all(
            _relative_error(assembly_native - value, rejected_total) <= tolerance
            for value in assembly_measured
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_jointed_assembly_heal_invariance_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "diagnosis": "persistent_solid_closure_loss"
        if not issues
        else "heal_sensitive_or_incomplete_evidence",
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "components": observed,
        "assembly": {
            "native_total_volume_mm3": assembly_native,
            "noheal_total_volume_mm3": assembly_measured[0],
            "heal_total_volume_mm3": assembly_measured[1],
            "portable_component_volume_mm3": portable_total,
            "rejected_component_volume_mm3": rejected_total,
        },
        "notes": [
            "A nominal external volume entity with measured volume zero is not a closed solver-ready solid.",
            "Do not assume an import heal option repairs a closure loss; measure both paths.",
            "Keep a portable component in the same assembly as a translation control.",
        ],
    }
def jointed_assembly_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate immutable source replay, joint identity, and headless CAD diagnosis."""

    components = summary.get("components") or []
    connections = summary.get("joint_connections") or []
    execution = summary.get("external_execution") or {}
    timing = summary.get("timing_breakdown_s") or {}
    replay_identity_value = summary.get("replay_identity")
    replay_identity_present = isinstance(replay_identity_value, Mapping)
    commit_identity_ok = True
    kernel_version_identity_ok = True
    export_follows_replay_ok = True
    kernel_session_identity_ok = True
    topology_replay_identity_ok = True
    unit_conversion_identity_ok = True
    boolean_clean_identity_ok = True
    tessellation_identity_ok = True
    step_export_tolerance_identity_ok = True
    assembly_replacement_identity_ok = True
    assembly_mass_property_coordinate_identity_ok = True
    boolean_final_shape_report_identity_ok = True
    step_geometry_unit_scale_identity_ok = True
    selector_cache_shape_identity_ok = True
    step_assembly_placement_unit_identity_ok = True
    brep_serialization_tolerance_kernel_identity_ok = True
    step_color_label_topology_identity_ok = True
    brep_surface_parameter_orientation_identity_ok = True
    step_assembly_child_parent_unit_identity_ok = True
    selector_normal_world_frame_identity_ok = True
    step_occurrence_hierarchy_identity_ok = True
    brep_edge_tolerance_shape_fix_identity_ok = True
    step_occurrence_color_material_inheritance_identity_ok = True
    stl_tolerance_model_length_unit_generation_identity_ok = True
    step_import_tolerance_unit_healing_generation_identity_ok = True
    tessellation_vertex_index_normal_transform_generation_identity_ok = True
    brep_serialization_shape_digest_occt_location_generation_identity_ok = True
    dxf_wire_plane_orientation_layer_generation_identity_ok = True
    step_assembly_product_color_location_generation_identity_ok = True
    sketch_constraint_entity_solver_order_generation_identity_ok = True
    step_ap242_component_identity_ok = True
    curved_mesh_export_identity_ok = True
    if replay_identity_value is not None and not replay_identity_present:
        commit_identity_ok = False
        kernel_version_identity_ok = False
        export_follows_replay_ok = False
        kernel_session_identity_ok = False
        topology_replay_identity_ok = False
        unit_conversion_identity_ok = False
        boolean_clean_identity_ok = False
        tessellation_identity_ok = False
        step_export_tolerance_identity_ok = False
        assembly_replacement_identity_ok = False
        assembly_mass_property_coordinate_identity_ok = False
        boolean_final_shape_report_identity_ok = False
        step_geometry_unit_scale_identity_ok = False
        selector_cache_shape_identity_ok = False
        step_assembly_placement_unit_identity_ok = False
        brep_serialization_tolerance_kernel_identity_ok = False
        step_color_label_topology_identity_ok = False
        brep_surface_parameter_orientation_identity_ok = False
        step_occurrence_hierarchy_identity_ok = False
        brep_edge_tolerance_shape_fix_identity_ok = False
        step_occurrence_color_material_inheritance_identity_ok = False
        stl_tolerance_model_length_unit_generation_identity_ok = False
        step_import_tolerance_unit_healing_generation_identity_ok = False
        tessellation_vertex_index_normal_transform_generation_identity_ok = False
        brep_serialization_shape_digest_occt_location_generation_identity_ok = False
        dxf_wire_plane_orientation_layer_generation_identity_ok = False
        step_assembly_product_color_location_generation_identity_ok = False
        sketch_constraint_entity_solver_order_generation_identity_ok = False
        step_ap242_component_identity_ok = False
        curved_mesh_export_identity_ok = False
    elif replay_identity_present:
        source_commit = str(replay_identity_value.get("source_commit", "")).lower()
        replayed_commit = str(
            replay_identity_value.get("replayed_source_commit", "")
        ).lower()
        artifacts = replay_identity_value.get("cad_artifacts") or []
        kernel = replay_identity_value.get("external_kernel") or {}
        artifact_rows_valid = isinstance(artifacts, list) and all(
            isinstance(row, Mapping) for row in artifacts
        )
        artifact_names = (
            [str(row.get("name", "")).strip() for row in artifacts]
            if artifact_rows_valid
            else []
        )
        valid_commit = (
            len(source_commit) == 40
            and all(character in "0123456789abcdef" for character in source_commit)
            and source_commit == replayed_commit
        )
        commit_identity_ok = (
            valid_commit
            and artifact_rows_valid
            and bool(artifacts)
            and all(artifact_names)
            and len(set(artifact_names)) == len(artifact_names)
            and all(
                row.get("fresh") is True
                and len(str(row.get("sha256", "")).lower()) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in str(row.get("sha256", "")).lower()
                )
                and str(row.get("source_commit", "")).lower() == source_commit
                for row in artifacts
            )
        )
        kernel_valid = isinstance(kernel, Mapping)
        claimed_version = (
            str(kernel.get("claimed_version", "")).strip() if kernel_valid else ""
        )
        replay_versions = [
            str(value).strip() for value in (kernel.get("replay_versions") or [])
        ] if kernel_valid else []
        kernel_version_identity_ok = (
            kernel_valid
            and bool(str(kernel.get("name", "")).strip())
            and bool(claimed_version)
            and len(replay_versions) >= 2
            and all(version == claimed_version for version in replay_versions)
        )
        replay_started = _timestamp(replay_identity_value.get("source_replay_started_utc"))
        export_times_present = artifact_rows_valid and any(
            row.get("export_completed_utc") is not None for row in artifacts
        )
        freshness_evidence_present = (
            replay_identity_value.get("source_replay_started_utc") is not None
            or export_times_present
        )
        export_times = (
            [_timestamp(row.get("export_completed_utc")) for row in artifacts]
            if artifact_rows_valid
            else []
        )
        export_follows_replay_ok = not freshness_evidence_present or (
            replay_started is not None
            and bool(export_times)
            and all(value is not None and value >= replay_started for value in export_times)
        )

        session_evidence_present = kernel_valid and (
            kernel.get("claimed_session_generation") is not None
            or kernel.get("replay_sessions") is not None
        )
        claimed_session = (
            str(kernel.get("claimed_session_generation", "")).strip()
            if kernel_valid
            else ""
        )
        replay_sessions = kernel.get("replay_sessions") if kernel_valid else None
        session_rows_valid = isinstance(replay_sessions, list) and all(
            isinstance(row, Mapping) for row in replay_sessions
        )
        session_starts = (
            {str(row.get("process_start_utc", "")).strip() for row in replay_sessions}
            if session_rows_valid
            else set()
        )
        kernel_session_identity_ok = not session_evidence_present or (
            bool(claimed_session)
            and session_rows_valid
            and len(replay_sessions) >= 2
            and all(
                str(row.get("session_generation", "")).strip() == claimed_session
                for row in replay_sessions
            )
            and len(session_starts) == 1
            and all(session_starts)
        )

        topology_value = replay_identity_value.get("topology_replay_identity")
        if topology_value is not None:
            topology = topology_value if isinstance(topology_value, Mapping) else {}
            source_topology = str(topology.get("source_topology_sha256", "")).lower()
            imports = topology.get("imports")
            import_rows_valid = isinstance(imports, list) and all(
                isinstance(row, Mapping) for row in imports
            )
            import_modes = (
                {str(row.get("mode", "")) for row in imports}
                if import_rows_valid
                else set()
            )
            topology_replay_identity_ok = (
                len(source_topology) == 64
                and import_rows_valid
                and import_modes == {"heal", "noheal"}
                and all(
                    str(row.get("topology_sha256", "")).lower()
                    == source_topology
                    for row in imports
                )
            )

        unit_value = replay_identity_value.get("unit_conversion_identity")
        if unit_value is not None:
            unit_identity = unit_value if isinstance(unit_value, Mapping) else {}
            try:
                length_scale = float(unit_identity.get("length_scale_to_target"))
                volume_scale = float(
                    unit_identity.get("declared_volume_scale_to_target")
                )
            except (TypeError, ValueError):
                length_scale = math.nan
                volume_scale = math.nan
            unit_conversion_identity_ok = (
                unit_identity.get("source_geometry_unit") == "mm"
                and unit_identity.get("target_geometry_unit") == "m"
                and unit_identity.get("external_measurement_stage")
                == "after_unit_conversion"
                and unit_identity.get("external_volume_unit") == "m^3"
                and math.isfinite(length_scale)
                and math.isfinite(volume_scale)
                and length_scale == 0.001
                and volume_scale == 1.0
            )

        boolean_value = replay_identity_value.get("boolean_clean_identity")
        if boolean_value is not None:
            boolean_identity = (
                boolean_value if isinstance(boolean_value, Mapping) else {}
            )
            boolean_digest = str(
                boolean_identity.get("boolean_result_sha256", "")
            ).lower()
            cleaned_digest = str(
                boolean_identity.get("cleaned_topology_sha256", "")
            ).lower()
            boolean_clean_identity_ok = (
                len(boolean_digest) == 64
                and boolean_identity.get("shape_clean_input_sha256")
                == boolean_digest
                and len(cleaned_digest) == 64
                and boolean_identity.get("export_topology_sha256")
                == cleaned_digest
                and bool(boolean_identity.get("shape_generation"))
                and boolean_identity.get("export_shape_generation")
                == boolean_identity.get("shape_generation")
            )

        tessellation_value = replay_identity_value.get("tessellation_identity")
        if tessellation_value is not None:
            tessellation = (
                tessellation_value if isinstance(tessellation_value, Mapping) else {}
            )
            try:
                linear_deflection = float(tessellation.get("linear_deflection"))
                angular_tolerance = float(
                    tessellation.get("angular_tolerance_rad")
                )
            except (TypeError, ValueError):
                linear_deflection = math.nan
                angular_tolerance = math.nan
            tessellation_identity_ok = (
                bool(tessellation.get("shape_generation"))
                and tessellation.get("tolerance_shape_generation")
                == tessellation.get("shape_generation")
                and math.isfinite(linear_deflection)
                and linear_deflection > 0.0
                and math.isfinite(angular_tolerance)
                and 0.0 < angular_tolerance <= math.pi
                and bool(tessellation.get("tessellation_generation"))
            )

        tolerance_value = replay_identity_value.get(
            "step_export_tolerance_identity"
        )
        if tolerance_value is not None:
            tolerance_identity = (
                tolerance_value if isinstance(tolerance_value, Mapping) else {}
            )
            try:
                sewing_tolerance = float(
                    tolerance_identity.get("sewing_tolerance")
                )
                brep_tolerance = float(tolerance_identity.get("brep_tolerance"))
            except (TypeError, ValueError):
                sewing_tolerance = math.nan
                brep_tolerance = math.nan
            current_session = str(
                kernel.get("claimed_session_generation", "")
            ).strip()
            boolean_identity = (
                boolean_value if isinstance(boolean_value, Mapping) else {}
            )
            current_shape = str(
                boolean_identity.get("shape_generation", "")
            ).strip()
            export_digest = str(
                tolerance_identity.get("export_artifact_sha256", "")
            ).lower()
            step_export_tolerance_identity_ok = (
                bool(current_session)
                and tolerance_identity.get("kernel_session_generation")
                == current_session
                and tolerance_identity.get("tolerance_kernel_session_generation")
                == current_session
                and bool(current_shape)
                and tolerance_identity.get("shape_generation") == current_shape
                and tolerance_identity.get("tolerance_shape_generation")
                == current_shape
                and math.isfinite(sewing_tolerance)
                and sewing_tolerance > 0.0
                and math.isfinite(brep_tolerance)
                and brep_tolerance > 0.0
                and len(export_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in export_digest
                )
            )

        replacement_value = replay_identity_value.get(
            "assembly_replacement_identity"
        )
        if replacement_value is not None:
            replacement = (
                replacement_value if isinstance(replacement_value, Mapping) else {}
            )
            assembly_generation = str(
                replacement.get("assembly_generation", "")
            ).strip()
            replacement_rows = replacement.get("components")
            rows_valid = isinstance(replacement_rows, list) and bool(replacement_rows) and all(
                isinstance(row, Mapping) for row in replacement_rows
            )
            slot_ids = (
                [str(row.get("slot_id", "")).strip() for row in replacement_rows]
                if rows_valid
                else []
            )

            def valid_digest(value: object) -> bool:
                digest = str(value or "").lower()
                return len(digest) == 64 and all(
                    character in "0123456789abcdef" for character in digest
                )

            assembly_replacement_identity_ok = (
                bool(assembly_generation)
                and rows_valid
                and all(slot_ids)
                and len(set(slot_ids)) == len(slot_ids)
                and all(
                    bool(str(row.get("replacement_generation", "")).strip())
                    and bool(str(row.get("removed_instance_uuid", "")).strip())
                    and bool(str(row.get("current_instance_uuid", "")).strip())
                    and row.get("removed_instance_uuid")
                    != row.get("current_instance_uuid")
                    and valid_digest(row.get("removed_shape_sha256"))
                    and valid_digest(row.get("current_shape_sha256"))
                    and row.get("removed_shape_sha256")
                    != row.get("current_shape_sha256")
                    and row.get("placement_shape_sha256")
                    == row.get("current_shape_sha256")
                    and row.get("placement_assembly_generation")
                    == assembly_generation
                    for row in replacement_rows
                )
            )

        coordinate_value = replay_identity_value.get(
            "assembly_mass_property_coordinate_identity"
        )
        if coordinate_value is not None:
            coordinate = coordinate_value if isinstance(coordinate_value, Mapping) else {}
            placement_generation = str(
                coordinate.get("placement_matrix_generation", "")
            ).strip()
            coordinate_frame = str(coordinate.get("coordinate_frame_id", "")).strip()
            placement_digest = str(
                coordinate.get("placement_matrix_sha256", "")
            ).lower()
            assembly_mass_property_coordinate_identity_ok = (
                bool(coordinate.get("assembly_generation"))
                and bool(placement_generation)
                and coordinate.get("centroid_transform_generation")
                == placement_generation
                and coordinate.get("inertia_transform_generation")
                == placement_generation
                and bool(coordinate_frame)
                and coordinate.get("centroid_coordinate_frame_id")
                == coordinate_frame
                and coordinate.get("inertia_coordinate_frame_id")
                == coordinate_frame
                and len(placement_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in placement_digest
                )
                and coordinate.get("centroid_placement_matrix_sha256")
                == placement_digest
                and coordinate.get("inertia_placement_matrix_sha256")
                == placement_digest
            )

        final_shape_value = replay_identity_value.get(
            "boolean_final_shape_report_identity"
        )
        if final_shape_value is not None:
            final_shape = (
                final_shape_value if isinstance(final_shape_value, Mapping) else {}
            )
            final_generation = str(
                final_shape.get("final_shape_generation", "")
            ).strip()
            pre_heal_digest = str(
                final_shape.get("pre_heal_brep_sha256", "")
            ).lower()
            final_digest = str(final_shape.get("final_brep_sha256", "")).lower()
            digests_valid = all(
                len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                for digest in (pre_heal_digest, final_digest)
            )
            boolean_final_shape_report_identity_ok = (
                bool(final_shape.get("boolean_result_generation"))
                and bool(final_shape.get("healing_generation"))
                and bool(final_generation)
                and digests_valid
                and pre_heal_digest != final_digest
                and final_shape.get("mass_property_shape_generation")
                == final_generation
                and final_shape.get("validity_shape_generation") == final_generation
                and final_shape.get("topology_shape_generation") == final_generation
                and final_shape.get("mass_property_brep_sha256") == final_digest
                and final_shape.get("validity_brep_sha256") == final_digest
                and final_shape.get("topology_brep_sha256") == final_digest
            )

        step_unit = replay_identity_value.get("step_geometry_unit_scale_identity")
        if step_unit is not None:
            step_unit = step_unit if isinstance(step_unit, Mapping) else {}
            units = {"m": 1.0, "mm": 1.0e-3, "um": 1.0e-6}
            geometry_unit = str(step_unit.get("geometry_length_unit", "")).strip()
            metadata_unit = str(step_unit.get("metadata_length_unit", "")).strip()
            try:
                geometry_scale = float(step_unit.get("geometry_scale_to_m"))
                metadata_scale = float(step_unit.get("metadata_scale_to_m"))
            except (TypeError, ValueError):
                geometry_scale = math.nan
                metadata_scale = math.nan
            generation = str(step_unit.get("step_import_generation", "")).strip()
            step_geometry_unit_scale_identity_ok = (
                geometry_unit in units and metadata_unit == geometry_unit
                and math.isclose(geometry_scale, units[geometry_unit], rel_tol=0.0, abs_tol=0.0)
                and math.isclose(metadata_scale, geometry_scale, rel_tol=0.0, abs_tol=0.0)
                and bool(generation) and step_unit.get("geometry_coordinate_generation") == generation
                and step_unit.get("metadata_generation") == generation
            )

        selector_cache = replay_identity_value.get("selector_cache_shape_identity")
        if selector_cache is not None:
            selector_cache = selector_cache if isinstance(selector_cache, Mapping) else {}
            generation = str(selector_cache.get("active_shape_generation", "")).strip()
            digest = str(selector_cache.get("selector_query_sha256", "")).lower()
            selector_cache_shape_identity_ok = (
                bool(generation) and selector_cache.get("selector_cache_shape_generation") == generation
                and selector_cache.get("selected_face_shape_generation") == generation
                and len(digest) == 64 and all(character in "0123456789abcdef" for character in digest)
                and selector_cache.get("cached_selector_query_sha256") == digest
                and list(selector_cache.get("selected_face_ids") or []) == list(selector_cache.get("live_face_ids") or [])
                and bool(selector_cache.get("selected_face_ids"))
            )

        placement_unit = replay_identity_value.get(
            "step_assembly_placement_unit_identity"
        )
        if placement_unit is not None:
            placement_unit = placement_unit if isinstance(placement_unit, Mapping) else {}
            units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            geometry_unit = str(placement_unit.get("part_geometry_length_unit", "")).strip()
            placement_translation_unit = str(
                placement_unit.get("placement_translation_unit", "")
            ).strip()
            generation = str(placement_unit.get("step_import_generation", "")).strip()
            transform_digest = str(
                placement_unit.get("placement_transform_sha256", "")
            ).lower()
            try:
                geometry_scale = float(placement_unit.get("part_geometry_scale_to_m"))
                placement_scale = float(
                    placement_unit.get("placement_translation_scale_to_m")
                )
            except (TypeError, ValueError):
                geometry_scale = math.nan
                placement_scale = math.nan
            expected_scale = units.get(geometry_unit)
            step_assembly_placement_unit_identity_ok = (
                expected_scale is not None
                and placement_translation_unit == geometry_unit
                and math.isclose(geometry_scale, expected_scale, rel_tol=0.0, abs_tol=0.0)
                and math.isclose(placement_scale, geometry_scale, rel_tol=0.0, abs_tol=0.0)
                and bool(generation)
                and placement_unit.get("part_geometry_generation") == generation
                and placement_unit.get("placement_metadata_generation") == generation
                and len(transform_digest) == 64
                and all(character in "0123456789abcdef" for character in transform_digest)
                and placement_unit.get("applied_transform_sha256") == transform_digest
            )

        brep_cache = replay_identity_value.get(
            "brep_serialization_tolerance_kernel_identity"
        )
        if brep_cache is not None:
            brep_cache = brep_cache if isinstance(brep_cache, Mapping) else {}
            kernel_generation = str(
                brep_cache.get("active_kernel_generation", "")
            ).strip()
            shape_digest = str(brep_cache.get("shape_sha256", "")).lower()
            try:
                modeling_tolerance = float(brep_cache.get("modeling_tolerance_value"))
                cache_tolerance = float(brep_cache.get("cache_tolerance_value"))
            except (TypeError, ValueError):
                modeling_tolerance = math.nan
                cache_tolerance = math.nan
            brep_serialization_tolerance_kernel_identity_ok = (
                bool(str(brep_cache.get("cache_generation", "")).strip())
                and bool(kernel_generation)
                and brep_cache.get("serialization_kernel_generation") == kernel_generation
                and brep_cache.get("modeling_tolerance_unit") == "m"
                and brep_cache.get("cache_tolerance_unit") == "m"
                and math.isfinite(modeling_tolerance)
                and modeling_tolerance > 0.0
                and math.isclose(cache_tolerance, modeling_tolerance, rel_tol=0.0, abs_tol=0.0)
                and len(shape_digest) == 64
                and all(character in "0123456789abcdef" for character in shape_digest)
                and brep_cache.get("cached_shape_sha256") == shape_digest
            )

        step_attributes = replay_identity_value.get(
            "step_color_label_topology_identity"
        )
        if step_attributes is not None:
            step_attributes = (
                step_attributes if isinstance(step_attributes, Mapping) else {}
            )
            topology_generation = str(
                step_attributes.get("topology_generation", "")
            ).strip()
            face_ids = list(step_attributes.get("face_ids") or [])
            attribute_digest = str(
                step_attributes.get("attribute_map_sha256", "")
            ).lower()
            step_color_label_topology_identity_ok = (
                bool(str(step_attributes.get("step_import_generation", "")).strip())
                and bool(topology_generation)
                and step_attributes.get("attribute_map_topology_generation")
                == topology_generation
                and bool(face_ids)
                and len(set(face_ids)) == len(face_ids)
                and list(step_attributes.get("attribute_face_ids") or []) == face_ids
                and len(list(step_attributes.get("labels") or [])) == len(face_ids)
                and len(list(step_attributes.get("colors_rgb") or []))
                == len(face_ids)
                and len(attribute_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in attribute_digest
                )
                and step_attributes.get("resolved_attribute_map_sha256")
                == attribute_digest
            )

        surface_parameters = replay_identity_value.get(
            "brep_surface_parameter_orientation_identity"
        )
        if surface_parameters is not None:
            surface_parameters = (
                surface_parameters
                if isinstance(surface_parameters, Mapping)
                else {}
            )
            serialization_generation = str(
                surface_parameters.get("serialization_generation", "")
            ).strip()
            surface_ids = list(surface_parameters.get("surface_ids") or [])
            parameter_digest = str(
                surface_parameters.get("parameter_range_sha256", "")
            ).lower()
            parameter_orientation = str(
                surface_parameters.get("parameter_orientation", "")
            ).strip()
            brep_surface_parameter_orientation_identity_ok = (
                bool(serialization_generation)
                and surface_parameters.get("surface_parameter_generation")
                == serialization_generation
                and bool(surface_ids)
                and len(set(surface_ids)) == len(surface_ids)
                and list(surface_parameters.get("exported_surface_ids") or [])
                == surface_ids
                and parameter_orientation == "u_cross_v_outward"
                and surface_parameters.get("exported_parameter_orientation")
                == parameter_orientation
                and len(parameter_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in parameter_digest
                )
                and surface_parameters.get("exported_parameter_range_sha256")
                == parameter_digest
            )

        assembly_units = replay_identity_value.get(
            "step_assembly_child_parent_unit_identity"
        )
        if assembly_units is not None:
            assembly_units = (
                assembly_units if isinstance(assembly_units, Mapping) else {}
            )
            import_generation = str(
                assembly_units.get("step_import_generation", "")
            ).strip()
            unit = str(assembly_units.get("assembly_length_unit", "")).strip()
            unit_scales = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            digest = str(
                assembly_units.get("assembly_placement_sha256", "")
            ).lower()
            try:
                assembly_scale = float(assembly_units.get("assembly_scale_to_m"))
                child_scale = float(
                    assembly_units.get("child_placement_scale_to_m")
                )
                parent_scale = float(
                    assembly_units.get("parent_placement_scale_to_m")
                )
            except (TypeError, ValueError):
                assembly_scale = child_scale = parent_scale = math.nan
            expected_scale = unit_scales.get(unit)
            step_assembly_child_parent_unit_identity_ok = (
                bool(import_generation)
                and assembly_units.get("child_placement_import_generation")
                == import_generation
                and assembly_units.get("parent_placement_import_generation")
                == import_generation
                and expected_scale is not None
                and assembly_units.get("child_placement_length_unit") == unit
                and assembly_units.get("parent_placement_length_unit") == unit
                and math.isclose(
                    assembly_scale, expected_scale, rel_tol=0.0, abs_tol=0.0
                )
                and child_scale == assembly_scale
                and parent_scale == assembly_scale
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and assembly_units.get("resolved_assembly_placement_sha256")
                == digest
            )

        selector_frame = replay_identity_value.get(
            "selector_normal_world_frame_identity"
        )
        if selector_frame is not None:
            selector_frame = (
                selector_frame if isinstance(selector_frame, Mapping) else {}
            )
            shape_generation = str(
                selector_frame.get("shape_generation", "")
            ).strip()
            placement_generation = str(
                selector_frame.get("placement_generation", "")
            ).strip()
            face_ids = list(selector_frame.get("selected_face_ids") or [])
            digest = str(selector_frame.get("normal_table_sha256", "")).lower()
            selector_normal_world_frame_identity_ok = (
                bool(shape_generation)
                and selector_frame.get("selector_shape_generation")
                == shape_generation
                and bool(placement_generation)
                and selector_frame.get("selector_placement_generation")
                == placement_generation
                and selector_frame.get("normal_predicate_frame") == "world"
                and selector_frame.get("evaluated_normal_frame") == "world"
                and bool(face_ids)
                and len(set(face_ids)) == len(face_ids)
                and list(selector_frame.get("resolved_face_ids") or []) == face_ids
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and selector_frame.get("evaluated_normal_table_sha256") == digest
            )

        occurrence_hierarchy = replay_identity_value.get(
            "step_occurrence_name_color_hierarchy_identity"
        )
        if occurrence_hierarchy is not None:
            occurrence_hierarchy = (
                occurrence_hierarchy
                if isinstance(occurrence_hierarchy, Mapping)
                else {}
            )
            import_generation = str(
                occurrence_hierarchy.get("step_import_generation", "")
            ).strip()
            hierarchy_generation = str(
                occurrence_hierarchy.get("assembly_hierarchy_generation", "")
            ).strip()
            occurrence_ids = list(occurrence_hierarchy.get("occurrence_ids") or [])
            names = list(occurrence_hierarchy.get("occurrence_names") or [])
            colors = list(
                occurrence_hierarchy.get("occurrence_colors_rgb") or []
            )
            parent_paths = list(
                occurrence_hierarchy.get("occurrence_parent_paths") or []
            )
            digest = str(
                occurrence_hierarchy.get("hierarchy_metadata_sha256", "")
            ).lower()
            try:
                colors_ok = all(
                    len(color) == 3
                    and all(
                        math.isfinite(float(channel))
                        and 0.0 <= float(channel) <= 1.0
                        for channel in color
                    )
                    for color in colors
                )
            except (TypeError, ValueError):
                colors_ok = False
            step_occurrence_hierarchy_identity_ok = (
                bool(import_generation)
                and occurrence_hierarchy.get(
                    "occurrence_metadata_import_generation"
                )
                == import_generation
                and bool(hierarchy_generation)
                and occurrence_hierarchy.get("occurrence_hierarchy_generation")
                == hierarchy_generation
                and bool(occurrence_ids)
                and len(set(occurrence_ids)) == len(occurrence_ids)
                and list(occurrence_hierarchy.get("metadata_occurrence_ids") or [])
                == occurrence_ids
                and len(names) == len(colors) == len(parent_paths) == len(occurrence_ids)
                and all(str(name).strip() for name in names)
                and all(str(path).strip() for path in parent_paths)
                and colors_ok
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and occurrence_hierarchy.get("imported_hierarchy_metadata_sha256")
                == digest
            )

        edge_tolerance = replay_identity_value.get(
            "brep_edge_tolerance_shape_fix_topology_identity"
        )
        if edge_tolerance is not None:
            edge_tolerance = (
                edge_tolerance if isinstance(edge_tolerance, Mapping) else {}
            )
            shape_fix_generation = str(
                edge_tolerance.get("shape_fix_generation", "")
            ).strip()
            topology_digest = str(
                edge_tolerance.get("topology_sha256", "")
            ).lower()
            try:
                edge_ids = [int(value) for value in edge_tolerance.get("edge_ids", [])]
                tolerance_edge_ids = [
                    int(value)
                    for value in edge_tolerance.get("edge_tolerance_edge_ids", [])
                ]
                tolerances = [
                    float(value)
                    for value in edge_tolerance.get("edge_tolerances_m", [])
                ]
                topology_edge_count = int(
                    edge_tolerance.get("topology_edge_count", -1)
                )
            except (TypeError, ValueError):
                edge_ids = []
                tolerance_edge_ids = []
                tolerances = []
                topology_edge_count = -1
            brep_edge_tolerance_shape_fix_identity_ok = (
                bool(shape_fix_generation)
                and edge_tolerance.get("edge_tolerance_shape_fix_generation")
                == shape_fix_generation
                and edge_tolerance.get("topology_digest_shape_fix_generation")
                == shape_fix_generation
                and bool(edge_ids)
                and len(set(edge_ids)) == len(edge_ids)
                and tolerance_edge_ids == edge_ids
                and len(tolerances) == len(edge_ids) == topology_edge_count
                and all(
                    math.isfinite(tolerance) and tolerance > 0.0
                    for tolerance in tolerances
                )
                and len(topology_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in topology_digest
                )
                and edge_tolerance.get("edge_tolerance_topology_sha256")
                == topology_digest
            )

        occurrence_inheritance = replay_identity_value.get(
            "step_occurrence_color_material_inheritance_identity"
        )
        if occurrence_inheritance is not None:
            occurrence_inheritance = (
                occurrence_inheritance
                if isinstance(occurrence_inheritance, Mapping)
                else {}
            )
            import_generation = str(
                occurrence_inheritance.get("step_import_generation", "")
            ).strip()
            assembly_generation = str(
                occurrence_inheritance.get("assembly_generation", "")
            ).strip()
            occurrence_ids = list(occurrence_inheritance.get("occurrence_ids") or [])
            metadata_ids = list(
                occurrence_inheritance.get("metadata_occurrence_ids") or []
            )
            parent_ids = list(
                occurrence_inheritance.get("parent_occurrence_ids") or []
            )
            metadata_parent_ids = list(
                occurrence_inheritance.get("metadata_parent_occurrence_ids") or []
            )
            inherited_colors = list(
                occurrence_inheritance.get("inherited_colors_rgb") or []
            )
            imported_colors = list(
                occurrence_inheritance.get("imported_colors_rgb") or []
            )
            inherited_materials = list(
                occurrence_inheritance.get("inherited_material_names") or []
            )
            imported_materials = list(
                occurrence_inheritance.get("imported_material_names") or []
            )
            digest = str(
                occurrence_inheritance.get("occurrence_metadata_sha256", "")
            ).lower()
            try:
                colors_ok = all(
                    len(color) == 3
                    and all(
                        math.isfinite(float(channel))
                        and 0.0 <= float(channel) <= 1.0
                        for channel in color
                    )
                    for color in inherited_colors
                )
            except (TypeError, ValueError):
                colors_ok = False
            step_occurrence_color_material_inheritance_identity_ok = (
                bool(import_generation)
                and occurrence_inheritance.get(
                    "occurrence_metadata_import_generation"
                )
                == import_generation
                and bool(assembly_generation)
                and occurrence_inheritance.get(
                    "occurrence_metadata_assembly_generation"
                )
                == assembly_generation
                and occurrence_inheritance.get(
                    "color_inheritance_assembly_generation"
                )
                == assembly_generation
                and occurrence_inheritance.get(
                    "material_inheritance_assembly_generation"
                )
                == assembly_generation
                and bool(occurrence_ids)
                and len(set(occurrence_ids)) == len(occurrence_ids)
                and metadata_ids == occurrence_ids
                and len(parent_ids) == len(occurrence_ids)
                and metadata_parent_ids == parent_ids
                and len(inherited_colors) == len(occurrence_ids)
                and imported_colors == inherited_colors
                and colors_ok
                and len(inherited_materials) == len(occurrence_ids)
                and all(str(material).strip() for material in inherited_materials)
                and imported_materials == inherited_materials
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and occurrence_inheritance.get(
                    "imported_occurrence_metadata_sha256"
                )
                == digest
            )

        stl_tolerance = replay_identity_value.get(
            "stl_tolerance_model_length_unit_generation_identity"
        )
        if stl_tolerance is not None:
            stl_tolerance = stl_tolerance if isinstance(stl_tolerance, Mapping) else {}
            model_unit_generation = str(
                stl_tolerance.get("model_length_unit_generation", "")
            ).strip()
            model_unit = str(stl_tolerance.get("model_length_unit", "")).strip()
            chordal_unit = str(
                stl_tolerance.get("chordal_tolerance_length_unit", "")
            ).strip()
            angular_unit = str(
                stl_tolerance.get("angular_tolerance_unit", "")
            ).strip()
            length_scales = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            angular_scales = {"rad": 1.0, "deg": math.pi / 180.0}
            digest = str(stl_tolerance.get("tolerance_contract_sha256", "")).lower()
            try:
                chordal_value = float(stl_tolerance.get("chordal_tolerance_value"))
                chordal_si = float(stl_tolerance.get("chordal_tolerance_si_m"))
                tessellator_chordal_si = float(
                    stl_tolerance.get("tessellator_chordal_tolerance_si_m")
                )
                angular_value = float(stl_tolerance.get("angular_tolerance_value"))
                angular_rad = float(stl_tolerance.get("angular_tolerance_rad"))
                tessellator_angular_rad = float(
                    stl_tolerance.get("tessellator_angular_tolerance_rad")
                )
            except (TypeError, ValueError):
                chordal_value = chordal_si = tessellator_chordal_si = math.nan
                angular_value = angular_rad = tessellator_angular_rad = math.nan
            length_scale = length_scales.get(chordal_unit)
            angular_scale = angular_scales.get(angular_unit)
            stl_tolerance_model_length_unit_generation_identity_ok = (
                bool(model_unit_generation)
                and stl_tolerance.get("tessellation_model_unit_generation")
                == model_unit_generation
                and stl_tolerance.get("tolerance_conversion_model_unit_generation")
                == model_unit_generation
                and model_unit in length_scales
                and chordal_unit == model_unit
                and length_scale is not None
                and angular_scale is not None
                and math.isfinite(chordal_value)
                and chordal_value > 0.0
                and math.isclose(
                    chordal_value * length_scale,
                    chordal_si,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-30,
                )
                and math.isclose(
                    tessellator_chordal_si,
                    chordal_si,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-30,
                )
                and math.isfinite(angular_value)
                and angular_value > 0.0
                and math.isclose(
                    angular_value * angular_scale,
                    angular_rad,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                and math.isclose(
                    tessellator_angular_rad,
                    angular_rad,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and stl_tolerance.get("tessellator_tolerance_contract_sha256")
                == digest
            )

        step_import_healing = replay_identity_value.get(
            "step_import_tolerance_unit_healing_generation_identity"
        )
        if step_import_healing is not None:
            step_import_healing = (
                step_import_healing
                if isinstance(step_import_healing, Mapping)
                else {}
            )
            import_generation = str(
                step_import_healing.get("step_import_generation", "")
            ).strip()
            healing_generation = str(
                step_import_healing.get("healing_generation", "")
            ).strip()
            source_unit = str(
                step_import_healing.get("source_length_unit", "")
            ).strip()
            tolerance_unit = str(
                step_import_healing.get("tolerance_length_unit", "")
            ).strip()
            length_scales = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            try:
                tolerance_value = float(step_import_healing.get("tolerance_value"))
                tolerance_si = float(step_import_healing.get("tolerance_si_m"))
                healed_edge_ids = [
                    int(item)
                    for item in step_import_healing.get("healed_edge_ids", [])
                ]
                imported_edge_ids = [
                    int(item)
                    for item in step_import_healing.get(
                        "imported_healed_edge_ids", []
                    )
                ]
            except (TypeError, ValueError):
                tolerance_value = tolerance_si = math.nan
                healed_edge_ids = []
                imported_edge_ids = []
            digest = str(
                step_import_healing.get("healed_edge_map_sha256", "")
            ).lower()
            scale = length_scales.get(tolerance_unit)
            step_import_tolerance_unit_healing_generation_identity_ok = (
                bool(import_generation)
                and step_import_healing.get("tolerance_import_generation")
                == import_generation
                and step_import_healing.get("healed_edge_map_import_generation")
                == import_generation
                and bool(healing_generation)
                and step_import_healing.get("tolerance_healing_generation")
                == healing_generation
                and step_import_healing.get("healed_edge_map_healing_generation")
                == healing_generation
                and source_unit in length_scales
                and tolerance_unit == source_unit
                and scale is not None
                and math.isfinite(tolerance_value)
                and tolerance_value > 0.0
                and math.isclose(
                    tolerance_value * scale,
                    tolerance_si,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-30,
                )
                and bool(healed_edge_ids)
                and len(set(healed_edge_ids)) == len(healed_edge_ids)
                and imported_edge_ids == healed_edge_ids
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and step_import_healing.get("imported_healed_edge_map_sha256")
                == digest
            )

        tessellation_transform = replay_identity_value.get(
            "tessellation_vertex_index_normal_transform_generation_identity"
        )
        if tessellation_transform is not None:
            tessellation_transform = (
                tessellation_transform
                if isinstance(tessellation_transform, Mapping)
                else {}
            )
            shape_generation = str(
                tessellation_transform.get("shape_generation", "")
            ).strip()
            transform_generation = str(
                tessellation_transform.get("location_transform_generation", "")
            ).strip()
            try:
                vertex_count = int(tessellation_transform.get("vertex_count", -1))
                normal_count = int(tessellation_transform.get("normal_count", -1))
                triangles = [
                    [int(item) for item in row]
                    for row in tessellation_transform.get("triangle_indices", [])
                ]
                transformed_triangles = [
                    [int(item) for item in row]
                    for row in tessellation_transform.get(
                        "transformed_triangle_indices", []
                    )
                ]
            except (TypeError, ValueError):
                vertex_count = normal_count = -1
                triangles = []
                transformed_triangles = []
            digest = str(
                tessellation_transform.get("tessellation_sha256", "")
            ).lower()
            tessellation_vertex_index_normal_transform_generation_identity_ok = (
                bool(shape_generation)
                and tessellation_transform.get("vertex_shape_generation")
                == shape_generation
                and tessellation_transform.get("index_shape_generation")
                == shape_generation
                and tessellation_transform.get("normal_shape_generation")
                == shape_generation
                and bool(transform_generation)
                and tessellation_transform.get(
                    "vertex_location_transform_generation"
                )
                == transform_generation
                and tessellation_transform.get(
                    "index_location_transform_generation"
                )
                == transform_generation
                and tessellation_transform.get(
                    "normal_location_transform_generation"
                )
                == transform_generation
                and vertex_count > 0
                and normal_count == vertex_count
                and bool(triangles)
                and all(
                    len(row) == 3
                    and len(set(row)) == 3
                    and all(0 <= item < vertex_count for item in row)
                    for row in triangles
                )
                and transformed_triangles == triangles
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and tessellation_transform.get("rendered_tessellation_sha256")
                == digest
            )

        brep_serialization = replay_identity_value.get(
            "brep_serialization_shape_digest_occt_location_generation_identity"
        )
        if brep_serialization is not None:
            brep_serialization = (
                brep_serialization
                if isinstance(brep_serialization, Mapping)
                else {}
            )
            serialization_generation = str(
                brep_serialization.get("serialization_generation", "")
            ).strip()
            shape_generation = str(
                brep_serialization.get("shape_generation", "")
            ).strip()
            kernel_version = str(
                brep_serialization.get("kernel_version", "")
            ).strip()
            location_generation = str(
                brep_serialization.get("location_generation", "")
            ).strip()
            digests = {
                key: str(brep_serialization.get(key, "")).lower()
                for key in (
                    "shape_sha256",
                    "serialized_shape_sha256",
                    "deserialized_shape_sha256",
                    "top_level_location_sha256",
                    "serialized_top_level_location_sha256",
                    "deserialized_top_level_location_sha256",
                    "brep_payload_sha256",
                    "deserialized_brep_payload_sha256",
                )
            }
            digests_valid = all(
                len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                for digest in digests.values()
            )
            brep_serialization_shape_digest_occt_location_generation_identity_ok = (
                bool(serialization_generation)
                and brep_serialization.get(
                    "deserialization_serialization_generation"
                )
                == serialization_generation
                and bool(shape_generation)
                and brep_serialization.get("serialized_shape_generation")
                == shape_generation
                and brep_serialization.get("deserialized_shape_generation")
                == shape_generation
                and bool(kernel_version)
                and brep_serialization.get("serialized_kernel_version")
                == kernel_version
                and brep_serialization.get("deserialized_kernel_version")
                == kernel_version
                and bool(location_generation)
                and brep_serialization.get("serialized_location_generation")
                == location_generation
                and brep_serialization.get("deserialized_location_generation")
                == location_generation
                and digests_valid
                and digests["serialized_shape_sha256"] == digests["shape_sha256"]
                and digests["deserialized_shape_sha256"] == digests["shape_sha256"]
                and digests["serialized_top_level_location_sha256"]
                == digests["top_level_location_sha256"]
                and digests["deserialized_top_level_location_sha256"]
                == digests["top_level_location_sha256"]
                and digests["deserialized_brep_payload_sha256"]
                == digests["brep_payload_sha256"]
            )

        dxf_wire = replay_identity_value.get(
            "dxf_wire_plane_orientation_layer_generation_identity"
        )
        if dxf_wire is not None:
            dxf_wire = dxf_wire if isinstance(dxf_wire, Mapping) else {}
            import_generation = str(
                dxf_wire.get("dxf_import_generation", "")
            ).strip()
            plane_generation = str(
                dxf_wire.get("plane_generation", "")
            ).strip()
            layer_generation = str(
                dxf_wire.get("layer_generation", "")
            ).strip()
            closure_generation = str(
                dxf_wire.get("closure_generation", "")
            ).strip()
            try:
                wire_ids = [int(item) for item in dxf_wire.get("wire_ids", [])]
                imported_wire_ids = [
                    int(item) for item in dxf_wire.get("imported_wire_ids", [])
                ]
            except (TypeError, ValueError):
                wire_ids = []
                imported_wire_ids = []
            layers = [str(item) for item in dxf_wire.get("wire_layers", [])]
            imported_layers = [
                str(item) for item in dxf_wire.get("imported_wire_layers", [])
            ]
            closed = list(dxf_wire.get("wire_closed", []))
            imported_closed = list(dxf_wire.get("imported_wire_closed", []))
            plane_digest = str(
                dxf_wire.get("plane_orientation_sha256", "")
            ).lower()
            table_digest = str(dxf_wire.get("wire_table_sha256", "")).lower()
            dxf_wire_plane_orientation_layer_generation_identity_ok = (
                bool(import_generation)
                and all(
                    dxf_wire.get(key) == import_generation
                    for key in (
                        "wire_import_generation",
                        "plane_import_generation",
                        "layer_import_generation",
                        "closure_import_generation",
                    )
                )
                and bool(plane_generation)
                and dxf_wire.get("wire_plane_generation") == plane_generation
                and dxf_wire.get("extrusion_plane_generation")
                == plane_generation
                and bool(layer_generation)
                and dxf_wire.get("wire_layer_generation") == layer_generation
                and bool(closure_generation)
                and dxf_wire.get("wire_closure_generation")
                == closure_generation
                and bool(wire_ids)
                and len(set(wire_ids)) == len(wire_ids)
                and imported_wire_ids == wire_ids
                and len(layers) == len(wire_ids)
                and all(layer.strip() for layer in layers)
                and imported_layers == layers
                and len(closed) == len(wire_ids)
                and all(isinstance(item, bool) for item in closed)
                and imported_closed == closed
                and len(plane_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in plane_digest
                )
                and dxf_wire.get("imported_plane_orientation_sha256")
                == plane_digest
                and len(table_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in table_digest
                )
                and dxf_wire.get("imported_wire_table_sha256") == table_digest
            )

        step_assembly = replay_identity_value.get(
            "step_assembly_product_id_color_location_generation_identity"
        )
        if step_assembly is not None:
            step_assembly = (
                step_assembly if isinstance(step_assembly, Mapping) else {}
            )
            export_generation = str(
                step_assembly.get("step_export_generation", "")
            ).strip()
            assembly_generation = str(
                step_assembly.get("assembly_generation", "")
            ).strip()
            product_ids = [
                str(item).strip() for item in step_assembly.get("product_ids", [])
            ]
            decoded_product_ids = [
                str(item).strip()
                for item in step_assembly.get("decoded_product_ids", [])
            ]
            parent_ids = [
                str(item).strip()
                for item in step_assembly.get("parent_product_ids", [])
            ]
            decoded_parent_ids = [
                str(item).strip()
                for item in step_assembly.get("decoded_parent_product_ids", [])
            ]
            try:
                colors = [
                    [int(channel) for channel in color]
                    for color in step_assembly.get("colors_rgb", [])
                ]
                decoded_colors = [
                    [int(channel) for channel in color]
                    for color in step_assembly.get("decoded_colors_rgb", [])
                ]
            except (TypeError, ValueError):
                colors = []
                decoded_colors = []
            locations = [
                str(item).lower()
                for item in step_assembly.get("component_location_sha256", [])
            ]
            decoded_locations = [
                str(item).lower()
                for item in step_assembly.get(
                    "decoded_component_location_sha256", []
                )
            ]
            metadata_digest = str(
                step_assembly.get("assembly_metadata_sha256", "")
            ).lower()
            step_assembly_product_color_location_generation_identity_ok = (
                bool(export_generation)
                and step_assembly.get("decoder_step_export_generation")
                == export_generation
                and bool(assembly_generation)
                and all(
                    step_assembly.get(key) == assembly_generation
                    for key in (
                        "product_id_assembly_generation",
                        "color_assembly_generation",
                        "hierarchy_assembly_generation",
                        "location_assembly_generation",
                    )
                )
                and bool(product_ids)
                and len(set(product_ids)) == len(product_ids)
                and all(product_ids)
                and decoded_product_ids == product_ids
                and len(parent_ids) == len(product_ids)
                and parent_ids[0] == ""
                and all(
                    not parent or parent in product_ids
                    for parent in parent_ids
                )
                and decoded_parent_ids == parent_ids
                and len(colors) == len(product_ids)
                and all(
                    len(color) == 3
                    and all(0 <= channel <= 255 for channel in color)
                    for color in colors
                )
                and decoded_colors == colors
                and len(locations) == len(product_ids)
                and all(
                    len(digest) == 64
                    and all(
                        character in "0123456789abcdef" for character in digest
                    )
                    for digest in locations
                )
                and decoded_locations == locations
                and len(metadata_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in metadata_digest
                )
                and step_assembly.get("decoded_assembly_metadata_sha256")
                == metadata_digest
            )

        sketch_constraint = replay_identity_value.get(
            "sketch_constraint_entity_id_solver_order_generation_identity"
        )
        if sketch_constraint is not None:
            sketch_constraint = (
                sketch_constraint
                if isinstance(sketch_constraint, Mapping)
                else {}
            )
            sketch_generation = str(
                sketch_constraint.get("sketch_generation", "")
            ).strip()
            try:
                entity_ids = [
                    int(item) for item in sketch_constraint.get("entity_ids", [])
                ]
                replay_entity_ids = [
                    int(item)
                    for item in sketch_constraint.get("replay_entity_ids", [])
                ]
                constraint_ids = [
                    int(item)
                    for item in sketch_constraint.get("constraint_ids", [])
                ]
                solver_order = [
                    int(item)
                    for item in sketch_constraint.get(
                        "solver_constraint_order", []
                    )
                ]
                replay_solver_order = [
                    int(item)
                    for item in sketch_constraint.get(
                        "replay_solver_constraint_order", []
                    )
                ]
                constraint_entities = [
                    [int(item) for item in row]
                    for row in sketch_constraint.get("constraint_entity_ids", [])
                ]
                replay_constraint_entities = [
                    [int(item) for item in row]
                    for row in sketch_constraint.get(
                        "replay_constraint_entity_ids", []
                    )
                ]
            except (TypeError, ValueError):
                entity_ids = []
                replay_entity_ids = []
                constraint_ids = []
                solver_order = []
                replay_solver_order = []
                constraint_entities = []
                replay_constraint_entities = []
            entity_digest = str(
                sketch_constraint.get("entity_table_sha256", "")
            ).lower()
            constraint_digest = str(
                sketch_constraint.get("constraint_table_sha256", "")
            ).lower()
            sketch_constraint_entity_solver_order_generation_identity_ok = (
                bool(sketch_generation)
                and all(
                    sketch_constraint.get(key) == sketch_generation
                    for key in (
                        "entity_table_sketch_generation",
                        "constraint_table_sketch_generation",
                        "solver_order_sketch_generation",
                    )
                )
                and bool(entity_ids)
                and len(set(entity_ids)) == len(entity_ids)
                and replay_entity_ids == entity_ids
                and bool(constraint_ids)
                and len(set(constraint_ids)) == len(constraint_ids)
                and solver_order == constraint_ids
                and replay_solver_order == solver_order
                and len(constraint_entities) == len(constraint_ids)
                and all(
                    row and set(row).issubset(entity_ids)
                    for row in constraint_entities
                )
                and replay_constraint_entities == constraint_entities
                and len(entity_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in entity_digest
                )
                and sketch_constraint.get("replay_entity_table_sha256")
                == entity_digest
                and len(constraint_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in constraint_digest
                )
                and sketch_constraint.get("replay_constraint_table_sha256")
                == constraint_digest
            )

        step_ap242 = replay_identity_value.get(
            "step_ap242_component_transform_name_product_generation_identity"
        )
        if step_ap242 is not None:
            step_ap242 = step_ap242 if isinstance(step_ap242, Mapping) else {}
            export_generation = str(
                step_ap242.get("step_export_generation", "")
            ).strip()
            assembly_generation = str(
                step_ap242.get("assembly_generation", "")
            ).strip()
            product_ids = [
                str(item).strip()
                for item in step_ap242.get("component_product_ids", [])
            ]
            decoded_product_ids = [
                str(item).strip()
                for item in step_ap242.get("decoded_component_product_ids", [])
            ]
            names = [str(item).strip() for item in step_ap242.get("component_names", [])]
            decoded_names = [
                str(item).strip()
                for item in step_ap242.get("decoded_component_names", [])
            ]
            transforms = [
                str(item).lower()
                for item in step_ap242.get("nested_transform_sha256", [])
            ]
            decoded_transforms = [
                str(item).lower()
                for item in step_ap242.get("decoded_nested_transform_sha256", [])
            ]
            product_digest = str(
                step_ap242.get("ap242_product_map_sha256", "")
            ).lower()
            step_ap242_component_identity_ok = (
                bool(export_generation)
                and step_ap242.get("decoder_step_export_generation")
                == export_generation
                and bool(assembly_generation)
                and all(
                    step_ap242.get(key) == assembly_generation
                    for key in (
                        "product_assembly_generation",
                        "name_assembly_generation",
                        "transform_assembly_generation",
                    )
                )
                and bool(product_ids)
                and len(set(product_ids)) == len(product_ids)
                and all(product_ids)
                and decoded_product_ids == product_ids
                and len(names) == len(product_ids)
                and all(names)
                and decoded_names == names
                and len(transforms) == len(product_ids)
                and all(
                    len(digest) == 64
                    and all(character in "0123456789abcdef" for character in digest)
                    for digest in transforms
                )
                and decoded_transforms == transforms
                and len(product_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in product_digest
                )
                and step_ap242.get("decoded_ap242_product_map_sha256")
                == product_digest
            )

        curved_mesh = replay_identity_value.get(
            "curved_mesh_export_edge_chord_surface_label_generation_identity"
        )
        if curved_mesh is not None:
            curved_mesh = curved_mesh if isinstance(curved_mesh, Mapping) else {}
            shape_generation = str(curved_mesh.get("shape_generation", "")).strip()
            export_generation = str(
                curved_mesh.get("mesh_export_generation", "")
            ).strip()
            edge_digests = [
                str(item).lower()
                for item in curved_mesh.get("edge_curve_sha256", [])
            ]
            exported_edge_digests = [
                str(item).lower()
                for item in curved_mesh.get("exported_edge_curve_sha256", [])
            ]
            labels = [
                str(item).strip()
                for item in curved_mesh.get("boundary_surface_labels", [])
            ]
            exported_labels = [
                str(item).strip()
                for item in curved_mesh.get("exported_boundary_surface_labels", [])
            ]
            mesh_digest = str(curved_mesh.get("curved_mesh_sha256", "")).lower()
            try:
                chord = float(curved_mesh.get("chordal_tolerance"))
                evaluated_chord = float(
                    curved_mesh.get("evaluated_chordal_tolerance")
                )
            except (TypeError, ValueError):
                chord = -1.0
                evaluated_chord = -2.0
            length_unit = str(curved_mesh.get("length_unit", "")).strip()
            curved_mesh_export_identity_ok = (
                bool(shape_generation)
                and all(
                    curved_mesh.get(key) == shape_generation
                    for key in (
                        "mesh_shape_generation",
                        "edge_curve_shape_generation",
                        "surface_label_shape_generation",
                    )
                )
                and bool(export_generation)
                and curved_mesh.get("metric_mesh_export_generation")
                == export_generation
                and bool(edge_digests)
                and all(
                    len(digest) == 64
                    and all(character in "0123456789abcdef" for character in digest)
                    for digest in edge_digests
                )
                and exported_edge_digests == edge_digests
                and math.isfinite(chord)
                and chord > 0.0
                and math.isclose(evaluated_chord, chord, rel_tol=0.0, abs_tol=1.0e-18)
                and length_unit in {"m", "cm", "mm"}
                and curved_mesh.get("evaluated_length_unit") == length_unit
                and bool(labels)
                and len(set(labels)) == len(labels)
                and all(labels)
                and exported_labels == labels
                and len(mesh_digest) == 64
                and all(character in "0123456789abcdef" for character in mesh_digest)
                and curved_mesh.get("exported_curved_mesh_sha256") == mesh_digest
            )
    joint_names = {
        str(name)
        for row in components
        for name in (row.get("joint_names") or [])
        if str(name)
    }
    connection_endpoints = {
        str(connection.get(side, ""))
        for connection in connections
        for side in ("from", "to")
    }
    checks = {
        "upstream_source_identity_recorded": str(summary.get("source_kind", "")).startswith(
            "upstream_source_native_example"
        )
        and len(str(summary.get("source_sha256", ""))) == 64
        and "/v0.10.0/" in str(summary.get("source_url", "")),
        "source_preserved_except_display_stub": summary.get("source_preserved") is True
        and summary.get("display_stubbed_only") is True,
        "component_inventory_and_joint_names_recorded": len(components) >= 2
        and all(str(row.get("name", "")).strip() and row.get("joint_names") for row in components),
        "joint_connection_endpoints_resolve": bool(connections)
        and bool(connection_endpoints)
        and connection_endpoints.issubset(joint_names),
        "headless_external_cad_replay": execution.get("mode")
        == "python_api_headless_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(set(execution.get("headless_flags") or []))
        and execution.get("gui_daemon_enabled") is False,
        "fresh_result_and_owned_process_cleanup": execution.get("result_artifact_fresh") is True
        and int(execution.get("owned_processes_remaining", -1)) == 0,
        "neutral_cad_artifacts_bind_current_source_commit": commit_identity_ok,
        "external_kernel_versions_are_replay_invariant": kernel_version_identity_ok,
        "neutral_cad_export_follows_source_replay": export_follows_replay_ok,
        "external_kernel_session_generation_is_continuous": kernel_session_identity_ok,
        "heal_and_noheal_imports_record_topology_identity": topology_replay_identity_ok,
        "external_volume_is_measured_after_unit_conversion": unit_conversion_identity_ok,
        "boolean_export_follows_shape_clean_identity": boolean_clean_identity_ok,
        "tessellation_tolerances_belong_to_current_shape": tessellation_identity_ok,
        "step_export_tolerances_belong_to_current_kernel_and_shape": (
            step_export_tolerance_identity_ok
        ),
        "replacement_rotates_instance_uuid_and_rebinds_placement": (
            assembly_replacement_identity_ok
        ),
        "assembly_mass_property_report_uses_current_placement_frame": (
            assembly_mass_property_coordinate_identity_ok
        ),
        "final_shape_report_uses_one_post_heal_generation": (
            boolean_final_shape_report_identity_ok
        ),
        "step_geometry_and_metadata_share_one_length_unit": step_geometry_unit_scale_identity_ok,
        "selector_cache_belongs_to_active_shape_generation": selector_cache_shape_identity_ok,
        "step_assembly_placement_uses_current_unit_transform_generation": (
            step_assembly_placement_unit_identity_ok
        ),
        "brep_cache_uses_current_kernel_tolerance_and_shape": (
            brep_serialization_tolerance_kernel_identity_ok
        ),
        "step_color_labels_follow_current_topology_generation": (
            step_color_label_topology_identity_ok
        ),
        "brep_surface_parameter_ranges_preserve_outward_orientation": (
            brep_surface_parameter_orientation_identity_ok
        ),
        "step_assembly_child_parent_placements_share_length_unit": (
            step_assembly_child_parent_unit_identity_ok
        ),
        "selector_normals_use_final_world_placement_frame": (
            selector_normal_world_frame_identity_ok
        ),
        "step_occurrence_metadata_uses_current_assembly_hierarchy": (
            step_occurrence_hierarchy_identity_ok
        ),
        "brep_edge_tolerances_follow_final_shape_fix_topology": (
            brep_edge_tolerance_shape_fix_identity_ok
        ),
        "step_occurrence_inheritance_uses_current_assembly_generation": (
            step_occurrence_color_material_inheritance_identity_ok
        ),
        "stl_tolerances_use_current_model_length_unit_generation": (
            stl_tolerance_model_length_unit_generation_identity_ok
        ),
        "step_import_tolerances_and_healed_edges_share_current_generation": (
            step_import_tolerance_unit_healing_generation_identity_ok
        ),
        "tessellation_vertices_indices_and_normals_use_final_transform": (
            tessellation_vertex_index_normal_transform_generation_identity_ok
        ),
        "brep_deserialization_uses_current_shape_kernel_and_location": (
            brep_serialization_shape_digest_occt_location_generation_identity_ok
        ),
        "dxf_wires_use_current_plane_layer_and_closure_generations": (
            dxf_wire_plane_orientation_layer_generation_identity_ok
        ),
        "step_assembly_uses_current_product_colors_hierarchy_and_locations": (
            step_assembly_product_color_location_generation_identity_ok
        ),
        "sketch_constraints_use_current_entity_ids_and_solver_order": (
            sketch_constraint_entity_solver_order_generation_identity_ok
        ),
        "step_ap242_components_use_current_products_names_and_nested_transforms": (
            step_ap242_component_identity_ok
        ),
        "curved_mesh_export_uses_current_edges_chord_and_surface_labels": (
            curved_mesh_export_identity_ok
        ),
        "component_closure_diagnosis_verified": summary.get("diagnosis_gate_status") == "ok"
        and summary.get("diagnosis") == "component_solid_closure_loss"
        and summary.get("solver_ready") is False,
        "four_dominant_timing_stages_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    warnings = [] if replay_identity_present else ["replay_identity_not_recorded"]
    return {
        "policy": "build123d_jointed_assembly_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "warnings": warnings,
        "joint_names": sorted(joint_names),
        "connection_endpoints": sorted(connection_endpoints),
        "notes": [
            "Build123d joint semantics are source metadata; persist their names and connection graph beside a neutral CAD export.",
            "A display stub is acceptable for headless replay, but geometry and joint-building statements must remain unchanged.",
            "Record component-level external entities so a zero-volume body cannot hide behind a surviving assembly member.",
        ],
    }
