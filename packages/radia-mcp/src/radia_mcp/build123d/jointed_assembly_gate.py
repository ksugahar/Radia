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


def _valid_sha256(value: object) -> bool:
    digest = str(value or "").lower()
    return len(digest) == 64 and all(
        character in "0123456789abcdef" for character in digest
    )


def _step_import_metadata_topology_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("import_generation", "")).strip()
    labels = [str(item).strip() for item in value.get("shape_labels", [])]
    decoded_labels = [
        str(item).strip() for item in value.get("decoded_shape_labels", [])
    ]
    colors = value.get("shape_colors_rgb", [])
    decoded_colors = value.get("decoded_shape_colors_rgb", [])
    if not isinstance(colors, list) or not isinstance(decoded_colors, list):
        return False
    try:
        colors = [[int(channel) for channel in row] for row in colors]
        decoded_colors = [
            [int(channel) for channel in row] for row in decoded_colors
        ]
        unit_scale = float(value.get("unit_scale_to_m"))
        decoded_unit_scale = float(value.get("decoded_unit_scale_to_m"))
    except (TypeError, ValueError):
        return False
    source_digest = str(value.get("source_content_sha256", "")).lower()
    topology_digest = str(value.get("brep_topology_sha256", "")).lower()
    unit = str(value.get("length_unit", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "source_content_import_generation",
                "label_import_generation",
                "color_import_generation",
                "unit_import_generation",
                "topology_import_generation",
            )
        )
        and _valid_sha256(source_digest)
        and value.get("imported_source_content_sha256") == source_digest
        and bool(labels)
        and all(labels)
        and len(set(labels)) == len(labels)
        and decoded_labels == labels
        and len(colors) == len(labels)
        and all(len(row) == 3 for row in colors)
        and all(0 <= channel <= 255 for row in colors for channel in row)
        and decoded_colors == colors
        and unit in {"m", "cm", "mm"}
        and value.get("decoded_length_unit") == unit
        and math.isfinite(unit_scale)
        and unit_scale > 0.0
        and math.isclose(decoded_unit_scale, unit_scale, rel_tol=0.0, abs_tol=1.0e-18)
        and _valid_sha256(topology_digest)
        and value.get("decoded_brep_topology_sha256") == topology_digest
    )


def _mesh_export_facet_normal_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("mesh_export_generation", "")).strip()
    try:
        facet_ids = [int(item) for item in value.get("facet_ids", [])]
        exported_facet_ids = [
            int(item) for item in value.get("exported_facet_ids", [])
        ]
        normals = [
            [float(component) for component in row]
            for row in value.get("facet_normals", [])
        ]
        exported_normals = [
            [float(component) for component in row]
            for row in value.get("exported_facet_normals", [])
        ]
        chord = float(value.get("chord_tolerance"))
        exported_chord = float(value.get("exported_chord_tolerance"))
        angle = float(value.get("angular_tolerance_deg"))
        exported_angle = float(value.get("exported_angular_tolerance_deg"))
    except (TypeError, ValueError):
        return False
    shape_digest = str(value.get("source_shape_sha256", "")).lower()
    topology_digest = str(value.get("facet_topology_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "shape_mesh_export_generation",
                "facet_mesh_export_generation",
                "normal_mesh_export_generation",
                "tolerance_mesh_export_generation",
            )
        )
        and _valid_sha256(shape_digest)
        and value.get("exported_source_shape_sha256") == shape_digest
        and bool(facet_ids)
        and all(item > 0 for item in facet_ids)
        and len(set(facet_ids)) == len(facet_ids)
        and exported_facet_ids == facet_ids
        and len(normals) == len(facet_ids)
        and all(len(row) == 3 for row in normals)
        and all(math.isfinite(component) for row in normals for component in row)
        and all(
            math.isclose(
                sum(component * component for component in row),
                1.0,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
            for row in normals
        )
        and exported_normals == normals
        and math.isfinite(chord)
        and chord > 0.0
        and math.isclose(exported_chord, chord, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isfinite(angle)
        and 0.0 < angle <= 180.0
        and math.isclose(exported_angle, angle, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(topology_digest)
        and value.get("exported_facet_topology_sha256") == topology_digest
    )


def _brep_step_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("roundtrip_generation", "")).strip()
    source_format = str(value.get("source_format", "")).strip()
    orientation = str(value.get("shell_orientation", "")).strip().lower()
    try:
        linear_tolerance = float(value.get("linear_tolerance"))
        decoded_linear_tolerance = float(value.get("decoded_linear_tolerance"))
        angular_tolerance = float(value.get("angular_tolerance_deg"))
        decoded_angular_tolerance = float(value.get("decoded_angular_tolerance_deg"))
        volume = float(value.get("volume"))
        decoded_volume = float(value.get("decoded_volume"))
    except (TypeError, ValueError):
        return False
    shape_digest = str(value.get("source_shape_sha256", "")).lower()
    topology_digest = str(value.get("topology_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "source_roundtrip_generation",
                "tolerance_roundtrip_generation",
                "orientation_roundtrip_generation",
                "volume_roundtrip_generation",
                "topology_roundtrip_generation",
                "result_roundtrip_generation",
            )
        )
        and bool(source_format)
        and value.get("decoded_source_format") == source_format
        and math.isfinite(linear_tolerance)
        and linear_tolerance > 0.0
        and math.isclose(
            decoded_linear_tolerance,
            linear_tolerance,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        )
        and math.isfinite(angular_tolerance)
        and 0.0 < angular_tolerance <= 180.0
        and math.isclose(
            decoded_angular_tolerance,
            angular_tolerance,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and orientation == "outward"
        and str(value.get("decoded_shell_orientation", "")).strip().lower()
        == orientation
        and math.isfinite(volume)
        and volume > 0.0
        and math.isclose(decoded_volume, volume, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and _valid_sha256(shape_digest)
        and value.get("decoded_source_shape_sha256") == shape_digest
        and _valid_sha256(topology_digest)
        and value.get("decoded_topology_sha256") == topology_digest
    )


def _fresh_subprocess_result_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("run_generation", "")).strip()
    script_digest = str(value.get("source_script_sha256", "")).lower()
    output_digest = str(value.get("output_shape_sha256", "")).lower()
    log_digest = str(value.get("process_log_sha256", "")).lower()
    try:
        remaining_processes = int(value.get("owned_process_count_after"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "process_run_generation",
                "interpreter_run_generation",
                "cache_run_generation",
                "temporary_output_run_generation",
                "result_run_generation",
            )
        )
        and value.get("fresh_interpreter") is True
        and value.get("timed_out") is False
        and value.get("exception_raised") is False
        and value.get("module_cache_preloaded") is False
        and value.get("temporary_directory_unique") is True
        and remaining_processes == 0
        and _valid_sha256(script_digest)
        and value.get("executed_source_script_sha256") == script_digest
        and _valid_sha256(output_digest)
        and value.get("accepted_output_shape_sha256") == output_digest
        and _valid_sha256(log_digest)
        and value.get("accepted_process_log_sha256") == log_digest
    )


def _step_label_color_unit_hierarchy_shape_roundtrip_identity_ok(
    value: object,
) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("roundtrip_generation", "")).strip()
    labels = [str(item).strip() for item in value.get("part_labels", [])]
    decoded_labels = [
        str(item).strip() for item in value.get("decoded_part_labels", [])
    ]
    colors = value.get("part_colors_rgb", [])
    decoded_colors = value.get("decoded_part_colors_rgb", [])
    hierarchy = value.get("assembly_hierarchy", [])
    decoded_hierarchy = value.get("decoded_assembly_hierarchy", [])
    if not all(
        isinstance(rows, list)
        for rows in (colors, decoded_colors, hierarchy, decoded_hierarchy)
    ):
        return False
    try:
        colors = [[float(channel) for channel in row] for row in colors]
        decoded_colors = [
            [float(channel) for channel in row] for row in decoded_colors
        ]
        hierarchy = [[str(item).strip() for item in row] for row in hierarchy]
        decoded_hierarchy = [
            [str(item).strip() for item in row] for row in decoded_hierarchy
        ]
    except (TypeError, ValueError):
        return False
    shape_digests = [
        str(item).lower() for item in value.get("part_shape_sha256", [])
    ]
    decoded_shape_digests = [
        str(item).lower() for item in value.get("decoded_part_shape_sha256", [])
    ]
    unit = str(value.get("length_unit", "")).strip()
    step_digest = str(value.get("step_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "label_roundtrip_generation",
                "color_roundtrip_generation",
                "unit_roundtrip_generation",
                "hierarchy_roundtrip_generation",
                "shape_roundtrip_generation",
                "result_roundtrip_generation",
            )
        )
        and bool(labels)
        and all(labels)
        and len(set(labels)) == len(labels)
        and decoded_labels == labels
        and len(colors) == len(labels)
        and all(len(row) == 3 for row in colors)
        and all(
            math.isfinite(channel) and 0.0 <= channel <= 1.0
            for row in colors
            for channel in row
        )
        and decoded_colors == colors
        and unit in {"m", "cm", "mm"}
        and value.get("decoded_length_unit") == unit
        and len(hierarchy) == len(labels)
        and all(len(row) == 2 and all(row) for row in hierarchy)
        and [row[1] for row in hierarchy] == labels
        and decoded_hierarchy == hierarchy
        and len(shape_digests) == len(labels)
        and all(_valid_sha256(digest) for digest in shape_digests)
        and decoded_shape_digests == shape_digests
        and _valid_sha256(step_digest)
        and value.get("decoded_step_sha256") == step_digest
    )


def _occ_version_tolerance_tessellation_cache_build_identity_ok(
    value: object,
) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("build_generation", "")).strip()
    occ_version = str(value.get("occ_version", "")).strip()
    try:
        tolerance = float(value.get("linear_tolerance"))
        result_tolerance = float(value.get("result_linear_tolerance"))
        linear_deflection = float(value.get("tessellation_linear_deflection"))
        result_linear_deflection = float(
            value.get("result_tessellation_linear_deflection")
        )
        angular_deflection = float(
            value.get("tessellation_angular_deflection_rad")
        )
        result_angular_deflection = float(
            value.get("result_tessellation_angular_deflection_rad")
        )
    except (TypeError, ValueError):
        return False
    cache_digest = str(value.get("module_cache_fingerprint_sha256", "")).lower()
    build_digest = str(value.get("build_fingerprint_sha256", "")).lower()
    tessellation_digest = str(value.get("tessellation_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "occ_build_generation",
                "tolerance_build_generation",
                "tessellation_build_generation",
                "cache_build_generation",
                "fingerprint_build_generation",
                "result_build_generation",
            )
        )
        and bool(occ_version)
        and value.get("result_occ_version") == occ_version
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and math.isclose(result_tolerance, tolerance, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isfinite(linear_deflection)
        and linear_deflection > 0.0
        and math.isclose(
            result_linear_deflection,
            linear_deflection,
            rel_tol=0.0,
            abs_tol=1.0e-18,
        )
        and math.isfinite(angular_deflection)
        and 0.0 < angular_deflection <= math.pi
        and math.isclose(
            result_angular_deflection,
            angular_deflection,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        and _valid_sha256(cache_digest)
        and value.get("result_module_cache_fingerprint_sha256") == cache_digest
        and _valid_sha256(build_digest)
        and value.get("result_build_fingerprint_sha256") == build_digest
        and _valid_sha256(tessellation_digest)
        and value.get("result_tessellation_sha256") == tessellation_digest
    )


def _stl_tessellation_component_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("stl_generation", "")).strip()
    try:
        linear = float(value.get("linear_deflection_m"))
        decoded_linear = float(value.get("decoded_linear_deflection_m"))
        angular = float(value.get("angular_deflection_rad"))
        decoded_angular = float(value.get("decoded_angular_deflection_rad"))
        triangles = int(value.get("triangle_count"))
        decoded_triangles = int(value.get("decoded_triangle_count"))
    except (TypeError, ValueError):
        return False
    components = [str(item).strip() for item in value.get("component_ids", [])]
    decoded_components = [str(item).strip() for item in value.get("decoded_component_ids", [])]
    source_digest = str(value.get("source_shape_sha256", "")).lower()
    stl_digest = str(value.get("stl_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "tolerance_stl_generation", "triangle_stl_generation", "normal_stl_generation",
            "component_stl_generation", "result_stl_generation"))
        and math.isfinite(linear) and linear > 0.0 and decoded_linear == linear
        and math.isfinite(angular) and 0.0 < angular <= math.pi and decoded_angular == angular
        and triangles > 0 and decoded_triangles == triangles
        and value.get("normal_orientation") == "outward"
        and value.get("decoded_normal_orientation") == "outward"
        and bool(components) and all(components) and len(set(components)) == len(components)
        and decoded_components == components
        and _valid_sha256(source_digest) and value.get("tessellated_source_shape_sha256") == source_digest
        and _valid_sha256(stl_digest) and value.get("decoded_stl_sha256") == stl_digest
    )


def _builder_context_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("context_generation", "")).strip()
    stack = [str(item).strip() for item in value.get("context_stack", [])]
    result_stack = [str(item).strip() for item in value.get("result_context_stack", [])]
    parts = [str(item).strip() for item in value.get("part_ids", [])]
    result_parts = [str(item).strip() for item in value.get("result_part_ids", [])]
    try:
        origin = [float(item) for item in value.get("workplane_origin_m", [])]
        result_origin = [float(item) for item in value.get("result_workplane_origin_m", [])]
        normal = [float(item) for item in value.get("workplane_normal", [])]
        result_normal = [float(item) for item in value.get("result_workplane_normal", [])]
        transform = [[float(item) for item in row] for row in value.get("local_frame_transform", [])]
        result_transform = [[float(item) for item in row] for row in value.get("result_local_frame_transform", [])]
    except (TypeError, ValueError):
        return False
    cache_digest = str(value.get("builder_cache_sha256", "")).lower()
    result_digest = str(value.get("builder_result_sha256", "")).lower()
    normal_norm = math.sqrt(sum(item * item for item in normal))
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "stack_context_generation", "workplane_context_generation", "frame_context_generation",
            "part_context_generation", "cache_context_generation", "result_context_generation"))
        and bool(stack) and all(stack) and result_stack == stack
        and len(origin) == 3 and all(math.isfinite(item) for item in origin) and result_origin == origin
        and len(normal) == 3 and all(math.isfinite(item) for item in normal)
        and math.isclose(normal_norm, 1.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_normal == normal
        and len(transform) == 3 and all(len(row) == 4 for row in transform)
        and all(math.isfinite(item) for row in transform for item in row)
        and result_transform == transform
        and bool(parts) and all(parts) and len(set(parts)) == len(parts) and result_parts == parts
        and _valid_sha256(cache_digest) and value.get("result_builder_cache_sha256") == cache_digest
        and _valid_sha256(result_digest) and value.get("result_builder_result_sha256") == result_digest
    )


def _step_import_hierarchy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation", "")).strip()
    hierarchy = [tuple(str(item).strip() for item in row) for row in value.get("product_hierarchy", [])]
    decoded_hierarchy = [tuple(str(item).strip() for item in row) for row in value.get("decoded_product_hierarchy", [])]
    placements = [tuple(row) for row in value.get("placements", [])]
    decoded_placements = [tuple(row) for row in value.get("decoded_placements", [])]
    colors = [tuple(row) for row in value.get("colors_rgb", [])]
    decoded_colors = [tuple(row) for row in value.get("decoded_colors_rgb", [])]
    shape_ids = [str(item).strip() for item in value.get("shape_ids", [])]
    decoded_shape_ids = [str(item).strip() for item in value.get("decoded_shape_ids", [])]
    step_digest = str(value.get("step_sha256", "")).lower()
    map_digest = str(value.get("shape_map_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "unit_step_generation", "hierarchy_step_generation", "placement_step_generation",
            "color_step_generation", "shape_step_generation", "result_step_generation"))
        and value.get("length_unit") in {"m", "cm", "mm"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and bool(hierarchy) and all(len(row) == 2 and all(row) for row in hierarchy)
        and decoded_hierarchy == hierarchy
        and len(placements) == len(shape_ids)
        and all(len(row) == 4 and str(row[0]).strip() in shape_ids for row in placements)
        and all(math.isfinite(float(item)) for row in placements for item in row[1:])
        and decoded_placements == placements
        and len(colors) == len(shape_ids)
        and all(len(row) == 4 and str(row[0]).strip() in shape_ids for row in colors)
        and all(0.0 <= float(item) <= 1.0 for row in colors for item in row[1:])
        and decoded_colors == colors
        and bool(shape_ids) and all(shape_ids) and len(set(shape_ids)) == len(shape_ids)
        and decoded_shape_ids == shape_ids
        and _valid_sha256(step_digest) and value.get("decoded_step_sha256") == step_digest
        and _valid_sha256(map_digest) and value.get("decoded_shape_map_sha256") == map_digest
    )


def _brep_cache_generation_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("brep_generation", "")).strip()
    kernel = str(value.get("kernel_version", "")).strip()
    shape_generation = str(value.get("shape_generation", "")).strip()
    try:
        tolerance = float(value.get("modeling_tolerance_m"))
        decoded_tolerance = float(value.get("decoded_modeling_tolerance_m"))
        location = [[float(item) for item in row] for row in value.get("location_transform", [])]
        decoded_location = [[float(item) for item in row] for row in value.get("decoded_location_transform", [])]
    except (TypeError, ValueError):
        return False
    shape_digest = str(value.get("shape_sha256", "")).lower()
    cache_digest = str(value.get("brep_cache_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "kernel_brep_generation", "tolerance_brep_generation", "location_brep_generation",
            "shape_brep_generation", "cache_brep_generation", "result_brep_generation"))
        and kernel.startswith("OCCT-") and value.get("decoded_kernel_version") == kernel
        and math.isfinite(tolerance) and tolerance > 0.0 and decoded_tolerance == tolerance
        and len(location) == 3 and all(len(row) == 4 for row in location)
        and all(math.isfinite(item) for row in location for item in row)
        and decoded_location == location
        and bool(shape_generation) and value.get("decoded_shape_generation") == shape_generation
        and _valid_sha256(shape_digest) and value.get("decoded_shape_sha256") == shape_digest
        and _valid_sha256(cache_digest) and value.get("decoded_brep_cache_sha256") == cache_digest
    )


def _step_ap242_import_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation", "")).strip()
    context = str(value.get("representation_context", "")).strip()
    products = [tuple(str(item).strip() for item in row) for row in value.get("product_uuid_map", [])]
    decoded_products = [tuple(str(item).strip() for item in row) for row in value.get("decoded_product_uuid_map", [])]
    colors = [tuple(row) for row in value.get("colors_rgb", [])]
    decoded_colors = [tuple(row) for row in value.get("decoded_colors_rgb", [])]
    placements = [tuple(row) for row in value.get("placements", [])]
    decoded_placements = [tuple(row) for row in value.get("decoded_placements", [])]
    owners = [tuple(str(item).strip() for item in row) for row in value.get("shape_owner_map", [])]
    decoded_owners = [tuple(str(item).strip() for item in row) for row in value.get("decoded_shape_owner_map", [])]
    names = [row[0] for row in products if len(row) == 2]
    uuids = [row[1] for row in products if len(row) == 2]
    step_digest = str(value.get("step_sha256", "")).lower()
    map_digest = str(value.get("shape_map_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "context_step_generation", "product_step_generation",
            "unit_step_generation", "color_step_generation",
            "placement_step_generation", "shape_step_generation",
            "file_step_generation", "result_step_generation"))
        and context == "mechanical_design_3d"
        and value.get("decoded_representation_context") == context
        and bool(products) and all(len(row) == 2 and all(row) for row in products)
        and len(set(names)) == len(names) and len(set(uuids)) == len(uuids)
        and all(len(uuid) == 36 and uuid.count("-") == 4 for uuid in uuids)
        and decoded_products == products
        and value.get("length_unit") in {"m", "cm", "mm"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and len(colors) == len(products)
        and all(len(row) == 4 and str(row[0]).strip() in names for row in colors)
        and all(0.0 <= float(item) <= 1.0 for row in colors for item in row[1:])
        and decoded_colors == colors
        and len(placements) == len(products)
        and all(len(row) == 4 and str(row[0]).strip() in names for row in placements)
        and all(math.isfinite(float(item)) for row in placements for item in row[1:])
        and decoded_placements == placements
        and bool(owners) and all(len(row) == 2 and all(row) and row[1] in names for row in owners)
        and len({row[0] for row in owners}) == len(owners)
        and decoded_owners == owners
        and _valid_sha256(step_digest) and value.get("decoded_step_sha256") == step_digest
        and _valid_sha256(map_digest) and value.get("decoded_shape_map_sha256") == map_digest
    )


def _occt_shape_cache_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("cache_generation", "")).strip()
    kernel = str(value.get("kernel_version", "")).strip()
    try:
        location = [[float(item) for item in row] for row in value.get("location_transform", [])]
        decoded_location = [[float(item) for item in row] for row in value.get("decoded_location_transform", [])]
        tolerance = float(value.get("modeling_tolerance_m"))
        decoded_tolerance = float(value.get("decoded_modeling_tolerance_m"))
        triangulation = [float(item) for item in value.get("triangulation_parameters", [])]
        decoded_triangulation = [float(item) for item in value.get("decoded_triangulation_parameters", [])]
    except (TypeError, ValueError):
        return False
    shape_digest = str(value.get("shape_sha256", "")).lower()
    triangulation_digest = str(value.get("triangulation_sha256", "")).lower()
    serialization_digest = str(value.get("serialization_sha256", "")).lower()
    cache_digest = str(value.get("cache_sha256", "")).lower()
    serialization_format = str(value.get("serialization_format", "")).strip()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "kernel_cache_generation", "shape_cache_generation",
            "location_cache_generation", "tolerance_cache_generation",
            "triangulation_cache_generation", "serialization_cache_generation",
            "result_cache_generation"))
        and kernel.startswith("OCCT-") and value.get("decoded_kernel_version") == kernel
        and _valid_sha256(shape_digest) and value.get("decoded_shape_sha256") == shape_digest
        and len(location) == 3 and all(len(row) == 4 for row in location)
        and all(math.isfinite(item) for row in location for item in row)
        and decoded_location == location
        and math.isfinite(tolerance) and tolerance > 0.0 and decoded_tolerance == tolerance
        and len(triangulation) == 3
        and triangulation[0] > 0.0 and 0.0 < triangulation[1] <= math.pi
        and triangulation[2] in {0.0, 1.0}
        and decoded_triangulation == triangulation
        and _valid_sha256(triangulation_digest)
        and value.get("decoded_triangulation_sha256") == triangulation_digest
        and serialization_format == "brep-binary-v1"
        and value.get("decoded_serialization_format") == serialization_format
        and _valid_sha256(serialization_digest)
        and value.get("decoded_serialization_sha256") == serialization_digest
        and _valid_sha256(cache_digest) and value.get("decoded_cache_sha256") == cache_digest
    )


def _dxf_face_reconstruction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("dxf_generation", "")).strip()
    arcs = [tuple(row) for row in value.get("arc_parameters", [])]
    decoded_arcs = [tuple(row) for row in value.get("decoded_arc_parameters", [])]
    splines = [tuple(row) for row in value.get("spline_parameters", [])]
    decoded_splines = [tuple(row) for row in value.get("decoded_spline_parameters", [])]
    layers = [tuple(str(item).strip() for item in row) for row in value.get("entity_layer_map", [])]
    decoded_layers = [
        tuple(str(item).strip() for item in row)
        for row in value.get("decoded_entity_layer_map", [])
    ]
    try:
        plane = [[float(item) for item in row] for row in value.get("workplane_matrix", [])]
        decoded_plane = [
            [float(item) for item in row]
            for row in value.get("decoded_workplane_matrix", [])
        ]
        area = float(value.get("face_area_mm2"))
        decoded_area = float(value.get("decoded_face_area_mm2"))
        arcs_numeric = [tuple(float(item) for item in row[1:]) for row in arcs]
        spline_degree_counts = [
            (int(row[1]), int(row[2])) for row in splines if len(row) == 4
        ]
    except (TypeError, ValueError):
        return False
    arc_names = [str(row[0]).strip() for row in arcs if len(row) == 6]
    spline_names = [str(row[0]).strip() for row in splines if len(row) == 4]
    entity_names = arc_names + spline_names
    dxf_digest = str(value.get("dxf_sha256", "")).lower()
    face_digest = str(value.get("face_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "arc_dxf_generation", "spline_dxf_generation", "layer_dxf_generation",
            "plane_dxf_generation", "unit_dxf_generation", "wire_dxf_generation",
            "face_dxf_generation", "result_dxf_generation"))
        and bool(arcs) and all(len(row) == 6 for row in arcs)
        and len(arc_names) == len(arcs) and all(arc_names)
        and len(set(arc_names)) == len(arc_names)
        and all(
            all(math.isfinite(item) for item in numeric)
            and numeric[2] > 0.0 and numeric[3] != numeric[4]
            for numeric in arcs_numeric
        )
        and decoded_arcs == arcs
        and bool(splines) and all(len(row) == 4 for row in splines)
        and len(spline_names) == len(splines) and all(spline_names)
        and len(set(spline_names)) == len(spline_names)
        and all(degree >= 1 and count > degree for degree, count in spline_degree_counts)
        and all(_valid_sha256(str(row[3]).lower()) for row in splines)
        and decoded_splines == splines
        and len(set(entity_names)) == len(entity_names)
        and len(layers) == len(entity_names)
        and {row[0] for row in layers if len(row) == 2} == set(entity_names)
        and all(len(row) == 2 and all(row) for row in layers)
        and decoded_layers == layers
        and len(plane) == 4
        and all(len(row) == 4 and all(math.isfinite(item) for item in row) for row in plane)
        and plane[3] == [0.0, 0.0, 0.0, 1.0]
        and decoded_plane == plane
        and value.get("length_unit") in {"m", "cm", "mm", "in"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and value.get("wire_closed") is True
        and value.get("decoded_wire_closed") is True
        and value.get("wire_orientation") == "counterclockwise"
        and value.get("decoded_wire_orientation") == value.get("wire_orientation")
        and math.isfinite(area) and area > 0.0 and decoded_area == area
        and _valid_sha256(dxf_digest) and value.get("decoded_dxf_sha256") == dxf_digest
        and _valid_sha256(face_digest) and value.get("decoded_face_sha256") == face_digest
    )


def _three_mf_reconstruction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("three_mf_generation", "")).strip()
    names = [str(item).strip() for item in value.get("component_names", [])]
    decoded_names = [str(item).strip() for item in value.get("decoded_component_names", [])]
    transforms = [tuple(str(item).strip() for item in row) for row in value.get("component_transform_sha256", [])]
    decoded_transforms = [
        tuple(str(item).strip() for item in row)
        for row in value.get("decoded_component_transform_sha256", [])
    ]
    try:
        materials = [(str(row[0]).strip(), int(row[1])) for row in value.get("material_id_map", [])]
        decoded_materials = [(str(row[0]).strip(), int(row[1])) for row in value.get("decoded_material_id_map", [])]
        volumes = [(str(row[0]).strip(), float(row[1])) for row in value.get("component_volumes_mm3", [])]
        decoded_volumes = [(str(row[0]).strip(), float(row[1])) for row in value.get("decoded_component_volumes_mm3", [])]
        total_volume = float(value.get("total_volume_mm3"))
        decoded_total_volume = float(value.get("decoded_total_volume_mm3"))
    except (TypeError, ValueError, IndexError):
        return False
    watertight = [str(item).strip() for item in value.get("watertight_component_names", [])]
    decoded_watertight = [
        str(item).strip() for item in value.get("decoded_watertight_component_names", [])
    ]
    file_digest = str(value.get("three_mf_sha256", "")).lower()
    mesh_digest = str(value.get("triangle_mesh_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "component_three_mf_generation", "transform_three_mf_generation",
            "winding_three_mf_generation", "material_three_mf_generation",
            "watertight_three_mf_generation", "volume_three_mf_generation",
            "file_three_mf_generation", "result_three_mf_generation"))
        and bool(names) and all(names) and len(set(names)) == len(names)
        and decoded_names == names
        and len(transforms) == len(names)
        and [row[0] for row in transforms if len(row) == 2] == names
        and all(len(row) == 2 and _valid_sha256(row[1].lower()) for row in transforms)
        and decoded_transforms == transforms
        and value.get("triangle_winding") == "outward_counterclockwise"
        and value.get("decoded_triangle_winding") == value.get("triangle_winding")
        and len(materials) == len(names) and [row[0] for row in materials] == names
        and all(row[1] > 0 for row in materials) and decoded_materials == materials
        and watertight == names and decoded_watertight == watertight
        and len(volumes) == len(names) and [row[0] for row in volumes] == names
        and all(math.isfinite(row[1]) and row[1] > 0.0 for row in volumes)
        and decoded_volumes == volumes
        and math.isfinite(total_volume) and total_volume > 0.0
        and math.isclose(total_volume, sum(row[1] for row in volumes), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and decoded_total_volume == total_volume
        and value.get("length_unit") in {"m", "cm", "mm", "in"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and _valid_sha256(file_digest) and value.get("decoded_three_mf_sha256") == file_digest
        and _valid_sha256(mesh_digest) and value.get("decoded_triangle_mesh_sha256") == mesh_digest
    )


def _step_assembly_reconstruction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation", "")).strip()
    names = [str(item).strip() for item in value.get("product_names", [])]
    decoded_names = [
        str(item).strip() for item in value.get("decoded_product_names", [])
    ]
    instances = [str(item).strip() for item in value.get("instance_order", [])]
    decoded_instances = [
        str(item).strip() for item in value.get("decoded_instance_order", [])
    ]
    transforms = [
        tuple(str(item).strip() for item in row)
        for row in value.get("instance_transform_sha256", [])
    ]
    decoded_transforms = [
        tuple(str(item).strip() for item in row)
        for row in value.get("decoded_instance_transform_sha256", [])
    ]
    hierarchy = [
        tuple(str(item).strip() for item in row)
        for row in value.get("assembly_hierarchy", [])
    ]
    decoded_hierarchy = [
        tuple(str(item).strip() for item in row)
        for row in value.get("decoded_assembly_hierarchy", [])
    ]
    try:
        colors = [
            (str(row[0]).strip(), *(float(item) for item in row[1:]))
            for row in value.get("product_colors_rgb", [])
        ]
        decoded_colors = [
            (str(row[0]).strip(), *(float(item) for item in row[1:]))
            for row in value.get("decoded_product_colors_rgb", [])
        ]
    except (TypeError, ValueError, IndexError):
        return False
    file_digest = str(value.get("step_sha256", "")).lower()
    shape_digest = str(value.get("assembly_shape_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "product_step_generation", "color_step_generation",
            "instance_step_generation", "transform_step_generation",
            "unit_step_generation", "hierarchy_step_generation",
            "file_step_generation", "result_step_generation"))
        and len(names) >= 2 and all(names) and len(set(names)) == len(names)
        and decoded_names == names
        and bool(colors)
        and all(len(row) == 4 and row[0] in names and all(math.isfinite(item) and 0.0 <= item <= 1.0 for item in row[1:]) for row in colors)
        and len({row[0] for row in colors}) == len(colors)
        and decoded_colors == colors
        and bool(instances) and all(instances) and len(set(instances)) == len(instances)
        and decoded_instances == instances
        and len(transforms) == len(instances)
        and [row[0] for row in transforms if len(row) == 2] == instances
        and all(len(row) == 2 and _valid_sha256(row[1].lower()) for row in transforms)
        and decoded_transforms == transforms
        and value.get("length_unit") in {"m", "cm", "mm", "in"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and len(hierarchy) == len(instances)
        and all(len(row) == 2 and row[0] in names and row[1] in instances for row in hierarchy)
        and [row[1] for row in hierarchy] == instances
        and decoded_hierarchy == hierarchy
        and _valid_sha256(file_digest) and value.get("decoded_step_sha256") == file_digest
        and _valid_sha256(shape_digest)
        and value.get("decoded_assembly_shape_sha256") == shape_digest
    )


def _stl_solid_reconstruction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("stl_generation", "")).strip()
    try:
        normals = [
            [float(item) for item in row] for row in value.get("facet_normals", [])
        ]
        decoded_normals = [
            [float(item) for item in row]
            for row in value.get("decoded_facet_normals", [])
        ]
        unmatched = int(value.get("unmatched_edge_count"))
        decoded_unmatched = int(value.get("decoded_unmatched_edge_count"))
        tolerance = float(value.get("merge_tolerance_m"))
        decoded_tolerance = float(value.get("decoded_merge_tolerance_m"))
        volume = float(value.get("signed_volume_m3"))
        decoded_volume = float(value.get("decoded_signed_volume_m3"))
    except (TypeError, ValueError):
        return False
    file_digest = str(value.get("stl_sha256", "")).lower()
    solid_digest = str(value.get("stl_solid_sha256", "")).lower()
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "normal_stl_generation", "winding_stl_generation",
            "watertight_stl_generation", "tolerance_stl_generation",
            "volume_stl_generation", "unit_stl_generation",
            "file_stl_generation", "result_stl_generation"))
        and bool(normals)
        and all(
            len(row) == 3
            and all(math.isfinite(item) for item in row)
            and math.isclose(sum(item * item for item in row), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for row in normals
        )
        and decoded_normals == normals
        and value.get("triangle_winding") == "outward_counterclockwise"
        and value.get("decoded_triangle_winding") == value.get("triangle_winding")
        and unmatched == 0 and decoded_unmatched == unmatched
        and math.isfinite(tolerance) and tolerance > 0.0
        and decoded_tolerance == tolerance
        and math.isfinite(volume) and volume > 0.0 and decoded_volume == volume
        and value.get("length_unit") in {"m", "cm", "mm", "in"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and _valid_sha256(file_digest) and value.get("decoded_stl_sha256") == file_digest
        and _valid_sha256(solid_digest)
        and value.get("decoded_stl_solid_sha256") == solid_digest
    )


def _step_ap242_context_owner_mass_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation") or "")
    try:
        contexts = [[str(item) for item in row] for row in value.get("representation_contexts", [])]
        decoded_contexts = [[str(item) for item in row] for row in value.get("decoded_representation_contexts", [])]
        owners = [[str(item) for item in row] for row in value.get("part_owners", [])]
        decoded_owners = [[str(item) for item in row] for row in value.get("decoded_part_owners", [])]
        occurrence_ids = [str(item) for item in value.get("occurrence_ids", [])]
        decoded_occurrence_ids = [str(item) for item in value.get("decoded_occurrence_ids", [])]
        transforms = [[str(item) for item in row] for row in value.get("occurrence_transform_sha256", [])]
        decoded_transforms = [[str(item) for item in row] for row in value.get("decoded_occurrence_transform_sha256", [])]
        mass_rows = [list(row) for row in value.get("occurrence_mass_properties", [])]
        decoded_mass_rows = [list(row) for row in value.get("decoded_occurrence_mass_properties", [])]
    except (TypeError, ValueError):
        return False
    valid_mass_rows = True
    for row in mass_rows:
        if len(row) != 5 or str(row[0]) not in occurrence_ids:
            valid_mass_rows = False
            break
        try:
            numbers = [float(item) for item in row[1:]]
        except (TypeError, ValueError):
            valid_mass_rows = False
            break
        if numbers[0] <= 0.0 or any(not math.isfinite(item) for item in numbers):
            valid_mass_rows = False
            break
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "schema_step_generation", "context_step_generation",
            "owner_step_generation", "occurrence_step_generation",
            "transform_step_generation", "mass_step_generation",
            "structure_step_generation", "file_step_generation",
            "result_step_generation"))
        and value.get("ap_schema") == "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF"
        and value.get("decoded_ap_schema") == value.get("ap_schema")
        and bool(contexts)
        and all(len(row) == 3 and all(row) for row in contexts)
        and decoded_contexts == contexts
        and bool(owners)
        and all(len(row) == 2 and all(row) for row in owners)
        and decoded_owners == owners
        and bool(occurrence_ids)
        and len(set(occurrence_ids)) == len(occurrence_ids)
        and decoded_occurrence_ids == occurrence_ids
        and len(transforms) == len(occurrence_ids)
        and all(
            len(row) == 2
            and row[0] in occurrence_ids
            and _valid_sha256(row[1])
            for row in transforms
        )
        and {row[0] for row in transforms} == set(occurrence_ids)
        and decoded_transforms == transforms
        and len(mass_rows) == len(occurrence_ids)
        and valid_mass_rows
        and {str(row[0]) for row in mass_rows} == set(occurrence_ids)
        and decoded_mass_rows == mass_rows
        and _valid_sha256(value.get("product_structure_sha256"))
        and value.get("decoded_product_structure_sha256")
        == value.get("product_structure_sha256")
        and _valid_sha256(value.get("step_file_sha256"))
        and value.get("decoded_step_file_sha256") == value.get("step_file_sha256")
    )


def _stl_brep_tessellation_error_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("stl_generation") or "")
    try:
        chord = float(value.get("chord_tolerance_m"))
        decoded_chord = float(value.get("decoded_chord_tolerance_m"))
        angle = float(value.get("angular_tolerance_deg"))
        decoded_angle = float(value.get("decoded_angular_tolerance_deg"))
        facets = int(value.get("facet_count"))
        decoded_facets = int(value.get("decoded_facet_count"))
        components = int(value.get("connected_component_count"))
        decoded_components = int(value.get("decoded_connected_component_count"))
        deviation = float(value.get("maximum_surface_deviation_m"))
        decoded_deviation = float(value.get("decoded_maximum_surface_deviation_m"))
        source_area = float(value.get("source_surface_area_m2"))
        decoded_area = float(value.get("decoded_surface_area_m2"))
        source_volume = float(value.get("source_volume_m3"))
        decoded_volume = float(value.get("decoded_volume_m3"))
        area_tolerance = float(value.get("relative_area_tolerance"))
        volume_tolerance = float(value.get("relative_volume_tolerance"))
    except (TypeError, ValueError):
        return False
    area_error = abs(decoded_area - source_area) / max(abs(source_area), 1.0e-300)
    volume_error = abs(decoded_volume - source_volume) / max(abs(source_volume), 1.0e-300)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "source_stl_generation", "tolerance_stl_generation",
            "facet_stl_generation", "component_stl_generation",
            "deviation_stl_generation", "area_stl_generation",
            "volume_stl_generation", "file_stl_generation",
            "result_stl_generation"))
        and _valid_sha256(value.get("source_brep_sha256"))
        and value.get("decoded_source_brep_sha256") == value.get("source_brep_sha256")
        and math.isfinite(chord) and chord > 0.0 and decoded_chord == chord
        and math.isfinite(angle) and 0.0 < angle < 180.0 and decoded_angle == angle
        and facets > 0 and decoded_facets == facets
        and components == 1 and decoded_components == components
        and math.isfinite(deviation) and 0.0 <= deviation <= chord
        and decoded_deviation == deviation
        and all(math.isfinite(item) and item > 0.0 for item in (source_area, decoded_area, source_volume, decoded_volume))
        and math.isfinite(area_tolerance) and 0.0 < area_tolerance < 1.0
        and math.isfinite(volume_tolerance) and 0.0 < volume_tolerance < 1.0
        and area_error <= area_tolerance
        and volume_error <= volume_tolerance
        and _valid_sha256(value.get("normal_table_sha256"))
        and value.get("decoded_normal_table_sha256") == value.get("normal_table_sha256")
        and _valid_sha256(value.get("stl_file_sha256"))
        and value.get("decoded_stl_file_sha256") == value.get("stl_file_sha256")
    )


def _step_assembly_instance_metadata_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation") or "")
    try:
        instances = [str(item) for item in value.get("instance_ids", [])]
        decoded_instances = [str(item) for item in value.get("decoded_instance_ids", [])]
        transforms = [[str(item) for item in row] for row in value.get("instance_transform_sha256", [])]
        decoded_transforms = [[str(item) for item in row] for row in value.get("decoded_instance_transform_sha256", [])]
        colors = [list(row) for row in value.get("instance_colors_rgb", [])]
        decoded_colors = [list(row) for row in value.get("decoded_instance_colors_rgb", [])]
        materials = [[str(item) for item in row] for row in value.get("material_labels", [])]
        decoded_materials = [[str(item) for item in row] for row in value.get("decoded_material_labels", [])]
        uuids = [[str(item) for item in row] for row in value.get("product_uuids", [])]
        decoded_uuids = [[str(item) for item in row] for row in value.get("decoded_product_uuids", [])]
        components = [[str(item) for item in row] for row in value.get("component_shape_sha256", [])]
        decoded_components = [[str(item) for item in row] for row in value.get("decoded_component_shape_sha256", [])]
        volume = float(value.get("total_volume_m3"))
        decoded_volume = float(value.get("decoded_total_volume_m3"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "instance_generation", "transform_generation", "unit_generation",
            "color_generation", "material_generation", "uuid_generation",
            "component_generation", "volume_generation", "file_generation",
            "result_generation"))
        and bool(instances)
        and len(set(instances)) == len(instances)
        and decoded_instances == instances
        and len(transforms) == len(colors) == len(materials) == len(uuids) == len(components) == len(instances)
        and all(len(row) == 2 and row[0] in instances and _valid_sha256(row[1]) for row in transforms + components)
        and {row[0] for row in transforms} == {row[0] for row in components} == set(instances)
        and decoded_transforms == transforms
        and all(
            len(row) == 4
            and str(row[0]) in instances
            and all(math.isfinite(float(item)) and 0.0 <= float(item) <= 1.0 for item in row[1:])
            for row in colors
        )
        and decoded_colors == colors
        and all(len(row) == 2 and row[0] in instances and row[1] for row in materials + uuids)
        and {row[0] for row in materials} == {row[0] for row in uuids} == set(instances)
        and decoded_materials == materials
        and decoded_uuids == uuids
        and decoded_components == components
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and math.isfinite(volume) and volume > 0.0 and decoded_volume == volume
        and _valid_sha256(value.get("step_file_sha256"))
        and value.get("decoded_step_file_sha256") == value.get("step_file_sha256")
    )


def _stl_repair_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("repair_generation") or "")
    try:
        tolerance = float(value.get("merge_tolerance_m"))
        decoded_tolerance = float(value.get("decoded_merge_tolerance_m"))
        duplicates = int(value.get("duplicate_vertex_count"))
        decoded_duplicates = int(value.get("decoded_duplicate_vertex_count"))
        boundaries = int(value.get("boundary_edge_count"))
        decoded_boundaries = int(value.get("decoded_boundary_edge_count"))
        components = int(value.get("watertight_component_count"))
        decoded_components = int(value.get("decoded_watertight_component_count"))
        volume = float(value.get("repaired_volume_m3"))
        decoded_volume = float(value.get("decoded_repaired_volume_m3"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "tolerance_generation", "normal_generation", "duplicate_generation",
            "boundary_generation", "watertight_generation", "volume_generation",
            "unit_generation", "file_generation", "result_generation"))
        and math.isfinite(tolerance) and tolerance > 0.0 and decoded_tolerance == tolerance
        and _valid_sha256(value.get("facet_normal_sha256"))
        and value.get("decoded_facet_normal_sha256") == value.get("facet_normal_sha256")
        and duplicates == 0 and decoded_duplicates == duplicates
        and boundaries == 0 and decoded_boundaries == boundaries
        and components == 1 and decoded_components == components
        and math.isfinite(volume) and volume > 0.0 and decoded_volume == volume
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and all(
            _valid_sha256(value.get(source))
            and value.get(target) == value.get(source)
            for source, target in (
                ("source_stl_sha256", "decoded_source_stl_sha256"),
                ("repaired_stl_sha256", "decoded_repaired_stl_sha256"),
            )
        )
    )


def _brep_semantic_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("brep_generation") or "")
    try:
        pcurves = [[int(row[0]), str(row[1])] for row in value.get("face_pcurve_sha256", [])]
        decoded_pcurves = [
            [int(row[0]), str(row[1])]
            for row in value.get("decoded_face_pcurve_sha256", [])
        ]
        orientations = [
            [int(row[0]), int(row[1])]
            for row in value.get("wire_orientation_signs", [])
        ]
        decoded_orientations = [
            [int(row[0]), int(row[1])]
            for row in value.get("decoded_wire_orientation_signs", [])
        ]
        locations = [
            [str(row[0]), str(row[1])]
            for row in value.get("nested_location_sha256", [])
        ]
        decoded_locations = [
            [str(row[0]), str(row[1])]
            for row in value.get("decoded_nested_location_sha256", [])
        ]
        tolerances = [
            [int(row[0]), float(row[1])]
            for row in value.get("edge_tolerances_m", [])
        ]
        decoded_tolerances = [
            [int(row[0]), float(row[1])]
            for row in value.get("decoded_edge_tolerances_m", [])
        ]
        surfaces = [
            [int(row[0]), str(row[1])] for row in value.get("surface_types", [])
        ]
        decoded_surfaces = [
            [int(row[0]), str(row[1])]
            for row in value.get("decoded_surface_types", [])
        ]
    except (IndexError, TypeError, ValueError):
        return False
    face_ids = {row[0] for row in pcurves}
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "pcurve_generation",
                "wire_generation",
                "location_generation",
                "tolerance_generation",
                "surface_generation",
                "serializer_generation",
                "shape_generation",
                "result_generation",
            )
        )
        and bool(pcurves)
        and len(face_ids) == len(pcurves)
        and all(row[0] > 0 and _valid_sha256(row[1]) for row in pcurves)
        and decoded_pcurves == pcurves
        and bool(orientations)
        and len({row[0] for row in orientations}) == len(orientations)
        and all(row[0] > 0 and row[1] in {-1, 1} for row in orientations)
        and decoded_orientations == orientations
        and bool(locations)
        and len({row[0] for row in locations}) == len(locations)
        and all(row[0] and _valid_sha256(row[1]) for row in locations)
        and decoded_locations == locations
        and bool(tolerances)
        and len({row[0] for row in tolerances}) == len(tolerances)
        and all(row[0] > 0 and math.isfinite(row[1]) and row[1] > 0.0 for row in tolerances)
        and decoded_tolerances == tolerances
        and bool(surfaces)
        and {row[0] for row in surfaces} == face_ids
        and all(row[0] > 0 and row[1] for row in surfaces)
        and decoded_surfaces == surfaces
        and str(value.get("serializer_version") or "").startswith("occt-brep-v")
        and value.get("decoded_serializer_version") == value.get("serializer_version")
        and _valid_sha256(value.get("brep_shape_sha256"))
        and value.get("decoded_brep_shape_sha256") == value.get("brep_shape_sha256")
    )


def _gltf_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("gltf_generation") or "")
    try:
        hierarchy = [
            [str(row[0]), str(row[1])] for row in value.get("node_hierarchy", [])
        ]
        decoded_hierarchy = [
            [str(row[0]), str(row[1])]
            for row in value.get("decoded_node_hierarchy", [])
        ]
        transforms = [
            [str(row[0]), str(row[1])]
            for row in value.get("instance_transform_sha256", [])
        ]
        decoded_transforms = [
            [str(row[0]), str(row[1])]
            for row in value.get("decoded_instance_transform_sha256", [])
        ]
        materials = [
            [str(row[0]), str(row[1])]
            for row in value.get("material_assignments", [])
        ]
        decoded_materials = [
            [str(row[0]), str(row[1])]
            for row in value.get("decoded_material_assignments", [])
        ]
        linear = float(value.get("linear_deflection_m"))
        decoded_linear = float(value.get("decoded_linear_deflection_m"))
        angular = float(value.get("angular_deflection_rad"))
        decoded_angular = float(value.get("decoded_angular_deflection_rad"))
        triangles = int(value.get("triangle_count"))
        decoded_triangles = int(value.get("decoded_triangle_count"))
        volume = float(value.get("enclosed_volume_m3"))
        decoded_volume = float(value.get("decoded_enclosed_volume_m3"))
    except (IndexError, TypeError, ValueError):
        return False
    nodes = {item for pair in hierarchy for item in pair}
    parents = {child: parent for parent, child in hierarchy}
    hierarchy_is_tree = bool(hierarchy) and len(parents) == len(hierarchy)
    if hierarchy_is_tree:
        for node in nodes:
            visited = set()
            current = node
            while current in parents:
                if current in visited:
                    hierarchy_is_tree = False
                    break
                visited.add(current)
                current = parents[current]
            if not hierarchy_is_tree:
                break
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "node_generation",
                "transform_generation",
                "winding_generation",
                "material_generation",
                "unit_generation",
                "tessellation_generation",
                "volume_generation",
                "file_generation",
                "result_generation",
            )
        )
        and hierarchy_is_tree
        and all(parent and child and parent != child for parent, child in hierarchy)
        and decoded_hierarchy == hierarchy
        and bool(transforms)
        and len({row[0] for row in transforms}) == len(transforms)
        and all(row[0] in nodes and _valid_sha256(row[1]) for row in transforms)
        and decoded_transforms == transforms
        and _valid_sha256(value.get("triangle_winding_sha256"))
        and value.get("decoded_triangle_winding_sha256")
        == value.get("triangle_winding_sha256")
        and bool(materials)
        and len({row[0] for row in materials}) == len(materials)
        and all(row[0] in nodes and row[1] for row in materials)
        and decoded_materials == materials
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and math.isfinite(linear)
        and linear > 0.0
        and decoded_linear == linear
        and math.isfinite(angular)
        and 0.0 < angular < math.pi
        and decoded_angular == angular
        and triangles > 0
        and decoded_triangles == triangles
        and math.isfinite(volume)
        and volume > 0.0
        and decoded_volume == volume
        and _valid_sha256(value.get("gltf_file_sha256"))
        and value.get("decoded_gltf_file_sha256") == value.get("gltf_file_sha256")
    )


def _step_semantic_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("step_generation") or "")
    try:
        products = [
            [int(row[0]), str(row[1]).strip()]
            for row in value.get("product_entities", [])
        ]
        decoded_products = [
            [int(row[0]), str(row[1]).strip()]
            for row in value.get("decoded_product_entities", [])
        ]
        colors = [
            [int(row[0]), [float(item) for item in row[1]]]
            for row in value.get("entity_colors_rgb", [])
        ]
        decoded_colors = [
            [int(row[0]), [float(item) for item in row[1]]]
            for row in value.get("decoded_entity_colors_rgb", [])
        ]
        transforms = [
            [int(row[0]), str(row[1])]
            for row in value.get("assembly_transform_sha256", [])
        ]
        decoded_transforms = [
            [int(row[0]), str(row[1])]
            for row in value.get("decoded_assembly_transform_sha256", [])
        ]
        shape_count = int(value.get("shape_count"))
        decoded_shape_count = int(value.get("decoded_shape_count"))
        validity = [
            [int(row[0]), bool(row[1])]
            for row in value.get("solid_validity", [])
        ]
        decoded_validity = [
            [int(row[0]), bool(row[1])]
            for row in value.get("decoded_solid_validity", [])
        ]
    except (IndexError, TypeError, ValueError):
        return False
    product_ids = [row[0] for row in products]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "unit_generation",
                "entity_generation",
                "color_generation",
                "transform_generation",
                "shape_generation",
                "validity_generation",
                "owner_generation",
                "file_generation",
                "result_generation",
            )
        )
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and bool(products)
        and len(set(product_ids)) == len(products)
        and all(row[0] > 0 and row[1] for row in products)
        and decoded_products == products
        and len(colors) == len(products)
        and {row[0] for row in colors} == set(product_ids)
        and all(
            len(row[1]) == 3
            and all(math.isfinite(channel) and 0.0 <= channel <= 1.0 for channel in row[1])
            for row in colors
        )
        and decoded_colors == colors
        and len(transforms) == len(products)
        and {row[0] for row in transforms} == set(product_ids)
        and all(_valid_sha256(row[1]) for row in transforms)
        and decoded_transforms == transforms
        and shape_count == len(products)
        and decoded_shape_count == shape_count
        and len(validity) == shape_count
        and {row[0] for row in validity} == set(product_ids)
        and all(row[1] is True for row in validity)
        and decoded_validity == validity
        and bool(str(value.get("source_owner") or "").strip())
        and value.get("decoded_source_owner") == value.get("source_owner")
        and _valid_sha256(value.get("step_file_sha256"))
        and value.get("decoded_step_file_sha256") == value.get("step_file_sha256")
    )


def _brep_periodic_seam_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("periodic_brep_generation") or "")
    try:
        face_ids = [int(item) for item in value.get("periodic_face_ids", [])]
        decoded_face_ids = [int(item) for item in value.get("decoded_periodic_face_ids", [])]
        seams = [
            [int(row[0]), int(row[1]), int(row[2])]
            for row in value.get("seam_edge_multiplicity", [])
        ]
        decoded_seams = [
            [int(row[0]), int(row[1]), int(row[2])]
            for row in value.get("decoded_seam_edge_multiplicity", [])
        ]
        orientations = [
            [int(row[0]), int(row[1])]
            for row in value.get("face_orientation_signs", [])
        ]
        decoded_orientations = [
            [int(row[0]), int(row[1])]
            for row in value.get("decoded_face_orientation_signs", [])
        ]
        pcurve_deviations = [
            [int(row[0]), float(row[1])]
            for row in value.get("edge_pcurve_max_deviation_m", [])
        ]
        decoded_pcurve_deviations = [
            [int(row[0]), float(row[1])]
            for row in value.get("decoded_edge_pcurve_max_deviation_m", [])
        ]
        vertex_tolerances = [
            [int(row[0]), float(row[1])]
            for row in value.get("vertex_tolerances_m", [])
        ]
        decoded_vertex_tolerances = [
            [int(row[0]), float(row[1])]
            for row in value.get("decoded_vertex_tolerances_m", [])
        ]
        incidence = [
            [int(row[0]), int(row[1])]
            for row in value.get("edge_face_incidence_counts", [])
        ]
        decoded_incidence = [
            [int(row[0]), int(row[1])]
            for row in value.get("decoded_edge_face_incidence_counts", [])
        ]
    except (IndexError, TypeError, ValueError):
        return False
    seam_edge_ids = {row[1] for row in seams}
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "face_generation",
                "seam_generation",
                "orientation_generation",
                "pcurve_generation",
                "tolerance_generation",
                "manifold_generation",
                "serializer_generation",
                "shape_generation",
                "result_generation",
            )
        )
        and bool(face_ids)
        and len(set(face_ids)) == len(face_ids)
        and all(face_id > 0 for face_id in face_ids)
        and decoded_face_ids == face_ids
        and len(seams) == len(face_ids)
        and {row[0] for row in seams} == set(face_ids)
        and len(seam_edge_ids) == len(seams)
        and all(row[1] > 0 and row[2] == 2 for row in seams)
        and decoded_seams == seams
        and len(orientations) == len(face_ids)
        and {row[0] for row in orientations} == set(face_ids)
        and all(row[1] in {-1, 1} for row in orientations)
        and decoded_orientations == orientations
        and {row[0] for row in pcurve_deviations} == seam_edge_ids
        and all(
            math.isfinite(row[1]) and 0.0 <= row[1] <= 1.0e-6
            for row in pcurve_deviations
        )
        and decoded_pcurve_deviations == pcurve_deviations
        and bool(vertex_tolerances)
        and len({row[0] for row in vertex_tolerances}) == len(vertex_tolerances)
        and all(
            row[0] > 0 and math.isfinite(row[1]) and 0.0 < row[1] <= 1.0e-5
            for row in vertex_tolerances
        )
        and decoded_vertex_tolerances == vertex_tolerances
        and {row[0] for row in incidence} == seam_edge_ids
        and all(row[1] == 2 for row in incidence)
        and decoded_incidence == incidence
        and str(value.get("serializer_version") or "").startswith("occt-brep-v")
        and value.get("decoded_serializer_version") == value.get("serializer_version")
        and _valid_sha256(value.get("periodic_brep_shape_sha256"))
        and value.get("decoded_periodic_brep_shape_sha256")
        == value.get("periodic_brep_shape_sha256")
    )


def _step_ap242_pmi_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("ap242_generation") or "")
    try:
        pmi = [[int(row[0]), str(row[1]), str(row[2]), float(row[3])] for row in value.get("pmi_annotations", [])]
        names = [[int(row[0]), str(row[1])] for row in value.get("product_names", [])]
        colors = [[int(row[0]), [float(item) for item in row[1]]] for row in value.get("colors_rgb", [])]
        occurrences = [[int(row[0]), str(row[1])] for row in value.get("occurrence_paths", [])]
        transforms = [[int(row[0]), str(row[1])] for row in value.get("transform_sha256", [])]
        validity = [[int(row[0]), bool(row[1])] for row in value.get("solid_validity", [])]
    except (IndexError, TypeError, ValueError):
        return False
    product_ids = {row[0] for row in names}
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("pmi_generation", "unit_generation", "name_generation", "color_generation", "occurrence_generation", "transform_generation", "validity_generation", "owner_generation", "file_generation", "result_generation"))
        and str(value.get("schema") or "").startswith("AP242_")
        and value.get("decoded_schema") == value.get("schema")
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and bool(pmi)
        and all(row[0] > 0 and row[1] and row[2] and math.isfinite(row[3]) for row in pmi)
        and value.get("decoded_pmi_annotations") == value.get("pmi_annotations")
        and bool(names)
        and len(product_ids) == len(names)
        and all(row[1] for row in names)
        and value.get("decoded_product_names") == value.get("product_names")
        and {row[0] for row in colors} == product_ids
        and all(len(row[1]) == 3 and all(0.0 <= channel <= 1.0 for channel in row[1]) for row in colors)
        and value.get("decoded_colors_rgb") == value.get("colors_rgb")
        and {row[0] for row in occurrences} == product_ids
        and all(row[1].startswith("root/") for row in occurrences)
        and value.get("decoded_occurrence_paths") == value.get("occurrence_paths")
        and {row[0] for row in transforms} == product_ids
        and all(_valid_sha256(row[1]) for row in transforms)
        and value.get("decoded_transform_sha256") == value.get("transform_sha256")
        and {row[0] for row in validity} == product_ids
        and all(row[1] is True for row in validity)
        and value.get("decoded_solid_validity") == value.get("solid_validity")
        and bool(str(value.get("source_owner") or ""))
        and value.get("decoded_source_owner") == value.get("source_owner")
        and _valid_sha256(value.get("ap242_file_sha256"))
        and value.get("decoded_ap242_file_sha256") == value.get("ap242_file_sha256")
    )


def _dxf_profile_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("dxf_generation") or "")
    try:
        layers = [str(item) for item in value.get("layer_names", [])]
        bulges = [[int(row[0]), float(row[1])] for row in value.get("arc_bulges", [])]
        winding = [[str(row[0]), int(row[1])] for row in value.get("loop_winding_signs", [])]
        loops = int(value.get("closed_loop_count"))
        height = float(value.get("extrusion_height_m"))
        area = float(value.get("profile_area_m2"))
        solids = int(value.get("solid_count"))
        volume = float(value.get("extruded_volume_m3"))
    except (IndexError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("unit_generation", "plane_generation", "layer_generation", "arc_generation", "loop_generation", "winding_generation", "extrusion_generation", "topology_generation", "owner_generation", "file_generation", "result_generation"))
        and value.get("length_unit") in {"mm", "m"}
        and value.get("decoded_length_unit") == value.get("length_unit")
        and value.get("sketch_plane") in {"XY", "XZ", "YZ"}
        and value.get("decoded_sketch_plane") == value.get("sketch_plane")
        and bool(layers)
        and len(set(layers)) == len(layers)
        and value.get("decoded_layer_names") == value.get("layer_names")
        and bool(bulges)
        and len({row[0] for row in bulges}) == len(bulges)
        and all(math.isfinite(row[1]) and abs(row[1]) <= 1.0 for row in bulges)
        and value.get("decoded_arc_bulges") == value.get("arc_bulges")
        and winding == [["outer", 1], ["hole", -1]]
        and value.get("decoded_loop_winding_signs") == value.get("loop_winding_signs")
        and loops == len(winding)
        and int(value.get("decoded_closed_loop_count", -1)) == loops
        and height > 0.0
        and float(value.get("decoded_extrusion_height_m")) == height
        and area > 0.0
        and float(value.get("decoded_profile_area_m2")) == area
        and solids == 1
        and int(value.get("decoded_solid_count", -1)) == solids
        and math.isclose(volume, area * height, rel_tol=1.0e-9, abs_tol=1.0e-12)
        and float(value.get("decoded_extruded_volume_m3")) == volume
        and bool(str(value.get("profile_owner") or ""))
        and value.get("decoded_profile_owner") == value.get("profile_owner")
        and _valid_sha256(value.get("dxf_file_sha256"))
        and value.get("decoded_dxf_file_sha256") == value.get("dxf_file_sha256")
    )


def _step_assembly_hierarchy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("assembly_generation") or "")
    try:
        paths = [str(item) for item in value.get("occurrence_paths", [])]
        parts = [str(item) for item in value.get("part_ids", [])]
        transforms = [str(item) for item in value.get("transform_sha256", [])]
        names = [str(item) for item in value.get("product_names", [])]
        colors = [[float(channel) for channel in row] for row in value.get("colors_rgb", [])]
    except (TypeError, ValueError):
        return False
    count = len(paths)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("hierarchy_generation", "occurrence_generation", "transform_generation", "part_generation", "name_generation", "color_generation", "unit_generation", "owner_generation", "file_generation", "result_generation"))
        and count >= 2
        and all(path.startswith("root/") for path in paths)
        and len(set(paths)) == count
        and value.get("decoded_occurrence_paths") == value.get("occurrence_paths")
        and len(parts) == count
        and all(parts)
        and len(set(parts)) < count
        and value.get("decoded_part_ids") == value.get("part_ids")
        and len(transforms) == count
        and all(_valid_sha256(item) for item in transforms)
        and value.get("decoded_transform_sha256") == value.get("transform_sha256")
        and len(names) == count
        and all(names)
        and value.get("decoded_product_names") == value.get("product_names")
        and len(colors) == count
        and all(len(row) == 3 and all(0.0 <= channel <= 1.0 for channel in row) for row in colors)
        and value.get("decoded_colors_rgb") == value.get("colors_rgb")
        and value.get("length_unit") == "m"
        and value.get("decoded_length_unit") == value.get("length_unit")
        and bool(str(value.get("assembly_owner") or ""))
        and value.get("decoded_assembly_owner") == value.get("assembly_owner")
        and _valid_sha256(value.get("step_file_sha256"))
        and value.get("decoded_step_file_sha256") == value.get("step_file_sha256")
    )


def _boolean_history_roundtrip_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("boolean_generation") or "")
    try:
        tolerance = float(value.get("fuzzy_tolerance_m"))
        slivers = int(value.get("sliver_face_count"))
        nonmanifold = int(value.get("nonmanifold_edge_count"))
        inputs = [str(item) for item in value.get("input_shape_ids", [])]
        history = [[str(item) for item in row] for row in value.get("operation_history", [])]
    except (TypeError, ValueError):
        return False
    output = str(value.get("output_shape_id") or "")
    healing = value.get("healing_actions")
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("tolerance_generation", "healing_generation", "sliver_generation", "manifold_generation", "history_generation", "input_generation", "output_generation", "owner_generation", "brep_generation", "result_generation"))
        and value.get("operation") in {"fuse", "cut", "intersect"}
        and value.get("decoded_operation") == value.get("operation")
        and math.isfinite(tolerance)
        and 0.0 < tolerance <= 1.0e-3
        and float(value.get("decoded_fuzzy_tolerance_m", math.nan)) == tolerance
        and isinstance(healing, list)
        and bool(healing)
        and all(isinstance(item, str) and item for item in healing)
        and len(set(healing)) == len(healing)
        and value.get("decoded_healing_actions") == healing
        and slivers == 0
        and int(value.get("decoded_sliver_face_count", -1)) == slivers
        and nonmanifold == 0
        and int(value.get("decoded_nonmanifold_edge_count", -1)) == nonmanifold
        and len(inputs) >= 2
        and all(inputs)
        and len(set(inputs)) == len(inputs)
        and value.get("decoded_input_shape_ids") == value.get("input_shape_ids")
        and bool(output)
        and value.get("decoded_output_shape_id") == output
        and len(history) == len(inputs)
        and {row[0] for row in history if len(row) == 2} == set(inputs)
        and all(len(row) == 2 and row[1] == output for row in history)
        and value.get("decoded_operation_history") == value.get("operation_history")
        and bool(str(value.get("input_owner") or ""))
        and value.get("decoded_input_owner") == value.get("input_owner")
        and bool(str(value.get("output_owner") or ""))
        and value.get("decoded_output_owner") == value.get("output_owner")
        and _valid_sha256(value.get("boolean_brep_sha256"))
        and value.get("decoded_boolean_brep_sha256") == value.get("boolean_brep_sha256")
    )


def _sketch_solve_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("sketch_generation") or "")
    constraints = value.get("constraint_ids")
    references = value.get("reference_geometry_ids")
    try:
        dof = int(value.get("remaining_dof"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("constraint_generation", "dof_generation", "solver_generation", "reference_generation", "unit_generation", "owner_generation", "source_generation", "result_generation"))
        and isinstance(constraints, list) and bool(constraints) and all(isinstance(item, str) and item for item in constraints) and len(set(constraints)) == len(constraints)
        and value.get("solved_constraint_ids") == constraints
        and dof == 0 and int(value.get("solved_remaining_dof", -1)) == dof
        and value.get("solver_status") == "fully_constrained" and value.get("solved_solver_status") == value.get("solver_status")
        and isinstance(references, list) and bool(references) and all(isinstance(item, str) and item for item in references) and len(set(references)) == len(references)
        and value.get("solved_reference_geometry_ids") == references
        and value.get("length_unit") == "m" and value.get("solved_length_unit") == value.get("length_unit")
        and bool(str(value.get("sketch_owner") or "")) and value.get("solved_sketch_owner") == value.get("sketch_owner")
        and _valid_sha256(value.get("sketch_source_sha256")) and value.get("solved_sketch_source_sha256") == value.get("sketch_source_sha256")
        and _valid_sha256(value.get("sketch_result_sha256")) and value.get("accepted_sketch_result_sha256") == value.get("sketch_result_sha256")
    )


def _topological_naming_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("toponame_generation") or "")
    edges = value.get("edge_names"); faces = value.get("face_names"); history = value.get("operation_history"); selected = value.get("selector_result")
    try:
        shape_generation = int(value.get("shape_generation_id"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in ("edge_generation", "face_generation", "history_generation", "ocp_generation", "selector_generation", "shape_generation", "owner_generation", "source_generation", "brep_generation", "result_generation"))
        and isinstance(edges, list) and bool(edges) and all(isinstance(item, str) and item.startswith("edge:") for item in edges) and len(set(edges)) == len(edges) and value.get("replayed_edge_names") == edges
        and isinstance(faces, list) and bool(faces) and all(isinstance(item, str) and item.startswith("face:") for item in faces) and len(set(faces)) == len(faces) and value.get("replayed_face_names") == faces
        and isinstance(history, list) and len(history) >= 2 and all(isinstance(item, str) and item for item in history) and value.get("replayed_operation_history") == history
        and bool(str(value.get("ocp_version") or "")) and value.get("replayed_ocp_version") == value.get("ocp_version")
        and isinstance(selected, list) and bool(selected) and set(selected).issubset(set(faces)) and value.get("replayed_selector_result") == selected
        and shape_generation >= 0 and int(value.get("replayed_shape_generation_id", -1)) == shape_generation
        and bool(str(value.get("feature_owner") or "")) and value.get("replayed_feature_owner") == value.get("feature_owner")
        and _valid_sha256(value.get("feature_source_sha256")) and value.get("replayed_feature_source_sha256") == value.get("feature_source_sha256")
        and _valid_sha256(value.get("feature_brep_sha256")) and value.get("accepted_feature_brep_sha256") == value.get("feature_brep_sha256")
    )


def _brep_roundtrip_topology_bounds_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("brep_generation") or "")
    try:
        tolerance = float(value.get("serialization_tolerance_m"))
        restored_tolerance = float(value.get("restored_serialization_tolerance_m"))
        counts = {str(key): int(item) for key, item in value.get("subshape_counts", {}).items()}
        restored_counts = {
            str(key): int(item) for key, item in value.get("restored_subshape_counts", {}).items()
        }
        bounds_min = [float(item) for item in value.get("bounding_box_min_m", [])]
        restored_bounds_min = [float(item) for item in value.get("restored_bounding_box_min_m", [])]
        bounds_max = [float(item) for item in value.get("bounding_box_max_m", [])]
        restored_bounds_max = [float(item) for item in value.get("restored_bounding_box_max_m", [])]
        volume = float(value.get("volume_m3"))
        restored_volume = float(value.get("restored_volume_m3"))
    except (AttributeError, TypeError, ValueError):
        return False
    expected_count_keys = {"solid", "shell", "face", "edge", "vertex"}
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "tolerance_generation",
                "ocp_generation",
                "topology_generation",
                "bounds_generation",
                "volume_generation",
                "owner_generation",
                "source_generation",
                "restored_generation",
                "result_generation",
            )
        )
        and math.isfinite(tolerance)
        and 0.0 < tolerance <= 1.0e-3
        and restored_tolerance == tolerance
        and bool(str(value.get("ocp_version") or ""))
        and value.get("restored_ocp_version") == value.get("ocp_version")
        and set(counts) == expected_count_keys
        and counts["solid"] >= 1
        and counts["shell"] >= 1
        and all(counts[key] > 0 for key in ("face", "edge", "vertex"))
        and counts["vertex"] - counts["edge"] + counts["face"] == 2 * counts["shell"]
        and restored_counts == counts
        and len(bounds_min) == len(bounds_max) == 3
        and all(math.isfinite(item) for item in bounds_min + bounds_max)
        and all(low < high for low, high in zip(bounds_min, bounds_max))
        and restored_bounds_min == bounds_min
        and restored_bounds_max == bounds_max
        and math.isfinite(volume)
        and volume > 0.0
        and math.isclose(restored_volume, volume, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and bool(str(value.get("shape_owner") or ""))
        and value.get("restored_shape_owner") == value.get("shape_owner")
        and _valid_sha256(value.get("source_brep_sha256"))
        and value.get("restored_source_brep_sha256") == value.get("source_brep_sha256")
        and _valid_sha256(value.get("restored_brep_sha256"))
        and value.get("accepted_restored_brep_sha256") == value.get("restored_brep_sha256")
    )


def _svg_extrusion_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("svg_generation") or "")
    commands = value.get("path_commands")
    transform = value.get("curve_transform")
    try:
        matrix = [[float(item) for item in row] for row in transform]
        scale = float(value.get("document_to_meter_scale"))
        replayed_scale = float(value.get("replayed_document_to_meter_scale"))
        area = float(value.get("profile_area_m2"))
        replayed_area = float(value.get("replayed_profile_area_m2"))
        distance = float(value.get("extrusion_distance_m"))
        replayed_distance = float(value.get("replayed_extrusion_distance_m"))
        volume = float(value.get("extrusion_volume_m3"))
        replayed_volume = float(value.get("replayed_extrusion_volume_m3"))
    except (TypeError, ValueError):
        return False
    unit_scales = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
    determinant = (
        matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        if len(matrix) == 3 and all(len(row) == 3 for row in matrix)
        else math.nan
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "path_generation",
                "fill_generation",
                "transform_generation",
                "unit_generation",
                "wire_generation",
                "face_generation",
                "volume_generation",
                "owner_generation",
                "digest_generation",
                "result_generation",
            )
        )
        and isinstance(commands, list)
        and len(commands) >= 3
        and commands[0] == "M"
        and commands[-1] == "Z"
        and all(
            command in {"M", "L", "H", "V", "C", "S", "Q", "T", "A", "Z"} for command in commands
        )
        and value.get("replayed_path_commands") == commands
        and value.get("fill_rule") in {"nonzero", "evenodd"}
        and value.get("replayed_fill_rule") == value.get("fill_rule")
        and len(matrix) == 3
        and all(len(row) == 3 and all(math.isfinite(item) for item in row) for row in matrix)
        and matrix[2] == [0.0, 0.0, 1.0]
        and math.isfinite(determinant)
        and determinant > 0.0
        and value.get("replayed_curve_transform") == value.get("curve_transform")
        and value.get("document_unit") in unit_scales
        and math.isclose(scale, unit_scales[value.get("document_unit")], rel_tol=0.0, abs_tol=0.0)
        and value.get("replayed_document_unit") == value.get("document_unit")
        and replayed_scale == scale
        and value.get("wire_closed") is True
        and value.get("replayed_wire_closed") is True
        and value.get("face_orientation") == "counterclockwise_positive"
        and value.get("replayed_face_orientation") == value.get("face_orientation")
        and all(math.isfinite(item) and item > 0.0 for item in (area, distance, volume))
        and replayed_area == area
        and replayed_distance == distance
        and math.isclose(volume, area * distance, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and replayed_volume == volume
        and bool(str(value.get("source_owner") or ""))
        and value.get("replayed_source_owner") == value.get("source_owner")
        and _valid_sha256(value.get("svg_source_sha256"))
        and value.get("replayed_svg_source_sha256") == value.get("svg_source_sha256")
        and _valid_sha256(value.get("svg_result_sha256"))
        and value.get("accepted_svg_result_sha256") == value.get("svg_result_sha256")
    )


def _occt_heal_replay_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("heal_generation") or "")
    try:
        tolerance = float(value.get("healing_tolerance_m"))
        replayed_tolerance = float(value.get("replayed_healing_tolerance_m"))
        shell_count = int(value.get("sewn_shell_count"))
        replayed_shell_count = int(value.get("replayed_sewn_shell_count"))
        solid_count = int(value.get("solid_count"))
        replayed_solid_count = int(value.get("replayed_solid_count"))
        orientations = [int(item) for item in value.get("face_orientations", [])]
        replayed_orientations = [int(item) for item in value.get("replayed_face_orientations", [])]
        removed = [int(item) for item in value.get("removed_degenerate_edge_ids", [])]
        replayed_removed = [
            int(item) for item in value.get("replayed_removed_degenerate_edge_ids", [])
        ]
        names = {str(key): int(item) for key, item in value.get("stable_subshape_names", {}).items()}
        replayed_names = {
            str(key): int(item)
            for key, item in value.get("replayed_stable_subshape_names", {}).items()
        }
        counts = {str(key): int(item) for key, item in value.get("roundtrip_subshape_counts", {}).items()}
        replayed_counts = {
            str(key): int(item)
            for key, item in value.get("replayed_roundtrip_subshape_counts", {}).items()
        }
    except (AttributeError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "tolerance_generation", "sew_generation", "orientation_generation",
                "degenerate_generation", "name_generation", "roundtrip_generation",
                "owner_generation", "result_generation",
            )
        )
        and math.isfinite(tolerance) and 0.0 < tolerance <= 1.0e-3
        and replayed_tolerance == tolerance
        and shell_count == 1 and replayed_shell_count == shell_count
        and solid_count == 1 and replayed_solid_count == solid_count
        and len(orientations) >= 4 and all(item == 1 for item in orientations)
        and replayed_orientations == orientations
        and bool(removed) and all(item > 0 for item in removed)
        and len(set(removed)) == len(removed) and replayed_removed == removed
        and bool(names)
        and all(key.startswith("face:") and item > 0 for key, item in names.items())
        and len(set(names.values())) == len(names)
        and replayed_names == names
        and set(counts) == {"solid", "shell", "face", "edge"}
        and counts["solid"] == 1 and counts["shell"] == 1
        and counts["face"] >= 4 and counts["edge"] >= counts["face"]
        and replayed_counts == counts
        and str(value.get("shape_owner") or "").startswith("headless:")
        and value.get("replayed_shape_owner") == value.get("shape_owner")
        and _valid_sha256(value.get("healed_brep_sha256"))
        and value.get("replayed_healed_brep_sha256") == value.get("healed_brep_sha256")
        and _valid_sha256(value.get("heal_result_sha256"))
        and value.get("accepted_heal_result_sha256") == value.get("heal_result_sha256")
    )


def _assembly_hierarchy_replay_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("assembly_generation") or "")
    hierarchy = value.get("hierarchy")
    replayed_hierarchy = value.get("replayed_hierarchy")
    local = value.get("local_locations_m")
    replayed_local = value.get("replayed_local_locations_m")
    global_locations = value.get("global_locations_m")
    replayed_global = value.get("replayed_global_locations_m")
    quantities = value.get("component_quantities")
    replayed_quantities = value.get("replayed_component_quantities")
    bom = value.get("bom_identity")
    replayed_bom = value.get("replayed_bom_identity")
    try:
        joint_axis = [float(item) for item in value.get("joint_axis", [])]
        replayed_joint_axis = [float(item) for item in value.get("replayed_joint_axis", [])]
        collisions = [[str(item) for item in pair] for pair in value.get("collision_pairs", [])]
        replayed_collisions = [
            [str(item) for item in pair] for pair in value.get("replayed_collision_pairs", [])
        ]
        local_vectors = {
            str(key): [float(item) for item in vector] for key, vector in local.items()
        }
        global_vectors = {
            str(key): [float(item) for item in vector] for key, vector in global_locations.items()
        }
    except (AttributeError, TypeError, ValueError):
        return False
    if not isinstance(hierarchy, Mapping) or not isinstance(quantities, Mapping) or not isinstance(bom, Mapping):
        return False
    parent_of = {
        str(child): str(parent)
        for parent, children in hierarchy.items()
        for child in children
    }
    components = set(local_vectors)
    locations_close = all(
        len(local_vectors[name]) == len(global_vectors.get(name, [])) == 3
        and all(math.isfinite(item) for item in local_vectors[name] + global_vectors[name])
        and all(
            math.isclose(
                global_vectors[name][axis],
                local_vectors[name][axis]
                + (global_vectors[parent_of[name]][axis] if parent_of.get(name) in components else 0.0),
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            for axis in range(3)
        )
        for name in components
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "hierarchy_generation", "location_generation", "joint_generation",
                "collision_generation", "quantity_generation", "bom_generation",
                "owner_generation", "result_generation",
            )
        )
        and bool(hierarchy) and replayed_hierarchy == hierarchy
        and components == set(global_vectors) == set(quantities) == set(bom)
        and components == set(parent_of)
        and len(parent_of) == sum(len(children) for children in hierarchy.values())
        and replayed_local == local and replayed_global == global_locations
        and locations_close
        and len(joint_axis) == 3
        and all(math.isfinite(item) for item in joint_axis)
        and math.isclose(sum(item * item for item in joint_axis), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and replayed_joint_axis == joint_axis
        and bool(collisions)
        and all(len(pair) == 2 and pair[0] != pair[1] and set(pair).issubset(components) for pair in collisions)
        and replayed_collisions == collisions
        and all(int(item) > 0 for item in quantities.values())
        and replayed_quantities == quantities
        and all(str(item).strip() for item in bom.values())
        and len(set(str(item) for item in bom.values())) == len(bom)
        and replayed_bom == bom
        and str(value.get("assembly_owner") or "").startswith("headless:")
        and value.get("replayed_assembly_owner") == value.get("assembly_owner")
        and _valid_sha256(value.get("assembly_result_sha256"))
        and value.get("accepted_assembly_result_sha256") == value.get("assembly_result_sha256")
    )


def _right_handed_orthonormal_frame(rows: list[list[float]]) -> bool:
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        return False
    if any(not math.isfinite(item) for row in rows for item in row):
        return False
    for left_index, left in enumerate(rows):
        for right_index, right in enumerate(rows):
            dot = sum(a * b for a, b in zip(left, right))
            expected = 1.0 if left_index == right_index else 0.0
            if not math.isclose(dot, expected, rel_tol=0.0, abs_tol=1.0e-12):
                return False
    determinant = (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )
    return math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1.0e-12)


def _oriented_bounding_box_inertia_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("obb_generation") or "")
    try:
        center = [float(item) for item in value.get("obb_center_m", [])]
        replayed_center = [float(item) for item in value.get("replayed_obb_center_m", [])]
        extents = [float(item) for item in value.get("obb_half_extents_m", [])]
        replayed_extents = [
            float(item) for item in value.get("replayed_obb_half_extents_m", [])
        ]
        axes = [[float(item) for item in row] for row in value.get("obb_axes", [])]
        replayed_axes = [
            [float(item) for item in row] for row in value.get("replayed_obb_axes", [])
        ]
        moments = [float(item) for item in value.get("principal_moments_kg_m2", [])]
        replayed_moments = [
            float(item) for item in value.get("replayed_principal_moments_kg_m2", [])
        ]
        principal_axes = [
            [float(item) for item in row] for row in value.get("principal_axes", [])
        ]
        replayed_principal_axes = [
            [float(item) for item in row]
            for row in value.get("replayed_principal_axes", [])
        ]
        local_com = [float(item) for item in value.get("center_of_mass_local_m", [])]
        replayed_local_com = [
            float(item) for item in value.get("replayed_center_of_mass_local_m", [])
        ]
        world_com = [float(item) for item in value.get("center_of_mass_world_m", [])]
        replayed_world_com = [
            float(item) for item in value.get("replayed_center_of_mass_world_m", [])
        ]
        transform = [
            [float(item) for item in row] for row in value.get("local_to_world_transform", [])
        ]
        replayed_transform = [
            [float(item) for item in row]
            for row in value.get("replayed_local_to_world_transform", [])
        ]
    except (TypeError, ValueError):
        return False
    transformed_com = (
        [
            sum(transform[row][column] * local_com[column] for column in range(3))
            + transform[row][3]
            for row in range(3)
        ]
        if len(transform) == 4
        and all(len(row) == 4 for row in transform)
        and len(local_com) == 3
        else []
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "box_generation",
                "inertia_generation",
                "axis_generation",
                "com_generation",
                "transform_generation",
                "unit_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and len(center) == 3
        and all(math.isfinite(item) for item in center)
        and replayed_center == center
        and len(extents) == 3
        and all(math.isfinite(item) and item > 0.0 for item in extents)
        and replayed_extents == extents
        and _right_handed_orthonormal_frame(axes)
        and replayed_axes == axes
        and len(moments) == 3
        and all(math.isfinite(item) and item > 0.0 for item in moments)
        and moments == sorted(moments)
        and moments[2] <= moments[0] + moments[1]
        and replayed_moments == moments
        and _right_handed_orthonormal_frame(principal_axes)
        and replayed_principal_axes == principal_axes
        and len(local_com) == len(world_com) == 3
        and all(math.isfinite(item) for item in local_com + world_com)
        and replayed_local_com == local_com
        and replayed_world_com == world_com
        and len(transform) == 4
        and all(len(row) == 4 and all(math.isfinite(item) for item in row) for row in transform)
        and transform[3] == [0.0, 0.0, 0.0, 1.0]
        and _right_handed_orthonormal_frame([row[:3] for row in transform[:3]])
        and replayed_transform == transform
        and all(
            math.isclose(transformed_com[index], world_com[index], rel_tol=0.0, abs_tol=1.0e-12)
            for index in range(3)
        )
        and value.get("length_unit") == "m"
        and value.get("replayed_length_unit") == "m"
        and value.get("inertia_unit") == "kg*m^2"
        and value.get("replayed_inertia_unit") == "kg*m^2"
        and str(value.get("shape_owner") or "").startswith("headless:")
        and value.get("replayed_shape_owner") == value.get("shape_owner")
        and _valid_sha256(value.get("obb_result_sha256"))
        and value.get("accepted_obb_result_sha256") == value.get("obb_result_sha256")
    )


def _tessellation_replay_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("tessellation_generation") or "")
    try:
        linear = float(value.get("linear_deflection_m"))
        replayed_linear = float(value.get("replayed_linear_deflection_m"))
        angular = float(value.get("angular_deflection_rad"))
        replayed_angular = float(value.get("replayed_angular_deflection_rad"))
        triangles = int(value.get("triangle_count"))
        replayed_triangles = int(value.get("replayed_triangle_count"))
        vertices = [
            [float(item) for item in vertex]
            for vertex in value.get("vertex_coordinates_m", [])
        ]
        replayed_vertices = [
            [float(item) for item in vertex]
            for vertex in value.get("replayed_vertex_coordinates_m", [])
        ]
        nonmanifold = int(value.get("nonmanifold_edge_count"))
        replayed_nonmanifold = int(value.get("replayed_nonmanifold_edge_count"))
    except (TypeError, ValueError):
        return False
    vertex_tuples = [tuple(vertex) for vertex in vertices]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "linear_generation",
                "angular_generation",
                "triangle_generation",
                "vertex_generation",
                "orientation_generation",
                "watertight_generation",
                "stl_generation",
                "brep_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and math.isfinite(linear)
        and linear > 0.0
        and replayed_linear == linear
        and math.isfinite(angular)
        and 0.0 < angular < math.pi
        and replayed_angular == angular
        and triangles >= 4
        and replayed_triangles == triangles
        and len(vertices) >= 4
        and all(len(vertex) == 3 and all(math.isfinite(item) for item in vertex) for vertex in vertices)
        and len(set(vertex_tuples)) == len(vertex_tuples)
        and triangles >= 2 * len(vertices) - 4
        and replayed_vertices == vertices
        and value.get("outward_orientation") is True
        and value.get("replayed_outward_orientation") is True
        and value.get("watertight") is True
        and value.get("replayed_watertight") is True
        and nonmanifold == 0
        and replayed_nonmanifold == 0
        and str(value.get("stl_owner") or "").startswith("headless:")
        and value.get("replayed_stl_owner") == value.get("stl_owner")
        and _valid_sha256(value.get("source_brep_sha256"))
        and value.get("replayed_source_brep_sha256") == value.get("source_brep_sha256")
        and _valid_sha256(value.get("stl_sha256"))
        and value.get("replayed_stl_sha256") == value.get("stl_sha256")
        and _valid_sha256(value.get("tessellation_result_sha256"))
        and value.get("accepted_tessellation_result_sha256")
        == value.get("tessellation_result_sha256")
    )


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
    step_import_metadata_topology_identity_ok = True
    mesh_export_facet_normal_identity_ok = True
    brep_step_roundtrip_identity_ok = True
    fresh_subprocess_result_identity_ok = True
    step_roundtrip_metadata_identity_ok = True
    occ_build_fingerprint_identity_ok = True
    stl_tessellation_component_identity_ok = True
    builder_context_identity_ok = True
    step_import_hierarchy_identity_ok = True
    brep_cache_generation_identity_ok = True
    step_ap242_import_identity_ok = True
    occt_shape_cache_identity_ok = True
    dxf_face_reconstruction_identity_ok = True
    three_mf_reconstruction_identity_ok = True
    step_assembly_reconstruction_identity_ok = True
    stl_solid_reconstruction_identity_ok = True
    step_ap242_context_owner_mass_identity_ok = True
    stl_brep_tessellation_error_identity_ok = True
    step_assembly_instance_metadata_identity_ok = True
    stl_repair_closure_identity_ok = True
    brep_semantic_roundtrip_identity_ok = True
    gltf_roundtrip_identity_ok = True
    step_semantic_roundtrip_identity_ok = True
    brep_periodic_seam_identity_ok = True
    step_ap242_pmi_roundtrip_identity_ok = True
    step_assembly_hierarchy_identity_ok = True
    boolean_history_roundtrip_identity_ok = True
    dxf_profile_roundtrip_identity_ok = True
    sketch_solve_identity_ok = True
    topological_naming_identity_ok = True
    brep_roundtrip_topology_bounds_identity_ok = True
    svg_extrusion_identity_ok = True
    occt_heal_replay_identity_ok = True
    assembly_hierarchy_replay_identity_ok = True
    oriented_bounding_box_inertia_identity_ok = True
    tessellation_replay_identity_ok = True
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
        step_import_metadata_topology_identity_ok = False
        mesh_export_facet_normal_identity_ok = False
        brep_step_roundtrip_identity_ok = False
        fresh_subprocess_result_identity_ok = False
        step_roundtrip_metadata_identity_ok = False
        occ_build_fingerprint_identity_ok = False
        stl_tessellation_component_identity_ok = False
        builder_context_identity_ok = False
        step_import_hierarchy_identity_ok = False
        brep_cache_generation_identity_ok = False
        step_ap242_import_identity_ok = False
        occt_shape_cache_identity_ok = False
        dxf_face_reconstruction_identity_ok = False
        three_mf_reconstruction_identity_ok = False
        step_ap242_context_owner_mass_identity_ok = False
        stl_brep_tessellation_error_identity_ok = False
        step_assembly_instance_metadata_identity_ok = False
        stl_repair_closure_identity_ok = False
        brep_semantic_roundtrip_identity_ok = False
        gltf_roundtrip_identity_ok = False
        step_semantic_roundtrip_identity_ok = False
        brep_periodic_seam_identity_ok = False
        sketch_solve_identity_ok = False
        topological_naming_identity_ok = False
        brep_roundtrip_topology_bounds_identity_ok = False
        svg_extrusion_identity_ok = False
        occt_heal_replay_identity_ok = False
        assembly_hierarchy_replay_identity_ok = False
        oriented_bounding_box_inertia_identity_ok = False
        tessellation_replay_identity_ok = False
    elif replay_identity_present:
        source_commit = str(replay_identity_value.get("source_commit", "")).lower()
        replayed_commit = str(replay_identity_value.get("replayed_source_commit", "")).lower()
        artifacts = replay_identity_value.get("cad_artifacts") or []
        kernel = replay_identity_value.get("external_kernel") or {}
        artifact_rows_valid = isinstance(artifacts, list) and all(
            isinstance(row, Mapping) for row in artifacts
        )
        artifact_names = (
            [str(row.get("name", "")).strip() for row in artifacts] if artifact_rows_valid else []
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
        step_import_metadata_topology_identity_ok = (
            _step_import_metadata_topology_identity_ok(
                replay_identity_value.get(
                    "step_import_label_color_unit_topology_generation_identity"
                )
            )
        )
        mesh_export_facet_normal_identity_ok = _mesh_export_facet_normal_identity_ok(
            replay_identity_value.get(
                "mesh_export_facet_normal_tolerance_shape_digest_generation_identity"
            )
        )
        brep_step_roundtrip_identity_ok = _brep_step_roundtrip_identity_ok(
            replay_identity_value.get(
                "brep_step_roundtrip_tolerance_orientation_volume_generation_identity"
            )
        )
        fresh_subprocess_result_identity_ok = _fresh_subprocess_result_identity_ok(
            replay_identity_value.get(
                "fresh_subprocess_timeout_exception_cache_output_generation_identity"
            )
        )
        step_roundtrip_metadata_identity_ok = (
            _step_label_color_unit_hierarchy_shape_roundtrip_identity_ok(
                replay_identity_value.get(
                    "step_label_color_unit_hierarchy_shape_roundtrip_identity"
                )
            )
        )
        occ_build_fingerprint_identity_ok = (
            _occ_version_tolerance_tessellation_cache_build_identity_ok(
                replay_identity_value.get(
                    "occ_version_tolerance_tessellation_cache_build_fingerprint_identity"
                )
            )
        )
        stl_tessellation_component_identity_ok = _stl_tessellation_component_identity_ok(
            replay_identity_value.get(
                "stl_chord_tolerance_triangle_normal_orientation_component_digest_generation_identity"
            )
        )
        builder_context_identity_ok = _builder_context_identity_ok(
            replay_identity_value.get(
                "builder_context_workplane_local_frame_part_identity_cache_generation_identity"
            )
        )
        step_import_hierarchy_identity_ok = _step_import_hierarchy_identity_ok(
            replay_identity_value.get(
                "step_import_unit_hierarchy_placement_color_shape_checksum_generation_identity"
            )
        )
        brep_cache_generation_identity_ok = _brep_cache_generation_identity_ok(
            replay_identity_value.get(
                "brep_serialization_kernel_tolerance_location_cache_generation_identity"
            )
        )
        step_ap242_import_identity_ok = _step_ap242_import_identity_ok(
            replay_identity_value.get(
                "step_ap242_context_product_uuid_unit_color_placement_shape_file_generation_identity"
            )
        )
        occt_shape_cache_identity_ok = _occt_shape_cache_identity_ok(
            replay_identity_value.get(
                "occt_kernel_shape_location_tolerance_triangulation_serialization_cache_generation_identity"
            )
        )
        dxf_face_reconstruction_identity_ok = _dxf_face_reconstruction_identity_ok(
            replay_identity_value.get(
                "dxf_arc_spline_layer_plane_unit_closed_wire_orientation_face_digest_generation_identity"
            )
        )
        three_mf_reconstruction_identity_ok = _three_mf_reconstruction_identity_ok(
            replay_identity_value.get(
                "three_mf_component_transform_triangle_winding_material_watertight_volume_unit_digest_generation_identity"
            )
        )
        step_assembly_reconstruction_identity_ok = (
            _step_assembly_reconstruction_identity_ok(
                replay_identity_value.get(
                    "step_assembly_product_color_instance_transform_unit_hierarchy_file_generation_identity"
                )
            )
        )
        stl_solid_reconstruction_identity_ok = _stl_solid_reconstruction_identity_ok(
            replay_identity_value.get(
                "stl_facet_normal_winding_watertight_tolerance_volume_unit_file_generation_identity"
            )
        )
        step_ap242_context_owner_mass_identity_ok = (
            _step_ap242_context_owner_mass_identity_ok(
                replay_identity_value.get(
                    "step_ap242_representation_context_external_owner_occurrence_transform_mass_product_structure_file_generation_identity"
                )
            )
        )
        stl_brep_tessellation_error_identity_ok = (
            _stl_brep_tessellation_error_identity_ok(
                replay_identity_value.get(
                    "stl_tessellation_source_brep_chord_angle_facet_component_deviation_area_volume_digest_generation_identity"
                )
            )
        )
        step_assembly_instance_metadata_identity_ok = (
            _step_assembly_instance_metadata_identity_ok(
                replay_identity_value.get(
                    "step_assembly_instance_transform_unit_color_material_uuid_component_volume_file_generation_identity"
                )
            )
        )
        stl_repair_closure_identity_ok = _stl_repair_closure_identity_ok(
            replay_identity_value.get(
                "stl_repair_merge_tolerance_normal_duplicate_boundary_watertight_volume_unit_file_generation_identity"
            )
        )
        brep_semantic_roundtrip_identity_ok = _brep_semantic_roundtrip_identity_ok(
            replay_identity_value.get(
                "brep_face_pcurve_wire_orientation_location_tolerance_surface_serializer_shape_generation_identity"
            )
        )
        gltf_roundtrip_identity_ok = _gltf_roundtrip_identity_ok(
            replay_identity_value.get(
                "gltf_node_hierarchy_transform_winding_material_unit_tessellation_volume_file_generation_identity"
            )
        )
        step_semantic_roundtrip_identity_ok = _step_semantic_roundtrip_identity_ok(
            replay_identity_value.get(
                "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
            )
        )
        brep_periodic_seam_identity_ok = _brep_periodic_seam_identity_ok(
            replay_identity_value.get(
                "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
            )
        )
        step_ap242_pmi_roundtrip_identity_ok = _step_ap242_pmi_roundtrip_identity_ok(
            replay_identity_value.get(
                "step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"
            )
        )
        dxf_profile_roundtrip_identity_ok = _dxf_profile_roundtrip_identity_ok(
            replay_identity_value.get(
                "dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"
            )
        )
        step_assembly_hierarchy_identity_ok = _step_assembly_hierarchy_identity_ok(
            replay_identity_value.get(
                "step_assembly_hierarchy_occurrence_transform_repeated_part_name_color_unit_owner_file_result_generation_identity"
            )
        )
        boolean_history_roundtrip_identity_ok = _boolean_history_roundtrip_identity_ok(
            replay_identity_value.get(
                "boolean_tolerance_healing_sliver_nonmanifold_operation_history_input_output_owner_brep_result_generation_identity"
            )
        )
        sketch_solve_identity_ok = _sketch_solve_identity_ok(
            replay_identity_value.get(
                "sketch_constraint_dof_solver_reference_unit_owner_source_result_generation_identity"
            )
        )
        topological_naming_identity_ok = _topological_naming_identity_ok(
            replay_identity_value.get(
                "topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"
            )
        )
        brep_roundtrip_topology_bounds_identity_ok = _brep_roundtrip_topology_bounds_identity_ok(
            replay_identity_value.get(
                "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
            )
        )
        svg_extrusion_identity_ok = _svg_extrusion_identity_ok(
            replay_identity_value.get(
                "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
            )
        )
        occt_heal_replay_identity_ok = _occt_heal_replay_identity_ok(
            replay_identity_value.get(
                "occt_heal_tolerance_sew_orientation_degenerate_stablename_roundtrip_owner_digest_generation_identity"
            )
        )
        assembly_hierarchy_replay_identity_ok = _assembly_hierarchy_replay_identity_ok(
            replay_identity_value.get(
                "assembly_hierarchy_location_joint_axis_collision_quantity_bom_owner_result_generation_identity"
            )
        )
        oriented_bounding_box_inertia_identity_ok = _oriented_bounding_box_inertia_identity_ok(
            replay_identity_value.get(
                "oriented_bounding_box_principal_inertia_axis_com_frame_transform_unit_owner_result_generation_identity"
            )
        )
        tessellation_replay_identity_ok = _tessellation_replay_identity_ok(
            replay_identity_value.get(
                "tessellation_linear_angular_deflection_triangle_vertex_orientation_watertight_stl_brep_owner_result_generation_identity"
            )
        )
    joint_names = {
        str(name) for row in components for name in (row.get("joint_names") or []) if str(name)
    }
    connection_endpoints = {
        str(connection.get(side, "")) for connection in connections for side in ("from", "to")
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
        "step_import_uses_current_source_labels_colors_units_and_topology": (
            step_import_metadata_topology_identity_ok
        ),
        "mesh_export_uses_current_shape_facets_normals_and_tolerances": (
            mesh_export_facet_normal_identity_ok
        ),
        "brep_step_roundtrip_uses_current_tolerances_orientation_volume_and_topology": (
            brep_step_roundtrip_identity_ok
        ),
        "fresh_subprocess_rejects_timeout_exception_cache_and_stale_output": (
            fresh_subprocess_result_identity_ok
        ),
        "step_roundtrip_preserves_labels_colors_units_hierarchy_and_shapes": (
            step_roundtrip_metadata_identity_ok
        ),
        "occ_build_uses_current_version_tolerances_tessellation_cache_and_fingerprint": (
            occ_build_fingerprint_identity_ok
        ),
        "stl_handoff_uses_current_tolerances_triangles_normals_components_and_digests": (
            stl_tessellation_component_identity_ok
        ),
        "builder_replay_uses_current_context_workplane_frame_parts_and_cache": (
            builder_context_identity_ok
        ),
        "step_import_uses_current_units_hierarchy_placements_colors_shapes_and_checksums": (
            step_import_hierarchy_identity_ok
        ),
        "brep_cache_uses_current_kernel_tolerance_location_shape_and_digest": (
            brep_cache_generation_identity_ok
        ),
        "step_ap242_import_uses_current_context_products_units_colors_placements_shapes_and_file": (
            step_ap242_import_identity_ok
        ),
        "occt_shape_cache_uses_current_kernel_shape_location_tolerance_triangulation_and_serialization": (
            occt_shape_cache_identity_ok
        ),
        "dxf_faces_use_current_arcs_splines_layers_plane_units_closure_orientation_and_digests": (
            dxf_face_reconstruction_identity_ok
        ),
        "three_mf_imports_use_current_components_transforms_winding_materials_watertight_volumes_units_and_digests": (
            three_mf_reconstruction_identity_ok
        ),
        "step_assemblies_use_current_products_colors_instances_transforms_units_hierarchy_and_digests": (
            step_assembly_reconstruction_identity_ok
        ),
        "stl_solids_use_current_normals_winding_watertight_edges_tolerance_volume_units_and_digests": (
            stl_solid_reconstruction_identity_ok
        ),
        "step_ap242_occurrences_use_current_schema_context_owners_transforms_mass_structure_and_file": (
            step_ap242_context_owner_mass_identity_ok
        ),
        "stl_tessellations_use_current_brep_chord_angle_facets_components_deviation_area_volume_and_digests": (
            stl_brep_tessellation_error_identity_ok
        ),
        "step_assemblies_use_current_instances_transforms_units_colors_materials_uuids_components_volume_and_file": (
            step_assembly_instance_metadata_identity_ok
        ),
        "stl_repairs_use_current_tolerance_normals_duplicates_boundaries_watertight_volume_unit_and_files": (
            stl_repair_closure_identity_ok
        ),
        "brep_roundtrips_use_current_pcurves_wire_orientations_locations_tolerances_surfaces_serializer_and_shape": (
            brep_semantic_roundtrip_identity_ok
        ),
        "gltf_roundtrips_use_current_hierarchy_transforms_winding_materials_units_tessellation_volume_and_file": (
            gltf_roundtrip_identity_ok
        ),
        "step_roundtrips_use_current_units_products_colors_transforms_shapes_validity_owner_and_file": (
            step_semantic_roundtrip_identity_ok
        ),
        "periodic_brep_roundtrips_use_current_seams_orientations_pcurves_tolerances_manifold_serializer_and_shape": (
            brep_periodic_seam_identity_ok
        ),
        "step_ap242_roundtrips_use_current_pmi_units_names_colors_occurrences_transforms_validity_owner_and_file": (
            step_ap242_pmi_roundtrip_identity_ok
        ),
        "dxf_profiles_use_current_units_plane_layers_arcs_bulges_winding_extrusion_topology_owner_and_file": (
            dxf_profile_roundtrip_identity_ok
        ),
        "step_assemblies_use_current_hierarchy_occurrences_repeated_parts_transforms_names_colors_units_owner_and_file": (
            step_assembly_hierarchy_identity_ok
        ),
        "boolean_roundtrips_use_current_operation_tolerance_healing_slivers_manifold_history_owners_and_brep": (
            boolean_history_roundtrip_identity_ok
        ),
        "sketch_solves_use_current_constraints_dof_status_references_units_owner_source_and_result": (
            sketch_solve_identity_ok
        ),
        "topological_names_use_current_edges_faces_history_ocp_selector_shape_owner_source_and_brep": (
            topological_naming_identity_ok
        ),
        "brep_roundtrips_use_current_tolerance_ocp_topology_bounds_volume_owner_source_and_restored_shape": (
            brep_roundtrip_topology_bounds_identity_ok
        ),
        "svg_extrusions_use_current_paths_fill_transform_units_wire_face_volume_owner_and_digests": (
            svg_extrusion_identity_ok
        ),
        "occt_heals_use_current_tolerance_sewing_orientation_degenerate_edges_names_roundtrip_owner_and_digests": (
            occt_heal_replay_identity_ok
        ),
        "assemblies_use_current_hierarchy_locations_joint_axis_collisions_quantities_bom_owner_and_result": (
            assembly_hierarchy_replay_identity_ok
        ),
        "mass_property_replays_use_current_obb_principal_inertia_axes_com_transform_units_owner_and_result": (
            oriented_bounding_box_inertia_identity_ok
        ),
        "tessellation_replays_use_current_deflections_triangles_vertices_orientation_watertight_stl_brep_owner_and_result": (
            tessellation_replay_identity_ok
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
