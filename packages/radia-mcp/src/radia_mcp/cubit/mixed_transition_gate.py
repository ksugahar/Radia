"""Gates for a conformal Cubit hex-pyramid-tet transition cell."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path


_VOLUME_FAMILIES = ("hex", "pyramid", "tet", "wedge")


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be a mapping")
    return value


def _rows(value: object, field: str) -> list[Mapping[str, object]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{field} must be a sequence of mappings")
    rows = list(value)
    if not all(isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{field} must contain mappings")
    return rows


def _finite(value: object, field: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a finite number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field} must be a finite number")
    return number


def _counts(value: object, field: str) -> dict[str, int]:
    row = _mapping(value, field)
    result = {family: int(row.get(family, 0)) for family in _VOLUME_FAMILIES}
    if any(count < 0 for count in result.values()):
        raise ValueError(f"{field} counts must be nonnegative")
    return result


def _opposed_transition_face_orientation_ok(identity: object) -> bool:
    if not isinstance(identity, Mapping):
        return False
    try:
        shared = [int(value) for value in identity.get("shared_face_node_ids", [])]
        pyramid = [int(value) for value in identity.get("pyramid_face_node_ids", [])]
        hex_face = [int(value) for value in identity.get("hex_face_node_ids", [])]
        normal_dot = float(identity.get("opposed_outward_normal_dot"))
    except (TypeError, ValueError):
        return False
    return (
        len(shared) >= 3
        and len(set(shared)) == len(shared)
        and set(pyramid) == set(shared) == set(hex_face)
        and pyramid == shared
        and hex_face == list(reversed(pyramid))
        and math.isclose(normal_dot, -1.0, rel_tol=0.0, abs_tol=1.0e-12)
    )


def _periodic_hex_node_pair_transform_frame_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    source_nodes = list(identity.get("source_node_ids") or [])
    target_nodes = list(identity.get("target_node_ids") or [])
    mesh_generation = str(identity.get("mesh_generation") or "")
    transform_generation = str(
        identity.get("periodic_transform_generation") or ""
    )
    transform_digest = str(identity.get("transform_matrix_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("node_pair_mesh_generation") == mesh_generation
        and identity.get("periodic_transform_frame_mesh_generation")
        == mesh_generation
        and bool(transform_generation)
        and identity.get("node_pair_periodic_transform_generation")
        == transform_generation
        and bool(source_nodes)
        and len(source_nodes) == len(target_nodes)
        and len(set(source_nodes)) == len(source_nodes)
        and len(set(target_nodes)) == len(target_nodes)
        and list(identity.get("paired_source_node_ids") or []) == source_nodes
        and list(identity.get("paired_target_node_ids") or []) == target_nodes
        and identity.get("coordinate_frame") == "global_cartesian"
        and identity.get("node_pair_coordinate_frame")
        == identity.get("coordinate_frame")
        and len(transform_digest) == 64
        and all(character in "0123456789abcdef" for character in transform_digest)
        and identity.get("applied_transform_matrix_sha256") == transform_digest
    )


def _pyramid_transition_face_diagonal_convention_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    pyramid_faces = list(identity.get("pyramid_face_ids") or [])
    convention = str(identity.get("diagonal_convention") or "")
    mesh_generation = str(identity.get("transition_mesh_generation") or "")
    connectivity_digest = str(
        identity.get("transition_face_connectivity_sha256") or ""
    )
    return (
        bool(mesh_generation)
        and identity.get("tet_neighbor_mesh_generation") == mesh_generation
        and identity.get("hex_neighbor_mesh_generation") == mesh_generation
        and bool(pyramid_faces)
        and len(set(pyramid_faces)) == len(pyramid_faces)
        and list(identity.get("tet_neighbor_face_ids") or []) == pyramid_faces
        and list(identity.get("hex_neighbor_face_ids") or []) == pyramid_faces
        and convention == "canonical_node_0_to_2"
        and identity.get("tet_neighbor_diagonal_convention") == convention
        and identity.get("hex_neighbor_diagonal_convention") == convention
        and len(connectivity_digest) == 64
        and all(
            character in "0123456789abcdef" for character in connectivity_digest
        )
        and identity.get("neighbor_face_connectivity_sha256")
        == connectivity_digest
    )


def _block_attribute_material_id_merge_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        block_ids = [int(value) for value in identity.get("block_ids", [])]
        material_ids = [int(value) for value in identity.get("material_ids", [])]
        attribute_block_ids = [
            int(value) for value in identity.get("block_attribute_block_ids", [])
        ]
        exported_material_ids = [
            int(value) for value in identity.get("exported_material_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    merge_generation = str(identity.get("final_merge_generation") or "")
    map_digest = str(identity.get("block_material_map_sha256") or "")
    return (
        bool(merge_generation)
        and identity.get("block_attribute_merge_generation") == merge_generation
        and identity.get("material_id_map_merge_generation") == merge_generation
        and bool(block_ids)
        and len(block_ids) == len(material_ids)
        and len(set(block_ids)) == len(block_ids)
        and len(set(material_ids)) == len(material_ids)
        and all(value > 0 for value in block_ids + material_ids)
        and attribute_block_ids == block_ids
        and exported_material_ids == material_ids
        and len(map_digest) == 64
        and all(character in "0123456789abcdef" for character in map_digest)
        and identity.get("exported_block_material_map_sha256") == map_digest
    )


def _high_order_exodus_node_permutation_export_order_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        element_order = int(identity.get("element_order", 0))
        source_order = [int(value) for value in identity.get("source_node_order", [])]
        permutation = [int(value) for value in identity.get("permutation_table", [])]
        written_order = [int(value) for value in identity.get("written_node_order", [])]
    except (TypeError, ValueError):
        return False
    export_generation = str(identity.get("export_order_generation") or "")
    permutation_digest = str(identity.get("node_permutation_sha256") or "")
    expected_written_order = (
        [source_order[index - 1] for index in permutation]
        if permutation
        and len(source_order) == len(permutation)
        and sorted(permutation) == list(range(1, len(permutation) + 1))
        else []
    )
    return (
        bool(export_generation)
        and identity.get("permutation_table_export_order_generation")
        == export_generation
        and identity.get("writer_export_order_generation") == export_generation
        and element_order >= 2
        and bool(source_order)
        and len(set(source_order)) == len(source_order)
        and written_order == expected_written_order
        and len(permutation_digest) == 64
        and all(
            character in "0123456789abcdef" for character in permutation_digest
        )
        and identity.get("written_node_permutation_sha256")
        == permutation_digest
    )


def _hex_sideset_outward_face_ordinal_volume_reorder_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        element_ids = [int(value) for value in identity.get("element_ids", [])]
        face_ordinals = [int(value) for value in identity.get("face_ordinals", [])]
        exported_element_ids = [
            int(value) for value in identity.get("exported_element_ids", [])
        ]
        exported_face_ordinals = [
            int(value) for value in identity.get("exported_face_ordinals", [])
        ]
        normal_signs = [
            int(value) for value in identity.get("outward_normal_signs", [])
        ]
        exported_normal_signs = [
            int(value)
            for value in identity.get("exported_outward_normal_signs", [])
        ]
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    reorder_generation = str(
        identity.get("volume_connectivity_reorder_generation") or ""
    )
    mapping_digest = str(identity.get("element_face_map_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("face_ordinal_mesh_generation") == mesh_generation
        and identity.get("normal_ownership_mesh_generation") == mesh_generation
        and bool(reorder_generation)
        and identity.get("face_ordinal_connectivity_reorder_generation")
        == reorder_generation
        and identity.get("normal_connectivity_reorder_generation")
        == reorder_generation
        and bool(element_ids)
        and len(set(element_ids)) == len(element_ids)
        and all(value > 0 for value in element_ids)
        and len(face_ordinals) == len(element_ids)
        and all(1 <= value <= 6 for value in face_ordinals)
        and exported_element_ids == element_ids
        and exported_face_ordinals == face_ordinals
        and len(normal_signs) == len(element_ids)
        and all(value in {-1, 1} for value in normal_signs)
        and exported_normal_signs == normal_signs
        and len(mapping_digest) == 64
        and all(character in "0123456789abcdef" for character in mapping_digest)
        and identity.get("exported_element_face_map_sha256") == mapping_digest
    )


def _sweep_layer_bias_source_curve_orientation_generation_ok(
    identity: object,
) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        curve_ids = [int(value) for value in identity.get("source_curve_ids", [])]
        curve_orientations = [
            int(value) for value in identity.get("source_curve_orientations", [])
        ]
        biased_curve_ids = [
            int(value) for value in identity.get("biased_curve_ids", [])
        ]
        biased_orientations = [
            int(value) for value in identity.get("biased_curve_orientations", [])
        ]
        interval_counts = [
            int(value) for value in identity.get("interval_counts", [])
        ]
        biased_interval_counts = [
            int(value) for value in identity.get("biased_interval_counts", [])
        ]
        bias_factors = [float(value) for value in identity.get("bias_factors", [])]
        applied_bias_factors = [
            float(value) for value in identity.get("applied_bias_factors", [])
        ]
    except (TypeError, ValueError):
        return False
    sweep_generation = str(identity.get("sweep_generation") or "")
    orientation_generation = str(
        identity.get("source_curve_orientation_generation") or ""
    )
    bias_digest = str(identity.get("curve_bias_map_sha256") or "")
    return (
        bool(sweep_generation)
        and identity.get("layer_bias_sweep_generation") == sweep_generation
        and bool(orientation_generation)
        and identity.get("layer_bias_curve_orientation_generation")
        == orientation_generation
        and bool(curve_ids)
        and len(set(curve_ids)) == len(curve_ids)
        and all(value > 0 for value in curve_ids)
        and len(curve_orientations) == len(curve_ids)
        and all(value in {-1, 1} for value in curve_orientations)
        and biased_curve_ids == curve_ids
        and biased_orientations == curve_orientations
        and len(interval_counts) == len(curve_ids)
        and all(value > 0 for value in interval_counts)
        and biased_interval_counts == interval_counts
        and len(bias_factors) == len(curve_ids)
        and all(math.isfinite(value) and value > 0.0 for value in bias_factors)
        and applied_bias_factors == bias_factors
        and len(bias_digest) == 64
        and all(character in "0123456789abcdef" for character in bias_digest)
        and identity.get("applied_curve_bias_map_sha256") == bias_digest
    )


def _exodus_sideset_element_face_topology_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        element_ids = [int(value) for value in identity.get("element_ids", [])]
        ordinals = [
            int(value)
            for value in identity.get("element_face_topology_ordinals", [])
        ]
        written_element_ids = [
            int(value) for value in identity.get("written_element_ids", [])
        ]
        written_ordinals = [
            int(value)
            for value in identity.get(
                "written_element_face_topology_ordinals", []
            )
        ]
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    export_generation = str(identity.get("exodus_export_generation") or "")
    topology_digest = str(identity.get("element_face_topology_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("topology_ordinal_mesh_generation") == mesh_generation
        and identity.get("writer_mesh_generation") == mesh_generation
        and bool(export_generation)
        and identity.get("topology_ordinal_export_generation") == export_generation
        and identity.get("writer_export_generation") == export_generation
        and bool(element_ids)
        and len(set(element_ids)) == len(element_ids)
        and all(value > 0 for value in element_ids)
        and len(ordinals) == len(element_ids)
        and all(1 <= value <= 6 for value in ordinals)
        and written_element_ids == element_ids
        and written_ordinals == ordinals
        and len(topology_digest) == 64
        and all(character in "0123456789abcdef" for character in topology_digest)
        and identity.get("written_element_face_topology_sha256") == topology_digest
    )


def _high_order_quality_reference_coordinate_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        element_order = int(identity.get("element_order", 0))
        reference_node_count = int(identity.get("reference_node_count", 0))
        quality_reference_node_count = int(
            identity.get("quality_reference_node_count", 0)
        )
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    order_generation = str(identity.get("element_order_generation") or "")
    sampling_rule = str(identity.get("jacobian_sampling_rule") or "")
    coordinate_digest = str(identity.get("reference_coordinates_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("reference_node_mesh_generation") == mesh_generation
        and identity.get("quality_mesh_generation") == mesh_generation
        and bool(order_generation)
        and identity.get("reference_node_element_order_generation")
        == order_generation
        and identity.get("quality_element_order_generation") == order_generation
        and element_order >= 2
        and reference_node_count > 0
        and quality_reference_node_count == reference_node_count
        and bool(sampling_rule)
        and identity.get("quality_jacobian_sampling_rule") == sampling_rule
        and len(coordinate_digest) == 64
        and all(character in "0123456789abcdef" for character in coordinate_digest)
        and identity.get("quality_reference_coordinates_sha256")
        == coordinate_digest
    )


def _periodic_hex_node_pair_transform_instance_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("source_node_ids"),
        identity.get("target_node_ids"),
        identity.get("paired_source_node_ids"),
        identity.get("paired_target_node_ids"),
        identity.get("transform_translation_m"),
        identity.get("paired_transform_translation_m"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        source_ids = [int(value) for value in fields[0]]
        target_ids = [int(value) for value in fields[1]]
        paired_source_ids = [int(value) for value in fields[2]]
        paired_target_ids = [int(value) for value in fields[3]]
        translation = [float(value) for value in fields[4]]
        paired_translation = [float(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    instance_generation = str(identity.get("volume_instance_generation") or "")
    digest = str(identity.get("node_pair_transform_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("node_pair_mesh_generation") == mesh_generation
        and bool(instance_generation)
        and identity.get("periodic_transform_volume_instance_generation")
        == instance_generation
        and identity.get("node_pair_volume_instance_generation")
        == instance_generation
        and bool(source_ids)
        and len(source_ids) == len(target_ids)
        and len(set(source_ids)) == len(source_ids)
        and len(set(target_ids)) == len(target_ids)
        and paired_source_ids == source_ids
        and paired_target_ids == target_ids
        and len(translation) == 3
        and all(math.isfinite(value) for value in translation)
        and paired_translation == translation
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("applied_node_pair_transform_sha256") == digest
    )


def _hex_boundary_layer_thickness_surface_normal_generation_ok(
    identity: object,
) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("surface_ids"),
        identity.get("surface_normal_signs"),
        identity.get("applied_surface_ids"),
        identity.get("applied_collapse_direction_signs"),
        identity.get("layer_thickness_m"),
        identity.get("applied_layer_thickness_m"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        surface_ids = [int(value) for value in fields[0]]
        normal_signs = [int(value) for value in fields[1]]
        applied_ids = [int(value) for value in fields[2]]
        applied_signs = [int(value) for value in fields[3]]
        thickness = [float(value) for value in fields[4]]
        applied_thickness = [float(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    geometry_generation = str(identity.get("geometry_generation") or "")
    layer_generation = str(identity.get("boundary_layer_generation") or "")
    digest = str(identity.get("surface_layer_map_sha256") or "")
    return (
        bool(geometry_generation)
        and identity.get("surface_normal_geometry_generation") == geometry_generation
        and identity.get("boundary_layer_geometry_generation") == geometry_generation
        and bool(layer_generation)
        and identity.get("thickness_boundary_layer_generation") == layer_generation
        and identity.get("collapse_direction_boundary_layer_generation")
        == layer_generation
        and bool(surface_ids)
        and len(set(surface_ids)) == len(surface_ids)
        and len(surface_ids) == len(normal_signs) == len(thickness)
        and all(sign in {-1, 1} for sign in normal_signs)
        and applied_ids == surface_ids
        and applied_signs == normal_signs
        and all(math.isfinite(value) and value > 0.0 for value in thickness)
        and applied_thickness == thickness
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("applied_surface_layer_map_sha256") == digest
    )


def _partition_ghost_owner_shared_node_map_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("partition_ids"),
        identity.get("element_ids"),
        identity.get("element_owner_partition_ids"),
        identity.get("ghost_element_ids"),
        identity.get("ghost_owner_partition_ids"),
        identity.get("shared_node_ids"),
        identity.get("shared_node_partition_pairs"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        partition_ids = [int(value) for value in fields[0]]
        element_ids = [int(value) for value in fields[1]]
        owners = [int(value) for value in fields[2]]
        ghost_ids = [int(value) for value in fields[3]]
        ghost_owners = [int(value) for value in fields[4]]
        shared_nodes = [int(value) for value in fields[5]]
        shared_pairs = [[int(value) for value in pair] for pair in fields[6]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("partition_generation") or "")
    digest = str(identity.get("partition_ownership_sha256") or "")
    owner_by_element = dict(zip(element_ids, owners))
    return (
        bool(generation)
        and identity.get("ghost_owner_partition_generation") == generation
        and identity.get("shared_node_partition_generation") == generation
        and bool(partition_ids)
        and len(set(partition_ids)) == len(partition_ids)
        and len(element_ids) == len(owners)
        and len(set(element_ids)) == len(element_ids)
        and all(owner in partition_ids for owner in owners)
        and len(ghost_ids) == len(ghost_owners)
        and all(owner_by_element.get(element) == owner for element, owner in zip(ghost_ids, ghost_owners))
        and len(shared_nodes) == len(shared_pairs)
        and len(set(shared_nodes)) == len(shared_nodes)
        and all(len(pair) == 2 and set(pair).issubset(partition_ids) for pair in shared_pairs)
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("exported_partition_ownership_sha256") == digest
    )


def _exodus_block_namespace_qa_mesh_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("block_ids"),
        identity.get("block_names"),
        identity.get("written_block_ids"),
        identity.get("written_block_names"),
        identity.get("qa_record"),
        identity.get("written_qa_record"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        block_ids = [int(value) for value in fields[0]]
        block_names = [str(value) for value in fields[1]]
        written_ids = [int(value) for value in fields[2]]
        written_names = [str(value) for value in fields[3]]
        qa_record = [str(value) for value in fields[4]]
        written_qa = [str(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    export_generation = str(identity.get("exodus_export_generation") or "")
    digest = str(identity.get("block_namespace_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("block_namespace_mesh_generation") == mesh_generation
        and identity.get("qa_record_mesh_generation") == mesh_generation
        and bool(export_generation)
        and identity.get("writer_export_generation") == export_generation
        and bool(block_ids)
        and len(block_ids) == len(block_names)
        and len(set(block_ids)) == len(block_ids)
        and len(set(block_names)) == len(block_names)
        and all(block_names)
        and written_ids == block_ids
        and written_names == block_names
        and len(qa_record) == 4
        and all(qa_record)
        and written_qa == qa_record
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("written_block_namespace_sha256") == digest
    )


def cubit_conformal_hex_pyramid_tet_interface_gate(
    summary: Mapping[str, object],
    *,
    mapped_volume_id: int = 1,
    transition_volume_id: int = 2,
    min_scaled_jacobian: float = 0.1,
    max_volume_relative_error: float = 1.0e-9,
) -> dict[str, object]:
    """Validate topology, adjacency, quality, and independent volume closure.

    The gate intentionally does not require hexes to outnumber tetrahedra.  A
    minimal conformal transition can contain one mapped hex, one pyramid, and
    several tetrahedra; ownership of the shared quad is the stronger contract.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    threshold = _finite(min_scaled_jacobian, "min_scaled_jacobian")
    tolerance = _finite(max_volume_relative_error, "max_volume_relative_error")
    if threshold <= 0.0 or tolerance < 0.0:
        raise ValueError("quality threshold must be positive and volume tolerance nonnegative")
    mapped_id = str(int(mapped_volume_id))
    transition_id = str(int(transition_volume_id))
    if mapped_id == transition_id:
        raise ValueError("mapped and transition volume IDs must differ")

    totals = _counts(summary.get("element_counts"), "element_counts")
    per_volume_raw = _mapping(summary.get("per_volume_element_counts"), "per_volume_element_counts")
    per_volume = {
        str(volume_id): _counts(row, f"per_volume_element_counts[{volume_id!r}]")
        for volume_id, row in per_volume_raw.items()
    }
    if mapped_id not in per_volume or transition_id not in per_volume:
        raise ValueError("mapped and transition volume inventories are required")

    quality_raw = _mapping(summary.get("quality"), "quality")
    quality_counts: dict[str, int] = {}
    quality_minima: dict[str, float] = {}
    for family in ("hex", "pyramid", "tet"):
        family_row = _mapping(quality_raw.get(family), f"quality.{family}")
        metric = _mapping(family_row.get("scaled_jacobian"), f"quality.{family}.scaled_jacobian")
        quality_counts[family] = int(metric.get("count", -1))
        quality_minima[family] = _finite(metric.get("min"), f"quality.{family}.scaled_jacobian.min")

    interfaces = _rows(summary.get("interface_surfaces"), "interface_surfaces")
    ownership = _rows(summary.get("interface_face_ownership"), "interface_face_ownership")
    interface_face_ids = {
        int(face_id)
        for row in interfaces
        for face_id in (row.get("face_ids") or [])
    }
    ownership_face_ids = {int(row.get("face_id", -1)) for row in ownership}
    all_interface_connectivity = [
        list(connectivity)
        for row in interfaces
        for connectivity in (row.get("face_connectivity") or [])
    ]
    interface_face_incidence_counts = [
        int(row.get("face_incidence_count", 2)) for row in interfaces
    ]

    geometry = _mapping(summary.get("geometry"), "geometry")
    cad_total = _finite(geometry.get("cad_total_volume_m3"), "geometry.cad_total_volume_m3")
    analytic_total = _finite(
        geometry.get("analytic_total_volume_m3"), "geometry.analytic_total_volume_m3"
    )
    gmsh_inventory = _mapping(summary.get("gmsh_inventory"), "gmsh_inventory")
    gmsh_counts = _counts(gmsh_inventory.get("volume_family_counts"), "gmsh_inventory.volume_family_counts")
    gmsh_volume = _mapping(summary.get("gmsh_volume_inventory"), "gmsh_volume_inventory")
    reconstructed_counts = _counts(gmsh_volume.get("family_counts"), "gmsh_volume_inventory.family_counts")
    reconstructed_total = _finite(
        gmsh_volume.get("total_volume_m3"), "gmsh_volume_inventory.total_volume_m3"
    )
    family_volumes = _mapping(gmsh_volume.get("family_volumes_m3"), "gmsh_volume_inventory.family_volumes_m3")
    reconstructed_family_sum = sum(
        _finite(family_volumes.get(family), f"gmsh_volume_inventory.family_volumes_m3.{family}")
        for family in ("hex", "pyramid", "tet")
    )

    denominator = max(abs(cad_total), 1.0e-300)
    analytic_error = abs(cad_total - analytic_total) / max(abs(analytic_total), 1.0e-300)
    reconstructed_error = abs(reconstructed_total - cad_total) / denominator
    reconstructed_sum_error = abs(reconstructed_family_sum - reconstructed_total) / max(
        abs(reconstructed_total), 1.0e-300
    )
    per_volume_sum = {
        family: sum(row[family] for row in per_volume.values()) for family in _VOLUME_FAMILIES
    }
    mapped = per_volume[mapped_id]
    transition = per_volume[transition_id]
    export = _mapping(summary.get("gmsh_export"), "gmsh_export")

    mesh_identity_value = summary.get("mesh_identity")
    quality_identity_value = summary.get("quality_report_identity")
    mesh_quality_identity_present = (
        mesh_identity_value is not None or quality_identity_value is not None
    )
    mesh_identity = (
        mesh_identity_value if isinstance(mesh_identity_value, Mapping) else {}
    )
    quality_identity = (
        quality_identity_value
        if isinstance(quality_identity_value, Mapping)
        else {}
    )
    mesh_quality_identity_ok = not mesh_quality_identity_present or (
        bool(mesh_identity.get("generation"))
        and len(str(mesh_identity.get("sha256") or "")) == 64
        and quality_identity.get("mesh_generation") == mesh_identity.get("generation")
        and quality_identity.get("mesh_sha256") == mesh_identity.get("sha256")
        and len(str(quality_identity.get("report_sha256") or "")) == 64
    )

    boundary_sets_value = summary.get("boundary_sets")
    boundary_sets_present = boundary_sets_value is not None
    boundary_sets_ok = True
    if boundary_sets_present:
        try:
            boundary_sets = _rows(boundary_sets_value, "boundary_sets")
        except ValueError:
            boundary_sets = []
        boundary_sets_ok = (
            bool(boundary_sets)
            and bool(mesh_identity.get("generation"))
            and all(
                bool(row.get("name"))
                and row.get("mesh_generation") == mesh_identity.get("generation")
                and row.get("mesh_sha256") == mesh_identity.get("sha256")
                and bool(list(row.get("entity_ids") or []))
                and len(str(row.get("connectivity_sha256") or "")) == 64
                for row in boundary_sets
            )
        )

    quality_scope_value = summary.get("quality_scope_identity")
    quality_scope_present = quality_scope_value is not None
    quality_scope = (
        quality_scope_value if isinstance(quality_scope_value, Mapping) else {}
    )
    quality_scope_ok = True
    if quality_scope_present:
        mesh_volume_ids = {str(value) for value in per_volume}
        minimum_volume_ids = {
            str(value) for value in (quality_scope.get("minimum_quality_volume_ids") or [])
        }
        histogram_volume_ids = {
            str(value) for value in (quality_scope.get("histogram_volume_ids") or [])
        }
        declared_mesh_volume_ids = {
            str(value) for value in (quality_scope.get("mesh_volume_ids") or [])
        }
        histogram_counts = quality_scope.get("histogram_owned_element_counts")
        quality_scope_ok = (
            isinstance(histogram_counts, Mapping)
            and mesh_volume_ids
            == minimum_volume_ids
            == histogram_volume_ids
            == declared_mesh_volume_ids
            and all(
                int(histogram_counts.get(family, -1)) == totals[family]
                for family in ("hex", "pyramid", "tet")
            )
        )

    partition_value = summary.get("partition_aggregation")
    partition_present = partition_value is not None
    partition = partition_value if isinstance(partition_value, Mapping) else {}
    partition_aggregation_ok = True
    if partition_present:
        try:
            partition_rows = _rows(partition.get("partitions"), "partition_aggregation.partitions")
        except ValueError:
            partition_rows = []
        families = ("hex", "pyramid", "tet")
        owned_sums = {
            family: sum(
                int(
                    (row.get("owned_counts") or {}).get(family, 0)
                    if isinstance(row.get("owned_counts"), Mapping)
                    else 0
                )
                for row in partition_rows
            )
            for family in families
        }
        reported = partition.get("reported_global_owned_counts")
        partition_aggregation_ok = (
            bool(partition_rows)
            and len({row.get("partition_id") for row in partition_rows})
            == len(partition_rows)
            and partition.get("aggregation_policy") == "owned_elements_only"
            and isinstance(reported, Mapping)
            and all(
                int(reported.get(family, -1))
                == owned_sums[family]
                == totals[family]
                for family in families
            )
        )

    signed_jacobian_value = summary.get("signed_jacobian_identity")
    signed_jacobian_present = signed_jacobian_value is not None
    signed_jacobian = (
        signed_jacobian_value
        if isinstance(signed_jacobian_value, Mapping)
        else {}
    )
    signed_jacobian_ok = True
    if signed_jacobian_present:
        try:
            minimum_signed = float(signed_jacobian.get("minimum_signed_jacobian"))
            maximum_signed = float(signed_jacobian.get("maximum_signed_jacobian"))
            sign_changes = int(signed_jacobian.get("interior_sign_change_count"))
        except (TypeError, ValueError):
            minimum_signed = math.nan
            maximum_signed = math.nan
            sign_changes = -1
        signed_jacobian_ok = (
            bool(signed_jacobian.get("mesh_generation"))
            and math.isfinite(minimum_signed)
            and math.isfinite(maximum_signed)
            and minimum_signed > 0.0
            and maximum_signed >= minimum_signed
            and sign_changes == 0
            and signed_jacobian.get("absolute_volume_matches_cad") is True
        )

    coordinate_scale_value = summary.get("coordinate_scale_identity")
    coordinate_scale_present = coordinate_scale_value is not None
    coordinate_scale = (
        coordinate_scale_value
        if isinstance(coordinate_scale_value, Mapping)
        else {}
    )
    coordinate_scale_ok = True
    if coordinate_scale_present:
        try:
            coordinate_to_si = float(coordinate_scale.get("coordinate_scale_to_si"))
            volume_to_si = float(coordinate_scale.get("volume_scale_to_si"))
        except (TypeError, ValueError):
            coordinate_to_si = math.nan
            volume_to_si = math.nan
        coordinate_scale_ok = (
            coordinate_scale.get("source_geometry_unit") == "mm"
            and coordinate_scale.get("export_coordinate_unit") == "m"
            and math.isclose(coordinate_to_si, 0.001, rel_tol=0.0, abs_tol=1.0e-15)
            and math.isclose(
                volume_to_si,
                coordinate_to_si**3,
                rel_tol=1.0e-12,
                abs_tol=1.0e-24,
            )
            and bool(coordinate_scale.get("coordinate_scale_generation"))
            and coordinate_scale.get("volume_scale_generation")
            == coordinate_scale.get("coordinate_scale_generation")
        )

    face_orientation_value = summary.get(
        "high_order_shared_face_orientation_identity"
    )
    face_orientation_present = face_orientation_value is not None
    face_orientation = (
        face_orientation_value if isinstance(face_orientation_value, Mapping) else {}
    )
    try:
        left_face = [int(value) for value in face_orientation.get("left_face_node_ids", [])]
        right_face = [
            int(value) for value in face_orientation.get("right_face_node_ids", [])
        ]
    except (TypeError, ValueError):
        left_face = []
        right_face = []
    reciprocal_face = (
        [
            left_face[0],
            left_face[3],
            left_face[2],
            left_face[1],
            left_face[7],
            left_face[6],
            left_face[5],
            left_face[4],
        ]
        if len(left_face) == 8
        else []
    )
    face_orientation_ok = not face_orientation_present or (
        bool(face_orientation.get("mesh_generation"))
        and face_orientation.get("left_element_generation")
        == face_orientation.get("mesh_generation")
        and face_orientation.get("right_element_generation")
        == face_orientation.get("mesh_generation")
        and len(left_face) == 8
        and len(set(left_face)) == 8
        and set(right_face) == set(left_face)
        and right_face == reciprocal_face
    )

    live_cad_value = summary.get("live_cad_mesh_identity")
    live_cad_present = live_cad_value is not None
    live_cad = live_cad_value if isinstance(live_cad_value, Mapping) else {}
    try:
        live_volume = float(live_cad.get("live_cad_volume"))
        mesh_reference_volume = float(live_cad.get("mesh_reference_cad_volume"))
    except (TypeError, ValueError):
        live_volume = math.nan
        mesh_reference_volume = math.nan
    live_digest = str(live_cad.get("live_cad_sha256") or "")
    live_cad_ok = not live_cad_present or (
        len(live_digest) == 64
        and live_cad.get("mesh_source_cad_sha256") == live_digest
        and bool(live_cad.get("live_cad_generation"))
        and live_cad.get("mesh_source_cad_generation")
        == live_cad.get("live_cad_generation")
        and math.isfinite(live_volume)
        and live_volume > 0.0
        and math.isclose(
            mesh_reference_volume, live_volume, rel_tol=1.0e-6, abs_tol=1.0e-12
        )
    )

    webcut_sideset_value = summary.get("webcut_sideset_topology_identity")
    webcut_sideset_present = webcut_sideset_value is not None
    webcut_sideset = (
        webcut_sideset_value if isinstance(webcut_sideset_value, Mapping) else {}
    )
    try:
        final_geometry_sequence = int(
            webcut_sideset.get("final_geometry_operation_sequence")
        )
        sideset_capture_sequence = int(
            webcut_sideset.get("sideset_capture_after_operation_sequence")
        )
        sideset_surface_ids = [
            int(value) for value in webcut_sideset.get("sideset_surface_ids", [])
        ]
        resolved_surface_ids = [
            int(value) for value in webcut_sideset.get("resolved_surface_ids", [])
        ]
    except (TypeError, ValueError):
        final_geometry_sequence = -1
        sideset_capture_sequence = -2
        sideset_surface_ids = []
        resolved_surface_ids = []
    sideset_connectivity_digest = str(
        webcut_sideset.get("sideset_connectivity_sha256") or ""
    )
    webcut_sideset_ok = not webcut_sideset_present or (
        final_geometry_sequence >= 0
        and sideset_capture_sequence >= final_geometry_sequence
        and bool(webcut_sideset.get("final_geometry_generation"))
        and webcut_sideset.get("sideset_geometry_generation")
        == webcut_sideset.get("final_geometry_generation")
        and bool(webcut_sideset.get("final_topology_generation"))
        and webcut_sideset.get("sideset_topology_generation")
        == webcut_sideset.get("final_topology_generation")
        and bool(sideset_surface_ids)
        and len(set(sideset_surface_ids)) == len(sideset_surface_ids)
        and resolved_surface_ids == sideset_surface_ids
        and len(sideset_connectivity_digest) == 64
        and webcut_sideset.get("resolved_connectivity_sha256")
        == sideset_connectivity_digest
    )

    curved_node_value = summary.get("high_order_curved_node_identity")
    curved_node_present = curved_node_value is not None
    curved_node = curved_node_value if isinstance(curved_node_value, Mapping) else {}
    try:
        export_order = int(curved_node.get("export_order"))
        edge_node_count = int(curved_node.get("high_order_edge_node_count"))
        face_node_count = int(curved_node.get("high_order_face_node_count"))
    except (TypeError, ValueError):
        export_order = 0
        edge_node_count = 0
        face_node_count = 0
    curved_node_ok = not curved_node_present or (
        export_order == 2
        and bool(curved_node.get("final_mesh_generation"))
        and curved_node.get("edge_node_mesh_generation")
        == curved_node.get("final_mesh_generation")
        and curved_node.get("face_node_mesh_generation")
        == curved_node.get("final_mesh_generation")
        and bool(curved_node.get("curving_generation"))
        and curved_node.get("export_curving_generation")
        == curved_node.get("curving_generation")
        and edge_node_count > 0
        and face_node_count > 0
    )

    block_material = summary.get("hex_block_material_topology_identity")
    block_material_ok = block_material is None or (
        isinstance(block_material, Mapping)
        and bool(block_material.get("final_imprint_generation"))
        and bool(block_material.get("final_topology_generation"))
        and block_material.get("block_topology_generation")
        == block_material.get("final_topology_generation")
        and block_material.get("material_assignment_topology_generation")
        == block_material.get("final_topology_generation")
        and list(block_material.get("block_volume_ids") or [])
        == list(block_material.get("material_assignment_volume_ids") or [])
        and len(str(block_material.get("block_material_map_sha256") or "")) == 64
        and block_material.get("resolved_material_map_sha256")
        == block_material.get("block_material_map_sha256")
    )
    transition_orientation = summary.get("pyramid_transition_face_orientation_identity")
    transition_orientation_ok = transition_orientation is None or (
        isinstance(transition_orientation, Mapping)
        and bool(transition_orientation.get("transition_generation"))
        and transition_orientation.get("pyramid_face_generation")
        == transition_orientation.get("transition_generation")
        and transition_orientation.get("hex_face_generation")
        == transition_orientation.get("transition_generation")
        and _opposed_transition_face_orientation_ok(transition_orientation)
    )
    high_order_hex_order = summary.get(
        "high_order_hex_curved_node_ordering_identity"
    )
    high_order_hex_order_ok = high_order_hex_order is None or (
        isinstance(high_order_hex_order, Mapping)
        and bool(high_order_hex_order.get("high_order_mesh_generation"))
        and bool(high_order_hex_order.get("curved_geometry_generation"))
        and high_order_hex_order.get("element_geometry_generation")
        == high_order_hex_order.get("curved_geometry_generation")
        and high_order_hex_order.get("element_type") == "hex20_serendipity"
        and high_order_hex_order.get("export_element_type")
        == high_order_hex_order.get("element_type")
        and high_order_hex_order.get("node_ordering_convention")
        == "cubit_hex20"
        and high_order_hex_order.get("export_node_ordering_convention")
        == high_order_hex_order.get("node_ordering_convention")
        and len(list(high_order_hex_order.get("canonical_node_ids") or [])) == 20
        and list(high_order_hex_order.get("export_node_ids") or [])
        == list(high_order_hex_order.get("canonical_node_ids") or [])
        and len(
            str(high_order_hex_order.get("canonical_node_order_sha256") or "")
        )
        == 64
        and high_order_hex_order.get("export_node_order_sha256")
        == high_order_hex_order.get("canonical_node_order_sha256")
    )
    sideset_normal = summary.get("sideset_outward_normal_merge_identity")
    sideset_face_ids = (
        list(sideset_normal.get("sideset_face_ids") or [])
        if isinstance(sideset_normal, Mapping)
        else []
    )
    outward_signs = (
        list(sideset_normal.get("outward_normal_signs") or [])
        if isinstance(sideset_normal, Mapping)
        else []
    )
    sideset_normal_ok = sideset_normal is None or (
        isinstance(sideset_normal, Mapping)
        and bool(sideset_normal.get("final_merge_generation"))
        and sideset_normal.get("sideset_topology_generation")
        == sideset_normal.get("final_merge_generation")
        and sideset_normal.get("normal_owner_topology_generation")
        == sideset_normal.get("final_merge_generation")
        and bool(sideset_face_ids)
        and len(set(sideset_face_ids)) == len(sideset_face_ids)
        and list(sideset_normal.get("normal_owner_face_ids") or [])
        == sideset_face_ids
        and list(sideset_normal.get("resolved_owner_volume_ids") or [])
        == list(sideset_normal.get("owner_volume_ids") or [])
        and len(outward_signs) == len(sideset_face_ids)
        and all(int(sign) == 1 for sign in outward_signs)
    )

    smoothing_orientation = summary.get(
        "mixed_interface_smoothing_orientation_identity"
    )
    try:
        orientation_products = [
            float(value)
            for value in smoothing_orientation.get(
                "paired_orientation_products", []
            )
        ]
    except (AttributeError, TypeError, ValueError):
        orientation_products = []
    smoothing_orientation_ok = smoothing_orientation is None or (
        isinstance(smoothing_orientation, Mapping)
        and bool(smoothing_orientation.get("smoothing_generation"))
        and smoothing_orientation.get("interface_face_orientation_generation")
        == smoothing_orientation.get("smoothing_generation")
        and bool(smoothing_orientation.get("interface_topology_generation"))
        and smoothing_orientation.get("orientation_topology_generation")
        == smoothing_orientation.get("interface_topology_generation")
        and bool(list(smoothing_orientation.get("hex_interface_face_ids") or []))
        and len(list(smoothing_orientation.get("hex_interface_face_ids") or []))
        == len(list(smoothing_orientation.get("pyramid_interface_face_ids") or []))
        == len(orientation_products)
        and all(product < 0.0 for product in orientation_products)
        and len(str(smoothing_orientation.get("interface_pair_sha256") or ""))
        == 64
        and smoothing_orientation.get("oriented_interface_pair_sha256")
        == smoothing_orientation.get("interface_pair_sha256")
    )

    jacobian_quadrature = summary.get("high_order_jacobian_quadrature_identity")
    try:
        element_order = int(jacobian_quadrature.get("element_order", 0))
        required_exactness = int(
            jacobian_quadrature.get("required_jacobian_exactness_degree", 0)
        )
        actual_exactness = int(
            jacobian_quadrature.get("jacobian_quadrature_exactness_degree", 0)
        )
    except (AttributeError, TypeError, ValueError):
        element_order = required_exactness = actual_exactness = 0
    jacobian_quadrature_ok = jacobian_quadrature is None or (
        isinstance(jacobian_quadrature, Mapping)
        and bool(jacobian_quadrature.get("high_order_mesh_generation"))
        and jacobian_quadrature.get("quality_evaluation_mesh_generation")
        == jacobian_quadrature.get("high_order_mesh_generation")
        and element_order >= 2
        and required_exactness >= 2 * element_order
        and actual_exactness >= required_exactness
        and bool(jacobian_quadrature.get("quadrature_rule_generation"))
        and jacobian_quadrature.get("quality_evaluation_quadrature_generation")
        == jacobian_quadrature.get("quadrature_rule_generation")
        and len(str(jacobian_quadrature.get("element_geometry_sha256") or ""))
        == 64
        and jacobian_quadrature.get("quality_evaluation_geometry_sha256")
        == jacobian_quadrature.get("element_geometry_sha256")
    )

    sweep_correspondence = summary.get(
        "hex_sweep_vertex_correspondence_heal_identity"
    )
    source_vertices = list(
        sweep_correspondence.get("source_vertex_ids") or []
    ) if isinstance(sweep_correspondence, Mapping) else []
    target_vertices = list(
        sweep_correspondence.get("target_vertex_ids") or []
    ) if isinstance(sweep_correspondence, Mapping) else []
    sweep_correspondence_ok = sweep_correspondence is None or (
        isinstance(sweep_correspondence, Mapping)
        and bool(sweep_correspondence.get("geometry_heal_generation"))
        and sweep_correspondence.get("sweep_geometry_heal_generation")
        == sweep_correspondence.get("geometry_heal_generation")
        and sweep_correspondence.get("source_vertex_map_heal_generation")
        == sweep_correspondence.get("geometry_heal_generation")
        and sweep_correspondence.get("target_vertex_map_heal_generation")
        == sweep_correspondence.get("geometry_heal_generation")
        and bool(source_vertices)
        and len(source_vertices) == len(target_vertices)
        and len(set(source_vertices)) == len(source_vertices)
        and len(set(target_vertices)) == len(target_vertices)
        and list(sweep_correspondence.get("sweep_source_vertex_ids") or [])
        == source_vertices
        and list(sweep_correspondence.get("sweep_target_vertex_ids") or [])
        == target_vertices
        and len(str(sweep_correspondence.get("vertex_correspondence_sha256") or ""))
        == 64
        and sweep_correspondence.get("sweep_vertex_correspondence_sha256")
        == sweep_correspondence.get("vertex_correspondence_sha256")
    )

    transition_jacobian = summary.get(
        "transition_jacobian_parent_orientation_identity"
    )
    try:
        signed_jacobian = float(
            transition_jacobian.get("minimum_signed_jacobian")
        )
        absolute_jacobian = float(
            transition_jacobian.get("minimum_absolute_jacobian")
        )
    except (AttributeError, TypeError, ValueError):
        signed_jacobian = absolute_jacobian = math.nan
    transition_jacobian_ok = transition_jacobian is None or (
        isinstance(transition_jacobian, Mapping)
        and bool(transition_jacobian.get("transition_mesh_generation"))
        and transition_jacobian.get("jacobian_mesh_generation")
        == transition_jacobian.get("transition_mesh_generation")
        and transition_jacobian.get("parent_orientation_generation")
        == transition_jacobian.get("transition_mesh_generation")
        and transition_jacobian.get("parent_orientation_convention")
        == "right_handed_positive"
        and transition_jacobian.get("jacobian_orientation_convention")
        == transition_jacobian.get("parent_orientation_convention")
        and math.isfinite(signed_jacobian)
        and signed_jacobian > 0.0
        and math.isclose(
            signed_jacobian,
            absolute_jacobian,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and len(str(transition_jacobian.get("parent_orientation_sha256") or ""))
        == 64
        and transition_jacobian.get("jacobian_parent_orientation_sha256")
        == transition_jacobian.get("parent_orientation_sha256")
    )
    periodic_hex_node_pair_transform_frame_ok = (
        _periodic_hex_node_pair_transform_frame_ok(
            summary.get("periodic_hex_node_pair_transform_frame_identity")
        )
    )
    pyramid_transition_face_diagonal_convention_ok = (
        _pyramid_transition_face_diagonal_convention_ok(
            summary.get("pyramid_transition_face_diagonal_convention_identity")
        )
    )
    hex_sideset_face_ordinal_reorder_ok = (
        _hex_sideset_outward_face_ordinal_volume_reorder_ok(
            summary.get("hex_sideset_outward_face_ordinal_volume_reorder_identity")
        )
    )
    sweep_layer_bias_orientation_ok = (
        _sweep_layer_bias_source_curve_orientation_generation_ok(
            summary.get(
                "sweep_layer_bias_source_curve_orientation_generation_identity"
            )
        )
    )

    checks = {
        "two_distinct_partition_volumes_recorded": set(per_volume) == {mapped_id, transition_id},
        "mapped_volume_is_hex_only": mapped["hex"] > 0
        and all(mapped[family] == 0 for family in ("pyramid", "tet", "wedge")),
        "transition_volume_is_tet_plus_pyramid": transition["tet"] > 0
        and transition["pyramid"] > 0
        and transition["hex"] == 0
        and transition["wedge"] == 0,
        "per_volume_inventory_matches_totals": per_volume_sum == totals,
        "quality_count_matches_topology": all(
            quality_counts[family] == totals[family] for family in ("hex", "pyramid", "tet")
        ),
        "quality_report_matches_current_mesh_generation": mesh_quality_identity_ok,
        "quality_histogram_covers_the_complete_mesh_scope": quality_scope_ok,
        "partition_aggregation_excludes_ghost_elements": partition_aggregation_ok,
        "signed_jacobians_remain_positive_inside_high_order_hexes": (
            signed_jacobian_ok
        ),
        "mesh_coordinates_and_volume_use_one_length_scale": coordinate_scale_ok,
        "high_order_hex_shared_faces_have_reciprocal_orientation": (
            face_orientation_ok
        ),
        "mesh_manifest_matches_live_cad_identity": live_cad_ok,
        "sidesets_follow_final_webcut_topology": webcut_sideset_ok,
        "high_order_export_nodes_match_current_mesh_and_curving_generation": (
            curved_node_ok
        ),
        "hex_block_materials_follow_final_imprint_topology": block_material_ok,
        "pyramid_hex_transition_faces_have_opposed_orientation": (
            transition_orientation_ok
        ),
        "curved_high_order_hex_uses_canonical_node_ordering": (
            high_order_hex_order_ok
        ),
        "merged_sideset_normals_follow_final_topology_owners": (
            sideset_normal_ok
        ),
        "smoothed_hex_pyramid_interface_faces_keep_opposed_orientation": (
            smoothing_orientation_ok
        ),
        "high_order_jacobian_uses_current_sufficient_quadrature": (
            jacobian_quadrature_ok
        ),
        "hex_sweep_uses_post_heal_vertex_correspondence": (
            sweep_correspondence_ok
        ),
        "transition_jacobian_uses_current_parent_orientation": (
            transition_jacobian_ok
        ),
        "periodic_hex_node_pairs_use_current_transform_frame": (
            periodic_hex_node_pair_transform_frame_ok
        ),
        "pyramid_transition_neighbors_share_one_face_diagonal_convention": (
            pyramid_transition_face_diagonal_convention_ok
        ),
        "hex_sideset_face_ordinals_and_normals_follow_connectivity_reorder": (
            hex_sideset_face_ordinal_reorder_ok
        ),
        "biased_sweep_layers_follow_current_source_curve_orientation": (
            sweep_layer_bias_orientation_ok
        ),
        "periodic_hex_pairs_follow_current_volume_instance_transform": (
            _periodic_hex_node_pair_transform_instance_generation_ok(
                summary.get("periodic_hex_node_pair_transform_instance_generation_identity")
            )
        ),
        "hex_boundary_layers_follow_current_healed_surface_normals": (
            _hex_boundary_layer_thickness_surface_normal_generation_ok(
                summary.get("hex_boundary_layer_thickness_surface_normal_generation_identity")
            )
        ),
        "boundary_sets_match_current_mesh_generation": boundary_sets_ok,
        "all_volume_families_above_quality_threshold": all(
            quality_minima[family] >= threshold for family in ("hex", "pyramid", "tet")
        ),
        "single_shared_interface_between_partitions": len(interfaces) == 1
        and {int(value) for value in (interfaces[0].get("adjacent_volumes") or [])}
        == {int(mapped_id), int(transition_id)},
        "interface_faces_are_quads": bool(interface_face_ids)
        and len(all_interface_connectivity) == len(interface_face_ids)
        and all(len(connectivity) == 4 for connectivity in all_interface_connectivity),
        "each_interface_quad_has_one_hex_and_one_pyramid_owner": bool(ownership)
        and ownership_face_ids == interface_face_ids
        and all(
            int(row.get("node_count", 0)) == 4
            and len(list(row.get("hex_owners") or [])) == 1
            and len(list(row.get("pyramid_owners") or [])) == 1
            for row in ownership
        ),
        "interface_quads_are_two_sided_manifold": bool(interfaces)
        and all(count == 2 for count in interface_face_incidence_counts)
        and all(
            not list(row.get("tet_owners") or [])
            and not list(row.get("wedge_owners") or [])
            and not list(row.get("other_owners") or [])
            for row in ownership
        ),
        "every_pyramid_serves_the_transition": int(summary.get("matched_pyramid_count", -1))
        == totals["pyramid"],
        "gmsh_export_is_fresh_nonempty_digest": int(export.get("bytes", 0)) > 0
        and len(str(export.get("sha256", ""))) == 64,
        "gmsh_is_ascii_v41_with_valid_connectivity": gmsh_inventory.get("status") == "ok"
        and str(gmsh_inventory.get("mesh_format", "")) == "4.1"
        and gmsh_inventory.get("binary") is False
        and not list(gmsh_inventory.get("connectivity_mismatches") or []),
        "gmsh_volume_families_match_live_inventory": gmsh_counts == totals
        and reconstructed_counts == totals,
        "cad_volume_matches_analytic": analytic_error <= tolerance,
        "independent_gmsh_volume_sum_matches_cad": reconstructed_error <= tolerance
        and reconstructed_sum_error <= tolerance
        and geometry.get("element_volume_source")
        == "independent_gmsh_v41_coordinate_reconstruction",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_conformal_hex_pyramid_tet_interface_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "element_counts": totals,
        "quality_minima": quality_minima,
        "interface_face_count": len(interface_face_ids),
        "matched_pyramid_count": int(summary.get("matched_pyramid_count", -1)),
        "cad_volume_relative_error": analytic_error,
        "gmsh_reconstructed_volume_relative_error": reconstructed_error,
        "notes": [
            "Do not require hex dominance: a valid minimal transition may contain more tetrahedra than hexes.",
            "A conformal transition is proven by quad ownership on both sides, not by family counts alone.",
            "Reconstruct volume from Gmsh 4.1 coordinates so CAD closure is independent of Cubit's quality API.",
        ],
    }


def cubit_mixed_transition_source_gate(
    summary: Mapping[str, object],
    *,
    mapped_volume_id: int = 1,
    transition_volume_id: int = 2,
) -> dict[str, object]:
    """Gate the source journal and a classified synchronous headless replay."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    commands_raw = summary.get("source_commands")
    if isinstance(commands_raw, (str, bytes)) or not isinstance(commands_raw, Sequence):
        raise ValueError("source_commands must be a sequence")
    commands = [str(command).strip().lower() for command in commands_raw if str(command).strip()]
    command_text = "\n".join(commands)

    def index(fragment: str) -> int:
        return next((offset for offset, command in enumerate(commands) if fragment in command), -1)

    source_sha = str(summary.get("source_sha256", "")).lower()
    quality_probe = _mapping(summary.get("quality_probe"), "quality_probe")
    process = _mapping(summary.get("process"), "process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    export_artifacts_value = summary.get("export_artifacts")
    export_artifacts_present = isinstance(export_artifacts_value, Mapping)
    export_artifacts_ok = True
    if export_artifacts_present:
        required_raw = export_artifacts_value.get("required")
        if isinstance(required_raw, (str, bytes)) or not isinstance(
            required_raw, Sequence
        ):
            raise ValueError("export_artifacts.required must be a sequence")
        required = [str(name) for name in required_raw]
        artifact_rows = _rows(
            export_artifacts_value.get("artifacts"), "export_artifacts.artifacts"
        )
        artifact_names = [str(row.get("name", "")) for row in artifact_rows]
        required_names = set(required)
        artifact_by_name = dict(zip(artifact_names, artifact_rows, strict=True))
        export_artifacts_ok = (
            bool(required_names)
            and all(required)
            and all(artifact_names)
            and len(required_names) == len(required)
            and len(artifact_by_name) == len(artifact_rows)
            and set(artifact_by_name) == required_names
        )
        export_artifacts_ok = export_artifacts_ok and all(
            row.get("fresh") is True
            and int(row.get("bytes", 0)) > 0
            and len(str(row.get("sha256", ""))) == 64
            for row in artifact_by_name.values()
        )

    replay_identity_value = summary.get("replay_identity")
    replay_identity_present = isinstance(replay_identity_value, Mapping)
    replay_identity_ok = True
    if replay_identity_present:
        pinned_journal = str(replay_identity_value.get("pinned_journal_sha256", ""))
        pinned_model = str(replay_identity_value.get("pinned_source_model_sha256", ""))
        replayed_journal = str(
            replay_identity_value.get("replayed_journal_sha256", "")
        )
        replayed_model = str(
            replay_identity_value.get("replayed_source_model_sha256", "")
        )
        replay_identity_ok = (
            all(
                len(value) == 64
                and all(character.lower() in "0123456789abcdef" for character in value)
                for value in (
                    pinned_journal,
                    pinned_model,
                    replayed_journal,
                    replayed_model,
                )
            )
            and pinned_journal == replayed_journal
            and pinned_model == replayed_model
            and pinned_journal == str(summary.get("source_sha256", ""))
        )

    export_manifest_value = summary.get("export_manifest")
    export_manifest_present = export_manifest_value is not None
    export_manifest_ok = True
    if export_manifest_present:
        export_manifest = (
            export_manifest_value
            if isinstance(export_manifest_value, Mapping)
            else {}
        )
        try:
            manifest_artifacts = _rows(
                export_manifest.get("artifacts"), "export_manifest.artifacts"
            )
        except ValueError:
            manifest_artifacts = []
        invocation_id = str(export_manifest.get("invocation_id") or "")
        model_generation = str(export_manifest.get("model_generation") or "")
        names = [str(row.get("name") or "") for row in manifest_artifacts]
        export_manifest_ok = (
            bool(invocation_id)
            and bool(model_generation)
            and bool(manifest_artifacts)
            and len(set(names)) == len(names)
            and all(names)
            and all(
                row.get("invocation_id") == invocation_id
                and row.get("model_generation") == model_generation
                and len(str(row.get("sha256") or "")) == 64
                for row in manifest_artifacts
            )
        )

    batch_invocation_value = summary.get("batch_invocation")
    batch_invocation_present = batch_invocation_value is not None
    batch_invocation_ok = True
    if batch_invocation_present:
        batch_invocation = (
            batch_invocation_value
            if isinstance(batch_invocation_value, Mapping)
            else {}
        )
        batch_log_value = batch_invocation.get("log")
        batch_log = batch_log_value if isinstance(batch_log_value, Mapping) else {}
        invocation_id = str(batch_invocation.get("invocation_id") or "")
        process_start = str(batch_invocation.get("process_start_utc") or "")
        batch_invocation_ok = (
            bool(invocation_id)
            and bool(process_start)
            and batch_log.get("invocation_id") == invocation_id
            and batch_log.get("process_start_utc") == process_start
            and len(str(batch_log.get("sha256") or "")) == 64
            and batch_invocation.get("exports_invocation_id") == invocation_id
        )

    operation_dag_value = summary.get("operation_dag_identity")
    operation_dag_present = operation_dag_value is not None
    operation_dag = (
        operation_dag_value if isinstance(operation_dag_value, Mapping) else {}
    )
    export_manifest_for_dag = (
        export_manifest_value if isinstance(export_manifest_value, Mapping) else {}
    )
    operation_dag_ok = True
    if operation_dag_present:
        try:
            final_sequence = int(operation_dag.get("final_operation_sequence"))
            export_sequence = int(operation_dag.get("export_after_operation_sequence"))
        except (TypeError, ValueError):
            final_sequence = -1
            export_sequence = -2
        operation_dag_ok = (
            bool(operation_dag.get("final_model_generation"))
            and operation_dag.get("export_model_generation")
            == operation_dag.get("final_model_generation")
            and export_manifest_for_dag.get("model_generation")
            == operation_dag.get("final_model_generation")
            and final_sequence >= 0
            and export_sequence >= final_sequence
        )

    length_scale_value = summary.get("length_scale_identity")
    length_scale_present = length_scale_value is not None
    length_scale = (
        length_scale_value if isinstance(length_scale_value, Mapping) else {}
    )
    length_scale_ok = True
    if length_scale_present:
        try:
            declared_scale = float(length_scale.get("declared_source_to_export_scale"))
            effective_scale = float(length_scale.get("effective_scale"))
        except (TypeError, ValueError):
            declared_scale = math.nan
            effective_scale = math.nan
        length_scale_ok = (
            length_scale.get("source_geometry_unit") == "mm"
            and length_scale.get("export_geometry_unit") == "m"
            and list(length_scale.get("scale_application_stages") or [])
            == ["source-command"]
            and math.isfinite(declared_scale)
            and math.isfinite(effective_scale)
            and declared_scale == 0.001
            and effective_scale == declared_scale
        )

    exodus_connectivity_value = summary.get("exodus_connectivity_identity")
    exodus_connectivity_present = exodus_connectivity_value is not None
    exodus_connectivity = (
        exodus_connectivity_value
        if isinstance(exodus_connectivity_value, Mapping)
        else {}
    )
    connectivity_digest = str(
        exodus_connectivity.get("permuted_connectivity_sha256") or ""
    )
    exodus_connectivity_ok = not exodus_connectivity_present or (
        bool(exodus_connectivity.get("connectivity_permutation_generation"))
        and exodus_connectivity.get("sideset_face_ordinal_generation")
        == exodus_connectivity.get("connectivity_permutation_generation")
        and len(connectivity_digest) == 64
        and exodus_connectivity.get("sideset_connectivity_sha256")
        == connectivity_digest
        and exodus_connectivity.get("target_ordering")
        == "solver-target-ordering-v1"
    )

    quality_generation_value = summary.get("quality_report_generation_identity")
    quality_generation_present = quality_generation_value is not None
    quality_generation = (
        quality_generation_value
        if isinstance(quality_generation_value, Mapping)
        else {}
    )
    try:
        final_smoothing_sequence = int(
            quality_generation.get("final_smoothing_sequence")
        )
        report_sequence = int(
            quality_generation.get("quality_report_after_operation_sequence")
        )
    except (TypeError, ValueError):
        final_smoothing_sequence = -1
        report_sequence = -2
    quality_generation_ok = not quality_generation_present or (
        bool(quality_generation.get("final_mesh_generation"))
        and quality_generation.get("quality_report_mesh_generation")
        == quality_generation.get("final_mesh_generation")
        and final_smoothing_sequence >= 0
        and report_sequence >= final_smoothing_sequence
    )

    block_material_value = summary.get("block_material_map_identity")
    block_material_present = block_material_value is not None
    block_material = (
        block_material_value if isinstance(block_material_value, Mapping) else {}
    )
    block_digest = str(block_material.get("block_table_sha256") or "")
    unmapped_blocks = block_material.get("unmapped_block_ids")
    block_material_ok = not block_material_present or (
        bool(block_material.get("final_mesh_generation"))
        and block_material.get("material_map_mesh_generation")
        == block_material.get("final_mesh_generation")
        and len(block_digest) == 64
        and block_material.get("material_map_block_table_sha256") == block_digest
        and isinstance(unmapped_blocks, list)
        and not unmapped_blocks
    )

    sculpt_value = summary.get("parallel_sculpt_completion_identity")
    sculpt_present = sculpt_value is not None
    sculpt = sculpt_value if isinstance(sculpt_value, Mapping) else {}
    try:
        expected_ranks = int(sculpt.get("expected_rank_count"))
        finalized_ranks = [int(value) for value in sculpt.get("finalized_rank_ids", [])]
    except (TypeError, ValueError):
        expected_ranks = 0
        finalized_ranks = []
    rank_digests = sculpt.get("rank_artifact_sha256")
    sculpt_ok = not sculpt_present or (
        expected_ranks > 0
        and finalized_ranks == list(range(expected_ranks))
        and isinstance(rank_digests, list)
        and len(rank_digests) == expected_ranks
        and all(len(str(value)) == 64 for value in rank_digests)
        and bool(sculpt.get("rank_manifest_generation"))
        and sculpt.get("global_aggregation_generation")
        == sculpt.get("rank_manifest_generation")
    )

    sideset_manifest_value = summary.get("headless_sideset_manifest_identity")
    sideset_manifest_present = sideset_manifest_value is not None
    sideset_manifest = (
        sideset_manifest_value
        if isinstance(sideset_manifest_value, Mapping)
        else {}
    )
    try:
        final_webcut_sequence = int(
            sideset_manifest.get("final_webcut_operation_sequence")
        )
        manifest_capture_sequence = int(
            sideset_manifest.get("manifest_capture_after_operation_sequence")
        )
        manifest_surface_ids = [
            int(value) for value in sideset_manifest.get("manifest_surface_ids", [])
        ]
        live_surface_ids = [
            int(value) for value in sideset_manifest.get("live_surface_ids", [])
        ]
    except (TypeError, ValueError):
        final_webcut_sequence = -1
        manifest_capture_sequence = -2
        manifest_surface_ids = []
        live_surface_ids = []
    manifest_connectivity_digest = str(
        sideset_manifest.get("manifest_connectivity_sha256") or ""
    )
    sideset_manifest_ok = not sideset_manifest_present or (
        bool(sideset_manifest.get("batch_invocation_id"))
        and sideset_manifest.get("manifest_invocation_id")
        == sideset_manifest.get("batch_invocation_id")
        and final_webcut_sequence >= 0
        and manifest_capture_sequence >= final_webcut_sequence
        and bool(sideset_manifest.get("final_topology_generation"))
        and sideset_manifest.get("manifest_topology_generation")
        == sideset_manifest.get("final_topology_generation")
        and bool(manifest_surface_ids)
        and len(set(manifest_surface_ids)) == len(manifest_surface_ids)
        and live_surface_ids == manifest_surface_ids
        and len(manifest_connectivity_digest) == 64
        and sideset_manifest.get("live_connectivity_sha256")
        == manifest_connectivity_digest
    )

    high_order_export_value = summary.get("netgen_high_order_export_identity")
    high_order_export_present = high_order_export_value is not None
    high_order_export = (
        high_order_export_value
        if isinstance(high_order_export_value, Mapping)
        else {}
    )
    try:
        netgen_export_order = int(high_order_export.get("export_order"))
        higher_order_node_count = int(
            high_order_export.get("higher_order_node_count")
        )
    except (TypeError, ValueError):
        netgen_export_order = 0
        higher_order_node_count = 0
    netgen_export_digest = str(
        high_order_export.get("netgen_export_sha256") or ""
    )
    high_order_export_ok = not high_order_export_present or (
        netgen_export_order == 2
        and bool(high_order_export.get("final_mesh_generation"))
        and high_order_export.get("higher_order_node_mesh_generation")
        == high_order_export.get("final_mesh_generation")
        and bool(high_order_export.get("active_model_generation"))
        and high_order_export.get("export_model_generation")
        == high_order_export.get("active_model_generation")
        and higher_order_node_count > 0
        and len(netgen_export_digest) == 64
    )
    headless_block = summary.get("headless_block_material_manifest_identity")
    headless_block_ok = headless_block is None or (
        isinstance(headless_block, Mapping)
        and bool(headless_block.get("batch_invocation_id"))
        and headless_block.get("manifest_invocation_id")
        == headless_block.get("batch_invocation_id")
        and bool(headless_block.get("final_imprint_generation"))
        and bool(headless_block.get("active_topology_generation"))
        and headless_block.get("manifest_topology_generation")
        == headless_block.get("active_topology_generation")
        and list(headless_block.get("manifest_volume_ids") or [])
        == list(headless_block.get("live_volume_ids") or [])
        and len(str(headless_block.get("manifest_material_map_sha256") or "")) == 64
        and headless_block.get("live_material_map_sha256")
        == headless_block.get("manifest_material_map_sha256")
    )
    export_orientation = summary.get("mesh_export_transition_orientation_identity")
    export_orientation_ok = export_orientation is None or (
        isinstance(export_orientation, Mapping)
        and bool(export_orientation.get("export_generation"))
        and export_orientation.get("pyramid_face_export_generation")
        == export_orientation.get("export_generation")
        and export_orientation.get("hex_face_export_generation")
        == export_orientation.get("export_generation")
        and _opposed_transition_face_orientation_ok(export_orientation)
    )
    journal_id_map = summary.get("journal_entity_id_map_reset_identity")
    requested_entity_ids = (
        list(journal_id_map.get("requested_entity_ids") or [])
        if isinstance(journal_id_map, Mapping)
        else []
    )
    journal_id_map_ok = journal_id_map is None or (
        isinstance(journal_id_map, Mapping)
        and bool(journal_id_map.get("reset_generation"))
        and journal_id_map.get("journal_replay_reset_generation")
        == journal_id_map.get("reset_generation")
        and journal_id_map.get("entity_id_map_reset_generation")
        == journal_id_map.get("reset_generation")
        and bool(requested_entity_ids)
        and len(list(journal_id_map.get("entity_kinds") or []))
        == len(requested_entity_ids)
        and list(journal_id_map.get("resolved_entity_ids") or [])
        == requested_entity_ids
        and len(str(journal_id_map.get("entity_id_map_sha256") or "")) == 64
        and journal_id_map.get("resolved_entity_id_map_sha256")
        == journal_id_map.get("entity_id_map_sha256")
    )
    exodus_id_width = summary.get("exodus_entity_id_width_identity")
    try:
        declared_id_width = int(
            exodus_id_width.get("declared_entity_id_width_bits", 0)
        ) if isinstance(exodus_id_width, Mapping) else 0
        decoder_id_width = int(
            exodus_id_width.get("decoder_entity_id_width_bits", 0)
        ) if isinstance(exodus_id_width, Mapping) else 0
        maximum_entity_id = int(
            exodus_id_width.get("maximum_entity_id", -1)
        ) if isinstance(exodus_id_width, Mapping) else -1
        decoded_maximum_entity_id = int(
            exodus_id_width.get("decoded_maximum_entity_id", -1)
        ) if isinstance(exodus_id_width, Mapping) else -1
    except (TypeError, ValueError):
        declared_id_width = decoder_id_width = 0
        maximum_entity_id = decoded_maximum_entity_id = -1
    exodus_id_width_ok = exodus_id_width is None or (
        isinstance(exodus_id_width, Mapping)
        and bool(exodus_id_width.get("export_generation"))
        and exodus_id_width.get("decoder_export_generation")
        == exodus_id_width.get("export_generation")
        and declared_id_width == 64
        and decoder_id_width == declared_id_width
        and exodus_id_width.get("integer_storage_type") == "int64"
        and maximum_entity_id > 2**31 - 1
        and decoded_maximum_entity_id == maximum_entity_id
        and len(str(exodus_id_width.get("entity_id_stream_sha256") or "")) == 64
        and exodus_id_width.get("decoded_entity_id_stream_sha256")
        == exodus_id_width.get("entity_id_stream_sha256")
    )
    tolerance_identity = summary.get("imprint_merge_tolerance_unit_identity")
    unit_scales = {"m": 1.0, "mm": 1.0e-3, "um": 1.0e-6}
    try:
        tolerance_si = float(tolerance_identity.get("tolerance_si_m"))
        imprint_tolerance_si = float(
            tolerance_identity.get("imprint_tolerance_value")
        ) * unit_scales[str(tolerance_identity.get("imprint_tolerance_unit"))]
        merge_tolerance_si = float(
            tolerance_identity.get("merge_tolerance_value")
        ) * unit_scales[str(tolerance_identity.get("merge_tolerance_unit"))]
    except (AttributeError, KeyError, TypeError, ValueError):
        tolerance_si = imprint_tolerance_si = merge_tolerance_si = math.nan
    tolerance_scale = max(abs(tolerance_si), 1.0e-30)
    tolerance_identity_ok = tolerance_identity is None or (
        isinstance(tolerance_identity, Mapping)
        and bool(tolerance_identity.get("geometry_generation"))
        and tolerance_identity.get("imprint_geometry_generation")
        == tolerance_identity.get("geometry_generation")
        and tolerance_identity.get("merge_geometry_generation")
        == tolerance_identity.get("geometry_generation")
        and bool(tolerance_identity.get("tolerance_generation"))
        and tolerance_identity.get("imprint_tolerance_generation")
        == tolerance_identity.get("tolerance_generation")
        and tolerance_identity.get("merge_tolerance_generation")
        == tolerance_identity.get("tolerance_generation")
        and tolerance_identity.get("model_length_unit") in unit_scales
        and math.isfinite(tolerance_si)
        and tolerance_si > 0.0
        and abs(imprint_tolerance_si - tolerance_si) / tolerance_scale <= 1.0e-12
        and abs(merge_tolerance_si - tolerance_si) / tolerance_scale <= 1.0e-12
    )

    exodus_renumber = summary.get("exodus_block_sideset_renumber_identity")
    exodus_renumber_ok = exodus_renumber is None or (
        isinstance(exodus_renumber, Mapping)
        and bool(exodus_renumber.get("renumber_generation"))
        and exodus_renumber.get("block_map_generation")
        == exodus_renumber.get("renumber_generation")
        and exodus_renumber.get("sideset_map_generation")
        == exodus_renumber.get("renumber_generation")
        and bool(list(exodus_renumber.get("block_ids") or []))
        and list(exodus_renumber.get("exported_block_ids") or [])
        == list(exodus_renumber.get("block_ids") or [])
        and bool(list(exodus_renumber.get("sideset_ids") or []))
        and list(exodus_renumber.get("exported_sideset_ids") or [])
        == list(exodus_renumber.get("sideset_ids") or [])
        and len(str(exodus_renumber.get("entity_map_sha256") or "")) == 64
        and exodus_renumber.get("exported_entity_map_sha256")
        == exodus_renumber.get("entity_map_sha256")
    )
    journal_imprint = summary.get("journal_entity_id_imprint_identity")
    journal_volume_ids = list(
        journal_imprint.get("journal_volume_ids") or []
    ) if isinstance(journal_imprint, Mapping) else []
    journal_surface_ids = list(
        journal_imprint.get("journal_surface_ids") or []
    ) if isinstance(journal_imprint, Mapping) else []
    journal_imprint_ok = journal_imprint is None or (
        isinstance(journal_imprint, Mapping)
        and bool(journal_imprint.get("imprint_generation"))
        and journal_imprint.get("journal_entity_generation")
        == journal_imprint.get("imprint_generation")
        and journal_imprint.get("resolved_entity_generation")
        == journal_imprint.get("imprint_generation")
        and bool(journal_volume_ids)
        and len(set(journal_volume_ids)) == len(journal_volume_ids)
        and list(journal_imprint.get("resolved_volume_ids") or [])
        == journal_volume_ids
        and bool(journal_surface_ids)
        and len(set(journal_surface_ids)) == len(journal_surface_ids)
        and list(journal_imprint.get("resolved_surface_ids") or [])
        == journal_surface_ids
        and len(str(journal_imprint.get("entity_table_sha256") or "")) == 64
        and journal_imprint.get("resolved_entity_table_sha256")
        == journal_imprint.get("entity_table_sha256")
    )

    sideset_normal = summary.get(
        "exodus_sideset_outward_normal_topology_identity"
    )
    sideset_ids = list(
        sideset_normal.get("sideset_ids") or []
    ) if isinstance(sideset_normal, Mapping) else []
    sideset_normal_ok = sideset_normal is None or (
        isinstance(sideset_normal, Mapping)
        and bool(sideset_normal.get("topology_generation"))
        and sideset_normal.get("sideset_map_topology_generation")
        == sideset_normal.get("topology_generation")
        and sideset_normal.get("normal_ownership_topology_generation")
        == sideset_normal.get("topology_generation")
        and bool(sideset_ids)
        and len(set(sideset_ids)) == len(sideset_ids)
        and list(sideset_normal.get("normal_ownership_sideset_ids") or [])
        == sideset_ids
        and sideset_normal.get("normal_orientation") == "outward"
        and sideset_normal.get("exported_normal_orientation") == "outward"
        and len(str(sideset_normal.get("normal_ownership_sha256") or "")) == 64
        and sideset_normal.get("exported_normal_ownership_sha256")
        == sideset_normal.get("normal_ownership_sha256")
    )
    block_attribute_material_id_merge_ok = _block_attribute_material_id_merge_ok(
        summary.get("block_attribute_material_id_merge_identity")
    )
    high_order_exodus_node_permutation_export_order_ok = (
        _high_order_exodus_node_permutation_export_order_ok(
            summary.get("high_order_exodus_node_permutation_export_order_identity")
        )
    )
    exodus_sideset_topology_generation_ok = (
        _exodus_sideset_element_face_topology_generation_ok(
            summary.get("exodus_sideset_element_face_topology_generation_identity")
        )
    )
    high_order_quality_reference_generation_ok = (
        _high_order_quality_reference_coordinate_generation_ok(
            summary.get("high_order_quality_reference_coordinate_generation_identity")
        )
    )
    public_gate = cubit_conformal_hex_pyramid_tet_interface_gate(
        summary,
        mapped_volume_id=mapped_volume_id,
        transition_volume_id=transition_volume_id,
    )
    block_families = ("hex", "pyramid", "tet", "tri", "face")
    ordered = [
        index("brick x 2 y 1 z 1"),
        index("webcut volume 1"),
        index("volume 1 scheme map"),
        index("mesh volume 1"),
        index("volume 2 scheme tetmesh"),
        index("mesh volume 2"),
    ]
    exit_code = int(process.get("exit_code", -1))
    unexpected = list(process.get("unexpected_error_lines") or [])
    known_diagnostics_only = process.get("known_headless_diagnostics_only") is True
    exit_explained = exit_code == 0 or (
        exit_code > 0
        and known_diagnostics_only
        and not unexpected
        and process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0
        and public_gate["status"] == "ok"
    )
    checks = {
        "source_native_journal_and_sha256_recorded": str(summary.get("source_kind", "")).startswith(
            "source_native_"
        )
        and Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou"
        and len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha),
        "source_builds_and_webcuts_brick": ordered[0] >= 0
        and ordered[1] > ordered[0]
        and "xplane" in commands[ordered[1]],
        "source_meshes_map_before_tetmesh": all(offset >= 0 for offset in ordered)
        and ordered == sorted(ordered),
        "source_registers_all_mixed_and_surface_families": all(
            f"add {family} all" in command_text for family in block_families
        ),
        "source_applies_meter_scale": any(
            "volume all scale 0.001" in command for command in commands
        ),
        "headless_batch_without_gui_daemon": summary.get("execution_mode")
        == "headless_combined_journal_then_python_inventory"
        and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or []))
        and summary.get("gui_daemon_enabled") is False,
        "cubit_version_recorded": bool(str(summary.get("version", "")).strip()),
        "unsupported_aggregate_quality_probe_is_diagnosed": quality_probe.get("command_supported")
        is False
        and "unknown metric name volume" in str(quality_probe.get("diagnostic", "")).lower()
        and quality_probe.get("failure_interpretation") == "unsupported_api_not_zero_quality",
        "per_element_scaled_jacobian_fallback_recorded": quality_probe.get("fallback")
        == "per_element_scaled_jacobian_by_family"
        and {str(value) for value in (quality_probe.get("families") or [])}
        == {"hex", "pyramid", "tet"},
        "nonzero_exit_is_semantically_classified": exit_explained,
        "fresh_artifact_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "all_required_export_artifacts_are_fresh": export_artifacts_ok,
        "export_manifest_uses_one_model_and_invocation_generation": export_manifest_ok,
        "batch_log_and_exports_share_invocation_identity": batch_invocation_ok,
        "exports_follow_the_final_geometry_operation": operation_dag_ok,
        "length_scale_is_applied_exactly_once": length_scale_ok,
        "exodus_sidesets_follow_connectivity_permutation": exodus_connectivity_ok,
        "quality_report_follows_final_smoothing_generation": quality_generation_ok,
        "block_material_map_matches_final_mesh_generation": block_material_ok,
        "parallel_sculpt_waits_for_every_rank_artifact": sculpt_ok,
        "headless_sideset_manifest_follows_final_webcut_topology": (
            sideset_manifest_ok
        ),
        "netgen_high_order_nodes_match_current_mesh_and_model_generation": (
            high_order_export_ok
        ),
        "headless_block_material_manifest_follows_final_imprint": headless_block_ok,
        "mesh_export_transition_faces_have_opposed_orientation": export_orientation_ok,
        "journal_entity_ids_follow_current_reset_generation": journal_id_map_ok,
        "exodus_decoder_preserves_declared_64bit_entity_ids": exodus_id_width_ok,
        "imprint_merge_tolerances_share_one_physical_length_basis": (
            tolerance_identity_ok
        ),
        "exodus_block_sideset_maps_follow_current_renumber_generation": (
            exodus_renumber_ok
        ),
        "journal_entity_ids_follow_final_imprint_generation": journal_imprint_ok,
        "exodus_sideset_normals_follow_current_topology": sideset_normal_ok,
        "block_material_attributes_follow_final_merge_generation": (
            block_attribute_material_id_merge_ok
        ),
        "high_order_exodus_nodes_use_current_export_order_permutation": (
            high_order_exodus_node_permutation_export_order_ok
        ),
        "exodus_sideset_ordinals_follow_current_mesh_and_export_topology": (
            exodus_sideset_topology_generation_ok
        ),
        "high_order_quality_uses_current_reference_coordinates_and_order": (
            high_order_quality_reference_generation_ok
        ),
        "partition_ghosts_and_shared_nodes_use_current_owner_map": (
            _partition_ghost_owner_shared_node_map_generation_ok(
                summary.get("partition_ghost_element_owner_shared_node_map_identity")
            )
        ),
        "exodus_blocks_and_qa_use_current_mesh_namespace": (
            _exodus_block_namespace_qa_mesh_generation_ok(
                summary.get("exodus_block_id_namespace_qa_record_mesh_generation_identity")
            )
        ),
        "journal_and_source_model_identity_match_replay": replay_identity_ok,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(_finite(value, f"timing_breakdown_s.{name}") >= 0.0 for name, value in timing.items()),
        "independent_interface_gate_passed": public_gate["status"] == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    warnings = []
    if not export_artifacts_present:
        warnings.append("per_artifact_export_freshness_not_recorded")
    if not replay_identity_present:
        warnings.append("journal_model_replay_identity_not_recorded")
    return {
        "policy": "cubit_mixed_transition_source_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "warnings": warnings,
        "source_journal": Path(str(summary.get("source_journal", ""))).name,
        "process_exit_code": exit_code,
        "public_gate_status": public_gate["status"],
        "notes": [
            "Replay source commands synchronously and headlessly; never infer completion from a queued GUI playback.",
            "An unsupported aggregate quality query is an API diagnostic, not evidence of zero elements or zero quality.",
            "Permit a nonzero batch exit only with classified diagnostics, a fresh passing artifact, and no leaked process.",
        ],
    }
