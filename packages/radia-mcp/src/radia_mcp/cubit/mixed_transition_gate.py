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


def _valid_sha256(value: object) -> bool:
    digest = str(value or "")
    return len(digest) == 64 and all(character in "0123456789abcdef" for character in digest)


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


def _high_order_hex_jacobian_node_order_scale_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("corner_node_ids"),
        identity.get("jacobian_corner_node_ids"),
        identity.get("reference_corner_coordinates"),
        identity.get("jacobian_reference_corner_coordinates"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        node_ids = [int(value) for value in fields[0]]
        jacobian_node_ids = [int(value) for value in fields[1]]
        reference = [[float(value) for value in row] for row in fields[2]]
        jacobian_reference = [[float(value) for value in row] for row in fields[3]]
        order = int(identity.get("element_order"))
        jacobian_order = int(identity.get("jacobian_element_order"))
        scale = float(identity.get("coordinate_scale_m"))
        jacobian_scale = float(identity.get("jacobian_coordinate_scale_m"))
        minimum_jacobian = float(identity.get("minimum_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    order_generation = str(identity.get("node_order_generation") or "")
    scale_generation = str(identity.get("coordinate_scale_generation") or "")
    digest = str(identity.get("jacobian_table_sha256") or "")
    return (
        bool(mesh_generation)
        and identity.get("curving_mesh_generation") == mesh_generation
        and identity.get("jacobian_mesh_generation") == mesh_generation
        and bool(order_generation)
        and identity.get("jacobian_node_order_generation") == order_generation
        and bool(scale_generation)
        and identity.get("jacobian_coordinate_scale_generation") == scale_generation
        and order >= 2
        and jacobian_order == order
        and len(node_ids) == 8
        and len(set(node_ids)) == 8
        and jacobian_node_ids == node_ids
        and len(reference) == 8
        and all(len(row) == 3 for row in reference)
        and all(value in {-1.0, 1.0} for row in reference for value in row)
        and len({tuple(row) for row in reference}) == 8
        and jacobian_reference == reference
        and math.isfinite(scale)
        and scale > 0.0
        and math.isclose(jacobian_scale, scale, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isfinite(minimum_jacobian)
        and minimum_jacobian > 0.0
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("evaluated_jacobian_table_sha256") == digest
    )


def _tet_hex_pyramid_interface_orientation_conformity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("quad_face_node_ids"),
        identity.get("interface_quad_face_node_ids"),
        identity.get("pyramid_base_node_ids"),
        identity.get("interface_pyramid_base_node_ids"),
        identity.get("pyramid_apex_node_ids"),
        identity.get("interface_pyramid_apex_node_ids"),
        identity.get("face_orientation_signs"),
        identity.get("interface_face_orientation_signs"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        quad_faces = [[int(value) for value in face] for face in fields[0]]
        interface_quads = [[int(value) for value in face] for face in fields[1]]
        pyramid_bases = [[int(value) for value in face] for face in fields[2]]
        interface_bases = [[int(value) for value in face] for face in fields[3]]
        apex_ids = [int(value) for value in fields[4]]
        interface_apex_ids = [int(value) for value in fields[5]]
        signs = [int(value) for value in fields[6]]
        interface_signs = [int(value) for value in fields[7]]
    except (TypeError, ValueError):
        return False
    tet_generation = str(identity.get("tet_mesh_generation") or "")
    hex_generation = str(identity.get("hex_mesh_generation") or "")
    pyramid_generation = str(identity.get("pyramid_transition_generation") or "")
    digest = str(identity.get("interface_conformity_sha256") or "")
    return (
        bool(str(identity.get("interface_generation") or ""))
        and bool(tet_generation)
        and identity.get("interface_tet_mesh_generation") == tet_generation
        and bool(hex_generation)
        and identity.get("interface_hex_mesh_generation") == hex_generation
        and bool(pyramid_generation)
        and identity.get("interface_pyramid_transition_generation") == pyramid_generation
        and bool(quad_faces)
        and all(len(face) == 4 and len(set(face)) == 4 for face in quad_faces)
        and interface_quads == quad_faces
        and pyramid_bases == quad_faces
        and interface_bases == pyramid_bases
        and len(apex_ids) == len(quad_faces)
        and all(
            apex not in face for apex, face in zip(apex_ids, pyramid_bases)
        )
        and interface_apex_ids == apex_ids
        and len(signs) == len(quad_faces)
        and all(sign in {-1, 1} for sign in signs)
        and interface_signs == signs
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("exported_interface_conformity_sha256") == digest
    )


def _journal_transaction_undo_entity_reuse_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("created_entity_ids"),
        identity.get("replay_created_entity_ids"),
        identity.get("group_entity_ids"),
        identity.get("replay_group_entity_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        created = [int(value) for value in fields[0]]
        replay_created = [int(value) for value in fields[1]]
        grouped = [int(value) for value in fields[2]]
        replay_grouped = [int(value) for value in fields[3]]
        reset_epoch = int(identity.get("reset_epoch"))
        replay_epoch = int(identity.get("replay_reset_epoch"))
        undo_depth = int(identity.get("undo_depth"))
        replay_undo = int(identity.get("replay_undo_depth"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("journal_generation") or "")
    transaction = str(identity.get("transaction_id") or "")
    digest = str(identity.get("transaction_entity_table_sha256") or "")
    return (
        bool(generation)
        and identity.get("transaction_journal_generation") == generation
        and identity.get("entity_table_journal_generation") == generation
        and identity.get("group_table_journal_generation") == generation
        and bool(transaction)
        and identity.get("replay_transaction_id") == transaction
        and reset_epoch >= 0
        and replay_epoch == reset_epoch
        and undo_depth >= 0
        and replay_undo == undo_depth
        and bool(created)
        and len(set(created)) == len(created)
        and replay_created == created
        and set(grouped).issubset(created)
        and replay_grouped == grouped
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("replay_transaction_entity_table_sha256") == digest
    )


def _netgen_vol_block_order_curving_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("element_block_ids"),
        identity.get("element_block_types"),
        identity.get("writer_element_block_types"),
        identity.get("element_orders"),
        identity.get("writer_element_orders"),
        identity.get("curving_node_counts"),
        identity.get("writer_curving_node_counts"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        block_ids = [int(value) for value in fields[0]]
        block_types = [str(value) for value in fields[1]]
        writer_types = [str(value) for value in fields[2]]
        orders = [int(value) for value in fields[3]]
        writer_orders = [int(value) for value in fields[4]]
        node_counts = [int(value) for value in fields[5]]
        writer_counts = [int(value) for value in fields[6]]
    except (TypeError, ValueError):
        return False
    mesh_generation = str(identity.get("mesh_generation") or "")
    export_generation = str(identity.get("export_generation") or "")
    curving_generation = str(identity.get("curving_generation") or "")
    digest = str(identity.get("element_block_table_sha256") or "")
    allowed_counts = {
        ("tet", 1): 4,
        ("tet", 2): 10,
        ("hex", 1): 8,
        ("hex", 2): 20,
        ("pyramid", 1): 5,
    }
    return (
        bool(mesh_generation)
        and identity.get("writer_mesh_generation") == mesh_generation
        and bool(export_generation)
        and identity.get("writer_export_generation") == export_generation
        and bool(curving_generation)
        and identity.get("writer_curving_generation") == curving_generation
        and bool(block_ids)
        and len(set(block_ids)) == len(block_ids)
        and len(block_ids) == len(block_types) == len(orders) == len(node_counts)
        and writer_types == block_types
        and writer_orders == orders
        and writer_counts == node_counts
        and all(
            allowed_counts.get((family, order)) == count
            for family, order, count in zip(block_types, orders, node_counts)
        )
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("writer_element_block_table_sha256") == digest
    )


def _hex_sweep_face_vertex_twist_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("source_vertex_ids"),
        identity.get("mapped_source_vertex_ids"),
        identity.get("target_vertex_ids"),
        identity.get("mapped_target_vertex_ids"),
        identity.get("twist_path_vertex_ids"),
        identity.get("mapped_twist_path_vertex_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        source_vertices = [int(value) for value in fields[0]]
        mapped_source_vertices = [int(value) for value in fields[1]]
        target_vertices = [int(value) for value in fields[2]]
        mapped_target_vertices = [int(value) for value in fields[3]]
        twist_path = [int(value) for value in fields[4]]
        mapped_twist_path = [int(value) for value in fields[5]]
        source_face = int(identity.get("source_face_id"))
        mapped_source_face = int(identity.get("mapped_source_face_id"))
        target_face = int(identity.get("target_face_id"))
        mapped_target_face = int(identity.get("mapped_target_face_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    face_digest = str(identity.get("face_vertex_map_sha256") or "")
    twist_digest = str(identity.get("twist_path_sha256") or "")
    return (
        bool(generation)
        and identity.get("source_face_sweep_generation") == generation
        and identity.get("target_face_sweep_generation") == generation
        and identity.get("vertex_map_sweep_generation") == generation
        and identity.get("twist_path_sweep_generation") == generation
        and source_face > 0
        and target_face > 0
        and source_face != target_face
        and mapped_source_face == source_face
        and mapped_target_face == target_face
        and len(source_vertices) == 4
        and len(set(source_vertices)) == 4
        and mapped_source_vertices == source_vertices
        and len(target_vertices) == 4
        and len(set(target_vertices)) == 4
        and mapped_target_vertices == target_vertices
        and not set(source_vertices).intersection(target_vertices)
        and len(twist_path) == 8
        and set(twist_path) == set(source_vertices + target_vertices)
        and mapped_twist_path == twist_path
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in (face_digest, twist_digest)
        )
        and identity.get("applied_face_vertex_map_sha256") == face_digest
        and identity.get("applied_twist_path_sha256") == twist_digest
    )


def _quality_histogram_metric_element_unit_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("element_ids"),
        identity.get("evaluated_element_ids"),
        identity.get("metric_values"),
        identity.get("evaluated_metric_values"),
        identity.get("histogram_bin_edges"),
        identity.get("histogram_counts"),
        identity.get("evaluated_histogram_counts"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        element_ids = [int(value) for value in fields[0]]
        evaluated_ids = [int(value) for value in fields[1]]
        metric_values = [float(value) for value in fields[2]]
        evaluated_values = [float(value) for value in fields[3]]
        bin_edges = [float(value) for value in fields[4]]
        counts = [int(value) for value in fields[5]]
        evaluated_counts = [int(value) for value in fields[6]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    digest = str(identity.get("quality_table_sha256") or "")
    computed_counts = [0] * max(0, len(bin_edges) - 1)
    for value in metric_values:
        for index, (lower, upper) in enumerate(zip(bin_edges, bin_edges[1:])):
            if lower <= value < upper or (
                index == len(computed_counts) - 1 and value == upper
            ):
                computed_counts[index] += 1
                break
    return (
        bool(generation)
        and identity.get("metric_mesh_generation") == generation
        and identity.get("element_set_mesh_generation") == generation
        and identity.get("coordinate_unit_mesh_generation") == generation
        and bool(str(identity.get("metric_name") or ""))
        and identity.get("evaluated_metric_name") == identity.get("metric_name")
        and bool(str(identity.get("coordinate_unit") or ""))
        and identity.get("evaluated_coordinate_unit")
        == identity.get("coordinate_unit")
        and bool(element_ids)
        and len(set(element_ids)) == len(element_ids)
        and evaluated_ids == element_ids
        and len(metric_values) == len(element_ids)
        and all(math.isfinite(value) for value in metric_values)
        and evaluated_values == metric_values
        and len(bin_edges) >= 2
        and all(math.isfinite(value) for value in bin_edges)
        and all(right > left for left, right in zip(bin_edges, bin_edges[1:]))
        and len(counts) == len(bin_edges) - 1
        and all(value >= 0 for value in counts)
        and counts == computed_counts
        and evaluated_counts == counts
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("evaluated_quality_table_sha256") == digest
    )


def _block_sideset_merge_renumber_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("block_ids"),
        identity.get("exported_block_ids"),
        identity.get("block_entity_ids"),
        identity.get("exported_block_entity_ids"),
        identity.get("sideset_ids"),
        identity.get("exported_sideset_ids"),
        identity.get("sideset_entity_ids"),
        identity.get("exported_sideset_entity_ids"),
        identity.get("group_entity_ids"),
        identity.get("exported_group_entity_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        block_ids = [int(value) for value in fields[0]]
        exported_block_ids = [int(value) for value in fields[1]]
        block_entities = [[int(value) for value in row] for row in fields[2]]
        exported_block_entities = [
            [int(value) for value in row] for row in fields[3]
        ]
        sideset_ids = [int(value) for value in fields[4]]
        exported_sideset_ids = [int(value) for value in fields[5]]
        sideset_entities = [[int(value) for value in row] for row in fields[6]]
        exported_sideset_entities = [
            [int(value) for value in row] for row in fields[7]
        ]
        group_entities = [int(value) for value in fields[8]]
        exported_group_entities = [int(value) for value in fields[9]]
    except (TypeError, ValueError):
        return False
    topology_generation = str(identity.get("topology_generation") or "")
    merge_generation = str(identity.get("merge_transaction_generation") or "")
    digest = str(identity.get("ownership_table_sha256") or "")
    owned_entities = [
        entity
        for membership in block_entities + sideset_entities
        for entity in membership
    ]
    return (
        bool(topology_generation)
        and identity.get("block_topology_generation") == topology_generation
        and identity.get("sideset_topology_generation") == topology_generation
        and identity.get("group_topology_generation") == topology_generation
        and identity.get("renumber_topology_generation") == topology_generation
        and bool(merge_generation)
        and identity.get("group_merge_transaction_generation")
        == merge_generation
        and bool(block_ids)
        and len(set(block_ids)) == len(block_ids)
        and len(block_entities) == len(block_ids)
        and all(membership for membership in block_entities)
        and exported_block_ids == block_ids
        and exported_block_entities == block_entities
        and bool(sideset_ids)
        and len(set(sideset_ids)) == len(sideset_ids)
        and len(sideset_entities) == len(sideset_ids)
        and all(membership for membership in sideset_entities)
        and exported_sideset_ids == sideset_ids
        and exported_sideset_entities == sideset_entities
        and len(set(group_entities)) == len(group_entities)
        and set(group_entities) == set(owned_entities)
        and exported_group_entities == group_entities
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("exported_ownership_table_sha256") == digest
    )


def _aprepro_include_variable_transaction_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("variable_names"),
        identity.get("expanded_variable_names"),
        identity.get("variable_values"),
        identity.get("expanded_variable_values"),
        identity.get("include_paths"),
        identity.get("expanded_include_paths"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        names = [str(value) for value in fields[0]]
        expanded_names = [str(value) for value in fields[1]]
        values = [float(value) for value in fields[2]]
        expanded_values = [float(value) for value in fields[3]]
        includes = [str(value) for value in fields[4]]
        expanded_includes = [str(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("journal_transaction_generation") or "")
    working_directory = str(identity.get("working_directory") or "")
    variable_digest = str(identity.get("variable_table_sha256") or "")
    include_digest = str(identity.get("include_tree_sha256") or "")
    return (
        bool(generation)
        and identity.get("variable_table_transaction_generation") == generation
        and identity.get("include_expansion_transaction_generation") == generation
        and identity.get("working_directory_transaction_generation") == generation
        and bool(working_directory)
        and not Path(working_directory).is_absolute()
        and ".." not in Path(working_directory).parts
        and identity.get("replay_working_directory") == working_directory
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and expanded_names == names
        and len(values) == len(names)
        and all(math.isfinite(value) for value in values)
        and expanded_values == values
        and bool(includes)
        and all(
            path
            and not Path(path).is_absolute()
            and ".." not in Path(path).parts
            for path in includes
        )
        and expanded_includes == includes
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in (variable_digest, include_digest)
        )
        and identity.get("expanded_variable_table_sha256") == variable_digest
        and identity.get("expanded_include_tree_sha256") == include_digest
    )


def _step_ap214_multibody_export_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("body_names"),
        identity.get("exported_body_names"),
        identity.get("body_placement_ids"),
        identity.get("exported_body_placement_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        names = [str(value) for value in fields[0]]
        exported_names = [str(value) for value in fields[1]]
        placements = [int(value) for value in fields[2]]
        exported_placements = [int(value) for value in fields[3]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("export_generation") or "")
    digest_pairs = (
        ("geometry_sha256", "exported_geometry_sha256"),
        ("placement_table_sha256", "exported_placement_table_sha256"),
        ("mass_property_table_sha256", "exported_mass_property_table_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "geometry_export_generation",
                "body_name_export_generation",
                "unit_export_generation",
                "placement_export_generation",
                "mass_property_export_generation",
            )
        )
        and identity.get("step_schema") == "AP214"
        and identity.get("exported_step_schema") == "AP214"
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and exported_names == names
        and bool(str(identity.get("length_unit") or ""))
        and identity.get("exported_length_unit") == identity.get("length_unit")
        and len(placements) == len(names)
        and all(value > 0 for value in placements)
        and len(set(placements)) == len(placements)
        and exported_placements == placements
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(exported) == identity.get(source)
            for source, exported in digest_pairs
        )
    )


def _hybrid_transition_topology_block_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("pyramid_face_orientations"),
        identity.get("exported_pyramid_face_orientations"),
        identity.get("shared_node_ids"),
        identity.get("exported_shared_node_ids"),
        identity.get("material_block_ids"),
        identity.get("exported_material_block_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        orientations = [int(value) for value in fields[0]]
        exported_orientations = [int(value) for value in fields[1]]
        nodes = [int(value) for value in fields[2]]
        exported_nodes = [int(value) for value in fields[3]]
        blocks = [int(value) for value in fields[4]]
        exported_blocks = [int(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    topology_digest = str(identity.get("transition_topology_sha256") or "")
    block_digest = str(identity.get("material_block_map_sha256") or "")
    return (
        bool(generation)
        and identity.get("pyramid_orientation_mesh_generation") == generation
        and identity.get("shared_node_topology_mesh_generation") == generation
        and identity.get("material_block_mesh_generation") == generation
        and bool(orientations)
        and all(value in {-1, 1} for value in orientations)
        and exported_orientations == orientations
        and bool(nodes)
        and all(value > 0 for value in nodes)
        and len(set(nodes)) == len(nodes)
        and exported_nodes == nodes
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and exported_blocks == blocks
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in (topology_digest, block_digest)
        )
        and identity.get("exported_transition_topology_sha256") == topology_digest
        and identity.get("exported_material_block_map_sha256") == block_digest
    )


def _headless_step_export_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("selected_body_ids"),
        identity.get("exported_body_ids"),
        identity.get("body_names"),
        identity.get("exported_body_names"),
        identity.get("transform_ids"),
        identity.get("exported_transform_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        selected = [int(value) for value in fields[0]]
        exported = [int(value) for value in fields[1]]
        names = [str(value) for value in fields[2]]
        exported_names = [str(value) for value in fields[3]]
        transforms = [int(value) for value in fields[4]]
        exported_transforms = [int(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("model_generation") or "")
    digest_pairs = (
        ("step_sha256", "exported_step_sha256"),
        ("transform_table_sha256", "exported_transform_table_sha256"),
        ("export_log_sha256", "recorded_export_log_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "selected_body_model_generation",
                "transform_model_generation",
                "name_model_generation",
                "export_log_model_generation",
            )
        )
        and bool(selected)
        and all(value > 0 for value in selected)
        and len(set(selected)) == len(selected)
        and exported == selected
        and len(names) == len(selected)
        and all(names)
        and len(set(names)) == len(names)
        and exported_names == names
        and len(transforms) == len(selected)
        and all(value > 0 for value in transforms)
        and len(set(transforms)) == len(transforms)
        and exported_transforms == transforms
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(recorded) == identity.get(source)
            for source, recorded in digest_pairs
        )
    )


def _geometry_heal_imprint_merge_ownership_generation_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("source_entity_ids"),
        identity.get("owned_entity_ids"),
        identity.get("imprint_pair_ids"),
        identity.get("owned_imprint_pair_ids"),
        identity.get("merge_survivor_ids"),
        identity.get("owned_merge_survivor_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        source_entities = [int(value) for value in fields[0]]
        owned_entities = [int(value) for value in fields[1]]
        imprint_pairs = [int(value) for value in fields[2]]
        owned_imprint_pairs = [int(value) for value in fields[3]]
        merge_survivors = [int(value) for value in fields[4]]
        owned_merge_survivors = [int(value) for value in fields[5]]
        tolerance = float(identity.get("heal_tolerance"))
        owned_tolerance = float(identity.get("ownership_heal_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("geometry_generation") or "")
    ownership_digest = str(identity.get("ownership_map_sha256") or "")
    operation_digest = str(identity.get("operation_log_sha256") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "heal_geometry_generation",
                "imprint_geometry_generation",
                "merge_geometry_generation",
                "ownership_geometry_generation",
            )
        )
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and math.isclose(owned_tolerance, tolerance, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and bool(source_entities)
        and all(value > 0 for value in source_entities)
        and len(set(source_entities)) == len(source_entities)
        and owned_entities == source_entities
        and bool(imprint_pairs)
        and all(value > 0 for value in imprint_pairs)
        and len(set(imprint_pairs)) == len(imprint_pairs)
        and owned_imprint_pairs == imprint_pairs
        and bool(merge_survivors)
        and all(value > 0 for value in merge_survivors)
        and len(set(merge_survivors)) == len(merge_survivors)
        and owned_merge_survivors == merge_survivors
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in (ownership_digest, operation_digest)
        )
        and identity.get("recorded_ownership_map_sha256") == ownership_digest
        and identity.get("recorded_operation_log_sha256") == operation_digest
    )


def _parallel_sculpt_partition_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("partition_ids"),
        identity.get("assembled_partition_ids"),
        identity.get("ghost_interface_ids"),
        identity.get("assembled_ghost_interface_ids"),
        identity.get("refinement_levels"),
        identity.get("assembled_refinement_levels"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        partitions = [int(value) for value in fields[0]]
        assembled_partitions = [int(value) for value in fields[1]]
        ghost_interfaces = [int(value) for value in fields[2]]
        assembled_ghost_interfaces = [int(value) for value in fields[3]]
        refinement_levels = [int(value) for value in fields[4]]
        assembled_refinement_levels = [int(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    digest_pairs = (
        ("partition_map_sha256", "assembled_partition_map_sha256"),
        ("ghost_interface_sha256", "assembled_ghost_interface_sha256"),
        ("sculpt_mesh_sha256", "assembled_sculpt_mesh_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "partition_mesh_generation",
                "ghost_interface_mesh_generation",
                "refinement_mesh_generation",
                "assembly_mesh_generation",
            )
        )
        and bool(partitions)
        and all(value > 0 for value in partitions)
        and len(set(partitions)) == len(partitions)
        and assembled_partitions == partitions
        and bool(ghost_interfaces)
        and all(value > 0 for value in ghost_interfaces)
        and len(set(ghost_interfaces)) == len(ghost_interfaces)
        and assembled_ghost_interfaces == ghost_interfaces
        and len(refinement_levels) == len(partitions)
        and all(value >= 0 for value in refinement_levels)
        and assembled_refinement_levels == refinement_levels
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(assembled) == identity.get(source)
            for source, assembled in digest_pairs
        )
    )


def _mixed_transition_interface_ownership_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("element_families"),
        identity.get("exported_element_families"),
        identity.get("interface_face_ids"),
        identity.get("exported_interface_face_ids"),
        identity.get("block_ids"),
        identity.get("exported_block_ids"),
        identity.get("sideset_ids"),
        identity.get("exported_sideset_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        families = [str(value) for value in fields[0]]
        exported_families = [str(value) for value in fields[1]]
        interface_faces = [int(value) for value in fields[2]]
        exported_interface_faces = [int(value) for value in fields[3]]
        blocks = [int(value) for value in fields[4]]
        exported_blocks = [int(value) for value in fields[5]]
        sidesets = [int(value) for value in fields[6]]
        exported_sidesets = [int(value) for value in fields[7]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    digest_pairs = (
        ("transition_table_sha256", "exported_transition_table_sha256"),
        ("ownership_table_sha256", "exported_ownership_table_sha256"),
    )
    id_pairs = (
        (interface_faces, exported_interface_faces),
        (blocks, exported_blocks),
        (sidesets, exported_sidesets),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "interface_mesh_generation",
                "element_family_mesh_generation",
                "block_mesh_generation",
                "sideset_mesh_generation",
            )
        )
        and families == ["tet", "pyramid", "hex"]
        and exported_families == families
        and all(
            bool(source)
            and all(value > 0 for value in source)
            and len(set(source)) == len(source)
            and exported == source
            for source, exported in id_pairs
        )
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(exported) == identity.get(source)
            for source, exported in digest_pairs
        )
    )


def _journal_replay_entity_version_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    fields = (
        identity.get("geometry_entity_ids"),
        identity.get("replay_geometry_entity_ids"),
        identity.get("entity_names"),
        identity.get("replay_entity_names"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in fields
    ):
        return False
    try:
        entity_ids = [int(value) for value in fields[0]]
        replay_entity_ids = [int(value) for value in fields[1]]
        names = [str(value) for value in fields[2]]
        replay_names = [str(value) for value in fields[3]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("replay_generation") or "")
    version = str(identity.get("application_version") or "")
    digest_pairs = (
        ("entity_map_sha256", "replay_entity_map_sha256"),
        ("journal_sha256", "replay_journal_sha256"),
        ("command_log_sha256", "replay_command_log_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "geometry_replay_generation",
                "entity_map_replay_generation",
                "application_version_replay_generation",
                "command_log_replay_generation",
            )
        )
        and bool(version)
        and identity.get("replay_application_version") == version
        and bool(entity_ids)
        and all(value > 0 for value in entity_ids)
        and len(set(entity_ids)) == len(entity_ids)
        and replay_entity_ids == entity_ids
        and len(names) == len(entity_ids)
        and all(names)
        and len(set(names)) == len(names)
        and replay_names == names
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(replay) == identity.get(source)
            for source, replay in digest_pairs
        )
    )


def _exodus64_mapping_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    flat_fields = (
        identity.get("entity_ids"),
        identity.get("decoded_entity_ids"),
        identity.get("sideset_ids"),
        identity.get("decoded_sideset_ids"),
    )
    nested_fields = (
        identity.get("sideset_entity_ids"),
        identity.get("decoded_sideset_entity_ids"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in (*flat_fields, *nested_fields)
    ):
        return False
    if not all(
        isinstance(row, Sequence) and not isinstance(row, (str, bytes))
        for values in nested_fields
        for row in values
    ):
        return False
    try:
        entity_ids = [int(value) for value in flat_fields[0]]
        decoded_entity_ids = [int(value) for value in flat_fields[1]]
        sideset_ids = [int(value) for value in flat_fields[2]]
        decoded_sideset_ids = [int(value) for value in flat_fields[3]]
        sideset_entities = [
            [int(value) for value in row] for row in nested_fields[0]
        ]
        decoded_sideset_entities = [
            [int(value) for value in row] for row in nested_fields[1]
        ]
        integer_width = int(identity.get("integer_width_bits"))
        decoded_integer_width = int(identity.get("decoded_integer_width_bits"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("export_generation") or "")
    digest_pairs = (
        ("element_map_sha256", "decoded_element_map_sha256"),
        ("schema_sha256", "decoded_schema_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "entity_id_export_generation",
                "sideset_export_generation",
                "element_map_export_generation",
                "schema_export_generation",
            )
        )
        and integer_width == 64
        and decoded_integer_width == integer_width
        and bool(entity_ids)
        and all(value > 2**32 for value in entity_ids)
        and len(set(entity_ids)) == len(entity_ids)
        and decoded_entity_ids == entity_ids
        and bool(sideset_ids)
        and all(value > 0 for value in sideset_ids)
        and len(set(sideset_ids)) == len(sideset_ids)
        and decoded_sideset_ids == sideset_ids
        and len(sideset_entities) == len(sideset_ids)
        and all(row for row in sideset_entities)
        and all(value in entity_ids for row in sideset_entities for value in row)
        and decoded_sideset_entities == sideset_entities
        and all(
            len(str(identity.get(source) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(source) or "")
            )
            and identity.get(decoded) == identity.get(source)
            for source, decoded in digest_pairs
        )
    )


def _high_order_hex_jacobian_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    sequence_keys = (
        "element_ids",
        "result_element_ids",
        "minimum_jacobian",
        "result_minimum_jacobian",
        "coordinate_transform",
        "result_coordinate_transform",
        "canonical_node_order",
        "exported_node_order",
        "block_ids",
        "result_block_ids",
    )
    if not all(
        isinstance(identity.get(key), Sequence)
        and not isinstance(identity.get(key), (str, bytes))
        for key in sequence_keys
    ):
        return False
    try:
        element_ids = [int(value) for value in identity["element_ids"]]
        result_element_ids = [int(value) for value in identity["result_element_ids"]]
        jacobians = [float(value) for value in identity["minimum_jacobian"]]
        result_jacobians = [
            float(value) for value in identity["result_minimum_jacobian"]
        ]
        transform = [
            [float(value) for value in row]
            for row in identity["coordinate_transform"]
        ]
        result_transform = [
            [float(value) for value in row]
            for row in identity["result_coordinate_transform"]
        ]
        node_order = [int(value) for value in identity["canonical_node_order"]]
        exported_node_order = [int(value) for value in identity["exported_node_order"]]
        block_ids = [int(value) for value in identity["block_ids"]]
        result_block_ids = [int(value) for value in identity["result_block_ids"]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "jacobian_mesh_generation",
                "transform_mesh_generation",
                "node_order_mesh_generation",
                "block_mesh_generation",
                "result_mesh_generation",
            )
        )
        and bool(element_ids)
        and len(set(element_ids)) == len(element_ids)
        and all(value > 0 for value in element_ids)
        and result_element_ids == element_ids
        and len(jacobians) == len(element_ids)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and len(transform) == 3
        and all(len(row) == 3 for row in transform)
        and all(math.isfinite(value) for row in transform for value in row)
        and result_transform == transform
        and len(node_order) >= 8
        and len(set(node_order)) == len(node_order)
        and exported_node_order == node_order
        and bool(block_ids)
        and len(set(block_ids)) == len(block_ids)
        and result_block_ids == block_ids
        and len(str(identity.get("mesh_sha256") or "")) == 64
        and all(
            character in "0123456789abcdef"
            for character in str(identity.get("mesh_sha256") or "")
        )
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
    )


def _mesh_cad_closure_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        cad_volume = float(identity.get("cad_volume"))
        mesh_volume = float(identity.get("mesh_volume"))
        cad_area = float(identity.get("cad_boundary_area"))
        mesh_area = float(identity.get("mesh_boundary_area"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("geometry_generation") or "")
    unit = str(identity.get("length_unit") or "")
    frame = str(identity.get("coordinate_frame") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "cad_geometry_generation",
                "mesh_geometry_generation",
                "boundary_geometry_generation",
                "unit_geometry_generation",
                "frame_geometry_generation",
                "result_geometry_generation",
            )
        )
        and bool(unit)
        and identity.get("result_length_unit") == unit
        and bool(frame)
        and identity.get("result_coordinate_frame") == frame
        and all(
            math.isfinite(value) and value > 0.0
            for value in (cad_volume, mesh_volume, cad_area, mesh_area)
        )
        and abs(mesh_volume - cad_volume) <= 1.0e-9 * cad_volume
        and abs(mesh_area - cad_area) <= 1.0e-9 * cad_area
        and identity.get("boundary_closed") is True
        and identity.get("result_boundary_closed") is True
        and len(str(identity.get("cad_shape_sha256") or "")) == 64
        and identity.get("mesh_source_shape_sha256") == identity.get("cad_shape_sha256")
        and len(str(identity.get("closure_table_sha256") or "")) == 64
        and identity.get("result_closure_table_sha256")
        == identity.get("closure_table_sha256")
    )


def _webcut_imprint_merge_topology_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        imprint_tolerance = float(identity.get("imprint_tolerance"))
        result_imprint_tolerance = float(identity.get("result_imprint_tolerance"))
        merge_tolerance = float(identity.get("merge_tolerance"))
        result_merge_tolerance = float(identity.get("result_merge_tolerance"))
        names = [str(value) for value in identity.get("entity_names", [])]
        result_names = [str(value) for value in identity.get("result_entity_names", [])]
    except (TypeError, ValueError):
        return False
    topology = identity.get("topology_counts")
    result_topology = identity.get("result_topology_counts")
    generation = str(identity.get("operation_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "webcut_operation_generation",
                "imprint_operation_generation",
                "merge_operation_generation",
                "topology_operation_generation",
                "result_operation_generation",
            )
        )
        and math.isfinite(imprint_tolerance)
        and imprint_tolerance > 0.0
        and result_imprint_tolerance == imprint_tolerance
        and math.isfinite(merge_tolerance)
        and merge_tolerance > 0.0
        and result_merge_tolerance == merge_tolerance
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and result_names == names
        and isinstance(topology, Mapping)
        and set(topology) == {"volume", "surface", "curve", "vertex"}
        and all(isinstance(value, int) and value > 0 for value in topology.values())
        and result_topology == topology
        and len(str(identity.get("topology_sha256") or "")) == 64
        and identity.get("result_topology_sha256") == identity.get("topology_sha256")
        and len(str(identity.get("command_log_sha256") or "")) == 64
        and identity.get("result_command_log_sha256")
        == identity.get("command_log_sha256")
    )


def _exodus_qa_coordinate_distribution_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    sequence_keys = (
        "qa_records",
        "decoded_qa_records",
        "coordinate_names",
        "decoded_coordinate_names",
        "sideset_distribution_factors",
        "decoded_sideset_distribution_factors",
    )
    if not all(
        isinstance(identity.get(key), Sequence)
        and not isinstance(identity.get(key), (str, bytes))
        for key in sequence_keys
    ):
        return False
    try:
        qa_records = [[str(value) for value in row] for row in identity["qa_records"]]
        decoded_qa = [
            [str(value) for value in row] for row in identity["decoded_qa_records"]
        ]
        coordinates = [str(value) for value in identity["coordinate_names"]]
        decoded_coordinates = [
            str(value) for value in identity["decoded_coordinate_names"]
        ]
        factors = [
            [float(value) for value in row]
            for row in identity["sideset_distribution_factors"]
        ]
        decoded_factors = [
            [float(value) for value in row]
            for row in identity["decoded_sideset_distribution_factors"]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("export_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "qa_export_generation",
                "coordinate_export_generation",
                "distribution_export_generation",
                "checksum_export_generation",
                "result_export_generation",
            )
        )
        and bool(qa_records)
        and all(len(row) == 4 and all(row) for row in qa_records)
        and decoded_qa == qa_records
        and coordinates == ["x", "y", "z"]
        and decoded_coordinates == coordinates
        and bool(factors)
        and all(row for row in factors)
        and all(math.isfinite(value) and value >= 0.0 for row in factors for value in row)
        and decoded_factors == factors
        and len(str(identity.get("payload_sha256") or "")) == 64
        and identity.get("decoded_payload_sha256") == identity.get("payload_sha256")
        and len(str(identity.get("qa_table_sha256") or "")) == 64
        and identity.get("decoded_qa_table_sha256")
        == identity.get("qa_table_sha256")
    )


def _hybrid_interface_conformity_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    counts = identity.get("element_counts")
    result_counts = identity.get("result_element_counts")
    try:
        faces = [
            [int(value) for value in row]
            for row in identity.get("interface_face_node_ids", [])
        ]
        result_faces = [
            [int(value) for value in row]
            for row in identity.get("result_interface_face_node_ids", [])
        ]
        orientations = [int(value) for value in identity.get("orientation_signs", [])]
        result_orientations = [
            int(value) for value in identity.get("result_orientation_signs", [])
        ]
        blocks = [int(value) for value in identity.get("block_ids", [])]
        result_blocks = [int(value) for value in identity.get("result_block_ids", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mesh_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "tet_mesh_generation",
                "hex_mesh_generation",
                "pyramid_mesh_generation",
                "interface_mesh_generation",
                "orientation_mesh_generation",
                "block_mesh_generation",
                "result_mesh_generation",
            )
        )
        and isinstance(counts, Mapping)
        and set(counts) == {"tet4", "hex8", "pyramid5"}
        and all(isinstance(value, int) and value > 0 for value in counts.values())
        and result_counts == counts
        and bool(faces)
        and all(len(row) == 4 and len(set(row)) == 4 and all(value > 0 for value in row) for row in faces)
        and result_faces == faces
        and identity.get("interface_conforming") is True
        and identity.get("result_interface_conforming") is True
        and orientations == [1, 1, 1]
        and result_orientations == orientations
        and len(blocks) == 3
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and result_blocks == blocks
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
    )


def _periodic_sideset_pairing_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        pairs = [[int(value) for value in row] for row in identity.get("node_pairs", [])]
        result_pairs = [
            [int(value) for value in row] for row in identity.get("result_node_pairs", [])
        ]
        transform = [
            [float(value) for value in row] for row in identity.get("rigid_transform", [])
        ]
        result_transform = [
            [float(value) for value in row]
            for row in identity.get("result_rigid_transform", [])
        ]
        tolerance = float(identity.get("pairing_tolerance"))
        result_tolerance = float(identity.get("result_pairing_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("periodic_generation") or "")
    master = str(identity.get("master_sideset") or "")
    slave = str(identity.get("slave_sideset") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "sideset_periodic_generation",
                "node_pair_periodic_generation",
                "transform_periodic_generation",
                "tolerance_periodic_generation",
                "geometry_periodic_generation",
                "result_periodic_generation",
            )
        )
        and bool(master)
        and bool(slave)
        and master != slave
        and identity.get("result_master_sideset") == master
        and identity.get("result_slave_sideset") == slave
        and bool(pairs)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 for row in pairs)
        and len({row[0] for row in pairs}) == len(pairs)
        and len({row[1] for row in pairs}) == len(pairs)
        and result_pairs == pairs
        and len(transform) == 3
        and all(len(row) == 4 for row in transform)
        and all(math.isfinite(value) for row in transform for value in row)
        and result_transform == transform
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and _valid_sha256(identity.get("geometry_sha256"))
        and identity.get("result_geometry_sha256") == identity.get("geometry_sha256")
        and _valid_sha256(identity.get("periodic_table_sha256"))
        and identity.get("result_periodic_table_sha256")
        == identity.get("periodic_table_sha256")
    )


def _journal_reset_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    commands = identity.get("command_sequence")
    result_commands = identity.get("result_command_sequence")
    entity_map = identity.get("entity_id_map")
    result_entity_map = identity.get("result_entity_id_map")
    if not all(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes))
        for value in (commands, result_commands, entity_map, result_entity_map)
    ):
        return False
    try:
        command_rows = [str(value) for value in commands]
        result_command_rows = [str(value) for value in result_commands]
        entities = [[str(row[0]), int(row[1])] for row in entity_map if len(row) == 2]
        result_entities = [
            [str(row[0]), int(row[1])] for row in result_entity_map if len(row) == 2
        ]
    except (TypeError, ValueError, IndexError):
        return False
    generation = str(identity.get("session_generation") or "")
    session = str(identity.get("session_id") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "reset_session_generation",
                "entity_session_generation",
                "undo_session_generation",
                "replay_session_generation",
                "topology_session_generation",
                "result_session_generation",
            )
        )
        and bool(session)
        and identity.get("result_session_id") == session
        and identity.get("reset_applied") is True
        and identity.get("result_reset_applied") is True
        and bool(command_rows)
        and command_rows[0].strip().lower() == "reset"
        and result_command_rows == command_rows
        and bool(entities)
        and all(row[0] and row[1] > 0 for row in entities)
        and len({row[0] for row in entities}) == len(entities)
        and result_entities == entities
        and _valid_sha256(identity.get("undo_checkpoint_sha256"))
        and identity.get("result_undo_checkpoint_sha256")
        == identity.get("undo_checkpoint_sha256")
        and _valid_sha256(identity.get("topology_sha256"))
        and identity.get("result_topology_sha256") == identity.get("topology_sha256")
        and _valid_sha256(identity.get("journal_sha256"))
        and identity.get("result_journal_sha256") == identity.get("journal_sha256")
    )


def _netgen_vol_export_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        order = int(identity.get("polynomial_order"))
        decoded_order = int(identity.get("decoded_polynomial_order"))
        boundaries = [str(value) for value in identity.get("boundary_names", [])]
        decoded_boundaries = [
            str(value) for value in identity.get("decoded_boundary_names", [])
        ]
        materials = [int(value) for value in identity.get("material_indices", [])]
        decoded_materials = [
            int(value) for value in identity.get("decoded_material_indices", [])
        ]
        entities = [int(value) for value in identity.get("source_entity_ids", [])]
        decoded_entities = [
            int(value) for value in identity.get("decoded_source_entity_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("export_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "family_export_generation",
                "order_export_generation",
                "boundary_export_generation",
                "material_export_generation",
                "entity_export_generation",
                "result_export_generation",
            )
        )
        and identity.get("volume_element_family") == "tet4"
        and identity.get("decoded_volume_element_family") == "tet4"
        and identity.get("surface_element_family") == "tri3"
        and identity.get("decoded_surface_element_family") == "tri3"
        and order == 1
        and decoded_order == order
        and bool(boundaries)
        and all(boundaries)
        and len(set(boundaries)) == len(boundaries)
        and decoded_boundaries == boundaries
        and bool(materials)
        and all(value > 0 for value in materials)
        and decoded_materials == materials
        and len(entities) == len(materials)
        and all(value > 0 for value in entities)
        and decoded_entities == entities
        and _valid_sha256(identity.get("source_mesh_sha256"))
        and identity.get("export_source_mesh_sha256") == identity.get("source_mesh_sha256")
        and _valid_sha256(identity.get("vol_sha256"))
        and identity.get("decoded_vol_sha256") == identity.get("vol_sha256")
    )


def _hex_sweep_layer_correspondence_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source = int(identity.get("source_face_id"))
        result_source = int(identity.get("result_source_face_id"))
        target = int(identity.get("target_face_id"))
        result_target = int(identity.get("result_target_face_id"))
        layers = int(identity.get("layer_count"))
        result_layers = int(identity.get("result_layer_count"))
        pairs = [[int(value) for value in row] for row in identity.get("source_target_vertex_pairs", [])]
        result_pairs = [[int(value) for value in row] for row in identity.get("result_source_target_vertex_pairs", [])]
        counts = [int(value) for value in identity.get("layer_element_counts", [])]
        result_counts = [int(value) for value in identity.get("result_layer_element_counts", [])]
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [float(value) for value in identity.get("result_scaled_jacobians", [])]
        block = int(identity.get("block_id"))
        result_block = int(identity.get("result_block_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_sweep_generation",
                "target_sweep_generation",
                "layer_sweep_generation",
                "correspondence_sweep_generation",
                "jacobian_sweep_generation",
                "block_sweep_generation",
                "result_sweep_generation",
            )
        )
        and source > 0
        and target > 0
        and source != target
        and result_source == source
        and result_target == target
        and layers > 0
        and result_layers == layers
        and len(pairs) >= 3
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 for row in pairs)
        and len({row[0] for row in pairs}) == len(pairs)
        and len({row[1] for row in pairs}) == len(pairs)
        and result_pairs == pairs
        and len(counts) == layers
        and all(value > 0 for value in counts)
        and result_counts == counts
        and len(jacobians) == layers
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and block > 0
        and result_block == block
        and _valid_sha256(identity.get("sweep_mesh_sha256"))
        and identity.get("result_sweep_mesh_sha256") == identity.get("sweep_mesh_sha256")
    )


def _high_order_hex_export_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        edges = [int(value) for value in identity.get("edge_node_order", [])]
        decoded_edges = [int(value) for value in identity.get("decoded_edge_node_order", [])]
        faces = [int(value) for value in identity.get("face_node_order", [])]
        decoded_faces = [int(value) for value in identity.get("decoded_face_node_order", [])]
        interior = [int(value) for value in identity.get("interior_node_order", [])]
        decoded_interior = [int(value) for value in identity.get("decoded_interior_node_order", [])]
        order = int(identity.get("curved_geometry_order"))
        decoded_order = int(identity.get("decoded_curved_geometry_order"))
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        decoded_jacobian = float(identity.get("decoded_minimum_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("export_generation") or "")
    all_nodes = edges + faces + interior
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "edge_export_generation",
                "face_export_generation",
                "interior_export_generation",
                "curvature_export_generation",
                "jacobian_export_generation",
                "result_export_generation",
            )
        )
        and identity.get("element_family") == "hex27"
        and identity.get("decoded_element_family") == "hex27"
        and len(edges) == 12
        and len(faces) == 6
        and len(interior) == 1
        and all(value > 0 for value in all_nodes)
        and len(set(all_nodes)) == len(all_nodes)
        and decoded_edges == edges
        and decoded_faces == faces
        and decoded_interior == interior
        and order == 2
        and decoded_order == order
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and decoded_jacobian == jacobian
        and _valid_sha256(identity.get("export_sha256"))
        and identity.get("decoded_export_sha256") == identity.get("export_sha256")
    )


def _sculpt_input_output_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        spacing = [float(value) for value in identity.get("voxel_spacing_m", [])]
        result_spacing = [float(value) for value in identity.get("result_voxel_spacing_m", [])]
        thresholds = [float(value) for value in identity.get("thresholds", [])]
        result_thresholds = [float(value) for value in identity.get("result_thresholds", [])]
        material_blocks = [[int(value) for value in row] for row in identity.get("material_to_block", [])]
        result_material_blocks = [[int(value) for value in row] for row in identity.get("result_material_to_block", [])]
        block_counts = [[int(value) for value in row] for row in identity.get("element_counts_by_block", [])]
        result_block_counts = [[int(value) for value in row] for row in identity.get("result_element_counts_by_block", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sculpt_generation") or "")
    session = str(identity.get("session_id") or "")
    block_ids = {row[1] for row in material_blocks if len(row) == 2}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "voxel_sculpt_generation",
                "threshold_sculpt_generation",
                "material_sculpt_generation",
                "block_sculpt_generation",
                "output_sculpt_generation",
                "result_sculpt_generation",
            )
        )
        and bool(session)
        and identity.get("result_session_id") == session
        and len(spacing) == 3
        and all(math.isfinite(value) and value > 0.0 for value in spacing)
        and result_spacing == spacing
        and bool(thresholds)
        and all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in thresholds)
        and all(right > left for left, right in zip(thresholds, thresholds[1:]))
        and result_thresholds == thresholds
        and bool(material_blocks)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 for row in material_blocks)
        and len({row[0] for row in material_blocks}) == len(material_blocks)
        and len(block_ids) == len(material_blocks)
        and result_material_blocks == material_blocks
        and {row[0] for row in block_counts if len(row) == 2} == block_ids
        and all(len(row) == 2 and row[1] > 0 for row in block_counts)
        and result_block_counts == block_counts
        and _valid_sha256(identity.get("input_volume_sha256"))
        and identity.get("result_input_volume_sha256") == identity.get("input_volume_sha256")
        and _valid_sha256(identity.get("output_mesh_sha256"))
        and identity.get("result_output_mesh_sha256") == identity.get("output_mesh_sha256")
    )


def _exodus_merge_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        tolerance = float(identity.get("merge_tolerance_m"))
        decoded_tolerance = float(identity.get("decoded_merge_tolerance_m"))
        pairs = [[int(value) for value in row] for row in identity.get("merged_node_pairs", [])]
        decoded_pairs = [[int(value) for value in row] for row in identity.get("decoded_merged_node_pairs", [])]
        global_ids = [int(value) for value in identity.get("global_node_ids", [])]
        decoded_global_ids = [int(value) for value in identity.get("decoded_global_node_ids", [])]
        blocks = [int(value) for value in identity.get("block_ids", [])]
        decoded_blocks = [int(value) for value in identity.get("decoded_block_ids", [])]
        sidesets = [int(value) for value in identity.get("sideset_ids", [])]
        decoded_sidesets = [int(value) for value in identity.get("decoded_sideset_ids", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("merge_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "tolerance_merge_generation",
                "node_merge_generation",
                "global_id_merge_generation",
                "block_merge_generation",
                "sideset_merge_generation",
                "result_merge_generation",
            )
        )
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and decoded_tolerance == tolerance
        and bool(pairs)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 for row in pairs)
        and len({row[0] for row in pairs}) == len(pairs)
        and len({row[1] for row in pairs}) == len(pairs)
        and decoded_pairs == pairs
        and bool(global_ids)
        and all(value > 0 for value in global_ids)
        and len(set(global_ids)) == len(global_ids)
        and decoded_global_ids == global_ids
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and decoded_blocks == blocks
        and bool(sidesets)
        and all(value > 0 for value in sidesets)
        and len(set(sidesets)) == len(sidesets)
        and decoded_sidesets == sidesets
        and _valid_sha256(identity.get("source_mesh_sha256"))
        and identity.get("merge_source_mesh_sha256") == identity.get("source_mesh_sha256")
        and _valid_sha256(identity.get("exodus_sha256"))
        and identity.get("decoded_exodus_sha256") == identity.get("exodus_sha256")
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
        "high_order_hex_jacobians_use_current_node_order_and_coordinate_scale": (
            _high_order_hex_jacobian_node_order_scale_generation_ok(
                summary.get("high_order_hex_jacobian_node_order_coordinate_scale_identity")
            )
        ),
        "tet_hex_pyramid_interface_uses_current_face_orientation_and_conformity": (
            _tet_hex_pyramid_interface_orientation_conformity_ok(
                summary.get("tet_hex_pyramid_interface_face_orientation_conformity_identity")
            )
        ),
        "hex_sweep_uses_current_source_target_vertex_map_and_twist_path": (
            _hex_sweep_face_vertex_twist_generation_ok(
                summary.get("hex_sweep_source_target_face_vertex_twist_generation_identity")
            )
        ),
        "quality_histogram_uses_current_metric_element_set_and_units": (
            _quality_histogram_metric_element_unit_generation_ok(
                summary.get("quality_histogram_metric_element_set_unit_generation_identity")
            )
        ),
        "step_ap214_multibody_export_uses_current_body_unit_frame_and_mass_properties": (
            _step_ap214_multibody_export_generation_ok(
                summary.get(
                    "step_ap214_body_name_unit_frame_mass_property_generation_identity"
                )
            )
        ),
        "hybrid_transition_uses_current_orientation_shared_nodes_and_blocks": (
            _hybrid_transition_topology_block_generation_ok(
                summary.get(
                    "hybrid_tet_hex_pyramid_transition_topology_block_generation_identity"
                )
            )
        ),
        "parallel_sculpt_assembly_uses_current_partitions_ghosts_and_refinement": (
            _parallel_sculpt_partition_identity_ok(
                summary.get(
                    "parallel_sculpt_partition_ghost_refinement_generation_identity"
                )
            )
        ),
        "mixed_transition_export_uses_current_interfaces_families_blocks_and_sidesets": (
            _mixed_transition_interface_ownership_ok(
                summary.get(
                    "mixed_transition_interface_block_sideset_generation_identity"
                )
            )
        ),
        "high_order_hex_uses_current_jacobians_transform_node_order_and_blocks": (
            _high_order_hex_jacobian_identity_ok(
                summary.get(
                    "high_order_hex_jacobian_transform_node_order_block_generation_identity"
                )
            )
        ),
        "mesh_boundary_closure_uses_current_cad_units_frame_and_geometry": (
            _mesh_cad_closure_identity_ok(
                summary.get(
                    "mesh_boundary_cad_volume_area_unit_frame_generation_identity"
                )
            )
        ),
        "hybrid_interfaces_use_current_families_nodes_orientation_and_blocks": (
            _hybrid_interface_conformity_identity_ok(
                summary.get(
                    "hybrid_tet_hex_pyramid_interface_conformity_orientation_block_identity"
                )
            )
        ),
        "periodic_sidesets_use_current_node_pairs_transform_tolerance_and_geometry": (
            _periodic_sideset_pairing_identity_ok(
                summary.get(
                    "periodic_sideset_node_pair_transform_tolerance_geometry_generation_identity"
                )
            )
        ),
        "hex_sweeps_use_current_faces_layers_correspondence_jacobians_and_blocks": (
            _hex_sweep_layer_correspondence_identity_ok(
                summary.get(
                    "hex_sweep_source_target_layer_correspondence_jacobian_block_generation_identity"
                )
            )
        ),
        "high_order_hex_exports_use_current_edge_face_interior_order_curvature_and_jacobian": (
            _high_order_hex_export_identity_ok(
                summary.get(
                    "high_order_hex_edge_face_interior_node_curvature_jacobian_export_generation_identity"
                )
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
        "journal_replay_uses_current_transaction_undo_and_entity_ids": (
            _journal_transaction_undo_entity_reuse_generation_ok(
                summary.get("journal_transaction_undo_entity_id_reuse_generation_identity")
            )
        ),
        "netgen_export_uses_current_block_order_and_curving_generation": (
            _netgen_vol_block_order_curving_generation_ok(
                summary.get("netgen_vol_element_block_order_curving_generation_identity")
            )
        ),
        "block_sideset_groups_use_current_merge_and_renumber_generation": (
            _block_sideset_merge_renumber_generation_ok(
                summary.get("block_sideset_group_entity_merge_renumber_generation_identity")
            )
        ),
        "aprepro_replay_uses_current_variables_includes_and_working_directory": (
            _aprepro_include_variable_transaction_generation_ok(
                summary.get("aprepro_include_variable_expansion_working_directory_generation_identity")
            )
        ),
        "headless_step_export_uses_current_bodies_transforms_names_and_log": (
            _headless_step_export_generation_ok(
                summary.get(
                    "headless_step_export_body_transform_name_generation_identity"
                )
            )
        ),
        "healed_geometry_ownership_uses_current_tolerance_imprint_and_merge": (
            _geometry_heal_imprint_merge_ownership_generation_ok(
                summary.get(
                    "geometry_heal_tolerance_imprint_merge_ownership_generation_identity"
                )
            )
        ),
        "journal_replay_uses_current_geometry_entity_map_version_and_command_log": (
            _journal_replay_entity_version_identity_ok(
                summary.get(
                    "journal_replay_geometry_entity_map_version_generation_identity"
                )
            )
        ),
        "exodus64_decode_uses_current_entity_sideset_element_map_and_schema": (
            _exodus64_mapping_identity_ok(
                summary.get(
                    "exodus64_entity_sideset_element_map_schema_generation_identity"
                )
            )
        ),
        "webcut_imprint_merge_uses_current_tolerances_topology_and_entities": (
            _webcut_imprint_merge_topology_identity_ok(
                summary.get(
                    "webcut_imprint_merge_tolerance_topology_entity_generation_identity"
                )
            )
        ),
        "exodus_qa_uses_current_coordinates_distribution_factors_and_checksum": (
            _exodus_qa_coordinate_distribution_identity_ok(
                summary.get(
                    "exodus_qa_coordinate_distribution_checksum_generation_identity"
                )
            )
        ),
        "journal_reset_replay_uses_current_session_commands_entities_and_topology": (
            _journal_reset_replay_identity_ok(
                summary.get(
                    "journal_reset_entity_id_reuse_undo_replay_session_generation_identity"
                )
            )
        ),
        "netgen_vol_export_uses_p1_tri_tet_boundaries_materials_and_source_mesh": (
            _netgen_vol_export_identity_ok(
                summary.get(
                    "netgen_vol_export_family_order_boundary_material_checksum_identity"
                )
            )
        ),
        "sculpt_outputs_use_current_voxels_thresholds_material_blocks_and_session": (
            _sculpt_input_output_identity_ok(
                summary.get(
                    "sculpt_voxel_spacing_threshold_material_block_output_session_generation_identity"
                )
            )
        ),
        "exodus_merges_use_current_tolerance_nodes_global_ids_blocks_sidesets_and_checksums": (
            _exodus_merge_identity_ok(
                summary.get(
                    "exodus_merge_node_tolerance_global_id_block_sideset_checksum_generation_identity"
                )
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
