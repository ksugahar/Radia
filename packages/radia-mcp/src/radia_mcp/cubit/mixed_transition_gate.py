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


def _sheet_pillow_layer_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        elements = [int(value) for value in identity.get("inserted_layer_element_ids", [])]
        result_elements = [int(value) for value in identity.get("result_inserted_layer_element_ids", [])]
        interfaces = [[int(value) for value in row] for row in identity.get("block_interface_pairs", [])]
        result_interfaces = [[int(value) for value in row] for row in identity.get("result_block_interface_pairs", [])]
        signs = [int(value) for value in identity.get("orientation_signs", [])]
        result_signs = [int(value) for value in identity.get("result_orientation_signs", [])]
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [float(value) for value in identity.get("result_scaled_jacobians", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("layer_generation") or "")
    operation = str(identity.get("operation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "topology_layer_generation", "block_layer_generation",
            "interface_layer_generation", "orientation_layer_generation",
            "jacobian_layer_generation", "result_layer_generation",
        ))
        and operation in {"sheet", "pillow"}
        and identity.get("result_operation") == operation
        and bool(elements)
        and all(value > 0 for value in elements)
        and len(set(elements)) == len(elements)
        and result_elements == elements
        and bool(interfaces)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 and row[0] != row[1] for row in interfaces)
        and result_interfaces == interfaces
        and len(signs) == len(elements)
        and all(value == 1 for value in signs)
        and result_signs == signs
        and len(jacobians) == len(elements)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and _valid_sha256(identity.get("layer_topology_sha256"))
        and identity.get("result_layer_topology_sha256") == identity.get("layer_topology_sha256")
        and _valid_sha256(identity.get("interface_map_sha256"))
        and identity.get("result_interface_map_sha256") == identity.get("interface_map_sha256")
    )


def _pyramid_transition_export_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        base = [int(value) for value in identity.get("pyramid_base_node_order", [])]
        result_base = [int(value) for value in identity.get("result_pyramid_base_node_order", [])]
        sides = [int(value) for value in identity.get("pyramid_side_orientations", [])]
        result_sides = [int(value) for value in identity.get("result_pyramid_side_orientations", [])]
        nodes = [int(value) for value in identity.get("interface_node_ids", [])]
        result_nodes = [int(value) for value in identity.get("result_interface_node_ids", [])]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        blocks = [int(value) for value in identity.get("block_ids", [])]
        result_blocks = [int(value) for value in identity.get("result_block_ids", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transition_generation") or "")
    families = [str(value) for value in identity.get("adjacent_element_families", [])]
    result_families = [str(value) for value in identity.get("result_adjacent_element_families", [])]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "base_transition_generation", "side_transition_generation",
            "interface_transition_generation", "jacobian_transition_generation",
            "block_transition_generation", "export_transition_generation",
            "result_transition_generation",
        ))
        and len(base) == 4
        and all(value > 0 for value in base)
        and len(set(base)) == 4
        and result_base == base
        and len(sides) == 4
        and all(value == 1 for value in sides)
        and result_sides == sides
        and len(nodes) == 5
        and all(value > 0 for value in nodes)
        and len(set(nodes)) == 5
        and result_nodes == nodes
        and families == ["hex8", "pyramid5", "tet4"]
        and result_families == families
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and result_jacobian == jacobian
        and len(blocks) == 3
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == 3
        and result_blocks == blocks
        and _valid_sha256(identity.get("transition_export_sha256"))
        and identity.get("result_transition_export_sha256") == identity.get("transition_export_sha256")
    )


def _journal_invocation_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("invocation_generation") or "")
    includes = [str(value) for value in identity.get("include_order", [])]
    scopes = [[str(value) for value in row] for row in identity.get("aprepro_scope", [])]
    workdir = str(identity.get("working_directory") or "")
    version = str(identity.get("cubit_version") or "")
    flags = [str(value) for value in identity.get("headless_flags", [])]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "include_invocation_generation", "aprepro_invocation_generation",
            "workdir_invocation_generation", "version_invocation_generation",
            "output_invocation_generation", "result_invocation_generation",
        ))
        and bool(includes)
        and len(set(includes)) == len(includes)
        and identity.get("result_include_order") == includes
        and bool(scopes)
        and all(len(row) == 2 and all(row) for row in scopes)
        and len({row[0] for row in scopes}) == len(scopes)
        and identity.get("result_aprepro_scope") == scopes
        and bool(workdir)
        and identity.get("result_working_directory") == workdir
        and bool(version)
        and identity.get("result_cubit_version") == version
        and set(flags) == {"-nographics", "-batch"}
        and identity.get("result_headless_flags") == flags
        and _valid_sha256(identity.get("journal_sha256"))
        and identity.get("result_journal_sha256") == identity.get("journal_sha256")
        and _valid_sha256(identity.get("output_sha256"))
        and identity.get("result_output_sha256") == identity.get("output_sha256")
    )


def _exodus_result_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        width = int(identity.get("integer_width_bits"))
        decoded_width = int(identity.get("decoded_integer_width_bits"))
        ids = [int(value) for value in identity.get("global_node_ids", [])]
        decoded_ids = [int(value) for value in identity.get("decoded_global_node_ids", [])]
        times = [float(value) for value in identity.get("time_steps_s", [])]
        decoded_times = [float(value) for value in identity.get("decoded_time_steps_s", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("exodus_generation") or "")
    qa = [[str(value) for value in row] for row in identity.get("qa_records", [])]
    variables = [str(value) for value in identity.get("nodal_variable_order", [])]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "id_exodus_generation", "qa_exodus_generation", "time_exodus_generation",
            "variable_exodus_generation", "mesh_exodus_generation",
            "result_exodus_generation",
        ))
        and width == 64
        and decoded_width == width
        and bool(ids)
        and all(value > 2**31 - 1 for value in ids)
        and len(set(ids)) == len(ids)
        and decoded_ids == ids
        and bool(qa)
        and all(len(row) == 4 and all(row) for row in qa)
        and identity.get("decoded_qa_records") == qa
        and len(times) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(right > left for left, right in zip(times, times[1:]))
        and decoded_times == times
        and bool(variables)
        and all(variables)
        and len(set(variables)) == len(variables)
        and identity.get("decoded_nodal_variable_order") == variables
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("decoded_mesh_sha256") == identity.get("mesh_sha256")
        and _valid_sha256(identity.get("exodus_sha256"))
        and identity.get("decoded_exodus_sha256") == identity.get("exodus_sha256")
    )


def _spine_sweep_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        spine_ids = [int(value) for value in identity.get("spine_curve_ids", [])]
        result_spine_ids = [int(value) for value in identity.get("result_spine_curve_ids", [])]
        twists = [float(value) for value in identity.get("twist_angles_deg", [])]
        result_twists = [float(value) for value in identity.get("result_twist_angles_deg", [])]
        intervals = int(identity.get("interval_count"))
        result_intervals = int(identity.get("result_interval_count"))
        source = int(identity.get("source_surface_id"))
        result_source = int(identity.get("result_source_surface_id"))
        target = int(identity.get("target_surface_id"))
        result_target = int(identity.get("result_target_surface_id"))
        block = int(identity.get("block_id"))
        result_block = int(identity.get("result_block_id"))
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [float(value) for value in identity.get("result_scaled_jacobians", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    frame = str(identity.get("frame_method") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "spine_sweep_generation", "frame_sweep_generation",
            "twist_sweep_generation", "interval_sweep_generation",
            "surface_sweep_generation", "block_sweep_generation",
            "quality_sweep_generation", "result_sweep_generation",
        ))
        and bool(spine_ids)
        and all(value > 0 for value in spine_ids)
        and len(set(spine_ids)) == len(spine_ids)
        and result_spine_ids == spine_ids
        and frame in {"parallel_transport", "frenet", "fixed"}
        and identity.get("result_frame_method") == frame
        and len(twists) >= 2
        and all(math.isfinite(value) for value in twists)
        and result_twists == twists
        and intervals > 0
        and result_intervals == intervals
        and source > 0
        and target > 0
        and source != target
        and result_source == source
        and result_target == target
        and block > 0
        and result_block == block
        and bool(jacobians)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and _valid_sha256(identity.get("sweep_mesh_sha256"))
        and identity.get("result_sweep_mesh_sha256") == identity.get("sweep_mesh_sha256")
    )


def _local_refinement_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        parents = [int(value) for value in identity.get("parent_element_ids", [])]
        result_parents = [int(value) for value in identity.get("result_parent_element_ids", [])]
        children = [[int(value) for value in row] for row in identity.get("child_parent_pairs", [])]
        result_children = [[int(value) for value in row] for row in identity.get("result_child_parent_pairs", [])]
        transitions = [[int(value) for value in row] for row in identity.get("transition_face_pairs", [])]
        result_transitions = [[int(value) for value in row] for row in identity.get("result_transition_face_pairs", [])]
        conformity = [[int(value) for value in row] for row in identity.get("conformity_node_pairs", [])]
        result_conformity = [[int(value) for value in row] for row in identity.get("result_conformity_node_pairs", [])]
        blocks = [int(value) for value in identity.get("block_ids", [])]
        result_blocks = [int(value) for value in identity.get("result_block_ids", [])]
        boundaries = [int(value) for value in identity.get("boundary_sideset_ids", [])]
        result_boundaries = [int(value) for value in identity.get("result_boundary_sideset_ids", [])]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("refinement_generation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "parent_refinement_generation", "child_refinement_generation",
            "transition_refinement_generation", "conformity_refinement_generation",
            "block_refinement_generation", "boundary_refinement_generation",
            "jacobian_refinement_generation", "export_refinement_generation",
            "result_refinement_generation",
        ))
        and bool(parents)
        and all(value > 0 for value in parents)
        and len(set(parents)) == len(parents)
        and result_parents == parents
        and bool(children)
        and all(len(row) == 2 and row[0] > 0 and row[1] in parents for row in children)
        and len({row[0] for row in children}) == len(children)
        and result_children == children
        and bool(transitions)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 and row[0] != row[1] for row in transitions)
        and result_transitions == transitions
        and bool(conformity)
        and all(len(row) == 2 and row[0] > 0 and row[0] == row[1] for row in conformity)
        and result_conformity == conformity
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and result_blocks == blocks
        and bool(boundaries)
        and all(value > 0 for value in boundaries)
        and len(set(boundaries)) == len(boundaries)
        and result_boundaries == boundaries
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and result_jacobian == jacobian
        and _valid_sha256(identity.get("refined_mesh_sha256"))
        and identity.get("result_refined_mesh_sha256") == identity.get("refined_mesh_sha256")
    )


def _sideset_export_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        dimension = int(identity.get("skin_entity_dimension"))
        decoded_dimension = int(identity.get("decoded_skin_entity_dimension"))
        faces = [int(value) for value in identity.get("skin_face_ids", [])]
        decoded_faces = [int(value) for value in identity.get("decoded_skin_face_ids", [])]
        signs = [int(value) for value in identity.get("orientation_signs", [])]
        decoded_signs = [int(value) for value in identity.get("decoded_orientation_signs", [])]
        sideset = int(identity.get("sideset_id"))
        decoded_sideset = int(identity.get("decoded_sideset_id"))
        owners = [int(value) for value in identity.get("geometric_owner_surface_ids", [])]
        decoded_owners = [int(value) for value in identity.get("decoded_geometric_owner_surface_ids", [])]
        exodus_pairs = [[int(value) for value in row] for row in identity.get("exodus_element_side_pairs", [])]
        decoded_exodus_pairs = [[int(value) for value in row] for row in identity.get("decoded_exodus_element_side_pairs", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sideset_generation") or "")
    namespace = str(identity.get("id_namespace") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "skin_sideset_generation", "dimension_sideset_generation",
            "orientation_sideset_generation", "namespace_sideset_generation",
            "owner_sideset_generation", "exodus_sideset_generation",
            "mesh_sideset_generation", "result_sideset_generation",
        ))
        and dimension == 2
        and decoded_dimension == dimension
        and bool(faces)
        and all(value > 0 for value in faces)
        and len(set(faces)) == len(faces)
        and decoded_faces == faces
        and len(signs) == len(faces)
        and all(value in {-1, 1} for value in signs)
        and decoded_signs == signs
        and namespace == "sideset"
        and identity.get("decoded_id_namespace") == namespace
        and sideset > 0
        and decoded_sideset == sideset
        and bool(owners)
        and all(value > 0 for value in owners)
        and len(set(owners)) == len(owners)
        and decoded_owners == owners
        and len(exodus_pairs) == len(faces)
        and all(len(row) == 2 and row[0] > 0 and 1 <= row[1] <= 6 for row in exodus_pairs)
        and decoded_exodus_pairs == exodus_pairs
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("decoded_mesh_sha256") == identity.get("mesh_sha256")
        and _valid_sha256(identity.get("exodus_sha256"))
        and identity.get("decoded_exodus_sha256") == identity.get("exodus_sha256")
    )


def _headless_python_invocation_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        undo_before = int(identity.get("undo_depth_before"))
        result_undo_before = int(identity.get("result_undo_depth_before"))
        undo_after = int(identity.get("undo_depth_after"))
        result_undo_after = int(identity.get("result_undo_depth_after"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("invocation_generation") or "")
    executable = str(identity.get("interpreter_executable") or "")
    python_version = str(identity.get("python_version") or "")
    module_version = str(identity.get("cubit_module_version") or "")
    flags = [str(value) for value in identity.get("headless_flags", [])]
    transaction = str(identity.get("command_transaction_id") or "")
    paths = [str(value) for value in identity.get("output_paths", [])]
    digests = [str(value) for value in identity.get("output_sha256", [])]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "interpreter_invocation_generation", "module_invocation_generation",
            "transaction_invocation_generation", "undo_invocation_generation",
            "output_invocation_generation", "result_invocation_generation",
        ))
        and bool(executable)
        and identity.get("result_interpreter_executable") == executable
        and bool(python_version)
        and identity.get("result_python_version") == python_version
        and bool(module_version)
        and identity.get("result_cubit_module_version") == module_version
        and set(flags) == {"-nographics", "-batch"}
        and identity.get("result_headless_flags") == flags
        and bool(transaction)
        and identity.get("result_command_transaction_id") == transaction
        and _valid_sha256(identity.get("command_log_sha256"))
        and identity.get("result_command_log_sha256") == identity.get("command_log_sha256")
        and undo_before >= 0
        and undo_after == undo_before
        and result_undo_before == undo_before
        and result_undo_after == undo_after
        and bool(paths)
        and len(set(paths)) == len(paths)
        and identity.get("result_output_paths") == paths
        and len(digests) == len(paths)
        and all(_valid_sha256(value) for value in digests)
        and identity.get("result_output_sha256") == digests
    )


def _sweep_hex_twist_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        twist = float(identity.get("twist_angle_deg"))
        result_twist = float(identity.get("result_twist_angle_deg"))
        source_face = int(identity.get("source_face_id"))
        result_source_face = int(identity.get("result_source_face_id"))
        target_face = int(identity.get("target_face_id"))
        result_target_face = int(identity.get("result_target_face_id"))
        intervals = int(identity.get("axial_interval_count"))
        result_intervals = int(identity.get("result_axial_interval_count"))
        source_edges = [int(value) for value in identity.get("source_edge_ids", [])]
        result_source_edges = [
            int(value) for value in identity.get("result_source_edge_ids", [])
        ]
        target_edges = [int(value) for value in identity.get("target_edge_ids", [])]
        result_target_edges = [
            int(value) for value in identity.get("result_target_edge_ids", [])
        ]
        edge_map = [
            [int(value) for value in pair] for pair in identity.get("edge_map", [])
        ]
        result_edge_map = [
            [int(value) for value in pair]
            for pair in identity.get("result_edge_map", [])
        ]
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [
            float(value) for value in identity.get("result_scaled_jacobians", [])
        ]
        block_id = int(identity.get("block_id"))
        result_block_id = int(identity.get("result_block_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "twist_sweep_generation",
                "pairing_sweep_generation",
                "interval_sweep_generation",
                "edge_map_sweep_generation",
                "orientation_sweep_generation",
                "quality_sweep_generation",
                "block_sweep_generation",
                "result_sweep_generation",
            )
        )
        and math.isfinite(twist)
        and result_twist == twist
        and source_face > 0
        and target_face > 0
        and source_face != target_face
        and result_source_face == source_face
        and result_target_face == target_face
        and intervals > 0
        and result_intervals == intervals
        and bool(source_edges)
        and len(source_edges) == len(target_edges)
        and len(set(source_edges)) == len(source_edges)
        and len(set(target_edges)) == len(target_edges)
        and all(value > 0 for value in source_edges + target_edges)
        and result_source_edges == source_edges
        and result_target_edges == target_edges
        and edge_map == [list(pair) for pair in zip(source_edges, target_edges)]
        and result_edge_map == edge_map
        and identity.get("sweep_orientation")
        in {"right_handed_source_to_target", "left_handed_source_to_target"}
        and identity.get("result_sweep_orientation")
        == identity.get("sweep_orientation")
        and bool(jacobians)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and block_id > 0
        and result_block_id == block_id
        and _valid_sha256(identity.get("sweep_mesh_sha256"))
        and identity.get("result_sweep_mesh_sha256")
        == identity.get("sweep_mesh_sha256")
    )


def _mixed_transition_face_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        hex_ids = [int(value) for value in identity.get("hex_element_ids", [])]
        result_hex_ids = [
            int(value) for value in identity.get("result_hex_element_ids", [])
        ]
        pyramid_ids = [
            int(value) for value in identity.get("pyramid_element_ids", [])
        ]
        result_pyramid_ids = [
            int(value) for value in identity.get("result_pyramid_element_ids", [])
        ]
        tet_ids = [int(value) for value in identity.get("tet_element_ids", [])]
        result_tet_ids = [
            int(value) for value in identity.get("result_tet_element_ids", [])
        ]
        faces = [
            [int(value) for value in face]
            for face in identity.get("transition_face_node_ids", [])
        ]
        result_faces = [
            [int(value) for value in face]
            for face in identity.get("result_transition_face_node_ids", [])
        ]
        owners = [
            [str(value) for value in pair]
            for pair in identity.get("face_owner_pairs", [])
        ]
        result_owners = [
            [str(value) for value in pair]
            for pair in identity.get("result_face_owner_pairs", [])
        ]
        orientations = [
            [int(value) for value in pair]
            for pair in identity.get("opposed_face_orientation_signs", [])
        ]
        result_orientations = [
            [int(value) for value in pair]
            for pair in identity.get("result_opposed_face_orientation_signs", [])
        ]
        unmatched = int(identity.get("unmatched_transition_face_count"))
        result_unmatched = int(identity.get("result_unmatched_transition_face_count"))
        quality = float(identity.get("minimum_scaled_jacobian"))
        result_quality = float(identity.get("result_minimum_scaled_jacobian"))
        block_ids = [int(value) for value in identity.get("block_ids", [])]
        result_block_ids = [int(value) for value in identity.get("result_block_ids", [])]
        sideset_ids = [int(value) for value in identity.get("sideset_ids", [])]
        result_sideset_ids = [
            int(value) for value in identity.get("result_sideset_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transition_generation") or "")
    all_elements = hex_ids + pyramid_ids + tet_ids
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "face_transition_generation",
                "node_transition_generation",
                "orientation_transition_generation",
                "conformity_transition_generation",
                "quality_transition_generation",
                "block_transition_generation",
                "sideset_transition_generation",
                "export_transition_generation",
                "result_transition_generation",
            )
        )
        and bool(hex_ids and pyramid_ids and tet_ids)
        and all(value > 0 for value in all_elements)
        and len(set(all_elements)) == len(all_elements)
        and result_hex_ids == hex_ids
        and result_pyramid_ids == pyramid_ids
        and result_tet_ids == tet_ids
        and bool(faces)
        and all(len(face) in {3, 4} and len(set(face)) == len(face) for face in faces)
        and result_faces == faces
        and len(owners) == len(faces)
        and all(len(pair) == 2 and all(pair) for pair in owners)
        and result_owners == owners
        and len(orientations) == len(faces)
        and all(len(pair) == 2 and pair[0] == -pair[1] for pair in orientations)
        and result_orientations == orientations
        and unmatched == 0
        and result_unmatched == unmatched
        and math.isfinite(quality)
        and quality > 0.0
        and result_quality == quality
        and bool(block_ids)
        and all(value > 0 for value in block_ids)
        and len(set(block_ids)) == len(block_ids)
        and result_block_ids == block_ids
        and bool(sideset_ids)
        and all(value > 0 for value in sideset_ids)
        and len(set(sideset_ids)) == len(sideset_ids)
        and result_sideset_ids == sideset_ids
        and _valid_sha256(identity.get("transition_mesh_sha256"))
        and identity.get("result_transition_mesh_sha256")
        == identity.get("transition_mesh_sha256")
        and _valid_sha256(identity.get("transition_export_sha256"))
        and identity.get("accepted_transition_export_sha256")
        == identity.get("transition_export_sha256")
    )


def _journal_transaction_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        commands = [str(value) for value in identity.get("command_sequence", [])]
        replayed_commands = [
            str(value) for value in identity.get("replayed_command_sequence", [])
        ]
        ordinals = [int(value) for value in identity.get("command_ordinals", [])]
        replayed_ordinals = [
            int(value) for value in identity.get("replayed_command_ordinals", [])
        ]
        entities = [int(value) for value in identity.get("entity_ids", [])]
        replayed_entities = [
            int(value) for value in identity.get("replayed_entity_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("journal_generation") or "")
    session = str(identity.get("active_session_id") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "undo_journal_generation",
                "command_journal_generation",
                "entity_journal_generation",
                "session_journal_generation",
                "model_journal_generation",
                "digest_journal_generation",
                "result_journal_generation",
            )
        )
        and bool(str(identity.get("undo_transaction_id") or ""))
        and identity.get("replayed_undo_transaction_id")
        == identity.get("undo_transaction_id")
        and bool(commands)
        and all(commands)
        and replayed_commands == commands
        and ordinals == list(range(1, len(commands) + 1))
        and replayed_ordinals == ordinals
        and bool(entities)
        and all(value > 0 for value in entities)
        and len(set(entities)) == len(entities)
        and replayed_entities == entities
        and session.startswith("headless-session-")
        and identity.get("replayed_active_session_id") == session
        and bool(str(identity.get("model_generation_id") or ""))
        and identity.get("replayed_model_generation_id")
        == identity.get("model_generation_id")
        and _valid_sha256(identity.get("journal_sha256"))
        and identity.get("loaded_journal_sha256") == identity.get("journal_sha256")
        and _valid_sha256(identity.get("journal_result_sha256"))
        and identity.get("accepted_journal_result_sha256")
        == identity.get("journal_result_sha256")
    )


def _exodus_sideset_distribution_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sidesets = [int(value) for value in identity.get("sideset_ids", [])]
        decoded_sidesets = [
            int(value) for value in identity.get("decoded_sideset_ids", [])
        ]
        topologies = [str(value) for value in identity.get("face_topologies", [])]
        decoded_topologies = [
            str(value) for value in identity.get("decoded_face_topologies", [])
        ]
        orientations = [
            int(value) for value in identity.get("face_orientation_signs", [])
        ]
        decoded_orientations = [
            int(value) for value in identity.get("decoded_face_orientation_signs", [])
        ]
        factors = [
            [float(value) for value in row]
            for row in identity.get("distribution_factors", [])
        ]
        decoded_factors = [
            [float(value) for value in row]
            for row in identity.get("decoded_distribution_factors", [])
        ]
        blocks = [int(value) for value in identity.get("block_owner_ids", [])]
        decoded_blocks = [
            int(value) for value in identity.get("decoded_block_owner_ids", [])
        ]
        time_index = int(identity.get("time_step_index"))
        decoded_time_index = int(identity.get("decoded_time_step_index"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sideset_generation") or "")
    nodes_per_face = {"tri3": 3, "quad4": 4, "tri6": 6, "quad8": 8}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "topology_sideset_generation",
                "orientation_sideset_generation",
                "factor_sideset_generation",
                "block_sideset_generation",
                "time_sideset_generation",
                "file_sideset_generation",
                "result_sideset_generation",
            )
        )
        and bool(sidesets)
        and all(value > 0 for value in sidesets)
        and len(set(sidesets)) == len(sidesets)
        and decoded_sidesets == sidesets
        and len(topologies) == len(sidesets)
        and all(value in nodes_per_face for value in topologies)
        and decoded_topologies == topologies
        and len(orientations) == len(sidesets)
        and all(value in {-1, 1} for value in orientations)
        and decoded_orientations == orientations
        and len(factors) == len(sidesets)
        and all(
            len(row) == nodes_per_face[topology]
            and all(math.isfinite(value) and value >= 0.0 for value in row)
            for row, topology in zip(factors, topologies)
        )
        and decoded_factors == factors
        and len(blocks) == len(sidesets)
        and all(value > 0 for value in blocks)
        and decoded_blocks == blocks
        and time_index >= 1
        and decoded_time_index == time_index
        and _valid_sha256(identity.get("exodus_file_sha256"))
        and identity.get("decoded_exodus_file_sha256")
        == identity.get("exodus_file_sha256")
        and _valid_sha256(identity.get("sideset_table_sha256"))
        and identity.get("decoded_sideset_table_sha256")
        == identity.get("sideset_table_sha256")
    )


def _hex_map_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        curves = [int(value) for value in identity.get("boundary_curve_ids", [])]
        result_curves = [int(value) for value in identity.get("result_boundary_curve_ids", [])]
        intervals = [int(value) for value in identity.get("curve_intervals", [])]
        result_intervals = [int(value) for value in identity.get("result_curve_intervals", [])]
        vertices = [int(value) for value in identity.get("corner_vertex_ids", [])]
        result_vertices = [int(value) for value in identity.get("result_corner_vertex_ids", [])]
        valences = [int(value) for value in identity.get("corner_valences", [])]
        result_valences = [int(value) for value in identity.get("result_corner_valences", [])]
        pairs = [[int(value) for value in row] for row in identity.get("source_target_face_pairs", [])]
        result_pairs = [[int(value) for value in row] for row in identity.get("result_source_target_face_pairs", [])]
        logical = [[int(value) for value in row] for row in identity.get("logical_corner_coordinates", [])]
        result_logical = [[int(value) for value in row] for row in identity.get("result_logical_corner_coordinates", [])]
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [float(value) for value in identity.get("result_scaled_jacobians", [])]
        block = int(identity.get("block_id"))
        result_block = int(identity.get("result_block_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("map_generation") or "")
    logical_cube = {(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)}
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "interval_map_generation", "valence_map_generation", "pairing_map_generation",
            "logical_map_generation", "jacobian_map_generation", "block_map_generation",
            "result_map_generation",
        ))
        and len(curves) == len(intervals) == 4
        and all(value > 0 for value in curves)
        and len(set(curves)) == 4
        and result_curves == curves
        and all(value > 0 and value % 2 == 0 for value in intervals)
        and intervals[0] == intervals[2]
        and intervals[1] == intervals[3]
        and result_intervals == intervals
        and len(vertices) == len(valences) == 8
        and len(set(vertices)) == 8
        and all(value > 0 for value in vertices)
        and result_vertices == vertices
        and valences == [3] * 8
        and result_valences == valences
        and bool(pairs)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 and row[0] != row[1] for row in pairs)
        and result_pairs == pairs
        and len(logical) == 8
        and {tuple(row) for row in logical} == logical_cube
        and result_logical == logical
        and bool(jacobians)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and block > 0
        and result_block == block
        and _valid_sha256(identity.get("map_mesh_sha256"))
        and identity.get("result_map_mesh_sha256") == identity.get("map_mesh_sha256")
    )


def _mesh_morph_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        nodes = [int(value) for value in identity.get("boundary_node_ids", [])]
        result_nodes = [int(value) for value in identity.get("result_boundary_node_ids", [])]
        displacements = [[float(value) for value in row] for row in identity.get("boundary_displacements", [])]
        result_displacements = [[float(value) for value in row] for row in identity.get("result_boundary_displacements", [])]
        fixed = [int(value) for value in identity.get("fixed_node_ids", [])]
        result_fixed = [int(value) for value in identity.get("result_fixed_node_ids", [])]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        sidesets = [int(value) for value in identity.get("sideset_ids", [])]
        result_sidesets = [int(value) for value in identity.get("result_sideset_ids", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("morph_generation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "boundary_morph_generation", "constraint_morph_generation", "frame_morph_generation",
            "smoothing_morph_generation", "jacobian_morph_generation", "sideset_morph_generation",
            "export_morph_generation", "result_morph_generation",
        ))
        and bool(nodes)
        and all(value > 0 for value in nodes)
        and len(set(nodes)) == len(nodes)
        and result_nodes == nodes
        and len(displacements) == len(nodes)
        and all(len(row) == 3 and all(math.isfinite(value) for value in row) for row in displacements)
        and result_displacements == displacements
        and bool(fixed)
        and all(value > 0 for value in fixed)
        and len(set(fixed)) == len(fixed)
        and not set(fixed).intersection(nodes)
        and result_fixed == fixed
        and identity.get("coordinate_frame") == "global_cartesian"
        and identity.get("result_coordinate_frame") == identity.get("coordinate_frame")
        and identity.get("interior_smoothing") == "winslow_volume"
        and identity.get("result_interior_smoothing") == identity.get("interior_smoothing")
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and result_jacobian == jacobian
        and bool(sidesets)
        and len(set(sidesets)) == len(sidesets)
        and result_sidesets == sidesets
        and _valid_sha256(identity.get("morphed_mesh_sha256"))
        and identity.get("result_morphed_mesh_sha256") == identity.get("morphed_mesh_sha256")
        and _valid_sha256(identity.get("export_sha256"))
        and identity.get("accepted_export_sha256") == identity.get("export_sha256")
    )


def _exodus_truth_table_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        variables = [str(value) for value in identity.get("element_variable_names", [])]
        decoded_variables = [str(value) for value in identity.get("decoded_element_variable_names", [])]
        blocks = [int(value) for value in identity.get("block_ids", [])]
        decoded_blocks = [int(value) for value in identity.get("decoded_block_ids", [])]
        table = [[int(value) for value in row] for row in identity.get("truth_table", [])]
        decoded_table = [[int(value) for value in row] for row in identity.get("decoded_truth_table", [])]
        times = [int(value) for value in identity.get("time_step_indices", [])]
        decoded_times = [int(value) for value in identity.get("decoded_time_step_indices", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("truth_generation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "variable_truth_generation", "block_truth_generation", "table_truth_generation",
            "timestep_truth_generation", "layout_truth_generation", "file_truth_generation",
            "result_truth_generation",
        ))
        and bool(variables)
        and all(variables)
        and len(set(variables)) == len(variables)
        and decoded_variables == variables
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and decoded_blocks == blocks
        and len(table) == len(blocks)
        and all(len(row) == len(variables) and set(row).issubset({0, 1}) for row in table)
        and decoded_table == table
        and bool(times)
        and all(value > 0 for value in times)
        and all(right > left for left, right in zip(times, times[1:]))
        and decoded_times == times
        and identity.get("value_layout") == "time_block_variable_element"
        and identity.get("decoded_value_layout") == identity.get("value_layout")
        and _valid_sha256(identity.get("exodus_sha256"))
        and identity.get("decoded_exodus_sha256") == identity.get("exodus_sha256")
        and _valid_sha256(identity.get("truth_table_sha256"))
        and identity.get("decoded_truth_table_sha256") == identity.get("truth_table_sha256")
    )


def _cad_import_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        bodies = [int(value) for value in identity.get("body_ids", [])]
        result_bodies = [int(value) for value in identity.get("result_body_ids", [])]
        sheets = [int(value) for value in identity.get("sheet_ids", [])]
        result_sheets = [int(value) for value in identity.get("result_sheet_ids", [])]
        lumps = [int(value) for value in identity.get("lump_ids", [])]
        result_lumps = [int(value) for value in identity.get("result_lump_ids", [])]
        matrix = [[float(value) for value in row] for row in identity.get("placement_matrix", [])]
        result_matrix = [[float(value) for value in row] for row in identity.get("result_placement_matrix", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("cad_generation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "classification_cad_generation", "transform_cad_generation", "unit_cad_generation",
            "heal_cad_generation", "topology_cad_generation", "session_cad_generation",
            "result_cad_generation",
        ))
        and bool(bodies)
        and all(value > 0 for value in bodies + sheets + lumps)
        and len(set(bodies + sheets)) == len(bodies) + len(sheets)
        and result_bodies == bodies
        and result_sheets == sheets
        and result_lumps == lumps
        and len(matrix) == 4
        and all(len(row) == 4 and all(math.isfinite(value) for value in row) for row in matrix)
        and matrix[3] == [0.0, 0.0, 0.0, 1.0]
        and result_matrix == matrix
        and identity.get("length_unit") in {"m", "mm", "cm", "in"}
        and identity.get("result_length_unit") == identity.get("length_unit")
        and all(bool(str(identity.get(key) or "")) for key in (
            "heal_transaction_id", "topology_revision", "cubit_session_id"
        ))
        and identity.get("result_heal_transaction_id") == identity.get("heal_transaction_id")
        and identity.get("result_topology_revision") == identity.get("topology_revision")
        and identity.get("result_cubit_session_id") == identity.get("cubit_session_id")
        and _valid_sha256(identity.get("cad_sha256"))
        and identity.get("loaded_cad_sha256") == identity.get("cad_sha256")
        and _valid_sha256(identity.get("import_result_sha256"))
        and identity.get("accepted_import_result_sha256") == identity.get("import_result_sha256")
    )


def _mixed_high_order_transition_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        topologies = [str(value) for value in identity.get("element_topologies", [])]
        result_topologies = [
            str(value) for value in identity.get("result_element_topologies", [])
        ]
        orders = [int(value) for value in identity.get("polynomial_orders", [])]
        result_orders = [
            int(value) for value in identity.get("result_polynomial_orders", [])
        ]
        corners = [
            [int(value) for value in row]
            for row in identity.get("interface_corner_node_ids", [])
        ]
        result_corners = [
            [int(value) for value in row]
            for row in identity.get("result_interface_corner_node_ids", [])
        ]
        midnodes = [
            [int(value) for value in row]
            for row in identity.get("interface_midnode_ids", [])
        ]
        result_midnodes = [
            [int(value) for value in row]
            for row in identity.get("result_interface_midnode_ids", [])
        ]
        permutations = [
            [int(value) for value in row]
            for row in identity.get("face_node_permutations", [])
        ]
        result_permutations = [
            [int(value) for value in row]
            for row in identity.get("result_face_node_permutations", [])
        ]
        jacobians = [
            [float(value) for value in row]
            for row in identity.get("quadrature_scaled_jacobians", [])
        ]
        result_jacobians = [
            [float(value) for value in row]
            for row in identity.get("result_quadrature_scaled_jacobians", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transition_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "order_transition_generation",
                "midnode_transition_generation",
                "permutation_transition_generation",
                "parametric_transition_generation",
                "jacobian_transition_generation",
                "export_transition_generation",
                "result_transition_generation",
            )
        )
        and topologies == ["hex27", "pyramid14", "tet10"]
        and result_topologies == topologies
        and orders == [2, 2, 2]
        and result_orders == orders
        and len(corners) == len(midnodes) == len(permutations) == len(jacobians) >= 1
        and all(
            len(row) in {3, 4}
            and len(set(row)) == len(row)
            and all(value > 0 for value in row)
            for row in corners
        )
        and result_corners == corners
        and all(
            len(nodes) == len(corner)
            and len(set(nodes)) == len(nodes)
            and all(value > 0 for value in nodes)
            for corner, nodes in zip(corners, midnodes, strict=True)
        )
        and result_midnodes == midnodes
        and all(
            sorted(permutation) == list(range(len(corner) + len(nodes)))
            for corner, nodes, permutation in zip(
                corners, midnodes, permutations, strict=True
            )
        )
        and result_permutations == permutations
        and all(
            bool(row)
            and all(math.isfinite(value) and value > 0.0 for value in row)
            for row in jacobians
        )
        and result_jacobians == jacobians
        and _valid_sha256(identity.get("parametric_face_coordinate_sha256"))
        and identity.get("result_parametric_face_coordinate_sha256")
        == identity.get("parametric_face_coordinate_sha256")
        and _valid_sha256(identity.get("export_connectivity_sha256"))
        and identity.get("accepted_export_connectivity_sha256")
        == identity.get("export_connectivity_sha256")
    )


def _periodic_high_order_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        order = int(identity.get("polynomial_order"))
        result_order = int(identity.get("result_polynomial_order"))
        transform = [
            [float(value) for value in row]
            for row in identity.get("affine_transform_4x4", [])
        ]
        result_transform = [
            [float(value) for value in row]
            for row in identity.get("result_affine_transform_4x4", [])
        ]
        pair_fields = []
        for source, target in (
            ("corner_node_pairs", "result_corner_node_pairs"),
            ("edge_node_pairs", "result_edge_node_pairs"),
            ("face_node_pairs", "result_face_node_pairs"),
        ):
            pairs = [[int(value) for value in row] for row in identity.get(source, [])]
            result_pairs = [
                [int(value) for value in row] for row in identity.get(target, [])
            ]
            pair_fields.append((pairs, result_pairs))
        orientations = [int(value) for value in identity.get("orientation_signs", [])]
        result_orientations = [
            int(value) for value in identity.get("result_orientation_signs", [])
        ]
        source_sideset = int(identity.get("source_sideset_id"))
        result_source_sideset = int(identity.get("result_source_sideset_id"))
        target_sideset = int(identity.get("target_sideset_id"))
        result_target_sideset = int(identity.get("result_target_sideset_id"))
        residual = float(identity.get("maximum_periodic_residual_m"))
        result_residual = float(identity.get("result_maximum_periodic_residual_m"))
        tolerance = float(identity.get("periodic_residual_tolerance_m"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("periodic_generation") or "")
    all_pairs = [pair for pairs, _ in pair_fields for pair in pairs]
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "affine_periodic_generation",
                "corner_periodic_generation",
                "edge_periodic_generation",
                "face_periodic_generation",
                "orientation_periodic_generation",
                "sideset_periodic_generation",
                "residual_periodic_generation",
                "result_periodic_generation",
            )
        )
        and order >= 2
        and result_order == order
        and len(transform) == 4
        and all(len(row) == 4 and all(math.isfinite(value) for value in row) for row in transform)
        and transform[3] == [0.0, 0.0, 0.0, 1.0]
        and result_transform == transform
        and all(bool(pairs) and result_pairs == pairs for pairs, result_pairs in pair_fields)
        and all(len(pair) == 2 and pair[0] > 0 and pair[1] > 0 for pair in all_pairs)
        and len({pair[0] for pair in all_pairs}) == len(all_pairs)
        and len({pair[1] for pair in all_pairs}) == len(all_pairs)
        and bool(orientations)
        and all(value in {-1, 1} for value in orientations)
        and result_orientations == orientations
        and source_sideset > 0
        and target_sideset > 0
        and source_sideset != target_sideset
        and result_source_sideset == source_sideset
        and result_target_sideset == target_sideset
        and math.isfinite(residual)
        and residual >= 0.0
        and result_residual == residual
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and residual <= tolerance
        and _valid_sha256(identity.get("periodic_mesh_sha256"))
        and identity.get("result_periodic_mesh_sha256")
        == identity.get("periodic_mesh_sha256")
    )


def _hex_boundary_layer_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        first = float(identity.get("first_layer_thickness_m"))
        result_first = float(identity.get("result_first_layer_thickness_m"))
        ratio = float(identity.get("growth_ratio"))
        result_ratio = float(identity.get("result_growth_ratio"))
        count = int(identity.get("layer_count"))
        result_count = int(identity.get("result_layer_count"))
        thicknesses = [
            float(value) for value in identity.get("layer_thicknesses_m", [])
        ]
        result_thicknesses = [
            float(value) for value in identity.get("result_layer_thicknesses_m", [])
        ]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        block = int(identity.get("block_id"))
        result_block = int(identity.get("result_block_id"))
        sideset = int(identity.get("wall_sideset_id"))
        result_sideset = int(identity.get("result_wall_sideset_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("boundary_layer_generation") or "")
    expected = [first * ratio**index for index in range(count)]
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "thickness_generation",
                "growth_generation",
                "collision_generation",
                "jacobian_generation",
                "block_generation",
                "sideset_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and math.isfinite(first)
        and first > 0.0
        and result_first == first
        and identity.get("growth_law") == "geometric"
        and identity.get("result_growth_law") == identity.get("growth_law")
        and math.isfinite(ratio)
        and ratio >= 1.0
        and result_ratio == ratio
        and count >= 1
        and result_count == count
        and len(thicknesses) == count
        and all(math.isfinite(value) and value > 0.0 for value in thicknesses)
        and all(
            math.isclose(actual, target, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, target in zip(thicknesses, expected, strict=True)
        )
        and result_thicknesses == thicknesses
        and identity.get("collision_handling") in {
            "truncate_and_rebalance",
            "stop_before_collision",
        }
        and identity.get("result_collision_handling")
        == identity.get("collision_handling")
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and result_jacobian == jacobian
        and block > 0
        and result_block == block
        and sideset > 0
        and result_sideset == sideset
        and _valid_sha256(identity.get("boundary_layer_mesh_sha256"))
        and identity.get("accepted_boundary_layer_mesh_sha256")
        == identity.get("boundary_layer_mesh_sha256")
    )


def _pyramid_transition_closure_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        diagonal = [int(value) for value in identity.get("interface_diagonal", [])]
        result_diagonal = [
            int(value) for value in identity.get("result_interface_diagonal", [])
        ]
        orientations = [
            int(value) for value in identity.get("pyramid_orientation_signs", [])
        ]
        result_orientations = [
            int(value)
            for value in identity.get("result_pyramid_orientation_signs", [])
        ]
        faces = [
            [int(value) for value in row]
            for row in identity.get("shared_face_connectivity", [])
        ]
        result_faces = [
            [int(value) for value in row]
            for row in identity.get("result_shared_face_connectivity", [])
        ]
        owners = [str(value) for value in identity.get("region_owners", [])]
        result_owners = [
            str(value) for value in identity.get("result_region_owners", [])
        ]
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [
            float(value) for value in identity.get("result_scaled_jacobians", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transition_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "diagonal_generation",
                "orientation_generation",
                "face_generation",
                "region_generation",
                "jacobian_generation",
                "export_generation",
                "result_generation",
            )
        )
        and len(diagonal) == 2
        and diagonal == sorted(set(diagonal))
        and all(value > 0 for value in diagonal)
        and result_diagonal == diagonal
        and bool(orientations)
        and all(value == 1 for value in orientations)
        and result_orientations == orientations
        and bool(faces)
        and all(
            len(face) in {3, 4}
            and len(set(face)) == len(face)
            and all(value > 0 for value in face)
            for face in faces
        )
        and result_faces == faces
        and owners == ["hex_region", "transition_region", "tet_region"]
        and result_owners == owners
        and len(jacobians) == len(owners)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and _valid_sha256(identity.get("transition_export_sha256"))
        and identity.get("accepted_transition_export_sha256")
        == identity.get("transition_export_sha256")
    )


def _boolean_entity_lineage_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        pre_ids = [int(value) for value in identity.get("pre_entity_ids", [])]
        resolved_pre_ids = [
            int(value) for value in identity.get("resolved_pre_entity_ids", [])
        ]
        post_ids = [int(value) for value in identity.get("post_entity_ids", [])]
        resolved_post_ids = [
            int(value) for value in identity.get("resolved_post_entity_ids", [])
        ]
        lineage = [
            [int(value) for value in row]
            for row in identity.get("entity_lineage_pairs", [])
        ]
        resolved_lineage = [
            [int(value) for value in row]
            for row in identity.get("resolved_entity_lineage_pairs", [])
        ]
        orientations = [
            int(value) for value in identity.get("surface_orientation_signs", [])
        ]
        resolved_orientations = [
            int(value)
            for value in identity.get("resolved_surface_orientation_signs", [])
        ]
        measures = [float(value) for value in identity.get("surface_measures_m2", [])]
        resolved_measures = [
            float(value) for value in identity.get("resolved_surface_measures_m2", [])
        ]
        adjacency = [
            [int(value) for value in row]
            for row in identity.get("block_adjacency", [])
        ]
        resolved_adjacency = [
            [int(value) for value in row]
            for row in identity.get("resolved_block_adjacency", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("lineage_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "imprint_lineage_generation",
                "merge_lineage_generation",
                "orientation_lineage_generation",
                "measure_lineage_generation",
                "adjacency_lineage_generation",
                "journal_lineage_generation",
                "result_lineage_generation",
            )
        )
        and bool(pre_ids and post_ids)
        and pre_ids == sorted(set(pre_ids))
        and post_ids == sorted(set(post_ids))
        and all(value > 0 for value in pre_ids + post_ids)
        and resolved_pre_ids == pre_ids
        and resolved_post_ids == post_ids
        and bool(lineage)
        and all(
            len(pair) == 2 and pair[0] in pre_ids and pair[1] in post_ids
            for pair in lineage
        )
        and {pair[1] for pair in lineage} == set(post_ids)
        and resolved_lineage == lineage
        and len(orientations) == len(measures) == len(adjacency) == len(post_ids)
        and all(value in {-1, 1} for value in orientations)
        and resolved_orientations == orientations
        and all(math.isfinite(value) and value > 0.0 for value in measures)
        and resolved_measures == measures
        and all(
            len(pair) == 2 and pair[0] in post_ids and pair[1] > 0
            for pair in adjacency
        )
        and {pair[0] for pair in adjacency} == set(post_ids)
        and resolved_adjacency == adjacency
        and bool(str(identity.get("journal_generation_id") or ""))
        and identity.get("resolved_journal_generation_id")
        == identity.get("journal_generation_id")
        and _valid_sha256(identity.get("lineage_table_sha256"))
        and identity.get("resolved_lineage_table_sha256")
        == identity.get("lineage_table_sha256")
        and _valid_sha256(identity.get("model_digest_sha256"))
        and identity.get("accepted_model_digest_sha256")
        == identity.get("model_digest_sha256")
    )


def _exodus_high_order_restart_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        block_ids = [int(value) for value in identity.get("block_ids", [])]
        decoded_block_ids = [
            int(value) for value in identity.get("decoded_block_ids", [])
        ]
        topologies = [str(value) for value in identity.get("element_topologies", [])]
        decoded_topologies = [
            str(value) for value in identity.get("decoded_element_topologies", [])
        ]
        midnodes = [
            [int(value) for value in row]
            for row in identity.get("midnode_orderings", [])
        ]
        decoded_midnodes = [
            [int(value) for value in row]
            for row in identity.get("decoded_midnode_orderings", [])
        ]
        word_size = int(identity.get("integer_word_size_bits"))
        decoded_word_size = int(identity.get("decoded_integer_word_size_bits"))
        qa_records = [[str(value) for value in row] for row in identity.get("qa_records", [])]
        decoded_qa_records = [
            [str(value) for value in row] for row in identity.get("decoded_qa_records", [])
        ]
        restart = int(identity.get("restart_step_index"))
        decoded_restart = int(identity.get("decoded_restart_step_index"))
        owners = [
            [int(value) for value in row]
            for row in identity.get("sideset_owner_block_ids", [])
        ]
        decoded_owners = [
            [int(value) for value in row]
            for row in identity.get("decoded_sideset_owner_block_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("exodus_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "block_exodus_generation",
                "midnode_exodus_generation",
                "id_exodus_generation",
                "qa_exodus_generation",
                "restart_exodus_generation",
                "sideset_exodus_generation",
                "digest_exodus_generation",
                "result_exodus_generation",
            )
        )
        and bool(block_ids)
        and all(value > 2**31 - 1 for value in block_ids)
        and len(set(block_ids)) == len(block_ids)
        and decoded_block_ids == block_ids
        and len(topologies) == len(block_ids)
        and all(value in {"HEX27", "PYRAMID14", "TET10"} for value in topologies)
        and decoded_topologies == topologies
        and len(midnodes) == len(block_ids)
        and all(bool(row) and len(set(row)) == len(row) and min(row) > 0 for row in midnodes)
        and decoded_midnodes == midnodes
        and word_size == 64
        and decoded_word_size == word_size
        and bool(qa_records)
        and all(len(row) == 4 and all(row) for row in qa_records)
        and decoded_qa_records == qa_records
        and restart >= 0
        and decoded_restart == restart
        and bool(owners)
        and all(len(pair) == 2 and pair[0] > 0 and pair[1] in block_ids for pair in owners)
        and decoded_owners == owners
        and _valid_sha256(identity.get("exodus_file_sha256"))
        and identity.get("decoded_exodus_file_sha256")
        == identity.get("exodus_file_sha256")
        and _valid_sha256(identity.get("connectivity_table_sha256"))
        and identity.get("decoded_connectivity_table_sha256")
        == identity.get("connectivity_table_sha256")
    )


def _journal_transaction_restore_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        undo_count = int(identity.get("undo_count"))
        replayed_undo_count = int(identity.get("replayed_undo_count"))
        redo_count = int(identity.get("redo_count"))
        replayed_redo_count = int(identity.get("replayed_redo_count"))
        entity_ids = [int(value) for value in identity.get("allocated_entity_ids", [])]
        replayed_entity_ids = [
            int(value) for value in identity.get("replayed_allocated_entity_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transaction_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "undo_generation",
                "allocation_generation",
                "idempotency_generation",
                "save_generation",
                "restore_generation",
                "session_generation",
                "result_generation",
            )
        )
        and all(
            bool(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("transaction_id", "replayed_transaction_id"),
                ("save_point_id", "restored_save_point_id"),
                ("session_owner", "restored_session_owner"),
            )
        )
        and undo_count >= 0
        and replayed_undo_count == undo_count
        and redo_count >= 0
        and replayed_redo_count == redo_count
        and bool(entity_ids)
        and entity_ids == sorted(set(entity_ids))
        and all(value > 0 for value in entity_ids)
        and replayed_entity_ids == entity_ids
        and all(
            _valid_sha256(identity.get(source))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("idempotency_sha256", "replayed_idempotency_sha256"),
                ("saved_model_sha256", "restored_model_sha256"),
                ("journal_result_sha256", "accepted_journal_result_sha256"),
            )
        )
    )


def _exodus_assembly_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        assembly_ids = [int(value) for value in identity.get("assembly_ids", [])]
        decoded_assembly_ids = [
            int(value) for value in identity.get("decoded_assembly_ids", [])
        ]
        members = [
            [int(value) for value in row]
            for row in identity.get("assembly_members", [])
        ]
        decoded_members = [
            [int(value) for value in row]
            for row in identity.get("decoded_assembly_members", [])
        ]
        qa = [str(value) for value in identity.get("qa_record", [])]
        decoded_qa = [str(value) for value in identity.get("decoded_qa_record", [])]
        times = [float(value) for value in identity.get("time_values_s", [])]
        decoded_times = [
            float(value) for value in identity.get("decoded_time_values_s", [])
        ]
        block_owners = [
            [int(value) for value in row]
            for row in identity.get("block_owners", [])
        ]
        decoded_block_owners = [
            [int(value) for value in row]
            for row in identity.get("decoded_block_owners", [])
        ]
        sideset_owners = [
            [int(value) for value in row]
            for row in identity.get("sideset_owners", [])
        ]
        decoded_sideset_owners = [
            [int(value) for value in row]
            for row in identity.get("decoded_sideset_owners", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("exodus_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "assembly_generation",
                "frame_generation",
                "qa_generation",
                "time_generation",
                "block_generation",
                "sideset_generation",
                "file_generation",
                "result_generation",
            )
        )
        and bool(assembly_ids)
        and assembly_ids == sorted(set(assembly_ids))
        and all(value > 0 for value in assembly_ids)
        and decoded_assembly_ids == assembly_ids
        and len(members) == len(assembly_ids)
        and all(bool(row) and len(set(row)) == len(row) and min(row) > 0 for row in members)
        and decoded_members == members
        and identity.get("coordinate_frame") == "global_cartesian_m"
        and identity.get("decoded_coordinate_frame") == identity.get("coordinate_frame")
        and len(qa) == 4
        and all(qa)
        and decoded_qa == qa
        and len(times) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(right > left for left, right in zip(times, times[1:]))
        and decoded_times == times
        and bool(block_owners and sideset_owners)
        and all(
            len(pair) == 2 and pair[0] > 0 and pair[1] in assembly_ids
            for pair in block_owners + sideset_owners
        )
        and decoded_block_owners == block_owners
        and decoded_sideset_owners == sideset_owners
        and all(
            _valid_sha256(identity.get(source))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("exodus_file_sha256", "decoded_exodus_file_sha256"),
                ("assembly_table_sha256", "decoded_assembly_table_sha256"),
            )
        )
    )


def _hex_sheet_pillow_topology_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        cells = [int(value) for value in identity.get("cell_ids", [])]
        result_cells = [int(value) for value in identity.get("result_cell_ids", [])]
        incidence = [
            [int(value) for value in row]
            for row in identity.get("cell_face_incidence", [])
        ]
        result_incidence = [
            [int(value) for value in row]
            for row in identity.get("result_cell_face_incidence", [])
        ]
        counts = {
            name: int(identity.get(name))
            for name in ("vertex_count", "edge_count", "face_count", "cell_count")
        }
        result_counts = {
            name: int(identity.get(f"result_{name}"))
            for name in ("vertex_count", "edge_count", "face_count", "cell_count")
        }
        euler = int(identity.get("euler_characteristic"))
        result_euler = int(identity.get("result_euler_characteristic"))
        shell_count = int(identity.get("boundary_shell_count"))
        result_shell_count = int(identity.get("result_boundary_shell_count"))
        orientations = [
            int(value) for value in identity.get("cell_orientation_signs", [])
        ]
        result_orientations = [
            int(value)
            for value in identity.get("result_cell_orientation_signs", [])
        ]
        block_id = int(identity.get("block_id"))
        result_block_id = int(identity.get("result_block_id"))
        jacobians = [float(value) for value in identity.get("scaled_jacobians", [])]
        result_jacobians = [
            float(value) for value in identity.get("result_scaled_jacobians", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sheet_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "incidence_generation",
                "euler_generation",
                "shell_generation",
                "orientation_generation",
                "block_generation",
                "jacobian_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and identity.get("operation") == "pillow"
        and identity.get("result_operation") == identity.get("operation")
        and bool(cells)
        and cells == sorted(set(cells))
        and all(value > 0 for value in cells)
        and result_cells == cells
        and bool(incidence)
        and len({tuple(row) for row in incidence}) == len(incidence)
        and all(
            len(row) == 2 and row[0] in cells and row[1] > 0 for row in incidence
        )
        and result_incidence == incidence
        and all(value > 0 for value in counts.values())
        and result_counts == counts
        and counts["cell_count"] == len(cells)
        and counts["vertex_count"]
        - counts["edge_count"]
        + counts["face_count"]
        - counts["cell_count"]
        == euler
        and euler == 1
        and result_euler == euler
        and shell_count == 1
        and result_shell_count == shell_count
        and len(orientations) == len(cells)
        and all(value == 1 for value in orientations)
        and result_orientations == orientations
        and block_id > 0
        and result_block_id == block_id
        and len(jacobians) == len(cells)
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and result_jacobians == jacobians
        and _valid_sha256(identity.get("sheet_mesh_sha256"))
        and identity.get("accepted_sheet_mesh_sha256")
        == identity.get("sheet_mesh_sha256")
    )


def _multiblock_interface_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        tolerance = float(identity.get("merge_tolerance_m"))
        result_tolerance = float(identity.get("result_merge_tolerance_m"))
        node_pairs = [
            [int(value) for value in row]
            for row in identity.get("interface_node_pairs", [])
        ]
        result_node_pairs = [
            [int(value) for value in row]
            for row in identity.get("result_interface_node_pairs", [])
        ]
        face_pairs = [
            [int(value) for value in row]
            for row in identity.get("coincident_face_pairs", [])
        ]
        result_face_pairs = [
            [int(value) for value in row]
            for row in identity.get("result_coincident_face_pairs", [])
        ]
        face_owners = [
            [int(value) for value in row] for row in identity.get("face_owners", [])
        ]
        result_face_owners = [
            [int(value) for value in row]
            for row in identity.get("result_face_owners", [])
        ]
        block_ids = [int(value) for value in identity.get("block_ids", [])]
        result_block_ids = [
            int(value) for value in identity.get("result_block_ids", [])
        ]
        sideset_ids = [int(value) for value in identity.get("sideset_ids", [])]
        result_sideset_ids = [
            int(value) for value in identity.get("result_sideset_ids", [])
        ]
        duplicate_count = int(identity.get("duplicate_cell_count"))
        result_duplicate_count = int(identity.get("result_duplicate_cell_count"))
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("interface_generation") or "")
    face_ids = {value for pair in face_pairs for value in pair}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "merge_generation",
                "node_generation",
                "face_generation",
                "owner_generation",
                "block_generation",
                "sideset_generation",
                "duplicate_generation",
                "jacobian_generation",
                "export_generation",
                "result_generation",
            )
        )
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and bool(node_pairs)
        and all(len(pair) == 2 and pair[0] > 0 and pair[1] > 0 for pair in node_pairs)
        and len({pair[0] for pair in node_pairs}) == len(node_pairs)
        and len({pair[1] for pair in node_pairs}) == len(node_pairs)
        and result_node_pairs == node_pairs
        and bool(face_pairs)
        and all(
            len(pair) == 2
            and pair[0] > 0
            and pair[1] > 0
            and pair[0] != pair[1]
            for pair in face_pairs
        )
        and len(face_ids) == 2 * len(face_pairs)
        and result_face_pairs == face_pairs
        and len(face_owners) == len(face_ids)
        and all(
            len(row) == 2 and row[0] in face_ids and row[1] in block_ids
            for row in face_owners
        )
        and {row[0] for row in face_owners} == face_ids
        and result_face_owners == face_owners
        and bool(block_ids)
        and block_ids == sorted(set(block_ids))
        and all(value > 0 for value in block_ids)
        and result_block_ids == block_ids
        and bool(sideset_ids)
        and sideset_ids == sorted(set(sideset_ids))
        and all(value > 0 for value in sideset_ids)
        and result_sideset_ids == sideset_ids
        and duplicate_count == 0
        and result_duplicate_count == duplicate_count
        and math.isfinite(jacobian)
        and jacobian > 0.0
        and result_jacobian == jacobian
        and _valid_sha256(identity.get("multiblock_export_sha256"))
        and identity.get("accepted_multiblock_export_sha256")
        == identity.get("multiblock_export_sha256")
    )


def _command_failure_atomic_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        failed_index = int(identity.get("failed_command_index"))
        reported_index = int(identity.get("reported_failed_command_index"))
        error_code = int(identity.get("error_code"))
        reported_error_code = int(identity.get("reported_error_code"))
        next_entity_id = int(identity.get("next_entity_id"))
        rolled_back_next_entity_id = int(identity.get("rolled_back_next_entity_id"))
        undo_depth = int(identity.get("undo_depth"))
        rolled_back_undo_depth = int(identity.get("rolled_back_undo_depth"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("failure_generation") or "")
    category = str(identity.get("failure_category") or "")
    session_owner = str(identity.get("session_owner") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "command_generation",
                "error_generation",
                "rollback_generation",
                "allocator_generation",
                "undo_generation",
                "session_generation",
                "result_generation",
            )
        )
        and failed_index >= 0
        and reported_index == failed_index
        and error_code != 0
        and reported_error_code == error_code
        and bool(category)
        and category != "success"
        and identity.get("reported_failure_category") == category
        and _valid_sha256(identity.get("pre_transaction_model_sha256"))
        and identity.get("rolled_back_model_sha256")
        == identity.get("pre_transaction_model_sha256")
        and next_entity_id > 0
        and rolled_back_next_entity_id == next_entity_id
        and undo_depth >= 0
        and rolled_back_undo_depth == undo_depth
        and session_owner.startswith("headless-")
        and identity.get("rolled_back_session_owner") == session_owner
        and _valid_sha256(identity.get("failure_result_sha256"))
        and identity.get("accepted_failure_result_sha256")
        == identity.get("failure_result_sha256")
    )


def _cub_roundtrip_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        names = [
            [str(row[0]), int(row[1]), str(row[2])]
            for row in identity.get("entity_names", [])
        ]
        reopened_names = [
            [str(row[0]), int(row[1]), str(row[2])]
            for row in identity.get("reopened_entity_names", [])
        ]
        attributes = [
            [str(row[0]), int(row[1]), str(row[2]), str(row[3])]
            for row in identity.get("entity_attributes", [])
        ]
        reopened_attributes = [
            [str(row[0]), int(row[1]), str(row[2]), str(row[3])]
            for row in identity.get("reopened_entity_attributes", [])
        ]
        groups = [
            [str(row[0]), str(row[1]), int(row[2])]
            for row in identity.get("group_memberships", [])
        ]
        reopened_groups = [
            [str(row[0]), str(row[1]), int(row[2])]
            for row in identity.get("reopened_group_memberships", [])
        ]
        mesh = dict(identity.get("mesh_state") or {})
        reopened_mesh = dict(identity.get("reopened_mesh_state") or {})
        hex_count = int(mesh.get("hex"))
        quad_count = int(mesh.get("quad"))
        block_ids = [int(value) for value in mesh.get("block_ids", [])]
        sideset_ids = [int(value) for value in mesh.get("sideset_ids", [])]
    except (IndexError, TypeError, ValueError):
        return False
    generation = str(identity.get("roundtrip_generation") or "")
    kernel_version = str(identity.get("kernel_version") or "")
    entity_keys = {(row[0], row[1]) for row in names}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "kernel_generation",
                "entity_generation",
                "attribute_generation",
                "group_generation",
                "mesh_generation",
                "file_generation",
                "result_generation",
            )
        )
        and bool(kernel_version)
        and identity.get("reopened_kernel_version") == kernel_version
        and bool(names)
        and len(entity_keys) == len(names)
        and all(row[0] and row[1] > 0 and row[2] for row in names)
        and reopened_names == names
        and all(
            row[0] and row[1] > 0 and row[2] and row[3]
            and (row[0], row[1]) in entity_keys
            for row in attributes
        )
        and reopened_attributes == attributes
        and all(
            row[0] and row[1] and row[2] > 0 and (row[1], row[2]) in entity_keys
            for row in groups
        )
        and reopened_groups == groups
        and hex_count >= 0
        and quad_count >= 0
        and bool(block_ids)
        and block_ids == sorted(set(block_ids))
        and all(value > 0 for value in block_ids)
        and bool(sideset_ids)
        and sideset_ids == sorted(set(sideset_ids))
        and all(value > 0 for value in sideset_ids)
        and reopened_mesh == mesh
        and bool(str(identity.get("model_generation") or ""))
        and identity.get("reopened_model_generation")
        == identity.get("model_generation")
        and _valid_sha256(identity.get("cub_file_sha256"))
        and identity.get("reopened_cub_file_sha256")
        == identity.get("cub_file_sha256")
        and _valid_sha256(identity.get("roundtrip_result_sha256"))
        and identity.get("accepted_roundtrip_result_sha256")
        == identity.get("roundtrip_result_sha256")
    )


def _high_order_hex_family_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        element_ids = [int(value) for value in identity.get("element_ids", [])]
        families = [str(value) for value in identity.get("element_families", [])]
        node_counts = [int(value) for value in identity.get("node_counts", [])]
        corners = [int(value) for value in identity.get("corner_node_order", [])]
        edges = [int(value) for value in identity.get("edge_node_order", [])]
        face_nodes = [int(value) for value in identity.get("hex27_face_node_order", [])]
        center = int(identity.get("hex27_center_node"))
        jacobians = [float(value) for value in identity.get("minimum_jacobians", [])]
        quadrature_volumes = [
            float(value) for value in identity.get("quadrature_volumes_m3", [])
        ]
        geometric_volumes = [
            float(value) for value in identity.get("geometric_volumes_m3", [])
        ]
        curved_face_owners = [
            [int(value) for value in row]
            for row in identity.get("curved_face_owners", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("high_order_hex_generation") or "")
    all_hex27_roles = corners + edges + face_nodes + [center]
    result_pairs = (
        ("result_element_ids", element_ids),
        ("result_element_families", families),
        ("result_node_counts", node_counts),
        ("result_corner_node_order", corners),
        ("result_edge_node_order", edges),
        ("result_hex27_face_node_order", face_nodes),
        ("result_hex27_center_node", center),
        ("result_minimum_jacobians", jacobians),
        ("result_quadrature_volumes_m3", quadrature_volumes),
        ("result_geometric_volumes_m3", geometric_volumes),
        ("result_curved_face_owners", curved_face_owners),
    )
    try:
        results_match = all(
            (
                [float(value) for value in identity.get(key, [])] == expected
                if key in {
                    "result_minimum_jacobians",
                    "result_quadrature_volumes_m3",
                    "result_geometric_volumes_m3",
                }
                else identity.get(key) == expected
            )
            for key, expected in result_pairs
        )
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "family_generation",
                "node_generation",
                "reference_generation",
                "jacobian_generation",
                "volume_generation",
                "face_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and element_ids == [20, 27]
        and families == ["hex20_serendipity", "hex27_lagrange"]
        and node_counts == [20, 27]
        and corners == list(range(8))
        and edges == list(range(8, 20))
        and face_nodes == list(range(20, 26))
        and center == 26
        and len(set(all_hex27_roles)) == 27
        and set(all_hex27_roles) == set(range(27))
        and len(jacobians) == 2
        and all(math.isfinite(value) and value > 0.0 for value in jacobians)
        and len(quadrature_volumes) == len(geometric_volumes) == 2
        and all(value > 0.0 and math.isfinite(value) for value in quadrature_volumes)
        and all(value > 0.0 and math.isfinite(value) for value in geometric_volumes)
        and all(
            math.isclose(qv, gv, rel_tol=1.0e-9, abs_tol=1.0e-12)
            for qv, gv in zip(quadrature_volumes, geometric_volumes)
        )
        and len(curved_face_owners) == 2
        and all(len(row) == 2 and row[1] > 0 for row in curved_face_owners)
        and [row[0] for row in curved_face_owners] == element_ids
        and len({row[1] for row in curved_face_owners}) == 2
        and results_match
        and _valid_sha256(identity.get("high_order_mesh_sha256"))
        and identity.get("accepted_high_order_mesh_sha256")
        == identity.get("high_order_mesh_sha256")
    )


def _sheet_midplane_mass_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source_volume_id = int(identity.get("source_volume_id"))
        offset = float(identity.get("midplane_offset_m"))
        thickness = float(identity.get("thickness_m"))
        block_id = int(identity.get("shell_block_id"))
        top_sideset_id = int(identity.get("top_sideset_id"))
        bottom_sideset_id = int(identity.get("bottom_sideset_id"))
        area = float(identity.get("midplane_area_m2"))
        source_volume = float(identity.get("source_volume_m3"))
        density = float(identity.get("density_kg_m3"))
        mass = float(identity.get("shell_mass_kg"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sheet_generation") or "")
    mirrored_fields = (
        "source_volume_id",
        "midplane_offset_m",
        "thickness_m",
        "normal_orientation",
        "shell_block_id",
        "top_sideset_id",
        "bottom_sideset_id",
        "midplane_area_m2",
        "source_volume_m3",
        "density_kg_m3",
        "shell_mass_kg",
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_generation",
                "midplane_generation",
                "thickness_generation",
                "normal_generation",
                "block_generation",
                "sideset_generation",
                "mass_generation",
                "geometry_generation",
                "result_generation",
            )
        )
        and source_volume_id > 0
        and math.isfinite(offset)
        and math.isfinite(thickness)
        and thickness > 0.0
        and identity.get("normal_orientation") == "outward_source_volume"
        and block_id > 0
        and top_sideset_id > 0
        and bottom_sideset_id > 0
        and top_sideset_id != bottom_sideset_id
        and all(math.isfinite(value) and value > 0.0 for value in (area, source_volume, density, mass))
        and math.isclose(source_volume, area * thickness, rel_tol=1.0e-9, abs_tol=1.0e-12)
        and math.isclose(mass, area * thickness * density, rel_tol=1.0e-9, abs_tol=1.0e-12)
        and all(identity.get(f"result_{key}") == identity.get(key) for key in mirrored_fields)
        and _valid_sha256(identity.get("sheet_geometry_sha256"))
        and identity.get("accepted_sheet_geometry_sha256")
        == identity.get("sheet_geometry_sha256")
    )


def _sweep_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        volume_id = int(identity.get("volume_id"))
        source_id = int(identity.get("source_surface_id"))
        target_id = int(identity.get("target_surface_id"))
        intervals = int(identity.get("interval_count"))
        bias_ratio = float(identity.get("bias_ratio"))
        node_layers = int(identity.get("node_layer_count"))
        match_pairs = [
            [int(value) for value in row] for row in identity.get("match_pairs", [])
        ]
        periodic_layers = [
            [int(value) for value in row]
            for row in identity.get("periodic_layer_map", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    mirrored_fields = (
        "volume_id",
        "source_surface_id",
        "target_surface_id",
        "interval_count",
        "bias_ratio",
        "bias_direction",
        "match_pairs",
        "periodic_layer_map",
        "node_layer_count",
        "volume_scheme",
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_generation",
                "target_generation",
                "interval_generation",
                "bias_generation",
                "match_generation",
                "periodic_generation",
                "layer_generation",
                "journal_generation",
                "result_generation",
            )
        )
        and volume_id > 0
        and source_id > 0
        and target_id > 0
        and source_id != target_id
        and intervals > 0
        and math.isfinite(bias_ratio)
        and bias_ratio > 0.0
        and identity.get("bias_direction") in {"source_to_target", "target_to_source"}
        and node_layers == intervals + 1
        and bool(match_pairs)
        and all(len(row) == 2 and row[0] > 0 and row[1] > 0 for row in match_pairs)
        and len({tuple(row) for row in match_pairs}) == len(match_pairs)
        and bool(periodic_layers)
        and all(
            len(row) == 2
            and 0 <= row[0] <= intervals
            and 0 <= row[1] <= intervals
            and row[0] + row[1] == intervals
            for row in periodic_layers
        )
        and identity.get("volume_scheme") == "sweep"
        and all(
            identity.get(f"replayed_{key}") == identity.get(key)
            for key in mirrored_fields
        )
        and _valid_sha256(identity.get("journal_sha256"))
        and identity.get("replayed_journal_sha256") == identity.get("journal_sha256")
        and _valid_sha256(identity.get("sweep_result_sha256"))
        and identity.get("accepted_sweep_result_sha256")
        == identity.get("sweep_result_sha256")
    )


def _checkpoint_partition_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False

    def partition_map(value: object) -> dict[int, list[int]]:
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise ValueError("partition rows must be a sequence")
        result: dict[int, list[int]] = {}
        for row in value:
            if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 2:
                raise ValueError("partition row must contain an id and entity ids")
            partition = int(row[0])
            entity_ids = row[1]
            if isinstance(entity_ids, (str, bytes)) or not isinstance(entity_ids, Sequence):
                raise ValueError("partition entity ids must be a sequence")
            if partition in result:
                raise ValueError("duplicate partition id")
            result[partition] = [int(entity_id) for entity_id in entity_ids]
        return result

    def membership_rows(value: object) -> list[list[object]]:
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise ValueError("membership rows must be a sequence")
        result: list[list[object]] = []
        for row in value:
            if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 2:
                raise ValueError("membership row must contain a set id and entity ids")
            set_id = int(row[0])
            members = row[1]
            if isinstance(members, (str, bytes)) or not isinstance(members, Sequence):
                raise ValueError("membership entity ids must be a sequence")
            result.append([set_id, [int(entity_id) for entity_id in members]])
        return result

    try:
        partition_count = int(identity.get("partition_count"))
        owned = partition_map(identity.get("owned_entity_ids"))
        ghost = partition_map(identity.get("ghost_entity_ids"))
        persistent = [int(value) for value in identity.get("persistent_entity_ids", [])]
        blocks = membership_rows(identity.get("block_membership"))
        sidesets = membership_rows(identity.get("sideset_membership"))
        quality = [
            float(value)
            for value in identity.get("partition_minimum_scaled_jacobian", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("checkpoint_generation") or "")
    partitions = set(range(partition_count))
    flattened_owned = [entity_id for ids in owned.values() for entity_id in ids]
    persistent_set = set(persistent)
    mirrored_fields = (
        "partition_count",
        "owned_entity_ids",
        "ghost_entity_ids",
        "persistent_entity_ids",
        "block_membership",
        "sideset_membership",
        "partition_minimum_scaled_jacobian",
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "partition_generation",
                "owned_generation",
                "ghost_generation",
                "persistent_generation",
                "block_generation",
                "sideset_generation",
                "quality_generation",
                "model_generation",
                "result_generation",
            )
        )
        and partition_count > 1
        and set(owned) == partitions
        and set(ghost) == partitions
        and all(owned[partition] for partition in partitions)
        and all(entity_id > 0 for entity_id in flattened_owned)
        and len(flattened_owned) == len(set(flattened_owned))
        and persistent == sorted(flattened_owned)
        and all(
            all(entity_id in persistent_set and entity_id not in owned[partition] for entity_id in ghost[partition])
            for partition in partitions
        )
        and all(
            set_id > 0
            and bool(members)
            and all(entity_id in persistent_set for entity_id in members)
            for set_id, members in blocks + sidesets
        )
        and len(quality) == partition_count
        and all(math.isfinite(value) and value > 0.0 for value in quality)
        and all(
            identity.get(f"restored_{key}") == identity.get(key)
            for key in mirrored_fields
        )
        and _valid_sha256(identity.get("checkpoint_model_sha256"))
        and identity.get("restored_checkpoint_model_sha256")
        == identity.get("checkpoint_model_sha256")
        and _valid_sha256(identity.get("checkpoint_result_sha256"))
        and identity.get("accepted_checkpoint_result_sha256")
        == identity.get("checkpoint_result_sha256")
    )


def _anisotropic_hex_metric_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        vectors = [[float(value) for value in row] for row in identity.get("metric_eigenvectors", [])]
        result_vectors = [[float(value) for value in row] for row in identity.get("result_metric_eigenvectors", [])]
        sizes = [float(value) for value in identity.get("principal_sizes_m", [])]
        result_sizes = [float(value) for value in identity.get("result_principal_sizes_m", [])]
        gradation = float(identity.get("gradation_ratio"))
        result_gradation = float(identity.get("result_gradation_ratio"))
        maximum_gradation = float(identity.get("maximum_gradation_ratio"))
        result_maximum_gradation = float(identity.get("result_maximum_gradation_ratio"))
        alignment = float(identity.get("alignment_error_deg"))
        result_alignment = float(identity.get("result_alignment_error_deg"))
        maximum_alignment = float(identity.get("maximum_alignment_error_deg"))
        result_maximum_alignment = float(identity.get("result_maximum_alignment_error_deg"))
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        minimum_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_minimum_jacobian = float(identity.get("result_minimum_allowed_scaled_jacobian"))
        blocks = [int(value) for value in identity.get("block_ids", [])]
        result_blocks = [int(value) for value in identity.get("result_block_ids", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("anisotropic_generation") or "")
    orthonormal = (
        len(vectors) == 3
        and all(len(vector) == 3 for vector in vectors)
        and all(
            math.isclose(sum(value * value for value in vector), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for vector in vectors
        )
        and all(
            math.isclose(
                sum(left * right for left, right in zip(vectors[i], vectors[j], strict=True)),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            for i in range(3)
            for j in range(i + 1, 3)
        )
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "metric_generation",
                "direction_generation",
                "size_generation",
                "gradation_generation",
                "quality_generation",
                "block_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and orthonormal
        and result_vectors == vectors
        and len(sizes) == 3
        and all(math.isfinite(value) and value > 0.0 for value in sizes)
        and result_sizes == sizes
        and 1.0 <= gradation <= maximum_gradation
        and result_gradation == gradation
        and result_maximum_gradation == maximum_gradation
        and 0.0 <= alignment <= maximum_alignment
        and result_alignment == alignment
        and result_maximum_alignment == maximum_alignment
        and jacobian >= minimum_jacobian > 0.0
        and result_jacobian == jacobian
        and result_minimum_jacobian == minimum_jacobian
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and result_blocks == blocks
        and _valid_sha256(identity.get("anisotropic_mesh_sha256"))
        and identity.get("accepted_anisotropic_mesh_sha256")
        == identity.get("anisotropic_mesh_sha256")
    )


def _curved_highorder_boundary_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        normal_dot = float(identity.get("minimum_normal_dot"))
        result_normal_dot = float(identity.get("result_minimum_normal_dot"))
        hausdorff = float(identity.get("hausdorff_error_m"))
        result_hausdorff = float(identity.get("result_hausdorff_error_m"))
        maximum_hausdorff = float(identity.get("maximum_hausdorff_error_m"))
        result_maximum_hausdorff = float(identity.get("result_maximum_hausdorff_error_m"))
        cad_area = float(identity.get("cad_surface_area_m2"))
        mesh_area = float(identity.get("mesh_surface_area_m2"))
        area_tolerance = float(identity.get("surface_measure_tolerance"))
        cad_volume = float(identity.get("cad_volume_m3"))
        mesh_volume = float(identity.get("mesh_volume_m3"))
        volume_tolerance = float(identity.get("volume_measure_tolerance"))
        jacobian = float(identity.get("minimum_curved_jacobian"))
        result_jacobian = float(identity.get("result_minimum_curved_jacobian"))
        order = int(identity.get("polynomial_order"))
        result_order = int(identity.get("result_polynomial_order"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("curved_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "normal_generation",
                "hausdorff_generation",
                "area_generation",
                "volume_generation",
                "jacobian_generation",
                "order_generation",
                "geometry_generation",
                "result_generation",
            )
        )
        and 0.95 <= normal_dot <= 1.0
        and result_normal_dot == normal_dot
        and 0.0 <= hausdorff <= maximum_hausdorff
        and result_hausdorff == hausdorff
        and result_maximum_hausdorff == maximum_hausdorff
        and cad_area > 0.0
        and mesh_area > 0.0
        and 0.0 <= area_tolerance < 1.0
        and abs(mesh_area - cad_area) / cad_area <= area_tolerance
        and cad_volume > 0.0
        and mesh_volume > 0.0
        and 0.0 <= volume_tolerance < 1.0
        and abs(mesh_volume - cad_volume) / cad_volume <= volume_tolerance
        and jacobian > 0.0
        and result_jacobian == jacobian
        and order >= 2
        and result_order == order
        and bool(str(identity.get("geometry_owner") or ""))
        and identity.get("result_geometry_owner") == identity.get("geometry_owner")
        and _valid_sha256(identity.get("curved_geometry_sha256"))
        and identity.get("accepted_curved_geometry_sha256")
        == identity.get("curved_geometry_sha256")
    )


def _sideset_skin_remesh_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sideset = int(identity.get("sideset_id"))
        restored_sideset = int(identity.get("restored_sideset_id"))
        blocks = [int(value) for value in identity.get("adjacent_block_ids", [])]
        restored_blocks = [int(value) for value in identity.get("restored_adjacent_block_ids", [])]
        normals = [[float(value) for value in row] for row in identity.get("outward_normals", [])]
        restored_normals = [[float(value) for value in row] for row in identity.get("restored_outward_normals", [])]
        faces = [int(value) for value in identity.get("face_ids", [])]
        restored_faces = [int(value) for value in identity.get("restored_face_ids", [])]
        multiplicities = [int(value) for value in identity.get("face_multiplicities", [])]
        restored_multiplicities = [int(value) for value in identity.get("restored_face_multiplicities", [])]
        keys = [str(value) for value in identity.get("source_entity_keys", [])]
        restored_keys = [str(value) for value in identity.get("restored_source_entity_keys", [])]
        revision = int(identity.get("remesh_revision"))
        restored_revision = int(identity.get("restored_remesh_revision"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sideset_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "skin_generation",
                "remesh_generation",
                "adjacency_generation",
                "normal_generation",
                "face_generation",
                "entity_generation",
                "journal_generation",
                "result_generation",
            )
        )
        and sideset > 0
        and restored_sideset == sideset
        and len(blocks) == 1
        and blocks[0] > 0
        and restored_blocks == blocks
        and bool(faces)
        and len(set(faces)) == len(faces)
        and all(value > 0 for value in faces)
        and restored_faces == faces
        and len(normals) == len(multiplicities) == len(keys) == len(faces)
        and all(len(normal) == 3 for normal in normals)
        and all(
            math.isclose(sum(value * value for value in normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for normal in normals
        )
        and restored_normals == normals
        and multiplicities == [1] * len(faces)
        and restored_multiplicities == multiplicities
        and all(keys)
        and len(set(keys)) == len(keys)
        and restored_keys == keys
        and revision >= 0
        and restored_revision == revision
        and _valid_sha256(identity.get("sideset_journal_sha256"))
        and identity.get("replayed_sideset_journal_sha256")
        == identity.get("sideset_journal_sha256")
        and _valid_sha256(identity.get("sideset_result_sha256"))
        and identity.get("accepted_sideset_result_sha256")
        == identity.get("sideset_result_sha256")
    )


def _parallel_sculpt_determinism_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        seed = int(identity.get("partition_seed"))
        replay_seed = int(identity.get("replay_partition_seed"))
        ranks = int(identity.get("rank_count"))
        replay_ranks = int(identity.get("replay_rank_count"))
        owned = [int(value) for value in identity.get("owned_cell_counts", [])]
        replay_owned = [int(value) for value in identity.get("replay_owned_cell_counts", [])]
        ghost = [int(value) for value in identity.get("ghost_cell_counts", [])]
        replay_ghost = [int(value) for value in identity.get("replay_ghost_cell_counts", [])]
        stitched = int(identity.get("stitched_interface_pair_count"))
        replay_stitched = int(identity.get("replay_stitched_interface_pair_count"))
        qa = [str(value) for value in identity.get("qa_record", [])]
        replay_qa = [str(value) for value in identity.get("replay_qa_record", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sculpt_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "seed_generation",
                "rank_generation",
                "partition_generation",
                "stitch_generation",
                "qa_generation",
                "connectivity_generation",
                "invocation_generation",
                "result_generation",
            )
        )
        and seed >= 0
        and replay_seed == seed
        and ranks > 1
        and replay_ranks == ranks
        and len(owned) == len(ghost) == ranks
        and all(value > 0 for value in owned)
        and all(value >= 0 for value in ghost)
        and replay_owned == owned
        and replay_ghost == ghost
        and stitched > 0
        and replay_stitched == stitched
        and len(qa) >= 3
        and all(qa)
        and replay_qa == qa
        and str(identity.get("invocation_owner") or "").startswith("headless:")
        and identity.get("replay_invocation_owner") == identity.get("invocation_owner")
        and _valid_sha256(identity.get("connectivity_sha256"))
        and identity.get("replay_connectivity_sha256") == identity.get("connectivity_sha256")
        and _valid_sha256(identity.get("sculpt_export_sha256"))
        and identity.get("accepted_sculpt_export_sha256")
        == identity.get("sculpt_export_sha256")
    )


def _hex_sweep_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source_nodes = int(identity.get("source_node_count"))
        result_source_nodes = int(identity.get("result_source_node_count"))
        target_nodes = int(identity.get("target_node_count"))
        result_target_nodes = int(identity.get("result_target_node_count"))
        layers = int(identity.get("layer_count"))
        result_layers = int(identity.get("result_layer_count"))
        bias = float(identity.get("interval_bias"))
        result_bias = float(identity.get("result_interval_bias"))
        hex_count = int(identity.get("hex_count"))
        result_hex_count = int(identity.get("result_hex_count"))
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        minimum_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_minimum_jacobian = float(identity.get("result_minimum_allowed_scaled_jacobian"))
        cad_volume = float(identity.get("cad_volume_m3"))
        mesh_volume = float(identity.get("mesh_volume_m3"))
        result_mesh_volume = float(identity.get("result_mesh_volume_m3"))
        tolerance = float(identity.get("volume_tolerance_m3"))
        result_tolerance = float(identity.get("result_volume_tolerance_m3"))
    except (TypeError, ValueError):
        return False
    owners = identity.get("boundary_owners")
    result_owners = identity.get("result_boundary_owners")
    generation = str(identity.get("sweep_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_generation",
                "target_generation",
                "layer_generation",
                "orientation_generation",
                "quality_generation",
                "boundary_generation",
                "volume_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("source_surface_topology") or ""))
        and identity.get("result_source_surface_topology")
        == identity.get("source_surface_topology")
        and identity.get("target_surface_topology")
        == identity.get("source_surface_topology")
        and identity.get("result_target_surface_topology")
        == identity.get("target_surface_topology")
        and source_nodes == target_nodes > 0
        and result_source_nodes == source_nodes
        and result_target_nodes == target_nodes
        and layers > 0
        and result_layers == layers
        and math.isfinite(bias)
        and bias > 0.0
        and result_bias == bias
        and hex_count > 0
        and result_hex_count == hex_count
        and math.isfinite(jacobian)
        and jacobian >= minimum_jacobian > 0.0
        and result_jacobian == jacobian
        and result_minimum_jacobian == minimum_jacobian
        and identity.get("orientation") == "source_to_target_positive"
        and identity.get("result_orientation") == identity.get("orientation")
        and isinstance(owners, list)
        and len(owners) >= 3
        and all(isinstance(value, str) and value for value in owners)
        and len(set(owners)) == len(owners)
        and result_owners == owners
        and all(math.isfinite(value) and value > 0.0 for value in (cad_volume, mesh_volume))
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and abs(mesh_volume - cad_volume) <= tolerance
        and result_mesh_volume == mesh_volume
        and _valid_sha256(identity.get("sweep_mesh_sha256"))
        and identity.get("accepted_sweep_mesh_sha256")
        == identity.get("sweep_mesh_sha256")
    )


def _mixed_transition_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        hex_count = int(identity.get("hex_count"))
        result_hex_count = int(identity.get("result_hex_count"))
        tet_count = int(identity.get("tet_count"))
        result_tet_count = int(identity.get("result_tet_count"))
        pyramid_count = int(identity.get("pyramid_count"))
        result_pyramid_count = int(identity.get("result_pyramid_count"))
        faces = [int(value) for value in identity.get("interface_face_ids", [])]
        result_faces = [int(value) for value in identity.get("result_interface_face_ids", [])]
        nodes = [int(value) for value in identity.get("shared_interface_node_ids", [])]
        result_nodes = [int(value) for value in identity.get("result_shared_interface_node_ids", [])]
        volumes = [float(value) for value in identity.get("signed_region_volumes_m3", [])]
        result_volumes = [float(value) for value in identity.get("result_signed_region_volumes_m3", [])]
        cad_volume = float(identity.get("cad_volume_m3"))
        mesh_volume = float(identity.get("mesh_volume_m3"))
        result_mesh_volume = float(identity.get("result_mesh_volume_m3"))
        tolerance = float(identity.get("volume_tolerance_m3"))
        result_tolerance = float(identity.get("result_volume_tolerance_m3"))
    except (TypeError, ValueError):
        return False
    labels = identity.get("region_labels")
    result_labels = identity.get("result_region_labels")
    generation = str(identity.get("transition_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "hex_generation",
                "tet_generation",
                "pyramid_generation",
                "interface_generation",
                "orientation_generation",
                "region_generation",
                "volume_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and hex_count > 0
        and tet_count > 0
        and pyramid_count > 0
        and result_hex_count == hex_count
        and result_tet_count == tet_count
        and result_pyramid_count == pyramid_count
        and len(faces) == pyramid_count
        and all(value > 0 for value in faces)
        and len(set(faces)) == len(faces)
        and result_faces == faces
        and bool(nodes)
        and all(value > 0 for value in nodes)
        and len(set(nodes)) == len(nodes)
        and result_nodes == nodes
        and identity.get("interface_face_orientation")
        == "hex_outward_pyramid_inward"
        and identity.get("result_interface_face_orientation")
        == identity.get("interface_face_orientation")
        and isinstance(labels, list)
        and len(labels) == len(volumes) >= 3
        and all(isinstance(value, str) and value for value in labels)
        and len(set(labels)) == len(labels)
        and result_labels == labels
        and all(math.isfinite(value) and value > 0.0 for value in volumes)
        and result_volumes == volumes
        and all(math.isfinite(value) and value > 0.0 for value in (cad_volume, mesh_volume))
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and abs(sum(volumes) - cad_volume) <= tolerance
        and abs(mesh_volume - cad_volume) <= tolerance
        and result_mesh_volume == mesh_volume
        and bool(str(identity.get("mixed_mesh_owner") or ""))
        and identity.get("result_mixed_mesh_owner") == identity.get("mixed_mesh_owner")
        and _valid_sha256(identity.get("mixed_mesh_sha256"))
        and identity.get("accepted_mixed_mesh_sha256")
        == identity.get("mixed_mesh_sha256")
    )


def _journal_replay_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        command_count = int(identity.get("command_count"))
        replay_command_count = int(identity.get("replay_command_count"))
        undo_groups = [[int(value) for value in row] for row in identity.get("undo_groups", [])]
        replay_undo_groups = [[int(value) for value in row] for row in identity.get("replay_undo_groups", [])]
        entities = [int(value) for value in identity.get("entity_ids_after_first_run", [])]
        replay_entities = [int(value) for value in identity.get("entity_ids_after_replay", [])]
        checkpoint = int(identity.get("checkpoint_generation_id"))
        replay_checkpoint = int(identity.get("replay_checkpoint_generation_id"))
    except (TypeError, ValueError):
        return False
    reset_sequence = identity.get("reset_sequence")
    replay_reset_sequence = identity.get("replay_reset_sequence")
    covered = [
        index
        for group in undo_groups
        if len(group) == 2 and group[0] <= group[1]
        for index in range(group[0], group[1] + 1)
    ]
    generation = str(identity.get("journal_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "undo_generation",
                "idempotence_generation",
                "entity_generation",
                "reset_generation",
                "checkpoint_generation",
                "replay_generation",
                "database_generation",
                "result_generation",
            )
        )
        and command_count > 0
        and replay_command_count == command_count
        and bool(undo_groups)
        and all(len(group) == 2 for group in undo_groups)
        and covered == list(range(command_count))
        and replay_undo_groups == undo_groups
        and bool(entities)
        and all(value > 0 for value in entities)
        and len(set(entities)) == len(entities)
        and replay_entities == entities
        and isinstance(reset_sequence, list)
        and len(reset_sequence) >= 2
        and reset_sequence[0] == "reset"
        and reset_sequence[-1] == "export"
        and replay_reset_sequence == reset_sequence
        and checkpoint >= 0
        and replay_checkpoint == checkpoint
        and all(
            _valid_sha256(identity.get(source))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("command_order_sha256", "replay_command_order_sha256"),
                ("first_database_sha256", "replay_database_sha256"),
                ("journal_result_sha256", "accepted_journal_result_sha256"),
            )
        )
        and str(identity.get("journal_owner") or "").startswith("headless:")
        and identity.get("replay_journal_owner") == identity.get("journal_owner")
    )


def _exodus_semantic_export_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        blocks = [int(value) for value in identity.get("block_ids", [])]
        exported_blocks = [int(value) for value in identity.get("exported_block_ids", [])]
        sidesets = [int(value) for value in identity.get("sideset_ids", [])]
        exported_sidesets = [int(value) for value in identity.get("exported_sideset_ids", [])]
        orientations = [int(value) for value in identity.get("sideset_orientation", [])]
        exported_orientations = [int(value) for value in identity.get("exported_sideset_orientation", [])]
        nodesets = [int(value) for value in identity.get("nodeset_ids", [])]
        exported_nodesets = [int(value) for value in identity.get("exported_nodeset_ids", [])]
        maximum_id = int(identity.get("maximum_entity_id"))
        exported_maximum_id = int(identity.get("exported_maximum_entity_id"))
        element_map = [int(value) for value in identity.get("element_map", [])]
        exported_element_map = [int(value) for value in identity.get("exported_element_map", [])]
        mesh_generation = int(identity.get("mesh_generation_id"))
        exported_mesh_generation = int(identity.get("exported_mesh_generation_id"))
    except (TypeError, ValueError):
        return False
    topologies = identity.get("block_topologies")
    exported_topologies = identity.get("exported_block_topologies")
    generation = str(identity.get("exodus_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "block_generation",
                "sideset_generation",
                "nodeset_generation",
                "integer_generation",
                "topology_generation",
                "map_generation",
                "owner_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(blocks)
        and all(value > 0 for value in blocks)
        and len(set(blocks)) == len(blocks)
        and exported_blocks == blocks
        and isinstance(topologies, list)
        and len(topologies) == len(blocks)
        and all(value in {"HEX8", "PYRAMID5", "TET4"} for value in topologies)
        and exported_topologies == topologies
        and bool(sidesets)
        and all(value > 0 for value in sidesets)
        and len(set(sidesets)) == len(sidesets)
        and exported_sidesets == sidesets
        and len(orientations) == len(sidesets)
        and all(value in {-1, 1} for value in orientations)
        and exported_orientations == orientations
        and bool(nodesets)
        and all(value > 0 for value in nodesets)
        and len(set(nodesets)) == len(nodesets)
        and exported_nodesets == nodesets
        and identity.get("int64_ids") is True
        and identity.get("exported_int64_ids") is True
        and maximum_id > 2_147_483_647
        and exported_maximum_id == maximum_id
        and bool(element_map)
        and all(value > 0 for value in element_map)
        and len(set(element_map)) == len(element_map)
        and exported_element_map == element_map
        and bool(str(identity.get("assembly_owner") or ""))
        and identity.get("exported_assembly_owner") == identity.get("assembly_owner")
        and mesh_generation >= 0
        and exported_mesh_generation == mesh_generation
        and all(
            _valid_sha256(identity.get(source))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("exodus_mesh_sha256", "exported_exodus_mesh_sha256"),
                ("exodus_file_sha256", "accepted_exodus_file_sha256"),
            )
        )
    )


def _periodic_hex_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        transform = [
            [float(value) for value in row]
            for row in identity.get("rigid_transform", [])
        ]
        result_transform = [
            [float(value) for value in row]
            for row in identity.get("result_rigid_transform", [])
        ]
        pairs = [
            [int(value) for value in row]
            for row in identity.get("periodic_node_pairs", [])
        ]
        result_pairs = [
            [int(value) for value in row]
            for row in identity.get("result_periodic_node_pairs", [])
        ]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        minimum_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_minimum_jacobian = float(
            identity.get("result_minimum_allowed_scaled_jacobian")
        )
    except (TypeError, ValueError):
        return False
    rotation = [row[:3] for row in transform[:3]]
    determinant = (
        rotation[0][0]
        * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1]
        * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2]
        * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    ) if len(rotation) == 3 and all(len(row) == 3 for row in rotation) else math.nan
    rotation_is_orthonormal = len(rotation) == 3 and all(
        math.isclose(
            sum(rotation[row][axis] * rotation[column][axis] for axis in range(3)),
            1.0 if row == column else 0.0,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        for row in range(3)
        for column in range(3)
    )
    generation = str(identity.get("periodic_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "pair_generation",
                "transform_generation",
                "node_order_generation",
                "orientation_generation",
                "jacobian_generation",
                "region_generation",
                "export_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("source_surface") or ""))
        and identity.get("result_source_surface") == identity.get("source_surface")
        and bool(str(identity.get("target_surface") or ""))
        and identity.get("target_surface") != identity.get("source_surface")
        and identity.get("result_target_surface") == identity.get("target_surface")
        and len(transform) == 4
        and all(len(row) == 4 for row in transform)
        and all(math.isfinite(value) for row in transform for value in row)
        and transform[3] == [0.0, 0.0, 0.0, 1.0]
        and rotation_is_orthonormal
        and math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_transform == transform
        and len(pairs) >= 4
        and all(len(pair) == 2 and all(value > 0 for value in pair) for pair in pairs)
        and len({pair[0] for pair in pairs}) == len(pairs)
        and len({pair[1] for pair in pairs}) == len(pairs)
        and result_pairs == pairs
        and identity.get("face_orientation") == "source_outward_target_inward"
        and identity.get("result_face_orientation") == identity.get("face_orientation")
        and math.isfinite(jacobian)
        and math.isfinite(minimum_jacobian)
        and jacobian >= minimum_jacobian > 0.0
        and result_jacobian == jacobian
        and result_minimum_jacobian == minimum_jacobian
        and bool(str(identity.get("region_owner") or ""))
        and identity.get("result_region_owner") == identity.get("region_owner")
        and str(identity.get("export_owner") or "").startswith("headless:")
        and identity.get("result_export_owner") == identity.get("export_owner")
        and _valid_sha256(identity.get("periodic_mesh_sha256"))
        and identity.get("accepted_periodic_mesh_sha256")
        == identity.get("periodic_mesh_sha256")
        and _valid_sha256(identity.get("periodic_result_sha256"))
        and identity.get("accepted_periodic_result_sha256")
        == identity.get("periodic_result_sha256")
    )


def _curved_hex_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        order = int(identity.get("polynomial_order"))
        result_order = int(identity.get("result_polynomial_order"))
        projection = float(identity.get("maximum_cad_projection_distance_m"))
        result_projection = float(
            identity.get("result_maximum_cad_projection_distance_m")
        )
        allowed_projection = float(identity.get("allowed_cad_projection_distance_m"))
        result_allowed_projection = float(
            identity.get("result_allowed_cad_projection_distance_m")
        )
        cad_volume = float(identity.get("cad_volume_m3"))
        volume = float(identity.get("curved_mesh_volume_m3"))
        result_volume = float(identity.get("result_curved_mesh_volume_m3"))
        tolerance = float(identity.get("volume_tolerance_m3"))
        result_tolerance = float(identity.get("result_volume_tolerance_m3"))
        jacobian = float(identity.get("minimum_high_order_jacobian"))
        result_jacobian = float(identity.get("result_minimum_high_order_jacobian"))
        minimum_jacobian = float(
            identity.get("minimum_allowed_high_order_jacobian")
        )
        result_minimum_jacobian = float(
            identity.get("result_minimum_allowed_high_order_jacobian")
        )
    except (TypeError, ValueError):
        return False
    edge_roles = identity.get("edge_midnode_roles")
    result_edge_roles = identity.get("result_edge_midnode_roles")
    face_roles = identity.get("face_midnode_roles")
    result_face_roles = identity.get("result_face_midnode_roles")
    generation = str(identity.get("curve_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "edge_generation",
                "face_generation",
                "projection_generation",
                "volume_generation",
                "jacobian_generation",
                "order_generation",
                "cad_generation",
                "mesh_generation",
                "export_generation",
                "result_generation",
            )
        )
        and identity.get("element_family") == "HEX27"
        and identity.get("result_element_family") == identity.get("element_family")
        and order == 2
        and result_order == order
        and isinstance(edge_roles, list)
        and edge_roles == ["edge_mid"] * 12
        and result_edge_roles == edge_roles
        and isinstance(face_roles, list)
        and face_roles == ["face_mid"] * 6
        and result_face_roles == face_roles
        and all(
            math.isfinite(value)
            for value in (
                projection,
                result_projection,
                allowed_projection,
                result_allowed_projection,
                cad_volume,
                volume,
                result_volume,
                tolerance,
                result_tolerance,
                jacobian,
                result_jacobian,
                minimum_jacobian,
                result_minimum_jacobian,
            )
        )
        and 0.0 <= projection <= allowed_projection
        and allowed_projection > 0.0
        and result_projection == projection
        and result_allowed_projection == allowed_projection
        and cad_volume > 0.0
        and volume > 0.0
        and tolerance > 0.0
        and abs(volume - cad_volume) <= tolerance
        and result_volume == volume
        and result_tolerance == tolerance
        and jacobian >= minimum_jacobian > 0.0
        and result_jacobian == jacobian
        and result_minimum_jacobian == minimum_jacobian
        and bool(str(identity.get("cad_owner") or ""))
        and identity.get("result_cad_owner") == identity.get("cad_owner")
        and str(identity.get("mesh_owner") or "").startswith("headless:")
        and identity.get("result_mesh_owner") == identity.get("mesh_owner")
        and _valid_sha256(identity.get("curved_export_sha256"))
        and identity.get("accepted_curved_export_sha256")
        == identity.get("curved_export_sha256")
    )


def _imprint_merge_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        tolerance = float(identity.get("merge_tolerance_m"))
        replay_tolerance = float(identity.get("replay_merge_tolerance_m"))
        count = int(identity.get("coincident_topology_count"))
        replay_count = int(identity.get("replay_coincident_topology_count"))
        lineage = [
            [int(value) for value in row]
            for row in identity.get("merged_entity_lineage", [])
        ]
        replay_lineage = [
            [int(value) for value in row]
            for row in identity.get("replay_merged_entity_lineage", [])
        ]
        model_generation = int(identity.get("model_generation_id"))
        replay_model_generation = int(identity.get("replay_model_generation_id"))
    except (TypeError, ValueError):
        return False
    commands = identity.get("command_sequence")
    replay_commands = identity.get("replay_command_sequence")
    topology = identity.get("final_topology_counts")
    replay_topology = identity.get("replay_final_topology_counts")
    generation = str(identity.get("imprint_merge_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "tolerance_generation",
                "topology_generation",
                "lineage_generation",
                "command_generation",
                "checkpoint_generation",
                "model_generation",
                "final_generation",
                "database_generation",
                "result_generation",
            )
        )
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and replay_tolerance == tolerance
        and count > 0
        and replay_count == count
        and len(lineage) == count
        and all(len(pair) == 2 and all(value > 0 for value in pair) for pair in lineage)
        and len({pair[0] for pair in lineage}) == count
        and len({pair[1] for pair in lineage}) == count
        and replay_lineage == lineage
        and commands == ["imprint all", "merge all", "compress all"]
        and replay_commands == commands
        and str(identity.get("checkpoint_owner") or "").startswith("headless:")
        and identity.get("replay_checkpoint_owner") == identity.get("checkpoint_owner")
        and model_generation >= 0
        and replay_model_generation == model_generation
        and isinstance(topology, Mapping)
        and set(topology) == {"volume", "surface", "curve"}
        and all(isinstance(value, int) and value > 0 for value in topology.values())
        and replay_topology == topology
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("imprint_result_sha256"))
        and identity.get("accepted_imprint_result_sha256")
        == identity.get("imprint_result_sha256")
    )


def _headless_batch_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        exit_code = int(identity.get("process_exit_code"))
        result_exit_code = int(identity.get("result_process_exit_code"))
        save_generation = int(identity.get("database_save_generation_id"))
        result_save_generation = int(
            identity.get("result_database_save_generation_id")
        )
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("batch_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "process_generation",
                "license_generation",
                "journal_generation",
                "log_generation",
                "database_generation",
                "command_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and identity.get("execution_mode") == "nographics_batch"
        and identity.get("result_execution_mode") == identity.get("execution_mode")
        and identity.get("gui_launched") is False
        and identity.get("result_gui_launched") is False
        and exit_code == 0
        and result_exit_code == exit_code
        and identity.get("result_license_warning") == identity.get("license_warning")
        and identity.get("license_fallback") == "limited_mode_batch_completed"
        and identity.get("result_license_fallback") == identity.get("license_fallback")
        and identity.get("journal_completion_marker") == "CAEAI_BATCH_COMPLETE"
        and identity.get("result_journal_completion_marker")
        == identity.get("journal_completion_marker")
        and str(identity.get("log_owner") or "").startswith("headless:")
        and identity.get("result_log_owner") == identity.get("log_owner")
        and save_generation >= 0
        and result_save_generation == save_generation
        and _valid_sha256(identity.get("command_sha256"))
        and identity.get("result_command_sha256") == identity.get("command_sha256")
        and str(identity.get("process_owner") or "").startswith("coreform_cubit:")
        and "-nographics" in str(identity.get("process_owner") or "")
        and "-batch" in str(identity.get("process_owner") or "")
        and identity.get("result_process_owner") == identity.get("process_owner")
        and _valid_sha256(identity.get("batch_result_sha256"))
        and identity.get("accepted_batch_result_sha256")
        == identity.get("batch_result_sha256")
    )


def _midsurface_shell_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        face_pairs = [
            [int(value) for value in pair] for pair in identity.get("paired_face_ids", [])
        ]
        result_face_pairs = [
            [int(value) for value in pair] for pair in identity.get("result_paired_face_ids", [])
        ]
        distances = [float(value) for value in identity.get("paired_face_distance_m", [])]
        result_distances = [
            float(value) for value in identity.get("result_paired_face_distance_m", [])
        ]
        thicknesses = [float(value) for value in identity.get("shell_thickness_m", [])]
        result_thicknesses = [
            float(value) for value in identity.get("result_shell_thickness_m", [])
        ]
        normal_dots = [float(value) for value in identity.get("paired_normal_dot", [])]
        result_normal_dots = [
            float(value) for value in identity.get("result_paired_normal_dot", [])
        ]
        areas = [float(value) for value in identity.get("midsurface_area_m2", [])]
        result_areas = [float(value) for value in identity.get("result_midsurface_area_m2", [])]
        reconstructed_volume = float(identity.get("reconstructed_volume_m3"))
        result_reconstructed_volume = float(identity.get("result_reconstructed_volume_m3"))
        source_volume = float(identity.get("source_volume_m3"))
        volume_tolerance = float(identity.get("volume_tolerance_m3"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("midsurface_generation") or "")
    flat_faces = [face_id for pair in face_pairs for face_id in pair]
    reconstructed_from_shells = sum(area * thickness for area, thickness in zip(areas, thicknesses))
    shell_sidesets = list(identity.get("shell_sidesets") or [])
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "facepair_generation",
                "thickness_generation",
                "normal_generation",
                "area_generation",
                "volume_generation",
                "block_generation",
                "geometry_generation",
                "export_generation",
                "result_generation",
            )
        )
        and bool(face_pairs)
        and all(
            len(pair) == 2 and pair[0] > 0 and pair[1] > 0 and pair[0] != pair[1]
            for pair in face_pairs
        )
        and len(set(flat_faces)) == len(flat_faces)
        and result_face_pairs == face_pairs
        and len(distances) == len(thicknesses) == len(normal_dots) == len(areas) == len(face_pairs)
        and all(math.isfinite(value) and value > 0.0 for value in distances + thicknesses + areas)
        and all(
            math.isclose(distance, thickness, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for distance, thickness in zip(distances, thicknesses)
        )
        and result_distances == distances
        and result_thicknesses == thicknesses
        and all(
            math.isfinite(value) and math.isclose(value, -1.0, rel_tol=0.0, abs_tol=1.0e-12)
            for value in normal_dots
        )
        and result_normal_dots == normal_dots
        and result_areas == areas
        and math.isfinite(volume_tolerance)
        and volume_tolerance > 0.0
        and math.isclose(
            reconstructed_volume,
            reconstructed_from_shells,
            rel_tol=1.0e-12,
            abs_tol=volume_tolerance,
        )
        and math.isclose(
            result_reconstructed_volume,
            reconstructed_volume,
            rel_tol=1.0e-12,
            abs_tol=volume_tolerance,
        )
        and math.isclose(
            source_volume, reconstructed_volume, rel_tol=1.0e-12, abs_tol=volume_tolerance
        )
        and bool(str(identity.get("shell_block") or ""))
        and identity.get("result_shell_block") == identity.get("shell_block")
        and bool(shell_sidesets)
        and len(set(shell_sidesets)) == len(shell_sidesets)
        and list(identity.get("result_shell_sidesets") or []) == shell_sidesets
        and bool(str(identity.get("geometry_owner") or ""))
        and identity.get("result_geometry_owner") == identity.get("geometry_owner")
        and _valid_sha256(identity.get("midsurface_export_sha256"))
        and identity.get("accepted_midsurface_export_sha256")
        == identity.get("midsurface_export_sha256")
        and _valid_sha256(identity.get("midsurface_result_sha256"))
        and identity.get("accepted_midsurface_result_sha256")
        == identity.get("midsurface_result_sha256")
    )


def _cohesive_crack_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        face_pairs = [
            [int(value) for value in pair] for pair in identity.get("duplicated_face_pairs", [])
        ]
        result_face_pairs = [
            [int(value) for value in pair]
            for pair in identity.get("result_duplicated_face_pairs", [])
        ]
        node_pairs = [
            [int(value) for value in pair] for pair in identity.get("duplicated_node_pairs", [])
        ]
        result_node_pairs = [
            [int(value) for value in pair]
            for pair in identity.get("result_duplicated_node_pairs", [])
        ]
        crack_front = [int(value) for value in identity.get("crack_front_nodes", [])]
        result_crack_front = [int(value) for value in identity.get("result_crack_front_nodes", [])]
        minimum_jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_minimum_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        allowed_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed_jacobian = float(identity.get("result_minimum_allowed_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("cohesive_generation") or "")
    source_faces = [pair[0] for pair in face_pairs if len(pair) == 2]
    target_faces = [pair[1] for pair in face_pairs if len(pair) == 2]
    source_nodes = [pair[0] for pair in node_pairs if len(pair) == 2]
    target_nodes = [pair[1] for pair in node_pairs if len(pair) == 2]
    cohesive_block = str(identity.get("cohesive_block") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "face_generation",
                "nodepair_generation",
                "front_generation",
                "normal_generation",
                "orientation_generation",
                "traction_generation",
                "block_generation",
                "jacobian_generation",
                "mesh_generation",
                "export_generation",
                "result_generation",
            )
        )
        and bool(face_pairs)
        and all(
            len(pair) == 2 and pair[0] > 0 and pair[1] > 0 and pair[0] != pair[1]
            for pair in face_pairs
        )
        and len(set(source_faces)) == len(source_faces)
        and len(set(target_faces)) == len(target_faces)
        and set(source_faces).isdisjoint(target_faces)
        and result_face_pairs == face_pairs
        and bool(node_pairs)
        and all(
            len(pair) == 2 and pair[0] > 0 and pair[1] > 0 and pair[0] != pair[1]
            for pair in node_pairs
        )
        and len(set(source_nodes)) == len(source_nodes)
        and len(set(target_nodes)) == len(target_nodes)
        and set(source_nodes).isdisjoint(target_nodes)
        and result_node_pairs == node_pairs
        and bool(crack_front)
        and len(set(crack_front)) == len(crack_front)
        and set(crack_front).issubset(source_nodes)
        and result_crack_front == crack_front
        and identity.get("interface_normal_relation") == "opposed_outward_normals"
        and identity.get("result_interface_normal_relation")
        == identity.get("interface_normal_relation")
        and identity.get("cohesive_orientation") == "positive_reference_orientation"
        and identity.get("result_cohesive_orientation") == identity.get("cohesive_orientation")
        and identity.get("traction_direction") == "source_to_target"
        and identity.get("result_traction_direction") == identity.get("traction_direction")
        and cohesive_block.startswith("block:")
        and "cohesive" in cohesive_block.lower()
        and identity.get("result_cohesive_block") == cohesive_block
        and math.isfinite(minimum_jacobian)
        and math.isfinite(allowed_jacobian)
        and allowed_jacobian > 0.0
        and minimum_jacobian >= allowed_jacobian
        and result_minimum_jacobian == minimum_jacobian
        and result_allowed_jacobian == allowed_jacobian
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("result_mesh_owner") == identity.get("mesh_owner")
        and _valid_sha256(identity.get("cohesive_export_sha256"))
        and identity.get("accepted_cohesive_export_sha256")
        == identity.get("cohesive_export_sha256")
        and _valid_sha256(identity.get("cohesive_result_sha256"))
        and identity.get("accepted_cohesive_result_sha256")
        == identity.get("cohesive_result_sha256")
    )


def _virtual_geometry_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        curves = [int(value) for value in identity.get("suppressed_curve_ids", [])]
        replay_curves = [int(value) for value in identity.get("replay_suppressed_curve_ids", [])]
        surfaces = [int(value) for value in identity.get("suppressed_surface_ids", [])]
        replay_surfaces = [
            int(value) for value in identity.get("replay_suppressed_surface_ids", [])
        ]
        topology_map = {
            str(key): int(value)
            for key, value in _mapping(
                identity.get("virtual_topology_map"), "virtual_topology_map"
            ).items()
        }
        replay_topology_map = {
            str(key): int(value)
            for key, value in _mapping(
                identity.get("replay_virtual_topology_map"), "replay_virtual_topology_map"
            ).items()
        }
        quality_before = float(identity.get("minimum_quality_before"))
        replay_quality_before = float(identity.get("replay_minimum_quality_before"))
        quality_after = float(identity.get("minimum_quality_after"))
        replay_quality_after = float(identity.get("replay_minimum_quality_after"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("virtual_generation") or "")
    suppressed_keys = {str(value) for value in curves + surfaces}
    block_inheritance = identity.get("block_inheritance")
    sideset_inheritance = identity.get("sideset_inheritance")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "suppression_generation",
                "topology_generation",
                "inheritance_generation",
                "quality_generation",
                "undo_generation",
                "checkpoint_generation",
                "database_generation",
                "result_generation",
            )
        )
        and bool(curves or surfaces)
        and all(value > 0 for value in curves + surfaces)
        and len(set(curves)) == len(curves)
        and len(set(surfaces)) == len(surfaces)
        and replay_curves == curves
        and replay_surfaces == surfaces
        and set(topology_map) == suppressed_keys
        and all(value > 0 for value in topology_map.values())
        and len(set(topology_map.values())) == len(topology_map)
        and replay_topology_map == topology_map
        and isinstance(block_inheritance, Mapping)
        and bool(block_inheritance)
        and all(str(value).startswith("block:") for value in block_inheritance.values())
        and identity.get("replay_block_inheritance") == block_inheritance
        and isinstance(sideset_inheritance, Mapping)
        and bool(sideset_inheritance)
        and all(str(value).startswith("sideset:") for value in sideset_inheritance.values())
        and identity.get("replay_sideset_inheritance") == sideset_inheritance
        and math.isfinite(quality_before)
        and math.isfinite(quality_after)
        and quality_before >= 0.0
        and quality_after > quality_before
        and replay_quality_before == quality_before
        and replay_quality_after == quality_after
        and identity.get("undo_restored_topology") is True
        and identity.get("replay_undo_restored_topology") is True
        and str(identity.get("checkpoint_owner") or "").startswith("headless:")
        and identity.get("replay_checkpoint_owner") == identity.get("checkpoint_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("virtual_result_sha256"))
        and identity.get("accepted_virtual_result_sha256") == identity.get("virtual_result_sha256")
    )


def _anisotropic_crack_contract_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        eigenvalues = [float(value) for value in identity.get("metric_eigenvalues", [])]
        result_eigenvalues = [
            float(value) for value in identity.get("result_metric_eigenvalues", [])
        ]
        eigenvectors = [
            [float(value) for value in row] for row in identity.get("metric_eigenvectors", [])
        ]
        result_eigenvectors = [
            [float(value) for value in row]
            for row in identity.get("result_metric_eigenvectors", [])
        ]
        sizes = [float(value) for value in identity.get("sizing_field_m", [])]
        result_sizes = [float(value) for value in identity.get("result_sizing_field_m", [])]
        tangent = [float(value) for value in identity.get("crack_front_tangent", [])]
        result_tangent = [float(value) for value in identity.get("result_crack_front_tangent", [])]
        aligned_axis = int(identity.get("aligned_metric_axis"))
        result_aligned_axis = int(identity.get("result_aligned_metric_axis"))
        transition_faces = int(identity.get("cohesive_transition_face_count"))
        result_transition_faces = int(identity.get("result_cohesive_transition_face_count"))
        minimum_quality = float(identity.get("minimum_quality"))
        result_minimum_quality = float(identity.get("result_minimum_quality"))
        allowed_quality = float(identity.get("minimum_allowed_quality"))
        result_allowed_quality = float(identity.get("result_minimum_allowed_quality"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("anisotropic_generation") or "")
    frame_ok = (
        len(eigenvectors) == 3
        and all(
            len(row) == 3 and all(math.isfinite(value) for value in row) for row in eigenvectors
        )
        and all(
            math.isclose(sum(value * value for value in row), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
            for row in eigenvectors
        )
        and all(
            math.isclose(
                sum(eigenvectors[i][k] * eigenvectors[j][k] for k in range(3)),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            for i in range(3)
            for j in range(i + 1, 3)
        )
    )
    metric_size_scales = (
        [size * math.sqrt(value) for size, value in zip(sizes, eigenvalues)]
        if len(sizes) == len(eigenvalues)
        else []
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "metric_generation",
                "eigenframe_generation",
                "sizing_generation",
                "alignment_generation",
                "cohesive_generation",
                "conformity_generation",
                "quality_generation",
                "region_generation",
                "export_generation",
                "result_generation",
            )
        )
        and len(eigenvalues) == 3
        and all(math.isfinite(value) and value > 0.0 for value in eigenvalues)
        and eigenvalues == sorted(eigenvalues)
        and len(set(eigenvalues)) == 3
        and result_eigenvalues == eigenvalues
        and frame_ok
        and result_eigenvectors == eigenvectors
        and len(sizes) == 3
        and all(math.isfinite(value) and value > 0.0 for value in sizes)
        and sizes[0] > sizes[1] > sizes[2]
        and result_sizes == sizes
        and len(metric_size_scales) == 3
        and all(
            math.isclose(value, metric_size_scales[0], rel_tol=1.0e-12, abs_tol=1.0e-15)
            for value in metric_size_scales[1:]
        )
        and len(tangent) == 3
        and all(math.isfinite(value) for value in tangent)
        and math.isclose(sum(value * value for value in tangent), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_tangent == tangent
        and 0 <= aligned_axis < 3
        and result_aligned_axis == aligned_axis
        and math.isclose(
            abs(sum(tangent[index] * eigenvectors[aligned_axis][index] for index in range(3))),
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and transition_faces > 0
        and result_transition_faces == transition_faces
        and identity.get("hex_interface_conformal") is True
        and identity.get("result_hex_interface_conformal") is True
        and math.isfinite(minimum_quality)
        and math.isfinite(allowed_quality)
        and allowed_quality > 0.0
        and minimum_quality >= allowed_quality
        and result_minimum_quality == minimum_quality
        and result_allowed_quality == allowed_quality
        and bool(str(identity.get("region_owner") or ""))
        and identity.get("result_region_owner") == identity.get("region_owner")
        and _valid_sha256(identity.get("anisotropic_export_sha256"))
        and identity.get("accepted_anisotropic_export_sha256")
        == identity.get("anisotropic_export_sha256")
        and _valid_sha256(identity.get("anisotropic_result_sha256"))
        and identity.get("accepted_anisotropic_result_sha256")
        == identity.get("anisotropic_result_sha256")
    )


def _periodic_hex_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        transform = [[float(value) for value in row] for row in identity.get("rigid_transform", [])]
        result_transform = [
            [float(value) for value in row] for row in identity.get("result_rigid_transform", [])
        ]
        pairs = [[int(value) for value in row] for row in identity.get("periodic_node_pairs", [])]
        result_pairs = [
            [int(value) for value in row] for row in identity.get("result_periodic_node_pairs", [])
        ]
        source_edges = [
            [int(value) for value in row] for row in identity.get("source_edge_node_order", [])
        ]
        result_source_edges = [
            [int(value) for value in row]
            for row in identity.get("result_source_edge_node_order", [])
        ]
        target_edges = [
            [int(value) for value in row] for row in identity.get("target_edge_node_order", [])
        ]
        result_target_edges = [
            [int(value) for value in row]
            for row in identity.get("result_target_edge_node_order", [])
        ]
        jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        minimum_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_minimum_jacobian = float(identity.get("result_minimum_allowed_scaled_jacobian"))
        database_generation = int(identity.get("database_generation_id"))
        result_database_generation = int(identity.get("result_database_generation_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("periodic_generation") or "")
    rotation = [row[:3] for row in transform[:3]]
    rotation_ok = len(rotation) == 3 and all(len(row) == 3 for row in rotation) and all(
        math.isclose(
            sum(rotation[row][axis] * rotation[column][axis] for axis in range(3)),
            1.0 if row == column else 0.0,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        for row in range(3)
        for column in range(3)
    )
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
        if rotation_ok
        else math.nan
    )
    source_nodes = [pair[0] for pair in pairs if len(pair) == 2]
    target_nodes = [pair[1] for pair in pairs if len(pair) == 2]
    node_map = dict(zip(source_nodes, target_nodes))
    sidesets = list(identity.get("sideset_owners") or [])
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "pair_generation",
                "transform_generation",
                "orientation_generation",
                "edge_order_generation",
                "jacobian_generation",
                "set_generation",
                "database_generation",
                "export_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("source_surface") or ""))
        and identity.get("result_source_surface") == identity.get("source_surface")
        and bool(str(identity.get("target_surface") or ""))
        and identity.get("target_surface") != identity.get("source_surface")
        and identity.get("result_target_surface") == identity.get("target_surface")
        and len(transform) == 4
        and all(len(row) == 4 and all(math.isfinite(value) for value in row) for row in transform)
        and transform[3] == [0.0, 0.0, 0.0, 1.0]
        and rotation_ok
        and math.isclose(determinant, 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_transform == transform
        and len(pairs) >= 4
        and all(len(pair) == 2 and pair[0] > 0 and pair[1] > 0 for pair in pairs)
        and len(set(source_nodes)) == len(source_nodes)
        and len(set(target_nodes)) == len(target_nodes)
        and result_pairs == pairs
        and len(source_edges) == len(target_edges) >= 4
        and all(len(edge) == 2 for edge in source_edges + target_edges)
        and result_source_edges == source_edges
        and result_target_edges == target_edges
        and all(
            source_edge[0] in node_map
            and source_edge[1] in node_map
            and target_edge == [node_map[source_edge[0]], node_map[source_edge[1]]]
            for source_edge, target_edge in zip(source_edges, target_edges)
        )
        and identity.get("face_orientation") == "source_outward_target_inward"
        and identity.get("result_face_orientation") == identity.get("face_orientation")
        and math.isfinite(jacobian)
        and jacobian >= minimum_jacobian > 0.0
        and result_jacobian == jacobian
        and result_minimum_jacobian == minimum_jacobian
        and str(identity.get("block_owner") or "").startswith("block:")
        and identity.get("result_block_owner") == identity.get("block_owner")
        and len(sidesets) == 2
        and len(set(sidesets)) == 2
        and all(str(value).startswith("sideset:") for value in sidesets)
        and list(identity.get("result_sideset_owners") or []) == sidesets
        and database_generation > 0
        and result_database_generation == database_generation
        and _valid_sha256(identity.get("periodic_export_sha256"))
        and identity.get("accepted_periodic_export_sha256") == identity.get("periodic_export_sha256")
    )


def _thin_sweep_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source = int(identity.get("source_surface_id"))
        result_source = int(identity.get("result_source_surface_id"))
        target = int(identity.get("target_surface_id"))
        result_target = int(identity.get("result_target_surface_id"))
        source_intervals = int(identity.get("source_interval_count"))
        result_source_intervals = int(identity.get("result_source_interval_count"))
        target_intervals = int(identity.get("target_interval_count"))
        result_target_intervals = int(identity.get("result_target_interval_count"))
        side_intervals = [int(value) for value in identity.get("propagated_side_intervals", [])]
        result_side_intervals = [
            int(value) for value in identity.get("result_propagated_side_intervals", [])
        ]
        layers = int(identity.get("layer_count"))
        result_layers = int(identity.get("result_layer_count"))
        total_thickness = float(identity.get("total_thickness_m"))
        result_total_thickness = float(identity.get("result_total_thickness_m"))
        thicknesses = [float(value) for value in identity.get("layer_thickness_m", [])]
        result_thicknesses = [
            float(value) for value in identity.get("result_layer_thickness_m", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("thin_sweep_generation") or "")
    topology = list(identity.get("side_topology") or [])
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_generation",
                "target_generation",
                "interval_generation",
                "layer_generation",
                "thickness_generation",
                "topology_generation",
                "orientation_generation",
                "volume_generation",
                "export_generation",
                "result_generation",
            )
        )
        and source > 0
        and target > 0
        and source != target
        and result_source == source
        and result_target == target
        and source_intervals == target_intervals > 0
        and result_source_intervals == source_intervals
        and result_target_intervals == target_intervals
        and len(side_intervals) == 4
        and all(value > 0 for value in side_intervals)
        and sum(side_intervals) == source_intervals
        and result_side_intervals == side_intervals
        and layers > 0
        and result_layers == layers
        and len(thicknesses) == layers
        and all(math.isfinite(value) and value > 0.0 for value in thicknesses)
        and math.isfinite(total_thickness)
        and total_thickness > 0.0
        and math.isclose(sum(thicknesses), total_thickness, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_total_thickness == total_thickness
        and result_thicknesses == thicknesses
        and len(topology) == 4
        and all(value == "structured_quad" for value in topology)
        and list(identity.get("result_side_topology") or []) == topology
        and identity.get("element_orientation") == "source_to_target_positive"
        and identity.get("result_element_orientation") == identity.get("element_orientation")
        and str(identity.get("volume_owner") or "").startswith("volume:")
        and identity.get("result_volume_owner") == identity.get("volume_owner")
        and _valid_sha256(identity.get("thin_sweep_export_sha256"))
        and identity.get("accepted_thin_sweep_export_sha256")
        == identity.get("thin_sweep_export_sha256")
    )


def _journal_recreate_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        deleted = [int(value) for value in identity.get("deleted_entity_ids", [])]
        replay_deleted = [int(value) for value in identity.get("replay_deleted_entity_ids", [])]
        recreated = [int(value) for value in identity.get("recreated_entity_ids", [])]
        replay_recreated = [int(value) for value in identity.get("replay_recreated_entity_ids", [])]
        groups = _mapping(identity.get("group_membership"), "group_membership")
        replay_groups = _mapping(identity.get("replay_group_membership"), "replay_group_membership")
        blocks = _mapping(identity.get("block_membership"), "block_membership")
        replay_blocks = _mapping(identity.get("replay_block_membership"), "replay_block_membership")
        sidesets = _mapping(identity.get("sideset_membership"), "sideset_membership")
        replay_sidesets = _mapping(identity.get("replay_sideset_membership"), "replay_sideset_membership")
        undo_depth = int(identity.get("undo_depth"))
        replay_undo_depth = int(identity.get("replay_undo_depth"))
        checkpoint = int(identity.get("checkpoint_generation_id"))
        replay_checkpoint = int(identity.get("replay_checkpoint_generation_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("journal_generation") or "")
    recreated_set = set(recreated)
    grouped_entities = {int(value) for values in groups.values() for value in values}
    blocked_entities = {int(value) for values in blocks.values() for value in values}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "delete_generation",
                "recreate_generation",
                "group_generation",
                "block_generation",
                "sideset_generation",
                "undo_generation",
                "checkpoint_generation",
                "database_generation",
                "result_generation",
            )
        )
        and bool(deleted)
        and len(deleted) == len(recreated)
        and all(value > 0 for value in deleted + recreated)
        and len(set(deleted)) == len(deleted)
        and len(recreated_set) == len(recreated)
        and set(deleted).isdisjoint(recreated_set)
        and replay_deleted == deleted
        and replay_recreated == recreated
        and bool(groups)
        and all(str(key).startswith("group:") and bool(value) for key, value in groups.items())
        and grouped_entities == recreated_set
        and replay_groups == groups
        and bool(blocks)
        and all(str(key).startswith("block:") and bool(value) for key, value in blocks.items())
        and blocked_entities == recreated_set
        and replay_blocks == blocks
        and bool(sidesets)
        and all(
            str(key).startswith("sideset:")
            and bool(value)
            and all(int(entity) > 0 for entity in value)
            for key, value in sidesets.items()
        )
        and replay_sidesets == sidesets
        and undo_depth > 0
        and replay_undo_depth == undo_depth
        and checkpoint > 0
        and replay_checkpoint == checkpoint
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("journal_result_sha256"))
        and identity.get("accepted_journal_result_sha256") == identity.get("journal_result_sha256")
    )


def _adaptive_size_field_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sizes = [float(value) for value in identity.get("size_field_samples_m", [])]
        result_sizes = [float(value) for value in identity.get("result_size_field_samples_m", [])]
        nodes = [int(value) for value in identity.get("projected_node_ids", [])]
        result_nodes = [int(value) for value in identity.get("result_projected_node_ids", [])]
        distances = [float(value) for value in identity.get("projection_distances_m", [])]
        result_distances = [
            float(value) for value in identity.get("result_projection_distances_m", [])
        ]
        maximum_distance = float(identity.get("maximum_projection_distance_m"))
        result_maximum_distance = float(identity.get("result_maximum_projection_distance_m"))
        boundaries = [int(value) for value in identity.get("region_boundary_ids", [])]
        result_boundaries = [int(value) for value in identity.get("result_region_boundary_ids", [])]
        refinement = int(identity.get("refinement_generation_id"))
        result_refinement = int(identity.get("result_refinement_generation_id"))
        histogram = [int(value) for value in identity.get("quality_histogram", [])]
        result_histogram = [int(value) for value in identity.get("result_quality_histogram", [])]
        block_map = _mapping(identity.get("block_map"), "block_map")
        result_block_map = _mapping(identity.get("result_block_map"), "result_block_map")
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("adaptive_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "field_generation",
                "projection_generation",
                "region_generation",
                "refinement_generation",
                "quality_generation",
                "block_generation",
                "session_generation",
                "export_generation",
                "result_generation",
            )
        )
        and len(sizes) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in sizes)
        and all(sizes[index] > sizes[index + 1] for index in range(len(sizes) - 1))
        and result_sizes == sizes
        and len(nodes) == len(distances) >= 3
        and all(value > 0 for value in nodes)
        and len(set(nodes)) == len(nodes)
        and result_nodes == nodes
        and math.isfinite(maximum_distance)
        and maximum_distance > 0.0
        and all(math.isfinite(value) and 0.0 <= value <= maximum_distance for value in distances)
        and result_distances == distances
        and result_maximum_distance == maximum_distance
        and len(boundaries) >= 3
        and all(value > 0 for value in boundaries)
        and len(set(boundaries)) == len(boundaries)
        and result_boundaries == boundaries
        and refinement > 0
        and result_refinement == refinement
        and len(histogram) >= 3
        and all(value >= 0 for value in histogram)
        and sum(histogram) > 0
        and result_histogram == histogram
        and bool(block_map)
        and all(str(key).startswith("volume:") for key in block_map)
        and all(str(value).startswith("block:") for value in block_map.values())
        and result_block_map == block_map
        and str(identity.get("session_owner") or "").startswith("headless:")
        and identity.get("result_session_owner") == identity.get("session_owner")
        and _valid_sha256(identity.get("adaptive_export_sha256"))
        and identity.get("accepted_adaptive_export_sha256") == identity.get("adaptive_export_sha256")
    )


def _medial_axis_hex_decomposition_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sheet_pairs = [
            [int(value) for value in pair]
            for pair in identity.get("paired_sheet_ids", [])
        ]
        result_sheet_pairs = [
            [int(value) for value in pair]
            for pair in identity.get("result_paired_sheet_ids", [])
        ]
        thicknesses = [float(value) for value in identity.get("local_thickness_m", [])]
        result_thicknesses = [
            float(value) for value in identity.get("result_local_thickness_m", [])
        ]
        cells = _mapping(identity.get("decomposition_cells"), "decomposition_cells")
        result_cells = _mapping(
            identity.get("result_decomposition_cells"), "result_decomposition_cells"
        )
        shared_faces = [int(value) for value in identity.get("shared_topology_faces", [])]
        result_shared_faces = [
            int(value) for value in identity.get("result_shared_topology_faces", [])
        ]
        intervals = {
            str(key): int(value)
            for key, value in _mapping(identity.get("interval_counts"), "interval_counts").items()
        }
        result_intervals = {
            str(key): int(value)
            for key, value in _mapping(
                identity.get("result_interval_counts"), "result_interval_counts"
            ).items()
        }
        minimum_jacobian = float(identity.get("minimum_scaled_jacobian"))
        result_minimum_jacobian = float(identity.get("result_minimum_scaled_jacobian"))
        allowed_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed_jacobian = float(
            identity.get("result_minimum_allowed_scaled_jacobian")
        )
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("medial_axis_generation") or "")
    flattened_sheets = [value for pair in sheet_pairs for value in pair]
    normalized_cells = {
        str(key): [int(value) for value in values] for key, values in cells.items()
    }
    normalized_result_cells = {
        str(key): [int(value) for value in values] for key, values in result_cells.items()
    }
    shared_face_set = set(shared_faces)
    face_incidence = {
        face: sum(face in values for values in normalized_cells.values())
        for face in shared_face_set
    }
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "sheet_generation",
                "thickness_generation",
                "decomposition_generation",
                "topology_generation",
                "interval_generation",
                "quality_generation",
                "block_generation",
                "export_generation",
                "result_generation",
            )
        )
        and bool(sheet_pairs)
        and all(len(pair) == 2 and pair[0] != pair[1] and min(pair) > 0 for pair in sheet_pairs)
        and len(flattened_sheets) == len(set(flattened_sheets))
        and result_sheet_pairs == sheet_pairs
        and len(thicknesses) == len(sheet_pairs)
        and all(math.isfinite(value) and value > 0.0 for value in thicknesses)
        and result_thicknesses == thicknesses
        and len(normalized_cells) >= 2
        and all(
            str(key).startswith("cell:")
            and len(values) >= 4
            and len(values) == len(set(values))
            and all(value > 0 for value in values)
            for key, values in normalized_cells.items()
        )
        and normalized_result_cells == normalized_cells
        and bool(shared_face_set)
        and len(shared_face_set) == len(shared_faces)
        and result_shared_faces == shared_faces
        and all(count == 2 for count in face_incidence.values())
        and bool(intervals)
        and all(key.startswith("curve:") and value > 0 and value % 2 == 0 for key, value in intervals.items())
        and result_intervals == intervals
        and identity.get("interval_parity") == "compatible_even"
        and identity.get("result_interval_parity") == identity.get("interval_parity")
        and math.isfinite(minimum_jacobian)
        and math.isfinite(allowed_jacobian)
        and minimum_jacobian >= allowed_jacobian > 0.0
        and result_minimum_jacobian == minimum_jacobian
        and result_allowed_jacobian == allowed_jacobian
        and str(identity.get("block_owner") or "").startswith("block:")
        and identity.get("result_block_owner") == identity.get("block_owner")
        and _valid_sha256(identity.get("medial_axis_export_sha256"))
        and identity.get("accepted_medial_axis_export_sha256")
        == identity.get("medial_axis_export_sha256")
    )


def _curve_chain_boundary_layer_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        chain = [int(value) for value in identity.get("curve_chain_order", [])]
        result_chain = [int(value) for value in identity.get("result_curve_chain_order", [])]
        intervals = [int(value) for value in identity.get("interval_counts", [])]
        result_intervals = [int(value) for value in identity.get("result_interval_counts", [])]
        corner_sums = [int(value) for value in identity.get("corner_interval_sums", [])]
        result_corner_sums = [
            int(value) for value in identity.get("result_corner_interval_sums", [])
        ]
        thicknesses = [
            float(value) for value in identity.get("boundary_layer_thickness_m", [])
        ]
        result_thicknesses = [
            float(value) for value in identity.get("result_boundary_layer_thickness_m", [])
        ]
        total_thickness = float(identity.get("total_boundary_layer_thickness_m"))
        result_total_thickness = float(
            identity.get("result_total_boundary_layer_thickness_m")
        )
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("curve_chain_generation") or "")
    biases = list(identity.get("bias_directions") or [])
    result_biases = list(identity.get("result_bias_directions") or [])
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "chain_generation",
                "interval_generation",
                "bias_generation",
                "corner_generation",
                "boundary_layer_generation",
                "orientation_generation",
                "sideset_generation",
                "export_generation",
                "result_generation",
            )
        )
        and len(chain) >= 4
        and all(value > 0 for value in chain)
        and len(chain) == len(set(chain))
        and result_chain == chain
        and identity.get("chain_orientation") == "counterclockwise"
        and identity.get("result_chain_orientation") == identity.get("chain_orientation")
        and len(intervals) == len(chain)
        and all(value > 0 for value in intervals)
        and result_intervals == intervals
        and len(biases) == len(chain)
        and all(value in {"forward", "reverse"} for value in biases)
        and result_biases == biases
        and len(corner_sums) == len(chain)
        and all(value > 0 and value % 2 == 0 for value in corner_sums)
        and result_corner_sums == corner_sums
        and len(thicknesses) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in thicknesses)
        and all(thicknesses[index] <= thicknesses[index + 1] for index in range(len(thicknesses) - 1))
        and math.isfinite(total_thickness)
        and math.isclose(sum(thicknesses), total_thickness, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_thicknesses == thicknesses
        and result_total_thickness == total_thickness
        and identity.get("element_orientation") == "outward_positive"
        and identity.get("result_element_orientation") == identity.get("element_orientation")
        and str(identity.get("sideset_owner") or "").startswith("sideset:")
        and identity.get("result_sideset_owner") == identity.get("sideset_owner")
        and _valid_sha256(identity.get("curve_chain_export_sha256"))
        and identity.get("accepted_curve_chain_export_sha256")
        == identity.get("curve_chain_export_sha256")
    )


def _smoothing_replay_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        iterations = int(identity.get("iteration_count"))
        replay_iterations = int(identity.get("replay_iteration_count"))
        fixed_nodes = [int(value) for value in identity.get("fixed_node_ids", [])]
        replay_fixed_nodes = [int(value) for value in identity.get("replay_fixed_node_ids", [])]
        fixed_displacements = [
            float(value) for value in identity.get("fixed_node_displacement_m", [])
        ]
        replay_fixed_displacements = [
            float(value) for value in identity.get("replay_fixed_node_displacement_m", [])
        ]
        moved_nodes = [int(value) for value in identity.get("moved_node_ids", [])]
        replay_moved_nodes = [int(value) for value in identity.get("replay_moved_node_ids", [])]
        displacements = [float(value) for value in identity.get("node_displacement_m", [])]
        replay_displacements = [
            float(value) for value in identity.get("replay_node_displacement_m", [])
        ]
        maximum_displacement = float(identity.get("maximum_allowed_displacement_m"))
        replay_maximum_displacement = float(
            identity.get("replay_maximum_allowed_displacement_m")
        )
        quality = [float(value) for value in identity.get("quality_history", [])]
        replay_quality = [float(value) for value in identity.get("replay_quality_history", [])]
        checkpoint = int(identity.get("checkpoint_generation_id"))
        replay_checkpoint = int(identity.get("replay_checkpoint_generation_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("smoothing_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "algorithm_generation",
                "iteration_generation",
                "constraint_generation",
                "motion_generation",
                "quality_generation",
                "checkpoint_generation",
                "database_generation",
                "result_generation",
            )
        )
        and str(identity.get("algorithm") or "") in {"smart_laplacian", "optimization"}
        and identity.get("replay_algorithm") == identity.get("algorithm")
        and iterations > 0
        and replay_iterations == iterations
        and bool(fixed_nodes)
        and len(fixed_nodes) == len(set(fixed_nodes))
        and all(value > 0 for value in fixed_nodes)
        and replay_fixed_nodes == fixed_nodes
        and len(fixed_displacements) == len(fixed_nodes)
        and all(math.isfinite(value) and abs(value) <= 1.0e-15 for value in fixed_displacements)
        and replay_fixed_displacements == fixed_displacements
        and bool(moved_nodes)
        and len(moved_nodes) == len(set(moved_nodes))
        and set(moved_nodes).isdisjoint(fixed_nodes)
        and replay_moved_nodes == moved_nodes
        and len(displacements) == len(moved_nodes)
        and math.isfinite(maximum_displacement)
        and maximum_displacement > 0.0
        and all(math.isfinite(value) and 0.0 <= value <= maximum_displacement for value in displacements)
        and replay_displacements == displacements
        and replay_maximum_displacement == maximum_displacement
        and len(quality) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in quality)
        and all(quality[index] <= quality[index + 1] for index in range(len(quality) - 1))
        and quality[-1] > quality[0]
        and replay_quality == quality
        and checkpoint > 0
        and replay_checkpoint == checkpoint
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("smoothing_result_sha256"))
        and identity.get("accepted_smoothing_result_sha256")
        == identity.get("smoothing_result_sha256")
    )


def _named_entity_transfer_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        names = _mapping(identity.get("entity_names"), "entity_names")
        replay_names = _mapping(identity.get("replay_entity_names"), "replay_entity_names")
        metadata = _mapping(identity.get("metadata_attributes"), "metadata_attributes")
        replay_metadata = _mapping(
            identity.get("replay_metadata_attributes"), "replay_metadata_attributes"
        )
        groups = _mapping(identity.get("group_membership"), "group_membership")
        replay_groups = _mapping(identity.get("replay_group_membership"), "replay_group_membership")
        blocks = _mapping(identity.get("block_membership"), "block_membership")
        replay_blocks = _mapping(identity.get("replay_block_membership"), "replay_block_membership")
        sidesets = _mapping(identity.get("sideset_membership"), "sideset_membership")
        replay_sidesets = _mapping(identity.get("replay_sideset_membership"), "replay_sideset_membership")
        save_generation = int(identity.get("save_generation_id"))
        replay_save_generation = int(identity.get("replay_save_generation_id"))
        open_generation = int(identity.get("open_generation_id"))
        replay_open_generation = int(identity.get("replay_open_generation_id"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("named_entity_generation") or "")
    entity_ids = set(names)
    grouped_ids = {str(value) for values in groups.values() for value in values}
    blocked_ids = {str(value) for values in blocks.values() for value in values}
    sideset_ids = {str(value) for values in sidesets.values() for value in values}
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "name_generation",
                "metadata_generation",
                "group_generation",
                "set_generation",
                "save_generation",
                "open_generation",
                "export_generation",
                "result_generation",
            )
        )
        and bool(names)
        and all(
            str(key).startswith(("volume:", "surface:", "curve:", "vertex:"))
            and bool(str(value))
            for key, value in names.items()
        )
        and len(set(str(value) for value in names.values())) == len(names)
        and replay_names == names
        and bool(metadata)
        and set(metadata).issubset(entity_ids)
        and all(isinstance(value, Mapping) and bool(value) for value in metadata.values())
        and replay_metadata == metadata
        and bool(groups)
        and all(str(key).startswith("group:") and bool(value) for key, value in groups.items())
        and grouped_ids.issubset(entity_ids)
        and replay_groups == groups
        and bool(blocks)
        and all(str(key).startswith("block:") and bool(value) for key, value in blocks.items())
        and blocked_ids.issubset(entity_ids)
        and all(value.startswith("volume:") for value in blocked_ids)
        and replay_blocks == blocks
        and bool(sidesets)
        and all(str(key).startswith("sideset:") and bool(value) for key, value in sidesets.items())
        and sideset_ids.issubset(entity_ids)
        and all(value.startswith("surface:") for value in sideset_ids)
        and replay_sidesets == sidesets
        and save_generation > 0
        and replay_save_generation == save_generation
        and open_generation == save_generation
        and replay_open_generation == open_generation
        and str(identity.get("export_owner") or "").startswith("headless:")
        and identity.get("replay_export_owner") == identity.get("export_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("named_entity_export_sha256"))
        and identity.get("accepted_named_entity_export_sha256")
        == identity.get("named_entity_export_sha256")
    )


def _multisweep_hex_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sources = [int(value) for value in identity.get("source_surface_ids", [])]
        targets = [int(value) for value in identity.get("target_surface_ids", [])]
        chain = [int(value) for value in identity.get("sweep_volume_chain", [])]
        pairs = [[int(value) for value in pair] for pair in identity.get("section_node_correspondence", [])]
        twists = [float(value) for value in identity.get("section_twist_deg", [])]
        intervals = [int(value) for value in identity.get("chain_interval_counts", [])]
        hex_count = int(identity.get("hex_element_count"))
        minimum_jacobian = float(identity.get("minimum_scaled_jacobian"))
        allowed_jacobian = float(identity.get("minimum_allowed_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("multisweep_generation") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "source_generation", "target_generation", "section_generation",
            "twist_generation", "interval_generation", "quality_generation",
            "block_generation", "export_generation", "result_generation",
        ))
        and len(sources) == len(targets) == len(chain) >= 2
        and len(set(sources)) == len(sources) and len(set(targets)) == len(targets)
        and set(sources).isdisjoint(targets) and len(set(chain)) == len(chain)
        and identity.get("result_source_surface_ids") == sources
        and identity.get("result_target_surface_ids") == targets
        and identity.get("result_sweep_volume_chain") == chain
        and len(pairs) >= 4 and all(len(pair) == 2 and min(pair) > 0 for pair in pairs)
        and len({pair[0] for pair in pairs}) == len(pairs)
        and len({pair[1] for pair in pairs}) == len(pairs)
        and identity.get("result_section_node_correspondence") == pairs
        and len(twists) == len(chain) and all(math.isfinite(value) for value in twists)
        and identity.get("result_section_twist_deg") == twists
        and len(intervals) == len(chain) and len(set(intervals)) == 1 and intervals[0] > 0
        and identity.get("result_chain_interval_counts") == intervals
        and hex_count > 0 and identity.get("result_hex_element_count") == hex_count
        and math.isfinite(minimum_jacobian) and minimum_jacobian >= allowed_jacobian > 0.0
        and identity.get("result_minimum_scaled_jacobian") == minimum_jacobian
        and identity.get("result_minimum_allowed_scaled_jacobian") == allowed_jacobian
        and str(identity.get("block_owner") or "").startswith("block:")
        and identity.get("result_block_owner") == identity.get("block_owner")
        and _valid_sha256(identity.get("multisweep_export_sha256"))
        and identity.get("accepted_multisweep_export_sha256") == identity.get("multisweep_export_sha256")
    )


def _imprint_merge_hex_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        tolerance = float(identity.get("merge_tolerance_m"))
        vertex_pairs = [[int(value) for value in pair] for pair in identity.get("coincident_vertex_pairs", [])]
        surface_pairs = [[int(value) for value in pair] for pair in identity.get("coincident_surface_pairs", [])]
        merged = _mapping(identity.get("merged_entity_map"), "merged_entity_map")
        volumes = _mapping(identity.get("volume_owners"), "volume_owners")
        blocks = _mapping(identity.get("block_membership"), "block_membership")
        sidesets = _mapping(identity.get("sideset_membership"), "sideset_membership")
        hex_count = int(identity.get("hex_element_count"))
        minimum_jacobian = float(identity.get("minimum_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("imprint_merge_generation") or "")
    block_entities = {str(value) for values in blocks.values() for value in values}
    sideset_entities = {str(value) for values in sidesets.values() for value in values}
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "tolerance_generation", "coincident_generation", "topology_generation",
            "volume_generation", "set_generation", "quality_generation",
            "export_generation", "result_generation",
        ))
        and math.isfinite(tolerance) and 0.0 < tolerance <= 1.0e-4
        and bool(vertex_pairs) and all(len(pair) == 2 and pair[0] != pair[1] and min(pair) > 0 for pair in vertex_pairs)
        and bool(surface_pairs) and all(len(pair) == 2 and pair[0] != pair[1] and min(pair) > 0 for pair in surface_pairs)
        and identity.get("result_merge_tolerance_m") == tolerance
        and identity.get("result_coincident_vertex_pairs") == vertex_pairs
        and identity.get("result_coincident_surface_pairs") == surface_pairs
        and bool(merged) and identity.get("result_merged_entity_map") == merged
        and bool(volumes) and all(str(key).startswith("volume:") and bool(str(value)) for key, value in volumes.items())
        and identity.get("result_volume_owners") == volumes
        and bool(blocks) and all(str(key).startswith("block:") for key in blocks)
        and block_entities == set(volumes)
        and identity.get("result_block_membership") == blocks
        and bool(sidesets) and all(str(key).startswith("sideset:") for key in sidesets)
        and all(value.startswith("surface:") for value in sideset_entities)
        and identity.get("result_sideset_membership") == sidesets
        and hex_count > 0 and identity.get("result_hex_element_count") == hex_count
        and math.isfinite(minimum_jacobian) and minimum_jacobian > 0.0
        and identity.get("result_minimum_scaled_jacobian") == minimum_jacobian
        and _valid_sha256(identity.get("imprint_export_sha256"))
        and identity.get("accepted_imprint_export_sha256") == identity.get("imprint_export_sha256")
    )


def _journal_undo_redo_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("journal_generation") or "")
    transactions = [str(value) for value in identity.get("transactions", [])]
    replay_transactions = [str(value) for value in identity.get("replay_transactions", [])]
    try:
        undo_index = int(identity.get("undo_boundary_index"))
        redo_index = int(identity.get("redo_boundary_index"))
        geometry_generation = int(identity.get("geometry_generation_id"))
        replay_geometry_generation = int(identity.get("replay_geometry_generation_id"))
        remap = _mapping(identity.get("entity_id_remap"), "entity_id_remap")
        replay_remap = _mapping(identity.get("replay_entity_id_remap"), "replay_entity_id_remap")
    except (TypeError, ValueError):
        return False
    selection = [str(value) for value in identity.get("active_selection", [])]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "transaction_generation", "undo_generation", "redo_generation",
            "entity_generation", "selection_generation", "database_generation",
            "owner_generation", "result_generation",
        ))
        and len(transactions) >= 4 and replay_transactions == transactions
        and 0 <= undo_index < redo_index < len(transactions)
        and transactions[undo_index].strip().lower() == "undo"
        and transactions[redo_index].strip().lower() == "redo"
        and identity.get("replay_undo_boundary_index") == undo_index
        and identity.get("replay_redo_boundary_index") == redo_index
        and bool(remap) and len(set(str(value) for value in remap.values())) == len(remap)
        and replay_remap == remap
        and bool(selection) and set(selection).issubset(str(value) for value in remap.values())
        and identity.get("replay_active_selection") == selection
        and geometry_generation > 0 and replay_geometry_generation == geometry_generation
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and str(identity.get("session_owner") or "").startswith("batch:")
        and identity.get("replay_session_owner") == identity.get("session_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("journal_result_sha256"))
        and identity.get("accepted_journal_result_sha256") == identity.get("journal_result_sha256")
    )


def _exodus_transient_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("exodus_generation") or "")
    try:
        times = [float(value) for value in identity.get("time_steps_s", [])]
        replay_times = [float(value) for value in identity.get("replay_time_steps_s", [])]
        nodal_values = [[float(value) for value in row] for row in identity.get("nodal_values", [])]
        replay_nodal_values = [[float(value) for value in row] for row in identity.get("replay_nodal_values", [])]
        element_values = [[float(value) for value in row] for row in identity.get("element_values", [])]
        replay_element_values = [[float(value) for value in row] for row in identity.get("replay_element_values", [])]
        mesh_generation = int(identity.get("mesh_generation_id"))
        replay_mesh_generation = int(identity.get("replay_mesh_generation_id"))
    except (TypeError, ValueError):
        return False
    nodal_names = [str(value) for value in identity.get("nodal_variable_names", [])]
    element_names = [str(value) for value in identity.get("element_variable_names", [])]
    truth = identity.get("block_truth_table")
    sidesets = identity.get("sideset_values")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "time_generation", "nodal_generation", "element_generation",
            "truth_generation", "sideset_generation", "mesh_generation",
            "database_generation", "result_generation",
        ))
        and len(times) >= 2 and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(times[index] < times[index + 1] for index in range(len(times) - 1))
        and replay_times == times
        and bool(nodal_names) and len(set(nodal_names)) == len(nodal_names)
        and identity.get("replay_nodal_variable_names") == nodal_names
        and len(nodal_values) == len(times) and all(row and all(math.isfinite(value) for value in row) for row in nodal_values)
        and replay_nodal_values == nodal_values
        and bool(element_names) and len(set(element_names)) == len(element_names)
        and identity.get("replay_element_variable_names") == element_names
        and len(element_values) == len(times) and all(row and all(math.isfinite(value) for value in row) for row in element_values)
        and replay_element_values == element_values
        and isinstance(truth, Mapping) and bool(truth)
        and identity.get("replay_block_truth_table") == truth
        and isinstance(sidesets, Mapping) and bool(sidesets)
        and all(len(values) == len(times) for values in sidesets.values())
        and identity.get("replay_sideset_values") == sidesets
        and mesh_generation > 0 and replay_mesh_generation == mesh_generation
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("exodus_result_sha256"))
        and identity.get("accepted_exodus_result_sha256") == identity.get("exodus_result_sha256")
    )


def _periodic_hex_surface_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source = int(identity.get("source_surface_id"))
        target = int(identity.get("target_surface_id"))
        transform = [[float(value) for value in row] for row in identity.get("rigid_transform", [])]
        result_transform = [[float(value) for value in row] for row in identity.get("result_rigid_transform", [])]
        node_map = [[int(value) for value in row] for row in identity.get("node_map", [])]
        result_node_map = [[int(value) for value in row] for row in identity.get("result_node_map", [])]
        source_normal = [float(value) for value in identity.get("source_face_normal", [])]
        result_source_normal = [float(value) for value in identity.get("result_source_face_normal", [])]
        target_normal = [float(value) for value in identity.get("target_face_normal", [])]
        result_target_normal = [float(value) for value in identity.get("result_target_face_normal", [])]
        source_intervals = [int(value) for value in identity.get("source_interval_counts", [])]
        result_source_intervals = [int(value) for value in identity.get("result_source_interval_counts", [])]
        target_intervals = [int(value) for value in identity.get("target_interval_counts", [])]
        result_target_intervals = [int(value) for value in identity.get("result_target_interval_counts", [])]
        quality = float(identity.get("minimum_scaled_jacobian"))
        result_quality = float(identity.get("result_minimum_scaled_jacobian"))
        allowed_quality = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed_quality = float(identity.get("result_minimum_allowed_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("periodic_generation") or "")
    if not (
        len(transform) == len(result_transform) == 4
        and all(len(row) == 4 for row in transform + result_transform)
        and result_transform == transform
        and len(source_normal) == len(target_normal) == 3
    ):
        return False
    rotation = [row[:3] for row in transform[:3]]
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    rotated_source = [
        sum(rotation[row][column] * source_normal[column] for column in range(3))
        for row in range(3)
    ]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "surface_generation", "transform_generation", "nodemap_generation",
            "orientation_generation", "interval_generation", "quality_generation",
            "block_generation", "export_generation", "result_generation",
        ))
        and source > 0 and target > 0 and source != target
        and identity.get("result_source_surface_id") == source
        and identity.get("result_target_surface_id") == target
        and all(math.isfinite(value) for row in transform for value in row)
        and all(math.isclose(sum(left * right for left, right in zip(rotation[i], rotation[j], strict=True)), 1.0 if i == j else 0.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for i in range(3) for j in range(3))
        and math.isclose(determinant, 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and transform[3] == [0.0, 0.0, 0.0, 1.0]
        and len(node_map) >= 4 and result_node_map == node_map
        and all(len(pair) == 2 and min(pair) > 0 for pair in node_map)
        and len({pair[0] for pair in node_map}) == len(node_map)
        and len({pair[1] for pair in node_map}) == len(node_map)
        and result_source_normal == source_normal and result_target_normal == target_normal
        and math.isclose(sum(value * value for value in source_normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(sum(value * value for value in target_normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(math.isclose(target_normal[index], -rotated_source[index], rel_tol=1.0e-12, abs_tol=1.0e-12) for index in range(3))
        and bool(source_intervals) and source_intervals == target_intervals
        and result_source_intervals == source_intervals
        and result_target_intervals == target_intervals
        and all(value > 0 for value in source_intervals)
        and math.isfinite(quality) and quality >= allowed_quality > 0.0
        and result_quality == quality and result_allowed_quality == allowed_quality
        and str(identity.get("block_owner") or "").startswith("block:")
        and identity.get("result_block_owner") == identity.get("block_owner")
        and _valid_sha256(identity.get("periodic_export_sha256"))
        and identity.get("accepted_periodic_export_sha256") == identity.get("periodic_export_sha256")
    )


def _boundary_layer_hex_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        first = float(identity.get("first_layer_thickness_m"))
        result_first = float(identity.get("result_first_layer_thickness_m"))
        growth = float(identity.get("growth_ratio"))
        result_growth = float(identity.get("result_growth_ratio"))
        count = int(identity.get("layer_count"))
        result_count = int(identity.get("result_layer_count"))
        total = float(identity.get("total_layer_thickness_m"))
        result_total = float(identity.get("result_total_layer_thickness_m"))
        normal = [float(value) for value in identity.get("wall_normal", [])]
        result_normal = [float(value) for value in identity.get("result_wall_normal", [])]
        direction = [float(value) for value in identity.get("layer_direction", [])]
        result_direction = [float(value) for value in identity.get("result_layer_direction", [])]
        quality = float(identity.get("minimum_scaled_jacobian"))
        result_quality = float(identity.get("result_minimum_scaled_jacobian"))
        allowed_quality = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed_quality = float(identity.get("result_minimum_allowed_scaled_jacobian"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("boundarylayer_generation") or "")
    topology = [str(value) for value in identity.get("transition_topology", [])]
    result_topology = [str(value) for value in identity.get("result_transition_topology", [])]
    expected_total = first * count if math.isclose(growth, 1.0) else first * (growth**count - 1.0) / (growth - 1.0)
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "thickness_generation", "growth_generation", "layer_generation",
            "transition_generation", "normal_generation", "quality_generation",
            "block_generation", "export_generation", "result_generation",
        ))
        and math.isfinite(first) and first > 0.0 and result_first == first
        and math.isfinite(growth) and growth >= 1.0 and result_growth == growth
        and count >= 2 and result_count == count
        and math.isfinite(total) and math.isclose(total, expected_total, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_total == total
        and len(topology) == count and result_topology == topology
        and topology[0] == "hex" and all(value in {"hex", "pyramid", "tet"} for value in topology)
        and len(normal) == len(direction) == 3
        and result_normal == normal and result_direction == direction
        and math.isclose(sum(value * value for value in normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(sum(value * value for value in direction), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and sum(left * right for left, right in zip(normal, direction, strict=True)) > 1.0 - 1.0e-12
        and math.isfinite(quality) and quality >= allowed_quality > 0.0
        and result_quality == quality and result_allowed_quality == allowed_quality
        and str(identity.get("block_owner") or "").startswith("block:")
        and identity.get("result_block_owner") == identity.get("block_owner")
        and _valid_sha256(identity.get("boundarylayer_export_sha256"))
        and identity.get("accepted_boundarylayer_export_sha256") == identity.get("boundarylayer_export_sha256")
    )


def _hybrid_tet_hex_pyramid_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("hybrid_generation") or "")
    try:
        counts = {
            name: int(identity.get(name))
            for name in (
                "hex_element_count", "tet_element_count", "pyramid_element_count",
                "hex_pyramid_interface_face_count",
                "pyramid_tet_interface_face_count",
            )
        }
        result_counts = {
            name: int(identity.get(f"result_{name}")) for name in counts
        }
        orientations = [
            float(value)
            for value in identity.get("interface_orientation_dot_products", [])
        ]
        result_orientations = [
            float(value)
            for value in identity.get(
                "result_interface_orientation_dot_products", []
            )
        ]
        allowed = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed = float(
            identity.get("result_minimum_allowed_scaled_jacobian")
        )
        qualities = {
            kind: float(identity.get(f"minimum_{kind}_scaled_jacobian"))
            for kind in ("hex", "pyramid", "tet")
        }
        result_qualities = {
            kind: float(identity.get(f"result_minimum_{kind}_scaled_jacobian"))
            for kind in qualities
        }
        blocks = _mapping(identity.get("block_membership"), "block_membership")
        result_blocks = _mapping(
            identity.get("result_block_membership"), "result_block_membership"
        )
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "transition_generation", "interface_generation",
            "orientation_generation", "quality_generation", "block_generation",
            "export_generation", "result_generation",
        ))
        and min(counts.values()) > 0
        and result_counts == counts
        and counts["hex_pyramid_interface_face_count"]
        == counts["pyramid_element_count"]
        and counts["pyramid_tet_interface_face_count"]
        == 4 * counts["pyramid_element_count"]
        and bool(orientations)
        and result_orientations == orientations
        and all(math.isfinite(value) and math.isclose(value, -1.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for value in orientations)
        and math.isfinite(allowed)
        and allowed > 0.0
        and result_allowed == allowed
        and all(math.isfinite(value) and value >= allowed for value in qualities.values())
        and result_qualities == qualities
        and bool(blocks)
        and result_blocks == blocks
        and all(str(key).startswith("block:") and bool(value) for key, value in blocks.items())
        and str(identity.get("mesh_owner") or "").startswith("headless:")
        and identity.get("result_mesh_owner") == identity.get("mesh_owner")
        and _valid_sha256(identity.get("hybrid_export_sha256"))
        and identity.get("accepted_hybrid_export_sha256")
        == identity.get("hybrid_export_sha256")
    )


def _webcut_multivolume_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("webcut_generation") or "")
    try:
        volumes = [int(value) for value in identity.get("volume_ids", [])]
        result_volumes = [
            int(value) for value in identity.get("result_volume_ids", [])
        ]
        pairs = [
            [int(value) for value in row]
            for row in identity.get("shared_face_pairs", [])
        ]
        result_pairs = [
            [int(value) for value in row]
            for row in identity.get("result_shared_face_pairs", [])
        ]
        blocks = _mapping(identity.get("block_connectivity"), "block_connectivity")
        result_blocks = _mapping(
            identity.get("result_block_connectivity"), "result_block_connectivity"
        )
        hex_count = int(identity.get("hex_element_count"))
        result_hex_count = int(identity.get("result_hex_element_count"))
        euler = int(identity.get("topology_euler_characteristic"))
        result_euler = int(identity.get("result_topology_euler_characteristic"))
        quality = float(identity.get("minimum_scaled_jacobian"))
        result_quality = float(identity.get("result_minimum_scaled_jacobian"))
        allowed = float(identity.get("minimum_allowed_scaled_jacobian"))
        result_allowed = float(
            identity.get("result_minimum_allowed_scaled_jacobian")
        )
    except (TypeError, ValueError):
        return False
    adjacency = {volume: set() for volume in volumes}
    valid_pairs = bool(pairs)
    for row in pairs:
        if len(row) != 3 or row[0] == row[1] or row[0] not in adjacency or row[1] not in adjacency or row[2] <= 0:
            valid_pairs = False
            continue
        adjacency[row[0]].add(row[1])
        adjacency[row[1]].add(row[0])
    visited: set[int] = set()
    pending = [volumes[0]] if volumes else []
    while pending:
        volume = pending.pop()
        if volume in visited:
            continue
        visited.add(volume)
        pending.extend(adjacency[volume] - visited)
    block_volumes = [int(value) for values in blocks.values() for value in values]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "topology_generation", "sharedface_generation",
            "connectivity_generation", "mesh_generation", "quality_generation",
            "block_generation", "export_generation", "result_generation",
        ))
        and len(volumes) >= 2
        and len(set(volumes)) == len(volumes)
        and min(volumes) > 0
        and result_volumes == volumes
        and valid_pairs
        and result_pairs == pairs
        and len(visited) == len(volumes)
        and len({row[2] for row in pairs}) == len(pairs)
        and bool(blocks)
        and result_blocks == blocks
        and all(str(key).startswith("block:") and bool(value) for key, value in blocks.items())
        and sorted(block_volumes) == sorted(volumes)
        and len(set(block_volumes)) == len(block_volumes)
        and hex_count > 0
        and result_hex_count == hex_count
        and euler == 1
        and result_euler == euler
        and math.isfinite(quality)
        and quality >= allowed > 0.0
        and result_quality == quality
        and result_allowed == allowed
        and str(identity.get("mesh_owner") or "").startswith("headless:")
        and identity.get("result_mesh_owner") == identity.get("mesh_owner")
        and _valid_sha256(identity.get("webcut_export_sha256"))
        and identity.get("accepted_webcut_export_sha256")
        == identity.get("webcut_export_sha256")
    )


def _sideset_normal_propagation_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("sideset_generation") or "")
    try:
        sidesets = _mapping(identity.get("sideset_membership"), "sideset_membership")
        replay_sidesets = _mapping(
            identity.get("replay_sideset_membership"), "replay_sideset_membership"
        )
        normals = _mapping(identity.get("outward_normals"), "outward_normals")
        replay_normals = _mapping(
            identity.get("replay_outward_normals"), "replay_outward_normals"
        )
        merge_map = _mapping(identity.get("merge_entity_map"), "merge_entity_map")
        replay_merge_map = _mapping(
            identity.get("replay_merge_entity_map"), "replay_merge_entity_map"
        )
        blocks = _mapping(identity.get("block_membership"), "block_membership")
        replay_blocks = _mapping(
            identity.get("replay_block_membership"), "replay_block_membership"
        )
        mesh_generation = int(identity.get("mesh_generation_id"))
        replay_mesh_generation = int(identity.get("replay_mesh_generation_id"))
        parsed_normals = {
            str(key): [float(value) for value in values]
            for key, values in normals.items()
        }
    except (TypeError, ValueError):
        return False
    sideset_surfaces = {
        f"surface:{int(value)}" for values in sidesets.values() for value in values
    }
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "normal_generation", "merge_generation", "mesh_generation",
            "entity_generation", "block_generation", "database_generation",
            "export_generation", "result_generation",
        ))
        and bool(sidesets)
        and replay_sidesets == sidesets
        and all(str(key).startswith("sideset:") and bool(value) for key, value in sidesets.items())
        and bool(parsed_normals)
        and replay_normals == normals
        and set(parsed_normals) == sideset_surfaces
        and all(len(vector) == 3 and all(math.isfinite(value) for value in vector) and math.isclose(sum(value * value for value in vector), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for vector in parsed_normals.values())
        and bool(merge_map)
        and replay_merge_map == merge_map
        and all(str(key).startswith("surface:") and str(value).startswith("surface:") for key, value in merge_map.items())
        and mesh_generation > 0
        and replay_mesh_generation == mesh_generation
        and bool(blocks)
        and replay_blocks == blocks
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and _valid_sha256(identity.get("mesh_export_sha256"))
        and identity.get("replay_mesh_export_sha256") == identity.get("mesh_export_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _cad_import_healing_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("cad_generation") or "")
    try:
        scale = float(identity.get("unit_scale_to_m"))
        replay_scale = float(identity.get("replay_unit_scale_to_m"))
        tolerance = float(identity.get("healing_tolerance_m"))
        replay_tolerance = float(identity.get("replay_healing_tolerance_m"))
        body_count = int(identity.get("body_count"))
        replay_body_count = int(identity.get("replay_body_count"))
        watertight = [int(value) for value in identity.get("watertight_body_ids", [])]
        replay_watertight = [
            int(value) for value in identity.get("replay_watertight_body_ids", [])
        ]
        topology = _mapping(identity.get("topology_counts"), "topology_counts")
        replay_topology = _mapping(
            identity.get("replay_topology_counts"), "replay_topology_counts"
        )
        entity_generation = int(identity.get("entity_generation_id"))
        replay_entity_generation = int(identity.get("replay_entity_generation_id"))
    except (TypeError, ValueError):
        return False
    expected_scales = {"m": 1.0, "cm": 0.01, "mm": 0.001, "in": 0.0254}
    source_unit = str(identity.get("source_unit") or "")
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "unit_generation", "healing_generation", "topology_generation",
            "entity_generation", "database_generation", "owner_generation",
            "result_generation",
        ))
        and source_unit in expected_scales
        and identity.get("replay_source_unit") == source_unit
        and identity.get("model_unit") == "m"
        and identity.get("replay_model_unit") == "m"
        and math.isclose(scale, expected_scales[source_unit], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and replay_scale == scale
        and math.isfinite(tolerance)
        and 0.0 < tolerance < 1.0e-3
        and replay_tolerance == tolerance
        and body_count > 0
        and replay_body_count == body_count
        and sorted(watertight) == list(range(1, body_count + 1))
        and replay_watertight == watertight
        and topology.get("volume") == body_count
        and all(int(value) > 0 for value in topology.values())
        and replay_topology == topology
        and entity_generation > 0
        and replay_entity_generation == entity_generation
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and _valid_sha256(identity.get("source_cad_sha256"))
        and identity.get("replay_source_cad_sha256") == identity.get("source_cad_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _python_batch_rollback_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("batch_generation") or "")
    commands = [str(value) for value in identity.get("commands", [])]
    replay_commands = [str(value) for value in identity.get("replay_commands", [])]
    try:
        exception_index = int(identity.get("exception_boundary_index"))
        replay_exception_index = int(identity.get("replay_exception_boundary_index"))
        before = _mapping(identity.get("entity_generations_before"), "entity_generations_before")
        replay_before = _mapping(identity.get("replay_entity_generations_before"), "replay_entity_generations_before")
        after = _mapping(identity.get("entity_generations_after"), "entity_generations_after")
        replay_after = _mapping(identity.get("replay_entity_generations_after"), "replay_entity_generations_after")
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "exception_generation", "rollback_generation", "entity_generation",
            "model_generation", "database_generation", "session_generation",
            "commandlog_generation", "result_generation",
        ))
        and len(commands) >= 2 and replay_commands == commands
        and 0 <= exception_index < len(commands) and replay_exception_index == exception_index
        and bool(str(identity.get("exception_type") or ""))
        and identity.get("replay_exception_type") == identity.get("exception_type")
        and identity.get("transaction_committed") is False
        and identity.get("replay_transaction_committed") is False
        and bool(before) and replay_before == before and after == before and replay_after == after
        and all(int(value) > 0 for value in before.values())
        and bool(str(identity.get("active_model") or ""))
        and identity.get("replay_active_model") == identity.get("active_model")
        and str(identity.get("database_owner") or "").startswith("headless:")
        and identity.get("replay_database_owner") == identity.get("database_owner")
        and str(identity.get("session_owner") or "").startswith("batch:")
        and identity.get("replay_session_owner") == identity.get("session_owner")
        and _valid_sha256(identity.get("command_log_sha256"))
        and identity.get("replay_command_log_sha256") == identity.get("command_log_sha256")
        and _valid_sha256(identity.get("database_sha256"))
        and identity.get("replay_database_sha256") == identity.get("database_sha256")
        and _valid_sha256(identity.get("batch_result_sha256"))
        and identity.get("accepted_batch_result_sha256") == identity.get("batch_result_sha256")
    )


def _mesh_recipe_dag_identity_ok(identity: object) -> bool:
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("recipe_generation") or "")
    try:
        dag = _mapping(identity.get("dependency_dag"), "dependency_dag")
        replay_dag = _mapping(identity.get("replay_dependency_dag"), "replay_dependency_dag")
        schemes = _mapping(identity.get("scheme_assignments"), "scheme_assignments")
        replay_schemes = _mapping(identity.get("replay_scheme_assignments"), "replay_scheme_assignments")
        blocks = _mapping(identity.get("block_owners"), "block_owners")
        replay_blocks = _mapping(identity.get("replay_block_owners"), "replay_block_owners")
        sidesets = _mapping(identity.get("sideset_owners"), "sideset_owners")
        replay_sidesets = _mapping(identity.get("replay_sideset_owners"), "replay_sideset_owners")
        export_generation = int(identity.get("export_generation_id"))
        replay_export_generation = int(identity.get("replay_export_generation_id"))
    except (TypeError, ValueError):
        return False
    order = [str(value) for value in identity.get("execution_order", [])]
    replay_order = [str(value) for value in identity.get("replay_execution_order", [])]
    nodes = [str(value) for value in dag]
    positions = {node: index for index, node in enumerate(order)}
    dependencies_valid = all(
        isinstance(dependencies, list)
        and all(str(dependency) in positions and positions[str(dependency)] < positions[str(node)] for dependency in dependencies)
        for node, dependencies in dag.items()
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "dag_generation", "parameter_generation", "execution_generation",
            "scheme_generation", "set_generation", "export_generation",
            "owner_generation", "result_generation",
        ))
        and bool(dag) and replay_dag == dag
        and len(order) == len(nodes) and set(order) == set(nodes) and replay_order == order
        and dependencies_valid
        and _valid_sha256(identity.get("parameter_sha256"))
        and identity.get("replay_parameter_sha256") == identity.get("parameter_sha256")
        and bool(schemes) and replay_schemes == schemes
        and all(str(key).startswith(("volume:", "surface:")) and bool(str(value)) for key, value in schemes.items())
        and bool(blocks) and replay_blocks == blocks
        and all(str(key).startswith("block:") and bool(value) for key, value in blocks.items())
        and bool(sidesets) and replay_sidesets == sidesets
        and all(str(key).startswith("sideset:") and bool(value) for key, value in sidesets.items())
        and export_generation > 0 and replay_export_generation == export_generation
        and str(identity.get("recipe_owner") or "").startswith("headless:")
        and identity.get("replay_recipe_owner") == identity.get("recipe_owner")
        and _valid_sha256(identity.get("recipe_export_sha256"))
        and identity.get("accepted_recipe_export_sha256") == identity.get("recipe_export_sha256")
    )


_V44_PUBLIC_KEY = (
    "hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_generation_identity"
)
_V44_SOURCE_JOURNAL_KEY = (
    "journal_replay_command_order_session_units_geometry_generation_database_owner_identity"
)
_V44_SOURCE_QUALITY_KEY = (
    "mesh_quality_metric_reference_element_dimension_block_export_owner_identity"
)


def _v44_positive_integer(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    if isinstance(value, float) and not value.is_integer():
        return None
    return parsed if parsed > 0 else None


def _v44_finite_number(value: object) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _v44_id_list(value: object) -> list[int] | None:
    if not isinstance(value, list) or not value:
        return None
    if not all(
        isinstance(item, int) and not isinstance(item, bool) and item > 0
        for item in value
    ):
        return None
    return value if len(set(value)) == len(value) else None


def _v44_numeric_matrix(value: object) -> list[list[float]] | None:
    if not isinstance(value, list) or not value:
        return None
    try:
        rows = [[float(item) for item in row] for row in value]
    except (TypeError, ValueError):
        return None
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        return None
    return rows if all(math.isfinite(item) for row in rows for item in row) else None


def _v44_result(policy: str, checks: Mapping[str, bool], generation: str) -> dict[str, object]:
    failed = sorted(name for name, ok in checks.items() if not ok)
    return {
        "policy": policy,
        "status": "ok" if not failed else "needs_attention",
        "checks": dict(checks),
        "issues": failed,
        "generation": generation,
    }


def _v44_public_gate(summary: Mapping[str, object]) -> dict[str, object]:
    identity = _mapping(summary.get(_V44_PUBLIC_KEY), _V44_PUBLIC_KEY)
    generation = str(identity.get("sweep_generation") or "")
    generation_names = (
        "periodic_generation", "interface_generation", "quality_generation",
        "jacobian_generation", "block_generation", "sideset_generation",
        "export_generation", "result_generation",
    )
    source_surface = _v44_positive_integer(identity.get("source_surface_id"))
    target_surface = _v44_positive_integer(identity.get("target_surface_id"))
    source_node_count = _v44_positive_integer(identity.get("source_node_count"))
    target_node_count = _v44_positive_integer(identity.get("target_node_count"))
    source_node_ids = _v44_id_list(identity.get("paired_source_node_ids"))
    target_node_ids = _v44_id_list(identity.get("paired_target_node_ids"))
    transform = _v44_numeric_matrix(identity.get("periodic_transform_matrix"))
    normal_dot = _v44_finite_number(identity.get("interface_normal_dot"))
    minimum_jacobian = _v44_finite_number(identity.get("minimum_scaled_jacobian"))
    allowed_jacobian = _v44_finite_number(
        identity.get("minimum_allowed_scaled_jacobian")
    )
    checks: dict[str, bool] = {
        "generation_lineage": bool(generation)
        and all(identity.get(name) == generation for name in generation_names),
        "periodic_surface_pair": (
            source_surface is not None
            and target_surface is not None
            and source_surface != target_surface
            and identity.get("source_surface_id")
            == identity.get("result_source_surface_id")
            and identity.get("target_surface_id") == identity.get("result_target_surface_id")
            and source_node_count is not None
            and target_node_count is not None
            and identity.get("result_source_node_count") == identity.get("source_node_count")
            and identity.get("result_target_node_count") == identity.get("target_node_count")
        ),
        "periodic_node_pairing": (
            source_node_ids is not None
            and target_node_ids is not None
            and source_node_count is not None
            and target_node_count is not None
            and identity.get("result_paired_source_node_ids") == source_node_ids
            and identity.get("result_paired_target_node_ids") == target_node_ids
            and len(source_node_ids) == len(target_node_ids)
            and len(source_node_ids) <= source_node_count
            and len(target_node_ids) <= target_node_count
        ),
        "coordinate_frame_and_transform": (
            identity.get("coordinate_frame") == "global_cartesian"
            and identity.get("result_coordinate_frame") == identity.get("coordinate_frame")
            and transform is not None
            and identity.get("result_periodic_transform_matrix")
            == identity.get("periodic_transform_matrix")
        ),
        "opposed_interface_normals": (
            normal_dot is not None
            and math.isclose(normal_dot, -1.0, abs_tol=1.0e-12)
            and identity.get("result_interface_normal_dot")
            == identity.get("interface_normal_dot")
        ),
        "positive_jacobian": (
            minimum_jacobian is not None
            and allowed_jacobian is not None
            and minimum_jacobian >= allowed_jacobian > 0.0
            and identity.get("result_minimum_scaled_jacobian")
            == identity.get("minimum_scaled_jacobian")
        ),
        "block_and_sideset_ownership": (
            identity.get("result_block_membership") == identity.get("block_membership")
            and identity.get("result_sideset_membership") == identity.get("sideset_membership")
            and bool(identity.get("block_membership"))
            and bool(identity.get("sideset_membership"))
        ),
        "headless_owner": (
            str(identity.get("mesh_owner") or "").startswith("headless:")
            and identity.get("result_mesh_owner") == identity.get("mesh_owner")
        ),
        "export_digest": (
            _valid_sha256(identity.get("mesh_export_sha256"))
            and identity.get("accepted_mesh_export_sha256")
            == identity.get("mesh_export_sha256")
        ),
    }
    return _v44_result("cubit_periodic_hex_transition_v44_gate_v1", checks, generation)


def _v44_source_journal_gate(summary: Mapping[str, object]) -> dict[str, object]:
    identity = _mapping(summary.get(_V44_SOURCE_JOURNAL_KEY), _V44_SOURCE_JOURNAL_KEY)
    generation = str(identity.get("journal_generation") or "")
    generation_names = (
        "command_order_generation", "session_units_generation", "geometry_generation",
        "status_generation", "database_generation", "mesh_generation", "result_generation",
    )
    commands = identity.get("journal_commands")
    replay_commands = identity.get("replay_journal_commands")
    statuses = identity.get("command_status")
    replay_statuses = identity.get("replay_command_status")
    geometry_generation_id = _v44_positive_integer(
        identity.get("geometry_generation_id")
    )
    checks = {
        "generation_lineage": bool(generation)
        and all(identity.get(name) == generation for name in generation_names),
        "journal_order_and_units": (
            isinstance(commands, list)
            and bool(commands)
            and all(isinstance(command, str) and bool(command.strip()) for command in commands)
            and replay_commands == commands
            and identity.get("session_units") == "mm"
            and identity.get("replay_session_units") == identity.get("session_units")
        ),
        "command_status": (
            isinstance(statuses, list)
            and statuses == ["success"] * len(statuses)
            and replay_statuses == statuses
            and bool(statuses)
        ),
        "geometry_generation": (
            geometry_generation_id is not None
            and identity.get("replay_geometry_generation_id")
            == identity.get("geometry_generation_id")
        ),
        "headless_database_owner": (
            str(identity.get("database_owner") or "").startswith("headless:")
            and identity.get("replay_database_owner") == identity.get("database_owner")
        ),
        "mesh_export_digest": (
            _valid_sha256(identity.get("mesh_export_sha256"))
            and identity.get("replay_mesh_export_sha256")
            == identity.get("mesh_export_sha256")
        ),
        "result_digest": (
            _valid_sha256(identity.get("result_sha256"))
            and identity.get("replay_result_sha256") == identity.get("result_sha256")
            and identity.get("accepted_result_sha256") == identity.get("result_sha256")
        ),
    }
    return _v44_result("cubit_v44_source_replay_quality_gate_v1", checks, generation)


def _v44_source_quality_gate(summary: Mapping[str, object]) -> dict[str, object]:
    identity = _mapping(summary.get(_V44_SOURCE_QUALITY_KEY), _V44_SOURCE_QUALITY_KEY)
    generation = str(identity.get("quality_generation") or "")
    generation_names = (
        "reference_element_generation", "dimension_generation", "metric_generation",
        "block_generation", "export_generation", "owner_generation", "result_generation",
    )
    dimension = _v44_positive_integer(identity.get("dimension"))
    minimum_jacobian = _v44_finite_number(identity.get("minimum_scaled_jacobian"))
    export_generation_id = _v44_positive_integer(
        identity.get("export_generation_id")
    )
    checks = {
        "generation_lineage": bool(generation)
        and all(identity.get(name) == generation for name in generation_names),
        "reference_element_and_dimension": (
            identity.get("reference_element") == "hex8"
            and identity.get("replay_reference_element") == identity.get("reference_element")
            and dimension == 3
            and identity.get("replay_dimension") == identity.get("dimension")
        ),
        "metric_definition_and_value": (
            identity.get("metric_definition") == "scaled_jacobian"
            and identity.get("replay_metric_definition") == identity.get("metric_definition")
            and minimum_jacobian is not None
            and minimum_jacobian > 0.0
            and identity.get("replay_minimum_scaled_jacobian")
            == identity.get("minimum_scaled_jacobian")
        ),
        "block_membership": (
            bool(identity.get("block_membership"))
            and identity.get("replay_block_membership") == identity.get("block_membership")
        ),
        "export_generation": (
            export_generation_id is not None
            and identity.get("replay_export_generation_id")
            == identity.get("export_generation_id")
        ),
        "headless_owner": (
            str(identity.get("database_owner") or "").startswith("headless:")
            and identity.get("replay_database_owner") == identity.get("database_owner")
        ),
        "quality_export_digest": (
            _valid_sha256(identity.get("quality_export_sha256"))
            and identity.get("accepted_quality_export_sha256")
            == identity.get("quality_export_sha256")
        ),
    }
    return _v44_result("cubit_v44_source_replay_quality_gate_v1", checks, generation)


def _v44_source_combined_gate(summary: Mapping[str, object]) -> dict[str, object]:
    results: list[tuple[str, dict[str, object]]] = []
    if _V44_SOURCE_JOURNAL_KEY in summary:
        results.append(("journal", _v44_source_journal_gate(summary)))
    if _V44_SOURCE_QUALITY_KEY in summary:
        results.append(("quality", _v44_source_quality_gate(summary)))
    checks: dict[str, bool] = {}
    generations: list[str] = []
    for prefix, result in results:
        checks.update(
            {
                f"{prefix}.{name}": ok
                for name, ok in result.get("checks", {}).items()
            }
        )
        generations.append(str(result.get("generation") or ""))
    return _v44_result(
        "cubit_v44_source_replay_quality_gate_v1",
        checks,
        "+".join(generations),
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
    if _V44_PUBLIC_KEY in summary:
        return _v44_public_gate(summary)
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
        "sheet_pillow_layers_use_current_topology_blocks_interfaces_orientation_and_jacobians": (
            _sheet_pillow_layer_identity_ok(
                summary.get("hex_sheet_pillow_layer_topology_block_interface_orientation_jacobian_identity")
            )
        ),
        "pyramid_transitions_use_current_base_sides_interface_jacobian_blocks_and_export": (
            _pyramid_transition_export_identity_ok(
                summary.get("pyramid_transition_orientation_interface_jacobian_block_export_identity")
            )
        ),
        "spine_sweeps_use_current_frame_twist_intervals_surfaces_block_and_quality": (
            _spine_sweep_identity_ok(
                summary.get("hex_sweep_spine_frame_twist_interval_surface_block_quality_identity")
            )
        ),
        "local_refinements_use_current_parent_child_transition_conformity_blocks_boundaries_and_export": (
            _local_refinement_identity_ok(
                summary.get("local_refinement_parent_child_transition_conformity_block_boundary_export_identity")
            )
        ),
        "mapped_hexes_use_current_intervals_valence_face_pairing_logical_coordinates_jacobians_and_blocks": (
            _hex_map_identity_ok(
                summary.get("hex_map_curve_interval_parity_vertex_valence_face_pairing_logical_jacobian_block_generation_identity")
            )
        ),
        "mesh_morphs_use_current_boundaries_constraints_frame_smoothing_jacobians_sidesets_and_export": (
            _mesh_morph_identity_ok(
                summary.get("mesh_morph_boundary_displacement_constraint_frame_smoothing_jacobian_sideset_export_generation_identity")
            )
        ),
        "swept_hexes_use_current_twist_pairing_intervals_edge_map_orientation_jacobians_block_and_mesh": (
            _sweep_hex_twist_identity_ok(
                summary.get("sweep_hex_twist_source_target_interval_edge_map_orientation_jacobian_block_mesh_generation_identity")
            )
        ),
        "mixed_transitions_use_current_elements_faces_nodes_orientation_conformity_quality_blocks_sidesets_and_export": (
            _mixed_transition_face_identity_ok(
                summary.get("hex_tet_pyramid_transition_face_nodes_orientation_conformity_quality_block_sideset_export_generation_identity")
            )
        ),
        "curved_mixed_transitions_use_current_orders_midnodes_permutations_parametric_faces_quadrature_jacobians_and_export": (
            _mixed_high_order_transition_identity_ok(
                summary.get("mixed_high_order_transition_midnode_permutation_parametric_face_quadrature_jacobian_export_generation_identity")
            )
        ),
        "periodic_high_order_meshes_use_current_affine_corner_edge_face_pairs_orientation_sidesets_and_residual": (
            _periodic_high_order_identity_ok(
                summary.get("periodic_high_order_affine_corner_edge_face_node_pair_orientation_sideset_residual_generation_identity")
            )
        ),
        "hex_boundary_layers_use_current_thickness_growth_count_collision_jacobian_block_sideset_and_mesh": (
            _hex_boundary_layer_identity_ok(
                summary.get("hex_boundary_layer_thickness_growth_count_collision_jacobian_block_sideset_mesh_generation_identity")
            )
        ),
        "pyramid_transitions_use_current_diagonal_orientation_faces_regions_jacobians_and_export": (
            _pyramid_transition_closure_identity_ok(
                summary.get("pyramid_transition_diagonal_orientation_shared_face_region_jacobian_export_generation_identity")
            )
        ),
        "hex_sheet_pillows_use_current_incidence_euler_shell_orientation_block_jacobian_and_mesh": (
            _hex_sheet_pillow_topology_identity_ok(
                summary.get("hex_sheet_pillow_incidence_euler_shell_orientation_block_jacobian_mesh_result_generation_identity")
            )
        ),
        "multiblock_interfaces_use_current_merge_nodes_faces_owners_sets_duplicates_jacobian_and_export": (
            _multiblock_interface_identity_ok(
                summary.get("multiblock_interface_merge_face_owner_block_sideset_duplicate_jacobian_export_generation_identity")
            )
        ),
        "high_order_hexes_use_current_family_node_roles_reference_order_jacobian_volume_faces_and_mesh": (
            _high_order_hex_family_identity_ok(
                summary.get("hex20_hex27_family_node_role_reference_order_jacobian_volume_face_mesh_result_generation_identity")
            )
        ),
        "sheet_midplanes_use_current_source_offset_thickness_normal_sets_area_mass_and_geometry": (
            _sheet_midplane_mass_identity_ok(
                summary.get("sheet_midplane_source_offset_thickness_normal_block_sideset_area_mass_geometry_result_generation_identity")
            )
        ),
        "anisotropic_hexes_use_current_metric_directions_sizes_gradation_alignment_jacobian_blocks_and_mesh": (
            _anisotropic_hex_metric_identity_ok(
                summary.get("anisotropic_hex_metric_direction_size_gradation_alignment_jacobian_block_mesh_result_generation_identity")
            )
        ),
        "curved_high_order_boundaries_use_current_normals_hausdorff_measures_jacobian_order_geometry_and_result": (
            _curved_highorder_boundary_identity_ok(
                summary.get("curved_highorder_boundary_normal_hausdorff_area_volume_jacobian_order_geometry_result_generation_identity")
            )
        ),
        "hex_sweeps_use_current_source_target_layers_bias_orientation_jacobian_boundaries_volume_and_mesh": (
            _hex_sweep_contract_identity_ok(
                summary.get("hex_sweep_source_target_topology_layer_interval_orientation_jacobian_boundary_volume_mesh_result_generation_identity")
            )
        ),
        "mixed_transitions_use_current_counts_interface_nodes_orientation_regions_volume_owner_and_mesh": (
            _mixed_transition_contract_identity_ok(
                summary.get("mixed_hex_tet_pyramid_count_interface_conformity_face_orientation_node_region_volume_mesh_result_generation_identity")
            )
        ),
        "periodic_hex_pairs_use_current_transform_node_order_orientation_jacobian_region_export_mesh_and_result": (
            _periodic_hex_contract_identity_ok(
                summary.get("periodic_hex_pair_transform_node_order_face_orientation_jacobian_region_export_mesh_result_generation_identity")
            )
        ),
        "high_order_curved_hexes_use_current_midnodes_projection_volume_jacobian_order_cad_mesh_and_export": (
            _curved_hex_contract_identity_ok(
                summary.get("high_order_hex_curve_midnode_projection_volume_jacobian_order_cad_mesh_export_generation_identity")
            )
        ),
        "midsurface_shells_use_current_facepairs_thickness_normals_area_volume_sets_geometry_and_result": (
            _midsurface_shell_contract_identity_ok(
                summary.get("midsurface_shell_facepair_thickness_normal_area_volume_block_sideset_geometry_export_result_generation_identity")
            )
        ),
        "cohesive_cracks_use_current_faces_nodes_front_normals_orientation_traction_quality_mesh_and_result": (
            _cohesive_crack_contract_identity_ok(
                summary.get("cohesive_crack_face_nodepair_front_normal_orientation_traction_block_jacobian_mesh_export_result_generation_identity")
            )
        ),
        "periodic_hex_replays_use_current_pairs_transform_edge_order_orientation_quality_sets_database_and_export": (
            _periodic_hex_replay_identity_ok(
                summary.get("periodic_hex_node_pair_transform_face_orientation_edge_order_jacobian_block_sideset_database_export_generation_identity")
            )
        ),
        "thin_sweeps_use_current_surfaces_intervals_layers_thickness_topology_orientation_owner_and_export": (
            _thin_sweep_replay_identity_ok(
                summary.get("thin_sweep_hex_source_target_interval_propagation_layer_thickness_side_topology_orientation_volume_export_generation_identity")
            )
        ),
        "medial_axis_hexes_use_current_sheet_pairs_thickness_cells_topology_intervals_quality_block_and_export": (
            _medial_axis_hex_decomposition_identity_ok(
                summary.get("medial_axis_hex_decomposition_sheet_pair_thickness_topology_interval_quality_block_export_generation_identity")
            )
        ),
        "curve_chain_hexes_use_current_order_intervals_bias_corner_parity_boundary_layers_orientation_sideset_and_export": (
            _curve_chain_boundary_layer_identity_ok(
                summary.get("curve_chain_interval_bias_corner_parity_boundary_layer_orientation_sideset_export_generation_identity")
            )
        ),
        "multisweep_hexes_use_current_sources_targets_chain_sections_twist_intervals_quality_block_and_export": (
            _multisweep_hex_identity_ok(
                summary.get("multisweep_source_target_chain_section_correspondence_twist_interval_quality_block_export_generation_identity")
            )
        ),
        "imprint_merge_hexes_use_current_tolerance_coincidence_topology_sets_quality_and_export": (
            _imprint_merge_hex_identity_ok(
                summary.get("imprint_merge_tolerance_coincident_topology_volume_block_sideset_quality_export_generation_identity")
            )
        ),
        "periodic_hex_surfaces_use_current_pair_transform_nodemap_orientation_intervals_quality_block_and_export": (
            _periodic_hex_surface_identity_ok(
                summary.get("periodic_hex_surface_pair_transform_nodemap_orientation_interval_quality_block_export_generation_identity")
            )
        ),
        "boundary_layer_hexes_use_current_thickness_growth_layers_transition_normal_quality_block_and_export": (
            _boundary_layer_hex_identity_ok(
                summary.get("boundarylayer_hex_thickness_growth_layers_transition_normal_quality_block_export_generation_identity")
            )
        ),
        "hybrid_tet_hex_pyramids_use_current_transition_interfaces_orientation_quality_blocks_owner_and_export": (
            _hybrid_tet_hex_pyramid_identity_ok(
                summary.get("hybrid_tet_hex_pyramid_transition_interface_orientation_quality_block_export_generation_identity")
            )
        ),
        "webcut_multivolumes_use_current_sharedfaces_connectivity_hexcount_euler_quality_owner_and_export": (
            _webcut_multivolume_identity_ok(
                summary.get("webcut_multivolume_sharedface_block_connectivity_hexcount_euler_quality_export_generation_identity")
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
    if _V44_SOURCE_JOURNAL_KEY in summary or _V44_SOURCE_QUALITY_KEY in summary:
        return _v44_source_combined_gate(summary)
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
        "journal_replays_use_current_includes_aprepro_scope_workdir_version_and_output": (
            _journal_invocation_identity_ok(
                summary.get("journal_include_aprepro_scope_workdir_version_output_invocation_identity")
            )
        ),
        "exodus_results_use_current_64bit_ids_qa_times_variables_mesh_and_checksum": (
            _exodus_result_identity_ok(
                summary.get("exodus_64bit_id_qa_time_nodal_variable_mesh_checksum_identity")
            )
        ),
        "sideset_exports_use_current_skin_dimension_orientation_namespace_owner_exodus_and_mesh": (
            _sideset_export_identity_ok(
                summary.get("sideset_skin_dimension_orientation_namespace_owner_exodus_mesh_identity")
            )
        ),
        "headless_python_runs_use_current_interpreter_module_transaction_undo_outputs_and_invocation": (
            _headless_python_invocation_identity_ok(
                summary.get("headless_python_interpreter_module_transaction_undo_output_invocation_identity")
            )
        ),
        "exodus_element_variables_use_current_truth_table_blocks_timesteps_layout_and_file": (
            _exodus_truth_table_identity_ok(
                summary.get("exodus_truth_table_element_variable_block_timestep_layout_file_generation_identity")
            )
        ),
        "cad_imports_use_current_classification_transform_units_heal_topology_session_and_result": (
            _cad_import_identity_ok(
                summary.get("cad_import_body_sheet_lump_transform_unit_heal_topology_session_result_generation_identity")
            )
        ),
        "journals_use_current_undo_transaction_command_order_entities_headless_session_model_and_digests": (
            _journal_transaction_identity_ok(
                summary.get("journal_undo_transaction_command_order_entity_session_model_digest_generation_identity")
            )
        ),
        "exodus_sidesets_use_current_ids_topology_orientation_distribution_blocks_time_and_digests": (
            _exodus_sideset_distribution_identity_ok(
                summary.get("exodus_sideset_distribution_factor_topology_orientation_block_time_file_generation_identity")
            )
        ),
        "boolean_imprint_merges_use_current_entity_lineage_orientation_measure_adjacency_journal_and_digest": (
            _boolean_entity_lineage_identity_ok(
                summary.get("boolean_imprint_merge_entity_lineage_surface_orientation_measure_adjacency_journal_digest_generation_identity")
            )
        ),
        "high_order_exodus_uses_current_blocks_midnodes_int64_qa_restart_sideset_owners_and_digests": (
            _exodus_high_order_restart_identity_ok(
                summary.get("exodus_high_order_block_midnode_order_int64_qa_restart_sideset_digest_generation_identity")
            )
        ),
        "journal_replays_use_current_transaction_undo_allocation_idempotency_save_restore_session_and_result": (
            _journal_transaction_restore_identity_ok(
                summary.get("journal_undo_transaction_idempotency_entity_allocation_save_restore_session_result_generation_identity")
            )
        ),
        "exodus_assemblies_use_current_membership_frame_qa_time_block_sideset_file_and_result": (
            _exodus_assembly_identity_ok(
                summary.get("exodus_assembly_membership_frame_qa_time_block_sideset_file_result_generation_identity")
            )
        ),
        "failed_commands_rollback_model_allocator_undo_session_and_result_atomically": (
            _command_failure_atomic_identity_ok(
                summary.get("command_failure_atomic_rollback_error_entity_allocator_undo_session_result_generation_identity")
            )
        ),
        "cub_roundtrips_preserve_kernel_entities_attributes_groups_mesh_model_file_and_result": (
            _cub_roundtrip_identity_ok(
                summary.get("cub_roundtrip_kernel_entity_name_attribute_group_mesh_model_file_result_generation_identity")
            )
        ),
        "sweep_replays_use_current_source_target_intervals_bias_matches_periodic_layers_journal_and_result": (
            _sweep_replay_identity_ok(
                summary.get("sweep_source_target_interval_bias_match_periodic_layer_scheme_journal_result_generation_identity")
            )
        ),
        "checkpoint_restores_use_current_partitions_owned_ghost_persistent_sets_quality_model_and_result": (
            _checkpoint_partition_identity_ok(
                summary.get("checkpoint_partition_owned_ghost_persistent_block_sideset_quality_model_result_generation_identity")
            )
        ),
        "skinned_sidesets_use_current_remesh_adjacency_normals_faces_entities_journal_and_result": (
            _sideset_skin_remesh_identity_ok(
                summary.get("sideset_skin_remesh_adjacent_block_normal_face_multiplicity_entity_journal_result_generation_identity")
            )
        ),
        "parallel_sculpt_uses_current_seed_ranks_owned_ghost_stitch_qa_connectivity_invocation_and_export": (
            _parallel_sculpt_determinism_identity_ok(
                summary.get("parallel_sculpt_seed_rank_owned_ghost_stitch_qa_connectivity_invocation_export_generation_identity")
            )
        ),
        "journal_replays_use_current_undo_groups_entities_reset_checkpoint_order_database_owner_and_result": (
            _journal_replay_contract_identity_ok(
                summary.get("journal_undo_idempotence_entity_reset_checkpoint_replay_order_database_result_generation_identity")
            )
        ),
        "exodus_exports_use_current_blocks_sets_int64_topology_map_owner_mesh_and_file": (
            _exodus_semantic_export_identity_ok(
                summary.get("exodus_block_sideset_nodeset_int64_topology_element_map_owner_mesh_export_generation_identity")
            )
        ),
        "imprint_merge_replays_use_current_tolerance_topology_lineage_order_checkpoint_model_database_and_result": (
            _imprint_merge_contract_identity_ok(
                summary.get("imprint_merge_tolerance_topology_count_entity_lineage_command_checkpoint_model_final_database_result_generation_identity")
            )
        ),
        "headless_batches_use_current_exit_license_fallback_journal_log_database_command_process_and_result": (
            _headless_batch_contract_identity_ok(
                summary.get("headless_batch_exit_license_fallback_journal_log_database_command_process_result_generation_identity")
            )
        ),
        "virtual_geometry_replays_use_current_suppression_topology_inheritance_quality_undo_checkpoint_database_and_result": (
            _virtual_geometry_contract_identity_ok(
                summary.get("virtual_geometry_suppression_topology_map_inheritance_quality_undo_checkpoint_database_result_generation_identity")
            )
        ),
        "anisotropic_crack_regions_use_spd_metric_eigenframe_sizing_alignment_cohesive_hex_quality_and_result": (
            _anisotropic_crack_contract_identity_ok(
                summary.get("anisotropic_crack_metric_eigenframe_sizing_alignment_cohesive_hexconformity_quality_region_export_result_generation_identity")
            )
        ),
        "journal_recreate_replays_use_current_entities_memberships_undo_checkpoint_database_and_result": (
            _journal_recreate_replay_identity_ok(
                summary.get("journal_delete_recreate_entity_group_block_sideset_undo_checkpoint_database_result_generation_identity")
            )
        ),
        "adaptive_mesh_replays_use_current_size_field_projection_region_quality_blocks_session_and_export": (
            _adaptive_size_field_replay_identity_ok(
                summary.get("adaptive_size_field_node_projection_region_boundary_quality_block_session_export_generation_identity")
            )
        ),
        "smoothing_replays_use_current_algorithm_iterations_constraints_motion_quality_checkpoint_database_and_result": (
            _smoothing_replay_identity_ok(
                summary.get("smoothing_iteration_constraint_node_motion_quality_history_checkpoint_database_result_generation_identity")
            )
        ),
        "named_entity_replays_use_current_names_metadata_groups_sets_save_open_owner_database_and_export": (
            _named_entity_transfer_identity_ok(
                summary.get("named_entity_metadata_group_transfer_save_open_export_owner_generation_identity")
            )
        ),
        "journal_undo_redo_replays_use_current_transactions_entities_selection_database_owner_and_result": (
            _journal_undo_redo_identity_ok(
                summary.get("journal_undo_redo_transaction_entityid_selection_database_owner_generation_identity")
            )
        ),
        "exodus_transients_use_current_times_variables_truth_sets_mesh_database_and_result": (
            _exodus_transient_identity_ok(
                summary.get("exodus_transient_timestep_nodal_element_variable_truth_sideset_mesh_database_result_generation_identity")
            )
        ),
        "python_batches_rollback_current_exception_entities_model_database_session_log_and_result": (
            _python_batch_rollback_identity_ok(
                summary.get("python_batch_exception_rollback_entity_generation_database_session_commandlog_result_generation_identity")
            )
        ),
        "mesh_recipes_use_current_dependency_dag_parameters_order_schemes_sets_export_owner_and_result": (
            _mesh_recipe_dag_identity_ok(
                summary.get("meshrecipe_dependency_dag_parameter_execution_scheme_set_export_owner_result_generation_identity")
            )
        ),
        "sidesets_use_current_normals_merge_mesh_blocks_database_owner_export_and_result": (
            _sideset_normal_propagation_identity_ok(
                summary.get("sideset_normal_propagation_merge_mesh_export_entity_owner_result_generation_identity")
            )
        ),
        "cad_imports_use_current_units_healing_bodies_topology_entities_database_owner_and_result": (
            _cad_import_healing_identity_ok(
                summary.get("cad_import_healing_tolerance_units_bodycount_topology_database_owner_result_generation_identity")
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
