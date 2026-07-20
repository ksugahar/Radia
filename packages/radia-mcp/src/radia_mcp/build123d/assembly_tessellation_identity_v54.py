"""Assembly, loft, AP242, and tessellation identity checks for v54."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .thread_sheet_identity_v55 import (
    validate_public_identity as validate_public_v55_identity,
    validate_source_identity as validate_source_v55_identity,
)


ASSEMBLY = "assembly_massproperty_density_location_center_inertia_owner_identity"
LOFT = "loft_section_orientation_parameter_seam_topology_owner_identity"
STEP = "step_ap242_unit_productstructure_color_layer_owner_identity"
TESSELLATION = "tessellation_deflection_angle_orientation_index_owner_identity"
_LENGTH_UNITS = {"m", "mm", "cm", "inch"}


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _vector(value: object, length: int) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == length and all(_finite(item) for item in value)


def _inertia_tensor(value: object) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3 or not all(_vector(row, 3) for row in value):
        return False
    matrix = [[float(item) for item in row] for row in value]
    if not all(math.isclose(matrix[i][j], matrix[j][i], rel_tol=0.0, abs_tol=1.0e-12) for i in range(3) for j in range(3)):
        return False
    a, b, c = matrix[0]
    _, d, e = matrix[1]
    _, _, f = matrix[2]
    determinant = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    return (
        a > 0.0
        and a * d - b * b > 0.0
        and determinant > 0.0
        and a + d >= f - 1.0e-12
        and a + f >= d - 1.0e-12
        and d + f >= a - 1.0e-12
    )


def _location(value: object) -> bool:
    if not isinstance(value, Mapping) or set(value) != {"translation_m", "quaternion_wxyz"}:
        return False
    quaternion = value["quaternion_wxyz"]
    return (
        _vector(value["translation_m"], 3)
        and _vector(quaternion, 4)
        and math.isclose(sum(float(item) ** 2 for item in quaternion), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    )


def _assembly_ok(row: Mapping[str, object]) -> bool:
    densities = row.get("solid_densities_kg_m3")
    locations = row.get("located_solids")
    records_ok = (
        isinstance(densities, Mapping)
        and bool(densities)
        and all(isinstance(name, str) and name.startswith("solid:") and _finite(value) and float(value) > 0.0 for name, value in densities.items())
        and isinstance(locations, Mapping)
        and set(locations) == set(densities)
        and all(_location(value) for value in locations.values())
    )
    return (
        _generation(row, "density_generation", "location_generation", "center_generation", "inertia_generation", "owner_generation", "result_generation")
        and records_ok
        and row.get("result_solid_densities_kg_m3") == densities
        and row.get("result_located_solids") == locations
        and _vector(row.get("center_of_mass_m"), 3)
        and row.get("result_center_of_mass_m") == row.get("center_of_mass_m")
        and _inertia_tensor(row.get("inertia_tensor_kg_m2"))
        and row.get("result_inertia_tensor_kg_m2") == row.get("inertia_tensor_kg_m2")
        and str(row.get("assembly_owner") or "").startswith("assembly:")
        and row.get("result_assembly_owner") == row.get("assembly_owner")
        and _result(row)
    )


def _loft_ok(row: Mapping[str, object]) -> bool:
    sections = row.get("section_correspondence")
    topology = row.get("resulting_topology")
    sections_ok = isinstance(sections, Sequence) and not isinstance(sections, (str, bytes)) and len(sections) >= 2
    if sections_ok:
        section_names: list[str] = []
        seams: list[str] = []
        for section in sections:
            if not isinstance(section, Mapping) or set(section) != {"section", "orientation", "parameter_start", "seam"}:
                sections_ok = False
                break
            name = section["section"]
            seam = section["seam"]
            parameter = section["parameter_start"]
            if not (
                isinstance(name, str)
                and name.startswith("wire:")
                and isinstance(section["orientation"], int)
                and not isinstance(section["orientation"], bool)
                and section["orientation"] in {-1, 1}
                and _finite(parameter)
                and 0.0 <= float(parameter) < 1.0
                and isinstance(seam, str)
                and seam.startswith("vertex:")
            ):
                sections_ok = False
                break
            section_names.append(name)
            seams.append(seam)
        sections_ok = sections_ok and len(section_names) == len(set(section_names)) and len(seams) == len(set(seams))
    topology_ok = isinstance(topology, Mapping) and set(topology) == {"solids", "shells", "faces", "edges", "vertices"}
    if topology_ok:
        topology_ok = (
            all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in topology.values())
            and topology["solids"] == 1
            and topology["shells"] == 1
            and topology["faces"] >= 3
            and topology["edges"] >= 3
            and topology["vertices"] >= 2
            and topology["vertices"] - topology["edges"] + topology["faces"] == 2
        )
    return (
        _generation(row, "section_generation", "parameter_generation", "seam_generation", "topology_generation", "owner_generation", "result_generation")
        and sections_ok
        and row.get("result_section_correspondence") == sections
        and topology_ok
        and row.get("result_resulting_topology") == topology
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _rgb(value: object) -> bool:
    return _vector(value, 3) and all(0.0 <= float(channel) <= 1.0 for channel in value)


def _product_structure(value: object) -> tuple[bool, set[str]]:
    if not isinstance(value, Mapping) or not value:
        return False, set()
    nodes: set[str] = set()
    for parent, children in value.items():
        if not isinstance(parent, str) or not parent.startswith(("product:", "part:")):
            return False, set()
        if not isinstance(children, Sequence) or isinstance(children, (str, bytes)) or not children:
            return False, set()
        if not all(isinstance(child, str) and child.startswith(("part:", "solid:")) for child in children):
            return False, set()
        if len(children) != len(set(children)) or parent in children:
            return False, set()
        nodes.add(parent)
        nodes.update(children)
    return True, nodes


def _step_ok(row: Mapping[str, object]) -> bool:
    structure = row.get("product_structure")
    structure_ok, nodes = _product_structure(structure)
    colors = row.get("component_colors")
    layers = row.get("component_layers")
    metadata_ok = (
        isinstance(colors, Mapping)
        and bool(colors)
        and all(isinstance(part, str) and part.startswith("part:") and part in nodes and _rgb(color) for part, color in colors.items())
        and isinstance(layers, Mapping)
        and set(layers) == set(colors)
        and all(isinstance(layer, str) and layer.startswith("layer:") for layer in layers.values())
    )
    return (
        _generation(row, "unit_generation", "structure_generation", "color_generation", "layer_generation", "owner_generation", "result_generation")
        and row.get("schema") == "AP242"
        and row.get("replayed_schema") == "AP242"
        and row.get("length_unit") in _LENGTH_UNITS
        and row.get("replayed_length_unit") == row.get("length_unit")
        and structure_ok
        and row.get("replayed_product_structure") == structure
        and metadata_ok
        and row.get("replayed_component_colors") == colors
        and row.get("replayed_component_layers") == layers
        and str(row.get("document_owner") or "").startswith("document:")
        and row.get("replayed_document_owner") == row.get("document_owner")
        and _result(row)
    )


def _tessellation_ok(row: Mapping[str, object]) -> bool:
    vertices = row.get("vertices_m")
    triangles = row.get("triangle_indices")
    orientations = row.get("triangle_orientations")
    vertices_ok = isinstance(vertices, Sequence) and not isinstance(vertices, (str, bytes)) and len(vertices) >= 3 and all(_vector(vertex, 3) for vertex in vertices)
    triangles_ok = isinstance(triangles, Sequence) and not isinstance(triangles, (str, bytes)) and bool(triangles) and vertices_ok
    canonical: list[tuple[int, int, int]] = []
    if triangles_ok:
        for triangle in triangles:
            if not (
                isinstance(triangle, Sequence)
                and not isinstance(triangle, (str, bytes))
                and len(triangle) == 3
                and all(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(vertices) for index in triangle)
                and len(set(triangle)) == 3
            ):
                triangles_ok = False
                break
            canonical.append(tuple(sorted(triangle)))
        triangles_ok = triangles_ok and len(canonical) == len(set(canonical))
    return (
        _generation(row, "deflection_generation", "angle_generation", "orientation_generation", "index_generation", "owner_generation", "result_generation")
        and _finite(row.get("linear_deflection_m"))
        and float(row["linear_deflection_m"]) > 0.0
        and row.get("replayed_linear_deflection_m") == row.get("linear_deflection_m")
        and _finite(row.get("angular_deflection_rad"))
        and 0.0 < float(row["angular_deflection_rad"]) <= math.pi
        and row.get("replayed_angular_deflection_rad") == row.get("angular_deflection_rad")
        and vertices_ok
        and row.get("replayed_vertices_m") == vertices
        and triangles_ok
        and row.get("replayed_triangle_indices") == triangles
        and isinstance(orientations, Sequence)
        and not isinstance(orientations, (str, bytes))
        and len(orientations) == len(triangles)
        and all(isinstance(orientation, int) and not isinstance(orientation, bool) and orientation in {-1, 1} for orientation in orientations)
        and row.get("replayed_triangle_orientations") == orientations
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("replayed_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    rows: list[Mapping[str, object]] = []
    reference = payload.get("reference")
    if isinstance(reference, Sequence) and not isinstance(reference, (str, bytes)):
        rows.extend(item for item in reference if isinstance(item, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for family in measured.values():
            if isinstance(family, Sequence) and not isinstance(family, (str, bytes)):
                rows.extend(item for item in family if isinstance(item, Mapping))
    return rows


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    checks: dict[str, bool] = {}
    v55 = validate_public_v55_identity(payload)
    if v55:
        checks.update(v55["checks"])
    assemblies = [row.get(ASSEMBLY) for row in rows if ASSEMBLY in row]
    lofts = [row.get(LOFT) for row in rows if LOFT in row]
    if assemblies:
        checks["v54_assembly_density_location_center_inertia_owner"] = len(assemblies) == len(rows) and all(isinstance(item, Mapping) and _assembly_ok(item) for item in assemblies)
    if lofts:
        checks["v54_loft_section_parameter_seam_topology_owner"] = len(lofts) == len(rows) and all(isinstance(item, Mapping) and _loft_ok(item) for item in lofts)
    return _report("build123d_v54_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    v55 = validate_source_v55_identity(payload)
    if v55:
        checks.update(v55["checks"])
    if identity.get(STEP) is not None:
        checks["v54_step_ap242_unit_structure_color_layer_owner"] = isinstance(identity[STEP], Mapping) and _step_ok(identity[STEP])
    if identity.get(TESSELLATION) is not None:
        checks["v54_tessellation_deflection_angle_orientation_index_owner"] = isinstance(identity[TESSELLATION], Mapping) and _tessellation_ok(identity[TESSELLATION])
    return _report("build123d_v54_source_identity_v1", checks) if checks else {}
