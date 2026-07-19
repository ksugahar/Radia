"""Semantic identity checks for configured CAD builds and source replays."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


SUPPRESSION = "feature_suppression_configuration_mass_cache_owner_identity"
LOFT = "loft_profile_orientation_wire_seam_correspondence_self_intersection_owner_identity"
IMPORT = "import_unit_inference_layer_color_subshape_mapping_owner_identity"
SELECTOR = "topology_selector_query_cardinality_witness_feature_history_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _result_ok(row: Mapping[str, object], *, owner: bool = False) -> bool:
    return (
        (not owner or bool(str(row.get("owner") or "")) and row.get("accepted_owner") == row.get("owner"))
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _generations_ok(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _finite_vector(value: object, *, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _suppression_ok(row: Mapping[str, object]) -> bool:
    suppressed = row.get("suppressed_features")
    configuration = row.get("configuration")
    mass = row.get("mass_kg")
    owner = str(row.get("owner") or "")
    return (
        _generations_ok(
            row,
            (
                "suppression_generation",
                "configuration_generation",
                "shape_generation",
                "cache_generation",
                "result_generation",
            ),
        )
        and isinstance(suppressed, list)
        and len(set(suppressed)) == len(suppressed)
        and all(bool(str(feature).strip()) for feature in suppressed)
        and row.get("result_suppressed_features") == suppressed
        and isinstance(configuration, Mapping)
        and bool(configuration)
        and row.get("result_configuration") == configuration
        and _digest(row.get("shape_sha256"))
        and row.get("cache_shape_sha256") == row.get("shape_sha256")
        and isinstance(mass, (int, float))
        and math.isfinite(float(mass))
        and float(mass) >= 0.0
        and row.get("cached_mass_kg") == mass
        and row.get("result_mass_kg") == mass
        and bool(owner)
        and row.get("cache_owner") == owner
        and _result_ok(row, owner=True)
    )


def _loft_ok(row: Mapping[str, object]) -> bool:
    profiles = row.get("profile_order")
    orientation = row.get("profile_orientation")
    seams = row.get("wire_seams")
    correspondence = row.get("profile_correspondence")
    valid_profiles = isinstance(profiles, list) and bool(profiles) and len(set(profiles)) == len(profiles)
    valid_correspondence = (
        isinstance(correspondence, list)
        and bool(correspondence)
        and all(isinstance(chain, list) and len(chain) == len(profiles) and len(set(chain)) == len(chain) for chain in correspondence)
    ) if valid_profiles else False
    return (
        _generations_ok(
            row,
            (
                "profile_generation",
                "orientation_generation",
                "seam_generation",
                "correspondence_generation",
                "diagnostic_generation",
                "result_generation",
            ),
        )
        and valid_profiles
        and row.get("result_profile_order") == profiles
        and isinstance(orientation, Mapping)
        and set(orientation) == set(profiles)
        and all(value in {-1, 1} for value in orientation.values())
        and row.get("result_profile_orientation") == orientation
        and isinstance(seams, Mapping)
        and set(seams) == set(profiles)
        and all(bool(str(value).strip()) for value in seams.values())
        and row.get("result_wire_seams") == seams
        and valid_correspondence
        and row.get("result_profile_correspondence") == correspondence
        and row.get("self_intersection_count") == row.get("result_self_intersection_count") == 0
        and _result_ok(row, owner=True)
    )


def _import_ok(row: Mapping[str, object]) -> bool:
    unit = str(row.get("inferred_unit") or "")
    scale = row.get("unit_scale_to_m")
    layers = row.get("layer_map")
    colors = row.get("color_map")
    subshapes = row.get("persistent_subshape_map")
    revision = str(row.get("source_revision") or "")
    valid_metadata = (
        isinstance(layers, Mapping)
        and bool(layers)
        and isinstance(colors, Mapping)
        and set(colors) == set(layers)
        and all(_finite_vector(color, length=3) and all(0.0 <= float(channel) <= 1.0 for channel in color) for color in colors.values())
    )
    return (
        _generations_ok(row, ("unit_generation", "metadata_generation", "subshape_generation", "result_generation"))
        and unit in {"m", "mm", "cm", "in"}
        and row.get("result_inferred_unit") == unit
        and isinstance(scale, (int, float))
        and math.isfinite(float(scale))
        and float(scale) > 0.0
        and row.get("result_unit_scale_to_m") == scale
        and valid_metadata
        and row.get("result_layer_map") == layers
        and row.get("result_color_map") == colors
        and isinstance(subshapes, Mapping)
        and bool(subshapes)
        and len(set(subshapes.values())) == len(subshapes)
        and row.get("result_persistent_subshape_map") == subshapes
        and bool(revision)
        and row.get("result_source_revision") == revision
        and row.get("owner") == f"import:{revision}"
        and _result_ok(row, owner=True)
    )


def _selector_ok(row: Mapping[str, object]) -> bool:
    query = str(row.get("query") or "")
    cardinality = row.get("expected_cardinality")
    entities = row.get("selected_entities")
    witnesses = row.get("witness_points")
    history_owner = str(row.get("feature_history_owner") or "")
    return (
        _generations_ok(row, ("query_generation", "witness_generation", "history_generation", "result_generation"))
        and bool(query)
        and row.get("result_query") == query
        and isinstance(cardinality, int)
        and cardinality >= 0
        and isinstance(entities, list)
        and len(entities) == cardinality
        and len(set(entities)) == len(entities)
        and row.get("result_cardinality") == cardinality
        and row.get("result_selected_entities") == entities
        and isinstance(witnesses, Mapping)
        and set(witnesses) == set(entities)
        and all(_finite_vector(point, length=3) for point in witnesses.values())
        and row.get("result_witness_points") == witnesses
        and history_owner.startswith("history:")
        and row.get("result_feature_history_owner") == history_owner
        and _result_ok(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows: list[Mapping[str, object]] = []
    if isinstance(payload.get("reference"), list):
        rows.extend(row for row in payload["reference"] if isinstance(row, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, Mapping))
    checks: dict[str, bool] = {}
    suppression = [row.get(SUPPRESSION) for row in rows if SUPPRESSION in row]
    loft = [row.get(LOFT) for row in rows if LOFT in row]
    if suppression:
        checks["v48_feature_configuration_mass_cache_owner"] = len(suppression) == len(rows) and all(
            isinstance(item, Mapping) and _suppression_ok(item) for item in suppression
        )
    if loft:
        checks["v48_loft_profile_seam_correspondence_owner"] = len(loft) == len(rows) and all(
            isinstance(item, Mapping) and _loft_ok(item) for item in loft
        )
    return _report("build123d_v48_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    imported = identity.get(IMPORT)
    selector = identity.get(SELECTOR)
    if imported is not None:
        checks["v48_import_unit_metadata_subshape_owner"] = isinstance(imported, Mapping) and _import_ok(imported)
    if selector is not None:
        checks["v48_selector_cardinality_witness_history"] = isinstance(selector, Mapping) and _selector_ok(selector)
    return _report("build123d_v48_source_identity_v1", checks) if checks else {}
