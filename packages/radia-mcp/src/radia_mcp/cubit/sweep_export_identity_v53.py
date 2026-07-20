"""Swept-hex, pyramid transition, Exodus sideset, and journal identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .topology_transaction_identity_v54 import (
    validate_public_identity as validate_public_v54_identity,
    validate_source_identity as validate_source_v54_identity,
)


HEX_SWEEP = "hex_sweep_source_target_face_layer_orientation_volume_owner_identity"
PYRAMID = "pyramid_transition_apex_basequad_conformity_jacobian_owner_identity"
EXODUS = "exodus_sideset_element_side_ordinal_topology_block_owner_identity"
JOURNAL = "journal_include_aprepro_expansion_workdir_owner_identity"
_SIDE_COUNTS = {"HEX8": 6, "TET4": 4, "PYRAMID5": 5}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _hex_sweep_ok(row: Mapping[str, object]) -> bool:
    source = str(row.get("source_surface") or "")
    target = str(row.get("target_surface") or "")
    source_count = row.get("source_quad_count")
    target_count = row.get("target_quad_count")
    layers = row.get("layer_count")
    return (
        _generation(row, "face_generation", "layer_generation", "orientation_generation", "volume_generation", "owner_generation", "result_generation")
        and source.startswith("surface:") and target.startswith("surface:") and source != target
        and isinstance(source_count, int) and not isinstance(source_count, bool) and source_count > 0
        and target_count == source_count
        and isinstance(layers, int) and not isinstance(layers, bool) and layers >= 1
        and row.get("orientation") == "source_to_target"
        and str(row.get("volume_id") or "").startswith("volume:")
        and all(row.get("result_" + name) == row.get(name) for name in ("source_surface", "target_surface", "source_quad_count", "target_quad_count", "layer_count", "orientation", "volume_id"))
        and str(row.get("volume_owner") or "").startswith("headless:")
        and row.get("result_volume_owner") == row.get("volume_owner") and _result(row)
    )


def _pyramid_ok(row: Mapping[str, object]) -> bool:
    apex = row.get("apex_node")
    base = row.get("base_quad_nodes")
    neighbors = row.get("conformal_neighbor_faces")
    jacobian = row.get("minimum_scaled_jacobian")
    return (
        _generation(row, "apex_generation", "base_generation", "conformity_generation", "jacobian_generation", "owner_generation", "result_generation")
        and isinstance(apex, int) and not isinstance(apex, bool) and apex > 0
        and isinstance(base, Sequence) and not isinstance(base, (str, bytes)) and len(base) == 4 and len(set(base)) == 4
        and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 and node != apex for node in base)
        and isinstance(neighbors, int) and not isinstance(neighbors, bool) and neighbors == 5
        and isinstance(jacobian, (int, float)) and not isinstance(jacobian, bool) and math.isfinite(float(jacobian)) and float(jacobian) > 0.0
        and all(row.get("result_" + name) == row.get(name) for name in ("apex_node", "base_quad_nodes", "conformal_neighbor_faces", "minimum_scaled_jacobian"))
        and str(row.get("transition_owner") or "").startswith("headless:")
        and row.get("result_transition_owner") == row.get("transition_owner") and _result(row)
    )


def _exodus_ok(row: Mapping[str, object]) -> bool:
    entries = row.get("sideset_entries")
    entries_ok = isinstance(entries, Sequence) and not isinstance(entries, (str, bytes)) and bool(entries)
    if entries_ok:
        seen: set[tuple[int, int]] = set()
        for entry in entries:
            if not isinstance(entry, Mapping) or set(entry) != {"element_id", "side_ordinal", "topology", "block_id"}:
                entries_ok = False; break
            topology = entry["topology"]
            key = (entry["element_id"], entry["side_ordinal"])
            if topology not in _SIDE_COUNTS or not all(isinstance(entry[name], int) and not isinstance(entry[name], bool) and entry[name] >= 1 for name in ("element_id", "side_ordinal", "block_id")) or entry["side_ordinal"] > _SIDE_COUNTS[topology] or key in seen:
                entries_ok = False; break
            seen.add(key)
    return (
        _generation(row, "sideset_generation", "element_generation", "ordinal_generation", "topology_generation", "block_generation", "owner_generation", "result_generation")
        and bool(str(row.get("sideset_name") or "")) and entries_ok
        and row.get("replayed_sideset_name") == row.get("sideset_name")
        and row.get("replayed_sideset_entries") == entries
        and str(row.get("export_owner") or "").startswith("headless:")
        and row.get("replayed_export_owner") == row.get("export_owner") and _result(row)
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    includes = row.get("include_order")
    hashes = row.get("include_sha256")
    variables = row.get("aprepro_variables")
    expanded = row.get("expanded_variables")
    workdir = str(row.get("working_directory") or "")
    include_ok = isinstance(includes, Sequence) and not isinstance(includes, (str, bytes)) and bool(includes) and len(includes) == len(set(includes)) and all(isinstance(path, str) and path and "/" not in path and "\\" not in path for path in includes)
    return (
        _generation(row, "include_generation", "variable_generation", "expansion_generation", "workdir_generation", "owner_generation", "result_generation")
        and include_ok and isinstance(hashes, Mapping) and set(hashes) == set(includes) and all(_digest(value) for value in hashes.values())
        and isinstance(variables, Mapping) and bool(variables) and isinstance(expanded, Mapping) and bool(expanded)
        and workdir.startswith("workspace:coreform/") and ".." not in workdir.replace("\\", "/").split("/")
        and row.get("replayed_include_order") == includes and row.get("replayed_include_sha256") == hashes
        and row.get("replayed_aprepro_variables") == variables and row.get("replayed_expanded_variables") == expanded
        and row.get("replayed_working_directory") == workdir
        and str(row.get("journal_owner") or "").startswith("headless:")
        and row.get("replayed_journal_owner") == row.get("journal_owner") and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping): return {}
    checks: dict[str, bool] = {}
    v54 = validate_public_v54_identity(payload)
    if v54: checks.update(v54["checks"])
    if payload.get(HEX_SWEEP) is not None: checks["v53_hex_sweep_face_layer_orientation_owner"] = isinstance(payload[HEX_SWEEP], Mapping) and _hex_sweep_ok(payload[HEX_SWEEP])
    if payload.get(PYRAMID) is not None: checks["v53_pyramid_apex_base_conformity_jacobian_owner"] = isinstance(payload[PYRAMID], Mapping) and _pyramid_ok(payload[PYRAMID])
    return _report("cubit_v53_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping): return {}
    checks: dict[str, bool] = {}
    v54 = validate_source_v54_identity(payload)
    if v54: checks.update(v54["checks"])
    if payload.get(EXODUS) is not None: checks["v53_exodus_sideset_element_side_topology_block_owner"] = isinstance(payload[EXODUS], Mapping) and _exodus_ok(payload[EXODUS])
    if payload.get(JOURNAL) is not None: checks["v53_journal_include_aprepro_workdir_owner"] = isinstance(payload[JOURNAL], Mapping) and _journal_ok(payload[JOURNAL])
    return _report("cubit_v53_source_identity_v1", checks) if checks else {}
