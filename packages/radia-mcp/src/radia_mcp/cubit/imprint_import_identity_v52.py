"""Imprint, scheme, ACIS-import, and Aprepro identity checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


SHEET = "sheet_imprint_merge_tolerance_curvesplit_topology_owner_identity"
SCHEME = "volume_scheme_autosmooth_curveinterval_seed_owner_identity"
ACIS = "acis_import_tolerance_healing_bodyname_layer_owner_identity"
APREPRO = "aprepro_scope_include_expression_unit_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _named_mapping(value: object, *, prefix: str) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(isinstance(name, str) and name.startswith(prefix) for name in value)


def _sheet_ok(row: Mapping[str, object]) -> bool:
    split_map = row.get("curve_split_map")
    topology = row.get("topology_counts")
    tolerance = row.get("merge_tolerance")
    return (
        _generation(row, "tolerance_generation", "split_generation", "topology_generation", "owner_generation", "result_generation")
        and isinstance(tolerance, (int, float))
        and not isinstance(tolerance, bool)
        and math.isfinite(float(tolerance))
        and 0.0 < float(tolerance) <= 1.0e-3
        and row.get("result_merge_tolerance") == tolerance
        and _named_mapping(split_map, prefix="curve:")
        and all(
            isinstance(children, Sequence)
            and not isinstance(children, (str, bytes))
            and len(children) >= 2
            and len(children) == len(set(children))
            and all(isinstance(child, str) and child.startswith("curve:") for child in children)
            for children in split_map.values()
        )
        and row.get("result_curve_split_map") == split_map
        and isinstance(topology, Mapping)
        and set(topology) == {"sheets", "surfaces", "curves", "vertices"}
        and all(isinstance(count, int) and not isinstance(count, bool) and count > 0 for count in topology.values())
        and topology["curves"] >= sum(len(children) for children in split_map.values())
        and row.get("result_topology_counts") == topology
        and str(row.get("body_owner") or "").startswith("headless:")
        and row.get("result_body_owner") == row.get("body_owner")
        and _result(row)
    )


def _scheme_ok(row: Mapping[str, object]) -> bool:
    schemes = row.get("volume_schemes")
    intervals = row.get("curve_intervals")
    seeds = row.get("volume_seeds")
    return (
        _generation(row, "scheme_generation", "smooth_generation", "interval_generation", "seed_generation", "owner_generation", "result_generation")
        and _named_mapping(schemes, prefix="volume:")
        and all(scheme in {"map", "sweep", "submap"} for scheme in schemes.values())
        and row.get("result_volume_schemes") == schemes
        and row.get("autosmooth") is True
        and row.get("result_autosmooth") is True
        and _named_mapping(intervals, prefix="curve:")
        and all(isinstance(count, int) and not isinstance(count, bool) and count >= 2 for count in intervals.values())
        and row.get("result_curve_intervals") == intervals
        and isinstance(seeds, Mapping)
        and set(seeds) == set(schemes)
        and all(isinstance(seed, int) and not isinstance(seed, bool) and seed > 0 for seed in seeds.values())
        and row.get("result_volume_seeds") == seeds
        and str(row.get("volume_owner") or "").startswith("headless:")
        and row.get("result_volume_owner") == row.get("volume_owner")
        and _result(row)
    )


def _acis_ok(row: Mapping[str, object]) -> bool:
    healing = row.get("healing_operations")
    names = row.get("body_names")
    layers = row.get("body_layers")
    tolerance = row.get("import_tolerance")
    return (
        _generation(row, "tolerance_generation", "healing_generation", "body_generation", "layer_generation", "owner_generation", "result_generation")
        and isinstance(tolerance, (int, float))
        and not isinstance(tolerance, bool)
        and math.isfinite(float(tolerance))
        and 0.0 < float(tolerance) <= 1.0e-3
        and row.get("replayed_import_tolerance") == tolerance
        and isinstance(healing, Sequence)
        and not isinstance(healing, (str, bytes))
        and bool(healing)
        and len(healing) == len(set(healing))
        and all(operation in {"stitch", "remove_sliver", "merge_coincident"} for operation in healing)
        and row.get("replayed_healing_operations") == healing
        and isinstance(names, Sequence)
        and not isinstance(names, (str, bytes))
        and bool(names)
        and len(names) == len(set(names))
        and all(isinstance(name, str) and name for name in names)
        and row.get("replayed_body_names") == names
        and isinstance(layers, Mapping)
        and set(layers) == set(names)
        and all(isinstance(layer, str) and layer for layer in layers.values())
        and row.get("replayed_body_layers") == layers
        and str(row.get("geometry_owner") or "").startswith("headless:")
        and row.get("replayed_geometry_owner") == row.get("geometry_owner")
        and _result(row)
    )


def _aprepro_ok(row: Mapping[str, object]) -> bool:
    scopes = row.get("variable_scopes")
    includes = row.get("include_order")
    expressions = row.get("expressions")
    try:
        gap = float(scopes["global"]["gap"]["value"])
        gap_unit = scopes["global"]["gap"]["unit"]
        expression = expressions["airgap_total"]
        resolved = float(expression["resolved"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        _generation(row, "scope_generation", "include_generation", "expression_generation", "unit_generation", "owner_generation", "result_generation")
        and isinstance(scopes, Mapping)
        and "global" in scopes
        and gap_unit in {"m", "mm", "um"}
        and math.isfinite(gap)
        and gap > 0.0
        and row.get("replayed_variable_scopes") == scopes
        and isinstance(includes, Sequence)
        and not isinstance(includes, (str, bytes))
        and len(includes) >= 2
        and len(includes) == len(set(includes))
        and all(isinstance(path, str) and path.endswith(".apr") and "/" not in path and "\\" not in path for path in includes)
        and row.get("replayed_include_order") == includes
        and isinstance(expressions, Mapping)
        and expression.get("expression") == "2*gap"
        and expression.get("unit") == gap_unit
        and math.isclose(resolved, 2.0 * gap, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and row.get("replayed_expressions") == expressions
        and str(row.get("journal_owner") or "").startswith("headless:")
        and row.get("replayed_journal_owner") == row.get("journal_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(SHEET) is not None:
        checks["v52_sheet_imprint_tolerance_split_topology_owner"] = isinstance(payload[SHEET], Mapping) and _sheet_ok(payload[SHEET])
    if payload.get(SCHEME) is not None:
        checks["v52_volume_scheme_smooth_interval_seed_owner"] = isinstance(payload[SCHEME], Mapping) and _scheme_ok(payload[SCHEME])
    return _report("cubit_v52_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(ACIS) is not None:
        checks["v52_acis_tolerance_healing_body_layer_owner"] = isinstance(payload[ACIS], Mapping) and _acis_ok(payload[ACIS])
    if payload.get(APREPRO) is not None:
        checks["v52_aprepro_scope_include_expression_unit_owner"] = isinstance(payload[APREPRO], Mapping) and _aprepro_ok(payload[APREPRO])
    return _report("cubit_v52_source_identity_v1", checks) if checks else {}
