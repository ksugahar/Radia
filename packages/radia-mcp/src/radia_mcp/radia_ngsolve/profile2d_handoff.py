"""Closed-world gate for portable 2-D CAD, mesh, and solver semantics."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any


PROFILE_SCHEMA = "radia.profile2d-handoff.v1"
RESULT_SCHEMA = "radia.profile2d-handoff-gate.v1"
ABI_SCHEMA = "radia.fixed-scalar-io.v1"
_HEX = re.compile(r"[0-9a-f]{64}")
_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,63}")
_TRANSPORTS = {"matlab_mex", "simulink_s_function", "python_batch"}


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha(value: Any) -> str:
    payload = value if isinstance(value, bytes) else _canonical(value).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be finite")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _identifier(value: Any, label: str, *, allow_empty: bool = False) -> str:
    result = str(value or "")
    if allow_empty and not result:
        return result
    if not _NAME.fullmatch(result):
        raise ValueError(f"{label} must match {_NAME.pattern}")
    return result


def _id_list(value: Any, label: str) -> list[str]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be an array")
    result = [_identifier(item, f"{label}[]") for item in value]
    if len(result) != len(set(result)):
        raise ValueError(f"{label} contains duplicates")
    return result


def _edge_semantics(
    row: Mapping[str, Any],
    label: str,
    boundary_ids: set[str],
    conductor_ids: set[str],
) -> tuple[str, str]:
    boundary = _identifier(row.get("boundary_id"), f"{label}.boundary_id", allow_empty=True)
    conductor = _identifier(row.get("conductor_id"), f"{label}.conductor_id", allow_empty=True)
    if boundary and boundary not in boundary_ids:
        raise ValueError(f"{label} references unknown boundary_id {boundary!r}")
    if conductor and conductor not in conductor_ids:
        raise ValueError(f"{label} references unknown conductor_id {conductor!r}")
    return boundary, conductor


def _arc_metrics(p0: tuple[float, float], p1: tuple[float, float], sweep_deg: float) -> tuple[float, float]:
    theta = math.radians(sweep_deg)
    chord_x = p1[0] - p0[0]
    chord_y = p1[1] - p0[1]
    chord = math.hypot(chord_x, chord_y)
    if chord <= 0.0:
        raise ValueError("arc endpoints must be distinct")
    if abs(theta) <= 1.0e-12 or abs(theta) > math.pi + 1.0e-12:
        raise ValueError("arc sweep_deg must satisfy 0 < abs(sweep) <= 180")
    radius = chord / (2.0 * math.sin(abs(theta) / 2.0))
    midpoint = ((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0)
    left = (-chord_y / chord, chord_x / chord)
    center_offset = chord / (2.0 * math.tan(theta / 2.0))
    center = (
        midpoint[0] + left[0] * center_offset,
        midpoint[1] + left[1] * center_offset,
    )
    area = 0.5 * (
        center[0] * (p1[1] - p0[1])
        - center[1] * (p1[0] - p0[0])
        + radius * radius * theta
    )
    return area, radius * abs(theta)


def _artifact_contract(
    raw: Any,
    semantic_sha256: str,
    profile_sha256: str,
    exact_area_m2: float | None,
    exact_perimeter_m: float | None,
) -> tuple[dict[str, Any], dict[str, bool]]:
    if raw is None:
        raw = {}
    if not isinstance(raw, Mapping):
        raise ValueError("artifacts must be an object")
    allowed = {"step", "vol", "semantic_sidecar", "cad_measurement"}
    extra = sorted(set(map(str, raw)) - allowed)
    if extra:
        raise ValueError(f"unsupported artifact keys: {extra}")
    normalized: dict[str, Any] = {}
    decoded_cad_measurement: dict[str, Any] | None = None
    for name in ("step", "vol", "semantic_sidecar", "cad_measurement"):
        if name not in raw:
            continue
        row = raw[name]
        if not isinstance(row, Mapping):
            raise ValueError(f"artifacts.{name} must be an object")
        content = row.get("content")
        if not isinstance(content, str) or not content:
            raise ValueError(f"artifacts.{name}.content must be nonempty inline text")
        content_sha = hashlib.sha256(content.encode("utf-8")).hexdigest()
        claimed = str(row.get("sha256") or content_sha)
        if not _HEX.fullmatch(claimed) or claimed != content_sha:
            raise ValueError(f"artifacts.{name}.sha256 does not match inline content")
        fmt = str(row.get("format") or "")
        if name == "step":
            if fmt not in {"STEP AP203", "STEP AP214", "STEP AP242"}:
                raise ValueError("STEP format must name AP203, AP214, or AP242")
            stripped = content.strip()
            if not stripped.startswith("ISO-10303-21;") or not stripped.endswith("END-ISO-10303-21;"):
                raise ValueError("STEP content is not an ISO-10303-21 exchange file")
        elif name == "vol":
            if fmt != "Netgen vol dimension 2" or "dimension\n2" not in content.lower():
                raise ValueError("vol content must be a dimension-2 Netgen .vol")
        elif name == "semantic_sidecar":
            if fmt != "radia semantic sidecar v1":
                raise ValueError("semantic sidecar format is invalid")
            try:
                sidecar = json.loads(content)
            except json.JSONDecodeError as exc:
                raise ValueError("semantic sidecar is not JSON") from exc
            if _sha(sidecar) != semantic_sha256:
                raise ValueError("semantic sidecar does not match semantics")
        else:
            if fmt != "radia CAD measurement v1":
                raise ValueError("CAD measurement sidecar format is invalid")
            try:
                measurement = json.loads(content)
            except json.JSONDecodeError as exc:
                raise ValueError("CAD measurement sidecar is not JSON") from exc
            if not isinstance(measurement, dict) or measurement.get("schema") != "radia.cad2d-measurement.v1":
                raise ValueError("CAD measurement sidecar schema is invalid")
            decoded_cad_measurement = measurement
        normalized[name] = {"format": fmt, "sha256": content_sha, "bytes": len(content.encode("utf-8"))}
    checks = {
        "step_has_semantic_sidecar": "step" not in normalized or "semantic_sidecar" in normalized,
        "step_has_cad_measurement": "step" not in normalized or "cad_measurement" in normalized,
        "vol_has_semantic_sidecar": "vol" not in normalized or "semantic_sidecar" in normalized,
    }
    if "step" in normalized and decoded_cad_measurement is not None:
        if exact_area_m2 is None or exact_perimeter_m is None:
            raise ValueError("STEP measurement crosscheck currently requires degree-2 closed loops")
        cad_area = _finite(decoded_cad_measurement.get("area_m2"), "cad_measurement.area_m2")
        cad_perimeter = _finite(decoded_cad_measurement.get("perimeter_m"), "cad_measurement.perimeter_m")
        checks.update({
            "cad_measurement_step_identity": decoded_cad_measurement.get("step_sha256") == normalized["step"]["sha256"],
            "cad_measurement_profile_identity": decoded_cad_measurement.get("profile_sha256") == profile_sha256,
            "cad_area_matches_profile": math.isclose(cad_area, exact_area_m2, rel_tol=1.0e-8, abs_tol=1.0e-12),
            "cad_perimeter_matches_profile": math.isclose(cad_perimeter, exact_perimeter_m, rel_tol=1.0e-8, abs_tol=1.0e-12),
        })
    if not all(checks.values()):
        raise ValueError("CAD/mesh geometry must travel with a matching semantic sidecar")
    return normalized, checks


def _abi_contract(raw: Any) -> tuple[dict[str, Any] | None, dict[str, bool]]:
    if raw is None:
        return None, {"fixed_width": True, "si_units": True, "no_dynamic_paths": True}
    if not isinstance(raw, Mapping):
        raise ValueError("execution_abi must be an object")
    if raw.get("schema") != ABI_SCHEMA:
        raise ValueError(f"execution_abi.schema must be {ABI_SCHEMA}")
    transport = str(raw.get("transport") or "")
    if transport not in _TRANSPORTS:
        raise ValueError(f"unsupported execution_abi.transport {transport!r}")
    inputs = raw.get("inputs")
    outputs = raw.get("outputs")
    if not isinstance(inputs, list) or not inputs or not isinstance(outputs, list) or not outputs:
        raise ValueError("execution_abi inputs and outputs must be nonempty arrays")

    def fields(rows: list[Any], label: str) -> list[dict[str, str]]:
        normalized = []
        names = []
        for index, row in enumerate(rows):
            if not isinstance(row, Mapping):
                raise ValueError(f"{label}[{index}] must be an object")
            name = _identifier(row.get("name"), f"{label}[{index}].name")
            unit = str(row.get("unit") or "")
            if not unit or unit.lower() in {"native", "model_unit", "unknown"}:
                raise ValueError(f"{label}[{index}].unit must be an explicit SI unit or 1")
            normalized.append({"name": name, "unit": unit})
            names.append(name)
        if len(names) != len(set(names)):
            raise ValueError(f"{label} contains duplicate field names")
        return normalized

    normalized_inputs = fields(inputs, "execution_abi.inputs")
    normalized_outputs = fields(outputs, "execution_abi.outputs")
    request_sha = str(raw.get("request_contract_sha256") or "")
    checks = {
        "fixed_width": raw.get("fixed_width") is True
        and int(raw.get("input_width", -1)) == len(normalized_inputs)
        and int(raw.get("output_width", -1)) == len(normalized_outputs),
        "si_units": True,
        "no_dynamic_paths": raw.get("dynamic_paths") is False,
        "request_digest": bool(_HEX.fullmatch(request_sha)),
    }
    if transport == "simulink_s_function":
        checks["positive_sample_period"] = _finite(raw.get("sample_period_s"), "sample_period_s") > 0.0
    if not all(checks.values()):
        raise ValueError(f"execution_abi failed fixed-runtime checks: {checks}")
    contract = {
        "schema": ABI_SCHEMA,
        "transport": transport,
        "inputs": normalized_inputs,
        "outputs": normalized_outputs,
        "input_width": len(normalized_inputs),
        "output_width": len(normalized_outputs),
        "fixed_width": True,
        "dynamic_paths": False,
        "request_contract_sha256": request_sha,
    }
    if transport == "simulink_s_function":
        contract["sample_period_s"] = _finite(raw["sample_period_s"], "sample_period_s")
    contract["contract_sha256"] = _sha(contract)
    return contract, checks


def profile2d_handoff_gate(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Validate a source-neutral 2-D profile and its execution handoff."""

    if not isinstance(packet, Mapping) or packet.get("schema") != PROFILE_SCHEMA:
        raise ValueError(f"packet.schema must be {PROFILE_SCHEMA}")
    if packet.get("length_unit") != "m":
        raise ValueError("profile coordinates must use SI metres")
    semantics = packet.get("semantics")
    profile = packet.get("profile")
    if not isinstance(semantics, Mapping) or not isinstance(profile, Mapping):
        raise ValueError("packet requires semantics and profile objects")
    material_ids = _id_list(semantics.get("material_ids"), "semantics.material_ids")
    boundary_ids = _id_list(semantics.get("boundary_ids"), "semantics.boundary_ids")
    conductor_ids = _id_list(semantics.get("conductor_ids"), "semantics.conductor_ids")
    semantic_contract = {
        "material_ids": material_ids,
        "boundary_ids": boundary_ids,
        "conductor_ids": conductor_ids,
    }
    semantic_sha = _sha(semantic_contract)

    raw_nodes = profile.get("nodes")
    raw_segments = profile.get("segments", [])
    raw_arcs = profile.get("arcs", [])
    raw_regions = profile.get("regions")
    if not isinstance(raw_nodes, list) or len(raw_nodes) < 2:
        raise ValueError("profile.nodes must contain at least two nodes")
    if not isinstance(raw_segments, list) or not isinstance(raw_arcs, list):
        raise ValueError("profile segments and arcs must be arrays")
    if not raw_segments and not raw_arcs:
        raise ValueError("profile needs at least one edge")
    if not isinstance(raw_regions, list) or not raw_regions:
        raise ValueError("profile.regions must be a nonempty array")

    nodes: dict[int, tuple[float, float]] = {}
    normalized_nodes = []
    for index, row in enumerate(raw_nodes):
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{index}] must be an object")
        node_id = int(row.get("id"))
        if node_id < 0 or node_id in nodes:
            raise ValueError("node ids must be unique nonnegative integers")
        point = (_finite(row.get("x_m"), f"nodes[{index}].x_m"), _finite(row.get("y_m"), f"nodes[{index}].y_m"))
        nodes[node_id] = point
        normalized_nodes.append({"id": node_id, "x_m": point[0], "y_m": point[1]})

    degree: Counter[int] = Counter()
    edge_keys: set[tuple[Any, ...]] = set()
    edge_ids: set[int] = set()
    normalized_segments = []
    normalized_arcs = []
    signed_area = 0.0
    perimeter = 0.0

    def endpoints(row: Mapping[str, Any], label: str) -> tuple[int, int, int]:
        edge_id = int(row.get("id"))
        start = int(row.get("start"))
        end = int(row.get("end"))
        if edge_id < 0 or edge_id in edge_ids:
            raise ValueError("edge ids must be unique nonnegative integers")
        if start == end or start not in nodes or end not in nodes:
            raise ValueError(f"{label} has invalid endpoints")
        edge_ids.add(edge_id)
        degree[start] += 1
        degree[end] += 1
        return edge_id, start, end

    for index, row in enumerate(raw_segments):
        if not isinstance(row, Mapping):
            raise ValueError(f"segments[{index}] must be an object")
        label = f"segments[{index}]"
        edge_id, start, end = endpoints(row, label)
        edge_key = ("segment", *sorted((start, end)))
        if edge_key in edge_keys:
            raise ValueError("duplicate line segment")
        edge_keys.add(edge_key)
        boundary, conductor = _edge_semantics(row, label, set(boundary_ids), set(conductor_ids))
        p0, p1 = nodes[start], nodes[end]
        length = math.dist(p0, p1)
        if length <= 0.0:
            raise ValueError(f"{label} has zero length")
        signed_area += 0.5 * (p0[0] * p1[1] - p1[0] * p0[1])
        perimeter += length
        normalized_segments.append({"id": edge_id, "start": start, "end": end, "boundary_id": boundary, "conductor_id": conductor})

    for index, row in enumerate(raw_arcs):
        if not isinstance(row, Mapping):
            raise ValueError(f"arcs[{index}] must be an object")
        label = f"arcs[{index}]"
        edge_id, start, end = endpoints(row, label)
        boundary, conductor = _edge_semantics(row, label, set(boundary_ids), set(conductor_ids))
        sweep = _finite(row.get("sweep_deg"), f"{label}.sweep_deg")
        max_segment = _finite(row.get("max_segment_deg"), f"{label}.max_segment_deg")
        edge_key = ("arc", start, end, sweep)
        if edge_key in edge_keys:
            raise ValueError("duplicate directed arc")
        edge_keys.add(edge_key)
        if max_segment <= 0.0 or max_segment > 180.0:
            raise ValueError(f"{label}.max_segment_deg must be in (0, 180]")
        area, length = _arc_metrics(nodes[start], nodes[end], sweep)
        signed_area += area
        perimeter += length
        normalized_arcs.append({"id": edge_id, "start": start, "end": end, "sweep_deg": sweep, "max_segment_deg": max_segment, "boundary_id": boundary, "conductor_id": conductor})

    if set(degree) != set(nodes) or any(value < 2 for value in degree.values()):
        raise ValueError("profile topology contains an unused node or dangling edge")
    degree_two_loops = all(value == 2 for value in degree.values())

    normalized_regions = []
    for index, row in enumerate(raw_regions):
        if not isinstance(row, Mapping):
            raise ValueError(f"regions[{index}] must be an object")
        region_id = int(row.get("id"))
        material_id = _identifier(row.get("material_id"), f"regions[{index}].material_id")
        if material_id not in material_ids:
            raise ValueError(f"regions[{index}] references unknown material_id {material_id!r}")
        normalized_regions.append({
            "id": region_id,
            "x_m": _finite(row.get("x_m"), f"regions[{index}].x_m"),
            "y_m": _finite(row.get("y_m"), f"regions[{index}].y_m"),
            "material_id": material_id,
        })
    if len({row["id"] for row in normalized_regions}) != len(normalized_regions):
        raise ValueError("region ids must be unique")

    profile_contract = {
        "nodes": sorted(normalized_nodes, key=lambda row: row["id"]),
        "segments": sorted(normalized_segments, key=lambda row: row["id"]),
        "arcs": sorted(normalized_arcs, key=lambda row: row["id"]),
        "regions": sorted(normalized_regions, key=lambda row: row["id"]),
    }
    profile_sha = _sha(profile_contract)
    for key, actual in (("expected_profile_sha256", profile_sha), ("expected_semantics_sha256", semantic_sha)):
        if packet.get(key) is not None and packet.get(key) != actual:
            raise ValueError(f"{key} does not match the normalized contract")
    artifacts, artifact_checks = _artifact_contract(
        packet.get("artifacts"),
        semantic_sha,
        profile_sha,
        abs(signed_area) if degree_two_loops else None,
        perimeter if degree_two_loops else None,
    )
    abi, abi_checks = _abi_contract(packet.get("execution_abi"))
    checks = {
        "closed_topology_without_dangling_edges": True,
        "positive_perimeter": perimeter > 0.0,
        "nonzero_green_area": abs(signed_area) > 0.0,
        "semantic_references_resolved": True,
        **artifact_checks,
        **{f"abi_{key}": value for key, value in abi_checks.items()},
    }
    passed = all(checks.values())
    return {
        "schema": RESULT_SCHEMA,
        "status": "verified" if passed else "needs_attention",
        "pass": passed,
        "profile_contract": profile_contract,
        "profile_sha256": profile_sha,
        "semantics_contract": semantic_contract,
        "semantics_sha256": semantic_sha,
        "measurements": {
            "signed_green_integral_m2": signed_area,
            "absolute_green_integral_m2": abs(signed_area),
            "total_edge_length_m": perimeter,
            "area_is_exact_for_degree2_loops": degree_two_loops,
            "exact_closed_loop_area_m2": abs(signed_area) if degree_two_loops else None,
            "exact_closed_loop_perimeter_m": perimeter if degree_two_loops else None,
            "node_count": len(nodes),
            "edge_count": len(edge_ids),
            "region_count": len(normalized_regions),
        },
        "artifacts": artifacts,
        "execution_abi": abi,
        "checks": checks,
        "cad_ready": "step" in artifacts,
        "solver_ready": "vol" in artifacts,
        "realtime_ready": abi is not None and abi["transport"] in {"matlab_mex", "simulink_s_function"},
        "policy": {
            "step_scope": "geometry_only",
            "semantic_sidecar_required": True,
            "mesh_scope": "dimension-2 Netgen vol",
            "runtime_io": "fixed-width SI fields; no dynamic paths",
        },
    }
