"""Readable motor-loss handoff gate for lumped and 3D thermal models."""
from __future__ import annotations

import json
import math
import re
from collections import Counter, deque
from typing import Any


HEX_CELL_TYPES = {"hex", "hex8", "hex20", "hex27", "hexahedron"}
ELECTROTHERMAL_STAGE_ORDER = (
    "electromagnetic",
    "stator_core_loss",
    "rotor_core_loss",
    "thermal",
)
ELECTROTHERMAL_SOURCE_OWNERS = {
    "rotor_conductor_joule": "electromagnetic",
    "phase_u_joule": "electromagnetic",
    "phase_v_joule": "electromagnetic",
    "phase_w_joule": "electromagnetic",
    "stator_core_loss": "stator_core_loss",
    "rotor_core_loss": "rotor_core_loss",
}
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _load_json(value: str, expected: type, name: str) -> Any:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} must be valid JSON: {exc.msg}") from exc
    if not isinstance(parsed, expected):
        raise ValueError(f"{name} must decode to {expected.__name__}")
    return parsed


def _close(a: float, b: float, relative_tolerance: float) -> bool:
    scale = max(abs(a), abs(b), 1.0)
    return abs(a - b) <= relative_tolerance * scale


def _reachable_nodes(nodes: set[str], branches: list[dict], ambient: str) -> set[str]:
    graph = {node: set() for node in nodes}
    for branch in branches:
        a = str(branch.get("from") or "")
        b = str(branch.get("to") or "")
        if a in graph and b in graph:
            graph[a].add(b)
            graph[b].add(a)
    seen = {ambient} if ambient in graph else set()
    queue = deque(seen)
    while queue:
        node = queue.popleft()
        for neighbour in graph[node] - seen:
            seen.add(neighbour)
            queue.append(neighbour)
    return seen


def evaluate_motor_thermal_handoff(
    loss_buckets: dict,
    network: dict,
    mesh_regions: list,
    relative_tolerance: float = 1.0e-9,
) -> dict:
    """Check one loss table against an LPTN and an all-hex thermal mesh.

    ``loss_buckets`` maps physical region names to non-negative watts.
    Network nodes map those regions with ``source_regions`` and connect to the
    fixed ``ambient_node`` through positive ``resistance_K_per_W`` branches.
    Mesh rows map the same regions to positive hexahedral cell counts and
    ``loss_W``.  The gate checks topology, one-to-one region ownership, and
    conservation of total heat without pretending to solve either model.
    """
    errors: list[str] = []
    if relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be positive")

    losses: dict[str, float] = {}
    for raw_name, raw_value in loss_buckets.items():
        name = str(raw_name).strip()
        if not name:
            errors.append("loss bucket names must be non-empty")
            continue
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            errors.append(f"loss bucket {name!r} is not numeric")
            continue
        if value < 0:
            errors.append(f"loss bucket {name!r} must be non-negative")
        losses[name] = value
    if not losses:
        errors.append("at least one loss bucket is required")

    nodes = network.get("nodes")
    branches = network.get("branches")
    ambient = str(network.get("ambient_node") or "").strip()
    if not isinstance(nodes, list) or not nodes:
        errors.append("network.nodes must be a non-empty list")
        nodes = []
    if not isinstance(branches, list) or not branches:
        errors.append("network.branches must be a non-empty list")
        branches = []

    node_ids = [str(row.get("id") or "").strip() for row in nodes if isinstance(row, dict)]
    node_set = {item for item in node_ids if item}
    if len(node_ids) != len(nodes) or len(node_set) != len(node_ids):
        errors.append("network node ids must be present and unique")
    if not ambient or ambient not in node_set:
        errors.append("network.ambient_node must name one network node")

    network_regions: list[str] = []
    for row in nodes:
        if not isinstance(row, dict):
            continue
        node_id = str(row.get("id") or "").strip()
        regions = row.get("source_regions", [])
        if not isinstance(regions, list):
            errors.append(f"node {node_id!r} source_regions must be a list")
            continue
        network_regions.extend(str(region).strip() for region in regions if str(region).strip())
        if node_id != ambient:
            try:
                capacitance = float(row.get("capacitance_J_per_K"))
            except (TypeError, ValueError):
                capacitance = -1.0
            if capacitance <= 0:
                errors.append(f"node {node_id!r} needs positive capacitance_J_per_K")

    for index, branch in enumerate(branches):
        if not isinstance(branch, dict):
            errors.append(f"network branch {index} must be an object")
            continue
        a = str(branch.get("from") or "").strip()
        b = str(branch.get("to") or "").strip()
        if a not in node_set or b not in node_set or a == b:
            errors.append(f"network branch {index} has invalid endpoints")
        try:
            resistance = float(branch.get("resistance_K_per_W"))
        except (TypeError, ValueError):
            resistance = -1.0
        if resistance <= 0:
            errors.append(f"network branch {index} needs positive resistance_K_per_W")

    reachable = _reachable_nodes(node_set, branches, ambient)
    if node_set and reachable != node_set:
        errors.append("every thermal-network node must have a path to ambient")

    network_counts = Counter(network_regions)
    missing_network = sorted(set(losses) - set(network_counts))
    extra_network = sorted(set(network_counts) - set(losses))
    duplicate_network = sorted(name for name, count in network_counts.items() if count != 1)
    if missing_network:
        errors.append(f"network is missing loss regions: {missing_network}")
    if extra_network:
        errors.append(f"network has unknown loss regions: {extra_network}")
    if duplicate_network:
        errors.append(f"network loss regions must be assigned once: {duplicate_network}")

    mesh_counts: Counter[str] = Counter()
    mesh_loss_by_region: dict[str, float] = {}
    for index, row in enumerate(mesh_regions):
        if not isinstance(row, dict):
            errors.append(f"mesh region {index} must be an object")
            continue
        region = str(row.get("region") or "").strip()
        cell_type = str(row.get("cell_type") or "").strip().lower()
        if not region:
            errors.append(f"mesh region {index} needs a region name")
            continue
        mesh_counts[region] += 1
        if cell_type not in HEX_CELL_TYPES:
            errors.append(f"mesh region {region!r} must use a hexahedral cell type")
        try:
            cell_count = int(row.get("cell_count"))
        except (TypeError, ValueError):
            cell_count = 0
        if cell_count <= 0:
            errors.append(f"mesh region {region!r} needs a positive cell_count")
        try:
            mesh_loss = float(row.get("loss_W"))
        except (TypeError, ValueError):
            mesh_loss = -1.0
        if mesh_loss < 0:
            errors.append(f"mesh region {region!r} needs non-negative loss_W")
        mesh_loss_by_region[region] = mesh_loss

    missing_mesh = sorted(set(losses) - set(mesh_counts))
    extra_mesh = sorted(set(mesh_counts) - set(losses))
    duplicate_mesh = sorted(name for name, count in mesh_counts.items() if count != 1)
    if missing_mesh:
        errors.append(f"mesh is missing loss regions: {missing_mesh}")
    if extra_mesh:
        errors.append(f"mesh has unknown loss regions: {extra_mesh}")
    if duplicate_mesh:
        errors.append(f"mesh loss regions must be assigned once: {duplicate_mesh}")

    mismatched_mesh = sorted(
        name for name, value in losses.items()
        if name in mesh_loss_by_region
        and not _close(value, mesh_loss_by_region[name], relative_tolerance)
    )
    if mismatched_mesh:
        errors.append(f"mesh loss_W does not match source buckets: {mismatched_mesh}")

    source_total = sum(losses.values())
    network_total = sum(losses.get(name, 0.0) for name in network_regions)
    mesh_total = sum(value for value in mesh_loss_by_region.values() if value >= 0)
    totals_match = (
        _close(source_total, network_total, relative_tolerance)
        and _close(source_total, mesh_total, relative_tolerance)
    )
    if not totals_match:
        errors.append("source, network, and mesh total heat must agree")

    return {
        "schema": "radia-motor-thermal-handoff/v1",
        "policy": "loss_to_lptn_and_hex_thermal_mesh",
        "status": "ok" if not errors else "needs_attention",
        "source_total_loss_W": source_total,
        "network_total_loss_W": network_total,
        "mesh_total_loss_W": mesh_total,
        "loss_region_count": len(losses),
        "network_node_count": len(node_set),
        "network_branch_count": len(branches),
        "mesh_region_count": len(mesh_counts),
        "checks": {
            "network_regions_one_to_one": not (missing_network or extra_network or duplicate_network),
            "network_connected_to_ambient": bool(node_set) and reachable == node_set,
            "mesh_regions_one_to_one": not (missing_mesh or extra_mesh or duplicate_mesh),
            "mesh_all_hexahedral": all(
                isinstance(row, dict)
                and str(row.get("cell_type") or "").strip().lower() in HEX_CELL_TYPES
                for row in mesh_regions
            ),
            "regional_loss_values_match": not mismatched_mesh,
            "total_heat_conserved": totals_match,
        },
        "errors": errors,
        "next_steps": [
            "Run the LPTN for fast parameter sweeps and controller-oriented studies.",
            "Run the 3D all-hex thermal model for spatial gradients and hotspot resolution.",
            "Compare region-average temperatures and total boundary heat flow before promotion.",
        ],
    }


def motor_thermal_handoff_gate(
    loss_buckets_json: str,
    network_json: str,
    mesh_regions_json: str,
    relative_tolerance: float = 1.0e-9,
) -> str:
    """JSON wrapper used by the MCP tool."""
    losses = _load_json(loss_buckets_json, dict, "loss_buckets_json")
    network = _load_json(network_json, dict, "network_json")
    mesh = _load_json(mesh_regions_json, list, "mesh_regions_json")
    result = evaluate_motor_thermal_handoff(losses, network, mesh, relative_tolerance)
    return json.dumps(result, indent=2, sort_keys=True)


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _close_power(
    actual: float,
    expected: float,
    *,
    absolute_tolerance_W: float,
    relative_tolerance: float,
) -> bool:
    return abs(actual - expected) <= max(
        absolute_tolerance_W,
        relative_tolerance * max(abs(actual), abs(expected), 1.0),
    )


def evaluate_motor_electrothermal_result_chain(
    chain: dict,
    absolute_tolerance_W: float = 1.1e-2,
    relative_tolerance: float = 1.0e-3,
) -> dict:
    """Gate a fresh four-stage motor electrothermal result handoff.

    The gate is intentionally solver-neutral.  It checks result identity and
    dependency digests before checking physics: three-phase winding and rotor
    Joule losses plus stator/rotor core losses must be owned exactly once, and
    the thermal input must equal the upstream full-model loss multiplied by an
    explicit symmetry fraction.  A steady thermal result must then show a
    finite temperature rise above ambient.
    """
    if absolute_tolerance_W <= 0 or relative_tolerance <= 0:
        raise ValueError("power tolerances must be positive")

    errors: list[str] = []
    stages = chain.get("stages")
    if not isinstance(stages, list):
        stages = []
        errors.append("stages must be a list")
    stage_names = [str(row.get("stage") or "") for row in stages if isinstance(row, dict)]
    stage_order_ok = stage_names == list(ELECTROTHERMAL_STAGE_ORDER)
    if not stage_order_ok:
        errors.append("stages must follow electromagnetic, stator loss, rotor loss, thermal order")

    stage_by_name = {
        str(row.get("stage")): row
        for row in stages
        if isinstance(row, dict) and str(row.get("stage") or "")
    }
    artifact_ids: list[str] = []
    result_digests: list[str] = []
    stage_evidence_ok = len(stages) == len(ELECTROTHERMAL_STAGE_ORDER)
    for stage_name in ELECTROTHERMAL_STAGE_ORDER:
        row = stage_by_name.get(stage_name)
        if not isinstance(row, dict):
            stage_evidence_ok = False
            continue
        artifact_id = str(row.get("artifact_id") or "").strip()
        digest = str(row.get("result_digest") or "").strip().lower()
        solve_s = _finite_float(row.get("solve_s"))
        artifact_ids.append(artifact_id)
        result_digests.append(digest)
        if (
            not artifact_id
            or not SHA256_PATTERN.fullmatch(digest)
            or row.get("completed") is not True
            or row.get("fresh") is not True
            or solve_s is None
            or solve_s <= 0
        ):
            stage_evidence_ok = False
    stage_identity_ok = (
        stage_evidence_ok
        and len(set(artifact_ids)) == len(ELECTROTHERMAL_STAGE_ORDER)
        and len(set(result_digests)) == len(ELECTROTHERMAL_STAGE_ORDER)
    )
    if not stage_identity_ok:
        errors.append("every stage needs a unique fresh result digest and positive solve time")

    dependency_ok = stage_identity_ok and stage_order_ok
    if dependency_ok:
        em = stage_by_name["electromagnetic"]
        stator = stage_by_name["stator_core_loss"]
        rotor = stage_by_name["rotor_core_loss"]
        thermal = stage_by_name["thermal"]
        em_dependency = {}
        stator_dependency = {em["artifact_id"]: em["result_digest"].lower()}
        rotor_dependency = {em["artifact_id"]: em["result_digest"].lower()}
        thermal_dependency = {
            em["artifact_id"]: em["result_digest"].lower(),
            stator["artifact_id"]: stator["result_digest"].lower(),
            rotor["artifact_id"]: rotor["result_digest"].lower(),
        }
        dependency_ok = (
            em.get("input_result_digests") == em_dependency
            and stator.get("input_result_digests") == stator_dependency
            and rotor.get("input_result_digests") == rotor_dependency
            and thermal.get("input_result_digests") == thermal_dependency
        )
    if not dependency_ok:
        errors.append("downstream stages must pin the exact upstream result digests")

    raw_sources = chain.get("source_buckets")
    if not isinstance(raw_sources, list):
        raw_sources = []
        errors.append("source_buckets must be a list")
    source_rows = [row for row in raw_sources if isinstance(row, dict)]
    channel_counts = Counter(str(row.get("channel") or "") for row in source_rows)
    source_channels_ok = (
        len(source_rows) == len(ELECTROTHERMAL_SOURCE_OWNERS)
        and set(channel_counts) == set(ELECTROTHERMAL_SOURCE_OWNERS)
        and all(count == 1 for count in channel_counts.values())
    )
    source_values_ok = source_channels_ok and stage_identity_ok
    source_total_W = 0.0
    for row in source_rows:
        channel = str(row.get("channel") or "")
        expected_owner = ELECTROTHERMAL_SOURCE_OWNERS.get(channel)
        power_W = _finite_float(row.get("power_W"))
        owner = stage_by_name.get(expected_owner or "")
        if (
            expected_owner is None
            or owner is None
            or row.get("upstream_stage") != expected_owner
            or row.get("upstream_artifact_id") != owner.get("artifact_id")
            or str(row.get("upstream_result_digest") or "").lower()
            != str(owner.get("result_digest") or "").lower()
            or power_W is None
            or power_W < 0
        ):
            source_values_ok = False
        elif power_W is not None:
            source_total_W += power_W
    if not source_channels_ok:
        errors.append("the six motor loss channels must each be owned exactly once")
    if not source_values_ok:
        errors.append("loss channels need non-negative power and matching upstream identity")

    symmetry_fraction = _finite_float(chain.get("symmetry_fraction"))
    symmetry_ok = symmetry_fraction is not None and 0 < symmetry_fraction <= 1
    if not symmetry_ok:
        errors.append("symmetry_fraction must be in (0, 1]")
    expected_thermal_input_W = (
        source_total_W * symmetry_fraction if symmetry_fraction is not None else math.nan
    )
    thermal_summary = chain.get("thermal_summary")
    if not isinstance(thermal_summary, dict):
        thermal_summary = {}
        errors.append("thermal_summary must be an object")
    thermal_input_W = _finite_float(thermal_summary.get("input_power_W"))
    power_closure_ok = (
        source_values_ok
        and symmetry_ok
        and thermal_input_W is not None
        and thermal_input_W >= 0
        and _close_power(
            thermal_input_W,
            expected_thermal_input_W,
            absolute_tolerance_W=absolute_tolerance_W,
            relative_tolerance=relative_tolerance,
        )
    )
    if not power_closure_ok:
        errors.append("thermal input must close to symmetry-scaled upstream losses")

    ambient_C = _finite_float(thermal_summary.get("ambient_temperature_C"))
    maximum_C = _finite_float(thermal_summary.get("maximum_temperature_C"))
    temperature_ok = (
        thermal_summary.get("steady_state") is True
        and ambient_C is not None
        and maximum_C is not None
        and maximum_C > ambient_C
    )
    if not temperature_ok:
        errors.append("steady thermal result must have a finite positive temperature rise")

    checks = {
        "schema_matches": chain.get("schema") == "motor-electrothermal-result-chain/v1",
        "four_stage_order": stage_order_ok,
        "fresh_unique_stage_results": stage_identity_ok,
        "upstream_result_digests_pinned": dependency_ok,
        "six_loss_channels_owned_once": source_channels_ok,
        "loss_channel_values_and_identity_valid": source_values_ok,
        "symmetry_fraction_explicit": symmetry_ok,
        "symmetry_scaled_power_closure": power_closure_ok,
        "steady_temperature_rise_present": temperature_ok,
    }
    if not checks["schema_matches"]:
        errors.append("schema must be motor-electrothermal-result-chain/v1")

    return {
        "schema": "radia-motor-electrothermal-result-chain/v1",
        "policy": "fresh_digest_pinned_symmetry_scaled_power_handoff",
        "status": "ok" if all(checks.values()) and not errors else "needs_attention",
        "checks": checks,
        "metrics": {
            "source_total_loss_W": source_total_W,
            "symmetry_fraction": symmetry_fraction,
            "expected_thermal_input_W": expected_thermal_input_W,
            "thermal_input_W": thermal_input_W,
            "power_closure_error_W": (
                abs(thermal_input_W - expected_thermal_input_W)
                if thermal_input_W is not None and math.isfinite(expected_thermal_input_W)
                else None
            ),
            "ambient_temperature_C": ambient_C,
            "maximum_temperature_C": maximum_C,
            "temperature_rise_C": (
                maximum_C - ambient_C
                if maximum_C is not None and ambient_C is not None
                else None
            ),
        },
        "errors": errors,
    }


def motor_electrothermal_result_chain_gate(
    chain_json: str,
    absolute_tolerance_W: float = 1.1e-2,
    relative_tolerance: float = 1.0e-3,
) -> str:
    """JSON wrapper for the fresh electrothermal result-chain gate."""
    chain = _load_json(chain_json, dict, "chain_json")
    result = evaluate_motor_electrothermal_result_chain(
        chain,
        absolute_tolerance_W=absolute_tolerance_W,
        relative_tolerance=relative_tolerance,
    )
    return json.dumps(result, indent=2, sort_keys=True)
