"""Solver-independent cross-method electromagnetic-force report gate."""
from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence


_FRAMES = {"global_cartesian", "global_cylindrical", "local_body"}


def _vector(value: object, field: str) -> tuple[float, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{field} must be a numeric vector")
    try:
        vector = tuple(float(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a numeric vector") from exc
    if len(vector) not in (2, 3) or not all(math.isfinite(item) for item in vector):
        raise ValueError(f"{field} must contain two or three finite components")
    return vector


def _norm(vector: Sequence[float]) -> float:
    return math.sqrt(sum(item * item for item in vector))


def _relative_difference(a: Sequence[float], b: Sequence[float]) -> float:
    return _norm(tuple(x - y for x, y in zip(a, b))) / max(_norm(a), _norm(b), 1.0e-30)


def evaluate_force_report_method_metadata(
    report: Mapping[str, object], relative_tolerance: float = 2.0e-2
) -> dict:
    """Require two independent force methods and an action-reaction closure.

    Each method record supplies ``family``, ``vector`` and ``domain``. The
    method vectors, action force and reaction force must use the same recorded
    frame and force unit. This catches the common false comparison between a
    local component, a global resultant and an incompletely integrated stress
    surface.
    """
    if relative_tolerance <= 0 or not math.isfinite(relative_tolerance):
        raise ValueError("relative_tolerance must be finite and positive")
    methods = report.get("methods")
    if not isinstance(methods, list) or len(methods) < 2:
        raise ValueError("methods must contain at least two force-method records")

    parsed = []
    dimension = None
    for index, item in enumerate(methods):
        if not isinstance(item, Mapping):
            raise ValueError(f"methods[{index}] must be an object")
        family = str(item.get("family", "")).strip().lower()
        domain = str(item.get("domain", "")).strip()
        vector = _vector(item.get("vector"), f"methods[{index}].vector")
        if dimension is None:
            dimension = len(vector)
        elif len(vector) != dimension:
            raise ValueError("all force vectors must have the same dimension")
        parsed.append({"family": family, "domain": domain, "vector": vector})

    action = _vector(report.get("action_force"), "action_force")
    reaction = _vector(report.get("reaction_force"), "reaction_force")
    if len(action) != dimension or len(reaction) != dimension:
        raise ValueError("action/reaction vectors must match the method-vector dimension")

    reference = parsed[0]["vector"]
    method_errors = [
        _relative_difference(reference, method["vector"]) for method in parsed[1:]
    ]
    closure = tuple(a + r for a, r in zip(action, reaction))
    action_reaction_error = _norm(closure) / max(_norm(action), _norm(reaction), 1.0e-30)
    action_reference_error = _relative_difference(action, reference)
    families = {method["family"] for method in parsed if method["family"]}
    checks = {
        "force_unit_is_newton": report.get("force_unit") == "N",
        "component_frame_is_explicit": report.get("component_frame") in _FRAMES,
        "independent_force_families": len(families) >= 2,
        "integration_domains_recorded": all(method["domain"] for method in parsed),
        "method_vectors_agree": max(method_errors, default=0.0) <= relative_tolerance,
        "action_force_matches_methods": action_reference_error <= relative_tolerance,
        "action_reaction_closes": action_reaction_error <= relative_tolerance,
    }
    return {
        "schema": "radia-motor-force-report-method-metadata/v1",
        "policy": "force_claims_require_independent_methods_frame_unit_and_action_reaction",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "force_unit": report.get("force_unit"),
        "component_frame": report.get("component_frame"),
        "method_families": sorted(families),
        "method_relative_errors": method_errors,
        "action_reference_relative_error": action_reference_error,
        "action_reaction_relative_error": action_reaction_error,
        "relative_tolerance": relative_tolerance,
        "checks": checks,
    }


def force_report_method_metadata_gate(
    report_json: str, relative_tolerance: float = 2.0e-2
) -> str:
    try:
        report = json.loads(report_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"report_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(report, dict):
        raise ValueError("report_json must decode to an object")
    return json.dumps(
        evaluate_force_report_method_metadata(report, relative_tolerance),
        indent=2,
        sort_keys=True,
    )
