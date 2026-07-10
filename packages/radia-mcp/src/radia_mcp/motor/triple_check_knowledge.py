"""ELF-seeded comparison planning for radia-motor.

The standard radia-motor workflow is an always-on primary comparison:

1. Use the public ELF/MAGIC MCP surface to choose a motor deck family and
   product-local handoff contract.
2. Verify the current supported finite-element rotating-machine path with the
   ``ngsolve_age`` lane contract.
3. Verify the ``hdiv_vim_reduced_fem`` lane with a solver-ready artifact.  The
   HDiv-VIM rotor source-field plus fixed-stator reduced-FEM coupling remains
   experimental, but a motor result is not accepted for MCP learning unless this
   lane is checked alongside ``ngsolve_age``.
Only the source deck family, public MCP call names, lane IDs, and reduced
engineering lessons belong in public radia-mcp.  Product solver outputs,
private paths, and raw commercial benchmark values stay in the private lane.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from .age_quality_knowledge import route_age_validation_plan
from .validation_lanes_knowledge import (
    lane_template,
    validate_motor_validation_artifact,
)


PRIMARY_COMPARISON_LANES = ("ngsolve_age", "hdiv_vim_reduced_fem")
OPTIONAL_AUXILIARY_LANES: tuple[str, ...] = ()


def validate_motor_source_deck_review_packet(packet: Mapping[str, Any] | str) -> dict[str, Any]:
    """Validate a public source-deck packet before dual-lane motor work."""

    data = json.loads(packet) if isinstance(packet, str) else dict(packet)
    selected = data.get("selected_decks", [])
    required_lanes = data.get("required_lanes", [])
    required_fields = {str(item) for item in data.get("required_result_fields", [])}
    mandatory_fields = {
        "observable_id",
        "observable_unit",
        "coordinate_frame",
        "sign_convention",
        "solver_version",
        "run_date_utc",
        "timing_breakdown_s",
    }
    checks = {
        "schema_matches": data.get("schema_version") == "motor-source-deck-review-packet/v1",
        "observable_id_recorded": bool(str(data.get("observable_id", "")).strip()),
        "selected_decks_present": isinstance(selected, list) and bool(selected),
        "selected_deck_identity_complete": isinstance(selected, list) and bool(selected) and all(
            isinstance(row, Mapping)
            and bool(str(row.get("family", "")).strip())
            and bool(str(row.get("case", "")).strip())
            and str(row.get("mai_path", "")).endswith(".mai")
            for row in selected
        ),
        "both_primary_lanes_required": list(required_lanes) == list(PRIMARY_COMPARISON_LANES),
        "result_contract_complete": mandatory_fields.issubset(required_fields),
        "publication_boundary_recorded": bool(str(data.get("publication_boundary", "")).strip()),
        "raw_solver_results_absent": not any(
            key in data for key in ("solver_results", "raw_results", "benchmark_values")
        ),
    }
    return {
        "policy": "motor_source_deck_review_packet_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "observable_id": str(data.get("observable_id", "")).strip(),
        "selected_deck_count": len(selected) if isinstance(selected, list) else 0,
        "required_lanes": list(required_lanes) if isinstance(required_lanes, list) else [],
        "required_result_fields": sorted(required_fields),
        "checks": checks,
    }


FAMILY_SEEDS: dict[str, dict[str, Any]] = {
    "spm": {
        "source_examples": (
            "application/motor/spm_surface_pm_10/spm001/spm001.mai",
            "application/motor/pm_cosine_pickup_72/pm001/pm001.mai",
            "application/motor/pm_square_2pole_pickup_100/pm001/pm001.mai",
        ),
        "hdiv_observables": ("pickup_flux", "flux_linkage", "demag_field"),
        "age_focus": ("back_emf", "cogging_torque", "ld_lq", "mtpa"),
    },
    "ipm": {
        "source_examples": (
            "application/motor/emdlab_ipm_hairpin_10/eip001/eip001.mai",
            "application/motor/ipm_interior_pm_10/ipm001/ipm001.mai",
            "application/motor/pm_cosine_pickup_72/pm001/pm001.mai",
        ),
        "hdiv_observables": ("pickup_flux", "flux_linkage", "demag_field"),
        "age_focus": ("ld_lq", "mtpa", "field_weakening", "demag_margin"),
    },
    "induction": {
        "source_examples": (
            "application/motor/emdlab_induction_bar_10/eim001/eim001.mai",
            "application/motor/induction_cage_10/im001/im001.mai",
        ),
        "hdiv_observables": ("flux_linkage", "force_or_torque_trend"),
        "age_focus": ("induction_machine", "airgap_eddy_machine", "deep_bar"),
    },
    "srm": {
        "source_examples": (
            "application/motor/emdlab_srm_pole_variants_10/esr001/esr001.mai",
            "application/motor/sr_motor_loop_10/sr001/sr001.mai",
        ),
        "hdiv_observables": ("coenergy", "force_or_torque_trend"),
        "age_focus": ("reluctance_torque", "saturating_inductance"),
    },
    "synrm": {
        "source_examples": (
            "application/motor/emdlab_synrm_flux_barrier_10/esy001/esy001.mai",
            "application/motor/reluctance_motor_10/rel001/rel001.mai",
        ),
        "hdiv_observables": ("coenergy", "force_or_torque_trend"),
        "age_focus": ("synchronous_power_angle", "mtpa", "cross_saturation"),
    },
    "hysteresis": {
        "source_examples": (
            "application/motor/hysteresis_motor_10/hys001/hys001.mai",
        ),
        "hdiv_observables": ("demag_field", "force_or_torque_trend"),
        "age_focus": ("hysteresis_motor_loss", "hysteresis_play"),
    },
}


def _infer_family(goal: str) -> str:
    g = f" {goal.lower()} "
    if any(term in g for term in (" induction", " cage", " im ", "slip", "deep bar")):
        return "induction"
    if any(term in g for term in (" srm", "switched reluctance", " sr motor")):
        return "srm"
    if any(term in g for term in ("synrm", "reluctance motor", "flux barrier")):
        return "synrm"
    if any(term in g for term in ("ipm", "interior", "hairpin", "buried")):
        return "ipm"
    if "hysteresis" in g:
        return "hysteresis"
    return "spm"


def route_motor_triple_check(goal: str) -> dict[str, Any]:
    """Return a structured ELF-seeded multi-lane triple-check plan."""
    family = _infer_family(goal)
    seed = FAMILY_SEEDS[family]
    age_plan = route_age_validation_plan(goal)
    return {
        "schema_version": "radia-motor-triple-check-plan/v1",
        "goal": goal,
        "inferred_family": family,
        "standard_comparison": {
            "primary_required_lanes": list(PRIMARY_COMPARISON_LANES),
            "optional_auxiliary_lanes": list(OPTIONAL_AUXILIARY_LANES),
            "rule": (
                "Every radia-motor validation that claims MCP learning must "
                "compare the NGSolve+AGE lane with the HDiv-VIM/reduced-FEM "
                "lane."
            ),
        },
        "source_mcp_seed": {
            "server": "mcp-server-elf",
            "calls": [
                f'elf_motor_hybrid_router("{goal}")',
                f'elf_sample_decks_route("{goal}", limit=3)',
                f'elf_local_simulation_handoff("{goal}")',
            ],
            "representative_public_decks": list(seed["source_examples"]),
            "public_boundary": (
                "Deck families and MCP call names may be public. Product-run "
                "outputs and raw benchmark values stay private."
            ),
        },
        "radia_lanes": {
            "ngsolve_age": {
                "role": "finite-element air-gap machine lane",
                "support_status": "supported_validation_path",
                "age_focus": list(seed["age_focus"]),
                "age_gate_ids": age_plan["required_gate_ids"],
                "pytest_targets": age_plan["pytest_targets"],
                "artifact_template": lane_template("ngsolve_age"),
            },
            "hdiv_vim_reduced_fem": {
                "role": "mandatory experimental reduced integral / VIM-to-reduced-FEM comparison lane",
                "support_status": "experimental_rfc",
                "observable_candidates": list(seed["hdiv_observables"]),
                "artifact_template": lane_template("hdiv_vim_reduced_fem"),
                "minimum_gate": (
                    "For radia-motor learning, this lane must include a "
                    "solver-ready artifact with a non-empty verification list. "
                    "The full coupling is still experimental until the rotor "
                    "VIM source operator, stator reduced basis, and interface "
                    "operator are regression-tested."
                ),
            },
        },
        "closure": {
            "primary_required_lanes": list(PRIMARY_COMPARISON_LANES),
            "optional_auxiliary_lanes": list(OPTIONAL_AUXILIARY_LANES),
            "required_artifacts": [
                "source_mcp_seed",
                "ngsolve_age artifact",
                "hdiv_vim_reduced_fem solver-ready artifact",
            ],
            "required_gates": [
                'motor_validation_artifact_gate(..., "ngsolve_age")',
                'motor_validation_artifact_gate(..., "hdiv_vim_reduced_fem")',
                "motor_triple_check_artifact_gate(...)",
            ],
            "learning_rule": (
                "radia-motor learned only after the supported AGE lane and "
                "the HDiv-VIM/reduced-FEM lane are both verified in the same "
                "combined artifact, with at least one public-safe MCP "
                "target/test changed and verified."
            ),
        },
    }


def format_motor_triple_check_plan(plan: Mapping[str, Any]) -> str:
    """Format a triple-check plan as Markdown."""
    src = plan["source_mcp_seed"]
    standard = plan["standard_comparison"]
    age = plan["radia_lanes"]["ngsolve_age"]
    hdiv = plan["radia_lanes"]["hdiv_vim_reduced_fem"]
    lines = [
        "# radia-motor triple-check plan",
        "",
        f"- schema: `{plan['schema_version']}`",
        f"- goal: {plan['goal']}",
        f"- inferred family: `{plan['inferred_family']}`",
        "- primary required lanes: "
        + ", ".join(f"`{lane}`" for lane in standard["primary_required_lanes"]),
        "- optional auxiliary lanes: "
        + (", ".join(f"`{lane}`" for lane in standard["optional_auxiliary_lanes"]) or "`none`"),
        f"- comparison rule: {standard['rule']}",
        "",
        "## Source MCP Seed",
        f"- server: `{src['server']}`",
        "- calls:",
    ]
    lines.extend(f"  - `{call}`" for call in src["calls"])
    lines.extend(["- representative public decks:"])
    lines.extend(f"  - `{deck}`" for deck in src["representative_public_decks"])
    lines.extend(
        [
            f"- boundary: {src['public_boundary']}",
            "",
            "## NGSolve+AGE Lane",
            "- role: " + age["role"],
            f"- support status: `{age['support_status']}`",
            "- AGE focus:",
        ]
    )
    lines.extend(f"  - `{item}`" for item in age["age_focus"])
    lines.extend(["- AGE gate IDs:"])
    lines.extend(f"  - `{item}`" for item in age["age_gate_ids"])
    lines.extend(["- pytest targets:"])
    lines.extend(f"  - `{item}`" for item in age["pytest_targets"])
    lines.extend(
        [
            "## HDiv-VIM + Reduced FEM RFC Lane",
            "- role: " + hdiv["role"],
            f"- support status: `{hdiv['support_status']}`",
            "- observable candidates:",
        ]
    )
    lines.extend(f"  - `{item}`" for item in hdiv["observable_candidates"])
    lines.extend(["- minimum gate: " + hdiv["minimum_gate"]])
    lines.extend(["", "## Closure"])
    lines.extend(f"- {item}" for item in plan["closure"]["required_artifacts"])
    lines.extend(["", "Required gates:"])
    lines.extend(f"- `{item}`" for item in plan["closure"]["required_gates"])
    lines.append("")
    lines.append(plan["closure"]["learning_rule"])
    return "\n".join(lines).rstrip()


def _loads_object(value: Mapping[str, Any] | str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    parsed = json.loads(value)
    if not isinstance(parsed, Mapping):
        msg = "artifact JSON must decode to an object"
        raise TypeError(msg)
    return parsed


def validate_motor_triple_check_artifact(
    artifact: Mapping[str, Any] | str,
) -> dict[str, Any]:
    """Validate a combined source + radia motor lane artifact."""
    try:
        data = _loads_object(artifact)
    except (TypeError, json.JSONDecodeError) as exc:
        return {
            "schema_version": "radia-motor-triple-check-artifact-gate/v1",
            "status": "fail",
            "errors": [str(exc)],
            "lane_results": {},
        }

    errors: list[str] = []
    warnings: list[str] = []
    if data.get("schema_version") != "radia-motor-triple-check-artifact/v1":
        errors.append("schema_version must be radia-motor-triple-check-artifact/v1")
    source = data.get("source_mcp_seed", {})
    if not isinstance(source, Mapping):
        errors.append("source_mcp_seed must be an object")
        source = {}
    if not source.get("representative_public_decks"):
        errors.append("source_mcp_seed.representative_public_decks is required")
    if not source.get("source_mcp_calls"):
        errors.append("source_mcp_seed.source_mcp_calls is required")

    lane_artifacts = data.get("lane_artifacts", {})
    if not isinstance(lane_artifacts, Mapping):
        errors.append("lane_artifacts must be an object")
        lane_artifacts = {}

    lane_results: dict[str, Any] = {}
    known_lanes = set(PRIMARY_COMPARISON_LANES)
    for lane_id in lane_artifacts:
        if lane_id not in known_lanes:
            errors.append(f"unknown lane artifact: {lane_id}")

    for lane_id in PRIMARY_COMPARISON_LANES:
        lane_data = lane_artifacts.get(lane_id)
        if lane_data is None:
            errors.append(f"missing lane artifact: {lane_id}")
            continue
        result = validate_motor_validation_artifact(lane_data, lane_id)
        lane_results[lane_id] = result
        if result["status"] != "pass":
            errors.append(f"{lane_id} artifact gate failed")

    feedback = data.get("mcp_feedback", {})
    if not isinstance(feedback, Mapping):
        errors.append("mcp_feedback must be an object")
    else:
        if str(feedback.get("public_status", "")).strip().lower() != "verified":
            warnings.append("mcp_feedback.public_status is not verified")
        if not str(feedback.get("public_summary", "")).strip():
            errors.append("mcp_feedback.public_summary is required")
        if not feedback.get("learning_targets"):
            errors.append("mcp_feedback.learning_targets is required")
        if not feedback.get("verification"):
            errors.append("mcp_feedback.verification is required")

    status = "pass" if not errors else "fail"
    research_triple_check_ready = (
        status == "pass"
        and lane_results.get("ngsolve_age", {}).get("validated_solver_path") is True
        and lane_results.get("hdiv_vim_reduced_fem", {}).get("support_status")
        == "experimental_rfc"
    )
    validated_supported_solver_check = (
        status == "pass"
        and lane_results.get("ngsolve_age", {}).get("validated_solver_path") is True
    )
    validated_dual_solver_check = (
        status == "pass"
        and lane_results.get("ngsolve_age", {}).get("validated_solver_path") is True
        and lane_results.get("hdiv_vim_reduced_fem", {}).get(
            "validated_experimental_solver_path"
        ) is True
    )
    accepted_for_supported_mcp_learning = (
        status == "pass"
        and not warnings
        and validated_supported_solver_check
        and lane_results.get("ngsolve_age", {}).get("accepted_for_mcp_learning")
        is True
    )
    accepted_for_mcp_rfc_learning = (
        status == "pass"
        and not warnings
        and research_triple_check_ready
        and lane_results.get("hdiv_vim_reduced_fem", {}).get(
            "accepted_for_mcp_rfc_learning"
        ) is True
    )
    accepted_for_primary_dual_learning = (
        status == "pass"
        and not warnings
        and validated_dual_solver_check
        and lane_results.get("ngsolve_age", {}).get("accepted_for_mcp_learning")
        is True
        and lane_results.get("hdiv_vim_reduced_fem", {}).get(
            "accepted_for_mcp_learning"
        ) is True
    )
    return {
        "schema_version": "radia-motor-triple-check-artifact-gate/v1",
        "status": status,
        "primary_required_lanes": list(PRIMARY_COMPARISON_LANES),
        "optional_auxiliary_lanes": list(OPTIONAL_AUXILIARY_LANES),
        "research_triple_check_ready": research_triple_check_ready,
        "validated_supported_solver_check": validated_supported_solver_check,
        "validated_dual_solver_check": validated_dual_solver_check,
        "accepted_for_supported_mcp_learning": accepted_for_supported_mcp_learning,
        "accepted_for_mcp_rfc_learning": accepted_for_mcp_rfc_learning,
        "accepted_for_primary_dual_learning": accepted_for_primary_dual_learning,
        "accepted_for_mcp_learning": accepted_for_primary_dual_learning,
        "errors": errors,
        "warnings": warnings,
        "lane_results": lane_results,
    }


def format_triple_check_gate_result(result: Mapping[str, Any]) -> str:
    """Format a triple-check gate result as Markdown."""
    lines = [
        "# Motor triple-check artifact gate",
        "",
        f"- schema: `{result.get('schema_version', '')}`",
        f"- status: `{result.get('status', '')}`",
        "- primary required lanes: "
        + ", ".join(f"`{lane}`" for lane in result.get("primary_required_lanes", ())),
        "- optional auxiliary lanes: "
        + (", ".join(f"`{lane}`" for lane in result.get("optional_auxiliary_lanes", ())) or "`none`"),
        f"- research triple check ready: `{result.get('research_triple_check_ready', False)}`",
        f"- validated supported solver check: `{result.get('validated_supported_solver_check', False)}`",
        f"- validated dual solver check: `{result.get('validated_dual_solver_check', False)}`",
        f"- accepted for supported MCP learning: `{result.get('accepted_for_supported_mcp_learning', False)}`",
        f"- accepted for MCP RFC learning: `{result.get('accepted_for_mcp_rfc_learning', False)}`",
        f"- accepted for primary dual learning: `{result.get('accepted_for_primary_dual_learning', False)}`",
        f"- accepted for MCP learning: `{result.get('accepted_for_mcp_learning', False)}`",
    ]
    lane_results = result.get("lane_results", {})
    if isinstance(lane_results, Mapping):
        lines.extend(["", "## Lane Results"])
        for lane, lane_result in lane_results.items():
            lines.append(
                f"- `{lane}`: `{lane_result.get('status', '')}`, "
                f"accepted=`{lane_result.get('accepted_for_mcp_learning', False)}`"
            )
    errors = list(result.get("errors", ()))
    warnings = list(result.get("warnings", ()))
    if errors:
        lines.extend(["", "## Errors"])
        lines.extend(f"- {item}" for item in errors)
    if warnings:
        lines.extend(["", "## Warnings"])
        lines.extend(f"- {item}" for item in warnings)
    return "\n".join(lines).rstrip()
