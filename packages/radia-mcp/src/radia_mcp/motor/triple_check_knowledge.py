"""ELF-seeded triple-check planning for radia-motor.

The triple-check workflow is:

1. Use the public ELF/MAGIC MCP surface to choose a motor deck family and
   product-local handoff contract.
2. Verify the current supported finite-element rotating-machine path with the
   ``ngsolve_age`` lane contract.
3. Verify the supported coarse/reduced planar path with the ``mmmm2d_coarse``
   lane contract.
4. Organize the new ``hdiv_vim_reduced_fem`` idea as an experimental RFC:
   HDiv-VIM rotor source field plus fixed-stator reduced FEM response.  This is
   not yet a supported radia-motor solver path.

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


FAMILY_SEEDS: dict[str, dict[str, Any]] = {
    "spm": {
        "source_examples": (
            "application/motor/spm_surface_pm_10/spm001/spm001.mai",
            "application/motor/pm_cosine_pickup_72/pm001/pm001.mai",
            "application/motor/pm_square_2pole_pickup_100/pm001/pm001.mai",
        ),
        "hdiv_observables": ("pickup_flux", "flux_linkage", "demag_field"),
        "mmmm_observables": ("torque", "per_region_magnetization", "demag_field"),
        "mmmm_focus": (
            "per-region soft-iron material dictionaries",
            "factor-once torque-angle sweep",
            "planar PM/iron sign and scale check before AGE",
        ),
        "age_focus": ("back_emf", "cogging_torque", "ld_lq", "mtpa"),
    },
    "ipm": {
        "source_examples": (
            "application/motor/emdlab_ipm_hairpin_10/eip001/eip001.mai",
            "application/motor/ipm_interior_pm_10/ipm001/ipm001.mai",
            "application/motor/pm_cosine_pickup_72/pm001/pm001.mai",
        ),
        "hdiv_observables": ("pickup_flux", "flux_linkage", "demag_field"),
        "mmmm_observables": ("torque", "per_region_magnetization", "demag_field"),
        "mmmm_focus": (
            "multi-grade rotor/stator soft-iron regions",
            "coarse torque-angle and demag sanity checks",
            "optimization-loop triage before AGE verification",
        ),
        "age_focus": ("ld_lq", "mtpa", "field_weakening", "demag_margin"),
    },
    "induction": {
        "source_examples": (
            "application/motor/emdlab_induction_bar_10/eim001/eim001.mai",
            "application/motor/induction_cage_10/im001/im001.mai",
        ),
        "hdiv_observables": ("flux_linkage", "force_or_torque_trend"),
        "mmmm_observables": ("force_or_torque_trend", "per_region_magnetization"),
        "mmmm_focus": (
            "coarse rotor/stator magnetization trend",
            "geometry/material sanity check before slip-frequency FE",
        ),
        "age_focus": ("induction_machine", "airgap_eddy_machine", "deep_bar"),
    },
    "srm": {
        "source_examples": (
            "application/motor/emdlab_srm_pole_variants_10/esr001/esr001.mai",
            "application/motor/sr_motor_loop_10/sr001/sr001.mai",
        ),
        "hdiv_observables": ("coenergy", "force_or_torque_trend"),
        "mmmm_observables": ("coenergy", "force_or_torque_trend", "per_region_magnetization"),
        "mmmm_focus": (
            "reluctance-torque sign and periodicity",
            "saturation-free geometry/material sanity check",
        ),
        "age_focus": ("reluctance_torque", "saturating_inductance"),
    },
    "synrm": {
        "source_examples": (
            "application/motor/emdlab_synrm_flux_barrier_10/esy001/esy001.mai",
            "application/motor/reluctance_motor_10/rel001/rel001.mai",
        ),
        "hdiv_observables": ("coenergy", "force_or_torque_trend"),
        "mmmm_observables": ("coenergy", "force_or_torque_trend", "per_region_magnetization"),
        "mmmm_focus": (
            "flux-barrier saliency trend before full FE",
            "multi-region soft-iron consistency check",
        ),
        "age_focus": ("synchronous_power_angle", "mtpa", "cross_saturation"),
    },
    "hysteresis": {
        "source_examples": (
            "application/motor/hysteresis_motor_10/hys001/hys001.mai",
        ),
        "hdiv_observables": ("demag_field", "force_or_torque_trend"),
        "mmmm_observables": ("demag_field", "force_or_torque_trend"),
        "mmmm_focus": (
            "demag-field sign and scale",
            "linearized precheck before hysteresis-specific validation",
        ),
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
            "mmmm2d_coarse": {
                "role": "supported coarse/reduced planar MMMM lane",
                "support_status": "supported_coarse_path",
                "observable_candidates": list(seed["mmmm_observables"]),
                "mmmm_focus": list(seed["mmmm_focus"]),
                "pytest_targets": ["validation_test/feec/test_moment2d_perregion.py"],
                "artifact_template": lane_template("mmmm2d_coarse"),
                "minimum_gate": (
                    "Use this as a verified fast planar sign/scale and "
                    "per-region material gate. It can train motor routing and "
                    "coarse optimization advice, but it does not replace AGE "
                    "for final rotating-machine validation."
                ),
            },
            "hdiv_vim_reduced_fem": {
                "role": "experimental reduced integral / VIM-to-reduced-FEM coupling RFC",
                "support_status": "experimental_rfc",
                "observable_candidates": list(seed["hdiv_observables"]),
                "artifact_template": lane_template("hdiv_vim_reduced_fem"),
                "minimum_gate": (
                    "A private source comparison may be collected as a research "
                    "artifact, but this lane is not a validated solver path "
                    "until the rotor VIM source operator, stator reduced basis, "
                    "and interface operator are implemented and regression-tested."
                ),
            },
        },
        "closure": {
            "required_artifacts": [
                "source_mcp_seed",
                "ngsolve_age artifact",
                "mmmm2d_coarse artifact",
                "hdiv_vim_reduced_fem RFC artifact",
            ],
            "required_gates": [
                'motor_validation_artifact_gate(..., "ngsolve_age")',
                'motor_validation_artifact_gate(..., "mmmm2d_coarse")',
                'motor_validation_artifact_gate(..., "hdiv_vim_reduced_fem")',
                "motor_triple_check_artifact_gate(...)",
            ],
            "learning_rule": (
                "radia-motor learned only after the supported AGE lane and "
                "the supported coarse MMMM lane are verified, the "
                "HDiv/reduced-FEM RFC is explicitly labeled as experimental, "
                "and at least one public-safe MCP target/test changed and "
                "verified."
            ),
        },
    }


def format_motor_triple_check_plan(plan: Mapping[str, Any]) -> str:
    """Format a triple-check plan as Markdown."""
    src = plan["source_mcp_seed"]
    age = plan["radia_lanes"]["ngsolve_age"]
    mmmm = plan["radia_lanes"]["mmmm2d_coarse"]
    hdiv = plan["radia_lanes"]["hdiv_vim_reduced_fem"]
    lines = [
        "# radia-motor triple-check plan",
        "",
        f"- schema: `{plan['schema_version']}`",
        f"- goal: {plan['goal']}",
        f"- inferred family: `{plan['inferred_family']}`",
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
            "",
            "## 2D MMMM Coarse Lane",
            "- role: " + mmmm["role"],
            f"- support status: `{mmmm['support_status']}`",
            "- observable candidates:",
        ]
    )
    lines.extend(f"  - `{item}`" for item in mmmm["observable_candidates"])
    lines.extend(["- MMMM focus:"])
    lines.extend(f"  - {item}" for item in mmmm["mmmm_focus"])
    lines.extend(["- pytest targets:"])
    lines.extend(f"  - `{item}`" for item in mmmm["pytest_targets"])
    lines.extend(
        [
            "- minimum gate: " + mmmm["minimum_gate"],
            "",
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
    required_lanes = ("ngsolve_age", "mmmm2d_coarse", "hdiv_vim_reduced_fem")
    for lane_id in required_lanes:
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
        and lane_results.get("mmmm2d_coarse", {}).get("validated_coarse_path") is True
        and lane_results.get("hdiv_vim_reduced_fem", {}).get("support_status")
        == "experimental_rfc"
    )
    validated_supported_solver_check = (
        status == "pass"
        and lane_results.get("ngsolve_age", {}).get("validated_solver_path") is True
        and lane_results.get("mmmm2d_coarse", {}).get("validated_coarse_path") is True
    )
    validated_dual_solver_check = (
        status == "pass"
        and all(
            lane_results.get(lane, {}).get("validated_solver_path") is True
            for lane in ("hdiv_vim_reduced_fem", "ngsolve_age")
        )
    )
    return {
        "schema_version": "radia-motor-triple-check-artifact-gate/v1",
        "status": status,
        "research_triple_check_ready": research_triple_check_ready,
        "validated_supported_solver_check": validated_supported_solver_check,
        "validated_dual_solver_check": validated_dual_solver_check,
        "accepted_for_mcp_learning": status == "pass"
        and not warnings
        and all(
            lane_results.get(lane, {}).get("accepted_for_mcp_learning") is True
            for lane in required_lanes
        ),
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
        f"- research triple check ready: `{result.get('research_triple_check_ready', False)}`",
        f"- validated supported solver check: `{result.get('validated_supported_solver_check', False)}`",
        f"- validated dual solver check: `{result.get('validated_dual_solver_check', False)}`",
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
