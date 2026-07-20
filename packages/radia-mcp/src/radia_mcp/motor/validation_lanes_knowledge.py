"""Public-safe validation-lane contract for radia-motor.

radia-motor has several motor-validation paths that should not be collapsed
into one score:

* NGSolve+AGE for the finite-element air-gap machine path: torque, dq,
  eddy, and nonlinear machine quantities.
* Radia HDiv-MMM coupled to the HCurl eddy-bubble basis for the independent
  material and eddy-current path.

This module gives MCP clients a small contract for naming the lane, checking
artifact metadata, and deciding which radia-motor knowledge should be updated
after a cross-validation slot.  It deliberately avoids private paths, product
solver outputs, and commercial benchmark numbers.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class MotorValidationLane:
    """One radia-motor validation lane."""

    lane_id: str
    label: str
    support_status: str
    support_note: str
    radia_path: str
    best_for: tuple[str, ...]
    observable_families: tuple[str, ...]
    required_fields: tuple[str, ...]
    required_metrics: tuple[str, ...]
    public_evidence: tuple[str, ...]
    private_reference_sources: tuple[str, ...]
    promotion_targets: tuple[str, ...]


COMMON_REQUIRED_FIELDS = (
    "schema_version",
    "timestamp_utc",
    "radia_version",
    "motor_validation_lane",
    "reference_source_class",
    "observable_family",
    "case_count",
    "status",
    "tolerances",
    "metrics",
    "timing_breakdown_s",
    "artifact_feedback",
    "shared_mesh_material_identity",
    "solver_ready_artifact",
)


LANES: dict[str, MotorValidationLane] = {
    "ngsolve_age": MotorValidationLane(
        lane_id="ngsolve_age",
        label="NGSolve+AGE",
        support_status="supported_validation_path",
        support_note=(
            "This is the current radia-motor supported validation path for "
            "2D rotating-machine finite-element studies."
        ),
        radia_path="radia-ngsolve air-gap element finite-element motor path",
        best_for=(
            "air-gap field and Maxwell-stress torque",
            "coenergy torque-angle and cogging periodicity",
            "dq quantities, MTPA, and field-weakening checks",
            "slip-frequency eddy-current and hysteresis-loss anchors",
            "linear motor air-gap thrust, flux-linkage, and end-effect reduced checks",
            "rotary motor family sweeps across SPM, BLDC, IPM, induction, SRM, SynRM, and AFPM",
            "nonlinear material iteration once the reduced invariant is clear",
        ),
        observable_families=(
            "airgap_flux",
            "torque",
            "motor_family_sweep",
            "rotary_flux_linkage",
            "reluctance_torque",
            "linear_thrust",
            "linear_pm_flux",
            "cogging_torque",
            "back_emf",
            "ld_lq",
            "mtpa",
            "field_weakening",
            "slip_loss",
            "hysteresis_loss",
        ),
        required_fields=COMMON_REQUIRED_FIELDS
        + (
            "age_gate_ids",
            "pytest_targets",
        ),
        required_metrics=(
            "field_relative_error",
            "torque_relative_error",
            "energy_or_coenergy_relative_error",
            "quantity_specific_residual",
        ),
        public_evidence=(
            "tests/test_airgap_element.py",
            "validation_test/radia_mcp/test_airgap_ngsolve_coupling.py",
            "validation_test/radia_mcp/test_airgap_two_region.py",
            "validation_test/radia_mcp/test_airgap_machine_rotation.py",
            "validation_test/radia_mcp/test_airgap_eddy_machine.py",
            "validation_test/radia_mcp/test_build123d_ipm_age_torque.py",
        ),
        private_reference_sources=(
            "product local reference",
            "lab-local finite-element reference",
            "open-source reference",
            "stored regression reference",
        ),
        promotion_targets=(
            "radia_mcp.motor.age_quality_knowledge",
            "radia_mcp.radia_ngsolve AGE / force recipes",
        ),
    ),
    "hdiv_mmm_hcurl_eddy_bubble": MotorValidationLane(
        lane_id="hdiv_mmm_hcurl_eddy_bubble",
        label="Radia HDiv-MMM + HCurl eddy-bubble",
        support_status="required_validation_path",
        support_note=(
            "The single-rotor planar reluctance path is implemented as "
            "radia.motor_hdiv.HDivReducedMotor and the radia_motor 'HDiv Reduced' "
            "study. It reuses one symmetric BDM1 or BDM2 charge Gram "
            "(BDM1/Q2 or BDM2/Q3 geometry) and checks torque "
            "through Maxwell stress, magnetization-volume coupling, and "
            "fixed-current coenergy. PlanarDemagBody.field_cf is the native "
            "rotating source/target-frame interface. A fixed-stator reduced-FEM "
            "basis and full AGE/transient coupling are still RFC work and must "
            "not be inferred from the validated reduced rotor path."
        ),
        radia_path=(
            "radia.vim.NgsolveHDivMMMResponseReduction + "
            "radia.vim.CoupleEddyBubbleHCurlBasisWithHDivMMM"
        ),
        best_for=(
            "passive pickup flux and signed flux-linkage sweeps",
            "linear PM motor thrust and eddy-current reaction checks",
            "rotary SPM/BLDC/IPM/IM/SRM/SynRM/AFPM flux, loss, and coenergy checks",
            "permanent-magnet demagnetizing-field anchors",
            "magnetic-material and conductor coupling on a shared mesh/material identity",
            "frequency-domain eddy-current and Joule-loss checks",
            "Maxwell, volume-force, coenergy, and magnetic-energy consistency",
        ),
        observable_families=(
            "pickup_flux",
            "flux_linkage",
            "motor_family_sweep",
            "rotary_flux_linkage",
            "linear_pm_flux",
            "linear_force_or_thrust",
            "demag_field",
            "coenergy",
            "force_or_torque_trend",
            "eddy_current",
            "joule_loss",
            "frequency_response",
        ),
        required_fields=COMMON_REQUIRED_FIELDS
        + (
            "hdiv_mmm_operator_contract",
            "hcurl_eddy_bubble_contract",
            "coupling_operator_contract",
        ),
        required_metrics=(
            "signed_agreement_count",
            "mean_abs_relative_error",
            "rms_abs_relative_error",
            "max_abs_relative_error",
            "mixed_block_residual",
            "magnetic_energy_closure",
            "eddy_power_nonnegative",
        ),
        public_evidence=(
            "analytic sign/scale checks",
            "tests/test_vim_eddy_hybrid.py::test_eddy_bubble_hcurl_basis_is_vim_and_hdiv_mmm_ready",
            "validation_test/cln/hcurl_vim_hdiv_mmm_end_to_end.py",
            "validation_test/cln/planar_hdiv_mmm_response_smoke.py",
            "stored public-safe regression artifacts",
        ),
        private_reference_sources=(
            "product local reference",
            "lab-local finite-element reference",
            "open-source reference",
            "stored regression reference",
        ),
        promotion_targets=(
            "radia_mcp.motor.validation_lanes_knowledge",
            "radia_mcp.radia_ngsolve force / flux recipes when applicable",
        ),
    ),
}


LEGACY_LANE_ALIASES = {
    "hdiv_vim_reduced_fem": "hdiv_mmm_hcurl_eddy_bubble",
}


def _canonical_lane_id(lane_id: str) -> str:
    normalized = lane_id.strip().lower()
    return LEGACY_LANE_ALIASES.get(normalized, normalized)


OVERVIEW = """\
# radia-motor validation lanes

radia-motor should keep independent cross-validation lanes, but its default
learning environment is an always-on two-lane comparison: every motor result
that claims radia-motor MCP learning must include both `ngsolve_age` and
`hdiv_mmm_hcurl_eddy_bubble` in the same combined artifact.  A single-lane artifact
can pass its own metadata gate, but it is not enough to say radia-motor learned.

- `ngsolve_age`: NGSolve+AGE.  This is the current supported radia-motor
  finite-element air-gap machine lane for torque, dq quantities, cogging,
  eddy/slip, hysteresis, and nonlinear machine studies.
- `hdiv_mmm_hcurl_eddy_bubble`: Radia HDiv-MMM + HCurl eddy-bubble.  The mixed
  operator is the independent material/eddy-current lane.  Every promoted
  motor artifact must record its shared mesh/material identity, solver-ready
  artifact, and concrete execution verification.

Private product, lab-local, open-source, analytic, and stored-regression
comparisons can all be reference sources.  The public MCP learning artifact
must still say which radia lane was exercised, which observable family was
checked, what metrics/tolerances were used, and which public-safe knowledge
target was updated.  Do not merge the lanes into one vague "motor passed"
result.
"""


SOURCE_POLICY = """\
# Reference-source policy

Commercial or lab-private tools are reference sources, not public evidence.
They may train radia-motor locally, but the public artifact must keep only the
generalized engineering lesson:

- Do record `reference_source_class`, such as `product_local_reference`,
  `opensource_reference`, `analytic_reference`, or `stored_regression`.
- Do use public-safe source classes such as `product_local_reference`,
  `lab_local_reference`, `opensource_reference`, `analytic_reference`, or
  `stored_regression` when a private run trained the lane.
- Do record the radia lane, observable family, tolerances, aggregate metrics,
  timing breakdown, and promotion target.
- Do not record private absolute paths, solver logs, product benchmark tables,
  license details, or raw commercial case files in public radia-mcp text.

This lets private product and lab-local references strengthen both radia lanes
while keeping the publication boundary clean.
"""


PROMOTION_POLICY = """\
# Artifact-to-MCP promotion policy

Each motor cross-validation slot should end with three decisions:

1. `which_lane`: `ngsolve_age` or `hdiv_mmm_hcurl_eddy_bubble`.
2. `which_observable`: one lane-supported observable family.
3. `which_promotion`: the exact public-safe knowledge or recipe that improved.

Promotion is allowed when:

- the artifact passes its own tolerance,
- the lane metadata is complete,
- the observable belongs to that lane,
- the timing breakdown has at least one named phase, and
- the lesson can be stated without private paths, product logs, or commercial
  benchmark numbers.

For radia-motor learning, promote through the combined comparison gate:
`ngsolve_age` and `hdiv_mmm_hcurl_eddy_bubble` must both be present. The Radia
lane must include a solver-ready mixed-system artifact with a non-empty
verification list.
If the slot only produced a useful private comparison, keep it in the private
cross-validation directory and mark `artifact_feedback.status = candidate`.
"""


RUNBOOK = """\
# Motor validation lane runbook

For a Radia HDiv-MMM + HCurl eddy-bubble slot:

```powershell
python validation_test\\cln\\hcurl_vim_hdiv_mmm_end_to_end.py
python -m pytest tests\\test_vim_eddy_hybrid.py -q
```

Use the first command as the solver-ready mixed-system gate.  It must exercise
HDiv-MMM material response, the HCurl eddy-bubble basis, their coupling block,
and the shared mesh/material registry.  Attach the comparison with
`motor_validation_lane = "hdiv_mmm_hcurl_eddy_bubble"`; a non-empty
`solver_ready_artifact.verification` list is mandatory for MCP learning.

For an NGSolve+AGE slot:

```powershell
python -m pytest packages\\radia-mcp\\tests\\test_airgap_element.py `
  validation_test\\radia_mcp\\test_airgap_ngsolve_coupling.py `
  validation_test\\radia_mcp\\test_airgap_two_region.py validation_test\\radia_mcp\\test_airgap_machine_rotation.py `
  validation_test\\radia_mcp\\test_airgap_eddy_machine.py
```

Then attach the private/local solver comparison as an artifact with
`motor_validation_lane = "ngsolve_age"` and the relevant `age_gate_ids`.
"""


SECTIONS = {
    "overview": OVERVIEW,
    "source_policy": SOURCE_POLICY,
    "promotion_policy": PROMOTION_POLICY,
    "runbook": RUNBOOK,
}


def _lane_lines(lane: MotorValidationLane) -> list[str]:
    return [
        f"## `{lane.lane_id}`: {lane.label}",
        "",
        f"- radia path: {lane.radia_path}",
        f"- support status: `{lane.support_status}`",
        f"- support note: {lane.support_note}",
        "- best for:",
        *[f"  - {item}" for item in lane.best_for],
        "- observable families:",
        *[f"  - `{item}`" for item in lane.observable_families],
        "- required artifact fields:",
        *[f"  - `{item}`" for item in lane.required_fields],
        "- required metrics:",
        *[f"  - `{item}`" for item in lane.required_metrics],
        "- public-safe evidence:",
        *[f"  - `{item}`" for item in lane.public_evidence],
        "- private reference sources:",
        *[f"  - {item}" for item in lane.private_reference_sources],
        "- MCP promotion targets:",
        *[f"  - `{item}`" for item in lane.promotion_targets],
    ]


def format_motor_validation_lanes(topic: str = "overview") -> str:
    """Return Markdown documentation for the motor validation-lane contract."""
    t = topic.strip().lower()
    if t == "all":
        parts = [SECTIONS[key] for key in SECTIONS]
        parts.append(format_motor_validation_lanes("lane_matrix"))
        return "\n\n---\n\n".join(parts)
    if t == "lane_matrix":
        lines: list[str] = ["# Motor validation lane matrix", ""]
        for lane in LANES.values():
            lines.extend(_lane_lines(lane))
            lines.append("")
        return "\n".join(lines).rstrip()
    if t not in SECTIONS:
        valid = ", ".join(sorted(tuple(SECTIONS) + ("all", "lane_matrix")))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]


def lane_template(lane_id: str = "all") -> dict[str, Any]:
    """Return a JSON-serializable artifact template for one lane or all lanes."""
    requested = _canonical_lane_id(lane_id)
    if requested == "all":
        return {
            "schema_version": "radia-motor-validation-lanes/v1",
            "lanes": {key: lane_template(key) for key in LANES},
        }
    if requested not in LANES:
        return {
            "schema_version": "radia-motor-validation-lanes/v1",
            "error": f"unknown lane: {lane_id}",
            "valid_lanes": sorted(LANES),
        }
    lane = LANES[requested]
    return {
        "schema_version": "radia-motor-validation-artifact/v1",
        "motor_validation_lane": lane.lane_id,
        "label": lane.label,
        "support_status": lane.support_status,
        "support_note": lane.support_note,
        "reference_source_class": "product_local_reference | opensource_reference | analytic_reference | stored_regression",
        "observable_family": list(lane.observable_families),
        "public_evidence": list(lane.public_evidence),
        "required_fields": list(lane.required_fields),
        "required_metrics": list(lane.required_metrics),
        "timing_breakdown_s": {
            "setup": 0.0,
            "solve": 0.0,
            "postprocess": 0.0,
            "artifact_write": 0.0,
        },
        "artifact_feedback": {
            "status": "candidate | promoted",
            "promotion_target": list(lane.promotion_targets),
            "public_lesson": "state the generalized engineering lesson here",
        },
    }


def _as_artifact(value: Mapping[str, Any] | str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    parsed = json.loads(value)
    if not isinstance(parsed, Mapping):
        msg = "artifact JSON must decode to an object"
        raise TypeError(msg)
    return parsed


def _valid_solver_ready_verification(value: Any) -> bool:
    """Return True for a concrete, non-empty solver verification list."""
    if not isinstance(value, (list, tuple)) or not value:
        return False
    for item in value:
        if isinstance(item, str):
            if not item.strip():
                return False
            continue
        if isinstance(item, Mapping):
            if not any(
                str(item.get(key, "")).strip()
                for key in ("command", "artifact_id", "result", "test", "path")
            ):
                return False
            continue
        return False
    return True


def validate_motor_validation_artifact(
    artifact: Mapping[str, Any] | str,
    expected_lane: str = "",
) -> dict[str, Any]:
    """Check whether a cross-validation artifact can train a motor lane."""
    try:
        data = _as_artifact(artifact)
    except (TypeError, json.JSONDecodeError) as exc:
        return {
            "schema_version": "radia-motor-validation-artifact-gate/v1",
            "status": "fail",
            "errors": [str(exc)],
            "warnings": [],
        }

    errors: list[str] = []
    warnings: list[str] = []
    declared_lane_id = str(data.get("motor_validation_lane", "")).strip().lower()
    lane_id = _canonical_lane_id(declared_lane_id)
    if declared_lane_id in LEGACY_LANE_ALIASES:
        warnings.append(
            f"deprecated motor_validation_lane {declared_lane_id!r}; use {lane_id!r}"
        )
    expected_canonical = _canonical_lane_id(expected_lane) if expected_lane else ""
    if expected_canonical and lane_id != expected_canonical:
        errors.append(
            f"expected lane {expected_lane!r}, artifact declares {lane_id!r}"
        )
    if lane_id not in LANES:
        errors.append(f"unknown or missing motor_validation_lane: {lane_id!r}")
        lane = None
    else:
        lane = LANES[lane_id]

    if lane is not None:
        missing = [field for field in lane.required_fields if field not in data]
        errors.extend(f"missing required field: {field}" for field in missing)

        observable = str(data.get("observable_family", "")).strip().lower()
        if observable and observable not in lane.observable_families:
            errors.append(
                f"observable_family {observable!r} is not valid for lane {lane_id!r}"
            )
        elif not observable:
            errors.append("missing observable_family")

        metrics = data.get("metrics", {})
        if not isinstance(metrics, Mapping):
            errors.append("metrics must be an object")
            metrics = {}
        present_metrics = [metric for metric in lane.required_metrics if metric in metrics]
        if not present_metrics:
            errors.append(
                "metrics must contain at least one lane-required metric: "
                + ", ".join(lane.required_metrics)
            )

    status_value = str(data.get("status", "")).strip().lower()
    if status_value not in {"pass", "warn", "fail"}:
        errors.append("status must be one of: pass, warn, fail")
    elif status_value != "pass":
        warnings.append(f"artifact status is {status_value!r}; keep as candidate")

    timing = data.get("timing_breakdown_s", {})
    if not isinstance(timing, Mapping) or not timing:
        errors.append("timing_breakdown_s must be a non-empty object")

    coupling_status = str(data.get("coupling_design_status", "")).strip().lower()
    solver_artifact = data.get("solver_ready_artifact", {})
    solver_verification = (
        solver_artifact.get("verification") if isinstance(solver_artifact, Mapping) else None
    )
    has_solver_ready_artifact = (
        isinstance(solver_artifact, Mapping)
        and bool(str(solver_artifact.get("artifact_id", "")).strip())
        and _valid_solver_ready_verification(solver_verification)
    )
    if coupling_status in {"solver_validated", "validated_solver_path"} and not has_solver_ready_artifact:
        errors.append(
            "solver_ready_artifact must include artifact_id and a non-empty "
            "verification list when coupling_design_status claims solver validation"
        )
    if lane_id in LANES and not has_solver_ready_artifact:
        errors.append(
            f"{lane_id} requires solver_ready_artifact with artifact_id and "
            "non-empty execution verification"
        )
    shared_identity = data.get("shared_mesh_material_identity")
    if lane_id in LANES and not isinstance(shared_identity, Mapping):
        errors.append("shared_mesh_material_identity must be an object")
    elif isinstance(shared_identity, Mapping):
        for key in ("geometry_sha256", "material_sha256", "excitation_sha256"):
            digest = str(shared_identity.get(key, "")).lower()
            if len(digest) != 64 or any(
                character not in "0123456789abcdef" for character in digest
            ):
                errors.append(
                    f"shared_mesh_material_identity.{key} must be a sha256 digest"
                )
    if lane_id == "hdiv_mmm_hcurl_eddy_bubble":
        for contract_name in (
            "hdiv_mmm_operator_contract",
            "hcurl_eddy_bubble_contract",
            "coupling_operator_contract",
        ):
            contract = data.get(contract_name)
            if not isinstance(contract, Mapping) or not contract:
                errors.append(f"{contract_name} must be a non-empty object")

    feedback = data.get("artifact_feedback", {})
    if not isinstance(feedback, Mapping):
        errors.append("artifact_feedback must be an object")
    else:
        fb_status = str(feedback.get("status", "")).strip().lower()
        if fb_status not in {"candidate", "promoted"}:
            errors.append("artifact_feedback.status must be candidate or promoted")
        if not str(feedback.get("public_lesson", "")).strip():
            errors.append("artifact_feedback.public_lesson is required")

    result_status = "pass" if not errors else "fail"
    support_status = LANES[lane_id].support_status if lane_id in LANES else "unknown"
    validated_solver_path = (
        result_status == "pass"
        and status_value == "pass"
        and support_status == "supported_validation_path"
        and has_solver_ready_artifact
    )
    accepted_for_mcp_rfc_learning = (
        result_status == "pass"
        and status_value == "pass"
        and support_status == "experimental_rfc"
    )
    validated_experimental_solver_path = (
        accepted_for_mcp_rfc_learning
        and coupling_status in {"solver_validated", "validated_solver_path"}
        and has_solver_ready_artifact
    )
    validated_required_solver_path = (
        result_status == "pass"
        and status_value == "pass"
        and support_status == "required_validation_path"
        and has_solver_ready_artifact
    )
    accepted_for_mcp_learning = (
        result_status == "pass"
        and status_value == "pass"
        and not warnings
        and (
            validated_solver_path
            or validated_experimental_solver_path
            or validated_required_solver_path
        )
    )
    return {
        "schema_version": "radia-motor-validation-artifact-gate/v1",
        "status": result_status,
        "lane": lane_id,
        "support_status": support_status,
        "validated_solver_path": validated_solver_path,
        "validated_experimental_solver_path": validated_experimental_solver_path,
        "validated_required_solver_path": validated_required_solver_path,
        "validated_supported_path": validated_solver_path,
        "accepted_for_mcp_learning": accepted_for_mcp_learning,
        "accepted_for_mcp_rfc_learning": accepted_for_mcp_rfc_learning,
        "errors": errors,
        "warnings": warnings,
    }


def format_artifact_gate_result(result: Mapping[str, Any]) -> str:
    """Format a validation-artifact gate result as Markdown."""
    lines = [
        "# Motor validation artifact gate",
        "",
        f"- schema: `{result.get('schema_version', '')}`",
        f"- status: `{result.get('status', '')}`",
        f"- lane: `{result.get('lane', '')}`",
        f"- support status: `{result.get('support_status', '')}`",
        f"- validated solver path: `{result.get('validated_solver_path', False)}`",
        f"- validated experimental solver path: `{result.get('validated_experimental_solver_path', False)}`",
        f"- validated required solver path: `{result.get('validated_required_solver_path', False)}`",
        f"- validated supported path: `{result.get('validated_supported_path', False)}`",
        f"- accepted for MCP learning: `{result.get('accepted_for_mcp_learning', False)}`",
        f"- accepted for MCP RFC learning: `{result.get('accepted_for_mcp_rfc_learning', False)}`",
    ]
    errors = list(result.get("errors", ()))
    warnings = list(result.get("warnings", ()))
    if errors:
        lines.extend(["", "## Errors"])
        lines.extend(f"- {item}" for item in errors)
    if warnings:
        lines.extend(["", "## Warnings"])
        lines.extend(f"- {item}" for item in warnings)
    return "\n".join(lines).rstrip()
