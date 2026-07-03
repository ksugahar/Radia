"""Public-safe validation-lane contract for radia-motor.

radia-motor has several motor-validation paths that should not be collapsed
into one score:

* NGSolve+AGE for the finite-element air-gap machine path: torque, dq,
  eddy, and nonlinear machine quantities.
* 2D collocation MMMM for fast planar per-region soft-iron and torque-sweep
  checks.
* HDiv-VIM + reduced FEM as an experimental RFC for future rotor/source-field
  plus fixed-stator reduced response coupling.

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
            "nonlinear material iteration once the reduced invariant is clear",
        ),
        observable_families=(
            "airgap_flux",
            "torque",
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
            "tests/test_airgap_ngsolve_coupling.py",
            "tests/test_airgap_two_region.py",
            "tests/test_airgap_machine_rotation.py",
            "tests/test_airgap_eddy_machine.py",
            "tests/test_build123d_ipm_age_torque.py",
        ),
        private_reference_sources=(
            "product local reference",
            "lab-local finite-element reference",
            "open-source reference",
            "stored regression reference",
        ),
        promotion_targets=(
            "radia_mcp.motor.age_quality_knowledge",
            "radia_mcp.motor.simple_mmm_2d validation routing",
            "radia_mcp.radia_ngsolve AGE / force recipes",
        ),
    ),
    "mmmm2d_coarse": MotorValidationLane(
        lane_id="mmmm2d_coarse",
        label="2D collocation MMMM coarse motor lane",
        support_status="supported_coarse_path",
        support_note=(
            "This is now a verified coarse/reduced radia path for planar "
            "multi-region soft-iron MMMM checks. It supports per-region "
            "mu_r/BH inputs and factor-once torque sweeps, but it is not the "
            "full AGE rotating-machine finite-element path."
        ),
        radia_path="radia.mmmm2d dense 2D collocation moment solver",
        best_for=(
            "fast planar soft-iron sanity checks before AGE",
            "multi-grade rotor/stator region experiments",
            "factor-once torque-angle sweeps for reduced motor studies",
            "coarse optimization loops where internal loop pollution is acceptable",
        ),
        observable_families=(
            "torque",
            "coenergy",
            "force_or_torque_trend",
            "per_region_magnetization",
            "demag_field",
        ),
        required_fields=COMMON_REQUIRED_FIELDS
        + (
            "mmmm2d_contract",
            "region_material_contract",
            "pytest_targets",
        ),
        required_metrics=(
            "torque_relative_error",
            "m_avg_relative_error",
            "region_magnetization_ratio",
            "quantity_specific_residual",
        ),
        public_evidence=(
            "validation_test/feec/test_moment2d_perregion.py",
        ),
        private_reference_sources=(
            "product local reference",
            "lab-local finite-element reference",
            "open-source reference",
            "stored regression reference",
        ),
        promotion_targets=(
            "radia_mcp.motor.validation_lanes_knowledge",
            "radia_mcp.motor.triple_check_knowledge",
            "radia_mcp.radia_ngsolve knowledge/mmm_core.py",
        ),
    ),
    "hdiv_vim_reduced_fem": MotorValidationLane(
        lane_id="hdiv_vim_reduced_fem",
        label="HDiv-VIM + reduced FEM (experimental RFC)",
        support_status="experimental_rfc",
        support_note=(
            "This is a new research coupling idea, not the historical "
            "radia-motor supported path. Treat HDiv-VIM rotor plus reduced-FEM "
            "stator as a design proposal until an interface operator, reduced "
            "basis, and regression artifact are implemented."
        ),
        radia_path="proposed radia.vim HDiv rotor source-field lane plus fixed-stator reduced FEM",
        best_for=(
            "passive pickup flux and signed flux-linkage sweeps",
            "permanent-magnet demagnetizing-field anchors",
            "source-field / surface-current intuition",
            "researching whether a rotor VIM source can drive a compact fixed-stator reduced FEM response",
        ),
        observable_families=(
            "pickup_flux",
            "flux_linkage",
            "demag_field",
            "coenergy",
            "force_or_torque_trend",
        ),
        required_fields=COMMON_REQUIRED_FIELDS
        + (
            "coupling_design_status",
            "interface_operator_contract",
            "reduced_fem_contract",
            "vim_operator_contract",
        ),
        required_metrics=(
            "signed_agreement_count",
            "mean_abs_relative_error",
            "rms_abs_relative_error",
            "max_abs_relative_error",
        ),
        public_evidence=(
            "analytic sign/scale checks",
            "stored public-safe regression artifacts",
            "reduced FEM consistency checks",
        ),
        private_reference_sources=(
            "product local reference",
            "lab-local finite-element reference",
            "open-source reference",
            "stored regression reference",
        ),
        promotion_targets=(
            "radia_mcp.motor.validation_lanes_knowledge",
            "radia_mcp.motor.simple_mmm_2d prompt triage text",
            "radia_mcp.radia_ngsolve force / flux recipes when applicable",
        ),
    ),
}


OVERVIEW = """\
# radia-motor validation lanes

radia-motor should keep independent cross-validation lanes:

- `ngsolve_age`: NGSolve+AGE.  This is the current supported radia-motor
  finite-element air-gap machine lane for torque, dq quantities, cogging,
  eddy/slip, hysteresis, and nonlinear machine studies.
- `mmmm2d_coarse`: 2D collocation MMMM.  This is a supported coarse/reduced
  lane for planar multi-region soft iron, per-region material dictionaries,
  and factor-once torque sweeps.  It is useful before AGE and for optimization
  triage, but it is not the full AGE moving-air-gap path.
- `hdiv_vim_reduced_fem`: HDiv-VIM plus reduced FEM.  This is an experimental
  RFC lane.  The idea of using HDiv-VIM for the rotor and a reduced FEM model
  for the fixed stator is new and intentionally unusual; do not describe it as
  supported until a coupling/interface operator and reduced-basis regression
  pass.

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

1. `which_lane`: `ngsolve_age`, `mmmm2d_coarse`, or `hdiv_vim_reduced_fem`.
2. `which_observable`: one lane-supported observable family.
3. `which_promotion`: the exact public-safe knowledge or recipe that improved.

Promotion is allowed when:

- the artifact passes its own tolerance,
- the lane metadata is complete,
- the observable belongs to that lane,
- the timing breakdown has at least one named phase, and
- the lesson can be stated without private paths, product logs, or commercial
  benchmark numbers.

If the slot only produced a useful private comparison, keep it in the private
cross-validation directory and mark `artifact_feedback.status = candidate`.
"""


RUNBOOK = """\
# Motor validation lane runbook

For a 2D MMMM coarse motor slot:

```powershell
python -m pytest validation_test\\feec\\test_moment2d_perregion.py -q
```

Then attach the comparison as `motor_validation_lane = "mmmm2d_coarse"` with
`mmmm2d_contract`, `region_material_contract`, and the pytest target. This lane
can be used as a verified coarse/reduced path, not as a replacement for AGE.

For an HDiv-VIM + reduced FEM research slot:

```powershell
python -m pytest tests\\test_loop_slot_gates.py -k hdiv
```

Then attach the private/local comparison as a research artifact with
`motor_validation_lane = "hdiv_vim_reduced_fem"` and
`coupling_design_status = "experimental_rfc"`.  Passing this metadata gate
means the idea is organized for learning; it does not mean radia-motor already
supports the coupled solver.

For an NGSolve+AGE slot:

```powershell
python -m pytest tests\\test_airgap_element.py tests\\test_airgap_ngsolve_coupling.py `
  tests\\test_airgap_two_region.py tests\\test_airgap_machine_rotation.py `
  tests\\test_airgap_eddy_machine.py
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
    requested = lane_id.strip().lower()
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
    lane_id = str(data.get("motor_validation_lane", "")).strip().lower()
    if expected_lane and lane_id != expected_lane.strip().lower():
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
    )
    validated_coarse_path = (
        result_status == "pass"
        and status_value == "pass"
        and support_status == "supported_coarse_path"
    )
    return {
        "schema_version": "radia-motor-validation-artifact-gate/v1",
        "status": result_status,
        "lane": lane_id,
        "support_status": support_status,
        "validated_solver_path": validated_solver_path,
        "validated_coarse_path": validated_coarse_path,
        "validated_supported_path": validated_solver_path or validated_coarse_path,
        "accepted_for_mcp_learning": result_status == "pass" and status_value == "pass",
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
        f"- validated coarse path: `{result.get('validated_coarse_path', False)}`",
        f"- validated supported path: `{result.get('validated_supported_path', False)}`",
        f"- accepted for MCP learning: `{result.get('accepted_for_mcp_learning', False)}`",
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
