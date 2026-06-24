"""Public-safe NGSolve AGE quality gates for radia-motor.

The goal of this module is not to claim turnkey commercial-solver parity.
It gives the MCP server a precise, repeatable definition of what "ready"
means for the NGSolve Air-Gap Harmonic Element (AGE) motor path:

* which physical quantities are reduced,
* which public tests gate them,
* which motor families they support, and
* which limitations must still be stated before publishing examples.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class AgeGate:
    """One public-safe quality gate for the AGE motor path."""

    gate_id: str
    level: str
    quantity: str
    evidence: tuple[str, ...]
    acceptance: tuple[str, ...]
    families: tuple[str, ...]


GATES: tuple[AgeGate, ...] = (
    AgeGate(
        gate_id="age_analytic_dtn",
        level="gold_age_invariant",
        quantity="annular DtN matrix, harmonic transfer, radius-independent torque",
        evidence=("tests/test_airgap_element.py",),
        acceptance=(
            "closed-form torque is zero for in-phase phasors",
            "torque sign flips with relative phase",
            "radius and mesh are not part of the torque readout",
        ),
        families=("spm", "ipm", "synrm", "srm", "induction", "hysteresis"),
    ),
    AgeGate(
        gate_id="age_ngsolve_robin",
        level="gold_age_invariant",
        quantity="single-boundary AGE Robin coupling against analytic annulus field",
        evidence=("tests/test_airgap_ngsolve_coupling.py",),
        acceptance=("field relative error below 5e-3 against full-annulus analytic solution",),
        families=("spm", "ipm", "synrm", "srm"),
    ),
    AgeGate(
        gate_id="age_two_region_coupling",
        level="gold_age_invariant",
        quantity="rotor FE region plus stator FE region coupled across an unmeshed gap",
        evidence=("tests/test_airgap_two_region.py",),
        acceptance=("two-region field relative error below 1e-3 against 3-region analytic solution",),
        families=("spm", "ipm", "synrm", "srm", "induction"),
    ),
    AgeGate(
        gate_id="age_multiharmonic_gap",
        level="gold_age_invariant",
        quantity="multi-harmonic dense gap block for slotting and non-sinusoidal fields",
        evidence=("tests/test_airgap_multiharmonic.py",),
        acceptance=("multi-harmonic field relative error below 1e-4 against superposed analytic modes",),
        families=("spm", "ipm", "synrm", "srm", "hysteresis"),
    ),
    AgeGate(
        gate_id="age_rotation_torque",
        level="gold_age_invariant",
        quantity="rotor rotation as harmonic phase plus closed-form air-gap torque",
        evidence=("tests/test_airgap_machine_rotation.py",),
        acceptance=(
            "torque relative error below 2e-3 against analytic rotating-field torque",
            "torque swings with rotor phase, proving this is not a static-only check",
        ),
        families=("spm", "ipm", "synrm", "srm", "hysteresis"),
    ),
    AgeGate(
        gate_id="age_eddy_gap",
        level="gold_age_invariant",
        quantity="complex AGE gap with conducting region and eddy phase",
        evidence=("tests/test_airgap_eddy_coupling.py",),
        acceptance=(
            "complex field relative error below 1e-2 against fully meshed gap reference",
            "imaginary field component is non-trivial, so the eddy phase is actually tested",
        ),
        families=("induction", "wound_field", "loss"),
    ),
    AgeGate(
        gate_id="age_eddy_machine",
        level="gold_age_invariant",
        quantity="complete rotating eddy machine: rotor, stator conductor, gap, rotation, torque",
        evidence=("tests/test_airgap_eddy_machine.py",),
        acceptance=(
            "conducting-stator field relative error below 1e-4 against fully meshed reference",
            "rotation torque relative error below 1e-2",
            "strong complex eddy phase is present",
        ),
        families=("induction", "loss"),
    ),
    AgeGate(
        gate_id="age_spm_field",
        level="silver_age_field",
        quantity="SPM magnetization labels to AGE air-gap harmonic spectrum",
        evidence=("tests/test_build123d_pmsm_field.py",),
        acceptance=("radial alternating N/S magnets create the requested pole harmonic in the AGE field",),
        families=("spm", "afpm"),
    ),
    AgeGate(
        gate_id="age_ipm_synchronous_torque",
        level="silver_age_torque",
        quantity="IPM/SPM synchronous torque and torque-angle law",
        evidence=("tests/test_build123d_ipm_age_torque.py",),
        acceptance=(
            "phase-locked torque ripple below 1e-2",
            "torque-angle residual below 5e-3 against a single-sine law for the baseline case",
        ),
        families=("spm", "ipm"),
    ),
    AgeGate(
        gate_id="dq_control_layer",
        level="gold_reduced_invariant",
        quantity="Ke/Kt, dq torque, Ld/Lq, MTPA, field-weakening operating regions",
        evidence=(
            "tests/test_pm_emf_constants.py",
            "tests/test_dq_torque.py",
            "tests/test_motor_mtpa.py",
            "tests/test_field_weakening.py",
            "tests/test_dq_operating_point.py",
        ),
        acceptance=(
            "back-EMF constants match the dq voltage equation",
            "MTPA closed form matches independent numeric current-angle maximization",
            "field-weakening region selection matches independent numeric argmax",
        ),
        families=("spm", "ipm", "synrm"),
    ),
    AgeGate(
        gate_id="cogging_reluctance_layer",
        level="silver_age_torque",
        quantity="cogging order, skew cancellation, reluctance torque sign and periodicity",
        evidence=(
            "tests/test_cogging_order.py",
            "tests/test_skew_factor.py",
            "tests/test_skew_average.py",
            "tests/test_motor_cogging_torque.py",
            "tests/test_machine_scaling.py",
        ),
        acceptance=(
            "cogging order is LCM(slots, poles)",
            "one-slot-pitch skew cancels the slot-passing/cogging harmonic",
            "reluctance torque has the expected zero mean and sin(2 theta) signature",
        ),
        families=("spm", "ipm", "synrm", "srm"),
    ),
    AgeGate(
        gate_id="induction_slip_layer",
        level="silver_age_eddy",
        quantity="induction slip coupling, torque-slip curve, deep-bar AC factors",
        evidence=(
            "tests/test_airgap_eddy_machine.py",
            "tests/test_motor_induction_coupling.py",
            "tests/test_induction_machine.py",
            "tests/test_deep_bar.py",
            "tests/test_dowell.py",
        ),
        acceptance=(
            "rotor flux linkage screens with increasing slip",
            "rotor eddy loss is positive and rises from near zero",
            "Thevenin torque-slip breakdown has a finite maximum",
            "deep-bar resistance/reactance factors satisfy the closed-form limits",
        ),
        families=("induction",),
    ),
    AgeGate(
        gate_id="hysteresis_loss_layer",
        level="silver_age_loss",
        quantity="stateful B-input Play hysteresis loss from an AGE rotor sweep",
        evidence=(
            "tests/test_hysteresis_play.py",
            "tests/test_hysteresis_fe_coupled.py",
            "tests/test_hysteresis_fe_variational.py",
            "tests/test_hysteresis_motor_loss.py",
        ),
        acceptance=(
            "per-point loop area is positive",
            "motor iron-loss wattage stays finite and physical",
            "rotating/elliptical B path exceeds scalar peak-B loss estimate",
        ),
        families=("hysteresis", "spm", "ipm", "loss"),
    ),
)


FAMILY_GATES: dict[str, tuple[str, ...]] = {
    "spm": (
        "age_analytic_dtn",
        "age_multiharmonic_gap",
        "age_rotation_torque",
        "age_spm_field",
        "age_ipm_synchronous_torque",
        "dq_control_layer",
        "cogging_reluctance_layer",
    ),
    "ipm": (
        "age_analytic_dtn",
        "age_multiharmonic_gap",
        "age_rotation_torque",
        "age_ipm_synchronous_torque",
        "dq_control_layer",
        "cogging_reluctance_layer",
    ),
    "induction": (
        "age_analytic_dtn",
        "age_two_region_coupling",
        "age_eddy_gap",
        "age_eddy_machine",
        "induction_slip_layer",
    ),
    "srm": (
        "age_analytic_dtn",
        "age_two_region_coupling",
        "age_multiharmonic_gap",
        "age_rotation_torque",
        "cogging_reluctance_layer",
    ),
    "synrm": (
        "age_analytic_dtn",
        "age_two_region_coupling",
        "age_multiharmonic_gap",
        "age_rotation_torque",
        "dq_control_layer",
        "cogging_reluctance_layer",
    ),
    "hysteresis": (
        "age_analytic_dtn",
        "age_multiharmonic_gap",
        "age_rotation_torque",
        "hysteresis_loss_layer",
    ),
}


FAMILY_QUANTITIES: dict[str, tuple[str, ...]] = {
    "spm": (
        "air-gap harmonic spectrum",
        "PM flux linkage and back-EMF constant",
        "cogging order and skew attenuation",
        "phase-locked torque-angle law",
        "Ld/Lq and MTPA control checks where saliency is present",
    ),
    "ipm": (
        "PM flux linkage and Ke/Kt",
        "Ld/Lq saliency",
        "PM plus reluctance torque split",
        "MTPA current angle",
        "field-weakening voltage ellipse",
        "cogging order and torque periodicity",
    ),
    "induction": (
        "slip-frequency complex flux linkage",
        "rotor eddy loss",
        "Thevenin torque-slip curve and breakdown point",
        "deep-bar resistance/reactance trend",
        "complex AGE field and torque against a fully meshed reference",
    ),
    "srm": (
        "reluctance torque sign",
        "angle-current torque periodicity",
        "nonlinear inductance trend",
        "sector/cogging symmetry",
    ),
    "synrm": (
        "Ld/Lq saliency",
        "pure-reluctance MTPA angle",
        "synchronous power-angle torque",
        "flux-barrier/cogging periodicity",
    ),
    "hysteresis": (
        "B(theta) path from AGE rotor sweep",
        "stateful loop area per point",
        "positive iron-loss wattage",
        "rotating-field loss excess over scalar peak-B estimate",
    ),
}


def _gate_by_id() -> dict[str, AgeGate]:
    return {gate.gate_id: gate for gate in GATES}


def _infer_family(goal: str) -> str:
    g = goal.lower()
    if any(term in g for term in ("induction", " cage", " im ", "slip", "deep bar", "rotor bar")):
        return "induction"
    if any(term in g for term in ("srm", "switched reluctance", "sr motor", "doubly salient")):
        return "srm"
    if any(term in g for term in ("synrm", "synchronous reluctance", "flux barrier", "reluctance motor")):
        return "synrm"
    if any(term in g for term in ("ipm", "interior", "buried magnet", "hairpin")):
        return "ipm"
    if "hysteresis" in g:
        return "hysteresis"
    return "spm"


def _format_gate(gate: AgeGate) -> str:
    evidence = ", ".join(f"`{item}`" for item in gate.evidence)
    accept = "; ".join(gate.acceptance)
    families = ", ".join(f"`{item}`" for item in gate.families)
    return (
        f"| `{gate.gate_id}` | `{gate.level}` | {gate.quantity} | "
        f"{evidence} | {accept} | {families} |"
    )


OVERVIEW = """\
# NGSolve AGE motor quality: readiness layer

The NGSolve Air-Gap Harmonic Element (AGE) is the main radia-motor solve
path for 2D rotating machines. It is stronger than a prompt-time magnetic
circuit check because it couples rotor and stator FE regions across an
unmeshed analytic air gap, rotates the rotor as a harmonic phase, supports
complex eddy-current regions, and reads torque from closed-form gap phasors.

The public quality policy is:

- Use AGE as the authoritative open validation path for 2D SPM/IPM/IM/SRM/SynRM
  reduced quantities.
- Use the small MMM quick check only for sign and scale triage.
- Publish examples only when their advertised quantity is covered by a passing
  public gate below.
- Do not publish commercial solver logs, private paths, or product-run numeric
  references. Public AGE gates are analytic, NGSolve self-consistency, or
  fully-meshed open reference comparisons.

This is the closest public-safe analogue to product-solver discipline: every
claim must name the physical quantity, the independent gate, and the limitation.
"""


GATE_MATRIX = """\
# AGE gate matrix

| Gate | Level | Physical quantity | Evidence | Acceptance | Families |
|---|---|---|---|---|---|
""" + "\n".join(_format_gate(gate) for gate in GATES)


FAMILY_MATRIX = """\
# Motor-family AGE validation matrix

| Family | Required AGE gates | Primary physical quantities |
|---|---|---|
""" + "\n".join(
    f"| `{family}` | "
    f"{', '.join(f'`{gate}`' for gate in gates)} | "
    f"{'; '.join(FAMILY_QUANTITIES[family])} |"
    for family, gates in FAMILY_GATES.items()
)


PUBLICATION_POLICY = """\
# Publication policy for radia-motor AGE examples

An example can be described as AGE-verified only when:

1. Its motor family maps to a concrete gate set in `motor_age_validation_plan`.
2. The advertised quantity is one of the gate quantities, not a vague "motor
   simulation works" claim.
3. The public test evidence is green locally.
4. The result text reports limitations: 2D model, no automatic end effects,
   no turnkey drive circuit, and no hidden commercial reference.

Suggested labels:

- `gold_age_invariant`: analytic or fully-meshed open reference with explicit
  tolerance.
- `silver_age_torque`: AGE field solve plus physical torque/sign/periodicity
  invariant.
- `silver_age_eddy`: complex eddy-current trend or torque-slip layer, with a
  public finite-element or closed-form consistency gate.
- `silver_age_loss`: hysteresis/core-loss invariant from a public stateful
  material model.
- `prompt_triage_only`: MMM quick-check output; not a publishable validation
  label by itself.
"""


RUNBOOK = """\
# AGE quality runbook

Fast algebra/control gates:

```powershell
python -m pytest tests\\test_airgap_element.py tests\\test_age_winding_factor.py `
  tests\\test_pm_emf_constants.py tests\\test_dq_torque.py tests\\test_motor_mtpa.py `
  tests\\test_field_weakening.py tests\\test_induction_machine.py tests\\test_deep_bar.py `
  tests\\test_cogging_order.py tests\\test_skew_factor.py tests\\test_skew_average.py
```

AGE finite-element gates:

```powershell
python -m pytest tests\\test_airgap_ngsolve_coupling.py tests\\test_airgap_two_region.py `
  tests\\test_airgap_multiharmonic.py tests\\test_airgap_machine_rotation.py `
  tests\\test_airgap_eddy_coupling.py tests\\test_airgap_eddy_machine.py
```

Motor archetype gates:

```powershell
python -m pytest tests\\test_build123d_pmsm_field.py `
  tests\\test_build123d_ipm_age_torque.py tests\\test_motor_cogging_torque.py `
  tests\\test_motor_induction_coupling.py tests\\test_hysteresis_motor_loss.py
```

MCP gates:

```powershell
$env:PYTHONPATH = ".\\src"
python -m radia_mcp.motor.server --selftest
python -m pytest tests\\test_each_server_selftest.py -k motor
```

Use the family plan first, then run only the relevant subset when iterating.
Run the broader list before tagging or publishing.
"""


LIMITATIONS = """\
# AGE limitations that must stay visible

AGE is a high-quality open 2D validation backbone, but it is not a complete
turnkey motor product by itself.

State these limitations when relevant:

- 2D cross-section first; end winding, end ring, skew, and 3D end-region effects
  need explicit correction, multi-slice, or 3D modeling.
- The current nonlinear workflow is Picard/fixed point first. A production
  Newton tangent path is not yet the default.
- Circuit/drive co-simulation is represented by reduced dq or equivalent-circuit
  gates unless a coupled study is explicitly built.
- Hysteresis and lamination losses are available as public models, but material
  calibration remains a separate validation task.
- MMM quick checks are prompt-time triage only; AGE or another open gate must
  carry any public validation claim.
"""


SECTIONS = {
    "overview": OVERVIEW,
    "gate_matrix": GATE_MATRIX,
    "family_matrix": FAMILY_MATRIX,
    "publication_policy": PUBLICATION_POLICY,
    "runbook": RUNBOOK,
    "limitations": LIMITATIONS,
}


def route_age_validation_plan(goal: str) -> dict[str, Any]:
    """Return a structured AGE validation plan for a motor prompt."""
    family = _infer_family(goal)
    by_id = _gate_by_id()
    gate_ids = FAMILY_GATES[family]
    gates = [by_id[gate_id] for gate_id in gate_ids]
    tests = []
    for gate in gates:
        for item in gate.evidence:
            if item not in tests:
                tests.append(item)
    quantities = FAMILY_QUANTITIES[family]
    return {
        "schema_version": "radia-motor-age-validation/v1",
        "goal": goal,
        "family": family,
        "required_gate_ids": list(gate_ids),
        "primary_quantities": list(quantities),
        "pytest_targets": tests,
        "quality_labels": sorted({gate.level for gate in gates}),
        "workflow": [
            "Use motor_mmm_quick_check only for prompt-time sign/scale triage.",
            "Run the AGE/public reduced gates for the advertised physical quantity.",
            "Upgrade the example label only to the strongest passed gate.",
            "Keep commercial solver provenance, private paths, and raw product references out of public text.",
        ],
    }


def format_age_validation_plan(plan: dict[str, Any]) -> str:
    """Format an AGE validation plan as Markdown."""
    lines = [
        "# NGSolve AGE validation plan",
        "",
        f"- schema: `{plan['schema_version']}`",
        f"- goal: {plan['goal']}",
        f"- inferred family: `{plan['family']}`",
        f"- quality labels: {', '.join(f'`{label}`' for label in plan['quality_labels'])}",
        "",
        "## Required Gates",
    ]
    lines.extend(f"- `{gate}`" for gate in plan["required_gate_ids"])
    lines.extend(["", "## Physical Quantities"])
    lines.extend(f"- {quantity}" for quantity in plan["primary_quantities"])
    lines.extend(["", "## Pytest Targets"])
    lines.extend(f"- `{target}`" for target in plan["pytest_targets"])
    lines.extend(["", "## Workflow"])
    for i, step in enumerate(plan["workflow"], 1):
        lines.append(f"{i}. {step}")
    return "\n".join(lines).rstrip()


def get_age_quality_report(topic: str = "overview") -> str:
    """Return AGE motor quality documentation."""
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[key] for key in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
