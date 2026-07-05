"""Public-safe dual-lane motor training catalog.

The catalog is fed by source-native motor examples and readable external
teaching material, but it keeps only scrubbed engineering lessons.  Public
clients should see which radia-motor lane to train, not where a private lesson
came from.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any


FORBIDDEN_PUBLIC_MARKERS = (
    "J" + "MAG",
    "J" + "AC",
    "J" + "FT",
    "C" + "OMSOL",
    "F" + "EMM",
    "C" + "ST",
    "A" + "NSYS",
    "S" + ":" + "\\",
    "W" + ":" + "\\",
    "_cross" + "val",
)


@dataclass(frozen=True)
class DualLaneTrainingCase:
    """One scrubbed motor lesson for AGE and VIM training."""

    case_id: str
    family: str
    title: str
    source_seed_class: str
    learning_axis: str
    age_targets: tuple[str, ...]
    vim_targets: tuple[str, ...]
    observable_family: str
    validation_focus: str
    teaching_gate: str

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["radia_motor_age"] = {
            "lane_id": "radia-motor-age",
            "validation_lane": "ngsolve_age",
            "targets": list(self.age_targets),
        }
        data["radia_motor_vim"] = {
            "lane_id": "radia-motor-vim",
            "validation_lane": "hdiv_vim_reduced_fem",
            "targets": list(self.vim_targets),
        }
        return data


CATALOG: tuple[DualLaneTrainingCase, ...] = (
    DualLaneTrainingCase(
        "nonlinear_c_core_flux_linkage",
        "magnetic_circuit",
        "Nonlinear C-core flux-linkage and coenergy closure",
        "external_readable_machine_fem_reference",
        "Newton residual, nonlinear B-H material, and FLUM-like pickup consistency",
        ("nonlinear_material_iteration", "flux_linkage", "coenergy"),
        ("pickup_flux", "coenergy", "source_field"),
        "flux_linkage",
        "compare current scaling, energy/coenergy ordering, and residual convergence",
        "nonlinear_flux_linkage_energy_gate",
    ),
    DualLaneTrainingCase(
        "srm_static_torque_curve",
        "srm",
        "SRM static torque curve from aligned to unaligned positions",
        "external_readable_machine_fem_reference",
        "reluctance torque from coenergy gradient and pole overlap",
        ("reluctance_torque", "saturating_inductance", "torque_angle"),
        ("force_or_torque_trend", "coenergy", "pickup_flux"),
        "torque",
        "check sign, periodicity, positive overlap torque, and near-zero absolute tolerance",
        "coenergy_static_torque_curve_gate",
    ),
    DualLaneTrainingCase(
        "induction_airgap_br_waveform",
        "induction",
        "Induction-machine radial air-gap flux-density waveform",
        "external_readable_machine_fem_reference",
        "three-phase winding excitation, cage screening, and air-gap sampling",
        ("induction_machine", "airgap_eddy_machine", "slip_loss"),
        ("pickup_flux", "source_field", "reduced_fem_response"),
        "airgap_flux",
        "preserve angle grid, component basis, and phase-current metadata before comparing Br",
        "airgap_flux_waveform_metadata_gate",
    ),
    DualLaneTrainingCase(
        "spmsm_pm_airgap_waveform",
        "spm",
        "Surface-PM air-gap waveform with magnet and gap refinement",
        "external_readable_machine_fem_reference",
        "permanent-magnet polarity, slot-opening ripple, and mesh-density sensitivity",
        ("back_emf", "cogging_torque", "pm_airgap_flux"),
        ("demag_field", "pickup_flux", "source_field"),
        "airgap_flux",
        "check magnet polarity reversal, air-gap Br waveform phase, and cogging order",
        "pm_airgap_waveform_refinement_gate",
    ),
    DualLaneTrainingCase(
        "spm_torque_angle_sweep",
        "spm",
        "SPM torque-angle sweep with current phase alignment",
        "private_source_native_motor_workflow",
        "torque table metadata before value-level comparison",
        ("age_rotation_torque", "torque_angle", "dq_control_layer"),
        ("pickup_flux", "force_or_torque_trend", "coenergy"),
        "torque",
        "lock angle unit, phase convention, symmetry factor, and endpoint policy",
        "torque_angle_table_export_health",
    ),
    DualLaneTrainingCase(
        "ipm_saliency_mtpa",
        "ipm",
        "IPM saliency, Ld/Lq, MTPA, and demagnetization margin",
        "private_source_native_motor_workflow",
        "separate magnet torque from reluctance torque before optimization",
        ("ld_lq", "mtpa", "field_weakening", "demag_margin"),
        ("demag_field", "flux_linkage", "pickup_flux"),
        "ld_lq",
        "check dq axis convention, current basis, voltage limit, and d-axis demag margin",
        "dq_saliency_mtpa_metadata_gate",
    ),
    DualLaneTrainingCase(
        "induction_skew_multislice",
        "induction",
        "Induction-machine skew as slice-averaged torque and loss",
        "private_source_native_motor_workflow",
        "multi-slice skew must average independent 2D slices without losing cage coupling metadata",
        ("induction_machine", "airgap_eddy_machine", "deep_bar"),
        ("source_field", "reduced_fem_response", "force_or_torque_trend"),
        "slip_loss",
        "check skew angle, slice count, slip frequency, cage-bar metadata, and torque averaging",
        "skew_slice_torque_average_gate",
    ),
    DualLaneTrainingCase(
        "force_balance_maxwell_nodal",
        "force",
        "Force balance between stress integration and virtual-work nodal force",
        "private_source_native_motor_workflow",
        "force method selection and component-frame metadata",
        ("maxwell_stress_force", "coenergy_force", "force_report"),
        ("force_or_torque_trend", "source_field", "coenergy"),
        "force_or_torque_trend",
        "check force unit, component frame, integration contour, and action-reaction balance",
        "force_report_method_metadata_gate",
    ),
    DualLaneTrainingCase(
        "harmonic_current_loss_control",
        "drive_loss",
        "Current-feedback harmonic loss and iron-loss bucket separation",
        "private_source_native_motor_workflow",
        "drive waveform metadata must survive before loss buckets are trusted",
        ("hysteresis_loss", "eddy_loss", "dq_control_layer"),
        ("pickup_flux", "reduced_fem_response", "coenergy"),
        "hysteresis_loss",
        "check carrier frequency, current basis, loss buckets, and power-balance residual",
        "loss_bucket_efficiency_gate",
    ),
    DualLaneTrainingCase(
        "variable_magnet_demag_margin",
        "pm_demag",
        "Variable permanent-magnet demagnetization margin",
        "private_source_native_motor_workflow",
        "irreversible Br update needs an operating-point package, not only a scalar field",
        ("demag_margin", "pm_operating_point", "nonlinear_material_iteration"),
        ("demag_field", "source_field", "pickup_flux"),
        "demag_field",
        "check H along magnetization, recoil slope, knee margin, and update provenance",
        "pm_demag_operating_point_gate",
    ),
    DualLaneTrainingCase(
        "parametric_rotor_template",
        "geometry",
        "Parametric rotor geometry template and label contract",
        "private_source_native_motor_workflow",
        "geometry parameters must produce stable region names and boundary tags",
        ("geometry_label_contract", "mesh_region_mapping", "age_rotation_torque"),
        ("source_field", "interface_operator_contract", "reduced_fem_response"),
        "airgap_flux",
        "check region labels, pole-pair symmetry, air-gap clearance, and mesh-zone identity",
        "parametric_motor_geometry_label_gate",
    ),
    DualLaneTrainingCase(
        "arbitrary_motion_table",
        "motion",
        "Arbitrary 6-DOF motion table for actuator or rotor trajectories",
        "private_source_native_motor_workflow",
        "motion interpolation and coordinate-frame metadata before transient comparison",
        ("motion_table", "torque_angle", "airgap_flux"),
        ("force_or_torque_trend", "source_field", "reduced_fem_response"),
        "force_or_torque_trend",
        "check time unit, coordinate axes, rotation order, interpolation, and endpoint policy",
        "motion_table_coordinate_gate",
    ),
    DualLaneTrainingCase(
        "stranded_litz_ac_loss",
        "winding_loss",
        "Stranded or litz conductor AC loss trend",
        "private_source_native_motor_workflow",
        "skin/proximity reduced model before fine-strand meshing",
        ("winding_ac_loss", "dowell_factor", "hysteresis_loss"),
        ("pickup_flux", "reduced_fem_response", "coenergy"),
        "slip_loss",
        "check conductor radius, strand count, frequency, resistivity, and loss sign",
        "stranded_ac_loss_scaling_gate",
    ),
    DualLaneTrainingCase(
        "one_turn_skin_effect",
        "winding_loss",
        "One-turn conductor skin-effect frequency sweep",
        "private_source_native_motor_workflow",
        "AC/DC resistance ratio and skin-depth scaling",
        ("skin_depth", "winding_ac_loss", "eddy_loss"),
        ("pickup_flux", "source_field", "reduced_fem_response"),
        "slip_loss",
        "check R_ac/R_dc monotonicity, frequency units, and conductor conductivity",
        "skin_depth_resistance_gate",
    ),
    DualLaneTrainingCase(
        "arago_disk_eddy_drag",
        "moving_conductor",
        "Moving-conductor eddy drag in a rotating disk",
        "private_source_native_motor_workflow",
        "velocity-dependent drag and phase-lag sanity before full transient solve",
        ("moving_conductor_eddy", "eddy_drag", "airgap_eddy_machine"),
        ("source_field", "force_or_torque_trend", "coenergy"),
        "force_or_torque_trend",
        "check low-speed linear drag, high-speed saturation trend, and power loss sign",
        "moving_conductor_eddy_drag_gate",
    ),
    DualLaneTrainingCase(
        "eddy_current_brake",
        "moving_conductor",
        "Eddy-current brake lift/drag and conductive-plate loss",
        "private_source_native_motor_workflow",
        "drag force and Joule-loss must share the same velocity and conductivity package",
        ("eddy_drag", "airgap_eddy_machine", "loss_bucket"),
        ("source_field", "force_or_torque_trend", "reduced_fem_response"),
        "force_or_torque_trend",
        "check drag sign, lift sign convention, loss positivity, and velocity metadata",
        "eddy_brake_force_loss_gate",
    ),
    DualLaneTrainingCase(
        "cylinder_magnet_demag",
        "pm_demag",
        "Cylindrical magnet self-demagnetizing-field benchmark",
        "private_source_native_motor_workflow",
        "demag factor as a source-field sanity gate",
        ("pm_operating_point", "demag_margin", "analytic_demag_factor"),
        ("demag_field", "source_field", "pickup_flux"),
        "demag_field",
        "check aspect ratio, magnetization axis, recoil permeability, and internal-field sign",
        "cylindrical_magnet_demag_factor_gate",
    ),
    DualLaneTrainingCase(
        "uniform_field_demag",
        "demag",
        "Permeable or conducting body in uniform external field",
        "private_source_native_motor_workflow",
        "uniform-field response as a scalar demag sanity check",
        ("uniform_field", "analytic_demag_factor", "nonlinear_material_iteration"),
        ("demag_field", "source_field", "coenergy"),
        "demag_field",
        "check applied-field unit, material slope, internal-field ratio, and symmetry",
        "uniform_field_demag_ratio_gate",
    ),
    DualLaneTrainingCase(
        "quadrupole_force_multipole",
        "multipole",
        "Four-pole electromagnet multipole and force consistency",
        "private_source_native_motor_workflow",
        "multipole expansion before force/torque interpretation",
        ("multipole_expansion", "maxwell_stress_force", "field_quality"),
        ("source_field", "force_or_torque_trend", "reduced_fem_response"),
        "airgap_flux",
        "check quadrupole harmonic dominance, center-field null, and force frame",
        "multipole_field_quality_gate",
    ),
    DualLaneTrainingCase(
        "cable_impedance_frequency",
        "winding_loss",
        "Cable impedance frequency sweep",
        "private_source_native_motor_workflow",
        "frequency-axis and impedance unit metadata for eddy/conductor studies",
        ("skin_depth", "winding_ac_loss", "frequency_response"),
        ("pickup_flux", "reduced_fem_response", "coenergy"),
        "slip_loss",
        "check frequency grid, R/L units, monotonic trend, and conductor material package",
        "cable_impedance_frequency_gate",
    ),
    DualLaneTrainingCase(
        "hysteresis_minor_loop",
        "hysteresis",
        "Minor-loop hysteresis response in a magnetic circuit",
        "private_source_native_motor_workflow",
        "stateful material history must be tied to the B/H path",
        ("hysteresis_loss", "play_hysteresis", "nonlinear_material_iteration"),
        ("demag_field", "source_field", "coenergy"),
        "hysteresis_loss",
        "check loop orientation, remanence, coercivity, and energy-loss positivity",
        "hysteresis_minor_loop_energy_gate",
    ),
    DualLaneTrainingCase(
        "phase_rotor_alignment_symmetry",
        "pmsm",
        "Phase convention, rotor initial angle, and symmetry scaling",
        "private_source_native_motor_workflow",
        "table agreement starts with phase and symmetry metadata, not torque values",
        ("dq_control_layer", "torque_angle", "winding_factor"),
        ("pickup_flux", "force_or_torque_trend", "coenergy"),
        "torque",
        "check electrical/mechanical angle basis, phase offset, symmetry factor, and step count",
        "phase_symmetry_alignment_gate",
    ),
    DualLaneTrainingCase(
        "switched_reluctance_commutation",
        "srm",
        "SRM sequential phase excitation and inductance map",
        "private_source_native_motor_workflow",
        "commutation table and L(theta,i) map for reluctance torque",
        ("reluctance_torque", "saturating_inductance", "current_table"),
        ("pickup_flux", "force_or_torque_trend", "reduced_fem_response"),
        "ld_lq",
        "check phase order, current plateau, turn count, and torque periodicity",
        "srm_commutation_inductance_gate",
    ),
    DualLaneTrainingCase(
        "magnet_static_airgap_refinement",
        "spm",
        "Static PM motor mesh refinement around magnets and tooth tips",
        "external_readable_machine_fem_reference",
        "air-gap and magnet refinement controls waveform fidelity",
        ("pm_airgap_flux", "cogging_torque", "mesh_refinement"),
        ("demag_field", "source_field", "pickup_flux"),
        "airgap_flux",
        "check refinement ratio, magnet region labels, and Br waveform convergence",
        "pm_mesh_refinement_waveform_gate",
    ),
    DualLaneTrainingCase(
        "topology_optimization_pm_motor",
        "optimization",
        "PM motor topology optimization handoff to validation gates",
        "private_source_native_motor_workflow",
        "candidate geometry must return to AGE/VIM gates before being trusted",
        ("topology_optimization", "torque_angle", "cogging_torque"),
        ("source_field", "reduced_fem_response", "force_or_torque_trend"),
        "torque",
        "check design-variable package, objective observable, and validation replay id",
        "motor_topology_candidate_replay_gate",
    ),
    DualLaneTrainingCase(
        "magnetothermal_loss_chain",
        "multiphysics",
        "Magnetic-to-thermal loss-chain handoff",
        "private_source_native_motor_workflow",
        "loss buckets need temperature and material-state ownership",
        ("loss_bucket", "thermal_handoff", "hysteresis_loss"),
        ("pickup_flux", "reduced_fem_response", "coenergy"),
        "hysteresis_loss",
        "check loss-bucket units, heat-source region labels, and temperature state",
        "magthermal_loss_handoff_gate",
    ),
    DualLaneTrainingCase(
        "external_mesh_import_quality",
        "mesh",
        "External mesh import quality and region-label preservation",
        "private_source_native_motor_workflow",
        "mesh quality is not enough without region and boundary ownership",
        ("mesh_region_mapping", "geometry_label_contract", "mesh_quality"),
        ("interface_operator_contract", "reduced_fem_response", "source_field"),
        "airgap_flux",
        "check element type, region labels, boundary sets, and air-gap band quality",
        "external_mesh_region_boundary_gate",
    ),
    DualLaneTrainingCase(
        "torque_angle_harmonics",
        "pmsm",
        "Torque-angle harmonic budget and endpoint policy",
        "private_source_native_motor_workflow",
        "ripple harmonics require uniform angle rows and no repeated endpoint",
        ("torque_angle", "cogging_torque", "harmonic_budget"),
        ("force_or_torque_trend", "coenergy", "pickup_flux"),
        "torque",
        "check uniform angle step, endpoint removal, mean torque, and dominant harmonic",
        "torque_angle_harmonic_budget_gate",
    ),
    DualLaneTrainingCase(
        "efficiency_map_loss_bucket",
        "drive_loss",
        "Efficiency-map loss buckets and torque-speed envelope",
        "private_source_native_motor_workflow",
        "efficiency rows must retain mechanical/electrical power balance",
        ("efficiency_map", "loss_bucket", "field_weakening"),
        ("pickup_flux", "reduced_fem_response", "coenergy"),
        "field_weakening",
        "check P_in = P_out + losses, speed/torque grid, and dominant loss bucket",
        "efficiency_map_power_balance_gate",
    ),
    DualLaneTrainingCase(
        "short_circuit_demag_fault",
        "fault",
        "PM short-circuit fault and demagnetization stress",
        "private_source_native_motor_workflow",
        "fault current and demag margin must share the same dq convention",
        ("short_circuit", "demag_margin", "dq_control_layer"),
        ("demag_field", "pickup_flux", "source_field"),
        "demag_field",
        "check terminal voltage residuals, current decay, and d-axis demag fraction",
        "pm_short_circuit_demag_gate",
    ),
)


def motor_dual_lane_training_catalog(query: str = "all") -> list[dict[str, Any]]:
    """Return scrubbed dual-lane motor training cases."""

    q = query.strip().lower()
    cases = [case.to_dict() for case in CATALOG]
    if q in ("", "all", "*"):
        return cases
    return [
        case
        for case in cases
        if q in case["case_id"].lower()
        or q in case["family"].lower()
        or q in case["title"].lower()
        or q in case["learning_axis"].lower()
    ]


def motor_dual_lane_training_catalog_gate() -> dict[str, Any]:
    """Check that the public catalog is complete and provenance-scrubbed."""

    cases = motor_dual_lane_training_catalog()
    text = json.dumps(cases, ensure_ascii=True, sort_keys=True)
    forbidden_hits = [marker for marker in FORBIDDEN_PUBLIC_MARKERS if marker in text]
    missing_age = [case["case_id"] for case in cases if "radia_motor_age" not in case]
    missing_vim = [case["case_id"] for case in cases if "radia_motor_vim" not in case]
    source_classes = Counter(case["source_seed_class"] for case in cases)
    families = Counter(case["family"] for case in cases)
    checks = {
        "case_count_is_30": len(cases) == 30,
        "no_forbidden_public_markers": not forbidden_hits,
        "all_cases_have_age_lane": not missing_age,
        "all_cases_have_vim_lane": not missing_vim,
        "has_external_readable_reference_gap_closure": (
            source_classes["external_readable_machine_fem_reference"] >= 4
        ),
        "has_private_source_native_motor_workflows": (
            source_classes["private_source_native_motor_workflow"] >= 20
        ),
        "covers_core_machine_families": all(
            family in families
            for family in ("spm", "ipm", "induction", "srm", "pm_demag")
        ),
    }
    return {
        "schema_version": "radia-motor-dual-lane-training-catalog-gate/v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "count": len(cases),
        "checks": checks,
        "forbidden_hits": forbidden_hits,
        "missing_age": missing_age,
        "missing_vim": missing_vim,
        "source_seed_classes": dict(source_classes),
        "families": dict(families),
        "promotion_targets": [
            "radia-motor-age",
            "radia-motor-vim",
            "radia_mcp.motor.validation_lanes_knowledge",
        ],
    }


def route_dual_lane_training_case(goal: str) -> dict[str, Any]:
    """Pick a catalog case and expose both radia-motor lanes."""

    matches = motor_dual_lane_training_catalog(goal)
    if matches:
        case = matches[0]
    else:
        tokens = set(re.findall(r"[a-z0-9]+", goal.lower()))
        all_cases = motor_dual_lane_training_catalog()
        case = max(
            all_cases,
            key=lambda item: len(tokens & set(re.findall(r"[a-z0-9]+", json.dumps(item).lower()))),
        )
    return {
        "schema_version": "radia-motor-dual-lane-training-route/v1",
        "goal": goal,
        "selected_case": case,
        "next_public_calls": [
            f'motor_age_validation_plan("{case["family"]} {case["observable_family"]}")',
            'motor_validation_lane_template("hdiv_vim_reduced_fem")',
            'motor_validation_artifact_gate(..., "ngsolve_age")',
            'motor_validation_artifact_gate(..., "hdiv_vim_reduced_fem")',
        ],
    }


def format_motor_dual_lane_training_catalog(query: str = "all") -> str:
    """Format the scrubbed training catalog for MCP clients."""

    cases = motor_dual_lane_training_catalog(query)
    gate = motor_dual_lane_training_catalog_gate()
    lines = [
        "# radia-motor dual-lane training catalog",
        "",
        f"- status: `{gate['status']}`",
        f"- cases: `{len(cases)}` shown / `{gate['count']}` total",
        "- lanes: `radia-motor-age` and `radia-motor-vim`",
        "- public boundary: source-native provenance is private; this catalog keeps only scrubbed engineering lessons.",
        "",
    ]
    for index, case in enumerate(cases, 1):
        lines.append(f"{index}. `{case['case_id']}` -- {case['title']}")
        lines.append(f"   family: `{case['family']}`, observable: `{case['observable_family']}`")
        lines.append(f"   AGE targets: {', '.join(case['radia_motor_age']['targets'])}")
        lines.append(f"   VIM targets: {', '.join(case['radia_motor_vim']['targets'])}")
        lines.append(f"   gate: `{case['teaching_gate']}`")
    return "\n".join(lines).rstrip()
