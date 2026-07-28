"""Motor MCP Server (radia_mcp.motor)

Rotating electric machine analysis: 2D magnetodynamic A-formulation,
moving-band air-gap coupling, motor-type vocabulary (PMSM/IM/SRM/WFSM),
SynRM topology optimization (Wakao 2025 autoencoder + LS), and the
Darwin model time-domain solver (Kaimori-Mifune-Kameari-Wakao 2024).

Distilled from:
- public-safe curated corpus   (Sabariego-Gyselinck-Geuzaine)
- public-safe curated corpus   (Liu Xinyao thesis)

Usage:
    mcp-server-motor              # Start MCP server (stdio transport)
    mcp-server-motor --selftest   # Run self-test
"""

import asyncio
import json
import sys

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations
from ..common import register_status_tool
from ..common.mcp_contract import apply_tool_contract

from .onelab_knowledge import get_onelab_knowledge
from .topology_opt_knowledge import get_topology_opt_knowledge
from .darwin_model_knowledge import get_darwin_knowledge
from .femm_transient_knowledge import get_femm_transient_knowledge
from .henrotte_lineage_knowledge import get_henrotte_lineage_knowledge
from .hollaus_eddy_knowledge import get_hollaus_eddy_knowledge
from .hollaus_genealogy_knowledge import get_hollaus_genealogy
from .tritool_cross_reference_knowledge import get_tritool_cross_reference
from .deck_bridge_knowledge import get_deck_bridge
from .age_quality_knowledge import (
    format_age_validation_plan,
    get_age_quality_report,
    route_age_validation_plan,
)
from ..radia_ngsolve.circuit_excitation import compile_circuit_age_application
from ..radia_ngsolve.circuit_system import analyze_circuit_field
from .validation_lanes_knowledge import (
    format_artifact_gate_result,
    format_motor_validation_lanes,
    lane_template,
    validate_motor_validation_artifact,
)
from .dual_lane_training_catalog import (
    format_motor_dual_lane_training_catalog,
    motor_dual_lane_training_catalog_gate as build_dual_lane_training_catalog_gate,
    route_dual_lane_training_case,
)
from .triple_check_knowledge import (
    format_motor_triple_check_plan,
    format_triple_check_gate_result,
    route_motor_triple_check,
    validate_motor_triple_check_artifact,
)
from .simple_field_2d import (
    FieldQuickInput,
    evaluate_field_quick_check,
    format_field_quick_check,
    format_motor_validation_route,
    route_motor_validation,
)
from .planar_coupling_knowledge import get_planar_coupling
from .angle_periodic_rom_knowledge import get_angle_periodic_rom_knowledge
from .thermal_handoff import (
    motor_electrothermal_result_chain_gate as build_motor_electrothermal_result_chain_gate,
    motor_thermal_handoff_gate as build_motor_thermal_handoff_gate,
)
from .force_covariance import force_rotation_covariance_gate as build_force_rotation_covariance_gate
from .force_report_gate import force_report_method_metadata_gate as build_force_report_method_metadata_gate
from .phase_flux_park_gate import phase_flux_park_alignment_gate as build_phase_flux_park_alignment_gate
from .two_run_ldlq_gate import ipm_two_run_ldlq_gate as build_ipm_two_run_ldlq_gate
from .periodic_torque_sampling_gate import periodic_torque_sampling_gate as build_periodic_torque_sampling_gate
from .rotating_circuit_transient_gate import rotating_circuit_transient_gate as build_rotating_circuit_transient_gate
from .motion_table_gate import motion_table_coordinate_gate as build_motion_table_coordinate_gate
from .magnet_model_handoff_gate import magnet_model_handoff_gate as build_magnet_model_handoff_gate
from .magnetization_group_symmetry_gate import (
    mirror_symmetric_three_magnet_handoff_gate as build_mirror_magnet_handoff_gate,
)
from .variable_magnet_gate import variable_magnet_material_parameter_gate as build_variable_magnet_material_parameter_gate
from .permanent_magnet_force_pair_gate import permanent_magnet_force_pair_gate as build_permanent_magnet_force_pair_gate
from .demagnetization_history_gate import permanent_magnet_demagnetization_history_gate as build_permanent_magnet_demagnetization_history_gate
from .dual_torque_curve_gate import dual_torque_method_curve_gate as build_dual_torque_method_curve_gate
from .virtual_work_width_gate import motor_virtual_work_width_ladder_gate as build_motor_virtual_work_width_ladder_gate
from .transient_no_load_load_gate import (
    motor_transient_no_load_load_cycle_gate as build_transient_no_load_load_cycle_gate,
)

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return "Motor bibliography index not yet generated."


mcp = FastMCP("mcp-server-motor")


def _decode_owned_worker_json(stdout: bytes) -> dict:
    """Decode the final JSON object while tolerating native stdout diagnostics."""

    text = stdout.decode("utf-8", errors="replace").strip()
    if not text:
        raise ValueError("owned worker returned no JSON output")
    for line in reversed(text.splitlines()):
        candidate = line.strip()
        if not candidate.startswith("{"):
            continue
        try:
            result = json.loads(candidate)
        except json.JSONDecodeError:
            continue
        if isinstance(result, dict):
            return result
    raise ValueError("owned worker output did not end with a JSON object")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def motor_transient_no_load_load_cycle_gate(summary_json: str) -> str:
    """Gate paired no-load and loaded three-phase transient cycles."""

    try:
        result = build_transient_no_load_load_cycle_gate(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "motor_transient_no_load_load_cycle_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_virtual_work_width_ladder_gate(summary_json: str) -> str:
    """Select a coenergy-difference angle width against independent torque."""
    try:
        result = build_motor_virtual_work_width_ladder_gate(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "motor_virtual_work_width_ladder_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_dual_torque_method_curve_gate(summary_json: str) -> str:
    """Gate two independently evaluated static-torque curves."""
    try:
        result = build_dual_torque_method_curve_gate(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "dual_torque_method_curve_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_permanent_magnet_demagnetization_history_gate(
    summary_json: str,
    state_tolerance: float = 1.0e-9,
    minimum_damage_fraction: float = 1.0e-3,
) -> str:
    """Gate irreversible PM state across one history or a replayed case family."""

    try:
        result = build_permanent_magnet_demagnetization_history_gate(
            summary_json,
            state_tolerance,
            minimum_damage_fraction,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "permanent_magnet_demagnetization_history_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_permanent_magnet_force_pair_gate(
    summary_json: str,
    magnitude_relative_tolerance: float = 2.0e-2,
    off_axis_relative_tolerance: float = 1.0e-3,
) -> str:
    """Gate attraction/repulsion reversal for a facing permanent-magnet pair."""

    try:
        result = build_permanent_magnet_force_pair_gate(
            summary_json,
            magnitude_relative_tolerance,
            off_axis_relative_tolerance,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "permanent_magnet_force_pair_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_variable_magnet_material_gate(
    parameters: dict,
    parameter_authority: str,
    study_label_is_parameter_authority: bool = False,
) -> str:
    """Gate variable-PM material parameters and their authoritative source."""

    try:
        result = build_variable_magnet_material_parameter_gate(
            parameters,
            parameter_authority=parameter_authority,
            study_label_is_parameter_authority=study_label_is_parameter_authority,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "variable_magnet_material_parameter_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_magnet_model_handoff_gate(
    residual_phases: list[list[float]],
    nonlinear_tolerance: float,
    source_result_artifact_id: str,
    source_result_digest: str,
    magnet_control_artifact_id: str,
    magnet_control_digest: str,
    magnet_geometry_artifact_id: str,
    magnet_geometry_digest: str,
    numbering_policy: str,
    element_id_offset: int,
    node_id_offset: int,
    material_mapping_count: int,
    geometry_transform: str,
) -> str:
    """Gate a converged source result and two-file downstream magnet model."""
    try:
        result = build_magnet_model_handoff_gate(
            residual_phases,
            nonlinear_tolerance=nonlinear_tolerance,
            source_result_artifact_id=source_result_artifact_id,
            source_result_digest=source_result_digest,
            magnet_control_artifact_id=magnet_control_artifact_id,
            magnet_control_digest=magnet_control_digest,
            magnet_geometry_artifact_id=magnet_geometry_artifact_id,
            magnet_geometry_digest=magnet_geometry_digest,
            numbering_policy=numbering_policy,
            element_id_offset=element_id_offset,
            node_id_offset=node_id_offset,
            material_mapping_count=material_mapping_count,
            geometry_transform=geometry_transform,
        )
    except (TypeError, ValueError) as exc:
        result = {"policy": "magnet_model_handoff_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def motor_mirror_symmetric_three_magnet_handoff_gate(summary_json: str) -> str:
    """Gate grouped magnetization vectors, mirror symmetry, and fresh replay."""

    try:
        result = build_mirror_magnet_handoff_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "mirror_symmetric_three_magnet_handoff_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_motion_table_coordinate_gate(
    translation_times_s: list[float],
    translation_vectors: list[list[float]],
    rotation_times_s: list[float],
    rotation_vectors: list[list[float]],
    coordinate_frame_id: str,
    translation_unit: str = "mm",
    rotation_unit: str = "deg",
    motion_semantics: str = "cumulative_displacement",
) -> str:
    """Validate independent 3D translation and rotation motion tables."""
    try:
        result = build_motion_table_coordinate_gate(
            translation_times_s,
            translation_vectors,
            rotation_times_s,
            rotation_vectors,
            coordinate_frame_id=coordinate_frame_id,
            translation_unit=translation_unit,
            rotation_unit=rotation_unit,
            motion_semantics=motion_semantics,
        )
    except (TypeError, ValueError) as exc:
        result = {"policy": "motion_table_coordinate_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_periodic_torque_sampling_gate(
    period_deg: float,
    sample_count: int,
    endpoint_included: bool,
    spectrum_excludes_duplicate_endpoint: bool,
    torque_min_Nm: float,
    torque_max_Nm: float,
    speed_rps: float,
    expected_step_deg: float | None = None,
    step_tolerance_deg: float = 1.0e-9,
) -> str:
    """Validate periodic torque sampling and FFT endpoint ownership."""
    try:
        result = build_periodic_torque_sampling_gate(
            period_deg=period_deg,
            sample_count=sample_count,
            endpoint_included=endpoint_included,
            spectrum_excludes_duplicate_endpoint=spectrum_excludes_duplicate_endpoint,
            torque_min_Nm=torque_min_Nm,
            torque_max_Nm=torque_max_Nm,
            speed_rps=speed_rps,
            expected_step_deg=expected_step_deg,
            step_tolerance_deg=step_tolerance_deg,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "periodic_torque_sampling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def motor_rotating_circuit_transient_gate(summary_json: str) -> str:
    """Gate rotating-circuit identities and endpoint state before FFT use."""
    try:
        result = build_rotating_circuit_transient_gate(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "rotating_circuit_transient_endpoint_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)

@mcp.tool()
def motor_onelab(topic: str = "overview") -> str:
    """
    ONELAB/GetDP electric-machine reference template knowledge.

    Args:
        topic: One of:
            "overview"        - Bundle layout + reference machine table
            "groups"          - Region / FunctionSpace / Constraint vocab
            "formulation"     - 2D A-formulation weak form + non-linear ν
            "magstadyn_source"- VERIFIED weak form/torque/slip/Park from .pro source
            "twod_corrections"- 2D→3D corrections: end-effects/skew/lamination/etc
            "motor_types"     - PMSM, IM, SRM, WFSM, shaded-pole table
            "ngsolve_xlate"   - GetDP → NGSolve translation patterns
            "analysis_modes"  - Static / Time-Domain / Frequency-Domain
            "circuit"         - External-circuit coupling (MNA)
            "post"            - Torque (Arkkio), iron-loss, FFT recipes
            "liu_thesis"      - 2025 SynRM thesis: ONELAB vs an independent ref
            "all"             - Everything
    """
    return get_onelab_knowledge(topic)


@mcp.tool()
def motor_topology_optimization(topic: str = "wakao_ae_ls") -> str:
    """
    SynRM topology optimization (Wakao 2025 AE + level-set), bridge and channel shaping.

    Args:
        topic: One of:
            "wakao_ae_ls"                - Autoencoder + LS hybrid method
            "liu_thesis_application"     - Liu Xinyao end-to-end recipe
            "pareto_navigation"          - 2D latent-space Pareto walk
            "saturable_bridge_hodograph" - IPM bridge/rib at the saturation cap:
                                           free boundary -> hodograph coordinate
                                           line, one linear solve (verified)
            "flux_channel_hodograph"     - SynRM channel: exact saturated sizing
                                           chart (pure turn = annulus, any
                                           material) + free-form collecting
                                           channel design (verified; naive
                                           +78% iron)
            "all"                        - Everything
    """
    return get_topology_opt_knowledge(topic)


@mcp.tool()
def motor_darwin_model(topic: str = "overview") -> str:
    """
    Darwin-model time-domain formulation (capacitive + inductive coupling).

    Args:
        topic: One of:
            "overview"     - When to use Darwin vs MQS / EQS / full Maxwell
            "formulation"  - Coulomb-gauge A-φ-χ symmetric system
            "vs_peec"      - Darwin (FE) vs PEEC (IE) practical comparison
            "all"          - Everything
    """
    return get_darwin_knowledge(topic)


@mcp.tool()
def motor_planar_coupling(topic: str = "overview") -> str:
    """
    2D PLANAR machine modelling in radia: HDiv-VIM soft-iron demag + the shared
    postprocessing (radia.planar_charges) + the staggered eddy-current coupling
    (radia.planar_eddy) for PM motors / induction machines / eddy-current brakes.

    Analytic-led + fully gated (Bessel / 2D dipole / monolithic FEM). Companion to
    the HDiv-VIM 2D and planar material validation tests.

    Args:
        topic: One of:
            "overview"      - the HDiv demag method + one shared coupling layer
            "eddy_coupling" - staggered HDiv-VIM <-> reduced-Az eddy FEM (maglev/IM/ECB)
            "pm_motor"      - permanent magnets: design A (magnets=) / B (pm=) / unified rotor
            "nonlinear"     - nonlinear soft iron + eddy (effective-chi AC)
            "anisotropic"   - GO-steel tensor chi (planar_aniso, direct-N) + embedded PM design-B
            "hysteresis"    - 2D play-hysteresis demag (planar_hysteresis, direct-N + Newton)
            "api"           - API quick reference
            "validation"    - the gated validation ladder
            "all"           - everything
    """
    return get_planar_coupling(topic)


@mcp.tool()
def motor_angle_periodic_rom(topic: str = "architecture") -> str:
    """HCurl Eddy Bubble + HDiv-MMM angle-periodic motor ROM knowledge.

    Args:
        topic: architecture, face_policy, angle_rom, time_domain, ports,
            mesh_gate, validation, limits, or all.
    """
    return get_angle_periodic_rom_knowledge(topic)


@mcp.tool()
def motor_femm_transient(topic: str = "lab_recommendation") -> str:
    """
    FEMM newbuild transient solver — Lange-Henrotte-Hameyer 2009
    incremental-permeability linearization with sliding-band air-gap BC.

    Source: David Meeker, https://www.femm.info/doku/doku.php?id=newbuild
    + E. Lange, F. Henrotte, K. Hameyer, IEEE Trans. Mag. 45(3):1258-1261,
    2009 (DOI: 10.1109/TMAG.2009.2012585).

    Distinguishing feature: factor the nonlinear FE solve from time
    integration.  At a chosen rotor angle, 1 nonlinear FE solve +
    N_phase linear incremental-permeability solves (shared factorization)
    → extract `(L_inc(i, theta), e_bemf, T_em)` → circuit simulator
    integrates `L_inc di/dt + R i = v - e_bemf` at a much smaller PWM
    time step without re-running the FE.  Rotor angle advances
    analytically between FE refreshes.  Much faster than fully-coupled
    nonlinear FE+circuit (ONELAB) for controller-in-the-loop simulations.

    Args:
        topic: One of:
            "lab_recommendation"     - ★ Lab-canonical path + Hameyer/Henrotte refs
            "lange_2009"             - Underlying 2009 paper formulation
            "henrotte_no_rotation"   - Why the rotor never rotates in the FE
            "femm_newbuild"          - Meeker's FEMM implementation + antunes.fem
            "sliding_vs_moving_band" - Sliding-band BC vs moving-band remesh
            "ngsolve_recipe"         - Port to NGSolve with shared factorization
            "simulink_coupling"      - Simulink S-function + HIL/RT notes
            "all"                    - Everything
    """
    return get_femm_transient_knowledge(topic)


@mcp.tool()
def motor_henrotte_lineage(topic: str = "research_arc") -> str:
    """
    The Henrotte–Hameyer–RWTH research arc (energy-consistent E&M FE).

    Consolidates 7 key papers from the 1993-2018 ULg/UCL/RWTH school
    that culminated in the Lange-Henrotte-Hameyer 2009 field-circuit
    coupling.  All available in the Sugahara lab library.

    Args:
        topic: One of:
            "axisym_1993"       - Henrotte basis s=r²/2 (foundation of axifem)
            "source_field_1997" - Dular-Henrotte cut-aware source field
            "energy_hys_2006"   - Energy-based vector hysteresis (friction)
            "variational_2013"  - Variational FE-friendly extension
            "jacques_2018"      - 254p PhD monograph (canonical reference)
            "carstensen_2007"   - 187p SRM winding eddy currents PhD
            "research_arc"      - Synthesis: how the papers connect + roadmap
            "all"               - Everything
    """
    return get_henrotte_lineage_knowledge(topic)


@mcp.tool()
def motor_hollaus_eddy(topic: str = "intro") -> str:
    """
    Karl Hollaus / TU Wien MSFEM for laminated-iron eddy currents.

    THE answer to "Lange-Henrotte-Hameyer doesn't handle eddy currents".
    MSFEM homogenizes the through-the-thickness lamination structure
    so a 3D-laminated machine reduces to a 2D-1D or single-cell
    problem with ~100x fewer DOFs at <1% loss error.

    Args:
        topic: One of:
            "intro"             - Why standard FE fails; MSFEM idea
            "msfem_2d1d"        - 2D-1D MSFEM T-formulation (rotating machines)
            "t_phi_open"        - T-Phi-Phi MSFEM for OPEN circuits (motors)
            "effective_material"- EM homogenization (practical workhorse)
            "error_estimator"   - Equilibrated error estimator (Prager-Synge)
            "msfem_3d"          - Full 3D MSFEM for end-region asymmetric cases
            "nonlinear_esi"     - Nonlinear ESI in magnetic scalar potential
                                  (cross-references radia_mcp.ih.esim_cell_problem)
            "msfem_plus_mor"    - Cascade with POD/snapshots for transient PWM
            "nonasymptotic_homogenization" - Schobinger-Hollaus-Tsukerman 2020:
                                             two mu_eff definitions, resonance
            "msfem_mor_deim_detail" - Hollaus-Schobinger 2024: MSFEM+MOR+DEIM
                                       with structure-preserving DEIM
            "hanser_2025_circuit_coupling" - Hanser 2025 COMPEL: A-formulation
                                              + Schur for circuit coupling
                                              (300x-5000x speedup, IMPLEMENTATION
                                              ROADMAP for calc_motor_lamination v3)
            "hollaus_2014_nonlinear_two_scale" - Foundational two-scale paper
                                                   (Hollaus-Hannukainen-Schoberl 2014)
            "hierarchical_error_estimator" - Schobinger-Hollaus 2021 hierarchical
                                              error estimator for 3D MSFEM
            "frljic_2026_perpendicular_flux" - Frljic et al. 2026: anisotropic
                                                 mu_eff for open-type cores with
                                                 perpendicular flux
            "motor_strategy"    - SYNTHESIS: hybrid stack for eddy-current-capable
                                  motor analysis combining Lange-HH + MSFEM + ESIM
            "all"               - Everything
    """
    return get_hollaus_eddy_knowledge(topic)


@mcp.tool()
def motor_em_force_recipe(topic: str = "method_choice") -> str:
    """
    Practical NGSolve EM-force recipe for motor analysis.

    Forwards to `differential_forms_em_force_recipe` (radia_mcp.
    differential_forms.em_force_ngsolve_recipe_knowledge).  Tells you
    which method to use for motor torque / Maxwell stress on stator
    teeth / Lorentz on coil ends / etc.

    Achieved 0.5% Newton-3 violation on iron+coil 2D benchmark via
    Mertens-Hameyer 2007 v_order matching + far-field outer Dirichlet
    box (10x better than Pile-2018 "production" 5%).

    Args:
        topic: One of:
            "method_choice"         - Decision tree
            "high_order_recipe"     - Mertens-Hameyer 2007 prescription
            "common_pitfalls"       - 6 traps that destroy accuracy
            "validation_protocol"   - Newton-3 + path-independence
            "implementation_status" - calc_em_force.py status
            "full_example"          - Working code template
            "all"                   - Everything
    """
    # Lazy import to keep cross-package dependency soft
    from radia_mcp.differential_forms.em_force_ngsolve_recipe_knowledge \
        import get_em_force_ngsolve_recipe
    return get_em_force_ngsolve_recipe(topic)


@mcp.tool()
def motor_em_force_extras(topic: str = "all") -> str:
    """
    Forward to `differential_forms_em_force_extras` -- advanced EM force
    topics beyond the 7-method catalog.

    Args:
        topic: One of:
            "permanent_magnet_force"  - PM force formulations (bound
                                          current K = M x n, Kelvin,
                                          scalar charge)
            "energy_method_classical" - constant flux vs current,
                                          Legendre W <-> W'
            "sensitivity_vwp"         - shape derivative (Pironneau-
                                          Allaire), proper VWP fix
            "lorentz_canonical"       - Lorentz as first-class method
            "meissner_force"          - Type-I superconductor mu_r -> 0
            "all"                     - Everything
    """
    from radia_mcp.differential_forms.em_force_extras_knowledge \
        import get_em_force_extras
    return get_em_force_extras(topic)


@mcp.tool()
def motor_hollaus_genealogy(view: str = "by_topic") -> str:
    """
    Visualize the Karl Hollaus / TU Wien MSFEM research genealogy
    (1992-2026 paper lineage with citation DAG).

    The MSFEM stack now spans 13+ papers (see motor_hollaus_eddy for the
    physics content). This tool maps the dependency graph: which paper
    builds on which, grouped by research thread.

    Args:
        view: One of:
            "by_topic"  - Grouped by research thread (default; best for
                          orienting a new reader)
            "chrono"    - Chronological flat list with DOIs and parents
            "ascii"     - ASCII art tree (roots -> descendants)
            "dot"       - Graphviz DOT format (paste into
                          https://dreampuf.github.io/GraphvizOnline/)
            "all"       - All four views
    """
    return get_hollaus_genealogy(view)


@mcp.tool()
def motor_bibliography(query: str = "") -> str:
    """Search the motor analysis bibliography catalog."""
    return get_bibliography_index(query)


@mcp.tool()
def motor_tritool_cross_reference(topic: str = "overview") -> str:
    """
    Tri-tool cross-reference: FEMM / JMAG / radia-ngsolve (相互学習).

    Cross-learning that ties the lab's three motor-FEA tools together so each
    strengthens the others. Mirrors the shared cross_ref.json carried by the
    two dedicated lab-internal knowledge servers (mcp-server-femm and
    mcp-server-jmag), whose commercial-tool specifics stay lab-private.

    Args:
        topic: One of:
            "overview"          - The ecosystem + the two new MCP servers
            "capability_matrix" - Per-capability FEMM/JMAG/radia table + strongest
            "radia_can_exceed"  - Where radia-ngsolve matches or beats the others
            "jmag_only"         - Genuine JMAG-only capabilities (gaps)
            "femm_role"         - FEMM as the open-source 2D yardstick + .fem strategy
            "roadmap"           - radia-motor strengthening roadmap (ranked)
            "all"               - Everything
    """
    return get_tritool_cross_reference(topic)


@mcp.tool()
def motor_deck_bridge(topic: str = "overview") -> str:
    """
    Public-safe motor deck corpus bridge for radia-motor.

    Explains what the external public motor-deck corpus already covers, where
    its validation depth is still thin, and which radia-motor / radia-ngsolve
    anchors should be strengthened next.

    Args:
        topic: One of:
            "overview"                  - public motor deck corpus status
            "coverage_matrix"           - deck archetypes -> radia targets
            "insufficiency_audit"       - What is broad vs still only proxy-validated
            "routing_playbook"          - How to use public decks and radia-motor together
            "radia_strengthening_queue" - Next radia-side upgrades
            "jmag_coverage_reality"     - Turnkey production motor coverage boundary
            "age_vs_field_strategy"       - Why AGE is the main 2D path; where field quick checks help
            "all"                       - Everything
    """
    return get_deck_bridge(topic)


@mcp.tool()
def motor_age_quality(topic: str = "overview") -> str:
    """
    NGSolve AGE quality gates for radia-motor.

    This is the public-safe readiness layer for treating AGE as the main
    radia-motor validation path. It lists physical quantities, test evidence,
    acceptance conditions, family coverage, publication labels, and limitations.

    Args:
        topic: One of:
            "overview"           - AGE readiness policy
            "gate_matrix"        - Gate -> quantity/test/acceptance table
            "family_matrix"      - Motor family -> required AGE gates
            "publication_policy" - What can be called AGE-verified
            "runbook"            - Targeted pytest commands
            "limitations"        - Required caveats for public claims
            "all"                - Everything
    """
    return get_age_quality_report(topic)


@mcp.tool()
def motor_age_validation_plan(goal: str) -> str:
    """
    Route a motor prompt to the required NGSolve AGE quality gates.

    Args:
        goal: Natural-language motor analysis request, e.g. "IPM hairpin MTPA"
            or "induction cage slip loss".
    """
    return format_age_validation_plan(route_age_validation_plan(goal))


@mcp.tool()
def motor_circuit_age_application_plan(application_json: str) -> str:
    """Compile planar/axisymmetric current circuits and optional AGE motion.

    Series regions receive explicit signed ampere-turn source densities.
    Parallel branches remain field-circuit unknowns with one common voltage
    and a total-current constraint. Rotary and linear AGE motion is represented
    by Fourier phase factors, so neither FE mesh is rebuilt.
    """

    try:
        payload = json.loads(application_json)
        result = compile_circuit_age_application(payload)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.circuit-age-application.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def motor_circuit_field_analysis(analysis_json: str) -> str:
    """Analyze a solved 2D field-circuit, AGE sweep, or MEX state-space request.

    The pure numerical kernel supports planar and axisymmetric source matrices.
    Parallel circuits are augmented with branch-current and common-voltage
    unknowns; AGE sweeps reuse one operator; linear RL reductions target the
    existing native Simulink state-space MEX/S-Function.
    """

    try:
        result = analyze_circuit_field(json.loads(analysis_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "schema": "radia.circuit-field-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
async def motor_vol2d_circuit_analysis(analysis_json: str) -> str:
    """Analyze a Netgen dimension-2 ``.vol`` and its coupled field circuit.

    NGSolve assembly runs in one owned worker process because its runtime is not
    thread-safe inside FastMCP's request executor.  A timeout terminates only
    that worker.  The mesh stays a replayable input artifact rather than a
    tracked fixture.
    Planar P/Q elements use NGSolve H1; axisymmetric P/Q elements use Radia's
    Henrotte space.  Signed turns are assembled on the named mesh materials
    before the v2 series/parallel circuit kernel is evaluated.
    """

    process = None
    try:
        json.loads(analysis_json)
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.vol2d_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(analysis_json.encode("utf-8")),
            timeout=120.0,
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "vol2d worker failed")
        result = _decode_owned_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.vol2d-circuit-analysis.v1",
            "status": "timeout",
            "error": "vol2d worker exceeded 120 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.vol2d-circuit-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
async def motor_vol2d_dynamic_analysis(analysis_json: str) -> str:
    """Analyze nonlinear, conductive-transient, or MEX-ready ``.vol`` data.

    The request must carry a complete SI material law for every mesh material.
    Supported operations are ``assemble``, ``nonlinear_static``, ``transient``,
    and ``state_space``.  NGSolve/axifem execution is isolated in one owned
    worker so timeout or cancellation never terminates a shared solver session.
    """

    process = None
    try:
        json.loads(analysis_json)
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.vol2d_dynamics_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(analysis_json.encode("utf-8")),
            timeout=120.0,
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "vol2d dynamics worker failed")
        result = _decode_owned_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.vol2d-dynamic-analysis.v1",
            "status": "timeout",
            "error": "vol2d dynamics worker exceeded 120 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.vol2d-dynamic-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
async def motor_vol2d_force_analysis(analysis_json: str) -> str:
    """Extract boundary-aware force from a solved dimension-2 ``.vol`` model.

    ``solve`` reconstructs the NGSolve/axifem field and evaluates an air-only
    weighted-stress band.  Planar conductor targets can request an independent
    Lorentz comparison; passive magnetic targets remain weighted-stress plus
    virtual-work problems.  ``virtual_work_gate`` and ``refinement_gate`` close
    the displacement and mesh-evidence contracts.  All NGSolve work runs in one
    owned worker that can be cancelled without touching shared solver sessions.
    """

    process = None
    try:
        json.loads(analysis_json)
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.vol2d_force_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(analysis_json.encode("utf-8")),
            timeout=120.0,
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "vol2d force worker failed")
        result = _decode_owned_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.vol2d-force-analysis.v1",
            "status": "timeout",
            "error": "vol2d force worker exceeded 120 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.vol2d-force-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool(
    annotations=ToolAnnotations(
        readOnlyHint=True,
        destructiveHint=False,
        idempotentHint=True,
        openWorldHint=False,
    )
)
async def motor_age_periodic_motion_analysis(analysis_json: str) -> str:
    """Run or validate a no-remesh AGE periodic-motion torque sweep.

    ``solve`` consumes a generated dimension-2 ``.vol`` with disconnected
    rotor/stator FE regions and named gap rings.  One AGE factorization is reused
    while rotor angle changes only the harmonic phase.  ``periodic_sector_gate``
    checks periodic/anti-periodic sign, sector scaling, and fixed execution
    identities.  The owned worker can be cancelled without touching shared CAE
    sessions.
    """

    process = None
    try:
        json.loads(analysis_json)
        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-m",
            "radia_mcp.radia_ngsolve.age_periodic_motion_worker",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(
            process.communicate(analysis_json.encode("utf-8")),
            timeout=180.0,
        )
        if process.returncode != 0:
            message = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(message[-1000:] or "AGE periodic-motion worker failed")
        result = _decode_owned_worker_json(stdout)
    except asyncio.TimeoutError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        result = {
            "schema": "radia.age-periodic-motion-analysis.v1",
            "status": "timeout",
            "error": "AGE periodic-motion worker exceeded 180 seconds",
        }
    except asyncio.CancelledError:
        if process is not None and process.returncode is None:
            process.kill()
            await process.wait()
        raise
    except (json.JSONDecodeError, OSError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": "radia.age-periodic-motion-analysis.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def motor_validation_lanes(topic: str = "overview") -> str:
    """
    Cross-validation lane policy for radia-motor.

    Use this before promoting a private product, lab-local, open-source,
    stored-regression, or analytic comparison into radia-motor knowledge.
    It keeps NGSolve+AGE and the coupled HDiv-MMM + HCurl eddy-bubble lane
    independent, so
    cross-validation artifacts train the correct solver path.

    Args:
        topic: One of:
            "overview"          - why the radia lanes are independent
            "lane_matrix"       - lane -> observable/metric/promotion table
            "source_policy"     - public-safe handling of private references
            "promotion_policy"  - artifact-to-MCP learning rules
            "runbook"           - targeted validation commands
            "all"               - Everything
    """
    return format_motor_validation_lanes(topic)


@mcp.tool()
def motor_validation_lane_template(lane_id: str = "all") -> str:
    """
    Return the JSON artifact template for a motor validation lane.

    Args:
        lane_id: "ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble", or "all".
    """
    return json.dumps(lane_template(lane_id), indent=2, sort_keys=True)


@mcp.tool()
def motor_validation_artifact_gate(artifact_json: str, expected_lane: str = "") -> str:
    """
    Check whether a motor cross-validation artifact can train radia-motor.

    Args:
        artifact_json: JSON object text containing the cross-validation summary.
        expected_lane: Optional lane id to enforce:
            "ngsolve_age" or "hdiv_mmm_hcurl_eddy_bubble".
    """
    result = validate_motor_validation_artifact(artifact_json, expected_lane)
    return format_artifact_gate_result(result)


@mcp.tool()
def motor_dual_lane_training_catalog(topic: str = "all") -> str:
    """
    Return the public-safe wide motor learning catalog.

    Every case routes to both `radia-motor-age` and `radia-motor-mmm-eddy`.
    Source-native provenance is deliberately scrubbed from this public surface.

    Args:
        topic: "all" or a family/case/search phrase.
    """
    return format_motor_dual_lane_training_catalog(topic)


@mcp.tool()
def motor_dual_lane_training_gate() -> str:
    """Check that the public motor catalog is complete and provenance-scrubbed."""
    return json.dumps(build_dual_lane_training_catalog_gate(), indent=2, sort_keys=True)


@mcp.tool()
def motor_dual_lane_training_route(goal: str) -> str:
    """Route a motor prompt to one catalog case and both radia-motor lanes."""
    return json.dumps(route_dual_lane_training_case(goal), indent=2, sort_keys=True)


@mcp.tool()
def motor_triple_check_plan(goal: str) -> str:
    """
    Plan the standard radia-motor comparison.

    The plan uses the public ELF/MAGIC MCP surface for motor examples, the
    `ngsolve_age` lane and the coupled `hdiv_mmm_hcurl_eddy_bubble` lane as the
    mandatory comparison pair. HDiv-MMM and HCurl eddy-bubble are verified as
    one mixed system with a solver-ready artifact and matching model identity.

    Args:
        goal: Natural-language motor goal, e.g.
            "IPM hairpin motor flux linkage and MTPA".
    """
    return format_motor_triple_check_plan(route_motor_triple_check(goal))


@mcp.tool()
def motor_triple_check_artifact_gate(artifact_json: str) -> str:
    """
    Validate a combined AGE and HDiv-MMM/HCurl eddy-bubble motor artifact.

    `ngsolve_age` and `hdiv_mmm_hcurl_eddy_bubble` are mandatory for radia-motor
    learning claims.

    Args:
        artifact_json: JSON object text with schema
            `radia-motor-triple-check-artifact/v1`.
    """
    result = validate_motor_triple_check_artifact(artifact_json)
    return format_triple_check_gate_result(result)


@mcp.tool()
def motor_field_quick_check(
    motor_type: str = "spm",
    pole_pairs: int = 4,
    airgap_radius_m: float = 0.05,
    stack_length_m: float = 0.05,
    airgap_m: float = 1.0e-3,
    turns_per_phase: float = 50.0,
    phase_current_a: float = 10.0,
    electrical_angle_deg: float = 0.0,
    magnet_br_t: float = 1.2,
    magnet_thickness_m: float = 3.0e-3,
    magnet_arc_fraction: float = 0.75,
    saliency_ratio_lq_over_ld: float = 1.5,
    slip_hz: float = 5.0,
) -> str:
    """
    First-order 2D magnetic-circuit/BEM-like motor quick check.

    This is a public-safe, approximate magnetic-circuit evaluator for prompt-time
    sanity checks. It estimates PM flux linkage, back-EMF constant, dq torque
    proxy, and an induction slip-loss proxy, then routes the result to the
    NGSolve AGE validation targets that should be used for real verification.

    Args:
        motor_type: "spm", "ipm", "induction", "srm", "synrm",
            "hysteresis", etc.
        pole_pairs: Number of pole pairs.
        airgap_radius_m: Air-gap radius in meters.
        stack_length_m: Active stack length in meters.
        airgap_m: Mechanical air gap in meters.
        turns_per_phase: Effective series turns per phase.
        phase_current_a: Peak phase current.
        electrical_angle_deg: Electrical current angle from q-axis convention.
        magnet_br_t: PM remanence in tesla.
        magnet_thickness_m: Magnet thickness in meters.
        magnet_arc_fraction: Fraction of pole pitch covered by magnet.
        saliency_ratio_lq_over_ld: Lq/Ld proxy for IPM/SynRM/SRM checks.
        slip_hz: Slip frequency for induction-machine proxy checks.
    """
    inp = FieldQuickInput(
        motor_type=motor_type,
        pole_pairs=pole_pairs,
        airgap_radius_m=airgap_radius_m,
        stack_length_m=stack_length_m,
        airgap_m=airgap_m,
        turns_per_phase=turns_per_phase,
        phase_current_a=phase_current_a,
        electrical_angle_deg=electrical_angle_deg,
        magnet_br_t=magnet_br_t,
        magnet_thickness_m=magnet_thickness_m,
        magnet_arc_fraction=magnet_arc_fraction,
        saliency_ratio_lq_over_ld=saliency_ratio_lq_over_ld,
        slip_hz=slip_hz,
    )
    return format_field_quick_check(evaluate_field_quick_check(inp))


@mcp.tool()
def motor_validation_router(goal: str) -> str:
    """
    Route a motor prompt to a public deck, field quick check, and NGSolve AGE validation.

    This is the dispatch layer for the hybrid workflow:
    public decks for input authoring, the lightweight 2D field quick check for
    first-order sign/scale checks, and NGSolve AGE / radia-ngsolve for
    independent motor-physics validation.
    """
    return format_motor_validation_route(route_motor_validation(goal))


@mcp.tool()
def motor_thermal_handoff_gate(
    loss_buckets_json: str,
    network_json: str,
    mesh_regions_json: str,
    relative_tolerance: float = 1.0e-9,
) -> str:
    """Validate one motor-loss table for both LPTN and 3D all-hex thermal paths.

    The gate does not solve the thermal problem. It verifies that the same
    non-negative regional losses are assigned exactly once to a connected
    lumped thermal network and exactly once to positive hexahedral mesh
    regions, with regional and total heat conservation.

    Args:
        loss_buckets_json: JSON object mapping region names to loss in watts.
        network_json: JSON object with ``ambient_node``, ``nodes`` and
            positive-resistance ``branches``. Non-ambient nodes require
            positive ``capacitance_J_per_K`` and may own ``source_regions``.
        mesh_regions_json: JSON list of ``region``, ``cell_type``,
            ``cell_count`` and ``loss_W`` records. Cell type must be hex.
        relative_tolerance: Positive relative tolerance for loss matching.
    """
    return build_motor_thermal_handoff_gate(
        loss_buckets_json,
        network_json,
        mesh_regions_json,
        relative_tolerance,
    )


@mcp.tool()
def motor_electrothermal_result_chain_gate(
    chain_json: str,
    absolute_tolerance_W: float = 1.1e-2,
    relative_tolerance: float = 1.0e-3,
) -> str:
    """Gate a four-stage motor electrothermal result handoff.

    The stage artifacts must be fresh, uniquely identified, and pin exact
    upstream result digests.  The six three-phase motor loss channels are
    owned once, then scaled by an explicit geometry-symmetry fraction before
    comparison with the steady thermal input and temperature rise.
    """
    return build_motor_electrothermal_result_chain_gate(
        chain_json,
        absolute_tolerance_W=absolute_tolerance_W,
        relative_tolerance=relative_tolerance,
    )


@mcp.tool()
def motor_force_rotation_covariance_gate(
    reference_force_json: str,
    rotated_force_json: str,
    rotation_deg: float,
    relative_tolerance: float = 1.0e-3,
) -> str:
    """Check that a planar force vector follows a rotated excitation/geometry.

    This is solver-independent and is useful for symmetric motors, magnetic
    bearings, and actuators. Force objects contain global ``Fx`` and ``Fy``.
    ``rotation_deg`` follows the standard counter-clockwise convention.
    """
    return build_force_rotation_covariance_gate(
        reference_force_json,
        rotated_force_json,
        rotation_deg,
        relative_tolerance,
    )


@mcp.tool()
def motor_force_report_method_metadata_gate(
    report_json: str,
    relative_tolerance: float = 2.0e-2,
) -> str:
    """Gate a force report using independent methods and action-reaction.

    The JSON report records ``force_unit``, ``component_frame``, at least two
    ``methods`` (each with ``family``, ``domain`` and ``vector``), plus
    ``action_force`` and ``reaction_force``. Method vectors and Newton's-third-
    law closure must agree within ``relative_tolerance``.
    """
    return build_force_report_method_metadata_gate(report_json, relative_tolerance)


@mcp.tool()
def motor_phase_flux_park_alignment_gate(
    mechanical_angles_deg_json: str,
    phase_flux_wb_json: str,
    pole_pairs: int,
    q_relative_tolerance: float = 3.0e-2,
    d_ripple_relative_tolerance: float = 2.0e-2,
) -> str:
    """Gate a PM-only three-phase flux sweep in the rotating d/q frame."""
    return build_phase_flux_park_alignment_gate(
        mechanical_angles_deg_json,
        phase_flux_wb_json,
        pole_pairs,
        q_relative_tolerance,
        d_ripple_relative_tolerance,
    )


@mcp.tool()
def motor_ipm_two_run_ldlq_gate(summary_json: str) -> str:
    """Gate same-angle PM-only/current-on runs and extract ``Ld``/``Lq``.

    Both runs must provide explicit angle grids and canonical phase order. The
    gate subtracts PM-only phase flux before Park projection and rejects stale
    total-flux shortcuts, phase-order mismatches, or implausible saliency.
    """
    return build_ipm_two_run_ldlq_gate(summary_json)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_motor_simulation(motor_type: str = "pmsm") -> str:
    """Set up a new electric-machine simulation."""
    t = motor_type.strip().lower()
    base = (
        f"Set up a {motor_type} simulation using the radia-mcp motor toolchain.\n\n"
        "Recommended path:\n"
        "1. motor_femm_transient('lab_recommendation') for the **lab-canonical**\n"
        "   transient approach (Lange-Henrotte-Hameyer 2009).\n"
        "2. motor_onelab('overview') for ONELAB ElectricMachines layout\n"
        "   (reference geometry templates).\n"
        "3. motor_onelab('motor_types') to confirm Region group names.\n"
        "4. motor_onelab('formulation') for the 2D weak form.\n"
        "5. motor_onelab('ngsolve_xlate') for the NGSolve port.\n"
        "6. motor_onelab('analysis_modes') to pick Static / TD / FD.\n\n"
    )
    if t == "pmsm":
        base += (
            "PMSM-specific:\n"
            "- Rotor: surface or interior permanent magnets (Br vector / region).\n"
            "- Anti-periodicity factor: NbrPolePairs (1 sector = 1 pole-pair).\n"
            "- Cogging torque map: motor_onelab('analysis_modes') Static, sweep θ.\n"
            "- d/q inductance: Static sweep I_d, I_q.\n"
        )
    elif t == "synrm":
        base += (
            "SynRM-specific:\n"
            "- No PM, no rotor coil — pure reluctance machine.\n"
            "- Dominant design objectives: T_ave, T_rip.\n"
            "- Use motor_topology_optimization('wakao_ae_ls') for shape design.\n"
            "- motor_topology_optimization('liu_thesis_application') for end-to-end recipe.\n"
        )
    elif t in ("im", "induction"):
        base += (
            "Induction-motor-specific:\n"
            "- Use Frequency-Domain analysis with slip-frequency for steady state.\n"
            "- Use Time-Domain for transient run-up.\n"
            "- Cage bars are RotorC region; end-rings via external circuit.\n"
            "- For broken-bar fault, see ONELAB `im` bundle.\n"
        )
    elif t == "srm":
        base += (
            "SRM-specific:\n"
            "- Doubly salient, no PM, no rotor coil, no cage.\n"
            "- Phase-pulsed excitation (square-wave or current-controlled).\n"
            "- Strongly non-linear ν(B) — use Newton with damping.\n"
        )
    elif t == "wfsm":
        base += (
            "WFSM-specific:\n"
            "- Rotor field winding (DC) + damper cage.\n"
            "- For sub-transient analysis, treat damper bars as RotorC.\n"
        )
    else:
        base += (
            f"Custom motor type '{motor_type}':\n"
            "- Map your rotor source to one of: Rotor_Magnets, RotorC, Rotor_Inds.\n"
            "- Choose anti-periodicity factor from pole count.\n"
        )
    return base


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-motor',
    description='Motor analysis: ONELAB transient, Hollaus effective material (lamination), Wakao autoencoder topology, Kaimori-Mifune Darwin TD',
    subpackage='radia_mcp.motor',
    related_servers=["electromagnet", "topology-optimization", "magnetic-materials"],
    optional_deps=["radia", "ngsolve"],
)

apply_tool_contract(
    mcp,
    server_name="mcp-server-motor",
    version="1.4.19",
)


def main():
    if "--selftest" in sys.argv:
        print("Motor MCP server self-test:")
        from .onelab_knowledge import SECTIONS as O_SEC
        from .topology_opt_knowledge import SECTIONS as T_SEC
        from .darwin_model_knowledge import SECTIONS as D_SEC
        from .femm_transient_knowledge import SECTIONS as F_SEC
        from .henrotte_lineage_knowledge import SECTIONS as H_SEC
        from .hollaus_eddy_knowledge import SECTIONS as E_SEC
        from .tritool_cross_reference_knowledge import SECTIONS as X_SEC
        from .deck_bridge_knowledge import SECTIONS as M_SEC
        from .age_quality_knowledge import SECTIONS as A_SEC
        from .validation_lanes_knowledge import SECTIONS as L_SEC
        from .angle_periodic_rom_knowledge import SECTIONS as R_SEC
        for k in O_SEC:
            r = motor_onelab(k)
            print(f"  motor_onelab({k!r}): {len(r)} chars")
            assert len(r) > 100, f"ONELAB topic {k} too short"
        for k in T_SEC:
            r = motor_topology_optimization(k)
            print(f"  motor_topology_optimization({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Topology-opt topic {k} too short"
        for k in D_SEC:
            r = motor_darwin_model(k)
            print(f"  motor_darwin_model({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Darwin topic {k} too short"
        for k in F_SEC:
            r = motor_femm_transient(k)
            print(f"  motor_femm_transient({k!r}): {len(r)} chars")
            assert len(r) > 100, f"FEMM transient topic {k} too short"
        for k in H_SEC:
            r = motor_henrotte_lineage(k)
            print(f"  motor_henrotte_lineage({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Henrotte lineage topic {k} too short"
        for k in E_SEC:
            r = motor_hollaus_eddy(k)
            print(f"  motor_hollaus_eddy({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Hollaus eddy topic {k} too short"
        for k in X_SEC:
            r = motor_tritool_cross_reference(k)
            print(f"  motor_tritool_cross_reference({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Tri-tool topic {k} too short"
        for k in M_SEC:
            r = motor_deck_bridge(k)
            print(f"  motor_deck_bridge({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Deck bridge topic {k} too short"
        for k in A_SEC:
            r = motor_age_quality(k)
            print(f"  motor_age_quality({k!r}): {len(r)} chars")
            assert len(r) > 100, f"AGE quality topic {k} too short"
        for k in L_SEC:
            r = motor_validation_lanes(k)
            print(f"  motor_validation_lanes({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Motor validation lane topic {k} too short"
        for k in R_SEC:
            r = motor_angle_periodic_rom(k)
            print(f"  motor_angle_periodic_rom({k!r}): {len(r)} chars")
            assert len(r) > 100, f"Motor angle-ROM topic {k} too short"
        assert "cycle basis" in motor_angle_periodic_rom("face_policy")
        assert "positive-real CLN" in motor_angle_periodic_rom("time_domain")
        ports_text = motor_angle_periodic_rom("ports")
        assert "FMI source boundary" in ports_text
        assert "packaged FMU" in ports_text
        bridge = motor_deck_bridge("insufficiency_audit")
        assert "gold_numeric_invariant" in bridge
        assert "radia-motor" in motor_deck_bridge("radia_strengthening_queue")
        assert "not a full" in motor_deck_bridge("jmag_coverage_reality")
        assert "NGSolve AGE" in motor_deck_bridge("age_vs_field_strategy")
        assert "gold_age_invariant" in motor_age_quality("publication_policy")
        assert "validation_test/radia_mcp/test_airgap_eddy_machine.py" in motor_age_quality("gate_matrix")
        assert "HDiv-MMM + HCurl eddy-bubble" in motor_validation_lanes("lane_matrix")
        assert "NGSolve+AGE" in motor_validation_lanes("overview")
        assert "product_local_reference" in motor_validation_lanes("source_policy")
        dual_catalog = motor_dual_lane_training_catalog("all")
        print(f"  motor_dual_lane_training_catalog('all'): {len(dual_catalog)} chars")
        assert "radia-motor-age" in dual_catalog
        assert "radia-motor-mmm-eddy" in dual_catalog
        dual_gate = json.loads(motor_dual_lane_training_gate())
        assert dual_gate["status"] == "PASS"
        assert dual_gate["count"] >= 50
        assert not dual_gate["forbidden_hits"]
        dual_route = motor_dual_lane_training_route("SRM static torque")
        assert "srm_static_torque_curve" in dual_route
        assert "hdiv_mmm_hcurl_eddy_bubble" in dual_route
        outer_route = motor_dual_lane_training_route("BLDC outer rotor polarity")
        assert "bldc_outer_rotor_polarity" in outer_route
        triple_plan = motor_triple_check_plan("IPM hairpin motor flux linkage and MTPA")
        print(f"  motor_triple_check_plan('IPM ...'): {len(triple_plan)} chars")
        assert "elf_motor_hybrid_router" in triple_plan
        assert "hdiv_mmm_hcurl_eddy_bubble" in triple_plan
        assert "ngsolve_age" in triple_plan
        assert "primary required lanes" in triple_plan
        lane_tpl = motor_validation_lane_template("hdiv_mmm_hcurl_eddy_bubble")
        assert "hdiv_mmm_operator_contract" in lane_tpl
        assert "hcurl_eddy_bubble_contract" in lane_tpl
        hdiv_selftest_artifact = {
            "schema_version": "radia-motor-validation-artifact/v1",
            "timestamp_utc": "2026-07-03T00:00:00Z",
            "radia_version": "selftest",
            "motor_validation_lane": "hdiv_mmm_hcurl_eddy_bubble",
            "reference_source_class": "analytic_reference",
            "observable_family": "pickup_flux",
            "case_count": 1,
            "status": "pass",
            "tolerances": {"max_abs_relative_error": 1.0e-2},
            "metrics": {"max_abs_relative_error": 1.0e-4},
            "timing_breakdown_s": {"solve": 0.01},
            "artifact_feedback": {
                "status": "candidate",
                "public_lesson": "HDiv-MMM/HCurl eddy-bubble artifact is complete.",
            },
            "coupling_design_status": "solver_validated",
            "hdiv_mmm_operator_contract": {
                "space": "HDiv",
                "observable": "pickup_flux",
            },
            "hcurl_eddy_bubble_contract": {
                "space": "HCurl",
                "basis": "eddy_bubble",
            },
            "coupling_operator_contract": {
                "blocks": ["hdiv_mmm", "hcurl_eddy_bubble"],
            },
            "shared_mesh_material_identity": {
                "geometry_sha256": "1" * 64,
                "material_sha256": "2" * 64,
                "excitation_sha256": "3" * 64,
            },
            "solver_ready_artifact": {
                "artifact_id": "hdiv_mmm_hcurl_eddy_bubble_selftest_v1",
                "verification": ["selftest mixed-system execution"],
            },
        }
        age_selftest_artifact = {
            "schema_version": "radia-motor-validation-artifact/v1",
            "timestamp_utc": "2026-07-03T00:00:00Z",
            "radia_version": "selftest",
            "motor_validation_lane": "ngsolve_age",
            "reference_source_class": "analytic_reference",
            "observable_family": "torque",
            "case_count": 1,
            "status": "pass",
            "tolerances": {"torque_relative_error": 1.0e-2},
            "metrics": {"torque_relative_error": 1.0e-4},
            "timing_breakdown_s": {"solve": 0.01},
            "artifact_feedback": {
                "status": "candidate",
                "public_lesson": "AGE torque lane selftest artifact is complete.",
            },
            "age_gate_ids": ["age_rotation_torque"],
            "pytest_targets": ["validation_test/radia_mcp/test_airgap_machine_rotation.py"],
            "shared_mesh_material_identity": {
                "geometry_sha256": "1" * 64,
                "material_sha256": "2" * 64,
                "excitation_sha256": "3" * 64,
            },
            "solver_ready_artifact": {
                "artifact_id": "ngsolve_age_selftest_v1",
                "verification": ["selftest AGE execution"],
            },
        }
        triple_gate = motor_triple_check_artifact_gate(
            json.dumps(
                {
                    "schema_version": "radia-motor-triple-check-artifact/v1",
                    "goal": "selftest",
                    "source_mcp_seed": {
                        "source_mcp_calls": ["elf_motor_hybrid_router('selftest')"],
                        "representative_public_decks": [
                            "application/motor/spm_surface_pm_10/spm001/spm001.mai"
                        ],
                    },
                    "lane_artifacts": {
                        "hdiv_mmm_hcurl_eddy_bubble": hdiv_selftest_artifact,
                        "ngsolve_age": age_selftest_artifact,
                    },
                    "mcp_feedback": {
                        "public_status": "verified",
                        "public_summary": "AGE and mixed-system metadata are complete.",
                        "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
                        "verification": ["selftest"],
                    },
                }
            )
        )
        assert "validated supported solver check: `True`" in triple_gate
        assert "validated dual solver check: `True`" in triple_gate
        assert "accepted for supported MCP learning: `True`" in triple_gate
        assert "accepted for MCP RFC learning: `False`" in triple_gate
        assert "accepted for MCP learning: `True`" in triple_gate
        gate = motor_validation_artifact_gate(
            json.dumps(
                age_selftest_artifact
            ),
            "ngsolve_age",
        )
        assert "accepted for MCP learning: `True`" in gate
        age_plan = motor_age_validation_plan("IPM hairpin MTPA field weakening")
        print(f"  motor_age_validation_plan('IPM ...'): {len(age_plan)} chars")
        assert "dq_control_layer" in age_plan
        assert "validation_test/radia_mcp/test_field_weakening.py" in age_plan
        im_plan = motor_age_validation_plan("induction cage slip loss")
        print(f"  motor_age_validation_plan('induction ...'): {len(im_plan)} chars")
        assert "age_eddy_machine" in im_plan
        assert "validation_test/radia_mcp/test_motor_induction_coupling.py" in im_plan
        field_check = motor_field_quick_check(motor_type="ipm", electrical_angle_deg=25)
        print(f"  motor_field_quick_check('ipm'): {len(field_check)} chars")
        assert "2D magnetic-circuit/BEM-like motor quick check" in field_check
        assert "ld_lq" in field_check
        assert "not a production solver" in field_check
        route = motor_validation_router("IPM hairpin motor flux linkage and MTPA")
        print(f"  motor_validation_router('IPM ...'): {len(route)} chars")
        assert "application/motor/emdlab_ipm_hairpin_10" in route
        assert "ngsolve_usage(\"mtpa\")" in route
        p = new_motor_simulation("synrm")
        print(f"  new_motor_simulation('synrm'): {len(p)} chars")
        assert "wakao_ae_ls" in p
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
