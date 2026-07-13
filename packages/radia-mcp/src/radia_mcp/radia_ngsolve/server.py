"""
Radia + NGSolve Unified MCP Server

Provides tools for both Radia (field/PEEC C++ core) and NGSolve (FEM/BEM):
- Unified linting (33 rules: Radia API + NGSolve FEM + BEM + PEEC)
- Radia C++ library usage (field computation, materials, solver)
- NGSolve FEM usage (22 topics: EM formulations, axisymmetric, materials)
- ngsolve.bem (BEM operators, inductance extraction)
- radia.sparsesolv_ngsolve (Compact AMS/COCR/ICCG preconditioners)
- Kelvin transformation for open boundary FEM
- md2html converter documentation

App-specific knowledge is in separate servers:
- mcp-server-ih: Induction heating (SIBC, ESIM, Karl iteration)
- mcp-server-cubit: Cubit scripting, mesh export
- mcp-server-gmsh: GMSH post-processing

Usage:
    mcp-server-radia-ngsolve              # Start MCP server (stdio transport)
    mcp-server-radia-ngsolve --selftest   # Run self-test
"""

import os
import json
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool
from ..common.learning_quality import build_balanced_learning_profile
from .rf_sweep_artifact_gate import rf_sweep_artifact_summary_gate as _rf_sweep_artifact_summary_gate
from .cq_urn import cq_response_reality_gate as _cq_response_reality_gate
from .cq_scattering_arrival_gate import cq_scattering_arrival_gate as _cq_scattering_arrival_gate
from .physics_result_preflight_gate import physics_result_preflight_gate as _physics_result_preflight_gate
from .field_profile_gate import (
    dual_formulation_symmetric_field_profile_gate as _dual_formulation_symmetric_field_profile_gate,
    symmetric_complex_field_curve_gate as _symmetric_complex_field_curve_gate,
    symmetric_axial_field_profile_gate as _symmetric_axial_field_profile_gate,
)
from .acoustic_kernel_gate import (
    helmholtz_double_layer_low_frequency_gate as _helmholtz_double_layer_low_frequency_gate,
)
from .terminal_source_sweep_gate import cyclic_terminal_source_sweep_gate as _cyclic_terminal_source_sweep_gate
from .terminal_phasor_balance_gate import (
    cyclic_terminal_phasor_balance_gate as _cyclic_terminal_phasor_balance_gate,
)
from .three_phase_winding_power_gate import (
    three_phase_winding_power_balance_gate as _three_phase_winding_power_balance_gate,
)
from .cogging_periodicity_gate import cogging_torque_periodicity_gate as _cogging_torque_periodicity_gate
from .nonlinear_actuator_gate import nonlinear_actuator_saturation_knee_gate as _nonlinear_actuator_saturation_knee_gate
from .source_free_static_gate import source_free_static_null_solution_gate as _source_free_static_null_solution_gate
from .harmonic_force_triplet_gate import (
    harmonic_magnetic_force_triplet_closure_gate as _harmonic_magnetic_force_triplet_closure_gate,
)
from .magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate as _magnetic_force_method_profile_gate,
)
from .harmonic_port_identity_gate import (
    harmonic_current_port_power_energy_identity_gate as _harmonic_current_port_power_energy_identity_gate,
)
from .periodic_pm_machine_gate import (
    periodic_unwrapped_pm_machine_replay_gate as _periodic_unwrapped_pm_machine_replay_gate,
)
from .permanent_magnet_recoil_state_gate import (
    permanent_magnet_recoil_state_gate as _permanent_magnet_recoil_state_gate,
)
from .eddy_levitation_force_gate import (
    linear_eddy_levitation_force_gate as _linear_eddy_levitation_force_gate,
)
from .motion_coupled_levitation_gate import (
    motion_coupled_eddy_levitation_transient_gate as _motion_coupled_eddy_levitation_transient_gate,
)
from .harmonic_zero_net_circuit_gate import (
    harmonic_zero_net_circuit_gate as _harmonic_zero_net_circuit_gate,
)
from .moving_conductor_brake_gate import moving_conductor_eddy_brake_gate as _moving_conductor_eddy_brake_gate
from .rotating_conductor_transient_gate import rotating_conductor_transient_gate as _rotating_conductor_transient_gate
from .linear_magnetization_scaling_gate import linear_magnetization_scaling_gate as _linear_magnetization_scaling_gate
from .linear_axisymmetric_circuit_energy_gate import (
    linear_axisymmetric_circuit_energy_gate as _linear_axisymmetric_circuit_energy_gate,
)
from .manual_auto_mixed_mesh_gate import (
    manual_auto_mixed_mesh_preservation_gate as _manual_auto_mixed_mesh_preservation_gate,
)
from .two_winding_frequency_gate import (
    two_winding_frequency_faraday_gate as _two_winding_frequency_faraday_gate,
)
from .conductive_shield_frequency_gate import (
    magnetic_conductive_shield_frequency_gate as _magnetic_conductive_shield_frequency_gate,
)
from .cylindrical_conductor_skin_gate import (
    cylindrical_conductor_skin_bessel_gate as _cylindrical_conductor_skin_bessel_gate,
)
from .conductive_network_monotonicity_gate import (
    conductive_network_resistance_monotonicity_gate as _conductive_network_resistance_monotonicity_gate,
)
from .autodiff_harmonic_balance_gate import (
    autodiff_harmonic_balance_convergence_gate as _autodiff_harmonic_balance_convergence_gate,
)
from .hall_effect_gate import (
    hall_effect_transverse_voltage_gate as _hall_effect_transverse_voltage_gate,
)
from .single_loop_normalized_field_gate import (
    single_loop_source_normalized_field_gate as _single_loop_source_normalized_field_gate,
)
from .coupled_cq_refinement_gate import (
    coupled_cq_refinement_gate as _coupled_cq_refinement_gate,
)
from .coil_self_resonance_gate import (
    coil_self_resonance_sweep_gate as _coil_self_resonance_sweep_gate,
)
from .passive_axial_bearing_gate import (
    passive_axial_bearing_stiffness_gate as _passive_axial_bearing_stiffness_gate,
)
from .radial_bearing_force_gate import (
    radial_bearing_force_symmetry_gate as _radial_bearing_force_symmetry_gate,
)
from .complex_field_maximum_gate import complex_vector_field_maximum_gate as _complex_vector_field_maximum_gate
from .one_port_vi_s_gate import one_port_vi_s_impedance_gate as _one_port_vi_s_impedance_gate
from .force_position_profile_gate import force_position_profile_gate as _force_position_profile_gate
from .force_coenergy_gate import force_coenergy_displacement_gate as _force_coenergy_displacement_gate
from .rotational_time_axis_gate import rotational_kinematics_time_axis_gate as _rotational_kinematics_time_axis_gate
from .inductance_matrix_gate import inductance_matrix_family_gate as _inductance_matrix_family_gate
from .sphere_mesh_convergence_gate import (
    linear_sphere_geometry_convergence_gate as _linear_sphere_geometry_convergence_gate,
)
from .leakage_inductance_closure_gate import (
    leakage_inductance_closure_gate as _leakage_inductance_closure_gate,
)
from .parallel_wire_force_refinement_gate import (
    parallel_wire_force_refinement_gate as _parallel_wire_force_refinement_gate,
)
from .capacitance_identity_gate import (
    two_conductor_capacitance_identity_gate as _two_conductor_capacitance_identity_gate,
)
from .capacitance_matrix_gate import two_conductor_capacitance_matrix_gate as _two_conductor_capacitance_matrix_gate
from .multiconductor_capacitance_gate import (
    multiconductor_capacitance_cross_formulation_gate as _multiconductor_capacitance_cross_formulation_gate,
)
from .impedance_sweep_gate import multiport_impedance_sweep_gate as _multiport_impedance_sweep_gate
from .radar_range_rcs_gate import radar_range_rcs_profile_gate as _radar_range_rcs_profile_gate
from .radar_range_angle_gate import radar_range_angle_localization_gate as _radar_range_angle_localization_gate
from .hmatrix_scaling_gate import hmatrix_compression_scaling_gate as _hmatrix_compression_scaling_gate
from .acoustic_duct_band_gap_gate import acoustic_duct_band_gap_gate as _acoustic_duct_band_gap_gate
from .force_error_convergence_gate import (
    dual_formulation_force_error_convergence_gate as _dual_formulation_force_error_convergence_gate,
)
from .voice_coil_gate import voice_coil_force_flux_sweep_gate as _voice_coil_force_flux_sweep_gate
from .linear_induction_gate import linear_induction_frequency_sweep_gate as _linear_induction_frequency_sweep_gate
from .conductor_frequency_gate import (
    homogenized_bundle_impedance_comparison_gate as _homogenized_bundle_impedance_comparison_gate,
    opposed_busbar_skin_force_gate as _opposed_busbar_skin_force_gate,
    twin_conductor_skin_effect_frequency_gate as _twin_conductor_skin_effect_frequency_gate,
)
from .transient_conductor_replay_gate import (
    transient_conductor_replay_identity_gate as _transient_conductor_replay_identity_gate,
)
from .adjoint_gradient_gate import adjoint_gradient_scaling_gate as _adjoint_gradient_scaling_gate
from .one_port_power_gate import one_port_power_balance_gate as _one_port_power_balance_gate
from .fsi_scattering_invariants_gate import (
    fsi_scattering_invariants_gate as _fsi_scattering_invariants_gate,
)
from .inductance_energy_gate import inductance_energy_mutual_gate as _inductance_energy_mutual_gate
from .loss_temperature_coupling_gate import loss_temperature_coupling_gate as _loss_temperature_coupling_gate
from .transient_coupled_coil_gate import (
    transient_coupled_coil_response_gate as _transient_coupled_coil_response_gate,
)
from .source_off_relaxation_gate import (
    source_off_linear_relaxation_gate as _source_off_linear_relaxation_gate,
)
from .nonlinear_bh_curve_gate import (
    nonlinear_bh_piecewise_material_gate as _nonlinear_bh_piecewise_material_gate,
)
from .skin_effect_adaptive_gate import skin_effect_adaptive_energy_loss_gate as _skin_effect_adaptive_energy_loss_gate
from .global_local_optimization_gate import global_local_optimization_replay_gate as _global_local_optimization_replay_gate
from .eddy_loss_formulation_gate import alternate_eddy_loss_formulation_gate as _alternate_eddy_loss_formulation_gate
from .open_boundary_magnetostatic_gate import (
    magnetostatic_open_boundary_equivalence_gate as _magnetostatic_open_boundary_equivalence_gate,
)
from .pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate as _pwm_controlled_motor_loss_gate,
)
from .material_contrast_force_gate import (
    material_contrast_force_gate as _material_contrast_force_gate,
)
from .two_body_force_magnitude_gate import (
    two_body_force_magnitude_replay_gate as _two_body_force_magnitude_replay_gate,
)
from .static_field_shim_family_gate import (
    static_field_shim_family_gate as _static_field_shim_family_gate,
)
from .energy_budgeted_trace_kkt_gate import (
    energy_budgeted_trace_kkt_gate as _energy_budgeted_trace_kkt_gate,
)
from .finite_solenoid_surface_current_gate import (
    finite_solenoid_surface_current_gate as _finite_solenoid_surface_current_gate,
)
from .linked_study_noop_gate import (
    linked_study_silent_noop_gate as _linked_study_silent_noop_gate,
)
from .two_port_power_gate import (
    reciprocal_two_port_power_sweep_gate as _reciprocal_two_port_power_sweep_gate,
)
from .fem_bem_capstone_gate import (
    fem_bem_capstone_suite_gate as _fem_bem_capstone_suite_gate,
)
from .helmholtz_dual_formulation_gate import (
    helmholtz_dual_formulation_axis_gate as _helmholtz_dual_formulation_axis_gate,
)
from .lossy_dielectric_power_gate import (
    lossy_dielectric_complex_power_refinement_gate as _lossy_power_refinement_gate,
)
from .heterogeneous_current_flow_gate import (
    heterogeneous_current_flow_p1_reintegration_gate as _heterogeneous_current_flow_gate,
)
from .thermal_robin_balance_gate import (
    thermal_robin_boundary_balance_gate as _thermal_robin_boundary_balance_gate,
)
from .hysteresis_minor_loop_gate import (
    hysteresis_minor_loop_replay_gate as _hysteresis_minor_loop_replay_gate,
)
from .heterogeneous_mesh_replay_gate import (
    heterogeneous_part_mesh_replay_gate as _heterogeneous_part_mesh_replay_gate,
)
from .two_terminal_dc_conduction_gate import (
    two_terminal_dc_conduction_power_gate as _two_terminal_dc_conduction_power_gate,
)
from .rwg_hcurl_trace_gate import (
    rwg_hcurl_trace_consistency_gate as _rwg_hcurl_trace_consistency_gate,
)
from .hartmann_profile_gate import hartmann_profile_gate as _hartmann_profile_gate

from .rules import ALL_RULES
from .knowledge.radia import get_radia_documentation
from .knowledge.md2html import get_md2html_documentation
from .knowledge.ngsolve import get_ngsolve_documentation
from .knowledge.sparsesolv import get_sparsesolv_documentation
from .knowledge.kelvin import get_kelvin_documentation
from .knowledge.kelvin_identify_post_hoc import (
    get_post_hoc_documentation as _get_kelvin_identify_post_hoc_doc)
from .knowledge.axifem import get_axifem_documentation
from .knowledge.ngsbem_inductance import get_ngsbem_inductance_documentation
from .knowledge.peec_inductance import get_peec_inductance_documentation
from .knowledge.esim import get_esim_documentation
from .knowledge.panel_gui_pitfalls import get_panel_gui_pitfalls
from .knowledge.analytical_formulas import get_analytical_formulas_documentation
from .knowledge.force_validation import get_force_validation_documentation
from .knowledge.install_deploy import get_install_deploy_documentation
from .knowledge.release_workflow import get_release_workflow_documentation
from .knowledge.standalone_panels import get_standalone_panels_documentation
from .knowledge.loop_learning import get_loop_learning_documentation
from .knowledge.basis_functions import get_basis_functions_documentation
from .knowledge.taskmanager import get_taskmanager_knowledge
from .knowledge.cln_sibc_orthogonal import (
    get_cln_sibc_orthogonal_documentation,
    get_cln_sibc_orthogonal_section,
)
from .knowledge.cln_3d import (
    get_cln_3d_documentation,
    get_cln_3d_notebook,
)
from .knowledge.bem_cln import get_bem_cln_documentation
from .knowledge.cln_sphere_dd import get_cln_sphere_dd_documentation
from .knowledge.hdiv_vim import get_hdiv_vim_documentation
from .knowledge.femm_parity import get_femm_parity_documentation
from .knowledge.fem_bem_schur import get_fem_bem_schur_documentation
from .knowledge.airgap_motor_workflow import get_airgap_motor_workflow_documentation
from .knowledge.dtn_coarse_mesh import get_dtn_coarse_mesh_documentation
from .knowledge.urn import get_urn_documentation, urn_fit_from_csv
from .acoustics import (
    acoustic_fembem_cross_learnings as _acoustic_fembem_cross_learnings,
)
from ..matlab_acoustic_fembem import (
    matlab_acoustic_fembem_agent_guide as _matlab_acoustic_fembem_agent_guide,
)
from .gmsh_post_spec import get_gmsh_post_spec
from .panel_describer import (
    find_panel_file as _find_panel_file,
    parse_panel_file as _parse_panel_file,
    describe_panel_jp as _describe_panel_jp,
    widget_locations as _widget_locations,
)

# NOTE: induction_heating_knowledge is in mcp-server-ih (not here)

mcp = FastMCP("mcp-server-radia-ngsolve")

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
REPO_ROOT = PACKAGE_ROOT.parent.parent if PACKAGE_ROOT.parent.name == "packages" else PACKAGE_ROOT
PROJECT_ROOT = PACKAGE_ROOT
MAX_LINT_FILE_BYTES = 2_000_000
MAX_LINT_DIRECTORY_FILES = 300


def _is_relative_to(path: Path, base: Path) -> bool:
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False


def _lint_allowed_roots() -> tuple[Path, ...]:
    """Roots that MCP lint tools may inspect.

    Public MCP operation must not become a local-file oracle.  By default the
    linter can inspect only this package.  Advanced local use can opt in to
    extra roots with RADIA_MCP_LINT_ROOTS, separated by os.pathsep.
    """
    roots = [PACKAGE_ROOT]
    extra = os.environ.get("RADIA_MCP_LINT_ROOTS", "")
    for raw in extra.split(os.pathsep):
        raw = raw.strip()
        if raw:
            roots.append(Path(raw).expanduser())
    resolved: list[Path] = []
    for root in roots:
        try:
            resolved.append(root.resolve())
        except OSError:
            continue
    return tuple(dict.fromkeys(resolved))


def _resolve_lint_path(value: str) -> Path:
    p = Path(value).expanduser()
    if p.is_absolute():
        return p.resolve(strict=False)

    parts_lower = tuple(part.lower() for part in p.parts)
    if parts_lower[:2] == ("packages", "radia-mcp"):
        return (REPO_ROOT / p).resolve(strict=False)
    if parts_lower[:1] == ("radia-mcp",):
        return (PACKAGE_ROOT.parent / p).resolve(strict=False)

    candidates = [PROJECT_ROOT / p, REPO_ROOT / p, PACKAGE_ROOT.parent / p]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve(strict=False)
    return candidates[0].resolve(strict=False)


def _lint_path_allowed(path: Path) -> bool:
    return any(_is_relative_to(path, root) for root in _lint_allowed_roots())


def _lint_display_path(path: Path) -> str:
    try:
        resolved = path.resolve(strict=False)
    except OSError:
        resolved = path
    if _is_relative_to(resolved, REPO_ROOT):
        return "repo:/" + resolved.relative_to(REPO_ROOT).as_posix()
    for idx, root in enumerate(_lint_allowed_roots()):
        if _is_relative_to(resolved, root):
            rel = resolved.relative_to(root).as_posix()
            prefix = "package-root" if root == PACKAGE_ROOT else f"allowed-root-{idx}"
            return f"{prefix}:/{rel}" if rel != "." else f"{prefix}:/"
    return "<outside allowed lint roots>"


def _lint_display_roots() -> str:
    labels = []
    for idx, root in enumerate(_lint_allowed_roots()):
        if _is_relative_to(root, REPO_ROOT):
            labels.append("repo:/" + root.relative_to(REPO_ROOT).as_posix())
        else:
            labels.append("package-root" if root == PACKAGE_ROOT else f"allowed-root-{idx}")
    return ", ".join(dict.fromkeys(labels))


def _lint_denied(path: Path) -> str:
    return (
        f"Error: Access denied: {_lint_display_path(path)}. "
        "Allowed lint roots are: "
        f"{_lint_display_roots()}. Set RADIA_MCP_LINT_ROOTS to opt in to another project root."
    )


def _lint_file(filepath: str) -> list[dict]:
    """Run all lint rules on a single file."""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except (OSError, IOError) as e:
        return [{'line': 0, 'severity': 'ERROR', 'rule': 'read-error',
                 'message': f'Cannot read file: {e}'}]

    findings = []
    for rule_fn in ALL_RULES:
        findings.extend(rule_fn(filepath, lines))

    severity_order = {'CRITICAL': 0, 'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'INFO': 4, 'ERROR': -1}
    findings.sort(key=lambda f: (severity_order.get(f['severity'], 9), f['line']))
    return findings


def _format_findings(filepath: str, findings: list[dict]) -> str:
    """Format findings for display."""
    if not findings:
        return f"[OK] {filepath}: No issues found."

    lines = [f"[{len(findings)} issue(s)] {filepath}:"]
    for f in findings:
        lines.append(
            f"  L{f['line']:>4d} [{f['severity']}] {f['rule']}: {f['message']}"
        )
    return '\n'.join(lines)


@mcp.tool()
def balanced_mcp_learning_profile() -> dict:
    """Return the ten-stage equal public/source MCP learning contract."""

    return build_balanced_learning_profile(
        "radia-mcp", "radia-mcp", "source-tool MCP or shared open-tool owner"
    )


@mcp.tool()
def lint_radia_script(filepath: str) -> str:
    """
    Lint a single Python script for Radia + NGSolve convention violations.

    Radia checks:
    - ObjBckg called with list instead of callable (CRITICAL)
    - Missing UtiDelAll cleanup (HIGH)
    - Removed APIs: FldUnits, FldBatch, old solver params (HIGH)

    NGSolve checks:
    - BEM on HDivSurface without .Trace() (CRITICAL)
    - Circular SIBC using jv instead of iv (CRITICAL)
    - HCurl magnetostatics without nograds=True (HIGH)
    - Eddy current FE space missing complex=True (HIGH)
    - EFIE V_LL term with wrong (minus) sign (HIGH)
    - PEEC P/(jw) low-frequency breakdown (HIGH)
    - BDDC preconditioner registered after assembly (MODERATE)
    - Overwriting x/y/z coordinate variables (MODERATE)
    - Direct .vec assignment without .data (MODERATE)
    - 2D OCC geometry without dim=2 (MODERATE)
    - CG on A-Omega saddle-point system (MODERATE)
    - Kelvin domain without bonus_intorder (MODERATE)
    - VectorH1 for electromagnetic fields (MODERATE)
    - PINVIT/LOBPCG without gradient projection (MODERATE)
    - Joule heat missing Conj() for complex fields (MODERATE)
    - PEEC n_seg too low for coupling accuracy (MODERATE)
    - Classical EFIE 1/kappa^2 low-frequency breakdown (MODERATE)
    - BEM GenerateMesh without curvaturesafety (MODERATE)
    - TaskManager with BEM non-determinism (LOW)

    Args:
        filepath: Absolute or relative path to the Python file to check.
    """
    p = _resolve_lint_path(filepath)

    if not _lint_path_allowed(p):
        return _lint_denied(p)
    if not p.exists():
        return f"Error: File not found: {_lint_display_path(p)}"
    if not p.suffix == '.py':
        return f"Error: Not a Python file: {_lint_display_path(p)}"
    try:
        if p.stat().st_size > MAX_LINT_FILE_BYTES:
            return (
                f"Error: File too large for MCP lint: {_lint_display_path(p)} "
                f"({p.stat().st_size} bytes > {MAX_LINT_FILE_BYTES})."
            )
    except OSError as exc:
        return f"Error: Cannot stat file: {_lint_display_path(p)} ({exc.__class__.__name__})."

    findings = _lint_file(str(p))
    return _format_findings(_lint_display_path(p), findings)


@mcp.tool()
def lint_radia_directory(directory: str = ".") -> str:
    """
    Lint all Python scripts in a directory for NGSolve convention violations.

    Recursively scans .py files and reports findings grouped by file.

    Args:
        directory: Directory path (default: current directory).
    """
    d = _resolve_lint_path(directory)
    if not _lint_path_allowed(d):
        return _lint_denied(d)
    if not d.exists():
        return f"Error: Directory not found: {_lint_display_path(d)}"

    py_files = []
    for py_file in sorted(d.rglob("*.py")):
        resolved = py_file.resolve(strict=False)
        if not _lint_path_allowed(resolved):
            return _lint_denied(resolved)
        try:
            if resolved.stat().st_size > MAX_LINT_FILE_BYTES:
                continue
        except OSError:
            continue
        py_files.append(resolved)
        if len(py_files) > MAX_LINT_DIRECTORY_FILES:
            return (
                f"Error: Too many Python files for MCP lint: more than "
                f"{MAX_LINT_DIRECTORY_FILES} under {_lint_display_path(d)}."
            )
    if not py_files:
        return f"No Python files found in {directory}."

    total_findings = 0
    file_results = []
    summary_by_severity = {'CRITICAL': 0, 'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'INFO': 0}

    for py_file in py_files:
        findings = _lint_file(str(py_file))
        if findings:
            total_findings += len(findings)
            file_results.append(_format_findings(_lint_display_path(py_file), findings))
            for f in findings:
                sev = f['severity']
                if sev in summary_by_severity:
                    summary_by_severity[sev] += 1

    output_parts = [
        f"NGSolve Lint Report: {len(py_files)} files scanned, {total_findings} issues found.",
        "",
        f"Summary: {summary_by_severity['CRITICAL']} CRITICAL, "
        f"{summary_by_severity['HIGH']} HIGH, "
        f"{summary_by_severity['MODERATE']} MODERATE, "
        f"{summary_by_severity['LOW']} LOW",
        "",
    ]

    if file_results:
        output_parts.append("=" * 70)
        output_parts.extend(file_results)
    else:
        output_parts.append("All files passed!")

    return '\n'.join(output_parts)


@mcp.tool()
def get_radia_lint_rules() -> str:
    """
    List all available NGSolve lint rules with descriptions.

    Returns a summary of each rule, its severity, and what it checks for.
    """
    rules_info = [
        {
            'rule': 'ngsolve-missing-trace-bem',
            'severity': 'CRITICAL',
            'description': (
                'BEM operators (LaplaceSL/HelmholtzSL) on HDivSurface require '
                '.Trace() on trial/test functions. Without it, boundary-edge '
                'DOFs get corrupted, causing wildly wrong results.'
            ),
            'fix': 'Use j_trial.Trace()*ds(...) instead of j_trial*ds(...).',
        },
        {
            'rule': 'bessel-jv-not-iv',
            'severity': 'CRITICAL',
            'description': (
                'Circular wire SIBC must use modified Bessel functions iv (I0, I1), '
                'NOT regular jv (J0, J1). jv gives correct R_ac/R_dc but wrong sign '
                'on internal inductance Im(Z).'
            ),
            'fix': 'from scipy.special import jv -> from scipy.special import iv',
        },
        {
            'rule': 'hcurl-missing-nograds',
            'severity': 'HIGH',
            'description': (
                'HCurl space for magnetostatics should use nograds=True to '
                'remove gradient null space. Without it, the curl-curl system '
                'is singular.'
            ),
            'fix': 'Add nograds=True: HCurl(mesh, order=2, nograds=True)',
        },
        {
            'rule': 'eddy-current-missing-complex',
            'severity': 'HIGH',
            'description': (
                'HCurl/H1 space in eddy current context without complex=True. '
                'Frequency-domain analysis requires complex-valued FE spaces.'
            ),
            'fix': 'Add complex=True: HCurl(mesh, order=2, complex=True)',
        },
        {
            'rule': 'efie-v-minus-sign',
            'severity': 'HIGH',
            'description': (
                'EFIE system (Zs*M_LL + jw*mu_0*V_LL)*I = -jw*b must use '
                'POSITIVE sign on V_LL term. A minus sign violates Lenz\'s law.'
            ),
            'fix': 'Change minus to plus: Zs*M_LL + jw*mu_0*V_LL (not minus).',
        },
        {
            'rule': 'peec-p-over-jw',
            'severity': 'HIGH',
            'description': (
                'PEEC Loop-Star P/(jw) causes low-frequency breakdown (40-340% error). '
                'Use reformulated Schur complement or stabilized EFIE.'
            ),
            'fix': 'Precompute P^{-1}@M_LS, multiply by jw. Or use mode="stabilized".',
        },
        {
            'rule': 'ngsolve-precond-after-assemble',
            'severity': 'MODERATE',
            'description': (
                'BDDC Preconditioner must be registered BEFORE .Assemble() '
                'to access element matrices.'
            ),
            'fix': 'Move Preconditioner(a, "bddc") BEFORE a.Assemble().',
        },
        {
            'rule': 'ngsolve-overwrite-xyz',
            'severity': 'MODERATE',
            'description': (
                'Loop variable x/y/z overwrites NGSolve coordinate '
                'CoefficientFunction. After the loop, the variable is a '
                'scalar, not a coordinate.'
            ),
            'fix': 'Use different loop variable: "for xi in ..." instead of "for x in ...".',
        },
        {
            'rule': 'ngsolve-vec-assign',
            'severity': 'MODERATE',
            'description': (
                'Direct .vec = assignment creates symbolic expression, '
                'not evaluated result. Must use .vec.data = to evaluate.'
            ),
            'fix': 'Use gfu.vec.data = ... instead of gfu.vec = ...',
        },
        {
            'rule': 'ngsolve-dim2-occ',
            'severity': 'MODERATE',
            'description': (
                '2D OCC geometry (Rectangle, Face) requires dim=2 parameter '
                'in OCCGeometry(). Without it, a 3D surface mesh is generated.'
            ),
            'fix': 'Add dim=2: OCCGeometry(shape, dim=2)',
        },
        {
            'rule': 'ngsolve-cg-on-saddle-point',
            'severity': 'MODERATE',
            'description': (
                'CG solver on A-Omega mixed formulation (saddle-point system). '
                'The system is indefinite, CG may diverge.'
            ),
            'fix': 'Replace solvers.CG() with solvers.GMRes() or MinRes().',
        },
        {
            'rule': 'ngsolve-kelvin-missing-bonus-intorder',
            'severity': 'MODERATE',
            'description': (
                'Kelvin domain integration without bonus_intorder. The varying '
                'Jacobian requires higher quadrature for accurate results.'
            ),
            'fix': 'Add bonus_intorder=4: dx("Kelvin", bonus_intorder=4)',
        },
        {
            'rule': 'ngsolve-vectorh1-for-em',
            'severity': 'MODERATE',
            'description': (
                'VectorH1 used in electromagnetic context. VectorH1 enforces full '
                'C^0 continuity on ALL components, which is wrong for EM fields.'
            ),
            'fix': 'Replace VectorH1 with HCurl (for E, A) or HDiv (for B, J).',
        },
        {
            'rule': 'ngsolve-pinvit-no-projection',
            'severity': 'MODERATE',
            'description': (
                'PINVIT/LOBPCG eigenvalue solver on HCurl without gradient '
                'projection. Curl-curl null space produces spurious zero eigenvalues.'
            ),
            'fix': 'Build gradient projection via fes.CreateGradient().',
        },
        {
            'rule': 'joule-heat-missing-conj',
            'severity': 'MODERATE',
            'description': (
                'Joule heat computed as InnerProduct(E, E) instead of '
                'InnerProduct(E, Conj(E)). Complex E*E != |E|^2.'
            ),
            'fix': 'Use: 0.5 * sigma * InnerProduct(E, Conj(E)).real',
        },
        {
            'rule': 'scattered-eddy-missing-a0',
            'severity': 'HIGH',
            'description': (
                'joule_loss_density() called WITHOUT A0= in a file that sets a '
                'background/applied field. The FE unknown gfA is the SCATTERED '
                'potential; total E = -jw (A0 + gfA + grad(Phi)), so dropping A0 '
                'overestimates the loss ~10x. (Coil/total-field sources where gfA '
                'is already total correctly use A0=None and are not flagged.)'
            ),
            'fix': 'Pass A0=<background>: joule_loss_density(gfA, gfPhi, sigma, omega, A0=A0).',
        },
        {
            'rule': 'ngsolve-set-definedon-in-loop',
            'severity': 'HIGH',
            'description': (
                'GridFunction.Set(..., definedon=...) called inside a loop. Set() '
                'zeros the whole vector first, so looping over boundaries/materials '
                'leaves only the LAST region (prior ones wiped -> e.g. conductor '
                'potentials collapse to 0, energy 0).'
            ),
            'fix': 'Set all regions in ONE call: gf.Set(mesh.BoundaryCF({...}, default=0), definedon=mesh.Boundaries("a|b")).',
        },
        {
            'rule': 'peec-low-nseg',
            'severity': 'MODERATE',
            'description': (
                'Circular coil PEEC with n_seg < 32 may give poor coupling accuracy.'
            ),
            'fix': 'Increase n_seg to 64 or higher.',
        },
        {
            'rule': 'classical-efie-breakdown',
            'severity': 'MODERATE',
            'description': (
                'Classical EFIE using 1/kappa^2 has O(kappa^{-2}) '
                'condition number blow-up at low frequency.'
            ),
            'fix': 'Use stabilized EFIE: [A_k, Q_k; Q_k^T, kappa^2*V_k].',
        },
    ]

    lines = ["NGSolve Lint Rules", "=" * 50, ""]
    for r in rules_info:
        lines.append(f"[{r['severity']}] {r['rule']}")
        lines.append(f"  {r['description']}")
        lines.append(f"  Fix: {r['fix']}")
        lines.append("")

    return '\n'.join(lines)


@mcp.tool()
def taskmanager(topic: str = "overview") -> str:
    """
    NGSolve TaskManager parallelism — usage, MKL interaction, audit, C++.

    ★ AGENTS.md policy: follow the NGSolve-native execution model. Use
    NGSolve TaskManager for Radia thread-level parallelisation (NOT raw
    OpenMP / std::thread / private pools).  External threaded kernels
    such as MKL dense LU are guarded with SuspendTaskManager.  This tool documents the
    `with TaskManager():` pattern, `SetNumThreads()` + `--nthreads`
    CLI convention, MKL nesting trap, and the C++ `ngcore::ParallelFor`
    equivalent.  Includes a 2026-05-27 audit of the IH panel solver
    scripts (`calc_inductance.py` / `calc_fem_kelvin.py` /
    `calc_fem_coilmesh.py` / `calc_heat*.py`) — all PASS.

    Args:
        topic: One of:
            "overview"          - Policy + what TaskManager does (DEFAULT)
            "usage"             - `with TaskManager():` pattern, when to wrap
            "set_num_threads"   - `SetNumThreads()` + `--nthreads` CLI convention
            "mkl_interaction"   - TaskManager vs MKL thread pool nesting
                                  (alias: "mkl", "pardiso")
            "cpp_kernels"       - `ngcore::ParallelFor` in custom C++ TUs
                                  (alias: "cpp", "parallel_for")
            "audit_radia_ih"    - 2026-05-27 IH panel solver audit
                                  (alias: "audit", "radia_ih")
            "common_mistakes"   - Wrap-outside-Assemble, missing
                                  --nthreads, OMP_NUM_THREADS inherited
                                  from Cubit, etc.
                                  (alias: "mistakes", "pitfalls")
            "all"               - Everything concatenated
    """
    return get_taskmanager_knowledge(topic)


@mcp.tool()
def ngsolve_usage(topic: str = "all") -> str:
    """
    Get NGSolve finite element library usage documentation.

    NGSolve is a high-performance FEM library for electromagnetic simulation.
    This tool provides API patterns, best practices, and common pitfalls
    gathered from official tutorials, documentation, and community forums.

    Sources:
      - https://docu.ngsolve.org/latest/i-tutorials/
      - https://forum.ngsolve.org/

    Args:
        topic: Documentation topic. Options:
            "all"              - Complete documentation
            "overview"         - Installation, workflow, direct solvers
            "spaces"           - FE spaces (H1, HCurl, HDiv, HDivSurface, SurfaceL2)
            "maxwell"          - Maxwell/magnetostatics (A-formulation, BDDC, materials)
            "solvers"          - Direct & iterative solver selection guide
            "preconditioners"  - BDDC, multigrid, Jacobi, AMG configuration
            "bem"              - Boundary element method (ngsolve.bem, LaplaceSL, FEM-BEM coupling)
            "ngsolve_bem_50"    - .vol visualization + 50-case NGSolve.BEM/MATLAB comparison lane
            "vol_double_click"  - Windows .vol/.sol double-click handler rules
            "vibroacoustic_drum" - struck-drum FEM/BEM teaching example target
            "curved_vol_geometry" - optional high-order geometry with P1 unknowns
            "mesh"             - Mesh generation (OCC geometry, STEP import, surface mesh)
            "nonlinear"        - Newton's method for nonlinear problems
            "pitfalls"         - Common mistakes and how to avoid them (40 items)
            "linalg"           - Vector/matrix operations, NumPy interop
            "formulations"     - EM formulations: A, Omega, A-Phi, T-Omega, Kelvin (EMPY)
            "adaptive"         - Adaptive mesh refinement with ZZ error estimator (EMPY)
            "darwin"           - Darwin approximation, Surface Impedance BC, Extended Darwin
            "esim"             - ESIM: nonlinear Zs(H,w) Robin BC for any FEM formulation
            "treecotree"       - Tree-cotree splitting, low-freq stability, field-circuit coupling
            "pml"              - Perfectly Matched Layers for open boundary (full-wave)
            "decomposition"    - Domain decomposition: FETI-DP, BDDC, DFDD, AWE/SSP
            "material"         - Material modeling: anisotropy, BH curves, Fixed-Point method
            "ironloss"         - Iron loss estimation: decomposition, FEM computation, steel grades
            "practical"        - Practical techniques: voltage source, force/torque, rotation, coupling
            "team7"            - TEAM Problem 7: eddy current benchmark (A-formulation, OCC geometry, BDDC/AMS solver)
            "multiphysics"     - COMSOL-class couplings: induction heating EM->thermal (joule_loss_density + solve_heat_steady), the scattered-field A0 gotcha
            "cross_validation_registry"
                               - Reusable validation scripts/summary JSONs and the
                                 public-safe MCP knowledge hooks that learned from them
    """
    return get_ngsolve_documentation(topic)


@mcp.tool()
def sparsesolv(topic: str = "all") -> str:
    """
    Get sparsesolv documentation and code examples.

    Since 2026-05-08, sparsesolv ships inside the radia wheel as the
    submodule `radia.sparsesolv_ngsolve` (built from src/ext/sparsesolv/).
    The legacy standalone `ngsolve-sparsesolv` PyPI package is retired.
    Import: from radia.sparsesolv_ngsolve import CompactAMSPreconditioner, COCRSolver, ...

    Source: src/ext/sparsesolv/ in the Radia monorepo.

    See also (theory + decision tree, not code-usage):
        - radia_mcp.matrix_solvers — solver theory + genealogy + which-to-pick
            * matrix_solvers_overview('decision_tree')
            * matrix_solvers_preconditioners('ams_hiptmair_xu') — CompactAMS theory
            * matrix_solvers_krylov('cocg_cocr') — COCR Sogabe-Zhang 2007 origin
            * matrix_solvers_em_specific('eddy_current_stab') — shifted-prec rationale

    Args:
        topic: Documentation topic. Options:
            "all"              - Complete documentation
            "overview"         - Library overview, add-on positioning, features
            "api"              - Python API reference (solvers, preconditioners)
            "examples"         - Usage examples (Poisson, curl-curl, complex, etc.)
            "abmc"             - ABMC ordering: parallel triangular solve optimization
            "compact_ams"      - Compact AMS: theory, benchmarks, COCR solver
            "best_practices"   - Preconditioner selection, complex systems, tips
            "build"            - Build and installation instructions
            "example_poisson"  - Ready-to-run: 2D Poisson with ICCG
            "example_curlcurl" - Ready-to-run: 3D curl-curl with auto-shift IC
            "example_eddy"     - Ready-to-run: Complex eddy current problem
            "example_precond"  - Ready-to-run: IC/SGS with NGSolve CGSolver
            "example_divergence" - Ready-to-run: Divergence detection
            "example_compact_ams" - Ready-to-run: Compact AMS + COCR eddy current
    """
    return get_sparsesolv_documentation(topic)


@mcp.tool()
def esim(topic: str = "all") -> str:
    """
    Get ESIM (Effective Surface Impedance Method) general documentation.

    ESIM extends linear SIBC to nonlinear magnetic materials by solving a 1D
    cell problem through the conductor depth.  Returns a field-dependent
    surface impedance Z_s(H_t) for use in BEM/FEM with surface impedance
    boundary conditions.

    This tool documents the GENERAL technique (cell problem mathematics,
    Karl iteration, module API) without coupling to any specific
    application.  For application-specific use of ESIM (induction heating
    workpieces with steel BH curves), see `mcp-server-ih.ih_sibc(topic="esim")`.

    Args:
        topic: Documentation topic. Options:
            "all"             - Complete documentation
            "overview"        - When ESIM vs linear SIBC; nonlinear conductors
            "cell_problem"    - 1D BVP, BCs, geometries (slab/cylinder/finite_slab)
            "karl_iteration"  - Picard relaxation, convergence pitfalls,
                                per-element vs per-node Z_s
            "module_api"      - radia.esim_cell_problem.ESIMFiniteSlabSolver
                                + BEM-SIBC / FEM-SIBC coupling examples
    """
    return get_esim_documentation(topic)


@mcp.tool()
def urn(topic: str = "all") -> str:
    """
    Universal Relaxation Network (URN): causal/passive rational fitting of a
    complex frequency response, with direct time-domain (relaxation-network /
    SPICE / auxiliary-ODE) synthesis.

    URN decomposes Z(omega) (impedance, dispersive eps/mu, or an open-boundary
    DtN symbol G_n(omega)) into a SPARSE sum of physical relaxation mechanisms
    (Debye, Cole-Cole, Cole-Davidson, Havriliak-Negami, CPE, Warburg, Gerischer,
    RLC, skin-effect; series + admittance branch) with KAN-style adaptive tau and
    attention.  Every basis is passive, so the fit is causal/passive BY
    CONSTRUCTION and maps to one first-order auxiliary ODE per pole (fractional
    terms -> short RC/RL ladder) -- the local-in-time operator an FETD /
    Newmark-beta solver needs.  Beats Vector Fitting on fractional/Cole-Cole data
    (avg ~22.8% lower NRMSE on NASA battery + TDK ferrite datasets).

    Use to turn a frequency-domain absorbing-BC / dispersive-layer response into
    a stable broadband time-domain model.  For transient FEM/BEM or Maxwell
    solvers, topic="cq" explains the URN H(s) -> convolution-quadrature bridge.
    Run the fit with the urn_fit tool.
    Ref: Sugahara & Sato, IEEE Access 2026; impl docs/universal_relaxation_network.

    Args:
        topic: all | overview | method | api | timedomain | cq | application
    """
    return get_urn_documentation(topic)


@mcp.tool()
def urn_fit(data_csv: str, freq_col: int = 0, real_col: int = 1,
            imag_col: int = 2, delimiter: str = ",", skip_rows: int = 0,
            n_debye: int = 3, n_cole_cole: int = 2, n_warburg: int = 1,
            n_cole_davidson: int = 0, sparsity_weight: float = 0.01,
            n_epochs: int = 2000, n_restarts: int = 3, spice_out: str = "") -> str:
    """
    Fit a complex frequency response with a Universal Relaxation Network and
    return the discovered relaxation mechanisms, the fit NRMSE, and a SPICE
    netlist (== the auxiliary-ODE ladder for a time-domain / Newmark-beta solver).

    Input is a CSV with columns (frequency_Hz, Re(Z), Im(Z)).  Requires torch;
    training is iterative -- lower n_epochs / n_restarts for a faster, rougher
    fit (defaults ~2000/3 are a responsive compromise; the paper uses 6000/10).

    Args:
        data_csv: path to a CSV of the frequency response.
        freq_col/real_col/imag_col: 0-based column indices.
        delimiter, skip_rows: CSV parsing.
        n_debye/n_cole_cole/n_warburg/n_cole_davidson: basis counts to try.
        sparsity_weight: L1-ish penalty that prunes unused bases.
        n_epochs/n_restarts: Adam iterations / random restarts.
        spice_out: if set, also write the SPICE netlist to this path.
    """
    return urn_fit_from_csv(
        data_csv, freq_col, real_col, imag_col, delimiter, skip_rows,
        n_debye, n_cole_cole, n_warburg, n_cole_davidson, sparsity_weight,
        n_epochs, n_restarts, spice_out)


@mcp.tool()
def cln_3d(topic: str = "all") -> str:
    """
    Get 3D Cauer Ladder Network (CLN) / Kameari-Tanimoto iteration
    documentation for eddy current analysis with NGSolve.

    Captures Tanimoto-Kameari iterative methods from the master's thesis
    + production code (public-safe curated corpus, ~25 notebooks):
      - A-T formulation (primary)
      - T-Ω formulation (H1 confined to conductor)
      - A-Φ formulation (HCurl + H1 mixed)
      - Constraint variants: penalty stabilization, explicit Coulomb gauge
      - Production solvers: SparseSolvPy ICCG, accICCG, NGSolve CG, direct

    Each formulation produces a Cauer-II ladder {R_n, L_n} via Schmidt
    orthogonalization on impressed J source. Validated against
    cylindrical TM-mode analytical R/L for n=0..9.

    Open research: Kameari + Kelvin transformation combination (3D
    HCurl A-formulation hits ~25× discrepancy with mpmath BEM Foster
    target due to A_ext gauge being unbounded at infinity; future
    work includes T-Ω with reduced-Ω = -H_0·z + Ω_r).

    Args:
        topic: Documentation topic. Options:
            "all"           - Complete documentation
            "overview"      - Mathematical foundation, three formulations
            "notebooks"     - Index of 修論 / 定式_誤差検証 / 静止器回転機用
            "formulas"      - Cauer-II synthesis, drift diagnostic,
                              bonus_intorder critical setting
    """
    full_doc = get_cln_3d_documentation()
    if topic == "all":
        return full_doc
    sections = {
        "overview": "OVERVIEW",
        "notebooks": "NOTEBOOK_INDEX",
        "formulas": "KEY_FORMULAS",
    }
    if topic in sections:
        # Split on H2 markdown headers as section breaks
        from .knowledge import cln_3d
        if topic == "overview":
            return cln_3d.CLN_3D_OVERVIEW
        if topic == "notebooks":
            return cln_3d.CLN_3D_NOTEBOOK_INDEX
        if topic == "formulas":
            return cln_3d.CLN_3D_KEY_FORMULAS
    return full_doc


@mcp.tool()
def bem_cln(topic: str = "all") -> str:
    """
    Get BEM-CLN (per-element multipole CLN with Schur-F termination)
    documentation: multi-conductor extension of single-conductor
    Schur-F CLN, using polarizability alpha(s) and integral-equation
    Green's function coupling.

    Backs Sugahara, Nagamine, Hane (2026) IEEE Trans Mag submission,
    sections V.G (DOF accounting) and V.H (verification).

    Key features:
      - polarizability alpha(s) = V - Y_cln(s) / sigma (DC = 0, PEC = V built in)
      - 2D coupling: 1/D^2, 3D coupling: mu_0 / (4 pi D^3)
      - bounded alpha -> no phenomenological saturation factor needed
      - per-element DOF = N_Cauer + 1; total = N (N_Cauer + 1)

    Args:
        topic: Documentation section. Options:
            "all"           - Complete documentation
            "overview"      - Framework summary, DOF accounting
            "2d_rigorous"   - Phase 2.5 canonical 2D rigorous
            "3d"            - Phase 3 B rigorous 3D cuboid extension
            "scripts"       - Index of Mathematica verification scripts
    """
    if topic == "all":
        return get_bem_cln_documentation()
    from .knowledge import bem_cln as bcln
    if topic == "overview":
        return bcln.BEM_CLN_OVERVIEW
    if topic == "2d_rigorous":
        return bcln.BEM_CLN_2D_RIGOROUS
    if topic == "3d":
        return bcln.BEM_CLN_3D
    if topic == "scripts":
        return bcln.BEM_CLN_NOTEBOOK_INDEX
    return get_bem_cln_documentation()


@mcp.tool()
def cln_sibc_orthogonal(section: str = "all") -> str:
    """
    Get CLN expansion-point + SIBC orthogonal-residual theory documentation.

    New theory established 2026-05-16 in IGTE 2026 Sugahara work:
    a small number N of Cauer ladder stages plus the SIBC analytical
    asymptote span the full eddy-current frequency response via
    Foster eigenmode L^2-orthogonality:

        L^2(Ω) = span{φ_0, ..., φ_{N-1}} ⊕ span{φ_N, φ_{N+1}, ...}
        Y_exact(jω) = Y_CLN^N(jω) + Y_SIBC^⊥(jω)
        Y_SIBC^⊥(jω) -> K_SIBC (jω)^{-1/2}  as  ω -> ∞

    K_SIBC = √(σ/μ) × geometric factor (slab=1, 2D prism=perimeter,
    3D body=surface integral with edge/corner corrections).  This
    avoids the precision frontier of high-stage Cauer extraction
    (Nagamine et al. 2026 verified interval arithmetic shows
    ~60×/stage interval growth even with 192-bit MPFR + affine).

    Verified analytically with Mathematica on:
      - 1D slab (c = 1 mm Cu)
      - 2D rectangular prism (17.72 × 2 mm, Dirichlet A_z = 0)

    Open: 2D square (Nagamine geometry, degeneracy bundling) and
    3D cuboid (HDivDivFreeHex + Green's theorem on surface integrals).

    Companion to Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune,
    Matsuo (JJIAM 2026 submitted) verified 2D square prism Cauer
    extraction — Nagamine = high-N verified rungs, this theory =
    small-N + SIBC asymptote.

    Args:
        section: Documentation section to return:
            "list"           - list available sections
            "overview"       - one-page summary of the construction
            "matsuo"         - relation to Matsuo SA-26-014 expansion point
            "kuriyama"       - Kuriyama 2019 multi-expansion canonical method
            "math"           - derivation, orthogonality, asymptote
            "verification"   - Mathematica results (1D slab + 2D rectangle)
            "nagamine"       - link to Nagamine 2026 verified extraction
            "outlook"        - 2D Nagamine square + 3D cuboid completion path
            "xfem_vs_sibc"   - XFEM (Hiruma 2023) / classical SIBC /
                               augmented CLN decision framework, with
                               port-driven scope and stacking strategy
                               for volume-source problems
            "all"            - full documentation (default)

    Returns:
        Markdown text of the requested section.
    """
    if section == "all":
        return get_cln_sibc_orthogonal_documentation()
    return get_cln_sibc_orthogonal_section(section)


@mcp.tool()
def cln_3d_notebook(name: str = "list") -> str:
    """
    Retrieve Tanimoto's raw 3D CLN notebook Python code.

    Provides direct access to the canonical Python code from Tanimoto's
    master's thesis + production notebooks at public-safe curated corpus
    Use this when you need to see the actual implementation details
    (HCurl space construction, Kameari iteration loop, ICCG solver
    invocation, output extraction, etc.).

    Args:
        name: Notebook identifier. Options:
            "list"       - List available notebooks with file sizes
            "AT"         - A-T formulation (primary 修論 reference,
                           cylinder, 10-stage Kameari, SparseSolvPy ICCG)
            "T_Omega"    - T-Ω formulation (HCurl × H1 with Ω confined
                           to conductor)
            "APhi"       - A-Φ formulation (HCurl + H1, body current
                           via σ∇Φ)
            "2D"         - 2D scalar reference (pedagogical, Kameari
                           formula validation)
            "production" - 2024-09-17 production: A + ICCG with inline
                           gauge correction, accICCG params, type1 HCurl

    Returns:
        Full Python script content (~3-9 KB each), or list of available.
    """
    return get_cln_3d_notebook(name)


@mcp.tool()
def cln_sphere_dd_pipeline() -> str:
    """
    Get the Sphere DD (double-double, ~32 digit) VIM Cauer Ladder Network
    extraction pipeline reference.

    Verified-arithmetic CLN extractor demonstrated on the canonical Cu
    sphere benchmark (R=10mm, sigma=5.8e7 S/m, uniform B_z=1T). Pure
    Python/CuPy implementation of the entire VIM-CLN chain in DD
    precision: kernel evaluation (mpmath elliptic K(m), E(m) at 40 digit
    via dd_axisym_kernel.py), DD K/M/b assembly with multiprocessing
    (dd_sphere_axisym_mp.py), mpmath Cholesky-based generalized eigh
    at dps=35-50, and verified-interval Hankel-Pade Cauer extraction
    (mpmath.iv at 80 digit).

    Reaches DD precision floor K_lo/K_hi = 1.1e-16 across all matrix
    entries (the entire 14 trailing decimal digits below FP64 are
    correctly captured), with verified-interval relative width < 1e-30
    on all extracted Cauer rungs.

    Returns:
        Markdown documentation covering algorithm steps, sigma-rescaling
        for canonical normalization, multiprocessing speedup, file
        inventory, production scaling estimates, key implementation
        lessons, and cross-references to Nagamine, Stoll, Hiruma,
        Sugahara TEAM 28.
    """
    return get_cln_sphere_dd_documentation()


@mcp.tool()
def hdiv_vim(topic: str = "overview") -> str:
    """
    HDiv-type VIM (Volume Integral Method) demag operator -- the lab's FEEC H(div) RT
    soft-iron demag route.  Canonical reference:
    docs/hdiv_vim/README.md.

    Key idea: SYMMETRIC demag operator N = B^T G B with the loop modes FIELD-NULL BY
    CONSTRUCTION (loops = ker B) -> mu_r-INDEPENDENT convergence + NO hand-crafted loop-star.
    VALIDATED: linear demag (sphere/spheroid/triaxial EXACT vs analytic), NONLINEAR (damped
    Newton; cube & C-yoke <1-3% vs shipped Radia in ~6 iters), distorted-mesh mu_r-independence,
    CURVED + high-order (demag exact; field accuracy-per-DOF ~10-30x vs flat Radia), and SYMMETRY
    models 1/2,1/4,1/8 (loops automatic + image-method demag).  Production entry
    radia.vim.Solve; the C++ _ChargeGramHMatrix kernel is the SOLE demag operator (the
    dense Python Gram path + the analytic_gram/wilton_surface kwargs were REMOVED 2026-06-23).
    TaskManager is assumed: shared HACApK build paths and long C++ solve loops stand up or reuse
    NGSolve RegionTaskManager; direct diagnostic `.matvec()` calls plus Python/NGSolve assembly
    follow caller-wraps `with ng.TaskManager():`.  The C++ HDiv CG/MINRES/Picard kernels use
    ParallelFor/ParallelForRange for charge gather, dot products, preconditioner/vector updates,
    and AtomicAdd for sparse face-vector scatters.

    Args:
        topic: One of:
            "overview"       - what it is + why (symmetric, loops field-null, mu_r-independent) [DEFAULT]
            "implementation" - C++/pybind/Python files + APIs (rad_hdiv_vim, _ChargeGramHMatrix, ...)
            "scaling"        - the C++ charge-Gram H-matrix (exact analytic near AND far)
            "verification"   - golden tests (tests/feec/) + the verify-first bug catches
            "nonlinear"      - damped Newton on the exact C++ charge Gram; C-yoke vs Radia;
                               fail-loud on non-convergence; the honest reference distinctions
            "curved"         - curved + high-order demag (ngsolve.bem single-layer): sphere/spheroid/
                               triaxial exact; accuracy-per-DOF ~10-30x vs flat Radia; curved x nonlinear
            "symmetry"       - 1/2, 1/4, 1/8 models: loops automatic (ker B) + image-method demag value
            "cross_method"   - the demag tensor cross-validated by 3 independent discretizations
                               (FEEC surface-charge VIM == volume-FE A-formulation == BEM == Osborn)
            "reference_audit"- HOW TO DEBUG a disagreement with a FEM cross-validation reference:
                               the audit ladder (solution scalar -> evaluator invariance ->
                               drive-equivalence probe -> split tests -> closed-form arbiter) +
                               the reference-side trap catalog (coil-polygon current deficit,
                               frozen-edge truncation, conjugate-potential sign/branch cuts, ...)
            "status"         - done / open summary
            "all"            - everything
    """
    return get_hdiv_vim_documentation(topic)


@mcp.tool()
def kelvin_identify_post_hoc(topic: str = "all", vol_path: str = "") -> str:
    """
    Add Kelvin Periodic Identifications to an existing NGSolve mesh
    AFTER load (no Cubit launcher / OCC Identify needed).

    Use this when you have a `.vol` file with `kelvin_int` /
    `kelvin_ext` boundary labels but no `Identifications` section yet
    -- typical when the Cubit C++ exporter's all-or-nothing vertex
    matching skipped due to a single tolerance miss (1/8-octant
    geometry), or when the mesh came from outside the Cubit panel
    pipeline.

    Public Python API (call this from user code, NOT from the AI):
        from radia import add_kelvin_identification
        info = add_kelvin_identification(mesh)
        # -> dict(n_pairs, n_unmatched, kelvin_offset, max_dist, ...)

    Args:
        topic: Documentation topic. Options:
            "all"           - Full documentation (overview + API + workflow + caveats)
            "overview"      - Why this exists, when to use, the 3-line solution
            "api"           - add_kelvin_identification() signature + return dict
            "workflow"      - 3-step usage: detect_offset -> add -> verify
            "snippet"       - Copy-paste minimal example code
            "caveats"       - 1:1 mesh requirement, tolerance, idnr handling
            "verify"        - If vol_path given: load .vol and run helper
                              (returns the helper's diagnostic dict as text)
        vol_path: optional path to a .vol file to actually run the helper
                  against (only used when topic="verify").

    Companion tool: `kelvin_transformation` provides the underlying
    theory (Kelvin inversion, material scaling, formulations).  This
    tool is specifically about the post-hoc identification helper.
    """
    if topic == "verify":
        if not vol_path:
            return ("topic='verify' requires vol_path argument.  "
                    "Pass the path to a .vol file with kelvin_int / "
                    "kelvin_ext labels.")
        import json
        import os
        if not os.path.exists(vol_path):
            return f"vol_path not found: {vol_path}"
        try:
            from ngsolve import Mesh
            from radia import (add_kelvin_identification,
                                detect_kelvin_offset,
                                has_kelvin_identification)
        except ImportError as e:
            return (f"Cannot import ngsolve / radia in this process: {e}.  "
                    f"This tool needs radia + ngsolve installed in the MCP "
                    f"server's Python environment.")
        mesh = Mesh(vol_path)
        report = {
            "vol_path": vol_path,
            "materials": list(mesh.GetMaterials()),
            "boundaries": list(mesh.GetBoundaries()),
            "has_existing_identifications": has_kelvin_identification(mesh),
            "detected_offset": detect_kelvin_offset(mesh),
        }
        try:
            info = add_kelvin_identification(mesh)
            report["add_kelvin_identification"] = info
            report["status"] = "ok"
        except ValueError as e:
            report["status"] = "error"
            report["error"] = str(e)
        return json.dumps(report, indent=2, default=str)
    return _get_kelvin_identify_post_hoc_doc(topic)


@mcp.tool()
def kelvin_transformation(topic: str = "all") -> str:
    """
    Get Kelvin transformation documentation for open boundary FEM problems.

    The Kelvin transformation maps an unbounded exterior domain to a bounded
    computational domain, enabling FEM solutions without artificial truncation.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - Mathematical foundation and key principles
            "h_formulation"  - H-field perturbation potential formulation
            "a_formulation"  - Vector potential formulation (coils)
            "3d"             - 3D sphere/solid examples (H1)
            "hcurl_3d"       - 3D HCurl A-formulation (calc_fem_kelvin.py)
            "verified_recipe" - Verified A-formulation + Periodic Kelvin
                                + BDDC cookbook (5 mandatory elements,
                                failure-mode -> root-cause table).
                                **Read this FIRST when debugging an
                                HCurl eddy-current script that is off
                                by O(10x) from analytical reference.**
            "adaptive"       - Adaptive mesh refinement with Kelvin
            "mesh_control"   - WHERE to spend elements (consolidated): the
                               Gamma-conforming constraint + six MEASURED
                               pillars (exterior volume is free, floor=Curve
                               order, p>=n & p-vs-h, optimal R/a~3, corner hp,
                               DoF-cost 1/45; holds in 2D on the -n/R ladder)
                               + the CLOSED-FORM mesh-adequacy criterion
                               (source order p*=ceil(ln eps/ln(d_max/R)),
                               geometry floor (h/R)^2k, eccentric/multi-body,
                               the A-form centre, the apparatus design calc)
                               + the non-separable build + DtN->CLN arc
                               (square C4v / cube O_h, 2D conformal disk)
            "identify"       - Periodic boundary Identify() best practices
            "tips"           - Common mistakes and performance tips
            "robustness"     - Robustness checklist (mesh copy, material scaling,
                               FreeDofs verification, symmetry models, GND)
            "verification"   - Numerical verification (single-domain approach)
            "periodic_wedge" - 1/n sector (symmetry model) with Periodic BC
            "material_exterior" - Provenance of material/conducting exterior
                               Kelvin: Freeman-Lowther (air exterior) vs
                               transformation-optics sigma (Ward-Pendry) vs
                               Sugahara 2022's validated sigma-conformal ECT
                               fusion (conductor crosses the truncation); the
                               formulation basis for radia.open_boundary.
                               kelvin_dtn's (a/r)^4 sigma / (a/r)^2 mu weights.
                               Read before claiming "material Kelvin" novel.
    """
    return get_kelvin_documentation(topic)


@mcp.tool()
def axifem_documentation(topic: str = "all") -> str:
    """
    Get radia-core axifem documentation: Henrotte axisymmetric FE
    add-on for NGSolve (registered FESpace name: "axihenrotte").

    Use this when designing or reviewing axisymmetric eddy-current /
    magnetostatic FEM problems with axis-touching elements, or when
    comparing axihenrotte to standard NGSolve H1 elements.

    Canonical API: ``FESpace("axihenrotte", mesh, order=k)`` or
    ``H1Henrotte(mesh, order=k)``.  order=1 dispatches to P1 triangles
    or Q1 quads; order=2 dispatches to P2 triangles or Q2 quads.  P2
    triangles are curved-mesh aware after ``mesh.Curve(2)``.  Curved Q2
    quads are production opt-in via ``H1Henrotte(mesh, order=2,
    curvedquad=True)``.

    Args:
        topic: Documentation section. Options:
            "all"             - all sections concatenated
            "overview"        - what it is and why NGSolve doesn't already have it
            "support_matrix"  - exact P1/P2/P2 curved/Q1/Q2/Q2 curved status
            "api"             - FESpace("axihenrotte", mesh, order=k) usage
            "hodge_geometry"  - differential-geometry/Hodge view of Henrotte
                                axisymmetric reduction
            "taskmanager"     - NGSolve-native parallel execution contract
            "basis_p1"        - order=1 P1 triangle + Q1 quad basis details
            "basis_p2"        - order=2 P2 triangle + Q2 quad basis details
            "curved_geometry" - P2 curved triangle and opt-in Q2 curved quad support
            "vs_standard_h1"  - 6-property comparison table vs H1 order=2
            "validation"      - cross-validation references; Hessian-of-W convention
            "kelvin"          - Phase B3 z-offset Kelvin recipe (Periodic + H1Henrotte,
                                sphere -0.001 % vs Stoll, mu-vs-nu factor & Curve(2) gotchas)
            "magnet"          - permanent-magnet source term (FEMM prob3big.cpp port):
                                weak-form RHS + magnetized-sphere validation (-0.05%)
            "file_layout"     - where each piece lives (C++, Mathematica, tests)
            "why_dropped_p3"  - why p=3 was attempted and reverted (Vandermonde cond ~ 1e30)
    """
    return get_axifem_documentation(topic)


@mcp.tool()
def femm_parity_documentation(topic: str = "all") -> str:
    """
    Get FEMM-parity documentation: which FEMM (Finite Element Method Magnetics,
    D. Meeker) analyses are reproduced as EXECUTABLE + TESTED NGSolve capability
    in radia-ngsolve, with the function to call and the analytical benchmark each
    was validated against (all <2%, most <0.2%).

    Use this when designing an FEM analysis that a FEMM user would run in FEMM 4.2,
    to find the equivalent radia-ngsolve function and its API.

    Args:
        topic: Documentation section. Options:
            "all"        - all sections concatenated
            "overview"   - design rule: build capability not a number
            "matrix"     - planar / axisymmetric capability matrix (21 analyses)
            "magnetics"  - magnetostatic / nonlinear / eddy / circuit + axi API
            "scalar"     - electrostatic / heat / current-flow (csolv/hsolv) API
            "lamination" - laminated steel (anisotropic + complex-mu), multi-
                           conductor proximity circuit, FEMM open-bdry cross-check
            "validation" - regression test list and per-test error bounds
    """
    return get_femm_parity_documentation(topic)


@mcp.tool()
def radia_usage(topic: str = "all") -> str:
    """
    Get Radia C++ library usage documentation.

    Covers: field computation, fixed-magnet materials, HDiv-backed soft iron,
    hysteresis, background fields, solver configuration (LU/BiCGSTAB/HACApK),
    NGSolve integration (RadiaField CF), memory management, IMA.

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - Architecture, Radia fields vs HDiv-VIM vs BEM
            "elements"       - ObjRecMag, ObjHexahedron, ObjTetrahedron, ObjWedge, ObjPyramid
            "materials"      - MatLin, MatSatIsoTab, hysteresis, permanent magnets
            "solver"         - rad.Solve, SolverConfig, LU/BiCGSTAB/HACApK
            "field"          - rad.Fld, batch evaluation, A field
            "removed_apis"   - Removed APIs reference (FldEnr/ObjDivMag/UtiDmp/...)
            "ngsolve"        - RadiaField CF, netgen_mesh_to_radia
            "background"     - ObjBckg, Biot-Savart source
            "ima"            - Image Method of Analysis
            "best_practices" - lab best practices, incl. #11 storing
                               numerical-experiment results as JSON
                               (reproducible) + delete-by-folder cleanup
    """
    return get_radia_documentation(topic)


@mcp.tool()
def md2html_usage() -> str:
    """Get md2html converter documentation (MathJax, reference links, styled HTML)."""
    return get_md2html_documentation()


@mcp.tool()
def analytical_formulas(topic: str = "all") -> str:
    """
    Get documentation for radia.analytical_formulas (closed-form reference layer).

    The package collects nine modules of closed-form expressions taken from
    the Wakao-Igarashi-Fujiwara-Kameari review series (IEE Japan, 2002-2004).
    These are the trusted-baseline analytical results that the rest of Radia
    (HDiv-VIM, fixed-magnet fields, PEEC, FEM panels, ngsolve.bem) is
    sanity-checked against.

    For any new analysis the FIRST QUESTION to ask is "is there a closed form
    here that I can validate against?". Use the validation_use_cases topic
    for the practical mapping from analysis type -> applicable formula.

    Args:
        topic: Documentation topic. Options:
            "all"                  - Complete documentation
            "overview"             - Module index, conventions, quick API surface
            "ellipsoid"            - Rotational ellipsoid demag factor + torque
                                     (Part 5 §5, eq 38-44; Osborn references included)
            "ac_locus"             - B_max / B_min of an AC vector phasor's locus
                                     (Part 5 §4, eq 29-37)
            "shielding"            - Magnetic-shell shielding factor S
                                     (Part 1 §5, eq 23-24)
            "rect_magnet_2d"       - 2D rectangular bar A_z, B_x, B_y
                                     (Part 2 §2, eq 2-3)
            "plate_eddy"           - Thin rectangular-plate eddy current
                                     (Part 1 §6.1, eq 26-27)
            "solenoid_central"     - Fabri form factor + closed-form axial field
                                     (Part 4 §4, eq 26-27)
            "three_phase_line"     - Triangle / planar / helical line field
                                     (Part 4 §5, Part 5 §3)
            "elliptic_integrals"   - K(k), E(k) Hastings polynomial approximation
                                     (Part 3 §3, Tables 1-2)
            "gauss_legendre"       - Gauss-Legendre nodes / weights to n=24
                                     (Part 3 §4, Table 3)
            "validation_use_cases" - "I have X, sanity-check it with Y" mapping

    Sources:
        - src/radia/analytical_formulas/   (Python modules, ~40 functions across Part 1-9)
        - tests/analytical_formulas/       (170+ pytest tests, < 1 s)
        - docs/analytical_formulas/        (analytical_formulas.ipynb: 12 demos + PNGs)
        - docs/analytical_formulas.md      (PDF -> code cross-reference)
        - PDFs themselves: lab-internal, not redistributed with the repo.
    """
    return get_analytical_formulas_documentation(topic)


@mcp.tool()
def acoustic_fembem_cross_learnings() -> str:
    """
    Method-selection and validation cross-learnings for radia-acoustic, distilled
    from the readable acoustic FEM/BEM teaching lane (Gypsilab-style solvers +
    their goldens). Tool-agnostic physics only -- no proprietary content.

    This is the acoustic companion to `analytical_formulas`: where that tool gives
    closed-form magnetostatic baselines, this one gives the *how to choose and
    validate* rules for exterior acoustic radiation / scattering, plus the closed-
    form acoustic references in `radia_mcp.radia_ngsolve.acoustics` (pulsating
    sphere, baffled piston, spherical / planar Helmholtz DtN, impedance-reflection).

    Read this when:
      * Deciding whether a scatterer needs a FEM interior or is pure BEM
        (rule: FEM interior ONLY when it carries physics BEM cannot -- elastic /
        inhomogeneous / lossy interior; rigid / soft = surface BC = pure BEM).
      * Standing up a time-domain acoustic solver -- the compact convolution-
        quadrature (Lubich CQ) lane = frequency-domain single-layer operator x a
        thin generic FFT wrapper, A-stability inherited from the Laplace-domain
        kernel, anchored by the imaginary-axis golden `V(-i c k) == Helmholtz
        single layer + coincident-node Delta`.
      * Validating a scatterer with NO analytic reference -- use the shape-
        independent invariants: far-field reciprocity `f(x_hat; d) == f(-d; -x_hat)`,
        Sommerfeld `1/r` radiation decay, radiation-force control-radius
        independence.
      * Coupling FEM pressure to an exterior radiation condition -- the spherical
        Helmholtz DtN eigenvalue is the exact "Kelvin operator on the sphere"
        fast exterior; the lab wave-boundary policy stays high-order Zs, no PML.
      * Asked whether the hodograph transform helps 3D acoustic BEM (it does not --
        acoustic BEM is already linear, 3D does not auto-linearise, and Helmholtz
        breaks the conformal / Kelvin structure; the tractable dual is the source-
        side surface-density linear inverse).

    Sources:
        - packages/radia-mcp/src/radia_mcp/radia_ngsolve/acoustics.py
          (closed-form helpers + this cross-learning document)
        - tests/test_acoustics.py (pytest coverage of the helpers + this document)
        - the Gypsilab `+acoustic_fembem` MCP `convolution_quadrature` knowledge
          topic (public teaching-lane counterpart; no proprietary content).
    """
    return _acoustic_fembem_cross_learnings()


@mcp.tool()
def matlab_acoustic_fembem_agent_guide() -> str:
    """
    Agent guide for MATLAB acoustic FEM-BEM / Gypsilab-style workflows.

    Use this when an agent needs to work on MATLAB acoustic FEM/BEM code while
    respecting the repository policy that MathWorks' official MATLAB MCP Server
    is the execution substrate.  The guide records the stack contract
    (toolbox detection, Code Analyzer, MATLAB unit tests, and mdx fallback),
    plus lab-specific domain gates for `.vol` meshes, P1 FEM/BEM assembly,
    convolution quadrature conventions, and NGSolve / `ngsolve.bem`
    cross-validation.
    """
    return _matlab_acoustic_fembem_agent_guide()


@mcp.tool()
def force_validation(topic: str = "all") -> str:
    """
    EM force extraction in NGSolve + independent <-> NGSolve cross-validation.

    Records how the radia-ngsolve FEM path chooses and computes electromagnetic
    force/torque: weighted Maxwell stress, Maxwell surface stress, Lorentz
    conductor force, air-gap pressure, energy/coenergy checks, dq torque, and
    the cross-validation results that make the NGSolve magnetostatic path
    trustworthy. Reference values are kept as stored regression references.

    Validated (linear magnetostatics, A-form, HCurl order 2):
      * uniformly magnetized sphere: reference == NGSolve to 0.11 %, both <0.5 %
        of the analytic (2/3)Br interior field.
      * coil + linear-iron sphere force: reference == NGSolve to ~3 %, both near
        the dipole-in-gradient analytic.

    Read this when: implementing/validating an EM force computation, deciding
    whether to trust an NGSolve magnetostatic result, or reviewing the
    reference-solve note (the reference_note topic summarises the shared
    geometry / material / source / force-method contract of the cross-check).

    Args:
        topic: all (default) | method_map | eggshell | cross_validation | reference_note
    """
    return get_force_validation_documentation(topic)


@mcp.tool()
def ngsbem_inductance(topic: str = "all") -> str:
    """
    Get ngsolve.bem boundary element method documentation for inductance extraction.

    ngsolve.bem is NGSolve's native boundary element module. Combined with Cubit
    mesh export and SetGeomInfo, it enables accurate inductance extraction
    on high-order curved surface elements.

    Key workflow:
      Cubit mesh -> export netgen "mesh.vol" order N -> LaplaceSL BEM -> L extraction

    Sources:
      - https://docu.ngsolve.org/latest/how_to/ngsbem.html
      - https://github.com/Weggler/docu-ngsbem/ (stabilized BEM)
      - https://github.com/ksugahar/Radia (cubit_mesh_export / export plugin)

    Args:
        topic: Documentation topic. Options:
            "all"            - Complete documentation
            "overview"       - What is ngsolve.bem, comparison with Radia PEEC
            "api"            - LaplaceSL operator, HDivSurface, matrix extraction
            "cubit_workflow" - Cubit -> SetGeomInfo -> Curve -> BEM pipeline
            "curve_order"    - Curve order convergence study for BEM accuracy
            "stabilized"     - Weggler's stabilized BEM for low-frequency
            "examples"       - Runnable examples (circular loop, Cubit torus)
            "best_practices" - Common pitfalls, validation, performance tips
            "known_limitations" - curvaturesafety, TaskManager, QUAD hang, grad(G) gap
    """
    return get_ngsbem_inductance_documentation(topic)


@mcp.tool()
def peec_inductance(topic: str = "all") -> str:
    """
    Get documentation for the Radia PEEC-inductance (coil only, STEP) panel mode.

    Lightest mode in the Radia panel family: STEP / Cubit .jou → perimeter
    filaments → Biot-Savart + Loop-bundle PEEC solve → L_coil, R_coil.
    No workpiece, no BEM, no FEM mesh.  Use when the question is "what
    is this coil's L / R at this frequency?" and you do not need
    heating or workpiece physics.

    Centerline extraction is **STEP-only** (v4.48.1, 2026-05-15):
    no caller-supplied JSON, no ``--path-points-m`` flag.
    Classification-based single dispatch picks ONE of 5 predicates;
    if none cover the conductor bbox, raises with a HINT pointing at
    CAD regeneration or BEM-A switch (NEVER a silent fallback).

    Five centerline predicates (positive-match, no try/except cascade):
      1. Loft of profiles (>= 5 planar end-caps)  -> NN-chain centroids
      2. United multi-turn pancake (CIRCLE-edge stations) -> NN-chain
         circle centres
      3. Single-loop revolution sweep (TORUS / CYLINDER / CONE /
         REVOLUTION + PLANE caps) -> analytical arc
      4. OPEN coil with caps -> longest open lateral rim edge
         (handles "arc + leads" e.g. keiko outsideline.step)
      5. CLOSED full revolution (no caps) -> coil_topology spine

    Sibling-.jou auto-preference: if the user picks ``foo.step`` and
    ``foo.jou`` (case-insensitive stem match) coexists in the same
    directory AND contains the PEEC explicit-centerline pattern,
    the calc switches to the .jou parser automatically.

    Unicode / Japanese path: verified end-to-end (Python + Cubit
    plugin via utf8_path.hpp + PySide6 QProcess).

    Args:
        topic: Options:
            "all"             - Complete documentation
            "overview"        - What it solves, when to use, perimeter placement
            "centerline"      - STEP centerline extraction (5 classification predicates)
            "filament_dispatch" - n_peri filament placement paths (single-dispatch)
            "step_authoring"  - Cubit / build123d recipes for auto-detect-friendly STEPs
            "jou"             - .jou explicit centerline parser
            "sibling_jou"     - Auto-prefer sibling .jou when co-located
            "japanese_path"   - Unicode / Japanese path support
    """
    return get_peec_inductance_documentation(topic)


@mcp.tool()
def fem_bem_schur(topic: str = "all") -> str:
    """
    Get FEM-BEM Schur coupling documentation -- exact open boundary for interior FEM.

    The exterior Laplace DtN assembled from ngsolve.bem replaces the unbounded
    exterior domain with a boundary operator on the coupling surface, giving an
    EXACT transparent BC without mesh truncation, PML, or Kelvin transformation.

    Coupled system:  (K_FEM + P^T Λ_ext P) u = f
    where Λ_ext = V⁻¹ (−½M + K)  [exterior BIE, NGSolve sign convention]
    and P is the H1 volume → SurfaceL2 boundary projection.

    Verification: spherical shell with inner Dirichlet u=cosθ, BEM DtN at outer
    sphere; exact solution u = R_inner² cosθ / r²; L2 rel_err < 2% at order=2.

    Companion code: radia_ngsolve.fem_bem_coupling
      - laplace_fem_bem_schur():  coupled solve (dense, for verification)
      - sphere_shell_fem_bem():   verification on spherical shell

    Args:
        topic: Documentation topic. Options:
            "all"          - Complete documentation
            "overview"     - Concept, DtN operator, coupling matrix
            "api"          - laplace_fem_bem_schur and sphere_shell_fem_bem usage
            "applications" - Open BC magnetostatics, motor models, AGE + BEM
    """
    return get_fem_bem_schur_documentation(topic)


@mcp.tool()
def airgap_motor_workflow(topic: str = "all") -> str:
    """
    Get AGE rotating machine workflow documentation -- nonlinear iron + AGE coupling.

    Wraps the linear AGE core (airgap_machine.py) in a Picard fixed-point iteration
    for machines with saturating iron B-H characteristics.  The air-gap element stays
    analytic (no gap mesh); only the iron FE stiffness is updated each iteration.

    Key validation:
      - nu_cf_fn=None (or constant nu) → 1 iteration, identical to direct airgap_solve
      - Froelich saturation model → converges in 5-15 iterations
      - age_motor_rotation_sweep → torque-angle curve, all angles converged

    Companion code: radia_ngsolve.airgap_motor_workflow
      - age_motor_nonlinear_solve(): Picard loop for one rotor position
      - age_motor_rotation_sweep():  torque–angle curve

    Args:
        topic: Documentation topic. Options:
            "all"               - Complete documentation
            "overview"          - AGE core, Picard iteration, linear limit
            "api"               - age_motor_nonlinear_solve and rotation_sweep usage
            "saturation_models" - Froelich, piecewise linear, NGSolve BSpline B-H
            "validated"         - Test patterns: linear limit, saturation, sweep
    """
    return get_airgap_motor_workflow_documentation(topic)


@mcp.tool()
def dtn_coarse_mesh(topic: str = "all") -> str:
    """
    Why open-boundary methods stay accurate on COARSE meshes -- a spectral
    (DtN-matrix) explanation of the Kelvin-transformation coarse-mesh accuracy.

    Kameari demonstrated the Kelvin transform's coarse-mesh accuracy empirically,
    by mesh refinement.  This reframes it as a PROPERTY OF THE DtN MATRIX across
    mesh sizes: every open-BC closure (Kelvin / BEM / PML / Robin) approximates
    the one exterior DtN operator Λ_ext, whose sphere eigenvalues are the
    mesh-independent ladder λ_n = −(n+1)/R.  The discrete matrix Λ_h reproduces
    the LOW-degree eigenvalues accurately and almost independently of h even on
    the coarsest mesh; the per-mode error grows with degree n.  Since a compact
    source's field is dominated by low multipoles, the coarse mesh already
    resolves everything that matters -- Kameari's result, read off the spectrum.

    Companion code (both sides MEASURED -- two discretisations of the one Λ_ext):
      - bem_integral.exterior_dtn_spectrum():  eigenvalues of the BEM Λ_h matched
                                  to −(n+1)/R (boundary operator spectrum)
      - bem_integral.dtn_spectrum_vs_mesh():   per-degree error vs mesh size,
                                  with coarse-low-mode and accurate-band summary
      - fem_bem_coupling.kelvin_dtn_eigenvalue(): the Kelvin closure's effective
                                  DtN (volume FEM); order≥n kills the polynomial
                                  error, then a curved-geometry floor (~5-6 digits
                                  in 3D = Kameari's result; the dominant dipole
                                  inverts to a linear field, accurate at order 1)

    Args:
        topic: Documentation topic. Options:
            "all"          - Complete documentation
            "overview"     - One operator behind every open BC; the spectral argument
            "numerics"     - Measured DtN spectrum vs mesh table + how to read it
            "api"          - exterior_dtn_spectrum / dtn_spectrum_vs_mesh usage
            "applications" - Air-box sizing, method choice, trusting coarse Kelvin, debugging
            "p_method"     - Kelvin is a p-method not an h-method (measured p-vs-h,
                             ~20-80x DOF gap); polyhedron faceting error scales with
                             multipole degree (dipole robust); (R,p) design rule
            "formulation"  - Differential-geometry view: Omega vs A as Hodge-dual,
                             complementary dual bracketing (certified bounds), DtN
                             gradient block is formulation-independent, conformal
                             pullback material (scalar 0-form / tensor 1-form),
                             infinity = one-point conformal compactification, FEEC
            "datasheet"    - Problem-INDEPENDENT performance: open-BC error factors
                             into source multipoles x method eigenvalue-defect; the
                             Kelvin closure has an analytic, universal "datasheet"
                             (certify once, predict any problem by multipole content).
                             Includes COST measured by DoF INCREMENT (not time): the
                             closure adds a Gamma-scale coarse ball (measured dDoF=58,
                             closure error ~1/45 of the interior FE error) -> Kelvin is
                             cheap; keep it SPARSE (condensing to a dense DtN BC is
                             N^4/3, 10-20x nnz).
                             Includes the C/L dual: capacitance C <- n=0 (monopole,
                             exact) and external inductance L_ext <- n=1 (dipole) are
                             two Steklov modes of the SAME exterior DtN (W_ext =
                             1/2 mu0 (n+1)/R oint phi^2); certified via the magnetic
                             POTENTIAL exterior (A no-cut / Omega+cut), NOT the ngsbem
                             vector single-layer L=mu0 J^T(LaplaceSL)J (a different
                             operator with no -(n+1)/R ladder)
            "method_map"   - UNIFIED open-boundary map: Kelvin / BEM / PML / IABC /
                             CLN on the three selection axes (frequency / geometry /
                             space-vs-time) + the no-free-lunch modal axis + the
                             selection table. Audit-verified anchors.
                             doc: docs/open_boundary/OPEN_BOUNDARY_MAP.md
    """
    return get_dtn_coarse_mesh_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_ngsolve_simulation(description: str, formulation: str = "magnetostatic") -> str:
    """Create a new NGSolve electromagnetic simulation script."""
    return (
        f"Create an NGSolve simulation script: {description}\n"
        f"Formulation: {formulation}\n\n"
        "Follow these conventions:\n"
        "1. Use appropriate FE spaces from the de Rham complex:\n"
        "   - HCurl for vector potential A, electric field E\n"
        "   - HDiv for magnetic flux B, current J\n"
        "   - H1 for scalar potential Phi, temperature T\n"
        "   - Do NOT use VectorH1 for EM fields\n"
        "2. Magnetostatics: use HCurl(mesh, order=2, nograds=True)\n"
        "3. Eddy current: use complex=True on HCurl/H1 spaces\n"
        "4. BDDC: register Preconditioner BEFORE .Assemble()\n"
        "5. Do not overwrite x/y/z variables in loops\n"
        "6. Use .vec.data = for vector assignment\n"
        "7. 2D OCC: OCCGeometry(shape, dim=2)\n"
        "8. Saddle-point systems: use GMRes/MinRes, not CG\n"
        "9. BEM with HDivSurface: use .Trace() on trial/test functions\n"
    )


@mcp.prompt()
def ngsolve_eddy_current(geometry: str) -> str:
    """Set up an NGSolve eddy current / induction heating simulation."""
    return (
        f"Set up an NGSolve eddy current simulation for: {geometry}\n\n"
        "Use the mcp-server-ih ih_sibc tool for SIBC method selection.\n"
        "Key points:\n"
        "1. A-Phi formulation: HCurl(complex=True) * H1(complex=True)\n"
        "2. P_total: use BEM (ScalarBIESIBCSolver), not FEM BND integral\n"
        "3. Thermal: transient theta-scheme with H1 space (real)\n"
        "4. Kelvin transform for open boundary (bonus_intorder=4)\n"
    )


# ============================================================
# MCP Resources
# ============================================================

@mcp.resource("ngsolve://spaces")
def ngsolve_spaces_reference() -> str:
    """NGSolve FE space selection quick reference."""
    return (
        "# NGSolve FE Space Selection\n\n"
        "## de Rham Complex\n"
        "```\n"
        "H1 --grad--> HCurl --curl--> HDiv --div--> L2\n"
        "```\n\n"
        "| Space | Continuity | Use For |\n"
        "|-------|-----------|----------|\n"
        "| H1 | Full C^0 | Scalar potential Phi, temperature T |\n"
        "| HCurl | Tangential | Vector potential A, electric field E |\n"
        "| HDiv | Normal | Magnetic flux B, current density J |\n"
        "| HDivSurface | Normal (surface) | BEM surface currents |\n"
        "| SurfaceL2 | None (surface) | BEM charges |\n"
        "| VectorH1 | Full C^0 (all) | Elasticity (NOT for EM!) |\n\n"
        "## Common Parameters\n"
        "- `order=2`: polynomial order (default 1)\n"
        "- `nograds=True`: remove gradient null space (magnetostatics)\n"
        "- `complex=True`: complex-valued (eddy current, time-harmonic)\n"
        "- `dirichlet='bnd'`: essential BC on named boundary\n"
    )


@mcp.resource("ngsolve://solvers")
def ngsolve_solvers_reference() -> str:
    """NGSolve solver selection quick reference."""
    return (
        "# NGSolve Solver Selection\n\n"
        "## Direct Solvers\n"
        "| Solver | Strengths | When to Use |\n"
        "|--------|-----------|------------|\n"
        "| UMFPACK | Default, robust | N < 50K DOFs |\n"
        "| PARDISO | Fast, parallel | N > 50K, Intel MKL available |\n"
        "| MUMPS | Distributed, out-of-core | Very large, MPI |\n\n"
        "## Iterative Solvers\n"
        "| Solver | System Type | Preconditioner |\n"
        "|--------|------------|----------------|\n"
        "| CG | SPD only | BDDC, Jacobi, AMG |\n"
        "| MinRes | Symmetric indefinite | BDDC |\n"
        "| GMRes | General | Any |\n\n"
        "## Preconditioners\n"
        "- BDDC: domain decomposition, must register BEFORE .Assemble()\n"
        "- local (Jacobi/GS): simple, good for smoothing\n"
        "- multigrid: geometric or algebraic, best for H1\n"
    )


# ============================================================
# Panel Registry Tools
# ============================================================

def _load_panel_registry():
    """Load panel_registry.json from panels directory."""
    import json
    registry_path = (Path(__file__).parent.parent.parent / "panels"
                     / "panel_registry.json")
    if not registry_path.exists():
        return None
    with open(registry_path, "r", encoding="utf-8") as f:
        return json.load(f)


@mcp.tool()
def panel_schema(panel_name: str = "") -> str:
    """
    Show Radia-NGSolve panel definitions with Japanese labels and physics.

    When called without arguments, lists all available panels.
    When called with a panel name, shows detailed parameter definitions
    including Japanese names, physical meaning, CLI flags, and defaults.

    This enables natural language <-> CLI parameter mapping:
      "周波数を50kHzに" -> --frequency 50000
      "銅のワークピース" -> --material copper --sigma 5.8e7

    Args:
        panel_name: Panel ID (e.g. "inductance", "fem_kelvin").
                    Empty string returns overview of all panels.
    """
    reg = _load_panel_registry()
    if reg is None:
        return ("Error: panel_registry.json not found. "
                "Run: python src/radia/panels/sync_registry.py")

    panels = reg.get("panels", {})

    if not panel_name:
        lines = ["# Radia-NGSolve Panels\n"]
        for pid, p in panels.items():
            n = len(p.get("params", []))
            lines.append(f"## {pid}: {p['ja_name']}")
            lines.append(f"  {p['ja_description']}")
            lines.append(f"  Script: {p['script']} | Method: {p['method']}")
            lines.append(f"  Parameters: {n}")
            lines.append("")
        lines.append("Use panel_schema(panel_name) for parameter details.")
        return "\n".join(lines)

    if panel_name not in panels:
        return (f"Unknown panel: {panel_name}. "
                f"Available: {', '.join(panels.keys())}")

    p = panels[panel_name]
    lines = [
        f"# {p['ja_name']} ({panel_name})",
        f"Script: `{p['script']}` | Function: `{p['function']}`",
        f"Method: {p['method']}",
        f"Description: {p['ja_description']}",
        "",
        "## Parameters",
        "",
        "| CLI | 日本語 | Type | Default | Physics |",
        "|-----|--------|------|---------|---------|",
    ]
    for param in p.get("params", []):
        cli = param.get("cli", "")
        ja = param.get("ja", "")
        typ = param.get("type", "str")
        default = param.get("default", "")
        if param.get("required"):
            default = "**required**"
        physics = param.get("physics", "")
        choices = param.get("choices", [])
        if choices:
            physics = f"{physics} [{'/'.join(str(c) for c in choices)}]"
        lines.append(f"| `{cli}` | {ja} | {typ} | {default} | {physics} |")

    if p.get("command_builder"):
        lines.append(f"\nCommand builder: `{p['command_builder']}`")

    return "\n".join(lines)


@mcp.tool()
def panel_add_param(panel_name: str, param_name: str, param_type: str = "float",
                    cli_flag: str = "", default: str = "",
                    ja: str = "", physics: str = "",
                    help_text: str = "") -> str:
    """
    Plan where to add a new parameter to a Radia-NGSolve panel.

    Does NOT modify code. Returns a checklist of files and locations
    that need to be updated, so the LLM can make precise edits.

    Args:
        panel_name: Panel ID (e.g. "fem_kelvin", "inductance")
        param_name: Python parameter name (e.g. "coil_sigma")
        param_type: "float", "int", "str", "bool"
        cli_flag: CLI flag (e.g. "--coil-sigma"). Auto-generated if empty.
        default: Default value as string
        ja: Japanese label (e.g. "コイル導電率")
        physics: Physics description (e.g. "R = L/(sigma*A)")
        help_text: English help text for argparse
    """
    reg = _load_panel_registry()
    if reg is None:
        return "Error: panel_registry.json not found."

    panels = reg.get("panels", {})
    if panel_name not in panels:
        return f"Unknown panel: {panel_name}. Available: {', '.join(panels.keys())}"

    p = panels[panel_name]
    if not cli_flag:
        cli_flag = "--" + param_name.replace("_", "-")

    # Check if param already exists
    existing = [x["cli"] for x in p.get("params", [])]
    if cli_flag in existing:
        return f"Parameter {cli_flag} already exists in {panel_name}."

    script = p["script"]
    function = p["function"]
    builder = p.get("command_builder", "")

    lines = [
        f"# Add `{param_name}` to {panel_name} ({p['ja_name']})",
        f"  日本語: {ja}",
        f"  Physics: {physics}",
        "",
        "## Checklist (4 locations):",
        "",
        f"### 1. `panels/{script}` — argparse",
        f"  Add: `parser.add_argument(\"{cli_flag}\", type={param_type}, "
        f"default={default}, help=\"{help_text}\")`",
        "",
        f"### 2. `panels/{script}` — function `{function}()`",
        f"  Add parameter: `{param_name}: {param_type} = {default}`",
        f"  Wire: `args.{param_name.replace('-', '_')}` -> function call",
        "",
    ]

    if builder:
        mod, method = builder.split(":")
        lines.extend([
            f"### 3. `{mod}` — `{method}()`",
            f"  Add: `cmd += [\"{cli_flag}\", self.val(\"{param_name}\")]`",
            "",
            f"### 4. `{mod}` — widget definition",
            f"  Add QLineEdit/QSpinBox for `{param_name}`",
            f"  Label: \"{ja}\" (displayed in Qt panel)",
            "",
        ])
    else:
        lines.extend([
            "### 3-4. No command builder (standalone script)",
            "",
        ])

    lines.extend([
        "### 5. Update registry",
        f"  Run: `python panels/sync_registry.py`",
        "",
        "### 6. Update MCP knowledge (if physics-relevant)",
        f"  File: mcp_server knowledge related to {panel_name}",
        "",
        "### 7. Avoid the GUI pitfalls",
        "  Before committing, call `panel_gui_pitfalls()` and check",
        "  that the new param does not regress any of the listed",
        "  bugs (combo state save/restore, hidden-widget read in",
        "  build_command, mode-switch widget visibility, GMSH viz,",
        "  subprocess argparse choices, Cubit .jou id capture, ...).",
    ])

    return "\n".join(lines)


@mcp.tool()
def panel_describe_jp(panel_name: str) -> str:
    """
    現在のパネルソースを AST 解析して日本語で詳細に説明する。

    Reads the actual ``radia_<panel_name>.py`` source file (NOT the
    cached panel_registry.json which can be stale) and returns a
    Japanese hierarchical description of:

      - all widgets (key, label, type, default, combo items)
      - mode-switch visibility logic per handler
      - subprocess command builders with their CLI flag mapping

    Use this to:
      1. Confirm what the panel actually looks like before editing
      2. Generate a "spec" for the user to confirm in plain Japanese
      3. Diff against panel_registry.json to find drift

    Args:
        panel_name: Panel id (e.g. "ih", "em", "pcb"). Resolved to
                    ``src/radia/radia_<panel_name>.py``.

    Returns markdown text. Combine with panel_gui_pitfalls() output
    when planning a panel modification — first describe the current
    state, then check the relevant pitfalls.
    """
    path = _find_panel_file(panel_name)
    if path is None:
        return (f"Panel file not found for {panel_name!r}. "
                f"Expected at src/radia/radia_{panel_name}.py.")
    try:
        info = _parse_panel_file(path)
    except SyntaxError as e:
        return f"SyntaxError in {path}:{e.lineno}: {e.msg}"
    return _describe_panel_jp(info)


@mcp.tool()
def panel_widget_locations(panel_name: str, widget_key: str) -> str:
    """
    Return file:line locations for everything that touches a widget.

    For a given widget key (e.g. ``"half_thickness"``), returns:

      - **Definition** location: which add_line/add_combo/add_spin
        call created the widget, with the line number, default
        value, and combo items.
      - **Visibility rules**: every ``self._set_row_visible(key, ...)``
        call across all _on_*_changed handlers, with the conditional
        branch and the visibility expression.
      - **Command builder uses**: every ``cmd += ["--flag",
        self.val("key")]`` line in _build_*_command methods.

    Use this BEFORE editing a widget so you can update every
    location that references it in one consistent commit. The MCP
    output is JSON-pretty so the LLM can structure follow-up
    edits programmatically.

    Args:
        panel_name:  Panel id (e.g. "ih")
        widget_key:  Internal widget key (e.g. "half_thickness",
                     "wp_sigma", "method")
    """
    import json
    path = _find_panel_file(panel_name)
    if path is None:
        return f"Panel file not found for {panel_name!r}."
    try:
        info = _parse_panel_file(path)
    except SyntaxError as e:
        return f"SyntaxError in {path}:{e.lineno}: {e.msg}"
    locs = _widget_locations(info, widget_key)
    return json.dumps(locs, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_spec() -> str:
    """
    GMSH post-processing specification for Radia panels.

    Returns the SINGLE SOURCE OF TRUTH for what the GMSH output must
    look like: file format, physical groups, mesh curving, display
    options (.msh.opt), NodeData requirements, and the "kirei"
    reference from v3.6.1.

    Read this BEFORE writing any GMSH export code in calc_*.py.
    Every item is mandatory — no exceptions.
    """
    return get_gmsh_post_spec()


@mcp.tool()
def panel_gui_pitfalls(topic: str = "") -> str:
    """
    Pitfalls and lessons learned from Radia GUI / Cubit panel development.

    Read this BEFORE adding a new parameter, mode, or method to a
    `radia_*.py` panel, BEFORE renaming a combo item, and BEFORE
    writing a new sample .jou. Each pitfall is paired with a "rule"
    that prevents it from coming back.

    Topics:
      combo_state             -- save/restore by text, not index
      mode_switch             -- hidden widgets must not feed build_command
      layout_unification      -- shared widget set across solver methods
      gmsh_viz                -- companion .geo, hide volume mesh, vector only
      gmsh_arrow_size         -- ArrowSizeMin/Max=20 — without this the
                                 field arrows are functionally invisible
      subprocess_args         -- calc_*.py choices must match GUI combos
      cubit_jou               -- subtract id semantics, surface id renumbering
      sample_jou              -- one .jou per (panel, method) pair
      silent_action           -- menu actions must produce visible feedback
      silent_except           -- never bare-except; always log type+traceback
                                 tail; always provide a fallback path
      result_keys             -- subprocess result dict is an API contract
      regression_blast_radius -- run BOTH panels after touching shared
                                 helpers; opaque casts (PointId) bite
      panel_qt_testing        -- retired PySide note; current gate is
                                 validation_test/panels/test_notebook_workbench.py
      learn_edition_cap       -- ignore the 50k warning, export bypasses it

    Args:
        topic: Empty for the full document, or one of the topic
               keywords above for a single section.
    """
    return get_panel_gui_pitfalls(topic)


@mcp.tool()
def install_deploy(topic: str = "") -> str:
    """
    Radia install / deploy policy and recipes — 2-tier configuration
    (LAB + 100号機 editable / mdx + hibino PyPI), reversible migration steps, and the
    non-obvious gotchas that cause silent breakage.

    Read this when:
      * Setting up a new lab machine.
      * Migrating a machine between editable / PyPI install.
      * Diagnosing "import works but pip says wrong version" or
        "DLL load failed" on a freshly-deployed machine.

    Topics:
      two_tier                    -- current LAB-editable / PyPI-consumer policy
      lab_editable                -- LAB editable install
      hyaku_editable              -- 100号機 NAS editable install
      mdx_pypi                    -- mdx PyPI install
      hibino_pypi                 -- hibino PyPI install via `ssh hibino`
      editable_to_pypi_migration  -- e.g. 100号機 NAS-editable -> PyPI
      pypi_to_editable_migration  -- e.g. mdx PyPI -> editable
      metadata_sync               -- pip metadata vs radia.__version__
      pyd_dll_bootstrap           -- cubit_mesh_curver requires `import radia` first
      cubit_plugin_layers         -- Cubit plugin lives in TWO places
      common_failure_modes        -- symptoms and fixes table

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_install_deploy_documentation(topic)


@mcp.tool()
def release_workflow(topic: str = "") -> str:
    """
    Release-QUD workflow for the Radia monorepo
    (3 packages / 4 machines: LAB, 100号機, mdx, hibino). Documents the 9-phase
    pipeline, the 4 pre-flight gates added 2026-05-03, the historical
    CI failure modes + their root causes, and the patch-bump recovery
    protocol when a tag CI fails.

    Read this when:
      * The user asks for a release / version bump / PyPI publish.
      * CI on a tag ref fails after `git push --tags`.
      * "Release" workflow shows skipped in `gh run list` (CI never
        went green).
      * A user reports "is X.Y.Z on PyPI yet?" and propagation is
        stuck.

    Topics:
      overview               -- what gets released and why atomically
      phases                 -- the 9-phase pipeline (table)
      preflight_gates        -- Phase 2.5 4-gate pre-push validation
      ci_failure_modes       -- known CI failures + cause + fix table
      recovery               -- when CI on a tag fails AFTER push
      patch_bump_protocol    -- exact steps for retry after CI failure
      lab_lock_release       -- pre-deploy: stop processes that hold .pyd
      monorepo_lockstep      -- 4-6 version files that must stay in sync
      ci_monitor_skill       -- companion skill for Phase 7

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_release_workflow_documentation(topic)


@mcp.tool()
def standalone_panels(topic: str = "") -> str:
    """
    Retired standalone PySide panel topic.  The canonical Radia panel surface
    is now the Jupyter notebook workbench (`radia_<app>.ipynb` +
    `radia.<app>_notebook`).  This tool remains as a compatibility redirect.

    Read this when:
      * A user asks about the old standalone panel entry-points.
      * An MCP client still calls the historical `standalone_panels` topic.
      * You need the post-migration notebook route and no-PySide boundary.

    Topics:
      quick_start      -- current notebook route
      four_panels      -- active notebook workbenches
      build_notebook_gui -- construction recipe for new notebook GUIs
      cubit_panels_migration -- examples/cubit_panels promotion route
      vol_sources      -- Cubit / Netgen-OCC / build123d / etc.
      vs_cubit         -- notebook route vs Cubit export/plugin boundary
      ih_methods       -- IH through `radia_ih.ipynb`
      troubleshooting  -- common post-migration issues

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_standalone_panels_documentation(topic)


@mcp.tool()
def loop_learning(topic: str = "overview") -> str:
    """
    Public-safe CAE loop learning rules distilled from repeated validation
    rotations. Use this after a multi-tool loop has produced artifacts and the
    user asks whether the MCP server has actually learned from them.

    Topics:
      overview                -- artifact -> MCP learning workflow
      dual_lane               -- split one artifact into public/source-tool lanes
      mesh_geometry_vol       -- .vol, tri/tet, block registration, geometry checks
      force_moment            -- Lorentz, Maxwell traction, coenergy, moments
      motor_airgap_torque     -- Br/Bt harmonic torque phase and sign gates
      rf_acoustic_passivity   -- acoustic impedance and two-port passivity gates
      artifact_feedback       -- JSON/notebook/result artifact -> MCP knowledge
      mcp_closure             -- collected/distilled/encoded/verified/learned labels
      all                     -- complete document

    Args:
        topic: Topic name, or "all".
    """
    return get_loop_learning_documentation(topic)


@mcp.tool()
def basis_functions(topic: str = "") -> str:
    """
    Finite-element basis function library — Mathematica-canonical
    reference for H1, HCurl, HDiv, L2 on triangle / tet / quad / hex
    / prism / pyramid, arbitrary order.

    Phase 1 (radia-mcp 0.38.0) covers triangle + tetrahedron with
    H1 Lagrange P1-P5 (triangle) / P1-P3 (tet), HDiv RT₀ (= RWG used
    by BEM-A), and L2 P0-P3.  Mathematica notation is the canonical
    form; SymPy / NumPy translations are CI-verified in
    ``tests/basis/test_basis_functions.py``.

    Read this when:
      * Implementing a new BEM / FEM panel and need the canonical
        basis formula + a verified NumPy translation.
      * Cross-checking NGSolve's hierarchical output against a
        textbook Lagrange definition.
      * Debugging RWG sign / divergence issues in BEM-A.
      * Generating a code-gen template (NumPy / MATLAB / SymPy /
        Maple) from the Mathematica canonical form.

    Topics:
      theory_overview               -- what each space (H1/HCurl/HDiv/L2) means
      code_gen_pattern              -- Mathematica → NumPy / SymPy translation rules
      triangle_h1_lagrange          -- Triangle P1, P2, P3, P4, P5 (Lagrange nodal)
      triangle_h1_hierarchical      -- Triangle hierarchical (NGSolve compat)
      triangle_rwg                  -- Triangle HDiv RT₀ (= RWG, BEM-A workhorse)
      triangle_l2                   -- Triangle L2 P0-P3, Dubiner orthogonal
      tet_h1_lagrange               -- Tetrahedron P1, P2, P3
      tet_rwg                       -- Tetrahedron HDiv RT₀ (3D RWG, future)
      verification_recipes          -- partition of unity, div=σ/A, etc.

    Companion files:
      packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m
        -- canonical Mathematica package
      tests/basis/test_basis_functions.py
        -- CI-verified NumPy ports

    Args:
        topic: Empty for the full document, or one of the topics above.
    """
    return get_basis_functions_documentation(topic)


def _selftest():
    """Run lint on built-in fixtures to verify rules work correctly."""
    print("=" * 70)
    print("NGSolve Lint Self-Test")
    print("=" * 70)
    print()

    fixtures_dir = Path(__file__).parent.parent.parent.parent.parent / "tests" / "mcp_server" / "fixtures"
    if not fixtures_dir.exists():
        fixtures_dir = Path(__file__).parent / "fixtures"

    if fixtures_dir.exists():
        print(f"Using fixtures: {fixtures_dir}")
        print()
        py_files = sorted(fixtures_dir.glob("*.py"))
        total_findings = 0
        total_files = 0
        for py_file in py_files:
            findings = _lint_file(str(py_file))
            if findings:
                print(_format_findings(str(py_file), findings))
                print()
            total_findings += len(findings)
            total_files += 1

        print("=" * 70)
        print(f"Summary: {total_findings} finding(s) in {total_files} fixture file(s)")

        # Verify: bad_ngsolve should have findings, clean_ngsolve should not
        ngsolve_bad = fixtures_dir / "bad_ngsolve_script.py"
        if ngsolve_bad.exists():
            findings = _lint_file(str(ngsolve_bad))
            if len(findings) == 0:
                print(f"  WARNING: {ngsolve_bad.name} expected findings but got none")

        ngsolve_clean = fixtures_dir / "clean_ngsolve_script.py"
        if ngsolve_clean.exists():
            findings = _lint_file(str(ngsolve_clean))
            if len(findings) > 0:
                print(f"  FAIL: {ngsolve_clean.name} should be clean but got {len(findings)} finding(s)")
                sys.exit(1)

        print("Self-test PASSED")
    else:
        print("SKIP: No fixtures found.")




@mcp.tool()
def rf_sweep_artifact_summary_gate(
    summary_json: str,
    passivity_tolerance: float = 1.0e-3,
    reciprocity_tolerance: float = 1.0e-3,
) -> str:
    """Gate a solved two-port sweep artifact and its process-neutral metadata."""
    return json.dumps(
        _rf_sweep_artifact_summary_gate(
            summary_json, passivity_tolerance, reciprocity_tolerance
        ), indent=2, sort_keys=True
    )


@mcp.tool()
def physics_result_preflight_gate(summary_json: str) -> str:
    """Gate physics namespace, selection, solution, and license metadata before result evaluation."""
    try: result=_physics_result_preflight_gate(summary_json)
    except (json.JSONDecodeError,TypeError,ValueError) as exc: result={"policy":"physics_result_preflight_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)


@mcp.tool()
def cq_scattering_arrival_gate(
    time_step_s: float,
    geometric_arrival_s: float,
    measured_peak_s: float,
    max_relative_residual: float,
    finite_response: bool,
    real_time_response: bool,
    max_peak_lag_steps: float = 3.0,
    max_residual: float = 1.0e-6,
) -> str:
    """Gate CQ scattered-field causality against a geometric ray arrival."""
    try:
        result = _cq_scattering_arrival_gate(
            time_step_s=time_step_s,
            geometric_arrival_s=geometric_arrival_s,
            measured_peak_s=measured_peak_s,
            max_relative_residual=max_relative_residual,
            finite_response=finite_response,
            real_time_response=real_time_response,
            max_peak_lag_steps=max_peak_lag_steps,
            max_residual=max_residual,
        )
    except (TypeError, ValueError) as exc:
        result = {"policy": "cq_scattering_arrival_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def fsi_scattering_invariants_gate(
    reciprocity_relative_error: float,
    optical_theorem_relative_error: float,
    bem_dtn_relative_error: float,
    max_solver_residual: float,
    lossless_material: bool,
    time_convention: str,
    max_invariant_error: float = 0.05,
    max_bem_dtn_error: float = 0.05,
    max_solver_residual_allowed: float = 1.0e-8,
) -> str:
    """Gate lossless FSI reciprocity, energy closure, and exterior-method agreement."""
    try:
        result = _fsi_scattering_invariants_gate(
            reciprocity_relative_error=reciprocity_relative_error,
            optical_theorem_relative_error=optical_theorem_relative_error,
            bem_dtn_relative_error=bem_dtn_relative_error,
            max_solver_residual=max_solver_residual,
            lossless_material=lossless_material,
            time_convention=time_convention,
            max_invariant_error=max_invariant_error,
            max_bem_dtn_error=max_bem_dtn_error,
            max_solver_residual_allowed=max_solver_residual_allowed,
        )
    except (TypeError, ValueError) as exc:
        result = {"policy": "fsi_scattering_invariants_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def inductance_energy_mutual_gate(
    self_inductance: float,
    energy_inductance: float,
    mutual_inductance: float,
    analytic_mutual_inductance: float,
    inductance_unit: str,
    max_energy_relative_error: float = 0.01,
    max_mutual_relative_error: float = 0.05,
) -> str:
    """Gate L=2W/I^2 and an analytic one-direction mutual inductance."""
    try:
        result = _inductance_energy_mutual_gate(
            self_inductance=self_inductance,
            energy_inductance=energy_inductance,
            mutual_inductance=mutual_inductance,
            analytic_mutual_inductance=analytic_mutual_inductance,
            inductance_unit=inductance_unit,
            max_energy_relative_error=max_energy_relative_error,
            max_mutual_relative_error=max_mutual_relative_error,
        )
    except (TypeError, ValueError) as exc:
        result = {"policy": "inductance_energy_mutual_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def cq_response_reality_gate(
    summary_json: str,
    residual_tolerance: float = 1.0e-10,
    imaginary_tolerance: float = 1.0e-10,
) -> str:
    """Gate a coupled CQ solve, including its real time-domain reconstruction."""
    return json.dumps(_cq_response_reality_gate(
        json.loads(summary_json), residual_tolerance=residual_tolerance,
        imaginary_tolerance=imaginary_tolerance), indent=2, sort_keys=True)


@mcp.tool()
def dual_formulation_symmetric_field_profile_gate(
    summary_json: str,
    max_profile_relative_difference: float = 0.01,
    max_center_relative_difference: float = 0.01,
    max_symmetry_relative: float = 0.01,
    min_sample_count: int = 21,
) -> str:
    """Gate full-profile agreement and symmetry for two field formulations."""

    return json.dumps(_dual_formulation_symmetric_field_profile_gate(
        json.loads(summary_json),
        max_profile_relative_difference=max_profile_relative_difference,
        max_center_relative_difference=max_center_relative_difference,
        max_symmetry_relative=max_symmetry_relative,
        min_sample_count=min_sample_count,
    ), indent=2, sort_keys=True)


@mcp.tool()
def symmetric_complex_field_curve_gate(
    axis_positions: list[float],
    field_real: list[float],
    log10_relative_residual: float,
    field_imag: list[float] | None = None,
    axis_unit: str = "m",
    field_unit: str = "A/m",
    min_sample_count: int = 9,
    max_axis_symmetry_relative: float = 1.0e-9,
    max_field_symmetry_relative: float = 2.0e-3,
    max_log10_relative_residual: float = -8.0,
) -> str:
    """Gate an even- or odd-sampled complex field curve by mirror symmetry."""

    try:
        result = _symmetric_complex_field_curve_gate(
            axis_positions,
            field_real,
            field_imag,
            axis_unit=axis_unit,
            field_unit=field_unit,
            log10_relative_residual=log10_relative_residual,
            min_sample_count=min_sample_count,
            max_axis_symmetry_relative=max_axis_symmetry_relative,
            max_field_symmetry_relative=max_field_symmetry_relative,
            max_log10_relative_residual=max_log10_relative_residual,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "symmetric_complex_field_curve_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def helmholtz_double_layer_low_frequency_gate(summary_json: str) -> str:
    """Gate the quadratic low-frequency correction of a Helmholtz double layer."""

    try:
        result = _helmholtz_double_layer_low_frequency_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "helmholtz_double_layer_low_frequency_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def cyclic_terminal_source_sweep_gate(summary_json: str) -> str:
    """Gate cyclic terminal charges without assuming formulations are identical."""
    try: result=_cyclic_terminal_source_sweep_gate(json.loads(summary_json))
    except (TypeError,ValueError) as exc: result={"policy":"cyclic_terminal_source_sweep_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)


@mcp.tool()
def cyclic_terminal_phasor_balance_gate(
    summary_json: str,
    max_magnitude_relative_spread: float = 1.0e-5,
    max_phase_step_error_deg: float = 1.0e-2,
    max_zero_sequence_residual: float = 1.0e-5,
    max_terminal_kcl_residual: float = 1.0e-5,
    max_reference_current_relative_error: float = 2.0e-2,
) -> str:
    """Gate cyclic voltage/current triplets and all-terminal phasor KCL."""
    try:
        result = _cyclic_terminal_phasor_balance_gate(
            json.loads(summary_json),
            max_magnitude_relative_spread=max_magnitude_relative_spread,
            max_phase_step_error_deg=max_phase_step_error_deg,
            max_zero_sequence_residual=max_zero_sequence_residual,
            max_terminal_kcl_residual=max_terminal_kcl_residual,
            max_reference_current_relative_error=max_reference_current_relative_error,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "cyclic_terminal_phasor_balance_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def three_phase_winding_power_balance_gate(
    summary_json: str,
    max_voltage_relative_spread: float = 1.0e-5,
    max_current_relative_spread: float = 5.0e-3,
    max_phase_step_error_deg: float = 1.0,
    max_star_kcl_residual: float = 1.0e-5,
    max_active_power_relative_residual: float = 1.0e-3,
) -> str:
    """Gate three-phase balance, STAR KCL, and coupled-winding copper power."""
    try:
        result = _three_phase_winding_power_balance_gate(
            json.loads(summary_json),
            max_voltage_relative_spread=max_voltage_relative_spread,
            max_current_relative_spread=max_current_relative_spread,
            max_phase_step_error_deg=max_phase_step_error_deg,
            max_star_kcl_residual=max_star_kcl_residual,
            max_active_power_relative_residual=max_active_power_relative_residual,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "three_phase_winding_power_balance_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def cogging_torque_periodicity_gate(summary_json: str) -> str:
    """Gate a zero-current torque sweep over one slot/pole LCM period."""
    try:
        result = _cogging_torque_periodicity_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "cogging_torque_periodicity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def nonlinear_actuator_saturation_knee_gate(summary_json: str) -> str:
    """Gate an axisymmetric nonlinear actuator by a shared L/F saturation knee."""
    try:
        result = _nonlinear_actuator_saturation_knee_gate(json.loads(summary_json))
    except (TypeError, ValueError, KeyError) as exc:
        result = {"policy": "nonlinear_actuator_saturation_knee_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def source_free_static_null_solution_gate(
    summary_json: str,
    absolute_tolerance: float = 1.0e-14,
) -> str:
    """Gate a source-free static Maxwell solve against the exact zero solution."""
    try:
        result = _source_free_static_null_solution_gate(
            json.loads(summary_json),
            absolute_tolerance=absolute_tolerance,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "source_free_static_null_solution_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def harmonic_magnetic_force_triplet_closure_gate(
    summary_json: str,
    maximum_body_method_relative_difference: float = 0.05,
    maximum_action_reaction_relative_residual: float = 0.01,
    maximum_transverse_relative: float = 1.0e-8,
) -> str:
    """Gate harmonic body-force methods and source/body action-reaction closure."""
    try:
        result = _harmonic_magnetic_force_triplet_closure_gate(
            json.loads(summary_json),
            maximum_body_method_relative_difference=maximum_body_method_relative_difference,
            maximum_action_reaction_relative_residual=maximum_action_reaction_relative_residual,
            maximum_transverse_relative=maximum_transverse_relative,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "harmonic_magnetic_force_triplet_closure_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def magnetic_force_method_profile_gate(
    summary_json: str,
    maximum_method_relative_difference: float = 0.05,
    maximum_independent_stress_relative_difference: float = 0.02,
    minimum_selection_scope_relative_difference: float = 0.25,
    maximum_all_body_to_target_magnitude_ratio: float = 0.75,
    maximum_work_relative_difference: float = 0.05,
    maximum_parsed_replay_absolute_difference: float = 1.0e-12,
    minimum_sample_count: int = 5,
) -> str:
    """Gate magnetic-force profiles with explicit body/surface selection scope."""
    try:
        result = _magnetic_force_method_profile_gate(
            json.loads(summary_json),
            maximum_method_relative_difference=maximum_method_relative_difference,
            maximum_independent_stress_relative_difference=maximum_independent_stress_relative_difference,
            minimum_selection_scope_relative_difference=minimum_selection_scope_relative_difference,
            maximum_all_body_to_target_magnitude_ratio=maximum_all_body_to_target_magnitude_ratio,
            maximum_work_relative_difference=maximum_work_relative_difference,
            maximum_parsed_replay_absolute_difference=maximum_parsed_replay_absolute_difference,
            minimum_sample_count=minimum_sample_count,
        )
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "magnetic_force_method_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def harmonic_current_port_power_energy_identity_gate(
    summary_json: str,
    maximum_identity_relative_error: float = 1.0e-9,
    maximum_cross_run_relative_error: float = 1.0e-9,
) -> str:
    """Gate peak-phasor port, loss, energy, flux, and profile identities."""
    try:
        result = _harmonic_current_port_power_energy_identity_gate(
            json.loads(summary_json),
            maximum_identity_relative_error=maximum_identity_relative_error,
            maximum_cross_run_relative_error=maximum_cross_run_relative_error,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "harmonic_current_port_power_energy_identity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def periodic_unwrapped_pm_machine_replay_gate(
    summary_json: str,
    maximum_field_symmetry_relative_error: float = 0.05,
    maximum_energy_replay_relative_error: float = 1.0e-3,
    maximum_field_replay_relative_error: float = 5.0e-3,
    maximum_mesh_cardinality_relative_difference: float = 0.05,
) -> str:
    """Gate topology-aware PM-machine field symmetry and replay stability."""
    try:
        result = _periodic_unwrapped_pm_machine_replay_gate(
            json.loads(summary_json),
            maximum_field_symmetry_relative_error=maximum_field_symmetry_relative_error,
            maximum_energy_replay_relative_error=maximum_energy_replay_relative_error,
            maximum_field_replay_relative_error=maximum_field_replay_relative_error,
            maximum_mesh_cardinality_relative_difference=maximum_mesh_cardinality_relative_difference,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "periodic_unwrapped_pm_machine_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def permanent_magnet_recoil_state_gate(
    summary_json: str,
    maximum_replay_relative_error: float = 1.0e-3,
    minimum_initial_axis_concentration: float = 100.0,
    maximum_open_axis_concentration: float = 2.0,
    minimum_recoil_axis_concentration: float = 10.0,
) -> str:
    """Gate nonlinear, open-circuit, and partial-recoil PM field states."""
    try:
        result = _permanent_magnet_recoil_state_gate(
            json.loads(summary_json),
            maximum_replay_relative_error=maximum_replay_relative_error,
            minimum_initial_axis_concentration=minimum_initial_axis_concentration,
            maximum_open_axis_concentration=maximum_open_axis_concentration,
            minimum_recoil_axis_concentration=minimum_recoil_axis_concentration,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "permanent_magnet_recoil_state_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linear_eddy_levitation_force_gate(summary_json: str) -> str:
    """Gate linear harmonic levitation force by dual extraction and I-squared laws."""
    try:
        result = _linear_eddy_levitation_force_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linear_eddy_levitation_force_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def motion_coupled_eddy_levitation_transient_gate(summary_json: str) -> str:
    """Gate motion-coupled lift while detecting aliased force output times."""
    try:
        result = _motion_coupled_eddy_levitation_transient_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "motion_coupled_eddy_levitation_transient_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def harmonic_zero_net_circuit_gate(summary_json: str) -> str:
    """Gate zero-net harmonic phasors, Faraday sign, loss, and force metadata."""
    try:
        result = _harmonic_zero_net_circuit_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "harmonic_zero_net_circuit_faraday_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def moving_conductor_eddy_brake_gate(summary_json: str) -> str:
    """Gate motion, Lorentz-force, and Joule-loss table identities."""
    try:
        result = _moving_conductor_eddy_brake_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "moving_conductor_eddy_brake_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def rotating_conductor_transient_gate(summary_json: str) -> str:
    """Gate moving-axis migration, rotational kinematics, and loss partition."""
    try:
        result = _rotating_conductor_transient_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "rotating_conductor_transient_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linear_magnetization_scaling_gate(summary_json: str) -> str:
    """Gate source scaling plus an independent refined P1 FEM reference."""
    try:
        result = _linear_magnetization_scaling_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linear_magnetization_scaling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linear_axisymmetric_circuit_energy_gate(summary_json: str) -> str:
    """Gate current, flux, field, and energy identities on one fixed mesh."""
    try:
        result = _linear_axisymmetric_circuit_energy_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linear_axisymmetric_circuit_energy_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def manual_auto_mixed_mesh_preservation_gate(summary_json: str) -> str:
    """Gate exact manual-region preservation and bounded automatic remeshing."""
    try:
        result = _manual_auto_mixed_mesh_preservation_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "manual_auto_mixed_mesh_preservation_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def two_winding_frequency_faraday_gate(summary_json: str) -> str:
    """Gate two-winding complex response against linked-flux Faraday identity."""
    try:
        result = _two_winding_frequency_faraday_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "two_winding_frequency_faraday_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def magnetic_conductive_shield_frequency_gate(summary_json: str) -> str:
    """Gate low-frequency magnetic loading and high-frequency eddy shielding."""
    try:
        result = _magnetic_conductive_shield_frequency_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "magnetic_conductive_shield_frequency_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def cylindrical_conductor_skin_bessel_gate(summary_json: str) -> str:
    """Gate cylindrical skin-effect identities and exact Bessel structure."""
    try:
        result = _cylindrical_conductor_skin_bessel_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "cylindrical_conductor_skin_bessel_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def conductive_network_resistance_monotonicity_gate(summary_json: str) -> str:
    """Gate Rayleigh resistance monotonicity for conductive contact networks."""
    try:
        result = _conductive_network_resistance_monotonicity_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "conductive_network_resistance_monotonicity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def autodiff_harmonic_balance_convergence_gate(summary_json: str) -> str:
    """Gate AD harmonic balance without mean-only false convergence."""
    try:
        result = _autodiff_harmonic_balance_convergence_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "autodiff_harmonic_balance_convergence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def hall_effect_transverse_voltage_gate(summary_json: str) -> str:
    """Gate Hall voltage by coefficient, drive, field, and replay controls."""
    try:
        result = _hall_effect_transverse_voltage_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "hall_effect_transverse_voltage_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def single_loop_source_normalized_field_gate(summary_json: str) -> str:
    """Gate a single-loop field transfer across two port formulations."""
    try:
        result = _single_loop_source_normalized_field_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "single_loop_source_normalized_field_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def coupled_cq_refinement_gate(summary_json: str) -> str:
    """Gate coupled FEM/BEM CQ symbols, contour balance, and refinement."""
    try:
        result = _coupled_cq_refinement_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "coupled_cq_refinement_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def coil_self_resonance_sweep_gate(summary_json: str) -> str:
    """Gate complex coil impedance, self-resonance, and sweep replay."""
    try:
        result = _coil_self_resonance_sweep_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "coil_self_resonance_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def passive_axial_bearing_stiffness_gate(summary_json: str) -> str:
    """Gate signed force, action-reaction, axial stability, and sweep replay."""
    try:
        result = _passive_axial_bearing_stiffness_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "passive_axial_bearing_stiffness_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def radial_bearing_force_symmetry_gate(summary_json: str) -> str:
    """Gate magnetic-body force with equal and mirrored excitation controls."""
    try:
        result = _radial_bearing_force_symmetry_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "radial_bearing_force_symmetry_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def complex_vector_field_maximum_gate(summary_json: str) -> str:
    """Gate complex vector-field magnitudes and per-material maxima."""
    try:
        result = _complex_vector_field_maximum_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "complex_vector_field_maximum_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def one_port_vi_s_impedance_gate(summary_json: str) -> str:
    """Gate one-port S, V/I, impedance-transform, and power identities."""
    try:
        result = _one_port_vi_s_impedance_gate(json.loads(summary_json))
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "one_port_vi_s_impedance_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def symmetric_axial_field_profile_gate(
    axis_positions: list[float],
    axial_field: list[float],
    expected_center_field: float,
    transverse_field_1: list[float] | None = None,
    transverse_field_2: list[float] | None = None,
    min_sample_count: int = 5,
    max_center_relative_error: float = 1.0e-6,
    max_symmetry_relative: float = 1.0e-9,
    max_transverse_relative: float = 1.0e-9,
    max_axis_symmetry_relative: float = 1.0e-9,
) -> str:
    """Gate an origin-centered axial profile by analytic value and symmetry."""

    try:
        result = _symmetric_axial_field_profile_gate(
            axis_positions,
            axial_field,
            expected_center_field=expected_center_field,
            transverse_field_1=transverse_field_1,
            transverse_field_2=transverse_field_2,
            min_sample_count=min_sample_count,
            max_center_relative_error=max_center_relative_error,
            max_symmetry_relative=max_symmetry_relative,
            max_transverse_relative=max_transverse_relative,
            max_axis_symmetry_relative=max_axis_symmetry_relative,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "symmetric_axial_field_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def finite_solenoid_surface_current_gate(summary_json: str) -> str:
    """Gate a finite-solenoid surface-current profile and signed linearity."""

    try:
        result = _finite_solenoid_surface_current_gate(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "finite_solenoid_surface_current_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def parallel_wire_force_refinement_gate(
    refinement_levels: list[float],
    force_wire1_rows: list[list[float]],
    force_wire2_rows: list[list[float]],
    expected_force_magnitude: float,
    separation_direction: list[float] | None = None,
    expected_wire2_radial_sign: int | None = None,
    min_sample_count: int = 3,
    max_final_relative_error: float = 0.01,
    max_final_pair_relative_residual: float = 0.01,
    max_final_transverse_relative_force: float = 0.01,
    min_initial_to_final_error_ratio: float = 1.2,
) -> str:
    """Gate a reciprocal two-wire force refinement sweep without requiring monotone error."""

    try:
        result = _parallel_wire_force_refinement_gate(
            refinement_levels,
            force_wire1_rows,
            force_wire2_rows,
            expected_force_magnitude=expected_force_magnitude,
            separation_direction=separation_direction or [1.0, 0.0],
            expected_wire2_radial_sign=expected_wire2_radial_sign,
            min_sample_count=min_sample_count,
            max_final_relative_error=max_final_relative_error,
            max_final_pair_relative_residual=max_final_pair_relative_residual,
            max_final_transverse_relative_force=max_final_transverse_relative_force,
            min_initial_to_final_error_ratio=min_initial_to_final_error_ratio,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "parallel_wire_force_refinement_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def two_conductor_capacitance_identity_gate(
    conductor_voltages_v: list[float],
    conductor_charges_c: list[float],
    stored_energy_j: float,
    driven_conductor_index: int = 1,
    planar_depth_m: float | None = None,
    max_capacitance_relative_error: float = 1.0e-5,
    max_charge_balance_relative_error: float = 1.0e-5,
) -> str:
    """Gate two-conductor capacitance using terminal charge and field energy."""

    try:
        result = _two_conductor_capacitance_identity_gate(
            conductor_voltages_v,
            conductor_charges_c,
            stored_energy_j,
            driven_conductor_index=driven_conductor_index,
            planar_depth_m=planar_depth_m,
            max_capacitance_relative_error=max_capacitance_relative_error,
            max_charge_balance_relative_error=max_charge_balance_relative_error,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "two_conductor_capacitance_identity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def two_conductor_capacitance_matrix_gate(summary_json: str) -> str:
    """Gate reciprocal Maxwell and mutual capacitance matrix representations."""
    try:
        result = _two_conductor_capacitance_matrix_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {"policy":"two_conductor_capacitance_matrix_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def multiconductor_capacitance_cross_formulation_gate(summary_json: str) -> str:
    """Gate N-conductor Maxwell matrices across volume and boundary formulations."""

    try:
        result = _multiconductor_capacitance_cross_formulation_gate(
            json.loads(summary_json)
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "multiconductor_capacitance_cross_formulation_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def force_position_profile_gate(
    positions: list[float],
    forces: list[float],
    node_counts: list[int] | None = None,
    element_counts: list[int] | None = None,
    min_sample_count: int = 5,
    max_mesh_count_relative_span: float = 0.02,
    require_interior_peak: bool = False,
    require_nonnegative: bool = False,
) -> str:
    """Gate a force-position sweep without assuming it is monotonic."""

    try:
        result = _force_position_profile_gate(
            positions,
            forces,
            node_counts=node_counts,
            element_counts=element_counts,
            min_sample_count=min_sample_count,
            max_mesh_count_relative_span=max_mesh_count_relative_span,
            require_interior_peak=require_interior_peak,
            require_nonnegative=require_nonnegative,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "force_position_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def force_coenergy_displacement_gate(
    positions_m: list[float],
    coenergy_j: list[float],
    forces_along_displacement_n: list[float],
    energy_kind: str = "constant_current_coenergy",
    max_central_relative_error: float = 0.02,
    min_sample_count: int = 5,
) -> str:
    """Gate direct force against the central derivative of magnetic coenergy."""

    try:
        result = _force_coenergy_displacement_gate(
            positions_m,
            coenergy_j,
            forces_along_displacement_n,
            energy_kind=energy_kind,
            max_central_relative_error=max_central_relative_error,
            min_sample_count=min_sample_count,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "force_coenergy_displacement_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def rotational_kinematics_time_axis_gate(
    time_values: list[float],
    angles_deg: list[float],
    speeds_rpm: list[float],
    reported_time_unit: str,
    time_value_basis: str = "si_seconds",
    max_central_relative_error: float = 1.0e-8,
    min_sample_count: int = 5,
) -> str:
    """Gate a result-table time axis using angle/speed kinematics."""
    try:
        result = _rotational_kinematics_time_axis_gate(
            time_values,
            angles_deg,
            speeds_rpm,
            reported_time_unit=reported_time_unit,
            time_value_basis=time_value_basis,
            max_central_relative_error=max_central_relative_error,
            min_sample_count=min_sample_count,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "rotational_kinematics_time_axis_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def inductance_matrix_family_gate(
    cases: list[dict],
    expected_strongest_coupling_case: str | None = None,
    max_reciprocity_relative_error: float = 0.02,
    psd_relative_tolerance: float = 1.0e-12,
    max_identity_relative_error: float = 1.0e-6,
    max_replay_relative_error: float = 1.0e-9,
    max_turn_scaling_relative_error: float = 0.02,
) -> str:
    """Gate two-winding matrices, identities, replay, and turn scaling."""
    try:
        result = _inductance_matrix_family_gate(
            cases,
            expected_strongest_coupling_case=expected_strongest_coupling_case,
            max_reciprocity_relative_error=max_reciprocity_relative_error,
            psd_relative_tolerance=psd_relative_tolerance,
            max_identity_relative_error=max_identity_relative_error,
            max_replay_relative_error=max_replay_relative_error,
            max_turn_scaling_relative_error=max_turn_scaling_relative_error,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "inductance_matrix_family_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linear_sphere_geometry_convergence_gate(
    rows: list[dict],
    analytic_volume: float,
    analytic_surface_area: float,
    replay: dict,
    max_reader_relative_error: float = 1.0e-12,
    max_surface_radius_error: float = 1.0e-12,
    max_final_geometry_relative_error: float = 3.0e-3,
    min_asymptotic_order: float = 1.8,
) -> str:
    """Gate first-order sphere tri/tet geometry convergence and replay."""
    try:
        result = _linear_sphere_geometry_convergence_gate(
            rows,
            analytic_volume=analytic_volume,
            analytic_surface_area=analytic_surface_area,
            replay=replay,
            max_reader_relative_error=max_reader_relative_error,
            max_surface_radius_error=max_surface_radius_error,
            max_final_geometry_relative_error=max_final_geometry_relative_error,
            min_asymptotic_order=min_asymptotic_order,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linear_sphere_geometry_convergence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def leakage_inductance_closure_gate(summary_json: str) -> str:
    """Gate compensated-energy and unit-current-matrix leakage inductance."""
    try:
        result = _leakage_inductance_closure_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "leakage_inductance_closure_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def multiport_impedance_sweep_gate(
    frequency_rows: list[list[float]],
    impedance_real_rows: list[list[float]],
    impedance_imag_rows: list[list[float]],
    port_ids: list[str] | None = None,
    min_sample_count: int = 5,
    min_frequency_decades: float = 1.0,
    passive_real_tolerance: float = 1.0e-9,
) -> str:
    """Gate common-grid, positive-real, nontrivial complex impedance sweeps."""

    try:
        result = _multiport_impedance_sweep_gate(
            frequency_rows,
            impedance_real_rows,
            impedance_imag_rows,
            port_ids=port_ids,
            min_sample_count=min_sample_count,
            min_frequency_decades=min_frequency_decades,
            passive_real_tolerance=passive_real_tolerance,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "multiport_impedance_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def one_port_power_balance_sweep_gate(
    summary_json: str,
    max_power_relative_residual: float = 1.0e-9,
    max_balance_abs_residual: float = 1.0e-9,
    max_reference_impedance_relative_drift: float = 1.0e-9,
) -> str:
    """Gate passive one-port accepted power against S11 and reference impedance."""

    try:
        result = _one_port_power_balance_gate(
            summary_json,
            max_power_relative_residual=max_power_relative_residual,
            max_balance_abs_residual=max_balance_abs_residual,
            max_reference_impedance_relative_drift=max_reference_impedance_relative_drift,
        )
    except (KeyError, TypeError, ValueError) as exc:
        result = {"policy": "one_port_power_balance_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def radar_range_angle_localization_gate(
    frequency_hz: list[float],
    targets_json: str,
    max_range_resolution_multiples: float = 1.0,
    max_angle_error_deg: float = 2.0,
    max_frequency_step_relative_drift: float = 1.0e-7,
) -> str:
    """Gate wideband range-angle localization of multiple targets."""

    try:
        targets = json.loads(targets_json)
        result = _radar_range_angle_localization_gate(
            frequency_hz,
            targets,
            max_range_resolution_multiples=max_range_resolution_multiples,
            max_angle_error_deg=max_angle_error_deg,
            max_frequency_step_relative_drift=max_frequency_step_relative_drift,
        )
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        result = {"policy": "radar_range_angle_localization_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def radar_range_rcs_profile_gate(
    frequency_hz: list[float],
    target_range_m: float,
    radar_peak_range_m: float,
    radar_peak_rcs_m2: float,
    generalized_peak_range_m: float,
    generalized_peak_rcs_m2: float,
    analytic_peak_rcs_m2: float,
    profile_relative_l2: float,
    max_frequency_step_relative_drift: float = 1.0e-7,
    max_peak_range_resolution_multiples: float = 1.0,
    max_profile_relative_l2: float = 1.0e-4,
    max_method_peak_relative_error: float = 1.0e-4,
    max_analytic_peak_relative_error: float = 0.05,
) -> str:
    """Gate wideband range-RCS localization, method agreement, and analytic amplitude."""

    try:
        result = _radar_range_rcs_profile_gate(
            frequency_hz,
            target_range_m=target_range_m,
            radar_peak_range_m=radar_peak_range_m,
            radar_peak_rcs_m2=radar_peak_rcs_m2,
            generalized_peak_range_m=generalized_peak_range_m,
            generalized_peak_rcs_m2=generalized_peak_rcs_m2,
            analytic_peak_rcs_m2=analytic_peak_rcs_m2,
            profile_relative_l2=profile_relative_l2,
            max_frequency_step_relative_drift=max_frequency_step_relative_drift,
            max_peak_range_resolution_multiples=max_peak_range_resolution_multiples,
            max_profile_relative_l2=max_profile_relative_l2,
            max_method_peak_relative_error=max_method_peak_relative_error,
            max_analytic_peak_relative_error=max_analytic_peak_relative_error,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "radar_range_rcs_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def hmatrix_compression_scaling_gate(
    rows_json: str,
    max_matvec_relative_error: float = 1.0e-8,
    max_rank: int = 20,
    max_storage_growth_exponent: float = 1.25,
    min_dense_growth_exponent: float = 1.9,
) -> str:
    """Gate H-matrix accuracy, bounded rank, and subquadratic storage scaling."""

    try:
        rows = json.loads(rows_json)
        result = _hmatrix_compression_scaling_gate(
            rows,
            max_matvec_relative_error=max_matvec_relative_error,
            max_rank=max_rank,
            max_storage_growth_exponent=max_storage_growth_exponent,
            min_dense_growth_exponent=min_dense_growth_exponent,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "hmatrix_compression_scaling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def acoustic_duct_band_gap_gate(summary_json: str) -> str:
    """Gate a confined acoustic band gap against empty and free-space controls."""

    try:
        result = _acoustic_duct_band_gap_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "acoustic_duct_band_gap_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def dual_formulation_force_error_convergence_gate(
    formulation_rows: list[dict],
    reference_force: float,
    max_final_relative_error: float = 0.02,
    min_initial_to_final_improvement: float = 1.1,
    max_final_to_best_error_ratio: float = 1.5,
    max_tail_relative_span: float = 0.005,
) -> str:
    """Gate force-error convergence envelopes across two or more formulations."""

    try:
        result = _dual_formulation_force_error_convergence_gate(
            formulation_rows,
            reference_force=reference_force,
            max_final_relative_error=max_final_relative_error,
            min_initial_to_final_improvement=min_initial_to_final_improvement,
            max_final_to_best_error_ratio=max_final_to_best_error_ratio,
            max_tail_relative_span=max_tail_relative_span,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "dual_formulation_force_error_convergence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def magnetostatic_open_boundary_equivalence_gate(
    formulation_rows: list[dict],
    physics_regime: str = "magnetostatic_open_boundary",
    axis_sample_indices: list[int] | None = None,
    max_dominant_b_relative_error: float = 0.01,
    max_axis_transverse_b_residual: float = 0.015,
    max_a_offset_relative_spread: float = 0.005,
    max_energy_coenergy_relative_error: float = 0.001,
    max_dominant_force_relative_error: float = 0.001,
    max_force_balance_relative: float = 0.002,
    max_transverse_force_difference_relative: float = 0.001,
) -> str:
    """Gate gauge-invariant equivalence of two magnetostatic open-boundary solutions."""

    try:
        result = _magnetostatic_open_boundary_equivalence_gate(
            formulation_rows,
            physics_regime=physics_regime,
            axis_sample_indices=axis_sample_indices,
            max_dominant_b_relative_error=max_dominant_b_relative_error,
            max_axis_transverse_b_residual=max_axis_transverse_b_residual,
            max_a_offset_relative_spread=max_a_offset_relative_spread,
            max_energy_coenergy_relative_error=max_energy_coenergy_relative_error,
            max_dominant_force_relative_error=max_dominant_force_relative_error,
            max_force_balance_relative=max_force_balance_relative,
            max_transverse_force_difference_relative=max_transverse_force_difference_relative,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "magnetostatic_open_boundary_equivalence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def pwm_controlled_motor_loss_gate(
    payload: dict,
    max_three_phase_kcl_relative_error: float = 1.0e-3,
    max_tail_control_tracking_rms_relative_error: float = 0.05,
    max_angle_speed_integral_relative_error: float = 2.0e-4,
    max_power_sum_relative_error: float = 1.0e-10,
    max_loss_identity_relative_error: float = 1.0e-10,
    max_frequency_step_relative_span: float = 1.0e-10,
) -> str:
    """Gate PWM current-control and aggregate/harmonic loss-table identities."""

    try:
        result = _pwm_controlled_motor_loss_gate(
            payload,
            max_three_phase_kcl_relative_error=max_three_phase_kcl_relative_error,
            max_tail_control_tracking_rms_relative_error=max_tail_control_tracking_rms_relative_error,
            max_angle_speed_integral_relative_error=max_angle_speed_integral_relative_error,
            max_power_sum_relative_error=max_power_sum_relative_error,
            max_loss_identity_relative_error=max_loss_identity_relative_error,
            max_frequency_step_relative_span=max_frequency_step_relative_span,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "pwm_controlled_motor_loss_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def material_contrast_force_gate(
    cases: list[dict],
    interaction_axis: str = "x",
    max_background_relative_force: float = 0.01,
    max_transverse_relative_force: float = 1.0e-6,
    min_stronger_repulsion_ratio: float = 1.5,
) -> str:
    """Gate null, attraction, and increasing-repulsion material-force cases."""

    try:
        result = _material_contrast_force_gate(
            cases,
            interaction_axis=interaction_axis,
            max_background_relative_force=max_background_relative_force,
            max_transverse_relative_force=max_transverse_relative_force,
            min_stronger_repulsion_ratio=min_stronger_repulsion_ratio,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "material_contrast_force_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def two_body_force_magnitude_replay_gate(summary_json: str) -> str:
    """Gate unsigned two-body force balance and two fresh solver replays."""
    try:
        result = _two_body_force_magnitude_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "two_body_force_magnitude_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def static_field_shim_family_gate(
    cases: list[dict],
    interaction_axis: str = "z",
    min_paired_source_field_ratio: float = 1.2,
    max_paired_source_uniformity_ratio: float = 0.5,
    min_shim_center_field_delta_relative: float = 0.01,
    max_center_transverse_relative: float = 1.0e-4,
    max_central_divergence_relative: float = 0.05,
) -> str:
    """Gate static-field scale, ROI uniformity, shim sensitivity, and map quality."""

    try:
        result = _static_field_shim_family_gate(
            cases,
            interaction_axis=interaction_axis,
            min_paired_source_field_ratio=min_paired_source_field_ratio,
            max_paired_source_uniformity_ratio=max_paired_source_uniformity_ratio,
            min_shim_center_field_delta_relative=min_shim_center_field_delta_relative,
            max_center_transverse_relative=max_center_transverse_relative,
            max_central_divergence_relative=max_central_divergence_relative,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "static_field_shim_family_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def energy_budgeted_trace_kkt_gate(
    payload: dict,
    max_gradient_relative_error: float = 1.0e-6,
    max_solution_relative_error: float = 1.0e-5,
    max_stationarity_inf: float = 1.0e-6,
    max_complementarity_abs: float = 1.0e-7,
    max_constraint_relative: float = 1.0e-8,
) -> str:
    """Gate KKT closure for an energy-budgeted FEM/BEM trace fit."""

    try:
        result = _energy_budgeted_trace_kkt_gate(
            payload,
            max_gradient_relative_error=max_gradient_relative_error,
            max_solution_relative_error=max_solution_relative_error,
            max_stationarity_inf=max_stationarity_inf,
            max_complementarity_abs=max_complementarity_abs,
            max_constraint_relative=max_constraint_relative,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "energy_budgeted_trace_kkt_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def voice_coil_force_flux_sweep_gate(
    rows: list[dict],
    max_zero_force_relative: float = 0.01,
    max_odd_residual_relative: float = 0.10,
    max_force_constant_relative_span: float = 0.12,
) -> str:
    """Gate a PM voice-coil current sweep by force, flux, symmetry, and mesh evidence."""

    try:
        result = _voice_coil_force_flux_sweep_gate(
            rows,
            max_zero_force_relative=max_zero_force_relative,
            max_odd_residual_relative=max_odd_residual_relative,
            max_force_constant_relative_span=max_force_constant_relative_span,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "voice_coil_force_flux_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linear_induction_frequency_sweep_gate(
    rows: list[dict],
    thrust_abs_tol_n: float = 0.75,
    thrust_rel_tol: float = 2.0e-3,
    phase_balance_atol_a: float = 1.0e-9,
) -> str:
    """Gate a linear-induction frequency sweep by thrust, loss, and phase balance."""

    try:
        result = _linear_induction_frequency_sweep_gate(
            rows,
            thrust_abs_tol_n=thrust_abs_tol_n,
            thrust_rel_tol=thrust_rel_tol,
            phase_balance_atol_a=phase_balance_atol_a,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linear_induction_frequency_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def twin_conductor_skin_effect_frequency_gate(
    frequencies_hz: list[float],
    resistance_ohm: list[list[float]],
    inductance_h: list[list[float]],
    symmetry_rtol: float = 5.0e-4,
) -> str:
    """Gate passive twin-conductor R/L and impedance trends over frequency."""

    try:
        result = _twin_conductor_skin_effect_frequency_gate(
            frequencies_hz,
            resistance_ohm,
            inductance_h,
            symmetry_rtol=symmetry_rtol,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "twin_conductor_skin_effect_frequency_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def opposed_busbar_skin_force_gate(
    rows: list[dict],
    conductor_thickness_mm: float,
    conductivity_s_per_m: float,
    commanded_current_a: float,
    replay_rtol: float = 1.0e-12,
    identity_rtol: float = 5.0e-8,
    force_balance_rtol: float = 5.0e-5,
) -> str:
    """Gate AC skin/proximity, phasor identities, and Lorentz action-reaction."""

    try:
        result = _opposed_busbar_skin_force_gate(
            rows,
            conductor_thickness_mm=conductor_thickness_mm,
            conductivity_s_per_m=conductivity_s_per_m,
            commanded_current_a=commanded_current_a,
            replay_rtol=replay_rtol,
            identity_rtol=identity_rtol,
            force_balance_rtol=force_balance_rtol,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "opposed_busbar_skin_force_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def transient_conductor_replay_identity_gate(
    summary_json: str,
    identity_rtol: float = 1.0e-10,
    equivalence_rtol: float = 1.0e-12,
) -> str:
    """Gate full transient conductor histories, identities, and independent replay."""

    try:
        result = _transient_conductor_replay_identity_gate(
            json.loads(summary_json),
            identity_rtol=identity_rtol,
            equivalence_rtol=equivalence_rtol,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "transient_conductor_replay_identity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def homogenized_bundle_impedance_comparison_gate(
    rows: list[dict],
    resistance_rtol: float = 0.03,
    inductance_rtol: float = 0.005,
    impedance_rtol: float = 0.01,
    observable_rtol: float = 1.0e-10,
    minimum_element_reduction: float = 5.0,
    minimum_speedup: float = 5.0,
) -> str:
    """Gate a stranded-bundle approximation against an explicit reference."""

    try:
        result = _homogenized_bundle_impedance_comparison_gate(
            rows,
            resistance_rtol=resistance_rtol,
            inductance_rtol=inductance_rtol,
            impedance_rtol=impedance_rtol,
            observable_rtol=observable_rtol,
            minimum_element_reduction=minimum_element_reduction,
            minimum_speedup=minimum_speedup,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "homogenized_bundle_impedance_comparison_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def loss_temperature_coupling_gate(
    magnetic_rows: list[dict],
    thermal_rows: list[dict],
    loss_to_heat_scale: float,
    coupling_rtol: float = 2.0e-5,
    decomposition_rtol: float = 1.0e-12,
    minimum_power_coverage: float = 0.90,
    initial_temperature_c: float = 20.0,
) -> str:
    """Gate an electromagnetic-loss to transient-temperature handoff."""

    try:
        result = _loss_temperature_coupling_gate(
            magnetic_rows,
            thermal_rows,
            loss_to_heat_scale=loss_to_heat_scale,
            coupling_rtol=coupling_rtol,
            decomposition_rtol=decomposition_rtol,
            minimum_power_coverage=minimum_power_coverage,
            initial_temperature_c=initial_temperature_c,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "loss_temperature_coupling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def linked_study_silent_noop_gate(
    summary: dict,
    maximum_noop_seconds: float = 1.0,
) -> str:
    """Verify a linked native run that returned without creating solver results."""

    try:
        result = _linked_study_silent_noop_gate(
            summary,
            maximum_noop_seconds=maximum_noop_seconds,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "linked_study_silent_noop_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def reciprocal_two_port_power_sweep_gate(
    rows: list[dict],
    adaptive_power_rows: list[dict],
    reference_impedance_ohm: float,
    reciprocity_atol: float = 5.0e-6,
    reflection_symmetry_atol: float = 2.0e-4,
    passivity_atol: float = 1.0e-10,
    power_balance_atol: float = 5.0e-8,
) -> str:
    """Gate complex two-port reciprocity, symmetry, passivity, and power closure."""

    try:
        result = _reciprocal_two_port_power_sweep_gate(
            rows,
            adaptive_power_rows,
            reference_impedance_ohm=reference_impedance_ohm,
            reciprocity_atol=reciprocity_atol,
            reflection_symmetry_atol=reflection_symmetry_atol,
            passivity_atol=passivity_atol,
            power_balance_atol=power_balance_atol,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "reciprocal_two_port_power_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def fem_bem_capstone_suite_gate(payload: dict) -> str:
    """Gate a ten-case first-order FEM/BEM reference capstone suite."""

    try:
        result = _fem_bem_capstone_suite_gate(payload)
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "fem_bem_capstone_suite_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def helmholtz_dual_formulation_axis_gate(summary_json: str) -> str:
    """Gate Helmholtz-coil axis symmetry, flatness, and formulation agreement."""

    try:
        result = _helmholtz_dual_formulation_axis_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "helmholtz_dual_formulation_axis_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def hartmann_profile_gate(summary_json: str) -> str:
    """Gate a Hartmann-number sweep against an independent channel profile."""

    try:
        result = _hartmann_profile_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "hartmann_profile_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def lossy_dielectric_complex_power_refinement_gate(summary_json: str) -> str:
    """Gate lossy-dielectric constitutive, energy, complex-power, and mesh closure."""

    try:
        result = _lossy_power_refinement_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "lossy_dielectric_complex_power_refinement_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def heterogeneous_current_flow_p1_reintegration_gate(summary_json: str) -> str:
    """Gate heterogeneous current-flow P1 reintegration and sign covariance."""

    try:
        result = _heterogeneous_current_flow_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "heterogeneous_current_flow_p1_reintegration_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def thermal_robin_boundary_balance_gate(summary_json: str) -> str:
    """Gate signed Robin heat balance, mesh plateau, replay, and reflection."""

    try:
        result = _thermal_robin_boundary_balance_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "thermal_robin_boundary_balance_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def hysteresis_minor_loop_replay_gate(summary_json: str) -> str:
    """Gate history, knot normalization, signed loss, and exact loop replay."""

    try:
        result = _hysteresis_minor_loop_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "hysteresis_minor_loop_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def heterogeneous_part_mesh_replay_gate(summary_json: str) -> str:
    """Diagnose deterministic heterogeneous part-mesh replay drift."""

    try:
        result = _heterogeneous_part_mesh_replay_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "heterogeneous_part_mesh_replay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def two_terminal_dc_conduction_power_gate(summary_json: str) -> str:
    """Gate current closure, Joule power, adaptive convergence, and replay."""

    try:
        result = _two_terminal_dc_conduction_power_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "two_terminal_dc_conduction_power_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def rwg_hcurl_trace_consistency_gate(summary_json: str) -> str:
    """Gate RWG/HCurl trace topology, de Rham closure, and reference matrices."""

    try:
        result = _rwg_hcurl_trace_consistency_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "rwg_hcurl_trace_consistency_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def transient_coupled_coil_response_gate(
    times_s: list[float],
    primary_current_a: list[float],
    secondary_current_a: list[float],
    secondary_resistance_ohm: float,
    secondary_turns: float,
    maximum_relative_residual: float = 1.0e-3,
) -> str:
    """Gate a passive shorted-secondary transient induced-current history."""

    try:
        result = _transient_coupled_coil_response_gate(
            times_s,
            primary_current_a,
            secondary_current_a,
            secondary_resistance_ohm=secondary_resistance_ohm,
            secondary_turns=secondary_turns,
            maximum_relative_residual=maximum_relative_residual,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "transient_coupled_coil_response_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def source_off_linear_relaxation_gate(
    summary: dict,
    max_decay_ratio_relative_span: float = 1.0e-3,
    max_field_current_scale_relative_span: float = 1.0e-3,
) -> str:
    """Gate a linear source-off RL relaxation using total current and field decay."""

    try:
        result = _source_off_linear_relaxation_gate(
            summary,
            max_decay_ratio_relative_span=max_decay_ratio_relative_span,
            max_field_current_scale_relative_span=max_field_current_scale_relative_span,
        )
    except (KeyError, TypeError, ValueError) as exc:
        result = {
            "policy": "linear_source_off_total_current_field_decay_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def nonlinear_bh_piecewise_material_gate(
    summary: dict,
    maximum_relative_identity_error: float = 1.0e-8,
) -> str:
    """Gate secant and left-interval differential permeability from B-H rows."""

    try:
        result = _nonlinear_bh_piecewise_material_gate(
            summary,
            maximum_relative_identity_error=maximum_relative_identity_error,
        )
    except (KeyError, TypeError, ValueError) as exc:
        result = {
            "policy": "piecewise_bh_secant_and_left_interval_differential_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def skin_effect_adaptive_energy_loss_gate(
    frequency_hz: float, current: dict, voltage: dict, impedance: dict, power: dict,
    flux_linkage: dict, total_energy_j: float, total_loss_w: float, adaptive_rows: list[dict],
) -> str:
    """Gate current-port identities and adaptive skin-effect loss convergence."""
    try:
        result = _skin_effect_adaptive_energy_loss_gate(
            frequency_hz=frequency_hz, current=current, voltage=voltage, impedance=impedance,
            power=power, flux_linkage=flux_linkage, total_energy_j=total_energy_j,
            total_loss_w=total_loss_w, adaptive_rows=adaptive_rows,
        )
    except (TypeError, ValueError, KeyError) as exc:
        result = {"policy": "skin_effect_adaptive_energy_loss_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def global_local_optimization_replay_gate(summary_json: str) -> str:
    """Gate a stochastic global-search to derivative-checked local-polish replay."""
    try:
        result = _global_local_optimization_replay_gate(summary_json)
    except (TypeError, ValueError, KeyError) as exc:
        result = {"policy":"global_local_optimization_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def alternate_eddy_loss_formulation_gate(summary_json: str) -> str:
    """Gate volume-resolved and surface-impedance losses as non-additive alternatives."""
    try:
        result = _alternate_eddy_loss_formulation_gate(json.loads(summary_json))
    except (TypeError, ValueError, KeyError, json.JSONDecodeError) as exc:
        result = {
            "policy": "alternate_eddy_loss_formulation_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def adjoint_gradient_scaling_gate(
    rows_json: str,
    max_gradient_relative_error: float = 1.0e-6,
    max_forward_affine_residual: float = 1.0e-10,
    min_final_objective_ratio: float = 1.0,
) -> str:
    """Gate reverse-mode solve scaling, FD agreement and ascent direction."""

    try:
        rows = json.loads(rows_json)
        result = _adjoint_gradient_scaling_gate(
            rows,
            max_gradient_relative_error=max_gradient_relative_error,
            max_forward_affine_residual=max_forward_affine_residual,
            min_final_objective_ratio=min_final_objective_ratio,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "adjoint_gradient_scaling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


register_status_tool(
    mcp,
    server_name='mcp-server-radia-ngsolve',
    description='Radia + NGSolve: Kelvin / sparsesolv / CLN / PEEC / analytical formulas / lint',
    subpackage='radia_mcp.radia_ngsolve',
    related_servers=["fem", "bem", "matrix-solvers"],
    optional_deps=["radia", "ngsolve"],
)


def main():
    """Entry point for mcp-server-ngsolve console script."""
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        from radia_mcp.common.utf8_stdout import use_utf8_stdout
        use_utf8_stdout()
        _selftest()
    else:
        mcp.run(transport="stdio")


if __name__ == '__main__':
    main()
