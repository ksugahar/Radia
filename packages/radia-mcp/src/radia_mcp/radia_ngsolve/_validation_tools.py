# Validation catalog tool wrappers extracted from server.py.

import json

from . import server as _srv
from .server import *  # noqa: F401, F403
globals().update({k: v for k, v in vars(_srv).items() if k.startswith('_') and not k.startswith('__')})

@_validation.tool()
def ellipsoid_delta_benchmark(
    order: int = 2, curve_order: int = 2, iron_maxh: float = 0.45,
    mu_r: float = 20.0, integration_order: int = 12,
    air_divisions: int = 5, assembly_bonus: int = 0,
) -> dict:
    """Compare reduced A and mixed Omega on a synthetic curved Kelvin ellipsoid.

    Returns the magnetic moment gap, field-energy defect, exact CAD moment
    errors, and timings. Hold curve_order and iron_maxh fixed for p convergence.
    The result is diagnostic and does not certify an exact-geometry bracket.
    """
    import ngsolve as ng
    from radia.dual_fem_delta import ellipsoid_delta_benchmark as benchmark

    with ng.TaskManager():
        return benchmark(order=order, curve_order=curve_order, iron_maxh=iron_maxh,
                         mu_r=mu_r, integration_order=integration_order,
                         air_divisions=air_divisions, assembly_bonus=assembly_bonus)


@_validation.tool()
def ellipsoid_delta_series(
    orders: list[int] | None = None, curve_order: int = 2, shape: str = "triaxial",
    iron_maxh: float = 0.45, mu_r: float = 20.0, integration_order: int = 12,
    air_divisions: int = 5, assembly_bonus: int = 0, vol_path: str = "",
) -> dict:
    """Compare A, mixed Omega, 2:1 and midpoint at several p on one verified mesh.

    Triaxial has a built-in OCC mesh; prolate requires a local affine Cubit
    .vol. An affine file must contain the full iron/air/Kelvin benchmark with
    the selected analytic dimensions and radius 2.75 m. Curve order must not
    exceed any requested p.
    The operation rejects folded elements and mismatched curved Kelvin seams
    before solving. Reported CAD errors apply to the source-conjugate moment.
    """
    import ngsolve as ng
    from radia.dual_fem_delta import ellipsoid_delta_series as series

    with ng.TaskManager():
        return series(orders=(2, 3, 4) if orders is None else orders,
                      curve_order=curve_order, shape=shape, iron_maxh=iron_maxh,
                      mu_r=mu_r, integration_order=integration_order,
                      air_divisions=air_divisions, assembly_bonus=assembly_bonus,
                      vol_path=vol_path)


@_validation.tool()
def dual_fem_reference_estimate(
    levels: list[dict], coarse: dict | None = None, convergence_fraction: float = 0.1,
) -> dict:
    """Freeze an estimated true value from converged A/Omega dual bounds.

    levels: coarse-to-fine list of {"q_lower": A bound, "q_upper": Omega bound}
    computed on one geometry with energy-functional (kappa-corrected) values.
    The finest level gives the reference midpoint, with its half-width as
    uncertainty, the last-step drift and a nesting check.
    The optional coarse level is judged only against that frozen reference:
    errors of A, Omega, midpoint and 2:1, the ratio r, and which blend is closer
    (midpoint iff r < 1.4; inconclusive inside the reference width or while the
    reference is not converged, with the would-be answer kept as
    provisional_closer). convergence_fraction must be finite and >= 0;
    non-finite derived values raise instead of returning a verdict.
    r depends on the case, mesh and order, so no fixed weight is assumed.
    Bounds certify the mapped mesh domain. Planar-faced bodies on affine
    meshes have no geometry error; curved CAD bodies do.
    """
    from radia.dual_fem_delta import dual_bracket_reference

    return dual_bracket_reference(levels, coarse=coarse,
                                  convergence_fraction=convergence_fraction)


@_validation.tool()
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


@_validation.tool()
def physics_result_preflight_gate(summary_json: str) -> str:
    """Gate physics namespace, selection, solution, and license metadata before result evaluation."""
    try: result=_physics_result_preflight_gate(summary_json)
    except (json.JSONDecodeError,TypeError,ValueError) as exc: result={"policy":"physics_result_preflight_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def cq_response_reality_gate(
    summary_json: str,
    residual_tolerance: float = 1.0e-10,
    imaginary_tolerance: float = 1.0e-10,
) -> str:
    """Gate a coupled CQ solve, including its real time-domain reconstruction."""
    return json.dumps(_cq_response_reality_gate(
        json.loads(summary_json), residual_tolerance=residual_tolerance,
        imaginary_tolerance=imaginary_tolerance), indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_magnetic_refinement_energy_gate(
    summary_json: str,
    max_refinement_growth_factor: float = 1.05,
    max_finest_pair_relative_change: float = 0.05,
    max_legendre_relative_residual: float = 1.0e-8,
    min_refinement_levels: int = 3,
) -> str:
    """Gate a nonlinear magnetic h ladder, material bounds, and field/energy identity."""

    try:
        result = _nonlinear_magnetic_refinement_energy_gate(
            json.loads(summary_json),
            max_refinement_growth_factor=max_refinement_growth_factor,
            max_finest_pair_relative_change=max_finest_pair_relative_change,
            max_legendre_relative_residual=max_legendre_relative_residual,
            min_refinement_levels=min_refinement_levels,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_magnetic_refinement_energy_gate_v3",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_magnetic_field_energy_parity_gate(
    summary_json: str,
    max_average_vector_relative_difference: float = 0.05,
    max_rms_magnitude_relative_difference: float = 0.05,
    max_energy_relative_difference: float = 0.05,
    max_coenergy_relative_difference: float = 0.05,
) -> str:
    """Compare nonlinear field/energy only after full physical identity matches."""

    try:
        result = _nonlinear_magnetic_field_energy_parity_gate(
            json.loads(summary_json),
            max_average_vector_relative_difference=(
                max_average_vector_relative_difference
            ),
            max_rms_magnitude_relative_difference=(
                max_rms_magnitude_relative_difference
            ),
            max_energy_relative_difference=max_energy_relative_difference,
            max_coenergy_relative_difference=max_coenergy_relative_difference,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_magnetic_field_energy_parity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_field_energy_identity_gate_v5(summary_json: str) -> str:
    """Gate canonical geometry/frame/refinement identity and nonlinear energy derivatives."""

    try:
        result = _nonlinear_field_energy_identity_gate_v5(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_field_energy_identity_gate_v5",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_field_energy_artifact_contract_gate_v6(summary_json: str) -> str:
    """Gate completed result schema, units, finite samples, and response digest."""

    try:
        result = _nonlinear_field_energy_artifact_contract_gate_v6(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_field_energy_artifact_contract_gate_v6",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_field_energy_lineage_gate_v7(summary_json: str) -> str:
    """Gate solver-output lineage, recomputable run identity, and cross-lane mesh identity."""

    try:
        result = _nonlinear_field_energy_lineage_gate_v7(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_field_energy_lineage_gate_v7",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_field_energy_physical_admissibility_gate_v8(summary_json: str) -> str:
    """Gate nonnegative energy, positive tangent permeability, and Legendre identity."""

    try:
        result = _nonlinear_field_energy_physical_admissibility_gate_v8(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_field_energy_physical_admissibility_gate_v8",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_field_energy_observable_comparison_gate_v9(summary_json: str) -> str:
    """Compare observables only after identity, units, and sample identities match."""

    try:
        result = _nonlinear_field_energy_observable_comparison_gate_v9(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_field_energy_observable_comparison_gate_v9",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_vector_observable_comparison_gate_v10(summary_json: str) -> str:
    """Compare field, force, and torque vectors after frame and unit checks."""

    try:
        result = _nonlinear_vector_observable_comparison_gate_v10(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_vector_observable_comparison_gate_v10",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_live_execution_identity_gate_v11(summary_json: str) -> str:
    """Gate completed live execution, physical identity, artifact digests, and cleanup evidence."""

    try:
        result = _nonlinear_live_execution_identity_gate_v11(json.loads(summary_json))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_live_execution_identity_gate_v11",
            "status": "invalid_input",
            "accepted": False,
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def build_constitutive_comparison_candidate(
    bh_table_json: str,
    h_values_json: str,
    constitutive_interpolation: str = "monotone_pchip",
) -> str:
    """Build a source-comparison B-H response without changing the Radia solver."""

    try:
        result = _build_constitutive_comparison_candidate(
            json.loads(bh_table_json),
            json.loads(h_values_json),
            constitutive_interpolation=constitutive_interpolation,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "schema": "radia.nonlinear-constitutive-response.v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_constitutive_response_parity_gate(
    summary_json: str,
    max_response_relative_difference: float = 1.0e-6,
    max_integrability_relative_residual: float = 1.0e-8,
    min_response_samples: int = 5,
) -> str:
    """Compare realized B-H response, tangent, energy, and coenergy on one SI grid."""

    try:
        result = _nonlinear_constitutive_response_parity_gate(
            json.loads(summary_json),
            max_response_relative_difference=max_response_relative_difference,
            max_integrability_relative_residual=max_integrability_relative_residual,
            min_response_samples=min_response_samples,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_constitutive_response_parity_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_constitutive_point_sample_gate(
    summary_json: str,
    max_response_relative_difference: float = 0.05,
    max_direction_sine: float = 0.02,
    max_vector_magnitude_relative_residual: float = 1.0e-9,
    min_response_samples: int = 5,
) -> str:
    """Admit only unsmoothed element-local B/H samples as constitutive evidence."""

    try:
        result = _nonlinear_constitutive_point_sample_gate(
            json.loads(summary_json),
            max_response_relative_difference=max_response_relative_difference,
            max_direction_sine=max_direction_sine,
            max_vector_magnitude_relative_residual=(
                max_vector_magnitude_relative_residual
            ),
            min_response_samples=min_response_samples,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_constitutive_point_sample_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def controlled_uniform_field_constitutive_sweep_gate(
    summary_json: str,
    max_response_relative_difference: float = 0.01,
    min_response_samples: int = 5,
) -> str:
    """Compare one identity-bound homogeneous B-H sweep with a candidate response."""

    try:
        result = _controlled_uniform_field_constitutive_sweep_gate(
            json.loads(summary_json),
            max_response_relative_difference=max_response_relative_difference,
            min_response_samples=min_response_samples,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "controlled_uniform_field_constitutive_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def nonlinear_magnetic_spatial_evidence_gate(
    summary_json: str,
    max_average_vector_relative_difference: float = 0.07,
    max_rms_magnitude_relative_difference: float = 0.10,
    min_tensor_gauss_samples: int = 27,
) -> str:
    """Gate nonlinear magnetic volume evidence and material/response order compatibility."""

    try:
        result = _nonlinear_magnetic_spatial_evidence_gate(
            json.loads(summary_json),
            max_average_vector_relative_difference=max_average_vector_relative_difference,
            max_rms_magnitude_relative_difference=max_rms_magnitude_relative_difference,
            min_tensor_gauss_samples=min_tensor_gauss_samples,
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        result = {
            "policy": "nonlinear_magnetic_spatial_evidence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def cyclic_terminal_source_sweep_gate(summary_json: str) -> str:
    """Gate cyclic terminal charges without assuming formulations are identical."""
    try: result=_cyclic_terminal_source_sweep_gate(json.loads(summary_json))
    except (TypeError,ValueError) as exc: result={"policy":"cyclic_terminal_source_sweep_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result,indent=2,sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def nonlinear_actuator_saturation_knee_gate(summary_json: str) -> str:
    """Gate an axisymmetric nonlinear actuator by a shared L/F saturation knee."""
    try:
        result = _nonlinear_actuator_saturation_knee_gate(json.loads(summary_json))
    except (TypeError, ValueError, KeyError) as exc:
        result = {"policy": "nonlinear_actuator_saturation_knee_gate_v1", "status": "invalid_input", "error": str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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
    try:
        payload = json.loads(summary_json)
    except (json.JSONDecodeError, TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        v44_checks = _validate_force_v44_identity(payload.get("artifact_identity"))
        if v44_checks:
            result.setdefault("checks", {}).update(v44_checks)
            result["force_v44_identity_checks"] = v44_checks
            if not all(v44_checks.values()):
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def rotational_eddy_brake_energy_gate(summary_json: str) -> str:
    """Gate free rotational braking with angular impulse and field energy."""
    try:
        result = _rotational_eddy_brake_energy_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
        result = {
            "policy": "rotational_eddy_brake_energy_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    try:
        payload = json.loads(summary_json)
    except (json.JSONDecodeError, TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        v46_checks = _validate_solver_state_v46_identity(payload)
        if v46_checks:
            result.setdefault("checks", {}).update(v46_checks["checks"])
            result["solver_state_v46_identity"] = v46_checks
            if v46_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v47_checks = _validate_eddy_brake_v47_identity(payload)
        if v47_checks:
            result.setdefault("checks", {}).update(v47_checks["checks"])
            result["comsol_v47_identity"] = v47_checks
            if v47_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v48_checks = _validate_eddy_brake_v48_identity(payload)
        if v48_checks:
            result.setdefault("checks", {}).update(v48_checks["checks"])
            result["transform_normalization_v48"] = v48_checks
            if v48_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v49_checks = _validate_eddy_brake_v49_identity(payload)
        if v49_checks:
            result.setdefault("checks", {}).update(v49_checks["checks"])
            result["solver_state_identity_v49"] = v49_checks
            if v49_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v50_checks = _validate_eddy_brake_v50_identity(payload)
        if v50_checks:
            result.setdefault("checks", {}).update(v50_checks["checks"])
            result["frequency_contact_identity_v50"] = v50_checks
            if v50_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v51_checks = _validate_eddy_brake_v51_identity(payload)
        if v51_checks:
            result.setdefault("checks", {}).update(v51_checks["checks"])
            result["modal_continuation_identity_v51"] = v51_checks
            if v51_checks["status"] != "ok":
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def grounded_sphere_capacitance_convergence_gate(
    summary_json: str,
    max_final_analytic_relative_error: float = 5.0e-4,
    max_corrected_energy_relative_error: float = 1.0e-5,
    max_replay_relative_error: float = 1.0e-13,
    max_sign_covariance_relative_error: float = 1.0e-13,
) -> str:
    """Gate grounded-sphere image-series convergence and mixed-boundary energy."""

    try:
        result = _grounded_sphere_capacitance_convergence_gate(
            json.loads(summary_json),
            max_final_analytic_relative_error=max_final_analytic_relative_error,
            max_corrected_energy_relative_error=max_corrected_energy_relative_error,
            max_replay_relative_error=max_replay_relative_error,
            max_sign_covariance_relative_error=max_sign_covariance_relative_error,
        )
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "grounded_sphere_capacitance_convergence_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def two_conductor_capacitance_matrix_gate(summary_json: str) -> str:
    """Gate reciprocal Maxwell and mutual capacitance matrix representations."""
    try:
        result = _two_conductor_capacitance_matrix_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {"policy":"two_conductor_capacitance_matrix_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def force_coenergy_displacement_gate(
    positions_m: list[float],
    coenergy_j: list[float],
    forces_along_displacement_n: list[float],
    energy_kind: str = "constant_current_coenergy",
    max_central_relative_error: float = 0.02,
    min_sample_count: int = 5,
    artifact_identity: dict | None = None,
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
            artifact_identity=artifact_identity,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "force_coenergy_displacement_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    if isinstance(artifact_identity, dict):
        v44_checks = _validate_axisymmetric_v44_identity(artifact_identity)
        if v44_checks:
            result.setdefault("checks", {}).update(v44_checks)
            result["axisymmetric_v44_identity_checks"] = v44_checks
            if not all(v44_checks.values()):
                result["status"] = "needs_attention"
        v46_checks = _validate_axisymmetric_v46_identity(artifact_identity)
        if v46_checks:
            result.setdefault("checks", {}).update(v46_checks)
            result["axisymmetric_v46_identity_checks"] = v46_checks
            if not all(v46_checks.values()):
                result["status"] = "needs_attention"
        v47_checks = _validate_coenergy_force_v47_identity(artifact_identity)
        if v47_checks:
            result.setdefault("checks", {}).update(v47_checks)
            result["magnetic_v47_identity_checks"] = v47_checks
            if not all(v47_checks.values()):
                result["status"] = "needs_attention"
        v48_checks = _validate_coenergy_force_v48_identity(artifact_identity)
        if v48_checks:
            result.setdefault("checks", {}).update(v48_checks)
            result["electromagnetic_v48_identity_checks"] = v48_checks
            if not all(v48_checks.values()):
                result["status"] = "needs_attention"
        v49_checks = _validate_coenergy_force_v49_identity(artifact_identity)
        if v49_checks:
            result.setdefault("checks", {}).update(v49_checks)
            result["electromagnetic_v49_identity_checks"] = v49_checks
            if not all(v49_checks.values()):
                result["status"] = "needs_attention"
        v50_checks = _validate_coenergy_force_v50_identity(artifact_identity)
        if v50_checks:
            result.setdefault("checks", {}).update(v50_checks)
            result["electromagnetic_v50_identity_checks"] = v50_checks
            if not all(v50_checks.values()):
                result["status"] = "needs_attention"
        v51_checks = _validate_coenergy_force_v51_identity(artifact_identity)
        if v51_checks:
            result.setdefault("checks", {}).update(v51_checks)
            result["electromagnetic_v51_identity_checks"] = v51_checks
            if not all(v51_checks.values()):
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def nonlinear_inductance_sweep_gate(summary_json: str) -> str:
    """Gate nonlinear apparent/incremental matrices, duality, and replay."""
    try:
        result = _nonlinear_inductance_sweep_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "nonlinear_inductance_sweep_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    try:
        payload = json.loads(summary_json)
    except (json.JSONDecodeError, TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        v44_checks = _validate_waveguide_emc_v44_identity(payload)
        if v44_checks:
            result.setdefault("checks", {}).update(v44_checks)
            result["waveguide_emc_v44_identity_checks"] = v44_checks
            if not all(v44_checks.values()):
                result["status"] = "needs_attention"
        v45_checks = _validate_network_v45_identity(payload)
        if v45_checks:
            result.setdefault("checks", {}).update(v45_checks)
            result["network_v45_identity_checks"] = v45_checks
            if not all(v45_checks.values()):
                result["status"] = "needs_attention"
        v46_checks = _validate_network_v46_identity(payload)
        if v46_checks:
            result.setdefault("checks", {}).update(v46_checks)
            result["network_v46_identity_checks"] = v46_checks
            if not all(v46_checks.values()):
                result["status"] = "needs_attention"
        v47_checks = _validate_network_v47_identity(payload)
        if v47_checks:
            result.setdefault("checks", {}).update(v47_checks)
            result["network_v47_identity_checks"] = v47_checks
            if not all(v47_checks.values()):
                result["status"] = "needs_attention"
        v48_checks = _validate_network_v48_identity(payload)
        if v48_checks:
            result.setdefault("checks", {}).update(v48_checks)
            result["network_v48_identity_checks"] = v48_checks
            if not all(v48_checks.values()):
                result["status"] = "needs_attention"
        v49_checks = _validate_network_v49_identity(payload)
        if v49_checks:
            result.setdefault("checks", {}).update(v49_checks)
            result["network_v49_identity_checks"] = v49_checks
            if not all(v49_checks.values()):
                result["status"] = "needs_attention"
        v50_checks = _validate_network_v50_identity(payload)
        if v50_checks:
            result.setdefault("checks", {}).update(v50_checks)
            result["network_v50_identity_checks"] = v50_checks
            if not all(v50_checks.values()):
                result["status"] = "needs_attention"
        v51_checks = _validate_network_v51_identity(payload)
        if v51_checks:
            result.setdefault("checks", {}).update(v51_checks)
            result["network_v51_identity_checks"] = v51_checks
            if not all(v51_checks.values()):
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def regularized_trace_inverse_path_gate(summary_json: str) -> str:
    """Gate a P1 trace Tikhonov path, L-curve, Morozov, and replay."""
    try:
        result = _regularized_trace_inverse_path_gate(json.loads(summary_json))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        result = {
            "policy": "regularized_trace_inverse_path_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    try:
        payload = json.loads(summary_json)
    except (json.JSONDecodeError, TypeError, ValueError):
        payload = None
    if isinstance(payload, dict):
        v44_checks = _validate_matlab_ml_rl_v44_identity(payload)
        if v44_checks:
            result.setdefault("checks", {}).update(v44_checks["checks"])
            result["matlab_ml_rl_v44_identity"] = v44_checks
            if v44_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v45_checks = _validate_matlab_ml_rl_v45_identity(payload)
        if v45_checks:
            result.setdefault("checks", {}).update(v45_checks["checks"])
            result["matlab_ml_rl_v45_identity"] = v45_checks
            if v45_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v46_checks = _validate_matlab_ml_rl_v46_identity(payload)
        if v46_checks:
            result.setdefault("checks", {}).update(v46_checks["checks"])
            result["matlab_ml_rl_v46_identity"] = v46_checks
            if v46_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v47_checks = _validate_matlab_ml_rl_v47_identity(payload)
        if v47_checks:
            result.setdefault("checks", {}).update(v47_checks["checks"])
            result["matlab_ml_rl_v47_identity"] = v47_checks
            if v47_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v48_checks = _validate_matlab_ml_rl_v48_identity(payload)
        if v48_checks:
            result.setdefault("checks", {}).update(v48_checks["checks"])
            result["matlab_ml_rl_v48_identity"] = v48_checks
            if v48_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v49_checks = _validate_matlab_ml_rl_v49_identity(payload)
        if v49_checks:
            result.setdefault("checks", {}).update(v49_checks["checks"])
            result["matlab_ml_rl_v49_identity"] = v49_checks
            if v49_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v50_checks = _validate_matlab_ml_rl_v50_identity(payload)
        if v50_checks:
            result.setdefault("checks", {}).update(v50_checks["checks"])
            result["matlab_ml_rl_v50_identity"] = v50_checks
            if v50_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v51_checks = _validate_matlab_ml_rl_v51_identity(payload)
        if v51_checks:
            result.setdefault("checks", {}).update(v51_checks["checks"])
            result["matlab_ml_rl_v51_identity"] = v51_checks
            if v51_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v52_checks = _validate_matlab_ml_rl_v52_identity(payload)
        if v52_checks:
            result.setdefault("checks", {}).update(v52_checks["checks"])
            result["matlab_ml_rl_v52_identity"] = v52_checks
            if v52_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v53_checks = _validate_matlab_ml_rl_v53_identity(payload)
        if v53_checks:
            result.setdefault("checks", {}).update(v53_checks["checks"])
            result["matlab_ml_rl_v53_identity"] = v53_checks
            if v53_checks["status"] != "ok":
                result["status"] = "needs_attention"
        v54_checks = _validate_matlab_ml_rl_v54_identity(payload)
        if v54_checks:
            result.setdefault("checks", {}).update(v54_checks["checks"])
            result["matlab_ml_rl_v54_identity"] = v54_checks
            if v54_checks["status"] != "ok":
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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
    if isinstance(payload, dict):
        v44_checks = _validate_motor_v44_identity(payload.get("artifact_identity"))
        if v44_checks:
            result.setdefault("checks", {}).update(v44_checks)
            result["motor_v44_identity_checks"] = v44_checks
            if not all(v44_checks.values()):
                result["status"] = "needs_attention"
        v46_checks = _validate_motor_v46_identity(payload.get("artifact_identity"))
        if v46_checks:
            result.setdefault("checks", {}).update(v46_checks)
            result["motor_v46_identity_checks"] = v46_checks
            if not all(v46_checks.values()):
                result["status"] = "needs_attention"
        v47_checks = _validate_motor_v47_identity(payload.get("artifact_identity"))
        if v47_checks:
            result.setdefault("checks", {}).update(v47_checks)
            result["motor_v47_identity_checks"] = v47_checks
            if not all(v47_checks.values()):
                result["status"] = "needs_attention"
        v48_checks = _validate_motor_v48_identity(payload.get("artifact_identity"))
        if v48_checks:
            result.setdefault("checks", {}).update(v48_checks)
            result["motor_v48_identity_checks"] = v48_checks
            if not all(v48_checks.values()):
                result["status"] = "needs_attention"
        v49_checks = _validate_motor_v49_identity(payload.get("artifact_identity"))
        if v49_checks:
            result.setdefault("checks", {}).update(v49_checks)
            result["motor_v49_identity_checks"] = v49_checks
            if not all(v49_checks.values()):
                result["status"] = "needs_attention"
        v50_checks = _validate_motor_v50_identity(payload.get("artifact_identity"))
        if v50_checks:
            result.setdefault("checks", {}).update(v50_checks)
            result["motor_v50_identity_checks"] = v50_checks
            if not all(v50_checks.values()):
                result["status"] = "needs_attention"
        v51_checks = _validate_motor_v51_identity(payload.get("artifact_identity"))
        if v51_checks:
            result.setdefault("checks", {}).update(v51_checks)
            result["motor_v51_identity_checks"] = v51_checks
            if not all(v51_checks.values()):
                result["status"] = "needs_attention"
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def finite_section_helmholtz_mixed_omega_gate(summary_json: str) -> str:
    """Gate a finite-section CoilBuilder source and mixed total/reduced Omega solve."""

    try:
        result = _finite_section_helmholtz_mixed_omega_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "finite_section_helmholtz_mixed_omega_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
def divergence_free_source_assembly_gate(summary_json: str) -> str:
    """Gate boundary-load acceleration using orientation and solved-field evidence."""

    try:
        result = _divergence_free_source_assembly_gate(json.loads(summary_json))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        result = {
            "policy": "divergence_free_source_assembly_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
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


@_validation.tool()
def nonlinear_bh_canonical_table_gate(contract: dict) -> str:
    """Gate explicit SI H,B data for a single-valued nonlinear material solve."""

    try:
        result = _nonlinear_bh_canonical_table_gate(contract)
    except (KeyError, TypeError, ValueError) as exc:
        result = {
            "policy": "explicit_si_h_b_single_valued_material_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
def global_local_optimization_replay_gate(summary_json: str) -> str:
    """Gate a stochastic global-search to derivative-checked local-polish replay."""
    try:
        result = _global_local_optimization_replay_gate(summary_json)
    except (TypeError, ValueError, KeyError) as exc:
        result = {"policy":"global_local_optimization_replay_gate_v1","status":"invalid_input","error":str(exc)}
    return json.dumps(result, indent=2, sort_keys=True)


@_validation.tool()
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


@_validation.tool()
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


