import numpy as np
import pytest

from radia.accelerator_abe_topopt import (
    ExactSectionEvaluation,
    bin_element_fill_to_interface_height,
    blended_interface_displacement,
    compose_specification_fill_response,
    contract_hdiv_element_fill_response,
    measured_element_fill_patterns,
    optimize_abe_section_contour,
    solve_abe_element_fill_plan,
)


def test_measured_fill_patterns_require_explicit_local_dof_compatibility():
    blocks = (np.array([0, 1]), np.array([2, 3]), np.array([4, 5]))
    state = np.array([1.0, 2.0, 3.0, 4.0, 0.0, 0.0])
    centroids = np.array([[0.0, 0.0, 0.0],
                          [1.0, 0.0, 0.0],
                          [1.1, 0.0, 0.0]])
    active = np.array([True, True, False])
    with pytest.raises(ValueError, match="pattern_transfer"):
        measured_element_fill_patterns(
            state, blocks, centroids, active, [1, 2])

    elements, patterns = measured_element_fill_patterns(
        state, blocks, centroids, active, [1, 2],
        assume_compatible_local_dofs=True)
    assert elements.tolist() == [1, 2]
    np.testing.assert_allclose(patterns[0], [3.0, 4.0])
    np.testing.assert_allclose(patterns[1], [3.0, 4.0])

    _, transferred = measured_element_fill_patterns(
        state, blocks, centroids, active, [2],
        pattern_transfer=lambda source, target, values: -values)
    np.testing.assert_allclose(transferred[0], [-3.0, -4.0])


def test_hdiv_rows_contract_to_one_column_per_element():
    rows = np.array([[1.0, 2.0, 3.0, 4.0],
                     [4.0, 3.0, 2.0, 1.0]])
    blocks = (np.array([0, 1]), np.array([2, 3]))
    response = contract_hdiv_element_fill_response(
        rows, blocks, [0, 1], ([2.0, -1.0], [0.5, 2.0]))
    expected = np.column_stack((
        rows[:, :2] @ np.array([2.0, -1.0]),
        rows[:, 2:] @ np.array([0.5, 2.0]),
    ))
    np.testing.assert_allclose(response, expected)
    assert response.flags.c_contiguous

    jacobian = np.array([[1.0, 2.0], [-1.0, 0.5], [3.0, -2.0]])
    composed = compose_specification_fill_response(
        jacobian, response, specification_rows=[0, 2])
    np.testing.assert_allclose(composed, jacobian[[0, 2]] @ response)
    assert composed.flags.c_contiguous


def test_abe_element_fill_plan_uses_signed_capacity_and_reports_volume():
    response = np.array([[1.0, 1.0], [1.0, -1.0]])
    plan = solve_abe_element_fill_plan(
        response, [0.75, 1.25], material_active=[False, True],
        element_volumes=[2.0, 3.0], element_ids=[7, 9],
        field_response=np.array([[2.0, 3.0]]),
        residual_rms=1.0e-11, method="dense", max_iterations=64)
    np.testing.assert_allclose(plan.fill_step, [1.0, -0.25], atol=2.0e-10)
    np.testing.assert_allclose(plan.delivered_specification,
                               [0.75, 1.25], atol=3.0e-10)
    np.testing.assert_allclose(plan.implied_field_difference,
                               [1.25], atol=3.0e-10)
    assert plan.element_ids.tolist() == [7, 9]
    assert plan.gross_material_volume == pytest.approx(2.75, abs=1.0e-9)
    assert plan.net_material_volume == pytest.approx(1.25, abs=1.0e-9)
    assert plan.numerical_rank == 2
    assert plan.bounded_solution.converged


def test_binned_fill_to_height_is_volume_conservative_and_interpolable():
    field = bin_element_fill_to_interface_height(
        [0.5, -0.25], [2.0, 4.0], [0.5, 1.5], [0.5, 0.5],
        [0.0, 1.0, 2.0], [0.0, 1.0], [[2.0], [4.0]])
    np.testing.assert_allclose(field.signed_volume[:, 0], [1.0, -1.0])
    np.testing.assert_allclose(field.height_change[:, 0], [-0.5, 0.25])
    assert field.element_count[:, 0].tolist() == [1, 1]
    np.testing.assert_allclose(
        field.sample([0.5, 1.5], [0.5, 0.5]), [-0.5, 0.25])
    assert -float(np.sum(field.height_change * field.bin_areas)) == \
        pytest.approx(float(np.sum(field.signed_volume)))


def test_blended_interface_displacement_preserves_fixed_surfaces():
    normal = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 4.0])
    displacement = blended_interface_displacement(
        normal, 1.0, 0.0, 3.0, 2.0)
    np.testing.assert_allclose(displacement, [0.0, 1.0, 2.0, 1.0, 0.0, 0.0])


def test_abe_contour_outer_loop_relinearizes_after_complete_solve():
    def realize(fill):
        return fill.copy()

    def evaluate(fill):
        return ExactSectionEvaluation(np.array([0.8 * fill[0]]), fill.copy())

    result = optimize_abe_section_contour(
        np.array([[1.0]]), target_specification=[0.6],
        response_band=[1.0e-5],
        initial_evaluation=ExactSectionEvaluation(np.array([0.0])),
        material_active=[False], element_volumes=[2.0],
        realize_fill=realize, evaluate_exact=evaluate,
        maximum_iterations=3, inner_residual_fraction=1.0e-6,
        relinearize=lambda exact, geometry: np.array([[0.8]]),
        solve_options={"method": "dense", "max_iterations": 64})
    assert result.converged
    assert result.stop_reason == "target_band_reached"
    assert len(result.history) == 2
    np.testing.assert_allclose(result.accumulated_fill, [0.75], atol=2.0e-6)
    np.testing.assert_allclose(result.final_evaluation.specification,
                               [0.6], atol=2.0e-6)
    assert result.final_max_band_ratio <= 1.0


def test_abe_contour_guard_rejects_invalid_exact_states_and_backtracks():
    result = optimize_abe_section_contour(
        np.array([[1.0]]), target_specification=[0.8],
        response_band=[1.0e-3],
        initial_evaluation=ExactSectionEvaluation(np.array([0.0])),
        material_active=[False], element_volumes=[1.0],
        realize_fill=lambda fill: fill.copy(),
        evaluate_exact=lambda fill: ExactSectionEvaluation(fill.copy()),
        maximum_iterations=2, inner_residual_fraction=1.0e-6,
        backtracking_scales=(1.0, 0.5),
        exact_guard=lambda exact: exact.specification[0] <= 0.5,
        solve_options={"method": "dense", "max_iterations": 64})
    assert not result.converged
    assert result.stop_reason == "exact_backtracking_exhausted"
    assert len(result.history) == 1
    assert result.history[0].backtracking_scale == 0.5
    np.testing.assert_allclose(result.accumulated_fill, [0.4], atol=2.0e-6)
