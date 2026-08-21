"""Abe/DUCAS material planning for smooth accelerator-magnet shape updates.

This module is the material half of the accelerator design cascade::

    section transfer specification
      -> material-realizable field/magnetization response
      -> ACA--QR--TSVD element fill fractions
      -> smooth interface-height change
      -> complete field solve and exact transfer-map acceptance

The design variable is one signed fill fraction per material cell.  It is
never an individual HDiv degree of freedom.  HDiv states enter only as a
measured, whole-element magnetization pattern that contracts the response
rows into one column per cell.  Air cells have capacity ``[0, 1]`` and iron
cells ``[-1, 0]`` relative to the reference geometry.

The outer driver is deliberately callback based.  Radia owns the material
inverse and acceptance bookkeeping; NGSolve owns meshes, transformations,
deformation, full reassembly, orbit tracking, and field evaluation.  The
geometry callback must realize an *absolute accumulated* fill vector from the
reference geometry so rejected/backtracked trials do not compound mesh
roundoff.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from radia.stream_function import (
    AbeBoundedCurrentPotentialSolution,
    solve_abe_bounded_current_potential,
)


def _finite_array(value, *, name, ndim=None):
    result = np.asarray(value, dtype=float)
    if ndim is not None and result.ndim != int(ndim):
        raise ValueError(f"{name} must have dimension {ndim}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _element_indices(value, element_count, *, name):
    supplied = np.asarray(value)
    if supplied.dtype == bool:
        supplied = supplied.reshape(-1)
        if supplied.shape != (int(element_count),):
            raise ValueError(f"{name} boolean mask has the wrong length")
        return np.flatnonzero(supplied).astype(np.int64)
    numeric = supplied.reshape(-1)
    if numeric.size == 0:
        raise ValueError(f"{name} must select at least one element")
    if not np.issubdtype(numeric.dtype, np.integer):
        as_float = np.asarray(numeric, dtype=float)
        if np.any(as_float != np.floor(as_float)):
            raise ValueError(f"{name} must contain integer element indices")
        numeric = as_float.astype(np.int64)
    result = np.asarray(numeric, dtype=np.int64)
    if (np.any(result < 0) or np.any(result >= int(element_count))
            or np.unique(result).size != result.size):
        raise ValueError(f"{name} contains invalid or duplicate elements")
    return result


def measured_element_fill_patterns(
        reference_state, element_dof_blocks, element_centroids,
        material_active, design_elements, *, pattern_source_elements=None,
        pattern_transfer=None, assume_compatible_local_dofs=False):
    """Return one measured full-material HDiv pattern per design element.

    Existing material uses its own solved HDiv block.  An inactive cell uses
    the nearest active source cell with the same local block size.  Copying a
    local HDiv block between cells is not generally orientation safe, so the
    caller must either provide ``pattern_transfer(source, target, values)`` or
    explicitly assert ``assume_compatible_local_dofs=True`` for a verified
    structured discontinuous mesh.  This function never reconstructs NGSolve
    basis functions or Piola transformations.
    """
    state = _finite_array(reference_state, name="reference_state").reshape(-1)
    blocks = tuple(np.asarray(block, dtype=np.int64).reshape(-1)
                   for block in element_dof_blocks)
    element_count = len(blocks)
    if element_count == 0 or any(block.size == 0 for block in blocks):
        raise ValueError("element_dof_blocks must contain non-empty blocks")
    if any(np.any(block < 0) or np.any(block >= state.size) for block in blocks):
        raise ValueError("element_dof_blocks contain out-of-range degrees of freedom")
    centroids = _finite_array(
        element_centroids, name="element_centroids", ndim=2)
    if centroids.shape != (element_count, 3):
        raise ValueError("element_centroids must have shape (ne, 3)")
    active = np.asarray(material_active, dtype=bool).reshape(-1)
    if active.shape != (element_count,):
        raise ValueError("material_active must have one value per element")
    design = _element_indices(
        design_elements, element_count, name="design_elements")
    if pattern_source_elements is None:
        source_mask = active.copy()
    else:
        source_mask = np.zeros(element_count, dtype=bool)
        source_mask[_element_indices(
            pattern_source_elements, element_count,
            name="pattern_source_elements")] = True
        source_mask &= active
    if not np.any(source_mask):
        raise ValueError("no active pattern-source element is available")
    if (np.any(~active[design]) and pattern_transfer is None
            and not bool(assume_compatible_local_dofs)):
        raise ValueError(
            "inactive element patterns require pattern_transfer or an explicit "
            "assume_compatible_local_dofs=True assertion")

    by_size = {}
    for size in {int(blocks[element].size) for element in design}:
        candidates = np.asarray([
            element for element in np.flatnonzero(source_mask)
            if blocks[element].size == size], dtype=np.int64)
        if candidates.size:
            by_size[size] = candidates
    patterns = []
    for target in design:
        target = int(target)
        if active[target]:
            pattern = state[blocks[target]].copy()
        else:
            candidates = by_size.get(int(blocks[target].size))
            if candidates is None or candidates.size == 0:
                raise ValueError(
                    "no active source has a compatible local HDiv block size")
            distances = np.linalg.norm(
                centroids[candidates] - centroids[target], axis=1)
            source = int(candidates[int(np.argmin(distances))])
            source_pattern = state[blocks[source]].copy()
            if pattern_transfer is None:
                pattern = source_pattern
            else:
                pattern = _finite_array(
                    pattern_transfer(source, target, source_pattern),
                    name="transferred element pattern").reshape(-1)
            if pattern.shape != (blocks[target].size,):
                raise ValueError(
                    "transferred element pattern has the wrong local size")
        patterns.append(np.ascontiguousarray(pattern))
    return design, tuple(patterns)


def contract_hdiv_element_fill_response(
        response_rows, element_dof_blocks, design_elements, fill_patterns):
    """Contract HDiv response rows into one fill-fraction column per cell."""
    rows = _finite_array(response_rows, name="response_rows", ndim=2)
    blocks = tuple(np.asarray(block, dtype=np.int64).reshape(-1)
                   for block in element_dof_blocks)
    design = _element_indices(
        design_elements, len(blocks), name="design_elements")
    patterns = tuple(fill_patterns)
    if len(patterns) != design.size:
        raise ValueError("fill_patterns must have one pattern per design element")
    if any(np.any(block < 0) or np.any(block >= rows.shape[1]) for block in blocks):
        raise ValueError("element_dof_blocks are incompatible with response_rows")
    result = np.empty((rows.shape[0], design.size), dtype=float)
    for column, (element, supplied) in enumerate(zip(design, patterns)):
        pattern = _finite_array(
            supplied, name="fill pattern").reshape(-1)
        dofs = blocks[int(element)]
        if pattern.shape != (dofs.size,):
            raise ValueError("a fill pattern does not match its element block")
        result[:, column] = rows[:, dofs] @ pattern
    return np.ascontiguousarray(result)


def compose_specification_fill_response(
        specification_field_jacobian, element_field_response, *,
        specification_rows=None):
    """Compose field-to-spec AD with whole-element fill response columns.

    This is the explicit Earlytimes-to-Abe junction.  The first matrix is the
    analytic/AD derivative of the physical section or Taylor-map
    specification with respect to the declared field coordinates.  The
    second matrix is the HDiv-MMM response of one full-material pattern per
    design cell.  No design finite difference or particle tracking is done in
    this contraction.
    """
    jacobian = _finite_array(
        specification_field_jacobian,
        name="specification_field_jacobian", ndim=2)
    field = _finite_array(
        element_field_response, name="element_field_response", ndim=2)
    if jacobian.shape[1] != field.shape[0]:
        raise ValueError(
            "specification Jacobian and field response dimensions do not match")
    if specification_rows is not None:
        rows = _element_indices(
            specification_rows, jacobian.shape[0],
            name="specification_rows")
        jacobian = jacobian[rows]
    return np.ascontiguousarray(jacobian @ field)


@dataclass(frozen=True)
class AbeElementFillPlan:
    """One bounded ACA--QR--TSVD element-fill proposal."""

    element_ids: np.ndarray
    fill_step: np.ndarray
    absolute_fill: np.ndarray
    lower_capacity: np.ndarray
    upper_capacity: np.ndarray
    delivered_specification: np.ndarray
    residual_specification: np.ndarray
    implied_field_difference: np.ndarray | None
    gross_material_volume: float
    net_material_volume: float
    singular_values: np.ndarray
    numerical_rank: int
    retained_condition: float
    bounded_solution: AbeBoundedCurrentPotentialSolution


def solve_abe_element_fill_plan(
        specification_response, requested_specification_difference, *,
        material_active, element_volumes, element_ids=None,
        current_fill=None, field_response=None,
        relative_singular_threshold=1.0e-12, **solve_options
        ) -> AbeElementFillPlan:
    """Solve a specification directly in physically bounded cell fills.

    ``specification_response`` has one column per design cell.  A column is
    normally ``J_spec_field @ element_field_response``.  Solving the compact
    specification rows directly lets the material inverse choose an
    inexpensive member of the field manifold; ``field_response`` is optional
    and is used only to report the implied field difference afterward.
    """
    response = _finite_array(
        specification_response, name="specification_response", ndim=2)
    row_count, cell_count = response.shape
    if min(response.shape) <= 0:
        raise ValueError("specification_response must be non-empty")
    requested = _finite_array(
        requested_specification_difference,
        name="requested_specification_difference").reshape(-1)
    if requested.shape != (row_count,):
        raise ValueError("requested specification has the wrong length")
    active = np.asarray(material_active, dtype=bool).reshape(-1)
    volumes = _finite_array(
        element_volumes, name="element_volumes").reshape(-1)
    if (active.shape != (cell_count,) or volumes.shape != (cell_count,)
            or np.any(volumes <= 0.0)):
        raise ValueError(
            "material_active and positive element_volumes must match columns")
    fill = (np.zeros(cell_count, dtype=float) if current_fill is None else
            _finite_array(current_fill, name="current_fill").reshape(-1))
    if fill.shape != (cell_count,):
        raise ValueError("current_fill must match the design cells")
    lower = np.where(active, -1.0, 0.0)
    upper = np.where(active, 0.0, 1.0)
    tolerance = 64.0 * np.finfo(float).eps
    if np.any(fill < lower - tolerance) or np.any(fill > upper + tolerance):
        raise ValueError("current_fill violates signed material capacity")
    remaining_lower = lower - fill
    remaining_upper = upper - fill
    if "lower_potential" in solve_options or "upper_potential" in solve_options:
        raise ValueError("element-fill bounds are owned by this function")
    if "relative_singular_threshold" in solve_options:
        raise ValueError("relative_singular_threshold was supplied twice")
    bounded = solve_abe_bounded_current_potential(
        response, requested, lower_potential=remaining_lower,
        upper_potential=remaining_upper,
        relative_singular_threshold=float(relative_singular_threshold),
        **solve_options)
    step = np.asarray(bounded.solution.potential, dtype=float).reshape(-1)
    absolute = np.clip(fill + step, lower, upper)
    step = absolute - fill
    delivered = response @ step
    residual = requested - delivered
    implied = None
    if field_response is not None:
        field = _finite_array(
            field_response, name="field_response", ndim=2)
        if field.shape[1] != cell_count:
            raise ValueError("field_response must use the same design columns")
        implied = np.ascontiguousarray(field @ step)
    singular = np.asarray(bounded.solution.factor.S, dtype=float).reshape(-1)
    if singular.size and singular[0] > 0.0:
        rank = int(np.count_nonzero(
            singular > float(relative_singular_threshold) * singular[0]))
    else:
        rank = 0
    condition = (float(singular[0] / singular[rank - 1])
                 if rank else np.inf)
    ids = (np.arange(cell_count, dtype=np.int64) if element_ids is None else
           np.asarray(element_ids, dtype=np.int64).reshape(-1))
    if ids.shape != (cell_count,):
        raise ValueError("element_ids must match the design cells")
    return AbeElementFillPlan(
        element_ids=ids.copy(), fill_step=np.ascontiguousarray(step),
        absolute_fill=np.ascontiguousarray(absolute),
        lower_capacity=np.ascontiguousarray(lower),
        upper_capacity=np.ascontiguousarray(upper),
        delivered_specification=np.ascontiguousarray(delivered),
        residual_specification=np.ascontiguousarray(residual),
        implied_field_difference=implied,
        gross_material_volume=float(np.sum(np.abs(step) * volumes)),
        net_material_volume=float(np.sum(step * volumes)),
        singular_values=np.ascontiguousarray(singular), numerical_rank=rank,
        retained_condition=condition, bounded_solution=bounded)


@dataclass(frozen=True)
class BinnedInterfaceHeight:
    """Conservative signed material-volume to interface-height conversion."""

    first_edges: np.ndarray
    second_edges: np.ndarray
    bin_areas: np.ndarray
    signed_volume: np.ndarray
    height_change: np.ndarray
    element_count: np.ndarray

    @property
    def first_centres(self):
        return 0.5 * (self.first_edges[:-1] + self.first_edges[1:])

    @property
    def second_centres(self):
        return 0.5 * (self.second_edges[:-1] + self.second_edges[1:])

    def sample(self, first_coordinate, second_coordinate):
        """Bilinearly sample the height field; return zero outside its bins."""
        from scipy.interpolate import RegularGridInterpolator

        first = np.asarray(first_coordinate, dtype=float)
        second = np.asarray(second_coordinate, dtype=float)
        first, second = np.broadcast_arrays(first, second)
        points = np.column_stack((first.reshape(-1), second.reshape(-1)))
        interpolator = RegularGridInterpolator(
            (self.first_centres, self.second_centres), self.height_change,
            bounds_error=False, fill_value=0.0)
        return np.asarray(interpolator(points), dtype=float).reshape(first.shape)


def bin_element_fill_to_interface_height(
        fill_fraction, element_volumes, first_coordinate, second_coordinate,
        first_edges, second_edges, bin_areas) -> BinnedInterfaceHeight:
    """Conservatively convert signed fill volume into a height field.

    The convention matches a pole face above an aperture: positive fill adds
    iron and therefore lowers the interface, so
    ``height_change = -signed_volume / area``.
    """
    fill = _finite_array(fill_fraction, name="fill_fraction").reshape(-1)
    volume = _finite_array(
        element_volumes, name="element_volumes").reshape(-1)
    first = _finite_array(
        first_coordinate, name="first_coordinate").reshape(-1)
    second = _finite_array(
        second_coordinate, name="second_coordinate").reshape(-1)
    if not (fill.shape == volume.shape == first.shape == second.shape):
        raise ValueError("fill, volume, and element coordinates must match")
    if np.any(volume <= 0.0):
        raise ValueError("element_volumes must be positive")
    edges1 = _finite_array(first_edges, name="first_edges").reshape(-1)
    edges2 = _finite_array(second_edges, name="second_edges").reshape(-1)
    if (edges1.size < 2 or edges2.size < 2
            or np.any(np.diff(edges1) <= 0.0)
            or np.any(np.diff(edges2) <= 0.0)):
        raise ValueError("interface bin edges must be strictly increasing")
    areas = _finite_array(bin_areas, name="bin_areas", ndim=2)
    shape = (edges1.size - 1, edges2.size - 1)
    if areas.shape != shape or np.any(areas <= 0.0):
        raise ValueError("bin_areas must be positive with one value per bin")
    i = np.searchsorted(edges1, first, side="right") - 1
    j = np.searchsorted(edges2, second, side="right") - 1
    i[first == edges1[-1]] = shape[0] - 1
    j[second == edges2[-1]] = shape[1] - 1
    valid = ((i >= 0) & (i < shape[0]) & (j >= 0) & (j < shape[1]))
    if not np.all(valid):
        raise ValueError("every design element must lie inside the interface bins")
    signed = np.zeros(shape, dtype=float)
    counts = np.zeros(shape, dtype=np.int64)
    np.add.at(signed, (i, j), fill * volume)
    np.add.at(counts, (i, j), 1)
    height = -signed / areas
    return BinnedInterfaceHeight(
        first_edges=np.ascontiguousarray(edges1),
        second_edges=np.ascontiguousarray(edges2),
        bin_areas=np.ascontiguousarray(areas),
        signed_volume=np.ascontiguousarray(signed),
        height_change=np.ascontiguousarray(height),
        element_count=np.ascontiguousarray(counts))


def blended_interface_displacement(
        normal_coordinate, interface_coordinate, lower_fixed_coordinate,
        upper_fixed_coordinate, interface_height_change):
    """Blend a face-height update to zero at fixed aperture/root surfaces."""
    normal, interface, lower, upper, change = np.broadcast_arrays(
        _finite_array(normal_coordinate, name="normal_coordinate"),
        _finite_array(interface_coordinate, name="interface_coordinate"),
        _finite_array(lower_fixed_coordinate, name="lower_fixed_coordinate"),
        _finite_array(upper_fixed_coordinate, name="upper_fixed_coordinate"),
        _finite_array(interface_height_change,
                      name="interface_height_change"))
    if np.any(interface <= lower) or np.any(upper <= interface):
        raise ValueError(
            "the reference interface must lie strictly between fixed surfaces")
    weight = np.zeros(normal.shape, dtype=float)
    below = (normal > lower) & (normal < interface)
    above = (normal >= interface) & (normal < upper)
    weight[below] = ((normal[below] - lower[below])
                     / (interface[below] - lower[below]))
    weight[above] = ((upper[above] - normal[above])
                     / (upper[above] - interface[above]))
    return np.ascontiguousarray(change * weight)


@dataclass(frozen=True)
class ExactSectionEvaluation:
    """Exact complete-solve result consumed by the outer material loop."""

    specification: np.ndarray
    payload: object = None


@dataclass(frozen=True)
class AbeContourIteration:
    iteration: int
    requested_difference: np.ndarray
    fill_step: np.ndarray
    accumulated_fill: np.ndarray
    backtracking_scale: float
    predicted_max_band_ratio: float
    exact_max_band_ratio: float
    gross_material_volume: float
    net_material_volume: float
    inner_converged: bool
    inner_stop_reason: str
    exact_evaluation: ExactSectionEvaluation
    realization: object


@dataclass(frozen=True)
class AbeContourOptimizationResult:
    target_specification: np.ndarray
    response_band: np.ndarray
    initial_evaluation: ExactSectionEvaluation
    final_evaluation: ExactSectionEvaluation
    accumulated_fill: np.ndarray
    history: tuple[AbeContourIteration, ...]
    converged: bool
    stop_reason: str

    @property
    def initial_max_band_ratio(self):
        return float(np.max(np.abs(
            (self.initial_evaluation.specification - self.target_specification)
            / self.response_band)))

    @property
    def final_max_band_ratio(self):
        return float(np.max(np.abs(
            (self.final_evaluation.specification - self.target_specification)
            / self.response_band)))


def optimize_abe_section_contour(
        specification_response, *, target_specification,
        response_band, initial_evaluation, material_active, element_volumes,
        realize_fill: Callable, evaluate_exact: Callable,
        element_ids=None, field_response=None, maximum_iterations=4,
        inner_residual_fraction=0.02,
        relative_singular_threshold=1.0e-12,
        backtracking_scales=(1.0, 0.5, 0.25, 0.125),
        minimum_exact_improvement=1.0e-12, exact_guard=None,
        relinearize=None, solve_options=None) -> AbeContourOptimizationResult:
    """Close the fill->smooth-contour->complete-solve loop.

    ``realize_fill(accumulated_fill)`` constructs geometry from the reference
    shape.  ``evaluate_exact(realization)`` must perform the complete field
    solve and native map evaluation and return :class:`ExactSectionEvaluation`.
    Optional ``relinearize(evaluation, realization)`` returns a refreshed
    specification-response matrix after an accepted step.  No design finite
    difference is performed here.
    """
    response = _finite_array(
        specification_response, name="specification_response", ndim=2)
    target = _finite_array(
        target_specification, name="target_specification").reshape(-1)
    band = _finite_array(response_band, name="response_band").reshape(-1)
    if (target.shape != (response.shape[0],) or band.shape != target.shape
            or np.any(band <= 0.0)):
        raise ValueError("target and positive response bands must match rows")
    if not isinstance(initial_evaluation, ExactSectionEvaluation):
        raise TypeError("initial_evaluation must be ExactSectionEvaluation")
    initial_spec = _finite_array(
        initial_evaluation.specification,
        name="initial exact specification").reshape(-1)
    if initial_spec.shape != target.shape:
        raise ValueError("initial exact specification has the wrong length")
    cell_count = response.shape[1]
    active = np.asarray(material_active, dtype=bool).reshape(-1)
    volumes = _finite_array(
        element_volumes, name="element_volumes").reshape(-1)
    if active.shape != (cell_count,) or volumes.shape != (cell_count,):
        raise ValueError("material metadata must match response columns")
    iterations = int(maximum_iterations)
    fraction = float(inner_residual_fraction)
    improvement = float(minimum_exact_improvement)
    scales = tuple(float(value) for value in backtracking_scales)
    if (iterations < 1 or not np.isfinite(fraction) or not 0.0 < fraction < 1.0
            or not np.isfinite(improvement) or improvement < 0.0
            or not scales or any(not np.isfinite(value) or not 0.0 < value <= 1.0
                              for value in scales)):
        raise ValueError("outer-loop controls are invalid")
    if exact_guard is not None and not bool(exact_guard(initial_evaluation)):
        raise ValueError("initial exact evaluation violates exact_guard")
    options = {} if solve_options is None else dict(solve_options)
    forbidden = {
        "residual_rms", "lower_potential", "upper_potential",
        "relative_singular_threshold", "initial_potential"}
    overlap = forbidden.intersection(options)
    if overlap:
        raise ValueError("solve_options contain outer-loop-owned keys: "
                         + ", ".join(sorted(overlap)))

    def ratio(specification):
        return float(np.max(np.abs(
            (np.asarray(specification, dtype=float) - target) / band)))

    current = ExactSectionEvaluation(initial_spec.copy(), initial_evaluation.payload)
    current_ratio = ratio(current.specification)
    accumulated = np.zeros(cell_count, dtype=float)
    history = []
    stop_reason = "maximum_iterations"
    for outer in range(iterations):
        if current_ratio <= 1.0:
            stop_reason = "target_band_reached"
            break
        requested = target - current.specification
        requested_rms = float(np.sqrt(np.mean(requested * requested)))
        plan = solve_abe_element_fill_plan(
            response, requested, material_active=active,
            element_volumes=volumes, element_ids=element_ids,
            current_fill=accumulated, field_response=field_response,
            residual_rms=max(
                fraction * requested_rms, np.finfo(float).tiny),
            relative_singular_threshold=relative_singular_threshold,
            **options)
        if np.max(np.abs(plan.fill_step), initial=0.0) <= 64.0 * np.finfo(float).eps:
            stop_reason = "material_inverse_proposed_no_fill_change"
            break
        accepted = None
        for scale in scales:
            trial_fill = accumulated + scale * plan.fill_step
            trial_fill = np.minimum(
                np.maximum(trial_fill, plan.lower_capacity),
                plan.upper_capacity)
            realization = realize_fill(np.ascontiguousarray(trial_fill.copy()))
            exact = evaluate_exact(realization)
            if not isinstance(exact, ExactSectionEvaluation):
                raise TypeError(
                    "evaluate_exact must return ExactSectionEvaluation")
            exact_spec = _finite_array(
                exact.specification,
                name="trial exact specification").reshape(-1)
            if exact_spec.shape != target.shape:
                raise ValueError("trial exact specification has the wrong length")
            exact = ExactSectionEvaluation(exact_spec.copy(), exact.payload)
            if exact_guard is not None and not bool(exact_guard(exact)):
                continue
            exact_ratio = ratio(exact_spec)
            if exact_ratio >= current_ratio - improvement:
                continue
            scaled_step = trial_fill - accumulated
            predicted = current.specification + response @ scaled_step
            accepted = (scale, trial_fill, scaled_step, realization, exact,
                        ratio(predicted), exact_ratio)
            break
        if accepted is None:
            stop_reason = "exact_backtracking_exhausted"
            break
        scale, accumulated, step, realization, current, predicted_ratio, current_ratio = accepted
        history.append(AbeContourIteration(
            iteration=outer, requested_difference=requested.copy(),
            fill_step=np.ascontiguousarray(step),
            accumulated_fill=np.ascontiguousarray(accumulated.copy()),
            backtracking_scale=float(scale),
            predicted_max_band_ratio=float(predicted_ratio),
            exact_max_band_ratio=float(current_ratio),
            gross_material_volume=float(np.sum(np.abs(step) * volumes)),
            net_material_volume=float(np.sum(step * volumes)),
            inner_converged=bool(plan.bounded_solution.converged),
            inner_stop_reason=str(plan.bounded_solution.stop_reason),
            exact_evaluation=current, realization=realization))
        if relinearize is not None:
            response = _finite_array(
                relinearize(current, realization),
                name="relinearized specification response", ndim=2)
            if response.shape != (target.size, cell_count):
                raise ValueError(
                    "relinearized response has an incompatible shape")
    converged = bool(current_ratio <= 1.0)
    if converged:
        stop_reason = "target_band_reached"
    return AbeContourOptimizationResult(
        target_specification=target.copy(), response_band=band.copy(),
        initial_evaluation=ExactSectionEvaluation(
            initial_spec.copy(), initial_evaluation.payload),
        final_evaluation=current,
        accumulated_fill=np.ascontiguousarray(accumulated),
        history=tuple(history), converged=converged,
        stop_reason=stop_reason)


__all__ = [
    "AbeContourIteration",
    "AbeContourOptimizationResult",
    "AbeElementFillPlan",
    "BinnedInterfaceHeight",
    "ExactSectionEvaluation",
    "bin_element_fill_to_interface_height",
    "blended_interface_displacement",
    "compose_specification_fill_response",
    "contract_hdiv_element_fill_response",
    "measured_element_fill_patterns",
    "optimize_abe_section_contour",
    "solve_abe_element_fill_plan",
]
