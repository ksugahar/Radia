"""Multi-momentum EarlyTimes Taylor-map fusion for FFAG HDiv-MMM.

This module is the fixed-design-orbit bridge between the single-orbit
second-order Taylor-map kernel and the binary HDiv-MMM material optimizer.
The raw electromagnetic coordinates are normal dipole, normal/skew
quadrupole, and normal/skew sextupole coefficients on every orbit segment.
Their forward-AD ``R/T`` Jacobian is contracted with all material candidates
before ACA--thin-QR--TSVD and connected graph-front selection.

Every accepted topology is a fully solved binary active set and is scored by
the native C++ variational-map value.  Finite differences, density variables,
and an air volume mesh are not used.  The caller owns ``ngsolve.TaskManager``.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .accelerator_magnet_topopt import (
    CoilBuilderHDivSource,
    PlanarDesignOrbit,
)
from .accelerator_taylor_topopt import (
    PlanarSecondOrderTaylorMapObjective,
    build_planar_orbit_multipole_response_matrix,
    planar_orbit_multipole_observations,
    second_order_taylor_map_from_multipoles,
)
from .topology_optimization import (
    GrowthTopologyReport,
    HDivMMMGenerationResult,
    grow_hdiv_mmm_by_superposition,
    ngsolve_growth_topology,
)


def _finite_array(value, *, name):
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class MultiMomentumSecondOrderTaylorMapObjective:
    """Block objective for several fixed design orbits/rigidities."""

    objectives: tuple[PlanarSecondOrderTaylorMapObjective, ...]

    def __post_init__(self):
        objectives = tuple(self.objectives)
        if (not objectives or any(
                not isinstance(item, PlanarSecondOrderTaylorMapObjective)
                for item in objectives)):
            raise TypeError(
                "objectives must contain PlanarSecondOrderTaylorMapObjective")
        object.__setattr__(self, "objectives", objectives)

    @property
    def orbits(self) -> tuple[PlanarDesignOrbit, ...]:
        return tuple(item.orbit for item in self.objectives)

    @property
    def derivative_backend(self) -> str:
        return "block-forward-mode-rk4-taylor-ad"

    @property
    def raw_field_response_size(self) -> int:
        return sum(item.raw_field_response_size for item in self.objectives)

    @property
    def raw_offsets(self) -> np.ndarray:
        sizes = [item.raw_field_response_size for item in self.objectives]
        return np.r_[0, np.cumsum(sizes)].astype(np.int64)

    @property
    def response_target(self) -> np.ndarray:
        return np.concatenate([item.response_target
                               for item in self.objectives])

    @property
    def response_band(self) -> np.ndarray:
        return np.concatenate([item.response_band
                               for item in self.objectives])

    def response_group_indices(self, groups) -> np.ndarray:
        """Return global design-response rows for selected metric groups."""
        requested = tuple(groups)
        allowed = ("normal_dipole", "R", "T")
        if (not requested or len(set(requested)) != len(requested)
                or any(name not in allowed for name in requested)):
            raise ValueError(
                "response groups must be a non-empty unique subset of "
                "('normal_dipole', 'R', 'T')")
        wanted = set(requested)
        rows = []
        offset = 0
        for item in self.objectives:
            for name, response_slice in item.response_slices:
                if name in wanted:
                    rows.extend(offset + np.arange(
                        response_slice.start, response_slice.stop,
                        dtype=np.int64))
            offset += item.response_target.size
        return np.asarray(rows, dtype=np.int64)

    def group_max_band_ratio(self, design_response, group) -> float:
        """Evaluate one response group's global multi-momentum minimax ratio."""
        values = _finite_array(
            design_response, name="multi-momentum design response").reshape(-1)
        if values.shape != self.response_target.shape:
            raise ValueError("design response does not match the objective")
        rows = self.response_group_indices((group,))
        return float(np.max(np.abs(
            (values[rows] - self.response_target[rows])
            / self.response_band[rows])))

    @property
    def source_calibration_rows(self) -> np.ndarray:
        rows = []
        for offset, item in zip(self.raw_offsets[:-1], self.objectives):
            count = len(item.orbit.segment_lengths)
            rows.extend(int(offset) + np.arange(count, dtype=np.int64))
        return np.asarray(rows, dtype=np.int64)

    @property
    def source_calibration_target(self) -> np.ndarray:
        return np.concatenate([
            item.required_normal_dipole for item in self.objectives])

    @property
    def source_calibration_band(self) -> np.ndarray:
        return np.concatenate([
            np.broadcast_to(
                item.normal_dipole_band,
                item.required_normal_dipole.shape).astype(float)
            for item in self.objectives])

    def split_raw_response(self, field_response) -> tuple[np.ndarray, ...]:
        values = _finite_array(
            field_response, name="multi-momentum multipole response"
        ).reshape(-1)
        if values.shape != (self.raw_field_response_size,):
            raise ValueError(
                "multipole response does not match the multi-momentum "
                "raw-row contract")
        offsets = self.raw_offsets
        return tuple(values[left:right] for left, right in zip(
            offsets[:-1], offsets[1:]))

    def transform(self, field_response) -> np.ndarray:
        return np.concatenate([
            item.transform(values)
            for item, values in zip(
                self.objectives, self.split_raw_response(field_response))])

    def transform_jacobian(self, field_response) -> np.ndarray:
        values = self.split_raw_response(field_response)
        blocks = [item.transform_jacobian(value)
                  for item, value in zip(self.objectives, values)]
        row_count = sum(block.shape[0] for block in blocks)
        result = np.zeros((row_count, self.raw_field_response_size))
        response_offset = 0
        for block, left, right in zip(
                blocks, self.raw_offsets[:-1], self.raw_offsets[1:]):
            next_offset = response_offset + block.shape[0]
            result[response_offset:next_offset, left:right] = block
            response_offset = next_offset
        return result


def build_multi_orbit_multipole_response_matrix(
        charge_gram, objective: MultiMomentumSecondOrderTaylorMapObjective,
        *, sample_radius, field_scale=None) -> np.ndarray:
    """Build HDiv rows for every momentum-indexed multipole observation."""
    if not isinstance(objective, MultiMomentumSecondOrderTaylorMapObjective):
        raise TypeError(
            "objective must be MultiMomentumSecondOrderTaylorMapObjective")
    radii = ((float(sample_radius),) * len(objective.objectives)
             if np.isscalar(sample_radius) else tuple(sample_radius))
    if len(radii) != len(objective.objectives):
        raise ValueError(
            "sample_radius must be scalar or have one value per momentum")
    options = {} if field_scale is None else {"field_scale": field_scale}
    rows = [build_planar_orbit_multipole_response_matrix(
        charge_gram, item.orbit, sample_radius=radius, **options)
        for item, radius in zip(objective.objectives, radii)]
    result = np.ascontiguousarray(np.vstack(rows))
    if result.shape[0] != objective.raw_field_response_size:
        raise RuntimeError("multi-orbit multipole row count mismatch")
    return result


def incident_multi_orbit_multipole_response(
        source: CoilBuilderHDivSource,
        objective: MultiMomentumSecondOrderTaylorMapObjective, *,
        sample_radius) -> np.ndarray:
    """Project the same CoilBuilder source onto all Taylor multipole rows."""
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    if not isinstance(objective, MultiMomentumSecondOrderTaylorMapObjective):
        raise TypeError(
            "objective must be MultiMomentumSecondOrderTaylorMapObjective")
    radii = ((float(sample_radius),) * len(objective.objectives)
             if np.isscalar(sample_radius) else tuple(sample_radius))
    if len(radii) != len(objective.objectives):
        raise ValueError(
            "sample_radius must be scalar or have one value per momentum")
    responses = []
    for item, radius in zip(objective.objectives, radii):
        points, weights = planar_orbit_multipole_observations(
            item.orbit, sample_radius=radius)
        responses.append(np.einsum(
            "rpc,pc->r", weights, source.b_field(points)))
    result = np.ascontiguousarray(np.concatenate(responses))
    if result.shape != (objective.raw_field_response_size,):
        raise RuntimeError("CoilBuilder multipole response size mismatch")
    return result


def build_ffag_second_order_taylor_objective(
        target_family, *, T_band=1.0e6, T_entries=None,
        maximum_step_m=1.0e-3
        ) -> MultiMomentumSecondOrderTaylorMapObjective:
    """Lift a fixed-orbit FFAG ``R`` target into an EarlyTimes ``R/T`` target.

    The requested first-order matrices remain exactly those in
    ``target_family``.  The otherwise unspecified second-order target is the
    Taylor map of the same ideal dipole/normal-gradient profile, while its
    band is intentionally caller controlled.  A loose band therefore gives a
    faithful EarlyTimes ``R`` comparison without silently demanding zero
    aberration; it can be tightened later as a continuation stage.
    """
    from .ffag_topopt import FFAGCellTargetFamily

    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    objectives = []
    for index, reference in enumerate(target_family.references):
        count = len(reference.orbit.segment_lengths)
        ideal = np.zeros((5, count), dtype=float)
        ideal[0] = reference.field_response[:count]
        ideal[1] = reference.field_response[count:]
        ideal_map = second_order_taylor_map_from_multipoles(
            ideal.reshape(-1), reference.orbit.segment_lengths,
            reference.orbit.magnetic_rigidity,
            curvature_sign=target_family.objective.curvature_sign,
            gradient_sign=target_family.objective.gradient_sign,
            maximum_step_m=maximum_step_m)
        options = {} if T_entries is None else {"T_entries": tuple(T_entries)}
        objectives.append(PlanarSecondOrderTaylorMapObjective(
            orbit=reference.orbit,
            target_R=target_family.objective.target_matrices[index],
            target_T=ideal_map.T,
            R_band=target_family.objective.transfer_matrix_band[index],
            T_band=T_band,
            normal_dipole_band=(
                target_family.objective.bend_field_band[index]),
            R_entries=target_family.objective.response_entries,
            curvature_sign=target_family.objective.curvature_sign,
            gradient_sign=target_family.objective.gradient_sign,
            maximum_step_m=maximum_step_m,
            **options))
    return MultiMomentumSecondOrderTaylorMapObjective(tuple(objectives))


@dataclass(frozen=True)
class FFAGSecondOrderTaylorTopologyResult:
    """Binary multi-momentum FFAG result scored by native ``R/T`` maps."""

    objective: MultiMomentumSecondOrderTaylorMapObjective
    generation: HDivMMMGenerationResult
    source_scale: float
    realized_multipole_responses: tuple[np.ndarray, ...]
    realized_R: np.ndarray
    realized_T: np.ndarray
    normal_dipole_max_band_ratios: np.ndarray
    R_max_band_ratios: np.ndarray
    T_max_band_ratios: np.ndarray
    topology: GrowthTopologyReport
    primary_response_groups: tuple[str, ...] = (
        "normal_dipole", "R", "T")
    primary_max_band_ratio: float = np.nan
    maximum_group_band_ratios: tuple[tuple[str, float], ...] = ()

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def max_band_ratio(self) -> float:
        return float(max(
            np.max(self.normal_dipole_max_band_ratios),
            np.max(self.R_max_band_ratios),
            np.max(self.T_max_band_ratios)))

    @property
    def converged(self) -> bool:
        return self.max_band_ratio <= 1.0


def optimize_ffag_hdiv_mmm_from_second_order_taylor_maps(
        objective: MultiMomentumSecondOrderTaylorMapObjective, *, source,
        charge_gram, fes, inv_chi, active_elements, element_volumes,
        volume_max, sample_radius, source_scale=1.0,
        optimize_source_scale=True, multipole_response_matrix=None,
        incident_multipole_response=None,
        primary_response_groups=None,
        maximum_group_band_ratios=None,
        **generation_options) -> FFAGSecondOrderTaylorTopologyResult:
    """Optimize one binary magnet directly against multi-momentum ``R/T``.

    Source amplitude is eliminated on every exact active solve when
    ``optimize_source_scale`` is true.  The derivative of that elimination is
    composed with the Taylor-map Jacobian inside the material contraction.
    ``primary_response_groups`` may select the rows minimized by ACA/QR/TSVD,
    LP, and exact acceptance.  ``maximum_group_band_ratios`` supplies absolute
    group-wise physics guards evaluated only after the calibrated complete
    active-set solve.  For example, primary ``('R',)`` with a
    ``normal_dipole`` limit minimizes the transfer matrix without allowing the
    bend-field profile to deteriorate.  Defaults preserve the original joint
    minimax objective.
    """
    if not isinstance(objective, MultiMomentumSecondOrderTaylorMapObjective):
        raise TypeError(
            "objective must be MultiMomentumSecondOrderTaylorMapObjective")
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    scale = float(source_scale)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("source_scale must be positive and finite")
    response_matrix = (
        build_multi_orbit_multipole_response_matrix(
            charge_gram, objective, sample_radius=sample_radius)
        if multipole_response_matrix is None else
        _finite_array(
            multipole_response_matrix,
            name="multipole_response_matrix"))
    expected = (objective.raw_field_response_size, int(fes.ndof))
    if response_matrix.shape != expected:
        raise ValueError(
            "multipole_response_matrix must have shape "
            "(objective.raw_field_response_size,fes.ndof)")
    incident = (
        incident_multi_orbit_multipole_response(
            source, objective, sample_radius=sample_radius)
        if incident_multipole_response is None else
        _finite_array(
            incident_multipole_response,
            name="incident_multipole_response").reshape(-1))
    if incident.shape != (objective.raw_field_response_size,):
        raise ValueError(
            "incident_multipole_response must match the objective")
    primary_groups = (("normal_dipole", "R", "T")
                      if primary_response_groups is None else
                      tuple(primary_response_groups))
    primary_rows = objective.response_group_indices(primary_groups)
    if (maximum_group_band_ratios is not None
            and not isinstance(maximum_group_band_ratios, Mapping)):
        raise TypeError("maximum_group_band_ratios must be a mapping")
    guards = ({} if maximum_group_band_ratios is None else
              dict(maximum_group_band_ratios))
    for name, value in guards.items():
        objective.response_group_indices((name,))
        limit = float(value)
        if not np.isfinite(limit) or limit < 0.0:
            raise ValueError("group band-ratio limits must be finite and nonnegative")
        guards[name] = limit
    full_target = objective.response_target
    full_band = objective.response_band

    def primary_transform(raw_response):
        return objective.transform(raw_response)[primary_rows]

    def primary_transform_jacobian(raw_response):
        return objective.transform_jacobian(raw_response)[primary_rows, :]

    def exact_group_guard(raw_response, primary_response):
        del primary_response
        design_response = objective.transform(raw_response)
        for name, limit in guards.items():
            rows = objective.response_group_indices((name,))
            ratio = float(np.max(np.abs(
                (design_response[rows] - full_target[rows])
                / full_band[rows])))
            tolerance = 1.0e-10 * max(1.0, limit)
            if ratio > limit + tolerance:
                return False
        return True
    source_rhs = source.assemble_hdiv_rhs(fes)
    rhs = scale * source_rhs
    incident = scale * incident
    reserved = {
        "response_matrix", "response_target", "response_band",
        "response_transform", "response_transform_jacobian",
        "incident_response", "source_calibration_rows",
        "source_calibration_target", "source_calibration_band",
        "source_calibration_norm",
        "exact_response_validator",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the FFAG Taylor contract: "
            + ", ".join(sorted(overlap)))
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=np.ascontiguousarray(response_matrix),
        active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=full_target[primary_rows],
        response_band=full_band[primary_rows],
        volume_max=volume_max,
        incident_response=np.ascontiguousarray(incident),
        response_transform=(objective.transform if
            primary_rows.size == full_target.size else primary_transform),
        response_transform_jacobian=(objective.transform_jacobian if
            primary_rows.size == full_target.size else
            primary_transform_jacobian),
        exact_response_validator=(exact_group_guard if guards else None),
        source_calibration_rows=(
            objective.source_calibration_rows
            if optimize_source_scale else None),
        source_calibration_target=(
            objective.source_calibration_target
            if optimize_source_scale else None),
        source_calibration_band=(
            objective.source_calibration_band
            if optimize_source_scale else None),
        source_calibration_norm=(
            "linf" if optimize_source_scale else "mean"),
        **generation_options)
    raw_by_orbit = objective.split_raw_response(generation.response)
    maps = [item.evaluate_taylor_map(raw)
            for item, raw in zip(objective.objectives, raw_by_orbit)]
    dipole_ratios = []
    R_ratios = []
    T_ratios = []
    for item, raw, transfer in zip(
            objective.objectives, raw_by_orbit, maps):
        count = len(item.orbit.segment_lengths)
        dipole_ratios.append(float(np.max(np.abs(
            (raw[:count] - item.required_normal_dipole)
            / item.normal_dipole_band))))
        R_ratios.append(float(max(abs(
            (transfer.R[index] - item.target_R[index])
            / item.R_band[index]) for index in item.R_entries)))
        T_ratios.append(float(max(abs(
            (transfer.T[index] - item.target_T[index])
            / item.T_band[index]) for index in item.T_entries)))
    topology = ngsolve_growth_topology(
        fes.mesh, generation.active_elements)
    return FFAGSecondOrderTaylorTopologyResult(
        objective=objective,
        generation=generation,
        source_scale=scale * generation.source_scale,
        realized_multipole_responses=tuple(
            np.asarray(raw, dtype=float).copy() for raw in raw_by_orbit),
        realized_R=np.asarray([item.R for item in maps], dtype=float),
        realized_T=np.asarray([item.T for item in maps], dtype=float),
        normal_dipole_max_band_ratios=np.asarray(dipole_ratios),
        R_max_band_ratios=np.asarray(R_ratios),
        T_max_band_ratios=np.asarray(T_ratios),
        topology=topology,
        primary_response_groups=tuple(primary_groups),
        primary_max_band_ratio=float(np.max(np.abs(
            (generation.objective_response - full_target[primary_rows])
            / full_band[primary_rows]))),
        maximum_group_band_ratios=tuple(
            (str(name), float(value)) for name, value in guards.items()))


__all__ = [
    "FFAGSecondOrderTaylorTopologyResult",
    "MultiMomentumSecondOrderTaylorMapObjective",
    "build_ffag_second_order_taylor_objective",
    "build_multi_orbit_multipole_response_matrix",
    "incident_multi_orbit_multipole_response",
    "optimize_ffag_hdiv_mmm_from_second_order_taylor_maps",
]
