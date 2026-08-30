"""Air-volume-mesh-free HDiv-MMM proof of concept for FFAG cells.

The module turns a momentum-indexed family of periodic design orbits and
first-order cell maps into the multi-orbit contract in
``accelerator_magnet_topopt``.  The Bell--Abell non-scaling FFAG parameters
provide a reproducible soft-edge fixture, not a claim to reproduce their PTC
ring: their complete placement and closed-orbit files are not published in the
paper.

The optimization path uses exact combined-function matrix-exponential Frechet
derivatives.  Enge ``I1``/``I2`` values are diagnostics for a reduced edge map;
they are never applied on top of a field response which already samples the
soft fringe.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .accelerator_magnet_topopt import (
    CoilBuilderHDivSource,
    CoilHDivTotalField,
    MultiMomentumAcceleratorMagnetTopologyResult,
    MultiMomentumTransferMatrixObjective,
    PlanarDesignOrbit,
    build_multi_orbit_field_response_matrix,
    optimize_hdiv_mmm_magnet_from_transfer_matrices,
    planar_orbit_field_observations,
    solve_transfer_matrix_field_correction,
    static_magnet_symplectic_residual,
    static_magnet_transfer_component_entries,
)
from .isochronous_topopt import (
    CombinedFunctionTransferMap,
    combined_function_transfer_map_from_field_response,
)


PROTON_REST_ENERGY_MEV = 938.27208816
GEV_C_PER_TESLA_METRE = 0.299792458


@dataclass(frozen=True)
class FFAGCyclicSectorContract:
    """Validated interpretation of one rotational FFAG sector.

    ``image_cyclic`` is sufficient when complete iron bodies fit inside one
    sector and their rotated copies are disjoint.  A continuous return yoke
    cut by the azimuthal sector planes needs a local closure in addition to
    the rotated nonlocal interaction.  Conforming FEM identifies the two HDiv
    normal traces.  Broken-HDiv VIM keeps the element unknowns independent and
    instead pairs the two surface-charge rows into one periodic jump.
    """

    fold: int
    field_antiperiodic: bool
    body_crosses_periodic_planes: bool
    formulation: str
    periodic_trace_identified: bool
    periodic_charge_paired: bool
    reduction_mode: str


@dataclass(frozen=True)
class FFAGCyclicDensityMap:
    """Independent material variables for a rotational FFAG quotient mesh.

    HDiv trace identification and broken-VIM charge pairing close the field
    space.  They do not by themselves make an element-wise topology density
    periodic.  This map identifies the volume elements adjacent to paired
    azimuthal faces, expands one reduced design vector to all elements, and
    contracts element gradients with the exact transpose map.
    """

    element_to_variable: np.ndarray
    groups: tuple[tuple[int, ...], ...]
    boundary_pair_count: int
    periodic_boundaries: tuple[str, str]

    def __post_init__(self):
        mapping = np.asarray(self.element_to_variable, dtype=np.int64).reshape(-1)
        groups = tuple(tuple(int(value) for value in group)
                       for group in self.groups)
        if (mapping.size == 0 or not groups
                or np.any(mapping < 0)
                or set(mapping.tolist()) != set(range(len(groups)))):
            raise ValueError("invalid FFAG cyclic density map")
        expected = tuple(tuple(np.flatnonzero(mapping == index).tolist())
                         for index in range(len(groups)))
        if groups != expected:
            raise ValueError(
                "FFAG cyclic density groups do not match element mapping")
        boundaries = tuple(str(name) for name in self.periodic_boundaries)
        if len(boundaries) != 2 or boundaries[0] == boundaries[1]:
            raise ValueError("FFAG cyclic density map needs two boundaries")
        mapping.setflags(write=False)
        object.__setattr__(self, "element_to_variable", mapping)
        object.__setattr__(self, "groups", groups)
        object.__setattr__(self, "periodic_boundaries", boundaries)
        object.__setattr__(self, "boundary_pair_count",
                           int(self.boundary_pair_count))

    @property
    def element_count(self) -> int:
        return int(self.element_to_variable.size)

    @property
    def variable_count(self) -> int:
        return len(self.groups)

    def expand(self, reduced_values) -> np.ndarray:
        """Expand one value per periodic equivalence class to elements."""
        values = np.asarray(reduced_values, dtype=float).reshape(-1)
        if values.size != self.variable_count:
            raise ValueError(
                f"FFAG cyclic reduced vector has {values.size} values, "
                f"expected {self.variable_count}")
        return np.ascontiguousarray(values[self.element_to_variable])

    def contract(self, element_values) -> np.ndarray:
        """Apply the transpose map, summing element gradients/volumes."""
        values = np.asarray(element_values, dtype=float).reshape(-1)
        if values.size != self.element_count:
            raise ValueError(
                f"FFAG cyclic element vector has {values.size} values, "
                f"expected {self.element_count}")
        reduced = np.zeros(self.variable_count, dtype=float)
        np.add.at(reduced, self.element_to_variable, values)
        return reduced

    def reduce(self, element_values, *, tolerance=1.0e-12) -> np.ndarray:
        """Reduce an already periodic element vector, failing on drift."""
        values = np.asarray(element_values, dtype=float).reshape(-1)
        if values.size != self.element_count:
            raise ValueError(
                f"FFAG cyclic element vector has {values.size} values, "
                f"expected {self.element_count}")
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("FFAG cyclic density tolerance must be non-negative")
        scale = max(float(np.max(np.abs(values))), 1.0)
        reduced = np.empty(self.variable_count, dtype=float)
        for index, group in enumerate(self.groups):
            local = values[np.asarray(group, dtype=int)]
            if float(np.max(local) - np.min(local)) > tolerance * scale:
                raise ValueError(
                    f"FFAG cyclic element densities are not equal in group {index}")
            reduced[index] = float(np.mean(local))
        return reduced


def _cyclic_boundary_facets(mesh, periodic_boundaries):
    """Return paired boundary facet keys and their owning volume elements."""
    import ngsolve as ng

    boundaries = tuple(str(name) for name in periodic_boundaries)
    if len(boundaries) != 2 or boundaries[0] == boundaries[1]:
        raise ValueError("FFAG cyclic boundaries must name two distinct faces")
    available = tuple(mesh.GetBoundaries())
    missing = sorted(set(boundaries) - set(available))
    if missing:
        raise ValueError(f"FFAG cyclic boundary labels are missing: {missing}")

    facet_key_by_number = {
        int(facet.nr): frozenset(int(vertex.nr) for vertex in facet.vertices)
        for facet in mesh.facets
    }
    owner_by_facet = {}
    for element in mesh.Elements(ng.VOL):
        for facet in element.facets:
            key = facet_key_by_number[int(facet.nr)]
            owner_by_facet.setdefault(key, []).append(int(element.nr))

    facets = {name: {} for name in boundaries}
    vertices = {name: set() for name in boundaries}
    for element in mesh.Elements(ng.BND):
        name = available[element.index]
        if name not in facets:
            continue
        key = frozenset(int(vertex.nr) for vertex in element.vertices)
        owners = owner_by_facet.get(key, ())
        if len(owners) != 1:
            raise ValueError(
                "FFAG cyclic boundary facet must have one volume owner")
        facets[name][key] = owners[0]
        vertices[name].update(key)
    if any(not facets[name] for name in boundaries):
        raise ValueError("FFAG cyclic boundary has no facets")
    if vertices[boundaries[0]] & vertices[boundaries[1]]:
        raise ValueError("FFAG cyclic boundary vertex sets must be disjoint")

    master_to_slave = {}
    for endpoint_a, endpoint_b in mesh.ngmesh.GetIdentifications():
        endpoint_a = int(getattr(endpoint_a, "nr", endpoint_a)) - 1
        endpoint_b = int(getattr(endpoint_b, "nr", endpoint_b)) - 1
        if (endpoint_a in vertices[boundaries[0]]
                and endpoint_b in vertices[boundaries[1]]):
            master_to_slave[endpoint_a] = endpoint_b
        elif (endpoint_b in vertices[boundaries[0]]
                and endpoint_a in vertices[boundaries[1]]):
            master_to_slave[endpoint_b] = endpoint_a
    if (set(master_to_slave) != vertices[boundaries[0]]
            or set(master_to_slave.values()) != vertices[boundaries[1]]
            or len(set(master_to_slave.values())) != len(master_to_slave)):
        raise ValueError(
            "NGSolve PERIODIC identifications must pair every FFAG cut vertex")

    pairs = []
    used_slave = set()
    for master_key, master_owner in facets[boundaries[0]].items():
        slave_key = frozenset(master_to_slave[vertex] for vertex in master_key)
        if slave_key not in facets[boundaries[1]]:
            raise ValueError("FFAG cyclic master facet has no slave facet")
        if slave_key in used_slave:
            raise ValueError("FFAG cyclic facet pairing is not one-to-one")
        used_slave.add(slave_key)
        pairs.append((master_key, slave_key, master_owner,
                      facets[boundaries[1]][slave_key]))
    if len(used_slave) != len(facets[boundaries[1]]):
        raise ValueError("not every FFAG cyclic slave facet was paired")
    return boundaries, tuple(pairs)


def identify_ffag_cyclic_sector_vertices(
        mesh, fold, *, periodic_boundaries=("periodic_min", "periodic_max"),
        relative_tolerance=1.0e-10):
    """Add rotational PERIODIC point identifications to a Cubit sector mesh.

    Cubit owns the curved HEX geometry and named cut faces.  Netgen owns the
    periodic point-identification contract consumed by ``Periodic(HDiv)`` and
    by the broken-VIM charge pullback.  The operation is intentionally
    explicit and fails if a cut vertex cannot be matched by the first positive
    cyclic rotation.
    """
    import ngsolve as ng
    from netgen.meshing import IdentificationType

    contract = validate_ffag_cyclic_sector_contract(fold)
    relative_tolerance = float(relative_tolerance)
    if not np.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise ValueError(
            "FFAG cyclic relative_tolerance must be finite and positive")
    if mesh.ngmesh.GetIdentifications():
        raise ValueError("FFAG sector mesh already has point identifications")
    boundaries = tuple(str(name) for name in periodic_boundaries)
    available = tuple(mesh.GetBoundaries())
    missing = sorted(set(boundaries) - set(available))
    if len(boundaries) != 2 or boundaries[0] == boundaries[1] or missing:
        raise ValueError(
            "FFAG sector needs two distinct named cyclic boundaries; missing=%s"
            % missing)
    boundary_vertices = {name: set() for name in boundaries}
    for element in mesh.Elements(ng.BND):
        name = available[element.index]
        if name in boundary_vertices:
            boundary_vertices[name].update(
                int(vertex.nr) for vertex in element.vertices)
    if any(not boundary_vertices[name] for name in boundaries):
        raise ValueError("FFAG cyclic boundary has no vertices")

    angle = 2.0 * np.pi / contract.fold
    rotation = np.array(((np.cos(angle), -np.sin(angle), 0.0),
                         (np.sin(angle), np.cos(angle), 0.0),
                         (0.0, 0.0, 1.0)))
    slave = sorted(boundary_vertices[boundaries[1]])
    slave_points = np.asarray([mesh[mesh.vertices[index]].point
                               for index in slave], dtype=float)
    all_points = np.asarray([vertex.point for vertex in mesh.vertices], dtype=float)
    scale = max(float(np.max(np.linalg.norm(
        all_points - np.mean(all_points, axis=0), axis=1))), 1.0e-300)
    limit = float(relative_tolerance) * scale
    used = set()
    pairs = []
    maximum_residual = 0.0
    for master in sorted(boundary_vertices[boundaries[0]]):
        point = np.asarray(mesh[mesh.vertices[master]].point, dtype=float)
        distances = np.linalg.norm(slave_points - rotation @ point, axis=1)
        index = int(np.argmin(distances))
        residual = float(distances[index])
        target = slave[index]
        if residual > limit or target in used:
            raise ValueError(
                "FFAG cyclic vertex rotation mismatch "
                f"(residual={residual:g}, limit={limit:g})")
        used.add(target)
        maximum_residual = max(maximum_residual, residual)
        mesh.ngmesh.AddPointIdentification(
            master + 1, target + 1, identnr=1,
            type=IdentificationType.PERIODIC)
        pairs.append((master, target))
    if len(used) != len(slave):
        raise ValueError("not every FFAG cyclic slave vertex was identified")
    return {
        "pair_count": len(pairs),
        "maximum_rotation_residual": maximum_residual,
        "relative_rotation_residual": maximum_residual / scale,
        "periodic_boundaries": boundaries,
    }


def build_ffag_cyclic_density_map(
        mesh, *, periodic_boundaries=("periodic_min", "periodic_max")):
    """Tie topology-density variables adjacent to identified cyclic faces."""
    boundaries, facet_pairs = _cyclic_boundary_facets(
        mesh, periodic_boundaries)
    parent = np.arange(int(mesh.ne), dtype=np.int64)

    def find(value):
        value = int(value)
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(first, second):
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[root_second] = root_first

    for _master, _slave, master_owner, slave_owner in facet_pairs:
        union(master_owner, slave_owner)
    roots = [find(index) for index in range(int(mesh.ne))]
    root_to_variable = {}
    mapping = np.empty(int(mesh.ne), dtype=np.int64)
    for element, root in enumerate(roots):
        mapping[element] = root_to_variable.setdefault(
            root, len(root_to_variable))
    groups = tuple(tuple(np.flatnonzero(mapping == index).tolist())
                   for index in range(len(root_to_variable)))
    return FFAGCyclicDensityMap(
        mapping, groups, len(facet_pairs), boundaries)


def validate_ffag_cyclic_sector_contract(
        fold, *, field_antiperiodic=False,
        body_crosses_periodic_planes=False,
        formulation="fem", periodic_trace_identified=False,
        periodic_charge_paired=False) -> FFAGCyclicSectorContract:
    """Validate the cyclic symmetry contract before an FFAG sector solve.

    FFAG F/D alternation normally changes the radial gradient within a cell;
    it does not make the vertical guide field antiperiodic from cell to cell.
    ``field_antiperiodic`` is therefore an explicit exceptional mode and, as
    required by closure around the ring, accepts only an even fold count.
    """
    numeric_fold = float(fold)
    if (isinstance(fold, (bool, np.bool_)) or not np.isfinite(numeric_fold)
            or not numeric_fold.is_integer() or numeric_fold < 2):
        raise ValueError("FFAG cyclic fold must be an integer >= 2")
    count = int(numeric_fold)
    antiperiodic = bool(field_antiperiodic)
    crosses = bool(body_crosses_periodic_planes)
    formulation = str(formulation).strip().lower().replace("_", "-")
    if formulation not in {"fem", "vim-broken"}:
        raise ValueError("FFAG cyclic formulation must be 'fem' or 'vim-broken'")
    identified = bool(periodic_trace_identified)
    paired = bool(periodic_charge_paired)
    if antiperiodic and count % 2:
        raise ValueError(
            "FFAG antiperiodic field symmetry requires an even fold count")
    if crosses and formulation == "fem" and not identified:
        raise ValueError(
            "continuous FFAG FEM iron crosses the azimuthal sector planes: "
            "identify the two rotation-related HDiv normal traces before "
            "using image_cyclic")
    if crosses and formulation == "vim-broken" and not paired:
        raise ValueError(
            "continuous FFAG broken-HDiv VIM iron crosses the azimuthal "
            "sector planes: pair the two rotation-related surface-charge "
            "rows before using image_cyclic")
    if identified and (not crosses or formulation != "fem"):
        raise ValueError(
            "periodic_trace_identified is only valid for connected FEM")
    if paired and (not crosses or formulation != "vim-broken"):
        raise ValueError(
            "periodic_charge_paired is only valid for connected broken-HDiv "
            "VIM")
    mode = (
        "connected-periodic-fem-sector"
        if crosses and formulation == "fem" else
        "connected-periodic-vim-sector"
        if crosses else "disjoint-cell-cyclic-images")
    return FFAGCyclicSectorContract(
        count, antiperiodic, crosses, formulation, identified, paired, mode)


def magnetic_rigidity_from_kinetic_energy(
        kinetic_energy_mev, *, rest_energy_mev=PROTON_REST_ENERGY_MEV,
        charge_number=1.0):
    """Return relativistic magnetic rigidity ``B rho`` in tesla-metre."""
    kinetic = np.asarray(kinetic_energy_mev, dtype=float)
    rest = float(rest_energy_mev)
    charge = abs(float(charge_number))
    if (not np.all(np.isfinite(kinetic)) or np.any(kinetic < 0.0)
            or not np.isfinite(rest) or rest <= 0.0
            or not np.isfinite(charge) or charge <= 0.0):
        raise ValueError(
            "kinetic energy must be nonnegative and rest energy/charge "
            "must be positive")
    momentum_gev_c = np.sqrt(
        kinetic * (kinetic + 2.0 * rest)) / 1000.0
    result = momentum_gev_c / (GEV_C_PER_TESLA_METRE * charge)
    return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class EngeFringeIntegrals:
    """Equal-integral boundary and dimensionless Enge form factors."""

    effective_boundary_m: float
    i1: float
    i2: float
    equal_integral_residual: float
    full_gap_m: float


def enge_fringe_integrals(
        coordinate_m, bending_field_t, *, body_field_t, full_gap_m
        ) -> EngeFringeIntegrals:
    """Evaluate exit-fringe ``I1`` and ``I2`` from a sampled field.

    Samples must run from the constant-field side to the zero-field side.
    ``sigma=(z-z_eff)/g`` uses the full gap ``g``.  The effective boundary
    ``z_eff`` equates the soft-field integral to a sharp unit step.  With
    ``f=B/B0`` and ``q=H(-sigma)-f``, the numerically stable first-moment form

    ``I1 = - integral sigma*q d sigma``

    is equivalent to the conventional nested integral when
    ``integral q d sigma=0``.  ``I2=integral f*(1-f) d sigma``.
    """
    coordinate = np.asarray(coordinate_m, dtype=float).reshape(-1)
    field = np.asarray(bending_field_t, dtype=float).reshape(-1)
    body = float(body_field_t)
    gap = float(full_gap_m)
    if (coordinate.size < 8 or field.shape != coordinate.shape
            or not np.all(np.isfinite(np.r_[coordinate, field, body, gap]))
            or np.any(np.diff(coordinate) <= 0.0) or body == 0.0
            or gap <= 0.0):
        raise ValueError(
            "Enge integration needs at least eight ordered finite samples, "
            "nonzero body field, and positive full gap")
    normalized = field / body
    end_tolerance = 5.0e-3
    if (abs(normalized[0] - 1.0) > end_tolerance
            or abs(normalized[-1]) > end_tolerance):
        raise ValueError(
            "samples must start in the body field and end in the zero field")
    if np.min(normalized) < -0.25 or np.max(normalized) > 1.25:
        raise ValueError(
            "normalized fringe field is too far outside the body/zero range")
    integral = float(np.trapezoid(normalized, coordinate))
    effective = float(coordinate[0] + integral)
    sigma = (coordinate - effective) / gap
    # Integrate the discontinuous sharp step analytically.  Sampling it with
    # a trapezoid would leave a grid-dependent half-cell in both I1 and the
    # equal-integral residual when z_eff happens to be a sample location.
    hard_moment = -0.5 * float(sigma[0] * sigma[0])
    soft_moment = float(np.trapezoid(sigma * normalized, sigma))
    i1 = soft_moment - hard_moment
    i2 = float(np.trapezoid(normalized * (1.0 - normalized), sigma))
    residual = float(
        ((effective - coordinate[0]) - integral) / gap)
    return EngeFringeIntegrals(
        effective, i1, i2, residual, gap)


def _tanh_window(coordinate, start, stop, epsilon):
    return 0.5 * (
        np.tanh((coordinate - start) / epsilon)
        - np.tanh((coordinate - stop) / epsilon))


@dataclass(frozen=True)
class FFAGSoftEdgeCellSpec:
    """Periodic doublet-cell parameters for a soft-edge FFAG PoC.

    The default constructor is generic.  :meth:`bell_abell` supplies Table 1
    of Bell and Abell, arXiv:1202.0805.  ``full_gap_m`` is deliberately a PoC
    input because that paper states that the 5 cm Enge scale is similar to the
    aperture but does not publish a full pole gap.
    """

    cell_count: int
    long_drift_m: float
    defocusing_length_m: float
    short_drift_m: float
    focusing_length_m: float
    defocusing_b0_t: float
    defocusing_gradient_t_per_m: float
    focusing_b0_t: float
    focusing_gradient_t_per_m: float
    fringe_epsilon_m: float
    full_gap_m: float

    def __post_init__(self):
        count = int(self.cell_count)
        lengths = np.asarray([
            self.long_drift_m, self.defocusing_length_m,
            self.short_drift_m, self.focusing_length_m,
            self.fringe_epsilon_m, self.full_gap_m], dtype=float)
        fields = np.asarray([
            self.defocusing_b0_t, self.defocusing_gradient_t_per_m,
            self.focusing_b0_t, self.focusing_gradient_t_per_m], dtype=float)
        if (count < 2 or not np.all(np.isfinite(lengths))
                or np.any(lengths <= 0.0) or not np.all(np.isfinite(fields))
                or (self.defocusing_gradient_t_per_m
                    + self.focusing_gradient_t_per_m) == 0.0):
            raise ValueError("invalid FFAG soft-edge cell specification")
        object.__setattr__(self, "cell_count", count)

    @classmethod
    def bell_abell(cls, *, full_gap_m=0.10):
        """Return the published 24-cell, 31--250 MeV proton fixture."""
        return cls(
            cell_count=24,
            long_drift_m=0.40,
            defocusing_length_m=0.22,
            short_drift_m=0.075,
            focusing_length_m=0.44,
            defocusing_b0_t=0.803952,
            defocusing_gradient_t_per_m=-12.8,
            focusing_b0_t=0.555057,
            focusing_gradient_t_per_m=8.0,
            fringe_epsilon_m=0.05,
            full_gap_m=full_gap_m)

    @property
    def cell_length_m(self) -> float:
        return float(
            self.long_drift_m + self.defocusing_length_m
            + self.short_drift_m + self.focusing_length_m)

    @property
    def cell_bend_angle_rad(self) -> float:
        return float(2.0 * np.pi / self.cell_count)

    @property
    def magnet_intervals_m(self):
        bd_start = 0.5 * self.long_drift_m
        bd_stop = bd_start + self.defocusing_length_m
        bf_start = bd_stop + self.short_drift_m
        bf_stop = bf_start + self.focusing_length_m
        return ((bd_start, bd_stop), (bf_start, bf_stop))

    def sampled_profiles(self, *, n_segments=256, periodic_images=2):
        """Return midpoint ``s,ds,B0(s),G(s)`` with overlapping fringes."""
        count = int(n_segments)
        images = int(periodic_images)
        if count < 16 or images < 1:
            raise ValueError(
                "FFAG profile needs at least 16 segments and one image")
        edges = np.linspace(0.0, self.cell_length_m, count + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        lengths = np.diff(edges)
        bd = np.zeros(count)
        bf = np.zeros(count)
        for image in range(-images, images + 1):
            shift = image * self.cell_length_m
            (bd_start, bd_stop), (bf_start, bf_stop) = (
                self.magnet_intervals_m)
            bd += _tanh_window(
                centers, bd_start + shift, bd_stop + shift,
                self.fringe_epsilon_m)
            bf += _tanh_window(
                centers, bf_start + shift, bf_stop + shift,
                self.fringe_epsilon_m)
        b0 = (self.defocusing_b0_t * bd
              + self.focusing_b0_t * bf)
        gradient = (self.defocusing_gradient_t_per_m * bd
                    + self.focusing_gradient_t_per_m * bf)
        return (np.ascontiguousarray(centers),
                np.ascontiguousarray(lengths),
                np.ascontiguousarray(b0),
                np.ascontiguousarray(gradient))

    def symmetric_tanh_fringe_integrals(self, *, sample_count=4001):
        """Return the isolated tanh edge diagnostic in the declared gap."""
        count = int(sample_count)
        if count < 101:
            raise ValueError("sample_count must be at least 101")
        extent = 12.0 * self.fringe_epsilon_m
        coordinate = np.linspace(-extent, extent, count)
        normalized = 0.5 * (
            1.0 - np.tanh(coordinate / self.fringe_epsilon_m))
        return enge_fringe_integrals(
            coordinate, normalized, body_field_t=1.0,
            full_gap_m=self.full_gap_m)


def _periodic_planar_orbit(curvature, segment_lengths, rigidity, bend_axis):
    curvature = np.asarray(curvature, dtype=float).reshape(-1)
    lengths = np.asarray(segment_lengths, dtype=float).reshape(-1)
    if curvature.shape != lengths.shape:
        raise ValueError("curvature and segment lengths must match")
    angle = np.r_[0.0, np.cumsum(curvature * lengths)]
    turning = np.diff(angle)
    midpoint_angle = angle[:-1] + 0.5 * turning
    steps = lengths[:, None] * np.column_stack((
        np.cos(midpoint_angle), np.sin(midpoint_angle)))
    relative = np.vstack((np.zeros(2), np.cumsum(steps, axis=0)))
    total_angle = float(angle[-1])
    rotation = np.array([
        [np.cos(total_angle), -np.sin(total_angle)],
        [np.sin(total_angle), np.cos(total_angle)],
    ])
    start = np.linalg.solve(rotation - np.eye(2), relative[-1])
    positions_2d = relative + start
    # All momenta must cross the same radial cell-boundary plane.  Rotate the
    # otherwise arbitrary local solution so its entrance lies on -y.  The
    # entrance tangent is allowed to vary with momentum, as it does in an FFAG.
    entrance_angle = float(np.arctan2(
        positions_2d[0, 1], positions_2d[0, 0]))
    alignment = -0.5 * np.pi - entrance_angle
    align_rotation = np.array([
        [np.cos(alignment), -np.sin(alignment)],
        [np.sin(alignment), np.cos(alignment)],
    ])
    positions_2d = positions_2d @ align_rotation.T
    tangent_2d = np.column_stack((np.cos(angle), np.sin(angle)))
    tangent_2d = tangent_2d @ align_rotation.T
    positions = np.column_stack((positions_2d, np.zeros(len(positions_2d))))
    tangents = np.column_stack((tangent_2d, np.zeros(len(angle))))
    return PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=float(rigidity),
        bend_axis=np.asarray(bend_axis, dtype=float),
        path_length_stations=np.r_[0.0, np.cumsum(lengths)])


@dataclass(frozen=True)
class FFAGCellReference:
    """One periodic reference orbit and its soft-edge first-order map."""

    kinetic_energy_mev: float
    magnetic_rigidity_tm: float
    transverse_offset_m: float
    orbit: PlanarDesignOrbit
    field_response: np.ndarray
    transfer: CombinedFunctionTransferMap
    bend_angle_rad: float
    periodic_position_residual_m: float
    periodic_tangent_residual: float


def build_ffag_cell_reference(
        spec: FFAGSoftEdgeCellSpec, kinetic_energy_mev, *,
        n_segments=256, response_entries=None) -> FFAGCellReference:
    """Construct the periodic reduced closed orbit for one kinetic energy.

    The Bell--Abell field law is linear in transverse displacement.  The
    single cell-wide displacement is therefore eliminated analytically from
    the total-bend condition.  This is a deterministic soft-edge target
    fixture.  A realized 3-D HDiv field must subsequently recover its own
    closed orbit; the reduced construction is not used as an acceptance
    substitute.
    """
    if not isinstance(spec, FFAGSoftEdgeCellSpec):
        raise TypeError("spec must be an FFAGSoftEdgeCellSpec")
    energy = float(kinetic_energy_mev)
    rigidity = magnetic_rigidity_from_kinetic_energy(energy)
    _, lengths, b0, gradient = spec.sampled_profiles(
        n_segments=n_segments)
    b0_integral = float(b0 @ lengths)
    gradient_integral = float(gradient @ lengths)
    scale = max(1.0, abs(b0_integral))
    if abs(gradient_integral) <= 1.0e-12 * scale:
        raise RuntimeError(
            "cell gradient integral cannot set the reference-orbit offset")
    offset = (
        rigidity * spec.cell_bend_angle_rad - b0_integral
    ) / gradient_integral
    field = b0 + gradient * offset
    curvature = field / rigidity
    orbit = _periodic_planar_orbit(
        curvature, lengths, rigidity, np.array([0.0, 0.0, 1.0]))
    raw = np.r_[field, gradient]
    transfer = combined_function_transfer_map_from_field_response(
        raw, lengths, rigidity, response_entries=response_entries)
    theta = spec.cell_bend_angle_rad
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    position_residual = float(np.linalg.norm(
        orbit.positions[-1] - rotation @ orbit.positions[0]))
    tangent_residual = float(np.linalg.norm(
        orbit.tangents[-1] - rotation @ orbit.tangents[0]))
    bend_angle = float(np.sum(curvature * lengths))
    return FFAGCellReference(
        energy, rigidity, float(offset), orbit,
        np.ascontiguousarray(raw), transfer, bend_angle,
        position_residual, tangent_residual)


@dataclass(frozen=True)
class FFAGCellTargetFamily:
    """Momentum-indexed reduced cell targets ready for HDiv-MMM fusion."""

    spec: FFAGSoftEdgeCellSpec
    references: tuple[FFAGCellReference, ...]
    objective: MultiMomentumTransferMatrixObjective
    fringe_integrals: EngeFringeIntegrals

    @property
    def kinetic_energies_mev(self) -> np.ndarray:
        return np.asarray([
            reference.kinetic_energy_mev for reference in self.references])


@dataclass(frozen=True)
class FFAGFixedDesignOrbitTargetFamily:
    """Caller-supplied one-pass design orbits and target transfer maps.

    This is the target contract for a beam line or one isolated magnet.  It has
    no soft-edge field fixture and no periodic-orbit reconstruction: the
    supplied :class:`PlanarDesignOrbit` objects define the observation paths,
    and ``objective`` owns the requested maps and engineering bands about those
    paths.
    """

    objective: MultiMomentumTransferMatrixObjective
    controlled_components: tuple[str, ...] = ()
    target_symplectic_residuals: np.ndarray | None = None

    def __post_init__(self):
        if not isinstance(
                self.objective, MultiMomentumTransferMatrixObjective):
            raise TypeError(
                "objective must be a MultiMomentumTransferMatrixObjective")
        components = tuple(str(value) for value in self.controlled_components)
        residuals = self.target_symplectic_residuals
        if residuals is None:
            residuals = np.asarray([
                static_magnet_symplectic_residual(matrix)
                for matrix in self.objective.target_matrices])
        else:
            residuals = np.asarray(residuals, dtype=float).reshape(-1)
        if (residuals.shape != (len(self.objective.orbits),)
                or not np.all(np.isfinite(residuals))
                or np.any(residuals < 0.0)):
            raise ValueError(
                "target_symplectic_residuals must contain one finite "
                "nonnegative value per design orbit")
        object.__setattr__(self, "controlled_components", components)
        object.__setattr__(
            self, "target_symplectic_residuals",
            np.ascontiguousarray(residuals))

    @property
    def design_orbits(self) -> tuple[PlanarDesignOrbit, ...]:
        return self.objective.orbits

    @property
    def magnetic_rigidities_tm(self) -> np.ndarray:
        return np.asarray([
            orbit.magnetic_rigidity for orbit in self.objective.orbits])


def build_ffag_cell_target_family(
        kinetic_energies_mev, *, spec=None, n_segments=256,
        transfer_matrix_band=1.0e-3, bend_field_band=1.0e-3,
        response_entries=None) -> FFAGCellTargetFamily:
    """Create a multi-momentum cell objective from soft-edge FFAG targets."""
    if spec is None:
        spec = FFAGSoftEdgeCellSpec.bell_abell()
    energies = np.asarray(kinetic_energies_mev, dtype=float).reshape(-1)
    if (energies.size < 2 or not np.all(np.isfinite(energies))
            or np.any(energies <= 0.0)
            or np.any(np.diff(energies) <= 0.0)):
        raise ValueError(
            "FFAG target family needs at least two increasing positive "
            "kinetic energies")
    references = tuple(build_ffag_cell_reference(
        spec, energy, n_segments=n_segments,
        response_entries=response_entries) for energy in energies)
    entries = (references[0].transfer.response_entries
               if response_entries is None else tuple(response_entries))
    objective = MultiMomentumTransferMatrixObjective(
        tuple(reference.orbit for reference in references),
        np.asarray([reference.transfer.matrix for reference in references]),
        transfer_matrix_band, bend_field_band, entries)
    return FFAGCellTargetFamily(
        spec, references, objective,
        spec.symmetric_tanh_fringe_integrals())


def build_ffag_fixed_design_orbit_target_family(
        design_orbits, target_transfer_matrices, *,
        transfer_matrix_band=1.0e-3, bend_field_band=1.0e-3,
        response_entries=None, controlled_components=None,
        require_symplectic=True, symplectic_tolerance=1.0e-9,
        curvature_sign=1.0, gradient_sign=1.0
        ) -> FFAGFixedDesignOrbitTargetFamily:
    """Build the direct ``design orbit + target map`` one-pass contract.

    ``target_transfer_matrices[i]`` is interpreted about
    ``design_orbits[i]``.  The orbit geometry, magnetic rigidity, target map,
    and bands are all caller-owned; no Enge profile, reference-field fixture,
    or closed-orbit solve is inserted by this constructor.  Named
    ``controlled_components`` select physically interpretable entries such as
    focusing and horizontal dispersion, but they do not make the unselected
    entries independent.  By default each complete target map must satisfy the
    static-magnet symplectic condition before any material optimization starts.
    """
    orbits = tuple(design_orbits)
    if response_entries is not None and controlled_components is not None:
        raise ValueError(
            "response_entries and controlled_components are mutually "
            "exclusive")
    components = ()
    if controlled_components is not None:
        if isinstance(controlled_components, str):
            controlled_components = (controlled_components,)
        components = tuple(
            str(value).strip().lower().replace("-", "_")
            for value in controlled_components)
        response_entries = static_magnet_transfer_component_entries(
            components)
    options = dict(
        transfer_matrix_band=transfer_matrix_band,
        bend_field_band=bend_field_band,
        curvature_sign=curvature_sign,
        gradient_sign=gradient_sign)
    if response_entries is not None:
        options["response_entries"] = response_entries
    objective = MultiMomentumTransferMatrixObjective(
        orbits, target_transfer_matrices, **options)
    residuals = np.asarray([
        static_magnet_symplectic_residual(matrix)
        for matrix in objective.target_matrices])
    tolerance = float(symplectic_tolerance)
    if (not np.isfinite(tolerance) or tolerance < 0.0):
        raise ValueError(
            "symplectic_tolerance must be finite and nonnegative")
    if require_symplectic and np.any(residuals > tolerance):
        index = int(np.argmax(residuals))
        raise ValueError(
            "target transfer matrix is not symplectic in static-magnet "
            f"coordinates: orbit {index}, residual {residuals[index]:.6g}, "
            f"tolerance {tolerance:.6g}")
    return FFAGFixedDesignOrbitTargetFamily(
        objective, components, residuals)


def _evaluate_b_field(field, points) -> np.ndarray:
    """Evaluate a callable or ``b_field`` provider in tesla."""
    values = np.asarray(points, dtype=float)
    single = values.shape == (3,)
    points_2d = np.ascontiguousarray(values.reshape(-1, 3))
    if not np.all(np.isfinite(points_2d)):
        raise ValueError("field-evaluation points must be finite")
    evaluator = getattr(field, "b_field", None)
    if evaluator is not None:
        result = np.asarray(evaluator(points_2d), dtype=float)
    else:
        try:
            result = np.asarray(field(points_2d), dtype=float)
        except (TypeError, ValueError):
            result = np.asarray([
                field(float(point[0]), float(point[1]), float(point[2]))
                for point in points_2d], dtype=float)
    if result.shape == (3,) and len(points_2d) == 1:
        result = result[None, :]
    if result.shape != points_2d.shape or not np.all(np.isfinite(result)):
        raise ValueError(
            "magnetic field provider must return one finite 3-vector per "
            "point")
    return result[0] if single else result


def sample_planar_orbit_field_response(
        field, orbit: PlanarDesignOrbit, *, gradient_offset) -> np.ndarray:
    """Sample ``[B_binormal, dB_binormal/dnormal]`` on an orbit.

    The centered normal stencil defines the physical field observable.  It is
    not a design finite difference: whole-element HDiv-MMM sensitivities still
    use the analytic Schur/adjoint contraction of the corresponding rows.
    """
    points, weights = planar_orbit_field_observations(
        orbit, gradient_offset=gradient_offset)
    response = np.einsum(
        "rpc,pc->r", weights, _evaluate_b_field(field, points))
    return np.ascontiguousarray(response, dtype=float)


@dataclass(frozen=True)
class FullFieldClosedOrbit:
    """Periodic planar reference orbit recovered from a realized 3-D field."""

    magnetic_rigidity_tm: float
    orbit: PlanarDesignOrbit
    path_length_m: float
    entrance_radius_m: float
    entrance_incidence_angle_rad: float
    periodic_position_residual_m: float
    periodic_tangent_residual: float
    vertical_position_residual_m: float
    vertical_tangent_residual: float
    root_evaluations: int
    field_response: np.ndarray
    transfer: CombinedFunctionTransferMap

    @property
    def closure_residual(self) -> float:
        return float(max(
            self.periodic_position_residual_m,
            self.periodic_tangent_residual,
            abs(self.vertical_position_residual_m),
            abs(self.vertical_tangent_residual)))


@dataclass(frozen=True)
class ReclosedOrbitShapeJacobian:
    """Total analytic GetTrafo derivative about a reclosed design orbit.

    The HDiv state/source derivative is analytic.  The trajectory tangent is
    propagated through the Lorentz variational equation and the entrance
    radius/incidence response follows from the 2-by-2 periodic-closure
    implicit function.  Hence the returned transfer derivative includes the
    moving observation path and segment lengths, not merely a lagged orbit.
    """

    recovered: FullFieldClosedOrbit
    entrance_parameter_jacobian: np.ndarray
    path_length_jacobian: np.ndarray
    position_jacobian: np.ndarray
    tangent_jacobian: np.ndarray
    field_response: np.ndarray
    field_response_jacobian: np.ndarray
    transfer: CombinedFunctionTransferMap
    closure_partial_jacobian: np.ndarray
    observation_points: np.ndarray
    fixed_field_shape_jacobian: np.ndarray
    observation_point_jacobian: np.ndarray


def differentiate_recovered_planar_orbit_shape_native(
        recovered: FullFieldClosedOrbit, *, charge_gram, state,
        state_jacobian, cell_vertex_velocity, face_vertex_velocity,
        iron_evaluator, iron_scale, constant_field_t=(0.0, 0.0, 0.0),
        cell_angle_rad, gradient_offset=1.0e-3,
        tracking_step_m=5.0e-3, integration_stations=257,
        maximum_path_m=None, curvature_sign=1.0,
        response_entries=None) -> ReclosedOrbitShapeJacobian:
    """Differentiate a periodic orbit and its map without design FD.

    This flat-TET production lane uses the exact native HDiv field spatial
    Jacobian and ``C*dm+dC*m`` source tangent.  A dense nominal native track
    supplies the coefficients of the Lorentz variational equation.  Its
    entrance radius/incidence derivatives are closed analytically with the
    implicit-function theorem.  Radia-object and extra tracker mirror terms
    are intentionally absent: every field term in this derivative must expose
    the same analytic spatial/design tangent contract.
    """
    from . import _radia_pybind as _native

    if not isinstance(recovered, FullFieldClosedOrbit):
        raise TypeError("recovered must be a FullFieldClosedOrbit")
    state = np.asarray(state, dtype=float).reshape(-1)
    state_jacobian = np.asarray(state_jacobian, dtype=float)
    cells = np.ascontiguousarray(cell_vertex_velocity, dtype=float)
    faces = np.ascontiguousarray(face_vertex_velocity, dtype=float)
    modes = state_jacobian.shape[0] if state_jacobian.ndim == 2 else 0
    constant_field = np.asarray(constant_field_t, dtype=float).reshape(-1)
    theta = float(cell_angle_rad)
    offset = float(gradient_offset)
    step = float(tracking_step_m)
    station_count = int(integration_stations)
    rigidity = float(recovered.magnetic_rigidity_tm)
    iron_scale = float(iron_scale)
    curvature_sign = float(curvature_sign)
    segment_count = len(recovered.orbit.segment_lengths)
    if (modes < 1 or state_jacobian.shape[1:] != state.shape
            or cells.shape[0] != modes or faces.shape[0] != modes
            or cells.ndim != 4 or cells.shape[2:] != (4, 3)
            or faces.ndim != 4 or faces.shape[2:] != (3, 3)
            or constant_field.shape != (3,)
            or not np.all(np.isfinite(np.r_[
                state, state_jacobian.ravel(), cells.ravel(), faces.ravel(),
                constant_field, theta, offset, step, rigidity, iron_scale,
                curvature_sign]))
            or theta <= 0.0 or theta >= np.pi or offset <= 0.0
            or step <= 0.0 or rigidity <= 0.0 or iron_scale == 0.0
            or curvature_sign == 0.0 or station_count < 17
            or (station_count-1) % segment_count != 0):
        raise ValueError("invalid reclosed-orbit shape derivative settings")
    if maximum_path_m is None:
        maximum_path = max(1.5*recovered.path_length_m,
                           recovered.path_length_m+4.0*step)
    else:
        maximum_path = float(maximum_path_m)
    if not np.isfinite(maximum_path) or maximum_path <= recovered.path_length_m:
        raise ValueError("maximum_path_m must exceed the recovered path length")

    axis = np.array([0.0, 0.0, 1.0])
    radial_0 = np.array([0.0, -1.0, 0.0])
    tangent_0 = np.array([1.0, 0.0, 0.0])
    cosine = float(np.cos(theta))
    sine = float(np.sin(theta))
    rotation = np.array([
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ])
    radial_1 = rotation@radial_0
    tangent_1 = np.ascontiguousarray(rotation@tangent_0)
    radius = float(recovered.entrance_radius_m)
    alpha = float(recovered.entrance_incidence_angle_rad)
    entrance_tangent = (
        np.cos(alpha)*tangent_0 + np.sin(alpha)*radial_0)
    entrance_alpha_tangent = (
        -np.sin(alpha)*tangent_0 + np.cos(alpha)*radial_0)
    native_rigidity = -rigidity/curvature_sign
    tracked = _native.track_reference_orbit_to_plane_native(
        iron_evaluator, iron_scale, -1, False, constant_field,
        native_rigidity, np.ascontiguousarray(radius*radial_0),
        np.ascontiguousarray(entrance_tangent), tangent_1, 0.0,
        step, maximum_path, 1.0e-6, station_count, "direct")
    positions = np.ascontiguousarray(tracked[0], dtype=float)
    tangents = np.ascontiguousarray(tracked[1], dtype=float)
    stations = np.ascontiguousarray(tracked[2], dtype=float)
    length = float(tracked[4])
    if abs(length-recovered.path_length_m) > 5.0e-10*max(1.0, length):
        raise RuntimeError(
            "dense shape track does not reproduce the recovered path length")

    raw_field = np.asarray(
        iron_evaluator.field(positions, "direct"), dtype=float)
    field = iron_scale*raw_field + constant_field[None, :]
    field_gradient = iron_scale*np.asarray(
        iron_evaluator.field_gradient(positions), dtype=float)
    # The configured-value API includes 1/(4*pi), whereas the persistent
    # evaluator is deliberately raw.  Convert both to the same B units.
    fixed_shape_field = (4.0*np.pi*iron_scale)*np.asarray(
        charge_gram.configured_field_values_shape_derivative(
            positions, np.ascontiguousarray(state),
            np.ascontiguousarray(state_jacobian), cells, faces),
        dtype=float)
    if (field.shape != (station_count, 3)
            or field_gradient.shape != (station_count, 3, 3)
            or fixed_shape_field.shape != (modes, station_count, 3)):
        raise RuntimeError("native field tangent API returned invalid shapes")

    partial_count = modes+2
    sensitivities = np.zeros((station_count, partial_count, 6))
    sensitivities[0, 0, :3] = radial_0
    sensitivities[0, 1, 3:] = entrance_alpha_tangent

    def sensitivity_rhs(values, tangent, b_value, gradient, shape_value):
        result = np.empty_like(values)
        result[:, :3] = values[:, 3:]
        displaced_field = values[:, :3]@gradient.T
        displaced_field[2:] += shape_value
        result[:, 3:] = (
            np.cross(values[:, 3:], b_value[None, :])
            + np.cross(tangent[None, :], displaced_field)
        )/native_rigidity
        return result

    for index in range(station_count-1):
        ds = stations[index+1]-stations[index]
        t_left = tangents[index]
        t_right = tangents[index+1]
        t_middle = t_left+t_right
        t_middle /= np.linalg.norm(t_middle)
        b_left = field[index]
        b_right = field[index+1]
        b_middle = 0.5*(b_left+b_right)
        g_left = field_gradient[index]
        g_right = field_gradient[index+1]
        g_middle = 0.5*(g_left+g_right)
        q_left = fixed_shape_field[:, index]
        q_right = fixed_shape_field[:, index+1]
        q_middle = 0.5*(q_left+q_right)
        value = sensitivities[index]
        k1 = sensitivity_rhs(value, t_left, b_left, g_left, q_left)
        k2 = sensitivity_rhs(
            value+0.5*ds*k1, t_middle, b_middle, g_middle, q_middle)
        k3 = sensitivity_rhs(
            value+0.5*ds*k2, t_middle, b_middle, g_middle, q_middle)
        k4 = sensitivity_rhs(
            value+ds*k3, t_right, b_right, g_right, q_right)
        sensitivities[index+1] = value + ds*(k1+2*k2+2*k3+k4)/6.0

    exit_tangent = tangents[-1]
    denominator = float(tangent_1@exit_tangent)
    if abs(denominator) <= 1.0e-8:
        raise RuntimeError("exit-plane event is tangent to the recovered orbit")
    partial_length = -(sensitivities[-1, :, :3]@tangent_1)/denominator
    exit_acceleration = np.cross(exit_tangent, field[-1])/native_rigidity
    event_position = (sensitivities[-1, :, :3]
                      + partial_length[:, None]*exit_tangent)
    event_tangent = (sensitivities[-1, :, 3:]
                     + partial_length[:, None]*exit_acceleration)
    back_position = event_position@rotation
    back_tangent = event_tangent@rotation
    entrance_partials = np.zeros((partial_count, 3))
    entrance_partials[1] = entrance_alpha_tangent
    cross_value = float(axis@np.cross(
        entrance_tangent, rotation.T@exit_tangent))
    dot_value = float(entrance_tangent@(rotation.T@exit_tangent))
    closure_jacobian = np.empty((2, partial_count))
    radius_scale = max(radius, 1.0e-6)
    closure_jacobian[0] = (back_position@radial_0)/radius_scale
    closure_jacobian[0, 0] -= 1.0/radius_scale
    for column in range(partial_count):
        dcross = float(axis@(
            np.cross(entrance_partials[column], rotation.T@exit_tangent)
            + np.cross(entrance_tangent, back_tangent[column])))
        ddot = float(
            entrance_partials[column]@(rotation.T@exit_tangent)
            + entrance_tangent@back_tangent[column])
        closure_jacobian[1, column] = (
            dot_value*dcross-cross_value*ddot)/(dot_value**2+cross_value**2)
    entrance_jacobian = -np.linalg.solve(
        closure_jacobian[:, :2], closure_jacobian[:, 2:])
    length_jacobian = (
        partial_length[2:] + partial_length[:2]@entrance_jacobian)

    fractions = stations/length
    acceleration = np.cross(tangents, field)/native_rigidity
    normalized_sensitivity = sensitivities.copy()
    normalized_sensitivity[:, :, :3] += (
        fractions[:, None, None]*partial_length[None, :, None]
        * tangents[:, None, :])
    normalized_sensitivity[:, :, 3:] += (
        fractions[:, None, None]*partial_length[None, :, None]
        * acceleration[:, None, :])
    total_sensitivity = (
        normalized_sensitivity[:, 2:]
        + np.einsum(
            "sui,uq->sqi", normalized_sensitivity[:, :2],
            entrance_jacobian))
    total_sensitivity[:, :, 3:] -= (
        tangents[:, None, :]
        * np.einsum("si,sqi->sq", tangents,
                    total_sensitivity[:, :, 3:])[:, :, None])

    stride = (station_count-1)//segment_count
    indices = np.arange(0, station_count, stride, dtype=np.int64)
    coarse_position = positions[indices]
    coarse_tangent = tangents[indices]
    position_jacobian = np.transpose(
        total_sensitivity[indices, :, :3], (1, 0, 2))
    tangent_jacobian = np.transpose(
        total_sensitivity[indices, :, 3:], (1, 0, 2))
    segment_lengths = np.full(segment_count, length/segment_count)
    segment_length_jacobian = np.tile(
        length_jacobian[None, :]/segment_count, (segment_count, 1))
    centers = (
        0.5*(coarse_position[:-1]+coarse_position[1:])
        + segment_lengths[:, None]
        *(coarse_tangent[:-1]-coarse_tangent[1:])/8.0)
    center_jacobian = (
        0.5*(position_jacobian[:, :-1]+position_jacobian[:, 1:])
        + segment_length_jacobian.T[:, :, None]
        *(coarse_tangent[:-1]-coarse_tangent[1:])[None, :, :]/8.0
        + segment_lengths[None, :, None]
        *(tangent_jacobian[:, :-1]-tangent_jacobian[:, 1:])/8.0)
    tangent_sum = coarse_tangent[:-1]+coarse_tangent[1:]
    midpoint_tangent = tangent_sum/np.linalg.norm(
        tangent_sum, axis=1)[:, None]
    tangent_sum_jacobian = tangent_jacobian[:, :-1]+tangent_jacobian[:, 1:]
    midpoint_tangent_jacobian = (
        tangent_sum_jacobian
        - midpoint_tangent[None, :, :]*np.einsum(
            "si,qsi->qs", midpoint_tangent, tangent_sum_jacobian)[:, :, None]
    )/np.linalg.norm(tangent_sum, axis=1)[None, :, None]
    normals = np.cross(axis[None, :], midpoint_tangent)
    normal_jacobian = np.cross(
        axis[None, None, :], midpoint_tangent_jacobian)
    observation_points = np.vstack((
        centers, centers+offset*normals, centers-offset*normals))
    observation_jacobian = np.concatenate((
        center_jacobian,
        center_jacobian+offset*normal_jacobian,
        center_jacobian-offset*normal_jacobian), axis=1)
    observation_field = (
        iron_scale*np.asarray(
            iron_evaluator.field(observation_points, "direct"), dtype=float)
        + constant_field[None, :])
    observation_gradient = iron_scale*np.asarray(
        iron_evaluator.field_gradient(observation_points), dtype=float)
    observation_fixed_shape = (4.0*np.pi*iron_scale)*np.asarray(
        charge_gram.configured_field_values_shape_derivative(
            np.ascontiguousarray(observation_points),
            np.ascontiguousarray(state),
            np.ascontiguousarray(state_jacobian), cells, faces), dtype=float)
    observation_total_shape = (
        observation_fixed_shape
        + np.einsum("pij,qpj->qpi", observation_gradient,
                    observation_jacobian))
    field_response = np.r_[
        observation_field[:segment_count, 2],
        (observation_field[segment_count:2*segment_count, 2]
         - observation_field[2*segment_count:, 2])/(2.0*offset)]
    field_response_jacobian = np.concatenate((
        observation_total_shape[:, :segment_count, 2],
        (observation_total_shape[:, segment_count:2*segment_count, 2]
         - observation_total_shape[:, 2*segment_count:, 2])/(2.0*offset),
    ), axis=1).T
    transfer = combined_function_transfer_map_from_field_response(
        field_response, segment_lengths, rigidity,
        field_response_jacobian=field_response_jacobian,
        curvature_sign=curvature_sign,
        segment_length_jacobian=segment_length_jacobian,
        response_entries=response_entries)
    return ReclosedOrbitShapeJacobian(
        recovered, entrance_jacobian, length_jacobian,
        position_jacobian, tangent_jacobian,
        field_response, field_response_jacobian, transfer,
        closure_jacobian, observation_points,
        observation_fixed_shape, observation_jacobian)


def recover_periodic_planar_closed_orbit(
        field, *, magnetic_rigidity, cell_angle_rad, initial_radius_m,
        initial_incidence_angle_rad=0.0, n_segments=128,
        gradient_offset=1.0e-3, max_path_length_m=None,
        curvature_sign=1.0, position_tolerance=1.0e-9,
        tangent_tolerance=1.0e-9, root_max_evaluations=80,
        response_entries=None) -> FullFieldClosedOrbit:
    """Recover one-cell periodic orbit and its local transfer map.

    Two geometric unknowns are solved: the entrance radius and incidence
    angle on the radial cell boundary.  A trajectory is integrated through
    the supplied total magnetic field until it reaches the next radial cell
    boundary.  The exit point and tangent, rotated back by one cell angle,
    must equal the entrance data.  Numerical differencing used internally by
    the two-variable root finder concerns orbit recovery only; no design
    derivative or topology sensitivity is approximated.

    ``curvature_sign=+1`` matches :class:`PlanarTransferMatrixObjective`: a
    positive binormal field produces positive signed curvature.  It therefore
    fixes the otherwise conventional charge/field orientation in the Lorentz
    equation used here.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import least_squares

    rigidity = float(magnetic_rigidity)
    theta = float(cell_angle_rad)
    radius_guess = float(initial_radius_m)
    alpha_guess = float(initial_incidence_angle_rad)
    segments = int(n_segments)
    gradient_offset = float(gradient_offset)
    curvature_sign = float(curvature_sign)
    position_tolerance = float(position_tolerance)
    tangent_tolerance = float(tangent_tolerance)
    max_evaluations = int(root_max_evaluations)
    if (not np.all(np.isfinite([
            rigidity, theta, radius_guess, alpha_guess, gradient_offset,
            curvature_sign, position_tolerance, tangent_tolerance]))
            or rigidity <= 0.0 or radius_guess <= 0.0
            or theta <= 0.0 or theta >= np.pi
            or segments < 8 or gradient_offset <= 0.0
            or curvature_sign == 0.0 or position_tolerance <= 0.0
            or tangent_tolerance <= 0.0 or max_evaluations < 4):
        raise ValueError("invalid periodic closed-orbit recovery settings")
    if max_path_length_m is None:
        max_path = max(
            4.0 * radius_guess * theta, 0.25 * radius_guess)
    else:
        max_path = float(max_path_length_m)
    if not np.isfinite(max_path) or max_path <= 0.0:
        raise ValueError("max_path_length_m must be positive")

    axis = np.array([0.0, 0.0, 1.0])
    radial_0 = np.array([0.0, -1.0, 0.0])
    tangent_0 = np.array([1.0, 0.0, 0.0])
    cosine = float(np.cos(theta))
    sine = float(np.sin(theta))
    rotation = np.array([
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ])
    radial_1 = rotation @ radial_0
    tangent_1 = rotation @ tangent_0
    radius_scale = max(radius_guess, 1.0e-6)

    def integrate(radius, alpha, *, dense_output=False):
        start_position = radius * radial_0
        start_tangent = (
            np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
        state_0 = np.r_[start_position, start_tangent]

        def ode(_path_length, state):
            tangent = state[3:]
            tangent_norm = float(np.linalg.norm(tangent))
            if tangent_norm <= 0.0:
                raise RuntimeError("particle tangent vanished")
            tangent = tangent / tangent_norm
            magnetic_field = _evaluate_b_field(field, state[:3])
            curvature = (-curvature_sign
                         * np.cross(tangent, magnetic_field) / rigidity)
            return np.r_[tangent, curvature]

        def next_boundary(_path_length, state):
            return float(state[:3] @ tangent_1)

        next_boundary.direction = 1.0
        next_boundary.terminal = True
        return solve_ivp(
            ode, (0.0, max_path), state_0, method="DOP853",
            events=next_boundary, dense_output=dense_output,
            rtol=min(1.0e-10, 0.1 * tangent_tolerance),
            atol=min(1.0e-12, 0.1 * position_tolerance),
            max_step=max_path / max(segments, 32))

    def terminal_state(radius, alpha, *, dense_output=False):
        solution = integrate(radius, alpha, dense_output=dense_output)
        if (not solution.success or len(solution.t_events) != 1
                or len(solution.t_events[0]) != 1):
            return None, solution
        terminal = np.asarray(solution.y_events[0][0], dtype=float)
        if terminal[:3] @ radial_1 <= 0.0:
            return None, solution
        terminal[3:] /= np.linalg.norm(terminal[3:])
        return terminal, solution

    def residual(parameters):
        radius, alpha = parameters
        terminal, _ = terminal_state(radius, alpha)
        if terminal is None:
            return np.array([1.0e3, 1.0e3])
        position_back = rotation.T @ terminal[:3]
        tangent_back = rotation.T @ terminal[3:]
        entrance_tangent = (
            np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
        tangent_cross = float(axis @ np.cross(
            entrance_tangent, tangent_back))
        tangent_dot = float(entrance_tangent @ tangent_back)
        return np.array([
            (position_back @ radial_0 - radius) / radius_scale,
            np.arctan2(tangent_cross, tangent_dot),
        ])

    root = least_squares(
        residual, np.array([radius_guess, alpha_guess]),
        bounds=(np.array([0.2 * radius_guess, -0.45 * np.pi]),
                np.array([5.0 * radius_guess, 0.45 * np.pi])),
        xtol=min(1.0e-12, position_tolerance / radius_scale),
        ftol=min(1.0e-12, tangent_tolerance),
        gtol=min(1.0e-12, tangent_tolerance),
        max_nfev=max_evaluations)
    radius, alpha = (float(value) for value in root.x)
    terminal, solution = terminal_state(radius, alpha, dense_output=True)
    if terminal is None:
        raise RuntimeError(
            "periodic orbit did not reach the next radial boundary")
    final_residual = residual((radius, alpha))
    position_residual = abs(float(final_residual[0])) * radius_scale
    tangent_residual = abs(float(final_residual[1]))
    position_back = rotation.T @ terminal[:3]
    tangent_back = rotation.T @ terminal[3:]
    if (not root.success or position_residual > position_tolerance
            or tangent_residual > tangent_tolerance):
        raise RuntimeError(
            "periodic orbit closure failed: position residual "
            f"{position_residual:.6e} m, tangent residual "
            f"{tangent_residual:.6e}")
    path_length = float(solution.t_events[0][0])
    path_stations = np.linspace(0.0, path_length, segments + 1)
    states = np.asarray(solution.sol(path_stations), dtype=float).T
    tangents = states[:, 3:]
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    orbit = PlanarDesignOrbit(
        states[:, :3], tangents, magnetic_rigidity=rigidity,
        bend_axis=axis, path_length_stations=path_stations)
    response = sample_planar_orbit_field_response(
        field, orbit, gradient_offset=gradient_offset)
    transfer = combined_function_transfer_map_from_field_response(
        response, orbit.segment_lengths, rigidity,
        curvature_sign=curvature_sign,
        response_entries=response_entries)
    return FullFieldClosedOrbit(
        rigidity, orbit, path_length, radius, alpha,
        position_residual, tangent_residual,
        float(position_back[2]), float(tangent_back[2]), int(root.nfev),
        response, transfer)


def recover_periodic_planar_closed_orbit_native(
        field, *, iron_evaluator=None, iron_scale=0.0,
        iron_algorithm="auto", radia_object=-1,
        mirror_z=False, constant_field_t=(0.0, 0.0, 0.0),
        magnetic_rigidity, cell_angle_rad, initial_radius_m,
        initial_incidence_angle_rad=0.0, n_segments=128,
        gradient_offset=1.0e-3, tracking_step_m=5.0e-4,
        max_path_length_m=None, planarity_tolerance_m=1.0e-7,
        curvature_sign=1.0, position_tolerance=1.0e-9,
        tangent_tolerance=1.0e-9, root_max_evaluations=80,
        response_entries=None) -> FullFieldClosedOrbit:
    """Recover a periodic orbit with the fixed-step C++ field tracker.

    Each complete trajectory integration and exit-plane interpolation is
    performed in C++.  Python owns only the two-variable radius/incidence
    Broyden iteration: it forms one physically scaled numerical root Jacobian
    at the initial seed, then updates that 2-by-2 Jacobian by rank one.
    ``field`` is the matching total-field provider used only for the final
    orbit/gradient samples.  The native field is the sum of
    ``iron_evaluator*iron_scale``, the optional Radia object, and
    ``constant_field_t``.

    Numerical differencing inside the orbit root finder is not an optimizer
    design derivative.  Shape/topology sensitivities must continue to use the
    analytic HDiv-MMM contractions; a candidate geometry is accepted only
    after this routine re-integrates and re-closes its realized orbit.
    """
    from . import _radia_pybind as _native

    rigidity = float(magnetic_rigidity)
    theta = float(cell_angle_rad)
    radius_guess = float(initial_radius_m)
    alpha_guess = float(initial_incidence_angle_rad)
    segments = int(n_segments)
    gradient_offset = float(gradient_offset)
    step = float(tracking_step_m)
    planarity_tolerance = float(planarity_tolerance_m)
    curvature_sign = float(curvature_sign)
    position_tolerance = float(position_tolerance)
    tangent_tolerance = float(tangent_tolerance)
    max_evaluations = int(root_max_evaluations)
    iron_scale = float(iron_scale)
    iron_algorithm = str(iron_algorithm)
    radia_object = int(radia_object)
    constant_field = np.ascontiguousarray(
        constant_field_t, dtype=float).reshape(-1)
    if (iron_algorithm not in ("auto", "direct", "tree")
            or constant_field.shape != (3,)
            or not np.all(np.isfinite(np.r_[
                constant_field, rigidity, theta, radius_guess, alpha_guess,
                gradient_offset, step, planarity_tolerance, curvature_sign,
                position_tolerance, tangent_tolerance, iron_scale]))
            or rigidity <= 0.0 or radius_guess <= 0.0
            or theta <= 0.0 or theta >= np.pi or segments < 8
            or gradient_offset <= 0.0 or step <= 0.0
            or planarity_tolerance <= 0.0 or curvature_sign == 0.0
            or position_tolerance <= 0.0 or tangent_tolerance <= 0.0
            or max_evaluations < 4):
        raise ValueError("invalid native periodic closed-orbit settings")
    if max_path_length_m is None:
        max_path = max(
            4.0 * radius_guess * theta, 0.25 * radius_guess)
    else:
        max_path = float(max_path_length_m)
    if (not np.isfinite(max_path) or max_path <= step):
        raise ValueError(
            "max_path_length_m must be finite and exceed tracking_step_m")

    axis = np.array([0.0, 0.0, 1.0])
    radial_0 = np.array([0.0, -1.0, 0.0])
    tangent_0 = np.array([1.0, 0.0, 0.0])
    cosine = float(np.cos(theta))
    sine = float(np.sin(theta))
    rotation = np.array([
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ])
    radial_1 = rotation @ radial_0
    tangent_1 = np.ascontiguousarray(rotation @ tangent_0)
    radius_scale = max(radius_guess, 1.0e-6)
    # C++ uses dt/ds=t x B/(B rho), whereas the public planar convention is
    # positive B_binormal -> positive signed curvature.
    native_rigidity = -rigidity / curvature_sign

    def track(radius, alpha, station_count):
        start_position = np.ascontiguousarray(radius * radial_0)
        start_tangent = np.ascontiguousarray(
            np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
        return _native.track_reference_orbit_to_plane_native(
            iron_evaluator, iron_scale, radia_object, bool(mirror_z),
            constant_field, native_rigidity, start_position, start_tangent,
            tangent_1, 0.0, step, max_path, planarity_tolerance,
            int(station_count), iron_algorithm)

    track_evaluations = 0

    def terminal_state(radius, alpha):
        nonlocal track_evaluations
        track_evaluations += 1
        try:
            tracked = track(radius, alpha, 2)
        except (RuntimeError, ValueError):
            return None
        position = np.asarray(tracked[0][-1], dtype=float)
        tangent = np.asarray(tracked[1][-1], dtype=float)
        if position @ radial_1 <= 0.0:
            return None
        return position, tangent

    cached_parameters = None
    cached_residual = None

    def residual(parameters):
        nonlocal cached_parameters, cached_residual
        parameters = np.asarray(parameters, dtype=float)
        if (cached_parameters is not None
                and np.array_equal(parameters, cached_parameters)):
            return cached_residual.copy()
        radius, alpha = parameters
        terminal = terminal_state(radius, alpha)
        if terminal is None:
            value = np.array([1.0e3, 1.0e3])
        else:
            position_back = rotation.T @ terminal[0]
            tangent_back = rotation.T @ terminal[1]
            entrance_tangent = (
                np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
            tangent_cross = float(axis @ np.cross(
                entrance_tangent, tangent_back))
            tangent_dot = float(entrance_tangent @ tangent_back)
            value = np.array([
                (position_back @ radial_0 - radius) / radius_scale,
                np.arctan2(tangent_cross, tangent_dot),
            ])
        cached_parameters = parameters.copy()
        cached_residual = value.copy()
        return value

    def root_jacobian(parameters):
        # Use physical absolute perturbations larger than the fixed-step RK4
        # interpolation floor.  A generic sqrt(eps) perturbation can
        # differentiate that floor and trigger many needless integrations.
        # This Jacobian belongs only to the two-variable orbit-recovery root;
        # it is never exposed as a magnet design derivative.
        nonlocal cached_parameters, cached_residual
        parameters = np.asarray(parameters, dtype=float)
        base = residual(parameters)
        steps = np.array([
            max(1.0e-5 * radius_scale, 100.0 * position_tolerance),
            max(1.0e-5, 100.0 * tangent_tolerance),
        ])
        jacobian = np.empty((2, 2))
        for column, difference in enumerate(steps):
            shifted = parameters.copy()
            shifted[column] += difference
            jacobian[:, column] = (
                residual(shifted) - base) / difference
        cached_parameters = parameters.copy()
        cached_residual = base.copy()
        return jacobian

    lower = np.array([0.2 * radius_guess, -0.45 * np.pi])
    upper = np.array([5.0 * radius_guess, 0.45 * np.pi])
    parameters = np.array([radius_guess, alpha_guess])
    root_residual = residual(parameters)
    jacobian = root_jacobian(parameters)
    position_root_tolerance = position_tolerance / radius_scale

    def root_merit(value):
        return float(max(
            abs(value[0]) / position_root_tolerance,
            abs(value[1]) / tangent_tolerance))

    root_success = root_merit(root_residual) <= 1.0
    while not root_success and track_evaluations < max_evaluations:
        try:
            update = np.linalg.solve(jacobian, -root_residual)
        except np.linalg.LinAlgError:
            update = np.linalg.lstsq(
                jacobian, -root_residual, rcond=1.0e-12)[0]
        # Keep a poor initial field seed inside the same local sector rather
        # than asking the plane-crossing tracker to globalize the problem.
        update[0] = np.clip(
            update[0], -0.25 * radius_guess, 0.25 * radius_guess)
        update[1] = np.clip(update[1], -0.15, 0.15)
        base_merit = root_merit(root_residual)
        accepted = False
        damping = 1.0
        while track_evaluations < max_evaluations:
            candidate = np.clip(
                parameters + damping * update, lower, upper)
            candidate_residual = residual(candidate)
            candidate_merit = root_merit(candidate_residual)
            if candidate_merit < base_merit:
                accepted = True
                break
            damping *= 0.5
            if damping < 1.0 / 128.0:
                break
        if not accepted:
            break
        step_parameters = candidate - parameters
        step_residual = candidate_residual - root_residual
        denominator = float(step_parameters @ step_parameters)
        if denominator > 1.0e-30:
            jacobian += np.outer(
                step_residual - jacobian @ step_parameters,
                step_parameters) / denominator
        parameters = candidate
        root_residual = candidate_residual
        root_success = root_merit(root_residual) <= 1.0
    radius, alpha = (float(value) for value in parameters)
    if not root_success:
        raise RuntimeError(
            "periodic orbit closure failed after "
            f"{track_evaluations} native tracks: position residual "
            f"{abs(root_residual[0]) * radius_scale:.6e} m, tangent "
            f"residual {abs(root_residual[1]):.6e}")
    try:
        tracked = track(radius, alpha, segments + 1)
    except (RuntimeError, ValueError) as error:
        raise RuntimeError(
            "periodic orbit did not reach the next radial boundary") from error
    positions = np.ascontiguousarray(tracked[0], dtype=float)
    tangents = np.ascontiguousarray(tracked[1], dtype=float)
    path_stations = np.ascontiguousarray(tracked[2], dtype=float)
    curvature = np.ascontiguousarray(tracked[3], dtype=float)
    path_length = float(tracked[4])
    final_residual = residual((radius, alpha))
    position_residual = abs(float(final_residual[0])) * radius_scale
    tangent_residual = abs(float(final_residual[1]))
    position_back = rotation.T @ positions[-1]
    tangent_back = rotation.T @ tangents[-1]
    if (position_residual > position_tolerance
            or tangent_residual > tangent_tolerance):
        raise RuntimeError(
            "periodic orbit closure failed: position residual "
            f"{position_residual:.6e} m, tangent residual "
            f"{tangent_residual:.6e}")
    orbit = PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=rigidity,
        bend_axis=axis, path_length_stations=path_stations,
        signed_curvature_per_m=curvature)
    response = sample_planar_orbit_field_response(
        field, orbit, gradient_offset=gradient_offset)
    transfer = combined_function_transfer_map_from_field_response(
        response, orbit.segment_lengths, rigidity,
        curvature_sign=curvature_sign,
        response_entries=response_entries)
    return FullFieldClosedOrbit(
        rigidity, orbit, path_length, radius, alpha,
        position_residual, tangent_residual,
        float(position_back[2]), float(tangent_back[2]), track_evaluations,
        response, transfer)


def recover_ffag_closed_orbit(
        field, spec: FFAGSoftEdgeCellSpec, kinetic_energy_mev, *,
        initial_reference=None, n_segments=128, gradient_offset=1.0e-3,
        **recovery_options) -> FullFieldClosedOrbit:
    """Recover one realized periodic FFAG orbit using a reduced target seed."""
    if not isinstance(spec, FFAGSoftEdgeCellSpec):
        raise TypeError("spec must be an FFAGSoftEdgeCellSpec")
    energy = float(kinetic_energy_mev)
    reference = (build_ffag_cell_reference(
        spec, energy, n_segments=max(64, int(n_segments)))
        if initial_reference is None else initial_reference)
    if isinstance(reference, FFAGCellReference):
        if (abs(reference.kinetic_energy_mev - energy)
                > 1.0e-12 * max(1.0, energy)):
            raise ValueError("initial_reference kinetic energy does not match")
        seed_orbit = reference.orbit
        rigidity = reference.magnetic_rigidity_tm
        radius = None
        alpha = None
    elif isinstance(reference, FullFieldClosedOrbit):
        rigidity = magnetic_rigidity_from_kinetic_energy(energy)
        if (abs(reference.magnetic_rigidity_tm - rigidity)
                > 1.0e-10 * max(1.0, rigidity)):
            raise ValueError("initial full-field orbit rigidity does not match")
        seed_orbit = reference.orbit
        radius = reference.entrance_radius_m
        alpha = reference.entrance_incidence_angle_rad
    else:
        raise TypeError(
            "initial_reference must be an FFAGCellReference or "
            "FullFieldClosedOrbit")
    if radius is None:
        radial_0 = np.array([0.0, -1.0, 0.0])
        tangent_0 = np.array([1.0, 0.0, 0.0])
        radius = float(seed_orbit.positions[0] @ radial_0)
        entrance_tangent = seed_orbit.tangents[0]
        alpha = float(np.arctan2(
            entrance_tangent @ radial_0, entrance_tangent @ tangent_0))
    recovery_options.setdefault(
        "max_path_length_m", 2.5 * spec.cell_length_m)
    return recover_periodic_planar_closed_orbit(
        field, magnetic_rigidity=rigidity,
        cell_angle_rad=spec.cell_bend_angle_rad,
        initial_radius_m=radius,
        initial_incidence_angle_rad=alpha,
        n_segments=n_segments, gradient_offset=gradient_offset,
        **recovery_options)


def recover_ffag_closed_orbit_native(
        field, spec: FFAGSoftEdgeCellSpec, kinetic_energy_mev, *,
        iron_evaluator=None, iron_scale=0.0, iron_algorithm="auto",
        radia_object=-1,
        mirror_z=False, constant_field_t=(0.0, 0.0, 0.0),
        initial_reference=None, n_segments=128, gradient_offset=1.0e-3,
        **recovery_options) -> FullFieldClosedOrbit:
    """Native-tracker counterpart of :func:`recover_ffag_closed_orbit`."""
    if not isinstance(spec, FFAGSoftEdgeCellSpec):
        raise TypeError("spec must be an FFAGSoftEdgeCellSpec")
    energy = float(kinetic_energy_mev)
    reference = (build_ffag_cell_reference(
        spec, energy, n_segments=max(64, int(n_segments)))
        if initial_reference is None else initial_reference)
    if isinstance(reference, FFAGCellReference):
        if (abs(reference.kinetic_energy_mev - energy)
                > 1.0e-12 * max(1.0, energy)):
            raise ValueError("initial_reference kinetic energy does not match")
        seed_orbit = reference.orbit
        rigidity = reference.magnetic_rigidity_tm
        radius = None
        alpha = None
    elif isinstance(reference, FullFieldClosedOrbit):
        rigidity = magnetic_rigidity_from_kinetic_energy(energy)
        if (abs(reference.magnetic_rigidity_tm - rigidity)
                > 1.0e-10 * max(1.0, rigidity)):
            raise ValueError("initial full-field orbit rigidity does not match")
        seed_orbit = reference.orbit
        radius = reference.entrance_radius_m
        alpha = reference.entrance_incidence_angle_rad
    else:
        raise TypeError(
            "initial_reference must be an FFAGCellReference or "
            "FullFieldClosedOrbit")
    if radius is None:
        radial_0 = np.array([0.0, -1.0, 0.0])
        tangent_0 = np.array([1.0, 0.0, 0.0])
        radius = float(seed_orbit.positions[0] @ radial_0)
        entrance_tangent = seed_orbit.tangents[0]
        alpha = float(np.arctan2(
            entrance_tangent @ radial_0, entrance_tangent @ tangent_0))
    recovery_options.setdefault(
        "max_path_length_m", 2.5 * spec.cell_length_m)
    return recover_periodic_planar_closed_orbit_native(
        field, iron_evaluator=iron_evaluator, iron_scale=iron_scale,
        iron_algorithm=iron_algorithm,
        radia_object=radia_object, mirror_z=mirror_z,
        constant_field_t=constant_field_t, magnetic_rigidity=rigidity,
        cell_angle_rad=spec.cell_bend_angle_rad,
        initial_radius_m=radius,
        initial_incidence_angle_rad=alpha,
        n_segments=n_segments, gradient_offset=gradient_offset,
        **recovery_options)


def recover_ffag_closed_orbit_family(
        field, target_family: FFAGCellTargetFamily, *, n_segments=128,
        gradient_offset=1.0e-3, initial_references=None,
        **recovery_options
        ) -> tuple[FullFieldClosedOrbit, ...]:
    """Recover all momentum orbits from one shared realized magnet field."""
    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    references = (target_family.references if initial_references is None
                  else tuple(initial_references))
    if len(references) != len(target_family.references):
        raise ValueError(
            "initial_references must have one seed per target momentum")
    return tuple(recover_ffag_closed_orbit(
        field, target_family.spec, target.kinetic_energy_mev,
        initial_reference=seed, n_segments=n_segments,
        gradient_offset=gradient_offset, **recovery_options)
        for target, seed in zip(target_family.references, references))


def recover_ffag_closed_orbit_family_native(
        field, target_family: FFAGCellTargetFamily, *,
        iron_evaluator=None, iron_scale=0.0, iron_algorithm="auto",
        radia_object=-1,
        mirror_z=False, constant_field_t=(0.0, 0.0, 0.0),
        n_segments=128, gradient_offset=1.0e-3,
        initial_references=None, **recovery_options
        ) -> tuple[FullFieldClosedOrbit, ...]:
    """Recover a momentum family using the C++ trajectory integrator."""
    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    references = (target_family.references if initial_references is None
                  else tuple(initial_references))
    if len(references) != len(target_family.references):
        raise ValueError(
            "initial_references must have one seed per target momentum")
    return tuple(recover_ffag_closed_orbit_native(
        field, target_family.spec, target.kinetic_energy_mev,
        iron_evaluator=iron_evaluator, iron_scale=iron_scale,
        iron_algorithm=iron_algorithm,
        radia_object=radia_object, mirror_z=mirror_z,
        constant_field_t=constant_field_t, initial_reference=seed,
        n_segments=n_segments, gradient_offset=gradient_offset,
        **recovery_options)
        for target, seed in zip(target_family.references, references))


def _realized_ffag_band_ratios(
        recovered_orbits, objective: MultiMomentumTransferMatrixObjective
        ) -> tuple[np.ndarray, np.ndarray]:
    field_ratios = []
    matrix_ratios = []
    for recovered, target in zip(recovered_orbits, objective.objectives):
        count = len(recovered.orbit.segment_lengths)
        if count != len(target.bend_field_band):
            raise ValueError(
                "recovered and target orbit segment counts must match")
        required_field = (
            recovered.orbit.magnetic_rigidity
            * recovered.orbit.signed_curvature / target.curvature_sign)
        field_ratios.append(float(np.max(np.abs(
            (recovered.field_response[:count] - required_field)
            / target.bend_field_band))))
        realized_values = np.asarray([
            recovered.transfer.matrix[row, column]
            for row, column in target.response_entries])
        target_values = np.asarray([
            target.target_matrix[row, column]
            for row, column in target.response_entries])
        bands = np.asarray([
            target.transfer_matrix_band[row, column]
            for row, column in target.response_entries])
        matrix_ratios.append(float(np.max(np.abs(
            (realized_values - target_values) / bands))))
    return np.asarray(field_ratios), np.asarray(matrix_ratios)


@dataclass(frozen=True)
class FFAGHDivMMMOuterIteration:
    """One Lego update followed by exact full-field orbit reconstruction."""

    index: int
    material_move_fraction: float
    active_count_before: int
    active_count_after: int
    source_scale_before: float
    source_scale_after: float
    max_band_ratio_before: float
    max_band_ratio_after: float
    max_position_closure_residual_m: float
    max_tangent_closure_residual: float
    accepted: bool
    reason: str
    topology_result: MultiMomentumAcceleratorMagnetTopologyResult


@dataclass(frozen=True)
class FFAGHDivMMMTopologyResult:
    """Full-field closed-orbit outer loop around binary HDiv-MMM growth."""

    target_family: FFAGCellTargetFamily
    source_scale: float
    active_elements: np.ndarray
    state: np.ndarray
    recovered_orbits: tuple[FullFieldClosedOrbit, ...]
    orbit_field_max_band_ratios: np.ndarray
    transfer_matrix_max_band_ratios: np.ndarray
    history: tuple[FFAGHDivMMMOuterIteration, ...]
    converged: bool
    stop_reason: str
    topology: object

    @property
    def realized_transfer_matrices(self) -> np.ndarray:
        return np.asarray([
            recovered.transfer.matrix for recovered in self.recovered_orbits])

    @property
    def max_band_ratio(self) -> float:
        return float(max(
            np.max(self.orbit_field_max_band_ratios),
            np.max(self.transfer_matrix_max_band_ratios)))


@dataclass(frozen=True)
class FFAGFixedOrbitMapTrial:
    """One exactly resolved material proposal checked against the one-pass map."""

    optics_iteration: int
    trial_index: int
    material_move_fraction: float | None
    active_count_before: int
    active_count_after: int
    max_band_ratio_before: float
    max_band_ratio_after: float
    accepted: bool
    reason: str
    proposal_model: str = "field-target"
    exact_search_trace: tuple = ()


@dataclass(frozen=True)
class FFAGFixedOrbitHDivMMMTopologyResult:
    """One-pass FFAG result about caller-supplied design orbits.

    Unlike :class:`FFAGHDivMMMTopologyResult`, this contract performs no
    periodic closed-orbit search.  The momentum-indexed design orbits in the
    target family are frozen observation paths from entrance to exit.  Their
    transfer matrices are therefore the one-pass maps that the material
    inverse must reproduce.
    """

    target_family: FFAGCellTargetFamily | FFAGFixedDesignOrbitTargetFamily
    source_scale: float
    topology_result: MultiMomentumAcceleratorMagnetTopologyResult
    optics_history: tuple[MultiMomentumAcceleratorMagnetTopologyResult, ...]
    termination_reason: str
    initial_max_band_ratio: float
    map_trust_history: tuple[FFAGFixedOrbitMapTrial, ...]

    @property
    def active_elements(self) -> np.ndarray:
        return self.topology_result.active_elements

    @property
    def state(self) -> np.ndarray:
        return self.topology_result.generation.state

    @property
    def realized_transfer_matrices(self) -> np.ndarray:
        return self.topology_result.realized_transfer_matrices

    @property
    def orbit_field_max_band_ratios(self) -> np.ndarray:
        return self.topology_result.orbit_field_max_band_ratios

    @property
    def transfer_matrix_max_band_ratios(self) -> np.ndarray:
        return self.topology_result.transfer_matrix_max_band_ratios

    @property
    def topology(self):
        return self.topology_result.topology

    @property
    def field_correction(self):
        return self.topology_result.field_correction

    @property
    def history(self):
        return tuple(
            item
            for result in self.optics_history
            for item in result.generation.history)

    @property
    def converged(self) -> bool:
        return self.max_band_ratio <= 1.0

    @property
    def stop_reason(self) -> str:
        return self.termination_reason

    @property
    def max_band_ratio(self) -> float:
        return float(max(
            np.max(self.orbit_field_max_band_ratios),
            np.max(self.transfer_matrix_max_band_ratios)))


def optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        target_family: (FFAGCellTargetFamily
                        | FFAGFixedDesignOrbitTargetFamily), *, source,
        charge_gram, fes, inv_chi, active_elements, element_volumes,
        volume_max, gradient_offset=1.0e-3, source_scale=1.0,
        optimize_source_scale=True, initial_state=None,
        max_optics_iterations=1, material_iterations_per_optics=1,
        field_inverse_relative_tolerance=1.0e-3,
        field_inverse_basis=None,
        field_inverse_maximum_step_scale=1.0,
        field_inverse_line_search_steps=8,
        map_trust_region_trials=3, map_ratio_tolerance=1.0e-8,
        direct_map_oracle_fallback=False,
        direct_map_oracle_exact_beam_width=0,
        direct_map_oracle_exact_beam_depth=0,
        direct_map_oracle_graph_front_proposal_limit=0,
        **generation_options) -> FFAGFixedOrbitHDivMMMTopologyResult:
    """Optimize one magnet about fixed entrance-to-exit design orbits.

    This is the production proof-of-concept path when the design orbit and
    desired one-pass transfer matrix are inputs.  It deliberately omits ring
    closure and periodic-orbit recovery.  The optics TSVD first maps the
    transfer-matrix error to a sampled field target on those frozen paths;
    the independent Abe--Murata ACA--QR--TSVD material inverse then selects
    binary HDiv-MMM element additions/removals.  The field-to-map Jacobian is
    propagated by forward-mode AD with an exact Frechet matrix-exponential
    primitive.

    The caller owns ``ngsolve.TaskManager``.  No design finite difference,
    density interpolation, air volume mesh, or tracking root solve is used.
    A caller may reuse an ``initial_state`` already solved for the identical
    active set, source scale, RHS, and material law.  The configured operator
    and inactive-DOF constraints are reapplied and the true residual is gated
    before that state can skip the otherwise redundant initial solve.
    """
    from .topology_optimization import (
        solve_hdiv_mmm_active_elements,
        validate_hdiv_mmm_active_state,
    )

    if not isinstance(target_family, (
            FFAGCellTargetFamily, FFAGFixedDesignOrbitTargetFamily)):
        raise TypeError(
            "target_family must be an FFAGCellTargetFamily or "
            "FFAGFixedDesignOrbitTargetFamily")
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    active = np.asarray(active_elements, dtype=bool).reshape(-1).copy()
    volumes = np.asarray(element_volumes, dtype=float).reshape(-1)
    scale = float(source_scale)
    optics_count = int(max_optics_iterations)
    material_count = int(material_iterations_per_optics)
    map_trial_count = int(map_trust_region_trials)
    map_ratio_tolerance = float(map_ratio_tolerance)
    oracle_beam_width = int(direct_map_oracle_exact_beam_width)
    oracle_beam_depth = int(direct_map_oracle_exact_beam_depth)
    oracle_graph_limit = int(
        direct_map_oracle_graph_front_proposal_limit)
    if (active.shape != volumes.shape or not np.any(active)
            or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0)
            or not np.isfinite(scale) or scale <= 0.0):
        raise ValueError("invalid fixed-design-orbit topology settings")
    if optics_count < 1 or material_count < 1 or map_trial_count < 1:
        raise ValueError(
            "fixed-design-orbit iteration counts must be positive")
    if not np.isfinite(map_ratio_tolerance) or map_ratio_tolerance < 0.0:
        raise ValueError(
            "map_ratio_tolerance must be nonnegative and finite")
    if (oracle_beam_width < 0 or oracle_beam_depth < 0
            or ((oracle_beam_width == 0)
                != (oracle_beam_depth == 0))):
        raise ValueError(
            "direct map oracle beam width and depth must both be zero or "
            "both be positive")
    if oracle_graph_limit < 0:
        raise ValueError(
            "direct map oracle graph-front proposal limit must be "
            "nonnegative")
    if "max_iterations" in generation_options:
        raise TypeError(
            "use material_iterations_per_optics instead of max_iterations")
    if "exact_state_cache" in generation_options:
        raise TypeError("the fixed-design-orbit loop owns exact_state_cache")

    objective = target_family.objective
    source_rhs = source.assemble_hdiv_rhs(fes)
    rhs = scale * source_rhs
    solve_tolerance = float(generation_options.get(
        "solve_tolerance", 1.0e-9))
    if initial_state is None:
        state = solve_hdiv_mmm_active_elements(
            charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
            response_matrix=np.zeros((1, int(fes.ndof))),
            active_elements=active,
            solve_tolerance=solve_tolerance,
            solve_max_iterations=int(generation_options.get(
                "solve_max_iterations", 5000)),
            mass_riesz=bool(generation_options.get("mass_riesz", True)),
            cluster_coarse_size=int(generation_options.get(
                "cluster_coarse_size", 0)),
            cluster_deflation_size=int(generation_options.get(
                "cluster_deflation_size", 0)),
            recycle_size=int(generation_options.get("recycle_size", 0)))[0]
    else:
        state = np.asarray(initial_state, dtype=float).reshape(-1).copy()
        validate_hdiv_mmm_active_state(
            charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
            active_elements=active, state=state,
            solve_tolerance=solve_tolerance)

    if optimize_source_scale:
        total_field = CoilHDivTotalField(
            source, charge_gram, state, source_scale=scale)
        sampled = []
        required = []
        bands = []
        for target in objective.objectives:
            response = sample_planar_orbit_field_response(
                total_field, target.orbit,
                gradient_offset=gradient_offset)
            count = len(target.orbit.segment_lengths)
            sampled.extend(response[:count])
            required.extend(target.required_bend_field)
            bands.extend(target.bend_field_band)
        sampled = np.asarray(sampled, dtype=float)
        required = np.asarray(required, dtype=float)
        weights = 1.0 / np.asarray(bands, dtype=float)
        denominator = float((weights * sampled) @ (weights * sampled))
        calibration = float(
            (weights * sampled) @ (weights * required) / denominator)
        if (not np.isfinite(calibration) or calibration <= 0.0
                or denominator <= np.finfo(float).tiny):
            raise RuntimeError(
                "CoilBuilder source cannot be positively calibrated to the "
                "fixed design-orbit bend fields")
        scale *= calibration
        rhs = scale * source_rhs
        state *= calibration

    response_matrix = build_multi_orbit_field_response_matrix(
        charge_gram, objective, gradient_offset=gradient_offset)
    incident = scale * source.incident_orbit_field_response(
        objective, gradient_offset=gradient_offset)
    optics_history = []
    map_trust_history = []
    exact_state_cache = {}
    accepted_result = None
    initial_max_band_ratio = None
    termination_reason = "maximum fixed-orbit optics iterations reached"
    for optics_iteration in range(optics_count):
        current_raw_field = np.asarray(
            response_matrix @ state + incident, dtype=float)
        current_ratio = float(np.max(np.abs(
            (objective.transform(current_raw_field)
             - objective.response_target) / objective.response_band)))
        if initial_max_band_ratio is None:
            initial_max_band_ratio = current_ratio
        field_correction = solve_transfer_matrix_field_correction(
            objective, current_raw_field,
            field_basis=field_inverse_basis,
            relative_tolerance=field_inverse_relative_tolerance,
            maximum_step_scale=field_inverse_maximum_step_scale,
            line_search_steps=field_inverse_line_search_steps)
        requested_initial = generation_options.get(
            "initial_material_move_fraction")
        requested_maximum = generation_options.get(
            "maximum_material_move_fraction")
        trial_fraction = (None if requested_initial is None else
                          float(requested_initial))
        accepted = False
        last_attempt = None
        for trial_index in range(map_trial_count):
            trial_options = dict(generation_options)
            if trial_fraction is not None:
                trial_options["initial_material_move_fraction"] = trial_fraction
                trial_options["maximum_material_move_fraction"] = min(
                    trial_fraction,
                    trial_fraction if requested_maximum is None else
                    float(requested_maximum))
            last_attempt = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                objective.orbits, objective.target_matrices,
                transfer_matrix_band=objective.transfer_matrix_band,
                bend_field_band=objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                field_correction=field_correction,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                initial_state=state,
                exact_state_cache=exact_state_cache,
                response_entries=objective.response_entries,
                curvature_sign=objective.curvature_sign,
                gradient_sign=objective.gradient_sign,
                max_iterations=material_count,
                **trial_options)
            next_active = np.asarray(
                last_attempt.active_elements, dtype=bool).copy()
            candidate_ratio = float(max(
                np.max(last_attempt.orbit_field_max_band_ratios),
                np.max(last_attempt.transfer_matrix_max_band_ratios)))
            changed = not np.array_equal(next_active, active)
            accepted = bool(
                candidate_ratio <= 1.0 + map_ratio_tolerance or
                (changed and candidate_ratio
                 < current_ratio - map_ratio_tolerance))
            if accepted:
                reason = "accepted by exact fixed one-pass map gate"
            elif not changed:
                reason = "material inverse proposed no active-set change"
            else:
                reason = "rejected by exact fixed one-pass map gate"
            map_trust_history.append(FFAGFixedOrbitMapTrial(
                optics_iteration, trial_index, trial_fraction,
                int(np.count_nonzero(active)),
                int(np.count_nonzero(next_active)), current_ratio,
                candidate_ratio, accepted, reason, "field-target",
                tuple(last_attempt.generation.exact_search_trace)))
            if accepted:
                accepted_result = last_attempt
                optics_history.append(last_attempt)
                active = next_active
                state = last_attempt.generation.state.copy()
                break
            if trial_fraction is None:
                break
            trial_fraction *= 0.5
        if not accepted and direct_map_oracle_fallback:
            # The two-stage transfer->field->material inverse can stall when
            # its reachable field target is a poor local surrogate for the
            # original map norm.  Reuse the same forward-mode AD field-to-map
            # Jacobian directly in the all-candidate material contraction as
            # a bounded fallback.  This is still an exact chain rule, not a
            # design finite difference or a density relaxation.
            oracle_options = dict(generation_options)
            oracle_fraction = (None if requested_initial is None else
                               float(requested_initial))
            if oracle_fraction is not None:
                oracle_options["initial_material_move_fraction"] = (
                    oracle_fraction)
                oracle_options["maximum_material_move_fraction"] = min(
                    oracle_fraction,
                    oracle_fraction if requested_maximum is None else
                    float(requested_maximum))
            # The default remains one bounded global all-candidate proposal.
            # A caller may explicitly enable shallow nonmonotone look-ahead
            # for the direct oracle after the primary field-target lane stalls.
            oracle_options["exact_beam_width"] = oracle_beam_width
            oracle_options["exact_beam_depth"] = oracle_beam_depth
            oracle_options["graph_front_proposal_limit"] = (
                oracle_graph_limit)
            last_attempt = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                objective.orbits, objective.target_matrices,
                transfer_matrix_band=objective.transfer_matrix_band,
                bend_field_band=objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                initial_state=state,
                exact_state_cache=exact_state_cache,
                response_entries=objective.response_entries,
                curvature_sign=objective.curvature_sign,
                gradient_sign=objective.gradient_sign,
                max_iterations=material_count,
                **oracle_options)
            next_active = np.asarray(
                last_attempt.active_elements, dtype=bool).copy()
            candidate_ratio = float(max(
                np.max(last_attempt.orbit_field_max_band_ratios),
                np.max(last_attempt.transfer_matrix_max_band_ratios)))
            changed = not np.array_equal(next_active, active)
            accepted = bool(
                candidate_ratio <= 1.0 + map_ratio_tolerance or
                (changed and candidate_ratio
                 < current_ratio - map_ratio_tolerance))
            reason = (
                "accepted by direct analytic map-Jacobian oracle"
                if accepted else
                ("direct map oracle proposed no active-set change"
                 if not changed else
                 "rejected by exact fixed one-pass map gate"))
            map_trust_history.append(FFAGFixedOrbitMapTrial(
                optics_iteration, map_trial_count, oracle_fraction,
                int(np.count_nonzero(active)),
                int(np.count_nonzero(next_active)), current_ratio,
                candidate_ratio, accepted, reason,
                "direct-map-jacobian",
                tuple(last_attempt.generation.exact_search_trace)))
            if accepted:
                accepted_result = last_attempt
                optics_history.append(last_attempt)
                active = next_active
                state = last_attempt.generation.state.copy()
        if not accepted:
            termination_reason = "map-level trust-region proposals rejected"
            break
        if candidate_ratio <= 1.0 + map_ratio_tolerance:
            termination_reason = "fixed one-pass transfer bands reached"
            break

    if accepted_result is None:
        # Preserve the incumbent when every material proposal is rejected.
        # This fallback performs one exact active solve only on that failure
        # path; it never returns the last rejected topology as the design.
        baseline_options = dict(generation_options)
        accepted_result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
            objective.orbits, objective.target_matrices,
            transfer_matrix_band=objective.transfer_matrix_band,
            bend_field_band=objective.bend_field_band,
            charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
            rhs=rhs, field_response_matrix=response_matrix,
            incident_field_response=incident,
            field_correction=field_correction,
            active_elements=active, element_volumes=volumes,
            volume_max=volume_max,
            initial_state=state,
            exact_state_cache=exact_state_cache,
            response_entries=objective.response_entries,
            curvature_sign=objective.curvature_sign,
            gradient_sign=objective.gradient_sign,
            max_iterations=0, **baseline_options)
    return FFAGFixedOrbitHDivMMMTopologyResult(
        target_family, scale, accepted_result, tuple(optics_history),
        termination_reason, float(initial_max_band_ratio),
        tuple(map_trust_history))


def optimize_ffag_hdiv_mmm_from_design_orbits(
        design_orbits, target_transfer_matrices, *,
        transfer_matrix_band=1.0e-3, bend_field_band=1.0e-3,
        response_entries=None, controlled_components=None,
        require_symplectic=True, symplectic_tolerance=1.0e-9,
        curvature_sign=1.0, gradient_sign=1.0,
        **optimization_options) -> FFAGFixedOrbitHDivMMMTopologyResult:
    """Optimize HDiv-MMM material for caller-supplied one-pass optics.

    This is the direct public PoC entry point: the caller supplies one or more
    design orbits and the 6-by-6 transfer matrix required about each orbit.
    It builds no reduced reference field and performs no periodic-orbit search.
    ``controlled_components`` can restrict the objective to selected pole-face
    observables while the complete caller-supplied map remains subject to the
    symplectic input gate.
    The numerical path delegates to
    :func:`optimize_ffag_hdiv_mmm_from_fixed_design_orbits`, preserving its
    analytic field-to-map AD, ACA--QR--TSVD material inverse, exact active-set
    re-solves, and caller-owned ``ngsolve.TaskManager`` contract.
    """
    target_family = build_ffag_fixed_design_orbit_target_family(
        design_orbits, target_transfer_matrices,
        transfer_matrix_band=transfer_matrix_band,
        bend_field_band=bend_field_band,
        response_entries=response_entries,
        controlled_components=controlled_components,
        require_symplectic=require_symplectic,
        symplectic_tolerance=symplectic_tolerance,
        curvature_sign=curvature_sign,
        gradient_sign=gradient_sign)
    return optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        target_family, **optimization_options)


def optimize_ffag_hdiv_mmm_from_transfer_matrices(
        target_family: FFAGCellTargetFamily, *, source,
        charge_gram, fes, inv_chi, active_elements, element_volumes,
        volume_max, gradient_offset=1.0e-3, source_scale=1.0,
        optimize_source_scale=True,
        hdiv_order=1, orbit_segments=None, max_outer_iterations=8,
        inner_iterations=1, outer_initial_material_move_fraction=0.10,
        outer_trust_region_trials=3, outer_ratio_tolerance=1.0e-8,
        field_inverse_relative_tolerance=1.0e-3,
        field_inverse_basis=None,
        field_inverse_maximum_step_scale=1.0,
        field_inverse_line_search_steps=8,
        recovery_options=None, **generation_options
        ) -> FFAGHDivMMMTopologyResult:
    """Optimize a binary FFAG magnet with orbit recovery after every move.

    The fixed CoilBuilder source is assembled once into the HDiv RHS.  At each
    outer iteration the realized coil-plus-magnet field is tracked to recover
    every momentum-dependent periodic orbit.  A small dense optics TSVD first
    converts transfer-matrix error to a target correction of the sampled orbit
    field using the forward-mode AD field-to-map Jacobian.  Native HDiv rows
    are then rebuilt on those orbits and the separate Abe--Murata DUCAS
    ACA--QR--TSVD material inverse proposes exactly one (by default) Lego
    update.  HDiv candidates never enter the optics inverse.
    The topology is accepted only if a complete
    active-set solve, full-field orbit recovery, and transfer-map rebuild
    improve the actual band-normalized objective.  A rejected update shrinks
    the whole-element material trust region; it is never passed to Trafo.

    The caller owns ``ngsolve.TaskManager`` for the whole call.  No design
    finite difference, density interpolation, or air volume mesh is used.
    """
    from .topology_optimization import (
        ngsolve_growth_topology,
        solve_hdiv_mmm_active_elements,
    )

    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    active = np.asarray(active_elements, dtype=bool).reshape(-1).copy()
    volumes = np.asarray(element_volumes, dtype=float).reshape(-1)
    source_scale = float(source_scale)
    optimize_source_scale = bool(optimize_source_scale)
    outer_count = int(max_outer_iterations)
    inner_count = int(inner_iterations)
    trial_count = int(outer_trust_region_trials)
    move_fraction = float(outer_initial_material_move_fraction)
    ratio_tolerance = float(outer_ratio_tolerance)
    if (active.shape != volumes.shape or not np.any(active)
            or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0)
            or not np.isfinite(source_scale) or source_scale <= 0.0
            or outer_count < 1 or inner_count < 1 or trial_count < 1
            or not np.isfinite(move_fraction)
            or move_fraction <= 0.0 or move_fraction > 1.0
            or not np.isfinite(ratio_tolerance) or ratio_tolerance < 0.0):
        raise ValueError("invalid FFAG HDiv-MMM outer-loop settings")
    reserved = {
        "max_iterations", "initial_material_move_fraction",
        "maximum_material_move_fraction", "source_calibration_rows",
        "source_calibration_target"}
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "outer loop owns " + ", ".join(sorted(overlap)))
    segment_count = (len(target_family.references[0].orbit.segment_lengths)
                     if orbit_segments is None else int(orbit_segments))
    if (segment_count < 8 or any(
            len(reference.orbit.segment_lengths) != segment_count
            for reference in target_family.references)):
        raise ValueError(
            "orbit_segments must match every target-family orbit")
    recovery = {} if recovery_options is None else dict(recovery_options)
    recovery.setdefault("response_entries", target_family.objective.response_entries)

    source_rhs = source.assemble_hdiv_rhs(fes)
    rhs = source_scale * source_rhs
    zero_response = np.zeros((1, int(fes.ndof)))
    state = solve_hdiv_mmm_active_elements(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=zero_response, active_elements=active,
        solve_tolerance=float(generation_options.get(
            "solve_tolerance", 1.0e-9)),
        solve_max_iterations=int(generation_options.get(
            "solve_max_iterations", 5000)),
        mass_riesz=bool(generation_options.get("mass_riesz", True)),
        cluster_coarse_size=int(generation_options.get(
            "cluster_coarse_size", 0)),
        cluster_deflation_size=int(generation_options.get(
            "cluster_deflation_size", 0)),
        recycle_size=int(generation_options.get("recycle_size", 0)))[0]

    total_field = CoilHDivTotalField(
        source, charge_gram, state, source_scale=source_scale,
        hdiv_order=hdiv_order)
    if optimize_source_scale:
        sampled = []
        required = []
        bands = []
        for reference, target in zip(
                target_family.references,
                target_family.objective.objectives):
            response = sample_planar_orbit_field_response(
                total_field, reference.orbit,
                gradient_offset=gradient_offset)
            count = len(reference.orbit.segment_lengths)
            sampled.extend(response[:count])
            required.extend(target.required_bend_field)
            bands.extend(target.bend_field_band)
        sampled = np.asarray(sampled)
        required = np.asarray(required)
        weights = 1.0 / np.asarray(bands)
        denominator = float((weights * sampled) @ (weights * sampled))
        calibration = float(
            (weights * sampled) @ (weights * required) / denominator)
        if (not np.isfinite(calibration) or calibration <= 0.0
                or denominator <= np.finfo(float).tiny):
            raise RuntimeError(
                "CoilBuilder source cannot be positively calibrated to the "
                "target bend fields")
        source_scale *= calibration
        rhs = source_scale * source_rhs
        state *= calibration
        total_field = CoilHDivTotalField(
            source, charge_gram, state, source_scale=source_scale,
            hdiv_order=hdiv_order)
    recovered = recover_ffag_closed_orbit_family(
        total_field, target_family, n_segments=segment_count,
        gradient_offset=gradient_offset, **recovery)
    field_ratios, matrix_ratios = _realized_ffag_band_ratios(
        recovered, target_family.objective)
    current_ratio = float(max(np.max(field_ratios), np.max(matrix_ratios)))
    history = []
    stop_reason = "maximum outer iterations reached"

    for outer_index in range(outer_count):
        if current_ratio <= 1.0 + ratio_tolerance:
            stop_reason = "full-field orbit and transfer bands reached"
            break
        dynamic_objective = MultiMomentumTransferMatrixObjective(
            tuple(value.orbit for value in recovered),
            target_family.objective.target_matrices,
            target_family.objective.transfer_matrix_band,
            target_family.objective.bend_field_band,
            target_family.objective.response_entries,
            target_family.objective.curvature_sign,
            target_family.objective.gradient_sign)
        response_matrix = build_multi_orbit_field_response_matrix(
            charge_gram, dynamic_objective,
            gradient_offset=gradient_offset)
        incident = source_scale * source.incident_orbit_field_response(
            dynamic_objective, gradient_offset=gradient_offset)
        current_raw_field = np.asarray(
            response_matrix @ state + incident, dtype=float)
        field_correction = solve_transfer_matrix_field_correction(
            dynamic_objective, current_raw_field,
            field_basis=field_inverse_basis,
            relative_tolerance=field_inverse_relative_tolerance,
            maximum_step_scale=field_inverse_maximum_step_scale,
            line_search_steps=field_inverse_line_search_steps)
        calibration_rows = []
        calibration_target = []
        if optimize_source_scale:
            offsets = dynamic_objective.raw_offsets
            for index, objective in enumerate(dynamic_objective.objectives):
                count = len(objective.orbit.segment_lengths)
                calibration_rows.extend(range(
                    int(offsets[index]), int(offsets[index]) + count))
                calibration_target.extend(objective.required_bend_field)
            calibration_rows = np.asarray(calibration_rows, dtype=np.int64)
            calibration_target = np.asarray(calibration_target, dtype=float)
        accepted = False
        attempted_result = None
        trial_fraction = move_fraction
        trial_reason = "no improving full-field update"
        for _ in range(trial_count):
            attempted_result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                dynamic_objective.orbits,
                dynamic_objective.target_matrices,
                transfer_matrix_band=dynamic_objective.transfer_matrix_band,
                bend_field_band=dynamic_objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                field_correction=field_correction,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                response_entries=dynamic_objective.response_entries,
                curvature_sign=dynamic_objective.curvature_sign,
                gradient_sign=dynamic_objective.gradient_sign,
                max_iterations=inner_count,
                initial_material_move_fraction=trial_fraction,
                maximum_material_move_fraction=trial_fraction,
                source_calibration_rows=(
                    calibration_rows if optimize_source_scale else None),
                source_calibration_target=(
                    calibration_target if optimize_source_scale else None),
                **generation_options)
            candidate_active = attempted_result.active_elements
            if np.array_equal(candidate_active, active):
                trial_reason = attempted_result.generation.stop_reason
                break
            candidate_source_scale = (
                source_scale * attempted_result.generation.source_scale)
            candidate_field = CoilHDivTotalField(
                source, charge_gram, attempted_result.generation.state,
                source_scale=candidate_source_scale, hdiv_order=hdiv_order)
            try:
                candidate_recovered = recover_ffag_closed_orbit_family(
                    candidate_field, target_family, n_segments=segment_count,
                    gradient_offset=gradient_offset,
                    initial_references=recovered, **recovery)
            except RuntimeError as exc:
                trial_reason = "orbit recovery rejected: " + str(exc)
                trial_fraction *= 0.5
                continue
            candidate_field_ratios, candidate_matrix_ratios = (
                _realized_ffag_band_ratios(
                    candidate_recovered, target_family.objective))
            candidate_ratio = float(max(
                np.max(candidate_field_ratios),
                np.max(candidate_matrix_ratios)))
            if (candidate_ratio < current_ratio - ratio_tolerance
                    or candidate_ratio <= 1.0 + ratio_tolerance):
                accepted = True
                trial_reason = "accepted after full-field orbit rebuild"
                break
            trial_reason = (
                "full-field objective did not improve "
                f"({current_ratio:.6e} -> {candidate_ratio:.6e})")
            trial_fraction *= 0.5

        if attempted_result is None:
            raise RuntimeError("FFAG outer loop did not evaluate a proposal")
        if accepted:
            before_count = int(np.count_nonzero(active))
            source_scale_before = source_scale
            active = np.asarray(candidate_active, dtype=bool).copy()
            state = attempted_result.generation.state.copy()
            source_scale = candidate_source_scale
            rhs = source_scale * source_rhs
            recovered = candidate_recovered
            previous_ratio = current_ratio
            field_ratios = candidate_field_ratios
            matrix_ratios = candidate_matrix_ratios
            current_ratio = candidate_ratio
            move_fraction = min(1.0, 1.5 * trial_fraction)
            after_count = int(np.count_nonzero(active))
        else:
            before_count = int(np.count_nonzero(active))
            after_count = before_count
            source_scale_before = source_scale
            previous_ratio = current_ratio
        history.append(FFAGHDivMMMOuterIteration(
            outer_index, trial_fraction, before_count, after_count,
            source_scale_before, source_scale,
            previous_ratio, current_ratio,
            max(value.periodic_position_residual_m for value in recovered),
            max(value.periodic_tangent_residual for value in recovered),
            accepted, trial_reason, attempted_result))
        if not accepted:
            stop_reason = trial_reason
            break
    else:
        stop_reason = "maximum outer iterations reached"

    converged = bool(current_ratio <= 1.0 + ratio_tolerance)
    if converged:
        stop_reason = "full-field orbit and transfer bands reached"
    return FFAGHDivMMMTopologyResult(
        target_family, source_scale, active, state,
        tuple(recovered), field_ratios, matrix_ratios, tuple(history),
        converged, stop_reason,
        ngsolve_growth_topology(fes.mesh, active))


__all__ = [
    "EngeFringeIntegrals",
    "FFAGCyclicDensityMap",
    "FFAGCyclicSectorContract",
    "FFAGCellReference",
    "FFAGCellTargetFamily",
    "FFAGFixedDesignOrbitTargetFamily",
    "FFAGHDivMMMOuterIteration",
    "FFAGHDivMMMTopologyResult",
    "FFAGFixedOrbitHDivMMMTopologyResult",
    "FFAGFixedOrbitMapTrial",
    "FFAGSoftEdgeCellSpec",
    "FullFieldClosedOrbit",
    "ReclosedOrbitShapeJacobian",
    "GEV_C_PER_TESLA_METRE",
    "PROTON_REST_ENERGY_MEV",
    "build_ffag_cell_reference",
    "build_ffag_cell_target_family",
    "build_ffag_cyclic_density_map",
    "build_ffag_fixed_design_orbit_target_family",
    "enge_fringe_integrals",
    "differentiate_recovered_planar_orbit_shape_native",
    "identify_ffag_cyclic_sector_vertices",
    "magnetic_rigidity_from_kinetic_energy",
    "optimize_ffag_hdiv_mmm_from_fixed_design_orbits",
    "optimize_ffag_hdiv_mmm_from_design_orbits",
    "optimize_ffag_hdiv_mmm_from_transfer_matrices",
    "recover_ffag_closed_orbit",
    "recover_ffag_closed_orbit_family",
    "recover_ffag_closed_orbit_family_native",
    "recover_ffag_closed_orbit_native",
    "recover_periodic_planar_closed_orbit",
    "recover_periodic_planar_closed_orbit_native",
    "sample_planar_orbit_field_response",
    "validate_ffag_cyclic_sector_contract",
]
