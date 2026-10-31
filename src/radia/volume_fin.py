"""Bounded finite-volume PEEC circuit and exact Cartesian regional measurements.

Linear nonmagnetic MQS, complex peak phasors, no SIBC. This development route
does not certify curved fins, a return-current geometry, or physical accuracy.
The explicit partial-potential convention must match any independent reference.
"""
from dataclasses import dataclass
import hashlib
import math
import numpy as np


@dataclass(frozen=True)
class VolumeFinSolution:
    topology_sha256: str
    frequency_hz: float
    drive_current_peak_A: complex
    branch_current_peak_A: np.ndarray
    inductance_H: np.ndarray
    impedance_ohm: complex
    resistance_ohm: float
    partial_inductance_H: float
    kcl_relative_residual: float
    kvl_relative_residual: float
    native_port_impedance_ohm: complex


def _topology_hash(topology):
    digest = hashlib.sha256()
    for name in topology.__dataclass_fields__:
        value = getattr(topology, name)
        digest.update(name.encode('ascii'))
        if isinstance(value, np.ndarray):
            digest.update(str((value.dtype.str, value.shape)).encode('ascii'))
            digest.update(value.tobytes())
        else:
            digest.update(str(value).encode('ascii'))
    return digest.hexdigest()


def _validate_topology(topology):
    """Validate the public DTO before entering unchecked native array access."""
    from .volume_peec import VolumePEECTopology
    if not isinstance(topology, VolumePEECTopology):
        raise ValueError('VolumePEECTopology from build_volume_peec required')
    n, nodes, cells = len(topology.branches), len(topology.nodes_m), len(topology.cell_indices)
    shapes = dict(nodes_m=(nodes, 3), branches=(n, 2), support_boxes_m=(n, 2, 3),
                  axes=(n,), signs=(n,), areas_m2=(n,), resistance_ohm=(n,),
                  support_cells=(n, 2), half_resistance_ohm=(n, 2),
                  cell_indices=(cells, 3), cell_boxes_m=(cells, 2, 3), sigma_S_per_m=(cells,))
    if not n or not cells or nodes != cells+2:
        raise ValueError('Nonempty connected volume with two contact nodes required')
    for name, shape in shapes.items():
        value = np.asarray(getattr(topology, name))
        if value.shape != shape or not np.isfinite(value).all():
            raise ValueError(f'Invalid finite topology field {name}: expected {shape}')
    for name in ('branches', 'axes', 'signs', 'support_cells', 'cell_indices'):
        if not np.issubdtype(getattr(topology, name).dtype, np.integer):
            raise ValueError(f'Integer topology field {name} required')
    if (np.any(topology.branches < 0) or np.any(topology.branches >= nodes)
            or np.any(topology.branches[:, 0] == topology.branches[:, 1])):
        raise ValueError('Invalid branch node index')
    if (np.any(topology.support_cells < -1) or np.any(topology.support_cells >= cells)
            or np.any(topology.axes < 0) or np.any(topology.axes > 2)
            or np.any(abs(topology.signs) != 1)):
        raise ValueError('Invalid support/material index or Cartesian orientation')
    if any(np.any(getattr(topology, name) <= 0) for name in ('areas_m2','resistance_ohm','sigma_S_per_m')):
        raise ValueError('Positive area, resistance and conductivity required')
    if (np.any(topology.half_resistance_ohm < 0)
            or np.any(topology.support_boxes_m[:, 1] <= topology.support_boxes_m[:, 0])
            or np.any(topology.cell_boxes_m[:, 1] <= topology.cell_boxes_m[:, 0])):
        raise ValueError('Positive support volume and nonnegative half resistance required')
    for port in (topology.port_plus, topology.port_minus):
        if isinstance(port, (bool, np.bool_)) or not isinstance(port, (int, np.integer)) or not 0 <= port < nodes:
            raise ValueError('Invalid contact node index')
    if topology.port_plus == topology.port_minus:
        raise ValueError('Distinct physical contacts required')


def solve_volume_fin(topology, *, frequency_hz, drive_current_peak_A,
                     quadrature_order):
    """Explicit MKL dense-MNA reference route, with original-system checks.

    Not a scalable backend: reject signed-index overflow before native casts.
    No star/capacitive terms and no added internal/surface impedance.
    """
    frequency = float(frequency_hz)
    drive = complex(drive_current_peak_A)
    if not math.isfinite(frequency) or frequency < 0:
        raise ValueError('Finite nonnegative frequency required')
    if not np.isfinite(drive) or abs(drive) == 0:
        raise ValueError('Finite nonzero peak drive current required')
    _validate_topology(topology)
    n, nodes = len(topology.branches), len(topology.nodes_m)
    if max(n, nodes, n+nodes+1)**2 > np.iinfo(np.int32).max:
        raise ValueError('Dense MNA signed-index budget exceeded; no automatic fallback')
    from . import peec_matrices
    builder = topology.native_builder(quadrature_order=quadrature_order)
    L, _, _, _ = builder.build(include_star=False)
    L = np.asarray(L, dtype=float)
    if not np.isfinite(L).all():
        raise RuntimeError('Nonfinite finite-volume inductance operator')
    solver = peec_matrices.MNASolver(
        L, topology.resistance_ohm, np.asarray(topology.branches, dtype=np.int32),
        nodes, [(topology.port_plus, topology.port_minus, 0)])
    solver.set_solver_method(0)  # explicit native MKL direct route
    current = np.asarray(solver.solve_branch_currents(frequency, np.array([drive])))
    port_z = complex(np.asarray(solver.solve_z_matrix(frequency))[0, 0])
    if current.shape != (n,) or not np.isfinite(current).all() or not np.isfinite(port_z):
        raise RuntimeError('Nonfinite/invalid native MNA result')
    drop = topology.resistance_ohm*current + 2j*np.pi*frequency*(L@current)
    load = np.zeros(nodes, dtype=complex)
    load[topology.port_plus] = -drive
    load[topology.port_minus] = drive
    kcl = float(np.linalg.norm(topology.incidence()@current-load)/np.linalg.norm(load))
    # Tree recovery avoids a second linear solver or least-squares projection.
    adjacency = [[] for _ in range(nodes)]
    for k, (a, b) in enumerate(topology.branches):
        adjacency[a].append((b, k, 1))
        adjacency[b].append((a, k, -1))
    voltage = np.zeros(nodes, dtype=complex)
    seen = {topology.port_minus}; stack = [topology.port_minus]
    while stack:
        a = stack.pop()
        for b, k, sign in adjacency[a]:
            if b not in seen:
                voltage[b] = voltage[a]-sign*drop[k]
                seen.add(b); stack.append(b)
    if len(seen) != nodes:
        raise RuntimeError('Disconnected circuit, no implicit component removal')
    expected_drop = voltage[topology.branches[:, 0]]-voltage[topology.branches[:, 1]]
    kvl = float(np.linalg.norm(drop-expected_drop)/max(np.linalg.norm(drop), np.finfo(float).tiny))
    if not np.isfinite([kcl, kvl]).all() or max(kcl, kvl) > 1e-6:
        raise RuntimeError(f'Original MNA residual exceeds 1e-6: KCL={kcl}, KVL={kvl}')
    norm = abs(drive)**2
    resistance = float(np.dot(topology.resistance_ohm, abs(current)**2)/norm)
    inductance = float(np.vdot(current, L@current).real/norm)
    impedance = complex((voltage[topology.port_plus]-voltage[topology.port_minus])/drive)
    mismatch = abs(impedance-port_z)/max(abs(port_z), np.finfo(float).tiny)
    energy_z = resistance+2j*np.pi*frequency*inductance
    energy_error = abs(impedance-energy_z)/max(abs(impedance), np.finfo(float).tiny)
    if max(mismatch, energy_error) > 1e-6:
        raise RuntimeError(f'Port/branch inconsistency: native={mismatch}, energy={energy_error}')
    current = np.array(current, copy=True); current.setflags(write=False)
    L = np.array(L, copy=True); L.setflags(write=False)
    return VolumeFinSolution(_topology_hash(topology), frequency, drive, current, L, impedance, resistance,
                             inductance, kcl, kvl, port_z)


def _intersection(a, b):
    box = np.array([np.maximum(a[0], b[0]), np.minimum(a[1], b[1])])
    return box if np.all(box[1] > box[0]) else None


def _partition(box, regions):
    # Rectilinear support restriction is exact; never approximate rotated planes.
    axes = []
    for vector in (regions.fin_axis, regions.span_axis):
        v = np.asarray(vector)
        nonzero = np.flatnonzero(v)
        if len(nonzero) != 1 or abs(v[nonzero[0]]) != 1:
            raise ValueError('Volume ROI currently requires signed Cartesian axes')
        axes.append((int(nonzero[0]), float(v[nonzero[0]])))
    cuts = [[float(box[0, a]), float(box[1, a])] for a in range(3)]
    for (axis, sign), offsets in zip(axes, (
            (-regions.root_half_width_m, regions.root_half_width_m, regions.tip_start_m),
            regions.axial_bounds_m)):
        for offset in offsets:
            plane = regions.origin_m[axis]+sign*offset
            if box[0, axis] < plane < box[1, axis]:
                cuts[axis].append(plane)
    cuts = [sorted(set(c)) for c in cuts]
    import itertools
    for index in itertools.product(*(range(len(c)-1) for c in cuts)):
        part = np.array([[cuts[a][index[a]] for a in range(3)],
                         [cuts[a][index[a]+1] for a in range(3)]])
        yield regions.point_region(part.mean(axis=0)) or 'outside_window', part


def measure_volume_fin(topology, solution, regions, *, quadrature_order, convention):
    """Original restricted Green pairings, with body mutual allocated by halves.

    Subdivide CURRENT SUPPORTS at actual measurement planes, not centroids or
    fractions of complete L rows. A cut's current scales by its transverse area;
    this preserves J and the original volume integral. Heat uses the actual
    material half-cell intersection, including conductivity jumps.
    """
    from .fin_observables import allocate_fin_partial_inductance
    from . import peec_matrices
    _validate_topology(topology)
    # Validate quadrature through the same public production guard.
    topology.native_builder(quadrature_order=quadrature_order)
    if solution.topology_sha256 != _topology_hash(topology):
        raise ValueError('Solution/topology identity mismatch')
    if not isinstance(convention, str) or not convention:
        raise ValueError('Explicit partial-potential/return convention required')
    names = tuple(name for name, _, _ in regions.intervals())+('outside_window',)
    pieces = []
    heat = dict.fromkeys(names, 0.)
    for i, support in enumerate(topology.support_boxes_m):
        axis = int(topology.axes[i])
        for name, box in _partition(support, regions):
            size = box[1]-box[0]
            area = float(np.prod(np.delete(size, axis)))
            current = solution.branch_current_peak_A[i]*area/topology.areas_m2[i]
            pieces.append((name, box, axis, int(topology.signs[i]), current))
            for cell in topology.support_cells[i]:
                if cell < 0:
                    continue
                half = _intersection(box, topology.cell_boxes_m[cell])
                if half is not None:
                    volume = float(np.prod(half[1]-half[0]))
                    heat[name] += .5*abs(solution.branch_current_peak_A[i]/topology.areas_m2[i])**2*volume/topology.sigma_S_per_m[cell]
    if len(pieces)**2 > np.iinfo(np.int32).max:
        raise ValueError('Regional Green pair budget exceeded; no automatic fallback')
    builder = peec_matrices.PEECBuilder(volume_order=quadrature_order)
    for _, box, axis, sign, _ in pieces:
        size = box[1]-box[0]; direction = np.zeros(3); direction[axis] = sign
        builder.add_segment(box.mean(axis=0), direction, size[axis],
                            size[1 if axis == 2 else 2], size[1 if axis == 0 else 0], 1.)
    q = np.zeros((len(names), len(names)), dtype=complex)
    for i, (name_i, _, _, _, ci) in enumerate(pieces):
        for j, (name_j, _, _, _, cj) in enumerate(pieces):
            q[names.index(name_i), names.index(name_j)] += ci.conjugate()*builder.inductance_entry(i, j)*cj
    allocation = allocate_fin_partial_inductance(
        names, q, drive_current_peak_A=solution.drive_current_peak_A, convention=convention)
    norm = abs(solution.drive_current_peak_A)**2
    return dict(status='diagnostic-not-validated', region_definition=regions.as_dict(),
                mean_power_W=heat, equivalent_resistance_ohm={k: 2*v/norm for k, v in heat.items()},
                partial_inductance=allocation, restricted_support_count=len(pieces),
                green_partition_difference_H=allocation['total_H']-solution.partial_inductance_H,
                heat_partition_difference_W=sum(heat.values())-sum(topology.mean_branch_power_W(solution.branch_current_peak_A)),
                physics_tolerances=None)
