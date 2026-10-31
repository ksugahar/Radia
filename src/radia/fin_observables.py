"""Common finite-region measurements for straight-prism fin comparisons.

No field solver lives here. PEEC loss is piecewise uniform over each branch's
dual surface support; this is a measurement of that discretization, not proof
of skin-effect accuracy. Surface/volume references must use the same planes
with their own clipped-cell quadrature, never unmatched centroid selections.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math


def _point(value):
    result = tuple(float(x) for x in value)
    if len(result) != 3 or not all(map(math.isfinite, result)):
        raise ValueError("Expected a finite SI xyz point/vector")
    return result


def _sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def _add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def _scale(a, factor):
    return tuple(x*factor for x in a)


def _dot(a, b):
    return math.fsum(x*y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _area(poly):
    if len(poly) < 3:
        return 0.0
    vectors = [_cross(_sub(poly[i], poly[0]), _sub(poly[i+1], poly[0]))
               for i in range(1, len(poly)-1)]
    return 0.5*math.sqrt(math.fsum(math.fsum(v[k] for v in vectors)**2
                                  for k in range(3)))


def _clip(poly, origin, axis, limit, lower):
    """Clip by a plane; upper faces are excluded (half-open ownership)."""
    if not poly:
        return []
    def distance(p):
        d = _dot(_sub(p, origin), axis)-limit
        return d if lower else -d
    def inside(d):
        return d >= 0 if lower else d > 0
    output = []
    previous = poly[-1]
    dp = distance(previous)
    for current in poly:
        dc = distance(current)
        if inside(dp) != inside(dc):
            output.append(_add(previous, _scale(_sub(current, previous), dp/(dp-dc))))
        if inside(dc):
            output.append(current)
        previous, dp = current, dc
    return output


@dataclass(frozen=True)
class FinMeasurementRegions:
    """Disjoint body/root/shaft/tip intervals in an explicit orthonormal frame.

    Origin is on the root plane; all lengths are metres. The span interval is
    relative to origin and excludes its upper plane. These are geometry
    settings, not accuracy tolerances. No default physical region is inferred.
    """
    origin_m: tuple
    fin_axis: tuple
    span_axis: tuple
    root_half_width_m: float
    tip_start_m: float
    axial_bounds_m: tuple

    def __post_init__(self):
        for name in ("origin_m", "fin_axis", "span_axis"):
            object.__setattr__(self, name, _point(getattr(self, name)))
        bounds = tuple(float(x) for x in self.axial_bounds_m)
        if (len(bounds) != 2 or not all(map(math.isfinite, bounds))
                or bounds[0] >= bounds[1]):
            raise ValueError("An explicit increasing finite axial interval is required")
        object.__setattr__(self, "axial_bounds_m", bounds)
        if (not math.isfinite(self.root_half_width_m) or not math.isfinite(self.tip_start_m)
                or not 0 < self.root_half_width_m < self.tip_start_m):
            raise ValueError("Require 0 < root half-width < tip start in metres")
        # Numerical frame consistency only, never a physical acceptance gate.
        if (not math.isclose(_dot(self.fin_axis, self.fin_axis), 1., rel_tol=0., abs_tol=1e-12)
                or not math.isclose(_dot(self.span_axis, self.span_axis), 1., rel_tol=0., abs_tol=1e-12)
                or abs(_dot(self.fin_axis, self.span_axis)) > 1e-12):
            raise ValueError("Fin/span axes must be orthonormal; no implicit normalization")

    def intervals(self):
        w = self.root_half_width_m
        return (("body", None, -w), ("root", -w, w),
                ("fin_shaft", w, self.tip_start_m), ("tip", self.tip_start_m, None))

    def as_dict(self):
        return {"schema": "radia.fin.measurement-regions.v1", **asdict(self),
                "ownership": "half-open intervals; upper boundary belongs to next region"}

    def point_region(self, point):
        """For reference quadrature points, not centroid-only cell assignment."""
        d = _sub(_point(point), self.origin_m)
        z, u = _dot(d, self.span_axis), _dot(d, self.fin_axis)
        if not self.axial_bounds_m[0] <= z < self.axial_bounds_m[1]:
            return None
        return next(name for name, lo, hi in self.intervals()
                    if (lo is None or lo <= u) and (hi is None or u < hi))

    def polygon_weights(self, vertices):
        """Exact planar clipping of a convex constant-density support polygon.

        Weights are relative to the FULL support, not renormalized to the
        axial window. The signed partition difference is reported, not hidden.
        """
        poly = [_point(p) for p in vertices]
        area = _area(poly)
        if len(poly) < 3 or area <= 0:
            raise ValueError("Positive-area convex planar support required")
        normal = _cross(_sub(poly[1], poly[0]), _sub(poly[2], poly[0]))
        norm = math.sqrt(_dot(normal, normal))
        length = max(math.dist(p, poly[0]) for p in poly)
        if norm == 0 or any(abs(_dot(_sub(p, poly[0]), normal)) > 1e-12*length*norm for p in poly):
            raise ValueError("Support must be planar with non-collinear leading vertices")
        if any(_dot(_cross(_sub(poly[(i+1) % len(poly)], poly[i]),
                          _sub(poly[(i+2) % len(poly)], poly[(i+1) % len(poly)])), normal)
               < -1e-12*length**2*norm for i in range(len(poly))):
            raise ValueError("Support must be ordered and convex")
        active = _clip(poly, self.origin_m, self.span_axis, self.axial_bounds_m[0], True)
        active = _clip(active, self.origin_m, self.span_axis, self.axial_bounds_m[1], False)
        weights = {}
        for name, lo, hi in self.intervals():
            clipped = active
            if lo is not None:
                clipped = _clip(clipped, self.origin_m, self.fin_axis, lo, True)
            if hi is not None:
                clipped = _clip(clipped, self.origin_m, self.fin_axis, hi, False)
            weights[name] = _area(clipped)/area
        window = _area(active)/area
        return dict(area_m2=area, window_weight=window, regions=weights,
                    partition_weight_defect=math.fsum(weights.values())-window)


def peec_fin_loss_observables(graph, resistance_diagonal_ohm, current_peak_A, regions):
    """Measure ALL longitudinal and transverse branch losses on a straight prism.

    Each branch's R|I_peak|²/2 is distributed uniformly over its dual surface
    support. This keeps the existing PEEC loss law; no BEM/FEM matching factor
    or dropped transverse modes are allowed. Dense/off-diagonal R is rejected:
    its cross terms have no uniquely defined branch-local loss allocation.
    """
    if not isinstance(regions, FinMeasurementRegions):
        raise TypeError("Explicit FinMeasurementRegions required")
    r = [float(x) for x in resistance_diagonal_ohm]
    currents = [complex(x) for x in current_peak_A]
    count = len(graph.branch_xyz)
    if len(r) != count or len(currents) != count or len(graph.branch_kind) != count:
        raise ValueError("R/current/kind arrays must include every branch")
    if (any(not math.isfinite(x) or x < 0 for x in r)
            or any(not math.isfinite(x.real) or not math.isfinite(x.imag) for x in currents)):
        raise ValueError("Finite passive branch R and finite complex peak currents required")
    ns, nl = int(graph.n_stations), int(graph.n_lanes)
    stations = tuple(graph.mesh_stations)
    if (ns < 3 or nl < 3 or stations != tuple(sorted(set(stations)))
            or any(not isinstance(s, int) or not 0 < s < ns-1 for s in stations)):
        raise ValueError("At least three stations/lanes and interior mesh stations required")
    if count != (ns-1+len(stations))*nl:
        raise ValueError("Unexpected branch topology; no implicit branch reduction")
    rings = [[_point(graph.branch_xyz[s*nl+k][0]) for k in range(nl)] for s in range(ns-1)]
    rings.append([_point(graph.branch_xyz[(ns-2)*nl+k][1]) for k in range(nl)])
    offsets = [_dot(_sub(ring[0], rings[0][0]), regions.span_axis) for ring in rings]
    if any(b <= a for a, b in zip(offsets, offsets[1:])):
        raise ValueError("Stations must increase along the declared span axis")
    for s, ring in enumerate(rings):
        for k, point in enumerate(ring):
            expected = _add(rings[0][k], _scale(regions.span_axis, offsets[s]))
            if math.dist(point, expected) > 1e-12*max(1e-3, offsets[-1]):
                raise ValueError("Only a constant-section straight prism is supported")
    totals = {name: 0. for name, _, _ in regions.intervals()}
    by_kind = {kind: dict(total_W=0., window_W=0.) for kind in ("longitudinal", "transverse")}
    rows = []
    for index, (kind, endpoints, resistance, current) in enumerate(zip(graph.branch_kind, graph.branch_xyz, r, currents)):
        p0, p1 = map(_point, endpoints)
        patches = []
        if index < (ns-1)*nl:
            s, k = divmod(index, nl)
            if kind != "longitudinal" or p0 != rings[s][k] or p1 != rings[s+1][k]:
                raise ValueError("Longitudinal branch order/geometry mismatch")
            for neighbor in ((k-1) % nl, (k+1) % nl):
                a = _scale(_add(p0, rings[s][neighbor]), .5)
                b = _scale(_add(p1, rings[s+1][neighbor]), .5)
                patches.append((p0, a, b, p1))
        else:
            j, k = divmod(index-(ns-1)*nl, nl)
            s = stations[j]
            if kind != "transverse" or p0 != rings[s][k] or p1 != rings[s][(k+1) % nl]:
                raise ValueError("Transverse branch order/geometry mismatch")
            low = _scale(regions.span_axis, -.5*(offsets[s]-offsets[s-1]))
            high = _scale(regions.span_axis, .5*(offsets[s+1]-offsets[s]))
            patches.append((_add(p0, low), _add(p1, low), _add(p1, high), _add(p0, high)))
        measures = [regions.polygon_weights(p) for p in patches]
        area = math.fsum(m['area_m2'] for m in measures)
        power = .5*resistance*abs(current)**2
        if not math.isfinite(power):
            raise ValueError("Nonfinite branch mean power")
        weights = {name: math.fsum(m['area_m2']*m['regions'][name] for m in measures)/area
                   for name in totals}
        window = math.fsum(m['area_m2']*m['window_weight'] for m in measures)/area
        for name, weight in weights.items():
            totals[name] += power*weight
        by_kind[kind]['total_W'] += power
        by_kind[kind]['window_W'] += power*window
        rows.append(dict(branch=index, kind=kind, mean_power_W=power,
                         support_area_m2=area, window_weight=window, region_weights=weights))
    window = math.fsum(x['window_W'] for x in by_kind.values())
    return dict(schema='radia.fin.peec-loss-observables.v1', regions=regions.as_dict(),
                convention='complex peak branch current; mean power=R*abs(I)^2/2',
                status='diagnostic-not-validated',
                assumption='uniform branch loss over straight-prism dual surface supports',
                all_branch_power_W=math.fsum(x['total_W'] for x in by_kind.values()),
                window_power_W=window, region_power_W=totals,
                region_fraction={name: value/window if window > 0 else None for name, value in totals.items()},
                partition_power_defect_W=math.fsum(totals.values())-window,
                by_branch_kind=by_kind, branches=rows)


def _field_region_setup(mesh, regions, domain, label, quadrature_orders):
    """One shared set of planes for power and reactive-form measurements."""
    orders = tuple(quadrature_orders)
    if (domain not in ('surface', 'volume') or not isinstance(label, str) or not label
            or not isinstance(regions, FinMeasurementRegions)
            or not orders or any(type(x) is not int or x < 1 for x in orders)
            or tuple(sorted(set(orders))) != orders):
        raise ValueError("Explicit domain, label, regions and increasing positive quadrature orders required")
    import re
    import ngsolve as ng
    names = mesh.GetBoundaries() if domain == 'surface' else mesh.GetMaterials()
    if label not in names:
        raise ValueError(f"Missing {domain} region label {label!r}; actual={list(names)!r}")
    selection = mesh.Boundaries(re.escape(label)) if domain == 'surface' else mesh.Materials(re.escape(label))
    u = sum(((ng.x, ng.y, ng.z)[i]-regions.origin_m[i])*regions.fin_axis[i] for i in range(3))
    z = sum(((ng.x, ng.y, ng.z)[i]-regions.origin_m[i])*regions.span_axis[i] for i in range(3))
    window = ng.IfPos(regions.axial_bounds_m[0]-z, 0., 1.) * ng.IfPos(regions.axial_bounds_m[1]-z, 1., 0.)
    masks = {}
    for name, lo, hi in regions.intervals():
        mask = window
        if lo is not None:
            mask = mask * ng.IfPos(lo-u, 0., 1.)
        if hi is not None:
            mask = mask * ng.IfPos(hi-u, 1., 0.)
        masks[name] = mask
    return orders, selection, ng.BND if domain == 'surface' else ng.VOL, window, masks


def integrate_ngsolve_fin_power(mesh, mean_loss_density, regions, *,
                                domain, label, quadrature_orders):
    """Measure a solved BEM surface or FEM volume field in the common planes.

    This uses indicator quadrature on actual (possibly curved) elements, NOT
    centroid assignment or exact polygon clipping. All requested quadrature
    levels are retained: cuts inside elements require a separate integration
    convergence study. No physical error threshold or convergence PASS is set.
    Caller owns TaskManager and supplies REAL mean loss density (peak factor
    already included), a nonempty exact region label and explicit orders.
    """
    orders, region, vb, window, masks = _field_region_setup(
        mesh, regions, domain, label, quadrature_orders)
    import ngsolve as ng

    def integrate(cf, order):
        value = complex(ng.Integrate(cf, mesh, vb, definedon=region, order=order))
        if not math.isfinite(value.real) or not math.isfinite(value.imag) or value.imag != 0:
            raise ValueError(f"Expected finite REAL integral, actual={value!r}")
        return value.real

    rows = []
    for order in orders:
        measure = integrate(ng.CF(1.), order)
        if measure <= 0:
            raise ValueError(f"Empty/nonpositive {domain} measure for {label!r}: {measure}")
        full = integrate(mean_loss_density, order)
        active = integrate(mean_loss_density*window, order)
        power = {name: integrate(mean_loss_density*mask, order) for name, mask in masks.items()}
        rows.append(dict(quadrature_order=order, domain_measure=measure,
                         full_domain_power_W=full, window_power_W=active, region_power_W=power,
                         region_fraction={name: value/active if active != 0 else None for name, value in power.items()},
                         partition_power_defect_W=math.fsum(power.values())-active))
    return dict(schema='radia.fin.field-power-observables.v1', status='diagnostic-not-validated',
                domain=domain, label=label, measure_unit='m2' if domain == 'surface' else 'm3',
                regions=regions.as_dict(), integration='indicator quadrature on actual elements; cut-integration error not certified',
                quadrature_levels=rows, physics_tolerances=None)


def allocate_fin_partial_inductance(region_names, original_pair_form_H_A2, *,
                                     drive_current_peak_A, convention):
    """Allocate measured ORIGINAL self/mutual forms, never infer missing entries.

    Q[g,h] is the regional current/potential pairing, with physical prefactors
    already included. Both directions are required; reciprocity is reported,
    not assumed or silently imposed. Half of each actual pair contribution is
    allocated to each participant. Signed regional contributions are retained.
    The caller includes the whole conductor (also outside the window), and
    explicitly records the partial-potential/return convention. No 1% verdict.
    """
    names = tuple(region_names)
    if (not names or any(not isinstance(x, str) or not x for x in names)
            or len(set(names)) != len(names) or not isinstance(convention, str) or not convention):
        raise ValueError("Unique nonempty region names and explicit convention required")
    current = complex(drive_current_peak_A)
    norm = abs(current)**2
    if not math.isfinite(norm) or norm <= 0:
        raise ValueError(f"Finite nonzero drive current required: {current}")
    q = tuple(tuple(complex(x) for x in row) for row in original_pair_form_H_A2)
    if len(q) != len(names) or any(len(row) != len(names) for row in q):
        raise ValueError("Both directions of every original regional pairing required")
    if any(not math.isfinite(x.real) or not math.isfinite(x.imag) for row in q for x in row):
        raise ValueError("Nonfinite original regional pairing")
    own = {name: q[i][i].real/norm for i, name in enumerate(names)}
    allocated = dict(own)
    pairs = []
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            value = (q[i][j].real+q[j][i].real)/norm
            allocated[names[i]] += .5*value
            allocated[names[j]] += .5*value
            pairs.append(dict(regions=[names[i], names[j]], pair_contribution_H=value,
                              reciprocity_defect_H_A2=abs(q[i][j]-q[j][i].conjugate())))
    total = math.fsum(x.real for row in q for x in row)/norm
    return dict(schema='radia.fin.partial-inductance-allocation.v1',
                status='diagnostic-not-validated', convention=convention,
                normalization='complex peak drive-current amplitude',
                region_names=names, self_H=own, mutual_pairs=pairs,
                allocated_H=allocated, total_H=total,
                allocation_defect_H=math.fsum(allocated.values())-total,
                diagonal_imag_H_A2={name: q[i][i].imag for i, name in enumerate(names)},
                physics_tolerances=None)


def integrate_ngsolve_fin_impedance(mesh, current_density, original_partial_A,
                                    mean_loss_density, regions, *, domain, label,
                                    quadrature_orders, drive_current_peak_A,
                                    reactive_model, internal_reactive_density,
                                    convention):
    """Measure R and allocated partial-L on the SAME planes as regional heat.

    Caller owns TaskManager and supplies original, normalized solved fields.
    Volume model uses its FULL Green A (internal energy already present) and
    requires internal_reactive_density=None. Surface model requires an explicit
    additional real density Im(Zs)*|K|^2/omega; zero must be explicit too.
    L rows are Re(J_g^* A_total)/|I_peak|^2, the row-sum allocation of reciprocal
    original self/mutual forms. Self/mutual separation is a separate measurement
    using allocate_fin_partial_inductance, not invented from A_total alone.
    Full-domain totals refer to the selected exact label only. Reconcile all
    conductor labels and the declared terminal L separately for a device total.
    """
    current = complex(drive_current_peak_A)
    norm = abs(current)**2
    if not math.isfinite(norm) or norm <= 0 or not isinstance(convention, str) or not convention:
        raise ValueError("Nonzero finite peak drive current and explicit convention required")
    if reactive_model == 'volume-full-green':
        if domain != 'volume' or internal_reactive_density is not None:
            raise ValueError("Volume full-Green model forbids an added internal reactive term")
    elif reactive_model == 'surface-external-plus-internal':
        if domain != 'surface' or internal_reactive_density is None:
            raise ValueError("Surface model requires explicit internal reactive density")
    else:
        raise ValueError("Explicit supported reactive model required; no inference")
    orders, selection, vb, _, masks = _field_region_setup(
        mesh, regions, domain, label, quadrature_orders)
    import ngsolve as ng
    masks['outside_window'] = ng.CF(1.)-sum(masks.values())
    reactive = (ng.Conj(current_density)*original_partial_A).real
    if internal_reactive_density is not None:
        reactive = reactive+internal_reactive_density
    def integral(value, order):
        result = complex(ng.Integrate(value, mesh, vb, definedon=selection, order=order))
        if not math.isfinite(result.real) or not math.isfinite(result.imag) or result.imag != 0:
            raise ValueError(f"Expected a finite real measured form: {result}")
        return result.real
    rows = []
    for order in orders:
        power = {name: integral(mean_loss_density*mask, order) for name, mask in masks.items()}
        inductance = {name: integral(reactive*mask, order)/norm for name, mask in masks.items()}
        full_power = integral(mean_loss_density, order)
        full_l = integral(reactive, order)/norm
        rows.append(dict(quadrature_order=order, region_power_W=power,
                         region_R_equivalent_ohm={name: 2*p/norm for name, p in power.items()},
                         region_L_allocated_H=inductance,
                         full_R_equivalent_ohm=2*full_power/norm, full_L_H=full_l,
                         partition_power_defect_W=math.fsum(power.values())-full_power,
                         partition_L_defect_H=math.fsum(inductance.values())-full_l))
    return dict(schema='radia.fin.field-impedance-observables.v1',
                status='diagnostic-not-validated', regions=regions.as_dict(),
                domain=domain, label=label, reactive_model=reactive_model,
                convention=convention, drive_current_peak_A=[current.real, current.imag],
                sum_scope='selected exact label; not automatically the whole device',
                integration='actual-element indicator quadrature; cut error not certified',
                allocation='original total-potential row sum; includes mutual coupling',
                quadrature_levels=rows, physics_tolerances=None)


def peec_fin_impedance_observables(graph, resistance_diagonal_ohm,
                                   current_peak_A, regions, *,
                                   inductance_action, drive_current_peak_A,
                                   convention, operator_identity):
    """Measure a branched PEEC response without materializing its L matrix.

    The supplied explicit operator acts on ALL branches and returns L @ I in
    H A. It must include the declared internal/external inductance split; this
    measurement never supplies a missing material law. Five regional actions
    retain BOTH directions of every pairing, including outside-window current.
    Region weights partition each branch's existing dual support. Thus this is
    a fractional branch-basis allocation, not a newly assembled clipped-volume
    Green operator. Refinement against field integration remains required.
    HACApK/native/GPU callers may supply their selected action explicitly; no
    automatic backend selection or fallback and no compression certification.
    """
    import numpy as np

    if not callable(inductance_action):
        raise TypeError("Explicit inductance_action required")
    if not isinstance(operator_identity, str) or not operator_identity.strip():
        raise ValueError("Nonempty operator identity required")
    drive = complex(drive_current_peak_A)
    norm = abs(drive)**2
    if not math.isfinite(norm) or norm <= 0:
        raise ValueError("Finite nonzero peak drive current required")
    currents = np.asarray(current_peak_A, dtype=complex)
    loss = peec_fin_loss_observables(
        graph, resistance_diagonal_ohm, currents, regions)
    names = tuple(loss['region_power_W'])+('outside_window',)
    weights = np.asarray([
        [row['region_weights'][name] for name in names[:-1]]
        + [1.-row['window_weight']] for row in loss['branches']], dtype=float)
    parts = weights.T*currents[None, :]
    pair = np.empty((len(names), len(names)), dtype=complex)
    for j, source in enumerate(parts):
        # Copy: even an in-place native action cannot overwrite recorded I_g.
        response = np.asarray(inductance_action(source.copy()), dtype=complex)
        if response.shape != currents.shape or not np.isfinite(response).all():
            raise ValueError("Inductance action must return finite all-branch H A")
        pair[:, j] = parts.conjugate() @ response
    allocation = allocate_fin_partial_inductance(
        names, pair, drive_current_peak_A=drive, convention=convention)
    powers = dict(loss['region_power_W'])
    powers['outside_window'] = loss['all_branch_power_W']-loss['window_power_W']
    return dict(schema='radia.fin.peec-impedance-observables.v1',
                status='diagnostic-not-validated', operator_identity=operator_identity,
                allocation_scope='fractional dual-support branch basis, not clipped Green reassembly',
                inductance_action_calls=len(names), mean_power_W=powers,
                equivalent_R_ohm={k: 2.*v/norm for k, v in powers.items()},
                total_R_ohm=2.*loss['all_branch_power_W']/norm,
                inductance=allocation, loss=loss, physics_tolerances=None)
