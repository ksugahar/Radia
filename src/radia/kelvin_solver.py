"""kelvin_solver.py

Layer 3 of the Kelvin helper API (api_plan.md): FEM drivers.

Two drivers are provided for the 3D HCurl A-formulation on a Sugahara
two-sphere Kelvin geometry built via
``kelvin_geometry.add_kelvin_exterior_domain``:

* ``solve_full_A_kelvin`` -- meshed, volume-J source (the source
  region is part of the mesh as a separate material, e.g. ``coil``).
  Solves ``curl(nu_kelvin curl A) = J_src`` on the whole domain.

* ``solve_reduced_A_kelvin`` -- external A_s (analytic / Biot-Savart)
  source. The total A is decomposed as ``A = A_r + A_s`` with
  ``A_s`` satisfying ``curl(nu_0 curl A_s) = J_s`` in free space.
  The reduced unknown A_r is solved with the Kelvin-weighted form,
  and the RHS picks up ``-(nu_kelvin - nu_0) curl(A_s_cf) . curl(v)``
  on the Kelvin exterior material (where nu != nu_0). A_s_cf must
  already be Kelvin-aware (e.g. built via
  ``kelvin_material.make_kelvin_aware_A_s_cf``) so that inside the
  Kelvin domain it is the pulled-back field in computational coords.
"""

from __future__ import annotations

import math
import hashlib
import json
import time
from collections.abc import Mapping

import numpy as np

from ngsolve import (H1, HCurl, BilinearForm, LinearForm, GridFunction,
                      Periodic, Compress, CoefficientFunction, TaskManager,
                      curl, dx, ds, grad, InnerProduct, Conj, Integrate)

from radia.kelvin_material import make_kelvin_mu_cf, make_kelvin_nu_cf, MU_0, NU_0, _is_kelvin_material


def _canonical_array_sha256(values):
    array = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    digest = hashlib.sha256()
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _mesh_topology_geometry_sha256(mesh):
    from ngsolve import VOL

    vertices = [
        [int(vertex.nr), *[float(value) for value in vertex.point]]
        for vertex in mesh.vertices
    ]
    elements = [
        [
            int(element.nr),
            int(element.index),
            *[int(vertex.nr) for vertex in element.vertices],
        ]
        for element in mesh.Elements(VOL)
    ]
    payload = {
        "vertices": vertices,
        "volume_elements": elements,
        "materials": [str(value) for value in mesh.GetMaterials()],
        "geometry_order": int(mesh.GetCurveOrder()),
        "deformation_sha256": None,
    }
    deformation = getattr(mesh, "deformation", None)
    if deformation is not None:
        try:
            payload["deformation_sha256"] = _canonical_array_sha256(
                deformation.vec.FV().NumPy()
            )
        except (AttributeError, TypeError) as exc:
            raise ValueError(
                "cannot establish restart identity for this deformed mesh"
            ) from exc
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _identity_sha256(identity):
    return hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _validate_restart_identity_digest(value):
    canonical = str(value or "").strip().lower()
    if len(canonical) != 64 or any(
        character not in "0123456789abcdef" for character in canonical
    ):
        raise ValueError("material restart state identity_sha256 is not a SHA-256 digest")
    return canonical


def _constrains_point_gauge(mesh, selector, dirichlet_bbbnd):
    """Does a space on ``selector`` own the ``dirichlet_bbbnd`` point gauge?

    ``definedon`` keeps unused degrees of freedom in the numbering and merely
    marks them unfree, so counting free degrees of freedom with and without the
    constraint is the unambiguous ownership test.  Order 1 answers it: a BBBND
    gauge sits on vertices.
    """
    if not dirichlet_bbbnd:
        return False
    unconstrained = H1(mesh, order=1, definedon=selector)
    constrained = H1(mesh, order=1, definedon=selector,
                     dirichlet_bbbnd=dirichlet_bbbnd)
    free_unconstrained = sum(1 for flag in unconstrained.FreeDofs() if flag)
    free_constrained = sum(1 for flag in constrained.FreeDofs() if flag)
    return free_unconstrained - free_constrained > 0


def _assemble_and_solve(a_bf, f_lf, fes, inverse="pardiso"):
    """Assemble and solve.  Caller MUST be inside `with TaskManager():`
    per CLAUDE.md "Caller Wraps, Helper Does NOT" (2026-05-27).
    """
    a_bf.Assemble()
    f_lf.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec
    return gfu


def project_source_interface_potential(
        mesh, H_s, interface_boundary, *, order=2, inverse="pardiso",
        gauge_epsilon=1.0e-12, relative_tolerance=None):
    """Project the scalar source-potential trace on a source/total interface.

    A current-linked Radia/CoilBuilder field has no globally single-valued
    scalar potential.  Its tangential trace on a simply connected source/total
    interface *does* have one: in the current-free interface neighbourhood,
    ``H_s,t = -grad_Gamma(Phi_s)``.  This surface Poisson projection constructs
    that trace directly from ``H_s`` without relying on ``rad.Fld(..., "phi")``
    or choosing an arbitrary branch sheet through the total-potential region.

    The returned potential is defined only on ``interface_boundary`` and is
    suitable as ``source_potential`` for
    :func:`solve_magnetostatic_mixed_total_reduced_omega_kelvin`.  Its additive
    constant is fixed by a negligible surface mass gauge.  A non-small residual
    means that the interface has nontrivial topology or intersects current; in
    that case create an explicit cut/cohomology representation rather than
    using this trace as though it were exact.

    Caller wraps this operation in :class:`ngsolve.TaskManager`.
    """
    from ngsolve import specialcf

    if interface_boundary not in mesh.GetBoundaries():
        raise ValueError(
            f"interface_boundary={interface_boundary!r} is not a mesh boundary")
    if int(order) < 1:
        raise ValueError("order must be positive")
    if relative_tolerance is not None and (
            not math.isfinite(float(relative_tolerance)) or float(relative_tolerance) <= 0):
        raise ValueError("relative_tolerance must be positive and finite, or None for diagnostics only")
    if gauge_epsilon <= 0.0 or not math.isfinite(gauge_epsilon):
        raise ValueError("gauge_epsilon must be positive and finite")

    interface_selector = mesh.Boundaries(interface_boundary)
    fes = Compress(H1(mesh, order=int(order), definedon=interface_selector))
    potential, test = fes.TnT()
    d_interface = ds(definedon=interface_selector)
    surface_gradient = grad(potential).Trace()
    test_surface_gradient = grad(test).Trace()
    normal = specialcf.normal(mesh.dim)
    H_tangential = H_s - InnerProduct(H_s, normal) * normal

    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += InnerProduct(surface_gradient, test_surface_gradient) * d_interface
    # The physical trace is determined only up to a constant.  The gauge term
    # acts only on that null mode; it is deliberately far below discretisation
    # accuracy of the source projection.
    a_bf += float(gauge_epsilon) * potential * test * d_interface
    f_lf = LinearForm(fes)
    f_lf += -InnerProduct(H_s, test_surface_gradient) * d_interface
    a_bf.Assemble()
    f_lf.Assemble()
    potential_gf = GridFunction(fes, name="source_interface_potential")
    potential_gf.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec

    residual = grad(potential_gf).Trace() + H_tangential
    residual_norm = float(math.sqrt(Integrate(
        InnerProduct(residual, residual) * d_interface, mesh)))
    source_norm = float(math.sqrt(Integrate(
        InnerProduct(H_tangential, H_tangential) * d_interface, mesh)))
    relative_residual = residual_norm / max(source_norm, 1.0e-30)
    if relative_tolerance is not None and relative_residual > float(relative_tolerance):
        raise RuntimeError(
            "source/total interface does not admit the requested scalar source "
            f"trace: relative tangential residual={relative_residual:.3e}, "
            f"tolerance={float(relative_tolerance):.3e}; supply an explicit "
            "cut/cohomology source representation")
    return {
        "potential": potential_gf,
        "fes": fes,
        "relative_tangential_residual": relative_residual,
        "tangential_residual_norm": residual_norm,
        "tangential_source_norm": source_norm,
        "gate_enabled": relative_tolerance is not None,
        "acceptance": "passed" if relative_tolerance is not None else "not_evaluated",
    }


def project_source_physical_potential(
        mesh, H_s, physical_materials, *, order=2, inverse="pardiso",
        gauge_epsilon=1.0e-12, relative_tolerance=None):
    """Project one globally consistent source scalar potential in physical space.

    This is the permanent-magnet counterpart to
    :func:`project_source_interface_potential`.  Outside fixed magnetization
    bodies, a prescribed-magnetization field has no free-current circulation,
    so ``H_s = -grad(Phi_s)`` holds throughout the connected physical
    air/iron domain.  Projecting it once in that volume preserves the relative
    constants between separate iron-air interfaces.  Projecting each surface
    independently would erase precisely those constants and is therefore not
    a valid mixed total/reduced-Omega source contract for segmented magnets.

    ``physical_materials`` must exclude Kelvin-transformed exterior materials:
    its source field is a physical-coordinate quantity.  Current-linked coil
    sources remain on the interface/cut-cohomology path because they generally
    do not admit one globally single-valued scalar potential.

    The returned GridFunction is suitable directly as both source-potential
    traces consumed by the mixed solver.  The caller owns the surrounding
    :class:`ngsolve.TaskManager` region.
    """
    names = tuple(str(name) for name in physical_materials)
    if not names or len(names) != len(set(names)) or any(not name for name in names):
        raise ValueError("physical_materials must contain unique non-empty names")
    if int(order) < 1:
        raise ValueError("order must be positive")
    if relative_tolerance is not None and (
            not math.isfinite(float(relative_tolerance)) or float(relative_tolerance) <= 0):
        raise ValueError("relative_tolerance must be positive and finite, or None for diagnostics only")
    if gauge_epsilon <= 0.0 or not math.isfinite(gauge_epsilon):
        raise ValueError("gauge_epsilon must be positive and finite")
    actual = {str(name) for name in mesh.GetMaterials()}
    unknown = sorted(set(names) - actual)
    if unknown:
        raise ValueError(
            "physical_materials contains mesh materials that are absent: "
            f"{unknown}"
        )

    physical_selector = mesh.Materials("|".join(names))
    fes = Compress(H1(mesh, order=int(order), definedon=physical_selector))
    potential, test = fes.TnT()
    d_physical = dx(definedon=physical_selector)
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += InnerProduct(grad(potential), grad(test)) * d_physical
    # A negligible mass gauge fixes the one physical scalar-potential constant
    # without modifying the source trace at discretisation accuracy.
    a_bf += float(gauge_epsilon) * potential * test * d_physical
    f_lf = LinearForm(fes)
    f_lf += -InnerProduct(H_s, grad(test)) * d_physical
    a_bf.Assemble()
    f_lf.Assemble()
    potential_gf = GridFunction(fes, name="physical_source_potential")
    potential_gf.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec

    residual = grad(potential_gf) + H_s
    residual_norm = float(math.sqrt(Integrate(
        InnerProduct(residual, residual) * d_physical, mesh)))
    source_norm = float(math.sqrt(Integrate(
        InnerProduct(H_s, H_s) * d_physical, mesh)))
    relative_residual = residual_norm / max(source_norm, 1.0e-300)
    if relative_tolerance is not None and relative_residual > float(relative_tolerance):
        raise RuntimeError(
            "source field is not one globally exact physical scalar potential; "
            f"relative_residual={relative_residual:.3e} exceeds "
            f"relative_tolerance={float(relative_tolerance):.3e}. "
            "Use the interface trace with an explicit cut/cohomology "
            "representation for current-linked sources."
        )
    return {
        "potential": potential_gf,
        "fes": fes,
        "relative_volume_residual": relative_residual,
        "volume_residual_norm": residual_norm,
        "volume_source_norm": source_norm,
        "physical_materials": names,
        "gate_enabled": relative_tolerance is not None,
        "acceptance": "passed" if relative_tolerance is not None else "not_evaluated",
    }


SOURCE_LOADS = ("volume", "surface_flux")


def _check_source_load(name, value):
    if value not in SOURCE_LOADS:
        raise ValueError(f"{name} must be one of {SOURCE_LOADS}; got {value!r}")


def material_set_boundary_signs(mesh, materials):
    """Outward-normal sign of every boundary label that encloses ``materials``.

    ``specialcf.normal`` on a boundary face points from the face's ``domin``
    material to its ``domout`` material, so a label is ``+1`` when the set lies
    on the ``domin`` side and ``-1`` when it lies on the ``domout`` side.
    Faces with the set on both sides are interior and are not returned.  A
    label that mixes orientations, or that also names interior or unrelated
    faces, cannot be integrated with one sign and is rejected.
    """
    names = tuple(mesh.GetMaterials())
    labels = tuple(mesh.GetBoundaries())
    inside = {str(name) for name in materials}
    unknown = sorted(inside - set(names))
    if not inside or unknown:
        raise ValueError(f"materials must name mesh materials; unknown={unknown}")
    signs = {}
    other = set()
    for face in mesh.ngmesh.FaceDescriptors():
        label = labels[face.bc - 1]
        side_in = face.domin > 0 and names[face.domin - 1] in inside
        side_out = face.domout > 0 and names[face.domout - 1] in inside
        if side_in == side_out:
            other.add(label)
            continue
        sign = 1.0 if side_in else -1.0
        if signs.setdefault(label, sign) != sign:
            raise ValueError(
                f"boundary {label!r} encloses {sorted(inside)} with both orientations; "
                "split the label before using a surface-flux source load")
    mixed = sorted(other & set(signs))
    if mixed:
        raise ValueError(
            f"boundary labels {mixed} also name faces that do not enclose "
            f"{sorted(inside)}; a surface-flux source load needs exclusive labels")
    if not signs:
        raise ValueError(f"materials {sorted(inside)} have no enclosing boundary")
    return signs


def outward_source_flux_form(mesh, H_s, signs, test_trace, *, bonus_intorder, scale=1.0):
    """``scale (H_s . n_out - q) v`` on the enclosing boundary in ``signs``.

    For a divergence-free source, ``int_Omega H_s . grad(v) dx`` equals
    ``oint_dOmega (H_s . n_out) v ds``; this is that boundary form, so a
    volume source load is replaced by source evaluations on faces only.

    The net outward flux of a divergence-free field through a closed surface
    is zero, and the volume form satisfies that exactly (``grad 1 = 0``).  The
    surface quadrature does not, and on a Neumann block a nonzero net load
    drives the constant mode through its gauge.  ``q`` is the mean outward
    flux, removed so the boundary load is compatible like the volume load.
    Returns ``(form, diagnostics)``; ``relative_flux_imbalance`` is the
    removed net flux over the absolute flux, a quadrature diagnostic.
    """
    from ngsolve import specialcf

    sign = mesh.BoundaryCF(dict(signs), default=0.0)
    measure = ds(definedon=mesh.Boundaries("|".join(sorted(signs))),
                 bonus_intorder=int(bonus_intorder))
    outward = sign * InnerProduct(H_s, specialcf.normal(mesh.dim))
    area = float(Integrate(CoefficientFunction(1.0) * measure, mesh))
    from ngsolve import IfPos

    # Same measure as the assembly, so the removed mean cancels its net load.
    net = float(Integrate(outward * measure, mesh))
    absolute = float(Integrate(IfPos(outward, outward, -outward) * measure, mesh))
    mean = net / area
    form = float(scale) * (outward - mean) * test_trace * measure
    return form, {"mean_outward_flux": mean, "net_outward_flux": net,
                  "relative_flux_imbalance": abs(net) / max(absolute, 1.0e-300)}


def _evaluate_cf_at_points(mesh, field, points):
    if hasattr(field, "PrepareCache"):
        # RadiaField: one parallel batch, then cache hits during point lookup.
        field.PrepareCache(points.tolist())
    located = mesh(points[:, 0], points[:, 1], points[:, 2])
    return np.asarray(field(located), dtype=float).reshape(len(points), -1)


def interpolate_source_field(mesh, H_s, region, *, order=2):
    """Nodal H1 interpolant of a smooth vector source on ``region``.

    ``region`` is a material or boundary region (``mesh.Materials(...)`` or
    ``mesh.Boundaries(...)``).  The source is evaluated once per interpolation
    node, at the vertices and, for ``order=2``, at the edge midpoints of the
    (possibly curved) geometry, instead of at every quadrature point of every
    element.  NGSolve's H1 basis is hierarchical, so an edge coefficient is
    ``(f(mid) - (f(a) + f(b)) / 2) / c`` with ``c`` the edge shape function
    at the midpoint, measured on a straight edge rather than assumed.

    The interpolant carries an ``O(h^{order+1})`` error that the exact source
    does not: use it for loads, and the exact source where a field is reported.
    Returns the vector coefficient, its space, the node count and timings.
    """
    from ngsolve import EDGE, VERTEX, NodeId, x, y, z

    if int(order) not in (1, 2):
        raise ValueError("order must be 1 or 2")
    started = time.perf_counter()
    scalar = H1(mesh, order=int(order), definedon=region)
    active = scalar.GetDofs(region)
    vertex_nodes, vertex_dofs = [], []
    for vertex in range(mesh.nv):
        dofs = scalar.GetDofNrs(NodeId(VERTEX, vertex))
        if len(dofs) and dofs[0] >= 0 and active[dofs[0]]:
            vertex_nodes.append(vertex)
            vertex_dofs.append(dofs[0])
    if not vertex_nodes:
        raise ValueError("region has no interpolation nodes")
    vertex_points = np.asarray([mesh.vertices[v].point for v in vertex_nodes], dtype=float)
    points = [vertex_points]
    edge_rows = []
    if int(order) == 2:
        coordinate = []
        for component in (x, y, z):
            gf = GridFunction(scalar)
            gf.Set(component, definedon=region)
            coordinate.append(gf.vec.FV().NumPy().copy())
        coordinate = np.stack(coordinate, axis=1)
        position = {vertex: index for index, vertex in enumerate(vertex_nodes)}
        for edge in range(mesh.nedge):
            dofs = scalar.GetDofNrs(NodeId(EDGE, edge))
            if len(dofs) and dofs[0] >= 0 and active[dofs[0]]:
                a, b = (v.nr for v in mesh.edges[edge].vertices)
                edge_rows.append((dofs[0], position[a], position[b]))
        edge_rows = np.asarray(edge_rows, dtype=np.int64)
        curvature = coordinate[edge_rows[:, 0]]
        scale = float(np.ptp(vertex_points, axis=0).max())
        straight = np.flatnonzero(np.abs(curvature).max(axis=1) <= 1.0e-12 * scale)
        if not len(straight):
            raise RuntimeError("no straight edge to measure the edge shape value on")
        from ngsolve import ElementId

        kind = region.VB()
        mask = region.Mask()
        probe = GridFunction(scalar)
        shape_mid = 0.0
        for candidate in straight[:200]:
            row = edge_rows[candidate]
            middle = 0.5 * (vertex_points[row[1]] + vertex_points[row[2]])
            located = mesh(*middle, VOL_or_BND=kind)
            if located.nr < 0 or not mask[mesh[ElementId(kind, located.nr)].index]:
                continue
            probe.vec[:] = 0.0
            probe.vec[int(row[0])] = 1.0
            shape_mid = float(probe(located))
            break
        if not shape_mid:
            raise RuntimeError("could not measure the edge shape value inside the region")
        endpoints = 0.5 * (vertex_points[edge_rows[:, 1]] + vertex_points[edge_rows[:, 2]])
        points.append(endpoints + shape_mid * curvature)
    points = np.concatenate(points, axis=0)
    geometry_seconds = time.perf_counter() - started

    started = time.perf_counter()
    values = _evaluate_cf_at_points(mesh, H_s, points)
    evaluation_seconds = time.perf_counter() - started

    components = []
    n_vertex = len(vertex_nodes)
    for axis in range(values.shape[1]):
        gf = GridFunction(scalar)
        gf.vec[:] = 0.0
        array = gf.vec.FV().NumPy()
        array[np.asarray(vertex_dofs)] = values[:n_vertex, axis]
        if len(edge_rows):
            linear = 0.5 * (values[edge_rows[:, 1], axis] + values[edge_rows[:, 2], axis])
            array[edge_rows[:, 0]] = (values[n_vertex:, axis] - linear) / shape_mid
        components.append(gf)
    return {
        "field": CoefficientFunction(tuple(components)),
        "space": scalar,
        "nodes": int(len(points)),
        "order": int(order),
        "timings_seconds": {"geometry": geometry_seconds, "source_evaluation": evaluation_seconds},
    }


def _material_centroid(mesh, material):
    from ngsolve import VOL

    for element in mesh.Elements(VOL):
        if element.mat == material:
            points = [mesh.vertices[vertex.nr].point for vertex in element.vertices]
            return tuple(sum(point[axis] for point in points) / len(points) for axis in range(3))
    raise ValueError(f"material {material!r} has no volume element")


def _uniform_permeability(mesh, mu_cf, materials, mu_r_by_material, *, role):
    """Constant ``mu`` shared by ``materials``, checked against ``mu_cf``.

    The surface-flux identity moves ``mu`` outside the integral, so it is
    valid only where ``mu`` is one constant.  The declared value is taken
    from ``mu_r_by_material`` (default 1) and verified at each material's
    element centroid against the permeability actually assembled.
    """
    values = {name: MU_0 * float((mu_r_by_material or {}).get(name, 1.0)) for name in materials}
    if len(set(values.values())) != 1:
        raise ValueError(f"surface-flux {role} load needs one permeability; got {values}")
    mu = next(iter(values.values()))
    for name in materials:
        point = _material_centroid(mesh, name)
        actual = float(mu_cf(mesh(*point)))
        if not math.isclose(actual, mu, rel_tol=1e-12):
            raise ValueError(
                f"surface-flux {role} load: assembled permeability in {name!r} is "
                f"{actual!r}, declared {mu!r}; the identity needs a constant mu")
    return mu


def project_source_total_hodge(
        mesh, H_s, total_source_materials, *, order=2, inverse="pardiso",
        gauge_epsilon=1.0e-12, bonus_intorder=4, source_load="volume",
        tangential_tolerance=None):
    """Split a linked source inside total-potential materials.

    On a multiply connected iron body a curl-free coil field need not be the
    gradient of one single-valued scalar.  This volume projection computes

    ``H_s = -grad(Phi_s) + H_harmonic``.

    ``Phi_s`` supplies the reduced/total interface jump and ``H_harmonic``
    carries the non-exact cohomology class in the total region.  The caller
    owns the surrounding :class:`ngsolve.TaskManager` region.
    ``bonus_intorder`` controls projection assembly and norm diagnostics;
    use the mixed solver's value when composing the two operations. Matching
    bonuses does not establish source quadrature convergence or orthogonality
    against a different test space.

    ``source_load="surface_flux"`` assembles the same right-hand side from
    the source normal flux on the enclosing boundary (``div H_s = 0``), so the
    source is evaluated on faces only.  That load determines the exact part
    alone: a harmonic (linked-current) remainder is not represented by
    boundary normal data.  It therefore requires ``tangential_tolerance`` and
    fails when the boundary tangential residual
    ``|n x (H_s + grad Phi_s)| / |n x H_s|`` exceeds it, which is what a
    linked source produces.  The volume norms are not evaluated in this mode.
    """
    _check_source_load("source_load", source_load)
    if source_load == "surface_flux" and (
            tangential_tolerance is None or not math.isfinite(float(tangential_tolerance))
            or float(tangential_tolerance) <= 0.0):
        raise ValueError(
            "surface_flux Hodge projection needs a positive tangential_tolerance; "
            "without it a linked-source harmonic remainder would be dropped silently")
    names = tuple(str(name) for name in total_source_materials)
    if not names or len(names) != len(set(names)) or any(not name for name in names):
        raise ValueError(
            "total_source_materials must contain unique non-empty names")
    if int(order) < 1:
        raise ValueError("order must be positive")
    if isinstance(bonus_intorder, bool) or int(bonus_intorder) != bonus_intorder or bonus_intorder < 0:
        raise ValueError("bonus_intorder must be a nonnegative integer")
    if gauge_epsilon <= 0.0 or not math.isfinite(gauge_epsilon):
        raise ValueError("gauge_epsilon must be positive and finite")
    actual = {str(name) for name in mesh.GetMaterials()}
    unknown = sorted(set(names) - actual)
    if unknown:
        raise ValueError(
            "total_source_materials contains mesh materials that are absent: "
            f"{unknown}"
        )

    selector = mesh.Materials("|".join(names))
    fes = Compress(H1(mesh, order=int(order), definedon=selector))
    potential, test = fes.TnT()
    d_total = dx(definedon=selector, bonus_intorder=int(bonus_intorder))
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += InnerProduct(grad(potential), grad(test)) * d_total
    a_bf += float(gauge_epsilon) * potential * test * d_total
    f_lf = LinearForm(fes)
    if source_load == "volume":
        f_lf += -InnerProduct(H_s, grad(test)) * d_total
    else:
        signs = material_set_boundary_signs(mesh, names)
        flux_form, flux_balance = outward_source_flux_form(
            mesh, H_s, signs, test.Trace(), bonus_intorder=bonus_intorder, scale=-1.0)
        f_lf += flux_form
    a_bf.Assemble()
    f_lf.Assemble()
    potential_gf = GridFunction(fes, name="total_source_potential")
    potential_gf.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec

    harmonic_field = H_s + grad(potential_gf)
    if source_load == "surface_flux":
        from ngsolve import specialcf

        normal = specialcf.normal(mesh.dim)
        d_boundary = ds(definedon=mesh.Boundaries("|".join(sorted(signs))),
                        bonus_intorder=int(bonus_intorder))
        source_t = H_s - InnerProduct(H_s, normal) * normal
        residual_t = source_t + grad(potential_gf).Trace()
        residual_t = residual_t - InnerProduct(residual_t, normal) * normal
        residual_norm = float(math.sqrt(Integrate(
            InnerProduct(residual_t, residual_t) * d_boundary, mesh)))
        source_t_norm = float(math.sqrt(Integrate(
            InnerProduct(source_t, source_t) * d_boundary, mesh)))
        relative_t = residual_norm / max(source_t_norm, 1.0e-300)
        if relative_t > float(tangential_tolerance):
            raise RuntimeError(
                "surface_flux Hodge projection: boundary tangential residual "
                f"{relative_t:.3e} exceeds {float(tangential_tolerance):.3e}; the "
                "source is not exact on this total region (linked current?). Use "
                "source_load='volume', which retains the harmonic remainder")
        return {
            "potential": potential_gf,
            "harmonic_field": harmonic_field,
            "fes": fes,
            "relative_harmonic_norm": None,
            "harmonic_norm": None,
            "source_norm": None,
            "relative_tangential_residual": relative_t,
            "tangential_tolerance": float(tangential_tolerance),
            "boundary_signs": dict(signs),
            "flux_balance": flux_balance,
            "total_source_materials": names,
            "bonus_intorder": int(bonus_intorder),
            "source_load": source_load,
        }
    harmonic_norm = float(math.sqrt(Integrate(
        InnerProduct(harmonic_field, harmonic_field) * d_total, mesh)))
    source_norm = float(math.sqrt(Integrate(
        InnerProduct(H_s, H_s) * d_total, mesh)))
    return {
        "potential": potential_gf,
        "harmonic_field": harmonic_field,
        "fes": fes,
        "relative_harmonic_norm": harmonic_norm / max(source_norm, 1.0e-300),
        "harmonic_norm": harmonic_norm,
        "source_norm": source_norm,
        "total_source_materials": names,
        "bonus_intorder": int(bonus_intorder),
        "source_load": source_load,
    }


def project_kelvin_A_source(mesh, A_s_cf, *, order=2):
    """Project a Kelvin-aware source potential into a periodic HCurl space.

    A Radia-backed coefficient provides exact values but deliberately has no
    symbolic derivative. The source curl used by the reduced-A and reduced
    Omega--Omega weak forms must therefore be the curl of this *same*
    conforming projection. Callers comparing formulations must create it once
    and pass the returned grid function to both solvers.

    Caller wraps this operation in :class:`ngsolve.TaskManager`.
    """
    source_fes = Periodic(HCurl(mesh, order=int(order)))
    source = GridFunction(source_fes, name="kelvin_source_A")
    source.Set(A_s_cf)
    return source


def solve_full_A_kelvin(mesh, J_source_cf, R_K, offset,
                         source_material="coil",
                         nu_0=NU_0,
                         order=1,
                         dirichlet_bbnd="GND",
                         gauge_eps=1e-6,
                         bonus_intorder=4,
                         kelvin_mats=("kelvin",),
                         inverse="pardiso"):
    """Full-A 3D HCurl FEM with Sugahara Kelvin convention.

    Solves ``curl(nu_kelvin curl A) = J_src`` on the two-sphere
    domain, with ``J_src`` supported on the mesh region
    ``source_material``.

    Args:
        mesh: NGSolve ``Mesh`` built from
            ``add_kelvin_exterior_domain`` output.
        J_source_cf: vector CF for J_src (A/m^2).
        R_K, offset: Kelvin sphere parameters (matching those used in
            ``add_kelvin_exterior_domain``).
        source_material: mesh material name carrying the volume current
            (default ``"coil"``).
        nu_0: vacuum reluctivity (default 1/mu_0).
        order: HCurl polynomial order.
        dirichlet_bbnd: BBND name for the Dirichlet GND vertex
            (A = 0 at the Kelvin sphere center).
        gauge_eps: mass regularization coefficient to fix the HCurl
            gauge (factor of nu_0 scaled).
        bonus_intorder: extra quadrature order for the Kelvin domain.
        kelvin_mats: substrings identifying Kelvin exterior materials
            (forwarded to ``make_kelvin_nu_cf``).
        inverse: NGSolve sparse solver name.

    Returns:
        dict with keys ``gfu`` (GridFunction, total A), ``fes``
        (Periodic HCurl space), ``nu_cf`` (Kelvin-modulated nu CF).
    """
    nu_cf = make_kelvin_nu_cf(mesh, R_K, offset, nu_0=nu_0,
                                kelvin_mats=kelvin_mats)
    fes = Periodic(HCurl(mesh, order=order, dirichlet_bbnd=dirichlet_bbnd))
    u, v = fes.TnT()
    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=bonus_intorder)
    a_bf += gauge_eps * nu_0 * u * v * dx
    f_lf = LinearForm(fes)
    f_lf += J_source_cf * v * dx(source_material)

    gfu = _assemble_and_solve(a_bf, f_lf, fes, inverse=inverse)
    return {"gfu": gfu, "fes": fes, "nu_cf": nu_cf}


def solve_reduced_A_kelvin(mesh, A_s_cf, R_K, offset,
                            nu_0=NU_0,
                            order=1,
                            dirichlet_bbnd="GND",
                            gauge_eps=1e-6,
                            bonus_intorder=4,
                            kelvin_mats=("kelvin",),
                            inverse="pardiso"):
    """Reduced-A 3D HCurl FEM: A = A_r + A_s with external A_s.

    The source A_s is supplied as a CF that is ALREADY Kelvin-aware
    (built via ``make_kelvin_aware_A_s_cf``): it returns the physical
    A in non-Kelvin materials and the pulled-back A' in Kelvin
    materials. The reduced unknown A_r satisfies (CORRECTED 2026-05-04)

        a(A_r, v) = - int_kext nu_kelvin * curl(A_s_cf) . curl(v) dV

    See docs/kelvin/KELVIN_TRANSFORMATION.md §7.5 for derivation.

    The previous form ``-(nu - nu_0) * curl(A_s) * curl(v) dx`` is INVALID
    when A_s is a Kelvin pullback because the pullback satisfies nu'
    Maxwell (NOT nu_0 Maxwell) in the Kelvin region; the (nu - nu_0)
    simplification requires nu_0 Maxwell globally. The bug previously
    caused +43% inductance error vs +6% with the correct form on a
    PEEC torus benchmark.

    The total vector potential is ``A_total = A_r + A_s``.

    Args:
        mesh, R_K, offset, nu_0, order, dirichlet_bbnd, gauge_eps,
        bonus_intorder, kelvin_mats, inverse: same as
        ``solve_full_A_kelvin``.
        A_s_cf: Kelvin-aware A_s CoefficientFunction.

    Returns:
        dict with keys ``gfu_r`` (reduced A_r GridFunction),
        ``fes``, ``nu_cf``, ``A_s_cf`` (passthrough for convenience).
    """
    nu_cf = make_kelvin_nu_cf(mesh, R_K, offset, nu_0=nu_0,
                                kelvin_mats=kelvin_mats)
    fes = Periodic(HCurl(mesh, order=order, dirichlet_bbnd=dirichlet_bbnd))
    u, v = fes.TnT()
    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=bonus_intorder)
    a_bf += gauge_eps * nu_0 * u * v * dx
    f_lf = LinearForm(fes)
    # RHS (CORRECTED 2026-05-04): - int_kext nu' curl(A_s) . curl(v) dV
    # The OLD form -(nu - nu_0) curl(A_s) . curl(v) dx is WRONG when A_s
    # is a Kelvin pullback (Convention A): the pullback satisfies nu'
    # Maxwell, not nu_0 Maxwell, breaking the (nu - nu_0) simplification.
    # See docs/kelvin/KELVIN_TRANSFORMATION.md §7.5.
    kelvin_mat_str = "|".join(kelvin_mats)
    f_lf += -nu_cf * curl(A_s_cf) * curl(v) \
        * dx(kelvin_mat_str, bonus_intorder=bonus_intorder)

    gfu_r = _assemble_and_solve(a_bf, f_lf, fes, inverse=inverse)
    return {"gfu_r": gfu_r, "fes": fes, "nu_cf": nu_cf, "A_s_cf": A_s_cf}


def solve_magnetostatic_reduced_A_kelvin(
        mesh, A_s, R_K, offset, *, mu_r_by_material,
        nu_0=NU_0, order=1, dirichlet_bbnd="GND", gauge_eps=1e-6,
        bonus_intorder=4, kelvin_mats=("kelvin",), inverse="pardiso"):
    """Solve the linear reduced-A magnetostatic problem with iron and Kelvin.

    This is the production three-dimensional route for a compact external
    current source and linear magnetic materials. ``A_s`` must be the
    :class:`ngsolve.GridFunction` returned by :func:`project_kelvin_A_source`.
    Construct its input coefficient with
    :func:`radia.kelvin_material.make_kelvin_aware_radia_A_s_cf` for a Radia
    source object.

    Let ``nu_source`` be vacuum reluctivity in the physical model and the
    Kelvin-transformed vacuum metric in the exterior. Since ``A_s`` solves the
    source-only problem with ``nu_source``, the reaction potential satisfies

    ``curl(nu curl(A_r)) = curl((nu_source - nu) curl(A_s))``.

    The RHS therefore exists only where the physical material differs from
    vacuum, while the Kelvin source is already carried by the pullback. This
    distinction is essential: treating Kelvin as an unassembled region or
    applying a physical-space source there gives a finite-domain surrogate,
    not an open-boundary reduced-A solve.

    Args:
        mesh: two-sphere Kelvin mesh with periodic point identifications.
        A_s_cf: Kelvin-aware source vector potential.
        R_K, offset: radius and translated centre of the Kelvin sphere.
        mu_r_by_material: mapping from physical mesh material name to positive
            scalar relative permeability. Kelvin materials are rejected by
            :func:`make_kelvin_nu_cf` because their metric is prescribed.
    """
    nu_source_cf = make_kelvin_nu_cf(
        mesh, R_K, offset, nu_0=nu_0, kelvin_mats=kelvin_mats)
    nu_cf = make_kelvin_nu_cf(
        mesh, R_K, offset, nu_0=nu_0, kelvin_mats=kelvin_mats,
        mu_r_by_material=mu_r_by_material)
    fes = Periodic(HCurl(mesh, order=order, dirichlet_bbnd=dirichlet_bbnd))
    u, v = fes.TnT()
    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=bonus_intorder)
    a_bf += gauge_eps * nu_0 * u * v * dx
    f_lf = LinearForm(fes)
    f_lf += ((nu_source_cf - nu_cf) * curl(A_s) * curl(v)
             * dx(bonus_intorder=bonus_intorder))

    gfu_r = _assemble_and_solve(a_bf, f_lf, fes, inverse=inverse)
    return {
        "gfu_r": gfu_r,
        "fes": fes,
        "nu_cf": nu_cf,
        "nu_source_cf": nu_source_cf,
        "A_s": A_s,
        # A_s and A_r intentionally live in separately constructed periodic
        # spaces (the source has no GND restriction). NGSolve cannot take the
        # curl of their symbolic sum, but curl is linear and each term has the
        # correct HCurl differential operator on its own space.
        "B_cf": curl(gfu_r) + curl(A_s),
    }


def solve_magnetostatic_reduced_omega_kelvin(
        mesh, H_s, R_K, offset, *, mu_r_by_material,
        order=1, dirichlet_bbbnd="GND", bonus_intorder=4,
        kelvin_mats=("kelvin",), inverse="pardiso"):
    """Retired plain reduced-Omega Kelvin entry point; always raises.

    The legacy exterior/interface convention is not validated. Use
    :func:`solve_magnetostatic_mixed_total_reduced_omega_kelvin` with explicit
    material partitions, interface traces and Kelvin source data instead.
    This is not a drop-in substitution: missing physical input must not be
    inferred. The retirement does not reject reduced-potential formulations
    in general. Historical evidence lives in
    ``validation_test/hdiv_vim/reduced_omega_retirement``.
    """
    raise NotImplementedError(
        "solve_magnetostatic_reduced_omega_kelvin is retired: the legacy "
        "Kelvin exterior/interface convention is not validated. Migrate to "
        "solve_magnetostatic_mixed_total_reduced_omega_kelvin with explicit "
        "material partitions, source-potential traces and Kelvin interface "
        "data; automatic substitution is not supported.")


def _report_direct_solve_residual(a_bf, f_lf, solution, fes, block_names=None):
    """Measure ``r = b - A x`` for a solve that has no iteration history.

    A direct factorisation is not exact arithmetic, and reporting nothing
    because there is no iteration count is the same gap as treating a CG that
    ran out of iterations as converged.  Norms are reported over the FREE dofs
    (the system that was solved) and, for a compound space, per block, so an
    interface-constraint defect cannot hide inside a healthy PDE residual.
    """
    import numpy as _np

    residual = solution.vec.CreateVector()
    residual.data = f_lf.vec - a_bf.mat * solution.vec
    r = _np.asarray(residual.FV(), dtype=float).copy()
    b = _np.asarray(f_lf.vec.FV(), dtype=float).copy()
    free = _np.fromiter((bool(bit) for bit in fes.FreeDofs()), bool, len(r))

    def _relative(mask):
        if not mask.any():
            return None
        numerator = float(_np.linalg.norm(r[mask]))
        denominator = float(_np.linalg.norm(b[mask]))
        return {"residual_l2": numerator,
                "rhs_l2": denominator,
                "relative": (numerator / denominator if denominator > 0.0
                             else None)}

    report = {"free_dofs": _relative(free),
              "all_dofs": _relative(_np.ones_like(free))}
    try:
        component_count = len(solution.components)
    except TypeError:                                   # not a compound space
        component_count = 0
    names = tuple(block_names or ()) or tuple(
        "block_%d" % index for index in range(component_count))
    blocks = {}
    for index in range(component_count):
        span = fes.Range(index)
        mask = _np.zeros_like(free)
        mask[span.start:span.stop] = True
        blocks[names[index] if index < len(names) else "block_%d" % index] = (
            _relative(mask & free))
    if blocks:
        report["blocks"] = blocks
    return report


def _assembled_primal_energy(a_bf, f_lf, solution, fes, primal_blocks):
    """The functional the discrete solve actually minimises, from the system.

    The solve is a saddle point of the Lagrangian; its primal part is
    ``J(x_P) = 1/2 x_P^T A_PP x_P - b_P^T x_P`` over the potential blocks,
    minimised subject to the interface constraint. This reports the ASSEMBLED
    functional, rather than an energy recomputed with different quadrature.
    Monotonicity across orders additionally requires nested admissible sets
    and a common functional, including consistent quadrature across spaces.

    The multiplier entries are zeroed before the product, so only primal rows
    and columns contribute and the constraint right-hand side drops out.
    """
    import numpy as _np

    x_p = solution.vec.CreateVector()
    x_p.data = solution.vec
    values = x_p.FV().NumPy()
    keep = _np.zeros(len(values), dtype=bool)
    for index in primal_blocks:
        span = fes.Range(index)
        keep[span.start:span.stop] = True
    values[~keep] = 0.0
    product = x_p.CreateVector()
    product.data = a_bf.mat * x_p
    half_quadratic = 0.5 * float(_np.dot(x_p.FV().NumPy(), product.FV().NumPy()))
    linear = float(_np.dot(_np.asarray(f_lf.vec.FV()), x_p.FV().NumPy()))
    return {"half_xAx": half_quadratic,
            "b_dot_x": linear,
            "energy": half_quadratic - linear,
            "primal_blocks": list(primal_blocks)}


def solve_magnetostatic_mixed_total_reduced_omega_kelvin(
        mesh, H_s, source_potential, R_K, offset, *, mu_r_by_material=None,
        reduced_materials, total_materials, interface_boundary,
        order=1, dirichlet_bbbnd="GND", bonus_intorder=4,
        kelvin_mats=("kelvin",), inverse="pardiso",
        interface_constraint_scale=None, total_dirichlet_cf=None,
        mu_cf=None, kelvin_interface_boundary=None,
        kelvin_source_potential=None, kelvin_source_h=None,
        total_source_h=None, total_source_materials=(), return_system=False,
        phase_callback=None, source_rhs_reduced=None,
        reduced_normal_flux=None, reduced_flux_boundary=None,
        total_normal_flux=None, total_flux_boundary=None,
        reduced_dirichlet_boundary=None, total_dirichlet_boundary=None,
        interface_multiplier_dirichlet_boundary=None,
        reduced_zero_normal_boundary=None, surface_dirichlet=None, _fixed_rhs_cache=None,
        _rhs_material=None, kelvin_match_exact=False,
        reduced_source_load="volume", total_source_load="volume",
        total_source_potential=None, load_source_h=None):
    """Solve the mixed total/reduced Omega formulation.

    The total/reduced scalar-potential split of Simkin and Trowbridge
    (``simkin1979use``, doi:10.1002/nme.1620140308; nonlinear
    three-dimensional treatment in ``simkin1980three``,
    doi:10.1049/ip-b.1980.0052): a total potential where the material is
    permeable, a reduced potential carrying ``H_s`` where the source lives,
    and an interface condition between them.

    ``return_system=True`` retains assembled forms for explicit diagnostics.
    The default returns ``system=None`` to avoid retaining their matrix storage.
    ``phase_callback`` optionally receives start/complete dictionaries for
    matrix assembly, source RHS assembly, factorization and backsolve.
    Completed subphase durations are returned in ``phase_timings_seconds``;
    these do not include mesh, space setup or postprocessing time.

    ``H_s`` is used in ``reduced_materials`` (the source enclosure), where
    ``H = H_s - grad(phi_reduced)``.  A linked coil may additionally supply
    ``total_source_h`` on ``total_source_materials``; this is the harmonic
    remainder of the iron-volume Hodge split and gives
    ``H = total_source_h - grad(phi_total)`` there.  Other total materials use
    ``H = -grad(phi_total)``.  On ``interface_boundary`` the potentials are
    coupled by the exact source trace

    ``phi_total - phi_reduced = source_potential``.

    When the reduced physical air reaches the inner Kelvin sphere, provide
    ``kelvin_interface_boundary`` and ``kelvin_source_potential`` as well; a
    Kelvin material in ``total_materials`` without them is rejected.  The
    source enclosure and the Kelvin exterior are then discretised by one
    periodic H1 space, which is what makes the exterior transmit.  The magnetic
    scalar potential is a twisted 0-form, so the orientation-reversing Kelvin
    inversion pulls it back as ``Omega_comp = -Omega_phys`` with no metric
    factor (``docs/kelvin/KELVIN_TRANSFORMATION.md`` 2.3).  Solving for
    ``-Omega_comp`` in the exterior leaves one unknown that is continuous
    across the identified spheres up to the source jump

    ``Omega_exterior - phi_reduced = kelvin_source_potential``,

    which enters as a lift of that trace rather than as a second Lagrange
    multiplier.  A space restricted to the total materials cannot express this
    coupling at all: the material touching the physical sphere is the reduced
    one, so such a space owns no degree of freedom there, the identification
    pairs nothing, and the exterior block is eliminated with an identically
    zero field.  ``H_cf`` in a Kelvin material is the computational-frame
    ``H_comp``; its sign follows the twisted convention, so it is ``+grad`` of
    the stored unknown, not ``-grad``.
    ``B_cf = mu_cf * H_cf`` is likewise in computational coordinates in the
    Kelvin exterior, not a directly sampled physical-space exterior field.

    This prevents the source field and the reduced correction from cancelling
    inside high-permeability material.  The normal flux condition is natural
    in the variational form; the scalar-trace constraint supplies the
    tangential-H condition.  It is the finite-element form of the total /
    reduced Omega split of Simkin and Trowbridge (1980), ``simkin1980three``.

    ``source_potential`` is intentionally mandatory.  It must satisfy
    ``H_s = -grad(source_potential)`` in the current-free neighbourhood of the
    source/total interface.  Do not reconstruct it from an HDiv-projected B
    field.  Linked filament coils require an explicit cut or a cohomology
    representative before they provide such a single-valued trace.

    The formulation has an interface Lagrange multiplier and is symmetric
    indefinite, so the default direct PARDISO solve is deliberate.  The
    returned field is continuous in the physical tangential/normal sense but
    not represented as one global H1 GridFunction.

    Args:
        mesh: Kelvin-periodic NGSolve mesh.
        H_s: physical-coordinate source H coefficient on
            ``reduced_materials``.  The production split keeps every Kelvin
            material in ``total_materials``, so no source evaluation is made
            in the Kelvin exterior.
        source_potential: physical scalar-potential trace on
            ``interface_boundary``.
        R_K, offset: Kelvin sphere radius and translated centre.
        mu_r_by_material: positive relative permeability by material name.
            Ignored when ``mu_cf`` is supplied by a nonlinear outer iteration.
        reduced_materials: exact material names for the source enclosure.
        total_materials: exact names for iron, ordinary air, and Kelvin.
        interface_boundary: named source/total internal boundary.
        kelvin_interface_boundary: inner physical Kelvin boundary that is
            periodically paired with the Kelvin exterior.  Supply this with
            ``kelvin_source_potential`` when a reduced material touches the
            inner Kelvin sphere.
        kelvin_source_potential: physical source-potential trace on
            ``kelvin_interface_boundary``.  It is lifted into the shared
            periodic space, so the exterior physical potential leaves the
            sphere as ``phi_reduced + kelvin_source_potential``.
        total_dirichlet_cf: optional value on the total-region BBBND point
            gauge named by ``dirichlet_bbbnd``. This is NOT an exterior-surface
            Dirichlet condition. Surface Dirichlet data on total and reduced
            potentials require different source-potential lifts and are not
            exposed by this parameter. Production Kelvin meshes use ``None``.
        surface_dirichlet: finite-domain mapping ``{"reduced": {label: value},
            "total": {label: value}}`` of potential values on exact exterior
            boundary labels. Replaces the point gauge; the caller must supply
            the correct source lift for reduced values. Not supported with a
            Kelvin exterior. Passed unchanged through P1/P2 Picard solves.
        mu_cf: optional fully Kelvin-aware permeability coefficient. This is
            the narrow extension point used by the nonlinear Picard driver;
            callers must not supply a physical-space coefficient in the Kelvin
            material.
        total_source_h: optional non-exact harmonic source field retained in
            the total-potential region.
        total_source_materials: total-region materials on which
            ``total_source_h`` is defined.  Supply both arguments together.
        source_rhs_reduced: optional preassembled reduced-region volume load
            in ``fes_reduced`` DOF order. The caller must establish the same
            source, permeability and quadrature contract as the omitted
            ``mu_cf * H_s * grad(test_reduced)`` volume term. Interface and
            exterior loads remain assembled normally.
        reduced_normal_flux: prescribed outward normal B on the named physical
            boundary of the reduced region. Without it the natural condition
            is zero normal flux. Supply ``reduced_flux_boundary`` with it.
        total_normal_flux: analogous outward normal B for the total region.
            Supply ``total_flux_boundary`` with it.
        reduced_dirichlet_boundary: physical boundary where the reduced
            correction potential is zero (not the total scalar potential).
        total_dirichlet_boundary: physical boundary where the total scalar
            potential is zero. These selectors do not identify the two fields.
        interface_multiplier_dirichlet_boundary: optional boundary on which
            the interface constraint is already fixed by essential traces.
            Use only after checking the intersecting interface DOFs.
        reduced_source_load: ``"volume"`` (default) assembles
            ``mu H_s . grad(v)`` over the reduced region.  ``"surface_flux"``
            assembles the equal boundary form ``mu (H_s . n_out) v`` on the
            region's enclosing labels, valid because ``div H_s = 0`` and the
            reduced permeability is one constant (checked).  It evaluates the
            source on faces only; the discrete loads differ by quadrature
            error, so compare before substituting one for the other.
        total_source_load: the same choice for ``total_source_h`` on
            ``total_source_materials``.  ``"surface_flux"`` needs
            ``total_source_h`` to be the raw source ``H_s`` and
            ``total_source_potential`` the Hodge potential ``Phi_s`` of
            ``H_s + grad(Phi_s)``, and a constant permeability per material
            (linear only; a supplied ``mu_cf`` must match it).
        load_source_h: optional representation of ``H_s`` used only in the
            reduced load (for example :func:`interpolate_source_field`).  The
            reported reduced field still uses the exact ``H_s``.
    """
    _check_source_load("reduced_source_load", reduced_source_load)
    _check_source_load("total_source_load", total_source_load)
    if (total_source_load == "surface_flux") != (total_source_potential is not None):
        raise ValueError(
            "total_source_potential is required by, and only used with, "
            "total_source_load='surface_flux'")
    if total_source_load == "surface_flux" and _rhs_material is not None:
        raise ValueError("the cached material RHS operator is a volume form; "
                         "total_source_load='surface_flux' is linear only")
    if reduced_source_load == "surface_flux" and source_rhs_reduced is not None:
        raise ValueError("source_rhs_reduced replaces the reduced load; do not combine it "
                         "with reduced_source_load='surface_flux'")
    reduced_materials = tuple(reduced_materials)
    total_materials = tuple(total_materials)
    actual_materials = set(mesh.GetMaterials())
    reduced_set = set(reduced_materials)
    total_set = set(total_materials)
    if not reduced_set:
        raise ValueError("reduced_materials must name the current-source enclosure")
    if not total_set:
        raise ValueError("total_materials must name the total-potential region")
    if reduced_set & total_set:
        raise ValueError("reduced_materials and total_materials must be disjoint")
    if reduced_set | total_set != actual_materials:
        missing = sorted(actual_materials - (reduced_set | total_set))
        unknown = sorted((reduced_set | total_set) - actual_materials)
        raise ValueError(
            "reduced_materials and total_materials must partition mesh materials; "
            f"missing={missing}, unknown={unknown}")
    if interface_boundary not in mesh.GetBoundaries():
        raise ValueError(
            f"interface_boundary={interface_boundary!r} is not a mesh boundary")
    if kelvin_source_h is not None and kelvin_source_potential is not None:
        raise ValueError(
            "supply either kelvin_source_h (the exact pulled-back exterior "
            "source) or kelvin_source_potential (its projected interface "
            "trace), not both")
    if (kelvin_interface_boundary is None) != (
            kelvin_source_potential is None and kelvin_source_h is None):
        raise ValueError(
            "kelvin_interface_boundary and one exterior source representation "
            "must be supplied together")
    if (kelvin_interface_boundary is not None
            and kelvin_interface_boundary not in mesh.GetBoundaries()):
        raise ValueError(
            f"kelvin_interface_boundary={kelvin_interface_boundary!r} is not "
            "a mesh boundary")
    total_source_materials = tuple(str(name) for name in total_source_materials)
    if (total_source_h is None) != (not total_source_materials):
        raise ValueError(
            "total_source_h and total_source_materials must be supplied together")
    if not set(total_source_materials) <= total_set:
        raise ValueError(
            "total_source_materials must be contained in total_materials")
    if interface_constraint_scale is None:
        # Stiffness entries scale as mu * L; interface trace entries scale as
        # L^2.  This balancing does not alter the constraint, only the saddle
        # matrix conditioning seen by the direct solver.
        interface_constraint_scale = (1.0 / NU_0) / float(R_K)
    if interface_constraint_scale <= 0.0 or not math.isfinite(interface_constraint_scale):
        raise ValueError("interface_constraint_scale must be positive and finite")
    if (reduced_normal_flux is None) != (reduced_flux_boundary is None):
        raise ValueError("reduced_normal_flux and reduced_flux_boundary must be supplied together")
    if reduced_flux_boundary is not None:
        flux_names = str(reduced_flux_boundary).split("|")
        if not flux_names or not set(flux_names) <= set(mesh.GetBoundaries()):
            raise ValueError("reduced_flux_boundary must name existing boundaries")
    if (total_normal_flux is None) != (total_flux_boundary is None):
        raise ValueError("total_normal_flux and total_flux_boundary must be supplied together")
    if total_flux_boundary is not None:
        flux_names = str(total_flux_boundary).split("|")
        if not flux_names or not set(flux_names) <= set(mesh.GetBoundaries()):
            raise ValueError("total_flux_boundary must name existing boundaries")
    for selector, label in ((reduced_dirichlet_boundary, "reduced_dirichlet_boundary"),
                            (total_dirichlet_boundary, "total_dirichlet_boundary"),
                            (interface_multiplier_dirichlet_boundary,
                             "interface_multiplier_dirichlet_boundary")):
        if selector is not None:
            names = str(selector).split("|")
            if not names or not set(names) <= set(mesh.GetBoundaries()):
                raise ValueError(f"{label} must name existing boundaries")

    kelvin_mats = tuple(kelvin_mats)
    reduced_kelvin = sorted(
        material for material in reduced_set
        if _is_kelvin_material(material, kelvin_mats, exact=kelvin_match_exact))
    if reduced_kelvin:
        raise ValueError(
            f"Kelvin selectors must not match reduced materials: {reduced_kelvin}")
    if mu_cf is None:
        mu_cf = make_kelvin_mu_cf(
            mesh, R_K, offset, kelvin_mats=kelvin_mats,
            mu_r_by_material=mu_r_by_material, kelvin_match_exact=kelvin_match_exact)
    kelvin_total_materials = tuple(
        material for material in total_materials
        if _is_kelvin_material(material, kelvin_mats, exact=kelvin_match_exact))
    core_total_materials = tuple(
        material for material in total_materials
        if material not in kelvin_total_materials)
    if kelvin_total_materials and kelvin_interface_boundary is None:
        raise ValueError(
            "a Kelvin exterior material is declared in total_materials without "
            "kelvin_interface_boundary; the reduced source enclosure is then "
            "the material that touches the identified physical sphere, the "
            "total space owns no degree of freedom there, and the Kelvin block "
            "stays uncoupled with an identically zero field: "
            f"kelvin_materials={list(kelvin_total_materials)}")
    if kelvin_total_materials and not core_total_materials:
        raise ValueError(
            "total_materials names only Kelvin exterior materials; the "
            "source/total interface then has no total-potential side")

    reduced_selector = mesh.Materials("|".join(reduced_materials))
    total_selector = mesh.Materials("|".join(total_materials))
    core_total_selector = (
        total_selector if not kelvin_total_materials
        else mesh.Materials("|".join(core_total_materials)))
    kelvin_selector = (
        None if not kelvin_total_materials
        else mesh.Materials("|".join(kelvin_total_materials)))
    # The source enclosure and the Kelvin exterior share one periodic space.
    # The magnetic scalar potential is a TWISTED 0-form, so the orientation
    # reversing Kelvin inversion pulls it back as ``Omega_comp = -Omega_phys``
    # with no metric factor (docs/kelvin/KELVIN_TRANSFORMATION.md 2.3).  Solving
    # for ``-Omega_comp`` in the exterior therefore leaves one unknown that is
    # continuous across the identified spheres up to the source jump, which a
    # periodic H1 space represents directly.  Restricting the total space to
    # the total materials instead leaves it without a single degree of freedom
    # on the physical sphere, because the material touching that sphere is the
    # reduced one; the identification then has nothing to pair and the whole
    # exterior block is silently eliminated.
    coupled_selector = (
        reduced_selector if kelvin_selector is None
        else mesh.Materials("|".join(tuple(reduced_materials)
                                     + kelvin_total_materials)))
    interface_selector = mesh.Boundaries(interface_boundary)
    kelvin_interface_selector = (
        None if kelvin_interface_boundary is None
        else mesh.Boundaries(kelvin_interface_boundary))

    # Surface values belong to the named potential, not automatically to H.
    if surface_dirichlet is not None and not isinstance(surface_dirichlet, Mapping):
        raise ValueError("surface_dirichlet must be a mapping")
    surface_values = {} if surface_dirichlet is None else dict(surface_dirichlet)
    if set(surface_values) - {"reduced", "total"}:
        raise ValueError("surface_dirichlet keys must be reduced or total")
    surface_labels = {}
    for region, materials in (("reduced", reduced_set), ("total", set(core_total_materials))):
        values = surface_values.get(region, {})
        if not isinstance(values, Mapping):
            raise ValueError("surface_dirichlet region values must be boundary/value mappings")
        labels = set(values)
        if not labels <= set(mesh.GetBoundaries()) or any(not label or "|" in label for label in labels):
            raise ValueError("surface_dirichlet must name exact mesh boundary labels")
        for face in mesh.ngmesh.FaceDescriptors():
            if mesh.GetBoundaries()[face.bc - 1] in labels:
                adjacent = [d for d in (face.domin, face.domout) if d != 0]
                if len(adjacent) != 1 or mesh.GetMaterials()[adjacent[0] - 1] not in materials:
                    raise ValueError("surface Dirichlet must be exterior to its named potential region")
        surface_labels[region] = "|".join(sorted(labels))
    has_surface = any(surface_labels.values())
    if has_surface:
        if kelvin_selector is not None:
            raise ValueError("surface Dirichlet with Kelvin exterior is not supported")
        if total_dirichlet_cf is not None:
            raise ValueError("surface Dirichlet fixes the gauge; do not also supply total_dirichlet_cf")
        if set(surface_values.get("reduced", {})) & set(str(reduced_zero_normal_boundary).split("|")):
            raise ValueError("Dirichlet and zero-normal conditions overlap")
        if all(surface_labels.values()):
            from ngsolve import BND
            nodes = {
                region: {vertex.nr for face in mesh.Elements(BND)
                         if str(face.mat) in surface_values.get(region, {})
                         for vertex in face.vertices}
                for region in ("reduced", "total")}
            if nodes["reduced"] & nodes["total"]:
                raise ValueError("surface Dirichlet on both sides of an interface junction is not supported")

    if kelvin_selector is None:
        def _dirichlet_labels(*parts):
            return "|".join(part for part in parts if part)

        fes_reduced = H1(mesh, order=int(order), definedon=coupled_selector,
                         dirichlet=_dirichlet_labels(surface_labels["reduced"],
                                                     reduced_dirichlet_boundary))
        fes_total = Compress(Periodic(H1(
            mesh, order=int(order), definedon=core_total_selector,
            dirichlet=_dirichlet_labels(surface_labels["total"],
                                        total_dirichlet_boundary),
            dirichlet_bbbnd=("" if has_surface else (dirichlet_bbbnd or "")))))
    else:
        # Exactly one of the two blocks carries the point gauge.  Pinning both
        # would fight the interface jump; pinning neither leaves the coupled
        # air/exterior block floating.
        coupled_gauged = _constrains_point_gauge(
            mesh, coupled_selector, dirichlet_bbbnd)
        total_gauged = _constrains_point_gauge(
            mesh, core_total_selector, dirichlet_bbbnd)
        if coupled_gauged and total_gauged:
            raise ValueError(
                f"dirichlet_bbbnd={dirichlet_bbbnd!r} constrains both the "
                "source/Kelvin block and the total block; the interface jump "
                "already ties their constants, so place the gauge inside one "
                "of them")
        if not (coupled_gauged or total_gauged):
            raise ValueError(
                f"dirichlet_bbbnd={dirichlet_bbbnd!r} constrains neither the "
                "source/Kelvin block nor the total block; the mixed system has "
                "no gauge")
        coupled_gauge = ({"dirichlet_bbbnd": dirichlet_bbbnd}
                         if coupled_gauged else {})
        total_gauge = ({"dirichlet_bbbnd": dirichlet_bbbnd}
                       if total_gauged else {})
        fes_reduced = Periodic(H1(
            mesh, order=int(order), definedon=coupled_selector,
            **coupled_gauge,
            **({"dirichlet": reduced_dirichlet_boundary}
               if reduced_dirichlet_boundary else {})))
        fes_total = Compress(Periodic(H1(
            mesh, order=int(order), definedon=core_total_selector,
            **total_gauge,
            **({"dirichlet": total_dirichlet_boundary}
               if total_dirichlet_boundary else {}))))
    fes_multiplier = Compress(H1(
        mesh, order=int(order), definedon=interface_selector,
        **({"dirichlet": interface_multiplier_dirichlet_boundary}
           if interface_multiplier_dirichlet_boundary else {})))
    fes = fes_reduced * fes_total * fes_multiplier
    source_values = None
    if source_rhs_reduced is not None:
        source_values = np.asarray(source_rhs_reduced, dtype=float)
        if (source_values.shape != (fes_reduced.ndof,)
                or not np.isfinite(source_values).all()):
            raise ValueError("source_rhs_reduced must be a finite vector in "
                             "fes_reduced DOF order")
    (phi_reduced, phi_total, multiplier), (
        test_reduced, test_total, test_multiplier) = fes.TnT()

    # The exterior needs the source 0-form that the identified spheres jump by.
    # Two explicit representations are supported and the caller picks one.
    #
    # ``kelvin_source_h`` is the exact pulled-back exterior source field, for
    # example ``radia.KelvinRadiaFieldStrength``.  It equals ``grad`` of the
    # exact lift, so the exterior unknown becomes the REDUCED potential there
    # and the discrete exterior field carries the analytic source itself.  No
    # interface trace is projected at all.
    #
    # Otherwise the projected trace ``kelvin_source_potential`` is lifted into
    # the shared periodic space.  Its periodic degrees of freedom are shared,
    # so setting it on the physical sphere gives the Kelvin sphere the same
    # values, and any discrete extension into the ball yields the same exterior
    # potential.  That representation inherits the projection residual of the
    # trace, which the exact-source form does not.
    kelvin_lift = None
    if kelvin_selector is not None and kelvin_source_h is None:
        kelvin_lift = GridFunction(fes_reduced, name="kelvin_source_lift")
        kelvin_lift.vec[:] = 0.0
        kelvin_lift.Set(kelvin_source_potential,
                        definedon=kelvin_interface_selector)

    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += mu_cf * grad(phi_reduced) * grad(test_reduced) * dx(
        definedon=coupled_selector, bonus_intorder=bonus_intorder)
    a_bf += mu_cf * grad(phi_total) * grad(test_total) * dx(
        definedon=core_total_selector, bonus_intorder=bonus_intorder)
    d_interface = ds(definedon=interface_selector, bonus_intorder=bonus_intorder)
    jump_trial = phi_total.Trace() - phi_reduced.Trace()
    jump_test = test_total.Trace() - test_reduced.Trace()
    a_bf += interface_constraint_scale * (
        multiplier * jump_test + test_multiplier * jump_trial) * d_interface

    fixed_rhs = LinearForm(fes)
    source_load_diagnostics = {}
    reduced_load_h = H_s if load_source_h is None else load_source_h
    if source_rhs_reduced is None and reduced_source_load == "volume":
        fixed_rhs += mu_cf * reduced_load_h * grad(test_reduced) * dx(
            definedon=reduced_selector, bonus_intorder=bonus_intorder)
    elif source_rhs_reduced is None and _fixed_rhs_cache is not None and "fixed" in _fixed_rhs_cache:
        # A Picard cache already holds this assembled load; building the form
        # again would re-evaluate the source for its flux balance only.
        source_load_diagnostics["reduced"] = _fixed_rhs_cache.get("reduced_flux_balance")
    elif source_rhs_reduced is None:
        mu_reduced = _uniform_permeability(
            mesh, mu_cf, reduced_materials, mu_r_by_material, role="reduced")
        flux_form, source_load_diagnostics["reduced"] = outward_source_flux_form(
            mesh, reduced_load_h, material_set_boundary_signs(mesh, reduced_materials),
            test_reduced.Trace(), bonus_intorder=bonus_intorder, scale=mu_reduced)
        fixed_rhs += flux_form
        if _fixed_rhs_cache is not None:
            _fixed_rhs_cache["reduced_flux_balance"] = source_load_diagnostics["reduced"]
    if reduced_zero_normal_boundary is not None:
        from ngsolve import specialcf, BoundaryFromVolumeCF
        labels = set(str(reduced_zero_normal_boundary).split("|"))
        if not labels or not labels <= set(mesh.GetBoundaries()):
            raise ValueError("reduced_zero_normal_boundary must name mesh boundaries")
        for face in mesh.ngmesh.FaceDescriptors():
            if mesh.GetBoundaries()[face.bc - 1] not in labels:
                continue
            domains = (face.domin, face.domout)
            adjacent = [domain for domain in domains if domain != 0]
            if len(adjacent) != 1 or mesh.GetMaterials()[adjacent[0] - 1] not in reduced_set:
                raise ValueError("zero-normal correction boundary must be exterior to a reduced region")
        # Zero normal correction field does not mean zero normal total field.
        fixed_rhs += -BoundaryFromVolumeCF(mu_cf) * InnerProduct(H_s, specialcf.normal(3)) * test_reduced.Trace() * ds(
            definedon=mesh.Boundaries(reduced_zero_normal_boundary),
            bonus_intorder=bonus_intorder)
    if kelvin_lift is not None:
        # Exterior equation for Omega_t = phi_reduced + lift, tested with the
        # same continuous space; the lift moves to the right-hand side.
        fixed_rhs += -mu_cf * grad(kelvin_lift) * grad(test_reduced) * dx(
            definedon=kelvin_selector, bonus_intorder=bonus_intorder)
    elif kelvin_selector is not None:
        # Same right-hand side with the exact source in place of the lift:
        # grad of the exact lift IS the pulled-back exterior source field.
        fixed_rhs += -mu_cf * kelvin_source_h * grad(test_reduced) * dx(
            definedon=kelvin_selector, bonus_intorder=bonus_intorder)
    f_lf = LinearForm(fes)
    if total_source_h is not None:
        total_source_selector = mesh.Materials("|".join(total_source_materials))
        from ngsolve import IntegrationRule, VOL
        # Linear and rectangular bilinear forms infer different default rules.
        # Bind one rule so caching cannot change the discrete excitation.
        material_rhs_measure = dx(
            definedon=total_source_selector,
            intrules={kind: IntegrationRule(kind, order=2 * int(order) + int(bonus_intorder))
                      for kind in {element.type for element in mesh.Elements(VOL)}})
        if total_source_load == "volume":
            f_lf += mu_cf * total_source_h * grad(test_total) * material_rhs_measure
        else:
            # mu (H_s + grad Phi_s) . grad v = mu [(H_s . n_out) v on the
            # boundary + grad Phi_s . grad v], material by material.
            for name in total_source_materials:
                mu_material = _uniform_permeability(
                    mesh, mu_cf, (name,), mu_r_by_material, role="total")
                flux_form, source_load_diagnostics[f"total:{name}"] = outward_source_flux_form(
                    mesh, total_source_h, material_set_boundary_signs(mesh, (name,)),
                    test_total.Trace(), bonus_intorder=bonus_intorder, scale=mu_material)
                f_lf += flux_form
                f_lf += mu_material * InnerProduct(
                    grad(total_source_potential), grad(test_total)) * dx(
                        definedon=mesh.Materials(name), bonus_intorder=bonus_intorder)
    fixed_rhs += interface_constraint_scale * test_multiplier * source_potential * d_interface
    # Prescribed normal flux is part of the source, so it belongs to the cached
    # fixed right-hand side rather than the material-dependent one.
    if reduced_normal_flux is not None:
        fixed_rhs += -reduced_normal_flux * test_reduced * ds(
            definedon=mesh.Boundaries(reduced_flux_boundary),
            bonus_intorder=bonus_intorder)
    if total_normal_flux is not None:
        fixed_rhs += -total_normal_flux * test_total * ds(
            definedon=mesh.Boundaries(total_flux_boundary),
            bonus_intorder=bonus_intorder)
    phase_timings = {}

    def timed_phase(name, action):
        started = time.perf_counter()
        if phase_callback is not None:
            phase_callback({"phase": name, "event": "start"})
        value = action()
        elapsed = time.perf_counter() - started
        phase_timings[name] = elapsed
        if phase_callback is not None:
            phase_callback({"phase": name, "event": "complete", "seconds": elapsed})
        return value

    timed_phase("matrix_assembly", a_bf.Assemble)

    def assemble_rhs():
        nonlocal f_lf
        if _fixed_rhs_cache is None:
            fixed_rhs.Assemble()
            f_lf.Assemble()
            f_lf.vec.data += fixed_rhs.vec
        else:
            # This cache belongs to one Picard call, with fixed source and spaces.
            layout = tuple(space.ndof for space in fes.components)
            if "fixed" not in _fixed_rhs_cache:
                fixed_rhs.Assemble()
                vector = fixed_rhs.vec.CreateVector()
                vector.data = fixed_rhs.vec
                _fixed_rhs_cache.update(fixed=vector, layout=layout)
            elif _fixed_rhs_cache["layout"] != layout:
                raise ValueError("Picard RHS cache space layout changed")
            if _rhs_material is None:
                f_lf.Assemble()
            else:
                if "material_operator" not in _fixed_rhs_cache:
                    operator = BilinearForm(trialspace=_rhs_material.space, testspace=fes)
                    if total_source_h is not None:
                        coefficient = _rhs_material.space.TrialFunction()
                        operator += coefficient * total_source_h * grad(test_total) * material_rhs_measure
                    operator.Assemble()
                    _fixed_rhs_cache["material_operator"] = operator
                f_lf = LinearForm(fes)
                f_lf.Assemble()
                f_lf.vec.data = _fixed_rhs_cache["material_operator"].mat * _rhs_material.vec
            f_lf.vec.data += _fixed_rhs_cache["fixed"]
        if source_values is not None:
            f_lf.vec.FV().NumPy()[:fes_reduced.ndof] += source_values

    timed_phase("source_rhs_assembly", assemble_rhs)

    solution = GridFunction(fes)
    inverse_mat = timed_phase(
        "factorization", lambda: a_bf.mat.Inverse(fes.FreeDofs(), inverse=inverse))
    if total_dirichlet_cf is None and not has_surface:
        def apply_inverse():
            solution.vec.data = inverse_mat * f_lf.vec
        timed_phase("backsolve", apply_inverse)
    else:
        if has_surface:
            for index, region in enumerate(("reduced", "total")):
                values = surface_values.get(region, {})
                if values:
                    solution.components[index].Set(
                        mesh.BoundaryCF(values, default=0),
                        definedon=mesh.Boundaries(surface_labels[region]))
        else:
            solution.components[1].Set(
                total_dirichlet_cf, definedon=mesh.BBBoundaries(dirichlet_bbbnd))
        residual = solution.vec.CreateVector()
        residual.data = f_lf.vec - a_bf.mat * solution.vec
        def apply_inverse():
            solution.vec.data += inverse_mat * residual
        timed_phase("backsolve", apply_inverse)

    # A direct solve has no iteration history, but it still has a residual.
    # r = b - A x on the system that was ACTUALLY solved -- non-zero Dirichlet
    # lift included, because `solution` already carries it.  Report the free
    # dofs separately from the constrained ones: a constrained row holds the
    # Dirichlet equation, not the system, so mixing them in hides the number
    # that matters.  The multiplier block is split out on its own because
    # those rows ARE the interface jump condition in weak form,
    # ``scale * int test_multiplier (phi_total - phi_reduced - Phi_s) ds``.
    linear_residual = _report_direct_solve_residual(
        a_bf, f_lf, solution, fes,
        block_names=("phi_reduced", "phi_total", "interface_constraint"))
    assembled_energy = _assembled_primal_energy(
        a_bf, f_lf, solution, fes, primal_blocks=(0, 1))

    phi_reduced_gf, phi_total_gf, multiplier_gf = solution.components[:3]
    H_reduced = H_s - grad(phi_reduced_gf)
    if total_source_load == "surface_flux":
        # The load used the raw source and its Hodge potential separately;
        # the field keeps the harmonic-remainder meaning of total_source_h.
        total_source_h = total_source_h + grad(total_source_potential)
    zero_h = CoefficientFunction((0.0, 0.0, 0.0))
    total_source_by_material = mesh.MaterialCF({
        material: total_source_h
        if total_source_h is not None and material in total_source_materials
        else zero_h
        for material in total_materials
    })
    H_total = total_source_by_material - grad(phi_total_gf)
    # Twisted 0-form: the stored exterior unknown is -Omega_comp, so the
    # computational-frame H is +grad of it rather than -grad.
    kelvin_exterior_source = (
        grad(kelvin_lift) if kelvin_lift is not None else kelvin_source_h)
    kelvin_potential_cf = (
        None if kelvin_lift is None else -(phi_reduced_gf + kelvin_lift))
    H_kelvin = (
        None if kelvin_selector is None
        else grad(phi_reduced_gf) + kelvin_exterior_source)
    kelvin_set = set(kelvin_total_materials)
    h_components = []
    for component in range(3):
        values = {}
        for material in mesh.GetMaterials():
            if material in reduced_set:
                values[material] = H_reduced[component]
            elif material in kelvin_set:
                values[material] = H_kelvin[component]
            else:
                values[material] = H_total[component]
        h_components.append(mesh.MaterialCF(values))
    H_cf = CoefficientFunction(tuple(h_components))
    return {
        "solution": solution,
        "linear_residual": linear_residual,
        "phase_timings_seconds": phase_timings,
        "assembled_energy": assembled_energy,
        # The assembled system, ONLY on request.  A caller testing whether
        # another order's solution is admissible HERE needs it: with the
        # multiplier entries set to zero, the multiplier rows of f - A x are the
        # interface constraint violation g - B x_P.  Returning it by default
        # would keep the assembled matrix alive in every production result the
        # wrappers hand back, where it used to be freed after the solve.
        "system": ({"bilinear_form": a_bf, "linear_form": f_lf}
                   if return_system else None),
        "phi_reduced": phi_reduced_gf,
        "phi_total": phi_total_gf,
        "interface_multiplier": multiplier_gf,
        "kelvin_source_lift": kelvin_lift,
        "kelvin_exterior_source": kelvin_exterior_source,
        "kelvin_total_potential": kelvin_potential_cf,
        "kelvin_materials": kelvin_total_materials,
        "fes": fes,
        "fes_reduced": fes_reduced,
        "fes_total": fes_total,
        "mu_cf": mu_cf,
        "H_s": H_s,
        "source_potential": source_potential,
        "kelvin_source_potential": kelvin_source_potential,
        "total_source_h": total_source_h,
        "total_source_materials": total_source_materials,
        "source_loads": {"reduced": reduced_source_load, "total": total_source_load,
                         "surface_flux_balance": source_load_diagnostics},
        "H_cf": H_cf,
        "B_cf": mu_cf * H_cf,
    }


def _map_matching_h1_load_to_global(
        mesh, source_space, target_space, source_values, *, order):
    """Map an uncompressed P1/P2 H1 load by shared mesh entities.

    The helper is deliberately narrower than a generic prolongation.  It is
    the numbering bridge needed when a reduced-region surface-source load was
    assembled before the matching-trace total/reduced spaces were condensed
    into one global H1 space.
    """
    import ngsolve as ng

    order = int(order)
    if order not in (1, 2):
        raise ValueError("matching H1 load mapping supports only order 1 or 2")
    values = np.asarray(source_values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError(
            "source_values must be a finite H1 load vector")

    # A producer may deliberately assemble directly in the global H1
    # numbering.  This is the preferred high-throughput contract.
    if values.shape == (target_space.ndof,):
        return values.copy()
    if values.shape != (source_space.ndof,):
        raise ValueError(
            "source_values must use source-space or global H1 DOF order")

    if source_space.ndof > target_space.ndof:
        raise ValueError("source H1 space is larger than the global H1 space")

    # A defined-on source space is compressed.  Keep the generic entity map
    # for compatibility; performance-sensitive producers should return the
    # global-numbered vector handled above.
    mapped = np.zeros(target_space.ndof, dtype=float)
    node_groups = [(ng.VERTEX, mesh.vertices)]
    if order == 2:
        node_groups.append((ng.EDGE, mesh.edges))
    seen_source = set()
    for node_type, entities in node_groups:
        for entity in entities:
            node = ng.NodeId(node_type, entity.nr)
            source_dofs = tuple(d for d in source_space.GetDofNrs(node) if d >= 0)
            if not source_dofs:
                continue
            target_dofs = tuple(d for d in target_space.GetDofNrs(node) if d >= 0)
            if len(source_dofs) != 1 or len(target_dofs) != 1:
                raise ValueError(
                    "matching H1 load mapping requires one uncompressed DOF "
                    "per active vertex/edge")
            source_dof, target_dof = source_dofs[0], target_dofs[0]
            seen_source.add(source_dof)
            mapped[target_dof] += values[source_dof]
    if not set(np.flatnonzero(np.abs(values) > 1.0e-30)).issubset(seen_source):
        raise ValueError("source load contains unmapped active H1 DOFs")
    return mapped


def _copy_matching_h1_trace_to_global(
        mesh, source_trace, target, interface_boundary, *, order):
    """Copy projected trace coefficients without re-projecting through Set."""
    import ngsolve as ng

    vertex_ids = set()
    edge_ids = set()
    for element in mesh.Elements(ng.BND):
        if element.mat != interface_boundary:
            continue
        vertex_ids.update(vertex.nr for vertex in element.vertices)
        if int(order) == 2:
            edge_ids.update(edge.nr for edge in element.edges)
    groups = [(ng.VERTEX, sorted(vertex_ids))]
    if int(order) == 2:
        groups.append((ng.EDGE, sorted(edge_ids)))
    target_values = target.vec.FV().NumPy()
    source_values = source_trace.vec.FV().NumPy()
    for node_type, entity_ids in groups:
        for entity_id in entity_ids:
            source_dofs = tuple(
                dof for dof in source_trace.space.GetDofNrs(
                    ng.NodeId(node_type, entity_id)) if dof >= 0)
            target_dofs = tuple(
                dof for dof in target.space.GetDofNrs(
                    ng.NodeId(node_type, entity_id)) if dof >= 0)
            if len(source_dofs) != 1 or len(target_dofs) != 1:
                raise ValueError(
                    "matching trace condensation requires one projected "
                    "P1/P2 DOF per interface vertex/edge")
            target_values[target_dofs[0]] = source_values[source_dofs[0]]


def _zero_h1_boundary_dofs(mesh, grid_function, boundary, *, order):
    """Zero P1/P2 boundary entities, including interface/Dirichlet junctions."""
    import ngsolve as ng

    boundary_names = set(str(boundary).split("|"))
    vertex_ids = set()
    edge_ids = set()
    for element in mesh.Elements(ng.BND):
        if element.mat not in boundary_names:
            continue
        vertex_ids.update(vertex.nr for vertex in element.vertices)
        if int(order) == 2:
            edge_ids.update(edge.nr for edge in element.edges)
    values = grid_function.vec.FV().NumPy()
    for node_type, entity_ids in ((ng.VERTEX, vertex_ids), (ng.EDGE, edge_ids)):
        if node_type == ng.EDGE and int(order) != 2:
            continue
        for entity_id in entity_ids:
            for dof in grid_function.space.GetDofNrs(
                    ng.NodeId(node_type, entity_id)):
                if dof >= 0:
                    values[dof] = 0.0


class LinearSolveNotConverged(RuntimeError):
    """The linear solve returned, but not with a solution.

    NGSolve's CGSolver stops at ``maxiter`` without raising, and with
    ``printrates=False`` without a word; a direct factorisation can return
    garbage on a singular or ill-scaled system just as quietly.  The residual
    of what came back is the only witness, so it is checked before the result
    is handed on, and this carries it.
    """

    def __init__(self, message, residual):
        super().__init__(message)
        self.residual = residual


def boundaries_touching_materials(mesh, names, materials):
    """Which of ``names`` has a face on an element of one of ``materials``?

    The mesh already knows this. A boundary region's ``Neighbours(VOL)`` is a
    region, and intersecting its material mask with the materials' mask is the
    question, in two calls and no Python iteration.

    Rebuilding the same adjacency by hand -- walking every volume element, then
    every face of every element, into a dictionary of sets -- is the obvious
    thing to write and costs 2.05 s on a 274k-element mesh, about a quarter of
    the solve it belongs to. The region form takes 0.012 s and returns the same
    set; ``tests/test_boundary_material_adjacency.py`` holds them to that.

    Args:
        mesh: an NGSolve mesh.
        names: boundary names to test. Names the mesh does not carry are
            ignored rather than raising, matching what a scan over the mesh's
            own boundary elements would have found.
        materials: material names whose elements count as touching.

    Returns:
        The subset of ``names`` that touches, as a set.
    """
    import ngsolve as ng

    wanted = set(str(name) for name in names)
    available = set(mesh.GetMaterials())
    selected = [str(m) for m in materials if str(m) in available]
    if not selected:
        return set()
    material_mask = mesh.Materials("|".join(sorted(selected))).Mask()
    touching = set()
    for name in sorted(wanted & set(mesh.GetBoundaries())):
        neighbours = mesh.Boundaries(name).Neighbours(ng.VOL).Mask()
        if (neighbours & material_mask).NumSet():
            touching.add(name)
    return touching


def solve_magnetostatic_matching_trace_total_reduced_omega(
        mesh, H_s, source_trace, *, mu_r_by_material,
        reduced_materials, total_materials, interface_boundary,
        order=1, dirichlet_boundary=None, bonus_intorder=4,
        source_rhs_reduced=None, reduced_normal_flux=None,
        reduced_flux_boundary=None, solver="direct", inverse="pardiso",
        cg_preconditioner="local", cg_tolerance=1.0e-10,
        cg_max_iterations=2000, return_system=False,
        dirichlet_bbbnd=None, total_source_h=None, total_source_materials=()):
    """Solve the finite-domain matching-trace mixed Omega problem as SPD.

    This is a strict acceleration path for a joined, non-periodic mesh whose
    total and reduced potentials have the same P1/P2 trace space.  With a
    projected interface lift ``L`` the constrained variables are
    ``phi_total = u`` and ``phi_reduced = u - L``.  The Lagrange multiplier is
    therefore eliminated exactly and the remaining global H1 system is SPD,
    allowing CG.  Kelvin exteriors, periodic/compressed spaces, nonmatching
    traces, nonlinear permeability and orders above two stay on the general
    saddle-point driver.

    ``source_trace`` must be an NGSolve GridFunction.  Accepting an arbitrary
    CoefficientFunction would silently change the projection used by the
    multiplier formulation and has been observed to change P1 fields.
    """
    import ngsolve as ng

    function_started = time.perf_counter()

    order = int(order)
    if order not in (1, 2):
        raise ValueError("matching-trace condensation supports only order 1 or 2")
    if not isinstance(source_trace, GridFunction):
        raise TypeError(
            "source_trace must be a projected NGSolve GridFunction, not an "
            "arbitrary CoefficientFunction")
    reduced_materials = tuple(str(name) for name in reduced_materials)
    total_materials = tuple(str(name) for name in total_materials)
    actual = set(mesh.GetMaterials())
    reduced_set = set(reduced_materials)
    total_set = set(total_materials)
    total_source_materials = tuple(total_source_materials)
    if not set(total_source_materials) <= total_set:
        raise ValueError("total_source_materials must be in total_materials")
    if bool(total_source_materials) != (total_source_h is not None):
        raise ValueError("total_source_h and total_source_materials must be supplied together")
    if (not reduced_set or not total_set or reduced_set & total_set
            or reduced_set | total_set != actual):
        raise ValueError(
            "reduced_materials and total_materials must be a disjoint, "
            "exhaustive mesh partition")
    if interface_boundary not in mesh.GetBoundaries():
        raise ValueError("interface_boundary must name an existing boundary")
    if (reduced_normal_flux is None) != (reduced_flux_boundary is None):
        raise ValueError(
            "reduced_normal_flux and reduced_flux_boundary must be supplied together")
    if solver not in ("direct", "cg"):
        raise ValueError("solver must be 'direct' or 'cg'")

    material_mu = {}
    supplied_mu = dict(mu_r_by_material or {})
    for material in mesh.GetMaterials():
        mu_r = float(supplied_mu.get(material, 1.0))
        if not math.isfinite(mu_r) or mu_r <= 0.0:
            raise ValueError(f"mu_r for {material!r} must be positive and finite")
        material_mu[material] = MU_0 * mu_r
    mu_cf = mesh.MaterialCF(material_mu)
    reduced_selector = mesh.Materials("|".join(reduced_materials))
    interface_selector = mesh.Boundaries(interface_boundary)
    fes = H1(
        mesh, order=order,
        **({"dirichlet_bbbnd": dirichlet_bbbnd} if dirichlet_bbbnd else {}),
        **({"dirichlet": dirichlet_boundary} if dirichlet_boundary else {}))
    if all(fes.FreeDofs()):
        raise ValueError(
            "matching-trace Omega requires a gauge: dirichlet_boundary or "
            "dirichlet_bbbnd must constrain at least one scalar degree of freedom")
    reduced_space = H1(mesh, order=order, definedon=reduced_selector)
    u, v = fes.TnT()

    lift = GridFunction(fes, name="source_interface_lift")
    lift.vec[:] = 0.0
    _copy_matching_h1_trace_to_global(
        mesh, source_trace, lift, interface_boundary, order=order)
    if dirichlet_boundary:
        # The multiplier formulation omits the jump equation where both
        # scalar potentials are already fixed.  Zeroing the lift on that
        # junction is its exact eliminated counterpart.
        _zero_h1_boundary_dofs(
            mesh, lift, dirichlet_boundary, order=order)
    if dirichlet_bbbnd:
        lift.vec.FV().NumPy()[~np.asarray(list(fes.FreeDofs()), dtype=bool)] = 0.0

    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += mu_cf * grad(u) * grad(v) * dx(
        bonus_intorder=bonus_intorder)
    preconditioner = None
    if solver == "cg":
        from ngsolve import Preconditioner
        preconditioner = Preconditioner(a_bf, cg_preconditioner)
    f_lf = LinearForm(fes)
    if source_rhs_reduced is None:
        f_lf += mu_cf * H_s * grad(v) * dx(
            definedon=reduced_selector, bonus_intorder=bonus_intorder)
    f_lf += mu_cf * grad(lift) * grad(v) * dx(
        definedon=reduced_selector, bonus_intorder=bonus_intorder)
    if total_source_h is not None:
        f_lf += mu_cf * total_source_h * grad(v) * dx(
            definedon=mesh.Materials("|".join(total_source_materials)),
            bonus_intorder=bonus_intorder)
    if reduced_normal_flux is not None:
        reduced_flux_names = boundaries_touching_materials(
            mesh, str(reduced_flux_boundary).split("|"), reduced_materials)
        if not reduced_flux_names:
            raise ValueError(
                "reduced_flux_boundary has no face adjacent to reduced materials")
        f_lf += -reduced_normal_flux * v * ds(
            definedon=mesh.Boundaries("|".join(sorted(reduced_flux_names))),
            bonus_intorder=bonus_intorder)

    setup_before_assembly = time.perf_counter() - function_started

    started = time.perf_counter()
    a_bf.Assemble()
    matrix_assembly = time.perf_counter() - started
    started = time.perf_counter()
    f_lf.Assemble()
    if source_rhs_reduced is not None:
        mapped = _map_matching_h1_load_to_global(
            mesh, reduced_space, fes, source_rhs_reduced, order=order)
        f_lf.vec.FV().NumPy()[:] += mapped
    rhs_assembly = time.perf_counter() - started

    solution = GridFunction(fes, name="phi_total_condensed")
    residual = solution.vec.CreateVector()
    residual.data = f_lf.vec - a_bf.mat * solution.vec
    started = time.perf_counter()
    iterations = None
    if solver == "direct":
        inv = a_bf.mat.Inverse(fes.FreeDofs(), inverse=inverse)
        setup_seconds = time.perf_counter() - started
        started = time.perf_counter()
        solution.vec.data += inv * residual
        solve_seconds = time.perf_counter() - started
    else:
        preconditioner.Update()
        setup_seconds = time.perf_counter() - started
        from ngsolve.krylovspace import CGSolver
        started = time.perf_counter()
        inv = CGSolver(
            mat=a_bf.mat, pre=preconditioner.mat,
            tol=float(cg_tolerance), maxiter=int(cg_max_iterations),
            printrates=False)
        solution.vec.data += inv * residual
        solve_seconds = time.perf_counter() - started
        iterations = getattr(inv, "iterations", None)

    postprocess_started = time.perf_counter()
    phi_reduced = solution - lift
    H_reduced = H_s - grad(solution) + grad(lift)
    zero = CoefficientFunction((0.0, 0.0, 0.0))
    total_source = mesh.MaterialCF({name: total_source_h for name in total_source_materials}, default=zero)
    H_total = total_source - grad(solution)
    components = []
    for component in range(3):
        components.append(mesh.MaterialCF({
            material: (H_reduced[component] if material in reduced_set
                       else H_total[component] if material in total_set
                       else zero[component])
            for material in mesh.GetMaterials()
        }))
    H_cf = CoefficientFunction(tuple(components))
    final_residual = f_lf.vec.CreateVector()
    final_residual.data = f_lf.vec - a_bf.mat * solution.vec
    free = np.asarray(list(fes.FreeDofs()), dtype=bool)
    residual_values = final_residual.FV().NumPy()[free]
    rhs_values = f_lf.vec.FV().NumPy()[free]
    relative_residual = float(
        np.linalg.norm(residual_values) / max(np.linalg.norm(rhs_values), 1.0e-30))
    linear_residual = {
        "free_dofs": {
            "l2": float(np.linalg.norm(residual_values)),
            "rhs_l2": float(np.linalg.norm(rhs_values)),
            "relative": relative_residual,
        }
    }
    # A solve that ran out of iterations, or came back non-finite, used to be
    # returned like any other; CGSolver does not raise at maxiter and says
    # nothing with printrates=False.  The true residual is checked here, on
    # the solution the caller will receive.  CG's own tolerance is on the
    # preconditioned residual, so the true one is allowed a factor over it.
    residual_limit = (100.0 * float(cg_tolerance) if solver == "cg" else 1.0e-8)
    linear_residual["accepted_below"] = residual_limit
    if not (math.isfinite(relative_residual) and relative_residual <= residual_limit):
        raise LinearSolveNotConverged(
            "matching-trace linear solve did not converge: "
            f"solver={solver}, iterations={locals().get('iterations')}, "
            f"relative_residual={relative_residual:.3e}, "
            f"accepted_below={residual_limit:.3e}",
            linear_residual)
    postprocess_seconds = time.perf_counter() - postprocess_started
    return {
        "solution": solution,
        "phi_total": solution,
        "phi_reduced": phi_reduced,
        "source_lift": lift,
        "fes": fes,
        "fes_reduced_source": reduced_space,
        "mu_cf": mu_cf,
        "H_cf": H_cf,
        "B_cf": mu_cf * H_cf,
        "linear_residual": linear_residual,
        "linear_residual_relative": relative_residual,
        "solver": solver,
        "iterations": iterations,
        "phase_timings_seconds": {
            "setup_before_assembly": setup_before_assembly,
            "matrix_assembly": matrix_assembly,
            "source_rhs_assembly": rhs_assembly,
            "preconditioner_or_factorization": setup_seconds,
            "backsolve": solve_seconds,
            "solver_postprocess": postprocess_seconds,
        },
        "system": ({"bilinear_form": a_bf, "linear_form": f_lf}
                   if return_system else None),
    }


class MixedOmegaPicardNotConverged(RuntimeError):
    """Fail-loud non-convergence of the mixed total/reduced Omega Picard loop.

    ``state`` carries the per-element permeability (``mu_r_elements`` in the order
    of ``element_numbers``) and the full ``nonlinear_stats`` so a caller can persist
    them as an explicitly partial state and warm-start a later run through
    ``mu_r_initial``.  No field result leaves a non-converged loop.
    """

    def __init__(self, message, state):
        super().__init__(message)
        self.state = state


def solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
        mesh, H_s, source_potential, R_K, offset, *, bh_table,
        nonlinear_materials, reduced_materials, total_materials,
        interface_boundary, order=1, dirichlet_bbbnd="GND",
        bonus_intorder=4, kelvin_mats=("kelvin",), inverse="pardiso",
        mu_r_initial=1000.0, tolerance=2.0e-5, max_iterations=80,
        relaxation=0.3, interface_constraint_scale=None,
        kelvin_interface_boundary=None, kelvin_source_potential=None,
        kelvin_source_h=None,
        total_source_h=None, total_source_materials=(),
        anderson_depth=0, anderson_transform="log", observation_points=None,
        material_update_order=None, material_log_state_initial=None,
        refinement_parent_identity=None, cache_fixed_rhs=True,
        reduced_zero_normal_boundary=None, surface_dirichlet=None,
        kelvin_match_exact=False, mu_r_by_material=None,
        material_sampling="element_centroid", bh_interpolation="pchip",
        progress_callback=None, reduced_source_load="volume"):
    """Picard solve for the mixed total/reduced Omega formulation.

    The source split and its interface trace stay fixed throughout the
    iteration.  Only the permeability in ``nonlinear_materials`` is updated
    from the shared monotone ``B(H)`` law.  This is deliberately separate from
    hysteretic state evolution: the latter needs its own committed material
    history and is not silently approximated by a memoryless Picard update.
    Linear physical materials in the total region require explicit entries in
    ``mu_r_by_material``. The B-H table must start at the soft-magnetic origin.
    P1 convergence checks both flux iteration change and secant consistency;
    returned H, mu and B belong to the same assembled iterate.
    ``constitutive_field_audit`` separately measures their spatial B(H) defect
    at quadrature points. Iterative convergence does not imply that this defect
    is small, nor that the returned field has passed a physical accuracy gate.

    ``mu_r_initial`` is one scalar or one value per nonlinear element in mesh
    element order (``nonlinear_stats["element_numbers"]`` of an earlier solve),
    which is the warm start of a resumed iteration.  ``anderson_depth`` > 0 mixes
    the per-element permeability with the constrained Anderson accelerator of
    :mod:`radia.picard_acceleration` (real-valued, projected onto the secant range
    of the B(H) law, extrapolated in ``anderson_transform`` space -- ``"log"`` by
    default because a shielded yoke spans decades of secant permeability -- and
    reverted whenever an accelerated iterate worsens the residual); depth 0 is the
    plain damped update.  ``observation_points`` (N, 3) records the iteration-to-
    iteration change of B at those points, so the stopping criterion can be judged
    where the result is consumed instead of only through the material step.

    The default material update is one order-0 permeability value sampled at
    each nonlinear element centroid and therefore supports ``order=1`` only.
    A genuinely nonlinear ``order=2`` solve must explicitly select
    ``material_update_order=1``. That path projects the logarithm of the B-H
    secant permeability into discontinuous L2 and exponentiates a bounded field,
    so the physical permeability remains positive while its spatial order is
    compatible with the response. Its restart input is the complete
    ``nonlinear_stats["material_restart_state"]`` mapping. Bare coefficient
    arrays are rejected because their mesh, material, order, and B-H identity
    cannot be verified.

    ``refinement_parent_identity`` is an optional non-empty mapping that names
    the common physical geometry, material, and excitation contract shared by
    an h-refinement family. Its canonical digest is attached to field and
    energy observables, allowing a convergence gate to reject a stale or
    unrelated mesh level before comparing numbers.

    The result has the same field keys as the linear mixed solve plus
    ``nonlinear_stats`` (with the per-iteration ``history``, a contraction-rate
    estimate, and its restart state). The projected high-order path also returns
    energy and coenergy evaluated from the final solution and the exact PCHIP
    H-potential under a shared identity. A loop that reaches
    ``max_iterations`` raises :class:`MixedOmegaPicardNotConverged` carrying that
    state.  Caller wraps the complete operation in :class:`ngsolve.TaskManager`.

    ``reduced_source_load="surface_flux"`` moves the fixed reduced load to the
    region's boundary as in the linear solve.  It is valid for B(H) iron
    because the reduced region's permeability stays constant; the iron volume
    term, whose permeability is updated, has no such form here.
    """
    from ngsolve import GridFunction, L2, VOL
    from radia.picard_acceleration import (
        ConstrainedAndersonAccelerator, estimate_contraction_rate)
    from radia.scalar_potential_solver import _build_bh_interpolator

    _check_source_load("reduced_source_load", reduced_source_load)
    nonlinear_materials = tuple(nonlinear_materials)
    nonlinear_set = set(nonlinear_materials)
    actual_materials = set(mesh.GetMaterials())
    if not nonlinear_set:
        raise ValueError("nonlinear_materials must name at least one material")
    if not nonlinear_set <= set(total_materials):
        raise ValueError("nonlinear_materials must be contained in total_materials")
    if nonlinear_set - actual_materials:
        raise ValueError(
            f"nonlinear_materials are not mesh materials: {sorted(nonlinear_set - actual_materials)}")
    kelvin_mats = tuple(kelvin_mats)
    kelvin_set = {name for name in actual_materials
                  if _is_kelvin_material(name, kelvin_mats, exact=kelvin_match_exact)}
    if kelvin_set & nonlinear_set:
        raise ValueError("nonlinear_materials must not include Kelvin materials")
    mu_r_by_material = dict(mu_r_by_material or {})
    if set(mu_r_by_material) & nonlinear_set:
        raise ValueError("mu_r_by_material must not override nonlinear materials")
    missing_linear = set(total_materials) - nonlinear_set - kelvin_set - set(mu_r_by_material)
    if missing_linear:
        raise ValueError(f"declare mu_r_by_material for linear total materials: {sorted(missing_linear)}")
    # Validate names, values and Kelvin overrides before any nonlinear work.
    make_kelvin_mu_cf(mesh, R_K, offset, kelvin_mats=kelvin_mats,
                      kelvin_match_exact=kelvin_match_exact,
                      mu_r_by_material=mu_r_by_material)
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive and finite")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be positive")
    if not 0.0 < relaxation <= 1.0:
        raise ValueError("relaxation must lie in (0, 1]")
    if int(anderson_depth) < 0:
        raise ValueError("anderson_depth must be non-negative")
    if refinement_parent_identity is None:
        normalized_refinement_parent_identity = None
        refinement_parent_identity_sha256 = None
    else:
        if not isinstance(refinement_parent_identity, Mapping):
            raise ValueError("refinement_parent_identity must be a mapping")
        normalized_refinement_parent_identity = dict(refinement_parent_identity)
        if not normalized_refinement_parent_identity:
            raise ValueError("refinement_parent_identity must not be empty")
        try:
            refinement_parent_identity_sha256 = _identity_sha256(
                normalized_refinement_parent_identity
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "refinement_parent_identity must be JSON serializable"
            ) from exc

    bh_array = np.asarray(bh_table, dtype=float)
    if bh_array.ndim != 2 or bh_array.shape[1] < 2:
        raise ValueError("bh_table must contain [H, B] rows")
    if len(bh_array) < 2 or not np.all(np.isfinite(bh_array[:, :2])):
        raise ValueError("bh_table requires at least two finite [H, B] rows")
    if not np.array_equal(bh_array[0, :2], [0.0, 0.0]):
        raise ValueError("memoryless soft-magnetic bh_table must start at [0, 0]")
    positive_rows = (bh_array[:, 0] > 0.0) & (bh_array[:, 1] > 0.0)
    secants = bh_array[positive_rows, 1] / bh_array[positive_rows, 0]
    genuinely_nonlinear = (
        secants.size > 1
        and float(np.ptp(secants))
        > 1.0e-12 * max(float(np.max(np.abs(secants))), MU_0)
    )
    requested_material_order = (
        None if material_update_order is None else int(material_update_order)
    )
    if material_sampling not in ("element_centroid", "integration_point"):
        raise ValueError("material_sampling must be element_centroid or integration_point")
    if bh_interpolation not in ("pchip", "linear_spline"):
        raise ValueError("bh_interpolation must be pchip or linear_spline")
    if material_sampling != "integration_point" and bh_interpolation != "pchip":
        raise ValueError("linear_spline is available only for integration_point sampling")
    if material_sampling == "integration_point":
        if progress_callback is not None and not callable(progress_callback):
            raise ValueError("progress_callback must be callable")
        if int(order) != 1 or requested_material_order not in (None, 0):
            raise ValueError("integration_point material sampling currently requires P1")
        if int(anderson_depth) != 0:
            raise ValueError("integration_point material sampling requires anderson_depth=0")
        if float(relaxation) != 1.0:
            raise ValueError("integration_point material sampling requires relaxation=1")
        if material_log_state_initial is not None:
            raise ValueError("integration_point material sampling has no projected restart state")
        return _solve_mixed_omega_pointwise_picard(
            mesh, H_s, source_potential, R_K, offset,
            bh_array=bh_array, nonlinear_materials=nonlinear_materials,
            reduced_materials=reduced_materials, total_materials=total_materials,
            interface_boundary=interface_boundary, order=int(order),
            dirichlet_bbbnd=dirichlet_bbbnd, bonus_intorder=bonus_intorder,
            kelvin_mats=kelvin_mats, inverse=inverse, mu_r_initial=mu_r_initial,
            tolerance=tolerance, max_iterations=max_iterations,
            interface_constraint_scale=interface_constraint_scale,
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            kelvin_source_h=kelvin_source_h, total_source_h=total_source_h,
            total_source_materials=total_source_materials,
            reduced_zero_normal_boundary=reduced_zero_normal_boundary,
            surface_dirichlet=surface_dirichlet, kelvin_match_exact=kelvin_match_exact,
            mu_r_by_material=mu_r_by_material,
            bh_interpolation=bh_interpolation,
            progress_callback=progress_callback,
            reduced_source_load=reduced_source_load)
    if genuinely_nonlinear and int(order) > 1 and requested_material_order is None:
        raise ValueError(
            "nonlinear mixed total/reduced Omega currently updates permeability "
            "as one order-0 centroid value per element; response order > 1 would "
            "use an unmatched material polynomial order. Use order=1 with an "
            "h-convergence and volume-observable gate, or explicitly select "
            "material_update_order=order-1."
        )
    if requested_material_order is not None and requested_material_order < 0:
        raise ValueError("material_update_order must be nonnegative")
    if genuinely_nonlinear and int(order) > 1:
        if reduced_zero_normal_boundary is not None:
            raise ValueError("reduced zero-normal boundary currently requires order=1")
        if requested_material_order != int(order) - 1:
            raise ValueError(
                "a high-order nonlinear response requires "
                "material_update_order=order-1"
            )
        if int(anderson_depth) != 0:
            raise ValueError(
                "anderson_depth must be 0 for the projected high-order material path"
            )
        return _solve_mixed_omega_projected_log_material(
            mesh,
            H_s,
            source_potential,
            R_K,
            offset,
            bh_array=bh_array,
            nonlinear_materials=nonlinear_materials,
            reduced_materials=reduced_materials,
            total_materials=total_materials,
            interface_boundary=interface_boundary,
            order=int(order),
            material_update_order=requested_material_order,
            dirichlet_bbbnd=dirichlet_bbbnd,
            bonus_intorder=bonus_intorder,
            kelvin_mats=kelvin_mats,
            inverse=inverse,
            mu_r_initial=mu_r_initial,
            kelvin_match_exact=kelvin_match_exact,
            mu_r_by_material=mu_r_by_material,
            tolerance=tolerance,
            max_iterations=max_iterations,
            relaxation=relaxation,
            interface_constraint_scale=interface_constraint_scale,
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            kelvin_source_h=kelvin_source_h,
            total_source_h=total_source_h,
            total_source_materials=total_source_materials,
            observation_points=observation_points,
            material_log_state_initial=material_log_state_initial,
            refinement_parent_identity=normalized_refinement_parent_identity,
            refinement_parent_identity_sha256=refinement_parent_identity_sha256,
            surface_dirichlet=surface_dirichlet,
            cache_fixed_rhs=cache_fixed_rhs,
            reduced_source_load=reduced_source_load,
        )
    if requested_material_order not in (None, 0):
        raise ValueError("order=1 supports material_update_order=0 only")
    if material_log_state_initial is not None:
        raise ValueError(
            "material_log_state_initial is only valid for the projected high-order path"
        )
    B_of_H = _build_bh_interpolator(bh_array[:, :2])
    B_scale = max(float(bh_array[:, 1].max()), MU_0)
    positive = (bh_array[:, 0] > 0.0) & (bh_array[:, 1] > 0.0)
    secant_upper = (float(np.max(bh_array[positive, 1] / (MU_0 * bh_array[positive, 0])))
                    if np.any(positive) else 1.0)

    mu_elements = GridFunction(L2(mesh, order=0), name="mixed_omega_mu")
    nonlinear_elements = []
    for element in mesh.Elements(VOL):
        material = str(element.mat)
        if material in nonlinear_set:
            coordinates = [mesh.vertices[vertex.nr].point for vertex in element.vertices]
            centroid = tuple(
                sum(float(point[component]) for point in coordinates) / len(coordinates)
                for component in range(mesh.dim))
            nonlinear_elements.append((element.nr, centroid))
        else:
            mu_elements.vec[element.nr] = MU_0
    element_numbers = [int(number) for number, _ in nonlinear_elements]
    initial = np.asarray(mu_r_initial, dtype=float)
    if initial.ndim == 0:
        if not math.isfinite(float(initial)) or float(initial) <= 0.0:
            raise ValueError("mu_r_initial must be positive and finite")
        mu_r_current = np.full(len(nonlinear_elements), float(initial))
    else:
        if initial.shape != (len(nonlinear_elements),):
            raise ValueError(
                "mu_r_initial must be a scalar or one value per nonlinear element "
                f"({len(nonlinear_elements)}); got shape {initial.shape}")
        if not np.all(np.isfinite(initial)) or np.any(initial < 1.0):
            raise ValueError("per-element mu_r_initial must be finite and >= 1")
        mu_r_current = initial.astype(float, copy=True)
    from scipy.interpolate import PchipInterpolator
    zero_field_mu_r = float(PchipInterpolator(
        bh_array[:, 0], bh_array[:, 1]).derivative()(0.0)) / MU_0
    if not math.isfinite(zero_field_mu_r) or zero_field_mu_r < 1.0:
        raise ValueError("bh_table origin tangent must be finite and >= vacuum permeability")
    mu_r_zero_field = np.full(len(nonlinear_elements), zero_field_mu_r)
    mu_r_upper = max(float(np.max(mu_r_current)) if mu_r_current.size else 1.0,
                     secant_upper, 1.0)
    for index, (element_nr, _) in enumerate(nonlinear_elements):
        mu_elements.vec[element_nr] = MU_0 * float(mu_r_current[index])

    kelvin_mu = make_kelvin_mu_cf(
        mesh, R_K, offset, kelvin_mats=kelvin_mats, mu_r_by_material={},
        kelvin_match_exact=kelvin_match_exact)

    def mixed_mu_cf():
        values = {}
        for material in mesh.GetMaterials():
            if material in nonlinear_set:
                values[material] = mu_elements
            elif _is_kelvin_material(material, kelvin_mats, exact=kelvin_match_exact):
                values[material] = kelvin_mu
            else:
                values[material] = MU_0 * float(mu_r_by_material.get(material, 1.0))
        return mesh.MaterialCF(values, default=MU_0)

    observation = None
    if observation_points is not None:
        observation = np.asarray(observation_points, dtype=float).reshape(-1, 3)
        for point in observation:
            if not mesh(*map(float, point)):
                raise ValueError(
                    f"observation point lies outside the mesh: {point.tolist()}")

    accelerator = ConstrainedAndersonAccelerator(
        depth=int(anderson_depth), relaxation=float(relaxation),
        lower=1.0, upper=mu_r_upper, transform=str(anderson_transform))
    history = []
    B_previous = np.zeros(len(nonlinear_elements))
    observed_previous = None
    observed_field = None
    result = None
    converged = False
    relative_change = float("inf")
    rhs_cache = {} if cache_fixed_rhs else None
    cache_material = (mu_elements if cache_fixed_rhs and total_source_h is not None and
                      set(total_source_materials) <= nonlinear_set else None)
    for iteration in range(1, int(max_iterations) + 1):
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, H_s, source_potential, R_K, offset,
            mu_r_by_material=None, reduced_materials=reduced_materials,
            total_materials=total_materials, interface_boundary=interface_boundary,
            order=order, dirichlet_bbbnd=dirichlet_bbbnd,
            bonus_intorder=bonus_intorder, kelvin_mats=kelvin_mats,
            inverse=inverse, interface_constraint_scale=interface_constraint_scale,
            kelvin_match_exact=kelvin_match_exact,
            mu_cf=mixed_mu_cf(),
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            kelvin_source_h=kelvin_source_h,
            total_source_h=total_source_h,
            total_source_materials=total_source_materials,
            reduced_zero_normal_boundary=reduced_zero_normal_boundary,
            surface_dirichlet=surface_dirichlet,
            _fixed_rhs_cache=rhs_cache, _rhs_material=cache_material,
            reduced_source_load=reduced_source_load)
        B_current = np.zeros(len(nonlinear_elements))
        mu_r_target = np.empty(len(nonlinear_elements))
        for index, (element_nr, centroid) in enumerate(nonlinear_elements):
            H_value = result["H_cf"](mesh(*centroid))
            H_magnitude = math.sqrt(sum(float(value) ** 2 for value in H_value))
            B_current[index] = B_of_H(H_magnitude)
            if H_magnitude <= 1.0e-12:
                mu_r_target[index] = mu_r_zero_field[index]
            else:
                mu_r_target[index] = max(1.0, B_current[index] / (MU_0 * H_magnitude))
        entry = {"iteration": int(iteration)}
        if observation is not None:
            # Evaluate before the permeability update so the recorded field is the
            # one this iteration actually solved for.
            observed_field = np.asarray(
                [[float(value) for value in result["B_cf"](mesh(*point))]
                 for point in observation], dtype=float)
            if observed_previous is not None:
                scale = max(float(np.max(np.linalg.norm(observed_field, axis=1))), 1.0e-30)
                entry["observation_relative_change"] = float(
                    np.max(np.linalg.norm(observed_field - observed_previous, axis=1)) / scale)
            observed_previous = observed_field
        constitutive_change = float(np.max(
            np.abs(mu_r_target - mu_r_current) / np.maximum(mu_r_target, 1.0)))
        entry["relative_constitutive_change"] = constitutive_change
        if iteration > 1:
            relative_change = float(
                np.max(np.abs(B_current - B_previous)) / max(B_scale, 1.0e-30))
            entry["relative_B_change"] = relative_change
        entry["mu_r_min"] = float(np.min(mu_r_current)) if mu_r_current.size else None
        entry["mu_r_max"] = float(np.max(mu_r_current)) if mu_r_current.size else None
        history.append(entry)
        # Keep the returned coefficient at the state used to assemble H.
        if iteration > 1 and max(relative_change, constitutive_change) <= tolerance:
            converged = True
            break
        # PCHIP secants can exceed the secants at the tabulated nodes.
        # A safeguard must not clip the constitutive law's own target.
        if mu_r_target.size:
            accelerator.upper = max(accelerator.upper, float(np.max(mu_r_target)))
        if iteration == 1 and initial.ndim == 0:
            # A scalar start is a guess: jump to the first material estimate
            # undamped.  A per-element warm start is an iterate of this very
            # loop, so it takes the damped / accelerated update from the start.
            mu_r_next = np.maximum(mu_r_target, 1.0)
        else:
            mu_r_next = accelerator.step(mu_r_current, mu_r_target)
        mu_r_current = np.asarray(mu_r_next, dtype=float)
        for index, (element_nr, _) in enumerate(nonlinear_elements):
            mu_elements.vec[element_nr] = MU_0 * float(mu_r_current[index])
        B_previous = B_current

    if result is None:  # pragma: no cover - guarded by max_iterations validation
        raise RuntimeError("mixed total/reduced Omega Picard iteration did not start")
    stats = {
        "method": "Picard",
        "rhs_reuse": "fixed_vector_and_material_operator" if cache_material is not None
                     else "fixed_vector" if rhs_cache is not None else "none",
        "iterations": iteration,
        "converged": converged,
        "relative_B_change": relative_change,
        "tolerance": float(tolerance),
        "relaxation": float(relaxation),
        "anderson_depth": int(anderson_depth),
        "response_order": int(order),
        "material_update_order": 0,
        "material_sampling": "element_centroid",
        "relative_constitutive_change": constitutive_change,
        "returned_material_state": "assembled_iterate",
        "warm_start": bool(initial.ndim > 0),
        "history": history,
        "contraction_rate_estimate": estimate_contraction_rate(
            [row["relative_B_change"] for row in history if "relative_B_change" in row]),
        "mu_r_elements": mu_r_current.tolist(),
        "element_numbers": element_numbers,
        "anderson": accelerator.stats(),
    }
    if observed_field is not None:
        stats["observation_points_m"] = observation.tolist()
        stats["observation_field_T"] = observed_field.tolist()
    if not converged:
        raise MixedOmegaPicardNotConverged(
            "mixed total/reduced Omega Picard iteration did not converge: "
            f"iterations={iteration}, relative_B_change={relative_change:.3e}, "
            f"tolerance={tolerance:.3e}",
            {"mu_r_elements": mu_r_current.tolist(),
             "element_numbers": element_numbers,
             "nonlinear_stats": stats})
    result["mu_cf"] = mixed_mu_cf()
    result["B_cf"] = result["mu_cf"] * result["H_cf"]
    result["nonlinear_stats"] = stats
    result["constitutive_field_audit"] = audit_mixed_omega_constitutive_field(
        mesh, result["H_cf"], result["B_cf"], bh_array, nonlinear_materials,
        integration_order=max(4, 2 * int(order) + int(bonus_intorder)))
    return result


def _solve_mixed_omega_projected_log_material(
        mesh, H_s, source_potential, R_K, offset, *, bh_array,
        nonlinear_materials, reduced_materials, total_materials,
        interface_boundary, order, material_update_order,
        dirichlet_bbbnd, bonus_intorder, kelvin_mats, inverse,
        mu_r_initial, tolerance, max_iterations, relaxation,
        interface_constraint_scale, kelvin_interface_boundary,
        kelvin_source_potential, kelvin_source_h,
        total_source_h, total_source_materials, observation_points,
        material_log_state_initial, refinement_parent_identity,
        refinement_parent_identity_sha256, surface_dirichlet=None,
        kelvin_match_exact=False, mu_r_by_material=None, cache_fixed_rhs=True,
        reduced_source_load="volume"):
    """Picard lane with a positive spatial L2 secant-permeability field."""

    from ngsolve import (
        GridFunction,
        IfPos,
        InnerProduct,
        Integrate,
        L2,
        exp,
        log,
        sqrt,
    )
    from radia.picard_acceleration import estimate_contraction_rate
    from radia.scalar_potential_solver import (
        _build_bh_coefficient_function,
        _build_bh_coenergy_coefficient_function,
    )

    initial = np.asarray(mu_r_initial, dtype=float)
    if initial.ndim != 0:
        raise ValueError(
            "projected high-order material updates currently require scalar "
            "mu_r_initial; resume states use nonlinear_stats.material_log_state_dofs"
        )
    initial_mu_r = float(initial)
    if not math.isfinite(initial_mu_r) or initial_mu_r < 1.0:
        raise ValueError("mu_r_initial must be finite and >= 1")

    # The admissible permeability ceiling must come from the constitutive law
    # as it is actually evaluated -- the PCHIP interpolant -- not from the
    # tabulated nodes alone.  Between and below the nodes the interpolated
    # secant B/(mu0 H) exceeds every node secant: on the Picard test table the
    # nodes top out at 2000 while the law reaches 2620, so a node-derived cap
    # silently rewrote the material by up to 31%, and by a different amount
    # for each mu_r_initial.  The P1 lane already widened its cap to the
    # law's own targets; this lane now does the same, and starts from the
    # law's dense maximum rather than from the nodes.
    from ngsolve import Parameter
    from radia.scalar_potential_solver import _build_bh_interpolator

    positive = (bh_array[:, 0] > 0.0) & (bh_array[:, 1] > 0.0)
    if np.any(positive):
        law = _build_bh_interpolator(bh_array)
        h_max = float(np.max(bh_array[positive, 0]))
        h_min = float(np.min(bh_array[positive, 0]))
        dense_h = np.geomspace(max(h_min * 1.0e-3, 1.0e-9), h_max, 4001)
        law_secant = np.asarray([law(h) for h in dense_h]) / (MU_0 * dense_h)
        node_secant = bh_array[positive, 1] / (MU_0 * bh_array[positive, 0])
        secant_upper = float(max(np.max(law_secant), np.max(node_secant)))
    else:
        secant_upper = 1.0
    mu_r_upper = max(initial_mu_r, secant_upper, 1.0)
    # A Parameter, so the loop can raise the ceiling whenever the law's own
    # target exceeds it; the coefficient function below reads it live.
    log_mu_upper = Parameter(math.log(mu_r_upper))

    nonlinear_materials = tuple(nonlinear_materials)
    nonlinear_set = set(nonlinear_materials)
    nonlinear_selector = mesh.Materials("|".join(nonlinear_materials))
    material_fes = L2(mesh, order=int(material_update_order))
    log_mu = GridFunction(material_fes, name="mixed_omega_log_mu")
    log_mu.Set(math.log(initial_mu_r), definedon=nonlinear_selector)
    target_log_mu = GridFunction(material_fes, name="mixed_omega_target_log_mu")
    active_mask = material_fes.GetDofs(nonlinear_selector)
    active_dofs = [index for index in range(material_fes.ndof) if active_mask[index]]
    if not active_dofs:
        raise ValueError("nonlinear material selector has no material-update degrees of freedom")
    material_state_identity = {
        "schema": "radia.mixed-omega-material-state.v1",
        "mesh_topology_geometry_sha256": _mesh_topology_geometry_sha256(mesh),
        "bh_table_sha256": _canonical_array_sha256(bh_array[:, :2]),
        "nonlinear_materials": sorted(str(value) for value in nonlinear_materials),
        "linear_mu_r_by_material": dict(mu_r_by_material or {}),
        "response_order": int(order),
        "material_update_order": int(material_update_order),
        "active_material_dof_numbers": [int(value) for value in active_dofs],
    }
    material_state_identity_sha256 = _identity_sha256(material_state_identity)
    warm_state = None
    if material_log_state_initial is not None:
        if not isinstance(material_log_state_initial, Mapping):
            raise ValueError(
                "material_log_state_initial must be the material_restart_state "
                "mapping from an earlier solve, not a bare coefficient array"
            )
        supplied_identity = material_log_state_initial.get("identity")
        supplied_identity_sha256 = _validate_restart_identity_digest(
            material_log_state_initial.get("identity_sha256")
        )
        if not isinstance(supplied_identity, Mapping):
            raise ValueError("material restart state is missing its identity mapping")
        if _identity_sha256(dict(supplied_identity)) != supplied_identity_sha256:
            raise ValueError("material restart state identity digest is internally inconsistent")
        if supplied_identity_sha256 != material_state_identity_sha256:
            raise ValueError(
                "material restart state does not match the current mesh, B-H table, "
                "material selector, response order, material order, and active DOFs"
            )
        warm_state = np.asarray(material_log_state_initial.get("values"), dtype=float)
        if warm_state.shape != (len(active_dofs),):
            raise ValueError(
                "material restart state must contain one value per active "
                f"material degree of freedom ({len(active_dofs)}); got {warm_state.shape}"
            )
        if not np.all(np.isfinite(warm_state)):
            raise ValueError("material_log_state_initial must be finite")
        log_mu.vec.FV().NumPy()[active_dofs] = warm_state

    bounded_log_mu = IfPos(
        log_mu,
        IfPos(log_mu_upper - log_mu, log_mu, log_mu_upper),
        0.0,
    )
    kelvin_mu = make_kelvin_mu_cf(
        mesh, R_K, offset, kelvin_mats=kelvin_mats, mu_r_by_material={},
        kelvin_match_exact=kelvin_match_exact
    )

    def mixed_mu_cf():
        values = {}
        for material in mesh.GetMaterials():
            if material in nonlinear_set:
                values[material] = MU_0 * exp(bounded_log_mu)
            elif _is_kelvin_material(material, kelvin_mats, exact=kelvin_match_exact):
                values[material] = kelvin_mu
            else:
                values[material] = MU_0 * float(mu_r_by_material.get(material, 1.0))
        return mesh.MaterialCF(values, default=MU_0)

    observation = None
    if observation_points is not None:
        observation = np.asarray(observation_points, dtype=float).reshape(-1, 3)
        for point in observation:
            if not mesh(*map(float, point)):
                raise ValueError(
                    f"observation point lies outside the mesh: {point.tolist()}"
                )

    b_fes = L2(mesh, order=int(material_update_order))
    b_projected = GridFunction(b_fes, name="mixed_omega_B_magnitude")
    b_previous = GridFunction(b_fes, name="mixed_omega_B_magnitude_previous")
    have_previous = False
    observed_previous = None
    observed_field = None
    history = []
    result = None
    converged = False
    relative_change = float("inf")
    constitutive_change = float("inf")
    final_state_constitutive_residual = None
    integration_order = max(4, 2 * int(material_update_order) + 2)
    rhs_cache = {} if cache_fixed_rhs else None

    def solve_current_material_state():
        return solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,
            H_s,
            source_potential,
            R_K,
            offset,
            mu_r_by_material=None,
            reduced_materials=reduced_materials,
            total_materials=total_materials,
            interface_boundary=interface_boundary,
            order=order,
            dirichlet_bbbnd=dirichlet_bbbnd,
            bonus_intorder=bonus_intorder,
            kelvin_mats=kelvin_mats,
            inverse=inverse,
            interface_constraint_scale=interface_constraint_scale,
            kelvin_match_exact=kelvin_match_exact,
            mu_cf=mixed_mu_cf(),
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            kelvin_source_h=kelvin_source_h,
            total_source_h=total_source_h,
            total_source_materials=total_source_materials,
            surface_dirichlet=surface_dirichlet,
            _fixed_rhs_cache=rhs_cache,
            reduced_source_load=reduced_source_load,
        )

    for iteration in range(1, int(max_iterations) + 1):
        result = solve_current_material_state()
        H_magnitude = sqrt(InnerProduct(result["H_cf"], result["H_cf"]) + 1.0e-24)
        B_target = _build_bh_coefficient_function(H_magnitude, bh_array)
        mu_r_target = B_target / (MU_0 * H_magnitude)
        admissible_mu_r_target = IfPos(mu_r_target - 1.0, mu_r_target, 1.0)
        # Evaluate the shared field once per quadrature rule, rather than once
        # per occurrence in every PCHIP branch. This compiles the expression
        # graph without changing its interpolation or quadrature.
        target_log_mu.Set(log(admissible_mu_r_target).Compile(), definedon=nonlinear_selector)
        # The ceiling must never clip the law's own target (P1 lane does the
        # same): raise it to whatever the projection asked for.
        target_max = float(np.max(target_log_mu.vec.FV().NumPy()[active_dofs]))
        if target_max > log_mu_upper.Get():
            log_mu_upper.Set(target_max)
        # Material fixed-point residual: how far the state that assembled
        # this solve sits from the law's answer to it, as an RMS relative
        # permeability change over the nonlinear region.  A B that stopped
        # moving because the update stalled is not a converged material.
        nonlinear_measure = float(
            Integrate(1.0, mesh, definedon=nonlinear_selector,
                      order=integration_order).real)
        constitutive_change = math.sqrt(max(float(Integrate(
            (exp(target_log_mu - log_mu) - 1.0) ** 2, mesh,
            definedon=nonlinear_selector, order=integration_order).real), 0.0)
            / max(nonlinear_measure, 1.0e-300))

        b_magnitude = sqrt(InnerProduct(result["B_cf"], result["B_cf"]) + 1.0e-30)
        b_projected.Set(b_magnitude, definedon=nonlinear_selector)
        entry = {"iteration": int(iteration)}
        if have_previous:
            difference_sq = float(
                Integrate(
                    (b_projected - b_previous) ** 2,
                    mesh,
                    definedon=nonlinear_selector,
                    order=integration_order,
                ).real
            )
            scale_sq = float(
                Integrate(
                    b_projected**2,
                    mesh,
                    definedon=nonlinear_selector,
                    order=integration_order,
                ).real
            )
            relative_change = math.sqrt(
                max(difference_sq, 0.0) / max(scale_sq, 1.0e-60)
            )
            entry["relative_B_l2_change"] = relative_change
        b_previous.vec.data = b_projected.vec
        have_previous = True

        if observation is not None:
            observed_field = np.asarray(
                [
                    [float(value) for value in result["B_cf"](mesh(*point))]
                    for point in observation
                ],
                dtype=float,
            )
            if observed_previous is not None:
                scale = max(
                    float(np.max(np.linalg.norm(observed_field, axis=1))), 1.0e-30
                )
                entry["observation_relative_change"] = float(
                    np.max(np.linalg.norm(observed_field - observed_previous, axis=1))
                    / scale
                )
            observed_previous = observed_field

        current_coefficients = log_mu.vec.FV().NumPy()
        target_coefficients = target_log_mu.vec.FV().NumPy()
        if iteration == 1:
            current_coefficients[active_dofs] = target_coefficients[active_dofs]
        else:
            current_coefficients[active_dofs] = (
                (1.0 - float(relaxation)) * current_coefficients[active_dofs]
                + float(relaxation) * target_coefficients[active_dofs]
            )
        entry["material_log_state_min"] = float(
            np.min(current_coefficients[active_dofs])
        )
        entry["material_log_state_max"] = float(
            np.max(current_coefficients[active_dofs])
        )
        entry["relative_constitutive_change"] = constitutive_change
        entry["mu_r_upper"] = float(math.exp(log_mu_upper.Get()))
        history.append(entry)
        # Both have to be small: the field between iterations AND the
        # material against the law.  Stalling under a small relaxation makes
        # the first tiny while the second stays large, and that was being
        # reported as convergence.
        if (iteration > 1 and relative_change <= tolerance
                and constitutive_change <= tolerance):
            converged = True
            break

    if result is None:  # pragma: no cover - guarded by validation
        raise RuntimeError("projected mixed Omega Picard iteration did not start")
    if converged:
        result = solve_current_material_state()
        # The returned field was assembled with the relaxed state, not the
        # state the loop tested.  Check the material against the law on THIS
        # solve; it is the one the caller receives, so it is the one that has
        # to be self-consistent.
        H_final = sqrt(InnerProduct(result["H_cf"], result["H_cf"]) + 1.0e-24)
        mu_final_target = _build_bh_coefficient_function(H_final, bh_array) / (
            MU_0 * H_final)
        target_log_mu.Set(log(IfPos(mu_final_target - 1.0, mu_final_target, 1.0)).Compile(),
                          definedon=nonlinear_selector)
        nonlinear_measure = float(
            Integrate(1.0, mesh, definedon=nonlinear_selector,
                      order=integration_order).real)
        final_state_constitutive_residual = math.sqrt(max(float(Integrate(
            (exp(target_log_mu - log_mu) - 1.0) ** 2, mesh,
            definedon=nonlinear_selector, order=integration_order).real), 0.0)
            / max(nonlinear_measure, 1.0e-300))
        if not (final_state_constitutive_residual <= tolerance):
            converged = False
        if observation is not None:
            observed_field = np.asarray(
                [
                    [float(value) for value in result["B_cf"](mesh(*point))]
                    for point in observation
                ],
                dtype=float,
            )
    state = log_mu.vec.FV().NumPy()[active_dofs].copy()
    material_restart_state = {
        "values": state.tolist(),
        "identity": material_state_identity,
        "identity_sha256": material_state_identity_sha256,
    }
    stats = {
        "method": "Picard projected log-permeability",
        "rhs_cache": "fixed_vector" if rhs_cache is not None else "none",
        "iterations": int(iteration),
        "converged": converged,
        "relative_B_change": relative_change,
        "relative_B_change_norm": "L2(nonlinear_materials)",
        "tolerance": float(tolerance),
        "relaxation": float(relaxation),
        "anderson_depth": 0,
        "anderson": None,
        "response_order": int(order),
        "material_update_order": int(material_update_order),
        "material_sampling": "L2_log_secant_projection",
        # The upper bound as it stands after the loop raised it to the law's
        # own targets; mu_r_upper alone is only where it started.
        "physical_permeability_bounds": [1.0, float(math.exp(log_mu_upper.Get()))],
        "physical_permeability_bounds_initial": [1.0, float(mu_r_upper)],
        "relative_constitutive_change": constitutive_change,
        "final_state_constitutive_residual": final_state_constitutive_residual,
        "history": history,
        "contraction_rate_estimate": estimate_contraction_rate(
            [
                row["relative_B_l2_change"]
                for row in history
                if "relative_B_l2_change" in row
            ]
        ),
        "material_log_state_dofs": state.tolist(),
        "material_state_dof_numbers": active_dofs,
        "material_state_identity": material_state_identity,
        "material_state_identity_sha256": material_state_identity_sha256,
        "material_restart_state": material_restart_state,
        "warm_start": warm_state is not None,
        "final_material_state_resolved": converged,
    }
    if observed_field is not None:
        stats["observation_points_m"] = observation.tolist()
        stats["observation_field_T"] = observed_field.tolist()
    if not converged:
        raise MixedOmegaPicardNotConverged(
            "projected mixed total/reduced Omega Picard iteration did not converge: "
            f"iterations={iteration}, relative_B_change={relative_change:.3e}, "
            f"relative_constitutive_change={constitutive_change:.3e}, "
            f"final_state_constitutive_residual={final_state_constitutive_residual}, "
            f"tolerance={tolerance:.3e}",
            {
                "material_log_state_dofs": state.tolist(),
                "material_state_dof_numbers": active_dofs,
                "material_restart_state": material_restart_state,
                "nonlinear_stats": stats,
            },
        )
    result["mu_cf"] = mixed_mu_cf()
    result["B_cf"] = result["mu_cf"] * result["H_cf"]
    H_magnitude = sqrt(InnerProduct(result["H_cf"], result["H_cf"]) + 1.0e-24)
    B_law = _build_bh_coefficient_function(H_magnitude, bh_array)
    coenergy_density = _build_bh_coenergy_coefficient_function(H_magnitude, bh_array)
    nonlinear_coenergy = float(
        Integrate(
            coenergy_density,
            mesh,
            definedon=nonlinear_selector,
            order=max(integration_order, 2 * int(order) + 2),
        ).real
    )
    nonlinear_energy = float(
        Integrate(
            H_magnitude * B_law - coenergy_density,
            mesh,
            definedon=nonlinear_selector,
            order=max(integration_order, 2 * int(order) + 2),
        ).real
    )
    nonlinear_h_dot_b = float(
        Integrate(
            H_magnitude * B_law,
            mesh,
            definedon=nonlinear_selector,
            order=max(integration_order, 2 * int(order) + 2),
        ).real
    )
    physical_linear_materials = [
        str(material)
        for material in mesh.GetMaterials()
        if material not in nonlinear_set
        and not _is_kelvin_material(material, kelvin_mats, exact=kelvin_match_exact)
    ]
    linear_energy = 0.0
    linear_h_dot_b = 0.0
    if physical_linear_materials:
        linear_selector = mesh.Materials("|".join(physical_linear_materials))
        linear_energy = float(
            Integrate(
                0.5 * result["mu_cf"] * H_magnitude**2,
                mesh,
                definedon=linear_selector,
                order=max(integration_order, 2 * int(order) + 2),
            ).real
        )
        linear_h_dot_b = float(
            Integrate(
                result["mu_cf"] * H_magnitude**2,
                mesh,
                definedon=linear_selector,
                order=max(integration_order, 2 * int(order) + 2),
            ).real
        )
    energy_total = nonlinear_energy + linear_energy
    coenergy_total = nonlinear_coenergy + linear_energy
    h_dot_b_total = nonlinear_h_dot_b + linear_h_dot_b
    legendre_residual = abs(energy_total + coenergy_total - h_dot_b_total) / max(
        abs(energy_total), abs(coenergy_total), abs(h_dot_b_total), 1.0e-300
    )
    solution_sha256 = _canonical_array_sha256(result["solution"].vec.FV().NumPy())
    energy_identity = {
        "schema": "radia.nonlinear-magnetic-energy-identity.v1",
        "mesh_topology_geometry_sha256": material_state_identity[
            "mesh_topology_geometry_sha256"
        ],
        "bh_table_sha256": material_state_identity["bh_table_sha256"],
        "constitutive_interpolation": "monotone_pchip",
        "constitutive_extrapolation": "vacuum_slope",
        "magnetic_anisotropy": "isotropic",
        "material_state_identity_sha256": material_state_identity_sha256,
        "solution_sha256": solution_sha256,
        "nonlinear_materials": sorted(str(value) for value in nonlinear_materials),
        "linear_physical_materials": sorted(physical_linear_materials),
        "material_domain": "|".join(
            sorted([*map(str, nonlinear_materials), *physical_linear_materials])
        ),
        "nonlinear_state_id": material_state_identity_sha256,
        "coordinate_system": "right-handed Cartesian",
        "unit_system": "SI",
        "unit": "J",
        "integration_order": max(integration_order, 2 * int(order) + 2),
        "refinement_parent_identity_sha256": refinement_parent_identity_sha256,
    }
    if refinement_parent_identity is not None:
        energy_identity["refinement_parent_identity"] = refinement_parent_identity
    field_identity = {
        key: value
        for key, value in energy_identity.items()
        if key not in {"schema", "integration_order", "unit"}
    }
    field_identity.update(
        {
            "schema": "radia.nonlinear-magnetic-field-identity.v1",
            "unit": "T",
        }
    )
    result["field_observable_identity"] = field_identity
    result["energy_observables"] = {
        "identity": energy_identity,
        "identity_sha256": _identity_sha256(energy_identity),
        "energy_J": energy_total,
        "coenergy_J": coenergy_total,
        "h_dot_b_integral_J": h_dot_b_total,
        "legendre_residual_relative": legendre_residual,
        "nonlinear_energy_J": nonlinear_energy,
        "nonlinear_coenergy_J": nonlinear_coenergy,
        "nonlinear_h_dot_b_integral_J": nonlinear_h_dot_b,
        "linear_physical_energy_J": linear_energy,
        "linear_physical_h_dot_b_integral_J": linear_h_dot_b,
        "nonlinear_constitutive_relation": "energy=H*B-integral_0^H_Bdh",
        "nonlinear_coenergy_relation": "coenergy=integral_0^H_Bdh",
        "nonlinear_B_evaluation": "reconstructed_BH_not_returned_B",
    }
    result["nonlinear_stats"] = stats
    result["constitutive_field_audit"] = audit_mixed_omega_constitutive_field(
        mesh, result["H_cf"], result["B_cf"], bh_array, nonlinear_materials,
        integration_order=max(4, 2 * int(order) + int(bonus_intorder)))
    return result


def _solve_mixed_omega_pointwise_picard(
        mesh, H_s, source_potential, R_K, offset, *, bh_array,
        nonlinear_materials, reduced_materials, total_materials,
        interface_boundary, order, dirichlet_bbbnd, bonus_intorder,
        kelvin_mats, inverse, mu_r_initial, tolerance, max_iterations,
        interface_constraint_scale, kelvin_interface_boundary,
        kelvin_source_potential, kelvin_source_h, total_source_h,
        total_source_materials, reduced_zero_normal_boundary,
        surface_dirichlet, kelvin_match_exact, mu_r_by_material,
        bh_interpolation, progress_callback, reduced_source_load="volume"):
    """P1 secant Picard with the coefficient evaluated by volume quadrature."""
    from ngsolve import IfPos, sqrt
    from scipy.interpolate import PchipInterpolator
    from radia.scalar_potential_solver import (
        _build_bh_coefficient_function,
        _build_bh_linear_spline_coefficient_function,
        _build_bh_coenergy_coefficient_function,
        _build_bh_linear_spline_coenergy_coefficient_function)
    build_b = (_build_bh_coefficient_function if bh_interpolation == "pchip"
               else _build_bh_linear_spline_coefficient_function)

    initial = np.asarray(mu_r_initial, dtype=float)
    if initial.ndim != 0 or not math.isfinite(float(initial)) or float(initial) <= 0:
        raise ValueError("integration_point mu_r_initial must be a positive scalar")
    if bh_interpolation == "pchip":
        origin_mu = float(PchipInterpolator(
            bh_array[:, 0], bh_array[:, 1]).derivative()(0.))
    else:
        origin_mu = float((bh_array[1, 1] - bh_array[0, 1])
                          / (bh_array[1, 0] - bh_array[0, 0]))
    if not math.isfinite(origin_mu) or origin_mu <= 0:
        raise ValueError("B(H) origin tangent must be positive and finite")
    material_names = tuple(nonlinear_materials)
    selector = mesh.Materials("|".join(material_names))
    current_mu = MU_0 * float(initial)
    previous_b = None
    history = []
    fixed_rhs_cache = {}
    integration_order = max(4, 2 * int(order) + int(bonus_intorder))

    for iteration in range(1, int(max_iterations) + 1):
        iteration_started = time.perf_counter()
        coefficients = dict(mu_r_by_material)
        mu = make_kelvin_mu_cf(
            mesh, R_K, offset, kelvin_mats=kelvin_mats,
            mu_r_by_material=coefficients, kelvin_match_exact=kelvin_match_exact)
        values = {name: current_mu for name in material_names}
        mu = mesh.MaterialCF(values, default=mu)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, H_s, source_potential, R_K, offset,
            mu_cf=mu, reduced_materials=reduced_materials,
            total_materials=total_materials, interface_boundary=interface_boundary,
            order=order, dirichlet_bbbnd=dirichlet_bbbnd,
            bonus_intorder=bonus_intorder, kelvin_mats=kelvin_mats,
            inverse=inverse, interface_constraint_scale=interface_constraint_scale,
            kelvin_match_exact=kelvin_match_exact,
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            kelvin_source_h=kelvin_source_h,
            total_source_h=total_source_h,
            total_source_materials=total_source_materials,
            reduced_zero_normal_boundary=reduced_zero_normal_boundary,
            surface_dirichlet=surface_dirichlet,
            _fixed_rhs_cache=fixed_rhs_cache,
            reduced_source_load=reduced_source_load)
        linear_solve_seconds = time.perf_counter() - iteration_started
        audit = audit_mixed_omega_constitutive_field(
            mesh, result["H_cf"], result["B_cf"], bh_array, material_names,
            integration_order=integration_order, interpolation=bh_interpolation)
        relative_change = float("inf")
        if previous_b is not None:
            difference = result["B_cf"] - previous_b
            numerator = float(Integrate(
                InnerProduct(difference, difference), mesh,
                definedon=selector, order=integration_order).real)
            denominator = float(Integrate(
                InnerProduct(result["B_cf"], result["B_cf"]), mesh,
                definedon=selector, order=integration_order).real)
            relative_change = math.sqrt(max(numerator, 0.) / max(denominator, 1.e-60))
        history.append({"iteration": iteration,
                        "linear_solve_seconds": linear_solve_seconds,
                        "iteration_seconds": time.perf_counter() - iteration_started,
                        "relative_B_change": relative_change,
                        "relative_B_constitutive_L2": audit["relative_B_constitutive_L2"]})
        if progress_callback is not None:
            progress_callback(dict(history[-1]))
        if (relative_change <= tolerance
                and audit["relative_B_constitutive_L2"] <= tolerance):
            magnitude = sqrt(InnerProduct(result["H_cf"], result["H_cf"]))
            build_energy = (
                _build_bh_coenergy_coefficient_function
                if bh_interpolation == "pchip"
                else _build_bh_linear_spline_coenergy_coefficient_function)
            coenergy_density = build_energy(magnitude, bh_array)
            coenergy = float(Integrate(
                coenergy_density, mesh, definedon=selector,
                order=integration_order).real)
            h_dot_b = float(Integrate(
                InnerProduct(result["H_cf"], result["B_cf"]), mesh,
                definedon=selector, order=integration_order).real)
            result["constitutive_field_audit"] = audit
            result["energy_observables"] = {
                "scope": "nonlinear_materials_only",
                "bh_interpolation": bh_interpolation,
                "coenergy_J": coenergy,
                "energy_J": h_dot_b - coenergy,
                "h_dot_b_integral_J": h_dot_b,
                "integration_order": integration_order,
                "accuracy_accepted": False,
            }
            result["nonlinear_stats"] = {
                "method": "pointwise_secant_picard",
                "material_sampling": "integration_point",
                "bh_interpolation": bh_interpolation,
                "iterations": iteration,
                "converged": True,
                "relative_B_change": relative_change,
                "relative_B_constitutive_L2": audit["relative_B_constitutive_L2"],
                "history": history,
            }
            return result
        previous_b = result["B_cf"]
        magnitude = sqrt(InnerProduct(result["H_cf"], result["H_cf"]))
        law_b = build_b(magnitude, bh_array)
        current_mu = IfPos(magnitude - 1.e-12,
                           law_b / (magnitude + 1.e-30), origin_mu)

    raise MixedOmegaPicardNotConverged(
        "pointwise mixed Omega Picard did not converge: "
        f"iterations={max_iterations}, relative_B_change={relative_change:.3e}, "
        f"relative_B_constitutive_L2={audit['relative_B_constitutive_L2']:.3e}",
        {"nonlinear_stats": {"method": "pointwise_secant_picard",
                             "material_sampling": "integration_point",
                             "bh_interpolation": bh_interpolation,
                             "converged": False, "history": history}})


def audit_mixed_omega_constitutive_field(mesh, H_cf, B_cf, bh_table,
                                       nonlinear_materials, *, integration_order=6,
                                       interpolation="pchip"):
    """Measure returned B against B(H), independently of iterate convergence.

    This quadrature diagnostic is not a Maxwell residual or an error bound.
    A Legendre identity evaluated using a reconstructed B(H) cannot replace it.
    """
    from ngsolve import Integrate, InnerProduct, sqrt
    from radia.scalar_potential_solver import (
        _build_bh_coefficient_function,
        _build_bh_linear_spline_coefficient_function)

    names = tuple(nonlinear_materials)
    if not names or not set(names) <= set(mesh.GetMaterials()):
        raise ValueError("nonlinear_materials must name existing mesh materials")
    if int(integration_order) != integration_order or integration_order < 1:
        raise ValueError("integration_order must be a positive integer")
    if interpolation not in ("pchip", "linear_spline"):
        raise ValueError("interpolation must be pchip or linear_spline")
    magnitude = sqrt(InnerProduct(H_cf, H_cf) + 1.e-30)
    build_b = (_build_bh_coefficient_function if interpolation == "pchip"
               else _build_bh_linear_spline_coefficient_function)
    law_b = build_b(magnitude, np.asarray(bh_table, dtype=float))
    target = law_b * H_cf / magnitude
    delta = B_cf - target
    region = mesh.Materials("|".join(names))
    def integral(value):
        return float(Integrate(value.Compile(), mesh, definedon=region,
                               order=int(integration_order)).real)
    defect = integral(InnerProduct(delta, delta))
    reference = integral(InnerProduct(target, target))
    if not math.isfinite(defect) or not math.isfinite(reference):
        raise ValueError("non-finite constitutive field audit")
    return {
        "relative_B_constitutive_L2": math.sqrt(max(defect, 0.) / max(reference, 1.e-60)),
        "absolute_B_constitutive_L2_T_m32": math.sqrt(max(defect, 0.)),
        "integration_order": int(integration_order),
        "interpolation": interpolation,
        "scope": "returned_B_against_pointwise_BH",
        "accuracy_accepted": False,
    }


def inductance_from_energy(gfu, nu_cf, mesh, I_total,
                            material_filter=None, order=10):
    """Extract self-inductance from the magnetic energy integral.

        L = 2 W / I^2,    W = 0.5 * int nu(r) |curl A|^2 dV

    Args:
        gfu: A GridFunction (total A).
        nu_cf: Kelvin-modulated reluctivity CF.
        mesh: NGSolve Mesh.
        I_total: coil total current (A).
        material_filter: optional NGSolve Materials selector
            (pipe-joined string, e.g. ``"air|coil"``) to restrict the
            energy integral. If None, integrates over the whole mesh.
        order: integration order.

    Returns:
        Tuple ``(L, W)`` in SI units.
    """
    integrand = 0.5 * nu_cf * InnerProduct(curl(gfu), Conj(curl(gfu)))
    if material_filter is None:
        W = Integrate(integrand, mesh, order=order).real
    else:
        W = Integrate(integrand, mesh,
                       definedon=mesh.Materials(material_filter),
                       order=order).real
    return 2.0 * W / (I_total ** 2), W
