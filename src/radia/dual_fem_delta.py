"""Complementary magnetostatic diagnostics for a uniform field at infinity."""

from hashlib import sha256
from math import isfinite, pi, prod, sqrt
from time import perf_counter


def evaluate_delta(mesh, h_a, h_omega, permeability, applied_h, susceptibility,
                   *, iron, air, kelvin, integration_order):
    """Evaluate the complementary moment gap and field-energy defect in joules.

    Both fields must use the same mesh and coefficients. In the Kelvin region
    pass only the reaction fields, in the twisted computational frame. The
    uniform background cancels in the difference and has a singular pullback.
    The defect is a consistency diagnostic, without a pass threshold.
    """
    from ngsolve import InnerProduct, Integrate
    from radia.kelvin_material import MU_0

    moment_a = float(Integrate(susceptibility * InnerProduct(h_a, applied_h),
                               mesh, definedon=iron, order=integration_order))
    moment_omega = float(Integrate(susceptibility * InnerProduct(h_omega, applied_h),
                                   mesh, definedon=iron, order=integration_order))
    difference = h_omega - h_a
    integrand = permeability * InnerProduct(difference, difference)
    energies = {
        name: float(Integrate(integrand, mesh, definedon=region, order=integration_order))
        for name, region in (("iron", iron), ("air", air), ("kelvin", kelvin))
    }
    # Reintegrate the linear load and physical energy of each reaction field.
    # Their differences expose assembly quadrature and A gauge contributions.
    b_reaction = permeability * h_a - MU_0 * applied_h
    h_reaction = h_omega - applied_h
    load_a = susceptibility / (1 + susceptibility) * InnerProduct(applied_h, b_reaction)
    load_omega = -MU_0 * susceptibility * InnerProduct(applied_h, h_reaction)
    kappa_a = float(Integrate(load_a, mesh, definedon=iron, order=integration_order))
    kappa_omega = float(Integrate(load_omega, mesh, definedon=iron, order=integration_order))
    for region in (iron, air):
        kappa_a -= float(Integrate(InnerProduct(b_reaction, b_reaction) / permeability,
                                   mesh, definedon=region, order=integration_order))
        kappa_omega -= float(Integrate(permeability * InnerProduct(h_reaction, h_reaction),
                                       mesh, definedon=region, order=integration_order))
    kappa_a -= float(Integrate(permeability * InnerProduct(h_a, h_a), mesh,
                               definedon=kelvin, order=integration_order))
    kappa_omega -= float(Integrate(permeability * InnerProduct(h_omega, h_omega), mesh,
                                   definedon=kelvin, order=integration_order))
    gap = MU_0 * (moment_omega - moment_a)
    delta = gap - sum(energies.values())
    if not all(isfinite(value) for value in (moment_a, moment_omega, gap, delta, kappa_a, kappa_omega)):
        raise ValueError("Nonfinite complementary FEM diagnostic")
    return {
        "moment_a": moment_a,
        "moment_omega": moment_omega,
        "moment_2_to_1": (2 * moment_a + moment_omega) / 3,
        "moment_midpoint": (moment_a + moment_omega) / 2,
        "moment_gap_j": gap,
        "difference_energy_j": energies,
        "delta_j": delta,
        "load_minus_energy_j": {"a": kappa_a, "omega": kappa_omega},
        "orthogonality_residual_j": delta - kappa_a - kappa_omega,
        "relative_delta": delta / gap if gap > 0 else None,
    }


def ellipsoid_delta_benchmark(*, order=2, curve_order=2, iron_maxh=0.45,
                              mu_r=20.0, integration_order=12,
                              air_divisions=5, assembly_bonus=0):
    """Run one paired FEM solve on the original triaxial benchmark."""
    return ellipsoid_delta_series(
        orders=(order,), curve_order=curve_order, iron_maxh=iron_maxh,
        mu_r=mu_r, integration_order=integration_order,
        air_divisions=air_divisions, assembly_bonus=assembly_bonus)["rows"][0]


def ellipsoid_delta_series(*, orders=(2, 3, 4), curve_order=2,
                           shape="triaxial", iron_maxh=0.45, mu_r=20.0,
                           integration_order=12, air_divisions=5, assembly_bonus=0, vol_path=""):
    """Compare paired FEM solves at several orders on one verified curved mesh.

    The CAD reference is analytic. The solver and consistency checks apply to
    the common polynomial mesh. The caller owns the NGSolve TaskManager.
    """
    orders = tuple(orders)
    if not orders or tuple(sorted(set(orders))) != orders:
        raise ValueError("orders must be nonempty, unique and increasing")
    if shape not in ("triaxial", "prolate"):
        raise ValueError("shape must be triaxial or prolate")
    if shape == "prolate" and not vol_path:
        raise ValueError("The prolate benchmark requires an affine Cubit vol_path")
    order = max(orders)
    if not 2 <= curve_order <= min(orders) <= order <= 4:
        raise ValueError("Require 2 <= curve_order <= order <= 4")
    if not 0.25 <= iron_maxh <= 0.75 or not 1 <= mu_r <= 10000:
        raise ValueError("Require 0.25 <= iron_maxh <= 0.75 and 1 <= mu_r <= 10000")
    if not 2 * order + 2 <= integration_order <= 24:
        raise ValueError("Integration order must be between 2*order+2 and 24")
    if not 4 <= air_divisions <= 6 or not 0 <= assembly_bonus <= 12:
        raise ValueError("Require 4 <= air_divisions <= 6 and 0 <= assembly_bonus <= 12")

    import ngsolve as ng
    from netgen import occ
    from radia.dual_fem_geometry import (audit_curved_kelvin_mesh, curve_analytic_ellipsoid,
                                         ellipsoid_boundary_residuals)
    from radia.kelvin_geometry import add_kelvin_exterior_domain, ground_kelvin_scalar_mesh

    started = perf_counter()
    axes = (1.0, 0.8, 0.6) if shape == "triaxial" else (1.0, 0.7, 0.7)
    direction = (2 / sqrt(14), -1 / sqrt(14), 3 / sqrt(14))
    radius = 2.75
    offset = (3 * radius, 0.0, 0.0)
    mesh_digest = None
    if vol_path:
        from pathlib import Path
        from radia.kelvin_identify_ngsolve import add_kelvin_identification, has_kelvin_identification

        path = Path(vol_path)
        if path.suffix.lower() != ".vol" or not path.is_file():
            raise ValueError("vol_path must be a local affine .vol file")
        with path.open("rb") as stream:
            mesh_digest = sha256(stream.read()).hexdigest()
        mesh = ng.Mesh(str(path))
        if (mesh.GetCurveOrder() != 1 or not 0 < mesh.ne <= 16000
                or set(mesh.GetMaterials()) != {"iron", "air", "kelvin"}):
            raise ValueError("Expected an affine iron/air/kelvin mesh with at most 16000 elements")
        if set(mesh.GetBoundaries()) != {"iron_air", "kelvin_int", "kelvin_ext"}:
            raise ValueError("Expected only iron_air, kelvin_int and kelvin_ext boundaries")
        if not has_kelvin_identification(mesh):
            pairs = add_kelvin_identification(mesh, offset, point_tolerance=1e-9)
            if pairs["n_unmatched"]:
                raise ValueError("Unmatched Kelvin vertices in affine mesh")
            mesh = ng.Mesh(mesh.ngmesh)
    else:
        iron = occ.Ellipsoid(occ.Axes(occ.Pnt(0, 0, 0), occ.Z, occ.X),
                             axes[2], axes[0], axes[1])
        iron.mat("iron")
        iron.maxh = iron_maxh
        iron.faces.name = "iron_air"
        ball = occ.Sphere(occ.Pnt(0, 0, 0), radius)
        ball.faces.name = "kelvin_int"
        air = ball - iron
        air.mat("air")
        air_maxh = radius / air_divisions
        air.maxh = air_maxh
        geometry, _ = add_kelvin_exterior_domain([air, iron], offset, radius,
                                                 inner_maxh=air_maxh, outer_maxh_factor=1.0)
        mesh = ng.Mesh(occ.OCCGeometry(geometry).GenerateMesh(maxh=air_maxh))

    mesh = ground_kelvin_scalar_mesh(mesh, offset, name="ellipsoid_gauge")
    mesh_seconds = perf_counter() - started

    started = perf_counter()
    if vol_path:
        curve_analytic_ellipsoid(mesh, axes, offset, radius, curve_order)
    else:
        mesh.Curve(curve_order)
    audit = audit_curved_kelvin_mesh(
        mesh, offset, radius, deformation_order=curve_order if vol_path else None)
    audit["boundary_levelset_rms"] = ellipsoid_boundary_residuals(
        mesh, axes, offset, radius, integration_order)
    if audit["failed_tetrahedra_by_material"]:
        raise ValueError(f"Folded or singular tetrahedra: {audit['failed_tetrahedra_by_material']}")
    if audit["kelvin_curved_match_relative"] > 1e-10:
        raise ValueError(f"Nonmatching curved Kelvin seam: {audit['kelvin_curved_match_max_m']:.6g} m")
    geometry_seconds = perf_counter() - started
    rows = []
    for p in orders:
        row = _solve_ellipsoid(mesh, axes, direction, radius, offset, p, mu_r,
                               integration_order, assembly_bonus)
        row.update(curve_order=curve_order, iron_maxh=None if vol_path else iron_maxh, shape=shape,
                   air_divisions=None if vol_path else air_divisions, assembly_bonus=assembly_bonus,
                   mesh_source="affine_vol" if vol_path else "occ",
                   geometry=audit)
        row["seconds"].update(mesh=mesh_seconds, geometry=geometry_seconds)
        rows.append(row)
    return {"shape": shape, "curve_order": curve_order, "geometry": audit,
            "affine_vol_sha256": mesh_digest, "rows": rows}


def _solve_ellipsoid(mesh, axes, direction, radius, offset, order, mu_r,
                     integration_order, assembly_bonus):
    import ngsolve as ng
    from scipy.special import elliprd
    from radia.kelvin_material import MU_0
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    from radia.vector_potential_solver import VectorPotentialSolver

    source = ng.CoefficientFunction(direction)
    source_potential = -sum(h * x for h, x in zip(direction, (ng.x, ng.y, ng.z)))
    started = perf_counter()
    vector = VectorPotentialSolver(mesh, iron_domains="iron", mu_r=mu_r, order=order,
                                   kelvin_region="kelvin", kelvin_radius=radius,
                                   kelvin_center=offset)
    vector.set_source_cf(MU_0 * source, kelvin_source_cf=ng.CoefficientFunction((0., 0., 0.)))
    vector.solve_linear(dirichlet="", solver="sparsecholesky", eps=1e-10,
                        kelvin_eps=1e-10, bonus_intorder=assembly_bonus)
    a_seconds = perf_counter() - started

    started = perf_counter()
    omega = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
        mesh, source, source_potential, radius, offset,
        mu_r_by_material={"iron": mu_r}, reduced_materials=("air",),
        total_materials=("iron", "kelvin"), interface_boundary="iron_air",
        kelvin_interface_boundary="kelvin_int",
        kelvin_source_potential=ng.CoefficientFunction(0.),
        reduced_zero_normal_boundary="kelvin_int",
        dirichlet_bbbnd="ellipsoid_gauge", order=order, bonus_intorder=4 + assembly_bonus)
    omega_seconds = perf_counter() - started

    started = perf_counter()
    iron_region = mesh.Materials("iron")
    result = evaluate_delta(
        mesh, vector.get_H(), omega["H_cf"], omega["mu_cf"], source, mu_r - 1,
        iron=iron_region, air=mesh.Materials("air"), kelvin=mesh.Materials("kelvin"),
        integration_order=integration_order)
    volume = float(ng.Integrate(1.0, mesh, definedon=iron_region, order=integration_order))
    trace = omega["phi_total"].Trace() - omega["phi_reduced"].Trace() - source_potential
    interface = mesh.Boundaries("iron_air")
    trace_l2 = sqrt(max(0.0, float(ng.Integrate(
        trace * trace, mesh, ng.BND, definedon=interface, order=integration_order))))
    exact_volume = 4 * pi * prod(axes) / 3
    factors = [
        prod(axes) * elliprd(axes[(i + 1) % 3]**2, axes[(i + 2) % 3]**2, a**2) / 3
        for i, a in enumerate(axes)
    ]
    exact_moment = float(exact_volume * sum(
        (mu_r - 1) * h*h / (1 + (mu_r - 1) * factor)
        for h, factor in zip(direction, factors)))
    error = {
        name: (result["moment_" + name] / exact_moment - 1) if exact_moment else None
        for name in ("a", "omega", "2_to_1", "midpoint")
    }
    return {
        **result,
        "order": order,
        "mu_r": mu_r, "integration_order": integration_order,
        "semiaxes_m": axes, "applied_h_a_per_m": direction,
        "elements": mesh.ne, "ndof_a": vector.get_A().space.ndof,
        "ndof_omega": omega["fes"].ndof,
        "volume_m3": volume, "exact_volume_m3": exact_volume,
        "iron_interface_potential_jump_l2": trace_l2,
        "omega_free_dof_residual": omega["linear_residual"]["free_dofs"]["relative"],
        "exact_moment": exact_moment, "relative_moment_error": error,
        "cad_moment_bracketed": result["moment_a"] <= exact_moment <= result["moment_omega"],
        "seconds": {"a": a_seconds, "omega": omega_seconds,
                    "diagnostic": perf_counter() - started},
    }


def _checked_bounds(level, name):
    """Return finite ``(q_lower, q_upper)`` with q_lower <= q_upper, else raise."""
    try:
        lower, upper = float(level["q_lower"]), float(level["q_upper"])
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"{name} needs numeric q_lower and q_upper") from None
    if not (isfinite(lower) and isfinite(upper)):
        raise ValueError(f"{name} needs finite bounds, got q_lower={lower}, q_upper={upper}")
    if lower > upper:
        raise ValueError(f"{name} needs q_lower <= q_upper, got q_lower={lower}, q_upper={upper}")
    return lower, upper


def dual_bracket_reference(levels, coarse=None, convergence_fraction=0.1):
    """Freeze a reference value from converged lower/upper dual bounds.

    ``levels`` lists refinement levels from coarse to fine, each a mapping with
    ``q_lower`` (A side) and ``q_upper`` (Omega side).  They must be the
    energy-functional (kappa-corrected) bounds of conforming solutions on one
    geometry.  The finest level gives the reference: its midpoint, with its
    half-width as uncertainty.  The previous level gives the last-step drift of
    each bound.

    ``coarse`` (optional, a level not used for the reference) is then judged
    against the frozen reference.  The midpoint error is e_A (r-1)/2 and the 2:1
    error is e_A (r-2)/3, with r = (Q_upper - Q*) / (Q* - Q_lower).  So the
    midpoint is closer iff r < 1.4.  r is a property of the case, mesh and
    order, not a constant, so no fixed weight is assumed.  The reference is
    marked converged when both last-step drifts are at most
    ``convergence_fraction`` of the coarse |midpoint - 2:1| difference (or of
    the final width when no coarse level is given).  A comparison is
    inconclusive when the two absolute errors differ by no more than the
    reference width.
    """
    if len(levels) < 2:
        raise ValueError("Need at least two refinement levels")
    rows = [_checked_bounds(level, f"levels[{i}]") for i, level in enumerate(levels)]
    if coarse is not None:
        coarse_bounds = _checked_bounds(coarse, "coarse")
    nesting = [i for i in range(1, len(rows))
               if rows[i][0] < rows[i - 1][0] or rows[i][1] > rows[i - 1][1]]
    (prev_lower, prev_upper), (lower, upper) = rows[-2], rows[-1]
    reference = 0.5 * (lower + upper)
    half_width = 0.5 * (upper - lower)
    drift = {"lower": lower - prev_lower, "upper": upper - prev_upper}
    result = {
        "q_ref": reference, "half_width": half_width,
        "relative_half_width": half_width / abs(reference) if reference else None,
        "last_step_drift": drift, "nesting_violations": nesting,
    }
    scale = 2 * half_width
    if coarse is not None:
        c_lower, c_upper = coarse_bounds
        estimates = {"lower": c_lower, "upper": c_upper,
                     "midpoint": 0.5 * (c_lower + c_upper),
                     "two_to_one": (2 * c_lower + c_upper) / 3}
        errors = {name: value - reference for name, value in estimates.items()}
        difference = abs(errors["two_to_one"]) - abs(errors["midpoint"])
        scale = abs(estimates["midpoint"] - estimates["two_to_one"])
        if abs(difference) <= 2 * half_width:
            closer = "inconclusive"
        else:
            closer = "midpoint" if difference > 0 else "two_to_one"
        below = reference - c_lower
        result["coarse"] = {
            "estimates": estimates, "errors": errors,
            "relative_errors": {k: v / reference for k, v in errors.items()} if reference else None,
            "r": (c_upper - reference) / below if below > 0 else None,
            "closer": closer,
        }
    result["converged"] = bool(not nesting and max(abs(drift["lower"]), abs(drift["upper"]))
                               <= convergence_fraction * scale)
    return result
