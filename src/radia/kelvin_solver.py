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

import numpy as np

from ngsolve import (H1, HCurl, BilinearForm, LinearForm, GridFunction,
                      Periodic, Compress, CoefficientFunction, TaskManager,
                      curl, dx, ds, grad, InnerProduct, Conj, Integrate)

from radia.kelvin_material import make_kelvin_mu_cf, make_kelvin_nu_cf, MU_0, NU_0


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
    }


def project_source_total_hodge(
        mesh, H_s, total_source_materials, *, order=2, inverse="pardiso",
        gauge_epsilon=1.0e-12):
    """Split a linked source inside total-potential materials.

    On a multiply connected iron body a curl-free coil field need not be the
    gradient of one single-valued scalar.  This volume projection computes

    ``H_s = -grad(Phi_s) + H_harmonic``.

    ``Phi_s`` supplies the reduced/total interface jump and ``H_harmonic``
    carries the non-exact cohomology class in the total region.  The caller
    owns the surrounding :class:`ngsolve.TaskManager` region.
    """
    names = tuple(str(name) for name in total_source_materials)
    if not names or len(names) != len(set(names)) or any(not name for name in names):
        raise ValueError(
            "total_source_materials must contain unique non-empty names")
    if int(order) < 1:
        raise ValueError("order must be positive")
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
    d_total = dx(definedon=selector)
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += InnerProduct(grad(potential), grad(test)) * d_total
    a_bf += float(gauge_epsilon) * potential * test * d_total
    f_lf = LinearForm(fes)
    f_lf += -InnerProduct(H_s, grad(test)) * d_total
    a_bf.Assemble()
    f_lf.Assemble()
    potential_gf = GridFunction(fes, name="total_source_potential")
    potential_gf.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec

    harmonic_field = H_s + grad(potential_gf)
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
    """Solve a periodic reduced Omega--Omega magnetostatic problem with Kelvin.

    ``H_s`` is the complete source field in computational coordinates. For a
    compact Radia current source, construct it with
    :func:`radia.kelvin_material.make_kelvin_aware_radia_H_s_cf`.

    The weak form consumes a twisted 1-form.  Do not pass a projected HDiv
    flux density multiplied by ``nu`` here: an HDiv projection preserves the
    normal trace/divergence of ``B``, but does not preserve the Kelvin Hodge
    relation or the curl-free ``H`` source contract required by this scalar
    potential formulation.
    """
    mu_cf = make_kelvin_mu_cf(
        mesh, R_K, offset, kelvin_mats=kelvin_mats,
        mu_r_by_material=mu_r_by_material)
    fes = Periodic(H1(mesh, order=order, dirichlet_bbbnd=dirichlet_bbbnd))
    phi, test = fes.TnT()
    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += mu_cf * grad(phi) * grad(test) * dx(
        bonus_intorder=bonus_intorder)
    f_lf = LinearForm(fes)
    f_lf += mu_cf * H_s * grad(test) * dx(
        bonus_intorder=bonus_intorder)
    phi_gf = _assemble_and_solve(a_bf, f_lf, fes, inverse=inverse)
    return {
        "phi": phi_gf,
        "fes": fes,
        "mu_cf": mu_cf,
        "H_s": H_s,
        "B_cf": mu_cf * (H_s - grad(phi_gf)),
    }


def solve_magnetostatic_mixed_total_reduced_omega_kelvin(
        mesh, H_s, source_potential, R_K, offset, *, mu_r_by_material=None,
        reduced_materials, total_materials, interface_boundary,
        order=1, dirichlet_bbbnd="GND", bonus_intorder=4,
        kelvin_mats=("kelvin",), inverse="pardiso",
        interface_constraint_scale=None, total_dirichlet_cf=None,
        mu_cf=None, kelvin_interface_boundary=None,
        kelvin_source_potential=None, total_source_h=None,
        total_source_materials=()):
    """Solve the TOSCA-style mixed total/reduced Omega formulation.

    ``H_s`` is used in ``reduced_materials`` (the source enclosure), where
    ``H = H_s - grad(phi_reduced)``.  A linked coil may additionally supply
    ``total_source_h`` on ``total_source_materials``; this is the harmonic
    remainder of the iron-volume Hodge split and gives
    ``H = total_source_h - grad(phi_total)`` there.  Other total materials use
    ``H = -grad(phi_total)``.  On ``interface_boundary`` the potentials are
    coupled by the exact source trace

    ``phi_total - phi_reduced = source_potential``.

    When the reduced physical air reaches the inner Kelvin sphere, provide
    ``kelvin_interface_boundary`` and ``kelvin_source_potential`` as well.
    This adds the periodic source jump between the reduced physical-air trace
    and the Kelvin total-potential trace.  Omitting it would silently remove
    the source 0-form from the Kelvin pair, even though the physical source
    reaches the Kelvin interface.

    This prevents the source field and the reduced correction from cancelling
    inside high-permeability material.  The normal flux condition is natural
    in the variational form; the scalar-trace constraint supplies the
    tangential-H condition.  It is the finite-element form of the total /
    reduced Omega split used by TOSCA-class magnetostatic solvers.

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
            ``kelvin_interface_boundary``.  Kelvin inversion is orientation
            reversing for this twisted 0-form, so the enforced jump there is
            ``phi_total - phi_reduced = -kelvin_source_potential``.
        total_dirichlet_cf: optional non-homogeneous total-potential lift for
            a finite-domain verification problem. Production Kelvin meshes use
            the default ``None`` and a point/edge ``GND`` constraint.
        mu_cf: optional fully Kelvin-aware permeability coefficient. This is
            the narrow extension point used by the nonlinear Picard driver;
            callers must not supply a physical-space coefficient in the Kelvin
            material.
        total_source_h: optional non-exact harmonic source field retained in
            the total-potential region.
        total_source_materials: total-region materials on which
            ``total_source_h`` is defined.  Supply both arguments together.
    """
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
    if (kelvin_interface_boundary is None) != (kelvin_source_potential is None):
        raise ValueError(
            "kelvin_interface_boundary and kelvin_source_potential must be "
            "supplied together")
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

    if mu_cf is None:
        mu_cf = make_kelvin_mu_cf(
            mesh, R_K, offset, kelvin_mats=kelvin_mats,
            mu_r_by_material=mu_r_by_material)
    reduced_selector = mesh.Materials("|".join(reduced_materials))
    total_selector = mesh.Materials("|".join(total_materials))
    interface_selector = mesh.Boundaries(interface_boundary)
    kelvin_interface_selector = (
        None if kelvin_interface_boundary is None
        else mesh.Boundaries(kelvin_interface_boundary))

    # Only the total region crosses the Kelvin periodic identification.  The
    # source enclosure is intentionally an independent H1 space.
    fes_reduced = H1(mesh, order=int(order), definedon=reduced_selector)
    fes_total = Compress(Periodic(H1(
        mesh, order=int(order), definedon=total_selector,
        dirichlet_bbbnd=dirichlet_bbbnd)))
    fes_multiplier = Compress(H1(
        mesh, order=int(order), definedon=interface_selector))
    fes_kelvin_multiplier = (
        None if kelvin_interface_selector is None
        else Compress(H1(mesh, order=int(order), definedon=kelvin_interface_selector)))
    fes = fes_reduced * fes_total * fes_multiplier
    if fes_kelvin_multiplier is None:
        (phi_reduced, phi_total, multiplier), (
            test_reduced, test_total, test_multiplier) = fes.TnT()
        kelvin_multiplier = test_kelvin_multiplier = None
    else:
        fes = fes * fes_kelvin_multiplier
        (phi_reduced, phi_total, multiplier, kelvin_multiplier), (
            test_reduced, test_total, test_multiplier, test_kelvin_multiplier) = fes.TnT()

    a_bf = BilinearForm(fes, symmetric=True)
    a_bf += mu_cf * grad(phi_reduced) * grad(test_reduced) * dx(
        definedon=reduced_selector, bonus_intorder=bonus_intorder)
    a_bf += mu_cf * grad(phi_total) * grad(test_total) * dx(
        definedon=total_selector, bonus_intorder=bonus_intorder)
    d_interface = ds(definedon=interface_selector, bonus_intorder=bonus_intorder)
    jump_trial = phi_total.Trace() - phi_reduced.Trace()
    jump_test = test_total.Trace() - test_reduced.Trace()
    a_bf += interface_constraint_scale * (
        multiplier * jump_test + test_multiplier * jump_trial) * d_interface
    if kelvin_interface_selector is not None:
        d_kelvin_interface = ds(
            definedon=kelvin_interface_selector, bonus_intorder=bonus_intorder)
        kelvin_jump_trial = phi_total.Trace() - phi_reduced.Trace()
        kelvin_jump_test = test_total.Trace() - test_reduced.Trace()
        a_bf += interface_constraint_scale * (
            kelvin_multiplier * kelvin_jump_test
            + test_kelvin_multiplier * kelvin_jump_trial
        ) * d_kelvin_interface

    f_lf = LinearForm(fes)
    f_lf += mu_cf * H_s * grad(test_reduced) * dx(
        definedon=reduced_selector, bonus_intorder=bonus_intorder)
    if total_source_h is not None:
        total_source_selector = mesh.Materials("|".join(total_source_materials))
        f_lf += mu_cf * total_source_h * grad(test_total) * dx(
            definedon=total_source_selector, bonus_intorder=bonus_intorder)
    f_lf += interface_constraint_scale * test_multiplier * source_potential * d_interface
    if kelvin_interface_selector is not None:
        f_lf += -interface_constraint_scale * test_kelvin_multiplier * (
            kelvin_source_potential) * d_kelvin_interface
    a_bf.Assemble()
    f_lf.Assemble()

    solution = GridFunction(fes)
    if total_dirichlet_cf is None:
        solution.vec.data = a_bf.mat.Inverse(
            fes.FreeDofs(), inverse=inverse) * f_lf.vec
    else:
        solution.components[1].Set(
            total_dirichlet_cf, definedon=mesh.BBBoundaries(dirichlet_bbbnd))
        residual = solution.vec.CreateVector()
        residual.data = f_lf.vec - a_bf.mat * solution.vec
        solution.vec.data += a_bf.mat.Inverse(
            fes.FreeDofs(), inverse=inverse) * residual

    phi_reduced_gf, phi_total_gf, multiplier_gf = solution.components[:3]
    kelvin_multiplier_gf = (
        None if fes_kelvin_multiplier is None else solution.components[3])
    H_reduced = H_s - grad(phi_reduced_gf)
    zero_h = CoefficientFunction((0.0, 0.0, 0.0))
    total_source_by_material = mesh.MaterialCF({
        material: total_source_h
        if total_source_h is not None and material in total_source_materials
        else zero_h
        for material in total_materials
    })
    H_total = total_source_by_material - grad(phi_total_gf)
    h_components = []
    for component in range(3):
        values = {
            material: H_reduced[component] if material in reduced_set
            else H_total[component]
            for material in mesh.GetMaterials()
        }
        h_components.append(mesh.MaterialCF(values))
    H_cf = CoefficientFunction(tuple(h_components))
    return {
        "solution": solution,
        "phi_reduced": phi_reduced_gf,
        "phi_total": phi_total_gf,
        "interface_multiplier": multiplier_gf,
        "kelvin_interface_multiplier": kelvin_multiplier_gf,
        "fes": fes,
        "fes_reduced": fes_reduced,
        "fes_total": fes_total,
        "mu_cf": mu_cf,
        "H_s": H_s,
        "source_potential": source_potential,
        "kelvin_source_potential": kelvin_source_potential,
        "total_source_h": total_source_h,
        "total_source_materials": total_source_materials,
        "H_cf": H_cf,
        "B_cf": mu_cf * H_cf,
    }


def solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
        mesh, H_s, source_potential, R_K, offset, *, bh_table,
        nonlinear_materials, reduced_materials, total_materials,
        interface_boundary, order=1, dirichlet_bbbnd="GND",
        bonus_intorder=4, kelvin_mats=("kelvin",), inverse="pardiso",
        mu_r_initial=1000.0, tolerance=2.0e-5, max_iterations=80,
        relaxation=0.3, interface_constraint_scale=None,
        kelvin_interface_boundary=None, kelvin_source_potential=None,
        total_source_h=None, total_source_materials=()):
    """Picard solve for the mixed total/reduced Omega formulation.

    The source split and its interface trace stay fixed throughout the
    iteration.  Only the permeability in ``nonlinear_materials`` is updated
    from the shared monotone ``B(H)`` law.  This is deliberately separate from
    hysteretic state evolution: the latter needs its own committed material
    history and is not silently approximated by a memoryless Picard update.

    The result has the same field keys as the linear mixed solve plus
    ``nonlinear_stats``.  Caller wraps the complete operation in
    :class:`ngsolve.TaskManager`.
    """
    from ngsolve import GridFunction, L2, VOL
    from radia.scalar_potential_solver import _build_bh_interpolator

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
    if not math.isfinite(mu_r_initial) or mu_r_initial <= 0.0:
        raise ValueError("mu_r_initial must be positive and finite")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive and finite")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be positive")
    if not 0.0 < relaxation <= 1.0:
        raise ValueError("relaxation must lie in (0, 1]")

    bh_array = np.asarray(bh_table, dtype=float)
    if bh_array.ndim != 2 or bh_array.shape[1] < 2:
        raise ValueError("bh_table must contain [H, B] rows")
    B_of_H = _build_bh_interpolator(bh_array[:, :2])
    B_scale = max(float(bh_array[:, 1].max()), MU_0)

    mu_elements = GridFunction(L2(mesh, order=0), name="mixed_omega_mu")
    nonlinear_elements = []
    for element in mesh.Elements(VOL):
        material = str(element.mat)
        if material in nonlinear_set:
            mu_elements.vec[element.nr] = MU_0 * float(mu_r_initial)
            coordinates = [mesh.vertices[vertex.nr].point for vertex in element.vertices]
            centroid = tuple(
                sum(float(point[component]) for point in coordinates) / len(coordinates)
                for component in range(mesh.dim))
            nonlinear_elements.append((element.nr, centroid))
        else:
            mu_elements.vec[element.nr] = MU_0

    kelvin_mu = make_kelvin_mu_cf(
        mesh, R_K, offset, kelvin_mats=kelvin_mats, mu_r_by_material={})

    def mixed_mu_cf():
        values = {}
        for material in mesh.GetMaterials():
            if material in nonlinear_set:
                values[material] = mu_elements
            elif any(key in material.lower() for key in kelvin_mats):
                values[material] = kelvin_mu
            else:
                values[material] = MU_0
        return mesh.MaterialCF(values, default=MU_0)

    B_previous = np.zeros(len(nonlinear_elements))
    result = None
    converged = False
    relative_change = float("inf")
    for iteration in range(1, int(max_iterations) + 1):
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, H_s, source_potential, R_K, offset,
            mu_r_by_material=None, reduced_materials=reduced_materials,
            total_materials=total_materials, interface_boundary=interface_boundary,
            order=order, dirichlet_bbbnd=dirichlet_bbbnd,
            bonus_intorder=bonus_intorder, kelvin_mats=kelvin_mats,
            inverse=inverse, interface_constraint_scale=interface_constraint_scale,
            mu_cf=mixed_mu_cf(),
            kelvin_interface_boundary=kelvin_interface_boundary,
            kelvin_source_potential=kelvin_source_potential,
            total_source_h=total_source_h,
            total_source_materials=total_source_materials)
        B_current = np.zeros(len(nonlinear_elements))
        for index, (element_nr, centroid) in enumerate(nonlinear_elements):
            H_value = result["H_cf"](mesh(*centroid))
            H_magnitude = math.sqrt(sum(float(value) ** 2 for value in H_value))
            B_current[index] = B_of_H(H_magnitude)
            if H_magnitude <= 1.0e-12:
                mu_r_next = float(mu_r_initial)
            else:
                mu_r_next = max(1.0, B_current[index] / (MU_0 * H_magnitude))
            if iteration > 1:
                mu_r_previous = float(mu_elements.vec[element_nr]) / MU_0
                mu_r_next = (
                    relaxation * mu_r_next + (1.0 - relaxation) * mu_r_previous)
            mu_elements.vec[element_nr] = MU_0 * mu_r_next

        if iteration > 1:
            relative_change = float(
                np.max(np.abs(B_current - B_previous))
                / max(B_scale, 1.0e-30))
            if relative_change <= tolerance:
                converged = True
                break
        B_previous = B_current

    if result is None:  # pragma: no cover - guarded by max_iterations validation
        raise RuntimeError("mixed total/reduced Omega Picard iteration did not start")
    result["mu_cf"] = mixed_mu_cf()
    result["B_cf"] = result["mu_cf"] * result["H_cf"]
    result["nonlinear_stats"] = {
        "method": "Picard",
        "iterations": iteration,
        "converged": converged,
        "relative_B_change": relative_change,
        "tolerance": float(tolerance),
        "relaxation": float(relaxation),
    }
    if not converged:
        raise RuntimeError(
            "mixed total/reduced Omega Picard iteration did not converge: "
            f"iterations={iteration}, relative_B_change={relative_change:.3e}, "
            f"tolerance={tolerance:.3e}")
    return result


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
