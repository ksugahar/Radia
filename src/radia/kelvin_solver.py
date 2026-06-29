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

from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                      Periodic, TaskManager, curl, dx, InnerProduct,
                      Conj, Integrate)

from radia.kelvin_material import make_kelvin_nu_cf, NU_0


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


def solve_full_A_kelvin(mesh, J_source_cf, R_K, offset,
                         source_material="coil",
                         nu_0=NU_0,
                         order=1,
                         dirichlet_bbnd="GND",
                         gauge_eps=1e-8,
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
                            gauge_eps=1e-8,
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
