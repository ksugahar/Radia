"""Maxwell-stress surface force on an NGSolve boundary region.

``radia.force`` turns *sampled* fields into force; this module is its symbolic
twin for callers that already hold an NGSolve ``CoefficientFunction`` and want
the boundary integral evaluated by NGSolve's own quadrature.  The stress
identity is the same one, so the two routes agree by construction instead of
being retyped per call site.

Following the lab TaskManager policy these are helpers: they never open a
``with TaskManager():`` region.  The caller owns that context.
"""

from __future__ import annotations

import math

from ngsolve import (BND, BoundaryFromVolumeCF, CoefficientFunction, Conj,
                     Integrate, dx, specialcf)
from ngsolve import sqrt as ng_sqrt
from ngsolve import x as ng_x
from ngsolve import y as ng_y

MU0 = 4.0e-7 * math.pi

__all__ = [
    "MU0",
    "air_gap_maxwell_torque_2d",
    "air_gap_maxwell_torque_arkkio_2d",
    "air_gap_maxwell_torque_line_2d",
    "maxwell_surface_force",
    "time_average_maxwell_surface_force",
]


def _amplitude_factor(amplitude: str) -> float:
    key = str(amplitude).strip().lower()
    if key == "peak":
        return 0.5
    if key == "rms":
        return 1.0
    raise ValueError("amplitude must be 'peak' or 'rms'")


def _orientation_sign(normal_points_out_of_body: bool) -> float:
    """+1 when the mesh normal already points out of the body, -1 otherwise.

    ``specialcf.normal`` points out of the element being integrated.  When the
    body is a HOLE in the mesh -- the SIBC workpiece of ``calc_fem_kelvin`` is
    the standard case -- that normal points INTO the body, and the traction
    integral must be negated to give the force on the body.
    """

    return 1.0 if bool(normal_points_out_of_body) else -1.0


def maxwell_surface_force(
    magnetic_flux_density,
    mesh,
    boundary,
    *,
    normal_points_out_of_body=True,
    permeability_H_per_m=MU0,
):
    """Static Maxwell-stress force in N over a closed boundary region in air.

        F_k = oint_S (1/mu) [ (B.n) B_k - 0.5 |B|^2 n_k ] dS

    ``magnetic_flux_density`` is a real three-component CoefficientFunction,
    typically ``curl(gfu)``.  ``boundary`` is an NGSolve boundary region, e.g.
    ``mesh.Boundaries("sibc")``.  Returns ``[Fx, Fy, Fz]``.
    """

    if permeability_H_per_m <= 0.0:
        raise ValueError("permeability_H_per_m must be > 0")
    sign = _orientation_sign(normal_points_out_of_body)
    normal = specialcf.normal(3)
    field = magnetic_flux_density
    normal_field = sum(field[k] * normal[k] for k in range(3))
    squared = sum(field[k] * field[k] for k in range(3))
    inverse_mu = 1.0 / permeability_H_per_m
    return [
        sign
        * float(
            Integrate(
                inverse_mu * field[k] * normal_field
                - 0.5 * inverse_mu * squared * normal[k],
                mesh,
                BND,
                definedon=boundary,
            )
        )
        for k in range(3)
    ]


def time_average_maxwell_surface_force(
    magnetic_flux_density_phasor,
    mesh,
    boundary,
    *,
    normal_points_out_of_body=True,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Cycle-averaged Maxwell-stress force in N from a complex phasor B.

        <F_k> = oint_S (f/mu) [ Re(B_k conj(B.n)) - 0.5 |B|^2 n_k ] dS

    with ``f = 0.5`` for peak phasors and ``f = 1`` for RMS phasors -- the
    time average of cos^2(wt).  For a field with real peak amplitude B0 this
    equals ``0.5 * maxwell_surface_force(B0)``.

    ``magnetic_flux_density_phasor`` is a complex three-component
    CoefficientFunction, typically ``curl(gfu)`` from a harmonic A-solve.
    Returns ``[Fx, Fy, Fz]``.
    """

    if permeability_H_per_m <= 0.0:
        raise ValueError("permeability_H_per_m must be > 0")
    factor = _amplitude_factor(amplitude)
    sign = _orientation_sign(normal_points_out_of_body)
    normal = specialcf.normal(3)
    field = magnetic_flux_density_phasor
    normal_field = sum(field[k] * normal[k] for k in range(3))
    squared = sum((field[k] * Conj(field[k])).real for k in range(3))
    scale = factor / permeability_H_per_m
    return [
        sign
        * float(
            Integrate(
                scale * (field[k] * Conj(normal_field)).real
                - 0.5 * scale * squared * normal[k],
                mesh,
                BND,
                definedon=boundary,
            ).real
        )
        for k in range(3)
    ]


def _polar_air_gap_components(magnetic_flux_density, center_xy):
    """Return (radius, B_radial, B_tangential) about ``center_xy`` in the plane."""

    try:
        center_x, center_y = (float(value) for value in center_xy)
    except (TypeError, ValueError) as exc:
        raise ValueError("center_xy must be two finite numbers") from exc
    dx_cf = ng_x - center_x
    dy_cf = ng_y - center_y
    radius = ng_sqrt(dx_cf * dx_cf + dy_cf * dy_cf)
    cos_phi = dx_cf / radius
    sin_phi = dy_cf / radius
    field_x = magnetic_flux_density[0]
    field_y = magnetic_flux_density[1]
    radial = field_x * cos_phi + field_y * sin_phi
    tangential = -field_x * sin_phi + field_y * cos_phi
    return radius, radial, tangential


def air_gap_maxwell_torque_line_2d(
    magnetic_flux_density,
    mesh,
    boundary,
    *,
    stack_length_m=1.0,
    center_xy=(0.0, 0.0),
    permeability_H_per_m=MU0,
):
    """Maxwell-stress torque in N m from a single air-gap contour.

        T = (L / mu) * oint_C r B_r B_phi dl

    ``boundary`` is a closed contour lying wholly in the air gap, e.g.
    ``mesh.Boundaries("airgap_mid")``.  ``magnetic_flux_density`` is the planar
    two-component field, ``CF((grad(A)[1], -grad(A)[0]))`` for an A_z solve.

    This is the classic single-contour Maxwell torque.  It is NOT Arkkio's
    method: Arkkio averages the same integrand over the radial thickness of the
    gap, which is what buys the reduced mesh sensitivity -- see
    :func:`air_gap_maxwell_torque_arkkio_2d`.
    """

    permeability = float(permeability_H_per_m)
    if not permeability > 0.0:
        raise ValueError("permeability_H_per_m must be > 0")
    length = float(stack_length_m)
    if not length > 0.0:
        raise ValueError("stack_length_m must be > 0")
    # The contour is normally an INTERNAL edge and the field is normally
    # grad() of a volume GridFunction.  Integrating such a CF directly over a
    # BND region evaluates it in the wrong space and returns a near-zero,
    # mesh-dependent number: measured 0.0116 -> 0.0007 against an exact 28.125
    # as the gap thinned.  BoundaryFromVolumeCF takes the volume trace and
    # recovers 28.110 -> 28.124.  It is a no-op for a purely analytic CF.
    radius, radial, tangential = _polar_air_gap_components(
        BoundaryFromVolumeCF(magnetic_flux_density), center_xy
    )
    value = Integrate(radius * radial * tangential, mesh, BND, definedon=boundary)
    return (length / permeability) * float(getattr(value, "real", value))


def air_gap_maxwell_torque_arkkio_2d(
    magnetic_flux_density,
    mesh,
    region,
    *,
    inner_radius_m,
    outer_radius_m,
    stack_length_m=1.0,
    center_xy=(0.0, 0.0),
    permeability_H_per_m=MU0,
):
    """Arkkio air-gap torque in N m, averaged over the gap thickness.

        T = L / (mu (r_out - r_in)) * int_S r B_r B_phi dS

    ``region`` is the meshed air-gap annulus, e.g. ``mesh.Materials("airgap")``,
    and ``inner_radius_m`` / ``outer_radius_m`` are its radial bounds.  Averaging
    over the annulus instead of reading one contour is the whole point of
    Arkkio's method (Arkkio 1987): it removes the sensitivity to which contour
    the mesh happens to place elements on.  As the annulus thins about a radius
    R this converges to :func:`air_gap_maxwell_torque_line_2d` at R.
    """

    permeability = float(permeability_H_per_m)
    if not permeability > 0.0:
        raise ValueError("permeability_H_per_m must be > 0")
    length = float(stack_length_m)
    if not length > 0.0:
        raise ValueError("stack_length_m must be > 0")
    inner = float(inner_radius_m)
    outer = float(outer_radius_m)
    if not 0.0 <= inner < outer:
        raise ValueError("require 0 <= inner_radius_m < outer_radius_m")
    radius, radial, tangential = _polar_air_gap_components(
        magnetic_flux_density, center_xy
    )
    value = Integrate(radius * radial * tangential, mesh, definedon=region)
    span = outer - inner
    return (length / (permeability * span)) * float(getattr(value, "real", value))


def air_gap_maxwell_torque_2d(
    magnetic_flux_density,
    mesh,
    domain,
    *,
    method="line",
    stack_length_m=1.0,
    center_xy=(0.0, 0.0),
    inner_radius_m=None,
    outer_radius_m=None,
    permeability_H_per_m=MU0,
):
    """Select between the single-contour and Arkkio air-gap torque routes.

    ``method="line"`` reads one contour (``domain`` is a boundary region);
    ``method="arkkio"`` averages over the gap annulus (``domain`` is a material
    region and both radii are required).  The two agree as the annulus thins,
    so a validation lane can quote both from one solve.
    """

    key = str(method).strip().lower()
    if key == "line":
        if inner_radius_m is not None or outer_radius_m is not None:
            raise ValueError("radii apply only to method='arkkio'")
        return air_gap_maxwell_torque_line_2d(
            magnetic_flux_density,
            mesh,
            domain,
            stack_length_m=stack_length_m,
            center_xy=center_xy,
            permeability_H_per_m=permeability_H_per_m,
        )
    if key == "arkkio":
        if inner_radius_m is None or outer_radius_m is None:
            raise ValueError(
                "method='arkkio' requires inner_radius_m and outer_radius_m"
            )
        return air_gap_maxwell_torque_arkkio_2d(
            magnetic_flux_density,
            mesh,
            domain,
            inner_radius_m=inner_radius_m,
            outer_radius_m=outer_radius_m,
            stack_length_m=stack_length_m,
            center_xy=center_xy,
            permeability_H_per_m=permeability_H_per_m,
        )
    raise ValueError("method must be 'line' or 'arkkio'")
