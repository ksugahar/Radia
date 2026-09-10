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

from ngsolve import BND, Conj, Integrate, specialcf

MU0 = 4.0e-7 * math.pi

__all__ = [
    "MU0",
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
