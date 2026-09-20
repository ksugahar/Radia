# Coil source precision

Separate direct source-field accuracy from finite-element response accuracy.
Start with a finite rectangular-section circular coil in vacuum and compare
its axial vector field to the closed-form volume-current integral. A filament
formula is not the reference for a finite conductor cross-section.

For inner/outer radii `a,b`, axial half-height `h`, uniform azimuthal current
density `J`, and axial observation coordinate `z`, the axial field is
`B(z) = mu0*J/2 * [G(z+h)-G(z-h)]`, where
`G(u) = u*(asinh(b/abs(u))-asinh(a/abs(u)))` and `G(0)=0`.
This integrates the circular-filament field over the radial and axial section;
the sign follows the winding direction. Subtraction can lose relative digits
far from the coil, so a double-precision formula is not an unlimited-precision
oracle.

The legacy `FldLenRndSw` name controls small length perturbations, not
rendering. Current repair offsets are deterministic. For precision tests at
nonsingular points, call `FldLenRndSw("off")` before constructing the coil
and retain that setting through field evaluation. Record this policy with
geometry, current, units, implementation version and observation points.
Switching it off only after construction cannot undo perturbed geometry.

The global default is unchanged. Disabling repair is not blanket validation
of conductor-boundary or singular evaluations. Do not mutate global precision
settings concurrently with other calculations.

`tests/test_coil_axis_closed_form.py` covers current reversal and three length
scales with a finite-section analytic reference. Its tolerance is a regression
gate, not a claim about arbitrary off-axis coils, material interfaces, nonlinear
response, or every solver observable.
