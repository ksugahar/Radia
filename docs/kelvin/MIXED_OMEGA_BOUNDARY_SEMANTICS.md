# Mixed potential boundary semantics

In a reduced region, `H = Hs - grad(phi_r)`. For an exact source lift
`Hs = -grad(phi_s)`, the physical total potential is `phi_t = phi_s + phi_r`,
up to a common gauge. Thus prescribed total potential `phi_t = g` requires
`phi_r = g - phi_s`, not `phi_r = g`. Conversely, `phi_r = 0` prescribes
`phi_t = phi_s`, not zero total potential. A globally single-valued source
lift must not be assumed across current-linked topology.

For isotropic permeability, zero normal correction field means
`d(phi_r)/dn = 0`; the total normal flux is then `mu * Hs.n` and need not
vanish. In the weak form `a(phi_r,v) = integral(mu*Hs.grad(v))`, this boundary
requires subtracting `integral_boundary(mu*Hs.n*v)` from the right-hand side.
Omitting that correction gives the natural zero-normal-total-flux condition.
Use `reduced_zero_normal_boundary` for the former condition. Boundary labels
alone do not establish identical physical conditions.

The low-level solver's legacy `total_dirichlet_cf` parameter sets only the
total-region BBBND point gauge. It does not impose a surface Dirichlet value.
Changing this gauge must shift the potentials without changing H or B.
General exterior-surface Dirichlet data are not exposed by this parameter;
do not claim they were tested through a nonzero point gauge.

Finite-domain linear and Picard solves accept `surface_dirichlet`, a mapping
from `reduced` / `total` to exact exterior boundary labels and scalar values.
For example, `{"reduced": {"outer": g - phi_s}}` imposes physical total
potential `g` when a valid source lift `phi_s` is known on that surface.
The caller supplies this lift; an iron-only source projection cannot silently
be extended to an exterior air boundary. A nonempty surface condition replaces
the point gauge. Simultaneous `total_dirichlet_cf`, internal/wrong-region labels,
and overlapping correction-Neumann labels are rejected. Kelvin-exterior plus
surface Dirichlet is explicitly unsupported.
Prescribing both sides where exterior faces meet an interface junction is
also rejected: that requires eliminating redundant multiplier constraints.

Regression covers uniform-field positive/negative controls, nonzero total
surface values, and nonlinear P1/P2 source-split invariance using full RHS
reassembly. The nonlinear comparison keeps identical physical H/B while
moving an exactly representable gradient between the source and potential.
This is a formulation test, not validation of arbitrary source projections,
external solver equivalence, or continuum accuracy.

Material-dependent RHS caching uses an explicit volume integration rule shared
by the linear form and rectangular material operator. Equal bonus orders alone
do not establish identical quadrature for those two form types. A smooth,
spatially varying source regression protects full-reassembly/cache parity;
constant sources alone would miss this error.

Nonlinear constitutive laws do not remove this distinction. Compare physical
H and B under the same boundary, source lift, material law and gauge convention,
not the raw total and reduced potential values. Those potentials generally
differ by the source lift. If a harmonic source remainder is retained in the
total region, reconstruct H with that remainder as well.

The total-source Hodge projection fixes its additive constant by setting one
actual mesh vertex DOF to zero in each connected projection component. It
does not add a mass penalty. `gauge_vertices` can select those zero-based
mesh vertex numbers explicitly; otherwise the smallest vertex number in
each component is used. The retired `gauge_epsilon` argument is no longer
accepted by `project_source_total_hodge`. A different gauge must preserve
the projected gradient and the final physical field. The strict true-residual
gate applies to the free rows. Returned gauge diagnostics also retain fixed-row
reactions and the full residual, so a flux compatibility defect cannot be
mistaken for a healthy unconstrained Neumann solve.

Surface-flux Hodge loads enforce Neumann compatibility on the actual assembled
H1 load, separately for every connected component, even when components share
a material name. The constant test is represented by its FE coefficient vector;
summing all hierarchical H1 coefficients does not represent a constant test.
A NumberSpace flux integral generally uses a different quadrature order and
cannot balance this load.

Before removing the component boundary mean, the absolute net flux divided by
the absolute boundary flux must be at most `RELATIVE_LIMIT` (1e-6). Both are
computed with the H1 load's assembly rule. Larger defects raise
`SurfaceFluxCompatibilityFailed`; they do not trigger an automatic volume-load
fallback. Resolve source quadrature or source closure explicitly. Diagnostics
retain raw/corrected compatibility, removed means, and fixed-row reactions.
The mixed-system source boundary loads use the same component compatibility
treatment in their scalar source H1 space. Interface constraints still leave
a common constant mode, so they do not justify absorbing raw flux defects at
the global point gauge. Prescribed physical flux boundary terms are separate
and are not changed by this source correction.
Compatibility, the tangential source gate, and physical-field accuracy are
independent requirements.
