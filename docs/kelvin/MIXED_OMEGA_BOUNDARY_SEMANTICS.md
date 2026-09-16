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

Nonlinear constitutive laws do not remove this distinction. Compare physical
H and B under the same boundary, source lift, material law and gauge convention,
not the raw total and reduced potential values. Those potentials generally
differ by the source lift. If a harmonic source remainder is retained in the
total region, reconstruct H with that remainder as well.
