# Fixed-source right-hand sides in mixed scalar potential solves

In a fixed-mesh magnetostatic Picard solve, the source field and interface
potential jump do not change. Re-evaluating those fields at quadrature points
on every iteration is unnecessary.

The order-one material path now separates:

- the fixed reduced-region, exterior and interface contributions;
- the total-region contribution weighted by the changing permeability.

When the total source acts only in nonlinear regions, its contribution is
assembled once as a rectangular operator from the piecewise-constant material
space to the mixed test space. Each iteration applies this operator to the
updated permeability coefficients. This retains the original quadrature and
does not replace the source field with a spatial interpolation.

The cache is private to one Picard call. It is never shared across meshes,
source changes, or runs. `cache_fixed_rhs=False` retains full reassembly for
regression comparisons. `nonlinear_stats['rhs_reuse']` identifies the path.
The projected higher-order material path is unchanged.

## Finite Exterior Boundaries

The default natural boundary specifies zero **total** normal flux.
`reduced_zero_normal_boundary` instead specifies zero normal **correction**
field, so that total normal H equals source normal H. These are different
physical conditions. The latter requires the fixed source-normal boundary
term in the weak right-hand side. Only named exterior faces of reduced
regions are accepted; internal interfaces are rejected.

A uniform-source control distinguishes these conditions analytically.
Cached/uncached nonlinear controls, both with and without a total-region
source, protect the discrete field and iteration-count equivalence.
