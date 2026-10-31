# Candidate independent CPU reference

This reference is under investigation; it is not an adopted validation oracle.
The scope is the constant tetrahedral and affine triangular self-diagonal
charges, not the complete HDiv-MMM operator or solve. No accuracy threshold
has been approved for this study.

For a charge m, compute

`G = integral_T integral_T m(x) m(y) / |x-y| dy dx / (4*pi)`.

For each outer point x inside a triangle, partition the source triangle into
three triangles with vertex x. On one part, write
`y = x + u*((1-v)*a + v*b)`, where a and b point from x to the edge vertices.
The surface Jacobian is `u*|a cross b|` and the distance is
`u*|(1-v)*a+v*b|`. For affine m, integrating the cancelled radial variable
exactly gives `(m(x)+m(edge(v)))/2`. The remaining edge integral uses open
Gauss nodes. Evaluate `(1, eta, xi)` directly from reference coordinates;
the implementation imports no Radia potential, moment or gradient routine.

For constant tetrahedral charge, partition the source tetrahedron into four
tetrahedra with vertex x. On each, write
`y=x+u*(a+xi*(b-a)+eta*(c-a))`. The volume Jacobian is
`u^2*|det(a,b,c)|`. Cancellation of distance leaves the radial integral
`integral_0^1 u du = 1/2`. Integrate the remaining reference triangle with
an independently selected Gauss rule. No singular node is removed and no
distance cutoff is used.

The outer triangle uses `(xi,eta)=(s,(1-s)*t)` with Jacobian `1-s`.
The outer tetrahedron uses `(xi,eta,zeta)=(s,(1-s)*t,(1-s)*(1-t)*v)`
with Jacobian `(1-s)^2*(1-t)`. Physical Jacobians are applied separately.
Absolute determinants handle arbitrary vertex winding. These constructions
assume straight, nondegenerate simplices and interior outer Gauss nodes.

General background: M. G. Duffy, *Quadrature Over a Pyramid or Cube of
Integrands with a Singularity at a Vertex*, SIAM Journal on Numerical
Analysis 19(6), 1260–1262 (1982),
<https://epubs.siam.org/doi/10.1137/0719090>.
Only its publisher metadata/abstract was consulted. This citation motivates
singularity transformations; it does not verify this derivation or code.

Before adoption, independently check the transformations symbolically,
refine inner and outer orders separately, compare identical charge inputs
to the production CPU entry, obtain reciprocal review and a human-approved
criterion. Finite differences between two quadrature levels are diagnostics,
not a rigorous error bound. CPU/GPU parity against the same analytic kernel
does not independently establish its physical accuracy. P1 quadrature
refinement is not FE p-convergence.
