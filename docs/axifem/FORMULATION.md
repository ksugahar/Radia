# axifem — Formulation

A self-contained derivation of the Henrotte / Meeker axisymmetric finite-
element formulation that `radia.axifem` implements.  Read this
before reading the C++ source under `src/ext/axifem/`.  Companion doc
[`AXIFEM.md`](AXIFEM.md) covers usage, API, and validation results.
The executed P1/Q1/P2/Q2/P2-curved/Q2-curved proof artifact is
[`AXIFEM_ELEMENT_EVIDENCE.ipynb`](../../validation_test/axifem/AXIFEM_ELEMENT_EVIDENCE.ipynb).

## 1. Problem class

Axisymmetric magnetostatic + linear quasi-static eddy current.  Cylindrical
coordinates `(r, φ, z)`, all fields independent of `φ`.  Source current and
magnetisation are circumferential only:

$$
\mathbf{J} = J_\varphi(r, z)\,\hat\varphi,
\qquad
\mathbf{M} = M_\varphi(r, z)\,\hat\varphi.
$$

The induced eddy current is also circumferential under axial-symmetric
excitation, so the magnetic vector potential reduces to a single scalar
component:

$$
\mathbf{A}(r, z) = A_\varphi(r, z)\,\hat\varphi.
$$

## 2. Magnetic field from the vector potential

In cylindrical coordinates the curl of a φ-only vector field is

$$
\mathbf{B} = \nabla \times \mathbf{A}
= -\frac{\partial A_\varphi}{\partial z}\,\hat r
+ \frac{1}{r}\frac{\partial (r A_\varphi)}{\partial r}\,\hat z.
$$

Thus the radial component carries an explicit `∂/∂z` and the axial
component carries an `(1/r) ∂/∂r` of `r A_φ`.  The latter is the source
of the axis (`r = 0`) singularity that motivates the Henrotte basis
(see § 6).

It is convenient to introduce the **flux function**

$$
\psi(r, z) \equiv 2\pi r A_\varphi(r, z).
$$

Then `2π·∫B_z dr = ∂ψ/∂r` (a useful relationship for visualising flux
tubes), and the field becomes

$$
B_z = \frac{1}{2\pi r}\frac{\partial \psi}{\partial r},
\qquad
B_r = -\frac{1}{2\pi r}\frac{\partial \psi}{\partial z}.
$$

Both components carry a `1/r` factor when expressed via `ψ`.  The
Henrotte trick is to choose a basis so that the polynomial part of `ψ`
contains an explicit `r²` factor near the axis, exactly cancelling the
`1/r`.

## 3. Maxwell equations and the energy functional

Magnetostatic constitutive law inside a region of permeability `μ` and
magnetisation `M`:

$$
\mathbf{H} = \frac{1}{\mu}(\mathbf{B} - \mu_0 \mathbf{M}),
\qquad
\nabla \times \mathbf{H} = \mathbf{J}.
$$

Combining with `B = ∇ × A` gives the vector elliptic equation

$$
\nabla \times \left(\frac{1}{\mu}\nabla \times \mathbf{A}\right)
= \mathbf{J} + \nabla \times \left(\frac{\mu_0 \mathbf{M}}{\mu}\right).
$$

For axisymmetric `J = J_φ ĵ_φ` and `A = A_φ ĵ_φ` this reduces to a
scalar PDE on the meridian plane `Ω ⊂ {(r, z) : r ≥ 0}`:

$$
-\frac{\partial}{\partial r}\!\left[\frac{1}{\mu r}\frac{\partial(r A_\varphi)}{\partial r}\right]
-\frac{\partial}{\partial z}\!\left[\frac{1}{\mu}\frac{\partial A_\varphi}{\partial z}\right]
= J_\varphi - \text{div}\,(\mu_0 \mathbf{M}/\mu)\big|_\varphi.
$$

The corresponding energy functional (per radian of φ) on the meridian
plane is

$$
W[A_\varphi]
= \int_\Omega \left[\frac{1}{2\mu}|\mathbf{B}|^2 - J_\varphi A_\varphi\right]
  \cdot 2\pi r \, dr\, dz,
$$

with the Jacobian `2π r` from the φ-integration of the cylindrical
volume element.  Substituting `B_r, B_z` from § 2 and choosing the flux
function `ψ` as the unknown:

$$
W[\psi]
= \int_\Omega \frac{1}{8\pi^2 \mu r}
  \!\left[\!\left(\frac{\partial \psi}{\partial r}\right)^{\!2}
        + \!\left(\frac{\partial \psi}{\partial z}\right)^{\!2}\right]
  \cdot 2\pi r \, dr\, dz
- \int_\Omega \frac{J_\varphi \psi}{2\pi r} \cdot 2\pi r\, dr\, dz.
$$

The `r` in the Jacobian cancels one of the `1/r` factors, leaving a
weighted Dirichlet energy in `ψ`:

$$
W[\psi]
= \int_\Omega \frac{1}{4\pi \mu r}
  \!\left[\!\left(\frac{\partial \psi}{\partial r}\right)^{\!2}
        + \!\left(\frac{\partial \psi}{\partial z}\right)^{\!2}\right] dr\, dz
- \int_\Omega J_\varphi \psi \, dr\, dz.
$$

There remains one `1/r` factor multiplying the gradient term.  This is
the integrand singularity that any axisymmetric finite-element scheme
must cope with on axis-touching elements (`r → 0`).

## 4. Eddy currents add a σ-mass term

For a linear quasi-static eddy current with conductivity `σ(r, z)`,
Faraday's law gives `E = -∂A/∂t`, so the induced volume current is
`J_eddy = σ E = -σ ∂A_φ/∂t ĵ_φ`.  The corresponding term in the
weak form (multiplied by the test function and the `2π r` Jacobian) is

$$
\int_\Omega \sigma A_\varphi^{(\text{trial})} A_\varphi^{(\text{test})}
\cdot 2\pi r\, dr\, dz
= \int_\Omega \frac{\sigma}{2\pi r}\,\psi^{(\text{trial})} \psi^{(\text{test})}
   \cdot 2\pi r \, dr\, dz
= \int_\Omega \frac{\sigma}{2\pi r}\,\psi \psi'\, dr\, dz \cdot 2\pi
$$

Up to overall normalisation:

$$
M_\sigma[\psi, \psi'] = \frac{\sigma}{4\pi^2}
\int_\Omega \frac{\psi \psi'}{r^2}\, r\, dr\, dz
= \frac{\sigma}{4\pi^2}
\int_\Omega \frac{\psi \psi'}{r}\, dr\, dz.
$$

Or in the variable `s = r²` (see § 5), `dr = ds/(2r)`, so
`(1/r) dr dz = ds dz / (2 s)`:

$$
M_\sigma[\psi, \psi'] = \frac{\sigma}{8\pi^2}
\int_\Omega \frac{\psi \psi'}{s} \, ds\, dz.
$$

This is the σ-mass operator that `AxiHenrotteSigmaMassBFI` assembles.

## 5. Henrotte / Meeker change of variable: s = r²

Both the stiffness energy of § 3 and the σ-mass of § 4 contain a `1/r`
factor (or equivalently a `1/s` factor under `s = r²`).  Standard
finite-element shape functions polynomial in `(r, z)` with degree `p`
make the integrands

$$
\frac{(\partial_r \psi)^2}{r}\quad\text{and}\quad\frac{\psi^2}{r}
$$

singular as `r → 0`, which Gauss quadrature struggles to integrate
accurately.  In the limit, axis-touching elements lose effective
convergence rate from `p` to `p_eff = 1` even for `p ≥ 2`.

The Henrotte / Meeker observation [^Henrotte93][^MeekerFEMM]: choose
shape functions polynomial in `(s, z)` with `s = r²` instead.  Then

$$
\partial_r \psi = 2 r \, \partial_s \psi,
\qquad
(\partial_r \psi)^2 / r = 4 r \, (\partial_s \psi)^2 = 4 s^{1/2} \, (\partial_s \psi)^2,
$$

but the stiffness integrand transforms cleanly:

$$
\frac{(\partial_r \psi)^2}{r}\, dr\, dz
= 4 r (\partial_s \psi)^2 \, \frac{ds}{2r} \, dz
= 2 (\partial_s \psi)^2 \, ds \, dz,
$$

— the `r` and `1/r` cancel exactly, leaving a polynomial integrand in
`(s, z)` that integrates in **closed form** via Mathematica.

Similarly

$$
\frac{(\partial_z \psi)^2}{r}\, dr\, dz
= \frac{(\partial_z \psi)^2}{r} \cdot \frac{ds}{2r}\, dz
= \frac{(\partial_z \psi)^2}{2 s} \, ds\, dz,
$$

which keeps a `1/s` factor.  This second piece must be handled with
care on axis elements (`s_a = 0`) — see § 7.

## 6. Variational form on a single element

For axis-aligned quads, let `Ω_e = [r_a, r_b] × [z_a, z_b]` be one
rectangle, with `s_a = r_a²`, `s_b = r_b²`.  Choose `ψ` to be polynomial
in `(s, z)` of total degree `p` (in NGSolve API:
`H1Henrotte(mesh, order=p)`):

$$
\psi(s, z) = \sum_{m=1}^{N_p} c_m \, m(s, z),
$$

where the monomial set is

| `p` | `N_p` | Monomials `m(s, z)`                                  |
|----:|------:|------------------------------------------------------|
| 1   | 4     | `{1, s, z, s·z}`                                     |
| 2   | 9     | `{1, s, s², z, s·z, s²·z, z², s·z², s²·z²}`          |

The element stiffness in **monomial coefficients** is the symmetric
`N_p × N_p` matrix

$$
K_\phi^{(\text{mono})}_{m,n}
= \int_{s_a}^{s_b}\!\!\int_{z_a}^{z_b}
  \!\left[\frac{2}{\mu_z}\,\partial_s m \cdot \partial_s n
       + \frac{1}{2\mu_r s}\,\partial_z m \cdot \partial_z n\right]\!ds\, dz,
$$

where `μ_r, μ_z` are the (anisotropic) per-element permeabilities.  The
`1/s` in the second term is the residual `1/r` after the variable
change.  Both integrals are polynomial in `s` (for the first) and
admit a `log(s_b/s_a)` or `log(s_b)` closed form (for the second);
all entries are derived once symbolically and emitted as C++ code in
`q2_henrotte_generated.hpp`; the generated values are locked by
`validation_test/axifem/research/validate_q2_codegen.py`.

The σ-mass element matrix in monomial coefficients is

$$
M_{\sigma,\phi}^{(\text{mono})}_{m,n}
= \frac{\sigma}{4\pi^2}
  \int_{s_a}^{s_b}\!\!\int_{z_a}^{z_b}
    \frac{m(s,z)\, n(s,z)}{s}\, ds\, dz.
$$

Same `1/s` issue, same closed-form treatment.

### Conversion from monomial coefficients to nodal DOFs

The DOFs `V_j = A_φ` at each node form the API surface.  The monomial
coefficients are recovered from the DOFs via a per-element Vandermonde
inversion in the chosen monomial basis.  The conversion respects the
flux-function transformation `ψ = 2π r A_φ`:

$$
K_V = T \cdot V^{-T} \cdot K_\phi^{(\text{mono})} \cdot V^{-1} \cdot T,
\qquad
T = \mathrm{diag}(2\pi\, r_{\text{node}}),
$$

so that the assembled global system is in `A_φ`-DOFs at each mesh
vertex / edge midnode / face centre.  This conversion is done inside the
finite elements and the BFI integrators, not by the user.

### One field convention for every evaluator

Every element exposes one basis, `AxiHenrotteBaseFE::CalcVdof`, and every
consumer uses it: the native stiffness and sigma-mass integrators, the
volume value (`u`, `GridFunction` evaluation, `Set`), the physical
gradient (`grad(u)` = `(∂A_φ/∂r, ∂A_φ/∂z)`) and the boundary trace
(`u.Trace()`).  The basis is the `A_φ` basis

$$
N_i(r, z) = \frac{r_i\, \psi_i(s, z)}{r}, \qquad s = r^2,
$$

where `ψ_i` is the element's Lagrange basis in `(s, z)` and `r_i` the
radius of node `i`, so `N_i(node_j) = δ_ij` and the DOFs are nodal values
of `A_φ`.  The stiffness and mass of § 6 are then exactly the symbolic
forms

$$
K = \frac{2\pi}{\mu} \int \big(\partial_z u\, \partial_z v
  + (\partial_r u + u/r)(\partial_r v + v/r)\big)\, r\, dr\, dz,
\qquad
M = 2\pi\sigma \int u\, v\, r\, dr\, dz,
$$

evaluated with the same `N_i`; `tests/axifem/test_vdof_consistency.py`
assembles both sides on axis-aligned, general and oblique-edge meshes
and requires agreement to `1e-6`.  Near the axis
(`s <= 1e-8 s_scale`) the basis is evaluated in its limit form
`N = r_i ψ_s r`, `∂N/∂r = r_i ψ_s`, so `A_φ` and `∂A_φ/∂r` stay finite
and `A_φ = 0` on the axis.

**Edge midnodes.**  One rule places every edge midnode, in the volume
element and in the boundary trace: an edge that is horizontal
(`z` constant) or has exactly one end on the axis takes its midnode at the
`s`-midpoint `s = (s_a + s_b)/2`; every other edge takes the mesh
midpoint.  The DOF of an axis-edge midnode is therefore `A_φ` at
`r = sqrt((r_a² + r_b²)/2)`, not at the chord midpoint.  On an oblique
axis edge the `z` of that node moves along the edge so that the node
stays on the straight `(r, z)` edge.  The boundary trace
(`AxiHenrotteFE_Edge`) evaluates `N_i` along the edge with the same rule,
so the trace equals the volume value on the edge
(`test_vdof_consistency.py`, `test_edge_trace_conformity.py`).

**Admission.**  The s-isoparametric elements accept straight elements
only (`mesh.Curve(1)`): every midnode must lie on its chord within
`1e-10` of the edge length and the quad centre at the corner average.
The `(s, z)` map must keep the mesh orientation at every quadrature
point and node; a folded or degenerate map raises.  The axis-aligned
closed-form quads additionally require the vertex order
`(r_min, z_min), (r_max, z_min), (r_max, z_max), (r_min, z_max)`.

## 7. Axis-touching elements (`s_a = 0`)

When `r_a = 0`, the bottom edge of the element sits on the axis.  The
`1/s` factor in the second stiffness term and the σ-mass diverges as
`s → 0`.  The Henrotte basis avoids this by **dropping the
axis-incompatible monomials** that are non-zero at `s = 0`:

| `p` | Interior monomials | **Axis monomials** (axis-restricted basis)        |
|----:|--------------------|---------------------------------------------------|
| 1   | `{1, s, z, s·z}` (4) | `{s, s·z}` (2 — drops `1, z`)                     |
| 2   | 9 (above)          | `{s, s², s·z, s²·z, s·z², s²·z²}` (6 — drops the 3 monomials with no `s` factor) |

Geometrically, the dropped axis-side DOFs are forced to zero — they
correspond to `A_φ` not having an `s = const` mode at the axis (which
would imply `B_z` constant down to the axis with non-zero r component,
inconsistent with the assumption `A_φ ∝ r` near axis).  After
removing those DOFs the reduced `M_σ,φ` integrand has factors of `s` in
both `m` and `n` that cancel the `1/s`, restoring closed-form
integrability.

`AxiHenrotteFE_Q1_AxisAligned` and `AxiHenrotteFE_Q2_AxisAligned`
internally toggle between the interior basis (4 / 9 DOFs) and the
axis basis (2 / 6 DOFs) based on whether `s_a < ε`.  The dropped DOFs
are treated as a Dirichlet zero in NGSolve's free-DOF bookkeeping.

## 8. Triangle paths: P1 and P2

For unstructured triangle meshes, `axifem` provides both P1 and P2 triangle
elements:

- `AxiHenrotteFE_P1_Triangle` (3 DOFs/cell) uses the `{1, r², z}` basis.
  The per-element Python reference under `tests/axifem/_reference_python/`
  preserves the FEMM `StaticAxisymmetric()` P1 formula lineage, including
  the axis-touching cases (1- or 2-vertex on axis).  The production C++
  path is not a line-for-line FEMM port: it is reached through
  `Mesh("model.vol")` / NGSolve and assembles the V-DOF stiffness lane
  used by `H1Henrotte`, with regression gates for the `(s, z)`-triangle
  quadrature and degenerate `(s, z)` triangles (see the last bullet).
- `AxiHenrotteFE_P2_Triangle` (6 DOFs/cell, 3 vertices plus 3 edge
  midnodes) is s-isoparametric: the element geometry is the quadratic
  Lagrange map from the reference triangle to `(s = r², z)` through the
  nodes, and the basis is the reference Lagrange basis.  The trace on a
  shared edge then depends only on that edge's three nodes, so the space is
  conforming on oblique edges.  The former basis, physical monomials
  `{1, s, z, s², s z, z²}` fitted at the nodes, was not: on an edge oblique
  in `(r, z)`, `s` is quadratic in the edge parameter and the monomial trace
  is quartic, which three edge DOFs cannot fix.  Corner coordinates come
  from NGSolve's element transformation and the midnodes follow the
  edge-midnode rule of § 6 (mesh midpoint, or `s`-midpoint on horizontal
  and axis edges); a straight edge is reproduced exactly because `s` is
  quadratic and `z` linear in the edge parameter, and the `s`-midpoint on
  axis edges keeps the `(s, z)` map regular and `ψ ~ s` at the axis.
  Curved elements (`mesh.Curve(p)`, `p > 1`) are refused.  Stiffness and
  sigma mass are integrated over the
  reference element through the `(s, z)` map
  (`∫ ∇ψ·∇ψ / r dr dz = ∫ (2 ψ_s² + ψ_z² / (2 s)) ds dz`).
- `AxiHenrotteFE_P1_Triangle` uses the linear Lagrange basis on the same
  quadratic `(s, z)` geometry as P2 (the six geometry nodes are sampled from
  the element transformation and placed by the § 6 rule).  The `(s, z)`
  triangle then follows the straight `(r, z)` triangle exactly, so the space
  is conforming on oblique edges and its trace on an edge depends only on
  that edge.  The price is that `ψ` is linear in the reference coordinates,
  not in `s`: P1 reproduces `A_φ = 1/r` and `z/r` exactly but no longer the
  uniform field `ψ = s`.  P1 stays out of the order-convergence gate until
  it is requalified.

P1 triangles have lower per-DOF accuracy than Q2 quads
(`disk_convergence` example shows ~5–6 % gap to BEM-Foster reference even
at very fine meshes).  P2 triangles serve unstructured regions (spheres,
Kelvin half-discs, OCC-generated regions) on straight meshes; Q2 straight
quads remain the preferred structured-mesh path for disks, cylinders, and
rectangular workpieces.

### General (non-axis-aligned) Q2 quads are opt-in production C++

The default `AxiHenrotteFE_Q2_AxisAligned` is still the straight,
axis-aligned closed-form quad.  It constructs the element from the four
corner coordinates and is the preferred path for structured rectangular
workpieces.

The 9-node general Q2 quad is available as an opt-in production path:
`H1Henrotte(mesh, order=2, curvedquad=True)` selects
`AxiHenrotteFE_Q2_Curved`.  This element accepts any straight quad
(skewed, trapezoidal, any vertex order), consumes the corner positions from
the element transformation and is s-isoparametric in the same way as the P2
triangle (biquadratic Lagrange map to `(s, z)`, reference Lagrange basis, the
same axis-edge midnode rule).  It is conforming with itself and with P2
triangles on oblique edges; the former physical-monomial version was not and
drifted under boundary-layer refinement.  `tests/axifem/test_q2_curved.py`
checks that curved Q2 and the axis-aligned closed form (a different space,
with `s`-midpoint edge nodes) converge together on straight quads and that
curved Q2 converges on skewed annular quads;
`tests/axifem/test_edge_trace_conformity.py` checks the shared-edge traces.

## 9. Boundary conditions

Standard FEM boundary conditions in cylindrical coordinates:

- **Axis (`r = 0`)**: `A_φ = 0` is mathematically required (the
  vector potential has no φ-component on the symmetry axis).  This is
  imposed by Dirichlet on the axis edge label.
- **Far-field truncation (`r → ∞`)**: `A_φ → 0`, imposed by Dirichlet
  on the outer edges.  For applications needing a true open-domain
  treatment (no truncation), use Kelvin transformation (see
  `docs/kelvin/`) instead.
- **Symmetry plane (`z = 0` for an even-parity problem)**: `∂A_φ/∂z = 0`,
  Neumann (no boundary integral).

## 10. Summary of the assembled system

The complete weak form on the meridian plane Ω is

$$
\underbrace{K_V \cdot \mathbf{V}}_{\text{stiffness}}
+ \underbrace{M_\sigma \cdot \frac{d\mathbf{V}}{dt}}_{\text{eddy current}}
= \underbrace{\mathbf{F}}_{\text{source}},
$$

where `V` is the global vector of nodal `A_φ` DOFs.  In the API:

- `AxiHenrotteStiffnessBFI(mu_cf)` assembles `K_V` (per § 6, with
  `μ` from the CoefficientFunction).
- `AxiHenrotteSigmaMassBFI(sigma_cf)` assembles `M_σ` (per § 4 + 6,
  with `σ` from the CoefficientFunction).
- `H1Henrotte(mesh, order=p)` provides the appropriate FESpace
  (P1 / P2 triangle, Q1 / Q2 quad, with axis-element bookkeeping).
- The user assembles `K`, `M` with standard NGSolve `BilinearForm`
  and feeds them into any solver (direct, iterative, eigenvalue,
  Hiruma 3-term recurrence).

## 10b. Temperature uses standard H1

The Henrotte space is the `A_φ` basis `N_i = r_i ψ_i / r` and carries no
temperature integrators: a scalar temperature is not of that form.
Axisymmetric heat
([`src/radia/panels/calc_heat_axisym.py`](../../src/radia/panels/calc_heat_axisym.py))
uses standard NGSolve `H1(order=2)` with the `2 pi r` weight, matching the
FEMM 4.2 heat solver.  Order 2 is the default because a piecewise-linear
radial profile cannot represent `dT/dr = 0` at `r = 0`; in the committed 9x9
uniform-flux cylinder validation it reduces the spurious radial axis slope
from about 139 C/m to below 0.3 C/m.  The former Henrotte heat BFIs used the
retired `ψ`-on-mesh basis and had no production consumer; they were removed
when the field convention was unified.

## 10c. Boundary trace + Neumann RHS

For boundary terms `LinearForm += q * v.Trace() * weight * ds(label)`
(current sheet, Neumann data, DtN couplings) the space provides the
boundary evaluator `AxiHenrotteDiffOpIdBnd`, wired into `evaluator[BND]`,
and the 1D element `AxiHenrotteFE_Edge(order, r0, z0, r1, z1)`.  The edge
element evaluates the same `A_φ` basis `N_i` as the volume element on that
edge: vertex DOFs at the edge ends and, at order 2, the midnode placed by
the § 6 rule.  Its local parameter follows NGSolve's segment convention
(`t = 1` at vertex 0).  `tests/axifem/test_vdof_consistency.py` checks that
the trace equals the volume value on horizontal, vertical, oblique and axis
edges.

### NGSolve Trace() requirement

The Henrotte FESpace requires the explicit form
`LinearForm += q * v.Trace() * weight * ds(label)`, also in boundary
bilinear terms such as `h * u.Trace() * v.Trace() * ds(label)`, because
the proxy-function machinery does not auto-inject the custom BND DiffOp at
the implicit-trace site.  The error
`NgException: Testfunction does not support BND-forms, maybe a Trace()
operator is missing` means the `.Trace()` is missing.

### Closed-form quads need structured axis-aligned meshes

The Q1 / Q2 closed-form quads accept only rectangles in `(r, z)` with the
vertex order of § 6.  Netgen's `quad_dominated=True` mesher produces
slightly skewed quads; use `curvedquad=True` for those, or a structured
grid generator (see
[`validation_test/panels/fixtures/generate_heat_cylinder_axisym.py`](../../validation_test/panels/fixtures/generate_heat_cylinder_axisym.py)).

## 11. Cross-validation

The closed-form C++ implementation is cross-checked against two
independent paths:

1. **Pure-Python reference** — `tests/axifem/_reference_python/`
   (axifem_core for P1 triangle, axifem_quad for Q1 quad,
   axifem_quad_q2 for Q2 quad with Gauss 8×8 numerical quadrature).
   Test `tests/axifem/test_python_reference_consistency.py` asserts
   the C++ stiffness eigenvalues match the Python references and
   smoke-tests the shipped P2 curved triangle assembly path.  The
   historical curved-Q2 Python prototype remains under
   `validation_test/axifem/research/prototypes/axifem_quad_q2_curved.py`;
   the production C++ curved-Q2 gate is
   `tests/axifem/test_q2_curved.py`.

2. **Mathematica derivation** —
   `validation_test/axifem/research/validate_q2_codegen.py` runs the
   closed-form C++ vs the same monomial integrals re-derived in
   sympy at 50-digit precision.  Should match to ~1e-13 relative.

## 12. Why this matters for induction heating

Induction heating coils + axisymmetric workpieces (cylinders, disks,
tubes) generate **circumferential eddy currents** — exactly the case
where `A_φ` is the only nontrivial vector-potential component.  The
axifem Henrotte / Meeker formulation:

- removes the standard P1 axis singularity in `B_r ∝ 1/r`,
- gives smooth per-element `B_z = const`, `B_r ∝ 1/r` — well suited to
  eddy-current time-constant extraction (the inductive transient of
  the workpiece is driven by these τ values),
- reduces the problem from 3-D to 2-D `(r, z)` without losing
  axisymmetric physics,
- runs at FEMM-grade accuracy per DOF, with current evidence recorded in
  [`AXIFEM_ELEMENT_EVIDENCE.ipynb`](../../validation_test/axifem/AXIFEM_ELEMENT_EVIDENCE.ipynb) and
  [`validation_test/axifem/axifem_element_evidence.json`](../../validation_test/axifem/axifem_element_evidence.json).

This is why `radia.axifem` is a critical module of the radia
distribution and not a separate research package.

## References

[^Henrotte93]: F. Henrotte, B. Meys, A. Genon, W. Legros, "A new method
    for axisymmetric linear and nonlinear problems," *IEEE Transactions
    on Magnetics* **29**(2):1352–1355, March 1993.
    DOI [10.1109/20.250664](https://doi.org/10.1109/20.250664).

[^MeekerFEMM]: D. Meeker, *Finite Element Method Magnetics 4.2 — User's
    Manual + axisymmetric formulation notes*.  Reference C++ source:
    [`prob3big.cpp::StaticAxisymmetric()`](https://www.femm.info/wiki/Documentation).
