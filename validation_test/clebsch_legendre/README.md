# Clebsch-Potential Magnetostatics: Legendre (Hodograph) Transformations to Flux Coordinates

Derivation of the magnetostatic equations with Clebsch potentials as
independent variables. Two transformations of the forward problem
(`B = grad(phi) x grad(psi)`, unknowns `phi(x,y,z)`, `psi(x,y,z)`):

1. **Full swap** `(x, y, z) -> (phi, psi, z)`: unknowns
   `x(phi,psi,z)`, `y(phi,psi,z)` — sections 1-5, verified by
   [`verify_clebsch_legendre_transform.py`](verify_clebsch_legendre_transform.py)
   (sympy, 27 checks).
2. **Partial swap** `(x, y, z) -> (x, y, psi)`: unknowns
   `z(x,y,psi)`, `phi(x,y,psi)` — section 6, verified by
   [`verify_clebsch_partial_legendre_transform.py`](verify_clebsch_partial_legendre_transform.py)
   (sympy, 26 checks). In the 2D limit `phi` reduces to the magnetic
   vector potential component `A_z`.

## 1. Forward formulation (independent variables x, y, z)

Introduce two Clebsch potentials and represent the flux density as

```
B = grad(phi) x grad(psi)
```

- `div B = 0` holds **identically** (no equation needed).
- Field lines are the intersections of the surfaces `phi = const` and
  `psi = const`; `dphi dpsi` is the flux element (`dPhi = dphi dpsi`).
- In a current-free region the governing equation is

```
curl( grad(phi) x grad(psi) ) = 0
```

Three components, of which only two are independent (the identity
`div(curl B) = 0` relates them) — consistent with the two unknowns
`phi`, `psi`.

## 2. The variable swap (Legendre / hodograph transformation)

Wherever `B_z != 0`, the triple `(phi, psi, z)` is a valid coordinate
system. Swap dependent and independent variables and take

```
x = x(phi, psi, z),    y = y(phi, psi, z)
```

as the new unknowns. With the Jacobian

```
J = x_phi y_psi - x_psi y_phi        (subscripts = partial derivatives)
```

the chain rule gives the exact relations

```
grad(phi) = ( y_psi, -x_psi, x_psi y_z - y_psi x_z ) / J
grad(psi) = ( -y_phi, x_phi, x_z y_phi - y_z x_phi ) / J

B = grad(phi) x grad(psi) = (1/J) ( x_z, y_z, 1 ),     B_z = 1/J
```

so the coordinate lines `phi, psi = const` ARE the field lines:
`x_z`, `y_z` are the field-line slopes and `1/J` is the flux density.
The transformation is admissible exactly where `J` is finite and
nonzero, i.e. where the field has a nonvanishing z-component.

## 3. Transformed magnetostatic equations

The covariant components of `B` (`b_i = B . dr/du^i` with
`u = (phi, psi, z)`, `r = (x, y, z)`) are

```
b_phi = ( x_phi x_z + y_phi y_z ) / J
b_psi = ( x_psi x_z + y_psi y_z ) / J
b_z   = ( 1 + x_z^2 + y_z^2 ) / J
```

`curl B = 0` is equivalent to the vanishing of the curl of the
covariant components, giving the **transformed governing equations**

```
(E1)  d/dphi [ (x_psi x_z + y_psi y_z)/J ] = d/dpsi [ (x_phi x_z + y_phi y_z)/J ]
(E2)  d/dphi [ (1 + x_z^2 + y_z^2)/J ]     = d/dz   [ (x_phi x_z + y_phi y_z)/J ]
(E3)  d/dpsi [ (1 + x_z^2 + y_z^2)/J ]     = d/dz   [ (x_psi x_z + y_psi y_z)/J ]
```

Only two of (E1)-(E3) are independent: the contravariant curl
components

```
C_phi = d(b_z)/dpsi - d(b_psi)/dz
C_psi = d(b_phi)/dz - d(b_z)/dphi
C_z   = d(b_psi)/dphi - d(b_phi)/dpsi
```

satisfy the differential identity
`d(C_phi)/dphi + d(C_psi)/dpsi + d(C_z)/dz = 0` for ANY `x`, `y`
(verified symbolically). Two equations, two unknowns `x`, `y` — the
system is closed.

Equivalently, vacuum means a scalar potential `chi(phi, psi, z)`
exists with `b_i = d(chi)/du^i`; (E1)-(E3) are its integrability
conditions.

## 4. Variational (Legendre) structure — how x and y are determined

The magnetic energy transforms with `dx dy dz = J dphi dpsi dz` and
`|B|^2 = (1 + x_z^2 + y_z^2)/J^2` into

```
W[x, y] = 1/(2 mu_0) Int (1 + x_z^2 + y_z^2) / J  dphi dpsi dz
```

The Lagrangian density contains only first derivatives of the
unknowns — the dependent/independent roles of `(x, y)` and
`(phi, psi)` are exactly exchanged relative to the forward energy
functional `W[phi, psi] = Int |grad phi x grad psi|^2/(2 mu_0) dV`.
The Euler-Lagrange equations (with `G = 1 + x_z^2 + y_z^2`) are

```
d/dz( 2 x_z / J ) - d/dphi( G y_psi / J^2 ) + d/dpsi( G y_phi / J^2 ) = 0
d/dz( 2 y_z / J ) + d/dphi( G x_psi / J^2 ) - d/dpsi( G x_phi / J^2 ) = 0
```

These are exact combinations of (E1)-(E3): symbolically

```
EL_x =  2 ( C_phi y_phi + C_psi y_psi ) / J
EL_y = -2 ( C_phi x_phi + C_psi x_psi ) / J
```

i.e. the stationarity conditions are the force-balance projections
`(curl B) x B = 0` (force-free); the current-free (vacuum) solutions
are the subset on which additionally (E1) (`C_z = 0`) holds.

**Solution procedure**:

1. Computational domain: a box in `(phi, psi, z)` — `phi`, `psi` are
   flux coordinates, so a curved flux tube in physical space becomes a
   rectangle (the main practical advantage of the transformation).
2. Boundary conditions: at the tube ends `z = z0, z1` prescribe the
   field-line footpoint maps `x(phi, psi, z0)`, `y(phi, psi, z0)`
   (e.g. from the flux distribution on a pole face); the lateral
   boundaries `phi, psi = const` are flux surfaces on which `x, y`
   trace the prescribed tube wall.
3. Initial guess: the uniform-field solution `x = phi/B0`, `y = psi`.
4. Discretize `W` (FEM) and iterate with Newton/Picard, damping to
   preserve `J > 0` (mesh non-degeneracy = `B_z > 0`).

## 5. Exact solutions (sanity checks)

**Uniform field** `B = (0, 0, B0)`: `phi = B0 x`, `psi = y`, inverse
`x = phi/B0`, `y = psi`, `J = 1/B0`. Satisfies (E1)-(E3) trivially.

**Hyperbolic field** `B = (y, x, B0)` (curl-free, div-free):

```
phi = (x + y) exp(-z/B0),     psi = -(B0/2)(x - y) exp(z/B0)

x = [ phi exp(z/B0) - (2 psi/B0) exp(-z/B0) ] / 2
y = [ phi exp(z/B0) + (2 psi/B0) exp(-z/B0) ] / 2
```

The inverse map satisfies (E1)-(E3) identically and `1/J = B0`
(verified symbolically).

## 6. Partial Legendre transformation: swap only z <-> psi

Independent variables `(x, y, psi)`, unknowns `z(x, y, psi)` and
`phi(x, y, psi)`. Admissible where `psi_z != 0`, i.e. the flux
surfaces `psi = const` are graphs `z = z(x, y, psi)` over the
horizontal plane (complementary to the full swap, which requires
`B_z != 0`).

The chain rule gives (`J = z_psi`)

```
grad(psi) = ( -z_x, -z_y, 1 ) / z_psi

B = grad(phi) x grad(psi)
  = (1/z_psi) ( phi_y, -phi_x, z_x phi_y - z_y phi_x )
```

Note `phi_psi` drops out of `B` entirely (verified): only the
in-surface gradient of `phi` carries the field, so `phi` is the
**stream function of the field within each flux surface**. (Indeed
`A = phi grad(psi)` is a vector potential for `B`; in the
`(x, y, psi)` coordinates its only covariant component is
`a_psi = phi`.)

**`psi` is the flux-function label** — the 3D counterpart of the 2D
vector-potential flux function. In the 2D cross-section limit (field
in the x-z plane, translation-invariant along y, `B = (B_x, 0, B_z)`)
the admissible Clebsch pair is `phi = y`, `psi = -A_y(x, z)`: the
independent variable `psi` IS the flux function whose level lines are
the field lines, `phi` is trivial, and the single remaining unknown
`z(x, psi)` is the **field-line shape as a graph over (x, psi)** with

```
B_x = 1/z_psi,    B_z = z_x/z_psi          (dz/dx along a line = z_x)

W = 1/(2 mu_0) Int (1 + z_x^2)/z_psi dx dpsi      (per unit y)

d/dx( 2 z_x / z_psi ) = d/dpsi( (1 + z_x^2) / z_psi^2 )
```

— the classical inverse (pole-face design) equation: prescribe flux
surfaces, solve for the geometry.

Covariant components of `B`:

```
b_x   = [ (1 + z_x^2) phi_y - z_x z_y phi_x ] / z_psi
b_y   = [ z_x z_y phi_y - (1 + z_y^2) phi_x ] / z_psi
b_psi = z_x phi_y - z_y phi_x
```

**Transformed governing equations** (`curl B = 0`):

```
(F1)  d(b_y)/dx   = d(b_x)/dy
(F2)  d(b_psi)/dy = d(b_y)/dpsi
(F3)  d(b_psi)/dx = d(b_x)/dpsi
```

Again only two are independent (same differential identity on
`C_x = d(b_psi)/dy - d(b_y)/dpsi`, `C_y = d(b_x)/dpsi - d(b_psi)/dx`,
`C_psi = d(b_y)/dx - d(b_x)/dy`) — two equations for the two unknowns
`z`, `phi`.

**Variational structure**: with `dV = z_psi dx dy dpsi` and
`D = z_x phi_y - z_y phi_x`,

```
W[z, phi] = 1/(2 mu_0) Int [ phi_x^2 + phi_y^2 + D^2 ] / z_psi  dx dy dpsi
```

whose Euler-Lagrange equations are (verified symbolically)

```
EL_phi:  d/dx[ (phi_x - D z_y)/z_psi ] + d/dy[ (phi_y + D z_x)/z_psi ] = 0
         == -2 C_psi        (exactly (F1):  (phi_x - D z_y)/z_psi = -b_y,
                                            (phi_y + D z_x)/z_psi =  b_x)
EL_z:    d/dx[ 2 D phi_y/z_psi ] - d/dy[ 2 D phi_x/z_psi ]
         - d/dpsi[ (phi_x^2 + phi_y^2 + D^2)/z_psi^2 ] = 0
         == -2 ( C_x phi_x + C_y phi_y ) / z_psi
```

i.e. stationarity of `W` is again the force-free condition
`(curl B) x B = 0`; vacuum is the subset where additionally (F2)/(F3)
hold individually.

**Solution procedure**: `psi` spans the flux-label interval, so a
stack of flux surfaces becomes a slab in `(x, y, psi)`. Prescribe the
bounding flux surfaces `z(x, y, psi_0)`, `z(x, y, psi_1)` (e.g. pole
faces of a magnet gap are exact flux surfaces of the ideal field) and
`phi` on the side boundary (field-line entry/exit positions within
each surface); start from the uniform-field solution `z = psi/...`,
`phi = B0 y` and iterate Newton/Picard on the discretized `W`,
keeping `z_psi > 0` (surfaces must not cross).

**Exact solutions (sanity checks)**: uniform transverse field
`B = (B0, 0, 0)`: `z = psi`, `phi = B0 y`. Hyperbolic field
`B = (y, x, B0)` with `psi = (B0/2)(y - x) exp(z/B0)`:

```
z   = B0 log( 2 psi / (B0 (y - x)) )
phi = B0 (y^2 - x^2) / (2 psi)
```

satisfies (F1)-(F3) identically and reconstructs `B = (y, x, B0)`.

## 7. Relation to prior work

This formulation is **not claimed as novel in concept**. The variable-swap
mechanism (flux/potential labels as coordinates, geometry as the unknown,
energy minimization) is classical across several fields; see
[`LITERATURE.md`](LITERATURE.md) for the full survey. Most directly:

- **The 2-D case is published**: A. Dervisha, A. Marjamaki, P. Rasilo,
  T. Tarhasaari, "Bidirectional Coordinate Transformation and Its Application
  to 2-D Magnetic Field Problems", CEFC 2026 (Tampere University; Bossavit
  school) -- the exact 2-D case of the full swap (`(x,y) <-> (A, phi)`
  potential coordinates), derived via exterior calculus, with nonlinear
  `mu(A,phi)`. Their conclusion explicitly announces the manifold/3-D
  generalization as ongoing work ("the introduction part for a long story").
- **Lineage**: Clebsch/Euler potentials (Stern 1970); inverse flux-coordinate
  equilibrium (VMEC, Hirshman & Whitson 1983; inverse Grad-Shafranov; BETA,
  Bauer-Betancourt-Garabedian); curl-free fields with position-as-unknown
  (Boozer 2019); the partial-swap twin in fluid mechanics (von Mises 1927;
  Stanitz); conformal pole design (Rogowski 1923 / Halbach).

What this directory adds is an **independent, complementary 3-D derivation**
via the Clebsch pair `B = grad(phi) x grad(psi)` (vs. exterior calculus), with
explicit 3-D equations (E1)-(E3)/(F1)-(F3), a variational form, and
sympy-verified exact solutions -- a verified repository capability, regardless
of publication priority.

## Run the verification

**Symbolic** (sympy only) -- the governing equations and exact solutions:
```
python verify_clebsch_legendre_transform.py          # full swap (27 checks)
python verify_clebsch_partial_legendre_transform.py  # partial swap (26 checks)
```
Prints PASS/FAIL for each identity and exits nonzero on any failure.

**Numerical solver** (NGSolve) -- the Variant-2 (partial-swap) **3-D vacuum
field-geometry solver**: minimises `W = Int (phi_x^2+phi_y^2+D^2)/z_psi` over
`z(x,y,psi)`, `phi(x,y,psi)` by Newton, on a box in `(x,y,psi)`:
```
python solve_clebsch_legendre_3d.py
```
Manufactured-solution verified (golden `tests/feec/test_clebsch_legendre_3d.py`):
the polynomial vacuum case (uniform field) is recovered to machine precision and
the non-polynomial case (hyperbolic field `B=(y,x,B0)`) at the FE convergence
rate (order ~3 in `h`) -- confirming the energy minimiser on Clebsch fields is
the **vacuum** field (not merely force-free). This is the 3-D field-side piece
that the Tampere 2-D map (Sec. 7) had only announced.

## Field-plane (Chaplygin) design verification

The transforms above use the **potential** pair as coordinates. Two further
drivers validate the **field-plane** hodograph `(q, theta)` (Chaplygin), where
a saturating material law `mu(q)` becomes a known coefficient and the design
problem is exactly linear. The Legendre potential `chi = H.r - Psi` collapses
the two coordinate unknowns into one scalar whose first derivatives ARE the
physical coordinates; reparameterizing the radius by `B = mu q` feeds a
measured secant curve `mu_s(B)` into the coefficients with no table inversion.

**Machinery vs closed forms** (numpy only):
```
python verify_chi_modesum_solver.py
```
Radial mode solver against the exact `q^s` / Froehlich closed forms
(rel err < 1e-5 at M=401, FD order 2.0), 2-D mode-sum assembly with exact
coordinate recovery, MMF-type Robin boundary `q chi_q - chi = g`
(manufactured solutions), orientation monitor `J` single-signed.
Writes `results_chi_modesum_solver.json`.

**End-to-end design check** (NGSolve; the headline result):
```
python verify_chaplygin_bend_design.py
```
A 90-degree saturable flux-guide bend is DESIGNED by one linear hodograph
solve (walls are flux lines; outer wall 1.00 T constant, inner wall tapered
1.30 -> 1.75 T; `mu_r(B) = 1 + 199/(1+(B/1 T)^2)`), then the designed shape is
meshed and an INDEPENDENT nonlinear FEM (`div(nu(|grad A|) grad A) = 0`,
damped Picard) is solved on it with the same flux data. 2026-07-23 baseline:

| check | result |
|---|---|
| constant-mu sanity (design must be an exact annulus) | circularity dev 1.3e-9 / 3.5e-9, radius-ratio err 1.8e-9 |
| inner wall \|B\| vs spec (5..85 deg core) | mean 0.29--0.31 %, max 0.86 % |
| outer wall \|B\| vs spec (5..85 deg core) | mean 0.27--0.43 % (corner boundary layer only above that) |
| MMF, independent global quantity | design 908.7 A vs FEM 909.8 A (0.12 %); flat case 0.03 % |
| mesh convergence | maxh w/8 -> w/16 moves core errors by < 0.02 points |
| orientation | `J` keeps one sign on every design (no folding) |

Both drivers assert these golden bands and exit nonzero on violation; the
committed `results_*.json` sidecars are the durable records. Note the forward
check REQUIRES damped Picard (omega = 0.35): the undamped iteration stalls at
a large residual and fakes a ~140 % "design error" (bug-pattern
`reference-secant-picard-oscillation` in radia-mcp). Scope: this validates the
design direction on a known hodograph domain; the analysis direction (given a
fixed pole shape, the hodograph image is unknown) remains open research.

## References

- A. Clebsch, "Ueber die Integration der hydrodynamischen
  Gleichungen", J. Reine Angew. Math. 56 (1859).
- D. P. Stern, "Euler potentials", Am. J. Phys. 38, 494 (1970).
- W. D. D'haeseleer, W. N. G. Hitchon, J. D. Callan, J. L. Shohet,
  "Flux Coordinates and Magnetic Field Structure", Springer (1991) —
  inverse (flux-coordinate) formulations of magnetostatic equilibria.
- A. Dervisha, A. Marjamaki, P. Rasilo, T. Tarhasaari, "Bidirectional
  Coordinate Transformation and Its Application to 2-D Magnetic Field
  Problems", CEFC 2026 — the 2-D case of the full swap (potential coordinates
  `(A, phi)`), exterior-calculus derivation, nonlinear `mu`; announces the 3-D
  generalization as ongoing work.
- A. Bossavit, "Computational Electromagnetism: Variational Formulations,
  Complementarity, Edge Elements", Academic Press (1998) — the exterior-calculus
  formalism underlying the Tampere derivation.
- S. P. Hirshman, J. C. Whitson, "Steepest-descent moment method for
  three-dimensional magnetohydrodynamic equilibria" (VMEC), Phys. Fluids 26,
  3553 (1983) — inverse-coordinate equilibrium `x = x(rho, theta, zeta)`;
  the toroidal `p=0` parent of Variant 1.
- A. H. Boozer, "Curl-free magnetic fields for stellarator optimization",
  Phys. Plasmas 26, 102504 (2019) — curl-free field with position as the
  unknown in flux/angle coordinates (closest single prior line).
- R. von Mises, "Bemerkungen zur Hydrodynamik", ZAMM 7, 425 (1927) — the
  `(x, psi)` partial-swap with `y(x, psi)` unknown; the 2-D fluid twin of
  Variant 2.

See [`LITERATURE.md`](LITERATURE.md) for the full prior-art survey and the
honest novelty assessment.
