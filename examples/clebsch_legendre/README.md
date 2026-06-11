# Clebsch-Potential Magnetostatics: Legendre (Hodograph) Transformation to Flux Coordinates

Derivation of the magnetostatic equations with two Clebsch potentials
`phi`, `psi` as independent variables: the roles of dependent and
independent variables are swapped from `(x, y, z)` with unknowns
`phi(x,y,z)`, `psi(x,y,z)` to `(phi, psi, z)` with unknowns
`x(phi,psi,z)`, `y(phi,psi,z)`.

All equations below are verified symbolically by
[`verify_clebsch_legendre_transform.py`](verify_clebsch_legendre_transform.py)
(sympy, 27 checks).

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

## Run the verification

```
python verify_clebsch_legendre_transform.py
```

Requires only `sympy`. Prints PASS/FAIL for each of the 27 identities
above and exits nonzero on any failure.

## References

- A. Clebsch, "Ueber die Integration der hydrodynamischen
  Gleichungen", J. Reine Angew. Math. 56 (1859).
- D. P. Stern, "Euler potentials", Am. J. Phys. 38, 494 (1970).
- W. D. D'haeseleer, W. N. G. Hitchon, J. D. Callan, J. L. Shohet,
  "Flux Coordinates and Magnetic Field Structure", Springer (1991) —
  inverse (flux-coordinate) formulations of magnetostatic equilibria.
