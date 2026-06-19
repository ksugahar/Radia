"""
Kelvin transformation knowledge base for Radia MCP server.

Covers the Kelvin inversion technique for open boundary FEM problems
using NGSolve, as implemented in examples/kelvin_transformation/.

Canonical convention (Nagamine, Yamaguchi, Sugahara, CEFC 2026, id 350,
"A Pullback-Based Formulation of Kelvin Transformation in EM Field
Analysis"; see also Sugahara 2022 IEEE TransMag 58(9) = ref [3] in
Nagamine):

    3D spherical (conformal) Kelvin:
      nu' = (rho'/R)^2 * nu_0     [A-formulation, HCurl]
      mu' = (R/rho')^2 * mu_0     [Omega / H-formulation, H1]
      (pointwise reciprocals: mu * nu = 1)

    2D cylindrical (non-conformal) Kelvin -- ANISOTROPIC tensor:
      nu' = diag(1, 1, (rho'/R)^4) * nu_0
      mu' = diag(1, 1, (R/rho')^4) * mu_0
      In-plane slots (xx, yy) are identity -- common case factor = 1.
      Axial slot (zz) carries (rho'/R)^4 or (R/rho')^4 -- rare cases
      where B or H has a z-component (e.g. 2D A-form with in-plane A).

Derivation: pullback of orthonormal 1-form basis + bilinear energy
functional. Validated numerically on a toroidal current loop
(analytical dipole exterior energy 3.333e-8 J vs FEM on Omega'
= 3.344e-8 J, +0.33%).

Reference docs (consolidated 2026-05-04):
  examples/kelvin_transformation/CONVENTION.md (canonical declaration, 1 page)
  docs/kelvin/KELVIN_TRANSFORMATION.md (comprehensive theory + API + workflow)
    - §2: 1-form / 2-form pullback derivation
    - §7: Reduced potential formulations + Kelvin
      including the (nu - nu_0) form pitfall (CRITICAL)
"""

KELVIN_OVERVIEW = """
# Kelvin Transformation for Open Boundary Problems

The Kelvin transformation maps an unbounded exterior domain to a bounded
computational domain, enabling FEM solutions for open boundary
magnetostatic/electromagnetic problems without artificial truncation.

## Mathematical Foundation

### 2D (Circle Inversion)
- Maps: r' = R^2 / r
- Physical exterior: r > R  ->  Computational interior: 0 < r' < R
- Metric factor: dr'/dr = -R^2/r^2

### 3D (Sphere Inversion)
- Maps: rho' = R^2 / rho
- Physical exterior: rho > R  ->  Computational interior: 0 < rho' < R
- Metric factor: drho'/drho = -R^2/rho^2

## Key Principles

1. **Two Identical Spheres, Offset in Space**:
   - Interior sphere: physical domain (coil + iron + surrounding air), radius R
   - Exterior sphere (Kelvin): same radius R, offset in space (e.g., offset_x = 3R)
   - Both spheres have **identical radius and identical surface mesh**
   - Periodic BC couples DOFs 1:1 between the two boundaries via TRANSLATION
   - **WRONG**: concentric shell (R_inner to R_outer) -- this is NOT Kelvin
   - **CORRECT**: two separate spheres of the SAME radius R, placed apart

2. **Identical DOF Structure**:
   - Periodic BC requires the two identified boundaries to have matching DOFs
   - Both sphere surfaces must have the same mesh topology (same nodes, same elements)
   - OCC workflow: `Identify()` + `GenerateMesh()` handles this automatically
   - Cubit workflow: `copy mesh surface` from interior sphere to exterior sphere
     (MANDATORY -- without it, tet mesh gives different triangulations)
   - The relationship between identified nodes is a pure TRANSLATION (offset)

3. **Material Modulation** (exterior domain) -- Nagamine CEFC 2026:

   Kelvin transformation is a *physical* coordinate transformation of
   the same material (not an ad-hoc numerical trick). Material
   modulations obey the reciprocal relation `mu * nu = 1` pointwise.
   The "image of infinity" is at rho'=0 (Kelvin sphere center), where
   mu diverges and nu = 1/mu vanishes:

   - **Direct material quantities** (mu, eps, sigma) -> **infinity** at rho'=0
   - **Reciprocal quantities** (nu = 1/mu, rho = 1/sigma) -> **0** at rho'=0

   For the 3D spherical (conformal) Kelvin inversion rho -> R^2/rho':

   | Quantity   | Formula                  | rho'=0 (infinity) | rho'=R (boundary) |
   |------------|--------------------------|-------------------|-------------------|
   | mu'        | (R / rho')^2 * mu_0      | infinity          | mu_0 (continuous) |
   | epsilon'   | (R / rho')^2 * eps_0     | infinity          | eps_0             |
   | sigma'     | (R / rho')^2 * sigma     | infinity          | sigma             |
   | nu' = 1/mu'| (rho' / R)^2 * nu_0      | **0**             | nu_0              |
   | rho_e'    | (rho' / R)^2 * rho_e     | **0**             | rho_e             |

   - rho' = distance from the CENTER of the exterior sphere
   - Use mu' for H1/Omega scalar potential, nu' for HCurl A-formulation
   - Derivation: Nagamine CEFC 2026 eq. 9 (pullback of 1-form basis +
     bilinear energy equivalence)

   For 2D cylindrical (non-conformal) Kelvin, the reluctivity is an
   ANISOTROPIC tensor (Nagamine eq. 12) -- the Kelvin factor is NOT
   scalar:

       nu' = diag(1, 1, (rho'/R)^4) * nu     (only axial component scaled)
       mu' = diag(1, 1, (R/rho')^4) * mu

   Which slot contracts into the bilinear form depends on which B/H
   component the problem uses. The common case (in-plane B and H, as in
   2D Omega-form with scalar phi, or 2D A-form with A = A_z z_hat) uses
   only the in-plane identity slots, so the effective Kelvin factor is
   **1**. The (rho'/R)^4 / (R/rho')^4 axial slot matters only when B_z
   or H_z is present (rare in standard 2D magnetostatics).

   **Accuracy note**: if a Kelvin setup gives wrong inductance /
   field / energy, DO NOT change the convention. The convention is
   mathematically / physically fixed by Nagamine's derivation. Instead
   debug the FEM setup separately: GND placement (Dirichlet at
   rho'=0), gauge regularization magnitude, mesh refinement near
   rho'=0 (where nu vanishes / mu diverges), integration order
   (bonus_intorder for rational nu), and periodic identification.

4. **Background Field Modulation** (exterior domain):
   - Physical Kelvin-transformed field: H'_i = -(R/r')^2 * H_i(R^2/r')
   - The sign flip comes from the negative Jacobian (orientation reversal)
   - For uniform dipole H_s=(0,0,1): exterior H'_s = (0, 0, -(R/r')^2)
   - For quadrupole H_s=(-z,0,-x): exterior H'_s = +(R^4/r'^4)*(z',0,x')
     (positive because the double negation: field already negative + Kelvin sign flip)
   - At the Kelvin boundary r'=R, the rule simplifies to H'(R) = -H(R)

5. **Periodic Boundary Conditions**:
   - Identify interior outer boundary with exterior outer boundary
   - Ensures field continuity at the inversion sphere/circle
   - Uses NGSolve `Periodic()` FE space wrapper
   - Boundary terms in the weak form AUTO-CANCEL with periodic BC

6. **GND Point (Dirichlet at Infinity)**:
   - Place vertex at exterior domain center (r'=0)
   - Corresponds to physical infinity
   - **H1 (scalar potential)**: Essential for uniqueness: `dirichlet_bbnd="GND"`
   - **HCurl (vector potential)**: Optional but recommended. Gauge
     regularization (`reg * nu0 * u * v * dx`) provides uniqueness without GND.
     GND improves convergence for iterative solvers.
   - Cubit: create vertex at Kelvin sphere center, assign as nodeset "GND"

## Cubit Workflow: Offset Spheres (3D with .vol export)

For Cubit-based meshes exported via `export netgen`.
Verified with NGSolve 6.2.2603, Coreform Cubit 2025.12.

### Step-by-step procedure

```
1. Create interior sphere (radius R) centered at origin
   - Contains coil + iron + air
2. Create exterior sphere (SAME radius R) at offset position (e.g., x = 3R)
   - Contains only air (Kelvin domain)
3. Webcut BOTH spheres with the SAME plane (e.g., zplane)
   - ACIS sphere has 0 curves by default -- webcut creates 1 curve per hemisphere
   - This matching topology is required for copy mesh surface
4. Merge BOTH pairs of half-spheres:
   - imprint volume <inner_top> <inner_bottom>; merge volume <inner_top> <inner_bottom>
   - imprint volume <outer_top> <outer_bottom>; merge volume <outer_top> <outer_bottom>
   - This ensures equator nodes are shared -> inner node count = outer node count
5. Imprint/merge coil with air (interior sphere only, NOT with exterior)
6. Mesh interior hemisphere SURFACES first:
   - The air outer surface must be triangulated before copy mesh
   - `mesh volume {air_top} {air_bot}` creates the surface mesh as a side effect
   - Alternatively: `surface <N> scheme trimesh; mesh surface <N>` explicitly
7. copy mesh surface: interior sphere -> exterior sphere (MANDATORY)
   - copy mesh surface <in_sid> onto surface <out_sid> source curve <ic>
     source vertex <iv> target curve <ec> target vertex <ev>
   - Identify surfaces by area: hemisphere = largest area surface per volume
   - Guarantees 1:1 node correspondence for Periodic BC
   - WITHOUT copy mesh, tet meshing produces different triangulations -> crash
8. Mesh all volumes (tet mesh, using existing surface mesh as boundary constraint)
   - Do NOT set volume size for Kelvin volumes after copy (overrides copied mesh)
9. Set blocks (coil, air, kelvin) and sidesets (source, sink, kelvin_int, kelvin_ext)
   - sideset "kelvin_int": interior sphere outer surface (air boundary)
   - sideset "kelvin_ext": exterior sphere surface (kelvin boundary)
   - These MUST be set explicitly in the .py -- offset spheres do not share
     a face, so DomainIn/DomainOut auto-detection does NOT work
10. export netgen: writes periodic identification as TRANSLATION
    - C++ reads bc names "kelvin_int" / "kelvin_ext" from sidesets
    - Computes offset = mean(kelvin_ext vertices) - mean(kelvin_int vertices)
    - Bijective nearest-neighbor matching (no duplicate targets)
    - Writes ident.Add() pairs to .vol (NGSolve Periodic() reads them)
```

### Copy mesh: dynamic surface identification (Cubit Python)

APREPRO cannot iterate over surfaces or compare areas. Use Cubit Python
for the mesh copy step:

```python
import math
A_hemi = 2.0 * math.pi * R**2  # hemisphere area

for inner_vid, outer_vid in [(air_top, kelvin_top), (air_bot, kelvin_bot)]:
    # Find hemisphere surface (largest area) in each volume
    inner_surfs = cubit.get_relatives("volume", inner_vid, "surface")
    outer_surfs = cubit.get_relatives("volume", outer_vid, "surface")
    src = max(inner_surfs, key=lambda s: cubit.surface(s).area())
    dst = max(outer_surfs, key=lambda s: cubit.surface(s).area())

    # Pick longest curve (equator) for orientation reference
    src_c = max(cubit.get_relatives("surface", src, "curve"),
                key=lambda c: cubit.curve(c).length())
    dst_c = max(cubit.get_relatives("surface", dst, "curve"),
                key=lambda c: cubit.curve(c).length())
    src_v = cubit.get_relatives("curve", src_c, "vertex")[0]
    dst_v = cubit.get_relatives("curve", dst_c, "vertex")[0]

    cubit.cmd(f"copy mesh surface {src} onto surface {dst} "
              f"source curve {src_c} source vertex {src_v} "
              f"target curve {dst_c} target vertex {dst_v}")
```

### Critical gotchas

| Issue | Symptom | Fix |
|-------|---------|-----|
| No mesh copy | different inner/outer vertex counts, crash | `copy mesh surface` is MANDATORY |
| Volume size after copy | Extra nodes on copied surface | Do NOT set kelvin volume size |
| ACIS sphere has 0 curves | copy mesh fails (no curve) | `webcut with plane zplane` |
| Coil imprinted with kelvin | Extra nodes on kelvin boundary | Only imprint coil with air, NOT kelvin |
| `mesh volume` without surface | copy mesh produces 0 tris | Mesh air volumes first |

### Verified result (2026-04-13)

```
copy mesh: surface 10 -> 14, surface 12 -> 16 (1:1 hemisphere copy)
Kelvin periodic: 270 inner, 270 outer verts (mesh copy guarantees match)
offset=(0.1990, 0.0000, -0.0001), max_dist=9.65e-04, 0 unmatched
Periodic HCurl: ndof=38333, FreeDofs: 38333 -> 37529 (804 coupled)
FEM-Kelvin L=80.46 nH (7000 Hz, copper cylinder workpiece)
```

**Do NOT use**:
- Concentric shells (R_inner to R_outer) -- this is NOT Kelvin transformation
- Radial scaling for node pairing -- use translation (offset)
- Independent tet mesh on both spheres without copy mesh -- different triangulations

## Material Modulation -- The Only Thing You Change

Kelvin transformation does not change the PDE -- you write the standard
weak form on the bounded computational domain (interior + Kelvin sphere)
and **only modify the material coefficients in the Kelvin region**.

| Formulation | Bilinear form           | Kelvin region material (3D sphere, Nagamine) |
|-------------|------------------------|-------------------------------|
| H1 / Omega  | int mu * grad u . grad v dx | mu = (R/r')^2 * mu_0     |
| HCurl A     | int nu * curl u . curl v dx | nu = (r'/R)^2 * nu_0     |
| H1 (eddy)   | + sigma * u * v dx          | sigma = (R/r')^2 * sigma_0|
| HCurl (eddy)| + (1/sigma) curl u * curl v | rho = (r'/R)^2 * rho_0   |

**The rule is the physical rule** (no coordinate substitution needed):

- mu, eps, sigma diverge to infinity at r'=0 (image of physical infinity)
- nu, rho_e (reciprocals) vanish to 0 at r'=0
- All factors are continuous (= material_0) at r'=R (kelvin boundary)

NGSolve example (Nagamine CEFC 2026 canonical):
```python
nu_dict = {}
for m in mesh.GetMaterials():
    if "kelvin" in m.lower():
        rp_sq = (x-kx)**2 + (y-ky)**2 + (z-kz)**2 + 1e-20
        nu_dict[m] = NU_0 * rp_sq / R**2   # (r'/R)^2 * nu_0  -- vanishes at r'=0
    else:
        nu_dict[m] = NU_0
nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)
```

Alternatively use the centralized helper (radia 4.5.16+):
```python
from radia.kelvin_source import kelvin_nu_factor_3d_cf, build_material_cf
nu_factor = kelvin_nu_factor_3d_cf(center=(kx, ky, kz), R=R)
nu_cf = build_material_cf(mesh, NU_0, nu_factor, outer_keyword="kelvin")
```

The periodic BC at r'=R ensures continuity of the field across the
inversion boundary, completing the open-boundary closure.

## Supported Formulations

| Formulation | Primary Variable | Use Case |
|-------------|-----------------|----------|
| H-formulation | Perturbation potential phi | Permanent magnets, soft iron |
| A-formulation | Vector potential A | Coils, eddy currents |
| Omega-ReducedOmega | Scalar potentials | Domain decomposition |

## Advantages over Truncation

- **Exact**: No artificial boundary approximation
- **No PML tuning**: Eliminates absorbing layer parameters
- **Compact mesh**: Bounded domain even for infinite problems
- **Analytical validation**: Interior fields match exact solutions
"""

KELVIN_H_FORMULATION = """
# H-Formulation with Kelvin Transformation

The perturbation potential formulation: H = H_s - grad(phi)
where H_s is the known source field and phi is the unknown perturbation.

## Setup Steps

### 1. Geometry (OCC-based)

```python
from netgen.occ import *

circle_radius = 0.5    # Material radius [m]
kelvin_radius = 1.0    # Inversion radius [m]
offset_x = 3.0         # Separate interior/exterior domains

# Interior domain
wp = WorkPlane()
mag_circle = wp.Circle(circle_radius).Face()
mag_circle.name = "magnetic"
inner_boundary = wp.Circle(kelvin_radius).Face()
inner_air = inner_boundary - mag_circle
inner_air.name = "air_inner"

# Exterior (Kelvin) domain
wp_ext = WorkPlane(Axes(Pnt(offset_x, 0, 0), n=Z, h=X))
outer_circle = wp_ext.Circle(kelvin_radius).Face()
outer_circle.name = "air_outer"

# GND vertex at exterior center (= physical infinity)
gnd = Vertex(Pnt(offset_x, 0, 0))
gnd.name = "GND"
```

### 2. Periodic Identification

```python
# Name Kelvin boundary edges BEFORE boolean operations
for edge in inner_boundary.edges:
    edge.name = "kelvin_int"
for edge in outer_circle.edges:
    edge.name = "kelvin_ext"

# Boolean subtraction preserves edge names
inner_air = inner_boundary - mag_circle
inner_air.name = "air_inner"

shape = Glue([inner_air, mag_circle, outer_circle, gnd])

# Find edges by name and identify
kelvin_int = [e for e in shape.edges if e.name == "kelvin_int"]
kelvin_ext = [e for e in shape.edges if e.name == "kelvin_ext"]
for ie, ee in zip(kelvin_int, kelvin_ext):
    ie.Identify(ee, "periodic", IdentificationType.PERIODIC)

geo = OCCGeometry(shape, dim=2)
```

See topic `identify` for detailed best practices and anti-patterns.

### 3. FE Space with Periodic BC

```python
from ngsolve import *

mesh = Mesh(geo.GenerateMesh(maxh=0.03, grading=0.7))

fes_before = H1(mesh, order=3, dirichlet_bbnd="GND")
fes = Periodic(fes_before)  # Couples identified DOFs

# Verify DOF reduction
n_before = sum(1 for d in fes_before.FreeDofs() if d)
n_after = sum(1 for d in fes.FreeDofs() if d)
print(f"FreeDofs: {n_before} -> {n_after} (coupled: {n_before - n_after})")

u, v = fes.TnT()
```

### 4. Material Properties

```python
mu0 = 4 * pi * 1e-7
mu_r = 100

mu_dict = {
    "air_inner": mu0,
    "magnetic": mu_r * mu0,
    "air_outer": mu0,   # Will be modulated by Kelvin metric
}
mu = CoefficientFunction([mu_dict[m] for m in mesh.GetMaterials()])
```

### 5. Background Field (with Kelvin modulation)

```python
# Interior: uniform applied field
Hs_y_inner = 1.0

# Exterior: modulated by Kelvin metric
# H_s_exterior = -(rho'/R)^2 * H_s_interior
x_local = x - offset_x
r_prime = sqrt(x_local**2 + y**2)
Hs_y_outer = -(r_prime / kelvin_radius)**2 * Hs_y_inner

# Domain switching
is_exterior = IfPos(x - offset_x / 2, 1.0, 0.0)
Hs_y = (1 - is_exterior) * Hs_y_inner + is_exterior * Hs_y_outer
Hs = CoefficientFunction((0, Hs_y))
```

### 6. Weak Form and Solve

```python
# Bilinear: integral(mu * grad(u) . grad(v))
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

# Linear: integral(mu * grad(v) . Hs)  -- NO boundary terms with Kelvin!
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs) * dx

a.Assemble()
f.Assemble()

gfu = GridFunction(fes)
c = Preconditioner(a, type="local")
solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
           tol=1e-8, printrates=True, maxsteps=10000)
```

### 7. Post-processing

```python
H_pert = -grad(gfu)   # Perturbation field

# Analytical solution for 2D cylinder in uniform field
H_interior_analytical = -1.0 + 2.0 / (mu_r + 1)

# Evaluate
H_numerical = H_pert[1](mesh(0, 0))
error = abs(H_numerical - H_interior_analytical) / abs(H_interior_analytical)
print(f"Error: {error*100:.3f}%")
```
"""

KELVIN_A_FORMULATION = """
# A-Formulation with Kelvin Transformation

Vector potential formulation, particularly useful for coil problems.

## Axisymmetric Coil Example (Z-offset Kelvin)

In axisymmetric problems, the z-offset strategy keeps the r-coordinate
unchanged between interior and exterior domains, simplifying the
transformation. The physical domain lives on x >= 0 (r = x, z = y)
and uses H1 scalar FES with u = r * A_phi substitution.

**Same recipe also works with `axihenrotte`** (Phase B3, 2026-05-12).
Substitute `H1Henrotte(...)` for `H1(...)`, wrap with `Periodic(...)`,
pass `kelvin_mu_factor_axisym_cf` to `AxiHenrotteStiffnessBFI`, and set
`check_unused=False` on the BilinearForm.  Sphere Cu R=10 mm hits
Stoll to -0.001 % with this configuration.  See
`axifemm_documentation(topic="kelvin")` for the canonical recipe and
the documented gotchas (`mu` vs `nu` factor convention, element-
centroid mu sampling, Curve(2) trade-off).

**Verified 2026-04-14**: L is independent of Kelvin radius R within
0.07% for R in [0.08, 0.25] m on a R_coil=30mm + R_wp=25mm test case,
confirming the reluctivity modulation formula nu_0 * (rho'/R)^2 is
correct.

**Helper function** (canonical:
`packages/cubit-mesh-export/src/cubit_mesh_export/cubit_helpers/add_kelvin.py`;
back-compat shim at `src/radia/panels/add_kelvin.py`):

```python
from cubit_mesh_export.cubit_helpers.add_kelvin import add_kelvin_2d_axisym
# interior_shape: your 2D OCC face on x>=0 with outer arc named "kelvin_int"
shape, info = add_kelvin_2d_axisym(interior_shape, R=0.10, z_offset=0.25)
# info['axis_labels']  -> "axis|axis_ext" (pass to H1 dirichlet=...)
# info['z_offset']     -> z_offset used
# Periodic pairs are established; GND vertex at (0, z_offset) is added.
```

The caller is responsible for:
  1. Building the interior half-circle on x >= 0 with physical inclusions
     subtracted (coil, workpiece, etc.)
  2. Naming the outer arc edge(s) "kelvin_int" BEFORE passing to this function
  3. Naming "axis" on x = 0 edges (interior)
  4. Passing `dirichlet="axis|axis_ext"` and `dirichlet_bbnd="GND"` to
     the H1 FES, wrapping with `Periodic()`

```python
from ngsolve import *
from netgen.occ import *

coil_radius = 0.6       # Coil center radius [m]
kelvin_radius = 1.0     # Inversion radius
z_offset = 2.5          # Vertical offset for Kelvin domain

# Geometry: half-plane (r >= 0)
# Interior at z=0, Exterior at z=z_offset

# FE Space: u = r * A_theta
fes = H1(mesh, order=3, dirichlet="axis|axis_ext",
         dirichlet_bbnd="GND")
fes = Periodic(fes)
u, v = fes.TnT()

# Reluctivity with Kelvin modulation
mu0 = 4 * pi * 1e-7
nu0 = 1 / mu0
rho_prime = sqrt(x**2 + (y - z_offset)**2)
nu_outer = nu0 * (rho_prime / kelvin_radius)**2

nu = CoefficientFunction([nu_dict[m] for m in mesh.GetMaterials()])

# Weak form for axisymmetric A-formulation
r_coord = IfPos(x - 1e-10, x, 1e-10)  # Avoid r=0 singularity

a = BilinearForm(fes)
a += nu / r_coord * grad(u) * grad(v) * dx

f = LinearForm(fes)
f += J0 * v * dx("coil")  # J_theta source

a.Assemble()
f.Assemble()

# Solve and recover B-field
gfu = GridFunction(fes)
solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
           tol=1e-8, maxsteps=10000)

# B-field: Br = -(1/r)*du/dz, Bz = (1/r)*du/dr
grad_u = grad(gfu)
Br = -grad_u[1] / r_coord
Bz = grad_u[0] / r_coord
```
"""

KELVIN_3D = """
# 3D Kelvin Transformation

## 3D Sphere in Uniform Field (H-formulation)

```python
from ngsolve import *
from netgen.occ import *

sphere_radius = 0.5
kelvin_radius = 1.0
offset_x = 3.0

# Interior: sphere + air shell
mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere.mat("magnetic")
inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_air = inner_sphere - mag_sphere
inner_air.mat("air_inner")

# Exterior: Kelvin domain
outer_sphere = Sphere(Pnt(offset_x, 0, 0), kelvin_radius)
outer_sphere.mat("air_outer")

# GND at center of exterior
gnd = Vertex(Pnt(offset_x, 0, 0))
gnd.name = "GND"

# Name Kelvin boundary faces BEFORE boolean operations
for face in inner_sphere.faces:
    face.name = "kelvin_int"
for face in outer_sphere.faces:
    face.name = "kelvin_ext"

inner_air = inner_sphere - mag_sphere
inner_air.mat("air_inner")

shape = Glue([inner_air, mag_sphere, outer_sphere, gnd])

# Identify by name (see topic "identify" for details)
int_face = next(f for s in shape.solids for f in s.faces if f.name == "kelvin_int")
ext_face = next(f for s in shape.solids for f in s.faces if f.name == "kelvin_ext")
int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)
geo = OCCGeometry(shape)
mesh = Mesh(geo.GenerateMesh(maxh=0.03))

# Periodic FE space
fes_before = H1(mesh, order=3, dirichlet="GND")
fes = Periodic(fes_before)
u, v = fes.TnT()

# 3D Kelvin permeability modulation
x_local = x - offset_x
r_prime_sq = x_local**2 + y**2 + z**2
mu_outer = kelvin_radius**2 / (r_prime_sq + 1e-20) * mu0

# 3D Kelvin background field modulation
r_prime = sqrt(r_prime_sq)
Hz_outer = -(r_prime / kelvin_radius)**2 * Hz_interior

# Solve (same pattern as 2D)
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs) * dx

# Analytical validation for 3D sphere
Hz_analytical = 3.0 / (mu_r + 2)  # Interior field
```

## Key Differences from 2D

| Aspect | 2D | 3D |
|--------|----|----|
| Inversion | Circle: r' = R^2/r | Sphere: rho' = R^2/rho |
| mu modulation (H1)   | mu = (R/r')^2 * mu_0 | mu = (R/r')^2 * mu_0 |
| nu modulation (HCurl)| nu = (r'/R)^2 * nu_0 | nu = (r'/R)^2 * nu_0 |
| Linear form Kelvin factor | (R/r')^2 | (R/r')^4 |
| Physical H' (uniform) | H' = -H (no spatial variation) | H'_z = -(R/r')^2 |
| Analytical (dipole) | 2/(mu_r+1) | 3/(mu_r+2) |

**Physical rule**: mu, eps, sigma -> infinity at r'=0; nu, rho_e -> 0 at r'=0.

**Note on 2D anisotropy**: In 2D, the metric-based permeability is
anisotropic: constant in the r-direction, (R/r')^2 in the theta-direction.
However, in Cartesian coordinates with isotropic materials, these factors
cancel in the bilinear form, so the effective isotropic mu' = mu0 suffices.
"""

KELVIN_HCURL_3D = """
# 3D HCurl A-Formulation with Kelvin (calc_fem_kelvin.py)

The HCurl formulation solves for vector potential A directly in 3D,
used for coil problems with SIBC workpieces (induction heating).
This is what `calc_fem_kelvin.py` implements.

## Key Difference from H1 Formulations

| Aspect | H1 (phi/Omega) | HCurl (A) |
|--------|----------------|-----------|
| Primary variable | Scalar potential phi | Vector potential A |
| Bilinear form | mu * grad(u) * grad(v) | **nu** * curl(u) * curl(v) |
| Kelvin factor in bilinear | mu' = (R/r')^2 * mu0 | **nu' = (r'/R)^2 * nu0** |
| GND requirement | Essential (uniqueness) | Optional (gauge reg suffices) |
| Source term | mu * grad(v) . Hs dx | J * v * dx("coil") |
| Source modulation | Hs modulated in exterior | **No modulation** (J=0 in Kelvin) |

## Why nu = (r'/R)^2 * nu0 (Nagamine CEFC 2026 canonical)

Derivation per Nagamine, Yamaguchi, Sugahara, CEFC 2026 (id 350):
pullback of orthonormal 1-form basis + bilinear energy functional
gives, for 3D spherical (conformal) Kelvin inversion:

    nu' = (rho'/R)^2 * nu          (eq. 9)

Physically: Kelvin is a coordinate transformation of the same
material, so mu' * nu' = 1 pointwise. Under the inversion
r -> R^2/r', the image of physical infinity is at r'=0. Material
quantities behave at r'=0 according to whether they are direct or
reciprocal:

| Quantity | At r'=0 (image of infinity) | Modulation                |
|----------|------------------------------|---------------------------|
| mu       | infinity                     | mu = (R/r')^2 * mu_0      |
| eps      | infinity                     | eps = (R/r')^2 * eps_0    |
| sigma    | infinity                     | sigma = (R/r')^2 * sigma_0|
| **nu = 1/mu**  | **0**                  | **nu = (r'/R)^2 * nu_0**  |
| rho_e = 1/sigma| 0                      | rho_e = (r'/R)^2 * rho_e_0|

The HCurl A-formulation uses nu in the bilinear form, hence:
```
nu_kelvin = (r'/R)^2 * nu_0     -- vanishes at the kelvin center (r'=0)
nu_kelvin = nu_0 at r'=R         -- continuous with the air domain
```

This is purely a **material distortion** -- the PDE is unchanged.
Kelvin transformation is among the most physical of numerical techniques
for open boundary problems precisely because it does not modify the weak
form, only the constitutive coefficient in the Kelvin domain.

## Implementation (from calc_fem_kelvin.py)

```python
from ngsolve import *

MU_0 = 4e-7 * pi
NU_0 = 1.0 / MU_0

mesh = Mesh("model.vol")  # Contains "coil", "air", "kelvin" materials

# Detect Kelvin sphere
kelvin_center = detect_kelvin_offset(mesh)  # centroid of "kelvin" elements
a_kelvin = estimate_kelvin_radius(mesh)     # min vertex distance from center

# Reluctivity with Kelvin modulation
kx, ky, kz = kelvin_center
nu_dict = {}
for m in mesh.GetMaterials():
    if "kelvin" in m.lower():
        dx_k = x - kx
        dy_k = y - ky
        dz_k = z - kz
        rp_sq = dx_k*dx_k + dy_k*dy_k + dz_k*dz_k + 1e-20
        kelvin_fac = rp_sq / a_kelvin**2  # (r'/R)^2
        nu_dict[m] = NU_0 * kelvin_fac
    else:
        nu_dict[m] = NU_0
nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

# FE space: HCurl + Periodic (no GND needed for HCurl)
dirichlet_bnd = "GND" if "GND" in mesh.GetBoundaries() else ""
base_fes = HCurl(mesh, order=1, complex=True, nograds=True, dirichlet=dirichlet_bnd)
fes = Periodic(base_fes)  # Couples Kelvin sphere DOFs
u, v = fes.TnT()

# Bilinear form: curl-curl with Kelvin-modulated nu
a = BilinearForm(fes)
a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
# Gauge regularization: nograds=True removes H1 gradient subspace,
# but residual harmonic cohomology can still leave curl-curl singular.
# Restrict to non-kelvin materials because nu_kelvin = (r'/R)^2 * nu0
# vanishes at r'=0 and a constant NU_0 mass term would dominate there
# and corrupt the magnetic energy integral.
non_kelvin_mats = [m for m in mesh.GetMaterials() if "kelvin" not in m.lower()]
if reg > 0 and non_kelvin_mats:
    a += reg * NU_0 * u * v * dx("|".join(non_kelvin_mats))

# SIBC = Robin BC on hole boundary (conductor NOT meshed, NOT solved)
# Validated 2026-04-14: 2D axisym Kelvin + SIBC (3D HCurl, hole approach)
#   matches within 2% in L, 4% in P vs 2D full-res with delta/5 mesh.
#   WARNING: an earlier claim of "L<1%, P<2%" compared 2D SIBC against a
#   COARSE 2D full-res (mesh/3, not delta/5); both were under-converged.
#   True full-res reference is reference_2d_axisym.solve_2d_axisym(kelvin=True)
#   which auto-refines WP to max(delta/5, 0.2mm).
# Z_s for solid cylinder: rho*gamma*I1(ga)/I0(ga), gamma=sqrt(jw*mu*sigma)
if has_sibc:
    robin = 1j * omega / Z_s
    a += robin * u.Trace() * v.Trace() * ds("sibc")

# Source: J in coil only (NO source modulation in Kelvin domain)
f = LinearForm(fes)
f += J_source * v * dx("coil")

a.Assemble()
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

# Inductance from magnetic energy (nu_cf includes Kelvin modulation)
W_mag = Integrate(0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)),
                  mesh, order=10).real
L = 2 * W_mag / I_total**2
```

## Why No Source Modulation is Needed

Unlike H-formulation where Hs must be modulated in the exterior domain,
the HCurl A-formulation has J_source = 0 in the Kelvin domain (no coil there).
The linear form `J * v * dx("coil")` only integrates over the coil region
in the interior sphere. The Kelvin domain contribution comes entirely from
the bilinear form (curl-curl stiffness with nu' modulation).

## GND Vertex in Cubit Workflow

For HCurl, GND is optional but improves iterative solver convergence:

```
# In Cubit journal (after creating exterior sphere):
create vertex X {kelvin_offset_x} Y 0 Z 0
nodeset 100 add vertex {last_vertex_id}
nodeset 100 name "GND"
```

The vertex at the Kelvin sphere center maps to physical infinity (r -> inf).
Setting A = 0 there is physically correct (no field at infinity).

## bonus_intorder for Kelvin Domain

The varying nu_cf = (r'/R)^2 * nu0 is a rational function of coordinates.
Standard quadrature rules are insufficient. Use `bonus_intorder=4`:

```python
a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
```

Without bonus_intorder, integration error in the Kelvin domain can reach
several percent, leading to incorrect inductance values.
"""

KELVIN_ADAPTIVE = """
# Adaptive Mesh Refinement with Kelvin

## Error Estimators

1. **Zienkiewicz-Zhu (ZZ)**: Patch-recovery based, automatic
2. **Metric-based**: Remeshing with element-size control
3. **Uniform refinement**: Simple h-refinement for convergence study

## Workflow

```python
from ngsolve import *

for level in range(max_refinements):
    # Solve
    a.Assemble()
    f.Assemble()
    solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
               tol=1e-8, maxsteps=10000)

    # Compute error estimator (ZZ)
    flux = -mu * grad(gfu)
    fes_flux = VectorL2(mesh, order=order - 1)
    gf_flux = GridFunction(fes_flux)
    gf_flux.Set(flux)

    # Element-wise error
    err = Integrate(
        InnerProduct(flux - gf_flux, flux - gf_flux),
        mesh, element_wise=True
    )

    # Mark elements (Doerfler marking)
    max_err = max(err)
    for el in mesh.Elements():
        if err[el.nr] > 0.5 * max_err:
            mesh.SetRefinementFlag(el, True)

    # Refine
    mesh.Refine()

    # Update spaces
    fes.Update()
    gfu.Update()
    a.Assemble()
    f.Assemble()
```

## Typical Orders and Accuracy

| Order | Elements | Error (vs analytical) |
|-------|----------|-----------------------|
| 2     | 5,000    | ~0.5%                 |
| 3     | 5,000    | ~0.05%                |
| 4     | 5,000    | ~0.005%               |
| 3     | 20,000   | ~0.01%                |
"""

KELVIN_IDENTIFY = """
# Periodic Boundary Identification (Identify Method)

The `Identify()` method pairs geometric edges (2D) or faces (3D) on the
Kelvin boundary so that NGSolve's `Periodic()` FE space couples their DOFs.
Correct identification is critical -- if it fails silently, the solver runs
but produces wrong results.

## Recommended Approaches (Best to Worst)

### 1. Named Edge/Face Approach (Best for H-formulation with offset)

Name edges/faces during geometry construction **before** Boolean operations
(subtraction, Glue). Names survive OCC Boolean operations.

```python
# 2D Example: Name edges before subtraction
for edge in inner_circle.edges:
    edge.name = "kelvin_int"

inner_air = inner_circle - mag_circle  # Name survives subtraction

for edge in outer_circle.edges:
    edge.name = "kelvin_ext"

# After Glue, find by name
shape = Glue([inner_air, mag_circle, outer_circle, gnd])

kelvin_int_edges = [e for e in shape.edges if e.name == "kelvin_int"]
kelvin_ext_edges = [e for e in shape.edges if e.name == "kelvin_ext"]

for int_edge, ext_edge in zip(kelvin_int_edges, kelvin_ext_edges):
    int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
```

```python
# 3D Example: Name faces before subtraction
for face in inner_sphere.faces:
    face.name = "kelvin_int"

inner_air = inner_sphere - mag_sphere  # Name survives

for face in outer_sphere.faces:
    face.name = "kelvin_ext"

shape = Glue([inner_air, mag_sphere, outer_sphere, gnd])

kelvin_int_face = None
kelvin_ext_face = None
for solid in shape.solids:
    for face in solid.faces:
        if face.name == "kelvin_int":
            kelvin_int_face = face
        elif face.name == "kelvin_ext":
            kelvin_ext_face = face

kelvin_int_face.Identify(kelvin_ext_face, "periodic",
                         IdentificationType.PERIODIC)
```

### 2. Geometric Distance Predicate (Best for AdaptiveMesh / distance-based)

Use `abs(dist - kelvin_radius) < tolerance` to identify edges/faces on the
Kelvin boundary. The symmetric band avoids false matches.

```python
# Correct: symmetric band around Kelvin radius
for edge in air_inner.edges:
    dist = sqrt(edge.center.x**2 + edge.center.y**2)
    if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
        edge.name = "kelvin_int"
```

### 3. Vertex Distance (Good for A-formulation)

For edges that are full circles (where `edge.center` returns the geometric
center, NOT a point on the circle), use `edge.start` or `edge.vertices`
instead.

```python
for edge in air_inner.edges:
    p = edge.start
    dist = sqrt(p.x**2 + p.y**2)
    if abs(dist - kelvin_radius) < kelvin_radius * 0.1:
        edge.name = "kelvin_int"
```

## Anti-Patterns (Do NOT Use)

### Bounding Box Heuristic
```python
# BAD: Fragile, depends on geometry ordering
max_bbox = 0
for i, edge in enumerate(shape.edges):
    bbox = edge.bounding_box
    size = max(bbox[3]-bbox[0], bbox[4]-bbox[1], bbox[5]-bbox[2])
    if size > max_bbox:
        max_bbox = size
        kelvin_edge = edge
```

### Direct Index Access
```python
# BAD: Index depends on OCC internal ordering, breaks silently
geo.solids[0].faces[0].Identify(geo.solids[1].faces[0], "periodic")
```

### One-Sided Threshold
```python
# BAD: Matches everything beyond 80% of R (too broad)
if dist > kelvin_radius * 0.8:
    face.name = "kelvin_int"

# GOOD: Symmetric band centered on R
if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
    face.name = "kelvin_int"
```

### WorkPlane.Arc for Azimuthal Sectors
```python
# BAD: WorkPlane.Arc creates a TANGENTIAL arc (center ≠ origin)
wp = WorkPlane(Axes(Pnt(0,0,-R-0.1), n=Vec(0,0,1), h=Vec(1,0,0)))
wp.MoveTo(0, 0); wp.LineTo(R, 0)
wp.Arc(R, 60)   # Arc center is at (R, R, ...) NOT at origin!
wp.LineTo(0, 0)

# The resulting shape is NOT a circular sector!
# Vertex verification: Arc(R=4, 60°) from (4,0) gives endpoint (7.46, 2.0),
# NOT the expected (2.0, 3.46) for a true 60° rotation.

# GOOD: Use HalfSpace intersection for azimuthal wedges
import math
dphi = 2*math.pi/6
big_ball = Sphere(Pnt(0,0,0), R_large + 0.1)
h1 = HalfSpace(Pnt(0,0,0), Vec(0, -1, 0))
h2 = HalfSpace(Pnt(0,0,0), Vec(-math.sin(dphi), math.cos(dphi), 0))
wedge_vol = big_ball * h1 * h2
```

### Identify AFTER Glue on Glued Geometry's Faces
```python
# BAD: After Glue(), identifying faces of the glued shape
#       → NrIdentifications = 0 (silently fails)
geo_obj = Glue([inner_sec, outer_sec])
geo_obj.faces[2].Identify(geo_obj.faces[1], 'periodic', IdentificationType.PERIODIC, rot)

# GOOD: Identify BEFORE Glue on individual solid faces
inner_sec.faces[2].Identify(inner_sec.faces[1], 'periodic', IdentificationType.PERIODIC, rot)
outer_sec.faces[2].Identify(outer_sec.faces[1], 'periodic', IdentificationType.PERIODIC, rot)
geo_obj = Glue([inner_sec, outer_sec])  # Glue preserves identifications
```

## Verification: FreeDofs Check

Always verify that `Periodic()` actually coupled DOFs:

```python
fes_before = H1(mesh, order=3, dirichlet_bbnd="GND")
freedof_before = sum(1 for d in fes_before.FreeDofs() if d)

fes = Periodic(fes_before)
freedof_after = sum(1 for d in fes.FreeDofs() if d)

print(f"FreeDofs: {freedof_before} -> {freedof_after}")
if freedof_before == freedof_after:
    raise RuntimeError("Periodic BC NOT working! Check Identify() setup.")
```

If FreeDofs does not decrease, the identification failed silently.
Common causes:
- Names lost during Boolean operations (name before subtraction/Glue)
- Wrong tolerance in distance predicate
- `edge.center` returns geometric center for full circles (use `edge.start`)

## Axisymmetric Identification (Named Edge + Z-Sign)

For axisymmetric problems with z-offset, match interior/exterior edges
by name and z-coordinate sign:

```python
int_edges = [(e, e.center.y) for e in shape.edges if e.name == "kelvin_int"]
ext_edges = [(e, e.center.y) for e in shape.edges if e.name == "kelvin_ext"]

for ie, iy in int_edges:
    for ee, ey in ext_edges:
        if iy * ey > 0:  # Same z-sign
            ie.Identify(ee, "periodic", IdentificationType.PERIODIC)
```
"""

KELVIN_VERIFICATION = """
# Kelvin変換の数値検証 (Single-Domain Approach)

## 概要

Kelvin 変換の理論的正しさを FEM で検証する最も簡潔な方法は、
**単一の有限領域** (大球) で Poisson 方程式を解き、その FEM 解に対して
Kelvin 変換恒等式を後処理として確認することである。

## 手法: 単一大球 + 厳密 Dirichlet BC

### 設定
- 電荷球: 半径 R_in=0.5, 一様電荷密度 ρ=3.0
- 計算領域: 半径 R_large=3.0 の大球 (外側境界に厳密値を与える)
- Kelvin 変換半径: R_K=1.0

### 理論 (モノポール解)
外部解: φ(r) = ρ·R_in³ / (3r) = 0.125/r
Kelvin 変換: φ*(r*) = (R_K/r*) · φ(R_K²/r*) = ρ·R_in³/(3·R_K) = 0.125 (定数)

### NGSolve 実装

```python
from ngsolve import *
import netgen.occ as occ
import numpy as np

R_in = 0.5; R_K = 1.0; R_large = 3.0; rho_val = 3.0

phi_BC         = rho_val * R_in**3 / (3.0 * R_large)   # 0.041667 (外側Dirichlet)
phi_star_exact = rho_val * R_in**3 / (3.0 * R_K)        # 0.125000 (Kelvin定数)

# ジオメトリ: 内球 + 外殻
inner = occ.Sphere(occ.Pnt(0,0,0), R_in); inner.mat('inner'); inner.maxh = 0.05
outer_ball = occ.Sphere(occ.Pnt(0,0,0), R_large); outer_ball.bc('dirichlet_outer')
outer = outer_ball - occ.Sphere(occ.Pnt(0,0,0), R_in); outer.mat('outer'); outer.maxh = 0.2

geo  = occ.OCCGeometry(occ.Glue([inner, outer]))
mesh = Mesh(geo.GenerateMesh(maxh=0.2))
# → ~87,421 要素, 16,564 節点, 123,700 DOF

fes = H1(mesh, order=2, dirichlet='dirichlet_outer')
u, v = fes.TnT()
a = BilinearForm(fes); a += grad(u)*grad(v)*dx; a.Assemble()
f = LinearForm(fes);   f += rho_val*v*dx('inner'); f.Assemble()

gfu = GridFunction(fes)
gfu.Set(phi_BC, definedon=mesh.Boundaries('dirichlet_outer'))
inv = a.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky')
gfu.vec.data += inv * (f.vec - a.mat * gfu.vec)
```

### 検証 1: 物理中間領域 FEM vs 解析解

```python
sample_r = np.linspace(R_in*1.05, R_K*0.95, 25)
phi_fem = np.array([gfu(mesh(r, 0, 0)) for r in sample_r])
phi_ana = np.array([rho_val*R_in**3/(3.0*r) for r in sample_r])
err_mid = np.abs(phi_fem - phi_ana) / np.abs(phi_ana)
print(f'最大相対誤差: {err_mid.max():.2e}')  # 期待値: < 0.5%
```

### 検証 2: Kelvin 変換恒等式

```python
sample_rstar = np.linspace(R_K**2/R_large*1.1, R_K*0.95, 25)
phi_star_fem = np.array([
    (R_K/rs) * gfu(mesh(R_K**2/rs, 0, 0))
    for rs in sample_rstar
])
err_kelvin = np.abs(phi_star_fem - phi_star_exact) / phi_star_exact
print(f'Kelvin 恒等式 最大誤差: {err_kelvin.max():.2e}')  # 期待値: < 0.4%
print(f'phi* 範囲: [{phi_star_fem.min():.5f}, {phi_star_fem.max():.5f}]')
# 期待: [0.12465, 0.12489] (理想: 0.12500)
```

## 確認済み精度 (NGSolve 6.2.2603+, maxh=0.2)

| 検証項目 | 最大相対誤差 | 平均相対誤差 |
|----------|-------------|-------------|
| 物理中間領域 FEM vs 解析解 | 3.53e-03 (0.35%) | 3.07e-03 |
| Kelvin 恒等式 φ*(r*) 定数性 | 2.81e-03 (0.28%) | 2.06e-03 |

## メッシュ精度の影響

精度を上げるには:
- `inner.maxh = 0.05` → 電荷球内を細かく (必須)
- `outer.maxh = 0.2`  → 外殻は粗くてよい
- `R_large` は大きいほど精度良いが、`R_large=3` で十分 (< 0.4% 誤差)
- `R_large=5` かつ `maxh=0.3` では誤差 ~9.7% (粗すぎ)

## 旧来の TnT() 呼び出し注意

```python
# NGSolve 6.2.2603+: キーワード引数不可
u, v = fes.TnT()           # OK
u, v = fes.TnT(u_='u', v_='v')  # TypeError: NGSolve 6.2.2603+ でも不可
```
"""

KELVIN_PERIODIC_WEDGE = """
# 周期境界条件: 扇形セクタの DOF 削減

## 概要

軸対称・回転対称問題では、全球の代わりに **1/n セクタ** のみを解くことで計算コストを削減できる。
NGSolve では `Identify` + `Periodic(FESpace)` でこれを実現する。

## 正しい扇形ウェッジの作成

`WorkPlane.Arc` は **接線方向弧**（中心が原点でない）を作るため扇形 (pizza slice) にならない。
`HalfSpace` で 2 枚の半空間を切り出すのが正確。

```python
from netgen.occ import HalfSpace, Pnt, Vec, Sphere, Glue, OCCGeometry, IdentificationType, Rotation
import math

n_sector = 6                     # 6 分割 → 60° セクタ
dphi     = 2 * math.pi / n_sector

# NG: 半径 R_large+0.1 の大球を 2 つの HalfSpace で切断
big_ball = Sphere(Pnt(0,0,0), R_large + 0.1)
h1 = HalfSpace(Pnt(0,0,0), Vec(0, -1, 0))                              # y >= 0 (φ >= 0°)
h2 = HalfSpace(Pnt(0,0,0), Vec(-math.sin(dphi), math.cos(dphi), 0))    # φ <= dphi
wedge_vol = big_ball * h1 * h2

# 球殻をウェッジで切断
inner_sec = Sphere(Pnt(0,0,0), R_in)                       * wedge_vol
outer_sec = (Sphere(Pnt(0,0,0), R_large) - Sphere(Pnt(0,0,0), R_in)) * wedge_vol
inner_sec.mat('inner'); outer_sec.mat('outer')
outer_sec.faces[0].name = 'dirichlet_outer'   # 最大面積の外球面
```

## 周期境界の識別 (Identify)

**重要**: Identify は `Glue()` より**前**に個別のソリッド上の face で呼び出す。

```python
# face[2]: φ=0° 面 (y≈0 側), face[1]: φ=dphi 面 (y>0 側)
# → 面積と重心 y 座標で確認可能:
#     face[2]: area≈π·R²/2, y_center≈0
#     face[1]: area≈π·R²/2, y_center>0

rot = Rotation(Pnt(0,0,0), Vec(0,0,1), math.degrees(dphi))
inner_sec.faces[2].Identify(inner_sec.faces[1], 'periodic', IdentificationType.PERIODIC, rot)
outer_sec.faces[2].Identify(outer_sec.faces[1], 'periodic', IdentificationType.PERIODIC, rot)
# ↑ Glue() の前に実行

geo_obj  = Glue([inner_sec, outer_sec])
geo      = OCCGeometry(geo_obj)
mesh_sec = Mesh(geo.GenerateMesh(maxh=0.2))
```

## DOF 削減の確認

```python
pairs = mesh_sec.GetPeriodicNodePairs(VERTEX)
print(f'周期頂点対数: {len(list(pairs))}')   # 0 なら Identify 失敗

fes_plain = H1(mesh_sec, order=2, dirichlet='dirichlet_outer')
fes_per   = Periodic(H1(mesh_sec, order=2, dirichlet='dirichlet_outer'))

free_plain = sum(fes_plain.FreeDofs())
free_per   = sum(fes_per.FreeDofs())

print(f'ndof (変化なし): {fes_plain.ndof}')
print(f'FreeDofs 前: {free_plain}')
print(f'FreeDofs 後: {free_per}  (削減: {free_plain - free_per})')
# 例: 12022 → 10555 (1467 DOF 削減, ~12%)
```

**注意**: `ndof` は同じまま、`FreeDofs()` が減る。
周期スレーブ DOF は拘束扱いとなり線形系から除外される。

## Identify が失敗する典型的な原因

| 現象 | 原因 | 対処 |
|------|------|------|
| 周期頂点対数 = 0 | Glue 後の face オブジェクトに対して Identify | Glue 前に呼び出す |
| 周期頂点対数 > 0 でも FreeDofs 変化なし | Rotation の方向が逆 | `dphi` の符号を確認 |
| NrIdentifications = 0 | WorkPlane.Arc で作った扇形 (幾何的に一致しない) | HalfSpace を使用 |

## 確認済み数値 (NGSolve 6.2.2603+, 1/6 セクタ, maxh=0.2)

| 項目 | 値 |
|------|----|
| 要素数 | 2,738 |
| 節点数 | 781 |
| 周期頂点対数 | 198 |
| FreeDofs (Periodic なし) | 12,022 |
| FreeDofs (Periodic あり) | 10,555 |
| 削減 DOF 数 | 1,467 (12.2%) |
| セクタ中央 最大誤差 | < 1.5% (maxh=0.2) / < 1% (maxh=0.1) |
"""

KELVIN_TIPS = """
# Kelvin Transformation: Tips and Pitfalls

## Common Mistakes

1. **Forgetting material modulation in exterior domain**
   - **Direct quantities** (mu, eps, sigma): diverge at r'=0
       mu = (R/rho')^2 * mu_0
   - **Reciprocal quantities** (nu, rho_e): vanish at r'=0
       nu = (rho'/R)^2 * nu_0
   - Use mu modulation for H1/Omega (scalar potential)
   - Use nu modulation for HCurl A-formulation
   - The PDE itself is unchanged -- only the material is distorted

2. **Wrong sign on background field modulation**
   - The Kelvin BC at r'=R requires: H'(R) = -H(R) (sign flip)
   - For uniform dipole field: sign flip needed (interior positive -> exterior negative)
   - For quadrupole H_s=(-z,0,-x): NO explicit sign flip needed in code
     (interior already negative -> Kelvin negation makes it positive)
   - Rule: the physical transformed field is H'_i = -(R/r')^2 * H_i(R^2/r')

3. **Missing GND point**
   - **H1 (phi, Omega)**: Without GND, solution is non-unique. Essential.
   - **HCurl (A)**: Gauge regularization provides uniqueness. GND is optional
     but improves iterative solver convergence.
   - Place vertex at center of exterior sphere (maps to physical infinity).
   - Cubit: `create vertex X {offset} Y 0 Z 0; nodeset N name "GND"`

4. **Boundary terms in linear form**
   - With Kelvin + Periodic BC, boundary terms CANCEL
   - Only volume integral remains: f(v) = int(mu * grad(v) . Hs) dx
   - Adding boundary terms leads to wrong results

5. **Domain offset too small**
   - Interior and exterior domains must NOT overlap in mesh
   - offset_x should be >= 2 * kelvin_radius

6. **Face identification order**
   - Identify interior outer face WITH exterior outer face
   - Order matters for orientation of periodic coupling

7. **HCurl without nograds=True**
   - HCurl basis functions include gradients of H1 basis which form the
     curl-curl kernel (curl(grad phi) = 0 = gauge freedom)
   - Use `HCurl(mesh, order=p, complex=..., nograds=True, dirichlet=...)`
     at ANY polynomial order p (including p=1).
   - Reference: NGSolve Maxwell tutorial unit-2.4
     `HCurl(mesh, order=3, dirichlet="outer", nograds=True)`
   - Without nograds, the curl-curl system is rank deficient by the
     gradient subspace dimension. PARDISO pseudo-inverses and returns a
     solution polluted by an arbitrary gradient field; inductance or
     magnetic energy values are then order-dependent and unreliable.

8. **Cubit coil: sweep, NOT webcut**
   - `create torus + webcut + delete` can silently produce a half-coil
     (BEM = 46 nH instead of 87 nH for R=30mm/a=3mm/gap=5 deg)
   - `webcut` also breaks p-convergence of curved surface elements
   - Correct: build a disk cross section and sweep about the axis:
     `create curve arc ... normal 0 1 0 full;
      create surface curve 1;
      sweep surface 1 axis 0 0 0 0 0 1 angle 355`
   - Verified p-convergence O1 → O5: -3.35% → -1.8e-05%

## Performance Tips

- Use `Preconditioner(a, type="local")` for small-medium problems
- Use `Preconditioner(a, type="multigrid")` for large problems
- Kelvin domain typically needs FEWER elements than interior
- Higher polynomial order (3-4) is more efficient than h-refinement

## When to Use Kelvin

| Problem Type | Kelvin? | Alternative |
|-------------|---------|-------------|
| Permanent magnet stray field | Yes | Truncation with large air box |
| Coil field in free space | Yes | Biot-Savart (analytical) |
| Soft iron in external field | Yes | Infinite element methods |
| Enclosed cavity | No | Standard FEM (bounded domain) |
| Waveguide (1D open) | No | PML (absorbing boundary) |
"""


KELVIN_ROBUSTNESS = """
# Kelvin Robustness Checklist (2026-04-14)

## POLICY: Kelvin Only -- No Dirichlet Truncation

Use Kelvin transformation for ALL FEM open boundary problems.
Do NOT use large Dirichlet truncation spheres. Kelvin provides exact
open boundary with zero truncation error. Dirichlet truncation introduces
uncontrolled error that depends on sphere size and cannot be bounded.

## Checklist: Steps That MUST Succeed

### 1. Mesh Copy (Interior -> Exterior Sphere)

Both sphere surfaces must have **identical mesh topology** (1:1 node
correspondence). Without this, Periodic BC cannot pair DOFs.

**OCC workflow**: `Identify()` before `GenerateMesh()` handles this
automatically. Netgen matches identified surfaces during mesh generation.

**Cubit workflow**: `copy mesh surface` is MANDATORY:
```python
# After meshing interior air volumes (which creates their surface mesh):
copy mesh surface <interior_hemi> onto surface <exterior_hemi> \\
    source curve <ic> source vertex <iv> \\
    target curve <ec> target vertex <ev>
```
Do NOT set volume size on the exterior Kelvin volumes after copy
(this adds extra nodes and breaks the 1:1 correspondence).

**Verification**:
```python
# After export, check vertex counts
kelvin_int_verts = set()
kelvin_ext_verts = set()
for el in mesh.Elements(BND):
    bnd = mesh.GetBoundaries()[el.index]
    for v in el.vertices:
        if "kelvin_int" in bnd: kelvin_int_verts.add(v.nr)
        elif "kelvin_ext" in bnd: kelvin_ext_verts.add(v.nr)
assert len(kelvin_int_verts) == len(kelvin_ext_verts), \\
    f"Mesh copy failed: {len(kelvin_int_verts)} int vs {len(kelvin_ext_verts)} ext"
```

### 2. Exterior Sphere Must Be Meshed

The exterior (Kelvin) sphere is a real computational domain -- it must
contain volume elements with material "kelvin". It is NOT just a boundary.

**Common mistake**: creating the exterior sphere but forgetting to mesh it
or forgetting to assign it to a block.

```python
# Cubit: exterior sphere MUST be meshed
volume {kelvin_top} {kelvin_bot} scheme tetmesh
mesh volume {kelvin_top} {kelvin_bot}

block N add volume {kelvin_top} {kelvin_bot}
block N name "kelvin"
```

### 3. Material Scaling in Kelvin Domain

Kelvin transformation **only distorts material properties** -- the PDE is
unchanged. The image of physical infinity is at r'=0 (Kelvin sphere center).

| Quantity | Physical (r→inf) | Kelvin (r'→0) | Modulation |
|----------|-----------------|---------------|------------|
| mu       | mu_0            | **infinity**  | mu = (R/r')^2 * mu_0 |
| eps      | eps_0           | **infinity**  | eps = (R/r')^2 * eps_0 |
| sigma    | 0 (air)         | 0             | sigma = (R/r')^2 * sigma (0 for air) |
| **nu = 1/mu** | nu_0      | **0**         | **nu = (r'/R)^2 * nu_0** |
| rho_e = 1/sigma | inf     | inf           | rho_e = (r'/R)^2 * rho_e |

**Implementation**:
```python
# HCurl A-formulation: uses nu (not mu)
for m in mesh.GetMaterials():
    if "kelvin" in m.lower():
        rp_sq = (x-kx)**2 + (y-ky)**2 + (z-kz)**2 + 1e-20
        nu_dict[m] = NU_0 * rp_sq / R**2   # -> 0 at r'=0
    else:
        nu_dict[m] = NU_0
```

**The 1e-20 epsilon prevents division by zero but is NOT a hack**:
at r'=0 exactly, nu=0 (correct), and the 1e-20 ensures no NaN.
The gauge regularization `reg * NU_0 * u * v * dx` in the bilinear form
provides the missing mass term where nu vanishes.

**IMPORTANT for gauge regularization with Kelvin**:
Do NOT apply `reg * NU_0 * u * v * dx` in the Kelvin domain -- the
nu_cf = (r'/R)^2 * nu0 vanishes near r'=0, so a constant NU_0 mass
term would dominate and corrupt the magnetic energy integral there.
Restrict regularization to non-kelvin materials:
```python
non_kelvin_mats = [m for m in materials if "kelvin" not in m.lower()]
if non_kelvin_mats:
    a += reg * NU_0 * u * v * dx("|".join(non_kelvin_mats))
```

### 4. kelvin_int / kelvin_ext Identification

The periodic identification pairs must be embedded in the .vol file.
Without it, the Kelvin sphere is just a disconnected air domain with
Neumann BC (silently wrong, looks like a small Dirichlet truncation).

**OCC**: Call `Identify()` before `GenerateMesh()`:
```python
int_face.Identify(ext_face, "kelvin", IdentificationType.PERIODIC)
```

**Cubit**: The C++ exporter (`export netgen`) reads sideset names
"kelvin_int" and "kelvin_ext", computes the translation offset, and writes
`ident.Add()` pairs to the .vol automatically.

**Detection in calc_fem_kelvin.py**:
```python
n_ident = mesh.ngmesh.GetNrIdentifications()
if n_ident == 0:
    raise RuntimeError("Kelvin without periodic identification!")
```

### 5. FreeDofs Decrease After Periodic()

This is the **definitive verification** that periodic identification is
working. If FreeDofs does not decrease, the identification failed silently.

```python
base_fes = HCurl(mesh, order=1, complex=True, nograds=True, dirichlet="GND")
n_before = sum(1 for d in base_fes.FreeDofs() if d)

fes = Periodic(base_fes)
n_after = sum(1 for d in fes.FreeDofs() if d)

n_coupled = n_before - n_after
print(f"FreeDofs: {n_before} -> {n_after} (coupled: {n_coupled})")
if n_coupled == 0:
    raise RuntimeError(
        "Periodic BC NOT working! FreeDofs unchanged. "
        "Check: (1) Identify() before GenerateMesh/Glue, "
        "(2) sideset names kelvin_int/kelvin_ext in .jou, "
        "(3) export netgen wrote identification pairs.")
```

Typical: ~2-5% of DOFs are coupled (boundary DOFs on the Kelvin spheres).
Example: ndof=49244, FreeDofs: 49244 -> 48440 (804 coupled).

### 6. Symmetry Models: 1/2, 1/4, 1/8

Kelvin + symmetry = periodic wedge + Kelvin open boundary.

**1/2 model** (one symmetry plane, e.g., z=0):
- Webcut interior AND exterior spheres with the symmetry plane
- Symmetry BC (Neumann for A_n=0, Dirichlet for A_t=0)
- Kelvin identification only between matching hemispheres

**1/4 model** (two symmetry planes):
- Webcut with both planes
- HalfSpace intersection for clean wedge geometry
- Identify BEFORE Glue (critical for OCC)

**1/8 model** (three symmetry planes):
- Use HalfSpace intersection to create octant
- Same Kelvin principles apply

**Critical**: the symmetry planes must cut BOTH interior and exterior
spheres identically. The mesh copy must respect the symmetry boundary.

**Example (1/2 model with z-symmetry)**:
```python
# In Cubit .jou:
webcut volume {air_vid} with plane zplane
webcut volume {kelvin_vid} with plane zplane
# Keep only z>0 halves for mesh copy + identification
# Apply symmetry BC (Neumann or Dirichlet) on z=0 faces
```

### 7. GND at Infinity

The vertex at the exterior sphere center maps to physical infinity.

| Formulation | GND requirement |
|-------------|----------------|
| **H1 (phi, Omega)** | **Essential** -- uniqueness of scalar potential |
| **HCurl (A)** | Optional -- gauge `reg * nu * u * v * dx` suffices |
| **HCurl + iterative solver** | Recommended -- improves convergence |
| **HCurl + PARDISO** | Optional -- PARDISO handles near-singular |

**Implementation**:
```python
# In Cubit .jou:
create vertex X {kelvin_offset_x} Y 0 Z 0
nodeset 100 add vertex {last_vertex_id}
nodeset 100 name "GND"

# In calc_fem_kelvin.py:
if "GND" in boundaries:
    dirichlet_bnd = "GND"  # A(infinity) = 0
```

**For OCC workflow**:
```python
gnd = Vertex(Pnt(kelvin_offset, 0, 0))
gnd.name = "GND"
# Include in Glue
shape = Glue([air, torus, ext_sphere, gnd])
```

### 8. Kelvin Domain: Tet Only (No Hex)

**POLICY**: Kelvin domain must use **tetrahedral** mesh, not hexahedral.

Spherical geometry is best approximated by triangular high-order elements.
Hex elements on a sphere surface introduce systematic geometry error that
worsens with mesh refinement (opposite of convergence).

- OCC workflow: always produces tets (correct by default)
- Cubit workflow: use `scheme tetmesh` for Kelvin volumes
- Do NOT use hex sweep or mapped mesh on the Kelvin sphere

This is why the Kelvin sphere is always a Netgen/OCC tet mesh, even when
the physical domain (coil, iron yoke) uses Cubit hex mesh.

## Diagnostic Summary

| Check | Expected | If Wrong |
|-------|----------|----------|
| Interior/exterior vertex count | Equal | Mesh copy failed |
| kelvin material in mesh | Present | Exterior sphere not meshed/blocked |
| nu at r'=0 | 0 | Wrong material scaling (mu not nu) |
| NrIdentifications > 0 | True | No periodic ID in .vol |
| FreeDofs decrease | Yes (2-5%) | Identification failed silently |
| L vs analytical | < 5% | Check all above + mesh size |
| Iterative solver converges | Yes | Check GND, reg, nograds |

## Validated Results (2026-04-14)

| Test case | L [nH] | Analytical | Error | Configuration |
|-----------|--------|------------|-------|---------------|
| Torus R=30mm a=3mm 5deg gap (3D) | 90.71 | 88.5 | +2.5% | Kelvin, Cu 50kHz, hole SIBC |
| Torus (2D axisym) | ~89 | 88.5 | <1% | Kelvin, static |
| C-magnet york.jou (3D) | 80.46 | -- | -- | Kelvin, 7kHz, mesh copy |
"""


KELVIN_HELPERS_RECIPE = """
# Student-friendly Kelvin helpers (radia.kelvin_source)

The canonical way to apply Kelvin material modulation in new code is
to import factor CFs from `radia.kelvin_source`. This avoids re-deriving
the factor inline and guarantees Nagamine CEFC 2026 consistency.

## Available factor helpers

Most Radia-NGSolve problems are **3D**, so the 3D helpers below are the
primary case. Axisymmetric is a dimensional reduction; 2D is special.

### 3D spherical Kelvin (primary case, isotropic)

| Formulation | Helper                       | Factor     | At rho'=0 |
|-------------|------------------------------|------------|-----------|
| 3D A-form   | kelvin_nu_factor_3d_cf(c, R) | (rho'/R)^2 | vanishes  |
| 3D Omega/H  | kelvin_mu_factor_3d_cf(c, R) | (R/rho')^2 | diverges  |

### Axisymmetric (r, z) -- 3D sphere viewed in meridional plane

| Formulation | Helper                                | Factor     |
|-------------|---------------------------------------|------------|
| Axisym A    | kelvin_nu_factor_axisym_cf(z_off, R)  | (rho'/R)^2 |
| Axisym Ω/H  | kelvin_mu_factor_axisym_cf(z_off, R)  | (R/rho')^2 |

**!! CRITICAL AXISYM CAVEAT:** NGSolve has NO built-in axisymmetric
mode. The caller MUST write the r-weight explicitly:

    a += nu_cf * grad(u) * grad(v) * r_coord * dx   # correct
    a += nu_cf * grad(u) * grad(v)           * dx   # WRONG, off by O(r)

The helper returns the sphere-inversion pullback only; it does NOT
absorb the `* r_coord` from the axisymmetric volume element. Dropping
this factor is one of the most common axisymmetric Kelvin bugs.

### 2D cylindrical Kelvin -- ANISOTROPIC tensor

Nagamine CEFC 2026 eq. 12: `nu' = diag(1, 1, (rho'/R)^4) nu`. The
Kelvin factor is NOT scalar -- it is a tensor, and which slot enters
the bilinear form depends on which **B / H component** the problem
actually uses (not on "A-form vs Omega-form").

| Case          | Helper                           | Factor     |
|---------------|----------------------------------|------------|
| In-plane B/H  | kelvin_factor_2d_inplane_cf()    | 1 (identity) |
| Axial B_z (nu)| kelvin_nu_factor_2d_axial_cf(o,R)| (rho'/R)^4 |
| Axial H_z (mu)| kelvin_mu_factor_2d_axial_cf(o,R)| (R/rho')^4 |

**In-plane is the common case** (most 2D magnetostatic setups): both
2D Omega-form with scalar phi (H in-plane) and 2D A-form with
A = A_z(x, y) z_hat (B in-plane, B_z = 0) land here -- Kelvin factor
is literally 1. Deprecated aliases `kelvin_mu_factor_2d_Hx_Hy_cf()`,
`kelvin_nu_factor_2d_Hx_Hy_cf()`, and `kelvin_nu_factor_2d_Az_cf(...)`
all resolve to `kelvin_factor_2d_inplane_cf()`.

The axial helpers are for the rare case where B or H genuinely has a
z-component (e.g. 2D A-form with in-plane A gives B_z z_hat).

## Material builder

`build_material_cf(mesh, base_value, kelvin_factor_cf, overrides=..., outer_keyword=...)`

Builds a `CoefficientFunction` indexed by `mesh.GetMaterials()`:
- Material name containing `outer_keyword` (case-insensitive) gets
  `base_value * kelvin_factor_cf`
- `overrides={material: value}` explicitly overrides specific materials
- Other materials get `base_value`

## Recipe: migrate inline Kelvin to the centralized helpers

### Before (inline, Ω-form 3D sphere with z-offset):

```python
r_prime_sq = x**2 + y**2 + (z - offset_z)**2
r_prime = sqrt(r_prime_sq + 1e-20)
mu_kelvin = (kelvin_radius / r_prime)**2 * mu0

mu_dict = {
    "air_inner": mu0,
    "air_outer": mu_kelvin,
    "magnetic": mu_r * mu0,
}
Mu = CoefficientFunction([mu_dict[m] for m in mesh.GetMaterials()])
```

### After (centralized helpers):

```python
from radia.kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf

mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=(0.0, 0.0, offset_z),
                                            R=kelvin_radius)
Mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)
# For post-processing (energy integration etc.):
mu_kelvin = mu0 * mu_kelvin_factor
```

### A-formulation axisym (Z-offset Kelvin):

```python
from radia.kelvin_source import kelvin_nu_factor_axisym_cf, build_material_cf

nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=z_offset, R=a)
nu_cf = build_material_cf(
    mesh, nu0, nu_kelvin_factor,
    overrides={m: nu0 / mu_r for m in mesh.GetMaterials() if "magnetic" in m.lower()},
)
```

### Rule

NEVER hand-write `(rho'/R)^2` or `(R/rho')^2` in new code. Always
import from `radia.kelvin_source`. This is mandatory for the
LLM-driven coil design workflow to avoid Kelvin-origin bugs.

If a Kelvin setup shows wrong inductance / field / energy, DO NOT
change the factor. Debug the FEM setup instead: GND, gauge
regularization, mesh near rho'=0, integration order, periodic ID.
"""


KELVIN_DIFFERENTIAL_FORMS = """
# Differential-form pullback under Kelvin inversion (0-, 1-, 2-form)

The Kelvin transformation is not just a change of variables; it is a
*coordinate transformation* that pulls back differential forms between
the physical exterior domain and the bounded computational outer
sphere. This makes the three field quantities — scalar potential Omega
(0-form), vector potential A (1-form), and magnetic induction B
(2-form pseudovector) — a very clean pedagogical example of how the
pullback factor depends only on the form degree.

## The three pullback rules

Let `r_phys` be a point in the physical exterior, `c` the Kelvin outer
sphere center, `R` the sphere radius, and `r' = c + R^2 (r_phys-c)/|r_phys-c|^2`
the Kelvin image in the computational domain. Define
`rho' = |r' - c|` and `n = (r' - c)/rho'`. Then:

| Form | Field | Pullback: F_comp(r') = ...                           |
|------|-------|-----------------------------------------------------|
| 0    | Ω     | Ω_phys(r_phys)                                      |
| 1    | A     | (R/rho')^2 * H(n) * A_phys(r_phys)                  |
| 2    | B     | -(R/rho')^4 * H(n) * B_phys(r_phys)                 |

Here `H(n) = I - 2 n n^T` is the Householder reflection. The magnitude
factor (R/rho')^{2k} for a k-form comes from k copies of the inverse
Kelvin Jacobian. The Householder comes from the non-identity part of
the Jacobian (for 1- and 2-forms only — 0-forms just transport). The
minus sign on B is from `det(H) = -1` via the Hodge identity for a
pseudovector.

The factors (R/rho')^{2k} DIVERGE as rho' -> 0 (Kelvin sphere center,
which is the image of physical infinity). That is consistent — the
computational frame compresses infinity to a point, so a finite
physical field must blow up in magnitude there (the integration
measure compensates).

## Practical eval helpers (radia.kelvin_source)

`eval_Omega_physical_from_gf(gf_Omega_comp, mesh, r_phys, c, R)`
    0-form. Evaluates the FEM scalar potential at a physical point
    r_phys outside the inner domain. Internally looks up the mesh at
    r' and returns the scalar value unchanged.

`eval_A_physical_from_gf(gf_A_comp, mesh, r_phys, c, R)`
    1-form. Does the (R/rho')^2 H pullback from computational A_comp
    back to physical A_phys. For use when the Kelvin exterior FEM
    result must be compared against a Biot-Savart / Radia reference
    at physical observation points.

`eval_B_physical_from_gf(gf_B_comp, mesh, r_phys, c, R)`
    2-form. -(R/rho')^4 H pullback. Same use case as above, at the
    B-field level; also useful for sanity-checking flux and energy.

Because the Kelvin map is involutive (phi(phi(r)) = r), the same
`kelvin_pullback_vector` / `kelvin_pullback_B_pseudovector` function
works for both directions: passing r_phys as the "r_prime" argument
yields the inverse pullback, which is exactly what these helpers do.

## Differential-geometry practice example

Given the above, the following are good self-checks for any reader
working through pullback derivations:

1. Verify on paper that `curl_r'(A_comp) = -(R/rho')^4 H curl_r(A_phys)`
   when A is defined as a 1-form and B = curl(A) as a 2-form. The
   factor (R/rho')^4 and the Householder with det -1 must match.

2. Check the 0-form case: `grad_r'(Omega_comp) = (R/rho')^2 H grad_r(Omega_phys)`
   — the gradient of a 0-form IS a 1-form, so the H + (R/rho')^2 factor
   appears here even though Omega itself has no pullback factor.

3. Verify the energy equivalence
     integral_{physical exterior} nu_0 |B|^2 dV
     = integral_{Kelvin sphere} (rho'/R)^2 nu_0 |B_comp|^2 dV'
   by combining the 2-form pullback (factor (R/rho')^4) with the
   volume element dV = (R/rho')^6 dV' and the nu material factor
   (rho'/R)^2 from Nagamine eq. 9. The (R/rho')^{4+6} / (R/rho')^2
   = (R/rho')^{12} on the LHS becomes (rho'/R)^{12} * (rho'/R)^2 on
   the RHS, which with |B_comp|^2 = (R/rho')^8 |B_phys|^2 gives back
   the physical LHS — confirming pullback + material + measure all
   agree.

These are the kind of exercises the Nagamine CEFC 2026 paper goes
through in its appendix and the Sugahara 2022 IEEE TransMag paper
uses to validate the two-sphere Kelvin implementation.

## Differential-geometry essence (H 1-form / B 2-form / mu Hodge)

The slide-level statement of the same fact, started from the field-
intensity side (H, B) rather than the potential side (Omega, A). Under
the sphere inversion r' = a^2/r:

| Object         | Invariant            | Pullback                       |
|----------------|----------------------|--------------------------------|
| H  (1-form) h  | r' h' = r h          | h' = (r/r') h = (a/r')^2 h     |
| B  (2-form) b  | r'^2 b' = r^2 b      | b' = (r/r')^2 b = (a/r')^4 b   |
| mu (Hodge)     | mu' = b'/h'          | mu' = (a/r')^2 mu              |
| energy 3-form  | invariant            | r'^3 mu' h'.h' = r^3 mu h.h    |

Two principles this language makes explicit -- the whole reason to write
Kelvin transformation in differential-form terms:

  1. **Differential forms are coordinate-transformation INVARIANT.**
     h and b are natural cochains; they pull back cleanly, and the field
     DOFs (Whitney edge/face) on the two identified spheres at r = a match
     1:1, with tangential-H / normal-B continuity automatic. The field
     representation needs no "fix-up" across the interface.

  2. **mu (the Hodge star) is the ONLY metric-dependent object.** All of
     the geometric distortion of the inversion is absorbed into the
     constitutive Hodge operator: mu' = (a/r')^2 mu (equivalently the
     reciprocal nu' = (r'/a)^2 nu). So "map infinity into a finite ball"
     reduces operationally to "solve the SAME problem with a position-
     dependent material mu'(r')" -- no far-field elements, no truncation.
     The energy 3-form invariance r'^3 mu' h'.h' = r^3 mu h.h guarantees
     the variational (FEM stiffness) problem is preserved exactly.

Consistency with the A/B/Omega table above: with a == R and r' == rho',
mu' = (a/r')^2 mu == (R/rho')^2 mu_0 (H1 / Omega side) and
nu' = (r'/a)^2 nu == (rho'/R)^2 nu_0 (HCurl / A side) -- the same statement.

References for this field-side view:
  - K. Sugahara, "Electromagnetic Analysis of Eddy Current Testing With
    Kelvin Transformation," IEEE Trans. Magn., 2022.  [lab origin; = Sugahara
    2022 IEEE TransMag 58(9) in the module header]
  - Q. Chen, "A Review of Finite Element Open Boundary Techniques for
    Static and Quasi-Static EM Field Problems," IEEE Trans. Magn., 1997.
  - Nagamine/Yamaguchi/Sugahara CEFC 2026 (pullback) + Sugahara 2022
    IEEE TransMag 58(9) -- see the convention header at top of module.
"""


KELVIN_SOURCE_IN_OMEGA_FORM = """
# Accelerator Magnet Analysis: Coil + Iron + Kelvin Open Boundary

## Quick Start (2 lines after mesh setup)

```python
from radia.scalar_potential_solver import ScalarPotentialSolver

solver = ScalarPotentialSolver(
    mesh,                              # NGSolve Mesh (iron + air + kelvin)
    iron_domains='iron',               # material name of iron yoke
    mu_r=1000,                         # relative permeability
    kelvin_region='air_outer',         # Kelvin exterior material name
    kelvin_radius=R,                   # inner sphere radius [m]
    kelvin_center=[cx, cy, cz])        # outer sphere center [m]

solver.set_source_from_radia(coil)     # rad.Fld -> H1 GridFunction (smooth)
solver.solve()                         # auto-selects total_reduced for Kelvin
B_cf = solver.get_B()                  # per-material B CoefficientFunction

# Evaluate B at any point (physical or Kelvin exterior)
mip = mesh(x, y, z)
Bx, By, Bz = B_cf(mip)[0], B_cf(mip)[1], B_cf(mip)[2]
```

What this solves:
- Coil-excited iron yoke (any mu_r) with open-boundary (Kelvin transform)
- Leakage field, field uniformity, shimming analysis
- Future: hysteresis (MatEnergyHysteresis) with same API

Accuracy: max 3%, RMS 1.5% vs Radia BEM reference (order=2, mu_r=100).
Validated by pytest: `tests/test_omega_reduced_omega.py` (5/5 PASS).

## Mesh Requirements

```
offset > 2 * R        # Kelvin spheres must NOT overlap geometrically
dirichlet = "GND"     # GND vertex at kelvin_center (image of infinity)
                       # Do NOT add iron_surf as Dirichlet (kills demag)
```

See `examples/kelvin_transformation/Omega_ReducedOmega/test_solver_integration.py`
for a complete working example with mesh construction.

## Method Details

When the Omega-reduced-Omega scheme is combined with the Kelvin outer
sphere, the coil source enters the FEM system ONLY through two boundary
terms on the total-reduced interface. The Kelvin (outer sphere) region
itself carries NO source term. The periodic BC between the inner and
outer sphere boundaries automatically propagates the solution.

Contrast with the A-formulation: there, the Biot-Savart source A_s
typically has to be injected (with a 1-form pullback) into the Kelvin
exterior. Omega (being a 0-form) vanishes at infinity, so the
Kelvin-image-of-infinity (rho' = 0) naturally gives GND Dirichlet.

**!! WARNING: the fully-reduced single-potential method
(`ScalarPotentialSolver.solve_single_potential`) does NOT work for
coil sources with Kelvin.** It requires H_s defined in the Kelvin
region via 1-form pullback, which `set_source_from_radia` does not
provide. Setting H_s=0 in Kelvin makes it worse (tested:
99% error vs 43% baseline). The CORRECT pattern is the EMPY-style
**total+reduced split** described below.

This pattern was validated via stone-bridge analysis (2026-04-17):

  1. mu_r=1 (air-only) 3D Kelvin vs Radia Biot-Savart: 0.39% max
  2. mu_r=100 3D Kelvin vs **2D axisymmetric FEM reference**: 1.15% max
     (proper h/p convergence, refined maxh=4mm order 3)
  3. 2D axisym reference vs Radia Biot-Savart (air): 0.02% agreement

**DO NOT** use Radia MMM/MSC as a convergence reference -- MMM showed
4-8% variation between mesh densities at iron endcap/surface probes,
masquerading as FEM error. The 2D axisym FEM reference
(`reference_2d_axisym.py`, u=r*A_phi formulation, order 4) is
sub-0.02% accuracy and is the authoritative reference.

**gap > 0 (gapped coil) is OUT OF SCOPE** for Omega-reduced-Omega.
ObjArcCur with arc < 360 deg has no return path, so div(J) != 0 at
the terminals and Omega_s is uncontrolled near the cut surface.  The
12% x-axis error at gap=5 deg is a formulation limitation, not
a discretization issue.

Scripts:
  `validate_omega_coil_source_v2.py` (3D Omega + Kelvin)
  `reference_2d_axisym.py` (2D axisym A_phi reference)

## Recipe (Radia + NGSolve)

Given: coil object `coil` = rad.ObjArcCur, and a Kelvin mesh with
three materials ("iron", "air_inner", "air_outer") plus
`iron_surf` = iron-air interface, Kelvin outer sphere offset
("kelvin_ext") identified periodically with "kelvin_int".

```python
from radia.kelvin_source import kelvin_mu_factor_3d_cf

# (1) Source fields from Radia → H1 GridFunction (smooth, C0).
#     rad.Fld('phi') returns Phi where H = -grad(Phi).
#     Omega_s = -Phi (EMPY convention: H = +grad(Omega_s)).
#     rad.Fld('h') gives H_s, B_s = mu0 * H_s.
#
#     HIGH-LEVEL API (recommended):
solver = ScalarPotentialSolver(mesh, iron_domains='iron', mu_r=100,
    kelvin_region='air_outer', kelvin_radius=R, kelvin_center=c)
solver.set_source_from_radia(coil)   # builds Omega_s + H_s as GridFunctions
solver.solve()                        # auto-selects total_reduced
B = solver.get_B()
#
#     MANUAL (for custom setups):
Omega_s_cf = ...  # H1 GridFunction from -rad.Fld('phi')
Bs_cf = ...       # mu0 * H1 GridFunction from rad.Fld('h')

# (2) Kelvin material (Nagamine canonical):
mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=c_out, R=R_K)
Mu = CF([mat_dict[m] for m in mesh.GetMaterials()])
# mat_dict = {"iron": mu_r*MU_0, "air_inner": MU_0,
#             "air_outer": MU_0 * mu_kelvin_factor}

# (3) Bilinear form — Mu * grad . grad over ALL volume:
#     CRITICAL: dirichlet="GND" ONLY.  Do NOT add "iron_surf"
#     as Dirichlet — that pins omega to Omega_s and kills
#     demagnetization (gives Bt = mu_iron * Hs, the long-cylinder
#     limit). The iron_surf BC must enter WEAKLY via the RHS lift.
fes_raw = H1(mesh, order=fe_order, dirichlet="GND")
fes = Periodic(fes_raw)
omega, psi = fes.TnT()
a = BilinearForm(fes)
a += Mu * grad(omega) * grad(psi) * dx
a.Assemble()

# (4) Dirichlet lift: Set Omega_s on iron_surf BND (NOT Dirichlet).
#     This creates a GridFunction with Omega_s at iron-air boundary
#     nodes and 0 elsewhere, used as RHS source contribution.
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s_cf, BND, mesh.Boundaries("iron_surf"))

# (5) RHS: reduced-region lift + Neumann B_s.n on iron surface.
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
normal = -specialcf.normal(mesh.dim)   # outward from iron
f += (normal * Bs_cf) * psi * ds("iron_surf")
f.Assemble()

# (6) Solve directly (NO inhomogeneous Dirichlet correction):
sol = GridFunction(fes)
sol.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

# (7) Post-process: B = Mu * grad(sol) in iron, Mu*(grad(sol)-grad(lift))+Bs in air
```

## The three source-handling rules for Kelvin

| Formulation            | Inner reduced | Inner total   | Kelvin       |
|------------------------|---------------|---------------|--------------|
| Omega-reduced-Omega    | Omega_s lift  | B_s.n (Neum.) | 0 (periodic) |
| A-reduced-A            | A_s (curl)    | A_s.n (Neum.) | pullback(A_s)|
| H-form with source     | H_s = B_s/mu0 | --            | pullback(H_s)|

The Omega scheme gets away WITHOUT a Kelvin-region source because it
reduces to solving a pure Laplace problem in air+Kelvin with the source
factored into BCs.

## Radia API cheat-sheet for source evaluation

```python
rad.Fld(obj, 'phi', pt)   # scalar magnetic potential — **MAGNETIZED BODIES ONLY**
rad.Fld(obj, 'b',   pt)   # magnetic flux density (any source)
rad.Fld(obj, 'h',   pt)   # H field               (any source)
rad.Fld(obj, 'a',   pt)   # vector potential      (any source)
```

**!! CRITICAL GOTCHA: rad.Fld(..., 'phi', pt) does NOT return the
magnetic scalar potential for current-driven coils (ObjArcCur).** The
C++ `rad_arc_current.cpp` uses an incorrect solid-angle integrand for
arc coils (not full-loop); the returned phi has `grad(phi) != H`
(ratio 19-32x empirically at multiple probe points). For magnetized
bodies (ObjRecMag, ObjCylMag with M), `H = -grad(phi)` holds.

Workarounds for a current-driven coil source:

1. **Line-integration of rad.Fld('h')**: compute Omega_s(p) =
   integral_ref^p H_s . dl on a voxel grid. Choose reference point
   below iron; straight path must not cross coil cut surface. Build
   VoxelCoefficient. See `validate_omega_coil_source_v2.py` for the
   working implementation.

2. **Path-integration helper**:
   `radia.bem_sibc_solver.compute_phi_inc_from_filaments(obs_points,
   filament_paths, currents)` computes Omega_s for arbitrary filament
   bundles. Supports complex per-filament currents.

3. **Equivalent-magnetization coil**: model a solenoid as hollow
   ObjRecMag blocks with M_z = N*I/h. Then rad.Fld('phi') works.

## Critical geometry constraints (5 hard-won bugs)

1. **Kelvin sphere offset > 2*R required.** Inner sphere @origin R
   and outer sphere @offset R: if offset < 2R, the spheres overlap
   geometrically and OCC Glue assigns iron tets to the air_outer
   material. This silently gives wrong Mu inside iron.

2. **Bs voxel bbox must cover entire air_inner region.** A bbox
   covering only the iron surface causes VoxelCoefficient boundary
   clamping: constant |B|_FEM beyond the bbox edge.

3. **dirichlet="GND" only, NOT "GND|iron_surf".** Strict Dirichlet
   on iron_surf pins omega to harmonic-extension(Omega_s) in iron,
   giving grad(omega) = Hs and Bt = mu_iron * Hs (long-cylinder
   limit with no demagnetization). Weak BC via RHS lift allows
   demagnetization.

4. **Ot.Set(sol, VOL, definedon="iron") may return near-zero** when
   projecting from full-mesh to iron-restricted H1. Use per-material
   CF dispatch `CF([Mu*grad(sol) if mat=="iron" else ...])` instead.

5. **Source must live in the INNER physical region** (not the Kelvin
   exterior).

## Gotchas

- **Omega_s is multi-valued through a coil-spanning surface.** For
  closed-loop coils (gap_deg = 0, the ONLY valid case for
  Omega-reduced-Omega), arrange the FEM domain so no mesh boundary
  crosses the branch-cut surface of Omega_s.
- **Gapped coils (gap_deg > 0) are NOT supported** by the Omega
  formulation (div(J) != 0 at terminals makes Omega_s
  path-dependent). Use the A-formulation (HCurl) for gapped coils.
- **Kelvin region gets ZERO source** but the periodic BC must
  actually reduce DOFs (check FreeDofs before/after Periodic).
- **GND Dirichlet** at the Kelvin outer sphere center is the image
  of physical infinity and enforces Omega(infinity) = 0.
"""

KELVIN_RADIA_SOURCE_HELPERS = """
# Radia source field evaluation with Kelvin pullback

Student-friendly helpers in ``radia.kelvin_source`` for evaluating
Radia source fields (coils, magnets) at arbitrary points, including
points inside the Kelvin exterior mesh region.

## Quick reference

```python
from radia.kelvin_source import (
    eval_H_from_radia_at_physical,     # H at physical point (no Kelvin)
    eval_B_from_radia_at_physical,     # B at physical point (no Kelvin)
    eval_H_from_radia_in_kelvin,       # H at Kelvin mesh point (with pullback)
    eval_B_from_radia_in_kelvin,       # B at Kelvin mesh point (with pullback)
    eval_field_from_radia_with_kelvin, # auto-dispatch (physical or Kelvin)
)

# Auto-dispatch: "just give me the field" (B or H)
B = eval_field_from_radia_with_kelvin(coil, points, center, R, field_type='b')

# Explicit Kelvin-region evaluation (e.g. for source injection)
H_kelvin = eval_H_from_radia_in_kelvin(coil, [0.4, 0, 0],
                                        center=[0.5, 0, 0], R=0.18)
```

## Differential-form pullback rules (applied automatically)

| Field | Form degree | Pullback factor |
|-------|-------------|-----------------|
| Omega | 0-form      | 1 (trivial)     |
| H, A  | 1-form      | (R/rho')^2 H    |
| B     | 2-form      | -(R/rho')^4 H   |

where H is the Householder reflection about the radial direction.

## FEM GridFunction evaluation at physical exterior points

```python
from radia.kelvin_source import (
    eval_Omega_physical_from_gf,  # 0-form: no Jacobian
    eval_A_physical_from_gf,      # 1-form: (R/rho')^2 H pullback
    eval_B_physical_from_gf,      # 2-form: -(R/rho')^4 H pullback
)

# After FEM solve, evaluate B at a physical point outside the inner sphere
B_phys = eval_B_physical_from_gf(gf_B, mesh, r_phys=[0, 0, 1.0],
                                  center=[0, 0, 4], R=1.5)
```
"""


KELVIN_BENCHMARK_PANEL = """
# Kelvin Benchmark Panel Mode (radia_em.py "Kelvin Benchmark")

End-user verification path for the Cubit-meshed Kelvin pipeline.
The EM panel ships with a 4th formulation -- "Kelvin Benchmark" --
that runs `calc_kelvin_benchmark.py` on a `.vol` with a magnetic
sphere in uniform external Hz and compares the interior field to
the analytical solution `Hz_inside = 3 / (mu_r + 2) * H0`.

## Required `.vol` structure

The `.vol` must declare:
- **blocks**: `magnetic`, `air`, `kelvin`
- **sidesets**: `sphere` (mag-air interface), `kelvin_int`, `kelvin_ext`
- **bbnodeset**: `GND` (Kelvin sphere centre, Dirichlet anchor)

## Bundled samples (verified 2026-04-26)

Two `.vol` samples ship with the wheel for the panel's Browse dialog:

  | sample (frac) | reduction                                | offset_dir | error @ p=2 |
  |---------------|------------------------------------------|------------|------------:|
  | `kelvin_benchmark_sphere_1_2.vol` | {"y":"bn=0"}                             | z          | +1.07%      |
  | `kelvin_benchmark_sphere_1_4.vol` | {"x":"bn=0","y":"bn=0"}                  | z          | +0.71%      |

The 1/8 variant (`{"x":"bn=0","y":"bn=0","z":"ht=0"}`, offset along
x) is intentionally NOT shipped: the kelvin_far Dirichlet plane
passes through r' = 0 (Kelvin centre), where the reluctivity
Mu = mu0*(R/r')^2 is singular, and the linear solve gives Hz ~ 0
(~-101% error) instead of analytical.  The build script
`panels/samples/kelvin_benchmark_sphere_1_8_build.py` is kept for
research; the .vol is omitted until the singularity is fixed.

Build scripts:
- `panels/samples/kelvin_benchmark_sphere_1_2_build.py` (thin wrapper)
- `panels/samples/kelvin_benchmark_sphere_1_4_build.py` (thin wrapper)
- `panels/samples/kelvin_benchmark_sphere_1_8_build.py` (research-only)
- `panels/samples/kelvin_benchmark_sphere_build.py` (parameterised core,
  `--frac 1_2|1_4|1_8`)

All four scripts mu_r=100, R_air=0.20 m, Kelvin offset 3*R = 0.60 m.

## Two non-obvious Cubit fixes embedded in the build core

1. `subtract A from B keep` does NOT carve B in Cubit 2025.3.
   Workaround: drop `keep`, then re-create A as a fresh primitive.

2. Copy-mesh anchor curve selection for 1/8 octant caps was
   non-deterministic (3 equal-length quarter-arcs). Fix in
   `_add_kelvin_cubit_reduction`: pick `min(curves,
   key=(centroid_z, centroid_y, centroid_x))` for translation-
   equivalent source / target anchors.  See
   `memory/feedback_kelvin_1_8_blocker.md` (mesh fixed; 1/8 numerical
   solver still has a kelvin_far singularity issue, deferred).

## CLI

    python calc_kelvin_benchmark.py \\
        --vol model.vol \\
        --fes-order 2 \\
        --mu-r 100 --H0 1.0 --field-axis z \\
        --R-kelvin 0.20

JSON output: `Hi_origin`, `Hi_analytical`, `error_pct`, `slaved_dofs`.

## p-convergence (1/4 sphere, validated 2026-04-25)

  | order | typical error |
  |-------|--------------:|
  | p=1   |        +7.79% |
  | p=2   |        +0.71% |
  | p=3   |        +0.54% |

Panel-level golden test:
`tests/panels/test_kelvin_benchmark_golden.py` parametrises over
both 1/2 and 1/4 samples, locking |error| < 1.5% at p=2.

## Why 1/8 is unsupported for the sphere benchmark (physics)

The benchmark BVP is "magnetic sphere in **uniform Hz**".  The z=0
mirror plane reverses the source vector `H0 z_hat -> -H0 z_hat`,
so the problem is **not 1/8-symmetric**.  No amount of BC tuning on
the z=0 cut can recover the analytical solution -- the source itself
breaks the third symmetry.  1/2 (split y) and 1/4 (split x and y) are
the maximum reductions for this BVP and are the only samples shipped.

Do NOT conflate this with the **EM panel (C-type magnet) 1/8 case**:
that BVP has no source vector; the iron yoke + coil ampere-turns are
fully octant-symmetric, so the 1/8 reduction is geometric only and
ships in the EM panel via `panels/samples/em/em_1-8_eighth.jou` +
`tests/panels/golden/em_eighth_mu1000.json` (ELF CEFC 2020 convention
`-x-y+z`, opposite sign pattern from the sphere reductions).

## Surprise: rho_min sweep up to R (2026-04-26)

The Kelvin reluctivity Mu = mu_0 * (R/rho')^2 is bounded by clamping
rho' >= rho_min in the build script.  Sweeping rho_min from a small
value all the way up to R (which collapses Mu to a uniform mu_0
everywhere in the Kelvin region) STILL gives 1/2 +0.34% / 1/4 -0.02%
at p=2.  i.e. for compact geometry (Kelvin offset = 3*R), Periodic
identification on (kelvin_int <-> kelvin_ext) plus the symmetric
Dirichlet/Neumann BCs do MOST of the open-boundary work; the
(R/rho')^2 reluctivity correction is a small refinement.

**Practical debugging shortcut**: if a Kelvin solver gives wildly
wrong numbers, re-run with rho_min = R (Mu collapses to mu_0).  If
the answer becomes correct, the bug is in the Mu coefficient (rho_min
clamp / Kelvin centre detection / domain assignment).  If the answer
is still wrong, the bug is in BCs / Periodic identification / mesh.
This isolates Kelvin-reluctivity issues from BC issues in one solve.

## Cubit-meshed Kelvin needs `-specialcf.normal` (Neumann correction)

The reduced-Omega Neumann BC `H_t.n` correction term has an extra
sign flip on Cubit-meshed Kelvin compared to the OCC reference.
Cubit assigns surface element normals with the opposite sign
convention to NGSolve's WorkPlane-based OCC builder, so the line:

    f += -H_s * specialcf.normal(3) * v.Trace() * ds("kelvin_int")

becomes correct only when the sign matches the OCC reference run.
On Cubit-meshed `.vol`, use `-specialcf.normal(3)` (i.e. the term
above as-written, NOT `+specialcf.normal(3)`).  Verified by a
sign-flip A/B test that drove the 1/4 sample from -7% to +0.71% at
p=2.  The benchmark calc script (`calc_kelvin_benchmark.py`) hard-
codes the Cubit-correct sign; the OCC examples in
`examples/kelvin_transformation/` use the same expression because
the sign convention is now unified across both mesh sources.
"""


KELVIN_KAMEARI_CANONICAL = """
# Kameari Canonical Pattern: Boundary-Integral Source Injection (2026-05-05)

Reference: Kameari (2025/10/14 slides), W:/00_CAE/NGSolve/亀有/4-05.pdf,
"Electromagnetic Analyses Using Higher Order Hierarchic Finite Elements".

## Key Insight

Kelvin transformation in NGSolve is a PROVEN, HIGH-ACCURACY technique when
properly formulated. The critical distinction:

- WRONG (historical v11-v14): "(nu - nu_0) bulk source" reduced-A form.
  The (nu - nu_0) simplification requires A_s to satisfy nu_0 Maxwell
  globally, but the Kelvin pullback satisfies the metric-dependent nu'
  Maxwell instead. This caused +43% inductance error.

- CORRECT (Kameari canonical): Source coupling via BOUNDARY INTEGRAL on the
  inner-Kelvin interface dOmega, NOT via volume integral with (nu - nu_0).

## Three Canonical Reductions (Kameari slide 4)

| Method   | Inner region      | Outer (Kelvin) region  | Use case                     |
|----------|-------------------|------------------------|------------------------------|
| Omega-Omega_r | H = -grad(Omega_t) | H = -grad(Omega_r) + H_s | Linear magnetic, scalar |
| A-Omega_r | B = curl(A_t)   | H = -grad(Omega_r) + H_s | ** Recommended: sphere/cuboid |
| A-A_r    | B = curl(A_t)     | B = curl(A_r + A_s)    | Vector source / coil         |
| A-phi-A_r | A_t, phi in cond | A_r in air/Kelvin      | Eddy current (TEAM Problem 7)|

## Weak Forms (Kameari slides 7-9)

```python
# A-Omega_r (slide 8): vector A in conductor/inner, scalar Omega_r in Kelvin
# Bulk terms:
a += 1/mu * curl(N) * curl(A) * dx           # inner + Kelvin regions
a += mu * grad(omega) * grad(Omega_r) * dx   # Kelvin region only

# Source coupling via BOUNDARY INTEGRAL at inner-Kelvin interface dOmega:
f += N.Trace() * Cross(H_s, n) * ds("kelvin_int")   # A side
f += -omega.Trace() * (B_s * n) * ds("kelvin_int")  # Omega_r side
```

The applied field (H_s or B_s) enters via a boundary integral on the interface
between inner physical region and Kelvin region -- NOT in the bulk volume.

## Performance (Kameari slide 17)

Magnetic sphere a=1m, mu_r=1000, B_0=1T, analytical Bz0=2.9940 T, COARSE mesh:

| Method     | Order | Bz0    | Error  |
|------------|-------|--------|--------|
| Omega-Omega_r | 2  | 3.3813 | 12.9%  |
| Omega-Omega_r | 3  | 2.9995 | 0.184% |
| Omega-Omega_r | 4  | 2.9985 | 0.029% |
| A-A_r      | 3     | 2.9928 | 0.040% |
| A-Omega_r  | 3     | 2.9940 | **0.001%** (best) |

A-Omega_r Order 3 achieves sub-0.01% error on a coarse mesh.

## Independence from Kelvin Radius rk (Kameari slide 18)

rk = 100 (very large Kelvin sphere) gives the SAME accuracy as small rk
with adaptive refinement. The Kelvin region size is a FREE parameter --
no tuning needed.

## Practical Recommendation

For "isolated conductor + applied uniform B" problems:
- Use Kameari A-Omega_r with Order >= 3 and boundary-integral source
- Do NOT use (nu - nu_0) bulk source (the "reduced-A" pattern from pre-2026 code)
- The boundary integral at kelvin_int replaces all volume-source terms

For eddy current with external conductors (TEAM Problem 7):
- Use A-phi-A_r (slide 25): A + phi inside conductor, A_r in Kelvin region
- Proved in TEAM Workshop Problem 7 computations

## NGSolve Implementation Sketch (A-Omega_r)

```python
from ngsolve import *

# Two FE spaces: HCurl inside + Kelvin, H1 (Kelvin only)
fesA = HCurl(mesh, order=3, dirichlet="<conductor_bcs>")
fesO = H1(mesh, order=3, definedon=mesh.Materials("kelvin"),
          dirichlet_bbnd="GND")
fes = fesA * fesO
(A, Omega_r), (N, omega) = fes.TnT()

mu_cf = mesh.MaterialCF({"inner": mu_0, "kelvin": (R/rho_p)**2 * mu_0})
nu_cf = mesh.MaterialCF({"inner": nu_0, "kelvin": (rho_p/R)**2 * nu_0})

a = BilinearForm(fes)
a += 1/mu_cf * curl(N) * curl(A) * dx
a += nu_cf * grad(omega) * grad(Omega_r) * dx("kelvin")

f = LinearForm(fes)
# Boundary-integral source (NOT bulk (nu - nu_0) term)
f += N.Trace() * Cross(H_s_vec, specialcf.normal(3)) * ds("kelvin_int")
f += -omega.Trace() * (B_s_vec * specialcf.normal(3)) * ds("kelvin_int")
```

See docs/kelvin/KELVIN_TRANSFORMATION.md §7.4.1 for full derivation.
See docs/cln/CAUER_LADDER_NETWORK.md §6.4.1 for CLN + Kelvin specifics.
"""

KELVIN_VERIFIED_RECIPE = """
# Verified Production Recipe: A-formulation + Periodic Kelvin + BDDC

This section consolidates the recipe that production code in
`src/radia/panels/calc_fem_kelvin.py` has used since 2026-04 for the
combined "HCurl + Periodic Kelvin + iterative solver" problem.  If
you are debugging a new HCurl eddy-current script and the result is
off by O(10x) from the analytical / BEM reference, you are almost
certainly missing one of the five elements below.

## The five mandatory elements

1. **`nograds=True`** on the HCurl FE space.
   - HCurl basis includes gradients of H1 basis (curl(grad phi) = 0).
   - These form the curl-curl kernel = gauge freedom.
   - Without `nograds=True`, BDDC factorization hits singular blocks
     (`ZGETRF::info = k` for some k > 0) and a direct solver returns
     an arbitrary gradient-polluted solution.
   - Required at ANY polynomial order (including p=1).

2. **`Periodic(base_fes)` after `nograds=True`**, not before.
   - The order is: `HCurl(..., nograds=True, ...)` THEN `Periodic(...)`.
   - Periodic identification couples slave DOFs to master DOFs; you
     want the master/slave pairing on the already-gradient-stripped
     basis.

3. **Gauge regularization restricted to non-Kelvin materials**.
   - `nograds=True` removes the H1 gradient subspace, but the
     residual harmonic cohomology (nontrivial topology, in particular
     the period gluing) can still leave curl-curl singular.
   - Add `reg * NU_0 * u * v * dx(non_kelvin_mats)` to fix the kernel.
   - **DO NOT** apply this term in the Kelvin domain: `nu_kelvin =
     (r'/R)^2 * nu0` vanishes near r'=0, so a constant `NU_0` mass
     term would dominate there and corrupt the magnetic-energy
     integral.
   - Typical magnitude: `reg = 1e-6`.

4. **`bonus_intorder=4` on the curl-curl integral**.
   - The Kelvin-modulated `nu_cf = (r'/R)^2 * nu_0` is a rational
     function of coordinates; standard quadrature undertcounts the
     energy in the exterior shell by several percent.
   - `dx(bonus_intorder=4)` raises the quadrature order to handle the
     rational integrand cleanly.

5. **BDDC preconditioner for order >= 2; CompactAMS only at order=1**.
   - `Preconditioner(a, "bddc")` works at any order with `nograds=True`.
   - `radia.sparsesolv_ngsolve.ComplexCompactAMSPreconditioner` is
     order-1-only AND incompatible with Periodic Kelvin (DOF
     dimensions don't match the auxiliary H1 space).
   - At p=1 without Periodic: prefer CompactAMS (cheaper).
   - At p>=2 OR with Periodic Kelvin: BDDC.

## Copy-paste skeleton

```python
from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, Integrate,
                     curl, x, y, z, dx, ds, InnerProduct,
                     Preconditioner, CGSolver)

MU_0 = 4e-7 * 3.141592653589793
NU_0 = 1.0 / MU_0

mesh = Mesh("model.vol")
omega = 2 * 3.141592653589793 * freq

# 1. Material-dependent nu_cf with Kelvin modulation
materials = mesh.GetMaterials()
non_kelvin_mats = [m for m in materials if "kelvin" not in m.lower()]
nu_list = []
sigma_list = []
for m in materials:
    if "kelvin" in m.lower():
        # Sugahara 2022 canonical: nu' = (r'/R)^2 * nu_0
        nu_list.append((x*x + y*y + z*z) / R_K**2 * NU_0)
        sigma_list.append(0.0)
    elif m == "conductor_name":
        nu_list.append(NU_0)
        sigma_list.append(sigma_cond)
    else:  # air
        nu_list.append(NU_0)
        sigma_list.append(0.0)
nu_cf    = CoefficientFunction(nu_list)
sigma_cf = CoefficientFunction(sigma_list)

# 2. HCurl FE space with nograds=True, then Periodic wrap
base_fes = HCurl(mesh, order=p, complex=True, nograds=True)
fes      = Periodic(base_fes)
u, v     = fes.TnT()

# 3. Bilinear form: curl-curl + AC mass (conductors only) + gauge reg
a = BilinearForm(fes, symmetric=True)
a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
a += 1j * omega * sigma_cf * InnerProduct(u, v) * dx("conductor_name")
if non_kelvin_mats:
    a += 1e-6 * NU_0 * InnerProduct(u, v) * dx("|".join(non_kelvin_mats))

# 4. RHS (your applied source / Biot-Savart background)
f = LinearForm(fes)
# ... f += -1j * omega * sigma_cf * InnerProduct(A_ext, v) * dx("conductor_name")

# 5. BDDC + CG
prec = Preconditioner(a, "bddc")
a.Assemble(); f.Assemble()
gfu = GridFunction(fes)
inv = CGSolver(a.mat, prec.mat, complex=True, maxsteps=2000,
               precision=1e-8, printrates=False)
gfu.vec.data = inv * f.vec
```

## Verify-First diagnostic (BEFORE you run a physics solve)

Per CLAUDE.md "Verify-First Policy", run these checks first.  If any
of them fails, the issue is in geometry / Identify / labels --
debugging at the physics-solve layer is wasted effort.

```python
# (a) Materials / boundaries spelled as expected
print(mesh.GetMaterials(), mesh.GetBoundaries())

# (b) Periodic actually constrains DOFs
fb = HCurl(mesh, order=p, complex=True, nograds=True)
fp = Periodic(fb)
print("slaved:", sum(fb.FreeDofs()) - sum(fp.FreeDofs()))  # must be > 0

# (c) Functional boundary test for Periodic
gfu_test = GridFunction(fp); gfu_test.vec[:] = 0
gfu_test.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
r = Integrate(gfu_test * gfu_test, mesh, definedon=mesh.Boundaries("kelvin_ext")) \\
  / Integrate(gfu_test * gfu_test, mesh, definedon=mesh.Boundaries("kelvin_int"))
print("ratio:", r)  # must be 1.0
```

## Failure-mode -> root-cause table

| Symptom                                          | Most likely cause           |
|--------------------------------------------------|-----------------------------|
| BDDC: `ZGETRF::info > 0`                         | Missing `nograds=True`      |
| CG diverges / stalls at iter ~ few               | Missing `nograds=True` OR missing gauge reg |
| Result O(10x) too small in magnitude             | `reg * NU_0 * u * v * dx` applied in Kelvin domain (`nu->0` there + constant mass = corrupted energy) |
| Inductance/energy off ~5%, mesh-independent      | Missing `bonus_intorder=4`  |
| `Periodic` ratio test fails (not 1.0)            | Identify() missing or wrong tolerance; fix the .vol, not the solver |
| Inductance grows as you refine the mesh          | Pure-gradient pollution surviving in solution; `nograds=True` + reg both required |

## Common AI-debugging mistake (recorded 2026-05-20)

Symptom: re-run the script with progressively finer mesh / higher
order / Kelvin-on / iterative solver, and the discrepancy gets
*worse*, not better.  The instinct is "more refinement should be
better"; but each refinement amplifies the gauge null-space
pollution if `nograds=True` is missing.  A naive HCurl direct-solve
on a coarse mesh can be accidentally close to the truth because the
pseudo-inverse picks a small-norm gradient pollution; refining the
mesh enlarges the kernel and the pollution.

**Fix order**: (1) `nograds=True`, (2) gauge reg restricted to
non-Kelvin, (3) bonus_intorder=4, (4) THEN refine mesh / raise order.

## References

- `src/radia/panels/calc_fem_kelvin.py` (production reference for IH)
- `src/radia/panels/calc_fem_coilmesh.py` (production reference for
  volumetric coil + workpiece SIBC + Kelvin)
- CLAUDE.md "Verify-First Policy" (2026-04-25): FES inspection
  before physics solve.
- NGSolve Maxwell tutorial unit-2.4: `HCurl(..., nograds=True)` is
  the documented gauge-fixing recipe at ANY polynomial order.
"""


KELVIN_MESH_CONTROL = """
# Kelvin Mesh Control -- where to spend elements (and where NOT)

Kelvin meshing is counterintuitive relative to standard FEM ("refine until the
field converges").  Two distinct things are under your control: (1) a HARD
constraint -- the two spheres must share an IDENTICAL surface mesh on Gamma --
and (2) the element BUDGET, where the verified result is the OPPOSITE of the
usual instinct: do NOT refine the exterior air.

Provenance: paper kelvin_dtn.tex S4.2-4.4; topic `dtn_coarse_mesh`
(subsec:gammaonly result + optimal-R analysis); topics `identify` / `adaptive`
below; CLAUDE.md "Verify-First Policy" and "AI-Driven Cubit: Probe, Don't Guess".

## 1. HARD constraint: identical surface mesh on Gamma (the gluing condition)

The inner physical ball and the Kelvin (inverted exterior) ball are glued at the
truncation sphere Gamma (radius R).  Their triangulations on Gamma must be
IDENTICAL point-to-point, or the periodic / DtN coupling does not form and the
solve crashes.

  - Cubit: `copy mesh surface <in> onto surface <out> ...` (MANDATORY).  The
    identified nodes are related by a pure TRANSLATION (offset), so source and
    target hemispheres must have matching topology.  Mesh the air VOLUMES first
    (so the surface is triangulated) BEFORE copy mesh.
  - OCC / NGSolve: `int_face.Identify(ext_face, "periodic",
    IdentificationType.PERIODIC)` then `GenerateMesh()` conforms automatically.
  - Pitfalls:
      * ACIS sphere has 0 curves -> copy mesh has no anchor; `webcut ... with
        plane zplane` creates 1 curve per hemisphere to anchor on.
      * Anchor selection: a HEMISPHERE has a unique longest curve (equator), so
        `max(curve, key=length)` is fine.  A 1/8 OCTANT has THREE equal-length
        arcs -> `max(key=length)` ties NON-deterministically and can pick a
        different arc on source vs target.  Use a deterministic geometric
        tiebreak `min(curves, key=(centroid_z, y, x))` (see
        _add_kelvin_cubit_reduction; feedback_kelvin_1_8_blocker).
      * Never tet-mesh the two balls independently -- you get different
        triangulations and the copy/identify fails.
  - Verify-First (BEFORE any solve, sub-second; CLAUDE.md):
        slaved = sum(fes.FreeDofs()) - sum(Periodic(fes).FreeDofs())   # > 0
        # Gamma functional test: Set(1) on kelvin_int, integrate on kelvin_ext;
        # the kelvin_ext/kelvin_int ratio must be 1.0
    If slaved == 0 or the ratio != 1.0 the fault is geometry / labels /
    Identify -- fix it before trusting any solved number.

## 2. The element BUDGET (verified, counterintuitive): do NOT refine the exterior

At p >= n the open-boundary accuracy is decided by the Gamma SURFACE mesh and
the Curve (geometry) order ONLY; the interior VOLUME mesh of BOTH balls is
irrelevant (||u_h - P_n|| ~ 1e-15; paper S4.2 subsec:gammaonly).  Measured:
Delta-DoF ~ 58 buys an open-boundary error ~ 1/45 of the interior FEM error.

  => The exterior Kelvin ball can be a handful of COARSE elements.  Shrinking
     maxh in the air is wasted budget -- it does not improve the open boundary.

## 3. The four knobs that DO matter

  | knob                | guidance                                             |
  |---------------------|------------------------------------------------------|
  | p (poly order)      | p >= n for the highest excited multipole n (the      |
  |                     | Kelvin image is a degree-n polynomial -> exact at    |
  |                     | p >= n; the eigenvalue ladder -(n+1)/R).             |
  | Curve order         | sets the accuracy FLOOR (5-6 digits on a coarse      |
  |  (geometry)         | mesh); the sphere must be curved -- raise it to      |
  |                     | lower the floor.                                     |
  | R / centering       | R/a ~ 3 optimal for COMPACT bodies (the spectral     |
  |                     | basis for the classic "air box 2-5x" rule); the      |
  |                     | sphere is R-independent (pure dipole).  Center on    |
  |                     | the device, take the minimal enclosing R.            |
  | hp at device corner | a re-entrant corner (L 270deg, lambda=2/3) caps       |
  |                     | algebraic convergence (alpha_h~0.357, alpha_p~0.661);|
  |                     | the rate limiter is INSIDE the device, not the air   |
  |                     | -- geometric grading TO the corner restores          |
  |                     | exponential convergence.                             |

## 4. Radial grading, the inversion, and the infinite point

  - Inner physical ball: refine toward the device, coarsen toward Gamma
    (`GenerateMesh(maxh=..., grading=0.7)`).
  - Kelvin ball: the inversion r -> R^2/r makes a uniform Kelvin-ball element
    physically huge in the far field -- which is FINE (the far field is smooth /
    low-multipole, exactly what coarse elements resolve).  Do not spend grading
    effort in the Kelvin ball.
  - Infinite point (Kelvin center rho'=0): the material weight ~(R/rho')^2
    diverges at the center but is INTEGRABLE and is pinned by the GND Dirichlet
    there, so a nodal SCALAR (Omega / phi) FEM is benign.  Only an EDGE-element
    A-formulation sees a singular tensor material at the center -- then use a
    center-less DtN (FEM-BEM) that never creates the center node.

## 5. Adaptive refinement

A Zienkiewicz-Zhu estimator + SetRefinementFlag / Refine loop is available
(topic `adaptive`), but per section 2 the refinement must target the device
INTERIOR / corners, never the exterior air.

## One-line rule

Make Gamma conforming (copy mesh / Identify), then spend the budget on the Gamma
SURFACE mesh + Curve order + p + device-corner hp -- and leave the exterior air
coarse.
"""


def get_kelvin_documentation(topic: str = "all") -> str:
    """Return Kelvin transformation documentation by topic."""
    topics = {
        "overview": KELVIN_OVERVIEW,
        "h_formulation": KELVIN_H_FORMULATION,
        "a_formulation": KELVIN_A_FORMULATION,
        "3d": KELVIN_3D,
        "hcurl_3d": KELVIN_HCURL_3D,
        "verified_recipe": KELVIN_VERIFIED_RECIPE,
        "adaptive": KELVIN_ADAPTIVE,
        "identify": KELVIN_IDENTIFY,
        "mesh_control": KELVIN_MESH_CONTROL,
        "tips": KELVIN_TIPS,
        "verification": KELVIN_VERIFICATION,
        "periodic_wedge": KELVIN_PERIODIC_WEDGE,
        "robustness": KELVIN_ROBUSTNESS,
        "helpers_recipe": KELVIN_HELPERS_RECIPE,
        "differential_forms": KELVIN_DIFFERENTIAL_FORMS,
        "radia_source_helpers": KELVIN_RADIA_SOURCE_HELPERS,
        "source_in_omega_form": KELVIN_SOURCE_IN_OMEGA_FORM,
        "benchmark_panel": KELVIN_BENCHMARK_PANEL,
        "kameari_canonical": KELVIN_KAMEARI_CANONICAL,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
