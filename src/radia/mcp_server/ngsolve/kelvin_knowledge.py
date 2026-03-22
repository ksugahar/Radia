"""
Kelvin transformation knowledge base for Radia MCP server.

Covers the Kelvin inversion technique for open boundary FEM problems
using NGSolve, as implemented in examples/KelvinTransformation/.
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

1. **Domain Duplication**: Create two copies of the boundary domain:
   - Interior: physical domain (materials + surrounding air)
   - Exterior (Kelvin): inverted image of the outer space

2. **Permeability Modulation** (exterior domain):
   - mu' = (R / rho')^2 * mu_0  (3D)
   - Accounts for the Jacobian of the Kelvin map

3. **Background Field Modulation** (exterior domain):
   - Hs' = -(rho' / R)^2 * Hs_interior
   - Sign flip + spatial scaling to maintain field continuity

4. **Periodic Boundary Conditions**:
   - Identify interior outer boundary with exterior outer boundary
   - Ensures field continuity at the inversion sphere/circle
   - Uses NGSolve `Periodic()` FE space wrapper

5. **GND Point (Dirichlet at Infinity)**:
   - Place vertex at exterior domain center (r'=0)
   - Corresponds to physical infinity
   - Essential for unique solution: `dirichlet_bbnd="GND"`

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
transformation.

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
| Metric | (R/r')^1 | (R/rho')^2 |
| mu modulation | R^2/r'^2 * mu0 | R^2/rho'^2 * mu0 |
| Hs modulation | -(r'/R)^2 | -(rho'/R)^2 |
| Analytical (dipole) | 2/(mu_r+1) | 3/(mu_r+2) |
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

## 確認済み精度 (NGSolve 6.2.2601, maxh=0.2)

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
# NGSolve 6.2.2601: キーワード引数不可
u, v = fes.TnT()           # OK
u, v = fes.TnT(u_='u', v_='v')  # TypeError: NGSolve 6.2.2601 では不可
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

## 確認済み数値 (NGSolve 6.2.2601, 1/6 セクタ, maxh=0.2)

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

1. **Forgetting permeability modulation in exterior domain**
   - mu' = (R/rho')^2 * mu0, NOT just mu0
   - This accounts for the Jacobian of the Kelvin map

2. **Wrong sign on background field modulation**
   - Exterior Hs has OPPOSITE sign to interior Hs
   - Hs_outer = -(rho'/R)^2 * Hs_inner

3. **Missing GND point**
   - Without Dirichlet at infinity (GND), solution is non-unique
   - Place vertex at center of exterior domain

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


def get_kelvin_documentation(topic: str = "all") -> str:
    """Return Kelvin transformation documentation by topic."""
    topics = {
        "overview": KELVIN_OVERVIEW,
        "h_formulation": KELVIN_H_FORMULATION,
        "a_formulation": KELVIN_A_FORMULATION,
        "3d": KELVIN_3D,
        "adaptive": KELVIN_ADAPTIVE,
        "identify": KELVIN_IDENTIFY,
        "tips": KELVIN_TIPS,
        "verification": KELVIN_VERIFICATION,
        "periodic_wedge": KELVIN_PERIODIC_WEDGE,
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
