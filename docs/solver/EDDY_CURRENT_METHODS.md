# Conductor Eddy Current Modeling: Method Comparison

Date: 2026-02-22

Three approaches for modeling eddy currents in conductors using NGSolve + ngbem.

## Summary Table

| | **FEM-BEM** (ngbem) | **BEM + SIBC** (ngbem) | **FEM + Kelvin** (NGSolve) |
|---|---|---|---|
| **Mesh** | Volume (conductor) | Surface only | Volume (conductor + air) |
| **Exterior domain** | BEM (exact) | BEM (exact) | Kelvin transform (FEM) |
| **Skin effect** | Exact (PDE solved) | SIBC approx (delta << size) | Exact (PDE solved) |
| **mu_r != 1** | Yes (arbitrary) | Yes | Yes (nonlinear B-H) |
| **Frequency sweep** | Slow (re-assembly) | Fast (operators cached) | Slow (re-assembly) |
| **DOF count** | Medium | Low | High (volume + air) |
| **PEEC coupling** | Possible | Direct (Biot-Savart) | Not available |
| **Best use case** | General 3D eddy current | High-freq shields | Nonlinear materials, multi-loop conductors |

## 1. ngbem FEM-BEM Coupling

### Overview

Interior FEM for the conductor volume + exterior BEM using Maxwell SLP. No air mesh needed. The BEM handles the unbounded exterior exactly.

```
+------------------+
|  Conductor       |  <-- FEM: curl(1/mu * curl A) + jw*sigma*A = source
|  (volume mesh)   |
+--------+---------+
         |
    tangential trace (coupling)
         |
+--------v---------+
|  Exterior (air)  |  <-- BEM: Maxwell SLP (no mesh needed)
|  unbounded       |
+------------------+
```

### Available Formulations

| Class | Formulation | Unknowns | BEM Operators |
|---|---|---|---|
| `EddyCurrentFEMBEM` | Scalar Hz (Costabel) | Hz, dHz/dn | V, K, D |
| `VectorEddyCurrentFEMBEM` | Vector A (Johnson-Nedelec + Weggler) | A, j_surf, rho_surf | HelmholtzSL only |

### Usage (Vector A-formulation)

```python
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'path/to/src/radia'))

from ngbem_eddy import VectorEddyCurrentFEMBEM, create_conductor_mesh

# Create conductor volume mesh
mesh = create_conductor_mesh(0.02, 0.02, 0.01, maxh=0.005)  # 20x20x10 mm

# Create solver (arbitrary mu_r)
solver = VectorEddyCurrentFEMBEM(mesh, sigma=3.7e7, mu_r=100.0, order=1)

# Assemble (expensive, done once)
solver.assemble(kappa=0.01, intorder=4)

# Solve at frequency
result = solver.solve(freq=100e3, B_ext=[0, 0, 1.0])
print(f"Condition number: {result['cond']:.2e}")

# Compute loss
P = solver.compute_loss()
print(f"Eddy current loss: {P:.4e} W")
```

### Usage (Scalar Hz formulation)

```python
from ngbem_eddy import EddyCurrentFEMBEM, create_conductor_mesh

mesh = create_conductor_mesh(0.02, 0.02, 0.01, maxh=0.003)

solver = EddyCurrentFEMBEM(mesh, sigma=3.7e7, mu_r=1.0, order=2)
solver.assemble_fembem(freq=100e3, intorder=12)
solver.solve(B_ext=[0, 0, 1.0], mode='fembem')
P = solver.compute_loss()
```

### Advantages

- Exact interior solution (full PDE solved, no skin-depth approximation)
- No air mesh needed (BEM handles exterior)
- Arbitrary mu_r (magnetic conductors: steel, ferrite)
- Vector formulation captures all 3D eddy current patterns
- Weggler stabilization: condition number O(1) for all frequencies

### Limitations

- Volume mesh required for conductor
- FEM block must be re-assembled for each frequency (not fast for sweeps)
- Higher DOF count than surface-only methods
- Mesh must resolve skin depth for accurate loss computation

### When to Use

- General 3D eddy current problems
- Magnetic conductors (mu_r >> 1)
- Thick skin depth (delta ~ conductor dimension)
- When accuracy of interior fields is required

---

## 2. ngbem BEM + SIBC

### Overview

Surface-only BEM with Surface Impedance Boundary Condition. The SIBC replaces the interior FEM by approximating the Dirichlet-to-Neumann map analytically. No volume mesh needed -- only the conductor surface.

```
+------------------+
|  Conductor       |  <-- SIBC: dHz/dn = gamma * Hz
|  (surface only)  |      gamma = sqrt(jw*mu*sigma) = (1+j)/delta
+--------+---------+
         |
    surface operators (BEM)
         |
+--------v---------+
|  Exterior (air)  |  <-- BEM: V, K, D operators (cached)
|  unbounded       |
+------------------+
```

### Available Formulations

| Class | Formulation | Unknowns | BEM Operators |
|---|---|---|---|
| `EddyCurrentBEMSIBC` | Scalar Hz + SIBC | Hz, dHz/dn on surface | V, K, D (cached) |
| `ShieldBEMSIBC` | Loop-only (div-free currents) | Loop currents I_loop | Maxwell SLP only |

### Usage (ShieldBEMSIBC -- recommended for shielding + PEEC)

```python
from ngbem_eddy import ShieldBEMSIBC, create_conductor_mesh, _biot_savart_A

# Create conductor mesh (volume mesh, but BEM uses surface only)
mesh = create_conductor_mesh(0.1, 0.1, 0.002, maxh=0.025)  # 100x100x2 mm plate

# Create solver
solver = ShieldBEMSIBC(mesh, sigma=3.7e7, mu_r=1.0)

# Assemble BEM operators (expensive, done ONCE)
solver.assemble(intorder=4)

# Define incident field from PEEC coil
import numpy as np
wire_center = np.array([0, 0, 0.011])  # 10mm above plate
wire_dir = np.array([1, 0, 0])

def A_inc(points):
    return 1.0 * _biot_savart_A(points, wire_center, wire_dir, 0.2)

# Solve (fast -- operators already cached)
solver.solve(freq=100e3, A_inc_func=A_inc)
P = solver.compute_loss()
print(f"Loss: {P:.4e} W")

# Fast frequency sweep
import numpy as np
for f in [1e3, 10e3, 100e3, 1e6]:
    solver.solve(f, A_inc_func=A_inc)
    print(f"f={f/1e3:.0f} kHz: P={solver.compute_loss():.4e} W")
```

### Usage (EddyCurrentBEMSIBC -- scalar Hz)

```python
from ngbem_eddy import EddyCurrentBEMSIBC, create_conductor_mesh

mesh = create_conductor_mesh(0.01, 0.01, 0.005, maxh=0.003)

solver = EddyCurrentBEMSIBC(mesh, sigma=5.8e7, mu_r=1.0, order=2)
solver.assemble(intorder=12)  # Done once

# Fast frequency sweep (BEM operators cached)
results = solver.frequency_sweep(
    freqs=np.logspace(3, 7, 50),
    B_ext=[0, 0, 1.0]
)
```

### Advantages

- **Fast frequency sweeps** (BEM operators assembled once, only gamma changes)
- **Low DOF count** (surface only, loop basis even fewer)
- **Direct PEEC coupling** (ShieldBEMSIBC accepts Biot-Savart A_inc_func)
- No volume mesh computation per frequency
- ShieldBEMSIBC: mathematically clean (divergence-free loop basis, V2 term vanishes)

### Limitations

- **SIBC valid only when delta << conductor size** (typically f > 1 kHz for metals)
- Less accurate at low frequency (DC breakdown: gamma -> 0)
- Cannot resolve interior field distribution (surface-only)
- ShieldBEMSIBC requires topologically closed surface (genus-0)

### When to Use

- Electromagnetic shielding analysis
- Frequency sweeps (1 kHz - 100 MHz)
- Conductor skin effect (delta << conductor dimensions)
- PEEC coil + shield coupling
- Fast impedance extraction

---

## 3. NGSolve FEM + Kelvin Transform (EMPY)

### Overview

Pure FEM approach using the Kelvin transformation to handle the unbounded exterior domain. The Kelvin transform maps the infinite exterior to a bounded region, eliminating the need for BEM or artificial boundary conditions (PML/ABC). **Supports both magnetostatic and eddy current analysis.**

Two formulations are implemented in EMPY (`S:\NGSolve\EMPY\EMPY_Analysis\EddyCurrent`):

| Formulation | Unknowns | FE Spaces | Multi-Loop |
|---|---|---|---|
| **A-Phi** | Vector potential A, scalar potential Phi | HCurl x H1 | No |
| **T-Omega** | Current potential T, magnetic potential Omega | HCurl x H1 + loop fields | **Yes** |

```
+------------------+     +------------------+
|  Conductor       |     |  Air region      |
|  (volume mesh)   |     |  (reduced pot.)  |
|  sigma > 0       |     |  sigma = 0       |
+--------+---------+     +--------+---------+
         |                         |
    conductor boundary             |
         |                         |
+--------v---------+     +--------v---------+
|  Outer air       |     |  Kelvin region   |
|  (FEM mesh)      | <-> |  r -> R^2/r      |
+------------------+     +------------------+
                          (bounded substitute
                           for infinity)
```

### Usage (A-Phi formulation -- simpler)

```python
from netgen.occ import Box, Pnt, Sphere, OCCGeometry, Glue
from ngsolve import *
import math

# Create geometry: conductor + air
conductor = Box(Pnt(-0.055, -0.03, -0.003175), Pnt(0.055, 0.03, 0.003175))
conductor.mat("conductor")
conductor.faces.name = "conductorBND"

air_sphere = Sphere(Pnt(0,0,0), 0.2)
air = air_sphere - conductor
air.mat("air")
air_sphere.faces.name = "outer_boundary"

geo = OCCGeometry(Glue([conductor, air]))
mesh = Mesh(geo.GenerateMesh(maxh=0.01))

# Material properties
sigma = 3.278e7  # Al conductivity [S/m]
mu0 = 4e-7 * math.pi
f = 50  # Hz
s = 2 * math.pi * f  # MQS: real frequency

Sigma = CoefficientFunction([{"conductor": sigma, "air": 0}[m] for m in mesh.GetMaterials()])
Mu = CoefficientFunction([mu0 for m in mesh.GetMaterials()])

# A-Phi FE spaces
fesA = HCurl(mesh, order=1, nograds=True, complex=False)
fesPhi = H1(mesh, order=1, definedon="conductor", complex=False)
fesAPhi = fesA * fesPhi
(A, phi), (N, psi) = fesAPhi.TnT()

# Bilinear form: s*sigma*(A+grad(phi))*(N+grad(psi)) + 1/mu*curl(A)*curl(N)
a = BilinearForm(fesAPhi)
a += s * Sigma * (A + grad(phi)) * (N + grad(psi)) * dx("conductor")
a += (1.0 / Mu) * curl(A) * curl(N) * dx
a.Assemble()

# Source: uniform Bz=1T via Neumann BC
normal = specialcf.normal(mesh.dim)
h_source = (0, 0, 1.0 / mu0)  # H = B/mu0
sr = LinearForm(fesAPhi)
sr += (1.0 / mu0) * Cross(N.Trace(), CoefficientFunction(h_source)) * normal * ds("outer_boundary")
sr.Assemble()

# Solve
gfAPhi = GridFunction(fesAPhi)
gfAPhi.vec.data = a.mat.Inverse(fesAPhi.FreeDofs()) * sr.vec
gfA, gfPhi = gfAPhi.components

# Extract fields
Bfield = curl(gfA)
Jfield = -Sigma * s * (gfA + grad(gfPhi))
```

### Usage (T-Omega formulation -- multi-loop conductors)

```python
# T-Omega method for conductor with holes (multi-loop)
# Based on EMPY T_Omega_Method class
#
# H = T + grad(Omega) + sum(amp_k * loopField_k)
# T: current potential in HCurl(conductor)
# Omega: magnetic scalar potential in H1(all)
# loopField_k: topological loop basis from genus computation

fesT = HCurl(mesh, order=1, nograds=True, definedon="conductor",
             dirichlet="conductorBND", complex=False)
fesOmega = H1(mesh, order=1, dirichlet="outer_boundary", complex=False)
fesTOmega = fesT * fesOmega
(T, omega), (W, psi) = fesTOmega.TnT()

# System matrix
a = BilinearForm(fesTOmega)
a += (1.0 / Sigma) * curl(T) * curl(W) * dx("conductor")
a += s * Mu * T * (W + grad(psi)) * dx("conductor")
a += s * Mu * grad(omega) * grad(psi) * dx
a.Assemble()

# Set Omega boundary value for external Bz
gfTOmega = GridFunction(fesTOmega)
gfT, gfOmega = gfTOmega.components
Omega0 = Bz0 * box_half_height * 2 / mu0
gfOmega.Set(Omega0, definedon=mesh.Boundaries("outer_boundary"))

# Solve (simplified: single-loop case without loop field coupling)
source = -a.mat * gfTOmega.vec
gfTOmega.vec.data += a.mat.Inverse(fesTOmega.FreeDofs()) * source

# Direct R, L extraction
R = Integrate((1.0 / Sigma) * curl(gfT)**2 * dx("conductor"), mesh)
L = Integrate(Mu * (gfT + grad(gfOmega))**2 * dx, mesh)
Z = R + s * L
```

### Advantages

- **No BEM needed** (pure FEM, all standard NGSolve tools available)
- **Nonlinear mu_r(H)** supported (Newton iteration)
- **Exact open boundary** (Kelvin transform, no PML/ABC approximation)
- **Multi-loop conductors** (T-Omega with topological loop fields)
- **Direct R, L extraction** (T-Omega CalcRL formula)
- Standard NGSolve solvers (BDDC, multigrid, ICCG, etc.)
- Natural handling of multi-material interfaces
- **MQS regime**: valid for power frequency to ~MHz

### Limitations

- **Air mesh required** (3D volume mesh extends to Kelvin/outer radius)
- **High DOF count** (conductor + air volume all meshed)
- Mesh must resolve skin depth for accurate loss computation
- No circuit extraction for SPICE (use Radia PEEC instead)
- No PEEC coupling
- Each frequency requires full re-assembly and re-solve

### When to Use

- **Nonlinear magnetic materials** (saturation effects, B-H curves)
- When BEM is not available (no ngbem dependency)
- **Multi-loop conductors** (plates with holes, T-Omega method)
- Problems where interior field distribution in air is needed
- Multi-physics coupling (thermal, mechanical) where volume data is required
- Magnetostatic + eddy current combined analysis
- **Direct impedance extraction** (R + jwL)

---

## Decision Guide

### Which Eddy Current Method to Use

```
                         Start
                           |
               Need eddy currents?
                      /       \
                    Yes        No -> Kelvin Transform (Omega-Reduced Omega)
                     |
            Nonlinear mu_r(H)?
                 /          \
               Yes           No
                |             |
       FEM-SIBC + ESIM       |
       (Karl iteration,      |
        calc_fem_kelvin.py)   |
                              |
                    Need air mesh free? (BEM vs FEM)
                         /         \
                       Yes          No
                        |            |
                   BEM methods   FEM-SIBC (interface)
                        |        (calc_fem_kelvin.py)
                        |
              +-------------------+
              |                   |
        Workpiece (SIBC)    Coil self-L
              |                   |
     Scalar BIE + SIBC    EFIE saddle point
     (all Z_s, <0.1%)     (LaplaceSL, PEC)
     bem_sibc_solver.py    bem_inductance.py
```

### Method Selection Table

| Scenario | Recommended Method | Why |
|----------|-------------------|-----|
| Coil self-inductance (L_air) | EFIE saddle point | PEC assumption OK for copper coils |
| Workpiece eddy current (linear) | **Scalar BIE + SIBC** | All Z_s, surface-only, fast freq sweep |
| Workpiece P + L + B distribution | **FEM-SIBC (interface)** | P (-2.3%), L, B in one solve |
| Workpiece eddy current (nonlinear mu) | FEM-SIBC + ESIM (interface) | H-dependent Z_s via Karl iteration |
| Coil + magnetic core | FEM + Kelvin | Nonlinear iron needs volume mesh |
| High-frequency (>1 MHz) | ngsolve.bem Helmholtz | Wave effects matter |
| Circuit extraction (L, R, C, M) | Radia PEEC pipeline | Direct SPICE netlist output |

### FEM-SIBC: Hole Approach (Correct)

**SIBC = Robin BC** on the conductor surface.  The conductor interior is
NOT solved -- it is replaced by the impedance condition.  Subtract the
workpiece from the mesh (hole) and apply Robin BC on the boundary:

```
dphi/dn + (jw*mu_0/Z_s) * phi = 0      (H1/phi formulation)
n x E   + Z_s * (n x H x n)  = 0       (HCurl formulation)
```

**Interface approach (wrong for H1)**: workpiece meshed as air + Robin
on internal interface.  Flux bypasses the Robin BC through the
transparent interior.  Verified in 2D axisym: interface gives 58-94% of
the correct L (steel worst), while hole gives > 99.6%.

**Note on 3D HCurl**: the earlier interface result (-2.9% H_t) was for a
specific 3D HCurl test case where the curl-curl operator constrains
the interior differently.  The general principle remains: SIBC replaces
the interior, so don't solve inside it.

**H_t extraction**: BND integral with tangential projection using `specialcf.normal(3)`:

```python
n = specialcf.normal(3)
A_sq = sum(gfu[i].real**2 + gfu[i].imag**2 for i in range(3))
Adn_re = sum(gfu[i].real * n[i] for i in range(3))
Adn_im = sum(gfu[i].imag * n[i] for i in range(3))
At_sq = A_sq - Adn_re**2 - Adn_im**2     # |A_t|^2
H_t_rms = abs(1j*omega/Z_s) * sqrt(Integrate(At_sq, mesh, BND, ...) / A_wp)
P = 0.5 * Re(Z_s) * H_t_rms**2 * A_wp    # time-averaged
```

**Without tangential projection**, `|A|^2` includes the normal component,
which is inflated by the Robin penalty. The previous implementation
(before 2026-04-09) lacked this projection and returned H_t ~ 0.

**Robin sign**: `+jw/Z_s` (positive). With negative sign or without Robin,
the solution diverges or gives PEC.

**GMSH |B| export**: Per-vertex `|B|` field written to `.msh v4.1` via `GmshPostExport`.

```python
# CRITICAL: unpack tuple for GridFunction evaluation
node_B = np.array([gf_B(mesh(*mesh.vertices[v.nr].point))
                   for v in mesh.vertices])
# gf(mesh.vertices[v].point) FAILS - tuple is not a MeshPoint
# gf(*mesh.vertices[v].point) also works (unpacks to x, y, z)
```

**Panel script**: `src/radia/panels/calc_fem_kelvin.py`

### 2D Axisymmetric SIBC Validation (2026-04-14)

Independent validation using 2D axisymmetric FEM with z-offset Kelvin
open boundary.  Full-resolution (eddy currents resolved at delta/5)
compared with SIBC (hole + Robin BC).

Script: `validation_test/eddy_current_analytical_validation/reference_2d_axisym.py`

**Formulation**: phi = r * A_phi (H1 space).
SIBC = Robin BC on hole boundary: `(jw/Z_s) / r * u * v * ds("wp_bnd")`.
Conductor interior NOT meshed, NOT solved.

**Surface impedance**:
- **Cylindrical Bessel** (exact for solid cylinder):
  `Z_s = rho * gamma * I1(gamma*a) / I0(gamma*a)`,
  `gamma = sqrt(jw * mu_r * mu_0 * sigma)`.
- Dowell flat-slab: acceptable for L, but P off 16-23% at low R/delta.

**Results (Kelvin + cylindrical Bessel Z_s)**:

| Material | R/delta | L error | P error |
|----------|---------|---------|---------|
| Copper | 3.8 - 84.6 | < 0.3% | < 1% |
| Steel (mu_r=100) | 7.0 - 157 | < 0.4% | < 2% |
| Aluminum | 2.9 - 65.7 | < 0.8% | < 2% |

Air-only L = 99.47 nH (Neumann = 99.23, 0.2% match).

**Conclusion**: SIBC Robin BC on hole is validated.  The same approach
(hole + Robin + cylindrical Bessel Z_s) applies to 3D.

| Material | R/delta | L error (cyl) | P error (cyl) |
|----------|---------|-------------|-------------|
| Copper | 3.8 - 84.6 | < 0.1% | < 1.0% |
| Steel (mu_r=100) | 7.0 - 157 | < 0.4% | < 3.3% |
| Aluminum | 2.9 - 65.7 | < 0.6% | < 2.8% |

**SIBC = Robin BC on the conductor surface; do NOT solve inside.**
The conductor interior is replaced by the surface impedance condition.
Use the **hole approach** (conductor subtracted from the mesh) and apply
the Robin BC `dphi/dn + (jw*mu_0/Z_s)*phi = 0` on the hole boundary.

The interface approach (conductor meshed as air + Robin on internal
interface) is **wrong** because flux bypasses the Robin BC through the
transparent interior. Verified in 2D axisymmetric H1/phi: interface
gives 58-94% of the correct L (steel worst), while hole gives > 99.6%.

### Why Scalar BIE + SIBC (Not EFIE or MFIE)

| Formulation | Unknown | Z_s Support | Error | Issue |
|-------------|---------|-------------|-------|-------|
| EFIE-SIBC | J (HDivSurface) | PEC only | 65% at Z_s/beta=10 | SL eigenvalue R/3 |
| MFIE | J (HDivSurface) | PEC only | 0% (PEC) | No Z_s dependence (M->H=0 in MQS) |
| PMCHWT | J + M | Impossible in MQS | N/A | jw*eps*SL(M) ~ 10^{-14} |
| **Scalar BIE** | **phi (H1)** | **All Z_s** | **<0.1%** | **None** |

The scalar potential BIE avoids all three failure modes because:
1. No SL*J term (avoids R/3 eigenvalue)
2. Uses H1 stiffness K (nonzero surface Laplacian, unlike RT0 which has zero curl_s)
3. Single-equation formulation (no need for M coupling)

## Performance Characteristics

### Assembly Time (20x20x10 mm Al cube, coarse mesh)

| Method | Assembly | Per-Frequency Solve | Notes |
|---|---|---|---|
| VectorEddyCurrentFEMBEM | ~6 s | 0.02 s | BEM assembly is the bottleneck |
| EddyCurrentFEMBEM | ~3 s | 0.01 s | Fewer BEM DOFs (scalar) |
| ShieldBEMSIBC | ~60 s | 0.01 s | Larger surface mesh (plate geometry) |
| EddyCurrentBEMSIBC | ~3 s | 0.001 s | Fastest per-freq (gamma only) |
| FEM + Kelvin (A-Phi) | ~1 s | ~1 s | Large system (volume + air mesh) |
| FEM + Kelvin (T-Omega) | ~2 s | ~1 s | + loop field computation |

### DOF Count Comparison (same geometry)

| Method | Typical DOFs | Breakdown |
|---|---|---|
| **Scalar BIE + SIBC** | **~128** | **H1 vertices on surface** |
| VectorEddyCurrentFEMBEM | ~164 | 94 H(curl) + 42 HDivSurf + 28 SurfL2 |
| EddyCurrentFEMBEM | ~100 | H1 interior + H1/2 boundary |
| ShieldBEMSIBC | ~174 loops | From V-1 (genus-0 surface) |
| FEM + Kelvin (A-Phi) | ~5000+ | HCurl(all) + H1(conductor) |
| FEM + Kelvin (T-Omega) | ~3000+ | HCurl(conductor) + H1(all) + loops |

## Required ngbem/ngsolve.bem Operators

| Method | Operators Used | Package |
|---|---|---|
| **Scalar BIE + SIBC** | **LaplaceSL, LaplaceDL (scalar H1)** | **ngsolve.bem** |
| EddyCurrentFEMBEM | V (SLP), K (DLP), D (Hypersingular) | ngbem or ngsolve.bem |
| VectorEddyCurrentFEMBEM | HelmholtzSL only | ngsolve.bem |
| EddyCurrentBEMSIBC | V, K, D (cached) | ngbem or ngsolve.bem |
| ShieldBEMSIBC | MaxwellSingleLayerPotentialOperator | ngbem |
| FEM + Kelvin | (none -- pure FEM) | NGSolve only |

## Cross-Method Validation

### Test Setup

- Geometry: 20x20x10 mm aluminum block (sigma = 3.7e7 S/m, mu_r = 1)
- Excitation: Uniform Bz = 1 T
- Mesh: maxh = 5 mm (volume mesh for FEM-BEM, surface from same mesh for BEM+SIBC)
- Analytical reference: half-space model (P_sides = 0.5 * Rs * H0^2 * side_area)

### Results

| f (kHz) | delta (mm) | VectorFEMBEM | ShieldBEMSIBC | ScalarFEM | Analytical |
|---------|-----------|-------------|--------------|----------|-----------|
| 1       | 2.62      | 2443        | 292          | 1340     | 2617      |
| 3       | 1.51      | 4588        | 1582         | 2918     | 4532      |
| 10      | 0.83      | 7036        | 6742         | 4808     | 8274      |
| 30      | 0.48      | 7989        | 13371        | 5565     | 14331     |
| 100     | 0.26      | 8157        | 19945        | 5714     | 26165     |

### Method Validity Ranges

Each method has two validity parameters:
- **delta/D_min** (skin depth / minimum conductor dimension): controls SIBC validity
- **maxh/delta** (mesh size / skin depth): controls FEM mesh resolution

**VectorEddyCurrentFEMBEM** (FEM interior + BEM exterior):
- Valid when mesh resolves skin depth: **maxh/delta < 3**
- At maxh=5mm: accurate up to ~10 kHz (delta=0.83mm, maxh/delta=6)
- Saturates at high freq because volume mesh cannot resolve thin skin layer
- No SIBC approximation: exact solution within mesh resolution limits

**ShieldBEMSIBC** (surface BEM + SIBC):
- Valid when skin depth << ALL conductor dimensions: **delta/D_min < 0.1**
- For 10mm-thick block: accurate above ~10 kHz (delta/D_min=0.08)
- **CRITICAL: Massively overestimates for thin conductors** (see below)
- Fast frequency sweeps (BEM operators assembled once)

**ScalarFEM (Hz only, FEM-only mode)**:
- Captures only the Hz component of magnetic field
- Loss ratio = ScalarFEM/Vector correlates with side face fraction
- Useful as quick estimate but NOT quantitatively accurate for 3D

**ScalarSIBC (Hz + SIBC)**:
- Compound error: scalar Hz limitation + SIBC limitation
- Loss DECREASES with frequency (physically wrong for all 3D geometries)
- **NOT suitable for any 3D eddy current analysis**

### Key Finding: Geometry-Dependent Agreement

The agreement between VectorFEMBEM and ShieldBEMSIBC depends strongly on geometry:

**Block 20x20x10mm** (D_min=10mm):

| Frequency | delta/D_min | maxh/delta | Shield/Vector | Verdict |
|-----------|-----------|----------|-------------|---------|
| 1 kHz     | 0.262     | 1.9      | 0.12        | SIBC invalid |
| 3 kHz     | 0.151     | 3.3      | 0.34        | SIBC marginal |
| **10 kHz** | **0.083** | **6.0**  | **0.96**    | **4.2% agreement** |
| 30 kHz    | 0.048     | 10.5     | 1.67        | FEM mesh too coarse |
| 100 kHz   | 0.026     | 19.1     | 2.45        | FEM mesh too coarse |

**Thin 20x20x5mm** (D_min=5mm):

| Frequency | delta/D_min | maxh/delta | Shield/Vector | Verdict |
|-----------|-----------|----------|-------------|---------|
| 1 kHz     | 0.523     | 1.5      | 0.46        | SIBC invalid |
| 3 kHz     | 0.302     | 2.6      | 1.35        | Both marginal |
| 10 kHz    | 0.165     | 4.8      | 3.52        | SIBC overestimates |
| 30 kHz    | 0.096     | 8.4      | 7.53        | SIBC massively wrong |
| 100 kHz   | 0.052     | 15.3     | 12.2        | SIBC 12x overestimate |

### ShieldBEMSIBC Thin-Conductor Failure

**CRITICAL**: ShieldBEMSIBC overestimates loss by 3-12x for the thin (5mm) block, even when delta/D_min < 0.1. This is because SIBC assumes each surface is a semi-infinite half-space. For thin conductors, the top and bottom surfaces are close enough that their fields interact -- the SIBC approximation counts this interaction twice.

At 100 kHz for the thin block:
- ShieldSIBC: 69,480 W
- Analytical (ALL 6 faces, half-space): 39,247 W
- ShieldSIBC exceeds even the all-faces upper bound by 1.77x

**ShieldBEMSIBC is valid ONLY when delta << ALL dimensions** (not just D_min). For thin plates, use VectorEddyCurrentFEMBEM with a fine mesh, or FEM+Kelvin.

### Scalar Hz Limitations: Quantitative

The scalar Hz ratio correlates with **side face fraction** (area where Hz is tangential for Bz excitation):

| Geometry | Side fraction | ScalarFEM/Vector | Explanation |
|----------|--------------|-----------------|-------------|
| Thin 20x20x5mm | 0.33 | 0.26 - 0.79 | Hz normal on 67% of area |
| Block 20x20x10mm | 0.50 | 0.55 - 0.70 | Hz normal on 50% of area |

Physics: For Bz excitation, Hz is tangential ONLY on side faces (x,y normals). On top/bottom (z-normal), Hz is the normal component and Hx,Hy are tangential. Scalar Hz misses all Hx,Hy-driven eddy currents.

ScalarSIBC: P < 50 W for all cases (correct value: thousands of W). Hz_total/Hz_inc decreases with frequency (0.58 -> 0.12), causing loss to decrease. Fundamentally broken for 3D.

### Method Selection Guide

```
               delta/D_min > 0.3    0.1 < delta/D_min < 0.3    delta/D_min < 0.1
              (thick skin)          (transitional)              (thin skin)

Thick body   VectorFEMBEM          VectorFEMBEM or             ShieldBEMSIBC
(D_min~W~H)  (if maxh < delta)     ShieldBEMSIBC               (fast sweeps)

Thin body    VectorFEMBEM          VectorFEMBEM ONLY           VectorFEMBEM ONLY
(D_min<<W)   (fine mesh required)  (SIBC overestimates)        (SIBC overestimates)
```

**Rule of thumb**: ShieldBEMSIBC is safe when the conductor is roughly equi-dimensional (cube-like) AND delta/D_min < 0.1. For thin plates/shells, always use VectorFEMBEM or FEM+Kelvin.

---

## VectorFEMBEM Cross-Validation Analysis

Detailed analysis of VectorEddyCurrentFEMBEM accuracy and bugs discovered during cross-validation against ShieldBEMSIBC on the same 20x20x10 mm aluminum block test case (see Cross-Method Validation above for setup).

### Observation: Frequency-Independent Loss

Cross-validation revealed that VectorFEMBEM loss is **frequency-independent** (~4.3 kW at all frequencies on a coarse mesh). Root cause analysis identified **two independent bugs** that interact paradoxically.

Coarse mesh parameters used for diagnosis: maxh = 12 mm (18 volume elements, 94 H(curl) DOFs), skin depth range 0.83 mm (10 kHz) to 0.08 mm (1 MHz).

### Bug 1: Missing curl-curl RHS Term

**File**: `src/radia/ngbem_eddy.py`, `VectorEddyCurrentFEMBEM.solve()`, line 2138

The interior FEM equation for the scattered field A_s is:

```
curl(1/mu * curl A_s) + j*w*sigma*A_s = -curl(1/mu * curl A_inc) - j*w*sigma*A_inc
```

In weak form:

```
(1/mu * curl A_s, curl v) + jws*(A_s, v) = -(1/mu * curl A_inc, curl v) - jws*(A_inc, v)
```

So the RHS should be:

```
f_1 = -(a_curl + jws*a_mass) @ A_inc = -a_FEM @ A_inc
```

**Buggy code**:
```python
rhs[:n1] = -1j * omega * sigma * self._a_mass_np @ A_inc_coeffs
```

**Missing**: `-self._a_curl_np @ A_inc_coeffs`

#### Magnitude of the Missing Term

| Frequency | |w*s*M*A_inc| | |a_curl*A_inc| | curl/mass ratio |
|-----------|---------------|----------------|-----------------|
| 10 kHz    | 6.17e+05      | 1.40e+04       | 2.28%           |
| 100 kHz   | 6.17e+06      | 1.40e+04       | 0.23%           |
| 1 MHz     | 6.17e+07      | 1.40e+04       | 0.023%          |

The curl-curl term is small relative to the mass term (2.3% at 10 kHz, negligible at higher frequencies). This is expected for high-conductivity materials where `omega*sigma*mu >> 1/mu`.

### Bug 2 (Fundamental): Coarse Mesh Cannot Resolve Skin Layer

Even with the correct RHS, the loss is **exactly zero** (P ~ 1e-26 W):

| Frequency | P_original (buggy RHS) | P_fixed (correct RHS) |
|-----------|------------------------|----------------------|
| 1 kHz     | 2.525e+03 W            | 1.169e-26 W          |
| 10 kHz    | 4.315e+03 W            | 1.040e-24 W          |
| 100 kHz   | 4.366e+03 W            | 6.581e-23 W          |
| 1 MHz     | 4.367e+03 W            | 1.064e-20 W          |

#### Why P_fixed = 0

With the correct RHS `f_1 = -a_FEM @ A_inc`, the FEM equation becomes:

```
a_FEM * (A_s + A_inc) + B^T * j = 0
a_FEM * A_total + B^T * j = 0
```

On a coarse mesh (maxh = 12 mm >> delta = 0.08-0.83 mm), the H(curl) basis cannot represent the exponential boundary layer `exp(-z/delta)`. The "best" FEM solution is:

```
A_total = 0  (everywhere inside)
j = 0        (surface current)
```

This is physically the **perfect conductor limit** (sigma -> infinity). The coarse mesh FEM cannot distinguish between sigma = 3.7e7 and sigma = infinity.

#### Shielding Diagnostic Confirms PEC Limit

| Frequency | |A_total|/|A_inc| | cos(A_scat, A_inc) |
|-----------|-------------------|--------------------|
| 1 kHz     | 3.05e-01          | -0.952             |
| 10 kHz    | 3.92e-02          | -0.999             |
| 100 kHz   | 3.94e-03          | -1.000             |
| 1 MHz     | 3.94e-04          | -1.000             |

A_scat -> -A_inc (perfect cancellation), A_total -> 0.

#### Why P_original ~ const (the Paradox)

Without the curl-curl RHS, the FEM equation is:

```
(a_curl + jws*M) * A_s + B^T * j = -jws * M * A_inc
```

At high omega, `jws*M` dominates, so `A_s ~ -A_inc`. But the curl-curl part `a_curl * A_s ~ -a_curl * A_inc` creates an unbalanced residual:

```
A_total ~ a_curl * A_inc / (jws * M)
|A_total|^2 ~ |a_curl*A_inc|^2 / (w^2 * s^2 * ||M||^2)
```

Loss becomes:
```
P = 0.5 * w^2 * s * |A_total|^2 ~ 0.5 * |a_curl*A_inc|^2 / (s * ||M||^2)
```

This is **frequency-independent** (the w^2 in the loss formula exactly cancels the 1/w^2 in |A_total|^2). The "nonzero loss" is entirely a numerical artifact from the incomplete RHS.

### Root Cause Interaction Summary

| Code state       | A_total            | Loss              | Correct? |
|------------------|--------------------|-------------------|----------|
| Original (buggy) | Small, ~1/omega    | ~4.3 kW (const)   | NO (artifact from RHS imbalance) |
| Fixed (curl-curl)| Zero (to precision)| ~0 W              | NO (coarse mesh = perfect conductor) |
| Physical truth   | exp(-z/delta) layer| ~sqrt(omega)       | YES (needs delta-scale mesh) |

**Neither code state gives correct loss on a coarse mesh.** The volume FEM loss formula `P = 0.5*w^2*s*|A_total|^2` fundamentally requires mesh elements finer than the skin depth to resolve the boundary layer where the physical loss occurs.

### Fixes Applied

#### Fix 1: Total Field Formulation

Switched from scattered field (A_scat as unknown) to total field (A_total as unknown) to avoid catastrophic cancellation:

**Old (scattered)**: A_scat + A_inc ~ 0 (both large, nearly cancel)
**New (total)**: A_total is direct unknown (no cancellation)

```python
# Old: Row 1 RHS = -jws*M*A_inc (buggy: missing curl-curl)
# New: Row 1 RHS = 0 (homogeneous interior)
# New: Row 2 RHS = +B @ A_inc (BEM drives system)
```

Result: A_total -> 0 (PEC limit on coarse mesh), j is frequency-independent and nonzero (correct for PEC surface current).

#### Fix 2: Surface Current j Not Physical

Investigation showed that j in the Johnson-Nedelec FEM-BEM coupling is an auxiliary SLP representation density, NOT the physical surface current K = n x H. The normalization differs by ~10^10 from the physical current.

**Why**: In VectorFEMBEM, the BEM equation relates j to the boundary trace of the vector potential A (not to H or E). ShieldBEMSIBC uses the EFIE which directly solves for the physical surface current K.

#### Fix 3: Analytical SIBC Loss Estimate

Added `compute_loss_sibc()` method that computes loss from the known incident field H_inc on each boundary triangle:

```python
P = sum_faces 0.5 * Re(Zs) * |H_inc_tangential|^2 * area
```

This is the half-space SIBC approximation applied face-by-face.

#### Cross-Validation: compute_loss_sibc() vs ShieldBEMSIBC

| Frequency | P_analSIBC (W) | P_shield (W) | Ratio |
|-----------|----------------|--------------|-------|
| 1 kHz     | 2,617          | 678          | 3.87  |
| 10 kHz    | 8,274          | 9,005        | 0.92  |
| 100 kHz   | 26,165         | 19,604       | 1.33  |
| 1 MHz     | 82,741         | 50,796       | 1.63  |

The analytical SIBC gives the correct frequency scaling (sqrt(f)) and is within a factor of 0.9-3.9x of ShieldBEMSIBC. The deviation comes from:
- Low freq (delta ~ thickness): half-space assumption breaks down
- High freq (delta << thickness): no edge/corner enhancement in flat approximation

### VectorFEMBEM Solver Selection by Regime

| Regime | delta vs thickness | Recommended Solver | Notes |
|--------|--------------------|--------------------|-------|
| Thick skin | delta > thickness | VectorFEMBEM (fine mesh) | Need mesh resolution < delta |
| Moderate | delta ~ thickness | VectorFEMBEM (fine mesh) | Boundary layer meshing |
| Thin skin | delta << thickness | **ShieldBEMSIBC** | Mesh-independent, most accurate |
| Quick estimate | Any | VectorFEMBEM.compute_loss_sibc() | 1-4x of ShieldBEMSIBC |

**VectorFEMBEM is most useful** for:
1. Thick-skin problems where skin depth is resolvable by the mesh
2. Magnetic materials (mu_r >> 1) where ShieldBEMSIBC is less accurate
3. Computing A_total and H fields inside the conductor (not just loss)

**ShieldBEMSIBC is preferred** for:
1. Thin-skin shielding analysis (delta << thickness)
2. Loss computation on coarse meshes
3. Quick frequency sweeps (surface-only, no volume DOFs)

### Diagnostic Archive

The old one-off `examples/ngbem_diagnostics/*.py` scripts were promoted and
pruned on 2026-06-27. Use the result-bearing notebook
`docs/ngbem_diagnostics/ngbem_diagnostics_archive.ipynb` and its synchronized
JSON `docs/ngbem_diagnostics/ngbem_diagnostics_archive_results.json` for the
historical source, SHA-256 hashes, and diagnostic notes. New executable checks
should be added under `validation_test/` with explicit assertions.

---

## 4. ESIM Surface Impedance Cross-Check: BEM vs FEM vs FEM-full

Date: 2026-03-27

### Overview

Three independent methods to compute eddy current loss in a workpiece (cylinder)
excited by a coil (torus), compared for the same geometry and material.

```
Method A: BEM-ESIM    -- BEM(coil J) -> Biot-Savart(H) -> ESIM(Z_s) -> P
Method B: FEM-ESIM    -- FEM+Kelvin(static) -> grad(u)(H) -> ESIM(Z_s) -> P
Method C: FEM-full    -- FEM+Kelvin(AC, sigma in workpiece) -> direct P
```

### Geometry

- Coil: torus R=30 mm, a=3 mm (1 turn, I=1 A)
- Workpiece: cylinder r=10 mm, h=20 mm, at center
- Material: copper (sigma=5.8e7 S/m, mu_r=1)
- Frequency: 1 kHz (skin depth delta=2.09 mm, xi=R/delta=4.8)

### Method D: FEM-SIBC (Karl Hollaus Iteration)

The correct SIBC approach solves FEM with Robin BC and iterates Z_s:

```
Method D: FEM-SIBC   -- FEM+Kelvin(complex, Robin BC) + Karl iteration -> P
```

Robin BC weak form (2D axi, u = r*A_theta):
```
int (nu/r) grad(u).grad(v) dx - (jw/Zs) int u*v/r ds = int J*v dx
```

Power (Poynting flux through surface):
```
P = pi * w^2 * Re(Zs)/|Zs|^2 * int |u|^2/r ds
```
Note: factor is **pi** (not 2*pi). Derivation: `P = (1/2) * 2*pi * int Re(E*H*) r ds`.

Karl iteration:
1. Solve FEM with current Z_s (Robin BC)
2. Sample H_t near workpiece boundary from grad(u)
3. Update Z_s from ESIM cell problem at H_t
4. Under-relaxation: Z_s = 0.5*Z_s_old + 0.5*Z_s_new
5. Repeat until |dZ_s/Z_s| < tol

### Results (I = 1 A)

**Copper (sigma=5.8e7, mu_r=1), 1 kHz, xi=4.8:**

| Method | P [W] | L [nH] | DOFs | Time [s] |
|--------|-------|--------|------|----------|
| **FEM-full** (reference) | 2.30e-6 | 92.0 | 2,878,208 | 103 |
| **FEM-SIBC** (Karl) | 3.24e-6 | 91.4 | 82,835 | 5 |
| **BEM-ESIM** (transparent) | 1.50e-6 | 86.7 | 5,064 | 24 |

**Steel (sigma=2e6, nonlinear BH), 1 kHz (Karl converges in 10 iterations):**

| Method | P [W] | Notes |
|--------|-------|-------|
| **FEM-SIBC** (Karl) | 1.53e-6 | Converged, 10 iter |
| **BEM-ESIM** (transparent) | 1.51e-6 | Good agreement |
| FEM-full | N/A | delta=0.16mm, mesh impractical |

### Analysis

1. **FEM-SIBC vs FEM-full (copper, xi=4.8)**: +40% error. Source: 0th-order SIBC (slab
   approximation for cylindrical geometry). Higher-order SIBC with curvature correction
   (1/R term in cell problem) would reduce this. At xi >> 10, error < 10%.

2. **FEM-SIBC vs BEM-ESIM (steel, xi >> 20)**: 1% agreement. Both use ESIM Z_s, and
   at strong skin effect the slab model is accurate. Steel is the target application.

3. **BEM-ESIM is 35% below FEM-full** (copper): the "transparent workpiece" assumption
   ignores eddy current feedback. FEM-SIBC includes this feedback via Robin BC.

4. **Power formula**: Must use P = pi (not 2*pi) in 2D axisymmetric Poynting flux.

### Curvature Correction Analysis

0th-order SIBC (flat slab) error depends on xi = R/delta:

| xi | delta/R | Curvature err | Typical case |
|----|---------|--------------|--------------|
| 5 | 0.20 | ~10% | Copper 1 kHz |
| 8 | 0.13 | ~6% | Steel 50 Hz (power transformer) |
| 36 | 0.03 | ~1.4% | Steel 1 kHz |
| 140 | 0.007 | **< 0.4%** | **Steel 7 kHz (induction heating target)** |
| 250 | 0.004 | ~0.2% | Steel 50 kHz |

**Target application**: Steel induction heating at 7 kHz (nonlinear BH).
At xi ~ 140, 0th-order SIBC is sufficient (< 0.4% curvature error).
ESIM is required because Dowell/linear SIBC cannot handle nonlinear B-H curves.

Mitzner 1st-order correction (Yuferev & Ida, 2010):
```
E_u = -Z_s * [1 + (1-j)/2 * delta * (kv - ku)] * H_v
```
where ku, kv are principal curvatures. For cylinder: ku=1/R, kv=0.

Higher-order SIBC is needed only for:
- Power transformers (steel, 50 Hz, xi ~ 8)
- Copper/aluminum conductors at low frequency (xi < 10)

### Remaining Issues

| Issue | Cause | Priority |
|-------|-------|----------|
| +40% P at xi=4.8 | 0th-order SIBC + 2D geometry | Low (not target regime) |
| Per-element Z_s | Average H_t used | Medium (improves P distribution) |
| High-order SIBC | Mitzner curvature 1/R | Low (< 0.4% at target xi=140) |

### Cubit Panel: EM -> q_surf -> Thermal

The induction-heating panel (`radia_ih`) runs the EM solve and the
thermal solve as separate, independent modes (see "Cubit Panel
Architecture" in CLAUDE.md):

1. **EM**: `calc_inductance.py` (PEEC / BEM-A weak coupling: L + Delta_L
   + P_wp), `calc_fem_kelvin.py` (PEEC + FEM-Kelvin + SIBC), or
   `calc_fem_coilmesh.py` (full FEM A-V).  Each emits a workpiece
   surface power density `q_surf` (`*_qsurf.sol`).
2. **Thermal**: `calc_heat.py` / `calc_heat_axisym.py` consume the
   workpiece `.vol` + `q_surf` and solve the heat equation.

All EM solvers read a pre-meshed `.vol` (no auto air-mesh generation);
the Cubit/NGSolve interface is `.vol`-only per CLAUDE.md.

### When Each Method is Best

| Method | Best for | P accuracy | Cost |
|--------|----------|------------|------|
| **FEM-full** | Reference, moderate xi | Exact (mesh-dependent) | Very high (skin mesh) |
| **FEM-SIBC** (Karl) | Any xi, nonlinear BH | Good (xi>10), fair (xi~5) | Low (air mesh only) |
| **FEM-ESIM** (panel) | Design, steel induction heating | Good (xi>10) | Lowest (auto mesh) |
| **BEM-ESIM** | Quick L+P, xi>10 | Good (no feedback) | Low (surface DOFs) |

### Scripts

| Script | Method | Location |
|--------|--------|----------|
| `impedance_esim.py` | BEM-ESIM | `examples/cubit_panels/inductance/` |
| `fem_esim_kelvin.py --mode esim` | FEM-SIBC (Karl) | `examples/cubit_panels/inductance/` |
| `fem_esim_kelvin.py --mode full` | FEM-full | `examples/cubit_panels/inductance/` |
| `verify_esim.py` | ESIM vs NGSolve FEM 1D | `examples/cubit_panels/inductance/` |

### ESIM 1D Cell Problem Verification

The ESIM 1D cell problem itself (independent of BEM/FEM exterior) is verified
against NGSolve H1 FEM (p=4) as independent method:

| Test | ESIM vs | Max error |
|------|---------|-----------|
| Linear Z_s | Analytical rho*gamma*tanh(gamma*a) | 0.25% |
| Linear Z_s | NGSolve H1 FEM (p=4) | 0.04% |
| Nonlinear Z_s (steel BH) | NGSolve H1 FEM + Picard | 0.78% |

### Scattered-Field FEM-SIBC: Missing RHS Term (2026-03-28)

In a **scattered-field** FEM with SIBC (e.g., `verify_sphere_sibc.py`), the RHS
must include BOTH surface terms:

```
f(v) = -(jw/Z_s) * <A_inc, v>_sibc   (surface impedance term)
     + (-1)      * <n x H_inc, v>_sibc (incident field boundary term)
```

The second term `<n x H_inc, v>` arises from the curl-curl integration by parts
and was initially missing, causing a **factor-of-3 error** on the sphere benchmark.

The **total-field formulation** (used in `fem_esim_3d.py` with a coil current source)
does NOT have this issue because the source is `J_source` in the volume:
```
int nu0 * curl(A) . curl(v) dx + (jw/Z_s) * int A_t . v_t ds(wp) = int J . v dx
```

**Rule**: Use total-field formulation for coil+workpiece problems.  Use scattered-field
only when A_inc is known analytically (e.g., uniform field).

### BEM EFIE-SIBC Limitation: SL Eigenvalue Mismatch (2026-03-28)

The EFIE `Z_s*J + jw*mu0*SL(J) = -jw*A_inc` is fundamentally incorrect for finite
Z_s in MQS. The Laplace SL eigenvalue for l=1 on a sphere of radius R is R/3 (not R),
so the BEM EFIE gives an effective denominator `(3*Z_s + jw*mu0*R)` instead of the
correct `(Z_s + jw*mu0*R)`.

This produces BEM/Analytical ratios from 0.35 (Z_s/(jw*mu0*R)=10) to 1.00 (PEC).
Only correct for PEC (Z_s -> 0).

**SOLVED**: Use the **Scalar Potential BIE + SIBC** formulation (Section 6 below),
which achieves <0.1% error for all Z_s using existing ngsolve.bem operators.

See: `examples/cubit_panels/inductance/scalar_bie_sibc.py`

### Literature

SIBC references (local library):

| File | Topic |
|------|-------|
| `Surface Impedance Boundary Conditions a Comprehensive Approach.pdf` | 0th/higher-order SIBC theory |
| `Course G2ELab 2018 - Derivation of High Order SIBCs - lesson 2.pdf` | Curvature correction derivation |
| `Course G2ELab 2018 - FEM Formulation - lesson 4.pdf` | FEM weak form with SIBC |
| `A_Nonlinear_effective_surface_impedance_...pdf` | Hollaus ESIM (nonlinear, Karl iteration) |
| `FEM and BEM implementations of a high order surface impedance...pdf` | High-order SIBC in FEM/BEM |

---

## 5. BEM Formulation Selection for MQS-SIBC (2026-03-28)

### Background

In MQS (Magneto-Quasi-Static) regime, BEM formulations for eddy currents with
Surface Impedance Boundary Condition (SIBC) behave very differently from their
full-wave counterparts. This section summarizes findings from systematic analysis
on a sphere benchmark (analytical solution available).

### EFIE-SIBC in MQS: Limited to PEC-Like Conductors

The MQS EFIE-SIBC formulation:
```
Z_s * M * J + jw * mu0 * SL(J) = -jw * A_inc
```

**Works ONLY when** `Z_s / (jw * mu0 * R) < 0.1` (copper at high frequency).

**Fails for** steel, ferrite, or any material where `Z_s / (jw * mu0 * R) >> 1`.

**Root cause**: In MQS, the scattered vector potential representation
`A_scat = mu0 * SL(J_s)` is incomplete when the workpiece has magnetization
M != 0 (mu_r > 1). The SL operator eigenvalue for l=1 on a sphere of radius R
is R/3, not R. This gives an effective denominator `(3*Z_s + jw*mu0*R)` instead
of the correct `(Z_s + jw*mu0*R)`, producing factor-of-3 errors.

**Why full-wave EFIE-SIBC works**: At GHz frequencies, `Z_s / eta_0 << 1` for
good conductors (eta_0 = sqrt(mu0/eps0) ~ 377 Ohm). The Z_s term is negligible
relative to the Maxwell SL operator, so the eigenvalue mismatch does not matter.
In MQS: `Z_s / (jw*mu0*R)` can be 1-100 for steel, breaking the formulation.

| Z_s/(jw*mu0*R) | BEM/Analytical | Regime |
|-----------------|----------------|--------|
| 0.01 | ~1.00 | PEC-like (copper at high f) |
| 0.1 | ~0.97 | Acceptable |
| 1.0 | ~0.75 | Factor-of-3 effect visible |
| 5.0 | ~0.50 | Significantly wrong |
| 10.0 | ~0.35 | Steel at 7 kHz |

### MFIE in MQS: Two Variants

#### Tangential MFIE (n x H)

The tangential MFIE uses the MFIE K operator:
```
(1/2) * J + K(J) = n x H_inc
```
where `K(J)(x) = n(x) x curl integral_S G(x,y) J(y) dS_y`.

**Key property**: Gives `J = J_pec` for ALL values of Z_s (no Z_s dependence).

In ngsolve.bem, this requires `LaplaceDL` + `HDivSurface` to form the K operator.
On a PEC sphere, this gives 99.8% accuracy for the l=1 mode.

**Why Z_s-independent**: In MQS, magnetization M does not contribute to the
magnetic field H (only to B). The tangential MFIE relates J to H_inc via the
free-space Green's function, and M's contribution to H vanishes in the
Laplace regime.

**Limitation**: Cannot incorporate finite Z_s correction -- always gives the
PEC solution.

#### Normal MFIE (n . H) -- Sugahara Formulation

Uses the normal component of B on the surface:
```
n . B_total = n . B_inc + n . B_scat(J)
```

The SIBC enters via the surface divergence relation:
```
n . B_total = -(Z_s / jw) * curl_s(J) = -(Z_s / jw) * div_s(n x J)
```

**Advantages**:
- n . B_scat computed via Biot-Savart (collocation, not Galerkin)
- SIBC naturally incorporated through surface divergence
- Works for arbitrary Z_s

**Requirements**:
- Biot-Savart matrix assembly (not available as ngsolve.bem operator)
- GRAPH constraint for multiply-connected bodies (loops/holes)
- Collocation method (not Galerkin -- different from ngsolve.bem framework)

**Reference**: Sugahara et al., IEEE TAP, 2006.

### PMCHWT-SIBC in MQS: Impossible

The PMCHWT formulation couples electric and magnetic surface currents (J, M)
via the Stratton-Chu representation. In MQS:

```
jw * eps * SL(M) / ((1/2) * M)  =  O((kR)^2)  ~  10^{-14}
```

The electric displacement current term `jw*eps*SL(M)` that couples M into the
system is 14 orders of magnitude smaller than the identity term. The
`grad*SL*div(M)` correction acts only on the irrotational subspace, orthogonal
to the divergence-free physical solution.

**Weggler stabilization** (low-frequency stabilized BEM): NOT implemented in
current ngsolve.bem. Only standard Laplace/Helmholtz/Maxwell operators are
available. Even with stabilization, the 10^{-14} physics gap on the solenoidal
subspace cannot be bridged.

### Recommended Approach Summary

| Problem | Method | Module/File |
|---------|--------|-------------|
| **Any Z_s (BEM)** | **Scalar BIE + SIBC** | **`scalar_bie_sibc.py`** |
| PEC (Z_s -> 0) | MFIE tangential (LaplaceDL + HDivSurface) | ngsolve.bem |
| PEC quick check | EFIE (LaplaceSL) | `efie_sibc.py` |
| Finite Z_s (nonlinear) | FEM-SIBC + ESIM (Karl iteration) | `fem_esim_3d.py` |
| Circuit extraction (L, R) | Radia PEEC pipeline | `peec_topology.py` |

**Decision tree**:
```
Is Z_s frequency-dependent only (linear material)?
  Yes -> Scalar BIE + SIBC (BEM, surface-only, fast frequency sweep)
  No  -> Is Z_s field-dependent (nonlinear mu)?
           Yes -> FEM-SIBC + ESIM (Karl iteration)
           No  -> Scalar BIE + SIBC
```

---

## 6. Scalar Potential BIE with SIBC (2026-03-28)

### Overview

The Scalar Potential BIE + SIBC is a BEM formulation that correctly handles
finite surface impedance Z_s in the MQS regime. Unlike EFIE-SIBC (which
fails for Z_s/(jw*mu0*R) > 0.1) and MFIE (which gives PEC-only results),
this formulation achieves <0.1% error for ALL Z_s values.

**Key insight**: Use the scalar magnetic potential phi (H = -grad phi) as
the unknown, not the surface current J. The SIBC enters naturally through
the surface Laplacian (Laplace-Beltrami operator), which is well-defined
in the H1 finite element space.

### Formulation

**Exterior problem**: Laplace equation for H = -grad(phi) outside conductor.

**BIE** (Green's representation, exterior limit):
```
(1/2) phi + DL(phi) - SL(g) = phi_inc
```
where g = dphi/dn (Neumann data), DL = scalar double layer, SL = scalar single layer.

**SIBC** (derived from Faraday's law on the surface):
```
E_t = Z_s * J_s = -Z_s * (n x grad_s phi)
Faraday: curl_s(E_t).n = -jw * n.B = jw * mu0 * dphi/dn
Surface identity: curl_s(n x grad_s f).n = Delta_s f

=> dphi/dn = -(Z_s / (jw*mu0)) * Delta_s(phi)
```

**Combined system** (single unknown phi in H1):
```
(1/2*M - DL + gamma * SL * M_inv * K) phi = rhs

M  = H1 surface mass matrix
K  = H1 surface stiffness (Laplace-Beltrami)
DL = LaplaceDL (scalar, ngsolve.bem)
SL = LaplaceSL (scalar, ngsolve.bem)
gamma = Z_s / (jw * mu0)
```

**Sign convention**: NGSolve `LaplaceDL` has (1/2*M - DL) for the exterior
Dirichlet problem. The SIBC correction enters with **positive** gamma.

**Gauge**: Lagrange multiplier for int(phi) dS = 0 fixes the constant null space.

### Why It Works (and Why EFIE/MFIE Failed)

| Issue | EFIE-SIBC | MFIE | Scalar BIE + SIBC |
|-------|-----------|------|-------------------|
| SL eigenvalue R/3 | Causes 3x Z_s error | N/A | N/A (no SL*J term) |
| RT0 curl_s = 0 | N/A | No Z_s dependence | Uses H1 stiffness K |
| M -> H = 0 in MQS | Incomplete A_scat | No coupling | N/A (scalar pot.) |

The scalar potential phi naturally encodes both tangential H (via grad_s phi)
and normal B (via dphi/dn). The SIBC couples them through the surface Laplacian K,
which is nonzero for H1 order >= 1.

### Verification: Sphere Benchmark

Sphere (R=10mm, B0=1mT, f=1kHz), mesh maxh=R/5 (722 elements, 363 H1 DOFs):

| Z_s/(jw*mu0*R) | Analytical H_rms | Scalar BIE H_rms | Error |
|-----------------|------------------|------------------|-------|
| 0 (PEC) | 974.62 | 975.47 | +0.09% |
| 0.001 | 973.93 | 974.78 | +0.09% |
| 0.01 | 967.75 | 968.58 | +0.09% |
| 0.1 | 908.28 | 908.95 | +0.07% |
| 1.0 | 527.46 | 527.46 | +0.00% |
| 10.0 | 90.83 | 90.76 | -0.07% |

Errors decrease with mesh refinement (R/3: 0.25%, R/5: 0.09%).

### Usage

```python
from radia.scalar_bie_sibc import ScalarBIE_SIBC

# 1. Create solver (assembles BEM operators once)
solver = ScalarBIE_SIBC(mesh, order=1)

# 2. Solve for given H_inc, Z_s, omega
result = solver.solve(H_inc_cf, Zs=0.01+0.01j, omega=2*pi*1e3)
print(f"H_rms = {result['H_rms']:.2f} A/m")
print(f"P_loss = {result['P_loss']:.4f} W")

# 3. Frequency sweep (Z_s recomputed per frequency)
results = solver.frequency_sweep(H_inc_cf, freqs=[1e3, 10e3, 100e3],
                                  sigma=5.8e7)  # copper
```

### Advantages

- **All Z_s**: Works from PEC to highly resistive (steel, ferrite)
- **Surface-only**: No volume mesh needed (BEM)
- **Fast frequency sweep**: DL and SL are frequency-independent;
  only the gamma*SL*M_inv*K term changes with frequency
- **Existing operators**: No new C++ code; uses ngsolve.bem LaplaceSL/DL
- **Small DOF count**: H1 on surface (~N_vertices vs ~N_edges for HDivSurface)

### Limitations

- **Linear Z_s only**: For nonlinear mu(H), use FEM-SIBC + ESIM Karl iteration
- **Dense matrices**: Current implementation is O(N^2) memory, O(N^3) solve.
  For large meshes, H-matrix compression (HACApK) could be applied to SL/DL.
- **Closed surfaces**: The gauge condition assumes a closed surface.
  Open conductors need a different formulation.
- **H_inc reconstruction**: When phi_inc is not known analytically,
  it must be recovered from H_inc via surface Poisson (adds one extra solve)

### Files

| File | Purpose |
|------|---------|
| `src/radia/scalar_bie_sibc.py` | Core solver module (ScalarBIE_SIBC class) |
| `examples/cubit_panels/inductance/scalar_bie_sibc.py` | Sphere verification script |

---

## References

- **EMPY T-Omega method**: `S:\NGSolve\EMPY\EMPY_Analysis\EddyCurrent` (T_Omega_Method.py)
- **EMPY A-Phi method**: `S:\NGSolve\EMPY\EMPY_Analysis\SolverRun.py`
- **Weggler stabilization**: [Maxwell_DtN_Stabilized.ipynb](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb)
- **Weggler thesis**: L. Weggler, "High order boundary element methods", PhD thesis, Saarland University (2011)
- **ngbem**: [github.com/Weggler/ngbem](https://github.com/Weggler/ngbem) (MIT license)
- **Johnson-Nedelec coupling**: C. Johnson & J.C. Nedelec, "On the coupling of BEM and FEM", Math. Comp., 1980
- **Costabel symmetric coupling**: M. Costabel, "Symmetric methods for the coupling of FEM and BEM", 1987
- **Kelvin transform**: T. Kuwahara & T. Takeda, "Unbounded FEM using Kelvin transformation", 1990
- **Hollaus ESIM**: K. Hollaus et al., "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation", IEEE Trans. Magnetics, 2025
- **Sugahara normal MFIE**: Sugahara et al., "Electromagnetic analysis using the GRAPH method", IEEE Trans. Antennas Propagat., 2006
