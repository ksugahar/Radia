# Conductor Eddy Current Modeling: Method Comparison

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

### MCP Server Tools

The NGSolve MCP server provides eddy current tools based on the EMPY implementation:

```
eddy_create_conductor_mesh  -- Create conductor + air mesh
eddy_solve_a_phi            -- A-Phi formulation solver
eddy_solve_t_omega          -- T-Omega formulation solver
eddy_compute_loss           -- Compute eddy current loss
eddy_extract_impedance      -- Extract R, L, Z
eddy_frequency_sweep        -- Frequency sweep analysis
eddy_export_vtk             -- Export B, J, H to VTK
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

```
                    Start
                      |
            Need eddy currents?
                   /       \
                 Yes        No
                  |          |
                  |     Kelvin Transform
                  |     (magnetostatic Omega-Reduced Omega)
                  |
         Nonlinear mu_r(H)?
              /          \
            Yes           No
             |             |
    FEM + Kelvin           |
    (A-Phi or T-Omega,     |
     Newton iteration)     |
                           |
                  Need air mesh free?
                      /         \
                    Yes          No
                     |            |
               ngbem BEM     FEM + Kelvin
                     |       (A-Phi: simple)
                     |       (T-Omega: multi-loop)
                     |
            Need frequency sweep?
                  /           \
                Yes            No
                 |              |
        BEM + SIBC         FEM-BEM
        (fast sweep)    (exact interior)
                 |
        PEEC coupling needed?
             /       \
           Yes        No
            |          |
     ShieldBEMSIBC   EddyCurrentBEMSIBC
     (loop basis,    (scalar Hz,
      Biot-Savart)    Costabel)
```

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
| VectorEddyCurrentFEMBEM | ~164 | 94 H(curl) + 42 HDivSurf + 28 SurfL2 |
| EddyCurrentFEMBEM | ~100 | H1 interior + H1/2 boundary |
| ShieldBEMSIBC | ~174 loops | From V-1 (genus-0 surface) |
| FEM + Kelvin (A-Phi) | ~5000+ | HCurl(all) + H1(conductor) |
| FEM + Kelvin (T-Omega) | ~3000+ | HCurl(conductor) + H1(all) + loops |

## Required ngbem/ngsolve.bem Operators

| Method | Operators Used | Package |
|---|---|---|
| EddyCurrentFEMBEM | V (SLP), K (DLP), D (Hypersingular) | ngbem or ngsolve.bem |
| VectorEddyCurrentFEMBEM | HelmholtzSL only | ngsolve.bem |
| EddyCurrentBEMSIBC | V, K, D (cached) | ngbem or ngsolve.bem |
| ShieldBEMSIBC | MaxwellSingleLayerPotentialOperator | ngbem |
| FEM + Kelvin | (none -- pure FEM) | NGSolve only |

## MCP Server Tool Summary

| Tool Category | Prefix | Eddy Current | Tools |
|---|---|---|---|
| **Mesh generation** | `ngsolve_mesh_` | N/A | Create/import mesh |
| **Kelvin (magnetostatic)** | `kelvin_` | No | Omega-Reduced Omega, energy |
| **Eddy current (FEM+Kelvin)** | `eddy_` | **Yes** | A-Phi, T-Omega, loss, Z, sweep |
| **ngbem (BEM)** | `ngbem_` | **Yes** | Maxwell stabilized, FEM-BEM coupling |
| **Radia coupling** | `ngsolve_radia_` | No | Magnetostatic field coupling |

## Cross-Method Validation (2026-02-15)

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

### Method Selection Guide (Updated)

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

## References

- **EMPY T-Omega method**: `S:\NGSolve\EMPY\EMPY_Analysis\EddyCurrent` (T_Omega_Method.py)
- **EMPY A-Phi method**: `S:\NGSolve\EMPY\EMPY_Analysis\SolverRun.py`
- **Weggler stabilization**: [Maxwell_DtN_Stabilized.ipynb](https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb)
- **ngbem**: [github.com/Weggler/ngbem](https://github.com/Weggler/ngbem) (MIT license)
- **Johnson-Nedelec coupling**: C. Johnson & J.C. Nedelec, "On the coupling of BEM and FEM", Math. Comp., 1980
- **Costabel symmetric coupling**: M. Costabel, "Symmetric methods for the coupling of FEM and BEM", 1987
- **Kelvin transform**: T. Kuwahara & T. Takeda, "Unbounded FEM using Kelvin transformation", 1990
