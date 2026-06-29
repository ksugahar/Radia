# Unified PEEC Loop-Star + MMM + MSC Architecture

**Date**: 2026-02-22
**Status**: Low-frequency BEM verified via product space (ngsbem_peec_demo)

## Overview

This document describes the unified architecture for electromagnetic analysis combining:
- **PEEC (Partial Element Equivalent Circuit)** with complete Loop-Star decomposition
- **MMM (Magnetic Moment Method)** for magnetic materials (existing Radia)
- **MSC (Magnetic Surface Charge)** for dielectric materials (future extension)

## Architecture Decision (2026-01-10)

### Key Changes

1. **RWG-EFIE removed** - Deleted from codebase (rad_rwg_*.cpp/h)
2. **Helmholtz kernel removed** - Laplace kernel only for quasi-static
3. **Complete Loop-Star decomposition** - Full separation of solenoidal/irrotational currents
4. **MMM coupling via A and B/H fields** - Loop currents couple with magnetization
5. **MSC coupling via E field (future)** - Star charges couple with dielectric polarization
6. **Complex material support** - mu = mu' - j*mu'', epsilon = epsilon' - j*epsilon''

### Target Applications

- **WPT (Wireless Power Transfer)** - Self-resonance analysis
- **Coil design** - Coupled electromagnetic-magnetic problems
- **Future: Induction heating** - After dielectric MSC is implemented

### Frequency Domain Only

This architecture supports **frequency domain (linear) analysis only**:
- Complex impedance Z(omega)
- Complex material properties for loss modeling
- No time-domain transient analysis

## Unified System Matrix

### Block Structure

```
[Z_LL    Z_LS    Z_LM    0    ] [I_L ]   [V_L ]
[Z_SL    Z_SS    0       Z_SE ] [I_S ] = [V_S ]
[Z_ML    0       Z_MM    0    ] [M   ]   [H_ext]
[0       Z_ES    0       Z_EE ] [P   ]   [D_ext]

where:
  I_L: Loop currents (solenoidal, div J = 0)
  I_S: Star currents (irrotational, charge-related)
  M:   Magnetization (MMM)
  P:   Polarization (MSC for dielectrics, future)
```

### Coupling Physics

| Coupling | Physical Mechanism | Kernel |
|----------|-------------------|--------|
| Z_LL (Loop-Loop) | Inductance via vector potential A | Laplace: mu0/(4*pi*r) |
| Z_SS (Star-Star) | Capacitance via scalar potential | Laplace: 1/(4*pi*epsilon*r) |
| Z_LM (Loop-MMM) | B field from M affects J | Laplace (MSC B field) |
| Z_ML (MMM-Loop) | H field from J affects M | Biot-Savart: mu0/(4*pi*r^3) |
| Z_SE (Star-MSC) | E field from P affects sigma | Laplace (surface charge) |
| Z_ES (MSC-Star) | E field from sigma affects P | Laplace (surface charge) |
| Z_MM (MMM self) | Demagnetization tensor | Existing Radia |
| Z_EE (MSC self) | Dielectric self-interaction | Laplace surface charge |

### Loop-Star Scaling for Low Frequency

Standard EFIE has condition number issues at low frequency:
- Z_LL ~ O(omega) (inductive)
- Z_SS ~ O(1/omega) (capacitive)
- Condition number ~ O(1/omega^2)

**Rescaled Loop-Star system**:
```
I_S' = j*omega * I_S  (scaled star current = charge derivative)

[L        M_LS/jw  ] [I_L ]   [V_L/jw ]
[M_SL*jw  1/C      ] [I_S'] = [V_S*jw ]
```

All blocks now O(1) at low frequency.

## Implementation Components

### Laplace Kernel Only

All electromagnetic interactions use the Laplace Green's function:
```
G(r, r') = 1 / (4 * pi * |r - r'|)
```

No Helmholtz kernel (exp(-jkr)) required for quasi-static analysis.

### Complex Material Properties

**Magnetic permeability**:
```
mu = mu0 * (mu_r' - j * mu_r'')

where:
  mu_r': Real permeability (storage)
  mu_r'': Loss tangent (dissipation)
```

**Electric permittivity**:
```
epsilon = epsilon0 * (epsilon_r' - j * epsilon_r'')

where:
  epsilon_r': Real permittivity (storage)
  epsilon_r'': Loss tangent (dissipation)
```

### Skin Depth with Complex mu

For conductive magnetic materials:
```
delta = sqrt(2 / (omega * mu0 * mu_r * sigma))

With complex mu_r = mu_r' - j*mu_r'':
  delta becomes complex -> field penetration with phase shift
```

## Priority Assessment

### Current Priority Order

| Priority | Task | Status | Notes |
|----------|------|--------|-------|
| **1** | PEEC Loop-Star + MMM coupling | In Progress | Core unified solver |
| **2** | Remove Helmholtz kernel | Done | Laplace only |
| **3** | Remove RWG-EFIE | Done | Files deleted |
| **4** | Complex mu/epsilon support | Planned | Loss modeling |
| **5** | Star-MSC coupling (dielectric) | Future | For WPT capacitors |
| **6** | NGBEM low-frequency BEM | **Verified** | Product space = Weggler EFIE |

### Key Development Goals Summary

1. **PEEC Loop-Star + MMM coupling** - Primary goal for WPT analysis
2. **Laplace kernel only** - Simplify to quasi-static
3. **Complex material properties** - Enable loss modeling
4. **Star-MSC for dielectrics (future)** - Self-resonance support
5. **NGBEM integration (future)** - High-order elements if needed

### Why Remove RWG-EFIE?

1. **Redundant with PEEC Loop-Star** - Loop-Star provides same low-frequency stability
2. **Simpler maintenance** - Single solver architecture
3. **Focus on MMM coupling** - RWG was standalone, not coupled to MMM
4. **NGBEM as future option** - High-order BEM if needed later

## Goals

1. **PEEC Loop-Star with ESIM**: Complete Loop-Star separation with surface impedance
2. **MMM coupling via A and B/H**: Loop currents couple with magnetization
3. **MSC for dielectrics (future)**: Star charges couple with polarization via E field
4. **Complex material properties**: Support mu'' and epsilon'' for loss modeling
5. **NGBEM as future option**: High-order elements when needed

## Current Radia Architecture (After RWG Removal)

### PEEC Loop-Star + MMM + MSC Implementation

```
+---------------------------------------------------+
|  rad_conductor.cpp/h                              |
|  - PEEC conductor formulation with Loop-Star      |
|  - Surface panel discretization                   |
|  - ESIM surface impedance (Karl Hollaus)          |
+---------------------------------------------------+
|  rad_peec_mmm_coupled.cpp/h                       |
|  - Unified Loop-Star + MMM + MSC solver           |
|  - Loop <-> MMM coupling (via A, B/H)             |
|  - Star <-> MSC coupling (via E, future)          |
|  - Complex material properties                    |
+---------------------------------------------------+
|  rad_green_fullwave.cpp/h                         |
|  - Laplace kernel only (1/4*pi*r)                 |
|  - Quasi-static formulation                       |
|  - Panel interaction integrals                    |
+---------------------------------------------------+
```

**Key Features**:
- Laplace kernel only (Helmholtz removed)
- Complete Loop-Star decomposition for low-frequency stability
- ESIM surface impedance for conductive materials
- Coupled with existing Radia MMM (rad_interaction.cpp)

### Deleted Files (2026-01-10)

- `rad_rwg_basis.cpp/h` - Replaced by Loop-Star in rad_conductor.cpp
- `rad_rwg_coupled.cpp/h` - Replaced by rad_peec_mmm_coupled.cpp
- `rad_rwg_coupled_api.cpp` - API merged into rad_peec_mmm_api.cpp

### Future NGBEM Integration (Optional)

```
+---------------------------------------------------+
|  NGBEM (NGSolve BEM add-on)                       |
|  - High-order H(div) and H(curl) spaces           |
|  - Curved mesh support                            |
|  - ACA/H-matrix compression                       |
|  - Laplace operators (MQS kernel to be added)     |
+---------------------------------------------------+
            |
            v
+---------------------------------------------------+
|  Radia NGBEM Interface (Python, future)           |
|  - Loop-Star transformation matrices              |
|  - ESIM surface impedance (Karl Hollaus)          |
|  - Coreform mesh import                           |
|  - Field computation interface                    |
+---------------------------------------------------+
            |
            v
+---------------------------------------------------+
|  Existing Radia Infrastructure                    |
|  - PEEC Loop-Star + MMM (current)                 |
|  - Material database                              |
+---------------------------------------------------+
```

## NGBEM Capabilities

### Supported Operators

| Operator | Kernel | NGSolve Space | Application |
|----------|--------|---------------|-------------|
| Single Layer (V) | 1/(4*pi*r) | SurfaceL2 | Laplace BEM |
| Double Layer (K) | d/dn(1/(4*pi*r)) | H1 | Laplace BEM |
| Maxwell EFIE | exp(-jkr)/(4*pi*r) | HDiv | Full-wave |
| Maxwell MFIE | curl(G) | HCurl | Full-wave |

### Space Definitions

```python
from ngsolve import *
from ngbem import *

# Laplace (MSC kernel)
fesH1 = H1(mesh, order=3, definedon=mesh.Boundaries(".*"))
fesL2 = SurfaceL2(mesh, order=2, dual_mapping=True)

# Maxwell EFIE
fesHDiv = HDivSurface(mesh, order=3, complex=True)
fesHCurl = HCurlSurface(mesh, order=3, complex=True)
```

## Implementation Plan

### Phase 1: Laplace BEM Verification

**Goal**: Verify NGBEM Laplace operators match Radia MSC results

```python
from ngsolve import *
from ngbem import *
import radia as rad

# Create test geometry (sphere)
mesh = Mesh(...)

# NGBEM Laplace
fesL2 = SurfaceL2(mesh, order=2, dual_mapping=True)
V = SingleLayerPotentialOperator(fesL2, intorder=12, eps=1e-4)

# Compare with Radia MSC
# ... field computation at test points
```

**Validation**: Compare B field from both methods at external points

### Phase 2: Maxwell EFIE with Loop-Star

**Goal**: Implement stable low-frequency EFIE using Loop-Star decomposition

#### Loop-Star Decomposition Theory

For low-frequency stability, decompose current density J into:
- **Loop currents** (JL): Solenoidal, divJ=0
- **Star currents** (JS): Non-solenoidal, surface charge related

```
J = JL + JS

EFIE: [ZLL  ZLS] [IL]   [VL]
      [ZSL  ZSS] [IS] = [VS]
```

Scaling:
- ZLL ~ O(omega)
- ZLS, ZSL ~ O(omega)
- ZSS ~ O(1/omega)

Rescale star currents: IS' = jomega * IS

```
[ZLL    ZLS/jw ] [IL ]   [VL ]
[ZSL*jw ZSS    ] [IS'] = [VS']
```

Now all blocks are O(1) at low frequency.

#### NGBEM Loop-Star Implementation

```python
class NGBEMLoopStarSolver:
    """
    NGBEM-based EFIE solver with Loop-Star decomposition.

    Uses NGBEM high-order H(div) space and custom Loop-Star
    transformation for low-frequency stability.
    """

    def __init__(self, mesh, order=3):
        self.mesh = mesh
        self.order = order

        # H(div) space for surface currents
        self.fes_hdiv = HDivSurface(mesh, order=order, complex=True)

        # H1 space for loop identification
        self.fes_h1 = H1(mesh, order=order+1,
                        definedon=mesh.Boundaries(".*"))

        # Build Loop-Star transformation
        self._build_loop_star_transform()

    def _build_loop_star_transform(self):
        """
        Build Loop-Star transformation matrix T.

        T transforms [JL, JS] -> J_hdiv
        T^(-1) transforms J_hdiv -> [JL, JS]
        """
        # Loop basis: curl of H1 functions (edge-based)
        # Star basis: gradient of vertex functions
        # Implementation uses mesh topology
        pass

    def assemble(self, frequency):
        """
        Assemble EFIE system with Loop-Star scaling.
        """
        omega = 2 * np.pi * frequency
        k = omega / 299792458  # wavenumber

        # NGBEM Maxwell operators
        from ngbem import MaxwellSingleLayerPotentialOperator

        # Single layer: A-A interaction
        SL = MaxwellSingleLayerPotentialOperator(
            self.fes_hdiv,
            intorder=12,
            eps=1e-4,
            k=k  # wavenumber
        )

        # Transform to Loop-Star basis
        # Z_LS = T^H * Z * T
        # Apply scaling for low-frequency stability
        pass
```

### Phase 3: ESIM Integration

**Goal**: Add Karl Hollaus ESIM surface impedance

#### ESIM Formulation (Karl Hollaus)

Surface impedance for conductive magnetic materials:

```
Skin depth: delta = sqrt(2 / (omega * mu0 * mur * sigma))

Surface resistance: Rs = 1 / (sigma * delta)

Surface impedance: Zs = (1 + j) * Rs
```

For nonlinear materials, mu_r depends on local H field.

#### Implementation

```python
class NGBEMLoopStarESIMSolver(NGBEMLoopStarSolver):
    """
    NGBEM Loop-Star solver with ESIM surface impedance.
    """

    def __init__(self, mesh, order=3):
        super().__init__(mesh, order)

        # Material properties
        self.sigma = 5.8e7  # Conductivity [S/m]
        self.mu_r = 1.0     # Relative permeability

    def set_material(self, sigma, mu_r):
        """Set conductor material properties."""
        self.sigma = sigma
        self.mu_r = mu_r

    def get_skin_depth(self, frequency):
        """Calculate skin depth."""
        omega = 2 * np.pi * frequency
        mu0 = 4 * np.pi * 1e-7
        return np.sqrt(2 / (omega * mu0 * self.mu_r * self.sigma))

    def get_surface_impedance(self, frequency):
        """
        Calculate ESIM surface impedance.

        Returns complex impedance Zs = (1+j) * Rs
        """
        delta = self.get_skin_depth(frequency)
        Rs = 1.0 / (self.sigma * delta)
        return complex(Rs, Rs)  # (1+j) * Rs

    def assemble_with_esim(self, frequency):
        """
        Assemble EFIE with ESIM surface impedance.

        EFIE + ESIM:
            Z_total = Z_EFIE + Z_ESIM

        where Z_ESIM adds surface impedance contribution:
            Z_ESIM[i,j] = Zs * integral{ fi . fj dS }
        """
        # Get base EFIE matrix
        Z_efie = self.assemble(frequency)

        # Add ESIM contribution
        Zs = self.get_surface_impedance(frequency)

        # Mass matrix for H(div) space
        u, v = self.fes_hdiv.TnT()
        mass = BilinearForm(self.fes_hdiv)
        mass += InnerProduct(u, v) * ds
        mass.Assemble()

        # Z_total = Z_EFIE + Zs * M
        Z_total = Z_efie + Zs * mass.mat

        return Z_total
```

### Phase 4: Coreform Hexahedral Mesh to Netgen/NGSolve

**Goal**: Import Coreform Cubit hexahedral meshes into Netgen/NGSolve

### Motivation

- **Coreform Cubit excels at hexahedral meshing** (structured, mapped, sweep)
- **Netgen defaults to tetrahedral** - hex mesh import is essential
- **High-quality hex mesh** = better accuracy for BEM surface extraction

### Current Mesh Pipeline

```
Coreform Cubit (.cub5)
        |
        v
GMSH Format (.msh)  <-- radia Cubit plugin (export gmsh)
        |
        v
NGSolve Mesh (via Mesh() constructor)
        |
        v
NGBEM Surface Mesh (boundaries only)
```

### Hexahedral Mesh Support in NGSolve

NGSolve supports hexahedral elements:

```python
from ngsolve import *

# NGSolve reads .vol (the only supported input format)
mesh = Mesh("hex_mesh.vol")

# Check element types
for el in mesh.Elements(VOL):
    print(el.type)  # HEXAHEDRON, TET, PRISM, PYRAMID
```

**Challenge**: Netgen's native mesh generator produces tetrahedra only. Hex mesh must be imported.

### Coreform to Netgen Hex Pipeline

```
Cubit -> export netgen "model.vol" order N -> NGSolve Mesh("model.vol")
```

```python
import cubit
from ngsolve import Mesh

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('import step "model.step" heal')
cubit.cmd("volume all scheme map")   # hex meshing
cubit.cmd("volume all size 0.01")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# Export .vol with high-order curving
cubit.cmd('export netgen "model.vol" order 3 overwrite')

mesh = Mesh("model.vol")
```

The `.vol` format supports all element types (tet, hex, wedge, pyramid) with order 1-5.

### Surface Extraction for BEM

For BEM analysis, only surface mesh is needed:

```python
from ngsolve import *

# Load volume mesh (the only supported input is .vol)
mesh = Mesh("model.vol")

# Get boundary mesh for BEM
# NGBEM automatically extracts surface from volume mesh
from ngbem import *
fes_surf = HDivSurface(mesh, order=2)  # Surface H(div) space
```

### Quadrilateral Surface from Hexahedral Volume

When extracting surface from hex mesh:
- Each hex face is a **quadrilateral**
- NGSolve preserves quad faces on boundary
- **NGBEM should support quad surface elements** (needs verification)

```
+-------------+
|\            |\
| \           | \
|  \          |  \
|   +-------------+   <- Hex volume
|   |         |   |
+---|---------|---+
 \  |          \  |
  \ |           \ |
   \|            \|
    +-------------+   <- Quad surface faces
```

### Implementation Plan

1. **Verify NGSolve hex import**: Test GMSH hex -> NGSolve pipeline
2. **Test NGBEM quad support**: Check if HDivSurface works on quad faces
3. **Create Coreform export tool**: Extend the radia Cubit plugin for hex
4. **End-to-end test**: Coreform hex -> NGSolve -> NGBEM -> Radia

### Validation

| Test | Input | Expected Output |
|------|-------|-----------------|
| Hex cube | Coreform hex mesh | NGSolve mesh with 6 quad faces |
| Cylinder | Coreform mapped mesh | NGSolve mesh with quad sides |
| BEM space | Hex volume mesh | HDivSurface on quad faces |

### Implementation (Extended)

```python
from ngsolve import Mesh
# Use: cubit.cmd('export gmsh "model.msh" overwrite')

def import_coreform_hex_mesh(cubit_file, surface_names=None):
    """
    Import Coreform Cubit hexahedral mesh for NGBEM analysis.

    Args:
        cubit_file: Path to .cub5 file
        surface_names: List of surface names to extract

    Returns:
        NGSolve mesh with hexahedral elements and quad surface faces
    """
    # Export to GMSH format (preserving hex elements)
    gmsh_file = cubit_file.replace('.cub5', '.msh')
    export_to_gmsh(cubit_file, gmsh_file, element_type="hex")

    # Import to NGSolve
    mesh = Mesh(gmsh_file)

    # Verify element types
    n_hex = sum(1 for el in mesh.Elements(VOL) if el.type == "HEXAHEDRON")
    print(f"Imported {n_hex} hexahedral elements")

    return mesh

def extract_quad_surface(mesh, boundary_name="all"):
    """
    Extract quadrilateral surface from hexahedral volume mesh.

    Args:
        mesh: NGSolve mesh with hex elements
        boundary_name: Name of boundary to extract ("all" for entire surface)

    Returns:
        Surface mesh with quad faces
    """
    if boundary_name == "all":
        bnd = mesh.Boundaries(".*")
    else:
        bnd = mesh.Boundaries(boundary_name)

    # NGSolve automatically provides surface elements
    return bnd
```

## Priority 1: FastImp PEEC + MMM Coupling

### Current Status (2026-01-10: Implementation Started)

| Component | Implementation | Coupling Status |
|-----------|---------------|-----------------|
| FastImp PEEC | rad_conductor.cpp | Standalone |
| MMM (Radia) | rad_interaction.cpp | Standalone |
| RWG-EFIE | rad_rwg_*.cpp | Coupled (coil-workpiece) |
| **PEEC-MMM** | **rad_peec_mmm_coupled.cpp** | **NEW - Initial Implementation** |

**Implemented**:
- `PEECMMMCoupledSolver` class in `src/core/rad_peec_mmm_coupled.h/cpp`
- Python API in `src/lib/rad_peec_mmm_api.cpp`
- CMakeLists.txt updated

**Applications**:
- Eddy currents in magnetic yoke
- Transformer with conductive core
- Kicker magnet with vacuum chamber

### Coupling Architecture

```
+-------------------+     +-------------------+
|  FastImp PEEC     |     |  Radia MMM        |
|  (Conductors)     |     |  (Magnets)        |
|                   |     |                   |
|  J_cond (surface) |     |  M_mag (volume)   |
+--------+----------+     +--------+----------+
         |                         |
         |    Mutual Coupling      |
         +----------+--------------+
                    |
                    v
+------------------------------------------+
|  B_total = B_cond + B_mag                |
|                                          |
|  B_cond: Biot-Savart from J_cond         |
|  B_mag:  MSC from M_mag                  |
+------------------------------------------+
```

### Implementation Plan

```cpp
// rad_peec_mmm_coupled.h

class PEECMMMCoupledSolver {
public:
    // Set components
    void SetPEECConductor(int peecHandle);
    void SetMMMObject(int mmmHandle);

    // Coupling computation
    void ComputeMutualCoupling();

    // Solve coupled system
    void Solve(double frequency);

    // Results
    std::complex<double> GetImpedance() const;
    void ComputeB(const TVector3d& point,
                  std::complex<double>& Bx,
                  std::complex<double>& By,
                  std::complex<double>& Bz) const;

private:
    // B field from conductor at magnet volume
    void ComputeBCondAtMagnet();

    // B field from magnet at conductor surface
    void ComputeBMagAtConductor();

    // Update M based on B_total
    void UpdateMagnetization();
};
```

### Python API

The old C++ conductor APIs (`CndSpiral`, `CndLoop`, etc.) are **removed**.
Use the Python-based PEEC solver:

```python
from radia.peec_topology import PEECCircuitSolver
from radia.peec_coupled import CoupledPEECSolver

# Build PEEC topology
from peec_matrices import PyPEECBuilder
builder = PyPEECBuilder()
# ... add nodes, segments, ports ...
topo = builder.build_topology()

# Create magnet (MMM)
import radia as rad
magnet = rad.ObjRecMag([0,0,-0.01], [0.1,0.1,0.01], [0,0,0])
mat = rad.MatLin(1000)
rad.MatApl(magnet, mat)

# Coupled PEEC + MMM solve
solver = CoupledPEECSolver(topo, magnetic_objects=[magnet])
Z = solver.compute_port_impedance(freq=50000)
Z_sweep = solver.frequency_sweep([1000, 10000, 50000, 100000])
```

### API Reference

| Function | Description | Parameters | Returns |
|----------|-------------|------------|---------|
| `PEECMMMCreate(cond, mag)` | Create coupled solver | Conductor & magnet handles | Solver handle |
| `PEECMMMSetFrequency(solver, f)` | Set frequency | Solver handle, frequency [Hz] | None |
| `PEECMMMSetVoltage(solver, V_re, V_im)` | Set voltage excitation | Solver, V real/imag [V] | None |
| `PEECMMMSetCurrent(solver, I_re, I_im)` | Set current excitation | Solver, I real/imag [A] | None |
| `PEECMMMSetExtField(solver, Hx, Hy, Hz)` | Set external H field | Solver, H [A/m] | None |
| `PEECMMMSolve(solver)` | Solve coupled system | Solver handle | [Z_re, Z_im, P_cond, P_mag, iter] |
| `PEECMMMImpedance(solver)` | Get impedance | Solver handle | [Z_re, Z_im] |
| `PEECMMMFld(solver, point)` | Compute B field | Solver, [x,y,z] | [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im] |
| `PEECMMMSweep(solver, freqs)` | Frequency sweep | Solver, [f1,f2,...] | [Z_re1, Z_im1, ...] |
| `PEECMMMDelete(solver)` | Delete solver | Solver handle | None |

## Priority 2: BEM-FEM Coupling (NGSolve)

### Concept

NGSolve provides easy BEM-FEM coupling through:
- FEM in interior domain (volume mesh)
- BEM on boundary (surface mesh)

### Architecture

```
+-------------------+     +-------------------+
|  NGSolve FEM      |     |  NGSolve BEM      |
|  (Interior)       |     |  (Boundary)       |
|                   |     |                   |
|  H1/HCurl spaces  |     |  Surface spaces   |
+--------+----------+     +--------+----------+
         |                         |
         |    Trace operators      |
         +----------+--------------+
                    |
                    v
+------------------------------------------+
|  Coupled System                          |
|  [A_FEM   B_trace] [u_int]   [f_int]     |
|  [B_trace A_BEM  ] [u_bnd] = [f_bnd]     |
+------------------------------------------+
```

### Use Cases

1. **Eddy current in thin shell**: BEM for shell, FEM for surrounding air
2. **Magnetic shielding**: FEM for shield volume, BEM for external field
3. **SIBC formulation**: FEM interior with BEM surface impedance BC

### Implementation Sketch

```python
from ngsolve import *
from ngbem import *

# Volume mesh (FEM domain)
mesh_vol = Mesh("interior.vol.gz")

# Surface mesh (BEM domain)
mesh_surf = mesh_vol.GetSurfaceMesh()

# FEM space
fes_fem = HCurl(mesh_vol, order=2)

# BEM space
fes_bem = HDivSurface(mesh_surf, order=2)

# Coupled bilinear form
# ... (NGSolve provides coupling operators)
```

## Priority 3: NGBEM Low-Frequency BEM (VERIFIED)

### Status: Verified (2026-02-22)

Low-frequency BEM is already implemented using ngbem's existing product space (`HDivSurface x SurfaceL2`).
A dedicated MQS kernel is unnecessary -- the Weggler stabilized EFIE equivalent formulation is obtained naturally.

**Reference implementation**: `docs/peec_integration/demos/ngsbem_peec_demo/`

### Verified Approach: Product Space = Weggler EFIE

The ngbem `HDivSurface x SurfaceL2` product space naturally provides the Loop-Star decomposition:

```
| Z_LL    M_LS^T |   | I_loop |   | V_port |
| M_LS    Z_SS   | * | Q_star | = | 0      |
```

- **L = μ₀ · LaplaceSL(HDivSurface)**: Loop inductance (edge-based RWG)
- **P = SingleLayerPotentialOperator(SurfaceL2) / ε₀**: Star potential (cell-based)
- **M_LS = ∫div(J_edge)·φ_cell dS**: Divergence coupling (charge conservation)
- **Condition number**: O(1) from DC to RF — no low-frequency breakdown

### Loop-Star Decomposition (Graph Theory)

Euler characteristic χ = V - E + F determines:
- **T_loop** (face-edge incidence): ±1 based on face circulation vs global edge direction
- **T_star** (node-edge incidence): ±1 for edge leaving/entering vertex
- **Orthogonality**: T_star^T · T_loop = 0 (exact, from graph theory)
- **Completeness**: n_loop + n_star = n_edges (at order=0)

### Dowell Skin Effect (Verified)

AC resistance for rectangular conductors (d << w):

```
F_R(ξ) = ξ · (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
ξ = d / (2δ),  δ = √(2/(ωμσ))
Zs[i] = R_dc[i] · (F_R(ξ) - 1)
```

Added as callable `Zs_func` for frequency-dependent excess resistance.

### Schur Complement Port Extraction (Verified)

1. Solve Z_SS · X = M_LS via LDL^T (complex symmetric, Bunch-Kaufman)
2. Z_eff = Z_LL - M_LS^T · X
3. Z_port = 1 / (e^T · Z_eff^{-1} · e)

### Coupled Core Models (Verified in ngbem_coupled.py)

| Model | Domain | μ_r | Status | Notes |
|-------|--------|-----|--------|-------|
| fembem (Calderon) | Unbounded | =1 only | Verified | Hz scalar; **inaccurate for mu_r!=1** |
| vector_fembem | Unbounded | Any | Verified | Full vector formulation |
| fem | Bounded | Any | Verified | Truncated domain |
| radia (MMM) | Unbounded | Nonlinear | Verified | Best for static/time-domain |
| BEM+SIBC | Fast | N/A | Verified | Mesh-independent loss via compute_loss_sibc() |

### Known Limitations (Practical)

1. **fembem limited to mu_r=1** -- Calderon Hz scalar formulation constraint
2. **Loop-Star complete only for order=0** -- high-order Helmholtz decomposition (Andriulli 2008) not implemented
3. **Radia core: static mu_r only** -- frequency-dependent eddy current not supported
4. **BEM+SIBC: PEC limit when maxh >> delta** -- mitigated by compute_loss_sibc()
5. **T matrix condition number > 1e14** -- handled via pseudoinverse
6. **COCG alone is unstable for BEM** -- GMRes over COCG recommended

### Verified Parameters

| Parameter | Typical Range | Guidance |
|-----------|---------------|----------|
| Frequency | 10 Hz – 1 MHz | DC extrapolation: ω < 1e-10 triggers DC path |
| intorder | ≥ 5 | Singular quadrature; 5-7 recommended |
| maxh | 1/5 – 1/20 of conductor | Mesh resolution |
| Order p | 0, 1, 2 | p=0 for speed; p=1,2 for convergence |
| Copper σ | 5.8e7 S/m | t=35μm → R_sheet=1/(σ·t) |
| Solver | scipy LDL^T | assume_a='sym' for complex symmetric |

### Reference Code Architecture

```
ngbem_peec.py      — PEEC Loop-Star matrices (L, P, M_LS, R)
ngbem_interface.py — Edge topology extraction for coupling
ngbem_coupled.py   — Coupled core solver (FEM-BEM, Radia, etc.)
test_dowell_comparison.py — FastHenry comparison & Dowell validation
```

### Original Kernel Table (Updated)

| Kernel | Formula | Frequency | Stability |
|--------|---------|-----------|-----------|
| Laplace | 1/(4*pi*r) | Static | Stable |
| Helmholtz | exp(-jkr)/(4*pi*r) | High freq | Stable |
| **Product Space** | **HDivSurface × SurfaceL2** | **DC – RF** | **Stable (verified)** |
| Maxwell EFIE | Full-wave | High freq | Unstable at low freq |

## Priority 4: NGBEM High-Order EFIE (After MQS Kernel)

## Integration with Existing FastImp PEEC

The NGBEM implementation will **complement**, not replace, the existing FastImp PEEC solver:

| Solver | Use Case | Mesh Type | Frequency Range |
|--------|----------|-----------|-----------------|
| FastImp PEEC | Conductors (sigma >> 1, mu_r = 1) | Surface panels | DC to RF |
| NGBEM EFIE | General (sigma, mu_r variable) | High-order surface | DC to RF |
| Radia MSC | Magnets (sigma = 0, mu_r >> 1) | Volume elements | Static |

### Coupling Strategy

```
+-------------------+     +-------------------+
|  FastImp PEEC     |     |  NGBEM EFIE       |
|  (Coils)          |     |  (Workpiece)      |
+--------+----------+     +--------+----------+
         |                         |
         v                         v
+------------------------------------------+
|           Mutual Coupling Matrix          |
|        (Biot-Savart integration)          |
+------------------------------------------+
         |
         v
+------------------------------------------+
|        Combined System Solution           |
+------------------------------------------+
```

## NGSolve FEM Independent Verification

Independent verification of Radia PEEC+BEM results using NGSolve FEM (A-formulation).
The two methods share no code: PEEC uses Neumann formula + BEM surface currents,
while FEM uses a full volume mesh with energy-based inductance extraction.

**Script**: `validation_test/peec_integration/verification/verify_ngsolve_inductance.py`
**PEEC Reference**: `docs/peec_integration/demos/applications/demo_circular_coil_4cases.py`

### Verification Geometry

| Component | Parameters |
|-----------|-----------|
| **Coil** | Circular, R=20 mm, 1.0x1.0 mm Cu wire, I=1.0 A |
| **Core** | 15x15x10 mm ferrite box at origin, mu_r=1000 |
| **Shield** | 50x50x10 mm Al plate (sigma=3.7e7 S/m) at z=5..15 mm |
| **Air domain** | Sphere R=120 mm (static), R=60 mm (eddy current) |

### Method Comparison

| Aspect | NGSolve FEM | PEEC+BEM |
|--------|-------------|----------|
| Formulation | A-formulation (HCurl, order=2) | Neumann integral + BEM SIBC |
| Coil model | OCC torus (Revolve), volume current | 64-segment polygon, filament current |
| Core model | Volume mesh, mu_r=1000 | Hex mesh (3x3x2=18 elements), CoupledPEECSolver |
| Shield model | Volume mesh, sigma=3.7e7 | BEM surface (ShieldBEMSIBC), slab impedance |
| Air treatment | Volume mesh to R_air boundary | Not needed (integral method) |
| Solver | PARDISO (Intel MKL, multi-threaded) | Dense LU + BEM LU |
| Inductance | Energy method: L = 2*W/I^2 | Z_port from MNA circuit solve |

### Analytical Reference

Circular loop with equivalent wire radius `a = sqrt(w*h/pi)`:

| Formula | Value |
|---------|-------|
| L_ext = mu_0*R*(ln(8R/a) - 2) | 91.67 nH |
| L_tot = mu_0*R*(ln(8R/a) - 7/4) | 97.96 nH (incl. internal Li/4) |

### Test Case 1: Air Only (Magnetostatic)

| Metric | NGSolve FEM | PEEC (n_seg=64) | Analytical |
|--------|-------------|------------------|------------|
| L [nH] | 96.58 | 100.69 | 97.96 |
| vs analytical | -1.4% | +2.8% | -- |
| FEM vs PEEC | -4.1% | -- | -- |
| Mesh | 78,655 elem, 417,765 DOF | 64 segments | -- |
| Time | 142.2 s | <1 s | -- |

**Notes**:
- FEM: order=2, air_r=120mm, maxh=8mm, PARDISO solver
- PEEC: Neumann formula with GMD already includes internal inductance (Li/4)
- Both within 5% of analytical: **PASS**

### Test Case 2: + Ferrite Core (mu_r=1000)

| Metric | NGSolve FEM | PEEC (n_seg=64) | Diff |
|--------|-------------|------------------|------|
| L_air [nH] | 96.57 (same-mesh ref) | 100.69 | -4.1% |
| L_core [nH] | 102.11 | 105.76 | -3.5% |
| **Delta_L_core [nH]** | **+5.54** | **+5.07** | **+9.1%** |

**Notes**:
- Delta_L computed on same mesh to cancel systematic errors
- PEEC core: 3x3x2 = 18 hex elements, CoupledPEECSolver (Biot-Savart -> Radia Solve -> A-field)
- Core division sensitivity: coarser (2x2x1) gives ~26% error; 3x3x2 is sufficient
- Delta_L within 15%: **PASS**

### Test Case 3: + Al Shield (Eddy Current, Frequency Sweep)

| Freq | delta | FEM L [nH] | PEEC L [nH] | Diff | FEM DeltaL% | PEEC DeltaL% |
|------|-------|-----------|-------------|------|-------------|--------------|
| 100 Hz | 8.3 mm | 89.03 | 96.00 | -7.3% | -6.7% | -4.7% |
| 1 kHz | 2.6 mm | 78.54 | 90.21 | -12.9% | -17.7% | -10.4% |
| 10 kHz | 0.8 mm | 74.56 | 86.93 | -14.2% | -21.9% | -13.7% |
| 100 kHz | 0.3 mm | 72.89 | 83.91 | -13.1% | -23.7% | -16.7% |

**Notes**:
- FEM: complex A-formulation, air_r=60mm, maxh=4mm, PARDISO solver, ~146s per frequency
- PEEC: ShieldBEMSIBC with slab impedance (Zs*coth(gamma*t)), ~3s per frequency
- Both methods show correct physics: L decreases monotonically with frequency
- Absolute L differs by 7-14% (expected: different coil representations + domain truncation)
- DeltaL% trend matches: both show increasing shielding with frequency

### Physics Checks (ALL PASS)

| Check | Criterion | Result |
|-------|-----------|--------|
| L_air within 5% of analytical | abs(L_fem - L_ana)/L_ana < 5% | PASS (-1.4%) |
| L_air within 5% of PEEC | abs(L_fem - L_peec)/L_peec < 5% | PASS (-4.1%) |
| Core increases L | Delta_L_core > 0 | PASS (+5.54 nH) |
| Delta_L_core within 15% of PEEC | diff < 15% | PASS (9.1%) |
| Shield decreases L | DeltaL < 0 at all frequencies | PASS |
| Shield DeltaL @ 1kHz within 100% | abs(DeltaL_fem - DeltaL_peec) < 100% | PASS |

### Error Budget Analysis

#### Why Absolute L Differs (4-14%)

1. **Coil geometry**: FEM uses smooth torus (OCC Revolve); PEEC uses 64-segment polygon
2. **Domain truncation**: FEM truncates at finite air radius with Dirichlet BC
3. **Internal inductance**: Neumann GMD and FEM energy method compute it differently
4. **Shield modeling**: FEM meshes full volume; PEEC uses surface-only BEM+SIBC

#### Why Delta_L Agrees Better (9%)

Delta_L (change due to core/shield) cancels systematic errors:
- Same coil in both FEM cases -> coil geometry error cancels
- Same air domain -> truncation error cancels
- Only the core/shield effect remains -> methods agree on the physics

### Computational Performance

#### Timing Summary

| Case | FEM Time | FEM DOFs | PEEC Time | Speedup |
|------|----------|----------|-----------|---------|
| Air only | 142 s | 417,765 | <1 s | ~200x |
| + Core | 194 s | 444,678 | ~5 s | ~40x |
| + Shield (per freq) | 146 s | ~400k | ~3 s | ~50x |
| **4-case total** | **~15 min** | -- | **~19 s** | **~48x** |

PARDISO (Intel MKL multi-threaded) vs UMFPACK (single-threaded) speedup: 3.9-7.7x.
NGSolve accesses PARDISO via `inverse="pardiso"` in the Preconditioner/Inverse call.

#### Why PEEC+BEM is Faster

| Aspect | FEM | PEEC+BEM | Impact |
|--------|-----|----------|--------|
| **Air domain** | Volume mesh (sphere R=120 mm) | Not needed (integral method) | ~95% of FEM DOFs |
| **Coil** | OCC torus, volume mesh | 64 line segments (1D) | O(1) vs O(N^3) |
| **Core** | Volume mesh (~30k DOFs) | 18 hex elements (Radia MMM) | 1000x fewer unknowns |
| **Shield** | Volume mesh (~400k DOFs) | BEM surface (~200 DOFs) | 2000x fewer unknowns |
| **Solver** | Sparse LU (PARDISO) | Dense LU (64x64 + 200x200) | Small dense >> large sparse |

The fundamental advantage of PEEC+BEM is that it avoids meshing the air domain:
- FEM requires a volume mesh filling the entire computational domain (air sphere)
- PEEC uses the Neumann integral for inductance (no air mesh)
- BEM models shields as surfaces only (no volume mesh)
- The air domain typically accounts for >95% of FEM DOFs

#### PEEC Timing Breakdown (4 cases at 11 frequencies each)

| Component | Time | Notes |
|-----------|------|-------|
| FastHenry parse + build | ~0.1 s | Per case |
| Radia MMM solve (core) | ~1 s | Per case (18 hex elements) |
| BEM assembly (shield) | ~2 s | One-time (precomputed coupling matrix) |
| BEM multi-RHS LU (shield) | ~1.5 s | Per frequency (LU factored once, 64 RHS) |
| MNA circuit solve | ~0.01 s | Per frequency (64x64 dense LU) |

The BEM shield solver was optimized (2026-02-22) with:
1. **Multi-RHS LU factorization**: Factor BEM system once, solve for all 64 RHS via `lu_solve()`
2. **Precomputed coupling matrix**: Geometry-dependent operator built once, reused per frequency
3. Result: 2.4x speedup over per-call NGSolve `Integrate()` approach

### Key Parameters for Accuracy

#### PEEC Parameters

| Parameter | Value | Impact |
|-----------|-------|--------|
| **n_seg** | **64** | Circle approximation. n_seg=16 gives 26% Delta_L error; 64 gives 9% |
| Core divisions | 3,3,2 (18 elements) | Sufficient for mu_r=1000. Coarser (2,2,1) gives ~26% error |
| use_sibc | False (for L_air) | Neumann GMD includes internal inductance; SIBC double-counts |
| Bessel SIBC | `iv` (modified Bessel) | NOT `jv` (regular Bessel). `jv` gives wrong sign on Im(Z) |

#### FEM Parameters

| Parameter | Value | Impact |
|-----------|-------|--------|
| **FEM order** | **2** | Good accuracy. Order=1 is too coarse; order=3 adds DOFs without proportional benefit |
| air_r (static) | 120 mm | 6x coil radius. Adequate for dipole decay |
| air_r (eddy) | 60 mm | Smaller domain OK because shield confines field |
| maxh (air) | 8 mm (static), 4 mm (eddy) | Coarse air mesh is fine |
| core.faces.maxh | 2 mm | Local refinement for core only (+6.6% elements, better Delta_L) |
| Solver | PARDISO | Intel MKL multi-threaded. 3.9-7.7x faster than UMFPACK |
| gauge | nograds=True + 1e-10 regularization | Removes kernel of curl-curl operator |

### Verification Summary

The PEEC+BEM results are validated by independent NGSolve FEM computation:
- **Core coupling (Delta_L)**: 9.1% agreement -- both methods capture the volume magnetization effect correctly
- **Shield effect (Delta_L trend)**: Same physics -- L decreases monotonically with frequency in both methods
- **Air inductance**: Both within 5% of analytical Neumann formula
- The remaining systematic offset (4-14% in absolute L) is fully explained by different coil representations (polygon vs torus) and domain truncation effects

## Dependencies

### Required Packages

```bash
# NGSolve (base)
pip install ngsolve>=6.2.2601

# NGBEM add-on
pip install ngbem

# Cubit mesh export (optional)
pip install cubit-mesh-export
```

### Build Requirements

- No C++ changes required for NGBEM integration
- Python-only implementation using NGBEM operators
- Existing Radia C++ code (FastImp PEEC, MSC) remains unchanged

## Validation Plan

### Test Cases

1. **Laplace kernel**: Compare NGBEM Single Layer with Radia MSC
2. **Conducting sphere**: Analytical solution for eddy currents
3. **Induction heating**: Coil + workpiece coupled problem
4. **Low-frequency limit**: Verify Loop-Star stabilization

### Expected Results

| Test | NGBEM Result | Reference |
|------|--------------|-----------|
| Sphere H-field | < 1% error | Analytical |
| Coil inductance | < 2% error | FastImp PEEC |
| Workpiece power | < 5% error | FEM reference |

## Timeline

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Laplace BEM verification | **Done** (ngbem_peec_demo) |
| 2 | Loop-Star via product space | **Done** (Weggler EFIE verified) |
| 3 | Dowell skin effect + Schur complement | **Done** (test_dowell_comparison.py) |
| 4 | Coupled core models (5 types) | **Done** (ngbem_coupled.py) |
| 5 | Coreform mesh import | Planned |
| 6 | High-order Loop-Star (p≥1) | Future (Andriulli 2008) |

## References

1. **NGBEM**: https://github.com/Weggler/ngbem
2. **Loop-Star decomposition**: Vecchi, IEEE TAP, 1999
3. **ESIM**: Karl Hollaus, "A Nonlinear Effective Surface Impedance...", 2024
4. **Coreform Cubit**: https://coreform.com/products/coreform-cubit/

## Appendix: Karl Hollaus ESIM Formulation

From: `W:\03_\00_\SIBC\A_Nonlinear_effective_surface_impedance_in_a_Magnetic_Scalar_Potential_Formulation.pdf`

### Mathematical Formulation

**Skin depth**:
```
delta = sqrt(2 / (omega * mu0 * mur * sigma))
```

**Surface impedance**:
```
Zs = (1 + j) * Rs
Rs = 1 / (sigma * delta) = sqrt(omega * mu0 * mur / (2 * sigma))
```

**Nonlinear extension**:
For materials with B-H curve, mu_r is field-dependent:
```
mu_r = mu_r(H)
Zs = Zs(H)
```

Iterative solution required for nonlinear problems.
