# Unified PEEC Loop-Star + MMM + MSC Architecture

**Date**: 2026-01-10
**Status**: Architecture Redesign Phase

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
| **6** | NGBEM low-frequency kernel | Future | Optional high-order elements |

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
|  - Visualization (VTK export)                     |
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
GMSH Format (.msh)  <-- cubit_mesh_export tool
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

# NGSolve can read hex mesh from GMSH
mesh = Mesh("hex_mesh.msh")

# Check element types
for el in mesh.Elements(VOL):
    print(el.type)  # HEXAHEDRON, TET, PRISM, PYRAMID
```

**Challenge**: Netgen's native mesh generator produces tetrahedra only. Hex mesh must be imported.

### Coreform to Netgen Hex Pipeline

#### Option 1: GMSH Intermediate Format

```
Coreform (.cub5) -> GMSH (.msh) -> NGSolve Mesh
```

```python
from coreform_cubit_mesh_export import export_to_gmsh
from ngsolve import Mesh

# Export hex mesh from Cubit to GMSH
export_to_gmsh("model.cub5", "model.msh", element_type="hex")

# Import to NGSolve
mesh = Mesh("model.msh")
```

#### Option 2: Exodus II Format (Native Cubit)

```
Coreform (.cub5) -> Exodus II (.exo) -> NGSolve Mesh
```

NGSolve may support Exodus format via VTK or custom reader.

#### Option 3: Direct Python API

```python
import cubit
from ngsolve import *
from ngsolve.meshes import Make3DMesh

# Get hex elements directly from Cubit
cubit.init([""])
cubit.cmd("import mesh 'model.cub5'")

# Extract hex connectivity
hex_elements = []
for hex_id in cubit.get_hex_conn():
    nodes = cubit.get_connectivity("hex", hex_id)
    hex_elements.append(nodes)

# Get node coordinates
nodes = {}
for node_id in cubit.get_nodeset_nodes(1):
    coords = cubit.get_nodal_coordinates(node_id)
    nodes[node_id] = coords

# Create Netgen mesh directly
# ... (requires netgen.mesh API)
```

### Surface Extraction for BEM

For BEM analysis, only surface mesh is needed:

```python
from ngsolve import *

# Load volume mesh (hex or tet)
mesh = Mesh("model.msh")

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
3. **Create Coreform export tool**: Extend cubit_mesh_export for hex
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
from coreform_cubit_mesh_export import export_to_gmsh

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

### Python API (Implemented)

```python
import radia as rad

# Create conductor (PEEC)
coil = rad.CndSpiral([0,0,0], 0.02, 0.05, 0.005, 3, [0,0,1], 'r', 0.003, 0.002, 5.8e7, 8)
rad.CndSetFrequency(coil, 50000)

# Create magnet (MMM)
magnet = rad.ObjRecMag([0,0,-0.01], [0.1,0.1,0.01], [0,0,0])
mat = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(magnet, mat)

# Create coupled solver
solver = rad.PEECMMMCreate(coil, magnet)
rad.PEECMMMSetFrequency(solver, 50000)
rad.PEECMMMSetVoltage(solver, 1.0, 0.0)  # 1V excitation

# Solve coupled system
result = rad.PEECMMMSolve(solver)
# result = [Z_real, Z_imag, P_cond, P_mag, n_iter]

# Get impedance
Z = rad.PEECMMMImpedance(solver)  # Returns [Z_real, Z_imag]

# Get field
B = rad.PEECMMMFld(solver, [0, 0, 0.05])
# B = [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]

# Frequency sweep
freqs = [1000, 10000, 50000, 100000]
sweep = rad.PEECMMMSweep(solver, freqs)
# sweep = [Z_re1, Z_im1, Z_re2, Z_im2, ...]

# Clean up
rad.PEECMMMDelete(solver)
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

## Priority 3: NGBEM Low-Frequency Kernel Extension

### Objective

**Implement low-frequency/MQS kernel in NGBEM** to enable stable quasi-static analysis.

This is a key development goal that would:
1. Enable high-order elements for induction heating (f < 1 MHz)
2. Avoid numerical issues with standard Maxwell kernel at low frequencies
3. Support quadrilateral meshes from Coreform Cubit

### Current NGBEM Kernels

| Kernel | Formula | Frequency | Stability |
|--------|---------|-----------|-----------|
| Laplace | 1/(4*pi*r) | Static | Stable |
| Helmholtz | exp(-jkr)/(4*pi*r) | High freq | Stable |
| Maxwell | Full-wave EFIE | High freq | **Unstable at low freq** |

### Proposed MQS Kernel

**Magneto-Quasi-Static (MQS) EFIE**:

```
Z_MQS = jw*L + R

where:
  L: Inductance matrix (Neumann kernel)
  R: Resistance matrix (ESIM surface impedance)
```

**Neumann kernel** (low-frequency vector potential):
```
A(r) = mu0/(4*pi) * integral{ J(r') / |r-r'| dV' }

G_Neumann = mu0 / (4*pi*|r-r'|)
```

This is the **same kernel as Laplace** scaled by mu0, so NGBEM's existing Laplace operator can be reused.

### Loop-Star Scaling for MQS

Standard EFIE at low frequency:
```
[ZLL  ZLS] [IL]   [VL]
[ZSL  ZSS] [IS] = [VS]
```

Problem: ZLL ~ O(w), ZSS ~ O(1/w) -> condition number ~ O(1/w^2)

**MQS Loop-Star scaling**:
```
ZLL_mqs = ZLL / (jw)  -> O(1) (pure inductance)
ZSS_mqs = ZSS * (jw)  -> O(1) (pure capacitance)
```

Rescaled system:
```
[L       M_LS  ] [IL ]   [VL/jw ]
[M_SL    1/C   ] [IS'] = [VS*jw ]
```

where IS' = jw*IS (charge derivative = current)

### Implementation Strategy

**Phase 3a: Use existing Laplace kernel for MQS**

```python
from ngbem import SingleLayerPotentialOperator

# Laplace SLP = G(r,r') = 1/(4*pi*|r-r'|)
# MQS inductance = mu0 * Laplace SLP
V_laplace = SingleLayerPotentialOperator(fes, intorder=12, eps=1e-4)

# Scale by mu0 to get inductance operator
mu0 = 4 * np.pi * 1e-7
L_operator = mu0 * V_laplace.mat
```

**Phase 3b: Contribute MQS kernel to NGBEM (upstream)**

If successful, contribute the MQS implementation back to NGBEM project:
- Add `MQSInductanceOperator` class
- Add Loop-Star transformation utilities
- Add ESIM surface impedance support

### Technical Challenges

1. **Loop-Star basis construction**: Need to identify loops and stars from mesh topology
2. **Quadrilateral support**: Verify NGBEM supports quad elements for H(div) space
3. **Singular integrals**: MQS kernel has same singularity as Laplace (manageable)
4. **ACA compression**: Verify ACA works for MQS kernel (should work, same smoothness)

### Validation Plan

| Test Case | Reference | Expected Accuracy |
|-----------|-----------|-------------------|
| Circular loop inductance | Analytical | < 1% |
| Mutual inductance | Neumann formula | < 1% |
| Coil impedance vs frequency | FastImp PEEC | < 5% |
| Induction heating power | RWG-EFIE | < 5% |

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

## Dependencies

### Required Packages

```bash
# NGSolve (base)
pip install ngsolve==6.2.2405

# NGBEM add-on
pip install ngbem

# Coreform mesh export (optional)
pip install coreform-cubit-mesh-export
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
| 1 | Laplace BEM verification | Planned |
| 2 | Loop-Star EFIE implementation | Planned |
| 3 | ESIM integration | Planned |
| 4 | Coreform mesh import | Planned |
| 5 | Validation and benchmarks | Planned |

## References

1. **NGBEM**: https://github.com/Weggler/ngbem
2. **Loop-Star decomposition**: Vecchi, IEEE TAP, 1999
3. **ESIM**: Karl Hollaus, "A Nonlinear Effective Surface Impedance...", 2024
4. **Coreform Cubit**: https://coreform.com/products/coreform-cubit/

## Appendix: Karl Hollaus ESIM Formulation

From: `W:\03_\00_\SIBC\A_Nonlinear_Effective_Surface_Impedance_in_a_Magnetic_Scalar_Potential_Formulation.pdf`

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
