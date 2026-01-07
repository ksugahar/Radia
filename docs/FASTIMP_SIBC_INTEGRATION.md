# FastImp + SIBC Integration Design Document

**Date**: 2026-01-08
**Status**: Phase 4 Python API Complete (CndWire, CndSpiral added)

## Overview

This document describes the design for integrating FastImp-based conductor modeling and Nonlocal SIBC (Surface Impedance Boundary Condition) into Radia for wide-band electromagnetic analysis.

## Goals

1. **Coil/Conductor modeling**: Import FastImp formulation for surface current analysis
2. **Magnetic material**: Use existing ELF/Radia MSC formulation
3. **Conductive magnetic material**: Implement Nonlocal SIBC for materials with both conductivity and permeability

## Target Applications

- **Accelerator magnets**: Insertion devices, undulators, wigglers
- **Kicker magnets**: Fast pulsed magnets for beam injection/extraction (eddy currents in vacuum chamber)
- **WPT (Wireless Power Transfer)**: Coil impedance and resonance analysis
- **On-chip inductors**: High-frequency parasitic extraction
- **Transformers**: Core losses and winding impedance
- **Electromagnetic shielding**: Eddy current shielding effectiveness

## Architecture

```
+-------------------------------------------------------------+
|                    Radia Unified Solver                      |
+-------------------------------------------------------------+
|                                                             |
|  +------------------+  +------------------+  +--------------+|
|  |  Coil/Conductor  |  | Magnetic (s=0)   |  | Magnetic     ||
|  |                  |  |                  |  | Conductor    ||
|  |  FastImp         |  |  ELF/Radia MSC   |  | (s!=0,ur!=1) ||
|  |  Surface K, s    |  |  Surface sm, M   |  | Nonlocal SIBC||
|  +--------+---------+  +--------+---------+  +------+-------+|
|           |                     |                   |        |
|           +----------+----------+-------------------+        |
|                      v                                       |
|              +-------------------+                            |
|              |  Coupled Solver   |                            |
|              +-------------------+                            |
|                      |                                       |
|           +----------+----------+                            |
|           v                     v                            |
|     +-----------+         +-----------+                      |
|     |  HACApK   |         | pFFT(MKL) |                      |
|     | Low freq  |         | High freq |                      |
|     +-----------+         +-----------+                      |
+-------------------------------------------------------------+
```

## Three Modules

### 1. Coil/Conductor (FastImp Formulation)

**Target**: Copper coils, aluminum conductors, wiring
**Properties**: s ~ 10^7 S/m, ur = 1

**Unknowns**:
- K: Surface current density [A/m]
- s: Surface charge density [C/m^2]

**Formulation**: FastImp Full-wave IE
```
A = u * integral{ g(r,r') * K dF' }
Phi = (1/e) * integral{ g(r,r') * s dF' }
g = exp(-jkr) / (4*pi*r)
```

**Acceleration**: pFFT with MKL FFT backend (GPL-free)

### 2. Magnetic Material (ELF/Radia MSC)

**Target**: Permanent magnets, ferrite (high resistivity), soft iron (low frequency)
**Properties**: s ~ 0, ur >> 1

**Unknowns**:
- sm: Magnetic surface charge density [Wb/m^2]
- M: Magnetization vector [A/m]

**Formulation**: MSC (existing Radia)
```
H = -(1/4pi) * integral{ sm * (r-r')/|r-r'|^3 dF' }
sm = M . n_hat
```

**Acceleration**: HACApK (ACA+)

### 3. Conductive Magnetic Material (Nonlocal SIBC)

**Target**: Electrical steel sheets, iron yoke (high frequency)
**Properties**: s ~ 10^6 S/m, ur ~ 1000-10000

**Unknowns**:
- K: Surface current density [A/m]
- Internal: Solved by 2D cross-section FEM

**Formulation**: Nonlocal SIBC (Bilicz et al., 2023)
```
Et = Z{n_hat x K}
Z: Nonlocal impedance operator from 2D FEM

2D cross-section problem:
nabla^2 Eu - jw*u(H)*k*Eu = 0  in Omega (cross-section)
dEu/dn = -jw*u*Hv              on Gamma (boundary)
```

**Reference**:
- Bilicz S, Badics Z, Pavo J. "Wide-band nonlocal impedance boundary condition model for high-conductivity regions in integral equation framework", December 2023.

**Acceleration**: Surface integrals with pFFT/HACApK, 2D FEM solved once per cross-section

## Interaction Matrix

```
[Z_cc  Z_cm  Z_cs] [K_coil ]   [V_ext]
[Z_mc  Z_mm  Z_ms] [sm     ] = [H_ext]
[Z_sc  Z_sm  Z_ss] [K_sibc ]   [E_ext]

Z_cc: Coil-Coil (FastImp)
Z_mm: Magnetic-Magnetic (MSC)
Z_ss: SIBC-SIBC (Nonlocal SIBC)
Z_cm, Z_mc: Coil-Magnetic (cross terms)
Z_cs, Z_sc: Coil-SIBC (cross terms)
Z_ms, Z_sm: Magnetic-SIBC (cross terms)
```

## Implementation Phases

### Phase 1a: FastImp Port (FFTW -> MKL)

**Tasks**:
1. Analyze FastImp source code (https://github.com/ediloren/FastImp)
2. Identify FFTW dependencies in pfft++
3. Replace FFTW calls with Intel MKL FFT (DftiCreateDescriptor, etc.)
4. Build with MSVC + MKL (same as Radia)
5. Validate against original FastImp results

**Key Files**:
- `pfft++/src/` - pFFT implementation
- `fastImp/src/surf/formulation.cc` - Core formulation

**Dependencies**: Intel MKL (already used by Radia)

### Phase 1b: MSC Verification

**Tasks**:
1. Verify existing Radia MSC implementation
2. Ensure compatibility with FastImp surface element format
3. Document interface requirements

### Phase 2: Coil + Magnetic Material Coupling

**Tasks**:
1. Define cross-term computation (Z_cm, Z_mc)
2. Implement coupled solver (iterative or direct)
3. Test with simple coil + iron core problem
4. Benchmark against reference solutions

**Cross-term Physics**:
```
Coil -> Magnetic: B field from coil currents induces magnetization
Magnetic -> Coil: H field from magnetization affects coil impedance
```

### Phase 3: Nonlocal SIBC Implementation

**Tasks**:
1. Implement 2D cross-section FEM solver
   - Option A: Use NGSolve (already integrated)
   - Option B: Simple FEM for Helmholtz equation
2. Implement nonlocal SIBC operator Z{.}
3. Couple with FastImp surface formulation
4. Test with conductive magnetic material
5. Validate against 3D FEM (Ansys HFSS or similar)

**2D FEM Choice**:
- Recommend NGSolve (same build configuration as radia_ngsolve)
- Cross-section problem is simple: Helmholtz equation with Robin BC

### Phase 4: Full Integration

**Tasks**:
1. Unified API design
2. Automatic material type detection
3. Frequency sweep support
4. Resonance finding
5. Documentation and examples

### Phase 5: Transient Analysis (Future Work)

**Goal**: Convert frequency-domain Z(s) to time-domain response

**Candidate Methods**:

1. **Cauer Ladder Network (CLN) Method**
   - Z(s) -> continued fraction expansion -> L-C ladder
   - Physically meaningful equivalent circuit
   - Passive and stable by construction
   - **Requirement: Symmetric impedance matrix**
   - Challenge: May fail for coupled systems with non-symmetric cross terms

2. **Arnoldi-based Model Order Reduction** (Likely choice)
   - Krylov subspace projection
   - PRIMA (Passive Reduced-order Interconnect Macromodeling Algorithm)
   - **Works with non-symmetric matrices**
   - Mathematically robust for coupled multi-physics systems
   - Challenge: Less physical interpretation

3. **Vector Fitting**
   - Rational function approximation of Z(s)
   - Widely used in signal integrity
   - Can be converted to state-space or SPICE model

**Open Questions**:
- Which method works best for electromagnetic systems with eddy currents?
- Where does each method break down?
- How to handle nonlinear magnetic materials in time domain?

**Matrix Symmetry Consideration**:
- Single-physics blocks (Z_cc, Z_mm, Z_ss): Symmetric (reciprocity)
- Cross-coupling blocks (Z_cm, Z_mc, etc.): May be non-symmetric
- For coupled multi-physics systems, **Arnoldi-based methods are likely required**

**Note**: This phase requires further research to determine the optimal approach.
The choice between CLN, Arnoldi, or other methods depends on the specific
characteristics of the impedance function Z(s) obtained from the IE solver.

## API Design (Preliminary)

```python
import radia as rad

# Set units
rad.FldUnits('m')

# Create coil (FastImp)
coil = rad.ObjCoil(vertices, conductivity=5.8e7)

# Create magnetic core (MSC)
core_vertices = [...]
core = rad.ObjHexahedron(core_vertices, [0, 0, 0])
mat_core = rad.MatLin(1000)  # ur = 1000
rad.MatApl(core, mat_core)

# Create conductive magnetic shield (SIBC)
shield = rad.ObjConductiveMagnetic(
    vertices,
    conductivity=1e6,
    permeability=1000,
    cross_section='rectangular'
)

# Assemble and solve
assembly = rad.ObjCnt([coil, core, shield])

# Frequency domain analysis
freq = 1e6  # 1 MHz
Z = rad.SolveImpedance(assembly, freq)

# Frequency sweep
freqs = np.logspace(3, 9, 100)  # 1 kHz to 1 GHz
Z_sweep = rad.ImpedanceSweep(assembly, freqs)

# Find resonances
resonances = rad.FindResonances(assembly, freq_range=[1e6, 1e9])
```

## Build Configuration

All components use MSVC + Intel MKL (same as current Radia):

```
Compiler: MSVC (Visual Studio 2022)
Libraries:
  - Intel MKL (BLAS/LAPACK + FFT)
  - NGSolve (for 2D FEM in SIBC, optional)

No new dependencies:
  - FFTW (GPL) is NOT used - replaced by MKL FFT
  - Intel Compiler is NOT used - only MKL library
```

## Performance Considerations

### pFFT vs HACApK Selection

| Frequency Range | Kernel | Recommended |
|-----------------|--------|-------------|
| DC / Low freq | 1/r | HACApK |
| MQS | 1/r, jw/r | HACApK or pFFT |
| Full-wave | exp(-jkr)/r | pFFT |

### Memory Estimates

For N surface elements:
- Dense matrix: O(N^2)
- HACApK: O(N log N)
- pFFT: O(N)

### Expected Performance (from Bilicz paper)

| Problem | IE + SIBC | 3D FEM |
|---------|-----------|--------|
| Loop (720 elem) | 1 s/freq | 46 s/freq |
| Spiral (1800 elem) | 11 s/freq | 65 s/freq |

## References

### Nonlocal SIBC (Primary Reference)

[1] S. Bilicz, Z. Badics, and J. Pávó, "Wide-band nonlocal impedance boundary condition model for high-conductivity regions in integral equation framework," presented at ISEM 2023 (International Symposium on Electromagnetic Fields in Mechatronics, Electrical and Electronic Engineering), December 2023.
- Affiliation: Budapest University of Technology and Economics, Hungary; Tensor Research, LLC, USA
- ORCID: S. Bilicz (0000-0003-4995-6698), Z. Badics (0000-0001-6176-3675), J. Pávó (0000-0002-9501-7176)
- Funding: Hungarian Scientific Research Fund, Grant K-135307

[2] S. Bilicz, Z. Badics, S. Gyimóthy, and J. Pávó, "A Full-Wave Integral Equation Method Including Accurate Wide-Frequency-Band Wire Models for WPT Coils," IEEE Transactions on Magnetics, vol. 54, no. 3, pp. 1-4, March 2018.
- DOI: 10.1109/TMAG.2017.2771366

### FastImp

[3] Z. Zhu, B. Song, and J. K. White, "Algorithms in FastImp: a fast and wide-band impedance extraction program for complicated 3-D geometries," IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. 24, no. 7, pp. 981-998, July 2005.
- DOI: 10.1109/TCAD.2005.850814
- Affiliation: Massachusetts Institute of Technology
- Abstract: This paper presents algorithms underlying FastImp, an efficient 3-D impedance extraction program using integral equations with pFFT acceleration.

[4] Z. Zhu, B. Song, and J. K. White, "FastImp: A Fast and Wide-Band Impedance Extraction Program for Complicated 3D Geometries," Research Laboratory of Electronics, MIT, 2003.
- Original technical report describing the FastImp formulation and pFFT algorithm

[5] FastImp source code: https://github.com/ediloren/FastImp (MIT License)
- Original implementation by MIT, forked by Enrico Di Lorenzo
- Key algorithms used in Radia: Surface panel discretization, pFFT for O(N log N) matrix-vector products

### pFFT (pre-corrected FFT)

[6] J. R. Phillips and J. K. White, "A Precorrected-FFT Method for Electrostatic Analysis of Complicated 3-D Structures," IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, vol. 16, no. 10, pp. 1059-1072, October 1997.
- DOI: 10.1109/43.662670
- Affiliation: Massachusetts Institute of Technology
- This is the foundational paper for pFFT acceleration used in FastImp

### Related Work

[7] M. Al-Qedra, J. Aronsson, and V. Okhmatovski, "A Novel Skin-Effect Based Surface Impedance Formulation for Broadband Modeling of 3-D Interconnects With Electric Field Integral Equation," IEEE Transactions on Microwave Theory and Techniques, vol. 58, no. 12, pp. 3872-3881, December 2010.

[8] W. C. Gibson, "The Method of Moments in Electromagnetics," Boca Raton: Chapman & Hall/CRC, 2008.

## License Considerations

- FastImp: MIT license (source code)
- FFTW: GPL (NOT used - replaced by MKL)
- Intel MKL: Intel EULA (redistributable, already used by Radia)
- NGSolve: LGPL (dynamic linking OK)

**Result**: No GPL contamination, compatible with Radia's license.

## Implementation Status

### Completed Files (Phase 1-3)

| File | Description | Status |
|------|-------------|--------|
| `src/core/rad_pfft.h` | pFFT acceleration header | Complete |
| `src/core/rad_pfft.cpp` | pFFT with MKL DFTI implementation | Complete |
| `src/core/rad_conductor.h` | Conductor element classes | Complete |
| `src/core/rad_conductor.cpp` | FastImp formulation implementation | Complete |
| `src/core/rad_green_fullwave.h` | Green's function header | Complete |
| `src/core/rad_green_fullwave.cpp` | DC/MQS/Full-wave Green's functions | Complete |
| `src/core/rad_coupled_solver.h` | Coupled solver header | Complete |
| `src/core/rad_coupled_solver.cpp` | Cross-term computation (Z_cm, Z_mc) | Complete |
| `src/core/rad_sibc.h` | Nonlocal SIBC header | Complete |
| `src/core/rad_sibc.cpp` | 2D FEM + SIBC implementation | Complete |

### Key Classes Implemented

1. **radTConductor**: Conductor element with Radia-compatible geometry input
   - CreateFromRecBlock(), CreateFromHexahedron()
   - CreateWire(), CreateLoop(), CreateSpiral()
   - Surface panel discretization

2. **radTConductorSolver**: FastImp-based IE solver
   - DC/MQS/EMQS/FullWave formulation selection
   - EFIE + charge continuity equation
   - pFFT acceleration for large problems

3. **radTPfft**: Precorrected FFT acceleration
   - MKL DFTI backend (GPL-free)
   - Toeplitz->Circulant embedding
   - Near-field correction

4. **radTGreenFunction**: Scalar/Vector Green's functions
   - DC: 1/(4πr)
   - Full-wave: exp(-jkr)/(4πr)
   - Dyadic Green's function for EFIE

5. **radTCrossTerms**: Conductor-Magnetic coupling
   - Z_cm: Conductor current → Magnetic H field
   - Z_mc: Magnetization → Conductor E field

6. **radTCoupledSolver**: Unified coupled solver
   - Direct block LU factorization
   - Iterative coupling option
   - Automatic subsystem detection

7. **radTCrossSection2DFEM**: 2D FEM for SIBC
   - Rectangle/Circle mesh generation
   - Helmholtz equation: ∇²E - jωμσE = 0
   - LU factorization for efficient multi-RHS solve

8. **radTNonlocalSIBC**: Nonlocal impedance operator
   - Z{.} operator from 2D FEM
   - Local vs nonlocal comparison
   - Skin depth calculation

### Remaining Work (Phase 4+)

1. **Python API binding** (radpy_pyapi.cpp) ✓ (completed 2026-01-08)
   - `CndRecBlock`, `CndHexahedron`, `CndWire`, `CndLoop`, `CndSpiral` - conductor creation
   - `CndSetFrequency`, `CndSolve`, `CndGetImpedance`, `CndImpedanceSweep` - analysis
   - `CndFld`, `CndNumPanels` - field computation and info
   - `MatSIBC` - SIBC material definition
2. **CMake integration** for new source files ✓ (completed 2026-01-08)
3. **Test cases and validation** ✓ (test_fastimp_conductor.py, test_fastimp_core.py)
4. **Documentation and examples**
5. **Phase 5: Transient analysis** (CLN/Arnoldi methods)

## Python API Design (Phase 4)

### Overview

The Python API follows existing Radia conventions:
- Functions return object handles (integers)
- Handles are passed to subsequent functions
- Error handling via RuntimeError exceptions

### Proposed API Functions

#### Conductor Creation

```python
# Create conductor from rectangular block
# Similar to rad.ObjRecMag() but for conductor analysis
cond = rad.CndRecBlock(center, dimensions, conductivity)
# Parameters:
#   center: [x, y, z] center coordinates (in current units)
#   dimensions: [Lx, Ly, Lz] block dimensions
#   conductivity: electrical conductivity [S/m]
# Returns: conductor handle (int)

# Create conductor from hexahedron vertices
# Similar to rad.ObjHexahedron()
cond = rad.CndHexahedron(vertices, conductivity)
# Parameters:
#   vertices: list of 8 vertex coordinates [[x1,y1,z1], ...]
#   conductivity: electrical conductivity [S/m]

# Create conductor from existing Radia magnetic object
# Converts magnetic object to conductor (for coupling)
cond = rad.CndFromObj(mag_obj, conductivity)

# Create wire conductor along a path
cond = rad.CndWire(path, cross_section, width, height=0, conductivity=5.8e7)
# Parameters:
#   path: list of points defining wire center line
#   cross_section: "circular" or "rectangular"
#   width: wire width (or diameter for circular)
#   height: wire height (ignored for circular)
#   conductivity: default copper (5.8e7 S/m)

# Create circular loop coil
cond = rad.CndLoop(center, radius, normal, cross_section, wire_width, wire_height=0, conductivity=5.8e7)

# Create spiral coil
cond = rad.CndSpiral(center, inner_radius, outer_radius, pitch, num_turns, axis,
                     cross_section, wire_width, wire_height=0, conductivity=5.8e7)
```

#### Conductor Container

```python
# Create conductor container (like rad.ObjCnt for magnets)
cond_cnt = rad.CndCnt([cond1, cond2, ...])

# Add conductor to container
rad.CndCntAdd(cond_cnt, cond3)
```

#### Analysis Configuration

```python
# Set analysis formulation
rad.CndSetFormulation(cond, formulation)
# formulation: "dc", "mqs", "emqs", "fullwave"

# Set analysis frequency
rad.CndSetFrequency(cond, frequency)
# frequency: analysis frequency in Hz (0 for DC)

# Set surface panel discretization
rad.CndSetPanelDensity(cond, num_panels_per_face)

# Enable/disable pFFT acceleration
rad.CndSetPfft(cond, enable=True)
```

#### Port Definition (for Impedance Extraction)

```python
# Define port between two terminals
rad.CndDefinePort(cond, terminal1_panels, terminal2_panels)
# terminal1_panels, terminal2_panels: list of panel indices or "auto"

# For simple wire/loop, auto-detect terminals
rad.CndDefinePortAuto(cond)
```

#### Solver

```python
# Solve at single frequency (conductor only)
rad.CndSolve(cond)

# Solve coupled system (conductor + magnetic)
rad.CoupledSolve(cond_cnt, mag_cnt, precision=1e-4, max_iter=1000)
# Solves the full coupled system:
#   [Z_c   Z_cm] [J]   [V]
#   [Z_mc  Z_m ] [M] = [H_ext]

# Get impedance after solve
Z = rad.CndGetImpedance(cond)
# Returns: complex impedance [Ohm]

# Frequency sweep
freqs = [1e3, 10e3, 100e3, 1e6]
Z_list = rad.CndImpedanceSweep(cond, freqs)
# Returns: list of complex impedances
```

#### Field Computation

```python
# Compute B field from conductor currents
B = rad.CndFld(cond, 'b', point)
# Returns: [Bx, By, Bz] complex for AC, real for DC

# Compute E field from conductor
E = rad.CndFld(cond, 'e', point)

# Batch field computation
points = [[x1,y1,z1], [x2,y2,z2], ...]
B_list = rad.CndFldBatch(cond, 'b', points)
```

#### Solution Access

```python
# Get surface current density K [A/m]
K = rad.CndGetSurfaceCurrent(cond)
# Returns: list of complex 3-vectors for each panel

# Get surface charge density sigma [C/m^2]
sigma = rad.CndGetSurfaceCharge(cond)
# Returns: list of complex values for each panel

# Get panel information
panels = rad.CndGetPanels(cond)
# Returns: list of dicts with 'center', 'normal', 'area', 'vertices'
```

#### SIBC Functions (for Conductive Magnetic Materials)

```python
# Create conductive magnetic material with SIBC
sibc_mat = rad.MatSIBC(conductivity, mu_r)
# Parameters:
#   conductivity: electrical conductivity [S/m]
#   mu_r: relative permeability

# Apply SIBC material to object
rad.MatApl(hex_obj, sibc_mat)

# Set SIBC type
rad.SIBCSetType(sibc_mat, sibc_type)
# sibc_type: "local" or "nonlocal"

# Set cross-section mesh for nonlocal SIBC
rad.SIBCSetCrossSection(sibc_mat, shape, params)
# shape: "rectangle" or "circle"
# params: [width, height] for rectangle, [radius] for circle
```

### Implementation Plan

#### File Structure

```
src/radia/
├── radpy_pyapi.cpp          # Add new functions here
└── ...

src/lib/
├── radentry.h               # Add C API declarations
└── radentry.cpp             # Add C API implementations

src/core/
├── rad_conductor.h          # Already implemented
├── rad_conductor.cpp        # Already implemented
├── rad_sibc.h               # Already implemented
└── rad_sibc.cpp             # Already implemented
```

#### Implementation Steps

1. **Add C API to radentry.h/cpp** (wrapper functions)
   ```cpp
   EXP int CALL RadCndRecBlock(int* handle, double* center, double* dims, double sigma);
   EXP int CALL RadCndSolve(int handle);
   EXP int CALL RadCndFld(double* field, int handle, char fieldType, double* point);
   // etc.
   ```

2. **Add Python bindings to radpy_pyapi.cpp**
   ```cpp
   static PyObject* radia_CndRecBlock(PyObject* self, PyObject* args);
   static PyObject* radia_CndSolve(PyObject* self, PyObject* args);
   static PyObject* radia_CndFld(PyObject* self, PyObject* args);
   // etc.
   ```

3. **Register in module method table**
   ```cpp
   static PyMethodDef radia_methods[] = {
       // ... existing methods ...
       {"CndRecBlock", radia_CndRecBlock, METH_VARARGS, "Create conductor from rectangular block"},
       {"CndSolve", radia_CndSolve, METH_VARARGS, "Solve conductor system"},
       {"CndFld", radia_CndFld, METH_VARARGS, "Compute field from conductor"},
       // etc.
   };
   ```

### Example Usage

```python
import radia as rad
import numpy as np

rad.FldUnits('m')

# ========== Example 1: Simple wire loop impedance ==========

# Create circular loop coil (10cm radius, 1mm wire)
loop = rad.CndLoop(
    center=[0, 0, 0],
    radius=0.1,
    normal=[0, 0, 1],
    cross_section='circular',
    wire_width=1e-3,  # 1mm diameter wire
    conductivity=5.8e7  # Copper
)

# Frequency sweep
freqs = np.logspace(3, 7, 50)  # 1kHz to 10MHz
Z = rad.CndImpedanceSweep(loop, freqs.tolist())

# Extract L and R
R = [z.real for z in Z]
L = [z.imag / (2 * np.pi * f) for z, f in zip(Z, freqs)]

# ========== Example 2: Coil with magnetic core ==========

# Create magnetic core (soft iron cube)
core_vertices = [
    [-0.05, -0.05, -0.1], [0.05, -0.05, -0.1],
    [0.05, 0.05, -0.1], [-0.05, 0.05, -0.1],
    [-0.05, -0.05, 0.1], [0.05, -0.05, 0.1],
    [0.05, 0.05, 0.1], [-0.05, 0.05, 0.1]
]
core = rad.ObjHexahedron(core_vertices, [0, 0, 0])

# For low-frequency analysis, use standard magnetic material
mat_low_freq = rad.MatLin(4000)  # mu_r = 4000
rad.MatApl(core, mat_low_freq)

# For high-frequency analysis with eddy currents, use SIBC
# mat_high_freq = rad.MatSIBC(2e6, 4000)  # σ=2MS/m, μr=4000
# rad.MatApl(core, mat_high_freq)

# Create coil around core
coil = rad.CndSpiral(
    center=[0, 0, 0],
    inner_radius=0.06,
    outer_radius=0.08,
    pitch=0.01,
    num_turns=20,
    axis=[0, 0, 1],
    cross_section='rectangular',
    wire_width=0.005,
    wire_height=0.002
)

# Create containers
mag_cnt = rad.ObjCnt([core])
cond_cnt = rad.CndCnt([coil])

# Solve coupled system
rad.CoupledSolve(cond_cnt, mag_cnt, precision=1e-4, max_iter=1000)

# Get impedance
Z = rad.CndGetImpedance(coil)
print(f"Impedance at DC: {Z}")

# ========== Example 3: Field computation ==========

# Compute B field at observation points
obs_points = [[0, 0, z] for z in np.linspace(-0.2, 0.2, 41)]
B_list = rad.CndFldBatch(coil, 'b', obs_points)

# For AC analysis, B is complex
B_magnitude = [np.abs(np.sqrt(b[0]**2 + b[1]**2 + b[2]**2)) for b in B_list]
```

### Naming Convention Rationale

| Prefix | Meaning | Examples |
|--------|---------|----------|
| `Cnd` | Conductor operations | `CndRecBlock`, `CndSolve`, `CndFld` |
| `SIBC` | Surface impedance BC | `SIBCSetType`, `SIBCSetCrossSection` |
| `Coupled` | Coupled analysis | `CoupledSolve` |
| `Mat` | Material (existing) | `MatSIBC` (new) |

This follows Radia conventions:
- `Obj` for magnetic objects → `Cnd` for conductors
- `Fld` for field computation → `CndFld` for conductor field
- `Solve` for solver → `CndSolve` for conductor solver

## Coil on Magnetic Core: Frequency-Dependent Characteristics

This section describes the physics of coils wound on magnetic cores with various materials,
demonstrating when each solver module (MSC, FastImp, SIBC) is appropriate.

### Physical Model

```
┌─────────────────────────────────────────────┐
│           Coil wound on magnetic core        │
│                                             │
│    ╭───╮  ╭───╮  ╭───╮                     │
│   ╭┤   ├──┤   ├──┤   ├╮  ← Coil winding    │
│   │╰───╯  ╰───╯  ╰───╯│    (copper, σ=5.8e7)│
│   │ ┌─────────────┐   │                     │
│   │ │  Magnetic   │   │  μr >> 1           │
│   │ │    Core     │   │  σ > 0 (conductive)│
│   │ │   (σ, μr)   │   │                     │
│   │ └─────────────┘   │                     │
│   │╭───╮  ╭───╮  ╭───╮│                     │
│   ╰┤   ├──┤   ├──┤   ├╯                     │
│    ╰───╯  ╰───╯  ╰───╯                      │
└─────────────────────────────────────────────┘
```

### Material Properties Comparison

| Core Material | μr | σ [S/m] | L_DC (100 turns, 10cm path, 1cm^2) | Q @ 1kHz |
|--------------|-----|---------|-----------------------------------|----------|
| Air (no core) | 1 | 0 | 0.013 mH | 0.3 |
| Ferrite (MnZn) | 2000 | 0.1 | 25.1 mH | 337 |
| Ferrite (NiZn) | 200 | 1e-4 | 2.5 mH | 48 |
| Silicon Steel | 4000 | 2e6 | 50.3 mH | ~0 |
| Pure Iron | 5000 | 1e7 | 62.8 mH | ~0 |

**Key Insight**: High permeability does NOT guarantee high Q-factor. Conductive cores
(silicon steel, iron) have severe eddy current losses that reduce Q to near zero at AC frequencies.

### Frequency-Dependent Phenomena

| Frequency Range | Dominant Effect | Model Required |
|-----------------|-----------------|----------------|
| DC | Magnetic circuit (flux concentration) | Radia MSC |
| Low freq (< 1 kHz) | Eddy current loss begins | MQS |
| Mid freq (1 kHz - 1 MHz) | Skin effect significant | MQS + SIBC |
| High freq (> 1 MHz) | Surface currents dominate | Full-wave / FastImp |

### Skin Depth and SIBC Selection Guide

For a 10mm diameter core:

| Material | 50 Hz | 1 kHz | 10 kHz | 100 kHz | 1 MHz |
|----------|-------|-------|--------|---------|-------|
| Ferrite (MnZn) | DC | DC | DC | DC | Nonlocal |
| Silicon Steel | Surf | Surf | Surf | Surf | Surf |
| Pure Iron | Surf | Surf | Surf | Surf | Surf |
| Copper | Local | Local | Surf | Surf | Surf |

**Legend**:
- **DC**: δ >> d → Full penetration, quasi-static → Use **Radia MSC**
- **Nonlocal**: δ ~ d → Internal distribution matters → Use **Nonlocal SIBC**
- **Local**: 0.1d < δ < d → Thin skin approximation OK → Use **Local SIBC** (Zs = (1+j)/(σδ))
- **Surf**: δ << d → Surface currents only → Use **FastImp**

### Solver Selection Algorithm

```python
def select_solver(material_sigma, material_mu_r, frequency, dimension):
    """
    Select appropriate solver based on skin depth vs characteristic dimension.

    Parameters:
        material_sigma: Conductivity [S/m]
        material_mu_r: Relative permeability
        frequency: Operating frequency [Hz]
        dimension: Characteristic dimension [m]

    Returns:
        Recommended solver module
    """
    # Calculate skin depth
    if frequency <= 0 or material_sigma <= 0:
        delta = float('inf')
    else:
        omega = 2 * pi * frequency
        mu = MU_0 * material_mu_r
        delta = sqrt(2 / (omega * mu * material_sigma))

    ratio = delta / dimension

    if ratio > 10:
        return "Radia MSC (quasi-static)"
    elif ratio > 1:
        return "Nonlocal SIBC (2D FEM cross-section)"
    elif ratio > 0.1:
        return "Local SIBC (Zs = (1+j)/(sigma*delta))"
    else:
        return "FastImp (surface current only)"
```

### Practical Application Guidelines

| Application | Recommended Core | Frequency Range | Solver |
|-------------|------------------|-----------------|--------|
| Power transformer | Laminated Si steel | 50-60 Hz | Radia MSC + loss factor |
| Choke coil | MnZn ferrite | 100 Hz - 1 MHz | Nonlocal SIBC |
| RF inductor | NiZn ferrite | 1 MHz - 100 MHz | Local SIBC |
| Air-core coil | None | All frequencies | FastImp |
| Eddy current probe | None (air) | 100 kHz - 10 MHz | FastImp |

### Impedance Formulas

**Inductance (frequency-dependent)**:
```
L(f) = L_ext + L_int(f)

L_ext = μ0 * μr_eff * N^2 * A / l   (external inductance)
L_int(f) ∝ 1/√f                      (internal inductance, decreases with skin effect)
```

**Resistance (frequency-dependent)**:
```
R(f) = R_DC + R_eddy(f) + R_hyst(f)

R_DC = ρ * l / A_wire               (wire DC resistance)
R_eddy ∝ f^2                        (eddy current loss in core)
R_hyst ∝ f                          (hysteresis loss in core)
```

**Quality Factor**:
```
Q(f) = ω * L(f) / R(f)
```

### Example Analysis Script

See `examples/fastimp_integration/coil_on_magnetic_core_analysis.py` for a complete
analysis example that generates:
- Inductance vs frequency plots
- Resistance vs frequency plots
- Q-factor vs frequency plots
- Skin depth analysis for different materials

### Key Takeaways

1. **Ferrite cores are optimal for AC applications** due to low conductivity (high Q)
2. **Iron/steel cores require lamination** to reduce eddy current losses at power frequencies
3. **Nonlocal SIBC is essential** when skin depth is comparable to core dimension
4. **Local SIBC is sufficient** when skin depth is small but not negligible
5. **FastImp (surface current)** is appropriate when skin depth << dimension
