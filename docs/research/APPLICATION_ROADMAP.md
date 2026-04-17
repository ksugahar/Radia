# Radia Application Roadmap

Research laboratory applications and collaboration roadmap.

## Target Applications

### 1. Wireless Power Transfer (WPT)
**Status**: Implemented

| Feature | Status | Module |
|---------|--------|--------|
| Coil inductance/resistance | OK | `analysis.py` |
| Mutual inductance | OK | `lanczos_reduction.py` |
| Coupling coefficient | OK | `demo_wpt.py` |
| S-parameter analysis | OK | `demo_wpt.py` |
| Skin effect (Dowell) | OK | `veriloga_generator.py` |
| SPICE export | OK | Verilog-A |
| Ferrite core coupling | OK | CplMag |

**Examples**: `examples/peec_integration/demo_wpt.py`

---

### 2. Induction Heating
**Status**: Implemented

| Feature | Status | Module |
|---------|--------|--------|
| Surface impedance Zs | OK | `esim_cell_problem.py` |
| Skin depth calculation | OK | `demo_induction_heating.py` |
| Power density | OK | `demo_induction_heating.py` |
| Nonlinear mu(H) - ESIM | OK | `esim_coupled_solver.py` |
| Complex permeability | OK | `esim_cell_problem.py` |
| Temperature rise estimate | OK | `demo_induction_heating.py` |

**Examples**: `examples/peec_integration/demo_induction_heating.py`

**Collaboration**: Karl Hollaus (TU Wien) - ESIM method

---

### 3. Accelerator / Beam Tracking
**Status**: Partial

| Feature | Status | Module |
|---------|--------|--------|
| Static field maps | OK | `rad.Fld()` |
| Quadrupole magnets | OK | Background field API |
| Dipole magnets | OK | MMM solver |
| Sextupole/higher | TODO | Need implementation |
| Symplectic integrator | TODO | Need implementation |
| Space charge | TODO | Future (with PEEC?) |


**GPU Acceleration (JHPCN application submitted)**:

| Resource | GPU | Memory | Target Use |
|----------|-----|--------|------------|
| MDX I | A100 x2 | 80GB x2 | Large-scale tracking |
| MDX II | H200 | 141GB | Massive particle count |

GPU acceleration targets:
- **Particle tracking**: Parallel integration of 10^6+ particles
- **Field interpolation**: GPU texture memory for field maps
- **Space charge**: N-body on GPU

**TODO**:
- [ ] Symplectic integrator (Boris, Yoshida)
- [ ] Transfer matrix extraction
- [ ] Twiss parameter computation
- [ ] MAD-X/Elegant format export
- [ ] GPU particle pusher (CUDA/OpenCL)
- [ ] GPU field map interpolation

---

### 4. NMR Magnet Design
**Status**: Partial

| Feature | Status | Module |
|---------|--------|--------|
| High homogeneity field | OK | MMM solver |
| Magnet array design | OK | `smco_array.py` |
| Shim coil optimization | Partial | `demo_nmr_magnet.py` |
| Spherical harmonic analysis | Partial | `demo_nmr_magnet.py` |
| Gradient coil design | TODO | Need implementation |
| Cryostat modeling | OK | NGSolve FEM |
| Eddy current in shields | TODO | PEEC + surface Z |

**Examples**:
- `examples/smco_magnet_array/smco_array.py` - SmCo hexagonal magnet array
- `examples/peec_integration/demo_nmr_magnet.py` - Spherical harmonic shimming

**TODO**:
- [x] Spherical harmonic expansion of B field
- [x] Field homogeneity metrics (ppm)
- [ ] Shim coil current optimization (physical coils)
- [ ] Gradient coil Biot-Savart
- [ ] Halbach ring permanent magnet NMR

---

### 5. Magnetic Levitation (MagLev)
**Status**: In Development (student project)

| Feature | Status | Module |
|---------|--------|--------|
| Permanent magnet force | OK | `rad.FldForce()` |
| Soft iron interaction | OK | MMM solver |
| Halbach array | OK | ObjRecMag |
| Eddy current damping | TODO | PEEC + ESIM |
| Stability analysis | TODO | Need implementation |
| Control system interface | TODO | Future |

**Examples**: Student project in progress

**TODO**:
- [ ] Force/torque calculation demos (student)
- [ ] Halbach array optimization
- [ ] Eddy current braking (moving conductor)
- [ ] Stability margin computation
- [ ] Earnshaw stability analysis

---

## Collaborations

### Karl Hollaus (TU Wien)
**Research**: ESIM (Effective Surface Impedance Method)

**ESIM Theory**:

The 1D Cell Problem for nonlinear surface impedance:
```
rho * d^2H/ds^2 + j*omega*mu(|H|)*H = 0    for s in [0, infinity)

Boundary conditions:
    H(0) = H0           (surface tangential field)
    H(infinity) = 0     (field vanishes at infinity)

Surface impedance:
    Z(H0) = 2*(P' + j*Q') / |H0|^2
```

**Integration Status**:
| Item | Status | Module | Notes |
|------|--------|--------|-------|
| ESIM cell problem | OK | `esim_cell_problem.py` | 1D BVP solver |
| Complex mu(H) | OK | `BHCurveInterpolator` | Nonlinear permeability |
| NGSolve 2D FEM | OK | NGSolve integration | Training data for PyKAN |
| Homogenization | OK | Cell problem | Effective Zs extraction |
| Complex permeability | OK | `esim_cell_problem.py` | mu' - j*mu" support |

**Reference**: K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics, 2025, DOI: 10.1109/TMAG.2025.3613932

**PyKAN Integration**:

PyKAN can learn the mapping H0 -> Z(H0) from ESIM cell problem solutions:

```
Training: {H0_i, Z(H0_i)} from ESIM cell problem
PyKAN learns: Z_KAN(H0) ≈ Z(H0)
Output: Symbolic formula or Verilog-A for SPICE
```

**Future Work**:
- [ ] 3D ESIM extension (anisotropic surface impedance tensor)
- [ ] Temperature-dependent mu(T, H) with thermal coupling
- [ ] Joint paper on PyKAN + ESIM for frequency-dependent materials
- [ ] Multi-physics: ESIM + thermal + mechanical stress

---

### Hane Lab (Toyo University)
**Research**: RNA (Reluctance Network Analysis) with NGSolve

**Potential Collaboration**:
| Area | Radia Contribution | Hane Lab Contribution |
|------|-------------------|----------------------|
| FEM mesh | NGSolve integration | RNA expertise |
| Circuit extraction | SPICE export | Network topology |
| Model reduction | PRIMA/Lanczos | RNA reduction |
| Magnetic cores | MMM solver | Reluctance modeling |

**Integration Ideas**:
- RNA network from Radia MMM elements
- Hybrid RNA-PEEC for transformers
- NGSolve mesh to RNA topology conversion
- Comparison: RNA vs MMM accuracy

**Research Hypothesis: MMM H-matrix = RNA**:

The H-matrix representation of the MMM interaction matrix may be mathematically equivalent to the RNA (Reluctance Network Analysis) formulation. This is a key research direction to explore:

| Concept | MMM with H-matrix | RNA |
|---------|-------------------|-----|
| Matrix structure | Hierarchical low-rank | Sparse network |
| Physical meaning | Dipole-dipole interaction | Magnetic reluctance path |
| Compression | ACA+ approximation | Network topology |
| Sparsity pattern | Distance-based admissibility | Circuit connectivity |

If this equivalence holds, it would provide:
1. **Physical interpretation** of H-matrix blocks as reluctance paths
2. **Automatic network extraction** from MMM solver
3. **Unified framework** for BEM and circuit-based magnetic analysis

**TODO**:
- [ ] Contact and discuss collaboration
- [ ] Identify common use cases
- [ ] Prototype RNA export from Radia
- [ ] Investigate MMM H-matrix ≈ RNA equivalence

---

## Analytical Formula Library

### Goal: Fast & Accurate EM Sources for NGSolve

Incorporate well-known analytical solutions as Radia field sources. These provide:
- **Speed**: O(1) evaluation vs O(N) numerical integration
- **Accuracy**: Exact solutions (no mesh dependency)
- **NGSolve integration**: CoefficientFunction for FEM coupling

### Priority Analytical Formulas

| Formula | Geometry | Status | Module |
|---------|----------|--------|--------|
| Uniformly magnetized sphere | Sphere | TODO | `analytical_sources.py` |
| Uniformly magnetized cylinder | Cylinder | TODO | `analytical_sources.py` |
| Halbach array (2D) | Rectangular | Partial | `ObjRecMag` |
| Solenoid (thin wire) | Coil | OK | Arc current |
| Helmholtz coil pair | Coil pair | TODO | `analytical_sources.py` |
| Dipole field | Point | OK | `rad.Fld()` |

### CLN (Continued Fraction / Ladder Network) Representations

**Background**: MoM-CLN approach converts frequency-dependent integral equations into SPICE-compatible ladder networks. This was explored in:

> Compumag 2019 Paris submission (#267): "MoM-CLN for electromagnetic analysis"

Analytical solutions can be expressed in CLN form for SPICE simulation:

```
Frequency-dependent EM problem
    |
    v
Method of Moments (MoM) formulation
    |
    v
Lanczos tridiagonalization (PRIMA)
    |
    v
[R1]--[L1]--[R2]--[L2]--...--[Rn]--[Ln]
    |
    v
SPICE netlist / Verilog-A
```

**Key insight**: CLN is mathematically equivalent to PRIMA Lanczos reduction. The tridiagonal Lanczos matrix directly maps to a series RL ladder network.

**Target formulas for CLN**:
- [ ] Cylindrical conductor skin effect (Bessel functions -> CLN)
- [ ] Spherical harmonic field expansion (NMR shimming)
- [ ] Eddy current in conducting sphere
- [ ] Mutual inductance between coaxial circles (elliptic integrals)
- [ ] MoM interaction matrix -> CLN via PRIMA

### Implementation Approach

```python
# Example: Analytical cylinder field as NGSolve CoefficientFunction
from radia import AnalyticalCylinderField

# Create magnetized cylinder source (analytical)
cyl_field = AnalyticalCylinderField(
    center=[0, 0, 0],
    axis=[0, 0, 1],
    radius=0.01,      # 10 mm
    length=0.05,      # 50 mm
    magnetization=1e6  # A/m
)

# Use in NGSolve FEM
B_cf = cyl_field.as_coefficient_function('b')
```

---

## Development Priorities

### Phase 1: Core Completion (Current)
1. ~~WPT examples~~ DONE
2. ~~Induction heating examples~~ DONE
3. ~~Verilog-A export~~ DONE
4. PAMELA integration for PEEC
5. Fix NotImplementedError stubs
6. Analytical formula library (sphere, cylinder)

### Phase 2: Application Examples
1. Beam tracking demos
2. NMR magnet examples
3. MagLev force calculation demos
4. Transformer with ferrite core

### Phase 3: Advanced Features
1. Spherical harmonic analysis (NMR)
2. Symplectic tracking (accelerator)
3. Moving conductor eddy currents (MagLev)
4. 3D ESIM (induction heating)

### Phase 4: Collaborations
1. Karl: 3D ESIM, temperature effects
2. Hane: RNA-Radia integration
3. Community: Additional applications

---

## Module Status Summary

| Module | WPT | IH | Accel | NMR | MagLev |
|--------|-----|-----|-------|-----|--------|
| MMM solver | - | OK | OK | OK | OK |
| PEEC | OK | - | - | - | - |
| CplMag | OK | OK | - | - | - |
| ESIM | - | OK | - | - | TODO |
| Beam tracking | - | - | Partial | - | - |
| SPICE export | OK | OK | - | - | - |
| NGSolve FEM | OK | OK | OK | OK | OK |

Legend: OK = Ready, Partial = Basic implementation, TODO = Not started, - = Not applicable

---

## External Dependencies

| Library | Purpose | Status |
|---------|---------|--------|
| **NGSolve** | FEM solver | Integrated |
| **HACApK** | H-matrix | Integrated |
| **Intel MKL** | BLAS/LAPACK | Integrated |
| **PAMELA** | PEEC | TODO |
| **PyKAN** | Material learning | Partial |

---

**Last Updated**: 2026-01-16
