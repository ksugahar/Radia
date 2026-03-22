# Radia - Electromagnetic Simulation Framework for Magnetic Levitation Systems

[![CI](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml)
[![Policy Lint](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-BSD%20%2B%20MIT-green.svg)](LICENSE)

**A Python-native electromagnetic source framework for [NGSolve](https://ngsolve.org/) and [ngsolve.bem](https://docu.ngsolve.org/latest/i-tutorials/unit-8.1-ngbem/ngbem.html).**

Radia provides **analytical field sources** ($B$, $H$, $A$, $\Phi$) as native C++ `CoefficientFunction` objects that NGSolve and ngsolve.bem consume directly during finite element assembly. Where FEM needs boundary conditions and source terms, Radia delivers them — analytically, without meshing the source region.

## 🚀 Mission: The Design Tool for Open-Space Magnetics

**Radia** is a specialized simulation framework developed as a **Design Tool** targeting:

*   **Magnetic Levitation (MagLev)**
*   **Wireless Power Transfer (WPT)**
*   **Induction Heating**
*   **Particle Accelerators & Beamlines**

Unlike general-purpose FEM tools optimized for motors (rotating machinery) with narrow gaps and sliding meshes, Radia addresses the unique challenges of **Open-Space Magnetics**:

*   **Large Air Gaps**: Solves open boundary problems exactly without meshing the air.
*   **Moving Permanent Magnets**: Dynamic simulation of moving magnets (levitation, undulators) is trivial and noise-free because there is no air mesh to distort or regenerate.
*   **Complex Source Geometries**: Models race-track coils, helical undulators, and Halbach arrays analytically with perfect geometric fidelity.
*   **System Level Simulation**: Designed for systems where the field source topology (coils/magnets) defines the performance.

**This is not just a solver; it is a Framework.** We provide the architecture to build specific solvers for your unique magnetic systems.

### Current Capabilities & Active Development

**Implemented**:
*   **PEEC Circuit Extraction**: C++ MNA solver with MKL LAPACK for L, R, C, M extraction from conductor geometry.
*   **Multi-filament Skin Effect**: `nwinc`/`nhinc` cross-section subdivision for skin/proximity effect modeling.
*   **FastHenry Compatibility**: Parse `.inp` files directly with one-step solve.
*   **Coupled PEEC+MMM**: Conductor-core coupling via Biot-Savart + Radia magnetostatic solver.
*   **ESIM (Effective Surface Impedance Method)**: Nonlinear surface impedance for induction heating analysis.
*   **Templated BiCGSTAB**: Shared between MSC (real) and PEEC (complex) solvers via `rad_bicgstab.h`.

**In Development**:
*   **HACApK for PEEC**: H-matrix acceleration for large PEEC systems (L matrix compression).
*   **Application Library**: Reference examples for MagLev, WPT, and Accelerator magnets.

---

## 💎 Strategic Value: Solving What FEM Cannot

**Closing the Gap in Computational Electromagnetics.**

Commercially available Finite Element Method (FEM) tools are powerful, but they face inherent limitations when dealing with open regions and moving parts. Radia provides a **Complementary Framework** based on **Integral Methods (Green's Functions / Kernels)** to solve these specific classes of problems effectively.

*   **The "Open Boundary" Problem**: FEM requires truncating the universe with artificial boundaries (or expensive infinite elements).
    *   *Our Solution*: Integral methods naturally satisfy the condition at infinity. No air mesh is needed.
*   **The "Moving Source" Problem**: Moving a coil or magnet in FEM requires complex re-meshing or sliding interfaces, introducing numerical noise.
    *   *Our Solution*: Sources are analytical objects. Moving them is a simple coordinate transformation, free of discretization error.

**We do not replace FEM; We are the Source Provider for FEM.**
Radia actively feeds NGSolve and ngsolve.bem with high-fidelity electromagnetic sources — coils, magnets, and conductor impedances — so that FEM can focus on what it does best: solving material-dominated PDEs (saturation, eddy currents, thermal coupling). The integration is not a loose coupling; Radia's fields are evaluated **inside** NGSolve's assembly loop as native `CoefficientFunction` objects.

---

## ⚡ Paradigm Shift: Surface-Based Physics

**Volume Meshing is Obsolete for Conductors.**

For high-frequency applications (WPT, Induction Heating, Accelerators), traditional FEM struggles with the **Multi-Scale Challenge**:
*   **Macro Scale**: Large air gaps (meters)
*   **Micro Scale**: Skin depth (microns)

Attempting to mesh both simultaneously results in massive element counts and slow convergence. **We reject this approach.**

**The Radia PEEC Solution: SIBC + MKL LAPACK**
We solve the physics exactly where it happens: **On the Surface.**

1.  **SIBC (Surface Impedance Boundary Condition)**: Mathematical modeling of skin effect physics directly on the boundary. No internal mesh is required inside the conductor.
2.  **PEEC + MKL LAPACK**: C++ MNA (Modified Nodal Analysis) solver with Intel MKL LAPACK (`zgesv_`, `zgetrf_`, `zgetrs_`) and templated BiCGSTAB for fast impedance extraction.
3.  **FastHenry Compatibility**: Parse `.inp` files directly, including multi-filament (`nwinc`/`nhinc`) for skin/proximity effect.

**Result**: Direct circuit parameter extraction (L, R, C, M) from conductor geometry, with SPICE-ready output.

## 🦁 Academic Heritage & Citations

Radia is not a new invention; it is the **Modern Evolution** of battle-tested scientific codes developed at world-leading research institutes.
We stand on the shoulders of giants:

*   **Radia (ESRF)**: Developed by **O. Chubar, P. Elleaume, et al.** at the European Synchrotron Radiation Facility. The standard for undulator design for decades.
    *   *Ref*: O. Chubar, P. Elleaume, J. Chavanne, "A 3D Magnetostatics Computer Code for Insertion Devices", J. Synchrotron Rad. (1998).
*   **FastImp (MIT)**: Developed by **J. White, et al.** at MIT. The pioneer of pFFT-accelerated Surface Integral Equation methods.
    *   *Ref*: Z. Zhu, B. Song, J. White, "Algorithms in FastImp: A Fast and Wide-Band Impedance Extraction Program", DAC (2003).
*   **[HACApK](https://github.com/RIKENGITHUB/ppOpen-HPC) (JAMSTEC/RIKEN)**: Developed by **A. Ida, et al.** at JAMSTEC. Hierarchical matrices with Adaptive Cross Approximation for Krylov solvers.
    *   *Ref*: A. Ida, T. Iwashita, T. Mifune, Y. Takahashi, "Parallel Hierarchical Matrices with Adaptive Cross Approximation on Symmetric Multiprocessing Clusters", J. Inf. Process., Vol. 22, No. 4, pp. 642–650 (2014).

---

## 📐 Mathematical Foundations: The Power of Analytical Kernels

The core advantage of **Integral Element Method (IEM)** is the use of **Analytical Integration** over source volumes and surfaces, eliminating discretization error.

### 1. Analytical Sources (Radia Kernels)
Instead of approximating a coil as a bundle of sticks, we analytically integrate the Bio-Savart law:

$$ \vec{B}(\vec{r}) = \frac{\mu_0 I}{4\pi} \int_{Volume} \vec{J}(\vec{r}') \times \frac{\vec{r} - \vec{r}'}{|\vec{r} - \vec{r}'|^3} dV' $$

For specific geometries, this yields **Exact Closed-Form Solutions**:
*   **Polygonal Coils**: Exact integration of straight segments.
*   **Arc Segments**: Exact integration of circular arcs.
*   **Cylindrical Magnets**: Exact field formulas involving elliptic integrals.
*   **Polyhedral Magnets**: Exact surface charge integration ($\sigma_m = \vec{M} \cdot \vec{n}$).

### 2. Method of Magnetized Methods (MMM) with MSC
For iron saturation, we employ the **Magnetic Surface Charge (MSC)** formulation. The magnetization $\vec{M}$ inside a volume $\Omega$ creates an equivalent surface charge density:

$$ \phi_m(\vec{r}) = \frac{1}{4\pi} \oint_{\partial \Omega} \frac{\vec{M} \cdot \vec{n}'}{|\vec{r} - \vec{r}'|} dS' - \frac{1}{4\pi} \int_{\Omega} \frac{\nabla' \cdot \vec{M}}{|\vec{r} - \vec{r}'|} dV' $$

### 3. Surface Impedance & FastImp Kernels (MQS/Darwin Regime)
For conductor analysis, we solve the Surface Integral Equation (SIE) using the **Laplace kernel**:

$$ G(\vec{r}, \vec{r}') = \frac{1}{4\pi|\vec{r} - \vec{r}'|} $$

**Supported Frequency Regime**: Magneto-Quasi-Static (MQS) to Darwin approximation.
- **MQS**: Ignores displacement current ($\partial D/\partial t \approx 0$). Valid when $\lambda >> L$ (wavelength >> problem size).
- **Darwin**: Includes inductive effects but ignores radiation. Valid for $kL << 1$ where $k = \omega/c$.

Combined with **SIBC (Surface Impedance Boundary Condition)**, this reduces the volumetric skin-effect problem to a purely surface-based boundary element problem.

> [!NOTE]
> **Full-wave Helmholtz kernel** ($e^{ikr}/r$) has been removed. Radia targets MQS/Darwin applications (MagLev, WPT, Induction Heating) where wavelength >> device size, making the quasi-static approximation highly accurate and computationally efficient.

---

## 🧘 Philosophy: The Source Provider for NGSolve

**Radia exists to give NGSolve and ngsolve.bem the best possible electromagnetic sources.**

The division of labor is clear:

| Role | Engine | What it computes |
| :--- | :--- | :--- |
| **Source Provider** | **Radia** | $H_s$, $B_s$, $A_s$, $\Phi_s$ from coils, magnets, conductors — analytically |
| **Material Solver** | **NGSolve** | Reaction fields in iron/dielectric via FEM ($\nabla \cdot \mu \nabla \phi = -\nabla \cdot \mu H_s$) |
| **Eddy Current Solver** | **ngsolve.bem** | Surface eddy currents via BEM (HDivSurface $\times$ SurfaceL2) |
| **Circuit Extraction** | **Radia PEEC** | Impedance (L, R, C, M) from conductor geometry for SPICE |

**How the integration works**:
1.  Radia computes the source field analytically (no mesh, no discretization error).
2.  `rad.RadiaField()` wraps the result as an NGSolve `CoefficientFunction` — a native C++ object evaluated directly inside NGSolve's element assembly loop. Since v2.5.0, `RadiaField` is integrated into the main `_radia_pybind.pyd` module; no separate `radia_ngsolve` module is needed.
3.  NGSolve / ngsolve.bem solve the reaction problem driven by this source.
4.  Result = Source Field (Radia) + Reaction Field (NGSolve).

We bridge the gap between distinct mathematical communities:
*   **Integral Codes**: Radia (ESRF) & FastImp (MIT) $\rightarrow$ *Analytical source fields.*
*   **Finite Element Codes**: NGSolve (TU Wien) & ngsolve.bem $\rightarrow$ *Material & eddy current solvers.*

---

## 🤖 LLM-Agent Ready & Python Native

**"No GUI? No Problem."**

We believe that **Natural Language is the ultimate User Interface** for complex design.
Instead of clicking through nested menus to find a "Halbach Array" button, you simply describe what you want.

*   **Code-First Modeling**: Geometry and physics are defined in pure, human-readable Python.
*   **The "Nanobanana" Vision**: By combining Radia with modern AI, we turn text prompts into rigorous engineering models.
    *   *Prompt*: "Create a Halbach array for a MagLev slider with 12 periods, optimized for 5mm levitation gap."
    *   *Result*: An Agent generates the complete executable Radia script, including geometric parameters and material definitions.

> [!TIP]
> **Why Python?** GUI-based tools are excellent for standard tasks, but they limit you to what the developer imagined. Python + Radia limits you only by Python's endless ecosystem.

*   **Ecosystem Integration**: Seamlessly integrates with the rich Python scientific stack (NumPy, SciPy, PyVista, NGSolve) and modern version control (Git).

---

## 💡 Architecture: Source (IEM) + Material (FEM) + Eddy Current (BEM)

We define our unique approach as a hybrid of **Integral Element Method (IEM)**, **Finite Element Method (FEM)**, and **Boundary Element Method (BEM)**.

**What is "Integral Element Method (IEM)"?**
Unlike FEM, which uses uniform element formulations, IEM allows the combination of elements with **different integration kernels** (e.g., $1/r$ for monopoles, $\vec{J} \times \vec{r}/r^3$ for Biot-Savart) into a single system. All kernels use the **Laplace form** ($1/r$) for quasi-static analysis.

| Layer | Method | Role & Kernels | Advantage |
| :--- | :--- | :--- | :--- |
| **Source Layer** | **IEM** (Radia) | **Laplace Kernels** ($1/r$): Volume Magnets, Coils, SIBC Surfaces. Provides $H_s$, $B_s$, $A_s$ as `CoefficientFunction`. | **Composable & Analytical.** No mesh needed for sources. |
| **Material Layer** | **FEM** (NGSolve) | **Differential Operators** ($\nabla \cdot \mu \nabla$): Saturation, Hysteresis, Thermal. | **Non-Linear & Multi-Physics.** |
| **Eddy Current Layer** | **BEM** (ngsolve.bem) | **Surface BEM** (HDivSurface $\times$ SurfaceL2): Eddy currents, shielding. | **No volume mesh in conductors.** |

**The Workflow:**
1.  **Radia**: Computes the source field ($H_s$ or $T_s$) analytically.
2.  **NGSolve**: Solves for the reaction potential ($\phi$) in the iron regions using FEM.
    *   $\nabla \cdot (\mu \nabla \phi) = -\nabla \cdot (\mu H_s)$
    *   **Frequency Range**: Primarily targets **Low Frequency** (Magnetostatics / Eddy Currents), shielding, and extending up to the **Darwin Regime** (ignoring radiation, but including displacement currents if needed).
3.  **Result**: Superposition of Source Field + Reaction Field.

> [!NOTE]
> **Design**: The coupling is intentionally one-way (Radia Sources $\rightarrow$ NGSolve/ngsolve.bem). Radia provides the source; NGSolve and ngsolve.bem solve the reaction. This clean separation keeps the architecture composable and each solver optimal for its role.

### NGSolve Integration Details (Weak Coupling Mechanism)
The `rad.RadiaField()` function (integrated into `_radia_pybind.pyd` since v2.5.0) implements a high-performance **Weak Coupling** bridge using a native C++ `CoefficientFunction`. This allows Radia fields to be evaluated directly during NGSolve's finite element assembly process.

**Implementation Architecture:**
*   **Native C++ Shim**: A `RadiaFieldCF` class (inheriting from `ngfem::CoefficientFunction`) sits between NGSolve and Radia.
*   **Three-Tier Evaluation Strategy**:
    1.  **Fast FMM (C++)**: For `B`, `H`, and `A` fields, dipoles are extracted from Radia and evaluated using a C++ Fast Multipole Method (FMM) solver. This **bypasses the Python Global Interpreter Lock (GIL)** entirely, enabling maximum performance during massive parallel FEM assembly.
    2.  **Cached Evaluation**: A coordinate-hash cache prevents redundant re-calculation of fields at the same integration points.
    3.  **Python Fallback**: For complex material responses (Magnetization `M`, Scalar Potential `Phi`), it safely acquires the GIL and calls the Radia Python kernel.

### NGSolve Primer for Radia Users
*   **CoefficientFunction (CF)**: A generic function that can be evaluated anywhere in the 3D domain. Radia provides the source Magnetic Field ($H_s$) as a C++ `CoefficientFunction`. This means NGSolve can "query" Radia for the field value at any coordinate during matrix assembly **without needing to store values on a mesh** or interpolate from a grid.
*   **GridFunction (GF)**: A field defined on the finite element mesh (stored as vectors of coefficients). This typically represents the *solution* (like the Magnetic Potential $\phi$) or the *material property distribution* (like Permeability $\mu$) in the FEM model.

---

## Key Capabilities

### 1. Integrated Field Sources
Instead of simple "boundary conditions", Radia provides rich physical sources:

*   **Permanent Magnets**: Analytical surface charge method (Polyhedrons, Extrusions).
*   **Moving Magnets & Coils**: Sources can have arbitrary position and orientation transformations applied dynamically.
    *   *Development Status*: Comprehensive dynamic simulation examples and animation workflows are currently being developed.
*   **Coils & Current Loops**: Biot-Savart integration for arbitrary paths.
*   **Distributed Currents**: Arc segments, race-tracks, and helical filaments.
*   **Analytical Precision**: To eliminate source errors, **fully analytical formulas** are used wherever possible (e.g., exact integration for straight/arc segments, analytical surface charges) rather than approximate numerical integration.
*   **Versatile Field Types**: Supports computation of **A** (Vector Potential), **Phi** (Scalar Potential), **B** (Flux Density), and **H** (Field Intensity) to drive various FEM formulations ($A$-formulation, Reduced-Scalar-Potential, etc.).


### 2. High-Performance Solvers & Acceleration
To handle complex field sources efficiently, the framework employs state-of-the-art acceleration algorithms based on the **Laplace kernel** ($1/r$):

*   **Solver Acceleration (Source Definition)**:
    *   **$\mathcal{H}$-Matrix ([HACApK](https://github.com/RIKENGITHUB/ppOpen-HPC) ACA+)**: Used for Magnetostatics (MMM). Compresses dense interaction matrices to $O(N \log N)$, enabling large-scale iron/magnet simulations.
    *   **PEEC + MKL LAPACK**: C++ PEEC solver with Intel MKL LAPACK/BLAS for circuit parameter extraction. SIBC models skin depth effects as surface properties. Templated BiCGSTAB shared between MSC (real) and PEEC (complex) solvers.
*   **Field Evaluation Acceleration**:
    *   **FMM (ExaFMM-t)**: Fast Multipole Method using Laplace kernel for rapidly computing fields ($B, H, A$) from massive numbers of source elements. This is critical for the `CoefficientFunction` interface to NGSolve.
*   **Hybrid FEM**: Reduced Potential coupling with NGSolve.

> [!NOTE]
> **All acceleration methods use Laplace kernel** ($1/r$). This ensures consistency across the framework and optimal performance for MQS/Darwin applications.


### 3. Visualization & Export
*   **PyVista Viewer**: Modern, interactive 3D visualization within Python/Jupyter.
*   **VTK Export**: Compatible with ParaView.
*   **GMSH/STEP**: Mesh import via GMSH, CAD interoperability via Coreform Cubit (integrated `cubit_mesh_export` module).

---

### 4. MagLev Specific Capabilities
We provide built-in formulations for the unique physics of magnetic levitation:

*   **EDS (Electrodynamic Suspension)**:
    *   **Drag & Lift Forces**: Accurate computation of velocity-dependent forces on moving magnets over conductive plates (using `rad.ObjMpl` or FastImp).
    *   **Inductrack**: Simulation of Halbach arrays moving over passive coils or litz-wire tracks.
*   **EMS (Electromagnetic Suspension)**:
    *   **Control Inductances**: Fast extraction of differential inductance matrices ($L_{ij}$) for control loop design (differentiate Flux $\Phi$ w.r.t current $I$).
    *   **Force-Gap Characteristics**: High-precision force vs. air-gap curves for nonlinear controller tuning.

## ⚖️ Workflow Comparison: Why Switch?

| Feature | Traditional FEM (Commercial) | **Radia Framework (IEM + FEM)** |
| :--- | :--- | :--- |
| **Air Mesh** | **Required.** Must mesh the "nothingness" around the device. | **None.** Air is handled analytically. |
| **Moving Parts** | **Hard.** Mesh deformation, sliding interfaces, re-meshing noise. | **Trivial.** Just apply a coordinate transform `rad.Trsf`. |
| **Coil Geometry** | **Approximated.** Step-files or coarse filaments. | **Exact.** Analytical arcs, straight segments, and volumes. |
| **Skin Effect** | **Heavy.** Requires dense volume mesh inside conductors. | **Light.** SIBC solves it on the surface only. |
| **Optimization** | **Blackbox.** Slow parameters sweeps via GUI. | **Transparent.** Fast, gradient-friendly Python execution. |

---

## Quick Start

### Installation

```bash
# Windows (Python 3.12)
pip install radia
```

*Prerequisites for FEM features: `pip install ngsolve`*

### Example 1: Magnetostatic Source Field

```python
import radia as rad

# Define a Race-Track Coil
coil = rad.ObjRaceTrk(
    [0,0,0],       # Center
    [10, 30],      # Inner Radii (R_min, R_max)
    [20, 100],     # Straight section lengths (Lx, Ly)
    10.0,          # Height
    3.0,           # Curvature radius
    1000.0,        # Current [A]
    'man'          # Manually defined rectangular cross-section
)

# Export field to VTS for ParaView visualization
           [-50, 50], [-150, 150], [-20, 30],  # x, y, z ranges [mm]
           21, 31, 11)  # grid points
```

### Example 2: PEEC Circuit Parameter Extraction

```python
from radia.peec_matrices import PyPEECBuilder
from radia.peec_topology import PEECCircuitSolver
import numpy as np

# Build a simple inductor: 4 segments in series
builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.05, 0, 0)
n3 = builder.add_node_at(0.05, 0.05, 0)
n4 = builder.add_node_at(0, 0.05, 0)
for na, nb in [(n1,n2), (n2,n3), (n3,n4), (n4,n1)]:
    builder.add_connected_segment(na, nb, w=1e-3, h=1e-3, sigma=5.8e7, nwinc=3, nhinc=3)
builder.add_port(n1, n1)  # Single-turn loop

topo = builder.build_topology()
solver = PEECCircuitSolver(topo)

# Extract impedance vs frequency
freqs = np.logspace(2, 6, 20)
Z = solver.frequency_sweep(freqs)
R = np.real(Z)
L = np.imag(Z) / (2 * np.pi * freqs)
print(f"DC: R={R[0]*1e3:.2f} mOhm, L={L[0]*1e9:.1f} nH")
```

### Example 3: FastHenry .inp Import

```python
from radia.fasthenry_parser import FastHenryParser

parser = FastHenryParser()
parser.parse_string("""
.Units mm
N1 x=0 y=0 z=0
N2 x=100 y=0 z=0
E1 N1 N2 w=1 h=1 sigma=5.8e7 nwinc=5 nhinc=5
.external N1 N2
.freq fmin=100 fmax=1e6 ndec=5
.end
""")
result = parser.solve()
print(f"DC: R={result['R'][0]*1e3:.3f} mOhm, L={result['L'][0]*1e9:.1f} nH")
```

---

## Documentation & Resources

*   **[Installation Guide](BUILD.md)**: Build from source (Windows/Linux/macOS).
*   **[API Reference](docs/api/API_REFERENCE.md)**: Full Python API documentation.
*   **[NGSolve Integration](docs/NGSOLVE_INTEGRATION.md)**: Theory and usage of the hybrid FEM-Integral method.
*   **[ngsolve.bem Integration](docs/solver/NGBEM_INTEGRATION_DESIGN.md)**: Eddy current solver via ngsolve.bem.
*   **[Original Radia](https://github.com/ochubar/Radia)**: The core physics engine developed at ESRF.

---

## Cubit Mesh Export

Radia includes a mesh export module for [Coreform Cubit](https://coreform.com/products/coreform-cubit/), supporting **Gmsh**, **Nastran BDF**, **VTK/VTU**, **Netgen/NGSolve (.vol)**, **Exodus II**, and **MEG** formats.

### Quick Usage

```python
import cubit
import cubit_mesh_export

# Export to Gmsh v2.2 (for Radia import) or v4.1 (for visualization)
cubit_mesh_export.export_gmesh(cubit, "model.msh", version="2.2")

# Export to NGSolve mesh with high-order curving
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry_file="geometry.step")
mesh = ngsolve.Mesh(ngmesh)
mesh.Curve(3)  # 3rd order curved elements
```

### Installation

The `cubit_mesh_export` module is included in Radia:

```bash
pip install radia
```

To install the Cubit toolbar panels for GUI access (optional, requires Cubit installed):

```bash
cubit-mesh-export-install-panels
```

This registers a toolbar in `~/.cubit` so that export buttons appear in Cubit's GUI on next startup. To uninstall:

```bash
cubit-mesh-export-install-panels --uninstall
```

### Setup

Cubit scripts require the `cubit` Python module. Either:
1. Add Cubit's `bin` directory to your system PATH, or
2. Set the `CUBIT_PATH` environment variable:

```bash
# Windows
set CUBIT_PATH=C:\Program Files\Coreform Cubit 2025.3\bin

# Linux/Mac
export CUBIT_PATH=/opt/coreform/cubit/bin
```

**Important**: When using both NGSolve and Cubit in the same script, import NGSolve **before** Cubit to avoid DLL conflicts.

### Documentation

- **[Function Reference](docs/cubit/Function_Reference.md)** -- Full API documentation
- **[Examples](examples/cubit/)** -- Export examples for all supported formats

## License

*   **Radia Core**: BSD-style (ESRF)
*   **$\mathcal{H}$-Matrix Library ([HACApK](https://github.com/RIKENGITHUB/ppOpen-HPC))**: MIT (ppOpen-HPC/JAMSTEC)

---
*Radia: Empowering the next generation of magnetic system design.*
