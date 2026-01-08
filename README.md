# Radia - Electromagnetic Simulation Framework for Magnetic Levitation Systems

**A Python-native environment designed for the age of AI-driven Engineering.**

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

### Future Scope & Active Development
We are actively expanding the framework to cover:
*   **ESIM (Equivalent Source Integral Method)**: Currently prioritizing the implementation of ESIM for advanced source modeling.

---

## ⚡ Paradigm Shift: Surface-Based Physics

**Volume Meshing is Obsolete for Conductors.**

For high-frequency applications (WPT, Induction Heating, Accelerators), traditional FEM struggles with the **Multi-Scale Challenge**:
*   **Macro Scale**: Large air gaps (meters)
*   **Micro Scale**: Skin depth (microns)

Attempting to mesh both simultaneously results in massive element counts and slow convergence. **We reject this approach.**

**The Radia/FastImp Solution: SIBC + pFFT**
We solve the physics exactly where it happens: **On the Surface.**

1.  **SIBC (Surface Impedance Boundary Condition)**: Mathematical modeling of skin effect physics directly on the boundary. No internal mesh is required inside the conductor.
2.  **pFFT (Precorrected-FFT)**: Accelerates the dense matrix interactions to $O(N \log N)$.

**Result**: Simulations that took hours with FEM finish in minutes, with perfect geometric fidelity for Litz wires and complex coils.

---

## 🧘 Philosophy: Physics-First in a Multi-Physics World

We believe in **"The Right Tool for the Right Physics"**.

While Integral Methods (Radia) are superior for open-boundary magnetics, we recognize that modern engineering requires **Multi-Physics** (Thermal, Structural, Fluid). The world is not just magnetic; iron saturates, coils heat up, and structures deform.

**Our Vision**: Radia does not try to do everything. Instead, it acts as the **Precision Field Generator** within a larger Multi-Physics workflow.

*   **Radia**: Generates the exact electromagnetic sources (Coils, Magnets) analytically.
*   **NGSolve (FEM)**: Handles the multi-physics material response (Heat, Stress, Non-linear Saturation).

By coupling these dedicated solvers, we achieve a system that is both **Faster** (efficient algorithms) and **More Accurate** (no air mesh errors) than monolithic FEM approaches.

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

## 💡 Architecture: The "Best of Both Worlds"

We combine two powerful mathematical approaches into a single cohesive framework:

| Layer | Technology | Role | Advantage |
| :--- | :--- | :--- | :--- |
| **Source Layer** | **Radia** (Integral Method) | Defines Coils, Magnets, Current paths. | **Zero Meshing of Air.** Infinite boundaries are handled analytically. Perfect representation of curves. |
| **Material Layer** | **NGSolve** (FEM) | Defines Iron Yokes, Shields, Conductors. | Handles non-linear saturation (B-H curves) and eddy currents using **Reduced Scalar Potential**. Combined with **Kelvin Transformation**, it efficiently handles open boundaries for the reaction field. |

**The Workflow:**
1.  **Radia**: Computes the source field ($H_s$ or $T_s$) analytically.
2.  **NGSolve**: Solves for the reaction potential ($\phi$) in the iron regions using FEM.
    *   $\nabla \cdot (\mu \nabla \phi) = -\nabla \cdot (\mu H_s)$
    *   **Frequency Range**: Primarily targets **Low Frequency** (Magnetostatics / Eddy Currents), shielding, and extending up to the **Darwin Regime** (ignoring radiation, but including displacement currents if needed).
3.  **Result**: Superposition of Source Field + Reaction Field.

> [!NOTE]
> **Limitation**: Strong coupling with FEM is **not currently supported**. The integration is presently one-Way (Radia Sources $\rightarrow$ NGSolve).

### NGSolve Integration Details (Weak Coupling Mechanism)
The `radia_ngsolve` module implements a high-performance **Weak Coupling** bridge using a native C++ `CoefficientFunction`. This allows Radia fields to be evaluated directly during NGSolve's finite element assembly process.

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
To handle complex field sources efficiently, the framework employs state-of-the-art acceleration algorithms:

*   **Solver Acceleration (Source Definition)**:
    *   **H-Matrix**: Used for Magnetostatics (MMM). Compresses dense interaction matrices to $O(N \log N)$, enabling large-scale iron/magnet simulations.
    *   **pFFT & SIBC**: Used for Conductor Analysis (FastImp). **Surface Impedance Boundary Conditions (SIBC)** combined with Precorrected-FFT allow extremely fast impedance extraction by modeling skin depth effects as surface properties.
*   **Field Evaluation Acceleration**:
    *   **FMM (Fast Multipole Method)**: Used for rapidly computing fields ($B, H, A$) from massive numbers of source elements. This is critical for the `CoefficientFunction` interface to NGSolve.
*   **Hybrid FEM**: Reduced Potential coupling with NGSolve.


### 3. Visualization & Export
*   **PyVista Viewer**: Modern, interactive 3D visualization within Python/Jupyter.
*   **VTK Export**: Compatibile with ParaView.
*   **Nastran/Step**: Interoperability with CAD tools via [Coreform Cubit integration](https://github.com/ksugahar/Coreform_Cubit_Mesh_Export).

---

## Quick Start

### Installation

```bash
# Windows (Python 3.12)
pip install radia-ngsolve
```

*Prerequisites for FEM features: `pip install ngsolve`*

### Example: The "Agentic" Way

Modeling a complex coil doesn't require a GUI. It requires expressive code:

```python
import radia as rad

# Define a Race-Track Coil automatically
# An LLM can easily tweak parameters like 'current', 'radius', 'turns'
coil = rad.ObjRaceTrk(
    [0,0,0],       # Center
    [10, 30],      # Inner Radii (R_min, R_max)
    [20, 100],     # Straight section lengths (Lx, Ly)
    10.0,          # Height
    3.0,           # Curvature radius
    1000.0,        # Current [A]
    'man'          # Manually defined rectangular cross-section
)

# Visualize clearly
rad.ObjDrwVTK(coil, "coil_geometry.vtk")
```

---

## Documentation & Resources

*   **[Installation Guide](BUILD.md)**: Build from source (Windows/Linux/macOS).
*   **[API Reference](docs/API_REFERENCE.md)**: Full Python API documentation.
*   **[NGSolve Integration](docs/NGSOLVE_INTEGRATION.md)**: Theory and usage of the hybrid FEM-Integral method.
*   **[Original Radia](https://github.com/ochubar/Radia)**: The core physics engine developed at ESRF.

## License

*   **Radia Core**: BSD-style (ESRF)
*   **H-Matrix Library**: MIT (ppOpen-HPC)

---
*Radia: Empowering the next generation of magnetic system design.*
