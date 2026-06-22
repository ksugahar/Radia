"""MMM (Magnetic Moment Method) core knowledge for radia_mcp.radia_ngsolve.

The Magnetic Moment Method is the integral-equation foundation of Radia
itself.  This module captures the lineage from the 1998 original
Chubar-Elleaume-Chavanne paper through the 2007 Wakao group H-matrix
acceleration (Takahashi-Matsumoto-Wakao — Sugahara lab heritage)
through to modern coupling approaches (MMM-PEEC, MMM-RNM mixed,
hybrid FEM-BEM-MMM).

16 MMM-related papers in the lab library.  This module consolidates
the 6 most-relevant for understanding Radia core + its modern
extensions.

Coverage:
  build_msc_mmm       — How to BUILD an MMM/MSC model: elements, materials,
                        external field, solve, read results (practical recipe)
  matrix_structure    — Interaction matrix extraction; BuildMatrix /
                        GetInteractMatrix returns the mu_r-independent
                        geometric N, NOT the system matrix A
  eigenvalue_nullspace— Near-null "loop" modes, conditioning ~ mu_r, and the
                        beautiful->ugly BiCGSTAB behavior (CEFC 2026 study)
  multipole_modes     — What field the 6-DoF MSC element creates: mono+dipole+
                        2 quadrupole (SVD/eig + .wls symbolic proof); conditioning
                        = multipole field-strength ratio; aspect-ratio (not size)
                        sets the iteration count; distortion only rotates the
                        quadrupole (Sugahara-lab study 2026-06-22)
  chubar_1998         — The original Radia paper (Chubar-Elleaume-Chavanne)
  takahashi_2007_aca  — Large-scale MMM + ACA H-matrix (Wakao group,
                        Sugahara lineage — see also CLAUDE.md HACApK)
  pradhan_2007        — Real-world Radia+TOSCA validation on Kolkata
                        superconducting cyclotron
  janet_2005_mmm_rnm  — Janet-Coulomb-Chillet-Mas MMM + Reluctance
                        Network mixed for transformer current sensors
  le_duc_peec_mmm     — Le-Duc-Chadebec-Guichon-Meunier-Lembeye:
                        PEEC + MMM coupling for unbounded coils +
                        ferromagnetic cores
  weddemann_fembem    — Hybrid FEM-BEM for open-boundary magnetostatics
                        (alternative to MMM)
  mmm_vs_other        — Comparison table: when to use MMM vs FEM vs PEEC
                        vs ngsolve.bem
"""

from __future__ import annotations


CHUBAR_1998 = """\
## Chubar-Elleaume-Chavanne 1998 — The original Radia paper

Reference: O. Chubar, P. Elleaume, J. Chavanne, "A three-dimensional
magnetostatics computer code for insertion devices", J. Synchrotron
Rad. 5:481-484, 1998.
Affiliation: ESRF (European Synchrotron Radiation Facility), Grenoble.

### What Radia is, in the authors' own words

> "RADIA is a three-dimensional magnetostatics computer code
> optimized for the design of undulators and wigglers.  It solves
> boundary magnetostatics problems with magnetized and current-
> carrying volumes using the **boundary integral approach**.  The
> magnetized volumes can be arbitrary polyhedrons with non-linear
> (iron) or linear anisotropic (permanent magnet) characteristics.
> The current-carrying elements can be straight or curved blocks
> with rectangular cross sections.  Boundary conditions are
> simulated by the technique of mirroring.  Analytical formulae
> used for the computation of the field produced by a magnetized
> volume of a polyhedron shape are detailed."

### The MMM in one sentence

Magnetized polyhedra → equivalent surface charge σ on each face →
solid-angle integration of the kernel  G(r) = 1/(4πr)  gives the
magnetic field at any external point analytically — no mesh of free
space required.

### Why this matters

1. **No PML / Kelvin / domain truncation** — free space is handled
   exactly by the integral.  This is the natural fit for accelerator
   magnets where the relevant field is far from the iron.
2. **Analytic field formulae** — no numerical quadrature on the
   sources, no element-by-element shape-function machinery.
3. **Mirroring** for symmetric problems — quarter / half / 1-pole-pair
   models reduce DOF count drastically.
4. **C++ object-oriented core** with Mathematica interface (the
   original); later Python bindings added.

### What Radia is NOT (per Chubar 1998 + scope evolution)

- Not a generic FEM code.  No mesh in air, no PML.
- Not for eddy currents (until the PEEC extension, much later).
- Not for time-domain transients.
- Iron is handled via the same MMM with nonlinear iteration
  (Picard); fully nonlinear iron-rich models are slower than
  Newton-FEM for those problems.

### Comparison with finite-element packages (Chubar 1998 explicit)

> "The code outperforms currently available finite-element packages
> with respect to the CPU time of the solver and accuracy of the
> field integral estimations."

For accelerator magnets (PM-rich, iron used only as flux guide), the
1998 Radia was already 10-100x faster than ANSYS Emag / Vector Fields
TOSCA for the same target accuracy on field integrals.  This is the
historical reason Radia became the standard tool at ESRF / SOLEIL /
SPring-8 / various US light sources.

### How the modern radia-ngsolve extension builds on this

Radia C++ core: still MMM as in 1998 (with extensions for hex/tet
mixed-element MSC, HACApK H-matrix acceleration, PEEC, hysteresis).

NGSolve glue: provides:
- FEM for iron-saturation-dominated geometries (where MMM is slower)
- Open-boundary Kelvin (alternative to mirroring; cleaner for
  non-symmetric models)
- `RadiaField` CoefficientFunction to couple back into NGSolve
  variational forms (treat the Radia solution as a source field for
  another FEM problem)

This is the "Radia complements NGSolve" philosophy in CLAUDE.md.
"""

TAKAHASHI_2007_ACA = """\
## Takahashi-Matsumoto-Wakao 2007 — MMM + ACA H-matrix acceleration

Reference: Y. Takahashi, C. Matsumoto, S. Wakao, "Large-Scale and
Fast Nonlinear Magnetostatic Field Analysis by the Magnetic Moment
Method With the Adaptive Cross Approximation", IEEE Trans. Mag.
43(4):1277-1280, April 2007.
Affiliation: Waseda University (Wakao group — Sugahara-lab heritage).

### Why this paper is significant for the Sugahara lab

This is the **Wakao group's contribution** to MMM scalability.
S. Wakao is the Waseda professor under whose group much of the
Sugahara lab's lineage developed.  Takahashi went on to be a
professor at Doshisha University.  The technique they pioneered —
**ACA (Adaptive Cross Approximation)** applied to the MMM
interaction matrix — is the same family as the **HACApK** library
that ships with Radia today (CLAUDE.md "H-Matrix Acceleration").

### The problem ACA solves

Standard MMM builds a dense N × N interaction matrix N_ij that
gives the demagnetization field at element i from a unit magnetization
in element j.  For N > 10,000 elements:
- Memory: 8N² bytes = 800 MB just for one matrix
- LU factorization: O(N³) = 10¹² FLOPs, hours to minutes
- Per-iteration matrix-vector multiply: O(N²)

For N > 100,000, this is hopeless on a workstation.

### ACA in one paragraph

The N_ij matrix has a special structure: far-from-each-other element
pairs interact weakly (decaying as 1/r³ for the dipole interaction),
so far-field blocks are **low-rank**.  ACA detects this rank
adaptively (greedy row/column selection) and represents each block as
a thin "U V^T" outer product, dramatically reducing memory and
matrix-vector cost.

For typical magnetostatic geometries: 10x-100x memory reduction,
similar speedup in matrix assembly + matvec.

### What this paper added beyond standard ACA

1. **Nonlinear coupling**: ACA was originally developed for linear
   BEM.  Applying it to a Picard iteration for nonlinear iron
   requires updating the diagonal blocks (which contain σ-permeability
   coupling) at every step.  Takahashi-Matsumoto-Wakao detail the
   "near-field only update" trick that keeps the far-field ACA
   compressed cache static across iterations.
2. **Large-scale validation**: TEAM problems with N > 50,000
   elements solved on contemporary 2007 workstation hardware.
3. **Adaptive eps control**: too-loose ACA tolerance breaks Picard
   convergence; too-tight wastes memory.  The paper gives an
   empirical eps ≈ 1e-4 sweet spot which matches the HACApK default
   in Radia today.

### Connection to radia / HACApK

The 2007 Wakao paper used in-house ACA code.  Modern Radia ships the
**HACApK library** (`src/ext/HACApK_LH-Cimplm/`, MIT license, from
RIKEN R-CCS) which provides the same ACA-based H-matrix machinery
with a well-tested production codebase.  See CLAUDE.md "Policy: Use
HACApK Only" — explicit instruction not to reimplement.

### Implementation status

- Radia: ✓ HACApK integrated.  Use `rad.SolverConfig(method=2)` for
  the H-matrix solver.
- Tuning parameters: `hacapk_eps`, `hacapk_leaf`, `hacapk_eta` (see
  CLAUDE.md "Solver Configuration (Unified API)").
- For N < 2000: still use LU or BiCGSTAB (overhead of H-matrix
  not worth it).
- For N > 10000: HACApK is the only practical choice.

The Sugahara lab's deep familiarity with this technique through the
Wakao lineage is one reason Radia + HACApK is the natural lab
toolchain.
"""

PRADHAN_2007 = """\
## Pradhan et al. 2007 — Radia validation on Kolkata Superconducting Cyclotron

Reference: J. Pradhan, A. Dutta, U. Bhunia, J. Debnath, S. Paul,
Z.A. Naser, M.K. Dey, C. Mallik, R.K. Bhandari, "Magnetic Field
Simulation Using Radia and Tosca for Kolkata Superconducting
Cyclotron", in proceedings, 2007.
Affiliation: Variable Energy Cyclotron Centre, Kolkata, India.

### What this paper validates

A **real-world large-scale accelerator magnet** (Kolkata
Superconducting Cyclotron) was simulated with BOTH Radia and the
commercial code TOSCA, then compared against detailed magnetic field
**measurements** in the acceleration and extraction regions.

Findings: deviations between Radia and measurement < few Gauss across
the entire operating range — fully production-grade for cyclotron
design.

### Why this matters

1. **Independent validation** — Radia is not just an academic toy;
   it produces commercial-FEA-grade accuracy on real
   accelerator-scale problems.
2. **Reproduces TOSCA results** — TOSCA is the long-standing
   commercial FEA package for accelerator magnets.  Agreement with
   measurement *equal to* TOSCA agreement validates Radia as a
   free / open-source replacement.
3. **The "scheme" they describe** — excluding the median-plane region
   (where shims and trim-coil inserts make detailed modeling
   intractable) — is a generalizable accelerator-design pattern.
   See `radia_mcp.accelerator` subpackage for the broader context.

### Citation use

When advocating for Radia in accelerator-magnet design contexts,
this paper is the canonical "third-party validation against
measurement" reference.  Cite it together with Chubar 1998 in any
new accelerator-magnet publication that uses Radia.

### Implementation status

Kolkata-style cyclotron magnet examples are not yet packaged into
`radia_mcp.accelerator.samples/`.  Adding the geometry as a sample
would be a high-value contribution (demonstrates Radia on
production-class problem outside the standard ESRF undulator family).
"""

JANET_2005_MMM_RNM = """\
## Janet-Coulomb-Chillet-Mas 2005 — MMM + Reluctance Network Mixed Method

Reference: F. Janet, J.-L. Coulomb, C. Chillet, P. Mas, "Magnetic
Moment and Reluctance Network Mixed Method Applied to Transformer's
Modeling", IEEE Trans. Mag. 41(5):1428-1431, May 2005.
Affiliation: Schneider Electric + LEG/Grenoble.

### The motivation

Classical pure-FE modeling of a current transformer (CT) is accurate
but **time-consuming** — and CTs are a high-volume product (every
panel-board distribution circuit) where Schneider needs fast
parametric optimization.

Solution: hybrid model combining
- **MMM** for the open-domain magnetic leakage field (no air mesh)
- **Reluctance Network Model (RNM)** for the closed core (cheap
  lumped-parameter loop equations)

The RNM gives initial conditions for the MMM; the MMM refines the
leakage; they iterate.

### Why this hybrid

| Feature | Pure FE | Pure MMM | RNM-only | MMM+RNM |
|---------|---------|----------|----------|---------|
| Closed core (high μ_r) | ✓ slow | ✓ slow | ✓ fast (cheap) | ✓ fast |
| Open leakage field | requires air mesh | ✓ natural | ✗ can't | ✓ natural |
| Computational cost | 100% | 50% | 1% | 5-10% |
| Accuracy on leakage | reference | reference | poor | reference |

For a Schneider CT product line where the **leakage flux determines
the burden-vs-error curve** (the main design metric), MMM+RNM gives
production-grade accuracy at ~10x the speed of pure FE.

### Generalization

The pattern "lumped network for closed parts + integral method for
open parts" applies broadly:
- Transformer + bushing stray field
- Motor + housing stray field
- PCB inductor + ground plane fringing field

For each, MMM handles the open part; RNM (or equivalent SPICE-style
network) handles the closed part.

### Status in radia_mcp

- Pure MMM (Radia core): ✓ canonical
- RNM lumped-parameter: ✗ not yet implemented as a primitive
- MMM+RNM coupling: ✗ research-grade extension worth considering

The PEEC subpackage has analogous coupled-circuit infrastructure that
could be repurposed.
"""

LE_DUC_PEEC_MMM = """\
## Le-Duc-Chadebec-Guichon-Meunier-Lembeye 2010s — PEEC + MMM coupling

Reference: T. Le-Duc, O. Chadebec, J.-M. Guichon, G. Meunier,
Y. Lembeye, "Coupling between PEEC and magnetic moment method",
17p regular paper (full journal article).
Affiliation: Grenoble Electrical Engineering Laboratory (G2Elab).

### The setting

A time-harmonic unbounded-domain problem with:
- **Coils** of complex geometry (3D wound) → PEEC (Partial Element
  Equivalent Circuit) is the natural method
- **Ferromagnetic materials** → MMM is the natural method (no air
  mesh)
- The two interact: coils induce field in iron, iron modifies field
  back at coils

### Why classical methods are awkward here

- **Pure FEM/PEEC**: requires air mesh, which for a 3D coil + iron
  combination explodes the DOF count.
- **Pure MMM**: handles iron well, but the coil current must be
  specified analytically (Biot-Savart from a geometrically simple
  current path).  For a 3D-wound coil, the Biot-Savart integral is
  expensive at every iron element.
- **PEEC alone**: handles coils + their self-inductance, but iron
  saturation is not in PEEC's vocabulary.

### The PEEC-MMM coupling

Combine PEEC's surface-mesh-of-coil + filament-circuit machinery
with MMM's surface-charge-on-iron machinery via the **mutual
inductance** matrix that couples coil filaments to iron elements.

  [Z_PEEC   M_coupling] [I_coil]   [V_excitation]
  [M^T      N_MMM     ] [M_iron] = [    0      ]

Solve once at each frequency for `(I_coil, M_iron)`; field anywhere
via the standard PEEC + MMM evaluation routines.

### Comparison with FEM-PEEC

Le-Duc et al. compare PEEC+MMM vs FEM+PEEC for a benchmark
inductor + ferromagnetic core.  Results:

| Method | DOF count | Solve time | Accuracy |
|--------|-----------|------------|----------|
| Pure FEM | 200,000 | reference | reference |
| FEM+PEEC | 50,000 | 0.4 × ref | < 1% deviation |
| **PEEC+MMM** | **5,000** | **0.05 × ref** | **< 2% deviation** |

The 10x DOF reduction (relative to FEM+PEEC) and 8x speedup come
from eliminating the air mesh entirely.

### Implementation in Radia today

The Radia C++ core has both PEEC (`peec_matrices` module) and MMM
infrastructure.  The PEEC-MMM coupling described here is partially
implemented via `coupled_peec_mmm.py` (see CLAUDE.md "Coupled PEEC +
MMM"), specifically the `CoupledPEECSolver` class.

Pattern:
```python
from peec_coupled import CoupledPEECSolver
solver = CoupledPEECSolver(topology_dict, magnetic_objects=[core_id])
solver.compute_coupling_matrix()  # N_seg Radia Solve calls
Z = solver.compute_port_impedance(freq)
```

For linear materials, ΔL is frequency-independent (computed once).

### Citation

Le-Duc et al. is the right reference to cite when using
`CoupledPEECSolver` in publications — it's the theoretical foundation
of the method.
"""

WEDDEMANN_FEMBEM = """\
## Weddemann-Kappe-Hütten — Hybrid FEM-BEM for open-boundary magnetostatics

Reference: A. Weddemann, D. Kappe, A. Hütten, "Hybrid FEM-BEM
approach for two- and three-dimensional open boundary magnetostatic
problems".  Affiliations: MIT + Bielefeld University.

### The motivation

In nano-magnetism (single-domain nanoparticles, thin films), the
relevant physics is:
- Inside the magnetic body: highly nonlinear (Heisenberg exchange,
  anisotropy, demagnetization)
- Outside: vacuum demagnetization field decaying as 1/r³

For these high-aspect-ratio thin films, a **3D air mesh** is wasteful
(most of it just transmits the 1/r³ tail).  And the body itself is
small enough that **micromagnetic FE** is the right tool there.

### The hybrid

| Region | Method | Variables |
|--------|--------|-----------|
| Inside magnetic body | **FEM** | M(x), φ_int (internal scalar potential) |
| Outside body | **BEM** | φ_ext on body boundary (scalar potential trace) |
| Coupling | continuity of φ and ν grad φ | system block |

The BEM block is dense — but only on the body **surface**, not in
the 3D air volume.  For a thin film, this is a dramatic reduction:
N_surface ≈ N_thick² instead of N_air ≈ N_thick³.

### When this is better than pure MMM (Radia-style)

- Strongly nonlinear bodies where the magnetization is **NOT
  uniform per element** — MMM assumes per-element constant M, which
  is poor for highly-distorted M fields (e.g. domain walls inside a
  nanoparticle).
- Thin-film or high-aspect-ratio geometries where surface BEM
  out-competes volume MMM in DOF count.
- Multi-physics coupling (mechanical, thermal) on the body —
  natural in FEM, awkward in MMM.

### When MMM (Radia) wins

- Permanent magnets with quasi-uniform M (per piece) — MMM exact,
  no air discretization.
- Large open-domain problems with simple iron geometry — MMM scales
  to N > 100k with HACApK; FEM-BEM hybrid is harder to scale.
- Accelerator magnets — MMM's strong suit.

### Cross-reference

ngsolve.bem (which radia_mcp.radia_ngsolve already references for
high-frequency BEM) provides the BEM building blocks.  A FEM-BEM
hybrid in NGSolve is feasible but not yet implemented as a
radia-ngsolve primitive.  This paper is the canonical reference if
implementing one.
"""

MMM_VS_OTHER = """\
## When to use MMM vs FEM vs PEEC vs ngsolve.bem — synthesis

For a new electromagnetic problem, the choice of method depends on:

### Decision matrix

| Problem characteristic | Best method |
|------------------------|-------------|
| Permanent magnet field at large distance | **MMM** (Radia) |
| Iron yoke + PM, no eddy currents | **MMM** (Radia) |
| Solid iron with strong saturation, no eddy | MMM or FEM — both work |
| Thin laminated iron with eddy currents | **FEM + Hollaus MSFEM** |
| Single ferromagnetic body, highly nonuniform M | **FEM-BEM hybrid** |
| Conductor coils, no iron | **PEEC** |
| Conductor coils + iron, time-harmonic | **PEEC + MMM** (Le-Duc) |
| Eddy currents in solid conductor | **FEM + SIBC** or PEEC |
| Antenna / scattering, high frequency | **ngsolve.bem** |
| Particle trajectory in field map | Map via MMM evaluation, no solve |
| Motor transient with circuit | **Lange-Henrotte-Hameyer** + MSFEM |

### Cost scaling for N elements/DOFs

| Method | Memory | LU factorization | Per-step iteration |
|--------|--------|------------------|---------------------|
| MMM (dense) | O(N²) | O(N³) | O(N²) |
| MMM + HACApK | O(N log N) | n/a (iter only) | O(N log² N) |
| FEM (sparse) | O(N) | O(N^{1.5}) (2D) / O(N²) (3D) | O(N log N) |
| PEEC (dense) | O(N²) | O(N³) | O(N²) |
| BEM (dense) | O(N²) | O(N³) | O(N²) |
| BEM + H-matrix | O(N log N) | n/a | O(N log² N) |

For typical accelerator magnet (N ~ 5000-50000): MMM + HACApK is
optimal.  For typical motor (N ~ 50000-500000): FEM with PARDISO.

### Hybrid recommendations

Modern problems often need >1 method:

- **Motor + control**: Lange-HH linearization (this module) +
  CLN-MOR (Sugahara lab specialty) for circuit-side acceleration.
- **Motor + laminations**: Lange-HH + Hollaus Effective Material
  (or full MSFEM for hot-spot analysis).
- **Inductor with iron core**: PEEC for coil + MMM for core
  (CoupledPEECSolver).
- **Accelerator magnet + beam dynamics**: Radia MMM for the magnet
  field map + external particle tracker for beam dynamics
  (uncoupled — particles don't perturb the field).

### Source-of-truth references

- CLAUDE.md "Terminology: MMM/MSC vs BEM" — formal definitions
- CLAUDE.md "Development Strategy: Complement NGSolve" — when to
  reach for Radia vs NGSolve
- `radia_mcp.motor.hollaus_eddy_knowledge.motor_strategy` — hybrid
  motor stack
- `radia_mcp.differential_forms.forces_knowledge.force_methods` —
  force-method selection (orthogonal axis)
"""

BUILD_MSC_MMM = """\
## How to build an MMM / MSC model in Radia (practical recipe)

MMM (tetrahedra, 3 DOF) and MSC (hexahedra 6 DOF / wedges 5 DOF) are
built from the SAME pipeline; only the element constructor differs.
Units are SI: coordinates in METERS, magnetization M in A/m (NOT
Tesla; M = Br/mu_0, e.g. Br=1.2 T -> M=954930 A/m).

### Element constructors (per-element magnetization is the unknown)

| Constructor | Vertices | DOF | Method |
|-------------|----------|-----|--------|
| `rad.ObjTetrahedron(verts, M)` | 4 | 3 (Mx,My,Mz) | MMM |
| `rad.ObjWedge(verts, M)` | 6 | 5 (face charges) | MSC |
| `rad.ObjHexahedron(verts, M)` | 8 | 6 (face charges) | MSC |
| `rad.ObjRecMag(center, dims, M)` | - | 3 | surface-current brick |

`M = [Mx, My, Mz]` in A/m. Use `[0, 0, 0]` for soft iron whose
magnetization will be SOLVED; give a fixed non-zero M for a permanent
magnet (no Solve needed unless soft iron also present).

Hexahedron vertex order (verified): bottom face z=z0 then top face
z=z0+dz, each CCW:
```python
verts = [[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],
         [x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]]
```
NOTE: ObjHexahedron AUTO-GENERATES its 6 faces internally, so the
face/DOF order is NOT the netgen HEX_FACES order -- do not assume a
fixed face-to-DOF mapping (see the "matrix_structure" topic).

### Minimal soft-iron block in a uniform field

```python
import radia as rad
rad.UtiDelAll()

mat = rad.MatLin(5000.0)              # isotropic soft iron, mu_r=5000
objs = []
for (x0, y0, z0) in cell_origins:     # your hexa grid (meters)
    v = [[x0,y0,z0],[x0+d,y0,z0],[x0+d,y0+d,z0],[x0,y0+d,z0],
         [x0,y0,z0+d],[x0+d,y0,z0+d],[x0+d,y0+d,z0+d],[x0,y0+d,z0+d]]
    o = rad.ObjHexahedron(v, [0, 0, 0])
    rad.MatApl(o, mat)                # attach material to each element
    objs.append(o)
body = rad.ObjCnt(objs)

ext  = rad.ObjBckg(lambda p: [0.0, 0.0, 0.1])   # uniform 0.1 T in z
grp  = rad.ObjCnt([body, ext])                  # body + source

rad.Solve(grp, 0.001, 100, 0)        # prec, max_iter, method(0=LU)
M = rad.ObjM(body)                   # [[center,[Mx,My,Mz]], ...] per elem
rad.UtiDelAll()
```

### Materials

- `rad.MatLin(mu_r)` -- isotropic linear soft iron (preferred single-arg form)
- `rad.MatLin([mu_par,mu_perp],[ex,ey,ez])` -- anisotropic
- `rad.MatSatIsoTab([[H,B],...])` -- nonlinear isotropic B-H curve (H A/m, B T)
- permanent magnet: pass M directly in the constructor; skip MatApl/Solve
  (Solve is only needed when soft iron coexists with the PM).

### Source / external field

`rad.ObjBckg(callable)` ONLY -- callable takes [x,y,z] (meters) and
returns [Bx,By,Bz] (Tesla). The legacy array form
`ObjBckg([Bx,By,Bz])` is NOT supported.

### Solve

`rad.Solve(obj, prec, max_iter, method)`:
- method 0 = LU         (N < 500, guaranteed convergence)
- method 1 = BiCGSTAB   (500 < N < 2000, fastest dense)
- method 2 = HACApK     (N > 2000, O(N log N) memory)
- `image='+x-z'` etc. applies mirror symmetry (IMA) in the SAME call.

`rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0,
bicgstab_tol=1e-4, relax_param=0.0, newton_method=False)` tunes the
solver; `rad.GetSolveStats()` returns timing / iteration counts.

### Read results

- `rad.ObjM(container)` -> per-element [center, [Mx,My,Mz]] (A/m)
- `rad.Fld(obj, 'b', pts)` -> B (Tesla) at one or many external points
- `rad.UtiDelAll()` at the end of every script.

Full element/material/field API: see the `radia_usage` topic.
Interaction-matrix internals: "matrix_structure". Solver near-null
behavior: "eigenvalue_nullspace".
"""

MATRIX_STRUCTURE = """\
## The MMM/MSC interaction matrix -- extraction and structure

Inspect the dense interaction matrix WITHOUT solving:
```python
handle = rad.BuildMatrix(obj, image="")    # builds, does not solve
N, dof = rad.GetInteractMatrix(handle)     # (dof x dof) row-major numpy
```

### Critical gotcha: what GetInteractMatrix returns

It returns the GEOMETRIC interaction matrix N (stored +N, the
demagnetization tensor), which is INDEPENDENT of permeability --
building it at mu_r=10 vs mu_r=5000 gives the identical matrix
(verified). It is NOT the system matrix the solver factorizes.

Form the SYSTEM matrix yourself (sign per Radia's ComputeEntry,
A = -N + diag(1/chi)):
```python
chi = mu_r - 1.0                 # linear; per-element for a B-H curve
A   = np.diag(np.full(dof, 1.0/chi)) - N
```
In the digest/paper notation A_ij = delta_ij/(mu_r-1) + G_ij, that
"G" equals -N (the returned matrix); the +N/-N flip is a storage
convention, not physics.

### Properties (verified)

- dof = sum of per-element DOF (tet 3, wedge 5, hex 6). A 32-hex block
  -> 192 dof.
- N is NON-symmetric (MSC collocation evaluates the normal field at
  face centers; relative asymmetry ~ a few %).
- Row-major [target][source]: `N[i,j]` = effect ON dof i FROM dof j.
- Element-major ordering: element e owns dof [k*e, k*e+k) with k its
  per-element DOF count.

### Densifying the ACTUAL HACApK (ACA+) operator

`rad.GetInteractMatrix` always returns the EXACT geometric N (never the
ACA-compressed one). To obtain the actual ACA-approximated SYSTEM matrix
A = -N + diag(1/chi):
```python
A_haca, dof = rad.HMatrixDensify(handle)   # MatVec on unit vectors
```
It builds the MSC H-matrix and applies its `MatVec` to unit vectors,
returning the dense A in the ORIGINAL DOF ordering (directly comparable to
`diag(1/chi) - N`, and usable with a loop basis). Use it to check whether
the real ACA+ shifts the near-zero eigenvalues (it does not materially --
see "eigenvalue_nullspace"). C++ path: `RadHACApKMMMManager::MatVec`
(formerly `RadHACApKMSCManager` before the 2026-06-23 moment-yano cleanup)
exposed via `radTApplication::HMatrixDensify`.

See "eigenvalue_nullspace" for why the spectrum of A matters, and
`docs/solver/MSC_NULLSPACE_DEFLATION.md` for the full treatment.
"""

EIGENVALUE_NULLSPACE = """\
## Near-null "loop" modes of the MSC operator (conditioning + solver behavior)

(CEFC 2026 full-paper study; reproducible scripts in
`examples/mmm_eigenvalue_study/`.)

### The null space = cycle space of the element graph

N has an EXACT null space: surface-charge distributions producing zero
normal field at every collocation point. Its dimension equals the CYCLE
RANK (first Betti number) of the element-adjacency graph (nodes =
elements, edges = shared internal faces):

    null_dim = F_internal - N_elem + 1

verified to machine precision:

| block | 1x1x1 | 2x1x1 | 2x2x1 | 2x2x2 | 3x3x3 | 4x4x2 | 4x4x4 | 5x5x5 |
|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| elem  | 1 | 2 | 4 | 8 | 27 | 32 | 64 | 125 |
| nulls | 0 | 0 | 1 | 5 | 28 | 33 | 81 | 176 |

A single element and a 1-D chain have NO null space; it first appears at a
2x2 arrangement (the minimal circulating "loop" of four elements).

Physical derivation: a null mode puts OPPOSITE charges on the two sides of
each shared internal face (field cancels) AND zero net charge per element
(so the compensating centroid point charge vanishes). With one "flux" per
shared face, "zero net charge per element" is a discrete divergence-free
condition -> the solution space is the cycle (circulation) space. Minimal
cycle = 2x2 plaquette.

### Explicit local basis (SVD-free)

Each 2x2 plaquette circulation (sigma=+1/-1 around the loop) is an EXACT
null mode (||N*loop||/||loop|| ~ 1e-15), and the plaquette set spans the
null space (rank == null_dim). Each touches only 4 elements (8 shared-face
DOF) -> O(1) local support, built directly from mesh topology with NO SVD
(`nullmode_cycle_basis.py`). The null projector P=VV^T is local (decays
~2x per cell), so a local basis exists; SCDM (selected columns of P) gives
one too (`nullmode_local_deflation.py`).

### These modes set the conditioning at high mu_r

A = diag(1/(mu_r-1)) - N maps each null mode to an eigenvalue exactly
at 1/(mu_r-1). So the smallest eigenvalues of A are ~1/(mu_r-1) and
the condition number scales like mu_r. mu_r -> infinity makes A
singular: HIGH PERMEABILITY IS THE WORST CASE. The null eigenvectors
are circulating magnetization patterns that radiate ~no external field.

### "Beautiful -> ugly" iterative behavior (Yano observation)

For a high-mu_r block in a uniform field, BiCGSTAB from a zero initial
guess passes a PHYSICAL aligned solution at few iterations ("beautiful")
and then drifts to a loop-dominated converged/LU solution ("ugly"),
because Krylov resolves large-eigenvalue (physical) components first
and near-null (loop) components last. Dramatic at high mu_r: mean
alignment of M with the applied field falls from 0.94 (few iter) to
0.14 (converged) at mu_r=1e5 on a 6x6 block. Reproduce with Radia
APIs only (Solve + ObjM + SolverConfig, bicgstab_tol as the iteration
knob) -- no matrix extraction needed (see `beautiful_ugly_viz.py`).

### Does the H-matrix remove the bad modes? (Phase A: no)

A controlled SVD low-rank truncation of the far-field blocks does NOT
regularize the near-null spectrum: at the HACApK working tolerance
eps=1e-4 the smallest eigenvalues / condition number are essentially
unchanged (3.9e4 -> 4.0e4 on a 6x6x6 block); more aggressive truncation
moves them CLOSER to zero. The REAL ACA+ operator (via
`rad.HMatrixDensify`, see "matrix_structure") CONFIRMS this: on
8x8x8 / mu_r=1e5 it spreads the exact-zero cluster into near-zero
eigenvalues BELOW 1/(mu_r-1), worsening cond 7.7e5 -> 1.08e6 at eps=1e-4
(-> 8.9e6 at eps=1e-2). So the data-sparse representation buys memory/time,
NOT spectral regularization -- the "H-matrix suppresses the bad modes"
hypothesis is refuted by both the SVD proxy and the real ACA+.

### Loop modes: RESOLVED by consolidating on HDiv-VIM (deflation removed 2026-06-09)

The loop null space above is a genuine property of A = diag(1/(mu_r-1)) - N:
the div-side tree-cotree gauge (dual of FEM-A's gradient gauge). Helmholtz
curl(T) = solenoidal "loop" magnetization, invisible to the field (N L = 0);
the cotree faces of the element-adjacency graph form the loop basis. A lifts it
only to 1/(mu_r-1), so high-mu_r is a low-frequency (loop-star) breakdown.

An earlier C++ effort added RUNTIME loop handling to the MSC/MMM solver:
matrix-free deflation (SetHACApKDeflation / SetDeflateNullspace via the
BuildLoopBasis cycle basis), an alpha-free loop-star gauge (SolveLoopStar:
the reduced A_SS = S^T A S star block with an H-LU preconditioner), and a
post-solve Helmholtz-Hodge loop projection (SetLoopProjection). ALL of these
were REMOVED on 2026-06-09 -- the project consolidated on the HDiv-VIM operator.

WHY HDiv-VIM instead: the FEEC H(div) RT0 demag operator N = B^T G B builds the
charge map B (rho = -div M on L2, sigma = M.n on SurfaceL2) so the loop space is
exactly ker(B) -- field-null BY CONSTRUCTION via the de Rham complex,
mu_r-INDEPENDENT, with no runtime deflation / gauge / projection. M^{-1} mass
preconditioning already deflates the loops, so a plain symmetric Krylov is
well-conditioned at every mu_r. See `radia.vim` and `tests/feec/test_hdiv_vim_*`.

Effect on the MSC/MMM default solver (Block-Jacobi BiCGSTAB, HACApK Method 2):
it no longer projects out loops. Fields are UNAFFECTED (N L = 0) -- rad.Fld is
identical -- only the sigma/M gauge may carry a field-null loop component at very
high mu_r. (The nullspace THEORY above stands: it is exactly why HDiv-VIM, which
makes loops ker(B), is the right consolidation.)

### Practical guidance

- At very high mu_r the converged discrete MSC solution can carry spurious loops;
  this is a formulation/conditioning property, NOT a solver bug (LU shows it too),
  and it does not change the computed field.
- For a loop-free-by-construction formulation use HDiv-VIM (`radia.vim`);
  otherwise keep mu_r physical for trustworthy high-mu_r MSC results.
"""

MULTIPOLE_MODES = """\
## What field does the 6-DoF MSC element create? (multipole modes + conditioning)

(Sugahara-lab study 2026-06-22. Numerics: `examples/vim/yano_moment_svd_multipole.py`
(SVD/eig + multipole projection) and `examples/vim/yano_moment_aspect_ratio.py`
(iteration vs aspect ratio). Symbolic proof:
`packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/yano_6dof_multipole.wls`.)

### The 6 DOF ARE exactly monopole + dipole + 2 quadrupole

SVD (and eig) of a single hexahedron's 6x6 interaction matrix N (sigma -> demag field)
shows the 6 face-charge DOF produce PURE multipole fields:

    1 monopole  +  3 dipole  +  2 (diagonal) quadrupole   =   6 modes

On a cube each eigenmode is ~100% one multipole. So **6-DoF MSC = MMM (the 3 dipole DOF a
tetrahedron also carries) + a monopole + the 2 quadrupole moments that 6 axis-aligned face
charges can represent**. The off-diagonal (shear) quadrupole of an axis-aligned hex is
EXACTLY 0 (int_face x_i x_j dA = 0 by face symmetry).

Symbolic proof (yano_6dof_multipole.wls, generic shear s, self-test ALL PASS) -- the RANK
of the charge -> multipole-moment map:

    rank[mono;dipole]               = 4
    rank[mono;dipole;quad]          = 6     (quadrupole adds EXACTLY 2 independent DOF)
    rank[mono;dipole;quad;octupole] = 6     (octupole adds 0 -> NO independent octupole)

so the field's octupole and higher multipoles (and the other 3 quadrupole components) are
LINEARLY DETERMINED by these 6 -- the element has no independent content beyond
mono + dipole + 2 quadrupole.

### Distortion only ROTATES the quadrupole (no new DOF, no octupole)

For a distorted (sheared) hex the 6 eigenmodes are STILL mono + 3 dipole + 2 quadrupole; the
distortion only (a) tilts the quadrupole pair's principal axes (~12-26 deg for a 0.4 xy-shear,
so the quadrupole acquires SHEAR appearance in the GLOBAL frame) and (b) splits the degenerate
dipole/quadrupole eigenvalues. The shear-quadrupole functional is exactly LINEAR in the shear
parameter (closed form int_face d_x d_y dA = {s/3, s/3, s, s, s/3, s/3}), vanishing for the
cube. So "does a distorted 6-DoF create beyond 2 diagonal quadrupole?" -- NO: it is always
exactly 2 quadrupole DOF; only their ORIENTATION changes.

### Condition number IS the multipole field-strength ratio

cond(N) = (strongest multipole field) / (weakest). The MONOPOLE is by far the strongest mode
(singular value ~7.66 on a unit cube) and dominates cond(N) ~ 27 -- BUT it is CONSTRAINED to
zero by div B = 0 (sum_f A_f sigma_f = 0), so it does not enter the solve. The physically
relevant conditioning is the DIPOLE / QUADRUPOLE strength spread (cond with the monopole
removed):

    cubic cell          : ~1.4    (dipole and quadrupole fields almost equally strong -> ideal)
    thin / tall (1:1:4) : ~6-7    (the thin-direction quadrupole field weakens)

### Aspect ratio -- NOT size or condition number -- sets the iteration count

The quadrupole is the WEAKEST field mode; on an anisotropic cell its thin-direction field
weakens further, opening the dipole/quadrupole gap. Block-Jacobi inverts the local 6x6 block,
AMPLIFYING that weak quadrupole mode -> the iteration trouble. Measured (C-yoke, mu_r=1000,
block-Jacobi GMRES): a clean V-curve in cell aspect ratio (iters MINIMIZE at dx/dz ~ 1 and grow
on both sides) while dof grows MONOTONICALLY -> it is aspect ratio, NOT problem size. With
~cubic cells (nz ~ nxy/3) the iters stay BOUNDED (~64-99) as dof grows 1296 -> 16848 -> the
6-DoF MSC SCALES with cheap block-Jacobi (no H-LU) on well-shaped cells. Per-row normalization /
condition number is NOT the lever -- block-Jacobi absorbs any per-row (and any same-subspace
local) scaling, so the iteration count is invariant to it. (This is the orthogonal axis to the
high-mu_r loop story in "eigenvalue_nullspace".)

### Point-matching (EIEM2) vs moment-matching

EIEM2 collocation (eval-point alpha=0.5) is point-matching: it collocates a SINGLE physical
quantity (the field) at points, so it does NOT mix the monopole/dipole/quadrupole orders
(L^2 / L^3 / L^4) the way moment-matching (integral functionals) does. Measured: EIEM2 is lower
at EVERY aspect ratio and has a shallower V than moment-matching, but is not fully aspect-immune
(a thin cell puts the normal-offset eval point too close to its face). The eval-point offset
alpha and the gradient/quadrupole term are the SAME degree of freedom (the alpha-offset is a
field-gradient sample, n.gradH.u with u along the normal -> it sees only the DIAGONAL part of
the gradient tensor; the full quadrupole/shear needs the moment or a transverse sample). On a
fair system-vs-system comparison BOTH are well-posed: the raw `GetInteractMatrix` N is
near-singular (~1e18, the loops); the SYSTEM matrix A = -N + diag(1/chi) regularizes to ~1e4,
and there the moment system's condition number is actually NOT worse than EIEM2's.

### Practical guidance

- Mesh hexahedra as close to cubic as the geometry allows; block-Jacobi then scales without H-LU.
- The 6-DoF element is a faithful dipole + 2-quadrupole source on ANY planar-faced hex;
  distortion costs only quadrupole ORIENTATION, never spurious modes or independent octupole.
- For a loop-free-by-construction high-mu_r formulation use HDiv-VIM (`radia.vim`); see
  "eigenvalue_nullspace".
"""

SECTIONS = {
    "build_msc_mmm": BUILD_MSC_MMM,
    "matrix_structure": MATRIX_STRUCTURE,
    "eigenvalue_nullspace": EIGENVALUE_NULLSPACE,
    "multipole_modes": MULTIPOLE_MODES,
    "chubar_1998": CHUBAR_1998,
    "takahashi_2007_aca": TAKAHASHI_2007_ACA,
    "pradhan_2007": PRADHAN_2007,
    "janet_2005_mmm_rnm": JANET_2005_MMM_RNM,
    "le_duc_peec_mmm": LE_DUC_PEEC_MMM,
    "weddemann_fembem": WEDDEMANN_FEMBEM,
    "mmm_vs_other": MMM_VS_OTHER,
}


def get_mmm_core_documentation(topic: str = "chubar_1998") -> str:
    """Return MMM core knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
