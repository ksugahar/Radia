"""
Induction heating simulation knowledge base for the Radia MCP server.

Covers: complete EM-to-thermal workflow for induction heating using NGSolve,
including Gmsh mesh loading, A-Phi eddy current formulation, Joule heat
computation, transient thermal analysis, rotating workpiece, GMSH output,
and post-processing patterns.

Public docs/notebooks:
  - docs/induction_heating/induction_heating_demo_showcase.ipynb -- executed
    ESIM/Bessel public showcase with embedded figures.
  - validation_test/ih_esim_benchmark/ -- nonlinear ESIM benchmark scripts,
    checked JSON, and publication figures.

Sources:
  - https://docu.ngsolve.org/latest/i-tutorials/
  - https://forum.ngsolve.org/

Student element-Zs workflow: Geometry Update accepts panel_zs_file, an absolute
mesh/frequency-bound JSON path (ohm, BND-element-order). Use
radia.simulink.writeIHPanelImpedance with a full ordered complex vector; centroid
selection is performed by the caller. P1 weak coupling, genus-0 or one supported
flux-linked handle with dense P1 BIE and flat-triangle loop work, constant during each
linear simulation, with content-hash-triggered rebuild. Independent FEM accuracy
is separate: independent_fem_validation.status starts as not-performed.
See docs/induction_heating/SIMULINK_ELEMENT_ZS.md.
"""

from radia_mcp.matrix_solvers.direct_solvers_knowledge import SOLVER_CAPACITY

INDUCTION_HEATING_OVERVIEW = """
# Induction Heating Simulation with NGSolve

## Physics Overview

### Linear numerical acceptance scope

The repository-owned linear acceptance is complete for the declared 1 kHz
copper solid/bored axisymmetric cylinders, prescribed circular current and
uniform scalar SIBC. It compares independent axisymmetric A-form FEM-SIBC
against P1 surface BEM (including the full loop field for the bore).
Worst errors are 0.691% surface loss, 1.179% reaction power and 1.558%
complex tangential-field relative L2 on the declared side-surface sampling
stencil, all below the IH-specific 2% limit. Exterior/mesh refinement and
conservative heat projection also pass. This is not a whole-surface field
error certificate or a general material/frequency/geometry envelope.

Evidence: validation_test/eddy_current_analytical_validation/
ih_sibc_closure_20260915.json and ih_linear_acceptance_audit_20260916.json.
The latter rechecks saved measurements and recovered source bytes (allowing
only CRLF/LF conversion); it is not a new solve or installed-binary acceptance.
Do not present this as validation of SIBC approximation error against a
volume-resolved conductor, nonlinear steel, temperature-dependent BH, or
complete real-file Simulink delivery. Nonlinear experimental validation is a
separate activity; known nonlinear implementation issues remain issues.

### Native Simulink preview limits

The geometry assembler currently freezes the full EM unit-current loss field
and scales it by peak current squared; it emits one Eddy unknown and rotation=none.
Sinusoidal RMS current must be converted to peak amplitude before this interface.
The native thermal preview is linear and does not implement the Python nonlinear
enthalpy, latent-heat and temperature-dependent material model.
This is not temperature-dependent nonlinear BH coupling. Hand-authored dense
temperature-slope operators are a small-system runtime capability, not a
validated production generator. heat_projection represents diagonal quadratic
loss only: general overlapping FEM/BEM coefficients require cross terms.

For periodic sample configurations, Thermal state and heat are in workpiece
material coordinates. Eddy alone maps temperature to the source frame and heat
back; Thermal must not rotate temperature again. Linear fractional interpolation
can smooth hotspots; weighted-integral conservation does not imply peak accuracy.
Unsafe nearly cancelling integral corrections fail fast. High-order FE coefficient
temperature rotation remains unsupported. Use standard NGSolve H1 for heat.

Induction heating uses alternating current in a coil to generate eddy currents
in a conductive workpiece. The eddy currents produce Joule heating (I^2 R losses).

**Three-Phase Simulation Chain:**

```
Phase 1: Electromagnetic (Eddy Current) Analysis
  ├─ A-Phi formulation (vector potential + scalar potential)
  ├─ Frequency-domain complex solve at f = 8 kHz (typical)
  └─ Output: B, E, J fields

Phase 2: Heat Source Computation
  └─ Q = 0.5 * sigma * |E|^2 [W/m^3] (time-averaged Joule heating)

Phase 3: Transient Thermal Analysis
  ├─ Unsteady heat equation with Joule heat source
  ├─ theta-scheme (implicit) time stepping
  └─ Output: Temperature evolution T(x, t)
```

## Typical Parameters (Industrial Induction Heating)

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Frequency | f | 8,000 | Hz |
| Coil current | I | 7,000 | A |
| Cu conductivity | sigma_Cu | 57e6 | S/m |
| Work conductivity | sigma_Fe | 10e6 | S/m |
| Work relative permeability | mu_r | 1,000 | - |
| Skin depth (work) | delta | 0.18-0.56 | mm |
| Work density | rho | 7,800 | kg/m^3 |
| Work specific heat | c | 467 | J/(kg*K) |
| Work thermal conductivity | kappa | 46.6 | W/(m*K) |
| Convection coefficient | h | 10 | W/(m^2*K) |
| Ambient temperature | T_ext | 20 | deg C |

## Skin Depth Formula

```python
import numpy as np
mu0 = 4e-7 * np.pi
f = 8000          # Hz
sigma = 10e6      # S/m (work material)
mu_r = 1000       # relative permeability
omega = 2 * np.pi * f

delta = np.sqrt(2 / (omega * mu0 * mu_r * sigma))
# delta ~ 0.18 mm for iron at 8 kHz
```

The mesh boundary layer near the work surface must resolve the skin depth
(typically 2-10 layers with growth ratio 1.2-1.6).
"""

INDUCTION_HEATING_GMSH_MESH = """
# Gmsh Mesh Loading for Induction Heating

## Loading Meshes (.vol files)

NGSolve loads meshes from .vol files (exported from Cubit or Netgen).

```python
from ngsolve import Mesh

# Load .vol mesh (exported from Cubit or Netgen)
mesh = Mesh("toymodel.vol")

# Verify material regions and boundaries
print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())
print("Elements:", mesh.ne)
print("Vertices:", mesh.nv)
```

## Physical Groups (Material Regions and Boundaries)

The checked solver `.vol` defines NGSolve material regions and boundaries.
Gmsh 5 / MSH 4.2 physical groups label the post-processing output; they are
not a solver-mesh import route. Preserve the solver region names in exported
physical groups when the exporter supports the selected scope.

### Typical Physical Group Setup for Induction Heating

**Volume groups (materials):**

| Physical Group | NGSolve Material Name | Description |
|----------------|----------------------|-------------|
| work_main | work_main | Workpiece core (iron/steel) |
| work_skin | work_skin | Workpiece skin (boundary layer) |
| coil_main | coil_main | Coil conductor core (copper) |
| coil_skin | coil_skin | Coil skin (boundary layer) |
| Sair | Sair | Shielding/surrounding air |
| air | air | Far-field air |

**Surface groups (boundaries):**

| Physical Group | NGSolve Boundary Name | BC Type |
|----------------|----------------------|---------|
| in_tri / in_quad | in_tri / in_quad | Coil input (phi=1) |
| out_tri / out_quad | out_tri / out_quad | Coil output (phi=0) |
| work_surface | work_surface | Convection BC (thermal) |
| coil_surface | coil_surface | Coil surface |
| outer | outer | Far-field (A=0) |

### Mesh from Cubit (export vol)

The recommended path is `export vol` (Cubit plugin command),
which writes the .vol directly with proper FaceDescriptors, materials,
and (optionally) Kelvin Periodic identifications. The legacy
`Cubit2gmsh.py` workflow is obsolete.

```python
# Cubit-side: build geometry, mesh, export
import cubit, math, os
cubit.init(["cubit", "-nographics", "-nojournal"])

# Coil: sweep (NOT webcut — see policy below)
cubit.cmd('create vertex 0.030 0 0')
cubit.cmd('create vertex 0.033 0 0')
cubit.cmd('create vertex 0.030 0 0.003')
cubit.cmd('create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full')
cubit.cmd('create surface curve 1')
cubit.cmd('sweep surface 1 axis 0 0 0 0 0 1 angle 355')
coil_vid = cubit.get_last_id("volume")

# Air sphere
cubit.cmd('create sphere radius 0.120')
air0 = cubit.get_last_id("volume")

# CRITICAL: subtract before imprint+merge for nested volumes
# (imprint+merge alone does NOT share surfaces between fully nested
#  volumes; the resulting .vol has wrong FaceDescriptors)
cubit.cmd(f'subtract volume {coil_vid} from volume {air0} keep_tool')
all_vols = list(cubit.parse_cubit_list("volume", "all"))
air_vid = [v for v in all_vols if v != coil_vid][0]
cubit.cmd(f'imprint volume {coil_vid} {air_vid}')
cubit.cmd(f'merge volume {coil_vid} {air_vid}')
# Verify: "Consolidated N pairs of surfaces" with N > 0

# Mesh
cubit.cmd(f'volume {coil_vid} scheme tetmesh')
cubit.cmd(f'volume {coil_vid} size 0.0024')
cubit.cmd(f'mesh volume {coil_vid}')
cubit.cmd(f'volume {air_vid} scheme tetmesh')
cubit.cmd(f'volume {air_vid} size 0.030')
cubit.cmd(f'mesh volume {air_vid}')

# Blocks (named, NOT by hardcoded ID)
cubit.cmd(f'block 1 add volume {coil_vid}')
cubit.cmd('block 1 name "coil"')
cubit.cmd(f'block 2 add volume {air_vid}')
cubit.cmd('block 2 name "air"')

# Source/sink: identify by area (not by hardcoded surface ID).
# A_circle = pi * a_coil^2 ≈ 2.83e-5 for a=3mm.
A_gap = math.pi * 0.003**2
gap_faces = []
for sid in cubit.parse_cubit_list("surface", f"in volume {coil_vid}"):
    a = cubit.surface(sid).area()
    if a < 1.1 * A_gap:
        gap_faces.append((sid, cubit.get_center_point("surface", sid)[1]))
gap_faces.sort(key=lambda t: t[1])
cubit.cmd(f'sideset 1 add surface {gap_faces[0][0]}')
cubit.cmd('sideset 1 name "source"')
cubit.cmd(f'sideset 2 add surface {gap_faces[-1][0]}')
cubit.cmd('sideset 2 name "sink"')

# Export
out = os.path.abspath("model.vol")
cubit.cmd(f'export vol "{out}" order 2 overwrite')
```

Then on the NGSolve side: `mesh = Mesh("model.vol")`.

**Cubit policy reminders** (full details in cubit_scripting_knowledge):
- Use **sweep** to build coils, NOT `create torus + webcut + delete`
  (webcut+delete can silently produce a half coil)
- Use **subtract first, then imprint+merge** for nested volumes
- **Identify entities by geometric properties** (area, centroid),
  NOT by hardcoded IDs — IDs change after imprint/merge/version

## Naming convention contract (.jou ↔ calc scripts)

**POLICY**: Once the .vol is exported, NO label is passed at the GUI or
calc-script level. The .jou file is the single point where labels are
defined and verified by the user.

| Entity | Required name | Purpose |
|--------|---------------|---------|
| sideset on coil terminal (current injection)  | `source` | calc_inductance.py default `--source source` |
| sideset on coil terminal (current extraction) | `sink`   | calc_inductance.py default `--sink sink` |
| block of coil volume                          | `coil`   | calc_inductance.py default `--coil coil` (filters BEM to coil surfaces) |
| block of workpiece volume (optional)          | `workpiece` | calc_inductance.py `--workpiece workpiece` |
| block of air volume (optional)                | `air`    | (computation skips it) |

If the .jou follows this convention, the student opens the GUI, picks
the .vol file, hits Run — no label input needed. Errors at this stage
mean the .jou does not follow convention; fix the .jou, not the GUI.

## Curve Order vs FES Order (CRITICAL for IH accuracy)

**CRITICAL**: Curve a mesh only while its generating geometry is live.  Both
`fem_esim_kelvin.py` (2D axi) and `fem_esim_3d.py` (3D) should call
`mesh.Curve(>=2)` on the mesh freshly generated from their in-memory
SplineGeometry/OCCGeometry when the geometry contains circles, a torus, a
sphere, or a cylinder.  Missing that construction-time curving causes a
**systematic L under-prediction of ~10 %**, because curved entities collapse to
polygons.

Never call `mesh.Curve()` merely because the FES order was raised after loading
an arbitrary `.vol`.  The imported file already owns its geometry order; its
CAD association may be absent or incomplete.  Post-load `Curve()` rebuilds the
curved nodes from whatever geometry Netgen associates with the mesh: a CAD-less
export is flattened in a fresh process, is curved onto another mesh's
process-global geometry once such a mesh was loaded, and a mesh rescaled with
`ngmesh.Scale` after meshing is projected onto its unscaled CAD.  When a higher
order is unavoidable, use `radia.mesh_curve.ensure_curve_order`, which checks
the embedded CAD before curving.  To change geometry
order, curve in the originating mesher and re-export the `.vol`.  To change only
the approximation order, use `H1(mesh, order=p)` or `HCurl(mesh, order=p)` and
leave the loaded mesh unchanged.  This is guarded repository-wide by the HIGH
lint rule `ngsolve-curve-after-vol-import`.

### Practical recommendation (Sugahara, 2026-04-14)

General rule of thumb: **`curve_order=2` with `FES order=1 or 2`
is sufficient** for most IH problems. Going higher rarely pays off.

**Exception for strongly-curved 3D geometries** (torus coil + sphere
Kelvin + cylinder workpiece, such as `ih_fem_kelvin_sample.jou`):
the default `curve_order=2` with HCurl order=1 gives L error ~10 %
vs 2D axisym full-res reference, because the sphere and torus need
finer geometry than Curve(2) can capture on a maxh=15 mm mesh.

**Tested at Cu 7 kHz** (coil R=30 mm, a=3 mm, WP cylinder R=25 mm):

| 3D setting                  | L vs full-res | P vs full-res | t_solve |
|-----------------------------|---------------|---------------|---------|
| HCurl order=1 + Curve(2)    | -9.91 %       | -4.0 %        | 23 s    |
| HCurl order=1 + Curve(3)    | **-5.37 %**   | **-0.46 %**   | 26 s    |
| HCurl order=2 + Curve(2)    | +3.2 % (but spurious P -26 %) | | 1371 s  |

**Bottom line for geometry-construction code** (`fem_esim_3d.py`,
`calc_fem_kelvin.py` when they still own the live geometry):

- Start with `HCurl order=1 + mesh.Curve(max(order+1, 2))`
- If L seems too low (-5 % or worse) and P looks OK, **raise Curve
  to 3** instead of raising FES order. Curve is cheap; FES doubles
  ndof and can introduce spurious modes.
- HCurl order=2 costs ~60x more on this test case (ndof 83k -> 333k,
  23s -> 23min) and can OVERSHOOT L and WORSEN P at the same Curve
  level. Don't raise FES order without raising Curve too.

### Why 2D is less sensitive than 3D

2D axisym has only planar curves (coil circle + Kelvin arc). H1 scalar
elements curved at mesh-construction time with order 2 or 3 match the geometry
accurately even at moderate mesh resolution. 3D adds sphere + torus in space --
strong Gaussian curvature that needs better surface approximation.  This does
not authorize re-curving a loaded `.vol`.

### Companion fix (2026-04-14)

`fem_esim_kelvin.build_geometry` previously NEVER called `mesh.Curve()`
despite using `H1 order=3` by default. Adding `mesh.Curve(order)`
improved L by +7 nH (-11.7 % -> -4.6 %) on the Cu 7 kHz test case.
Users running old checkouts must pull this fix to get trustworthy L.

## Boundary Layer Meshing for Skin Depth

The skin depth at the work surface must be resolved by the mesh.
Cubit/Gmsh boundary layers provide graded elements near the surface.

```python
# Skin depth parameters
delta = 0.18e-3  # m (skin depth at 8 kHz for iron, mu_r=1000)

# Boundary layer mesh parameters
n_layers = 6           # Number of layers
growth_ratio = 1.2     # Element growth ratio
first_layer = delta/3  # First layer thickness ~ delta/3
```

### Multi-Parameter Skin Depth Study

Different configurations for parameter sweeps:

| Skin Depth (mm) | Layers | Growth | Application |
|------------------|--------|--------|-------------|
| 0.060 | 10 | 1.2 | High mu_r, high freq |
| 0.180 | 4-6 | 1.2 | Standard iron at 8 kHz |
| 0.246 | 2-4 | 1.6 | Lower mu_r material |
"""

INDUCTION_HEATING_EDDY_CURRENT = """
# A-Phi Eddy Current Formulation for Induction Heating

## Complete Production Implementation

This is a practical, production-quality A-Phi implementation for induction
heating with Gmsh meshes, based on real simulation workflows.

```python
from ngsolve import *
import numpy as np

# ============================================================
# 1. Load mesh (.vol)
# ============================================================
mesh = Mesh("toymodel.vol")

print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())

# ============================================================
# 2. Physical parameters
# ============================================================
mu0 = 4e-7 * np.pi
freq = 8000  # Hz
s = 2j * np.pi * freq  # Complex angular frequency

# Conductivity [S/m]
sigma_d = {
    "work_main": 10e6, "work_skin": 10e6,     # Iron/steel
    "coil_main": 57e6, "coil_skin": 57e6,     # Copper
    "Sair": 0, "air": 0,                       # Air (non-conductive)
}
sigma = CoefficientFunction(
    [sigma_d.get(mat, 0) for mat in mesh.GetMaterials()]
)

# Relative permeability
mur_d = {
    "work_main": 1000, "work_skin": 1000,     # Iron/steel
    "coil_main": 1, "coil_skin": 1,           # Copper
    "Sair": 1, "air": 1,                       # Air
}
nu = CoefficientFunction(
    [1 / (mu0 * mur_d.get(mat, 1)) for mat in mesh.GetMaterials()]
)

# ============================================================
# 3. Finite element spaces (A-Phi product space)
# ============================================================
order = 1  # Polynomial order

# HCurl for magnetic vector potential A
#   - nograds=True: removes gradient null space (gauge via Phi coupling)
#   - complex=True: frequency-domain analysis
#   - dirichlet="outer": A=0 on far-field boundary
fesA = HCurl(mesh, order=order, nograds=True,
             dirichlet="outer", complex=True)

# H1 for scalar electric potential Phi (conductor regions only)
#   - definedon: Phi exists only where sigma > 0
#   - dirichlet: fixed potential on coil terminals
#   - dirichlet_bbbnd: edge Dirichlet for 3D problems
fesPhi = H1(mesh, order=order,
            definedon="coil_main|coil_skin",
            dirichlet="in_tri|in_quad|out_tri|out_quad",
            dirichlet_bbbnd="in_tri|in_quad|out_tri|out_quad",
            complex=True)

# Product space
fesAPhi = fesA * fesPhi
(A, phi), (N, psi) = fesAPhi.TnT()

print(f"Total DOFs: {fesAPhi.ndof:,}")

# ============================================================
# 4. Dirichlet boundary conditions for Phi
# ============================================================
gfAPhi = GridFunction(fesAPhi)
gfA, gfPhi = gfAPhi.components

# Set coil terminal potentials: phi=1 at input, phi=0 at output
gfPhi.Set(1, definedon=mesh.Boundaries("in_tri|in_quad"))
gfPhi.Set(0, definedon=mesh.Boundaries("out_tri|out_quad"))

# ============================================================
# 5. Bilinear form (weak form)
# ============================================================
a = BilinearForm(fesAPhi)

# Curl-curl term (magnetic energy): integral(nu * curl(A) . curl(N))
a += nu * curl(A) * curl(N) * dx

# Stabilization term (gauge regularization)
a += 1e-6 * nu * A * N * dx

# Eddy current term in coil: s*sigma*(A+grad(phi))*(N+grad(psi))
a += s * sigma * (A + grad(phi)) * (N + grad(psi)) * dx("coil_main|coil_skin")

# Eddy current term in work: s*sigma*A*N (no phi in work)
a += s * sigma * A * N * dx("work_main|work_skin")

# ============================================================
# 6. Preconditioner and assembly
# ============================================================
# BDDC preconditioner: register BEFORE assembly
c = Preconditioner(a, "bddc")

with TaskManager():
    a.Assemble()
    c.Update()

# ============================================================
# 7. Solve with preconditioned CG or GMRes
# ============================================================
from ngsolve.krylovspace import CGSolver, GMRes

# RHS from Dirichlet BC
r = a.mat.CreateColVector()
r.data = -a.mat * gfAPhi.vec

with TaskManager():
    # Option A: CG solver (if system is close to symmetric)
    solver = CGSolver(mat=a.mat, pre=c.mat,
                      maxiter=1000, tol=1e-10, printrates=True)
    gfAPhi.vec.data += solver * r

    # Option B: GMRes solver (more robust for non-symmetric systems)
    # GMRes(mat=a.mat, pre=c.mat, rhs=r,
    #        sol=gfAPhi.vec, maxsteps=1500, tol=1e-10, printrates=True)

print("Solver converged.")
```

## Current Normalization

In practice, the coil current is normalized to a target value (e.g., 7000 A).
The current is computed from the solution and used as a scaling factor.

```python
# Compute actual current from solution
# J = -s * sigma * (A + grad(phi)) in coil, J = -s * sigma * A in work
# Current = integral of J . grad(psi) over coil cross-section

# Use a dummy psi to measure current
gfA_sol, gfPhi_sol = gfAPhi.components

# Direct current extraction via J . n on coil output face
J_coil = s * sigma * (gfA_sol + grad(gfPhi_sol))
I_out = Integrate(J_coil * grad(psi) * dx("coil_main|coil_skin"), mesh)
print(f"Computed current: {abs(I_out):.1f} A")

# Normalize to target current
I_target = 7000  # A
scale = I_target / abs(I_out)

# Scaled physical fields
B = curl(gfA_sol) * scale           # Magnetic flux density [T]
E = -s * (gfA_sol + grad(gfPhi_sol)) * scale  # Electric field [V/m]
J = sigma * E                        # Current density [A/m^2]
```

## Joule Heat Source Computation

```python
# Time-averaged Joule heating density [W/m^3]
# Q = 0.5 * Re(J . conj(E)) = 0.5 * sigma * |E|^2
Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real

# Alternative formulation using J:
# Q = 0.5 * (1/sigma) * |J|^2  (where sigma > 0)
# Q = 0.5 * (J.real * J.real + J.imag * J.imag) / sigma
```

## Solver Selection

The supported complex HCurl workpiece route uses AMS+COCR at order 1 and
BDDC with an AMS coarse solver at orders 2 and 3. A complex system must not
be treated as a real SPD system merely because its matrix is symmetric.
The compound A-V coil route has separate space and preconditioning requirements;
its supported direct solve explicitly selects SparseCholesky. Scalar thermal
and axisymmetric H1 equations also retain their own solver contracts.
Measure setup, solve, memory and true residual on the intended problem before
making a performance claim; small problems can favour a direct solve.
"""

INDUCTION_HEATING_THERMAL = """
# Transient Thermal Analysis for Induction Heating

## Production path: native IH Simulink blocks

The tracked production model is ``packages/radia-simulink/matlab/radia_ih.slx``.  Its readable Level-2
MATLAB S-Functions expose two independent lifecycle boundaries:

* ``radia.simulink.ihEddySFunction`` receives coil current, workpiece angle,
  and the accepted temperature distribution, then emits distributed heat.
* ``radia.simulink.ihThermalSFunction`` receives distributed heat, ambient
  temperature, and workpiece angle, then emits the accepted temperature field.

Both wrappers own checked ``uint64`` handles to the standalone ``radia_mex``
Eddy and Thermal kernels.  Python is not called per Simulink time step.  The
masked parameter block loads a checked MAT/JSON configuration; geometry-driven
updates first validate each ``.vol`` with ``check-vol`` and then assemble a new
configuration explicitly.

Workflow:

1. Run an EM solve method that emits ``qsurf.sol`` (PEEC+BEM,
   BEM-A+BEM, PEEC+FEM+Kelvin, or FEM-full).
   For the corrected weak BEM route, ``qsurf_method`` must be
   ``solved-total-field-lumped-P1``: heat is evaluated from the solved
   total magnetic field, including the loop carrier when present.
   Complete magnetic-plus-electric SIBC reciprocity is checked independently.
   The old incident-field distribution rescaled to BIE total power is NOT
   a validated local heat field, even if its integral matches the loss.
   Weak EM postprocessing currently requires P1 and uniform Zs; P2 BEM and
   per-panel ESIM fail fast. This does not restrict thermal H1 order 2.
2. Select or generate the native IH configuration from the masked
   ``IH Parameters`` / geometry-update blocks.
3. Select the workpiece thermal ``wp.vol`` (a SEPARATE mesh from the
   EM ``.vol``; the EM mesh treats the workpiece as a hole with a
   SIBC face, the thermal mesh treats it as a real solid).
4. Set the material, convection, time-step, and geometry parameters.
5. Drive current and angle with ordinary Simulink source blocks and run.

The headless ``calc_heat.py`` and ``calc_heat_axisym.py`` paths remain reusable
batch and validation entry points.  They are not desktop interfaces.

## Thermal dimensionality: 3D EM does not require 3D heat

Native operator assembly accepts ``--axisymmetric-thermal-vol`` for a separate
2D (r,z) workpiece mesh; the EM workpiece remains the 3D first argument.
The MATLAB ``assembleIHOperatorsFromGeometry`` options use
``axisymmetric_thermal_vol`` and ``n_phi_samples`` (default 128).
New Geometry Update blocks expose those settings and fingerprint the thermal
mesh. Existing blocks need rebuilding to expose the new fields.
The 2D label contract requires material ``workpiece``, physical boundary
``sibc`` and optional ``axis``. Do not include r=0 in ``sibc``.
Mass, conduction, surface heat and convection use 2*pi*r. The native route
supports thermal_order=1 or 2 independently of the P1 surface heat handoff.
The executed public reproduction is
``docs/induction_heating/axisymmetric_p2_thermal.ipynb``. It verifies an
in-memory quadratic transient and displays separately recorded MEX/Simulink
results; it does not certify production CAD/VOL coupling or release acceptance.
Order 2 uses NGSolve H1 coefficients (P2 on triangles, Q2 on quads), not
nodal kelvin values. The configuration stores the constant-function vector
and a sparse mapped evaluation operator; initialization, ambient convection,
and monitor temperature statistics must use them. The raw temperature port
contains these coefficients in this mode. Periodic coefficient transport is
not supported and fails explicitly. The monitor extrema are sampled at mapped
volume quadrature points, not certified global polynomial extrema.
Transfer uses the exact circumferential average of
``radia.ih_thermal`` (see "Heat-source transfer contract" below) and rejects
an integrated power difference over 2%; it does not silently rescale the
source. The transfer difference is
part of the error budget, not proof of 2% total application accuracy.
Convection in this initial native axisymmetric route acts on ``sibc`` only;
independently selected inner/outer convection coefficients are not supported.

Use ``calc_heat_axisym.py`` when the workpiece shape, thermal properties,
thermal boundary conditions, and the circumferentially averaged heat input are
rotation invariant.  The electromagnetic solve may still be 3D: transfer its
surface-loss field by averaging around each ``(r,z)`` circle, then solve the
scalar temperature field on the 2D meridian with standard NGSolve H1 and the
``2*pi*r`` Jacobian.  This is the production choice for an axisymmetric
workpiece after the EM hotspot has been rotation/time averaged.

Use the 3D thermal solver only when angular temperature structure is physically
required, for example a non-axisymmetric workpiece or material, angularly
different cooling/contact boundaries, or a transient rotating hotspot whose
period is not short compared with the thermal response.  Do not select 3D heat
merely because the supplied EM mesh or ``qsurf.sol`` is 3D.

A study that supplies a fixed full-workpiece 3D ``.vol`` and varies only H1
order may use 3D heat; that is an experiment constraint, not a formulation
requirement.  An independent 2D meridian cross-check (H1 order 2, 128 azimuth
samples) against the 3D order-2 solve is the recommended confirmation.

Axisymmetric heat requires a separately generated 2D ``(r,z)`` workpiece mesh.
Loading a 3D ``.vol`` and calling ``Curve()`` cannot turn it into a 2D mesh:
``Curve()`` changes geometry order, not dimension.  Moreover, post-load
``Curve()`` is forbidden by the thermal geometry contract below.  Generate the
2D mesh in memory or export it at the intended geometry order before solving.

## ``.sol + .vol`` strict contract (4.58.0+)

NGSolve ``.sol`` files are coefficient vectors only -- no embedded
mesh, no FES-order header.  Both files MUST be passed:

```bash
# Direct CLI (calc_heat.py):
python -m radia.panels.calc_heat \\
    --wp-vol workpiece_thermal.vol \\
    --qsurf-sol  <stem>_qsurf.sol \\  # REQUIRED for spatial mode
    --em-vol     <stem>_fem.vol   \\  # REQUIRED (no auto-locate)
    --qsurf-order 1                \\  # cross-mesh transfer is P1 only
    --heat-flux-boundaries heated_outer \\
    --convection-boundaries 'outer|top|bottom' \\
    --material steel --dt 0.5 --t-end 5.0 \\
    --rotation-rpm 12                  # see ``rotating`` topic
```

The boundary roles are independent.  ``--heat-flux-boundaries`` is always
required; ``--convection-boundaries`` is required when ``--h-conv`` is nonzero;
and ``--radiation-boundaries`` is required when emissivity is nonzero.  A role
may use an NGSolve expression such as ``outer|top`` and roles may overlap
intentionally.  Empty selectors never expand silently to every boundary.  The
removed ``--surface-label`` fails with migration guidance.

In a 2D axisymmetric thermal mesh, the true ``r=0`` center axis is the natural
symmetry boundary, not a physical heat-transfer surface.  Its revolved area is
zero because the weak form uses ``2*pi*r*ds``.  The axisymmetric solver fails
fast when any active heat-flux, convection, or radiation selector includes an
axis boundary element, even when a broad boundary label also covers non-axis
edges.  Exclude the axis and label each coil-facing or exposed physical surface
separately.  An inner cylindrical surface at ``r>0`` is not the center axis: it
has nonzero revolved area and may receive heat when explicitly selected.

Cross-mesh q_surf transfer currently accepts ``--qsurf-order 1`` only.  The
loader rejects higher order because its surface transfer is vertex-sampled and
higher-order H1 coefficients are hierarchical; accepting order 2 would silently
distort the heat source.  ``calc_fem_kelvin.py`` therefore saves this handoff
as P1 even when the electromagnetic solve itself uses a higher order.

## Heat-source transfer contract (radia 5.1)

``radia.ih_thermal`` owns the EM-to-thermal handoff for ``calc_heat.py``,
``calc_heat_axisym.py`` and the Simulink operator assembly:

* **Sidecars are required.** EM runs write ``<stem>_qsurf.sol.json`` (and
  a P1 ``<stem>_Ht.sol`` + sidecar): mesh file digest, H1 order, DOF count,
  quantity, heated boundaries, ``P_wp_W`` and (|H_t|) ``frequency_Hz``.  The
  heat solvers verify the pair, the quantity and the EM power.  A ``.sol``
  without a sidecar is an error; after confirming the pairing, adopt it
  explicitly:

  ```bash
  python -m radia.ih_thermal sidecar --sol q.sol --mesh em.vol --order 1 \
      --quantity surface_loss_density --boundaries "wp_outer|wp_top" --p-wp 1234.5
  ```

  (``--frequency`` for ``surface_tangential_H_amplitude``; ``--geometry 3d|axisymmetric-rz`` and
  ``--definedon`` for temperature.)  Temperature fields carry the same
  sidecar; ``ih_thermal.load_field(T.sol, quantity="temperature")`` rebuilds
  them exactly (order, ``definedon`` region) and never calls ``Mesh.Curve``.
* **Pointwise transfer.** Thermal heat-flux vertices are projected onto the
  nearest EM heated triangle; a vertex farther than ``--transfer-tolerance``
  (default a quarter of the median EM edge, at most one median edge) is an
  error; there is no zero-flux fallback.
* **Circumferential average** (``--q-phi-average``; always in the
  axisymmetric solver).  Every EM triangle is cut at meridian arc-length bin
  edges and integrated exactly, so the average conserves the EM power and has
  no azimuth sampling error; each thermal vertex takes the value at its
  meridian position.  Triangles pierced by the axis are split at the
  piercing point (exact for P1), so the bins at ``r = 0`` are integrated, not
  interpolated; an empty bin is an error.  A surface that is not a body of
  revolution about the axis is rejected.  ``--q-phi-bin`` sets the bin
  width.
  ``--q-phi-average-n`` was removed; in the axisymmetric solver
  ``--n-phi-samples`` (default 128) now drives an independent pointwise ring
  check reported as ``qsurf_projection.pointwise_ring_check``.
* **Power gate.** The thermal heat input must match the EM power over
  ``--em-heat-boundaries`` (default: the boundaries recorded in the sidecar;
  thermal boundary names are never borrowed) within ``--power-tolerance``
  (default 2%).  Every transfer is gated, including each rotation resample.
  A mismatch usually means the thermal model heats only part of the EM
  heated surface.
* ``--q-uniform`` is a constant flux; any option that shapes a spatial
  source (``--qsurf-sol``, ``--em-table``, ``--q-scale`` ...) given with it is
  an error.
* ``--q-scale`` multiplies q (current sweeps, ``(I/I_ref)^2``) and is
  recorded.

Result JSON ``qsurf_projection`` records the mode, provenance, EM boundaries
and rule, transfer distances, the power balance and (phi-average) the
profile audit.

The headless runner and Simulink initialization both fail before solving when
the ``.sol`` / ``.vol`` pair is incomplete or incompatible.  They do not infer
a missing companion file from a filename convention.

## ``.sol`` data model: H1 GridFunction with boundary-only DOFs

`calc_fem_kelvin.py` saves q_surf as:

```python
from radia.ih_thermal import lumped_surface_p1
gf_q, on_boundary = lumped_surface_p1(mesh, q_surf_cf, wp_region)
# Positive lumped projection: fixed P1 cross-mesh handoff; interior DOFs stay zero.
gf_q.Save(qsurf_sol_path)
```

The thermal solver re-projects this onto the wp-mesh surface vertices by
surface point evaluation.  The axisymmetric route explicitly requests a
boundary mapped point so equivalent independently faceted 2D/3D surfaces are
handled as a surface projection and audited.  See
``radia_mcp.radia_ngsolve.ngsolve`` Section 18c for the broader pattern
(surface-restricted H1 GF save/load).

## Output: T .sol re-loadable for later evaluation (radia 4.59.0+)

calc_heat / calc_heat_axisym now ALWAYS save the final
temperature GridFunction as a `.sol` file, symmetric with the
qsurf.sol contract on the EM side.  This means a thermal run is
fully re-loadable: you can come back later and sample T at
arbitrary points, feed T(x) back into a temperature-dependent
sigma for a second EM solve, or feed it into an in-house
post-processor without re-running the heat equation.

JSON output keys (calc_heat / calc_heat_axisym):

| key | content |
|---|---|
| ``T_mean_C`` | physical volume mean; the axisymmetric route uses ``int(T*2*pi*r)/int(2*pi*r)`` |
| ``T_sol_file`` | absolute path to the final-T NGSolve `.sol` |
| ``heat_vol_file`` | absolute path to the companion `.vol`.  Empty when no separate companion was written (then re-use the `--wp-vol` input as the companion). |
| ``msh_file`` | GMSH `.msh v4.2` (T_C + q_surf fields) when `--msh-output` was set |
| ``csv_file`` | probe history CSV when both `--probe-point` + `--csv-output` were set |
| ``boundary_audit`` | Concrete matched heat-flux, convection, and radiation boundaries; per-boundary area; and per-boundary heat input |
| ``qsurf_projection`` | transfer mode, provenance, EM boundaries, distances, power balance (or ``mode=uniform``) |
| ``thermal_exposure`` | volume above each ``--exposure-thresholds`` / ``--temperature-limit`` value, its bounding box and meridian range, the peak and its location, and (3D) the hottest-ring spread |
| ``temperature_limit`` | constraint record for optimisers: exceeded, excess over the limit, volume above it |
| ``T_max_history_C`` | nodal peak temperature per step |
| ``thermal_material`` / ``nonlinear_transient`` | material audit; Newton counts, halvings and energy balance when temperature-dependent |

Naming convention:

* With `--msh-output FILE.msh`: writes `FILE_T.sol` +
  `FILE_heat.vol` (fresh companion mesh) alongside `FILE.msh`.
* Without `--msh-output`: writes `<wp-stem>_heat_T.sol` next to
  `wp.vol`.  No separate companion mesh; reload directly against
  the same `wp.vol` the solve consumed.

Reload pattern (radia 5.1): read the field through its sidecar, which
verifies the mesh digest and order and never curves the mesh:

```python
from radia import ih_thermal
mesh, gfT, audit = ih_thermal.load_field("workpiece_thermal_heat_T.sol",
                                         quantity="temperature")
T_at = float(gfT(mesh(x, y, z)))
```

Do not call ``Mesh.Curve`` on the reader side "to match the writer": loading
already applies the stored order, and a reported thermal workpiece mesh had a
6.9e6 m^2 surface after ``Curve(2)``.  That mesh had been generated in
millimetres and rescaled with ``ngmesh.Scale(1e-3)``, so its saved CAD stayed in
millimetres and ``Curve`` projected the nodes onto it; every later point
evaluation was wrong.  Scale the CAD shape before meshing instead.

## Post-processing: exposure and case depth (radia 5.1)

``radia.ih_thermal_post`` (MATLAB ``radia.python.ihThermalPost``):

* ``thermal_exposure(mesh, gfT, thresholds, axisymmetric=..., region=...)``
  -- volume above each temperature from integration points (a molten skin
  thinner than an element is counted), bounding boxes, peak location, and for
  3D the temperature spread on the ring through the peak; a large spread
  under a supposedly axisymmetric source means the source was not
  axisymmetric.
* ``case_depth(mesh, gfT, 850.0, span=..., boundary_names=... | origins=,
  normals=, n_phi=...)`` -- ``span`` (probe length) is required.  Station
  origins are projected onto the material surface (a station farther than a
  quarter edge is an error; distances are reported), then marched inward;
  every ray is classified ``ok``, ``not_reached`` (surface below threshold),
  ``through`` (leaves the material still hot), ``beyond_span`` (probe too
  short; lower bound only) or ``edge`` (automatic station on an edge or
  corner, no depth).  Meridian stations on a 3D body report min/mean/max over
  the azimuth; a 2D (r, z) mesh takes (r, z) stations.
* CLI: ``python -m radia.ih_thermal_post {exposure,depth} --temperature T.sol
  ...`` reads the field, its ``definedon`` region and its geometry from the
  sidecar and writes JSON/CSV.  Boundaries are ``|``-separated.  A stations
  CSV has exactly the header ``x_m,y_m,z_m,nx,ny,nz`` or ``r_m,z_m,n_r,n_z``
  (``_mm`` variants; optional ``s_m``/``s_mm`` arc-length column).

## Temperature-dependent heating and the Curie band (radia 5.1)

* ``--material-table T_C,k_W_mK,cp_J_kgK`` (+ ``--latent-heat``,
  ``--latent-range``) switches the heat solvers to an enthalpy
  backward-Euler Newton integrator.  The energy balance is checked at the
  end of the run (relative tolerance 1e-4; it closes to rounding when Newton
  converges).  A step that does not converge is an error; ``--max-halvings N``
  allows N explicit step halvings (default 0), and the halvings are
  reported.  Leaving the table is an error unless
  ``--allow-table-extrapolation``.
* **Do not freeze |H_t| for ferromagnetic workpieces.**  For a
  current-driven coil the surface |H_t| falls as sigma(T) falls and rises
  when the surface passes the Curie band.  In the validation case
  (``validation_test/induction_heating/coupled_curie_cylinder.py``) the true
  power stays within +5% until the surface reaches the Curie band, peaks
  at +34% as the band passes and is back at the initial power 1 s later;
  a frozen-|H_t| table instead rises to +53% at once and then stalls the
  surface at the Curie temperature (767 C against 1066 C), a false
  self-limit.  ``--em-table`` therefore requires
  ``--allow-frozen-ht`` (non-magnetic parts, or between EM re-solves), and
  ``calc_heat_with_em_table.py`` holds for non-magnetic parts only.
* For bodies of revolution with coaxial coils use the coupled solver
  ``calc_ih_axisym_coupled.py`` (``radia.ih_axisym_coupled``): volumetric
  A_phi eddy currents with sigma(T), mu_r(T) (``T_C,sigma_S_m,mu_r``)
  re-solved from the temperature every ``--em-every`` steps, Joule heat in
  the volume, exact coil/Joule power identity checked every solve.  A
  skin-depth gate refuses a mesh whose workpiece-interface edges exceed
  ``--skin-resolution`` (default 1) skin depths at any EM solve: at order 2
  the Bessel error is 0.3 % at h/delta = 1 and 6 % at 2.  It resolves the skin
  depth, which above the
  Curie band reaches millimetres (6 mm at 7 kHz for hot steel) -- the regime
  where a surface-impedance model is itself no longer valid.  A rotating part
  under a coil with short leads is usually well described by the ring coil
  of the same cross-section.
* 3D temperature-dependent EM through the BEM-SIBC route needs a spatially
  varying Z_s with a consistent weighted operator; the per-panel path is not
  yet accepted by the verified heat post-processing.

Same three contracts as the qsurf side (see "Strict .sol + .vol
contract" above): the `.sol` is a raw coefficient vector, the
`.vol` carries the mesh, the FES order must match.

## Thermal `.vol` geometry contract: preserve after load

The thermal field order and the serialized mesh geometry order are independent:

```python
wp_mesh = Mesh("workpiece_thermal.vol")
fes_T = H1(wp_mesh, order=2)  # field P2/Q2; do NOT call wp_mesh.Curve(2)
```

`calc_heat.py`, `calc_heat_axisym.py`, and `calc_heat_with_em_table.py` preserve
the loaded `.vol` exactly.  They never infer a geometry order from
`--fes-order`.  If a curved thermal mesh is required, the mesher must create it
at the desired order before saving the `.vol`; re-export rather than attempting
to curve it after load.  This is a fail-safe geometry policy, not a flat-mesh
fallback.

The result JSON records this decision in `mesh_geometry`: policy
`preserve-input-vol-geometry`, `post_load_curve_applied=false`, the field and
input curve orders, and the consumed domain and boundary measures.  On the
reported curved CAD `.vol`, the removed post-load `Curve(2)` call returned
without an exception but inflated the boundary measure from about 0.0233 square
metres to 6.91 million square metres (about 297 million times).  Therefore
`GetCurveOrder()` or lack of an exception is not evidence that post-load curving
was valid.  The cause was the mesh-in-millimetres-then-`ngmesh.Scale(1e-3)`
workflow, which keeps the saved CAD in millimetres (lint rule
`netgen-scale-after-generate`; see `docs/ngsolve_integration/curve_order.md`).

## Axisymmetric discretization contract: Henrotte for EM, NGSolve H1 for heat

Do not reuse the electromagnetic space for temperature:

| Field | Required space | Axisymmetric operator |
|---|---|---|
| electromagnetic `A_phi` | `radia.axifem.H1Henrotte(...)` | Henrotte magnetic BFIs |
| temperature `T` | standard `ngsolve.H1(mesh, order=2)` | symbolic forms weighted by `2*pi*x` |

For standard NGSolve H1, `order=2` selects **Q2 on quadrilateral meshes** and
**P2 on triangular meshes**.  It is an approximation order, not a request for
one fixed element topology.  The Henrotte space carries the `A_phi` basis
`N_i = r_i psi_i / r`; it has no temperature integrators, because a scalar
temperature is not of that form (the former Henrotte heat BFIs were removed
with the V-DOF convention unification).

## Axisymmetric thermal FE order: default 2 (near-axis cusp at order 1)

`calc_heat_axisym.py` defaults to `--fes-order 2` (2026-09-03 near-axis
study).  Standard P1/Q1 cannot reproduce a nonconstant even quadratic
radial profile while also representing dT/dr = 0 at the r = 0 axis, so
its near-axis profile shows a piecewise-linear cusp.  Order 2 contains
r^2 and strongly suppresses that artifact: the committed 9x9
uniform-flux cylinder validation reduces the spurious radial axis slope
from about 139 C/m to below 0.3 C/m and compares both axis and surface
temperatures with an independent Bessel-series solution.

Field-report symptom map: "the axisym thermal result looks wrong near
the axis" == this order-1 cusp.  Fix: keep the default order 2 (or pass
`--fes-order 2`).  Do NOT switch the scalar solve to the Henrotte
basis: the reference standard-H1 policy stands, order 2 already
strongly suppresses the axis artifact.

Two order>=2 implementation notes (both fixed 2026-09-03; they apply
to `calc_heat.py` as well):

* the uniform initial state must be `gfT.Set(CF(T0))`, NOT
  `gfT.vec[:] = T0` -- hierarchical H1 edge coefficients are not nodal
  values, so a constant coefficient vector is not a constant field;
* `T_max_C` / `T_min_C` come from physical GridFunction evaluations
  over vertices plus deterministic volume and boundary samples.  Raw
  coefficient extrema report garbage at order 2, while vertex-only
  extrema can miss higher-order edge, face, or cell modes.

## Heat Equation with Joule Heat Source

Solves the unsteady heat equation in the workpiece:

    rho * c * dT/dt = div(kappa * grad(T)) + Q

with convection and/or radiation boundary conditions on the work surface:

    -kappa * dT/dn = h * (T - T_ext)                    (convection)
                   + eps * sigma_SB * (T^4 - T_ext^4)   (radiation, T in KELVIN)

Radiation matters at IH temperatures (steel hardening ~900-1100 deg C, where
eps*sigma*T^4 rivals or exceeds convection). The Heat panel exposes an
``Emissivity [0..1]`` field (0 = off, the default -> convection-only, unchanged);
``calc_heat.py`` / ``calc_heat_axisym.py`` take ``--emissivity`` and add the
radiation flux EXPLICITLY (previous-step T, internally converted to Kelvin) to
the per-step RHS, so the cached ``mstar = M + theta*dt*K`` factorisation is kept:

    if emissivity > 0:                       # T in Kelvin internally
        f_T += -emissivity * 5.670374419e-8 \\
               * ((T+273.15)**4 - (T_ext+273.15)**4) * v * ds("work_surface")

Validated by the eps-sigma steady-state balance
T_ss = (q/(eps*sigma) + T_ext_K^4)^(1/4)  (validation_test/panels/test_heat_radiation.py,
0.02 %).  Radiation ambient = ``--t-ext`` (shared with convection).

## Complete Thermal Implementation

```python
from ngsolve import *
import numpy as np

# ============================================================
# 1. Material properties (workpiece only)
# ============================================================
rho_steel = 7800      # density [kg/m^3]
c_steel = 467         # specific heat [J/(kg*K)]
kappa_steel = 46.6    # thermal conductivity [W/(m*K)]
h_conv = 10           # convection coefficient [W/(m^2*K)]
T_ext = 20            # ambient temperature [deg C]

# Material CoefficientFunctions (defined on work regions)
K = CoefficientFunction(
    [kappa_steel if mat in ("work_main", "work_skin") else 0
     for mat in mesh.GetMaterials()]
)
rho_c = CoefficientFunction(
    [rho_steel * c_steel if mat in ("work_main", "work_skin") else 0
     for mat in mesh.GetMaterials()]
)

# ============================================================
# 2. FE space (H1, work regions only)
# ============================================================
order = 1
fes_T = H1(mesh, order=order,
           definedon="work_main|work_skin")
u, v = fes_T.TnT()

# ============================================================
# 3. Bilinear forms
# ============================================================
# Stiffness: conduction + convection
a_T = BilinearForm(fes_T)
a_T += K * grad(u) * grad(v) * dx("work_main|work_skin")
a_T += h_conv * u * v * ds("work_surface")

# Mass: rho * c * u * v
m_T = BilinearForm(fes_T)
m_T += rho_c * u * v * dx("work_main|work_skin")

# Source: Joule heat + convection ambient
f_T = LinearForm(fes_T)
f_T += Q * v * dx("work_main|work_skin")       # Q from EM analysis
f_T += h_conv * T_ext * v * ds("work_surface")  # Convection ambient

with TaskManager():
    a_T.Assemble()
    m_T.Assemble()
    f_T.Assemble()

# ============================================================
# 4. Time stepping (theta-scheme, implicit)
# ============================================================
dt = 0.5       # time step [s]
t_end = 5.0    # total simulation time [s]

# Combined matrix: M + dt * K
mstar = m_T.mat.CreateMatrix()
mstar.AsVector().data = m_T.mat.AsVector() + dt * a_T.mat.AsVector()
mstar_inv = mstar.Inverse(fes_T.FreeDofs(), inverse="sparsecholesky")

# Initial condition
gfT = GridFunction(fes_T)
gfT.Set(T_ext)  # Initial temperature = ambient

# Time stepping loop
times = []
temperatures = []
t = 0.0

while t < t_end:
    # RHS: M * T_n + dt * f
    rhs = gfT.vec.CreateVector()
    rhs.data = m_T.mat * gfT.vec + dt * f_T.vec

    # Solve: (M + dt*K) * T_{n+1} = M * T_n + dt * f
    gfT.vec.data = mstar_inv * rhs

    t += dt

    # Monitor temperature at a point
    T_monitor = gfT(mesh(0.0305, 0, 0))  # Work surface point
    times.append(t)
    temperatures.append(T_monitor)
    print(f"t={t:.1f}s, T={T_monitor:.1f} deg C")

print(f"Final temperature: {temperatures[-1]:.1f} deg C")
```

## Important Notes

### Why theta=1 (Backward Euler)?

- **Unconditionally stable**: No CFL condition on dt
- **First-order accurate**: Sufficient for engineering (temperature rise is smooth)
- Crank-Nicolson (theta=0.5) is second-order but may oscillate with sharp sources

### Time Step Selection

```python
# Characteristic thermal diffusion time
alpha = kappa_steel / (rho_steel * c_steel)  # thermal diffusivity [m^2/s]
h_elem = 0.001  # typical element size [m]
dt_diffusion = h_elem**2 / (2 * alpha)  # ~10^-5 s (CFL for explicit)

# For implicit scheme: dt can be much larger
# Practical guideline: dt ~ 0.1-1.0 s for engineering accuracy
dt = 0.5  # s
```

### Temperature Monitoring

```python
# Point evaluation
T_surface = gfT(mesh(0.0305, 0, 0))

# Volume average
T_avg = Integrate(gfT * dx("work_main|work_skin"), mesh) / \
        Integrate(1 * dx("work_main|work_skin"), mesh)
```

## Coupled EM-Thermal with Pickle Serialization

For large EM problems, save the electromagnetic solution and reuse it
for multiple thermal analyses (parameter sweeps, time studies):

```python
import pickle

# ============================================================
# Save EM results (after electromagnetic solve)
# ============================================================
em_data = {
    'gfA_vec': gfA_sol.vec.FV().NumPy().copy(),
    'gfPhi_vec': gfPhi_sol.vec.FV().NumPy().copy(),
    'freq': freq,
    'sigma_d': sigma_d,
    'mur_d': mur_d,
    'I_target': I_target,
    'I_out': I_out,
}
with open('em_solution.pkl', 'wb') as f:
    pickle.dump(em_data, f)

# ============================================================
# Load EM results (for thermal analysis)
# ============================================================
with open('em_solution.pkl', 'rb') as f:
    em_data = pickle.load(f)

# Reconstruct GridFunction
gfA_sol.vec.FV().NumPy()[:] = em_data['gfA_vec']
gfPhi_sol.vec.FV().NumPy()[:] = em_data['gfPhi_vec']
scale = em_data['I_target'] / abs(em_data['I_out'])

# Recompute Q from loaded fields
E = -s * (gfA_sol + grad(gfPhi_sol)) * scale
Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real
```
"""

INDUCTION_HEATING_ROTATING = """
# Rotating Workpiece Simulation (radia 4.58.0+)

## Production path: angle is an explicit Simulink signal

The workpiece angle is generated on the Simulink side. Temperature stays in
workpiece material coordinates. Eddy maps temperature to the source frame and
heat back; Thermal must not rotate the accepted temperature a second time.
The generated native preview currently fixes rotation=none; a supplied angle
does not enable an unimplemented geometry-dependent EM solve. The headless path
can generate the equivalent angle history from ``--rotation-rpm``:

```bash
python -m radia.panels.calc_heat \\
    --wp-vol     workpiece_thermal.vol \\
    --heat-flux-boundaries heated_outer \\
    --convection-boundaries 'outer|top|bottom' \\
    --qsurf-sol  ih_em_qsurf.sol \\
    --em-vol     ih_em_fem.vol \\
    --material   steel \\
    --dt 0.5 --t-end 5.0 \\
    --rotation-rpm 12        # <-- workpiece spins around z at 12 rpm
```

## How calc_heat.py implements rotation (radia 5.1)

The part turns under a coil fixed in the world; a point ``p`` of the part
sits at ``R(theta) p``.  The thermal mesh, FE space and matrices stay fixed
in the part; only the heat-flux source changes, as
``q(p, theta) = q_EM(theta)(R(theta) p)``.

``ih_thermal.RotatingSurfaceSource`` tabulates that source at knots over one
revolution, each knot a gated transfer (off-surface vertices and a power
mismatch are errors), and takes it piecewise linear in the angle:

* **Body of revolution** (``--qsurf-sol`` + ``--em-vol``): one EM solution,
  turned with the part; knots every half EM edge along the outermost path.
  A part that is not of revolution leaves the EM surface when turned and is
  refused with a pointer to ``--rotor-states``.
* **Any other part** (cross holes, flats, keyways): ``--rotor-states
  rotor.json`` lists one EM solution per rotor angle over the symmetry
  period (schema ``radia.ih-rotor-states/1``: ``axis``, ``period_rad`` =
  2 pi / n, ``frame`` ``world`` (part turned in each EM mesh) or ``body``
  (coil turned by -angle), ``states`` of ``angle_rad``, ``qsurf_sol``,
  ``em_vol``).  An n-fold part repeats its EM solution every period, but a
  point returns only after a revolution, so the table is expanded to 2 pi by
  reading state k on the part turned by whole periods -- a declared symmetry
  the part does not have fails the transfer.  ``--angle-step-tolerance``
  (default 25 %) refuses states whose neighbours differ too much.

Each time step applies the **exact average over the angles it sweeps**
(``omega (t - dt) .. omega t``).  Sampling one angle per step is
stroboscopic: at 120 rpm and dt = 0.5 s every sample sees the part where it
started.  ``qsurf_projection.rotation`` records the knots, the largest
neighbour change, the worst knot power error and the angle swept per step.

## Why NOT mesh deformation

`mesh.SetDeformation(...)` deforms the FE space at each step, which
invalidates the cached mass / stiffness / factorisation.  For a 5-second
simulation with dt=0.5 (10 steps), that's 10x the full assemble + LU
factor cost (~30s each on a 200k-DOF cylinder mesh = ~5 minutes wasted).
Re-projection of q_surf is ~10ms per step on the same mesh.

## .sol contract for rotation users

The rotation path consumes the spatial qsurf path (`--qsurf-sol +
--em-vol`).  Uniform `--q-uniform` is rotation-invariant; configuration
validation rejects or diagnoses combinations that cannot express the requested
spatial rotation.

Three CONTRACTS that must hold between the EM solve (calc_fem_kelvin.py)
and the thermal solve (calc_heat.py):

1. **`.sol` is mesh-free**.  NGSolve's `GridFunction.Save()` writes
   the raw coefficient vector (no mesh, no header).  The thermal
   solver MUST be told the EM mesh path via `--em-vol`.  Pass the
   `*_fem.vol` file that `calc_fem_kelvin.py` saves alongside
   `*_qsurf.sol`.

2. **`qsurf_order` is currently fixed to 1 for cross-mesh transfer**.
   The thermal solver rebuilds `H1(em_mesh, order=qsurf_order)` and loads
   the .sol into it, then transfers surface vertex values.  It rejects
   any order other than 1 because higher-order H1 coefficients are
   hierarchical and vertex-only reconstruction is not a valid projection.
   ``calc_fem_kelvin.py`` always produces this handoff at order 1,
   independently of the electromagnetic solve order.

3. **q_surf is a fixed-P1 volume H1 GridFunction with non-zero values only
   on the workpiece boundary**. `lumped_surface_p1` performs the positive
   lumped surface projection independently of the EM polynomial order.
   Transfer uses boundary-mapped evaluation and explicit mesh/provenance checks.
   Unmapped points fail loudly; they are not silently replaced by zero.

4. **Thermal boundary roles are explicit and independently audited**.
   Heat flux, convection, and radiation use separate selectors.  The
   result JSON records each selector, every concrete matched boundary,
   its area, and the boundary-wise heat input.  Inspect ``boundary_audit``
   before interpreting a temperature hotspot.

## Test coverage

* `validation_test/panels/test_heat_rotation.py` -- validation tests for the rotation
  projection math on a synthetic unit-cube setup with q_em(x,y,z)=x.
  At theta=0 / pi / pi/2 the body point near (+1,0,0) reads
  q ~ 1 / -1 / 0 respectively.

* `validation_test/panels/test_heat_chain_golden.py` -- end-to-end chain test
  (EM solve -> qsurf.sol -> thermal solve).  Slow-marked. Export accepted
  spatial fields as indexed GMSH `.msh v4.2` artifacts.

5. **Full rotation**: At 12 rpm, one rotation = 5 seconds. Choose t_end
   to capture the desired number of rotations.
"""

INDUCTION_HEATING_POSTPROCESS = """
# Post-Processing for Induction Heating

## GMSH `.msh v4.2` Output

```python
from ngsolve import *
from radia.gmsh_post_export import GmshPostExport

# ============================================================
# Basic checked field output
# ============================================================
post = GmshPostExport(mesh)
post.add_vector_field("B_re", B.real)
post.add_vector_field("B_im", B.imag)
post.add_vector_field("J_re", J.real)
post.add_vector_field("J_im", J.imag)
post.add_scalar_field("Q", Q)
post.write("eddy_current_result.msh")

# ============================================================
# Output with material ID
# ============================================================
mat_names = mesh.GetMaterials()
mat_id_dict = {mat: i+1 for i, mat in enumerate(mat_names)}
mat_id = mesh.MaterialCF(mat_id_dict)

post = GmshPostExport(mesh)
post.add_scalar_field("MaterialID", mat_id)
post.add_vector_field("B_real", B.real)
post.add_scalar_field("Q", Q)
post.write("result_with_id.msh")

# ============================================================
# Thermal time-step files
# ============================================================
for step in range(n_steps):
    # ... solve time step ...
    post = GmshPostExport(mesh)
    post.add_scalar_field("Temperature", gfT)
    post.write(f"thermal_series_{step:04d}.msh", time=t, timestep=step)

# Index generated files from result.json for reproducible playback.
```

## Point and Line Evaluation

```python
import numpy as np

# ============================================================
# Single point evaluation
# ============================================================
# B_z component at point (0, 0, 0.03)
Bz_value = B[2](mesh(0, 0, 0.03))
print(f"B_z = {Bz_value:.6f} T")

# Temperature at work surface
T_surface = gfT(mesh(0.0305, 0, 0))
print(f"T = {T_surface:.1f} deg C")

# ============================================================
# Line evaluation (e.g., B_z along z-axis)
# ============================================================
z_points = np.linspace(-0.02, 0.08, 101)
Bz_line = np.array([B[2](mesh(0, 0, z)) for z in z_points], dtype=complex)

# Magnitude
Bz_mag = np.abs(Bz_line)

# Real and imaginary parts
Bz_real = Bz_line.real
Bz_imag = Bz_line.imag

# ============================================================
# Radial evaluation (e.g., B along r-direction)
# ============================================================
r_points = np.linspace(0.001, 0.05, 50)
Br_radial = np.array([B[0](mesh(r, 0, 0)) for r in r_points], dtype=complex)
```

## Data Export

### CSV Export

```python
import csv

# Temperature vs time
with open('temperature_history.csv', 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.writer(f)
    writer.writerow(['Time [s]', 'Temperature [deg C]'])
    for t_val, T_val in zip(times, temperatures):
        writer.writerow([f'{t_val:.2f}', f'{T_val:.2f}'])

# Field along a line
with open('Bz_along_z.csv', 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.writer(f)
    writer.writerow(['Z [m]', 'Bz_real [T]', 'Bz_imag [T]', 'Bz_abs [T]'])
    for z, bz in zip(z_points, Bz_line):
        writer.writerow([f'{z:.6f}', f'{bz.real:.6e}', f'{bz.imag:.6e}',
                        f'{abs(bz):.6e}'])
```

### MATLAB Format (.mat)

```python
from scipy.io import savemat

# Save field data for MATLAB post-processing
savemat('field_results.mat', {
    'Z': z_points,
    'Bz_real': Bz_line.real,
    'Bz_imag': Bz_line.imag,
    'Bz_abs': np.abs(Bz_line),
    'freq': freq,
    'I_target': I_target,
})
```

### Animated GIF from GMSH Field Frames

```python
# animated_gif.py: Create animated GIF from PNG frames
from PIL import Image
import glob

# Assuming the GMSH render step has exported PNG frames
png_files = sorted(glob.glob("thermal_step_*.png"))
frames = [Image.open(f) for f in png_files]

frames[0].save(
    'temperature_animation.gif',
    save_all=True,
    append_images=frames[1:],
    duration=200,    # ms per frame
    loop=0           # infinite loop
)
```

## Power and Energy Computation

```python
# Total Joule power dissipated in work [W]
P_work = Integrate(Q * dx("work_main|work_skin"), mesh)
print(f"Power in work: {P_work:.1f} W")

# Total Joule power in coil [W] (ohmic loss)
P_coil = Integrate(Q * dx("coil_main|coil_skin"), mesh)
print(f"Power in coil: {P_coil:.1f} W")

# Efficiency
eta = P_work / (P_work + P_coil) * 100
print(f"Heating efficiency: {eta:.1f}%")

# Total energy delivered over time [J]
E_total = P_work * t_end
print(f"Total energy: {E_total:.0f} J")

# Expected temperature rise (rough estimate)
V_work = Integrate(1 * dx("work_main|work_skin"), mesh)
m_work = rho_steel * V_work
dT_expected = E_total / (m_work * c_steel)
print(f"Expected dT (uniform): {dT_expected:.1f} deg C")
```
"""

INDUCTION_HEATING_PITFALLS = """
# Common Pitfalls in Induction Heating Simulation

## 1. CRITICAL: nograds with A-Phi Formulation

**Correct**: Use `nograds=True` in A-Phi formulation.

In A-Phi, the scalar potential Phi provides the gauge fixing in the conductor.
The `nograds=True` flag removes gradient basis functions from HCurl, which is
still needed because the curl-curl system in air regions has a gradient null space.

```python
# CORRECT: nograds=True even for eddy current A-Phi
fesA = HCurl(mesh, order=order, nograds=True,
             dirichlet="outer", complex=True)

# WRONG: nograds=True in T-Omega (T already gauged by curl)
# This would also be wrong for a pure time-harmonic A formulation
# without the Phi coupling.
```

## 2. CRITICAL: Missing complex=True for Frequency Domain

```python
# WRONG: real-valued spaces for eddy current analysis
fesA = HCurl(mesh, order=order, dirichlet="outer")  # Missing complex=True!

# CORRECT: complex-valued for frequency domain
fesA = HCurl(mesh, order=order, dirichlet="outer", complex=True)
fesPhi = H1(mesh, order=order, definedon="coil", complex=True)
```

## 3. HIGH: Phi definedon Must Match Conductive Regions

```python
# WRONG: Phi defined everywhere
fesPhi = H1(mesh, order=order, complex=True)  # Phi in air = singular!

# CORRECT: Phi only in conductive regions
fesPhi = H1(mesh, order=order,
            definedon="coil_main|coil_skin",
            complex=True)
```

## 4. HIGH: Stabilization Term for A-Phi

The A-Phi system requires a small stabilization (regularization) term
to improve conditioning:

```python
# Without stabilization: poor conditioning, slow convergence
a += nu * curl(A) * curl(N) * dx
a += s * sigma * (A+grad(phi)) * (N+grad(psi)) * dx("coil")

# With stabilization: well-conditioned
a += nu * curl(A) * curl(N) * dx
a += 1e-6 * nu * A * N * dx          # <-- Stabilization
a += s * sigma * (A+grad(phi)) * (N+grad(psi)) * dx("coil")
```

The factor `1e-6` is typically sufficient; larger values may affect accuracy.

## 5. HIGH: Preconditioner Registration Order

```python
# WRONG: Preconditioner after assembly
a.Assemble()
c = Preconditioner(a, "bddc")  # Too late!

# CORRECT: Register BEFORE assembly
c = Preconditioner(a, "bddc")
a.Assemble()
c.Update()
```

## 6. MODERATE: Q Computation Must Use Conj(E)

```python
# WRONG: Q = 0.5 * sigma * E * E (complex dot product != power)
Q_wrong = 0.5 * sigma * InnerProduct(E, E)  # Complex, not real power!

# CORRECT: Q = 0.5 * sigma * |E|^2 = 0.5 * sigma * E . conj(E)
Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real
```

## 7. MODERATE: Thermal definedon Must Match Work Regions

```python
# WRONG: Heat equation over entire mesh
fes_T = H1(mesh, order=order)  # Includes air = meaningless

# CORRECT: Only work regions
fes_T = H1(mesh, order=order, definedon="work_main|work_skin")
```

## 8. MODERATE: Mass Matrix Needed for Transient Thermal

```python
# WRONG: Steady-state solve for transient problem
a.Assemble()
gfT.vec.data = a.mat.Inverse(freedofs) * f.vec  # No time derivative!

# CORRECT: theta-scheme with mass matrix
mstar = m.mat + dt * a.mat  # M + dt*K
mstar_inv = mstar.Inverse(freedofs)
gfT.vec.data = mstar_inv * (m.mat * gfT.vec + dt * f.vec)
```

## 9. HIGH: Do Not Reuse the Henrotte A_phi Space for Temperature

```python
# WRONG: the Henrotte space is the A_phi basis N_i = r_i psi_i / r
fes_T = radia.axifem.H1Henrotte(mesh, order=2)

# CORRECT: standard scalar H1 with the cylindrical 2*pi*r measure
fes_T = H1(mesh, order=2)  # Q2 on quads, P2 on triangles
u, v = fes_T.TnT()
a += k * grad(u) * grad(v) * 2*pi*x * dx
m += rho_c * u * v * 2*pi*x * dx
```

The Henrotte space has no temperature integrators; select the standard H1
thermal formulation explicitly.

## 10. MODERATE: Missing TaskManager for Parallel Assembly

```python
# SLOW: Serial assembly
a.Assemble()  # Uses single thread

# FAST: Parallel assembly
with TaskManager():
    a.Assemble()  # Uses all CPU cores
```

## 11. LOW: Point Evaluation Syntax

```python
# WRONG: Direct coordinate (only works for scalar fields)
T = gfT(0.03, 0, 0)  # May fail for complex meshes

# CORRECT: mesh() mapping first
T = gfT(mesh(0.03, 0, 0))

# For vector field components
Bz = B[2](mesh(0, 0, 0.03))  # B_z component
```

## 12. LOW: Bypassing the Radia Post-Processing Contract

```python
# WRONG for a Radia-owned workflow: direct VTK-family output

# CORRECT: checked GMSH .msh v4.2 output
from radia.gmsh_post_export import GmshPostExport

post = GmshPostExport(mesh)
post.add_scalar_field("Temperature", gfT)
post.write("temperature.msh")
```
"""


INDUCTION_HEATING_ESIM_KELVIN = """
# ESIM + Kelvin Integration for Induction Heating

## Motivation: Three Approaches to Workpiece Heating

| Method | Workpiece mesh | Open BC | Skin mesh | Nonlinear BH |
|--------|---------------|---------|-----------|--------------|
| **FEM-full** (current) | Volume (skin layer) | Dirichlet/PML | Required | Via mu(H) |
| **FEM-ESIM + Kelvin** (target) | Surface only (SIBC) | Kelvin (exact) | NOT needed | Via ESIM Z_s(H) |
| **BEM-ESIM** | None (panels) | Exact (BEM) | NOT needed | Via ESIM Z_s(H) |

FEM-ESIM + Kelvin combines the best: no skin mesh, exact open BC, nonlinear BH.

## Reference: Go-Tech Toymodel (FEM-full, production code)

Source: internal production Go-Tech toymodel notes (machine-local path omitted).

### A-Phi Formulation (3D, frequency domain)

```python
# FE spaces: HCurl (A) + H1 (phi on coil only)
fesA = HCurl(mesh, order=order, nograds=True, dirichlet="...", complex=True)
fesPhi = H1(mesh, order=order, definedon="coil_main|coil_skin", dirichlet="...", complex=True)
fesAPhi = fesA * fesPhi
(A, phi), (N, psi) = fesAPhi.TnT()

# Bilinear form
s = 1j * 2 * pi * freq
a += nu * curl(A) * curl(N) * dx                           # curl-curl
a += 1e-6 * nu * A * N * dx                                # gauge regularization
a += s * sigma * (A + grad(phi)) * (N + grad(psi)) * dx("coil_...")  # eddy current

# Current scaling: solve for arbitrary I, then scale to target
I_out = Integrate(J * grad(psi) * dx("coil_..."), mesh)  # computed current
B = curl(gfA) * I_target / I_out                          # scale to I_target
E = -s * (gfA + grad(gfPhi)) * I_target / I_out
Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real                 # Joule heat [W/m^3]
```

### Mesh Structure (Gmsh)

| Region | Material | Role |
|--------|----------|------|
| `coil_main` | Copper (sigma=57e6) | Coil bulk |
| `coil_skin` | Copper | Coil boundary layer |
| `work_main` | Steel (sigma=10e6, mu_r=1000) | Workpiece bulk |
| `work_skin` | Steel | Workpiece boundary layer |
| `Sair` | Air | Near-field air |
| `air` | Air | Far-field air |

Boundary layer mesh (`work_skin`) resolves skin depth delta ~ 0.18 mm at 8 kHz.

### Thermal Coupling

```python
Q = 0.5 * sigma * InnerProduct(E, Conj(E)).real  # from EM solve
# theta-scheme: rho*c * dT/dt = div(k*grad(T)) + Q
# dt = 0.5s, t_end = 5s, h = 10 W/m^2K convection
```

## FEM-ESIM + Kelvin Design (Replacing FEM-full)

### Changes from FEM-full to FEM-ESIM

1. **Keep** workpiece as separate material with air/vacuum properties (NOT a hole!)
   ```python
   # OCC: CORRECT pattern (use boolean subtract)
   # workpiece meshed as air, internal interface for SIBC
   air = air_sphere - torus - wp_cyl
   air.name = "air"
   wp_cyl.name = "workpiece"   # same nu0 as air, separate material
   shape = Glue([air, wp_cyl, torus])  # THREE volumes (interface approach)
   # NOTE: for SIBC, hole approach (wp NOT in Glue) is also valid
   ```

   **Cubit equivalent** (same pattern, different syntax):
   ```python
   # Coil: sweep (NOT webcut, see cubit_scripting_knowledge sweep policy)
   cubit.cmd('sweep surface 1 axis 0 0 0 0 0 1 angle 355')
   coil_vid = cubit.get_last_id("volume")
   # Workpiece cylinder
   cubit.cmd('create cylinder height ... radius ...')
   wp_vid = cubit.get_last_id("volume")
   # Air sphere
   cubit.cmd('create sphere radius 0.120')
   air_vid_orig = cubit.get_last_id("volume")
   # CRITICAL: subtract BOTH coil and workpiece from air sphere.
   # imprint+merge alone does NOT share surfaces between fully nested
   # volumes — Cubit reports "Consolidated 0 pair of surfaces" and the
   # FEM treats coil/wp walls as PEC (wrong).
   cubit.cmd(f'subtract volume {coil_vid} {wp_vid} from volume {air_vid_orig} keep_tool')
   all_vols = list(cubit.parse_cubit_list("volume", "all"))
   air_vid = [v for v in all_vols if v not in (coil_vid, wp_vid)][0]
   cubit.cmd(f'imprint volume {coil_vid} {wp_vid} {air_vid}')
   cubit.cmd(f'merge volume {coil_vid} {wp_vid} {air_vid}')
   # Verify: "Consolidated N pairs of surfaces" with N > 0
   # blocks
   cubit.cmd(f'block 1 add volume {coil_vid}'); cubit.cmd('block 1 name "coil"')
   cubit.cmd(f'block 2 add volume {air_vid}');  cubit.cmd('block 2 name "air"')
   cubit.cmd(f'block 3 add volume {wp_vid}');   cubit.cmd('block 3 name "workpiece"')
   ```

   **Why this matters**: Without `subtract` first, Cubit's
   `imprint+merge` does nothing for fully nested volumes. The exported
   .vol then has FaceDescriptor `domin=coil, domout=0` (PEC) instead of
   `domin=coil, domout=air` (interface). FEM L is wrong by ~+150% for
   coil-only models, and the SIBC Robin BC has nothing to attach to for
   workpiece-bearing models.
2. **Add** SIBC Robin BC on internal workpiece surface:
   ```python
   # Robin BC: (jw/Z_s) * A_t * v_t on internal wp_surface
   a += (1j * omega / Z_s) * u.Trace() * v.Trace() * ds("sibc")
   # SIBC = Robin BC. Conductor interior not solved (hole approach).
   # Validated 2026-04-14: L<1%, P<2% for Cu/Steel/Al (2D axisym Kelvin).
   ```
3. **Add** Kelvin exterior sphere (replaces `air` Dirichlet BC):
   ```python
   # Kelvin: nu' = nu0 * (r'/R)^2 in exterior sphere (3D HCurl A;
   # radia.kelvin_source.kelvin_nu_factor_3d_cf). (r'/R)^4 is only the
   # 2D cylindrical nu_zz slot, never the 3D sphere.
   # Periodic BC: inner sphere <-> outer sphere
   ```
4. **Compute** H_t from SIBC relation (NOT from curl(A)):
   ```python
   # CORRECT: H_t = |J_s| = |jw/Z_s| * |A_t| (surface current for ESIM)
   At_sq = sum(gfu[i].real**2 + gfu[i].imag**2 for i in range(3))
   At_rms = sqrt(Integrate(At_sq, mesh, BND, definedon=wp_region) / A_wp)
   H_t = abs(1j * omega / Z_s) * At_rms
   # WRONG: curl(A)/mu0 gives total H (incident + scattered), NOT surface current
   sol = esim_solver.solve(H_t)
   Q_surface = sol['P_prime']  # W/m^2 (surface density, not volume)
   ```

### SIBC = Robin BC; Conductor Interior Not Solved (2026-04-14)

**SIBC is a Robin BC on the conductor surface.  Do NOT solve inside.**

The correct approach is the **hole approach**: subtract the workpiece
from the mesh and apply the Robin BC on the hole boundary.

For 3D HCurl: `a += (jw/Z_s) * u.Trace() * v.Trace() * ds("sibc")`
For 2D H1/phi: `a += (jw/Z_s) / r * u * v * ds("wp_bnd")`

**Validated (2026-04-14, 2D axisymmetric Kelvin)**:
Full-resolution (eddy currents resolved at delta/5) vs SIBC (hole + Robin):

| Material | R/delta | L error | P error |
|----------|---------|---------|---------|
| Copper   | 3.8-84.6 | < 0.3% | < 1% |
| Steel    | 7.0-157 | < 0.4% | < 2% |
| Aluminum | 2.9-65.7 | < 0.8% | < 2% |

Z_s for solid cylinder: `rho * gamma * I1(ga)/I0(ga)`,
`gamma = sqrt(jw * mu_r * mu_0 * sigma)` (cylindrical Bessel, not Dowell).

Script: `validation_test/eddy_current_analytical_validation/reference_2d_axisym.py`

### Karl Iteration (for nonlinear BH)

For steel with nonlinear B-H curve, Z_s depends on H_t. Linear materials
(copper, aluminum) have constant Z_s and converge in 1 iteration (no need
for Karl iteration).

```
1. Initial Z_s from ESIM at estimated H_t
2. Solve FEM with Robin BC -> get A on workpiece boundary
3. H_t = |jw/Z_s| * |A_t|_rms  (from SIBC relation, NOT curl(A))
4. Update Z_s from ESIM cell problem at new H_t
5. Under-relaxation: Z_s = 0.5*Z_s_old + 0.5*Z_s_new
6. Repeat until |dZ_s/Z_s| < tol (typically 4-5 iterations for steel)
```

### ESIM Geometry: Slab vs Cylinder

The ESIM cell problem ODE depends on geometry:

| Geometry | ODE | Analytical (linear) | Use |
|----------|-----|---------------------|-----|
| `slab` | rho*d^2H/dz^2 + jw*mu*H = 0 | cosh(gamma*(a-z))/cosh(gamma*a) | Flat plates |
| `cylinder` | (rho/r)*d/dr(r*dH/dr) + jw*mu*H = 0 | I0(gamma*r)/I0(gamma*R) | Cylindrical workpieces |

Impact on Z_s (mu_r=100, sigma=2e6, R=10mm):

| delta/R | |Z_s| diff | P' diff | When slab approx OK |
|---------|-----------|---------|---------------------|
| 0.04 (steel 7kHz) | -1% | -2% | Yes (thin skin) |
| 0.21 (copper 1kHz) | -5% | -11% | No (curvature matters) |

Rule: when delta/R < 0.1, slab and cylinder give <2% difference.
For delta/R > 0.1, use the correct geometry.

Both `fem_esim_3d.py` and `efie_sibc.py` support `--geometry cylinder|slab`.

### Per-Element Curvature (Future Extension)

For general curved surfaces (not just cylinders), each surface element has
local principal curvatures kappa_1, kappa_2. The ESIM cell problem can use
the local effective radius R_eff = 1/kappa_max:

```python
# Future: per-element ESIM with local curvature
for el in mesh.Elements(BND):
    # Get local curvature from NGSolve GetTrafo (high-order element mapping)
    trafo = mesh.GetTrafo(el)
    # Compute principal curvatures from Jacobian -> second fundamental form
    kappa = compute_principal_curvatures(trafo)
    R_local = 1.0 / max(abs(kappa[0]), abs(kappa[1]), 1e-10)
    geometry = 'slab' if R_local > 100 * delta else 'cylinder'
    esim_local = ESIMFiniteSlabSolver(half_thickness=R_local, ..., geometry=geometry)
```

**NGSolve advantage**: a mesh curved while its valid generating geometry is
attached provides high-order geometry via `GetTrafo`. From the element Jacobian,
the second fundamental form (Weingarten map) gives principal curvatures at any
point. For an imported `.vol`, consume the baked-in curved elements as-is; never
invoke post-load `Curve(3)` as a substitute for re-exporting the mesh.

### Accuracy at Target Conditions

Steel 7 kHz induction heating (R_wp=10mm):
- mu_r ~ 500 -> delta ~ 0.07 mm -> xi ~ 140
- 0th-order SIBC (flat slab) curvature error: **< 0.4%**
- ESIM handles nonlinear BH (Dowell/linear SIBC cannot)
- Boundary layer mesh for skin effect: **NOT needed**

### EM -> q_surf -> Thermal Workflow

The Radia Induction Heating Simulink application runs the EM and thermal
solvers as separate blocks:

- **EM** (`calc_inductance.py` / `calc_fem_kelvin.py` /
  `calc_fem_coilmesh.py`) emits a workpiece surface power density
  `q_surf` (`*_qsurf.sol`).
- **Thermal** (`calc_heat.py` / `calc_heat_axisym.py`) consumes the
  workpiece `.vol` + `q_surf` and solves the heat equation.

All EM solvers read a pre-meshed `.vol` (no auto air-mesh generation).
"""


INDUCTION_HEATING_AV_COIL_SIGMA = """
# A-V formulation for coil eddy currents (--coil-sigma)

## When to use

- **FEM truth / static benchmark** — validating PEEC against resolved skin.
- Coil is part of the .vol with source/sink sidesets.
- Skin-depth-resolved coil mesh (wire surface <= delta/3, e.g. 0.3 mm at 7kHz Cu).

## When NOT to use

- Rotating/moving workpiece (requires re-meshing coil + WP each step — expensive).
- Design sweeps where coil geometry varies (PEEC is faster).

## Formulation

A (vector potential) + phi (scalar potential in conductor) compound space.
The coil field is `E = -s * (A + grad(phi))` with `s = 1j * omega`.
The volume term uses the coupled combination `A + grad(phi)` and the
workpiece uses the tangential A trace in the SIBC term.

Use `python -m radia.panels.calc_fem_coilmesh --help` for the installed
entry point. It applies the port lift, computes terminal current from
the conductor integral of `J_coil * grad(psi_n)`, and scales the solution
to the requested current. The direct solve uses `inverse="sparsecholesky"`
with the free degrees of freedom of the formulation.

For the spanning-tree gauge, the separate port nodes, gap voltage
convention, and SIBC trace constraints, see
`ih_sibc(topic="ngsolve_recipes")`. A unit scaled-potential lift has
voltage magnitude `abs(s)` when the A line integral across the designated
gap path is fixed to zero.

## Gotchas

1. **Lowest-order HCurl gradients remain with `nograds=True`.** The
   coupled A–V space still needs a gauge consistent with the ports and
   SIBC traces.
2. **phi defined only on conductor** (`definedon="coil"`). If defined globally,
   the unconstrained phi in air breaks the solution (A becomes pure gradient,
   curl(A)=0 trivial solution).
3. **Periodic Kelvin requires separate validation**. Pure-HCurl AMS results
   do not certify this compound A-V system. Check finite values and true
   residuals for the chosen direct solve; never silently change the operator.
4. **A GND nodeset constrains a scalar point value.** HCurl edge degrees of
   freedom require the formulation's own gauge constraints. A small mass
   penalty can leave a nearly singular system and change normalized results.
5. **Mesh requirement**: coil surface elements < delta/3 (typically 0.3 mm
   for Cu at 7 kHz). Use `surface with length < X interval N` in Cubit.
6. **Check the port drive before normalizing.** Inspect the true linear
   residual and the raw terminal current. Verify inductance and both
   losses across admissible trees with the gap path fixed.

## Reference: existing A-V implementation
`radia.panels.calc_fem_coilmesh` implements the port-driven coil route.
"""

INDUCTION_HEATING_FAILED_APPROACHES = """
# Failed approaches for PEEC coupling to FEM workpiece (DO NOT retry)

These were tried during 2026-04-18 skin-depth comparison and shown to fail.
Documented so future sessions do not waste effort on broken formulations.

## 1. PEEC reduced-A FEM (HCurl projection of Biot-Savart A_s)

**Strategy**: sample A_s at mesh vertices via Biot-Savart, project into HCurl,
use `-(nu - nu_0)*curl(A_s_proj)*curl(v)` as RHS for scattered field A_r.
A_total = A_r + A_s on SIBC drives H_t.

**Why it fails**: HCurl projection of vertex-interpolated A_s has curl errors
up to 5x in the Kelvin domain. This is a "variational crime" — the
cancellation of J_s in the reduced-A weak form assumes
`curl(nu_0 curl A_s_proj) = J_s` exactly, which projection breaks.

**Symptoms**: L ~ 10x too high (W_mag blown up by Kelvin curl error).
H_t ~ 20x too small (Robin RHS `-robin*A_s` drives A_r to cancel A_s on SIBC).

Even computing B_s directly via `h_segments_batch` at vertices (avoiding
curl of projected A_s) did not fix H_t — the Robin surface term structure
itself drives A_total -> 0 on SIBC when A_s is weaker than the volume-current
reference.

## 2. PEEC line-integral source FEM (delta-function line source in RHS)

**Strategy**: assemble `f[dof_j] = sum_k I_k * int_path phi_j(x).dl` via
per-quadrature-point HCurl basis evaluation. Each filament quadrature point
finds containing element via `mesh(*pt).nr`, evaluates all 6 edge basis
functions via scratch GridFunction, adds weighted contribution.

**Why it fails partially**: Delta-function line source misses near-field
contribution between mesh vertices. For geometries where coil-WP gap is
similar to filament spacing (2mm gap, 1.5mm between filaments), the
near-field at workpiece is under-represented.

**Symptoms**: H_t = 7.77 A/m vs A-V truth 10.36 (-25%). Worse with more
filaments (nwinc=5: H_t=7.46) because phase cancellation between filaments
amplifies the miss. L is mesh-dependent (ln(1/h) divergence) — must use
PEEC Z_port for L, not W_mag.

**Implementation detail**: Works (no NaN, sensible order of magnitude),
just not accurate enough for validation-grade physics.

## 3. T0-modulated PEEC volume source (works but requires coil mesh)

**Strategy**: use T0 spatial distribution `J_T0 = (I/Phi)*grad(T0)` as
volume current pattern, multiply by PEEC filament current ratio
`I_k_nearest / I_avg` in each coil element. Solve FEM with modulated
volume source.

**Result**: H_t -0.9%, P -1.8% vs A-V truth — **accurate**. BUT requires
coil volume mesh. Defeats the purpose of PEEC (mesh-free coil for rotating
WP). Noted as a valid static-benchmark alternative to A-V but not
production-useful.

## 4. PEEC-modulated surface current -> compute_phi_inc_from_surface_J

**Strategy**: EFIE DC surface J (spatial pattern) multiplied by filament
`I_k/I_avg` (AC redistribution) at each panel. Feed modulated complex J
into `compute_phi_inc_from_surface_J`.

**Why it fails**: `compute_phi_inc_from_surface_J` is real-only.  It
historically cast complex J_vecs to real silently (discarding phase --
the modulation was destroyed); since 2026-07-02 it raises TypeError on
complex input instead.  Callers must bridge Re and Im in two separate
calls and combine `phi_re + 1j*phi_im` (as calc_inductance does for the
complex impedance-EFIE coil current).

**Resolution**: Direct filament -> `compute_phi_inc_from_filaments` is
simpler and already handles complex currents correctly. This IS the PEEC+BEM
production path. No modulation of surface panels needed.
"""


INDUCTION_HEATING_SOLVER_CAPACITY = """
# IH solver selection and capacity

For nonperiodic HCurl FEM-SIBC, calc_fem_kelvin --solver auto selects AMS at p=1
and BDDC+AMS at p=2/3, both with COCR; the BDDC route uses an edge-only wirebasket
and three AMS coarse cycles. Periodic Kelvin auto selects sparsecholesky;
explicit `ams` or `bddc` requests fail. Compound A-V requires explicit
sparsecholesky. Scalar thermal H1 does not inherit HCurl AMS support.
The optional IH iccg branch uses the complex-symmetric setting, not the native
Hermitian ICCG contract: it must converge and pass the original true residual.
It is not an automatic fallback or a blanket recommendation for eddy current.
""" + SOLVER_CAPACITY

def get_induction_heating_documentation(topic: str = "all") -> str:
    """Return induction heating simulation documentation by topic."""
    topics = {
        "overview": INDUCTION_HEATING_OVERVIEW,
        "gmsh_mesh": INDUCTION_HEATING_GMSH_MESH,
        "eddy_current": INDUCTION_HEATING_EDDY_CURRENT,
        "thermal": INDUCTION_HEATING_THERMAL,
        "rotating": INDUCTION_HEATING_ROTATING,
        "postprocess": INDUCTION_HEATING_POSTPROCESS,
        "pitfalls": INDUCTION_HEATING_PITFALLS,
        "esim_kelvin": INDUCTION_HEATING_ESIM_KELVIN,
        "av_coil_sigma": INDUCTION_HEATING_AV_COIL_SIGMA,
        "failed_approaches": INDUCTION_HEATING_FAILED_APPROACHES,
    }

    topic = topic.lower().strip()
    if topic in ("solvers", "capacity", "memory", "diagnostics", "bddc"):
        return INDUCTION_HEATING_SOLVER_CAPACITY
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
