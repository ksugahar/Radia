"""
Induction heating simulation knowledge base for the Radia MCP server.

Covers: complete EM-to-thermal workflow for induction heating using NGSolve,
including Gmsh mesh loading, A-Phi eddy current formulation, Joule heat
computation, transient thermal analysis, rotating workpiece, VTK output,
and post-processing patterns.

Public docs/notebooks:
  - docs/induction_heating/public_demo.ipynb -- human-facing ESIM/BEM/WPT
    induction-heating entry point with saved outputs and JSON sidecar.
  - docs/induction_heating/induction_heating_demo_showcase.ipynb -- saved
    promotion notebook for the closed public ESIM/WPT/RWG demo scripts.
  - docs/induction_heating/induction_heating_examples_catalog.ipynb -- full
    source/hash catalog for IH docs/API/validation migration lanes.
  - docs/ih_esim_benchmark/esim_showcase.ipynb -- nonlinear ESIM benchmark
    figures/results.

Sources:
  - Internal production induction-heating toymodel notes
    (K. Sugahara's production simulations; machine-local path omitted)
  - https://docu.ngsolve.org/latest/i-tutorials/
  - https://forum.ngsolve.org/
"""

INDUCTION_HEATING_OVERVIEW = """
# Induction Heating Simulation with NGSolve

## Physics Overview

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

Gmsh meshes use "physical groups" to define material regions (volumes)
and boundary conditions (surfaces). These map to NGSolve materials/boundaries.

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

### Mesh from Cubit (export netgen)

The recommended path is `export netgen` (Cubit plugin command),
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
cubit.cmd(f'export netgen "{out}" order 2 overwrite')
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

**CRITICAL**: Both `fem_esim_kelvin.py` (2D axi) and `fem_esim_3d.py` (3D)
MUST call `mesh.Curve(>=2)` after `Mesh()` if FES order >= 1 and the
geometry contains curved entities (circles, torus, sphere, cylinder).
Missing `mesh.Curve()` causes a **systematic L under-prediction of
~10 %**, because circular coil cross-sections and Kelvin arcs collapse
to polygons.

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

**Bottom line for current 3D IH panel code** (`fem_esim_3d.py`,
`calc_fem_kelvin.py`):

- Start with `HCurl order=1 + mesh.Curve(max(order+1, 2))`
- If L seems too low (-5 % or worse) and P looks OK, **raise Curve
  to 3** instead of raising FES order. Curve is cheap; FES doubles
  ndof and can introduce spurious modes.
- HCurl order=2 costs ~60x more on this test case (ndof 83k -> 333k,
  23s -> 23min) and can OVERSHOOT L and WORSEN P at the same Curve
  level. Don't raise FES order without raising Curve too.

### Why 2D is less sensitive than 3D

2D axisym has only planar curves (coil circle + Kelvin arc). H1 scalar
elements with `Curve(order=2 or 3)` match the geometry accurately
even at moderate mesh resolution. 3D adds sphere + torus in space --
strong Gaussian curvature that needs better surface approximation.

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

## Solver Comparison for A-Phi

| Preconditioner | Solver | Use Case | Notes |
|----------------|--------|----------|-------|
| BDDC | CG | Standard, symmetric-like | Fastest for well-conditioned |
| BDDC | GMRes | Non-symmetric systems | More robust, slightly slower |
| local | CG | Small problems | Simple but slow for large DOF |
| (none) | Direct (PARDISO) | Debug/small problems | Exact but O(n^2) memory |

**Typical performance** (mesh with ~222k elements, ~1.7M DOFs):
- Matrix assembly: ~400 sec
- BDDC+CG solve: ~90 sec (tol=1e-10)
- VTK output: ~16 sec
"""

INDUCTION_HEATING_THERMAL = """
# Transient Thermal Analysis for Induction Heating

## Production path: IH Simulink block plus temporary notebook comparison

Heat analysis is a **Method choice** in the Radia Simulink Induction Heating
block and `IHDesignSpec`; the temporary comparison notebook exposes the same contract,
alongside the existing PEEC-Inductance / PEEC-BEM / FEM-Kelvin /
FEM-coilmesh methods.  Pick one of three Thermal Method choices
(v4.63.0+):

  - ``Thermal: 3D static (no rotation)`` -- 3D solver with
    rotation_rpm = 0.  Use for static one-shot heat-up.
  - ``Thermal: 3D + rotation (q_surf re-sampled per step)`` --
    3D solver with workpiece rotation (v4.58.0 qsurf re-projection
    each timestep).  Use for non-axisymmetric workpieces / coils.
  - ``Thermal: 2D axisymmetric (rotation implicit)`` -- axisym
    solver (10-100x faster); rotation is implicit in the axisym
    assumption.  Use for rotationally-symmetric workpieces.

The Method dropdown encodes the (mesh_type, rotation_rpm) pair;
the embedded HeatPanel's individual fields are hidden.

Pre-v4.63.0 single-entry path (kept for context): pick **Method
→ "Thermal (heat transfer from
saved q_surf .sol)"**.

Workflow (single window):

1. Run an EM solve method that emits ``qsurf.sol`` (PEEC+BEM,
   BEM-A+BEM, PEEC+FEM+Kelvin, or FEM-full).  Note: PEEC+BEM and
   BEM-A+BEM gained the ``qsurf.sol`` output in **radia 4.65.0**.
   The path has been hardened across four follow-up releases:

     - **4.66.0**: parent-vol_mesh em_vol (was 2D surface mesh -->
       cross-mesh transfer mirror-flipped) + signed energy fix.
     - **4.67.0**: surface-path phi_inc + triangle-wise gradient.
       (Energy magnitude regressed to 1.92 W on 3turnCoil_work,
       below the weak-coupling lower bound 5.26 W -- the BIE was
       under-coupling when fed a different-gauge phi_inc.)
     - **4.68.0**: BIE-calibrated direct Biot-Savart for spatial.
       The qsurf spatial pattern is computed from
       q(x) = (P_wp / ∫|H_t_inc|² dS) * |H_t_inc(x)|² where
       H_t_inc is the tangential Biot-Savart field at each wp
       surface vertex (curl-free, no topological multivalued
       branch cut).  The BIE provides only P_wp (its trusted
       global integral).  Result: ∫q dS = P_wp exactly (energy
       preserved) AND spatial peak matches the true |H_t|² peak.
     - **4.69.0**: P2 BIE outward-orientation guard.  The Lagrange-
       P2 path (``--h1-order 2``, curved-Tri6 geometry) was
       producing a ~1500x P_wp discrepancy vs the P1 reference
       because ``bem/sibc_hacapk.py::extract_surface_p2_lagrange``
       used the vol_mesh's natural BND triangle orientation -- on
       a workpiece-as-hole sibc this points INWARD toward the hole,
       flipping the sign of the Galerkin trace(DL) integral.  Fixed
       by mirroring the P1 path's outward-flip logic against the wp
       centroid (commit 57381c78, "Release v4.69.0 / radia-mcp-
       v0.64.0 -- P2 BIE outward orientation fix").  After the fix,
       P1 = 6.1386 W vs P2 = 6.1565 W on 3turnCoil_work (0.29%
       agreement); the curved-P2 path is now production-ready for
       basis_order=2 BEM workpieces.

   Use 4.68.0+ for thermal analysis on BEM-derived qsurf; use
   4.69.0+ if you want to enable ``--h1-order 2`` curved BEM.
   Verification on 3turnCoil_work (PEEC-BEM, 150 kHz, 100 A):
     z of max q   = +11.25 mm (matches Biot-Savart peak exactly)
     ∫ q dS       = 6.1386 W (= BIE P_wp)
     peak / end   = 1465x (strong concentration under coil)
2. Click the **"Run thermal..."** chain button on the action row.
   This SWITCHES the method dropdown to "Thermal" and pre-fills the
   embedded thermal panel's ``qsurf_sol`` + ``em_vol`` fields.
3. Pick the workpiece thermal ``wp.vol`` (a SEPARATE mesh from the
   EM ``.vol``; the EM mesh treats the workpiece as a hole with a
   SIBC face, the thermal mesh treats it as a real solid).
4. Set material / convection / time scheme / rotation_rpm.
5. Click Run.

The pre-4.59.0 standalone ``radia-heat`` (HeatWindow) launcher was
REMOVED in radia 4.62.0.  Heat analysis is the Method = "Thermal"
choice in ``radia-ih`` exclusively.  The HeatPanel sub-widget lives
at ``radia._heat_panel`` as an internal implementation detail of
the IH panel; no public CLI replacement for ``radia-heat`` exists
(launch ``radia-ih`` and pick the Thermal method instead).
Programmatic imports: replace ``from radia.radia_heat import ...``
with ``from radia._heat_panel import ...``.

## ``.sol + .vol`` strict contract (4.58.0+)

NGSolve ``.sol`` files are coefficient vectors only -- no embedded
mesh, no FES-order header.  Both files MUST be passed:

```bash
# Direct CLI (calc_heat.py):
python -m radia.panels.calc_heat \\
    --wp-vol workpiece_thermal.vol \\
    --qsurf-sol  <stem>_qsurf.sol \\  # REQUIRED for spatial mode
    --em-vol     <stem>_fem.vol   \\  # REQUIRED (no auto-locate)
    --qsurf-order 1                \\  # MUST equal EM fes_order
    --material steel --dt 0.5 --t-end 5.0 \\
    --rotation-rpm 12                  # see ``rotating`` topic
    # --surface-label is OPTIONAL since v4.74.0: empty = all BND.
    # Pass a specific name only when the workpiece has MULTIPLE BND
    # sidesets and heating + convection should be restricted to a
    # subset.  The panel auto-fills the sole BND label when the .vol
    # has exactly one (e.g. 'sibc'), so users typically don't need
    # to type it.
```

`--qsurf-order` must match the EM solve's `--fes-order` exactly;
mismatch reads garbage silently (no NGSolve error).  Defaults: both 1.

The HeatPanel's Run button stays disabled until both .sol AND .vol
exist on disk.  Auto-locate by filename convention was REMOVED in
4.58.0 per CLAUDE.md "No Fallbacks" (was a silent footgun when the
user renamed or moved files).

## ``.sol`` data model: H1 GridFunction with boundary-only DOFs

`calc_fem_kelvin.py` saves q_surf as:

```python
fes_q = H1(mesh, order=fes_order)        # VOLUME H1 on the EM mesh
gf_q = GridFunction(fes_q)
gf_q.vec[:] = 0                          # interior DOFs stay 0
gf_q.Set(q_surf_cf, definedon=wp_region) # only workpiece-boundary DOFs touched
gf_q.Save(qsurf_sol_path)
```

The thermal solver re-projects this onto the wp-mesh surface
vertices by point evaluation.  When wp-mesh and em-mesh share the
same physical workpiece geometry, vertices land on the EM boundary
faces and H1 continuity guarantees the boundary node values are
recovered exactly.  See ``radia_mcp.radia_ngsolve.ngsolve`` Section
18c for the broader pattern (surface-restricted H1 GF save/load).

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
| ``T_sol_file`` | absolute path to the final-T NGSolve `.sol` |
| ``heat_vol_file`` | absolute path to the companion `.vol`.  Empty when no separate companion was written (then re-use the `--wp-vol` input as the companion). |
| ``msh_file`` | GMSH `.msh v4.1` (T_C + q_surf fields) when `--msh-output` was set |
| ``vtu_files`` | per-step `.vtu` paths when `--vtu-prefix` was set |
| ``csv_file`` | probe history CSV when both `--probe-point` + `--csv-output` were set |

Naming convention:

* With `--msh-output FILE.msh`: writes `FILE_T.sol` +
  `FILE_heat.vol` (fresh companion mesh) alongside `FILE.msh`.
* Without `--msh-output`: writes `<wp-stem>_heat_T.sol` next to
  `wp.vol`.  No separate companion mesh; reload directly against
  the same `wp.vol` the solve consumed.

Reload pattern:

```python
from ngsolve import Mesh, H1, GridFunction
wp_mesh = Mesh("workpiece_thermal.vol")  # or <msh-stem>_heat.vol
fes_T   = H1(wp_mesh, order=1)           # MUST match solve's --fes-order
gfT     = GridFunction(fes_T)
gfT.Load("workpiece_thermal_heat_T.sol") # or <msh-stem>_T.sol
T_at = float(gfT(wp_mesh(x, y, z)))      # sample at body point
```

Same three contracts as the qsurf side (see "Strict .sol + .vol
contract" above): the `.sol` is a raw coefficient vector, the
`.vol` carries the mesh, the FES order must match.

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
T_ss = (q/(eps*sigma) + T_ext_K^4)^(1/4)  (tests/panels/test_heat_radiation.py,
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

## Production path: pass `--rotation-rpm` to `calc_heat.py`

Since radia v4.58.0 the 3D thermal solver implements workpiece rotation
natively.  Users do NOT write rotation code by hand; they just supply
a non-zero `--rotation-rpm`:

```bash
python -m radia.panels.calc_heat \\
    --wp-vol     workpiece_thermal.vol \\
    --surface-label sibc \\
    --qsurf-sol  ih_em_qsurf.sol \\
    --em-vol     ih_em_fem.vol \\
    --material   steel \\
    --dt 0.5 --t-end 5.0 \\
    --rotation-rpm 12        # <-- workpiece spins around z at 12 rpm
```

Equivalent panel knob: the `Rotation [rpm]` field in `radia-ih`
(Layer 3 Cubit panel, Method = "Thermal" section; the field is
implemented in `radia._heat_panel`).  Default 0 (stationary).

## How calc_heat.py implements rotation

The implementation is a **per-step re-projection of q_surf on the body
frame** (NOT mesh deformation).  Cheaper than `mesh.SetDeformation`
because:

* The thermal mesh / FES / mass matrix / stiffness matrix are held FIXED
  for the entire simulation -- they're built once outside the time loop
  and the linear solver factorisation (`mstar = M + theta*dt*K`) is
  cached.  Only the LinearForm RHS reassembles each step (it already
  has to, for the convection term).
* Re-projection cost = N surface vertex point-evaluations of the EM-frame
  qsurf GridFunction.  Identical to the one-shot projection done at
  startup, so per-step overhead = O(one initial projection).

The body-frame projection at angle `theta = omega * t`:

```python
# At each timestep, sample q_em at the world-frame coord that
# corresponds to body coord (xb, yb, zb) after rotation by +theta:
c, s = math.cos(theta), math.sin(theta)
for vnr, (xb, yb, zb) in zip(surf_vnrs, surf_xyz):
    xw = xb*c - yb*s
    yw = xb*s + yb*c
    em_mip = em_mesh(xw, yw, zb)        # world-frame lookup
    val = gf_q_em(em_mip)               # q_em evaluated there
    gf_wp_q.vec.FV()[vnr] = float(val)  # update body-frame GF in place
```

The rebuilt `gf_wp_q` is the same GridFunction returned by
`_build_qsurf_cf(...)`; it's updated **in place** so q_cf in the
LinearForm picks up the new values automatically.

## Why NOT mesh deformation

`mesh.SetDeformation(...)` deforms the FE space at each step, which
invalidates the cached mass / stiffness / factorisation.  For a 5-second
simulation with dt=0.5 (10 steps), that's 10x the full assemble + LU
factor cost (~30s each on a 200k-DOF cylinder mesh = ~5 minutes wasted).
Re-projection of q_surf is ~10ms per step on the same mesh.

## .sol contract for rotation users

The rotation path consumes the spatial qsurf path (`--qsurf-sol +
--em-vol`).  Uniform `--q-uniform` is rotation-invariant and the panel
warns when `--rotation-rpm > 0` is paired with it.

Three CONTRACTS that must hold between the EM solve (calc_fem_kelvin.py)
and the thermal solve (calc_heat.py):

1. **`.sol` is mesh-free**.  NGSolve's `GridFunction.Save()` writes
   the raw coefficient vector (no mesh, no header).  The thermal
   solver MUST be told the EM mesh path via `--em-vol`.  Pass the
   `*_fem.vol` file that `calc_fem_kelvin.py` saves alongside
   `*_qsurf.sol`.

2. **`qsurf_order` MUST equal the EM `fes_order`**.  The thermal solver
   rebuilds `H1(em_mesh, order=qsurf_order)` and loads the .sol into
   it.  If the orders disagree the coefficient vector lands in a
   space with a different DOF count and the loaded field is garbage
   (no NGSolve-level error, silent corruption).  Default for both is
   1.  When the EM solve uses `--fes-order 2`, pass `--qsurf-order 2`
   to calc_heat too.

3. **q_surf is a volume H1 GF with non-zero values ONLY on the
   workpiece boundary**.  calc_fem_kelvin does
   `gf_q = GridFunction(H1(mesh, order=fes_order));
    gf_q.Set(q_surf_cf, definedon=wp_region)`.
   The rest of the volume is 0.  The thermal solver point-samples
   this GF at workpiece surface vertices -- H1 continuity guarantees
   the boundary node values are recovered exactly.  Wp surface
   vertices that fall OUTSIDE the EM mesh (mesh mismatch) are set to
   0 and a count is reported in the run log.

## Test coverage

* `tests/panels/test_heat_rotation.py` -- unit tests for the rotation
  projection math on a synthetic unit-cube setup with q_em(x,y,z)=x.
  At theta=0 / pi / pi/2 the body point near (+1,0,0) reads
  q ~ 1 / -1 / 0 respectively.

* `tests/panels/test_heat_chain_golden.py` -- end-to-end chain test
  (EM solve -> qsurf.sol -> thermal solve).  Slow-marked.
   Use ParaView to animate the sequence.

5. **Full rotation**: At 12 rpm, one rotation = 5 seconds. Choose t_end
   to capture the desired number of rotations.
"""

INDUCTION_HEATING_POSTPROCESS = """
# Post-Processing for Induction Heating

## VTK Output for ParaView

```python
from ngsolve import *

# ============================================================
# Basic VTK output (all fields)
# ============================================================
vtk = VTKOutput(mesh,
                coefs=[B.real, B.imag, J.real, J.imag, Q],
                names=["B_re", "B_im", "J_re", "J_im", "Q"],
                filename="eddy_current_result",
                subdivision=1,    # Subdivide elements for smooth plot
                legacy=False)     # Use VTK XML format (.vtu)
vtk.Do()

# ============================================================
# VTK output with material ID
# ============================================================
# Create integer field for material identification
mat_names = mesh.GetMaterials()
mat_id_dict = {mat: i+1 for i, mat in enumerate(mat_names)}
mat_id = CoefficientFunction(
    [mat_id_dict[mat] for mat in mat_names]
)

vtk = VTKOutput(mesh,
                coefs=[mat_id, B.real, Q],
                names=["MaterialID", "B_real", "Q"],
                filename="result_with_id",
                subdivision=1, legacy=False)
vtk.Do()

# ============================================================
# VTK output for thermal time steps (.pvd animation)
# ============================================================
# Create VTKOutput ONCE, call Do(time=t) for .pvd time series
vtk = VTKOutput(mesh,
                coefs=[gfT],
                names=["Temperature"],
                filename="thermal_series",
                subdivision=1, legacy=False)

for step in range(n_steps):
    # ... solve time step ...
    vtk.Do(time=t)  # Writes thermal_series_N.vtu + updates .pvd

# Open thermal_series.pvd in ParaView for animation playback
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

### Animated GIF from VTK Outputs

```python
# animated_gif.py: Create animated GIF from PNG frames
from PIL import Image
import glob

# Assuming ParaView has exported PNG frames
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

## 9. MODERATE: Missing TaskManager for Parallel Assembly

```python
# SLOW: Serial assembly
a.Assemble()  # Uses single thread

# FAST: Parallel assembly
with TaskManager():
    a.Assemble()  # Uses all CPU cores
```

## 10. LOW: Point Evaluation Syntax

```python
# WRONG: Direct coordinate (only works for scalar fields)
T = gfT(0.03, 0, 0)  # May fail for complex meshes

# CORRECT: mesh() mapping first
T = gfT(mesh(0.03, 0, 0))

# For vector field components
Bz = B[2](mesh(0, 0, 0.03))  # B_z component
```

## 11. LOW: VTK Legacy Format

```python
# OLD: Legacy VTK format (.vtk)
vtk = VTKOutput(mesh, ..., legacy=True)  # Large files, old format

# RECOMMENDED: XML VTK format (.vtu)
vtk = VTKOutput(mesh, ..., legacy=False)  # Smaller, modern format
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
Q = 0.5 * sigma * InnerProduct(E, E).real                 # Joule heat [W/m^3]
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
Q = 0.5 * sigma * InnerProduct(E, E).real  # from EM solve
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
   # Kelvin: nu' = nu0 * (r'/a)^4 in exterior sphere (3D)
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

**NGSolve advantage**: `mesh.Curve(order)` provides exact high-order geometry
via `GetTrafo`. From the element Jacobian, the second fundamental form
(Weingarten map) gives exact principal curvatures at any point. This is
unavailable with flat (order=1) elements — another reason to use Curve(3)+.

### Accuracy at Target Conditions

Steel 7 kHz induction heating (R_wp=10mm):
- mu_r ~ 500 -> delta ~ 0.07 mm -> xi ~ 140
- 0th-order SIBC (flat slab) curvature error: **< 0.4%**
- ESIM handles nonlinear BH (Dowell/linear SIBC cannot)
- Boundary layer mesh for skin effect: **NOT needed**

### EM -> q_surf -> Thermal Workflow

The radia_ih panel runs the EM solve and the thermal solve as separate
modes:

- **EM** (`calc_inductance.py` / `calc_fem_kelvin.py` /
  `calc_fem_coilmesh.py`) emits a workpiece surface power density
  `q_surf` (`*_qsurf.sol`).
- **Thermal** (`calc_heat.py` / `calc_heat_axisym.py`) consumes the
  workpiece `.vol` + `q_surf` and solves the heat equation.

All EM solvers read a pre-meshed `.vol` (no auto air-mesh generation).
"""


INDUCTION_HEATING_PEEC_BEM_SIBC = """
# PEEC (coil filament) + BEM-SIBC (workpiece) — production path for rotating WP

## POLICY: BEM-A = impedance-EFIE, solved by COCR or HACApK-COCR

**POLICY (2026-07-03, Sugahara):** the BEM-A coil solver
(``--coil-solver bem-a``, ``radia.bem.coil_inductance_ngsolve``) uses the
**impedance-EFIE** formulation and ONLY that, and it is solved by
**COCR or HACApK-COCR** -- the two ``--coil-saddle-solver`` choices:

- ``cocr``        : div-free loop reduction + complex-symmetric COCR, dense
                    SL matvec (``auto`` picks this for large systems).
- ``hacapk_cocr`` : the same COCR with the HACApKBEMManager-compressed
                    O(N log N) SL matvec (accuracy ~3e-7; identical R/L).

(Dense ``lu`` stays the small-N direct option; ``gmres`` / ``minres`` are kept
only for comparison -- unpreconditioned GMRES stalls on the indefinite AC
saddle.)  The Leontovich surface impedance ``Z_s = (1+j)/(sigma*delta)`` sits
INSIDE the saddle system::

    [ jw*mu0*SL + Z_s*M   D^T ] [J]   [0]
    [ D                    0  ] [p] = [g]      (complex, omega > 0)

so the recovered J is the FINITE-IMPEDANCE surface current, and
``R = Re(Z_s) * (J^H M J)`` (Leontovich SIBC dissipation),
``L = mu_0 * (J^H SL J)`` (EXTERNAL inductance -- the internal surface
reactance ``Im(Z_s)*(J^H M J)/omega`` is deliberately NOT folded into L).
At ``omega == 0`` it reduces to the real vacuum-L saddle (R = 0).

**The historical PEC post-hoc formulation is REMOVED and must not return.**
It solved a perfect-conductor saddle then integrated ``R = Re(Z_s)*|J|^2 dS``
afterwards; the PEC J concentrates singularly at near-contact turn gaps /
edges where it varies below the skin depth and the Leontovich integral
breaks down, over-estimating R ~3x on tightly-wound coils (kubota 3-turn
pancake: PEC 15.14 mOhm vs the physical 4.63, the latter confirmed by
volume PEEC 3.7 / perimeter PEEC 4.5 / analytic proximity 4.8).  On smooth
geometry both agreed (isolated straight wire = closed-form Bessel to <1%),
so the removal loses nothing.

**Scalable solver = COCR** (two ``--coil-saddle-solver`` choices: ``cocr``
[dense SL matvec] and ``hacapk_cocr`` [HACApK-compressed SL matvec]; ``auto``
selects ``cocr`` for large systems: > 5000 tris AC / 10000 DC).  It reduces
the saddle to the divergence-free (loop / stream-function) subspace via the
sparse projector onto ker(D), giving the complex-symmetric operator
``Pi A11 Pi`` that COCR (Sogabe-Zhang, UNCONJUGATED inner products) solves in
~24 MESH-INDEPENDENT iterations -- EXACT vs the dense LU (dR/dL < 0.01% on the
gapped torus) and replacing BOTH the O(N^3) LU and the unpreconditioned GMRES
that stalled (~1e5 matvecs) on the indefinite AC saddle.  ``hacapk_cocr`` uses
the HACApKBEMManager-compressed O(N log N) SL matvec (accuracy ~3e-7, the
bem_sibc_solver pattern); ``cocr`` and ``hacapk_cocr`` give identical R/L.
COCR does NOT fit the raw saddle
(structural breakdown -- rhs = [0; g] lives in the constraint block so the
initial r^T A r = 0) nor the Schur complement (diverges); the div-free
reduction is what makes it fit.

## L differs from PEEC by ~4% -- BENIGN, and curve order does NOT close it

R is the accepted metric and it AGREES (impedance-EFIE 4.63 mOhm ~= perimeter
PEEC 4.5 / analytic 4.8, all in the ~4.5-5 band).  The residual L gap --
BEM-A EXTERNAL L 411.6 nH vs PEEC 430 nH (~4.3%) on the kubota 3-turn coil --
is a KNOWN, BENIGN geometry/convention difference, LEFT AS-IS (Sugahara
2026-07-03: signed off because R agrees).  It is NOT a defect, and raising the
mesh curve order does NOT close it:

- MEASURED (torus 150 kHz, fixed RT0 DOF n_J=5745): curve_order 1/2/3 give
  L=87.2702 nH and R=1.04531 mOhm BIT-IDENTICAL.  ``mesh.Curve(p)`` is a NO-OP
  in this path because the in-memory surface extraction
  (``_extract_surface_mesh_filtered``) strips the CAD/curvedelements
  association, so there is no geometry to project the mid-side nodes onto.
  Raising geometry order would require a Cubit ``export netgen ... order N``
  surface ``.vol`` loaded directly.
- Even properly applied, curve order is NOT the lever: EXTERNAL L is
  loop-dominated (turn radius, enclosed area, N^2 -- already resolved by flat
  elements), and round-wire faceting enters L only logarithmically
  (ln(8R/a)) = a sub-1% effect.  The ~4% gap is the "surface mesh != PEEC
  filament geometry" mismatch + PEEC finite n_peri + external-L(BEM) vs
  bundle-L(PEEC) convention -- none touched by curve order.
- 430 nH is NOT validated ground truth (it is PEEC's own approximation; the
  FEM A-V 3-D reference is still TBD).  If L ever needs reconciling: match the
  SAME STEP geometry, h-refine the surface mesh, then anchor to FEM A-V
  (``calc_fem_coilmesh.py``) or measurement -- NOT curve order.

Refs: ``docs/peec/VOLUME_PEEC_DESIGN.md``,
``docs/esim/R_MISMATCH_PEEC_VS_BEMA.md``,
``validation_test/bem/test_coil_bem_a_impedance_efie.py``.

## Strong coupling and genus-1 workpieces

Both coil solvers expose an iterative self-consistent path through
``--coupling-mode strong``: BEM-A uses
``radia.bem_coupled_solver.CoupledBEMSolver`` and PEEC uses
``radia.peec_coupled_bem_solver.CoupledPEECBEMSolver``.  The BEM-A
solver is normalized to unit terminal current, so the driver scales
fields by ``I`` and dissipation by ``I^2``.  The scaling contract is
locked by ``test_strong_output_scales_with_terminal_current``.

Surface winding must be globally consistent before assembling the
double-layer operator.  ``surface_mesh_extract.orient_surface_triangles``
uses face-BFS propagation and signed-volume normalization, including
the inner wall of a genus-1 tube where a centroid heuristic is invalid.

A single-valued scalar potential cannot carry net current through a
surface cut.  For a flux-linked genus-1 workpiece, the loop extension
adds the harmonic shorted-turn current.  The extension is closed by
Faraday's law and is locked by the analytic
shorted-ring golden plus a frozen-subsystem equivalence check.
Theoretical anchor: K. Sugahara, "Investigation of a Boundary Integral
Equation n x H = J_s on Torus-Shaped Perfect Conductors," IEEE Trans.
Antennas Propag. 56(3), pp. 722-727 (2008) -- the DUAL defect of the
same H^1(S) != 0 topology (there the BIE admits a spurious harmonic
null-space mode violating B.n=0, closed by one virtual-magnetic-current
DOF + a one-point constraint; here the representation LACKS the
harmonic mode, closed by one loop DOF + the Faraday row).  Per handle:
one extra DOF, one extra condition, in both formulations.

Related work (cite BOTH, the digest regardless of whether a journal
version appears): M. Schoebinger and K. Hollaus, "An Effective
Interface Approach for Multiply Connected Electromagnetic Shields,"
IEEE CEFC 2026 conference digest (Thessaloniki, June 2026) -- the FEM
thin-shell sibling of the same topology treatment: the shield becomes
a 2-D effective interface (1-D through-thickness analytic solution,
nonlinear mu_eff lookup -- the ESIM idea) and each HOLE gets one
cohomology jump unknown (constant T = s_i e_z in hole i), replacing
the earlier non-physical auxiliary conductivity inside the holes.
P. Dlotko, B. Kapidani, S. Pitassi, R. Specogna, "Fake Conductivity
or Cohomology: Which to Use When Solving Eddy Current Problems With
h-Formulations?", IEEE Trans. Magn. 55(6), 1-4 (2019) -- established
that the cohomology treatment (the route ``radia.cohomology``
implements) is the stable, physical choice.  Differentiators of the
radia route: surface-only BIE (no volume mesh, exact open boundary) +
Leontovich SIBC + the surface-H^1 period-matrix engine with
class-pure cut selection + the Faraday closure and per-run screening
diagnostics.

**Why the MISSING mode makes the workpiece heat MORE (not less).**
Counter-intuitive on first sight ("an extra current should add Joule
heat"), resolved by the PHASE.  The tube is a shorted one-turn
secondary: ``I_2 = -j omega M I_1 / (R_loop + j omega L_loop)``, whose
phase lies between quadrature (-90 deg, resistance-dominated) and
anti-phase (-180 deg, inductance-dominated).  In the
INDUCTANCE-DOMINATED regime the shorted turn approaches a lossless flux
canceller.  Since
``P = 1/2 Re(Z_s) int |H_t|^2`` and
``|H_inc + H_alpha|^2 = |H_inc|^2 + 2 Re(H_inc . H_alpha*) +
|H_alpha|^2`` with a large NEGATIVE cross term, the Lenz screening
can remove more surface |H_t|^2 than the mode's own dissipation adds.
The emitted ``P/P_frozen`` ratio and alpha phase expose whether a run
is in that screening regime.  Clamping the mode to zero (single-valued phi)
is the physics of a tube with an insulating slit around its section:
the full coil flux swings through the bore unopposed and the surface
sees the unscreened field.  The
machine-design intuition "a shorted turn overheats" belongs to the
OPPOSITE regime (omega L << R_loop, phase near -90 deg, weak screening,
the turn mostly self-heats); the alpha phase is the regime
discriminator.  The cut is oriented to positive toroidal winding before
alpha is reported, so this phase label is independent of tree orientation.

Every loop-DOF run now emits the screening diagnostics in the JSON:
``wp_loop_P_frozen_W`` / ``wp_loop_H_t_frozen_A_per_m`` (the no-mode,
"slitted-tube" values from the frozen alpha=0 sub-solve),
``wp_loop_screening_ratio`` (= P_wp / P_frozen; < 1 in the screening
regime), and ``wp_loop_regime`` ("inductive-screening" for |phase| >=
135 deg / "mixed" / "resistive-dissipative" below 105 deg), so the
physics of the correction is auditable per run.  It
requires linear SIBC, ``--wp-bem-backend intree-dense``, and (on the
weak path) ``--h1-order 1``; BOTH coupling modes take it -- the weak
path applies it to its single solve, the strong drivers apply it ONCE
on the converged Picard state (the Picard loop itself keeps the plain
solve whose L_total / Delta_L convention is retained; the alpha
back-reaction onto the coil current is not iterated).  The output
records ``wp_loop_alpha_A`` and replaces the genus caveat with
``P_wp_note``.

On every P1 route -- weak AND the strong Picard iterations -- the
incident potential is basis-determined: the solver always uses the
Laplace-Beltrami projection of the exact vertex incident field
(``SurfacePoissonPhiInc``, stiffness factorized once and reused per
Picard iteration).  It performs one batched field evaluation per
rebuild, is winding invariant, and fails loud when
``||grad_S psi + H_t,inc|| / ||H_t,inc||`` exceeds 10 percent.  The
legacy selectable P1 path-integration route and ``--wp-phi-inc`` flag
are removed; the strong solvers' former per-iteration path integration
and former ngsolve.bem dense workpiece assembly are removed with it --
the dense wp backend is now the same in-tree
Sauter-Schwab Galerkin configuration as the weak path.  Only the
Lagrange-P2 edge-node route retains path integration, as its sole
implemented reconstruction.

``--wp-loop-dof`` accepts ``auto`` and ``on`` only.  ``auto`` applies
the extension when its prerequisites hold and records
``wp_loop_dof_skip_reason`` otherwise; ``on`` keeps fail-fast
prerequisite checks.  The known-invalid ``off`` route is not selectable.
The Simulink configuration and temporary notebook workbench expose this control.

**Part 1 (DOMINANT, FIXED 2026-07-17): inconsistent surface winding.**
The hole extractor's per-triangle "centroid-outward" flip is wrong on a
genus-1 tube (the bore-wall outward normal points TOWARD the centroid),
so it flipped the entire inner wall -- 199 directed-edge conflicts --
corrupting the double-layer operator.  Fixed by
``surface_mesh_extract.orient_surface_triangles`` (face-BFS flip
propagation + per-component signed-volume outward), wired into BOTH
extractors.  No-op on consistent meshes (sphere benchmark unchanged).

**Part 2 (was +27-32 % on P_wp, +11 % on H_t): genus-1 missing loop
current -- SOLVED by ``radia.bem_loop_extension`` (2026-07-17).**
chi = V - E + F = 0 -> genus 1; the coil flux links the bore, and the
BIE's ``J_s = n x (-grad phi)`` with single-valued phi carries ZERO net
current through any cut -- the shorted-turn eddy current and its Lenz
screening are unrepresentable.  The extension adds ONE DOF alpha (the
net toroidal current): ``phi = phi_u + alpha Theta`` with Theta the
potential of a unit mid-wall ring (ray-cast from the homology cut, same
class -> single-valued on the cut-open mesh, verified +-1 jump), alpha
column = ``SL(gamma M^-1 K(Theta) - q_Theta)`` (the membrane term
cancels against Theta's own identity), closed by Faraday on the cut
loop.  BEM operators stay on the CLOSED mesh (assembling on the open
mesh poisons the regular quadrature via coincident duplicated vertices).
The analytic shorted-ring golden and frozen-subsystem equivalence test
lock the added mode and its disabled-limit behavior.  Entry points:
``radia.bem_loop_extension.solve_loop_extended(solver, phi_inc, Z_s,
omega, A_inc_fn)`` and the CLI flag ``calc_inductance.py --wp-loop-dof``
(weak AND strong coupling, linear SIBC, ``--wp-bem-backend
intree-dense``, ``--h1-order 1`` on weak; works with BOTH coil sources
-- surface panels or PEEC filaments via the exact ``A_from_filaments``,
and both strong solvers take ``loop_dof=True`` applied on the converged
Picard state).  With the flag, P_wp / H_t are replaced by the
loop-extended values, the Telegen delta_L (weak) / L_total (strong)
keep the plain-phi convention, the genus ``P_wp_caveat`` becomes a
``P_wp_note``, and ``wp_loop_alpha_A`` reports the shorted-turn current;
a built-in frozen-vs-plain cross-check refuses to report on operator
mismatch.

**Weak-vs-strong method choice: WEAK FIRST (Sugahara, 2026-07-17).**
Run the weak path as the production route -- on Takahashi it is ~10x
faster on the workpiece stage (~25 s vs 228 s coupled; the coil solve
is common) and agrees with strong to ~1% on P_wp / H_t.  The weak
run's OWN output reports the coupling strength: ``|delta_L_nH| /
L_coil_nH`` measures the coil back-reaction.  MEASURED on Takahashi
(2026-07-17 corrected; an earlier note said 1.4%, which was wrong):
delta_L = -15.5 nH -> |delta_L|/L_coil = 15% -- and weak-vs-strong
P_wp STILL agrees to ~1%, so a double-digit back-reaction ratio does
NOT by itself mandate strong for heating.  Escalate to
``--coupling-mode strong`` when the self-consistent L_total including
the workpiece magnetic-energy term (the term weak's Telegen form
drops) is itself the target quantity, or in close-coupled / high-mu /
small-gap regimes where the coil-current REDISTRIBUTION (not just its
lumped delta_L) plausibly changes the incident field shape -- then
verify by one weak-vs-strong pair before batch runs.  Both routes
carry the identical genus-1 correction, so this is purely a
cost/physics-regime choice, not an accuracy fallback.

**PEEC-coil weak coupling: verified against BEM-A on Takahashi
(2026-07-17).**  The fastest forward route (``--coil-solver peec
--coil-step <coil.step>`` + weak + intree-dense, loop auto) was
cross-checked against the BEM-A coil on the same workpiece: P_wp
18.33 vs 18.41 kW (0.55%), H_t 46.22 vs 46.32 kA/m (0.3%), alpha
5110 A @ -174.1 deg vs 5264 A @ -173.9 deg, delta_L -15.9 vs
-15.5 nH -- with the coil stage collapsing from ~82 s
(impedance-EFIE, 4136 tris) to ~1.1 s (filament extraction 0.4 s +
proximity bundle solve 0.7 s); whole weak run ~23 s vs ~109 s.  The
PEEC coil there was a RECONSTRUCTED gapped rect-torus STEP (342 deg
arc, true 18.1 x 16.3 mm section, leads omitted -- Kubota's original
CAD was unavailable), so the 0.55% also bounds the lead contribution
to P_wp.  PEEC + nonlinear ESIM composes too (Karl converged in 15
iterations, P_wp 18.58 kW, loop auto-skipped with the recorded ESIM
skip reason).  L_coil differs from BEM-A by the leads (66 vs 102 nH):
use BEM-A when coil L / R themselves are deliverables.

**psi-Poisson incident -- the P1 weak route, basis-determined
(2026-07-17).**  Replaces the axis-ray + horizontal-ray path
integration of phi_inc with a surface-Poisson (Laplace-Beltrami)
projection of the EXACT vertex H_inc
(``radia.bem_sibc_solver.compute_phi_inc_surface_poisson``): mean-zero
psi with ``int grad_S psi . grad_S v = -int H_t,inc . grad_S v``, so
the branch-cut wall a spanning-path integral drags across the surface
disappears and the reconstruction is L2-optimal.  It also replaces
per-vertex quadrature rays with one batched field evaluation.  The
fail-loud gate raises when
``||grad_S psi + H_t,inc|| / ||H_t,inc|| > 10%``; no silent fallback is
used when the projected scalar potential does not exist.
Works with both coil sources (BEM-A panels / PEEC filaments).
Goldens: ``tests/test_phi_inc_poisson.py`` (icosphere: uniform field
recovers psi = -H0 z to <2%, rotational Killing field fires the gate,
winding-invariance, complex linearity).

**THE FIX IS THE DEFAULT; THE BUGGY ROUTES ARE DELETED (2026-07-17,
Sugahara: a known-buggy path must not stay selectable).**

- phi_inc is NOT a knob: the weak P1 route is surface-Poisson ALWAYS.
  The legacy P1 path-integration option and the transient
  ``--wp-phi-inc`` flag were REMOVED the same day they were added;
  path integration survives ONLY where poisson has no implementation
  yet (Lagrange-P2 edge-node DOFs, and inside the strong driver) as
  the sole route there, with the wall limitation documented.  The
  dead ``compute_phi_inc_from_filaments_surface_path`` (reverted
  2026-05-21) was deleted outright.
- ``--wp-loop-dof {auto,on}`` (default auto; bare flag = on): auto
  applies the loop DOF on a genus-1 workpiece when the prerequisites
  hold (weak + linear SIBC + ``--wp-bem-backend intree-dense`` + P1)
  and otherwise SKIPS with a recorded ``wp_loop_dof_skip_reason``
  (genus >= 1 skips keep the ``P_wp_caveat``).  "on" fails loud on
  unmet prerequisites.  There is deliberately NO "off": the
  un-extended genus-1 solve is a known +25-30% over-estimate.

With ``--wp-bem-backend intree-dense`` the defaults enable the complete
genus-1 correction.  With the HACApK backend the projected incident
potential still applies, while the loop-DoF skip reason and caveat are
emitted so the limitation is explicit. The Simulink block and comparison
notebook share one ``auto``-defaulted loop-DoF setting and no
incident-potential knob.

Detection is built in: ``bem_sibc_solver.surface_euler_characteristic``
+ ``calc_inductance._wp_genus_check`` -- every weak/strong run logs a
WARNING for genus >= 1 and emits ``wp_euler_chi`` / ``wp_genus`` +
``P_wp_caveat`` in the JSON (genus-conditional; genus-0 gets NO caveat,
backed by the sphere benchmark).

For a genus-1 workpiece, use ``--wp-bem-backend intree-dense`` when the
mesh fits so the dense SL/DL loop extension can be assembled.  The
HACApK path currently reports ``wp_loop_dof_skip_reason`` and
``P_wp_caveat`` instead of silently claiming the missing cohomology
mode.  Genus-0 workpieces carry no topology caveat.

## When to use this

- Workpiece moves/rotates relative to coil (different relative position per step)
- Coil geometry is fixed
- Linear or nonlinear SIBC on workpiece (delta << WP size)
- Frequency range where filament model captures coil AC redistribution

## Advantages

- **No coil volume mesh needed** — coil is a filament bundle from
  `CoilBuilder.to_filaments()` or `coil_from_build123d_sweep()`.
- Only workpiece surface mesh + per-step rebuild as WP moves.
- BEM solver construction cached once per WP mesh (`IHWorkpieceContext`).
- L from PEEC `Z_port` is mesh-independent (matches analytical to <1%).
- Per-coil-variant eval ~3 seconds (MNA + Biot-Savart + BEM linear solve).

## Pipeline (validated 2026-04-18)

```python
from coil_builder import CoilBuilder
from peec_bundle import build_bundle_solver
from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                              compute_phi_inc_from_filaments)

# 1. Coil geometry via CAD-first API (NO mesh)
coil = (CoilBuilder(current=I_total)
        .set_start([R_major, 0, 0])
        .set_cross_section(width=SIDE, height=SIDE)
        .add_arc(radius=R_major, arc_angle=355))
paths, _ = coil.to_filaments(nw=3, nh=3, frequency=f,
                              sigma=5.8e7, mu_r=1.0, n_arc=40)
info = coil._last_filament_info
dw, dh = info['w']/3, info['h']/3

# 2. PEEC MNA -> per-filament complex AC currents
solver_peec, seg_of_fil, _, _ = build_bundle_solver(paths, dw, dh, 5.8e7)
I_branch = solver_peec.compute_branch_currents(f, [I_total+0j])
I_peec = [np.mean([I_branch[s] for s in segs]) for segs in seg_of_fil]

# 3. Inductance from port impedance (MESH-INDEPENDENT)
Z_port = solver_peec.compute_port_impedance(f)
L = Z_port.imag / (2*pi*f)
R_ac = Z_port.real

# 4. phi_inc at WP surface via Biot-Savart line integrals
wp_pts = np.array([wp_mesh.vertices[v].point for v in range(wp_mesh.nv)])
phi_inc = compute_phi_inc_from_filaments(wp_pts, paths, I_peec, n_quad=20)

# 5. BEM-SIBC solve on workpiece (CACHE THE SOLVER)
bem = ScalarBIESIBCSolver(wp_mesh)  # one-time ~10s
sol = bem.solve(phi_inc, Z_s, omega)
H_t_rms = sol['H_t_rms']
P_total = sol['P_density'] * sol['area']
```

## Validation vs A-V FEM truth (2026-04-18)

Geometry: gapped torus R=30mm, a=3mm, 5 deg gap. Steel WP R=25mm H=25mm.
Frequency 7 kHz (Cu delta = 0.79 mm, delta/a = 0.26 — weak skin effect).

| Method | L [nH] | H_t [A/m] | P [W] |
|--------|--------|-----------|-------|
| A-V FEM (coil-sigma) — truth | 100.25 | 10.36 | 4.95e-4 |
| PEEC + BEM-SIBC | 88.57 | 10.67 | 5.24e-4 |
| Agreement | -11.6% | **+3.0%** | +4.4% |

L difference reflects FEM truncation vs analytical; H_t/P agree within method
differences (FEM volume vs BEM surface on different WP mesh).

## mu_r support

Linear SIBC in the workpiece: `Z_s = (1+j)/(sigma*delta)` with
`delta = sqrt(2/(omega*mu_0*mu_r*sigma))`. Pass mu_r via the material setup
in `IHWorkpieceContext(mesh, freq, sigma, mu_r=100, ...)`.

Nonlinear BH (saturation): use ESIM cell problem to compute per-panel Z_s,
then `ScalarBIESIBCSolver.solve(phi_inc, Z_s_array, omega)` with per-node
Z_s array. Karl iteration alternates ESIM and BEM solves.

## Coil side: copper (linear, mu_r=1) is the default assumption

PEEC filament MNA uses sigma_coil for AC R and M. Copper at typical IH
frequencies is non-magnetic (mu_r=1). For ferromagnetic coil materials,
PEEC filament model would need revision (not currently supported).

## Panel integration: calc_inductance.py (v4.25.0+ unified Layer 4)

The `src/radia/panels/calc_inductance.py` script is the Layer-4
subprocess for the IH panel's "PEEC + BEM weak coupling" method.  It
unified three legacy scripts (`calc_peec_inductance.py`,
`calc_peec_bem.py`, `calc_coil_bem_a_workpiece.py`) into one CLI with a
``--coil-solver {peec,bem-a}`` switch.  Invocation:

    python calc_inductance.py \\
        --coil-solver peec --coil-step <coil.step> \\
        --vol <wp.vol> --wp-label sibc \\
        --sigma <S/m> --mu-r <...> --half-thickness <m> \\
        --frequency <Hz> --current <A>

It wraps the demo pipeline above into a .vol-based / JSON-stdout
interface.

### BND extraction from wp-as-hole geometry (CRITICAL)

When the workpiece is modeled as a HOLE in the air volume with a sibc
sideset (closed-torus / skin sample convention), multiple FaceDescriptors
may share the same 'sibc' name (e.g. air_top + air_bot halves after
webcut).  Each FD's BND element winding gives an outward normal relative
to its own domin, which CAN BE OPPOSITE between FDs.  A blanket flip is
WRONG — 77% of triangles needed flipping in the closed-torus sample,
23% did not.

Robust fix (implemented in `_extract_bnd_only_inline`, formerly in the
deleted `calc_peec_bem.py`):
1. Collect the wp-surface vertex cloud.
2. Compute wp_centroid = mean(coords).
3. For each triangle, check sign of dot(normal, centroid_to_face).
4. Flip winding if inward.

Without this, BEM P output is silently wrong by 2-20x (not a failed
solve — a plausible-looking but incorrect number).

### Golden regression test (policy)

`tests/panels/test_peec_bem_golden.py` runs `calc_inductance.py` on
`ih_peec_bem_coarse.vol` at 7 kHz Cu and asserts:
- L_coil within 1% of captured 88.57 nH (geometry-only, tight)
- P_wp within 15% of 2D SIBC reference 6.63e-5 W (mesh-sensitive)

Expected values in `tests/panels/golden/peec_bem_coarse_7kHz_Cu.json`.
Re-generate the sample via Cubit headless:
`coreform_cubit -batch -nojournal -nographics -input
src/radia/panels/samples/ih_peec_bem_coarse.jou`.

### Sibling: calc_fem_coilmesh.py (2026-04-19 FINAL, A-V + wp SIBC)

Gapped-torus A-V formulation following the MCP-documented
`INDUCTION_HEATING_AV_COIL_SIGMA` pattern.  Coil is volumetrically
meshed with sigma; phi (scalar potential) is defined on coil material
only, Dirichlet 1 on 'source' / 0 on 'sink' port sidesets.  Volume-
integral current extraction `I_out = int J . grad(psi_n) dV`.  WP
SIBC Robin as usual, Kelvin exterior.

**Design journey that got us here (retired options)**:
- Neumann-K on coil surface alone (2026-04-19 am): L collapsed 10x
  because Robin SIBC cancels K_imposed exactly.  **Retired.**
- Neumann-K without Robin (pure Dirichlet-source): got P_wp within
  2D SIBC ref, but P_wp matched PEEC+BEM only to 24% — the uniform-K
  surface source misses coil proximity effects.  **Retired.**
- Dowell thin-wire analytical P_coil: good for research, but not a
  panel-level primitive (hides mesh-dependence from the user).
  **Demoted to examples/ on 2026-04-19.**

### Validation (2026-04-19, ih_fem_kelvin_skin_fine, Cu 7 kHz)

| Method | P_wp [W] | L [nH] | P_coil [W] |
|--------|----------|--------|------------|
| 2D axisym (ref) | 6.63e-5 (SIBC) | 50.54 | ~8.4e-5 (implied) |
| PEEC+BEM | 6.48e-5 (-2.2%) | — | — |
| FEM A-V   | **6.54e-5** (-1.3%) | 47.32 (-6.4%) | 1.48e-4 (+76% mesh-sens.) |

**P_wp agreement PEEC+BEM vs FEM A-V: 1%** (!) — the primary IH
engineering metric.  L is within 10% (coil mesh-sensitive).  P_coil
is volumetric `|J|^2/sigma` and mesh-sensitive; captured value with
tolerance 15% in golden.  For exact P_coil use Dowell analytical
in examples/.

### .vol / .jou requirements

- 'coil' material, VOLUMETRIC-MESHED with `volume size <= delta/3`
  (ideally; skin_fine sample uses 0.16mm surface + 0.5mm volume for
  Cu 7kHz delta=0.79mm).
- 'source', 'sink' sidesets on the two gap-face circles.  Real IH
  coils have physical port terminations; use those.  Closed-torus
  topology cannot drive A-V.
- 'sibc' sideset on wp hole boundary.
- 'kelvin' material (+ Periodic identification) for open boundary.
- Canonical .jou: `src/radia/panels/samples/ih_fem_kelvin_skin_fine.jou`.

### Application interface

`IHDesignSpec`, the Induction Heating Simulink block, and the temporary IH
notebook comparison present the above through
settings sections (Method / Drive / Coil material / Coil geometry /
Workpiece material / Workpiece impedance / Linear solver / Advanced).

Non-obvious features:

1. **Material presets** (Copper / Aluminum / Brass / Steel mu_r=100 /
   Steel mu_r=500 / Stainless 304 / Custom).  Selecting a preset
   auto-fills sigma+mu_r and disables their lineedits; `Custom`
   re-enables manual entry.  **Coil and wp materials are INDEPENDENT**
   — setting one does not change the other.

2. **Impedance model** (Linear SIBC / Nonlinear ESIM-WIP).  Selecting
   ESIM reveals bh_file, max_iter, tol widgets.  Calc scripts accept
   `--impedance-model esim` but return `{"error":"ESIM WIP"}` — the
   CLI path is plumbed for future extension.

3. **Linear solver** per-method:
   - PEEC+BEM: Dense LU (small) / HACApK (large, O(N log N))
   - FEM A-V: pardiso / shifted AMS (p=1) / BDDC (p>=2) / iccg
   `shifted_ams` and `hacapk` are CLI-plumbed but WIP (calc returns
   error).  `pardiso` + `dense` are the tested paths.

4. **.vol label validation** on load via `inspect_vol_labels`: status
   label shows 'OK' (green), 'warn' (amber, e.g. missing kelvin),
   'ERROR' (red, e.g. missing source/sink for FEM).  Run button
   disabled on errors.

5. **Physics sanity** shown under method:
   - wp delta and R/delta ratio; warns when R/delta < 3
     (SIBC approximation breaks, volumetric FEM preferred)
   - Post-solve: coil h_max vs delta (flags under-resolved coil mesh)

6. **Output** formatted in IH Summary section:
   - `L_coil (vacuum)` vs `L_total (with wp)` distinguished
   - `P_workpiece` (primary) / `P_coil (mesh-sensitive ±15%)` /
     `P_total` / `Heating efficiency = P_wp/P_total`
   - `H_t_rms`, `wp_area`, `coil delta`, `coil h_max` with OK/WARN tags

### Implementation notes (2026-04-19)

1. phi H1 DEFINED-ON coil only (`definedon=mesh.Materials('coil')`).
2. `A + grad(phi)` compound term in the eddy bilinear form on coil.
3. Dirichlet lift: set gf_phi = 1 on source, then solve
   `A x^(-1) (- A . gfu) + gfu = gfu_new`.
4. Post-solve volume current extraction uses a SEPARATE H1 scalar
   `psi_n` with dirichlet=sink, Set(1) on source, then
   `I_out = Integrate(J_coil * grad(psi_n) * dx('coil'))`.  Scale
   the entire gfu by I_target / I_out to normalize.
5. WP H_t from SIBC Robin:
   `H_t_rms = |jw/Z_s_wp| * sqrt(int|A_t|^2 dS / A_wp)`.
6. P_coil volumetric: `0.5 / sigma * int |J|^2 dV` over coil material.
   Warn if coil h_max > delta — mesh-sensitive accuracy.
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
Driving force: phi Dirichlet lift (phi=1 on source, phi=0 on sink).
Total current I_out computed post-solve, solution scaled to I_total.

```python
# FES: HCurl(A) x H1(phi, conductor-only)
fesA = HCurl(mesh, order=1, nograds=True, complex=True, dirichlet=dirichlet_bnd)
fesA_p = Periodic(fesA) if has_kelvin else fesA
fesPhi = H1(mesh, order=1, complex=True, definedon="coil", dirichlet="source|sink")
fes = fesA_p * fesPhi
(A, phi), (N, psi) = fes.TnT()

# Bilinear form (KEY: unified s*sigma*(A+grad(phi))*(N+grad(psi)) term)
a = BilinearForm(fes, symmetric=True)
a += nu_cf*curl(A)*curl(N)*dx(bonus_intorder=4)
a += reg*NU_0*A*N*dx("<non-kelvin materials>")
a += robin*A.Trace()*N.Trace()*ds("sibc")
a += s*sigma_coil*(A+grad(phi))*(N+grad(psi))*dx("coil")
a.Assemble()

# Dirichlet lift: phi=1 on source, phi=0 on sink
gfu = GridFunction(fes)
gf_A, gf_phi = gfu.components
gf_phi.Set(CF(1), definedon=mesh.Boundaries("source"))
r = gfu.vec.CreateVector()
r.data = -a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * r

# Post-normalization to I_total
s_eddy = 1j * omega
J_coil = -s_eddy * sigma_coil * (gf_A + grad(gf_phi))
fes_psi_n = H1(mesh, order=1, complex=True, dirichlet="sink", definedon="coil")
gf_psi_n = GridFunction(fes_psi_n)
gf_psi_n.Set(CF(1), definedon=mesh.Boundaries("source"))
I_out = Integrate(J_coil * grad(gf_psi_n) * dx("coil"), mesh)
scale = I_total / I_out
gfu.vec.data = complex(scale) * gfu.vec
```

## Gotchas

1. **nograds=True on HCurl** — standard NGSolve policy; phi provides the
   gauge for the gradient part of A in the conductor. Both together give
   a well-posed A-V system.
2. **phi defined only on conductor** (`definedon="coil"`). If defined globally,
   the unconstrained phi in air breaks the solution (A becomes pure gradient,
   curl(A)=0 trivial solution).
3. **Periodic Kelvin + pardiso works**. BDDC+CG diverges without additional
   regularization (needs careful preconditioner design for saddle-point system).
4. **GND nodeset does NOT constrain HCurl** — Dirichlet on HCurl works only
   on sidesets (boundary faces), not nodesets (vertices). Regularization term
   `reg * NU_0 * A * N * dx` provides uniqueness instead.
5. **Mesh requirement**: coil surface elements < delta/3 (typically 0.3 mm
   for Cu at 7 kHz). Use `surface with length < X interval N` in Cubit.
6. **Post-normalize**: solve with V_applied=1 (phi Dirichlet), compute I_out
   as integral of J.grad(psi_n) over coil, scale gfu by I_total/I_out.

## Reference: existing A-V implementation
Internal A-V implementation notebook by K. Sugahara (machine-local path
omitted). Production code pattern for A-V + BDDC+CG (non-Kelvin outer).
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
        "peec_bem_sibc": INDUCTION_HEATING_PEEC_BEM_SIBC,
        "av_coil_sigma": INDUCTION_HEATING_AV_COIL_SIGMA,
        "failed_approaches": INDUCTION_HEATING_FAILED_APPROACHES,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
