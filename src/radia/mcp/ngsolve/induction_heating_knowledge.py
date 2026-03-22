"""
Induction heating simulation knowledge base for the Radia MCP server.

Covers: complete EM-to-thermal workflow for induction heating using NGSolve,
including Gmsh mesh loading, A-Phi eddy current formulation, Joule heat
computation, transient thermal analysis, rotating workpiece, VTK output,
and post-processing patterns.

Sources:
  - W:\\31_Go-Tech\\10_誘導加熱の解析\\2025_09_20_toymodel_Gmsh\\NGSolve
    (K. Sugahara's production induction heating simulations)
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

## ReadGmsh: Loading External .msh Files

NGSolve can load Gmsh v2 format meshes via `netgen.read_gmsh.ReadGmsh`.
This is the standard workflow when using Cubit or Gmsh for mesh generation.

```python
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

# Load Gmsh v2 mesh (1st or 2nd order)
ngmesh = ReadGmsh("toymodel.msh")
mesh = Mesh(ngmesh)

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

### Gmsh Mesh from Cubit (Cubit2gmsh Workflow)

```python
# Cubit2gmsh.py: Generate Gmsh mesh from Cubit journal file
# 1. Create geometry in Cubit (.jou file)
# 2. Run Cubit2gmsh.py to export to Gmsh v2 format
# 3. Load in NGSolve with ReadGmsh

# Typical Cubit export command (inside Cubit2gmsh.py):
# cubit.cmd(f'export mesh "{output_msh}" overwrite gmsh_ver2')
```

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
from netgen.read_gmsh import ReadGmsh
import numpy as np

# ============================================================
# 1. Load Gmsh mesh
# ============================================================
ngmesh = ReadGmsh("toymodel.msh")
mesh = Mesh(ngmesh)

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

## Heat Equation with Joule Heat Source

Solves the unsteady heat equation in the workpiece:

    rho * c * dT/dt = div(kappa * grad(T)) + Q

with convection boundary condition on work surface:

    -kappa * dT/dn = h * (T - T_ext)

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
# Rotating Workpiece Simulation

## Mesh Deformation for Workpiece Rotation

For induction heating with a rotating workpiece, the Joule heat source Q
rotates with the workpiece. This is implemented by rotating the mesh
at each thermal time step.

```python
from ngsolve import *
import numpy as np

# ============================================================
# Rotation parameters
# ============================================================
rpm = 12            # rotations per minute
omega = 2 * np.pi * rpm / 60  # angular velocity [rad/s]

# ============================================================
# Rotation matrix as CoefficientFunction
# ============================================================
def get_rotation_matrix(angle):
    \"\"\"3D rotation matrix around Z-axis.\"\"\"
    return CF((
        cos(angle), -sin(angle), 0,
        sin(angle),  cos(angle), 0,
        0,           0,          1
    )).Reshape((3, 3))

# ============================================================
# Mesh deformation field
# ============================================================
# Create deformation = rotated_position - original_position
fes_def = VectorH1(mesh, order=1)
gf_deform = GridFunction(fes_def)

def apply_rotation(mesh, gf_deform, angle):
    \"\"\"Apply rotation to work regions via mesh deformation.\"\"\"
    rotmat = get_rotation_matrix(angle)

    # Deformation = R*x - x  (where x is original position)
    pos = CF((x, y, z))
    rotated = rotmat * pos
    deformation = rotated - pos

    # Only deform work regions (coil stays fixed)
    # Use definedon to restrict deformation
    gf_deform.Set(deformation)
    mesh.SetDeformation(gf_deform)

# ============================================================
# Time stepping with rotation
# ============================================================
dt = 0.5       # time step [s]
t_end = 5.0    # total time (= 1 full rotation at 12 rpm)

t = 0.0
step = 0

while t < t_end:
    # 1. Compute rotation angle
    angle = omega * t

    # 2. Apply mesh deformation (rotate workpiece)
    apply_rotation(mesh, gf_deform, angle)

    # 3. Re-project Q onto rotated mesh
    #    (Q was computed in the original orientation)
    #    The mesh deformation handles the coordinate transformation.

    # 4. Solve thermal time step (same as non-rotating case)
    rhs = gfT.vec.CreateVector()
    rhs.data = m_T.mat * gfT.vec + dt * f_T.vec
    gfT.vec.data = mstar_inv * rhs

    # 5. Output VTK at this time step
    vtk = VTKOutput(mesh,
                    coefs=[gfT, Q_gf],
                    names=["Temperature", "Q"],
                    filename=f"rotating_step_{step:03d}",
                    subdivision=1, legacy=False)
    vtk.Do()

    # 6. Remove deformation before next step
    mesh.UnsetDeformation()

    t += dt
    step += 1
    print(f"Step {step}: t={t:.1f}s, angle={np.degrees(angle):.1f} deg")
```

## Key Considerations for Rotation

1. **Mesh deformation vs re-meshing**: `SetDeformation` is much faster
   than re-meshing at each step. It deforms the existing mesh without
   changing topology.

2. **Heat source rotation**: The Joule heat Q is computed once (EM solve)
   and then the mesh is rotated so Q effectively rotates with the workpiece.

3. **Coil stays fixed**: Only the workpiece rotates. The heat source pattern
   (fixed in the lab frame) sweeps over the rotating workpiece surface.

4. **VTK output per step**: Each time step produces a separate .vtu file.
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
