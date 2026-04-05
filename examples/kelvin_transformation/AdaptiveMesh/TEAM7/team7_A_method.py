"""
TEAM Problem 7: 3D Eddy Current Analysis using A-method

Problem description:
- Asymmetric conductor with hole
- Excited by a racetrack coil
- Calculate eddy currents and magnetic field

Method:
- A-phi formulation (vector potential)
- Coil current given as CoefficientFunction (source current density J0)
- curl(1/mu * curl(A)) + sigma * dA/dt = J0  (in conductor)
- curl(1/mu * curl(A)) = J0  (in air/coil)

For static case (or time-harmonic):
- curl(1/mu * curl(A)) + j*omega*sigma*A = J0
"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

import math
from numpy import pi, sqrt
import numpy as np
from ngsolve import *
from netgen.occ import *

print("=" * 70)
print("TEAM Problem 7: 3D A-method with Direct Coil Current")
print("=" * 70)

# ============================================================
# Physical parameters
# ============================================================
mu0 = 4 * pi * 1e-7
sigma_al = 3.526e7  # Aluminum conductivity [S/m]
freq = 50  # Hz (for time-harmonic analysis)
omega = 2 * pi * freq
current = 2742  # Coil current [A]

# ============================================================
# Geometry: TEAM Problem 7
# ============================================================
print("\n1. Creating geometry...")

# Conductor plate with hole
pc0 = gp_Pnt(0, 0, 0)
pc1 = gp_Pnt(294e-3, 294e-3, 19e-3)
conductor = Box(pc0, pc1)

# Hole in conductor
ph0 = gp_Pnt(18e-3, 18e-3, 0)
ph1 = gp_Pnt(126e-3, 126e-3, 19e-3)
hole = Box(ph0, ph1)

conductor = conductor - hole
conductor.mat("conductor")
conductor.faces.name = "conductor_bnd"

# Coil parameters (racetrack coil)
coil_a = 12.5e-3  # half width of coil cross-section (radial)
coil_b = 50e-3    # half height of coil cross-section (axial)
coil_rc = 25e-3 + coil_a  # radius of corner arcs
coil_w1 = 50e-3   # half-length of straight sections
coil_w2 = coil_w1 + coil_rc
coil_w3 = coil_w2 + coil_a

# Coil center position
coil_x0 = 294e-3 - coil_w3
coil_y0 = coil_w3
coil_z0 = 19e-3 + 30e-3 + coil_b  # Above conductor

# Create coil geometry (simplified as a torus-like shape)
# For simplicity, we'll define the coil region and assign current density

# Coil cross-section dimensions
coil_inner_r = coil_rc - coil_a
coil_outer_r = coil_rc + coil_a
coil_height = 2 * coil_b

# Create a simplified rectangular coil region
coil_box = Box(
    gp_Pnt(coil_x0 - coil_w2 - coil_a, coil_y0 - coil_w2 - coil_a, coil_z0 - coil_b),
    gp_Pnt(coil_x0 + coil_w2 + coil_a, coil_y0 + coil_w2 + coil_a, coil_z0 + coil_b)
)
coil_inner = Box(
    gp_Pnt(coil_x0 - coil_w2 + coil_a, coil_y0 - coil_w2 + coil_a, coil_z0 - coil_b - 0.001),
    gp_Pnt(coil_x0 + coil_w2 - coil_a, coil_y0 + coil_w2 - coil_a, coil_z0 + coil_b + 0.001)
)
coil = coil_box - coil_inner
coil.mat("coil")
coil.faces.name = "coil_bnd"

# Air region (bounding box)
air_margin = 0.3  # 30cm margin
air_box = Box(
    gp_Pnt(-air_margin, -air_margin, -air_margin),
    gp_Pnt(294e-3 + air_margin, 294e-3 + air_margin, coil_z0 + coil_b + air_margin)
)

# Hole region (inside conductor hole, filled with air)
hole_air = Box(ph0, ph1)
hole_air.mat("hole")

# Air = bounding box - conductor - coil - hole
air = air_box - conductor - coil - hole_air
air.mat("air")

# Mark outer boundary
for face in air_box.faces:
    face.name = "outer"

# Combine all regions
geo = Glue([conductor, coil, hole_air, air])

# Generate mesh
print("  Generating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.05))
print(f"  Mesh elements: {mesh.ne}")
print(f"  Materials: {mesh.GetMaterials()}")

# ============================================================
# Define coil current density as CoefficientFunction
# ============================================================
print("\n2. Defining coil current density...")

# Coil cross-section area
coil_area = 4 * coil_a * coil_b  # Approximate cross-section area
J0_magnitude = current / coil_area  # Current density magnitude

print(f"  Coil current: {current} A")
print(f"  Coil cross-section area: {coil_area*1e6:.2f} mm^2")
print(f"  Current density magnitude: {J0_magnitude:.2e} A/m^2")

# Current density direction (follows coil path)
# For a rectangular coil centered at (coil_x0, coil_y0), the current flows:
# - In +y direction on the right side (x > coil_x0 + threshold)
# - In -y direction on the left side (x < coil_x0 - threshold)
# - In -x direction on the top (y > coil_y0 + threshold)
# - In +x direction on the bottom (y < coil_y0 - threshold)

# Simplified: define current density based on position relative to coil center
# We'll use a smooth function to determine direction

def make_coil_current_cf():
    """Create CoefficientFunction for coil current density"""
    # Relative position from coil center
    dx = x - coil_x0
    dy = y - coil_y0

    # Determine which segment of the coil we're in
    # Use smooth transition with IfPos

    # Right side: +y direction
    # Left side: -y direction
    # Top: -x direction
    # Bottom: +x direction

    # Simple approach: current flows tangentially around the coil center
    # J = J0 * (-dy, dx, 0) / r  (normalized tangential direction)

    r = sqrt(dx**2 + dy**2 + 1e-10)  # Avoid division by zero

    # Tangential direction (counterclockwise when viewed from +z)
    Jx = -J0_magnitude * dy / r
    Jy = J0_magnitude * dx / r
    Jz = CF(0)

    J_coil = CF((Jx, Jy, Jz))

    return J_coil

J0_cf = make_coil_current_cf()

# Apply only in coil region
J0 = CoefficientFunction([
    CF((0, 0, 0)),  # conductor
    J0_cf,          # coil
    CF((0, 0, 0)),  # hole
    CF((0, 0, 0))   # air
])

# ============================================================
# A-method formulation (static case first)
# ============================================================
print("\n3. Setting up A-method formulation...")

order = 2

# H(curl) space for vector potential A
fes = HCurl(mesh, order=order, dirichlet="outer")
A, v = fes.TnT()

print(f"  DOFs: {fes.ndof}")

# Material properties
mu_dict = {"conductor": mu0, "coil": mu0, "hole": mu0, "air": mu0}
sigma_dict = {"conductor": sigma_al, "coil": 0, "hole": 0, "air": 0}

mu_cf = CoefficientFunction([mu_dict[mat] for mat in mesh.GetMaterials()])
sigma_cf = CoefficientFunction([sigma_dict[mat] for mat in mesh.GetMaterials()])
inv_mu_cf = CoefficientFunction([1/mu_dict[mat] for mat in mesh.GetMaterials()])

# ============================================================
# Static magnetostatic problem (no eddy currents)
# curl(1/mu * curl(A)) = J0
# ============================================================
print("\n4. Solving static magnetostatic problem...")

a_static = BilinearForm(fes)
a_static += inv_mu_cf * InnerProduct(curl(A), curl(v)) * dx
# Gauge condition (small regularization)
a_static += 1e-10 * InnerProduct(A, v) * dx
a_static.Assemble()

f_static = LinearForm(fes)
f_static += InnerProduct(J0, v) * dx
f_static.Assemble()

gfA_static = GridFunction(fes)
gfA_static.vec.data = a_static.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f_static.vec

print("  Static solution computed.")

# Compute B = curl(A)
B_static = curl(gfA_static)

# ============================================================
# Time-harmonic problem (with eddy currents)
# curl(1/mu * curl(A)) + j*omega*sigma*A = J0
# ============================================================
print("\n5. Solving time-harmonic problem (f = 50 Hz)...")

# Complex-valued space
fes_c = HCurl(mesh, order=order, dirichlet="outer", complex=True)
A_c, v_c = fes_c.TnT()

a_harmonic = BilinearForm(fes_c)
a_harmonic += inv_mu_cf * InnerProduct(curl(A_c), curl(v_c)) * dx
a_harmonic += 1j * omega * sigma_cf * InnerProduct(A_c, v_c) * dx
# Gauge condition
a_harmonic += 1e-10 * InnerProduct(A_c, v_c) * dx
a_harmonic.Assemble()

f_harmonic = LinearForm(fes_c)
f_harmonic += InnerProduct(J0, v_c) * dx
f_harmonic.Assemble()

gfA_harmonic = GridFunction(fes_c)
gfA_harmonic.vec.data = a_harmonic.mat.Inverse(fes_c.FreeDofs(), inverse="sparsecholesky") * f_harmonic.vec

print("  Time-harmonic solution computed.")

# Compute B and J (eddy current)
B_harmonic = curl(gfA_harmonic)
J_eddy = -1j * omega * sigma_cf * gfA_harmonic  # Eddy current density

# ============================================================
# Post-processing
# ============================================================
print("\n6. Post-processing...")

# Magnetic energy (static)
W_static = 0.5 * Integrate(inv_mu_cf * InnerProduct(B_static, B_static), mesh)
print(f"  Static magnetic energy: {W_static:.6e} J")

# Magnetic energy (harmonic, time-averaged)
W_harmonic = 0.25 * Integrate(inv_mu_cf * InnerProduct(B_harmonic, Conj(B_harmonic)), mesh).real
print(f"  Time-harmonic magnetic energy (avg): {W_harmonic:.6e} J")

# Joule loss (harmonic)
P_joule = 0.5 * Integrate(InnerProduct(J_eddy, Conj(J_eddy)) / (sigma_cf + 1e-10),
                          mesh, definedon=mesh.Materials("conductor")).real
print(f"  Joule loss in conductor: {P_joule:.6e} W")

# Sample B field at specific points (TEAM 7 measurement points)
# Points along y = 72mm, z = 34mm
print("\n  B field at measurement line (y=72mm, z=34mm):")
print(f"  {'x [mm]':>10s} {'Bz [mT]':>12s}")
print(f"  {'-'*25}")

for x_mm in [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180]:
    pnt = mesh(x_mm * 1e-3, 72e-3, 34e-3)
    try:
        Bz = B_static(pnt)[2]
        print(f"  {x_mm:10.0f} {Bz*1000:12.4f}")
    except:
        print(f"  {x_mm:10.0f} {'N/A':>12s}")

# ============================================================
# Visualization
# ============================================================
print("\n7. Visualization setup...")

try:
    from ngsolve.webgui import Draw
    # Draw(B_static, mesh, order=2, draw_surf=False, name="B_static")
    # Draw(Norm(B_harmonic), mesh, order=2, draw_surf=False, name="|B_harmonic|")
    print("  Use Draw(B_static, mesh) to visualize")
except:
    print("  WebGUI not available")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("Summary: TEAM Problem 7 - A-method")
print("=" * 70)
print(f"""
Problem setup:
  - Conductor: Aluminum plate (294x294x19 mm) with hole
  - Coil: Racetrack coil, I = {current} A
  - Frequency: {freq} Hz

Results:
  - Static magnetic energy: {W_static:.6e} J
  - Time-harmonic energy:   {W_harmonic:.6e} J
  - Joule loss:             {P_joule:.6e} W

Note:
  This is a simplified implementation. For accurate TEAM 7 results:
  - Coil geometry should match exact racetrack shape
  - Mesh refinement near conductor surface needed
  - Compare with measurement data at specified points
""")
