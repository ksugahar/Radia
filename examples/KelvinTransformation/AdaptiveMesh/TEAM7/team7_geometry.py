"""
TEAM Problem 7: Asymmetrical Conductor with a Hole - Geometry Definition

Reference: https://modules.freefem.org/modules/team7/
         https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/

Problem Description:
- Thick aluminum plate with an eccentric square hole
- Multi-turn racetrack coil placed asymmetrically above the plate
- Calculate eddy currents and magnetic field at 50 Hz and 200 Hz

Geometry (all dimensions in mm):
- Aluminum plate: 294 x 294 x 19 mm
- Hole: 108 x 108 mm at position (18, 18) to (126, 126)
- Coil: Racetrack shape above the plate
- Coil current: 2742 A-turns (2742 turns x 1 A)

"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ksugahar\Lib\site-packages')

from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

# ============================================================
# Physical Parameters
# ============================================================
mu0 = 4 * pi * 1e-7  # Vacuum permeability [H/m]
sigma_al = 3.526e7   # Aluminum conductivity [S/m]

# Coil parameters
n_turns = 2742       # Number of turns
I_coil = 1.0         # Current per turn [A]
NI = n_turns * I_coil  # Total ampere-turns [A]

# Frequencies for analysis
freq_50 = 50.0   # Hz
freq_200 = 200.0 # Hz

# Skin depth at 50 Hz: delta = sqrt(2/(omega*mu*sigma)) ≈ 12 mm
omega_50 = 2 * pi * freq_50
delta_50 = sqrt(2 / (omega_50 * mu0 * sigma_al))
print(f"Skin depth at 50 Hz: {delta_50*1000:.2f} mm")

# ============================================================
# Geometry Parameters (mm -> m)
# ============================================================
# All coordinates are converted from mm to m

# Aluminum plate
plate_x1, plate_y1, plate_z1 = 0, 0, 0
plate_x2, plate_y2, plate_z2 = 294e-3, 294e-3, 19e-3

# Hole in the plate
hole_x1, hole_y1 = 18e-3, 18e-3
hole_x2, hole_y2 = 126e-3, 126e-3

# Coil parameters (racetrack coil)
# Based on TEAM 7 specification:
# - Coil cross-section: 25mm (radial) x 100mm (axial)
# - Inner corner radius: 25mm
# - Straight section length: 50mm each side
# - Coil positioned asymmetrically

# Coil cross-section
coil_a = 12.5e-3   # Half-width of coil cross-section (radial) [m]
coil_b = 50e-3     # Half-height of coil cross-section (axial) [m]

# Coil path parameters
coil_r_inner = 25e-3   # Inner corner radius [m]
coil_r_center = coil_r_inner + coil_a  # Center radius of coil cross-section
coil_r_outer = coil_r_inner + 2*coil_a  # Outer corner radius

# Straight section half-lengths
coil_Lx = 50e-3  # Half-length in x-direction [m]
coil_Ly = 50e-3  # Half-length in y-direction [m]

# Coil center position (in the xy-plane)
# Positioned asymmetrically above the plate
coil_x0 = 294e-3 - (coil_Lx + coil_r_center)  # 294 - 50 - 37.5 = 206.5 mm
coil_y0 = coil_Lx + coil_r_center              # 50 + 37.5 = 87.5 mm

# Coil vertical position
# Gap between plate top and coil bottom: 30 mm
coil_gap = 30e-3
coil_z_bottom = plate_z2 + coil_gap  # 19 + 30 = 49 mm
coil_z_center = coil_z_bottom + coil_b  # 49 + 50 = 99 mm
coil_z_top = coil_z_bottom + 2*coil_b  # 49 + 100 = 149 mm

print(f"\nCoil center: ({coil_x0*1000:.1f}, {coil_y0*1000:.1f}, {coil_z_center*1000:.1f}) mm")
print(f"Coil z-range: {coil_z_bottom*1000:.1f} to {coil_z_top*1000:.1f} mm")

# ============================================================
# Bounding Box (Air Region)
# ============================================================
# Large enough to capture far-field decay
margin = 0.3  # 300 mm margin (reduced for smaller problem)
air_x1, air_y1, air_z1 = -margin, -margin, -margin
air_x2, air_y2, air_z2 = plate_x2 + margin, plate_y2 + margin, coil_z_top + margin


def create_team7_geometry(maxh=0.03):
    """
    Create TEAM Problem 7 geometry using OCC

    Returns:
        mesh: NGSolve mesh with materials 'plate', 'coil', 'air'
    """
    print("\nCreating TEAM 7 geometry...")

    # ============================================================
    # 1. Aluminum Plate with Hole
    # ============================================================
    plate = Box(Pnt(plate_x1, plate_y1, plate_z1),
                Pnt(plate_x2, plate_y2, plate_z2))

    hole = Box(Pnt(hole_x1, hole_y1, plate_z1 - 0.001),
               Pnt(hole_x2, hole_y2, plate_z2 + 0.001))

    plate = plate - hole
    plate.mat("plate")

    # ============================================================
    # 2. Racetrack Coil
    # ============================================================
    # Create coil as a hollow rectangular prism with rounded corners

    # Outer boundary of coil path
    outer_x1 = coil_x0 - coil_Lx - coil_r_outer
    outer_x2 = coil_x0 + coil_Lx + coil_r_outer
    outer_y1 = coil_y0 - coil_Ly - coil_r_outer
    outer_y2 = coil_y0 + coil_Ly + coil_r_outer

    # Inner boundary of coil path
    inner_x1 = coil_x0 - coil_Lx - coil_r_inner
    inner_x2 = coil_x0 + coil_Lx + coil_r_inner
    inner_y1 = coil_y0 - coil_Ly - coil_r_inner
    inner_y2 = coil_y0 + coil_Ly + coil_r_inner

    # Simplified coil as rectangular ring (without rounded corners for now)
    coil_outer = Box(Pnt(outer_x1, outer_y1, coil_z_bottom),
                     Pnt(outer_x2, outer_y2, coil_z_top))
    coil_inner = Box(Pnt(inner_x1, inner_y1, coil_z_bottom - 0.001),
                     Pnt(inner_x2, inner_y2, coil_z_top + 0.001))

    coil = coil_outer - coil_inner
    coil.mat("coil")

    # ============================================================
    # 3. Air Region
    # ============================================================
    air_box = Box(Pnt(air_x1, air_y1, air_z1),
                  Pnt(air_x2, air_y2, air_z2))

    # Mark outer boundary
    for face in air_box.faces:
        face.name = "outer"

    # Hole region (air inside the plate hole)
    hole_air = Box(Pnt(hole_x1, hole_y1, plate_z1),
                   Pnt(hole_x2, hole_y2, plate_z2))
    hole_air.mat("hole_air")

    # Air = bounding box - plate - coil - hole_air
    air = air_box - plate - coil - hole_air
    air.mat("air")

    # ============================================================
    # 4. Combine and Mesh
    # ============================================================
    geo = Glue([plate, coil, hole_air, air])

    print(f"  Generating mesh (maxh={maxh*1000:.0f} mm)...")
    mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh))
    # Important: Curve order should match FES polynomial order for good convergence
    # Default to order=1 here; caller can adjust with mesh.Curve(order)

    print(f"  Mesh elements: {mesh.ne}")
    print(f"  Materials: {mesh.GetMaterials()}")

    return mesh


def get_coil_current_density(mesh):
    """
    Create CoefficientFunction for coil current density

    The coil carries 2742 A-turns. The current flows tangentially
    around the rectangular coil path.

    For a rectangular racetrack coil:
    - Right side (x > coil_x0 + coil_Lx): current in +y direction
    - Left side (x < coil_x0 - coil_Lx): current in -y direction
    - Top (y > coil_y0 + coil_Ly): current in -x direction
    - Bottom (y < coil_y0 - coil_Ly): current in +x direction
    - Corners: blend between adjacent directions

    Returns:
        J0: CoefficientFunction (vector) for source current density
    """
    # Coil cross-section area
    coil_area = (2 * coil_a) * (2 * coil_b)  # 25mm x 100mm

    # Current density magnitude
    J0_mag = NI / coil_area  # [A/m^2]

    print(f"\nCoil current density:")
    print(f"  Ampere-turns: {NI:.0f} A")
    print(f"  Cross-section: {coil_area*1e6:.1f} mm^2")
    print(f"  J0 magnitude: {J0_mag:.2e} A/m^2")

    # Relative position from coil center
    dx = x - coil_x0
    dy = y - coil_y0

    # Determine which part of the coil we're in
    # Use IfPos for piecewise definition

    # For a rectangular coil:
    # Right side (dx > coil_Lx): Jy = +J0_mag
    # Left side (dx < -coil_Lx): Jy = -J0_mag
    # Top (dy > coil_Ly): Jx = -J0_mag
    # Bottom (dy < -coil_Ly): Jx = +J0_mag

    # Simplified approach: use the sign of dx/dy to determine direction
    # Right/Left sides (|dx| > |dy| scaled by aspect ratio)
    # Top/Bottom sides (|dy| > |dx| scaled by aspect ratio)

    # Normalized coordinates
    dx_norm = dx / (coil_Lx + coil_r_center)
    dy_norm = dy / (coil_Ly + coil_r_center)

    # Determine dominant direction
    # If |dx_norm| > |dy_norm|: we're on left/right side -> Jy dominates
    # If |dy_norm| > |dx_norm|: we're on top/bottom side -> Jx dominates

    # Use smooth blending with atan2-like approach
    angle = atan2(dy_norm, dx_norm)

    # Current direction tangent to the coil path
    # For counterclockwise flow:
    # angle ~ 0 (right): Jy > 0
    # angle ~ pi/2 (top): Jx < 0
    # angle ~ pi or -pi (left): Jy < 0
    # angle ~ -pi/2 (bottom): Jx > 0

    Jx = -J0_mag * sin(angle)
    Jy = J0_mag * cos(angle)
    Jz = CF(0)

    J_coil = CF((Jx, Jy, Jz))

    # Apply only in coil region
    mats = mesh.GetMaterials()
    J_list = []
    for mat in mats:
        if mat == "coil":
            J_list.append(J_coil)
        else:
            J_list.append(CF((0, 0, 0)))

    J0 = CoefficientFunction(J_list)

    return J0


def get_material_properties(mesh):
    """
    Create CoefficientFunctions for material properties

    Returns:
        mu: Permeability CoefficientFunction
        sigma: Conductivity CoefficientFunction
        inv_mu: 1/mu CoefficientFunction
    """
    mats = mesh.GetMaterials()

    # All regions have mu = mu0 (no magnetic materials)
    mu_dict = {mat: mu0 for mat in mats}

    # Only plate has conductivity
    sigma_dict = {mat: 0 for mat in mats}
    sigma_dict["plate"] = sigma_al

    mu_list = [mu_dict[mat] for mat in mats]
    sigma_list = [sigma_dict[mat] for mat in mats]
    inv_mu_list = [1/mu_dict[mat] for mat in mats]

    mu = CoefficientFunction(mu_list)
    sigma = CoefficientFunction(sigma_list)
    inv_mu = CoefficientFunction(inv_mu_list)

    return mu, sigma, inv_mu


# ============================================================
# Measurement Points for Validation
# ============================================================
# Line A1-B1: y = 72mm, z = 34mm, x varies
# Line A2-B2: y = 144mm, z = 34mm, x varies
# Line A3-B3: y = 72mm, z = 0 (inside plate surface)

measurement_line_A1B1 = {
    'y': 72e-3,
    'z': 34e-3,
    'x_values': [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]  # mm
}

measurement_line_A2B2 = {
    'y': 144e-3,
    'z': 34e-3,
    'x_values': [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]
}

measurement_line_A3B3 = {
    'y': 72e-3,
    'z': 0,  # Inside plate
    'x_values': [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]
}


if __name__ == "__main__":
    print("=" * 70)
    print("TEAM Problem 7: Geometry Test")
    print("=" * 70)

    # Create geometry with moderate mesh
    mesh = create_team7_geometry(maxh=0.03)

    # Get material properties
    mu, sigma, inv_mu = get_material_properties(mesh)

    # Get coil current density
    J0 = get_coil_current_density(mesh)

    print("\n" + "=" * 70)
    print("Geometry created successfully!")
    print("=" * 70)
    print(f"""
Summary:
  Plate: 294 x 294 x 19 mm aluminum (σ = {sigma_al:.3e} S/m)
  Hole:  108 x 108 mm at (18,18) to (126,126) mm
  Coil:  Racetrack, {n_turns} turns, {I_coil:.1f} A
  Frequencies: {freq_50} Hz and {freq_200} Hz

Measurement lines:
  A1-B1: y=72mm, z=34mm (above plate)
  A2-B2: y=144mm, z=34mm (above plate)
  A3-B3: y=72mm, z=0 (plate surface)
""")
