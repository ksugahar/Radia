"""
3D A-formulation (Vector Potential) for Circular Coil with Kelvin Transformation
1/8 MODEL (Octant: x>=0, y>=0, z>=0)

Problem: Circular coil (cross-section 10x10mm, inner radius 50mm) in air

Coil parameters:
  - Inner radius: r_in = 50 mm = 0.05 m
  - Outer radius: r_out = 60 mm = 0.06 m
  - Height: 10 mm = 0.01 m (z = -5mm to +5mm, centered at z=0)
  - Cross-section: 10mm x 10mm (square)

Current density (analytical):
  J = J0 * (-y/rho, x/rho, 0)  where rho = sqrt(x^2 + y^2)
  This gives uniform current in the azimuthal direction.

A-formulation:
  curl(1/mu * curl(A)) = J
  H = (1/mu) * curl(A)
  B = curl(A)

Kelvin transformation for 3D:
  Maps infinite exterior domain to finite sphere interior
  mu'(r') = (R/r')^2 * mu0
  nu_kelvin = (r'/R)^2 * nu0

Note: Uses nograds=True to eliminate gradient kernel in HCurl space.
"""
import os
import sys
import glob

# Set environment for ksugahar's NGSolve build
os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ksugahar\Lib\site-packages')

# Delete existing output files
script_dir = os.path.dirname(os.path.abspath(__file__))
print("Deleting existing output files...")
deleted_count = 0
for ext in ['*.png', '*.vtu', '*.mat']:
    for f in glob.glob(os.path.join(script_dir, ext)):
        try:
            os.remove(f)
            print(f"  Deleted: {f}")
            deleted_count += 1
        except Exception as e:
            print(f"  Failed to delete {f}: {e}")
print(f"Deleted {deleted_count} files.")

from numpy import pi, sqrt, linspace, zeros, array, log10, isnan, nan
from numpy import cos, sin, meshgrid
from scipy.special import ellipk, ellipe
from ngsolve import *
from netgen.occ import *
import scipy.io as sio

# Import matplotlib for plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})
matplotlib.rcParams['font.family'] = 'Times New Roman'

print("=" * 60)
print("3D A-formulation (HCurl) with Kelvin Transform")
print("Circular Coil - 1/8 Model (Octant)")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
# Coil geometry (in meters)
r_coil_inner = 0.05    # Inner radius: 50 mm
r_coil_outer = 0.06    # Outer radius: 60 mm (10mm thickness)
coil_height = 0.01     # Height: 10 mm
coil_z_min = -coil_height / 2  # -5 mm
coil_z_max = coil_height / 2   # +5 mm

# Kelvin transformation radius
kelvin_radius = 0.07   # 70 mm

# Physical constants
mu0 = 4 * pi * 1e-7    # Vacuum permeability [H/m]
nu0 = 1.0 / mu0        # Vacuum reluctivity [m/H]
J0 = 1e6               # Current density magnitude [A/m^2]

# Mesh parameters
maxh_initial = 0.008   # Initial mesh size [m]
order = 2              # Finite element order

# Offset for exterior domain (x-direction)
offset_x = 2 * kelvin_radius

print(f"\nCoil parameters:")
print(f"  Inner radius: {r_coil_inner*1000:.1f} mm")
print(f"  Outer radius: {r_coil_outer*1000:.1f} mm")
print(f"  Height: {coil_height*1000:.1f} mm (z = {coil_z_min*1000:.1f} to {coil_z_max*1000:.1f} mm)")
print(f"  Cross-section: {(r_coil_outer-r_coil_inner)*1000:.1f} x {coil_height*1000:.1f} mm")
print(f"\nKelvin radius: {kelvin_radius*1000:.1f} mm")
print(f"Current density: J0 = {J0:.2e} A/m^2")


# ============================================================
# Analytical Solution (Biot-Savart for circular coil on axis)
# ============================================================
def Bz_analytical_on_axis(z, R, I):
    """
    Bz on axis for thin circular coil.
    R: coil radius, I: total current
    """
    return mu0 * I * R**2 / (2 * (R**2 + z**2)**(3/2))


# ============================================================
# Geometry Definition for 3D - 1/8 Model
# ============================================================
def create_geometry():
    """Create 1/8 geometry with Kelvin exterior (single domain approach)."""
    print("\nCreating 1/8 model geometry...")

    # ===== Interior domain: 1/8 sphere at origin =====
    air_inner_full = Sphere(Pnt(0, 0, 0), kelvin_radius)
    box_cut = Box(Pnt(0, 0, 0), Pnt(kelvin_radius, kelvin_radius, kelvin_radius))
    air_inner = air_inner_full * box_cut
    air_inner.mat("air_inner")

    # Name faces for air_inner
    for face in air_inner.faces:
        fc = face.center
        if abs(fc.x) < 1e-6:
            face.name = "sym_x"
        elif abs(fc.y) < 1e-6:
            face.name = "sym_y"
        elif abs(fc.z) < 1e-6:
            face.name = "sym_z"

    # ===== Exterior domain (Kelvin transformed): 1/8 sphere =====
    air_outer_full = Sphere(Pnt(offset_x, 0, 0), kelvin_radius)
    box_cut_ext = Box(Pnt(offset_x, 0, 0), Pnt(offset_x + kelvin_radius, kelvin_radius, kelvin_radius))
    air_outer = air_outer_full * box_cut_ext
    air_outer.mat("air_outer")

    # Name faces for air_outer
    for face in air_outer.faces:
        fc = face.center
        if abs(fc.x - offset_x) < 1e-6:
            face.name = "sym_x_ext"
        elif abs(fc.y) < 1e-6:
            face.name = "sym_y"
        elif abs(fc.z) < 1e-6:
            face.name = "sym_z"

    # ===== Identify periodic faces =====
    print("Identifying periodic boundaries...")

    kelvin_int_face = None
    kelvin_ext_face = None

    for face in air_inner.faces:
        fc = face.center
        dist = sqrt(fc.x**2 + fc.y**2 + fc.z**2)
        if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
            kelvin_int_face = face
            face.name = "kelvin_int"
            print(f"  Found kelvin_int face at center ({fc.x:.4f}, {fc.y:.4f}, {fc.z:.4f})")
            break

    for face in air_outer.faces:
        fc = face.center
        dist = sqrt((fc.x - offset_x)**2 + fc.y**2 + fc.z**2)
        if abs(dist - kelvin_radius) < kelvin_radius * 0.2:
            kelvin_ext_face = face
            face.name = "kelvin_ext"
            print(f"  Found kelvin_ext face at center ({fc.x:.4f}, {fc.y:.4f}, {fc.z:.4f})")
            break

    if kelvin_int_face is not None and kelvin_ext_face is not None:
        kelvin_ext_face.Identify(kelvin_int_face, "periodic", IdentificationType.PERIODIC)
        print("  Periodic identification applied!")
    else:
        print("  WARNING: Could not find periodic faces!")

    # GND vertex at center of exterior domain
    vertex = Vertex(Pnt(offset_x, 0, 0))
    vertex.name = "GND"

    # Glue all domains
    geo = Glue([air_inner, air_outer, vertex])

    return OCCGeometry(geo)


# ============================================================
# Solve A-formulation (HCurl)
# ============================================================
def solve_A_formulation(mesh, fe_order):
    """Solve A-formulation using HCurl elements with nograds=True."""

    # HCurl space with nograds=True to eliminate gradient kernel
    # Dirichlet BC on symmetry planes: n x A = 0
    fes_before = HCurl(mesh, order=fe_order, dirichlet="sym_x|sym_y|sym_x_ext", nograds=True)
    fes = Periodic(fes_before)

    print(f"  DOFs before Periodic: {fes_before.ndof}")
    print(f"  DOFs after Periodic: {fes.ndof}")
    print(f"  FreeDofs: {sum(fes.FreeDofs())}")

    A = fes.TrialFunction()
    v = fes.TestFunction()

    # Kelvin coefficient for exterior domain
    r_prime_sq = (x - offset_x)**2 + y**2 + z**2
    r_prime = sqrt(r_prime_sq + 1e-20)

    # nu_kelvin = (r'/R)^2 * nu0
    nu_kelvin = (r_prime / kelvin_radius)**2 * nu0

    # Material coefficients
    nu_dict = {
        "air_inner": nu0,
        "air_outer": nu_kelvin
    }
    Nu = CoefficientFunction([nu_dict[mat] for mat in mesh.GetMaterials()])

    # Bilinear form: (nu * curl(A), curl(v)) + regularization
    a = BilinearForm(fes)
    a += Nu * InnerProduct(curl(A), curl(v)) * dx
    a += 1e-6 * Nu * InnerProduct(A, v) * dx  # Small regularization
    c = Preconditioner(a, 'bddc')
    a.Assemble()

    # Current density: J = J0 * (-y/rho, x/rho, 0) in coil region
    # Define coil region using IfPos (only in air_inner domain)
    rho = sqrt(x**2 + y**2 + 1e-20)
    r_cyl = sqrt(x**2 + y**2)

    # Coil region: r_coil_inner < r < r_coil_outer and coil_z_min < z < coil_z_max
    in_coil_r = IfPos(r_cyl - r_coil_inner, 1, 0) * IfPos(r_coil_outer - r_cyl, 1, 0)
    in_coil_z = IfPos(z - coil_z_min, 1, 0) * IfPos(coil_z_max - z, 1, 0)
    in_coil = in_coil_r * in_coil_z

    Jx = -J0 * y / rho
    Jy = J0 * x / rho
    Jz = CF(0.0)
    J_cf = in_coil * CoefficientFunction((Jx, Jy, Jz))

    # Linear form: (J, v) only in air_inner region (where coil is)
    f = LinearForm(fes)
    f += InnerProduct(J_cf, v) * dx("air_inner")
    f.Assemble()

    # Solve using CG with BDDC preconditioner
    gfA = GridFunction(fes)

    from ngsolve.krylovspace import CGSolver
    inv = CGSolver(a.mat, c.mat, maxiter=1000, printrates=True, tol=1e-10)
    gfA.vec.data = inv * f.vec

    # Compute B = curl(A)
    B_cf = curl(gfA)

    # H field
    H_dict = {
        "air_inner": nu0 * curl(gfA),
        "air_outer": nu_kelvin * curl(gfA)
    }
    H_cf = CoefficientFunction([H_dict[mat] for mat in mesh.GetMaterials()])

    fields = {
        'B_cf': B_cf,
        'H_cf': H_cf,
        'J_cf': J_cf,
        'Nu': Nu,
        'nu0': nu0,
        'nu_kelvin': nu_kelvin
    }

    return fes, gfA, fields


# ============================================================
# VTK Output
# ============================================================
def output_vtk(mesh, gfA, fields, filename_base):
    """Output mesh and solution to VTK file."""
    import shutil
    import tempfile
    import gc

    temp_dir = tempfile.gettempdir()
    temp_vtk_path = os.path.join(temp_dir, filename_base)
    final_vtk_path = os.path.join(script_dir, filename_base + ".vtu")

    coefs = [gfA, fields['B_cf'], fields['H_cf']]
    names = ["A", "B", "H"]

    vtk = VTKOutput(mesh, coefs=coefs, names=names,
                    filename=temp_vtk_path, subdivision=0)
    vtk.Do()

    del vtk
    gc.collect()

    temp_vtu_file = temp_vtk_path + ".vtu"
    if os.path.exists(temp_vtu_file):
        shutil.copy2(temp_vtu_file, final_vtk_path)
        try:
            os.remove(temp_vtu_file)
        except:
            pass

    return final_vtk_path


# ============================================================
# Visualization
# ============================================================
def plot_results(mesh, gfA, fields):
    """Plot results on z=0 and y=0 cross-sections."""
    B_cf = fields['B_cf']

    fig, axes = plt.subplots(2, 2, figsize=(14, 12), dpi=150)

    # ===== Upper left: |B| on z=0 plane in interior =====
    ax1 = axes[0, 0]

    # Sample on grid
    x_grid = linspace(0.001, kelvin_radius * 0.95, 41)
    y_grid = linspace(0.001, kelvin_radius * 0.95, 41)
    xx, yy = meshgrid(x_grid, y_grid)

    B_mag = zeros(xx.shape)
    for iy in range(len(y_grid)):
        for ix in range(len(x_grid)):
            r = sqrt(x_grid[ix]**2 + y_grid[iy]**2)
            if r < kelvin_radius:
                try:
                    mip = mesh(x_grid[ix], y_grid[iy], 0.001)
                    B_val = B_cf(mip)
                    B_mag[iy, ix] = sqrt(B_val[0]**2 + B_val[1]**2 + B_val[2]**2)
                except:
                    B_mag[iy, ix] = nan

    B_mag_plot = B_mag.copy()
    B_mag_plot[isnan(B_mag_plot)] = 0

    levels = linspace(0, B_mag_plot.max() * 0.9, 20)
    if levels[-1] > 0:
        contour1 = ax1.contourf(xx * 1000, yy * 1000, B_mag_plot * 1000, levels=levels * 1000, cmap='jet', extend='max')
        cbar1 = plt.colorbar(contour1, ax=ax1, shrink=0.8)
        cbar1.set_label('$|\\mathbf{B}|$ [mT]')

    # Draw coil region
    theta = linspace(0, pi/2, 50)
    ax1.plot(r_coil_inner * cos(theta) * 1000, r_coil_inner * sin(theta) * 1000, 'w-', linewidth=2)
    ax1.plot(r_coil_outer * cos(theta) * 1000, r_coil_outer * sin(theta) * 1000, 'w-', linewidth=2)

    # Kelvin boundary
    ax1.plot(kelvin_radius * cos(theta) * 1000, kelvin_radius * sin(theta) * 1000, 'g--', linewidth=1.5, label='Kelvin boundary')

    ax1.set_xlim(0, kelvin_radius * 1000 * 1.05)
    ax1.set_ylim(0, kelvin_radius * 1000 * 1.05)
    ax1.set_aspect('equal')
    ax1.set_xlabel('$x$ [mm]')
    ax1.set_ylabel('$y$ [mm]')
    ax1.set_title('$|\\mathbf{B}|$ on $z=0$ plane (Interior)')
    ax1.legend(loc='upper right', fontsize=8)

    # ===== Upper right: Bz on z-axis =====
    ax2 = axes[0, 1]

    z_axis = linspace(0.001, kelvin_radius * 0.9, 51)
    Bz_numerical = zeros(len(z_axis))
    Bz_analytical = zeros(len(z_axis))

    # Effective current and average radius for analytical solution
    R_avg = (r_coil_inner + r_coil_outer) / 2
    A_coil = (r_coil_outer - r_coil_inner) * coil_height
    I_eff = J0 * A_coil

    r_sample = 0.002  # Slightly off axis
    for i, z_val in enumerate(z_axis):
        try:
            mip = mesh(r_sample, r_sample, z_val)
            B_val = B_cf(mip)
            Bz_numerical[i] = B_val[2]
        except:
            Bz_numerical[i] = nan

        Bz_analytical[i] = Bz_analytical_on_axis(z_val, R_avg, I_eff)

    # Scale for 1/8 model -> full model (factor of 1 since B is local)
    valid = ~isnan(Bz_numerical)
    ax2.plot(z_axis[valid] * 1000, Bz_numerical[valid] * 1000, 'b-',
             linewidth=2, label='Numerical (FEM)')
    ax2.plot(z_axis * 1000, Bz_analytical * 1000, 'r--',
             linewidth=2, label='Analytical (thin coil)')

    ax2.axvline(x=coil_z_max * 1000, color='gray', linestyle=':', linewidth=1,
                label='Coil boundary')

    ax2.set_xlabel('$z$ [mm]')
    ax2.set_ylabel('$B_z$ [mT]')
    ax2.set_title('$B_z$ on $z$-axis')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    # ===== Lower left: |B| on y=0 plane in interior =====
    ax3 = axes[1, 0]

    x_grid2 = linspace(0.001, kelvin_radius * 0.95, 41)
    z_grid2 = linspace(0.001, kelvin_radius * 0.95, 41)
    xx2, zz2 = meshgrid(x_grid2, z_grid2)

    B_mag2 = zeros(xx2.shape)
    for iz in range(len(z_grid2)):
        for ix in range(len(x_grid2)):
            r = sqrt(x_grid2[ix]**2 + z_grid2[iz]**2)
            if r < kelvin_radius:
                try:
                    mip = mesh(x_grid2[ix], 0.001, z_grid2[iz])
                    B_val = B_cf(mip)
                    B_mag2[iz, ix] = sqrt(B_val[0]**2 + B_val[1]**2 + B_val[2]**2)
                except:
                    B_mag2[iz, ix] = nan

    B_mag2_plot = B_mag2.copy()
    B_mag2_plot[isnan(B_mag2_plot)] = 0

    levels2 = linspace(0, B_mag2_plot.max() * 0.9, 20)
    if levels2[-1] > 0:
        contour3 = ax3.contourf(xx2 * 1000, zz2 * 1000, B_mag2_plot * 1000, levels=levels2 * 1000, cmap='jet', extend='max')
        cbar3 = plt.colorbar(contour3, ax=ax3, shrink=0.8)
        cbar3.set_label('$|\\mathbf{B}|$ [mT]')

    # Coil region on y=0 (visible as line at z=0, r_inner to r_outer)
    ax3.plot([r_coil_inner * 1000, r_coil_outer * 1000], [coil_z_max * 1000, coil_z_max * 1000], 'w-', linewidth=2)

    ax3.set_xlim(0, kelvin_radius * 1000 * 1.05)
    ax3.set_ylim(0, kelvin_radius * 1000 * 1.05)
    ax3.set_aspect('equal')
    ax3.set_xlabel('$x$ [mm]')
    ax3.set_ylabel('$z$ [mm]')
    ax3.set_title('$|\\mathbf{B}|$ on $y=0$ plane (Interior)')

    # ===== Lower right: Exterior domain =====
    ax4 = axes[1, 1]

    x_grid_ext = linspace(offset_x + 0.001, offset_x + kelvin_radius * 0.95, 31)
    z_grid_ext = linspace(0.001, kelvin_radius * 0.95, 31)
    xx_ext, zz_ext = meshgrid(x_grid_ext, z_grid_ext)

    B_mag_ext = zeros(xx_ext.shape)
    for iz in range(len(z_grid_ext)):
        for ix in range(len(x_grid_ext)):
            r = sqrt((x_grid_ext[ix] - offset_x)**2 + z_grid_ext[iz]**2)
            if r < kelvin_radius:
                try:
                    mip = mesh(x_grid_ext[ix], 0.001, z_grid_ext[iz])
                    B_val = B_cf(mip)
                    B_mag_ext[iz, ix] = sqrt(B_val[0]**2 + B_val[1]**2 + B_val[2]**2)
                except:
                    B_mag_ext[iz, ix] = nan

    B_mag_ext_plot = B_mag_ext.copy()
    B_mag_ext_plot[isnan(B_mag_ext_plot)] = 0

    levels_ext = linspace(0, B_mag_ext_plot.max() * 0.9, 20) if B_mag_ext_plot.max() > 0 else linspace(0, 1, 20)
    if levels_ext[-1] > 0:
        contour4 = ax4.contourf(xx_ext * 1000, zz_ext * 1000, B_mag_ext_plot * 1000, levels=levels_ext * 1000, cmap='jet', extend='max')
        cbar4 = plt.colorbar(contour4, ax=ax4, shrink=0.8)
        cbar4.set_label('$|\\mathbf{B}|$ [mT]')

    # Kelvin boundary
    theta2 = linspace(0, pi/2, 50)
    ax4.plot((offset_x + kelvin_radius * cos(theta2)) * 1000, kelvin_radius * sin(theta2) * 1000,
             'g--', linewidth=1.5, label='Kelvin boundary')

    ax4.set_xlim(offset_x * 1000, (offset_x + kelvin_radius) * 1000 * 1.05)
    ax4.set_ylim(0, kelvin_radius * 1000 * 1.05)
    ax4.set_aspect('equal')
    ax4.set_xlabel('$x$ [mm]')
    ax4.set_ylabel('$z$ [mm]')
    ax4.set_title('Exterior domain (Kelvin, $y=0$ plane)')
    ax4.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    png_file = os.path.join(script_dir, "CircularCoil_A_formulation_result.png")
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    print(f"  PNG saved: {png_file}")
    plt.close()

    return z_axis, Bz_numerical, Bz_analytical


# ============================================================
# Main
# ============================================================
print("\n" + "=" * 60)
print("Creating geometry...")
print("=" * 60)

geo = create_geometry()
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial, grading=0.3))
mesh.Curve(order)

print(f"\nMesh statistics:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {set(mesh.GetBoundaries())}")

# Count elements by region
for mat in set(mesh.GetMaterials()):
    count = sum(1 for el in mesh.Elements(VOL) if el.mat == mat)
    print(f"    {mat}: {count} elements")

print("\n" + "=" * 60)
print("Solving A-formulation...")
print("=" * 60)

fes, gfA, fields = solve_A_formulation(mesh, order)

print(f"\nSolution computed.")
print(f"  DOFs: {fes.ndof}")

# Compute energy
B_cf = fields['B_cf']

# Energy = (1/2) * integral( nu * |B|^2 dx )
energy_inner = Integrate(0.5 * nu0 * InnerProduct(B_cf, B_cf) * dx("air_inner"), mesh)

# For Kelvin region
r_prime_sq = (x - offset_x)**2 + y**2 + z**2
r_prime = sqrt(r_prime_sq + 1e-20)
nu_kelvin_cf = (r_prime / kelvin_radius)**2 * nu0
energy_outer = Integrate(0.5 * nu_kelvin_cf * InnerProduct(B_cf, B_cf) * dx("air_outer"), mesh)

energy = energy_inner + energy_outer
energy_full = 8 * energy  # Full model (8 octants)

print(f"\nMagnetic energy (1/8 model): {energy:.6e} J")
print(f"  air_inner: {energy_inner:.6e} J")
print(f"  air_outer: {energy_outer:.6e} J")
print(f"Magnetic energy (full model): {energy_full:.6e} J")

# Check Bz on axis
print("\nBz on z-axis:")
r_sample = 0.002
for z_val in [0.001, 0.01, 0.03, 0.05]:
    try:
        mip = mesh(r_sample, r_sample, z_val)
        B_val = B_cf(mip)
        Bz = B_val[2]
        print(f"  z={z_val*1000:.0f}mm: Bz = {Bz*1000:.4f} mT")
    except Exception as e:
        print(f"  z={z_val*1000:.0f}mm: Error - {e}")

# Analytical reference
R_avg = (r_coil_inner + r_coil_outer) / 2
A_coil = (r_coil_outer - r_coil_inner) * coil_height
I_eff = J0 * A_coil
Bz_center_ana = mu0 * I_eff * R_avg**2 / (2 * R_avg**3)
print(f"\nAnalytical Bz at center (thin coil): {Bz_center_ana*1000:.4f} mT")

# Output VTK
print("\n" + "=" * 60)
print("Saving outputs...")
print("=" * 60)

vtk_file = output_vtk(mesh, gfA, fields, "CircularCoil_A_formulation")
print(f"  VTK saved: {vtk_file}")

# Plot results
z_axis, Bz_numerical, Bz_analytical = plot_results(mesh, gfA, fields)

# Save data to MAT file
mat_data = {
    'r_coil_inner': r_coil_inner,
    'r_coil_outer': r_coil_outer,
    'coil_height': coil_height,
    'kelvin_radius': kelvin_radius,
    'J0': J0,
    'mu0': mu0,
    'order': order,
    'ndof': fes.ndof,
    'n_elements': mesh.ne,
    'energy_1_8': energy,
    'energy_full': energy_full,
    'z_axis': z_axis,
    'Bz_numerical': array(Bz_numerical),
    'Bz_analytical': Bz_analytical
}

mat_file = os.path.join(script_dir, "CircularCoil_A_formulation.mat")
sio.savemat(mat_file, mat_data)
print(f"  MAT saved: {mat_file}")

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
