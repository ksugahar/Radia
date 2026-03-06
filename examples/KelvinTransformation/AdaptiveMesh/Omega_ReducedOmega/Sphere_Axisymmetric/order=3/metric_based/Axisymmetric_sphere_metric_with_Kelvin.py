"""
Axisymmetric Omega Method for Magnetostatics with Kelvin Transformation
WITH METRIC-BASED REMESHING

Problem: Magnetic sphere (mu_r=100) in uniform z-directed background field

Formulation (Total Scalar Potential):
- H = -grad(Omega)
- From div(B) = 0: div(mu * grad(Omega)) = 0

Kelvin transformation:
- Maps infinite exterior domain to finite half-circle
- Permeability transformation: mu'(rho') = (R/rho')^2 * mu0
- Periodic BC couples interior (rho=R) with exterior (rho'=R)

Metric-based remeshing algorithm:
  1. Solve on current mesh
  2. Compute ZZ error estimator per element (H(div) recovery)
  3. For each element, compute ideal mesh size: h_ideal = h_current * (eta_target / eta)^(1/(p+1))
  4. Write local size field to meshsizefile and regenerate mesh
  5. Repeat until DOF limit reached

Stop condition: DOF >= max_dof
"""
import os
import glob
import tempfile

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

from numpy import pi, sqrt, linspace, zeros, nan, isnan, meshgrid, array, log, log10, sin, cos
from ngsolve import *
from netgen.occ import *
from netgen.meshing import MeshingParameters
import scipy.io as sio

# Import matplotlib for plotting
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend to avoid tkinter errors
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})
matplotlib.rcParams['font.family'] = 'Times New Roman'

print("=" * 60)
print("Axisymmetric Omega Method with Kelvin Transform")
print("METRIC-BASED REMESHING")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
sphere_radius = 0.5      # Magnetic sphere radius [m]
kelvin_radius = 1.0      # Kelvin transformation radius [m]
mu_r = 100               # Relative permeability
mu0 = 4 * pi * 1e-7      # Vacuum permeability [H/m]

# Source field: H_s = (0, 0, H0) uniform in z-direction
H0 = 1.0  # [A/m]

# Mesh parameters
maxh_initial = 0.5       # Initial mesh size (adjusted for ~100 elements)
order = 3                # Finite element order

# Offset for exterior domain (z-direction for axisymmetric)
offset_z = 3.0

# Metric-based remeshing parameters
max_iterations = 20      # Stop after 20 iterations
h_min = 0.001            # Minimum mesh size [m]
h_max = 0.5              # Maximum mesh size [m]
grading = 0.3            # Mesh grading parameter

# DOF growth target: 2x per iteration
# For 2D: DOF ~ h^(-2), so h_ratio = dof_growth^(-1/2) = 0.707 for 2x
dof_growth_target = 2.0
h_ratio_target = dof_growth_target ** (-0.5)  # ≈ 0.707 for 2x
eta_target_factor = h_ratio_target ** (order + 1)  # Computed based on order

print(f"\nProblem parameters:")
print(f"  Sphere radius: {sphere_radius} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")
print(f"\nMetric-based remeshing parameters:")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Stop condition: {max_iterations} iterations or DOF >= 1e6")
print(f"  h_min: {h_min} m, h_max: {h_max} m")
print(f"  Grading: {grading}")
print(f"  DOF growth target: {dof_growth_target}x per iteration")
print(f"  h_ratio_target: {h_ratio_target:.4f}")
print(f"  eta_target_factor: {eta_target_factor:.4f}")


# ============================================================
# Geometry Definition using HALF-CIRCLES (axisymmetric, r >= 0)
# ============================================================
def create_base_geometry():
    """Create geometry with periodic boundary conditions."""

    # ===== INTERIOR DOMAIN (half-circle, r >= 0) =====
    outer_half_int = Circle((0, 0), kelvin_radius).Face()
    cutter_int = MoveTo(-kelvin_radius-1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
    outer_half_int = outer_half_int - cutter_int  # Keep x >= 0 part

    inner_half_int = Circle((0, 0), sphere_radius).Face()
    inner_half_int = inner_half_int - cutter_int  # Keep x >= 0 part

    # Air region (half-annulus)
    air_inner = outer_half_int - inner_half_int

    # Name boundaries
    for edge in air_inner.edges:
        x_center = edge.center.x
        dist = sqrt(edge.center.x**2 + edge.center.y**2)
        if x_center < 1e-6:  # On z-axis (r = 0)
            edge.name = "axis_int"
        elif dist > (kelvin_radius + sphere_radius) / 2:
            edge.name = "kelvin_int"
        else:
            edge.name = "sphere"
    air_inner.faces.name = "air_inner"

    # Magnetic material (half-circle)
    for edge in inner_half_int.edges:
        if edge.center.x < 1e-6:
            edge.name = "axis_int"
        else:
            edge.name = "sphere"
    inner_half_int.faces.name = "magnetic"

    # ===== EXTERIOR DOMAIN (half-circle, z-offset for axisymmetric) =====
    outer_half_ext = Circle((0, offset_z), kelvin_radius).Face()
    cutter_ext = MoveTo(-kelvin_radius - 1, offset_z - kelvin_radius - 1).Rectangle(kelvin_radius + 1, 2 * kelvin_radius + 2).Face()
    outer_half_ext = outer_half_ext - cutter_ext  # Keep x >= 0 part (r >= 0)

    for edge in outer_half_ext.edges:
        r_center = edge.center.x  # r coordinate
        if r_center < 1e-6:  # On axis (r = 0)
            edge.name = "axis_ext"
        else:
            edge.name = "kelvin_ext"
    outer_half_ext.faces.name = "air_outer"

    # ===== GND VERTEX (center of exterior domain) =====
    vertex = Vertex(Pnt(0, offset_z, 0))
    vertex.name = "GND"

    # Glue all domains
    shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

    # ===== IDENTIFY ALL PERIODIC BOUNDARY EDGE PAIRS =====
    # Find all kelvin edges
    kelvin_int_edges = []
    kelvin_ext_edges = []
    for edge in shape.edges:
        if edge.name == "kelvin_int":
            kelvin_int_edges.append(edge)
        elif edge.name == "kelvin_ext":
            kelvin_ext_edges.append(edge)

    # Match edges by r-coordinate (x) of center
    # Interior: z > 0 or z < 0, Exterior: (z - offset_z) > 0 or < 0
    if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
        for int_edge in kelvin_int_edges:
            int_z = int_edge.center.y  # z in interior
            for ext_edge in kelvin_ext_edges:
                ext_z_local = ext_edge.center.y - offset_z  # z' in exterior (local coord)
                # Match: both positive z or both negative z
                if (int_z > 0 and ext_z_local > 0) or (int_z < 0 and ext_z_local < 0):
                    int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                    break

    return shape


# ============================================================
# Generate mesh with local size field
# ============================================================
def generate_mesh_with_local_sizes(shape, local_sizes=None, default_maxh=0.3):
    """
    Generate mesh with local mesh size specifications.
    """
    geo = OCCGeometry(shape, dim=2)
    mp = MeshingParameters(maxh=default_maxh, grading=grading)

    if local_sizes is not None and len(local_sizes) > 0:
        # Write mesh size file in Netgen format
        temp_dir = tempfile.gettempdir()
        meshsize_file = os.path.join(temp_dir, "meshsize_field_omega.txt")
        with open(meshsize_file, 'w') as f:
            f.write(f"{len(local_sizes)}\n")
            for (px, py, h) in local_sizes:
                h_clamped = max(h_min, min(h_max, h))
                f.write(f"{px} {py} 0 {h_clamped}\n")
            f.write("0\n")  # nr_edges = 0

        ngmesh = geo.GenerateMesh(mp=mp, meshsizefilename=meshsize_file)
    else:
        ngmesh = geo.GenerateMesh(mp=mp)

    mesh = Mesh(ngmesh)
    mesh.Curve(order)  # Curve the mesh for higher-order elements
    return mesh


# ============================================================
# Solve Omega formulation (Omega-Reduced Omega with Dirichlet Lifting)
# ============================================================
def solve_omega_formulation(mesh, fe_order):
    """
    Solve Omega-Reduced Omega formulation on given mesh.
    Following EMPY_Analysis/Omega_ReducedOmega.py convention.

    - Total region (magnetic): H = grad(Omega)
    - Reduced region (air): H = H_s - grad(Omega), Omega is perturbation
    """
    import numpy as np_arr

    # H1 space with Dirichlet BC at GND (infinity)
    fes_before = H1(mesh, order=fe_order, dirichlet_bbnd="GND")
    fes = Periodic(fes_before)  # Apply periodic BC

    # Trial and test functions
    Omega = fes.TrialFunction()
    psi = fes.TestFunction()

    # Coordinate functions (x = r, y = z)
    r_coord = x  # Radial coordinate
    z_coord = y  # Axial coordinate

    # Distance squared from exterior domain center (z-offset)
    rho_prime_sq = r_coord**2 + (z_coord - offset_z)**2

    # Transformed permeability for exterior domain (Kelvin)
    mu_kelvin = kelvin_radius**2 / (rho_prime_sq + 1e-20) * mu0

    Mu_dict = {
        "air_inner": mu0,
        "air_outer": mu_kelvin,
        "magnetic": mu_r * mu0
    }
    Mu = CoefficientFunction([Mu_dict[mat] for mat in mesh.GetMaterials()])

    # r-weight for axisymmetric formulation (same r for all domains with z-offset)
    r_weight = r_coord
    r_weight_inner = r_coord  # For boundary integrals

    # ========================================
    # Omega-Reduced Omega Formulation
    # (Following Omega_ReducedOmega.py)
    # ========================================

    # Source potential: Omega_s = H0 * z (so grad(Omega_s) = (0, H0) = H_s)
    Omega_s = H0 * z_coord

    # Source magnetic field H_s = (0, H0) and B_s = mu0 * H_s
    Hs = CoefficientFunction((0.0, H0))
    Bs = CoefficientFunction((0.0, mu0 * H0))

    # For Kelvin-transformed exterior domain:
    rho_prime = sqrt(rho_prime_sq + 1e-20)
    Hz_exterior = -(rho_prime / kelvin_radius)**2 * H0
    Bs_exterior = CoefficientFunction((0.0, mu0 * Hz_exterior))

    # ===== Bilinear Form =====
    a = BilinearForm(fes)
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("magnetic")
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_inner")
    a += Mu * grad(Omega) * grad(psi) * r_weight * dx("air_outer")
    a.Assemble()

    # ===== Dirichlet BC on Total/Reduced interface =====
    gfOmega = GridFunction(fes)
    gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # ===== Linear Form =====
    f = LinearForm(fes)
    # Source term in Reduced region only (NOT in Kelvin exterior)
    f += Mu * grad(gfOmega) * grad(psi) * r_weight * dx("air_inner")
    # NOTE: Do NOT include air_outer (Kelvin region) in source term - this is critical!
    f.Assemble()

    # Extract FreeDofs part of f
    fcut = np_arr.array(f.vec.FV())[fes.FreeDofs()]
    np_arr.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

    # Neumann boundary condition on sphere
    normal = specialcf.normal(mesh.dim)
    f += (normal * Bs) * psi * r_weight_inner * ds("sphere")
    f.Assemble()

    # Solve
    gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # ===== Compute full B field for error estimation =====
    # Following Omega_ReducedOmega.py post-processing
    fesOt = H1(mesh, order=fe_order, definedon="magnetic")
    fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")

    Ot = GridFunction(fesOt)
    Orr = GridFunction(fesOr)
    Oxr = GridFunction(fesOr)

    Ot.Set(gfOmega, VOL, definedon="magnetic")
    Orr.Set(gfOmega, VOL, definedon="air_inner|air_outer")
    Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # B field: Bt in Total, Br + Bs in Reduced
    Bt = grad(Ot) * Mu  # Total region
    Br = (grad(Orr) - grad(Oxr)) * mu0  # Reduced region perturbation

    # Source field Bs per region
    Bs_dict = {
        "air_inner": Bs,
        "air_outer": Bs_exterior,
        "magnetic": CoefficientFunction((0.0, 0.0))
    }
    Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

    # Total B field (for error estimation)
    BField = Bt + Br + Bs_cf

    return fes, gfOmega, Mu, r_weight, BField


# ============================================================
# Compute error estimator
# ============================================================
def compute_error_estimator(mesh, fes, BField):
    """
    Compute ZZ-type error estimator per element.
    Uses full B field (BField) for error estimation.
    """
    flux = BField  # Use full B field for error estimation

    # Total magnetic energy for normalization (using B field norm)
    r_weight = x  # Axisymmetric r-weight
    W = Integrate(InnerProduct(flux, flux) * r_weight, mesh)
    if abs(W) < 1e-20:
        W = 1.0

    # H(div) recovery (order - 1 for ZZ estimator)
    recovery_order = max(1, fes.globalorder - 1)
    fes_flux = HDiv(mesh, order=recovery_order)
    gf_flux = GridFunction(fes_flux)
    gf_flux.Set(flux)

    # Error = difference between direct and recovered flux
    err = (flux - gf_flux) * (flux - gf_flux) * r_weight / abs(W)
    element_errors = Integrate(err, mesh, element_wise=True)

    return element_errors


def compute_error_and_local_sizes(mesh, element_errors, iteration=0):
    """
    Compute ideal local mesh size for each element based on error distribution.

    Uses the a posteriori error estimate relationship:
      error ~ h^(p+1) for smooth solutions

    For optimal mesh (equidistributed errors):
      h_new = h_current * (eta_target / eta)^(1/(p+1))

    To ensure continuous DOF growth, we apply a global scaling factor
    based on the iteration number: h_scale = h_ratio_target^iteration

    iteration: current iteration number (0-indexed). Used to gradually increase refinement.
    """
    # Compute statistics
    nonzero_errors = [e for e in element_errors if e > 1e-20]
    if nonzero_errors:
        eta_mean = sum(nonzero_errors) / len(nonzero_errors)
        eta_max = max(nonzero_errors)
    else:
        eta_mean = 1.0
        eta_max = 1.0

    # Target error per element
    eta_target = eta_mean * eta_target_factor

    local_sizes = []
    exponent = 1.0 / (order + 1)

    # Global scale factor: shrink mesh by h_ratio_target each iteration
    # This ensures DOF grows by ~2x each iteration
    global_scale = h_ratio_target ** iteration

    for el_idx, el in enumerate(mesh.Elements()):
        verts = el.vertices
        pts = [mesh[v].point for v in verts]
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)

        # Compute current element size from area
        if len(pts) >= 3:
            x0, y0 = pts[0][0], pts[0][1]
            x1, y1 = pts[1][0], pts[1][1]
            x2, y2 = pts[2][0], pts[2][1]
            area = 0.5 * abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))
            h_current = sqrt(area)
        else:
            h_current = 0.1

        if h_current < 1e-10:
            h_current = 0.1

        eta = element_errors[el_idx]

        # Compute error-based relative size (for distribution within mesh)
        if eta > 1e-20:
            ratio = eta_target / eta
            # Wider range to allow more aggressive refinement
            ratio = max(0.1, min(5.0, ratio))
            h_relative = ratio ** exponent  # Relative size factor (< 1 for high error)
        else:
            # Very small error - keep same relative size
            h_relative = 1.0

        # Apply global scaling to achieve target DOF growth
        # h_base shrinks by h_ratio_target each iteration
        h_base = maxh_initial * global_scale

        # Combine: global scale determines overall size, error determines distribution
        # High error elements get smaller size, low error elements get larger size
        h_ideal = h_base * h_relative
        h_ideal = max(h_min, min(h_max, h_ideal))
        local_sizes.append((cx, cy, h_ideal))

    # Add boundary size control gradually
    # Start with few points and increase over iterations to avoid sudden mesh jumps
    if iteration >= 1:
        # Gradually increase number of points and decrease boundary size
        n_pts = min(4 + iteration, 12)  # 5, 6, 7, ... up to 12 points
        h_boundary = max(h_min, h_min * 2 * (0.9 ** iteration))  # Gradually decrease

        for theta in linspace(0, pi, n_pts, endpoint=False):
            r = sphere_radius * sin(theta)
            z = sphere_radius * cos(theta)
            if r > 0.01:
                local_sizes.append((r * 1.05, z, h_boundary))
                local_sizes.append((r * 0.95, z, h_boundary))

    total_error = sqrt(sum(element_errors))

    return local_sizes, total_error


# ============================================================
# VTK Output
# ============================================================
def output_vtk(mesh, gfu, iteration):
    """Output mesh and solution to VTK file."""
    vtk_filename = os.path.join(script_dir, f"axi_omega_metric_iter_{iteration:02d}")
    vtk = VTKOutput(mesh, coefs=[gfu], names=["Omega_pert"],
                    filename=vtk_filename, subdivision=0)
    vtk.Do()
    return vtk_filename + ".vtu"


# ============================================================
# Main: Metric-based Remeshing Loop
# ============================================================
print("\n" + "=" * 60)
print("Creating base geometry...")
print("=" * 60)

base_shape = create_base_geometry()

# Initial mesh
print("\nGenerating initial mesh...")
mesh = generate_mesh_with_local_sizes(base_shape, local_sizes=None, default_maxh=maxh_initial)

print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# History tracking
history = {
    'ndof': [],
    'elements': [],
    'error': [],
    'Hz_error_percent': [],
    'energy': []
}

# Analytical solution for interior (total field)
Hz_total_analytical = 3.0 / (mu_r + 2) * H0

# Analytical perturbation field energy (3D full energy, NOT divided by 2*pi)
#
# 3D analytical energy:
# - Inside sphere: W_in = (1/2) * mu_r * mu0 * |H_pert|^2 * (4*pi/3) * a^3
# - Outside sphere (dipole): W_out = mu0 * m^2 / (12*pi*a^3) where m = 4*pi*a^3*(mu_r-1)/(mu_r+2)*H0
#
# H_pert inside sphere: -(mu_r-1)/(mu_r+2) * H0
H_pert_analytical = -(mu_r - 1) / (mu_r + 2) * H0

# 3D energies (NOT divided by 2*pi)
V_sphere = (4.0/3.0) * pi * sphere_radius**3
W_in_analytical = 0.5 * mu_r * mu0 * H_pert_analytical**2 * V_sphere
m_dipole = 4 * pi * sphere_radius**3 * (mu_r - 1) / (mu_r + 2) * H0
W_out_analytical = mu0 * m_dipole**2 / (12 * pi * sphere_radius**3)
W_analytical = W_in_analytical + W_out_analytical

# Theta for circles
theta_circle = linspace(0, pi, 100)


def draw_sphere_and_kelvin(ax, z_offset=0, show_legend=True):
    """Draw sphere and Kelvin boundary."""
    r_sphere = sphere_radius * sin(theta_circle)
    z_sphere = sphere_radius * cos(theta_circle)
    if z_offset == 0:
        ax.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
        label1 = 'Magnetic sphere' if show_legend else None
        ax.plot(r_sphere, z_sphere, 'r-', linewidth=2, label=label1)

    r_kelvin = kelvin_radius * sin(theta_circle)
    z_kelvin = kelvin_radius * cos(theta_circle) + z_offset
    label2 = 'Kelvin boundary' if show_legend else None
    ax.plot(r_kelvin, z_kelvin, 'g--', linewidth=1.5, label=label2)


def generate_iteration_plot(iter_num, mesh, gfu, element_errors, history, n_elements, n_vertices):
    """Generate 2x2 plot for current iteration and save MAT file with plotting data."""
    fig = plt.figure(figsize=(14, 14), dpi=150)

    # Separate inner and outer domain elements with error information
    inner_error_map = []  # List of (polygon, log_err)
    outer_error_map = []
    inner_elements = []
    outer_elements = []
    inner_el_idx = []
    outer_el_idx = []

    for el_idx, el in enumerate(mesh.Elements()):
        el_verts = [mesh[v].point for v in el.vertices]
        cx = sum(v[0] for v in el_verts) / len(el_verts)
        cy = sum(v[1] for v in el_verts) / len(el_verts)

        # Get element error
        err_val = element_errors[el_idx]
        if err_val > 1e-20:
            log_err = log10(err_val)
        else:
            log_err = -20

        poly = [(v[0], v[1]) for v in el_verts]

        if cy < offset_z / 2:
            inner_elements.append(el_verts)
            inner_el_idx.append(el_idx)
            inner_error_map.append((poly, log_err))
        else:
            outer_elements.append(el_verts)
            outer_el_idx.append(el_idx)
            outer_error_map.append((poly, log_err))

    # ===== Compute B field and stream function Psi for flux lines =====
    # In axisymmetric (r, z) coordinates:
    # B = mu * H = mu * (-grad(Omega) + H_s)
    # Stream function Psi: B_r = (1/r) dPsi/dz, B_z = -(1/r) dPsi/dr
    # Flux lines are contours of Psi = r * A_phi

    n_grid = 100
    r_inner = linspace(0.01, kelvin_radius - 0.01, n_grid)
    z_inner = linspace(-kelvin_radius + 0.01, kelvin_radius - 0.01, n_grid)
    R_grid, Z_grid = meshgrid(r_inner, z_inner)
    Psi_inner = zeros((n_grid, n_grid))
    Omega_inner = zeros((n_grid, n_grid))

    dr = r_inner[1] - r_inner[0]
    dz = z_inner[1] - z_inner[0]

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_grid[i, j], Z_grid[i, j]
            if r_pt**2 + z_pt**2 < kelvin_radius**2 and r_pt > 0:
                try:
                    mip = mesh(r_pt, z_pt)
                    Omega_inner[i, j] = gfu(mip)
                except:
                    Omega_inner[i, j] = nan
            else:
                Omega_inner[i, j] = nan

    # Compute gradients
    dOmega_dr = zeros((n_grid, n_grid))
    dOmega_dz = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            if not isnan(Omega_inner[i, j]):
                if j > 0 and j < n_grid - 1 and not isnan(Omega_inner[i, j-1]) and not isnan(Omega_inner[i, j+1]):
                    dOmega_dr[i, j] = (Omega_inner[i, j+1] - Omega_inner[i, j-1]) / (2 * dr)
                elif j > 0 and not isnan(Omega_inner[i, j-1]):
                    dOmega_dr[i, j] = (Omega_inner[i, j] - Omega_inner[i, j-1]) / dr
                elif j < n_grid - 1 and not isnan(Omega_inner[i, j+1]):
                    dOmega_dr[i, j] = (Omega_inner[i, j+1] - Omega_inner[i, j]) / dr
                if i > 0 and i < n_grid - 1 and not isnan(Omega_inner[i-1, j]) and not isnan(Omega_inner[i+1, j]):
                    dOmega_dz[i, j] = (Omega_inner[i+1, j] - Omega_inner[i-1, j]) / (2 * dz)
                elif i > 0 and not isnan(Omega_inner[i-1, j]):
                    dOmega_dz[i, j] = (Omega_inner[i, j] - Omega_inner[i-1, j]) / dz
                elif i < n_grid - 1 and not isnan(Omega_inner[i+1, j]):
                    dOmega_dz[i, j] = (Omega_inner[i+1, j] - Omega_inner[i, j]) / dz

    # Compute B field
    B_r_inner = zeros((n_grid, n_grid))
    B_z_inner = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_grid[i, j], Z_grid[i, j]
            if not isnan(Omega_inner[i, j]):
                rho = sqrt(r_pt**2 + z_pt**2)
                if rho < sphere_radius:
                    mu_val = mu_r * mu0
                else:
                    mu_val = mu0
                H_r = -dOmega_dr[i, j]
                H_z = -dOmega_dz[i, j] + H0
                B_r_inner[i, j] = mu_val * H_r
                B_z_inner[i, j] = mu_val * H_z
            else:
                B_r_inner[i, j] = nan
                B_z_inner[i, j] = nan

    # Compute stream function Psi by integrating B_r in z-direction
    # Psi satisfies: B_r = (1/r) dPsi/dz, B_z = -(1/r) dPsi/dr
    # Integrate along z (row index i): Psi(r,z) = Psi(r,z0) + integral(r * B_r dz)
    for j in range(n_grid):  # For each r value (column)
        # Start from bottom (i=0, z=-kelvin_radius)
        Psi_inner[0, j] = 0.5 * R_grid[0, j]**2 * B_z_inner[0, j] if not isnan(B_z_inner[0, j]) else 0
        for i in range(1, n_grid):  # Integrate along z
            if not isnan(B_r_inner[i, j]) and not isnan(B_r_inner[i-1, j]):
                B_r_avg = 0.5 * (B_r_inner[i, j] + B_r_inner[i-1, j])
                r_val = R_grid[i, j]
                Psi_inner[i, j] = Psi_inner[i-1, j] + r_val * B_r_avg * dz
            elif not isnan(Psi_inner[i-1, j]):
                Psi_inner[i, j] = Psi_inner[i-1, j]
            else:
                Psi_inner[i, j] = nan

    for i in range(n_grid):
        for j in range(n_grid):
            if isnan(Omega_inner[i, j]):
                Psi_inner[i, j] = nan

    r_outer = linspace(0.01, kelvin_radius - 0.01, n_grid)
    z_outer = linspace(offset_z - kelvin_radius + 0.01, offset_z + kelvin_radius - 0.01, n_grid)
    R_out_grid, Z_out_grid = meshgrid(r_outer, z_outer)
    Psi_outer = zeros((n_grid, n_grid))
    Omega_outer = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_out_grid[i, j], Z_out_grid[i, j]
            if r_pt**2 + (z_pt - offset_z)**2 < kelvin_radius**2 and r_pt > 0:
                try:
                    mip = mesh(r_pt, z_pt)
                    Omega_outer[i, j] = gfu(mip)
                except:
                    Omega_outer[i, j] = nan
            else:
                Omega_outer[i, j] = nan

    # Compute B field for outer domain
    dOmega_dr_out = zeros((n_grid, n_grid))
    dOmega_dz_out = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            if not isnan(Omega_outer[i, j]):
                if j > 0 and j < n_grid - 1 and not isnan(Omega_outer[i, j-1]) and not isnan(Omega_outer[i, j+1]):
                    dOmega_dr_out[i, j] = (Omega_outer[i, j+1] - Omega_outer[i, j-1]) / (2 * dr)
                elif j > 0 and not isnan(Omega_outer[i, j-1]):
                    dOmega_dr_out[i, j] = (Omega_outer[i, j] - Omega_outer[i, j-1]) / dr
                elif j < n_grid - 1 and not isnan(Omega_outer[i, j+1]):
                    dOmega_dr_out[i, j] = (Omega_outer[i, j+1] - Omega_outer[i, j]) / dr
                if i > 0 and i < n_grid - 1 and not isnan(Omega_outer[i-1, j]) and not isnan(Omega_outer[i+1, j]):
                    dOmega_dz_out[i, j] = (Omega_outer[i+1, j] - Omega_outer[i-1, j]) / (2 * dz)
                elif i > 0 and not isnan(Omega_outer[i-1, j]):
                    dOmega_dz_out[i, j] = (Omega_outer[i, j] - Omega_outer[i-1, j]) / dz
                elif i < n_grid - 1 and not isnan(Omega_outer[i+1, j]):
                    dOmega_dz_out[i, j] = (Omega_outer[i+1, j] - Omega_outer[i, j]) / dz

    B_r_outer = zeros((n_grid, n_grid))
    B_z_outer = zeros((n_grid, n_grid))

    for i in range(n_grid):
        for j in range(n_grid):
            r_pt, z_pt = R_out_grid[i, j], Z_out_grid[i, j]
            z_local = z_pt - offset_z
            if not isnan(Omega_outer[i, j]):
                rho_prime = sqrt(r_pt**2 + z_local**2)
                mu_val = kelvin_radius**2 / (rho_prime**2 + 1e-20) * mu0
                H_s_z = -(rho_prime / kelvin_radius)**2 * H0
                H_r = -dOmega_dr_out[i, j]
                H_z = -dOmega_dz_out[i, j] + H_s_z
                B_r_outer[i, j] = mu_val * H_r
                B_z_outer[i, j] = mu_val * H_z
            else:
                B_r_outer[i, j] = nan
                B_z_outer[i, j] = nan

    # Compute stream function Psi for outer domain by integrating B_r in z-direction
    for j in range(n_grid):  # For each r value (column)
        r_val = R_out_grid[0, j]
        Psi_outer[0, j] = 0.5 * r_val**2 * B_z_outer[0, j] if not isnan(B_z_outer[0, j]) else 0
        for i in range(1, n_grid):  # Integrate along z
            if not isnan(B_r_outer[i, j]) and not isnan(B_r_outer[i-1, j]):
                B_r_avg = 0.5 * (B_r_outer[i, j] + B_r_outer[i-1, j])
                r_val = R_out_grid[i, j]
                Psi_outer[i, j] = Psi_outer[i-1, j] + r_val * B_r_avg * dz
            elif not isnan(Psi_outer[i-1, j]):
                Psi_outer[i, j] = Psi_outer[i-1, j]
            else:
                Psi_outer[i, j] = nan

    for i in range(n_grid):
        for j in range(n_grid):
            if isnan(Omega_outer[i, j]):
                Psi_outer[i, j] = nan

    # ===== Save MAT file with all plotting data =====
    inner_el_verts_r = []
    inner_el_verts_z = []
    for el_verts in inner_elements:
        inner_el_verts_r.append([v[0] for v in el_verts])
        inner_el_verts_z.append([v[1] for v in el_verts])

    outer_el_verts_r = []
    outer_el_verts_z = []
    for el_verts in outer_elements:
        outer_el_verts_r.append([v[0] for v in el_verts])
        outer_el_verts_z.append([v[1] for v in el_verts])

    inner_errors = array([element_errors[i] for i in inner_el_idx]) if len(inner_el_idx) > 0 else array([])

    mat_iter_data = {
        'iter_num': iter_num,
        'n_elements': n_elements,
        'n_vertices': n_vertices,
        'ndof': history['ndof'][-1],
        'error': history['error'][-1],
        'Hz_error_percent': history['Hz_error_percent'][-1],
        'history_ndof': array(history['ndof']),
        'history_error': array(history['error']),
        'history_Hz_error_percent': array(history['Hz_error_percent']),
        'inner_n_elements': len(inner_elements),
        'inner_el_verts_r': array(inner_el_verts_r, dtype=object),
        'inner_el_verts_z': array(inner_el_verts_z, dtype=object),
        'inner_element_errors': inner_errors,
        'outer_n_elements': len(outer_elements),
        'outer_el_verts_r': array(outer_el_verts_r, dtype=object),
        'outer_el_verts_z': array(outer_el_verts_z, dtype=object),
        'R_grid_inner': R_grid,
        'Z_grid_inner': Z_grid,
        'Omega_inner': Omega_inner,
        'Psi_inner': Psi_inner,
        'B_r_inner': B_r_inner,
        'B_z_inner': B_z_inner,
        'R_grid_outer': R_out_grid,
        'Z_grid_outer': Z_out_grid,
        'Omega_outer': Omega_outer,
        'Psi_outer': Psi_outer,
        'B_r_outer': B_r_outer,
        'B_z_outer': B_z_outer,
        'sphere_radius': sphere_radius,
        'kelvin_radius': kelvin_radius,
        'offset_z': offset_z,
        'order': order
    }
    mat_iter_file = os.path.join(script_dir, f"axi_omega_metric_iter_{iter_num:02d}.mat")
    sio.savemat(mat_iter_file, mat_iter_data)
    print(f"  MAT saved: {mat_iter_file}")

    # Colormap settings for error visualization
    norm = Normalize(vmin=-12, vmax=-3)
    cmap_obj = plt.colormaps['jet']

    # ===== Top-left: Inner domain ZZ error map (no flux lines, white mesh edges) =====
    ax1 = plt.subplot(2, 2, 1)

    # Draw error colormap using PolyCollection with white edges
    if len(inner_error_map) > 0:
        polygons = [item[0] for item in inner_error_map]
        err_values = array([item[1] for item in inner_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='white', linewidths=0.3, alpha=0.9)
        ax1.add_collection(pc)

        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar1 = plt.colorbar(sm, ax=ax1)
        cbar1.set_label('$\\log_{10}$(ZZ Error)')

    draw_sphere_and_kelvin(ax1, 0, show_legend=True)
    ax1.set_xlabel('$r$ (m)')
    ax1.set_ylabel('$z$ (m)')
    ax1.set_title(f'Inner Domain: ZZ Error ({len(inner_elements)} elements)')
    ax1.set_aspect('equal')
    ax1.set_xlim(-0.05, kelvin_radius + 0.05)
    ax1.set_ylim(-kelvin_radius - 0.05, kelvin_radius + 0.05)

    # ===== Top-right: Outer domain ZZ error map (no flux lines, white mesh edges) =====
    ax2 = plt.subplot(2, 2, 2)

    # Draw error colormap using PolyCollection with white edges
    if len(outer_error_map) > 0:
        polygons = [item[0] for item in outer_error_map]
        err_values = array([item[1] for item in outer_error_map])
        colors = cmap_obj(norm(err_values))
        pc = PolyCollection(polygons, facecolors=colors, edgecolors='white', linewidths=0.3, alpha=0.9)
        ax2.add_collection(pc)

        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar2 = plt.colorbar(sm, ax=ax2)
        cbar2.set_label('$\\log_{10}$(ZZ Error)')

    draw_sphere_and_kelvin(ax2, offset_z, show_legend=False)
    ax2.set_xlabel("$r'$ (m)")
    ax2.set_ylabel("$z'$ (m)")
    ax2.set_title(f'Outer Domain: ZZ Error ({len(outer_elements)} elements)')
    ax2.set_aspect('equal')
    ax2.set_xlim(-0.05, kelvin_radius + 0.05)
    ax2.set_ylim(offset_z - kelvin_radius - 0.05, offset_z + kelvin_radius + 0.05)

    # ===== Bottom-left: DOF vs Error convergence curve =====
    ax3 = plt.subplot(2, 2, 3)
    ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Metric-based')
    ax3.loglog(history['ndof'][-1], history['error'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)

    ndof_line = array([1e2, 1e6])
    err_ref_point = history['error'][0]
    N_ref = history['ndof'][0]
    err_line = err_ref_point * (N_ref / ndof_line) ** (order / 2)
    ax3.loglog(ndof_line, err_line, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/2}})$')

    ax3.set_xlabel('DOFs')
    ax3.set_ylabel('Error Estimator')
    ax3.set_title(f'DOF vs Error (iter {iter_num})')
    ax3.set_xlim(1e2, 1e6)
    ax3.set_ylim(1e-4, 1)
    ax3.legend(loc='lower left')
    ax3.grid(True, alpha=0.3)
    ax3.tick_params(direction='in')

    # ===== Bottom-right: DOF vs Magnetic Energy =====
    ax4 = plt.subplot(2, 2, 4)
    ax4.semilogx(history['ndof'], history['energy'], 'gs-', linewidth=2, markersize=6, label='Numerical')
    ax4.axhline(y=W_analytical, color='r', linestyle='--', linewidth=2, label=f'Analytical: {W_analytical:.4e} J')
    ax4.semilogx(history['ndof'][-1], history['energy'][-1], 'ro', markersize=12, markerfacecolor='none', markeredgewidth=2)
    ax4.set_xlabel('DOFs')
    ax4.set_ylabel('Magnetic Energy (J)')
    ax4.set_title('DOF vs Perturbation Field Energy')
    ax4.set_xlim(1e2, 1e6)
    ax4.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0), useOffset=False)
    ax4.legend(loc='lower right')
    ax4.grid(True, alpha=0.3)
    ax4.tick_params(direction='in')

    plt.tight_layout()

    fig.suptitle(f'Iteration {iter_num}: {n_elements} elements, DOFs={history["ndof"][-1]}, Error={history["error"][-1]:.2e}',
                 fontsize=14, fontweight='bold', y=1.02)

    png_file = os.path.join(script_dir, f"axi_omega_metric_iter_{iter_num:02d}.png")
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    print(f"  PNG saved: {png_file}")
    plt.close()

print("\n" + "=" * 60)
print("Starting Metric-Based Remeshing")
print("=" * 60)
print(f"  Analytical energy (3D): {W_analytical:.6e} J")

iteration = 0
prev_ndof = 0  # Track DOF from previous iteration
while True:
    # Check DOF limit at loop start (stop if previous DOF >= 1e6)
    if prev_ndof >= 1e6:
        print(f"\n  DOF limit reached ({prev_ndof} >= 1e6), stopping without computing.")
        break

    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}")
    print("=" * 60)

    # Solve
    fes, gfu, mu, r_weight, BField = solve_omega_formulation(mesh, order)

    # Compute error estimator
    element_errors = compute_error_estimator(mesh, fes, BField)
    total_error = sqrt(sum(element_errors))

    # Compute Hz error inside sphere
    # In Total region (magnetic): H = grad(Omega), so Hz = grad(gfu)[1] directly
    try:
        Hz_total_val = grad(gfu)[1](mesh(0.1, 0))
        Hz_error_percent = abs(Hz_total_val - Hz_total_analytical) / abs(Hz_total_analytical) * 100
    except:
        Hz_error_percent = nan

    # Compute magnetic energy of perturbation field
    # Following 3D_sphere_with_Kelvin.py implementation
    #
    # Total region (magnetic): H_pert = grad(Omega) - H_s, B_pert = mu * H_pert
    # Reduced region (air_inner): H_pert = grad(Omega_r) = grad(Orr) - grad(Oxr)
    # Kelvin region (air_outer): H_pert = grad(Orr)
    #
    # Energy: W = (1/2) * integral(B_pert * H_pert) dV

    # Source field H_s = grad(Omega_s) = (0, H0)
    H_s_cf = CoefficientFunction((0.0, H0))

    # Need separate GridFunctions for accurate perturbation field in Reduced region
    # Following 3D_sphere_with_Kelvin.py lines 360-373
    from ngsolve import VOL
    Omega_s = H0 * y  # y is z in axisymmetric

    fesOt_energy = H1(mesh, order=fes.globalorder, definedon="magnetic")
    fesOr_energy = H1(mesh, order=fes.globalorder, definedon="air_inner|air_outer")

    Ot_energy = GridFunction(fesOt_energy)
    Orr_energy = GridFunction(fesOr_energy)
    Oxr_energy = GridFunction(fesOr_energy)

    Ot_energy.Set(gfu, VOL, definedon="magnetic")
    Orr_energy.Set(gfu, VOL, definedon="air_inner|air_outer")
    Oxr_energy.Set(Omega_s, BND, mesh.Boundaries("sphere"))

    # --- Total region (magnetic sphere) ---
    # H_pert = grad(Omega) - H_s
    # B_pert = mu_r * mu0 * H_pert
    # Energy = (1/2) * mu_r * mu0 * |H_pert|^2
    # Axisymmetric integration with r-weight gives W_2D, multiply by 2*pi for 3D energy
    H_pert_total = grad(Ot_energy) - H_s_cf
    W_total = 2 * pi * 0.5 * (mu_r * mu0) * Integrate(InnerProduct(H_pert_total, H_pert_total) * r_weight,
                                                       mesh, definedon=mesh.Materials("magnetic"))

    # --- Reduced region (air_inner) ---
    # H_pert = grad(Orr) - grad(Oxr) (perturbation from source potential)
    H_pert_reduced = grad(Orr_energy) - grad(Oxr_energy)
    W_reduced = 2 * pi * 0.5 * mu0 * Integrate(InnerProduct(H_pert_reduced, H_pert_reduced) * r_weight,
                                                mesh, definedon=mesh.Materials("air_inner"))

    # --- Kelvin region (air_outer) ---
    # H_pert = grad(Orr)
    # Using mu_kelvin for Kelvin-transformed region
    rho_prime_sq_energy = x**2 + (y - offset_z)**2
    mu_kelvin_energy = kelvin_radius**2 / (rho_prime_sq_energy + 1e-20) * mu0
    H_pert_kelvin = grad(Orr_energy)
    W_kelvin = 2 * pi * 0.5 * Integrate(mu_kelvin_energy * InnerProduct(H_pert_kelvin, H_pert_kelvin) * r_weight,
                                         mesh, definedon=mesh.Materials("air_outer"))

    # Total perturbation energy (3D, not divided by 2*pi)
    W_energy = W_total + W_reduced + W_kelvin

    # Record history
    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)
    history['Hz_error_percent'].append(Hz_error_percent)
    history['energy'].append(W_energy)

    # Compute local sizes from error distribution
    local_sizes, _ = compute_error_and_local_sizes(mesh, element_errors, iteration)
    h_values = [h for (_, _, h) in local_sizes if h < h_max]
    if h_values:
        h_min_actual = min(h_values)
        h_max_actual = max(h_values)
        h_mean = sum(h_values) / len(h_values)
    else:
        h_min_actual, h_max_actual, h_mean = h_min, h_max, (h_min + h_max) / 2

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")
    print(f"  Hz error inside sphere: {Hz_error_percent:.4f}%")
    print(f"  Energy (3D): {W_energy:.6e} J (analytical: {W_analytical:.6e} J, ratio: {W_energy/W_analytical:.4f})")
    print(f"  Size field: h_min={h_min_actual:.5f}, h_max={h_max_actual:.5f}, h_mean={h_mean:.5f}")

    # Output VTK
    vtk_file = output_vtk(mesh, gfu, iteration)
    print(f"  VTK saved: {vtk_file}")

    # Generate PNG and MAT for this iteration (with current mesh and solution)
    generate_iteration_plot(iteration + 1, mesh, gfu, element_errors, history, mesh.ne, mesh.nv)

    # Check iteration limit
    if iteration + 1 >= max_iterations:
        print(f"\n  Iteration limit reached ({iteration + 1} >= {max_iterations}), stopping.")
        break

    # Update prev_ndof for next iteration's DOF check
    prev_ndof = fes.ndof

    # Regenerate mesh with error-based local sizes
    new_maxh = h_max_actual
    new_maxh = max(h_min, min(h_max, new_maxh))
    base_shape = create_base_geometry()
    mesh = generate_mesh_with_local_sizes(base_shape, local_sizes, default_maxh=new_maxh)
    print(f"  New mesh: {mesh.ne} elements, {mesh.nv} vertices (default maxh={new_maxh:.4f})")

    iteration += 1


# ============================================================
# Final Validation
# ============================================================
print("\n" + "=" * 60)
print("FINAL VALIDATION")
print("=" * 60)

print(f"\nAnalytical Hz_total inside sphere: {Hz_total_analytical:.6f} A/m")

print("\n" + "=" * 60)
print("INTERIOR (magnetic sphere):")
print("In Total region: H = grad(Omega)")
print("=" * 60)

for r_val in [0.1, 0.2, 0.3, 0.4]:
    try:
        # In Total region: H = grad(Omega), so Hz = grad(gfu)[1] directly
        Hz_total_val = grad(gfu)[1](mesh(r_val, 0))
        err = abs(Hz_total_val - Hz_total_analytical) / abs(Hz_total_analytical) * 100
        print(f"  r={r_val}: Hz_total = {Hz_total_val:.6f}, error = {err:.4f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")

print("\n" + "=" * 60)
print("EXTERIOR (air region):")
print("In Reduced region: H = H_s + grad(Omega_r) - grad(Omega_s)")
print("=" * 60)

# Post-processing for exterior - need to compute H = H_s + grad(Omega_r) - grad(Omega_s)
fesOr_val = H1(mesh, order=order, definedon="air_inner|air_outer")
Orr_val = GridFunction(fesOr_val)
Oxr_val = GridFunction(fesOr_val)
Omega_s_val = H0 * y  # y is z in axisymmetric

Orr_val.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr_val.Set(Omega_s_val, BND, mesh.Boundaries("sphere"))

for r_val in [0.6, 0.7, 0.8, 0.9]:
    Hz_ext_analytical = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r_val)**3)
    try:
        # In Reduced region: Hz = H0 + (grad(Orr) - grad(Oxr))[1]
        Hz_pert = (grad(Orr_val) - grad(Oxr_val))[1](mesh(r_val, 0))
        Hz_total_val = H0 + Hz_pert
        err = abs(Hz_total_val - Hz_ext_analytical) / abs(Hz_ext_analytical) * 100
        print(f"  r={r_val}: Hz_total = {Hz_total_val:.6f} (analytical: {Hz_ext_analytical:.6f}), error = {err:.4f}%")
    except Exception as e:
        print(f"  r={r_val}: Error - {e}")


# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Convergence History")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error Est':<12} {'Hz Err(%)':<10}")
print("-" * 48)
for i in range(len(history['ndof'])):
    print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
          f"{history['error'][i]:<12.4e} {history['Hz_error_percent'][i]:<10.4f}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
if history['error'][-1] > 0:
    print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
          f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")
print(f"  Hz error: {history['Hz_error_percent'][0]:.4f}% -> {history['Hz_error_percent'][-1]:.4f}%")

# Save final summary to .mat file
mat_filename = os.path.join(script_dir, os.path.splitext(os.path.basename(__file__))[0] + ".mat")
mat_data = {
    'ndof': array(history['ndof']),
    'elements': array(history['elements']),
    'error': array(history['error']),
    'Hz_error_percent': array(history['Hz_error_percent']),
    'energy': array(history['energy']),
    'order': order
}
sio.savemat(mat_filename, mat_data)
print(f"\nFinal summary saved to: {mat_filename}")

# Open the final iteration plot
final_png = os.path.join(script_dir, f"axi_omega_metric_iter_{max_iterations:02d}.png")
if os.path.exists(final_png):
    try:
        os.startfile(final_png)
    except:
        pass

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
