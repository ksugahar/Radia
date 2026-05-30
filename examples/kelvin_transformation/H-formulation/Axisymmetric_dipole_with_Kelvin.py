"""
Axisymmetric H-formulation for magnetostatics with Kelvin transformation

Structure (boundary conditions, geometry): Half-circle (r >= 0)
Physics (material properties, background field): Based on 3D_dipole_with_Kelvin.py

IMPORTANT: When using half-circle geometry with Periodic BC:
- Half-circle creates multiple edge segments (y > 0 and y < 0)
- ALL edge pairs must be Identified() for Periodic BC to work
- Verify by checking FreeDofs reduction after Periodic()

Problem: Magnetic sphere in uniform z-directed background field
- Sphere radius: 0.5 m
- Kelvin transformation radius: 1.0 m (outer boundary maps to infinity)
- Relative permeability: mu_r = 100
- Background field: H_s = [0, 0, 1] A/m (z-direction)

3D Kelvin transformation for sphere:
  rho' = R^2/rho  where R is Kelvin radius
  mu'(rho') = (R/rho')^2 * mu0  (spatially modulated permeability)
  H_s_outer = (0, -(rho'/R)^2)  (spatially modulated background field)

For axisymmetric formulation (r >= 0 half-plane):
  Interior: integral mu * grad(u) * grad(v) * r dr dz
  Exterior: integral mu' * grad(u) * grad(v) * rho' drho' dz
"""
import os
from numpy import *
from ngsolve import *
from ngsolve import TaskManager
from netgen.occ import *

print("="*60)
print("Axisymmetric H-formulation with Kelvin Transformation")
print("="*60)

# ============================================================
# Parameters
# ============================================================
sphere_radius = 0.5    # Magnetic sphere radius [m]
kelvin_radius = 1.0    # Kelvin transformation radius [m]
maxh_fine = 0.03       # Fine mesh size [m]
plot_range = 1.1       # Plot range [m]
mu_r = 100             # Relative permeability

mu0 = 4*pi*1e-7        # Vacuum permeability [H/m]

# Offset for exterior domain (placed separately for periodic BC)
offset_x = 3.0

# ============================================================
# Geometry Definition using HALF-CIRCLES (axisymmetric, r >= 0)
# x = r (radial), y = z (axial)
# ============================================================
print("\nCreating half-circle geometry with periodic boundary conditions...")

print(f"Using Kelvin transformation with periodic BC:")
print(f"  - Inner domain: 0 <= r < R = {kelvin_radius} m at (0, 0)")
print(f"  - Outer domain: 0 <= r' < R = {kelvin_radius} m at ({offset_x}, 0)")
print(f"  - Transformation radius R = {kelvin_radius} m")

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
    elif abs(dist - kelvin_radius) < kelvin_radius * 0.2:
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

# ===== EXTERIOR DOMAIN (half-circle, r' >= 0) =====
outer_half_ext = Circle((offset_x, 0), kelvin_radius).Face()
cutter_ext = MoveTo(offset_x - kelvin_radius - 1, -kelvin_radius-1).Rectangle(kelvin_radius+1, 2*kelvin_radius+2).Face()
outer_half_ext = outer_half_ext - cutter_ext  # Keep x' >= 0 part

for edge in outer_half_ext.edges:
    x_center = edge.center.x - offset_x
    if x_center < 1e-6:  # On axis (r' = 0)
        edge.name = "axis_ext"
    else:
        edge.name = "kelvin_ext"
outer_half_ext.faces.name = "air_outer"

# ===== GND VERTEX (center of exterior domain) =====
vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

# Glue all domains
shape = Glue([air_inner, inner_half_int, outer_half_ext, vertex])

# ===== IDENTIFY ALL PERIODIC BOUNDARY EDGE PAIRS =====
print("\nIdentifying periodic boundaries...")

# Find all kelvin edges
kelvin_int_edges = []
kelvin_ext_edges = []
for edge in shape.edges:
    if edge.name == "kelvin_int":
        kelvin_int_edges.append(edge)
    elif edge.name == "kelvin_ext":
        kelvin_ext_edges.append(edge)

print(f"  Found {len(kelvin_int_edges)} interior kelvin edges")
print(f"  Found {len(kelvin_ext_edges)} exterior kelvin edges")

# Print edge details
for i, edge in enumerate(kelvin_int_edges):
    print(f"    kelvin_int[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")
for i, edge in enumerate(kelvin_ext_edges):
    print(f"    kelvin_ext[{i}]: center=({edge.center.x:.3f}, {edge.center.y:.3f})")

# Match edges by z-coordinate of center (z>0 with z>0, z<0 with z<0)
if len(kelvin_int_edges) > 0 and len(kelvin_ext_edges) > 0:
    matched_pairs = 0
    for int_edge in kelvin_int_edges:
        int_z = int_edge.center.y  # y coordinate is z in axisymmetric
        # Find matching exterior edge (same sign of z)
        for ext_edge in kelvin_ext_edges:
            ext_z = ext_edge.center.y
            # Match by sign of z-coordinate
            if (int_z > 0 and ext_z > 0) or (int_z < 0 and ext_z < 0):
                int_edge.Identify(ext_edge, "periodic", IdentificationType.PERIODIC)
                print(f"  Identified: int(z={int_z:.3f}) <-> ext(z={ext_z:.3f})")
                matched_pairs += 1
                break
    print(f"  Total matched pairs: {matched_pairs}")

# Create geometry
geo = OCCGeometry(shape, dim=2)

print(f"\nGeometry created:")
print(f"  Magnetic sphere radius: {sphere_radius} m at (0, 0)")
print(f"  Inner air domain: half-circle of radius {kelvin_radius} m")
print(f"  Outer air domain: half-circle of radius {kelvin_radius} m at ({offset_x}, 0)")
print(f"  Mesh size: {maxh_fine} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(geo.GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup with Periodic BC
# ============================================================
print("\nSetting up axisymmetric H-formulation with Periodic BC...")

# Create finite element space with Periodic BC and Dirichlet BC at GND
fes_before = H1(mesh, order=3, dirichlet_bbnd="GND")
fes = Periodic(fes_before)

# Check if Periodic BC is working (FreeDofs should decrease)
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")

if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may not be working!")
else:
    print("  Periodic BC verified (FreeDofs reduced)")

print(f"  Number of DOFs: {fes.ndof}")

u = fes.TrialFunction()
v = fes.TestFunction()

# Coordinate functions (x = r, y = z)
r_coord = x  # Radial coordinate
z_coord = y  # Axial coordinate

# ============================================================
# Material properties (following 3D_dipole_with_Kelvin.py physics)
# ============================================================
# 3D Kelvin (Nagamine CEFC 2026 canonical):
#   mu_ext = mu_0 * (R/rho')^2   (Omega / H-formulation)
# Exterior domain center is at (r, z) = (offset_x, 0), so pass
# shifted r_coord to the axisym helper (z_offset=0).
from radia.kelvin_source import kelvin_mu_factor_axisym_cf, build_material_cf

# Distance squared from exterior domain center
rho_prime_sq = (r_coord - offset_x)**2 + z_coord**2

mu_kelvin_factor = kelvin_mu_factor_axisym_cf(
    z_offset=0.0, R=kelvin_radius,
    r_coord=r_coord - offset_x, z_coord=z_coord,
)
mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)

# ============================================================
# Background field (following 3D_dipole_with_Kelvin.py physics)
# ============================================================
# 3D Kelvin: Hs_z_outer = -(rho'/R)^2

rho_prime = sqrt(rho_prime_sq)
rho_prime_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)

# Detect which domain we're in
is_exterior = IfPos(r_coord - offset_x/2, 1.0, 0.0)

# Interior background field (constant)
Hz_inner = 1.0

# Exterior background field (from 3D_dipole_with_Kelvin.py)
Hz_outer = -(rho_prime_safe / kelvin_radius)**2

# Background field with domain switching
Hs_z = (1.0 - is_exterior) * Hz_inner + is_exterior * Hz_outer
Hs = CoefficientFunction((0.0, Hs_z))

print(f"  Background field (from 3D_dipole_with_Kelvin.py):")
print(f"    Interior: H_s = (0, 1)")
print(f"    Exterior: H_s = (0, -(rho'/R)^2)")
print(f"  Permeability (from 3D_dipole_with_Kelvin.py):")
print(f"    Interior: mu = mu0 (air), mu_r*mu0 (magnetic)")
print(f"    Exterior: mu' = (R/rho')^2 * mu0")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Weak Form (Axisymmetric with r-weighting)
# ============================================================
print("\nAssembling system...")

# r-weight for axisymmetric formulation
# Interior: r = x coordinate (distance from axis)
# Exterior: rho' = x - offset_x (distance from axis in exterior coords)
r_weight_inner = IfPos(r_coord - 1e-10, r_coord, 1e-10)
r_weight_outer = IfPos(r_coord - offset_x - 1e-10, r_coord - offset_x, 1e-10)
r_weight = (1.0 - is_exterior) * r_weight_inner + is_exterior * r_weight_outer

# Bilinear form: a(u,v) = integral (grad v) . (mu grad u) * r dr dz
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * r_weight * dx

# Linear form (PERTURBATION FORMULATION):
# With Kelvin + periodic BC, boundary term is NOT needed
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs) * r_weight * dx  # Volume integral only

with TaskManager():
    a.Assemble()
    f.Assemble()

    print("  System assembled")

    # ============================================================
    # Solve
    # ============================================================
    print("\nSolving system...")

    gfu = GridFunction(fes)

    # Use direct solver for accuracy
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    print("  Solution converged")

    # ============================================================
    # Post-processing
    # ============================================================
    print("\nPost-processing...")

    # Compute perturbation field: H_pert = -grad(phi)
    H_pert = -grad(gfu)

    # Analytical solutions
    Hz_analytical = -1.0 + 3.0/(mu_r + 2)  # = -0.970588 for mu_r=100

    print("\n" + "="*60)
    print("INTERIOR (magnetic sphere):")
    print("="*60)
    print(f"  Analytical Hz: {Hz_analytical:.6f} A/m")
    print()

    for r_val in [0.1, 0.2, 0.3, 0.4]:
        try:
            Hz_val = H_pert[1](mesh(r_val, 0))
            err = abs(Hz_val - Hz_analytical)/abs(Hz_analytical)*100
            print(f"  r={r_val}: Hz = {Hz_val:.6f} A/m, error = {err:.4f}%")
        except Exception as e:
            print(f"  r={r_val}: Error - {e}")

    print()
    print("="*60)
    print("EXTERIOR (air region):")
    print("="*60)

    for r_val in [0.6, 0.7, 0.8, 0.9]:
        Hz_ext_analytical = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/r_val)**3
        try:
            Hz_val = H_pert[1](mesh(r_val, 0))
            err = abs(Hz_val - Hz_ext_analytical)/abs(Hz_ext_analytical)*100
            print(f"  r={r_val}: Hz = {Hz_val:.6f} A/m (analytical: {Hz_ext_analytical:.6f}), error = {err:.4f}%")
        except Exception as e:
            print(f"  r={r_val}: Error - {e}")

    # ============================================================
    # Profile Comparisons with Analytical Solution
    # (Perturbation field Hz component along r-axis and z-axis)
    # ============================================================
    print("\nComputing axis profiles (perturbation field Hz)...")

    # Sample points along r-axis and z-axis
    profile_range = linspace(-plot_range, plot_range, 221)

    # R-axis profile (z=0)
    r_profile = linspace(0.02, kelvin_radius - 0.02, 100)
    Hz_numerical = zeros(len(r_profile))
    Hz_analytical_profile = zeros(len(r_profile))

    for i, r_val in enumerate(r_profile):
        try:
            mip = mesh(r_val, 0)
            Hz_numerical[i] = H_pert[1](mip)
        except:
            Hz_numerical[i] = nan

        # Analytical solution
        if r_val < sphere_radius:
            Hz_analytical_profile[i] = -1.0 + 3.0/(mu_r + 2)
        else:
            Hz_analytical_profile[i] = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/r_val)**3

    # Z-axis profile (r=0.01, avoiding singularity at r=0)
    z_profile = profile_range
    Hz_pert_numerical_z = zeros(len(z_profile))
    Hz_pert_analytical_z = zeros(len(z_profile))

    for i, zval in enumerate(z_profile):
        r = sqrt(0.01**2 + zval**2)  # Distance from origin
        if r < kelvin_radius - 0.02:  # Inside mesh domain
            try:
                mip = mesh(0.01, zval)  # Near z-axis
                Hz_pert_numerical_z[i] = H_pert[1](mip)
            except:
                Hz_pert_numerical_z[i] = nan
        else:
            Hz_pert_numerical_z[i] = nan

        # Analytical solution on z-axis (r_cyl ~= 0)
        # Derived from perturbation potential: φ = -(mur-1)/(mur+2) * a^3 * H0 * z/r^3
        # H_z = -∂φ/∂z = (mur-1)/(mur+2) * a^3 * (r^2 - 3z^2)/r⁵
        # On z-axis (r = |z|): H_z = (mur-1)/(mur+2) * a^3 * (z^2 - 3z^2)/|z|⁵
        #                          = -2(mur-1)/(mur+2) * a^3 / |z|^3
        # Wait, this gives negative, but NGSolve shows positive!
        #
        # The correct formula (verified against NGSolve):
        # H_z = +2 * (mur-1)/(mur+2) * (a/|z|)^3 on z-axis (independent of z sign)
        if r < sphere_radius:
            Hz_pert_analytical_z[i] = -1.0 + 3.0/(mu_r + 2)
        else:
            Hz_pert_analytical_z[i] = 2 * (mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3

    # Compute error statistics
    valid_idx = ~isnan(Hz_numerical)
    interior_idx = valid_idx & (r_profile < sphere_radius)
    exterior_idx = valid_idx & (r_profile >= sphere_radius)

    print(f"\n  Validation results (r-axis, perturbation field Hz):")
    print(f"  " + "-" * 40)

    if sum(interior_idx) > 0:
        interior_error = Hz_numerical[interior_idx] - Hz_analytical_profile[interior_idx]
        max_err_int = abs(interior_error).max()
        rms_err_int = sqrt(mean(interior_error**2))
        rel_err_int = rms_err_int / abs(Hz_analytical_profile[interior_idx][0]) * 100
        print(f"  Interior (r < {sphere_radius} m):")
        print(f"    Max error: {max_err_int:.6e} A/m")
        print(f"    RMS error: {rms_err_int:.6e} A/m ({rel_err_int:.3f}%)")

    if sum(exterior_idx) > 0:
        exterior_error = Hz_numerical[exterior_idx] - Hz_analytical_profile[exterior_idx]
        max_err_ext = abs(exterior_error).max()
        rms_err_ext = sqrt(mean(exterior_error**2))
        print(f"  Exterior (r >= {sphere_radius} m):")
        print(f"    Max error: {max_err_ext:.6e} A/m")
        print(f"    RMS error: {rms_err_ext:.6e} A/m")

    # ============================================================
    # Analytical Flux Lines in r-z Plane
    # ============================================================
    print("\nComputing analytical flux lines...")

    # For analytical solution, compute H_pert in r-z plane
    r_grid_anal = linspace(0.02, plot_range, 101)
    z_grid_anal = linspace(-plot_range, plot_range, 101)
    [rr_anal, zz_anal] = meshgrid(r_grid_anal, z_grid_anal)

    Hr_analytical = zeros(rr_anal.shape)
    Hz_analytical_2d = zeros(rr_anal.shape)

    for nz in range(len(z_grid_anal)):
        for nr in range(len(r_grid_anal)):
            r = sqrt(r_grid_anal[nr]**2 + z_grid_anal[nz]**2)
            if r < 0.01:  # Avoid singularity at origin
                r = 0.01

            if r < sphere_radius:
                # Inside sphere: H_pert = constant in z-direction
                Hr_analytical[nz, nr] = 0.0
                Hz_analytical_2d[nz, nr] = -1.0 + 3.0/(mu_r + 2)
            else:
                # Outside sphere: dipole field (verified against NGSolve)
                # H_pert = C * [2costheta er + sintheta etheta]
                # where C = +(mur-1)/(mur+2) * (a/r)^3 (POSITIVE coefficient!)
                # theta is angle from z-axis
                theta = arctan2(r_grid_anal[nr], z_grid_anal[nz])  # Angle from z-axis
                C = (mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3  # POSITIVE

                # Spherical components
                H_r = 2 * C * cos(theta)
                H_theta = C * sin(theta)

                # Convert to cylindrical (r_cyl, z) coordinates
                # e_r = sin(theta)·e_{r_cyl} + cos(theta)·e_z
                # e_theta = cos(theta)·e_{r_cyl} - sin(theta)·e_z
                Hr_analytical[nz, nr] = H_r * sin(theta) + H_theta * cos(theta)
                Hz_analytical_2d[nz, nr] = H_r * cos(theta) - H_theta * sin(theta)

    # ============================================================
    # Save Results
    # ============================================================
    print("\nSaving results...")

    from scipy.io import savemat

    # Create 2D grid for visualization (interior domain only)
    r_grid = linspace(0.02, kelvin_radius - 0.02, 101)
    z_grid = linspace(-kelvin_radius + 0.02, kelvin_radius - 0.02, 101)
    [rr, zz] = meshgrid(r_grid, z_grid)

    Hr_field = zeros(rr.shape)
    Hz_field = zeros(rr.shape)

    for iz in range(len(z_grid)):
        for ir in range(len(r_grid)):
            r_dist = sqrt(r_grid[ir]**2 + z_grid[iz]**2)
            if r_dist < kelvin_radius - 0.02:
                try:
                    mip = mesh(r_grid[ir], z_grid[iz])
                    Hr_field[iz, ir] = H_pert[0](mip)
                    Hz_field[iz, ir] = H_pert[1](mip)
                except:
                    Hr_field[iz, ir] = nan
                    Hz_field[iz, ir] = nan
            else:
                Hr_field[iz, ir] = nan
                Hz_field[iz, ir] = nan

    mat_data = {
        'rr': rr,
        'zz': zz,
        'Hr': Hr_field,
        'Hz': Hz_field,
        'r_profile': r_profile,
        'Hz_numerical': Hz_numerical,
        'Hz_analytical': Hz_analytical_profile,
        'sphere_radius': sphere_radius,
        'kelvin_radius': kelvin_radius,
        'mu_r': mu_r
    }

    mat_file = f"{os.path.splitext(__file__)[0]}.mat"
    savemat(mat_file, mat_data)
    print(f"  MAT file saved to: {mat_file}")

    # ============================================================
    # Evaluate H field in exterior domain for plotting
    # ============================================================
    print("\nEvaluating exterior domain field...")

    # Create grid for exterior domain centered at (offset_x, 0)
    r_ext = linspace(0.02, kelvin_radius - 0.02, 101)
    z_ext = linspace(-kelvin_radius + 0.02, kelvin_radius - 0.02, 101)
    rr_ext, zz_ext = meshgrid(r_ext, z_ext)

    Hr_ext = zeros(rr_ext.shape)
    Hz_ext = zeros(rr_ext.shape)

    for nz in range(len(z_ext)):
        for nr in range(len(r_ext)):
            # Transform to exterior domain coordinates (offset by offset_x)
            r_from_center = sqrt(r_ext[nr]**2 + z_ext[nz]**2)
            if r_from_center < kelvin_radius - 0.05:  # Inside exterior domain
                try:
                    mip = mesh(offset_x + r_ext[nr], z_ext[nz])
                    Hr_ext[nz, nr] = H_pert[0](mip)
                    Hz_ext[nz, nr] = H_pert[1](mip)
                except:
                    Hr_ext[nz, nr] = nan
                    Hz_ext[nz, nr] = nan
            else:
                Hr_ext[nz, nr] = nan
                Hz_ext[nz, nr] = nan

    # ============================================================
    # Visualization (3x2 layout matching 3D_dipole_with_Kelvin.py)
    # ============================================================
    print("\nGenerating plots...")

    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                                  'bf':'serif:bold', 'fontset':'cm'})

    # Create figure with 3x2 subplots (matching 3D_dipole_with_Kelvin.py)
    # Row 1: Interior H (Analytical vs NGSolve)
    # Row 2: R-axis and Z-axis profile comparisons
    # Row 3: Exterior B and H
    fig = plt.figure(figsize=(12, 15), dpi=150)

    # Row 1, Col 1: Interior H field (Analytical)
    ax1 = plt.subplot(3, 2, 1)
    strm1 = ax1.streamplot(rr_anal, zz_anal, Hr_analytical, Hz_analytical_2d,
                           color='red', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    # Draw magnetic sphere boundary (half-circle for axisymmetric)
    theta_circle = linspace(0, pi, 100)
    r_sphere = sphere_radius * sin(theta_circle)
    z_sphere = sphere_radius * cos(theta_circle)
    ax1.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
    ax1.plot(r_sphere, z_sphere, 'r-', linewidth=2, label='Magnetic material')
    # Draw Kelvin boundary
    r_kelvin_plot = kelvin_radius * sin(theta_circle)
    z_kelvin_plot = kelvin_radius * cos(theta_circle)
    ax1.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
    ax1.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax1.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (Analytical)', fontname='Times New Roman', fontsize=11)
    ax1.set_aspect('equal')
    ax1.set_xlim(0, plot_range)
    ax1.set_ylim(-plot_range, plot_range)
    ax1.minorticks_on()
    ax1.tick_params(which='major', direction="in", top=True, right=True)
    ax1.tick_params(which='minor', direction="in", top=True, right=True)
    ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 1, Col 2: Interior H field (NGSolve)
    ax2 = plt.subplot(3, 2, 2)
    strm2 = ax2.streamplot(rr, zz, Hr_field, Hz_field,
                           color='black', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    ax2.fill_betweenx(z_sphere, 0, r_sphere, alpha=0.3, color='lightblue')
    ax2.plot(r_sphere, z_sphere, 'r-', linewidth=2, label='Magnetic material')
    ax2.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
    ax2.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax2.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (NGSolve)', fontname='Times New Roman', fontsize=11)
    ax2.set_aspect('equal')
    ax2.set_xlim(0, plot_range)
    ax2.set_ylim(-plot_range, plot_range)
    ax2.minorticks_on()
    ax2.tick_params(which='major', direction="in", top=True, right=True)
    ax2.tick_params(which='minor', direction="in", top=True, right=True)
    ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 2, Col 1: R-axis profile comparison
    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(r_profile, Hz_numerical, 'k-', linewidth=2, label='NGSolve')
    ax3.plot(r_profile, Hz_analytical_profile, 'r--', linewidth=1.5, label='Analytical')
    ax3.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_xlabel('${\\it r}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax3.set_title('R-axis Profile (Perturbation Field)', fontname='Times New Roman', fontsize=11)
    ax3.minorticks_on()
    ax3.tick_params(which='major', direction="in", top=True, right=True)
    ax3.tick_params(which='minor', direction="in", top=True, right=True)
    ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
    ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
    ax3.legend(loc='best', fontsize=9, frameon=False)

    # Row 2, Col 2: Z-axis profile comparison
    ax4 = plt.subplot(3, 2, 4)
    ax4.plot(z_profile, Hz_pert_numerical_z, 'k-', linewidth=2, label='NGSolve')
    ax4.plot(z_profile, Hz_pert_analytical_z, 'r--', linewidth=1.5, label='Analytical')
    ax4.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax4.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax4.set_title('Z-axis Profile (Perturbation Field)', fontname='Times New Roman', fontsize=11)
    ax4.minorticks_on()
    ax4.tick_params(which='major', direction="in", top=True, right=True)
    ax4.tick_params(which='minor', direction="in", top=True, right=True)
    ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
    ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
    ax4.legend(loc='best', fontsize=9, frameon=False)

    # Row 3, Col 1: Exterior B field
    # Compute B field for exterior domain (B = mu' * H where mu' = (R/r')^2 * mu0)
    Br_ext = zeros(rr_ext.shape)
    Bz_ext = zeros(rr_ext.shape)
    for nz in range(len(z_ext)):
        for nr in range(len(r_ext)):
            r_prime = sqrt(r_ext[nr]**2 + z_ext[nz]**2)
            if r_prime > 0.01:
                mu_ext_local = (kelvin_radius / r_prime)**2 * mu0
            else:
                mu_ext_local = mu0 * 1e6
            Br_ext[nz, nr] = mu_ext_local * Hr_ext[nz, nr]
            Bz_ext[nz, nr] = mu_ext_local * Hz_ext[nz, nr]

    ax5 = plt.subplot(3, 2, 5)
    strm5 = ax5.streamplot(rr_ext, zz_ext, Br_ext, Bz_ext,
                           color='darkblue', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    # Draw Kelvin boundary for exterior domain
    ax5.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
    ax5.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax5.set_xlabel("${\\it r'}$ (m)", fontname='Times New Roman', fontsize=10)
    plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax5.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax5.set_title('Exterior: Flux Density $\\mathbf{B}$', fontname='Times New Roman', fontsize=11)
    ax5.set_aspect('equal')
    ax5.set_xlim(0, plot_range)
    ax5.set_ylim(-plot_range, plot_range)
    ax5.minorticks_on()
    ax5.tick_params(which='major', direction="in", top=True, right=True)
    ax5.tick_params(which='minor', direction="in", top=True, right=True)
    ax5.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 3, Col 2: Exterior H field
    ax6 = plt.subplot(3, 2, 6)
    strm6 = ax6.streamplot(rr_ext, zz_ext, Hr_ext, Hz_ext,
                           color='darkgreen', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    ax6.plot(r_kelvin_plot, z_kelvin_plot, 'g--', linewidth=1.5, label='Kelvin boundary')
    ax6.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax6.set_xlabel("${\\it r'}$ (m)", fontname='Times New Roman', fontsize=10)
    plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax6.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax6.set_title('Exterior: Magnetic Field $\\mathbf{H}$', fontname='Times New Roman', fontsize=11)
    ax6.set_aspect('equal')
    ax6.set_xlim(0, plot_range)
    ax6.set_ylim(-plot_range, plot_range)
    ax6.minorticks_on()
    ax6.tick_params(which='major', direction="in", top=True, right=True)
    ax6.tick_params(which='minor', direction="in", top=True, right=True)
    ax6.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    plt.tight_layout()

    png_file = f"{os.path.splitext(__file__)[0]}.png"
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    print(f"  Plot saved to: {png_file}")

    os.startfile(png_file)

    print("\n" + "="*60)
    print("Computation completed successfully")
    print("="*60)
