"""
H-formulation for magnetostatics with perturbation potential
Geometry created internally using OCC
Updated: 2025-11-22
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve
from ngsolve import TaskManager

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 3D - OCC Geometry")
print("="*60)

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry...")

# Parameters
sphere_radius = 0.5  # Magnetic sphere radius [m]
air_inner_radius = 1.0  # Outer boundary radius [m]
maxh_fine = 0.03     # Fine mesh size [m]
plot_range = 1.1    # Plot range [m]

# Create magnetic sphere (finest mesh)
mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere.mat("magnetic")
mag_sphere.maxh = maxh_fine

# Create air sphere and name its boundary
air_sphere = Sphere(Pnt(0, 0, 0), air_inner_radius)
for face in air_sphere.faces:
    face.name = "outer"  # Name outer boundary before boolean operation
air_sphere.maxh = maxh_fine

# Boolean operations to create two regions
# Air = air sphere - magnetic sphere
air = air_sphere - mag_sphere
air.mat("air")

# Combine into single geometry
geo = Glue([air, mag_sphere])

print(f"Geometry created with two regions:")
print(f"  Magnetic sphere radius: {sphere_radius} m")
print(f"  Outer boundary radius: {air_inner_radius} m")
print(f"  Fine mesh size: {maxh_fine} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup
# ============================================================
print("\nSetting up H-formulation...")

n = specialcf.normal(mesh.dim)
fes = H1(mesh, order=3)
print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability
mu_d = {"air": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field: H_s = [0, 0, 1] A/m (z-direction)
Hs = CoefficientFunction((0, 0, 1))
Hsb = BoundaryFromVolumeCF(Hs)

print(f"  Background field: H_s = [0, 0, 1] A/m (z-direction)")
print(f"  Relative permeability: mu_r = {mu_r}")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = ∫(∇v)·(mu∇u)dOmega
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (PERTURBATION FORMULATION):
# f(v) = ∫(∇v)·(muH_s)dOmega - ∫v(n·muH_s)dΓ
# 注意: Kelvin変換なし（有限領域）では外部境界での境界項が必要
#       外部境界を通る磁束を考慮するため
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                    # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer"))  # 境界項（必須）

with TaskManager():
    a.Assemble()
    f.Assemble()

    print("  System assembled")

    # ============================================================
    # Solve
    # ============================================================
    print("\nSolving system...")

    gfu = GridFunction(fes)
    c = Preconditioner(a, type="local")

    solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat,
               tol=1e-5, printrates=True, maxsteps=10000)

    print("  Solution converged")

    # ============================================================
    # Post-processing
    # ============================================================
    print("\nPost-processing...")

    # Create evaluation grids for x-y and x-z planes
    x = linspace(-plot_range, plot_range, 221)
    y = linspace(-plot_range, plot_range, 221)
    z = linspace(-plot_range, plot_range, 221)
    [xx, yy] = meshgrid(x, y)
    [xx_xz, zz_xz] = meshgrid(x, z)

    # Compute perturbation field: H_pert = -grad(phi)
    H = -grad(gfu)

    # X-Y plane (z=0) - for potential distribution
    Hx = zeros((shape(xx)))
    Hy = zeros((shape(xx)))
    phi = zeros((shape(xx)))

    for ny in range(len(y)):
        for nx in range(len(x)):
            r = sqrt(x[nx]**2 + y[ny]**2)
            if r < air_inner_radius - 0.01:  # Inside mesh domain
                try:
                    mip = mesh(x[nx], y[ny], 0)
                    phi[ny, nx] = gfu(mip)
                    Hx[ny, nx] = H[0](mip)
                    Hy[ny, nx] = H[1](mip)
                except:
                    phi[ny, nx] = nan
                    Hx[ny, nx] = nan
                    Hy[ny, nx] = nan
            else:
                phi[ny, nx] = nan
                Hx[ny, nx] = nan
                Hy[ny, nx] = nan

    # Center potential at origin
    center_value = gfu(mesh(0, 0, 0))
    phi = phi - center_value

    # X-Z plane (y=0) - for flux lines
    Hx_xz = zeros((shape(xx_xz)))
    Hz_xz = zeros((shape(xx_xz)))

    for nz in range(len(z)):
        for nx in range(len(x)):
            r = sqrt(x[nx]**2 + z[nz]**2)
            if r < air_inner_radius - 0.01:  # Inside mesh domain
                try:
                    mip = mesh(x[nx], 0, z[nz])
                    Hx_xz[nz, nx] = H[0](mip)
                    Hz_xz[nz, nx] = H[2](mip)
                except:
                    Hx_xz[nz, nx] = nan
                    Hz_xz[nz, nx] = nan
            else:
                Hx_xz[nz, nx] = nan
                Hz_xz[nz, nx] = nan

    print(f"  Potential at origin: {center_value:.6e}")
    print(f"  Field at origin: Hz = {H[2](mesh(0,0,0)):.6f} A/m")

    # Expected analytical value (perturbation field interior, z-component)
    Hz_analytical = -1.0 + 3.0/(mu_r + 2)  # = -0.75 for mur=10
    print(f"  Analytical (interior): Hz = {Hz_analytical:.6f} A/m")
    print(f"  Relative error: {abs(H[2](mesh(0,0,0)) - Hz_analytical)/abs(Hz_analytical)*100:.3f}%")

    # ============================================================
    # Profile Comparisons with Analytical Solution
    # (Perturbation field Hz component along x-axis and y-axis)
    # Background field: H_s = [0, 0, 1] (z-direction)
    # ============================================================
    print("\nComputing axis profiles (perturbation field Hz)...")

    # Sample points along x-axis and y-axis
    profile_range = linspace(-plot_range, plot_range, 221)

    # X-axis profile
    x_profile = profile_range
    Hz_pert_numerical_x = zeros(len(x_profile))
    Hz_pert_analytical_x = zeros(len(x_profile))

    for i, xval in enumerate(x_profile):
        r = abs(xval)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(xval, 0, 0)  # X-axis: (x, 0, 0)
                Hz_pert_numerical_x[i] = H[2](mip)  # Hz component
            except:
                Hz_pert_numerical_x[i] = nan
        else:
            Hz_pert_numerical_x[i] = nan

        # Analytical solution for PERTURBATION field Hz component
        if r < sphere_radius:
            Hz_pert_analytical_x[i] = -1.0 + 3.0/(mu_r + 2)  # = -0.75 for mur=10
        else:
            Hz_pert_analytical_x[i] = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3

    # Y-axis profile
    y_profile = profile_range
    Hz_pert_numerical_y = zeros(len(y_profile))
    Hz_pert_analytical_y = zeros(len(y_profile))

    for i, yval in enumerate(y_profile):
        r = abs(yval)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(0, yval, 0)  # Y-axis: (0, y, 0)
                Hz_pert_numerical_y[i] = H[2](mip)  # Hz component
            except:
                Hz_pert_numerical_y[i] = nan
        else:
            Hz_pert_numerical_y[i] = nan

        # Analytical solution
        if r < sphere_radius:
            Hz_pert_analytical_y[i] = -1.0 + 3.0/(mu_r + 2)
        else:
            Hz_pert_analytical_y[i] = -(mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3

    # Compute error statistics for x-axis
    valid_idx_x = ~isnan(Hz_pert_numerical_x)
    interior_idx_x = valid_idx_x & (abs(x_profile) < sphere_radius)
    exterior_idx_x = valid_idx_x & (abs(x_profile) >= sphere_radius)

    print(f"\n  Validation results (X-axis, perturbation field Hz):")
    print(f"  -" * 30)

    if sum(interior_idx_x) > 0:
        interior_error = Hz_pert_numerical_x[interior_idx_x] - Hz_pert_analytical_x[interior_idx_x]
        max_err_int = max(abs(interior_error))
        rms_err_int = sqrt(mean(interior_error**2))
        rel_err_int = rms_err_int / abs(Hz_pert_analytical_x[interior_idx_x][0]) * 100
        print(f"  Interior (|x| < {sphere_radius} m):")
        print(f"    Max error: {max_err_int:.6e} A/m")
        print(f"    RMS error: {rms_err_int:.6e} A/m ({rel_err_int:.3f}%)")

    if sum(exterior_idx_x) > 0:
        exterior_error = Hz_pert_numerical_x[exterior_idx_x] - Hz_pert_analytical_x[exterior_idx_x]
        max_err_ext = max(abs(exterior_error))
        rms_err_ext = sqrt(mean(exterior_error**2))
        print(f"  Exterior (|x| >= {sphere_radius} m):")
        print(f"    Max error: {max_err_ext:.6e} A/m")
        print(f"    RMS error: {rms_err_ext:.6e} A/m")

    # ============================================================
    # Analytical Flux Lines in X-Z Plane
    # ============================================================
    print("\nComputing analytical flux lines...")

    # For analytical solution, compute H_pert in x-z plane
    Hx_xz_analytical = zeros((shape(xx_xz)))
    Hz_xz_analytical = zeros((shape(xx_xz)))

    for nz in range(len(z)):
        for nx in range(len(x)):
            r = sqrt(x[nx]**2 + z[nz]**2)
            if r < 0.01:  # Avoid singularity at origin
                r = 0.01

            if r < sphere_radius:
                # Inside sphere: H_pert = constant in z-direction
                Hx_xz_analytical[nz, nx] = 0.0
                Hz_xz_analytical[nz, nx] = -1.0 + 3.0/(mu_r + 2)
            else:
                # Outside sphere: dipole field
                # H_pert = C * [2costheta er + sintheta etheta]
                # where C = +(mur-1)/(mur+2) * (a/r)^3 (POSITIVE coefficient)
                # In Cartesian: Hx = H_r * sintheta + H_theta * costheta
                #               Hz = H_r * costheta - H_theta * sintheta
                theta = arctan2(x[nx], z[nz])  # Angle from z-axis
                C = (mu_r - 1)/(mu_r + 2) * (sphere_radius/r)**3

                # Radial and tangential components
                H_r = 2 * C * cos(theta)
                H_theta = C * sin(theta)

                # Convert to Cartesian
                Hx_xz_analytical[nz, nx] = H_r * sin(theta) + H_theta * cos(theta)
                Hz_xz_analytical[nz, nx] = H_r * cos(theta) - H_theta * sin(theta)

    # ============================================================
    # Save Results
    # ============================================================
    print("\nSaving results...")

    # Save to .mat file
    from scipy.io import savemat
    mat_file = f"{os.path.splitext(__file__)[0]}.mat"
    savemat(mat_file, {
        'xx': xx_xz,
        'zz': zz_xz,
        'Hx_analytical': Hx_xz_analytical,
        'Hz_analytical': Hz_xz_analytical,
        'Hx': Hx_xz,
        'Hz': Hz_xz
    })
    print(f"  MAT file saved to: {mat_file}")

    # ============================================================
    # Visualization
    # ============================================================
    print("\nGenerating plots...")

    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                                  'bf':'serif:bold', 'fontset':'cm'})

    # Create figure with 2x2 subplots (no exterior domain without Kelvin)
    # Row 1: Analytical H vs NGSolve H (streamlines)
    # Row 2: X-axis and Y-axis profile comparisons
    fig = plt.figure(figsize=(12, 10), dpi=150)

    # Row 1, Col 1: Interior H field (Analytical)
    ax1 = plt.subplot(2, 2, 1)
    strm1 = ax1.streamplot(xx_xz, zz_xz, Hx_xz_analytical, Hz_xz_analytical,
                           color='red', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    circle1 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                         alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
    ax1.add_patch(circle1)
    ax1.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax1.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (Analytical)', fontname='Times New Roman', fontsize=11)
    ax1.set_aspect('equal')
    ax1.set_xlim(-plot_range, plot_range)
    ax1.set_ylim(-plot_range, plot_range)
    ax1.minorticks_on()
    ax1.tick_params(which='major', direction="in", top=True, right=True)
    ax1.tick_params(which='minor', direction="in", top=True, right=True)
    ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 1, Col 2: Interior H field (NGSolve)
    ax2 = plt.subplot(2, 2, 2)
    strm2 = ax2.streamplot(xx_xz, zz_xz, Hx_xz, Hz_xz,
                           color='black', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    circle2 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                         alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
    ax2.add_patch(circle2)
    ax2.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    ax2.set_title('Interior: $\\mathbf{H}_{\\mathrm{pert}}$ (NGSolve)', fontname='Times New Roman', fontsize=11)
    ax2.set_aspect('equal')
    ax2.set_xlim(-plot_range, plot_range)
    ax2.set_ylim(-plot_range, plot_range)
    ax2.minorticks_on()
    ax2.tick_params(which='major', direction="in", top=True, right=True)
    ax2.tick_params(which='minor', direction="in", top=True, right=True)
    ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 2, Col 1: X-axis profile comparison
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(x_profile, Hz_pert_numerical_x, 'k-', linewidth=2, label='NGSolve')
    ax3.plot(x_profile, Hz_pert_analytical_x, 'r--', linewidth=1.5, label='Analytical')
    ax3.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax3.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax3.set_title('X-axis Profile ($H_z$ component)', fontname='Times New Roman', fontsize=11)
    ax3.minorticks_on()
    ax3.tick_params(which='major', direction="in", top=True, right=True)
    ax3.tick_params(which='minor', direction="in", top=True, right=True)
    ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
    ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
    ax3.legend(loc='best', fontsize=9, frameon=False)

    # Row 2, Col 2: Y-axis profile comparison
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(y_profile, Hz_pert_numerical_y, 'k-', linewidth=2, label='NGSolve')
    ax4.plot(y_profile, Hz_pert_analytical_y, 'r--', linewidth=1.5, label='Analytical')
    ax4.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax4.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_xlabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_ylabel('$H_{z,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax4.set_title('Y-axis Profile ($H_z$ component)', fontname='Times New Roman', fontsize=11)
    ax4.minorticks_on()
    ax4.tick_params(which='major', direction="in", top=True, right=True)
    ax4.tick_params(which='minor', direction="in", top=True, right=True)
    ax4.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
    ax4.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
    ax4.legend(loc='best', fontsize=9, frameon=False)

    plt.tight_layout()

    png_file = f"{os.path.splitext(__file__)[0]}.png"
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    print(f"  Plot saved to: {png_file}")

    os.startfile(png_file)

    print("\n" + "="*60)
    print("Computation completed successfully")
    print("="*60)
