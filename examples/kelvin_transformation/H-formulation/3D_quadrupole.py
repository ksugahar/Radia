"""
H-formulation for 3D magnetostatics with perturbation potential
QUADRUPOLE background field: H_s = (-z, 0, -x) in X-Z plane
Potential: φ_s = xz
Geometry created internally using OCC
Updated: 2025-11-27
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve
from ngsolve import TaskManager

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 3D QUADRUPOLE - OCC Geometry")
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

# Background field: H_s = (-z, 0, -x) A/m (QUADRUPOLE field in X-Z plane)
# This satisfies div(H_s) = 0 and rot(H_s) = 0
# Potential: φ_s = xz, so H_s = -∇φ_s = (-z, 0, -x)
Hs = CoefficientFunction((-z, 0, -x))
Hsb = BoundaryFromVolumeCF(Hs)

print(f"  Background field: H_s = (-z, 0, -x) A/m (quadrupole in X-Z plane)")
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

    # ============================================================
    # Analytical Solution for 3D Quadrupole
    # ============================================================
    # For H_s = (-z, 0, -x), the potential is φ_s = xz (quadrupole in X-Z plane)
    #
    # Interior solution (r < a):
    #   φ_pert = B * xz
    #   where B = -2(mur-1)/(2mur+5) (from boundary conditions)
    #   H_pert = -∇φ_pert = (-B*z, 0, -B*x)
    #
    # Exterior solution (r > a):
    #   φ_pert = A * (a⁵/r⁵) * xz
    #   where A = -2(mur-1)/(2mur+5)
    #   H_pert = -∇(A*a⁵*xz/r⁵)

    # Coefficient for interior perturbation potential
    B_coeff = -2.0 * (mu_r - 1.0) / (2.0 * mu_r + 5.0)

    # Test at (0.3, 0, 0.3) where both x and z are non-zero (inside sphere)
    x_test2, z_test2 = 0.3, 0.3
    try:
        Hz_test2 = H[2](mesh(x_test2, 0, z_test2))
        Hx_test2 = H[0](mesh(x_test2, 0, z_test2))
        print(f"  Field at ({x_test2}, 0, {z_test2}): Hx = {Hx_test2:.6f}, Hz = {Hz_test2:.6f} A/m")

        # Interior analytical: H_pert,z = -∂φ_pert/∂z = -B*x, H_pert,x = -B*z
        Hx_analytical_test2 = -B_coeff * z_test2
        Hz_analytical_test2 = -B_coeff * x_test2
        print(f"  Analytical (interior): Hx = {Hx_analytical_test2:.6f}, Hz = {Hz_analytical_test2:.6f} A/m")
    except:
        print(f"  Could not evaluate at ({x_test2}, 0, {z_test2})")

    # ============================================================
    # Profile Comparisons with Analytical Solution
    # (Hz along x-axis and Hx along z-axis)
    # Background field: H_s = (-z, 0, -x) (quadrupole field in X-Z plane)
    # ============================================================
    print("\nComputing axis profiles...")

    # Sample points along x-axis and z-axis
    profile_range = linspace(-plot_range, plot_range, 221)

    # X-axis profile: (x, 0, 0) - evaluate Hz component
    # On x-axis (z=0): φ_pert = B*x*0 = 0, but Hz = -∂φ/∂z = -B*x != 0
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

        # Analytical solution for PERTURBATION field Hz component on x-axis
        # Interior (r < a): φ_pert = B·xz, Hz_pert = -∂φ_pert/∂z = -B·x
        # Exterior (r > a): φ_pert = A·xz/r^5, Hz_pert = -A·x/r^5 (at z=0)
        # Coefficients: B = -2(mu_r-1)/(2mu_r+5), A = B x a^5
        r = abs(xval)
        if r < sphere_radius:
            # Inside sphere
            B = -2*(mu_r - 1)/(2*mu_r + 5)
            Hz_pert_analytical_x[i] = -B * xval
        elif r > 1e-10:
            # Outside sphere
            A = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)
            Hz_pert_analytical_x[i] = -A * xval / r**5
        else:
            Hz_pert_analytical_x[i] = 0.0

    # Z-axis profile: (0, 0, z) - evaluate Hx component
    # On z-axis (x=0): φ_pert = B*0*z = 0, but Hx = -∂φ/∂x = -B*z != 0
    z_profile = profile_range
    Hx_pert_numerical_z = zeros(len(z_profile))
    Hx_pert_analytical_z = zeros(len(z_profile))

    for i, zval in enumerate(z_profile):
        r = abs(zval)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(0, 0, zval)  # Z-axis: (0, 0, z)
                Hx_pert_numerical_z[i] = H[0](mip)  # Hx component
            except:
                Hx_pert_numerical_z[i] = nan
        else:
            Hx_pert_numerical_z[i] = nan

        # Analytical solution for PERTURBATION field Hx component on z-axis
        # Interior (r < a): φ_pert = B·xz, Hx_pert = -∂φ_pert/∂x = -B·z
        # Exterior (r > a): φ_pert = A·xz/r^5, Hx_pert = -A·z/r^5 (at x=0)
        # Coefficients: B = -2(mu_r-1)/(2mu_r+5), A = B x a^5
        r = abs(zval)
        if r < sphere_radius:
            # Inside sphere
            B = -2*(mu_r - 1)/(2*mu_r + 5)
            Hx_pert_analytical_z[i] = -B * zval
        elif r > 1e-10:
            # Outside sphere
            A = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)
            Hx_pert_analytical_z[i] = -A * zval / r**5
        else:
            Hx_pert_analytical_z[i] = 0.0

    # Error statistics for x-axis
    valid_idx_x = ~isnan(Hz_pert_numerical_x)
    interior_idx_x = valid_idx_x & (abs(x_profile) < sphere_radius)

    print(f"\n  Validation results (X-axis, perturbation field Hz):")
    print(f"  -" * 30)

    if sum(interior_idx_x) > 0:
        interior_error = Hz_pert_numerical_x[interior_idx_x] - Hz_pert_analytical_x[interior_idx_x]
        max_err_int = max(abs(interior_error))
        rms_err_int = sqrt(mean(interior_error**2))
        print(f"  Interior (|x| < {sphere_radius} m):")
        print(f"    Max error: {max_err_int:.6e} A/m")
        print(f"    RMS error: {rms_err_int:.6e} A/m")
        print(f"    (Expected ~0 on x-axis due to symmetry)")

    # ============================================================
    # Analytical Flux Lines in X-Z Plane
    # ============================================================
    print("\nComputing analytical flux lines...")

    # For analytical solution, compute H_pert in x-z plane
    # φ_pert = B * xz (interior), φ_pert = A * (a⁵/r⁵) * xz (exterior)
    # H_pert = -grad(φ_pert)
    Hx_xz_analytical = zeros((shape(xx_xz)))
    Hz_xz_analytical = zeros((shape(xx_xz)))

    for nz in range(len(z)):
        for nx in range(len(x)):
            xval = x[nx]
            zval = z[nz]
            r = sqrt(xval**2 + zval**2)
            if r < 0.01:  # Avoid singularity at origin
                r = 0.01

            if r < sphere_radius:
                # Inside sphere: φ_pert = B * xz
                # H_pert,x = -∂φ/∂x = -B*z
                # H_pert,z = -∂φ/∂z = -B*x
                Hx_xz_analytical[nz, nx] = -B_coeff * zval
                Hz_xz_analytical[nz, nx] = -B_coeff * xval
            else:
                # Outside sphere: φ_pert = A * (a⁵/r⁵) * xz where A = B
                # Let f = a⁵/r⁵ = a⁵ * (x^2 + z^2)^(-5/2)
                # ∂f/∂x = -5 * a⁵ * x / r⁷
                # ∂(f*xz)/∂x = z*f + xz*∂f/∂x = a⁵*z/r⁷ * (r^2 - 5x^2)
                # ∂(f*xz)/∂z = a⁵*x/r⁷ * (r^2 - 5z^2)
                A = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)
                Hx_xz_analytical[nz, nx] = -A * zval / r**7 * (r**2 - 5*xval**2)
                Hz_xz_analytical[nz, nx] = -A * xval / r**7 * (r**2 - 5*zval**2)

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
    # Row 1: Analytical H streamline vs NGSolve H streamline
    # Row 2: X-axis and Z-axis profile comparisons
    fig = plt.figure(figsize=(12, 10), dpi=150)

    # Row 1, Col 1: Analytical H field streamline
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
    ax1.set_title('Analytical $\\mathbf{H}_{\\mathrm{pert}}$ Streamline (X-Z plane)', fontname='Times New Roman', fontsize=11)
    ax1.set_aspect('equal')
    ax1.set_xlim(-plot_range, plot_range)
    ax1.set_ylim(-plot_range, plot_range)
    ax1.minorticks_on()
    ax1.tick_params(which='major', direction="in", top=True, right=True)
    ax1.tick_params(which='minor', direction="in", top=True, right=True)
    ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 1, Col 2: NGSolve H field streamline
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
    ax2.set_title('NGSolve $\\mathbf{H}_{\\mathrm{pert}}$ Streamline (X-Z plane)', fontname='Times New Roman', fontsize=11)
    ax2.set_aspect('equal')
    ax2.set_xlim(-plot_range, plot_range)
    ax2.set_ylim(-plot_range, plot_range)
    ax2.minorticks_on()
    ax2.tick_params(which='major', direction="in", top=True, right=True)
    ax2.tick_params(which='minor', direction="in", top=True, right=True)
    ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 2, Col 1: X-axis profile comparison (Hz)
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

    # Row 2, Col 2: Z-axis profile comparison (Hx)
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(z_profile, Hx_pert_numerical_z, 'k-', linewidth=2, label='NGSolve')
    ax4.plot(z_profile, Hx_pert_analytical_z, 'r--', linewidth=1.5, label='Analytical')
    ax4.axvline(-sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax4.axvline(sphere_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_xlabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_ylabel('$H_{x,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax4.set_title('Z-axis Profile ($H_x$ component)', fontname='Times New Roman', fontsize=11)
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
