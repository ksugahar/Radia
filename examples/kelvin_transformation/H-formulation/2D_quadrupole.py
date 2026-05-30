"""
H-formulation for 2D magnetostatics with perturbation potential
Geometry created internally using OCC (2D circle)
Updated: 2025-11-23

NOTE: This code attempts to use H_s = (x, -y) as quadrupole background field.
However, simple polynomial fields like (y, -x) or (x, -y) produce negligible
perturbations for circular geometries due to symmetry constraints.

For proper quadrupole field simulation, use the Kelvin transformation approach:
  - 2D_quadrupole_with_Kelvin.py

The Kelvin transformation maps the infinite domain to a finite one, allowing
a uniform field in the transformed space to become a quadrupole in physical space.
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve
from ngsolve import TaskManager

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 2D - OCC Geometry")
print("="*60)

# ============================================================
# Geometry Definition (OCC 2D)
# ============================================================
print("\nCreating geometry...")

# Parameters
circle_radius = 0.5      # Magnetic circle radius [m]
air_inner_radius = 1.0   # Outer boundary radius [m]
maxh_fine = 0.03         # Fine mesh size [m]

# Create magnetic circle (finest mesh)
wp = WorkPlane()
mag_circle_shape = wp.Circle(circle_radius).Face()
mag_circle_shape.maxh = maxh_fine

# Create air circle and name its boundary
air_circle_shape = wp.Circle(air_inner_radius).Face()
for edge in air_circle_shape.edges:
    edge.name = "outer"  # Name outer boundary before boolean operation
air_circle_shape.maxh = maxh_fine

# Boolean operations to create two regions
# Air = air circle - magnetic circle
air_shape = air_circle_shape - mag_circle_shape
air_shape.name = "air"

# Set magnetic material
mag_circle_shape.name = "magnetic"

# Add grounding vertex at origin to fix potential uniqueness
vertex = Vertex(Pnt(0, 0, 0))
vertex.name = "GND"

# Combine into single geometry
geo = Glue([air_shape, mag_circle_shape, vertex])

print(f"Geometry created with two regions:")
print(f"  Magnetic circle radius: {circle_radius} m")
print(f"  Outer boundary radius: {air_inner_radius} m")
print(f"  Fine mesh size: {maxh_fine} m")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Problem Setup
# ============================================================
print("\nSetting up H-formulation...")

n = specialcf.normal(mesh.dim)

# Define relative permeability
mu_r = 10  # Relative permeability

# Use boundary-based named grounding (GND vertex)
fes = H1(mesh, order=3, dirichlet_bbnd="GND")
print(f"  Number of DOFs: {fes.ndof}")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties (mu_r already defined above for BC)
# Note: OCC geometry uses 'default' for air region
mu_d = {"default": 1*mu0, "air": 1*mu0, "magnetic": mu_r*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# Background field: H_s = (x, -y) A/m (QUADRUPOLE field)
# This satisfies div(H_s) = 0, making it physically valid
# In polar: H_r = r cos(2theta), H_theta = -r sin(2theta)
# Corresponds to potential φ_s = -(1/2)r^2 cos(2theta)
Hs = CoefficientFunction((x, -y))
Hsb = BoundaryFromVolumeCF(Hs)  # Boundary values for integration

print(f"  Background field: H_s = (x, -y) A/m (quadrupole)")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Note: div(H_s) = 0 (solenoidal field)")

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
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項（必須）

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
               tol=1e-8, printrates=True, maxsteps=10000)

    print("  Solution converged")

    # ============================================================
    # Post-processing
    # ============================================================
    print("\nPost-processing...")

    # Create evaluation grid (adjusted for outer boundary at r=1.0)
    plot_range = 1.1  # Plot range [m]
    x = linspace(-plot_range, plot_range, 221)
    y = linspace(-plot_range, plot_range, 221)
    [xx, yy] = meshgrid(x, y)

    # Compute perturbation field: H_pert = -grad(phi)
    H = -grad(gfu)
    Hx = zeros((shape(xx)))
    Hy = zeros((shape(xx)))
    phi = zeros((shape(xx)))

    for ny in range(len(y)):
        for nx in range(len(x)):
            r = sqrt(x[nx]**2 + y[ny]**2)
            if r < air_inner_radius - 0.01:  # Inside mesh domain
                try:
                    mip = mesh(x[nx], y[ny])
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
    center_value = gfu(mesh(0, 0))
    phi = phi - center_value

    print(f"  Potential at origin: {center_value:.6e}")

    # Test points for field comparison
    # For H_s = (x, -y), analytical solution:
    # Interior: φ_pert = A r^2 cos2theta, where A = (mu_r-1)/(2(mu_r+1))
    # H_pert = -∇φ, H_r = -2Ar cos2theta, H_theta = 2Ar sin2theta
    # H_pert = -2A(x cos2theta cos - sin2theta(-sin), -y cos2theta sin - sin2theta cos)

    A_coeff = (mu_r - 1.0)/(2.0*(mu_r + 1.0))
    B_coeff = A_coeff  # Same coefficient

    # Test at (x, y) = (0.4, 0): H_s = (0.4, 0), theta=0, cos2theta=1
    x_test1 = 0.4
    y_test1 = 0.0
    r_test1 = sqrt(x_test1**2 + y_test1**2)
    theta_test1 = arctan2(y_test1, x_test1)  # theta=0

    Hs_x_test1 = x_test1  # H_s = (x, -y)
    Hs_y_test1 = -y_test1
    Hx_pert_test1 = H[0](mesh(x_test1, y_test1))
    Hy_pert_test1 = H[1](mesh(x_test1, y_test1))

    # Interior: H_pert,r = -2A r cos2theta = -2A x 0.4 x 1 = -0.8A
    # H_pert,theta = 2A r sin2theta = 2A x 0.4 x 0 = 0
    # At theta=0: H_x = H_r = -0.8A, H_y = H_theta = 0
    Hx_pert_analytical1 = -2.0 * A_coeff * r_test1 * cos(2*theta_test1)
    Hy_pert_analytical1 = 2.0 * A_coeff * r_test1 * sin(2*theta_test1)

    print(f"\n  Field at (x,y) = ({x_test1}, {y_test1}):")
    print(f"    H_s = ({Hs_x_test1:.3f}, {Hs_y_test1:.3f}) A/m")
    print(f"    H_pert (NGSolve) = ({Hx_pert_test1:.6f}, {Hy_pert_test1:.6f}) A/m")
    print(f"    H_pert (Analytical) = ({Hx_pert_analytical1:.6f}, {Hy_pert_analytical1:.6f}) A/m")
    print(f"    Error: ({abs(Hx_pert_test1-Hx_pert_analytical1):.6e}, {abs(Hy_pert_test1-Hy_pert_analytical1):.6e}) A/m")

    # Test at (x, y) = (0, 0.4): H_s = (0, -0.4), theta=pi/2, cos2theta=-1
    x_test2 = 0.0
    y_test2 = 0.4
    r_test2 = sqrt(x_test2**2 + y_test2**2)
    theta_test2 = arctan2(y_test2, x_test2)

    Hs_x_test2 = x_test2
    Hs_y_test2 = -y_test2
    Hx_pert_test2 = H[0](mesh(x_test2, y_test2))
    Hy_pert_test2 = H[1](mesh(x_test2, y_test2))

    # Interior: H_pert,r = -2A r cos2theta = -2A x 0.4 x (-1) = 0.8A
    # H_pert,theta = 2A r sin2theta = 2A x 0.4 x 0 = 0
    # At theta=pi/2: H_x = 0, H_y = H_r = 0.8A
    Hx_pert_analytical2 = -2.0 * A_coeff * r_test2 * cos(2*theta_test2) * cos(theta_test2) + 2.0 * A_coeff * r_test2 * sin(2*theta_test2) * sin(theta_test2)
    Hy_pert_analytical2 = -2.0 * A_coeff * r_test2 * cos(2*theta_test2) * sin(theta_test2) - 2.0 * A_coeff * r_test2 * sin(2*theta_test2) * cos(theta_test2)

    print(f"\n  Field at (x,y) = ({x_test2}, {y_test2}):")
    print(f"    H_s = ({Hs_x_test2:.3f}, {Hs_y_test2:.3f}) A/m")
    print(f"    H_pert (NGSolve) = ({Hx_pert_test2:.6f}, {Hy_pert_test2:.6f}) A/m")
    print(f"    H_pert (Analytical) = ({Hx_pert_analytical2:.6f}, {Hy_pert_analytical2:.6f}) A/m")
    print(f"    Error: ({abs(Hx_pert_test2-Hx_pert_analytical2):.6e}, {abs(Hy_pert_test2-Hy_pert_analytical2):.6e}) A/m")

    # Test at (x, y) = (0.7, 0) - exterior point
    x_test3 = 0.7
    y_test3 = 0.0
    r_test3 = sqrt(x_test3**2 + y_test3**2)
    theta_test3 = arctan2(y_test3, x_test3)

    Hx_pert_test3 = H[0](mesh(x_test3, y_test3))
    Hy_pert_test3 = H[1](mesh(x_test3, y_test3))

    # Exterior: H_r = 2B(a⁴/r^3)cos2theta, H_theta = -2B(a⁴/r^3)sin2theta
    # At theta=0: H_r = 2B(a⁴/r^3), H_theta = 0 -> H_x = H_r, H_y = 0
    Hx_pert_analytical3 = 2.0 * B_coeff * (circle_radius**4 / r_test3**3) * cos(2*theta_test3)
    Hy_pert_analytical3 = -2.0 * B_coeff * (circle_radius**4 / r_test3**3) * sin(2*theta_test3)

    print(f"\n  Field at (x,y) = ({x_test3}, {y_test3}) [exterior]:")
    print(f"    H_pert (NGSolve) = ({Hx_pert_test3:.6f}, {Hy_pert_test3:.6f}) A/m")
    print(f"    H_pert (Analytical) = ({Hx_pert_analytical3:.6f}, {Hy_pert_analytical3:.6f}) A/m")
    print(f"    Error: ({abs(Hx_pert_test3-Hx_pert_analytical3):.6e}, {abs(Hy_pert_test3-Hy_pert_analytical3):.6e}) A/m")

    # Additional debug points
    print(f"\nDebug: NGSolve solution at specific points:")
    test_pts = [(0.3, 0.0), (0.0, 0.3), (0.8, 0.0), (0.0, 0.8)]
    for (xp, yp) in test_pts:
        try:
            Hx_ngs = H[0](mesh(xp, yp))
            Hy_ngs = H[1](mesh(xp, yp))
            r_pt = sqrt(xp**2 + yp**2)
            region = "interior" if r_pt < circle_radius else "exterior"
            print(f"  ({xp:.1f}, {yp:.1f}) [{region}]: H_NGSolve = ({Hx_ngs:.6f}, {Hy_ngs:.6f})")
        except:
            pass

    # ============================================================
    # Profile Comparisons with Analytical Solution
    # ============================================================
    print("\nComputing axis profiles (perturbation field)...")

    profile_range = linspace(-plot_range, plot_range, 221)

    # X-axis profile: Hx vs x at y=0 (quadrupole field)
    x_profile = profile_range
    Hx_pert_numerical_x = zeros(len(x_profile))
    Hx_pert_analytical_x = zeros(len(x_profile))

    for i, xval in enumerate(x_profile):
        r = abs(xval)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(xval, 0)
                Hx_pert_numerical_x[i] = H[0](mip)  # Hx component at (x, 0)
            except:
                Hx_pert_numerical_x[i] = nan
        else:
            Hx_pert_numerical_x[i] = nan

        # Analytical solution for H_s = (x, -y) at (x, 0): theta=0 or pi
        r = abs(xval)
        theta = arctan2(0, xval)  # theta=0 or pi
        if r < circle_radius:
            # Interior: H_r = -2Ar cos2theta, at y=0: H_x = H_r
            Hx_pert_analytical_x[i] = -2.0 * A_coeff * r * cos(2*theta)
        else:
            # Exterior: H_r = 2B(a⁴/r^3)cos2theta, at y=0: H_x = H_r
            Hx_pert_analytical_x[i] = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)

    # Y-axis profile: Hy vs y at x=0 (quadrupole field)
    y_profile = profile_range
    Hy_pert_numerical_y = zeros(len(y_profile))
    Hy_pert_analytical_y = zeros(len(y_profile))

    for i, yval in enumerate(y_profile):
        r = abs(yval)
        if r < air_inner_radius - 0.01:  # Inside mesh domain
            try:
                mip = mesh(0, yval)
                Hy_pert_numerical_y[i] = H[1](mip)  # Hy component at (0, y)
            except:
                Hy_pert_numerical_y[i] = nan
        else:
            Hy_pert_numerical_y[i] = nan

        # Analytical solution for H_s = (x, -y) at (0, y): theta=pi/2 or -pi/2
        r = abs(yval)
        theta = arctan2(yval, 0)  # theta=±pi/2
        if r < circle_radius:
            # Interior: H_r = -2Ar cos2theta, H_theta = 2Ar sin2theta
            # At x=0: H_y = H_r sintheta + H_theta costheta
            Hr = -2.0 * A_coeff * r * cos(2*theta)
            Htheta = 2.0 * A_coeff * r * sin(2*theta)
            Hy_pert_analytical_y[i] = Hr * sin(theta) + Htheta * cos(theta)
        else:
            # Exterior: H_r = 2B(a⁴/r^3)cos2theta, H_theta = -2B(a⁴/r^3)sin2theta
            Hr = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)
            Htheta = -2.0 * B_coeff * (circle_radius**4 / r**3) * sin(2*theta)
            Hy_pert_analytical_y[i] = Hr * sin(theta) + Htheta * cos(theta)

    # Error statistics
    valid_idx_x = ~isnan(Hx_pert_numerical_x)
    valid_idx_y = ~isnan(Hy_pert_numerical_y)

    if sum(valid_idx_x) > 0:
        error_x = Hx_pert_numerical_x[valid_idx_x] - Hx_pert_analytical_x[valid_idx_x]
        max_err_x = max(abs(error_x))
        rms_err_x = sqrt(mean(error_x**2))
        print(f"\n  X-axis profile errors (Hx at y=0):")
        print(f"    Max error: {max_err_x:.6e} A/m")
        print(f"    RMS error: {rms_err_x:.6e} A/m")

    if sum(valid_idx_y) > 0:
        error_y = Hy_pert_numerical_y[valid_idx_y] - Hy_pert_analytical_y[valid_idx_y]
        max_err_y = max(abs(error_y))
        rms_err_y = sqrt(mean(error_y**2))
        # For relative error, use exterior points only
        exterior_idx_y = valid_idx_y & (abs(y_profile) > circle_radius)
        if sum(exterior_idx_y) > 0:
            avg_abs_val = mean(abs(Hy_pert_analytical_y[exterior_idx_y]))
            rel_err_y = rms_err_y / avg_abs_val * 100 if avg_abs_val > 0 else 0
            print(f"\n  Y-axis profile errors (Hy at x=0):")
            print(f"    Max error: {max_err_y:.6e} A/m")
            print(f"    RMS error: {rms_err_y:.6e} A/m")
            print(f"    Relative error (exterior): {rel_err_y:.3f}%")

    # ============================================================
    # Analytical Flux Lines (2D) - H_s = (x, -y)
    # ============================================================
    print("\nComputing analytical flux lines...")

    Hx_analytical = zeros((shape(xx)))
    Hy_analytical = zeros((shape(xx)))

    for ny in range(len(y)):
        for nx in range(len(x)):
            r = sqrt(x[nx]**2 + y[ny]**2)
            if r < 0.01:
                r = 0.01
            theta = arctan2(y[ny], x[nx])

            if r < circle_radius:
                # Interior: φ_pert = Ar^2cos2theta
                # H_r = -∂φ/∂r = -2Ar cos2theta
                # H_theta = -(1/r)∂φ/∂theta = 2Ar sin2theta
                Hr = -2.0 * A_coeff * r * cos(2*theta)
                Htheta = 2.0 * A_coeff * r * sin(2*theta)
                # Convert to Cartesian
                Hx_analytical[ny, nx] = Hr * cos(theta) - Htheta * sin(theta)
                Hy_analytical[ny, nx] = Hr * sin(theta) + Htheta * cos(theta)
            else:
                # Exterior: φ = B(a⁴/r^2)cos2theta
                # H_r = -∂φ/∂r = 2B(a⁴/r^3)cos2theta
                # H_theta = -(1/r)∂φ/∂theta = 2B(a⁴/r^3)sin2theta
                Hr = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)
                Htheta = 2.0 * B_coeff * (circle_radius**4 / r**3) * sin(2*theta)
                # Convert to Cartesian
                Hx_analytical[ny, nx] = Hr * cos(theta) - Htheta * sin(theta)
                Hy_analytical[ny, nx] = Hr * sin(theta) + Htheta * cos(theta)

    # ============================================================
    # Save Results
    # ============================================================
    print("\nSaving results...")

    # Save to .mat file
    from scipy.io import savemat
    mat_file = f"{os.path.splitext(__file__)[0]}.mat"
    savemat(mat_file, {
        'xx': xx,
        'yy': yy,
        'Hx_analytical': Hx_analytical,
        'Hy_analytical': Hy_analytical,
        'Hx': Hx,
        'Hy': Hy
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

    # Create figure with 2x2 subplots
    # Row 1: Analytical H streamline vs NGSolve H streamline
    # Row 2: X-axis and Y-axis profile comparisons
    fig = plt.figure(figsize=(12, 10), dpi=150)

    # Row 1, Col 1: Analytical H field streamline
    ax1 = plt.subplot(2, 2, 1)
    strm1 = ax1.streamplot(xx, yy, Hx_analytical, Hy_analytical,
                           color='red', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    circle1 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
                         alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
    ax1.add_patch(circle1)
    ax1.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax1.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax1.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
    ax1.set_title('Analytical $\\mathbf{H}_{\\mathrm{pert}}$ Streamline', fontname='Times New Roman', fontsize=11)
    ax1.set_aspect('equal')
    ax1.set_xlim(-plot_range, plot_range)
    ax1.set_ylim(-plot_range, plot_range)
    ax1.minorticks_on()
    ax1.tick_params(which='major', direction="in", top=True, right=True)
    ax1.tick_params(which='minor', direction="in", top=True, right=True)
    ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Row 1, Col 2: NGSolve H field streamline
    ax2 = plt.subplot(2, 2, 2)
    strm2 = ax2.streamplot(xx, yy, Hx, Hy,
                           color='black', linewidth=1.0, density=1.5,
                           arrowsize=0.8, arrowstyle='->')
    circle2 = plt.Circle((0, 0), circle_radius, fill=True, facecolor='lightblue',
                         alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
    ax2.add_patch(circle2)
    ax2.legend(loc='upper right', fontsize=8, frameon=False)
    plt.setp(ax2.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax2.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax2.set_ylabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
    ax2.set_title('NGSolve $\\mathbf{H}_{\\mathrm{pert}}$ Streamline', fontname='Times New Roman', fontsize=11)
    ax2.set_aspect('equal')
    ax2.set_xlim(-plot_range, plot_range)
    ax2.set_ylim(-plot_range, plot_range)
    ax2.minorticks_on()
    ax2.tick_params(which='major', direction="in", top=True, right=True)
    ax2.tick_params(which='minor', direction="in", top=True, right=True)
    ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

    # Apply sign(x) correction for X-axis profile analytical solution only
    Hx_pert_analytical_x_corrected = Hx_pert_analytical_x * sign(x_profile)
    # Handle x=0 case
    Hx_pert_analytical_x_corrected[x_profile == 0] = Hx_pert_analytical_x[x_profile == 0]

    # Row 2, Col 1: X-axis profile comparison
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(x_profile, Hx_pert_numerical_x, 'k-', linewidth=2, label='NGSolve')
    ax3.plot(x_profile, Hx_pert_analytical_x_corrected, 'r--', linewidth=1.5, label='Analytical')
    ax3.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax3.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax3.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax3.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax3.set_ylabel('$H_{x,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax3.set_title('X-axis Profile ($H_x$ component)', fontname='Times New Roman', fontsize=11)
    ax3.minorticks_on()
    ax3.tick_params(which='major', direction="in", top=True, right=True)
    ax3.tick_params(which='minor', direction="in", top=True, right=True)
    ax3.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)
    ax3.grid(axis='both', which='minor', c='gainsboro', linestyle='--', linewidth=0.1)
    ax3.legend(loc='best', fontsize=9, frameon=False)

    # Row 2, Col 2: Y-axis profile comparison
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(y_profile, Hy_pert_numerical_y, 'k-', linewidth=2, label='NGSolve')
    ax4.plot(y_profile, Hy_pert_analytical_y, 'r--', linewidth=1.5, label='Analytical')
    ax4.axvline(-circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax4.axvline(circle_radius, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    plt.setp(ax4.get_xticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_xlabel('${\\it y}$ (m)', fontname='Times New Roman', fontsize=10)
    plt.setp(ax4.get_yticklabels(), fontname='Times New Roman', fontsize=10)
    ax4.set_ylabel('$H_{y,\\mathrm{pert}}$ (A/m)', fontname='Times New Roman', fontsize=10)
    ax4.set_title('Y-axis Profile ($H_y$ component)', fontname='Times New Roman', fontsize=11)
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
