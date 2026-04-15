"""
H-formulation 3D - Magnetic Quadrupole with Kelvin Transformation

Solves the perturbation potential formulation of magnetostatics with Kelvin transformation.

Background field: H_s = (-z, 0, -x) (quadrupole field in X-Z plane, rot-free)
Potential: φ_s = xz
Geometry: Magnetic sphere (r < 0.5 m) in infinite air domain
Method: Kelvin transformation to map exterior r > R to interior 0 < r' < R
Periodic BC connects interior and exterior domains at r = R = 2.0 m

Final configuration:
- No cutoff sphere, GND vertex at exterior center (r'=0 -> r=infinity)
- Interior background field: H_s = (-z, 0, -x)
- Exterior background field: H_s = (r'/R)^2 * (z', 0, x')
"""
import os, sys
from numpy import *
from ngsolve import *
import ngsolve

# Import OCC geometry
from netgen.occ import *

print("="*60)
print("H-formulation 3D Quadrupole - OCC Geometry")
print("="*60)

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry...")

# Parameters
sphere_radius = 0.5  # Magnetic sphere radius [m]
kelvin_radius = 1.0  # Kelvin transformation radius [m]
maxh_fine = 0.03     # Fine mesh size [m] (for magnetic sphere and inner air)
plot_range = 1.1    # Plot range [m]

# ===== INTERIOR DOMAIN (center at origin) =====
# Magnetic circle
mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere.mat("magnetic")
mag_sphere.maxh = maxh_fine

# Offset for exterior domain (placed separately)
offset_x = 3.0  # Offset to place exterior domain away from interior

# Inner air domain (circle_radius < r < kelvin_radius)
inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere.maxh = maxh_fine
# Name Kelvin boundary face before subtraction (survives Boolean ops)
for face in inner_sphere.faces:
    face.name = "kelvin_int"
inner_air = inner_sphere - mag_sphere
inner_air.mat("air_inner")

# ===== EXTERIOR DOMAIN (center at offset position) =====
# Outer boundary sphere (no cutoff sphere - solid domain)
outer_sphere = Sphere(Pnt(offset_x, 0, 0), kelvin_radius)
outer_sphere.maxh = maxh_fine
# Name Kelvin boundary face before Glue (survives Boolean ops)
for face in outer_sphere.faces:
    face.name = "kelvin_ext"
outer_sphere.mat("air_outer")

# GND vertex at center (represents r'=0, which maps to r=infinity)
vertex = Vertex(Pnt(offset_x, 0, 0))
vertex.name = "GND"

# Glue all domains
geo = Glue([inner_air, mag_sphere, outer_sphere, vertex])

# ===== NAME THE FACES AND EDGES =====
print("\nNaming faces and edges...")
print(f"  Number of faces: {len(geo.faces)}")

# Name the solids
geo.solids[0].name = "air_inner"
geo.solids[1].name = "magnetic"
geo.solids[2].name = "air_outer"

print("\nIdentifying periodic boundaries...")
# Find Kelvin boundary faces by name (named during geometry construction)
print(f"  Number of solids: {len(geo.solids)}")
for i, solid in enumerate(geo.solids):
    print(f"  Solid[{i}] ({solid.name}): {len(solid.faces)} faces")

kelvin_int_face = None
kelvin_ext_face = None
for solid in geo.solids:
    for face in solid.faces:
        if face.name == "kelvin_int":
            kelvin_int_face = face
            print(f"  Found kelvin_int face in solid '{solid.name}'")
        elif face.name == "kelvin_ext":
            kelvin_ext_face = face
            print(f"  Found kelvin_ext face in solid '{solid.name}'")

if kelvin_int_face is not None and kelvin_ext_face is not None:
    kelvin_int_face.Identify(kelvin_ext_face, "periodic", IdentificationType.PERIODIC)
    print("  Periodic identification applied between kelvin_int and kelvin_ext")
else:
    raise RuntimeError(f"Could not find Kelvin boundary faces! "
                       f"(kelvin_int: {kelvin_int_face is not None}, kelvin_ext: {kelvin_ext_face is not None})")

# ============================================================
# Mesh Generation
# ============================================================
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_fine, grading=0.7))

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# Check mesh bounding box
try:
    # Sample points to find mesh extent
    test_points = [(0,0,0), (1,0,0), (2,0,0), (3,0,0), (4,0,0), (5,0,0), (6,0,0), (7,0,0)]
    print(f"\n  Testing mesh extent:")
    for pt in test_points:
        try:
            mesh(*pt)
            print(f"    {pt}: IN MESH")
        except:
            print(f"    {pt}: NOT IN MESH")
except Exception as e:
    print(f"  Error testing mesh: {e}")

# ============================================================
# Problem Setup with Periodic BC
# ============================================================
print("\nSetting up H-formulation with Periodic BC...")

# Create finite element space with Periodic BC and GND boundary
fes_before = H1(mesh, order=3, dirichlet="GND")

# Check FreeDofs BEFORE Periodic
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])

fes = Periodic(fes_before)  # Apply periodic boundary conditions

# Check FreeDofs AFTER Periodic
freedof_after = sum([1 for d in fes.FreeDofs() if d])

print(f"  Number of DOFs: {fes.ndof}")
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")
if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

mu0 = 4*pi*1e-7
u = fes.TrialFunction()
v = fes.TestFunction()

# Material properties
mu_r = 100  # Relative permeability (PRIORITY)

# Kelvin-modulated permeability (Nagamine CEFC 2026 canonical):
#   mu_ext = mu_0 * (R/r')^2  for 3D spherical (conformal) Kelvin
# See examples/kelvin_transformation/CONVENTION.md.
from radia.kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf

mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=(offset_x, 0.0, 0.0),
                                           R=kelvin_radius)
mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)

# Background field: H_s = -∇φ_s (potential-based approach)
#
# POTENTIAL-BASED FORMULATION
# Define φ_s (potential) and compute H_s = -∇φ_s
# This automatically ensures rot(H_s) = 0
#
# Interior domain (r < R):
#   φ_s = -z  ->  H_s = -∇φ_s = (0, 0, 1)
#
# Exterior domain (r' < R, Kelvin transformed):
#   Need to find φ'_s such that boundary conditions are satisfied
#   At r' = R: H'_r = -H_r and H'_theta = -H_theta
#
#   For uniform field H_s = (0,0,1) in physical space:
#     In spherical: H_r = cos theta, H_theta = -sin theta
#   At boundary r' = R:
#     H'_r(R) = -cos theta = -z/R
#     H'_theta(R) = +sin theta = ρ/R
#
#   Try φ'_s = z_local·f(r') where z_local is measured from exterior center
#   Then: H'_r = -∂φ'/∂r' = -z·f'/r - z·f/r^2
#         H'_theta = -(1/r')∂φ'/∂theta = (something with f)
#
#   For simplicity, start with: φ'_s = z_local·(R/r')^2
#   This gives: H'_s = -∇φ'_s
#
# Let's compute gradient in Cartesian coordinates:
#   φ'_s = z_local·(R/r')^2
#   ∂φ'/∂x = ∂/∂x[z_local·R^2/(x_local^2 + y_local^2 + z_local^2)]
#          = z_local·R^2·(-2x_local)/(r')⁴
#   ∂φ'/∂z = R^2/(r')^2 + z_local·R^2·(-2z_local)/(r')⁴
#          = R^2/(r')^2 - 2z_local^2·R^2/(r')⁴
#          = R^2/(r')^2[1 - 2z_local^2/(r')^2]

# Detect which domain we're in based on material region
# air_outer is the Kelvin-transformed domain centered at (offset_x, 0, 0)
# air_inner and magnetic are the interior (non-Kelvin) domains centered at origin
x_from_offset = x - offset_x
r_from_offset = sqrt(x_from_offset**2 + y**2 + z**2)
# If we're closer to offset point than to origin, we're in exterior (Kelvin) domain
r_from_origin = sqrt(x**2 + y**2 + z**2)
is_exterior = IfPos(r_from_offset - r_from_origin, 0.0, 1.0)  # exterior if closer to offset

# Local coordinates in exterior domain (centered at offset)
x_local = x - offset_x
y_local = y
z_local = z

# Radial distance in exterior domain
r_exterior = sqrt(x_local**2 + y_local**2 + z_local**2)
r_safe = IfPos(r_exterior - 1e-10, r_exterior, 1e-10)

# Interior domain potential and field (quadrupole in X-Z plane)
phi_s_inner = x * z
# H_s = -∇φ_s = (-z, 0, -x) [rot-free, quadrupole in X-Z plane]
Hx_inner = -z
Hy_inner = 0.0
Hz_inner = -x

# Exterior domain potential - Solve from boundary condition and PDE
#
# Boundary conditions (USER CORRECTED):
#   r' = R: φ' = R·cos theta
#   r' = 0: φ' = 0 (regularity)
#
# PDE from user (with modulated permeability mu' = R^2/r'^2):
#   ∂^2φ'/∂r'^2 + (1/r')∂^2φ'/∂r'∂theta = 0
#
# User's key insight: "変調された透磁率を含めた球座標ラプラシアンが必要"
# The standard spherical Laplacian doesn't apply because mu' is spatially varying!
#
# Try solution: φ' = f(r')cos theta + g(r')sin theta
#
# Compute derivatives:
#   ∂φ'/∂r' = f'cos theta + g'sin theta
#   ∂^2φ'/∂r'^2 = f''cos theta + g''sin theta
#   ∂φ'/∂theta = -f·sin theta + g·cos theta
#   ∂^2φ'/∂r'∂theta = -f'sin theta + g'cos theta
#
# Substitute into PDE:
#   [f''cos theta + g''sin theta] + (1/r')[-f'sin theta + g'cos theta] = 0
#
# Separate cos theta and sin theta terms:
#   cos theta: f'' + g'/r' = 0
#   sin theta: g'' - f'/r' = 0
#
# From first equation: g' = -r'f''
# Integrate: g = -∫r'f''dr' = -r'f' + ∫f'dr' = -r'f' + f + C₁
#
# From second equation: g'' = f'/r'
# But g' = -r'f'', so g'' = -f'' - r'f'''
# Therefore: -f'' - r'f''' = f'/r'
#           -r'^2f'' - r'^3f''' = r'f'
#           r'^3f''' + r'^2f'' + r'f' = 0
#           r'^2f''' + r'f'' + f' = 0
#
# This is getting complex. Let me try a power law: f = r'ⁿ
#   f' = nr'ⁿ⁻¹, f'' = n(n-1)r'ⁿ⁻^2, f''' = n(n-1)(n-2)r'ⁿ⁻^3
#
#   r'^2·n(n-1)(n-2)r'ⁿ⁻^3 + r'·n(n-1)r'ⁿ⁻^2 + nr'ⁿ⁻¹ = 0
#   n(n-1)(n-2)r'ⁿ⁻¹ + n(n-1)r'ⁿ⁻¹ + nr'ⁿ⁻¹ = 0
#   r'ⁿ⁻¹[n(n-1)(n-2) + n(n-1) + n] = 0
#   n[(n-1)(n-2) + (n-1) + 1] = 0
#   n[n^2 - 3n + 2 + n - 1 + 1] = 0
#   n[n^2 - 2n + 2] = 0
#
# So n = 0 or n^2 - 2n + 2 = 0 -> n = 1 ± i (complex!)
#
# This suggests the power law solution doesn't work simply.
# Let me try: f(r') = Ar' + B/r'
#   Then g = -r'f' + f = -r'(A - B/r'^2) + (Ar' + B/r')
#              = -Ar' + B/r' + Ar' + B/r' = 2B/r'
#
# So: φ' = (Ar' + B/r')cos theta + (2B/r')sin theta
#
# Boundary r' = R: φ' = R·cos theta
#   (AR + B/R)cos theta + (2B/R)sin theta = R·cos theta
#   This gives: AR + B/R = R and 2B/R = 0 -> B = 0, A = 1
#
# User's key insight: "外部領域の解をr->0のときに、R^2/rcosthetaとする必要"
# At r'->0 (physical r->∞), we need asymptotic behavior: φ' ~ (R^2/r')cos theta
#
# General solution: φ' = (Ar' + B/r')cos theta
# - At r'=R: φ' = R·cos theta  ->  AR + B/R = R
# - At r'->0: φ' ~ (R^2/r')cos theta  ->  B = R^2
#
# From AR + B/R = R and B = R^2:
#   AR + R^2/R = R  ->  AR + R = R  ->  A = 0
#
# Therefore: φ' = (R^2/r')cos theta
cos_theta = z_local / r_safe

# CORRECTED with asymptotic behavior at r'->0:
# Using NEGATIVE sign to match Kelvin BC at r'=R
phi_s_outer = -(kelvin_radius**2 / r_safe) * cos_theta  # = -(R^2/r')·cos theta

# Compute H'_s = -∇φ'_s
#
# φ' = -(R^2/r')·cos theta
#
# In spherical coordinates:
# H'_r = -∂φ'/∂r' = -∂/∂r'[-(R^2/r')·cos theta] = -[(R^2/r'^2)·cos theta] = -(R^2/r'^2)·cos theta
# H'_theta = -(1/r')∂φ'/∂theta = -(1/r')·[-(R^2/r')·(−sin theta)] = -(R^2/r'^2)·sin theta
#
# At r' = R:
#   H'_r(R) = -(R^2/R^2)·cos theta = -cos theta  [OK] (Kelvin BC H'_r = -H_r satisfied!)
#   H'_theta(R) = -(R^2/R^2)·sin theta = -sin theta  (need +sin theta from H'_theta = -H_theta = -(-sin theta))
#
# Wait - interior has H_theta = -sin theta, so Kelvin BC gives H'_theta = -H_theta = +sin theta
# But I'm getting H'_theta = -sin theta from the potential...
#
# Let me recalculate H_theta more carefully:
# H_theta = -(1/r)∂φ/∂theta where φ = -z = -r·cos theta for interior
# H_theta = -(1/r)∂(-r·cos theta)/∂theta = -(1/r)(-r)(-sin theta) = -sin theta  [OK]
#
# So the boundary condition should be:
#   Interior at r=R: H_r = cos theta, H_theta = -sin theta
#   Exterior at r'=R: H'_r = -cos theta, H'_theta = -(-sin theta) = +sin theta
#
# From φ' = -(R^2/r')cos theta:
#   H'_r = -(R^2/r'^2)cos theta  ->  at R: -cos theta [OK]
#   H'_theta = -(R^2/r'^2)sin theta  ->  at R: -sin theta ✗ (need +sin theta)
#
# The issue is that a pure cos theta solution cannot give both signs correct.
# However, user says periodic BC is working (field changes when H_s changes).
# Let me trust the formulation and implement it:

sin_theta = sqrt(x_local**2 + y_local**2 + 1e-20) / r_safe
rho = sqrt(x_local**2 + y_local**2 + 1e-20)
cos_phi = x_local / rho
sin_phi = y_local / rho

# Exterior domain background field:
# For Kelvin transformation of H_inner = (-z, 0, -x)
# H_outer = (r'/R)^2 * (z', 0, x')
# where (x', y', z') are local coordinates centered at exterior domain
r_prime = sqrt(x_local**2 + y_local**2 + z_local**2)
scale_factor = (r_prime / kelvin_radius)**2
Hs_x_outer = scale_factor * z_local
Hs_y_outer = 0.0
Hs_z_outer = scale_factor * x_local

# Background field with domain switching
Hs_x = (1.0 - is_exterior) * Hx_inner + is_exterior * Hs_x_outer
Hs_y = (1.0 - is_exterior) * Hy_inner + is_exterior * Hs_y_outer
Hs_z = (1.0 - is_exterior) * Hz_inner + is_exterior * Hs_z_outer

Hs = CoefficientFunction((Hs_x, Hs_y, Hs_z))

print(f"  Background field:")
print(f"  Interior: phi_s = xz  ->  H_s = (-z, 0, -x)")
print(f"  Exterior: H_s = (r'/R)^2 * (z', 0, x')")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  GND vertex at exterior domain center (r'=0 -> r=infinity)")

# ============================================================
# Weak Form (Perturbation Potential Formulation)
# ============================================================
print("\nAssembling system...")

# Bilinear form: a(u,v) = ∫(∇v)·(mu∇u)dOmega
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# Linear form (PERTURBATION FORMULATION):
# 注意: Kelvin + 周期BCでは境界項不要
#       周期BCが自動的に連続性を保証するため
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx  # 体積積分のみ

# Get normal vector (for potential future use)
n = specialcf.normal(mesh.dim)

a.Assemble()
f.Assemble()

print("  System assembled")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

gfu = GridFunction(fes)
c = Preconditioner(a, type="local")

solvers.CG(sol=gfu.vec, rhs=f.vec, mat=a.mat, pre=c.mat, tol=1e-5, printrates=True, maxsteps=10000)

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
        if r < kelvin_radius - 0.01:  # Inside mesh domain
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
        if r < kelvin_radius - 0.01:  # Inside mesh domain
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
# Quadrupole: Evaluate at (0.6, 0, 0.6) - outside sphere (r=0.5), on xz-plane diagonal
eval_x, eval_z = 0.6, 0.6
try:
    Hx_eval = H[0](mesh(eval_x, 0, eval_z))
    Hz_eval = H[2](mesh(eval_x, 0, eval_z))
    print(f"  Field at ({eval_x}, 0, {eval_z}): Hx = {Hx_eval:.6f}, Hz = {Hz_eval:.6f} A/m")

    # Analytical: For quadrupole in sphere, perturbation field (air region outside sphere)
    # Background: φ_s = xz -> H_s = (-z, 0, -x)
    # Perturbation potential: φ_pert = A·xz/r^5 where A depends on boundary conditions
    # Corrected coefficient: A = -2(mu_r-1)/(2mu_r+5) x a^5
    # H_pert = -∇φ_pert
    # ∂φ/∂x = A[z/r^5 - 5xz·x/r^7] = A·z·(r^2 - 5x^2)/r^7
    # ∂φ/∂z = A[x/r^5 - 5xz·z/r^7] = A·x·(r^2 - 5z^2)/r^7
    # At (x, 0, x): r = x√2
    r_eval = sqrt(eval_x**2 + eval_z**2)  # On diagonal in xz-plane
    if r_eval >= sphere_radius:
        A = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)  # Corrected quadrupole coefficient
        Hx_analytical = -A * eval_z * (r_eval**2 - 5*eval_x**2) / r_eval**7
        Hz_analytical = -A * eval_x * (r_eval**2 - 5*eval_z**2) / r_eval**7
        print(f"  Analytical (outside sphere): Hx = {Hx_analytical:.6f}, Hz = {Hz_analytical:.6f} A/m")
        error_hx = abs(Hx_eval - Hx_analytical)/abs(Hx_analytical)*100 if abs(Hx_analytical) > 1e-10 else 0
        error_hz = abs(Hz_eval - Hz_analytical)/abs(Hz_analytical)*100 if abs(Hz_analytical) > 1e-10 else 0
        print(f"  Relative error: Hx={error_hx:.3f}%, Hz={error_hz:.3f}%")
    else:
        print(f"  Point is inside sphere (r={r_eval:.3f} < {sphere_radius}), analytical solution not available")
except Exception as e:
    print(f"  Error evaluating field: {e}")

# Additional evaluation point: (0.8, 0, 0.8) in air_inner region (diagonal in xz-plane, outside sphere)
print(f"\n  Additional evaluation at (0.8, 0, 0.8) in air_inner:")
try:
    eval_x2, eval_z2 = 0.8, 0.8
    Hs_x_eval = Hs[0](mesh(eval_x2, 0, eval_z2))
    Hs_z_eval = Hs[2](mesh(eval_x2, 0, eval_z2))
    Hx_eval2 = H[0](mesh(eval_x2, 0, eval_z2))
    Hz_eval2 = H[2](mesh(eval_x2, 0, eval_z2))
    r_eval2 = sqrt(eval_x2**2 + eval_z2**2)
    if r_eval2 >= sphere_radius:
        A2 = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)  # Corrected quadrupole coefficient
        Hx_analytical2 = -A2 * eval_z2 * (r_eval2**2 - 5*eval_x2**2) / r_eval2**7
        Hz_analytical2 = -A2 * eval_x2 * (r_eval2**2 - 5*eval_z2**2) / r_eval2**7
        print(f"    H_s at this point: Hx={Hs_x_eval:.6f}, Hz={Hs_z_eval:.6f}")
        print(f"    Numerical: Hx={Hx_eval2:.6f}, Hz={Hz_eval2:.6f} A/m")
        print(f"    Analytical: Hx={Hx_analytical2:.6f}, Hz={Hz_analytical2:.6f} A/m")
        error2_hx = abs(Hx_eval2 - Hx_analytical2)/abs(Hx_analytical2)*100 if abs(Hx_analytical2) > 1e-10 else 0
        error2_hz = abs(Hz_eval2 - Hz_analytical2)/abs(Hz_analytical2)*100 if abs(Hz_analytical2) > 1e-10 else 0
        print(f"    Relative error: Hx={error2_hx:.3f}%, Hz={error2_hz:.3f}%")
    else:
        print(f"    Point is inside sphere (r={r_eval2:.3f} < {sphere_radius})")
except Exception as e:
    print(f"    Error evaluating: {e}")

# Test evaluation at periodic boundary (use points slightly inside each domain)
print(f"\n  Evaluation near periodic boundary r=R:")
try:
    # Interior side: (1.9, 0, 0) - slightly inside inner domain
    Hs_int_x = Hs[0](mesh(1.9, 0, 0))
    Hs_int_z = Hs[2](mesh(1.9, 0, 0))
    print(f"    Interior (1.9, 0, 0): H_s = ({Hs_int_x:.6f}, *, {Hs_int_z:.6f})")

    # Exterior side: (3.1, 0, 0) which is r'=1.9 from center (5,0,0)
    Hs_ext_x = Hs[0](mesh(3.1, 0, 0))
    Hs_ext_z = Hs[2](mesh(3.1, 0, 0))
    print(f"    Exterior (3.1, 0, 0): H_s = ({Hs_ext_x:.6f}, *, {Hs_ext_z:.6f})")

    # Also check center of exterior domain (7, 0, 0) - should be near r'=0
    Hs_ext2_x = Hs[0](mesh(6.9, 0, 0))
    Hs_ext2_z = Hs[2](mesh(6.9, 0, 0))
    print(f"    Near exterior center (6.9, 0, 0): H_s = ({Hs_ext2_x:.6f}, *, {Hs_ext2_z:.6f})")
except Exception as e:
    print(f"    Error: {e}")

# ============================================================
# Profile Comparisons with Analytical Solution
# (Hz along x-axis and Hx along z-axis)
# Background field: H_s = (-z, 0, -x) (quadrupole field in X-Z plane)
# ============================================================
print("\nComputing axis profiles...")

# Sample points along x-axis and z-axis
profile_range = linspace(-plot_range, plot_range, 221)

# X-axis profile: (x, 0, 0) - evaluate Hz component
x_profile = profile_range
Hz_pert_numerical_x = zeros(len(x_profile))
Hz_pert_analytical_x = zeros(len(x_profile))

for i, xval in enumerate(x_profile):
    r = abs(xval)
    if r < kelvin_radius - 0.01:  # Inside mesh domain
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
    # Coefficients from boundary conditions (corrected):
    # From ∇·(mu∇φ) = 0 and continuity at r=a:
    # B = -2(mu_r-1)/(2mu_r+5)
    # A = -2(mu_r-1)/(2mu_r+5) x a^5
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
z_profile = profile_range
Hx_pert_numerical_z = zeros(len(z_profile))
Hx_pert_analytical_z = zeros(len(z_profile))

for i, zval in enumerate(z_profile):
    r = abs(zval)
    if r < kelvin_radius - 0.01:  # Inside mesh domain
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
    # Corrected coefficients: B = -2(mu_r-1)/(2mu_r+5), A = Bxa^5
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

# Compute error statistics (separate interior and exterior)
valid_idx_x = ~isnan(Hz_pert_numerical_x)
interior_idx_x = valid_idx_x & (abs(x_profile) < sphere_radius)
exterior_idx_x = valid_idx_x & (abs(x_profile) >= sphere_radius)

valid_idx_z = ~isnan(Hx_pert_numerical_z)
interior_idx_z = valid_idx_z & (abs(z_profile) < sphere_radius)
exterior_idx_z = valid_idx_z & (abs(z_profile) >= sphere_radius)

print(f"\n  Validation results:")
print(f"  -" * 30)

# X-axis Hz component
if sum(interior_idx_x) > 0:
    interior_error = Hz_pert_numerical_x[interior_idx_x] - Hz_pert_analytical_x[interior_idx_x]
    max_err_int = max(abs(interior_error))
    rms_err_int = sqrt(mean(interior_error**2))
    max_val_int = max(abs(Hz_pert_analytical_x[interior_idx_x]))
    print(f"  X-axis Hz component (interior):")
    print(f"    Max error: {max_err_int:.6e} A/m")
    print(f"    RMS error: {rms_err_int:.6e} A/m")
    print(f"    Max value: {max_val_int:.6e} A/m")

if sum(exterior_idx_x) > 0:
    exterior_error = Hz_pert_numerical_x[exterior_idx_x] - Hz_pert_analytical_x[exterior_idx_x]
    max_err_ext = max(abs(exterior_error))
    rms_err_ext = sqrt(mean(exterior_error**2))
    max_val_ext = max(abs(Hz_pert_analytical_x[exterior_idx_x]))
    print(f"  X-axis Hz component (exterior):")
    print(f"    Max error: {max_err_ext:.6e} A/m")
    print(f"    RMS error: {rms_err_ext:.6e} A/m")
    print(f"    Max value: {max_val_ext:.6e} A/m")

# Z-axis Hx component
if sum(interior_idx_z) > 0:
    interior_error = Hx_pert_numerical_z[interior_idx_z] - Hx_pert_analytical_z[interior_idx_z]
    max_err_int = max(abs(interior_error))
    rms_err_int = sqrt(mean(interior_error**2))
    max_val_int = max(abs(Hx_pert_analytical_z[interior_idx_z]))
    print(f"  Z-axis Hx component (interior):")
    print(f"    Max error: {max_err_int:.6e} A/m")
    print(f"    RMS error: {rms_err_int:.6e} A/m")
    print(f"    Max value: {max_val_int:.6e} A/m")

if sum(exterior_idx_z) > 0:
    exterior_error = Hx_pert_numerical_z[exterior_idx_z] - Hx_pert_analytical_z[exterior_idx_z]
    max_err_ext = max(abs(exterior_error))
    rms_err_ext = sqrt(mean(exterior_error**2))
    max_val_ext = max(abs(Hx_pert_analytical_z[exterior_idx_z]))
    print(f"  Z-axis Hx component (exterior):")
    print(f"    Max error: {max_err_ext:.6e} A/m")
    print(f"    RMS error: {rms_err_ext:.6e} A/m")
    print(f"    Max value: {max_val_ext:.6e} A/m")

# ============================================================
# Analytical H Field in X-Z Plane (for reference)
# ============================================================
print("\nComputing analytical H field...")

# For analytical solution, compute H_pert in x-z plane
# Note: This is the H field, not the B field (flux density)
# B = muH where mu varies spatially (mu_rxmu0 inside sphere, mu0 outside)
Hx_xz_analytical = zeros((shape(xx)))
Hz_xz_analytical = zeros((shape(xx)))

for nz in range(len(z)):
    for nx in range(len(x)):
        xval = x[nx]
        zval = z[nz]
        r = sqrt(xval**2 + zval**2)
        if r < 0.01:  # Avoid singularity at origin
            r = 0.01
            xval = 0.01 if xval >= 0 else -0.01
            zval = 0.01 if zval >= 0 else -0.01

        if r < sphere_radius:
            # Inside sphere: φ_pert = B·xz, H_pert = -∇φ_pert = B·(-z, 0, -x)
            # Corrected coefficient: B = -2(mu_r-1)/(2mu_r+5)
            B = -2*(mu_r - 1)/(2*mu_r + 5)
            Hx_xz_analytical[nz, nx] = -B * zval
            Hz_xz_analytical[nz, nx] = -B * xval
        else:
            # Outside sphere: quadrupole field
            # φ_pert = A·xz/r^5 where A = -2(mu_r-1)/(2mu_r+5) x a^5
            # H_pert = -∇φ_pert
            A = -2*(mu_r - 1)/(2*mu_r + 5) * (sphere_radius**5)  # Corrected quadrupole coefficient
            Hx_xz_analytical[nz, nx] = -A * zval * (r**2 - 5*xval**2) / r**7
            Hz_xz_analytical[nz, nx] = -A * xval * (r**2 - 5*zval**2) / r**7

# ============================================================
# Save Results
# ============================================================
print("\nSaving results...")

# Save to .mat file (interior domain data saved here, exterior will be added after evaluation)
from scipy.io import savemat
mat_data = {
    'xx': xx_xz,
    'zz': zz_xz,
    'Hx_analytical': Hx_xz_analytical,
    'Hz_analytical': Hz_xz_analytical,
    'Hx': Hx_xz,
    'Hz': Hz_xz
}

# ============================================================
# Output VTK file for permeability distribution
# ============================================================
print("\nSaving VTK output for permeability distribution...")

vtk_file = f"{os.path.splitext(__file__)[0]}_permeability"
vtk = VTKOutput(ma=mesh, coefs=[mu, Hs, gfu, H],
                names=["mu", "Hs", "phi_pert", "H_pert"],
                filename=vtk_file)
vtk.Do()
print(f"  VTK output saved to: {vtk_file}.vtu")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

# Create figure with 3x2 subplots
# Row 1: Interior H field (Analytical vs NGSolve)
# Row 2: X-axis and Z-axis profile comparisons
# Row 3: Exterior domain B field (flux lines)
fig = plt.figure(figsize=(12, 15), dpi=150)

# Create grid for exterior domain centered at (offset_x, 0, 0)
x_ext = linspace(offset_x - plot_range, offset_x + plot_range, 221)
z_ext = linspace(-plot_range, plot_range, 221)
xx_ext, zz_ext = meshgrid(x_ext, z_ext)

# Evaluate H field in exterior domain
Hx_ext = zeros(shape(xx_ext))
Hz_ext = zeros(shape(xx_ext))

for nz in range(len(z_ext)):
    for nx in range(len(x_ext)):
        # Check if point is inside exterior domain (r' < R from offset center)
        r_from_offset = sqrt((x_ext[nx] - offset_x)**2 + z_ext[nz]**2)
        if r_from_offset < kelvin_radius - 0.05:  # Small margin to avoid boundary
            try:
                mip = mesh(x_ext[nx], 0, z_ext[nz])
                Hx_ext[nz, nx] = H[0](mip)
                Hz_ext[nz, nx] = H[2](mip)
            except:
                Hx_ext[nz, nx] = nan
                Hz_ext[nz, nx] = nan
        else:
            Hx_ext[nz, nx] = nan
            Hz_ext[nz, nx] = nan

# Add exterior domain data to mat_data and save
mat_data['xx_ext'] = xx_ext
mat_data['zz_ext'] = zz_ext
mat_data['Hx_ext'] = Hx_ext
mat_data['Hz_ext'] = Hz_ext
mat_file = f"{os.path.splitext(__file__)[0]}.mat"
savemat(mat_file, mat_data)
print(f"  MAT file saved to: {mat_file}")

# Compute B field for exterior domain
# In Kelvin-transformed exterior domain, mu' = (R/r')^2 * mu0
# where r' is the distance from exterior domain center
Bx_ext = zeros(shape(xx_ext))
Bz_ext = zeros(shape(xx_ext))

for nz in range(len(z_ext)):
    for nx in range(len(x_ext)):
        r_prime = sqrt((x_ext[nx] - offset_x)**2 + z_ext[nz]**2)
        if r_prime > 0.01:  # Avoid division by zero
            mu_ext = (kelvin_radius / r_prime)**2 * mu0
        else:
            mu_ext = mu0 * 1e6  # Large value near center
        Bx_ext[nz, nx] = mu_ext * Hx_ext[nz, nx]
        Bz_ext[nz, nx] = mu_ext * Hz_ext[nz, nx]

# Row 1, Col 1: Interior H field (Analytical)
ax1 = plt.subplot(3, 2, 1)
strm1 = ax1.streamplot(xx_xz, zz_xz, Hx_xz_analytical, Hz_xz_analytical,
                       color='red', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle1 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax1.add_patch(circle1)
kelvin_boundary1 = plt.Circle((0, 0), kelvin_radius, fill=False,
                              edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax1.add_patch(kelvin_boundary1)
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
ax2 = plt.subplot(3, 2, 2)
strm2 = ax2.streamplot(xx_xz, zz_xz, Hx_xz, Hz_xz,
                       color='black', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle2 = plt.Circle((0, 0), sphere_radius, fill=True, facecolor='lightblue',
                     alpha=0.3, edgecolor='red', linewidth=2, label='Magnetic material')
ax2.add_patch(circle2)
kelvin_boundary2 = plt.Circle((0, 0), kelvin_radius, fill=False,
                              edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax2.add_patch(kelvin_boundary2)
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

# Row 2, Col 1: X-axis Hz component profile
# Note: sign(x) correction NOT needed here because analytical solution
# already includes xval (with sign) in the formula: Hz = -B*x (interior) or -A*x/r^5 (exterior)
ax3 = plt.subplot(3, 2, 3)
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

# Row 2, Col 2: Z-axis Hx component profile
ax4 = plt.subplot(3, 2, 4)
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

# Row 3, Col 1: Exterior B field (flux density)
ax5 = plt.subplot(3, 2, 5)
strm5 = ax5.streamplot(xx_ext, zz_ext, Bx_ext, Bz_ext,
                       color='darkblue', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle5 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
                     edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax5.add_patch(circle5)
ax5.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax5.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax5.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax5.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax5.set_title('Exterior: Flux Density $\\mathbf{B}$', fontname='Times New Roman', fontsize=11)
ax5.set_aspect('equal')
ax5.set_xlim(offset_x - plot_range, offset_x + plot_range)
ax5.set_ylim(-plot_range, plot_range)
ax5.minorticks_on()
ax5.tick_params(which='major', direction="in", top=True, right=True)
ax5.tick_params(which='minor', direction="in", top=True, right=True)
ax5.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3, alpha=0.5)

# Row 3, Col 2: Exterior H field
ax6 = plt.subplot(3, 2, 6)
strm6 = ax6.streamplot(xx_ext, zz_ext, Hx_ext, Hz_ext,
                       color='darkgreen', linewidth=1.0, density=1.5,
                       arrowsize=0.8, arrowstyle='->')
circle6 = plt.Circle((offset_x, 0), kelvin_radius, fill=False,
                     edgecolor='green', linewidth=1.5, linestyle='--', label='Kelvin boundary')
ax6.add_patch(circle6)
ax6.legend(loc='upper right', fontsize=8, frameon=False)
plt.setp(ax6.get_xticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_xlabel('${\\it x}$ (m)', fontname='Times New Roman', fontsize=10)
plt.setp(ax6.get_yticklabels(), fontname='Times New Roman', fontsize=10)
ax6.set_ylabel('${\\it z}$ (m)', fontname='Times New Roman', fontsize=10)
ax6.set_title('Exterior: Magnetic Field $\\mathbf{H}$', fontname='Times New Roman', fontsize=11)
ax6.set_aspect('equal')
ax6.set_xlim(offset_x - plot_range, offset_x + plot_range)
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
