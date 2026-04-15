"""
Omega-Reduced Omega Method for 3D Magnetostatics with Kelvin Transformation
Problem: Magnetic sphere (mu_r=100) in uniform z-directed background field

Based on Omega_ReducedOmega.py implementation:
- Sign convention: B = mu * grad(Omega), H = grad(Omega)
- Source potential: Omega_s = H0 * z (so grad(Omega_s) = H_s)
- Total region (magnetic sphere): No source term in weak form
- Reduced region (air): Source term from Omega_s

Kelvin transformation:
- Maps infinite exterior domain to finite sphere
- Permeability transformation: mu'(r') = (R/r')^2 * mu0
- Periodic BC couples interior (r=R) with exterior (r'=R)
- Dirichlet BC at exterior center (r'=0 -> r=infinity): Omega = 0

IMPORTANT: 3D Normal Direction
- In NGSolve 3D, specialcf.normal(mesh.dim) on the sphere boundary points INWARD
  (from air into magnetic sphere), but Omega-Reduced Omega requires OUTWARD normal.
- Therefore, we NEGATE the normal: normal = -specialcf.normal(mesh.dim)
- This is in contrast to 2D axisymmetric where the normal already points outward.

Analytical solution:
- Inside sphere: Hz = 3/(mu_r + 2) * H0
- Outside sphere: dipole + uniform field

Perturbation field energy:
- Interior: W_in = (1/2) * mu_r * mu0 * [(mu_r-1)/(mu_r+2)]^2 * H0^2 * V_sphere
- Exterior: W_out = mu0 * m^2 / (12*pi*a^3), where m = 4*pi*a^3*(mu_r-1)/(mu_r+2)*H0

Author: Claude Code
Date: 2025-12-23
"""
import os
from numpy import *
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("Omega-Reduced Omega Method - 3D Magnetic Sphere with Kelvin Transform")
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
maxh_sphere = 0.05       # Finer mesh for magnetic sphere
maxh_air = 0.08          # Moderate mesh for air region
fe_order = 2             # Finite element order

# Offset for exterior domain (z-direction to match axisymmetric)
offset_z = 3.0

print(f"\nProblem parameters:")
print(f"  Sphere radius: {sphere_radius} m")
print(f"  Kelvin radius: {kelvin_radius} m")
print(f"  Relative permeability: mu_r = {mu_r}")
print(f"  Source field: H_s = (0, 0, {H0}) A/m")

# ============================================================
# Geometry Definition (OCC)
# ============================================================
print("\nCreating geometry using OCC...")

# Interior domain: magnetic sphere at origin
mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
mag_sphere.mat("magnetic")
mag_sphere.maxh = maxh_sphere
# Name the magnetic sphere surface
for face in mag_sphere.faces:
    face.name = "sphere"

# Inner air domain (annulus between magnetic sphere and Kelvin boundary)
inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
inner_sphere.maxh = maxh_air
# Name the Kelvin boundary (interior side)
for face in inner_sphere.faces:
    face.name = "kelvin_int"

inner_air = inner_sphere - mag_sphere
inner_air.mat("air_inner")

# Exterior domain (Kelvin-transformed, centered at z-offset position)
outer_sphere = Sphere(Pnt(0, 0, offset_z), kelvin_radius)
outer_sphere.maxh = maxh_air
outer_sphere.mat("air_outer")
# Name the Kelvin boundary (exterior side)
for face in outer_sphere.faces:
    face.name = "kelvin_ext"

# GND vertex at center of exterior domain (represents r'=0 -> r=infinity)
vertex = Vertex(Pnt(0, 0, offset_z))
vertex.name = "GND"

# Glue all domains
geo = Glue([inner_air, mag_sphere, outer_sphere, vertex])

# Name the solids
geo.solids[0].name = "air_inner"
geo.solids[1].name = "magnetic"
geo.solids[2].name = "air_outer"

# ===== IDENTIFY PERIODIC FACES =====
print("\nIdentifying periodic boundaries...")

# Print solid and face information
print(f"  Number of solids: {len(geo.solids)}")
for i, solid in enumerate(geo.solids):
    print(f"  Solid[{i}] ({solid.name}): {len(solid.faces)} faces")
    for j, face in enumerate(solid.faces):
        print(f"    Face[{j}]: name='{face.name}'")

# Find kelvin_int and kelvin_ext faces by name
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
    print(f"  WARNING: Could not find periodic faces!")
    print(f"    kelvin_int_face: {kelvin_int_face}")
    print(f"    kelvin_ext_face: {kelvin_ext_face}")

# Generate mesh
print("\nGenerating mesh...")
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5))

# Apply curved elements for high-order accuracy
mesh.Curve(fe_order)

print(f"  Number of elements: {mesh.ne}")
print(f"  Number of vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {mesh.GetBoundaries()}")
print(f"  Curved elements: order {fe_order}")


# ============================================================
# Finite Element Space
# ============================================================
print("\nSetting up finite element space...")

# H1 space with Dirichlet BC at GND (infinity)
fes_before = H1(mesh, order=fe_order, dirichlet="GND")
fes = Periodic(fes_before)  # Apply periodic BC

# Check if Periodic BC is working
freedof_before = sum([1 for d in fes_before.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  Finite element order: {fe_order}")
print(f"  Number of DOFs: {fes.ndof}")
print(f"  FreeDofs: {freedof_before} -> {freedof_after} (diff: {freedof_before - freedof_after})")
if freedof_before == freedof_after:
    print("  WARNING: Periodic BC may NOT be working!")
else:
    print("  Periodic BC is working (FreeDofs reduced)")

# Trial and test functions
Omega = fes.TrialFunction()
psi = fes.TestFunction()

# ============================================================
# Material Properties
# ============================================================
print("\nSetting up material properties...")

# Kelvin-transformed permeability (Nagamine CEFC 2026 canonical):
#   mu_ext = mu_0 * (R/r')^2  for 3D spherical (conformal) Kelvin
# See examples/kelvin_transformation/CONVENTION.md.
from radia.kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf

mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=(0.0, 0.0, offset_z),
                                           R=kelvin_radius)
Mu = build_material_cf(
    mesh, mu0, mu_kelvin_factor,
    outer_keyword="air_outer",
    overrides={"magnetic": mu_r * mu0},
)
mu_kelvin = mu0 * mu_kelvin_factor   # CF alias for downstream energy integration

print(f"  air_inner: mu = mu0")
print(f"  air_outer: mu = (R/r')^2 * mu0 [Nagamine CEFC 2026]")
print(f"  magnetic: mu = {mu_r} * mu0")

# ============================================================
# Following Omega_ReducedOmega.py implementation
# ============================================================
print("\nSetting up weak form (following Omega_ReducedOmega.py)...")

# Source potential Omega_s = H0 * z (so grad(Omega_s) = (0, 0, H0) = H_s)
Omega_s = H0 * z

# Source magnetic field and B field
Hs = CoefficientFunction((0.0, 0.0, H0))
Bs = CoefficientFunction((0.0, 0.0, mu0 * H0))

# For Kelvin-transformed exterior domain (z-offset):
r_exterior = sqrt(x**2 + y**2 + (z - offset_z)**2 + 1e-20)
Hz_exterior = -(r_exterior / kelvin_radius)**2 * H0
Hs_exterior = CoefficientFunction((0.0, 0.0, Hz_exterior))
Bs_exterior = CoefficientFunction((0.0, 0.0, mu0 * Hz_exterior))

# ===== Omega-Reduced Omega Method (following Omega_ReducedOmega.py) =====
# Total region (magnetic): No source term, B = mu * grad(Omega)
# Reduced region (air): Source term from Omega_s, H = grad(Omega)

# Detect exterior domain for Kelvin (z-offset)
z_from_offset = z - offset_z
r_from_offset = sqrt(x**2 + y**2 + z_from_offset**2)
r_from_origin = sqrt(x**2 + y**2 + z**2)
is_exterior = IfPos(z - offset_z/2, 1.0, 0.0)  # Simple z-based detection

# Bilinear form: mu * grad(Omega) * grad(psi) for all regions
a = BilinearForm(fes)
a += Mu * grad(Omega) * grad(psi) * dx("magnetic")
a += Mu * grad(Omega) * grad(psi) * dx("air_inner")
a += Mu * grad(Omega) * grad(psi) * dx("air_outer")
a.Assemble()

# Create GridFunction for solution
# Set Omega_s on sphere boundary (same as Omega_ReducedOmega.py line 96)
gfOmega = GridFunction(fes)
gfOmega.Set(Omega_s, BND, mesh.Boundaries("sphere"))

# Linear form: source term in Reduced regions + Neumann BC
# NOTE: Omega_ReducedOmega.py does NOT include Kelvin region in source term
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
f.Assemble()

# Normal for Neumann BC
# IMPORTANT: In 3D NGSolve, the normal on the sphere boundary points INWARD
# (from air into magnetic), but we need it to point OUTWARD for the
# Omega-Reduced Omega formulation. So we NEGATE the normal.
# This is in contrast to 2D axisymmetric where the normal already points outward.
normal = -specialcf.normal(mesh.dim)  # NEGATED for correct direction

# Neumann boundary condition on sphere surface (Total/Reduced interface)
f += (normal * Bs) * psi * ds("sphere")
f.Assemble()

print("  Omega-Reduced Omega formulation:")
print("  Total region (magnetic): H = grad(Omega), no source")
print("  Reduced region (air): Source from Omega_s")
print("  Bilinear: mu*grad(Omega)*grad(psi) dx (all regions)")
print("  Linear: mu*grad(Omega_s)*grad(psi) dx (Reduced only)")

# ============================================================
# Solve
# ============================================================
print("\nSolving system...")

# Use the same gfOmega that has Dirichlet BC set
gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
gfu = gfOmega

print("  Solution converged")

# ============================================================
# Post-processing (following Axisymmetric lines 345-373)
# ============================================================
print("\nPost-processing (following Axisymmetric)...")

# Following Axisymmetric lines 345-373:
# fesOt = H1(mesh, order=fe_order, definedon="magnetic")
# fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")
# Ot.Set(gfu, VOL, definedon="magnetic")
# Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
# Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))
# Bt = grad(Ot) * Mu
# Br = (grad(Orr) - grad(Oxr)) * mu0
# BField = Bt + Br + Bs_cf

fesOt = H1(mesh, order=fe_order, definedon="magnetic")
fesOr = H1(mesh, order=fe_order, definedon="air_inner|air_outer")

Ot = GridFunction(fesOt)
Orr = GridFunction(fesOr)
Oxr = GridFunction(fesOr)

Ot.Set(gfu, VOL, definedon="magnetic")
Orr.Set(gfu, VOL, definedon="air_inner|air_outer")
Oxr.Set(Omega_s, BND, mesh.Boundaries("sphere"))

# B field computation (same as Axisymmetric lines 361-373)
Bt = grad(Ot) * Mu  # Total region
Br = (grad(Orr) - grad(Oxr)) * mu0  # Reduced region perturbation

# Source field Bs is only in Reduced region
Bs_dict = {
    "air_inner": Bs,
    "air_outer": Bs_exterior,
    "magnetic": CoefficientFunction((0.0, 0.0, 0.0))
}
Bs_cf = CoefficientFunction([Bs_dict[mat] for mat in mesh.GetMaterials()])

# Total B field
BField = Bt + Br + Bs_cf

# Also compute grad_Omega for validation
grad_Omega = grad(gfu)

# ============================================================
# Analytical Solution and Validation
# ============================================================
Hz_analytical_interior = 3.0 / (mu_r + 2) * H0  # = 0.029412 A/m

print("\n" + "=" * 60)
print("VALIDATION RESULTS")
print("=" * 60)
print(f"\nAnalytical solution:")
print(f"  Interior Hz = 3/(mu_r+2) * H0 = {Hz_analytical_interior:.6f} A/m")

# Evaluate at origin (center of magnetic sphere)
print(f"\nInterior (magnetic sphere - Total potential region):")
print(f"  In Total region: H = grad(Omega)")
print()


print()
for point_name, coords in [("origin", (0, 0, 0)), ("(0.2, 0, 0)", (0.2, 0, 0)), ("(0, 0, 0.3)", (0, 0, 0.3))]:
    try:
        mip = mesh(coords[0], coords[1], coords[2])
        # In Total region: H = grad(Omega)
        Hz_numerical = grad_Omega[2](mip)
        error = abs(Hz_numerical - Hz_analytical_interior) / abs(Hz_analytical_interior) * 100
        print(f"  {point_name}: Hz={Hz_numerical:.6f}, analytical={Hz_analytical_interior:.6f}, error={error:.2f}%")
    except Exception as e:
        print(f"  {point_name}: Error - {e}")

# Evaluate at exterior points
print(f"\nExterior (air - Reduced potential region):")
print(f"  In Reduced region: H = grad(Omega) - grad(Omega_s) + H_s")
print()

for r_val in [0.6, 0.7, 0.8, 0.9]:
    # On x-axis: Hz = H0 * (1 - (mu_r-1)/(mu_r+2) * (a/r)^3)
    Hz_analytical_ext = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r_val)**3)
    try:
        mip = mesh(r_val, 0, 0)
        # In Reduced region: total H = grad(Omega) - grad(Omega_s) + H_s
        # Since H_s = grad(Omega_s), this simplifies to H = grad(Omega)
        # But we need to use BField / mu0 for consistency
        B_val = BField(mip)
        Hz_numerical = B_val[2] / mu0
        error = abs(Hz_numerical - Hz_analytical_ext) / abs(Hz_analytical_ext) * 100
        print(f"  x={r_val}: Hz={Hz_numerical:.6f}, analytical={Hz_analytical_ext:.6f}, error={error:.2f}%")
    except Exception as e:
        print(f"  x={r_val}: Error - {e}")

# ============================================================
# Profile along z-axis
# ============================================================
print("\nComputing z-axis profile...")

z_vals = linspace(-0.9, 0.9, 51)
Hz_numerical_z = []
Hz_analytical_z = []

for zv in z_vals:
    r = abs(zv)
    try:
        mip = mesh(0, 0, zv)
        if r < sphere_radius:
            # Total region: H = grad(Omega)
            Hz_numerical_z.append(grad_Omega[2](mip))
        else:
            # Reduced region: use BField / mu0
            B_val = BField(mip)
            Hz_numerical_z.append(B_val[2] / mu0)
    except:
        Hz_numerical_z.append(nan)

    # Analytical (total field)
    if r < sphere_radius:
        Hz_ana = Hz_analytical_interior
    elif r > 0.01:
        # On z-axis (theta=0): H_z = H0 * [1 + 2*(mu_r-1)/(mu_r+2)*(a/r)^3]
        Hz_ana = H0 * (1.0 + 2 * (mu_r - 1) / (mu_r + 2) * (sphere_radius / r)**3)
    else:
        Hz_ana = nan
    Hz_analytical_z.append(Hz_ana)

Hz_numerical_z = array(Hz_numerical_z)
Hz_analytical_z = array(Hz_analytical_z)

# ============================================================
# Profile along x-axis
# ============================================================
print("Computing x-axis profile...")

x_vals = linspace(-0.9, 0.9, 51)
Hz_numerical_x = []
Hz_analytical_x = []

for xv in x_vals:
    r = abs(xv)
    try:
        mip = mesh(xv, 0, 0)
        if r < sphere_radius:
            # Total region: H = grad(Omega)
            Hz_numerical_x.append(grad_Omega[2](mip))
        else:
            # Reduced region: use BField / mu0
            B_val = BField(mip)
            Hz_numerical_x.append(B_val[2] / mu0)
    except:
        Hz_numerical_x.append(nan)

    # Analytical (total field)
    if r < sphere_radius:
        Hz_ana = Hz_analytical_interior
    elif r > 0.01:
        # On x-axis (theta=90deg): H_z = H0 * [1 - (mu_r-1)/(mu_r+2)*(a/r)^3]
        Hz_ana = H0 * (1.0 - (mu_r - 1) / (mu_r + 2) * (sphere_radius / r)**3)
    else:
        Hz_ana = nan
    Hz_analytical_x.append(Hz_ana)

Hz_numerical_x = array(Hz_numerical_x)
Hz_analytical_x = array(Hz_analytical_x)

# ============================================================
# VTK Output
# ============================================================
print("\nSaving VTK output...")

vtk_file = os.path.splitext(__file__)[0] + "_result"
VTKOutput(ma=mesh, coefs=[Mu, gfu, grad(gfu)],
          names=["mu", "Omega", "grad_Omega"],
          filename=vtk_file).Do()
print(f"  VTK file saved: {vtk_file}.vtu")

# ============================================================
# Visualization (will be done after perturbation calculation)
# ============================================================
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                             'bf': 'serif:bold', 'fontset': 'cm'})

# ============================================================
# Perturbation Field Energy Calculation
# ============================================================
print("\nCalculating perturbation field energy...")

# Perturbation field definition:
# Total region (magnetic): B_pert = mu*grad(Omega) - Bs, H_pert = grad(Omega) - Hs
# Reduced region (air): B_pert = mu0*grad(Omega_r), H_pert = grad(Omega_r)
#
# Energy: W = (1/2) * integral(B_pert * H_pert) dV
#
# In NGSolve with Omega-Reduced Omega:
# - Total region: grad(Omega) is the total H, perturbation is grad(Omega) - Hs
# - Reduced region: grad(Omega) already gives the perturbation (since source is subtracted)

# --- Total region (magnetic sphere) ---
# H_pert = grad(Omega) - Hs = grad(Omega) - (0, 0, H0)
# B_pert = mu * H_pert
grad_Omega_total = grad(gfu)
H_pert_total = grad_Omega_total - Hs
B_pert_total = (mu_r * mu0) * H_pert_total

# Energy in total region: (1/2) * B_pert * H_pert = (1/2) * mu * |H_pert|^2
energy_total = Integrate(0.5 * (mu_r * mu0) * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)

# --- Reduced region (air_inner) ---
# In Reduced region, the perturbation field is: H_pert = grad(Orr) - grad(Oxr)
H_pert_reduced = grad(Orr) - grad(Oxr)

# Energy: (1/2) * mu0 * |H_pert|^2
energy_reduced = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)

# --- Kelvin region (air_outer) ---
# In Kelvin region, need to account for the transformation
# The perturbation in Kelvin region is grad(Orr)
H_pert_kelvin = grad(Orr)
# Kelvin変換後のエネルギー: (1/2) * mu'(r') * |H'|^2 dV'
# mu'(r') = (R/r')^2 * mu0 は既に mu_kelvin として定義済み
energy_kelvin = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)

# Total perturbation energy
energy_total_pert = energy_total + energy_reduced + energy_kelvin

print(f"\nPerturbation field energy:")
print(f"  Total region (magnetic):   W_t = {energy_total:.6e} J")
print(f"  Reduced region (air_inner): W_r = {energy_reduced:.6e} J")
print(f"  Kelvin region (air_outer):  W_k = {energy_kelvin:.6e} J")
print(f"  Total perturbation energy:  W   = {energy_total_pert:.6e} J")

# ============================================================
# Analytical Energy for Comparison
# ============================================================
# 磁性体球の摂動場エネルギー解析解
#
# 球内部の摂動磁場: H_pert = H_in - H_s = (3/(mu_r+2) - 1) * H0 = -(mu_r-1)/(mu_r+2) * H0
# 球内部の摂動磁束密度: B_pert = mu_r * mu0 * H_pert
#
# 球内部のエネルギー:
# W_in = (1/2) * mu_r * mu0 * |H_pert|^2 * V_sphere
#      = (1/2) * mu_r * mu0 * [(mu_r-1)/(mu_r+2)]^2 * H0^2 * (4*pi/3) * a^3
#
# 球外部（空気中）の双極子場エネルギー:
# W_out = integral of (1/2) * mu0 * |H_dipole|^2 dV (from a to infinity)
# 双極子場: H_dipole = (m/4*pi) * {...} where m = (4*pi/3)*a^3 * chi * H0
# 積分結果: W_out = (1/6) * mu0 * chi^2 * H0^2 * V_sphere
#
# 総エネルギー: W_total = W_in + W_out

chi_sphere = 3 * (mu_r - 1) / (mu_r + 2)
V_sphere = (4.0/3.0) * pi * sphere_radius**3
H_pert_in = -(mu_r - 1) / (mu_r + 2) * H0

# 球内部エネルギー解析解
W_in_analytical = 0.5 * mu_r * mu0 * H_pert_in**2 * V_sphere

# 球外部エネルギー解析解 (双極子場)
# 双極子モーメント: m = 4*pi*a^3 * (mu_r-1)/(mu_r+2) * H0
# 双極子場エネルギー: W_out = mu0 * m^2 / (12*pi*a^3)
m_dipole = 4 * pi * sphere_radius**3 * (mu_r - 1) / (mu_r + 2) * H0
W_out_analytical = mu0 * m_dipole**2 / (12 * pi * sphere_radius**3)

W_total_analytical = W_in_analytical + W_out_analytical

print(f"\nAnalytical comparison:")
print(f"  Susceptibility chi = 3*(mu_r-1)/(mu_r+2) = {chi_sphere:.6f}")
print(f"  Sphere volume = {V_sphere:.6e} m^3")
print(f"  H_pert inside = -{(mu_r-1)/(mu_r+2):.6f} * H0 = {H_pert_in:.6f} A/m")
print(f"\nAnalytical energies:")
print(f"  W_in (magnetic sphere):  {W_in_analytical:.6e} J")
print(f"  W_out (exterior dipole): {W_out_analytical:.6e} J")
print(f"  W_total (analytical):    {W_total_analytical:.6e} J")
print(f"\nComparison:")
print(f"  Numerical W_t / Analytical W_in = {energy_total / W_in_analytical:.4f}")
print(f"  Numerical (W_r + W_k) / Analytical W_out = {(energy_reduced + energy_kelvin) / W_out_analytical:.4f}")

# ============================================================
# Perturbation Field Visualization
# ============================================================
print("\nComputing perturbation field profiles...")

# Z-axis profile for perturbation field
Hz_pert_numerical_z = []
Hz_pert_analytical_z = []

for zv in z_vals:
    r = abs(zv)
    try:
        mip = mesh(0, 0, zv)
        if r < sphere_radius:
            # Total region: H_pert = grad(Omega) - Hs
            Hz_num = grad_Omega[2](mip) - H0
        else:
            # Reduced region: H_pert = grad(Omega_r) = Br/mu0
            Hz_num = Br[2](mip) / mu0
        Hz_pert_numerical_z.append(Hz_num)
    except:
        Hz_pert_numerical_z.append(nan)

    # Analytical perturbation field
    if r < sphere_radius:
        # Inside: H_pert = H_in - H0 = 3/(mu_r+2)*H0 - H0 = -(mu_r-1)/(mu_r+2)*H0
        Hz_pert_ana = -(mu_r - 1) / (mu_r + 2) * H0
    elif r > 0.01:
        # On z-axis: H_z = H0 + dipole field = H0 * [1 + 2*(mu_r-1)/(mu_r+2)*(a/r)^3]
        # Perturbation = H_z - H0 = 2*(mu_r-1)/(mu_r+2)*(a/r)^3 * H0
        Hz_pert_ana = 2 * (mu_r - 1) / (mu_r + 2) * (sphere_radius / r)**3 * H0
    else:
        Hz_pert_ana = nan
    Hz_pert_analytical_z.append(Hz_pert_ana)

Hz_pert_numerical_z = array(Hz_pert_numerical_z)
Hz_pert_analytical_z = array(Hz_pert_analytical_z)

# X-axis profile for perturbation field
Hz_pert_numerical_x = []
Hz_pert_analytical_x = []

for xv in x_vals:
    r = abs(xv)
    try:
        mip = mesh(xv, 0, 0)
        if r < sphere_radius:
            Hz_num = grad_Omega[2](mip) - H0
        else:
            Hz_num = Br[2](mip) / mu0
        Hz_pert_numerical_x.append(Hz_num)
    except:
        Hz_pert_numerical_x.append(nan)

    if r < sphere_radius:
        Hz_pert_ana = -(mu_r - 1) / (mu_r + 2) * H0
    elif r > 0.01:
        # On x-axis (equatorial): H_z = H0 * [1 - (mu_r-1)/(mu_r+2)*(a/r)^3]
        # Perturbation = H_z - H0 = -(mu_r-1)/(mu_r+2)*(a/r)^3 * H0
        Hz_pert_ana = -(mu_r - 1) / (mu_r + 2) * (sphere_radius / r)**3 * H0
    else:
        Hz_pert_ana = nan
    Hz_pert_analytical_x.append(Hz_pert_ana)

Hz_pert_numerical_x = array(Hz_pert_numerical_x)
Hz_pert_analytical_x = array(Hz_pert_analytical_x)

# ============================================================
# Generate Combined Plot (Total Field + Perturbation Field)
# ============================================================
print("\nGenerating plots...")

fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=150)

# Row 1: Total field profiles
ax1 = axes[0, 0]
ax1.plot(z_vals, Hz_numerical_z, 'b-', linewidth=2, label='NGSolve')
ax1.plot(z_vals, Hz_analytical_z, 'r--', linewidth=1.5, label='Analytical')
ax1.axvline(-sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax1.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax1.set_xlabel('$z$ (m)', fontsize=11)
ax1.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax1.set_title('Z-axis (Total Field)', fontsize=12)
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

ax2 = axes[0, 1]
ax2.plot(x_vals, Hz_numerical_x, 'b-', linewidth=2, label='NGSolve')
ax2.plot(x_vals, Hz_analytical_x, 'r--', linewidth=1.5, label='Analytical')
ax2.axvline(-sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax2.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax2.set_xlabel('$x$ (m)', fontsize=11)
ax2.set_ylabel('$H_z$ (A/m)', fontsize=11)
ax2.set_title('X-axis (Total Field)', fontsize=12)
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

# Energy comparison text
ax3 = axes[0, 2]
ax3.axis('off')
energy_text = f"""Perturbation Field Energy

Numerical:
  W_in  = {energy_total:.4e} J
  W_out = {energy_reduced + energy_kelvin:.4e} J
  W_tot = {energy_total_pert:.4e} J

Analytical:
  W_in  = {W_in_analytical:.4e} J
  W_out = {W_out_analytical:.4e} J
  W_tot = {W_total_analytical:.4e} J

Ratio (Num/Ana):
  W_in:  {energy_total / W_in_analytical:.4f}
  W_out: {(energy_reduced + energy_kelvin) / W_out_analytical:.4f}
"""
ax3.text(0.1, 0.5, energy_text, transform=ax3.transAxes, fontsize=10,
         verticalalignment='center', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Row 2: Perturbation field profiles
ax4 = axes[1, 0]
ax4.plot(z_vals, Hz_pert_numerical_z, 'b-', linewidth=2, label='NGSolve')
ax4.plot(z_vals, Hz_pert_analytical_z, 'r--', linewidth=1.5, label='Analytical')
ax4.axvline(-sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax4.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax4.axhline(0, color='k', linestyle='-', alpha=0.3)
ax4.set_xlabel('$z$ (m)', fontsize=11)
ax4.set_ylabel('$H_{z,pert}$ (A/m)', fontsize=11)
ax4.set_title('Z-axis (Perturbation Field)', fontsize=12)
ax4.legend(loc='best', fontsize=9)
ax4.grid(True, alpha=0.3)

ax5 = axes[1, 1]
ax5.plot(x_vals, Hz_pert_numerical_x, 'b-', linewidth=2, label='NGSolve')
ax5.plot(x_vals, Hz_pert_analytical_x, 'r--', linewidth=1.5, label='Analytical')
ax5.axvline(-sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax5.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax5.axhline(0, color='k', linestyle='-', alpha=0.3)
ax5.set_xlabel('$x$ (m)', fontsize=11)
ax5.set_ylabel('$H_{z,pert}$ (A/m)', fontsize=11)
ax5.set_title('X-axis (Perturbation Field)', fontsize=12)
ax5.legend(loc='best', fontsize=9)
ax5.grid(True, alpha=0.3)

# Error plot
ax6 = axes[1, 2]
valid_z = ~isnan(Hz_numerical_z) & ~isnan(Hz_analytical_z) & (abs(Hz_analytical_z) > 1e-10)
valid_x = ~isnan(Hz_numerical_x) & ~isnan(Hz_analytical_x) & (abs(Hz_analytical_x) > 1e-10)
if sum(valid_z) > 0:
    with errstate(divide='ignore', invalid='ignore'):
        error_z = abs(Hz_numerical_z[valid_z] - Hz_analytical_z[valid_z]) / abs(Hz_analytical_z[valid_z]) * 100
    error_z = where(isfinite(error_z), error_z, 0)
    ax6.semilogy(z_vals[valid_z], error_z + 1e-6, 'b-', linewidth=2, label='Z-axis')
if sum(valid_x) > 0:
    with errstate(divide='ignore', invalid='ignore'):
        error_x = abs(Hz_numerical_x[valid_x] - Hz_analytical_x[valid_x]) / abs(Hz_analytical_x[valid_x]) * 100
    error_x = where(isfinite(error_x), error_x, 0)
    ax6.semilogy(x_vals[valid_x], error_x + 1e-6, 'g--', linewidth=2, label='X-axis')
ax6.axvline(-sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax6.axvline(sphere_radius, color='gray', linestyle=':', alpha=0.7)
ax6.set_xlabel('Position (m)', fontsize=11)
ax6.set_ylabel('Relative Error (%)', fontsize=11)
ax6.set_title('Error Distribution', fontsize=12)
ax6.legend(loc='best', fontsize=9)
ax6.grid(True, alpha=0.3)

plt.tight_layout()

# Save figure
png_file = os.path.splitext(__file__)[0] + ".png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved: {png_file}")

# Open the plot
try:
    os.startfile(png_file)
except:
    pass

print("\n" + "=" * 60)
print("Computation completed")
print("=" * 60)
