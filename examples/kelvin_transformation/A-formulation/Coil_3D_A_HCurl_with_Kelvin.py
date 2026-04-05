"""
3D A-formulation (HCurl) with Periodic Kelvin transformation.

First 3D HCurl Periodic Kelvin implementation.

Formulation:
  Variable: A (vector potential, HCurl edge elements)
  Equation: curl(nu * curl(A)) = J (source current)

  Interior: nu = nu0 (standard curl-curl)
  Exterior (Kelvin-mapped): nu = nu0 * (r'/R_K)^2
    where R_K = Kelvin radius, r' = distance from exterior center

  Derivation (Sugahara, IEEE Trans. Magn. 2022):
    Conformal symmetry of Maxwell's equations under sphere inversion:
      ALL material constants (mu, epsilon, sigma, nu) scale by (a/r')^2
      mu' = mu0 * (a/r')^2, nu' = nu0 * (a/r')^2
    This is formulation-independent (works for A, H, phi, E).
    Consistent with H-formulation: mu_outer = (a/r')^2 * mu0.

Geometry:
  - Interior sphere (origin): air + coil, radius R_K
  - Exterior sphere (offset): Kelvin-mapped exterior, radius R_K
  - Periodic identification on the two sphere surfaces
  - GND at exterior center (maps to physical infinity)

Usage:
    python Coil_3D_A_HCurl_with_Kelvin.py
"""
import math
import time
from numpy import pi, sqrt as np_sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 60)
print("3D HCurl A-formulation with Periodic Kelvin")
print("Gapped Torus Coil")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
R_coil = 0.030       # Coil center radius [m]
a_coil = 0.003       # Coil wire radius [m]
gap_deg = 5          # Gap angle [degrees]
R_K = 0.06           # Kelvin boundary radius [m]
offset = 0.15        # Offset for exterior domain [m]
maxh = 0.012         # Mesh size [m]
maxh_coil = 0.005    # Coil mesh size [m]
order = 1            # FE order

MU_0 = 4e-7 * pi
NU_0 = 1.0 / MU_0
I_total = 1.0
J0 = I_total / (pi * a_coil**2)

L_analytical = MU_0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) * \
    (360 - gap_deg) / 360

print(f"Coil: R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
print(f"Kelvin: R_K={R_K*1e3:.0f}mm, offset={offset*1e3:.0f}mm")
print(f"L_analytical = {L_analytical*1e9:.2f} nH")
print()

# ============================================================
# Geometry: Interior sphere + Exterior sphere (offset)
# ============================================================
print("[1/4] Building geometry...")
t0 = time.perf_counter()

# Gapped torus (coil volume)
wp = WorkPlane(Axes(p=Pnt(R_coil, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
circle = wp.Circle(a_coil).Face()
torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)
torus.name = "coil"
torus.maxh = maxh_coil

# Interior sphere at origin
inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
inner_sphere.maxh = maxh
for f in inner_sphere.faces:
    f.name = "kelvin_int"

# Subtract coil from interior
inner_air = inner_sphere - torus
inner_air.name = "air"

# Exterior sphere at offset (Kelvin-mapped exterior)
outer_sphere = Sphere(Pnt(offset, 0, 0), R_K)
outer_sphere.maxh = maxh * 2
for f in outer_sphere.faces:
    f.name = "kelvin_ext"
outer_sphere.name = "kelvin"

# GND vertex at exterior center (maps to physical infinity)
gnd = Vertex(Pnt(offset, 0, 0))
gnd.name = "GND"

# Identify periodic faces BEFORE Glue
kelvin_int_face = None
kelvin_ext_face = None
for f in inner_air.faces:
    if f.name == "kelvin_int":
        kelvin_int_face = f
        break
for f in outer_sphere.faces:
    if f.name == "kelvin_ext":
        kelvin_ext_face = f
        break

if kelvin_int_face and kelvin_ext_face:
    kelvin_int_face.Identify(kelvin_ext_face, "periodic",
                              IdentificationType.PERIODIC)
    print("  Periodic identification: kelvin_int <-> kelvin_ext")
else:
    raise RuntimeError("Kelvin faces not found!")

# Combine
shape = Glue([inner_air, torus, outer_sphere, gnd])
geo = OCCGeometry(shape)

with TaskManager():
    ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
mesh = Mesh(ngmesh)
mesh.Curve(order + 1)

t_mesh = time.perf_counter() - t0
mats = mesh.GetMaterials()
bnds = mesh.GetBoundaries()
print(f"  Mesh: {mesh.ne} tets, {mesh.nv} vertices ({t_mesh:.1f}s)")
print(f"  Materials: {mats}")
print(f"  Boundaries: {bnds}")

# ============================================================
# FE Space: HCurl with Periodic Kelvin
# ============================================================
print("[2/4] Setting up HCurl + Periodic...")

# GND: Dirichlet A=0 at mapped infinity
# Note: for HCurl, "dirichlet_bbnd" applies to BBND (edges on boundary)
fes_base = HCurl(mesh, order=order, dirichlet_bbnd="GND")
fes = Periodic(fes_base)

freedof_before = sum([1 for d in fes_base.FreeDofs() if d])
freedof_after = sum([1 for d in fes.FreeDofs() if d])
print(f"  FreeDofs: {freedof_before} -> {freedof_after} "
      f"(coupled: {freedof_before - freedof_after})")
if freedof_before == freedof_after:
    print("  WARNING: Periodic BC NOT working!")
else:
    print("  Periodic BC working (DOFs reduced)")

print(f"  DOFs: {fes.ndof}")
u, v = fes.TnT()

# ============================================================
# Bilinear/Linear forms with Kelvin weight
# ============================================================
print("[3/4] Assembling...")
t0 = time.perf_counter()

# Kelvin weight for 3D curl-curl (A-formulation):
# From Sugahara, IEEE Trans. Magn. 2022:
#   ALL material constants scale by (a/r')^2 (conformal symmetry)
#   mu' = mu0 * (a/r')^2, nu' = nu0 * (a/r')^2, sigma' = sigma0 * (a/r')^2
# Verified: L error +1.2% (coarse) to +1.8% (fine) against torus analytical.
# r' = distance from exterior center (offset, 0, 0)
r_prime_sq = (x - offset)**2 + y**2 + z**2
kelvin_fac = R_K**2 / (r_prime_sq + 1e-20)  # = (a/r')^2

# Build nu per material
nu_dict = {}
for m in mats:
    if "kelvin" in m.lower() or "outer" in m.lower():
        nu_dict[m] = NU_0 * kelvin_fac
    else:
        nu_dict[m] = NU_0
nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

# Source current: J = J0 * e_theta in coil
r_xy = sqrt(x * x + y * y)
r_xy_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
J_source = J0 * CF((-y / r_xy_safe, x / r_xy_safe, 0))

# Bilinear form (bonus_intorder for Kelvin domain accuracy)
a_bf = BilinearForm(fes)
a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
a_bf += 1e-8 * NU_0 * u * v * dx  # gauge regularization
a_bf.Assemble()

# Linear form
f_lf = LinearForm(fes)
f_lf += J_source * v * dx("coil")
f_lf.Assemble()

t_asm = time.perf_counter() - t0
print(f"  Assembly: {t_asm:.1f}s")

# ============================================================
# Solve
# ============================================================
print("[4/4] Solving...")
t0 = time.perf_counter()

gfu = GridFunction(fes)
with TaskManager():
    gfu.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse="pardiso") * f_lf.vec

t_solve = time.perf_counter() - t0

# Energy-based inductance
W = Integrate(0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)), mesh,
              order=10).real
L = 2 * W / I_total**2

print(f"  Solve time: {t_solve:.1f}s")
print()
print("=" * 60)
print("Results")
print("=" * 60)
print(f"  L (HCurl Periodic Kelvin) = {L*1e9:.2f} nH")
print(f"  L (analytical torus)      = {L_analytical*1e9:.2f} nH")
print(f"  Error: {(L/L_analytical - 1)*100:+.2f}%")
print(f"  DOFs: {fes.ndof}, Elements: {mesh.ne}")
print(f"  Time: mesh={t_mesh:.1f}s, solve={t_solve:.1f}s")
print("=" * 60)
