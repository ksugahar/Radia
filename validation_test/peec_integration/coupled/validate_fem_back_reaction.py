"""validate_fem_back_reaction.py

Quantitative validation of PEEC+FEM back-reaction (Delta_L).

Three-way comparison for the ferrite case (mu_r=100, sigma=0):

  Method A: Direct FEM (coil torus + workpiece sphere + air)
    L_total_A = 2 * int nu |curl(A)|^2 / I^2
    L_air_A   = same mesh without workpiece
    Delta_L_A = L_total_A - L_air_A

  Method B: Reduced-A energy integral
    A_s = Biot-Savart (analytical); solve for A_r
    L_total_B = 2 * int nu |curl(A_s + A_r)|^2 / I^2
    L_air_B   = 2 * int nu_0 |curl(A_s)|^2 / I^2
    Delta_L_B = L_total_B - L_air_B

  Method C: Reduced-A line integral (new code under test)
    Delta_L_C = int A_r . dl / I  along coil filaments

If all three agree, the line integral implementation is validated.

Sign checks:
    Cu  SIBC (sigma=5.8e7):  Delta_L < 0  (Lenz)
    Ferrite   (mu_r=100):     Delta_L > 0  (flux concentration)
"""

import math
import os
import sys
import time

import numpy as np

_repo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
_src = os.path.join(_repo, 'src', 'radia')
if _src not in sys.path:
    sys.path.insert(0, _src)
_panels = os.path.join(_src, 'panels')
if _panels not in sys.path:
    sys.path.insert(0, _panels)


MU_0 = 4 * math.pi * 1e-7
NU_0 = 1.0 / MU_0

# Shared geometry
R_COIL = 0.03       # coil major radius [m]
A_WIRE = 0.0015     # coil wire radius [m] (1.5 mm)
Z_COIL = 0.0        # coil z-position
R_WP = 0.01         # workpiece sphere radius [m]
CENTER_WP = (0, 0, 0.018)
R_AIR = 0.10
I_TOTAL = 1.0
N_SEG = 48           # polyline segments for Biot-Savart


def make_coil_polyline():
    """Circular coil centerline as segment tuples."""
    theta = np.linspace(0, 2 * np.pi, N_SEG + 1)
    pts = np.column_stack([R_COIL * np.cos(theta),
                           R_COIL * np.sin(theta),
                           np.full(N_SEG + 1, Z_COIL)])
    return [(tuple(pts[i]), tuple(pts[i + 1])) for i in range(N_SEG)]


# ============================================================
#  Mesh builders
# ============================================================

def _make_direct_fem_mesh(include_workpiece, mu_r_wp=1.0,
                           maxh_coil=0.0012, maxh_wp=0.004,
                           maxh_air=0.025):
    """OCC mesh with coil torus + optional workpiece sphere + air."""
    from netgen.occ import (Sphere, Pnt, Glue, OCCGeometry,
                             WorkPlane, Axes, X, Y, Z, Axis)
    from ngsolve import Mesh

    # Coil torus: revolve circle cross-section around Z axis
    face = WorkPlane(Axes(Pnt(R_COIL, 0, Z_COIL), Y, X)).Circle(A_WIRE).Face()
    coil = face.Revolve(Axis(Pnt(0, 0, 0), Z), 360)
    coil.mat("coil")
    coil.maxh = maxh_coil

    # Air sphere
    air_full = Sphere(Pnt(0, 0, 0), R_AIR)
    air_full.bc("outer")

    if include_workpiece:
        cx, cy, cz = CENTER_WP
        wp = Sphere(Pnt(cx, cy, cz), R_WP)
        wp.mat("workpiece")
        wp.maxh = maxh_wp
        air = air_full - coil - wp
        air.mat("air")
        geo = OCCGeometry(Glue([air, coil, wp]))
    else:
        air = air_full - coil
        air.mat("air")
        geo = OCCGeometry(Glue([air, coil]))

    ngmesh = geo.GenerateMesh(maxh=maxh_air)
    ngmesh.Curve(2)
    return Mesh(ngmesh)


def make_hole_mesh(maxh_wp=0.004, maxh_air=0.025):
    """SIBC hole mesh (no coil volume, no workpiece volume)."""
    from netgen.occ import Sphere, Pnt, OCCGeometry
    from ngsolve import Mesh

    cx, cy, cz = CENTER_WP
    wp = Sphere(Pnt(cx, cy, cz), R_WP)
    wp.bc("sibc")
    for f in wp.faces:
        f.maxh = maxh_wp

    air = Sphere(Pnt(0, 0, 0), R_AIR)
    air.bc("outer")
    shape = air - wp
    shape.mat("air")

    geo = OCCGeometry(shape)
    ngmesh = geo.GenerateMesh(maxh=maxh_air)
    ngmesh.Curve(2)
    return Mesh(ngmesh)


def make_volume_mesh(maxh_wp=0.004, maxh_air=0.025):
    """Workpiece volume mesh (no coil volume)."""
    from netgen.occ import Sphere, Pnt, Glue, OCCGeometry
    from ngsolve import Mesh

    cx, cy, cz = CENTER_WP
    wp = Sphere(Pnt(cx, cy, cz), R_WP)
    wp.mat("workpiece")
    wp.maxh = maxh_wp

    air_full = Sphere(Pnt(0, 0, 0), R_AIR)
    air_full.bc("outer")
    air = air_full - wp
    air.mat("air")

    geo = OCCGeometry(Glue([air, wp]))
    ngmesh = geo.GenerateMesh(maxh=maxh_air)
    ngmesh.Curve(2)
    return Mesh(ngmesh)


# ============================================================
#  Helpers
# ============================================================

def project_A_s(mesh, fil_paths, fil_currents):
    """Biot-Savart A_s -> HCurl GridFunction (real-only thin wrapper)."""
    from kelvin_source import project_A_s_to_hcurl
    return project_A_s_to_hcurl(mesh, fil_paths, fil_currents,
                                is_complex=False)


def L_energy(mesh, B_cf, nu_cf, I):
    """L from volume magnetic energy."""
    from ngsolve import Integrate, Conj
    W = Integrate(0.5 * nu_cf * B_cf * Conj(B_cf), mesh, order=6).real
    return 2 * W / I**2


# ============================================================
#  Method A: Direct FEM A-formulation (coil volume mesh)
# ============================================================

def solve_direct_fem(include_workpiece, mu_r_wp=1.0):
    """Direct FEM: curl(nu * curl(A)) = J in coil volume.

    Returns L from energy integral.
    """
    from ngsolve import (HCurl, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, dx, CF,
                         TaskManager, VOL)
    from ngsolve import x, y, sqrt

    mesh = _make_direct_fem_mesh(include_workpiece, mu_r_wp)
    ne = mesh.GetNE(VOL)

    # Material nu
    mats = set(mesh.GetMaterials())
    nu_dict = {}
    for m in mats:
        if m == "workpiece":
            nu_dict[m] = NU_0 / mu_r_wp
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Current density: J = I / A_cross * e_theta in coil
    A_cross = math.pi * A_WIRE**2
    J_mag = I_TOTAL / A_cross
    r_xy = sqrt(x * x + y * y + 1e-30)
    J_cf = J_mag * CF((-y / r_xy, x / r_xy, 0))

    fes = HCurl(mesh, order=1, complex=False, nograds=True, dirichlet="outer")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=2)
    a += 1e-6 * NU_0 * u * v * dx
    a.Assemble()

    f = LinearForm(fes)
    f += J_cf * v * dx("coil")
    f.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

    L = L_energy(mesh, curl(gfu), nu_cf, I_TOTAL)
    return L, ne, fes.ndof


# ============================================================
#  Method B+C: Reduced-A (Biot-Savart source)
# ============================================================

def solve_reduced_A_ferrite(mu_r_wp):
    """Reduced-A for ferrite: returns (gfu, gf_As, mesh, nu_cf)."""
    from ngsolve import (HCurl, BilinearForm, LinearForm,
                         GridFunction, curl, dx, TaskManager, VOL)

    fil_paths = [make_coil_polyline()]
    fil_currents = [complex(I_TOTAL)]

    mesh = make_volume_mesh()
    ne = mesh.GetNE(VOL)

    gf_As = project_A_s(mesh, fil_paths, fil_currents)

    nu_wp = NU_0 / mu_r_wp
    nu_cf = mesh.MaterialCF({"air": NU_0, "workpiece": nu_wp}, default=NU_0)

    fes = HCurl(mesh, order=1, complex=False, nograds=True, dirichlet="outer")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=2)
    a += 1e-6 * NU_0 * u * v * dx
    a.Assemble()

    f = LinearForm(fes)
    f += -(nu_cf - NU_0) * curl(gf_As) * curl(v) * dx(bonus_intorder=4)
    f.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

    return gfu, gf_As, mesh, nu_cf, ne, fes.ndof, fil_paths, fil_currents


# ============================================================
#  Case 1: Ferrite -- three-way comparison
# ============================================================

def run_ferrite(mu_r_wp=100.0):
    from ngsolve import curl, CF
    from kelvin_source import (line_integral_A_filaments,
                                compute_back_reaction)

    print(f"\n{'='*60}")
    print(f"  FERRITE: mu_r={mu_r_wp}, sigma=0")
    print(f"  Three-way: Direct FEM / Reduced-A energy / Line integral")
    print(f"{'='*60}")

    # --- Method A: Direct FEM (independent reference) ---
    print("\n  [A] Direct FEM (coil torus + workpiece sphere + air)")
    t0 = time.perf_counter()
    L_total_A, ne_A, ndof_A = solve_direct_fem(True, mu_r_wp)
    t_A = time.perf_counter() - t0
    print(f"      ne={ne_A}, ndof={ndof_A}, t={t_A:.1f}s")
    print(f"      L_total = {L_total_A*1e9:.4f} nH")

    print("  [A] Direct FEM (coil torus + air, no workpiece)")
    t0 = time.perf_counter()
    L_air_A, ne_A0, ndof_A0 = solve_direct_fem(False)
    t_A0 = time.perf_counter() - t0
    print(f"      ne={ne_A0}, ndof={ndof_A0}, t={t_A0:.1f}s")
    print(f"      L_air   = {L_air_A*1e9:.4f} nH")

    Delta_L_A = L_total_A - L_air_A
    print(f"      Delta_L = {Delta_L_A*1e9:.4f} nH")

    # --- Methods B+C: Reduced-A ---
    print("\n  [B+C] Reduced-A (Biot-Savart source)")
    t0 = time.perf_counter()
    (gfu, gf_As, mesh, nu_cf,
     ne_BC, ndof_BC, fil_paths, fil_currents) = solve_reduced_A_ferrite(mu_r_wp)
    t_BC = time.perf_counter() - t0
    print(f"      ne={ne_BC}, ndof={ndof_BC}, t={t_BC:.1f}s")

    # Method B: energy integral
    B_total = curl(gfu) + curl(gf_As)
    B_air = curl(gf_As)
    L_total_B = L_energy(mesh, B_total, nu_cf, I_TOTAL)
    L_air_B = L_energy(mesh, B_air, CF(NU_0), I_TOTAL)
    Delta_L_B = L_total_B - L_air_B
    print(f"  [B] Energy:  L_total={L_total_B*1e9:.4f}, "
          f"L_air={L_air_B*1e9:.4f}, Delta_L={Delta_L_B*1e9:.4f} nH")

    # Method C: line integral
    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=4)
    back = compute_back_reaction(fil_currents, delta_phi, I_TOTAL)
    Delta_L_C = back['Delta_L']
    print(f"  [C] Line:    Delta_L={Delta_L_C*1e9:.4f} nH")

    # --- Comparison table ---
    print(f"\n  {'Method':25s} {'Delta_L [nH]':>14s}  {'vs A':>10s}")
    print(f"  {'-'*55}")
    print(f"  {'A: Direct FEM':25s} {Delta_L_A*1e9:14.4f}  {'(ref)':>10s}")

    err_B = (Delta_L_B - Delta_L_A) / abs(Delta_L_A) * 100 if abs(Delta_L_A) > 0 else float('nan')
    print(f"  {'B: Reduced-A energy':25s} {Delta_L_B*1e9:14.4f}  {err_B:+9.1f}%")

    err_C = (Delta_L_C - Delta_L_A) / abs(Delta_L_A) * 100 if abs(Delta_L_A) > 0 else float('nan')
    print(f"  {'C: Line integral':25s} {Delta_L_C*1e9:14.4f}  {err_C:+9.1f}%")

    ok_sign = Delta_L_C > 0
    ok_quant_BC = abs(err_B - err_C) < 5  # B and C should nearly match
    ok_quant_AC = abs(err_C) < 30  # A vs C within 30% (different meshes)

    print(f"\n  Sign check (C > 0):      {'PASS' if ok_sign else 'FAIL'}")
    print(f"  B vs C (same mesh):      {abs(err_B - err_C):.1f}% "
          f"({'PASS' if ok_quant_BC else 'FAIL'}, <5%)")
    print(f"  A vs C (cross-mesh):     {abs(err_C):.1f}% "
          f"({'PASS' if ok_quant_AC else 'FAIL'}, <30%)")

    return ok_sign, ok_quant_BC, ok_quant_AC, Delta_L_A, Delta_L_B, Delta_L_C


# ============================================================
#  Case 2: Cu SIBC -- sign check
# ============================================================

def run_copper(frequency=50000, sigma_wp=5.8e7):
    from ngsolve import (HCurl, BilinearForm, LinearForm,
                         GridFunction, curl, dx, ds, CF,
                         TaskManager, VOL)
    from kelvin_source import (line_integral_A_filaments,
                                compute_back_reaction)

    print(f"\n{'='*60}")
    print(f"  CU SIBC: sigma={sigma_wp:.2e}, f={frequency} Hz")
    print(f"{'='*60}")

    fil_paths = [make_coil_polyline()]
    fil_currents = [complex(I_TOTAL)]

    omega = 2 * math.pi * frequency
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma_wp))
    Z_s = (1 + 1j) / (sigma_wp * delta)
    robin = 1j * omega / Z_s

    mesh = make_hole_mesh()
    ne = mesh.GetNE(VOL)
    print(f"  Mesh: {ne} elems, delta={delta*1e3:.3f} mm")

    gf_As = project_A_s(mesh, fil_paths, fil_currents)

    fes = HCurl(mesh, order=1, complex=True, nograds=True, dirichlet="outer")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += NU_0 * curl(u) * curl(v) * dx(bonus_intorder=2)
    a += 1e-6 * NU_0 * u * v * dx
    a += robin * u.Trace() * v.Trace() * ds("sibc")
    a.Assemble()

    f = LinearForm(fes)
    f += -robin * gf_As * v.Trace() * ds("sibc")
    f.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=4)
    back = compute_back_reaction(fil_currents, delta_phi, I_TOTAL)
    Delta_L = back['Delta_L']
    Delta_R_ow = back['Delta_R_over_omega']

    print(f"  Delta_L      = {Delta_L*1e9:.4f} nH")
    print(f"  Delta_R/w    = {Delta_R_ow:.4e}")

    ok_sign = Delta_L < 0
    print(f"  Sign check:  {'PASS' if ok_sign else 'FAIL'} (< 0 Lenz)")

    return ok_sign, Delta_L


# ============================================================
#  Main
# ============================================================

# ============================================================
#  Analytical reference: magnetizable sphere in loop field
# ============================================================

def delta_L_analytical_dipole(R_coil, z_wp, R_wp, mu_r):
    """Analytical Delta_L for a permeable sphere near a circular loop.

    Dipole approximation: the sphere sees a nearly uniform field H_0 at
    its center, develops an induced magnetic moment

        m = 3*(mu_r - 1)/(mu_r + 2) * V * H_0

    and the mutual flux linkage back through the loop is computed from
    the dipole A-field.

    This OVERESTIMATES Delta_L because the real field is non-uniform
    (higher multipoles reduce coupling). Expected accuracy: ~10-20%
    for R_wp / distance ~ 0.5.
    """
    # On-axis H from a circular current loop (I=1A) at height z_wp
    d2 = R_coil**2 + z_wp**2
    H_z = R_coil**2 / (2 * d2**1.5)  # A/m for I=1A

    # Induced dipole moment
    V = (4.0 / 3.0) * math.pi * R_wp**3
    chi_eff = 3.0 * (mu_r - 1) / (mu_r + 2)
    m_z = chi_eff * V * H_z  # A*m^2

    # Flux from dipole m at (0,0,z_wp) through loop at z=0:
    # Phi = mu_0/(4*pi) * m_z * 2*pi*R^2 / (R^2 + h^2)^(3/2)
    Phi = (MU_0 / (4 * math.pi)) * m_z * 2 * math.pi * R_coil**2 / d2**1.5

    # Delta_L = Phi / I (I=1A)
    return Phi


def main():
    print("PEEC+FEM Back-Reaction: Quantitative Validation")
    print("=" * 60)

    # Analytical reference (dipole approximation)
    dL_ana = delta_L_analytical_dipole(
        R_COIL, CENTER_WP[2], R_WP, mu_r=100.0)
    print(f"Analytical (dipole approx): Delta_L = {dL_ana*1e9:.4f} nH")
    print(f"  (overestimates by ~10-20%% due to field non-uniformity)")

    # Ferrite: four-way comparison (analytical + A + B + C)
    (ok_sign_f, ok_BC_f, ok_AC_f,
     dL_A, dL_B, dL_C) = run_ferrite(mu_r_wp=100.0)

    # Cu SIBC: sign check
    ok_sign_cu, dL_cu = run_copper()

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    # vs analytical (dipole overestimates, so FEM < analytical)
    err_A_ana = (dL_A - dL_ana) / dL_ana * 100
    err_C_ana = (dL_C - dL_ana) / dL_ana * 100

    print(f"\n  Ferrite Delta_L comparison:")
    print(f"  {'Method':30s} {'Delta_L [pH]':>14s}  {'vs Dipole':>10s}")
    print(f"  {'-'*60}")
    print(f"  {'Analytical (dipole)':30s} {dL_ana*1e12:14.1f}  {'(ref)':>10s}")
    print(f"  {'A: Direct FEM':30s} {dL_A*1e12:14.1f}  {err_A_ana:+9.1f}%")
    print(f"  {'B: Reduced-A energy':30s} {dL_B*1e12:14.1f}")
    print(f"  {'C: Line integral':30s} {dL_C*1e12:14.1f}  {err_C_ana:+9.1f}%")
    print(f"\n  Cu Delta_L (line integral): {dL_cu*1e12:+.1f} pH")

    checks = [
        ("Ferrite sign (>0)", ok_sign_f),
        ("Ferrite B~C same-mesh (<5%%)", ok_BC_f),
        ("Ferrite A~C cross-mesh (<30%%)", ok_AC_f),
        ("Cu sign (<0)", ok_sign_cu),
        ("A vs dipole (-10~-20%%)",
         -25 < err_A_ana < 0),  # FEM should be below dipole
        ("C vs dipole (same range as A)",
         abs(err_C_ana - err_A_ana) < 10),
    ]

    all_pass = True
    for name, ok in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}")
        if not ok:
            all_pass = False

    if all_pass:
        print("\nAll checks PASSED.")
    else:
        print("\nSome checks FAILED!")
        sys.exit(1)


if __name__ == "__main__":
    main()
