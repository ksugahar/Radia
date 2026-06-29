"""test_delta_L_robustness.py

Robustness tests for Delta_L with workpiece-only RHS + Kelvin.

Test 1: Air-only (no workpiece) -> Delta_L must be zero
Test 2: mu_r scan -> Delta_L must follow dipole chi_eff scaling
Test 3: Distance scan -> Delta_L must follow 1/d^6 dipole scaling
Test 4: h-convergence (repeat, for reference)
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
R_COIL = 0.03
I_TOTAL = 1.0
N_SEG = 48


def make_coil_polyline():
    theta = np.linspace(0, 2 * np.pi, N_SEG + 1)
    pts = np.column_stack([R_COIL * np.cos(theta),
                           R_COIL * np.sin(theta),
                           np.zeros(N_SEG + 1)])
    return [(tuple(pts[i]), tuple(pts[i + 1])) for i in range(N_SEG)]


def delta_L_dipole(R_wp, z_wp, mu_r):
    """Analytical dipole Delta_L."""
    d2 = R_COIL**2 + z_wp**2
    H_z = R_COIL**2 / (2 * d2**1.5)
    V = (4.0 / 3.0) * math.pi * R_wp**3
    chi_eff = 3.0 * (mu_r - 1) / (mu_r + 2)
    m_z = chi_eff * V * H_z
    return (MU_0 / (4 * math.pi)) * m_z * 2 * math.pi * R_COIL**2 / d2**1.5


def make_kelvin_mesh(R_wp, center_wp, R_kelvin, maxh_wp, maxh_air,
                      include_workpiece=True):
    """Kelvin mesh with optional workpiece."""
    from netgen.occ import (Sphere, Pnt, Glue, OCCGeometry,
                             Vertex, IdentificationType)
    from ngsolve import Mesh

    cx, cy, cz = center_wp
    offset_z = 3.0 * R_kelvin

    inner = Sphere(Pnt(0, 0, 0), R_kelvin)
    inner.maxh = maxh_air
    for f in inner.faces:
        f.name = "kelvin_int"

    if include_workpiece:
        wp = Sphere(Pnt(cx, cy, cz), R_wp)
        wp.mat("workpiece")
        wp.maxh = maxh_wp
        inner_air = inner - wp
        inner_air.mat("air")
        parts = [inner_air, wp]
    else:
        inner.mat("air")
        parts = [inner]

    outer = Sphere(Pnt(0, 0, offset_z), R_kelvin)
    outer.maxh = maxh_air
    outer.mat("kelvin")
    for f in outer.faces:
        f.name = "kelvin_ext"

    gnd = Vertex(Pnt(0, 0, offset_z))
    gnd.name = "GND"

    geo = Glue(parts + [outer, gnd])

    kelvin_int_face = kelvin_ext_face = None
    for solid in geo.solids:
        for face in solid.faces:
            if face.name == "kelvin_int":
                kelvin_int_face = face
            elif face.name == "kelvin_ext":
                kelvin_ext_face = face
    kelvin_int_face.Identify(
        kelvin_ext_face, "periodic", IdentificationType.PERIODIC)

    ngmesh = OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5)
    ngmesh.Curve(2)
    return Mesh(ngmesh), R_kelvin, offset_z


def solve_reduced_A_kelvin(R_wp, center_wp, mu_r,
                            maxh_wp=0.002, maxh_air=0.012,
                            R_kelvin=0.06,
                            include_workpiece=True):
    """Solve reduced-A with workpiece-only RHS + Kelvin."""
    from ngsolve import (HCurl, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, curl, dx, TaskManager, VOL)
    from ngsolve import x, y, z
    from kelvin_source import (biot_savart_A_at_points,
                                kelvin_map_3d, kelvin_pullback_vector,
                                is_in_kelvin_exterior_domain,
                                line_integral_A_filaments,
                                compute_back_reaction)

    fil_paths = [make_coil_polyline()]
    fil_currents = [complex(I_TOTAL)]

    mesh, R_K, offset_z = make_kelvin_mesh(
        R_wp, center_wp, R_kelvin, maxh_wp, maxh_air,
        include_workpiece=include_workpiece)
    ne = mesh.GetNE(VOL)
    nv = mesh.nv
    kc = np.array([0.0, 0.0, offset_z])

    # A_s projection (for line integral of A_r later)
    all_pts = np.array([mesh.vertices[v].point for v in range(nv)])
    in_kelv = is_in_kelvin_exterior_domain(all_pts, kc, R_K)
    A_all = np.zeros((nv, 3), dtype=float)
    phys_idx = np.where(~in_kelv)[0]
    if len(phys_idx) > 0:
        A_all[phys_idx] = biot_savart_A_at_points(
            all_pts[phys_idx], fil_paths, fil_currents).real
    kelv_idx = np.where(in_kelv)[0]
    if len(kelv_idx) > 0:
        r_phys = kelvin_map_3d(all_pts[kelv_idx], kc, R_K)
        A_phys = biot_savart_A_at_points(r_phys, fil_paths, fil_currents)
        A_all[kelv_idx] = kelvin_pullback_vector(
            A_phys, all_pts[kelv_idx], kc, R_K).real

    fes_h1 = H1(mesh, order=1, dim=3)
    gf_h1 = GridFunction(fes_h1)
    fv = gf_h1.vec.FV()
    for vi in range(nv):
        for k in range(3):
            fv[vi * 3 + k] = float(A_all[vi, k])
    fes_As = HCurl(mesh, order=1, complex=False)
    gf_As = GridFunction(fes_As)
    gf_As.Set(gf_h1)

    # nu
    materials = mesh.GetMaterials()
    nu_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            dx_k = x - kc[0]
            dy_k = y - kc[1]
            dz_k = z - kc[2]
            rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
            nu_dict[m] = NU_0 * rp_sq / R_K**2
        elif m == "workpiece":
            nu_dict[m] = NU_0 / mu_r
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # FEM
    base_fes = HCurl(mesh, order=1, complex=False, nograds=True,
                      dirichlet="GND")
    fes = Periodic(base_fes)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=2)
    non_kelvin = [m for m in materials if "kelvin" not in m.lower()]
    if non_kelvin:
        a += 1e-6 * NU_0 * u * v * dx("|".join(non_kelvin))
    a.Assemble()

    # RHS: workpiece-only (analytically exact)
    f = LinearForm(fes)
    if include_workpiece and "workpiece" in materials:
        f += -(nu_cf - NU_0) * curl(gf_As) * curl(v) * dx(
            "workpiece", bonus_intorder=4)
    f.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(
            fes.FreeDofs(), inverse="pardiso") * f.vec

    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=6)
    back = compute_back_reaction(fil_currents, delta_phi, I_TOTAL)

    return back['Delta_L'], ne, fes.ndof


def main():
    print("Delta_L Robustness Tests (workpiece-only RHS + Kelvin)")
    print("=" * 60)
    all_pass = True

    # ============================================================
    # Test 1: Air-only (no workpiece) -> Delta_L must be ~0
    # ============================================================
    print("\nTest 1: Air-only (no workpiece)")
    dL_air, ne, ndof = solve_reduced_A_kelvin(
        R_wp=0.01, center_wp=(0, 0, 0.018), mu_r=1.0,
        include_workpiece=False)
    # Compare with a typical workpiece Delta_L (~1500 pH)
    threshold = 1.0e-12  # 1 pH
    ok1 = abs(dL_air) < threshold
    print(f"  Delta_L = {dL_air*1e12:.4f} pH  (threshold: {threshold*1e12:.1f} pH)")
    print(f"  [{'PASS' if ok1 else 'FAIL'}] Delta_L ~ 0 for air-only")
    if not ok1:
        all_pass = False

    # ============================================================
    # Test 2: mu_r scan -> chi_eff scaling
    # ============================================================
    print("\nTest 2: mu_r scan (R_wp=0.01, z_wp=0.018)")
    print(f"  {'mu_r':>8s} {'FEM [pH]':>10s} {'Dipole [pH]':>12s} "
          f"{'ratio':>8s} {'chi_eff':>8s}")
    print(f"  {'-'*52}")

    R_wp = 0.01
    z_wp = 0.018
    center_wp = (0, 0, z_wp)
    mu_r_list = [2, 5, 10, 50, 100, 500, 1000]
    dL_fem_list = []
    dL_ana_list = []

    for mu_r in mu_r_list:
        dL_fem, _, _ = solve_reduced_A_kelvin(
            R_wp, center_wp, mu_r, maxh_wp=0.0015)
        dL_ana = delta_L_dipole(R_wp, z_wp, mu_r)
        chi_eff = 3.0 * (mu_r - 1) / (mu_r + 2)
        ratio = dL_fem / dL_ana if abs(dL_ana) > 1e-30 else float('nan')
        dL_fem_list.append(dL_fem)
        dL_ana_list.append(dL_ana)
        print(f"  {mu_r:8d} {dL_fem*1e12:10.1f} {dL_ana*1e12:12.1f} "
              f"{ratio:8.3f} {chi_eff:8.3f}")

    # Check: FEM/dipole ratio should be roughly constant across mu_r
    ratios = [f / a for f, a in zip(dL_fem_list, dL_ana_list) if abs(a) > 1e-30]
    ratio_spread = max(ratios) - min(ratios)
    ok2 = ratio_spread < 0.15  # within 15% spread
    print(f"\n  Ratio spread: {ratio_spread:.3f} "
          f"({'PASS' if ok2 else 'FAIL'}, <0.15)")
    if not ok2:
        all_pass = False

    # ============================================================
    # Test 3: Distance scan -> 1/d^6 scaling
    # ============================================================
    print("\nTest 3: Distance scan (R_wp=0.01, mu_r=100)")
    print(f"  {'z_wp':>8s} {'FEM [pH]':>10s} {'Dipole [pH]':>12s} "
          f"{'FEM/Dip':>8s}")
    print(f"  {'-'*42}")

    z_list = [0.025, 0.030, 0.040, 0.050]
    for z_wp in z_list:
        dL_fem, _, _ = solve_reduced_A_kelvin(
            0.01, (0, 0, z_wp), 100.0, maxh_wp=0.0015)
        dL_ana = delta_L_dipole(0.01, z_wp, 100.0)
        ratio = dL_fem / dL_ana if abs(dL_ana) > 1e-30 else float('nan')
        print(f"  {z_wp:8.3f} {dL_fem*1e12:10.2f} {dL_ana*1e12:12.2f} "
              f"{ratio:8.3f}")

    # Check dipole scaling: Delta_L ~ 1/(R^2+z^2)^3
    # Ratio of first/last should match dipole scaling
    dL_first, _, _ = solve_reduced_A_kelvin(
        0.01, (0, 0, z_list[0]), 100.0, maxh_wp=0.0015)
    dL_last, _, _ = solve_reduced_A_kelvin(
        0.01, (0, 0, z_list[-1]), 100.0, maxh_wp=0.0015)
    d2_first = R_COIL**2 + z_list[0]**2
    d2_last = R_COIL**2 + z_list[-1]**2
    expected_ratio = (d2_last / d2_first)**3  # 1/d^6 scaling
    actual_ratio = dL_first / dL_last if abs(dL_last) > 1e-30 else float('nan')
    ok3 = abs(actual_ratio / expected_ratio - 1) < 0.10
    print(f"\n  Scaling z={z_list[0]}->z={z_list[-1]}: "
          f"FEM ratio={actual_ratio:.2f}, dipole={expected_ratio:.2f} "
          f"({'PASS' if ok3 else 'FAIL'}, <10%)")
    if not ok3:
        all_pass = False

    # ============================================================
    # Summary
    # ============================================================
    print(f"\n{'='*60}")
    checks = [
        ("Air-only Delta_L ~ 0", ok1),
        ("mu_r scaling consistent", ok2),
        ("Distance 1/d^6 scaling", ok3),
    ]
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    if all_pass:
        print("\nAll robustness tests PASSED.")
    else:
        print("\nSome tests FAILED!")
        sys.exit(1)


if __name__ == "__main__":
    main()
