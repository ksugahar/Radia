"""convergence_delta_L.py

Mesh convergence study for Delta_L (ferrite sphere, mu_r=100).

Uses Kelvin transformation (NO domain truncation error).
The only remaining error sources are h-refinement and p-refinement.

Reference: dipole analytical Delta_L (l=1 only).
True Delta_L > dipole because higher multipoles (l>=2) contribute positively.
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
R_WP = 0.01
CENTER_WP = (0, 0, 0.018)
I_TOTAL = 1.0
N_SEG = 48
MU_R = 100.0


def make_coil_polyline():
    theta = np.linspace(0, 2 * np.pi, N_SEG + 1)
    pts = np.column_stack([R_COIL * np.cos(theta),
                           R_COIL * np.sin(theta),
                           np.zeros(N_SEG + 1)])
    return [(tuple(pts[i]), tuple(pts[i + 1])) for i in range(N_SEG)]


def delta_L_analytical():
    """Dipole approximation (l=1 term only, underestimates true value)."""
    h = CENTER_WP[2]
    d2 = R_COIL**2 + h**2
    H_z = R_COIL**2 / (2 * d2**1.5)
    V = (4.0 / 3.0) * math.pi * R_WP**3
    chi_eff = 3.0 * (MU_R - 1) / (MU_R + 2)
    m_z = chi_eff * V * H_z
    Phi = (MU_0 / (4 * math.pi)) * m_z * 2 * math.pi * R_COIL**2 / d2**1.5
    return Phi


def make_kelvin_mesh(maxh_wp, maxh_air, R_kelvin, fes_order=1):
    """Build Kelvin-transformed mesh: workpiece + air + Kelvin exterior.

    Two-sphere model (Sugahara 2022):
      inner sphere (R=R_kelvin): air + workpiece
      outer sphere (R=R_kelvin, offset): Kelvin exterior domain
      Periodic BC between inner/outer sphere surfaces
      GND vertex at outer sphere center (= physical infinity)
    """
    from netgen.occ import (Sphere, Pnt, Glue, OCCGeometry,
                             Vertex, IdentificationType)
    from ngsolve import Mesh

    cx, cy, cz = CENTER_WP
    offset_z = 3.0 * R_kelvin  # standard offset

    # Workpiece sphere
    wp = Sphere(Pnt(cx, cy, cz), R_WP)
    wp.mat("workpiece")
    wp.maxh = maxh_wp

    # Inner physical sphere
    inner = Sphere(Pnt(0, 0, 0), R_kelvin)
    inner.maxh = maxh_air
    for f in inner.faces:
        f.name = "kelvin_int"

    inner_air = inner - wp
    inner_air.mat("air")

    # Outer Kelvin sphere (offset in z)
    outer = Sphere(Pnt(0, 0, offset_z), R_kelvin)
    outer.maxh = maxh_air
    outer.mat("kelvin")
    for f in outer.faces:
        f.name = "kelvin_ext"

    # GND vertex at exterior center (r'=0 = physical infinity)
    gnd = Vertex(Pnt(0, 0, offset_z))
    gnd.name = "GND"

    geo = Glue([inner_air, wp, outer, gnd])

    # Periodic identification
    kelvin_int_face = None
    kelvin_ext_face = None
    for solid in geo.solids:
        for face in solid.faces:
            if face.name == "kelvin_int":
                kelvin_int_face = face
            elif face.name == "kelvin_ext":
                kelvin_ext_face = face

    if kelvin_int_face is None or kelvin_ext_face is None:
        raise RuntimeError("Could not find kelvin_int/kelvin_ext faces")

    kelvin_int_face.Identify(
        kelvin_ext_face, "periodic", IdentificationType.PERIODIC)

    ngmesh = OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5)
    ngmesh.Curve(max(2, fes_order))
    return Mesh(ngmesh), R_kelvin, offset_z


def run_kelvin(maxh_wp, maxh_air, R_kelvin, fes_order=1):
    """Reduced-A + line integral with Kelvin transformation."""
    from ngsolve import (HCurl, H1, Periodic, BilinearForm, LinearForm,
                         GridFunction, curl, dx, TaskManager, VOL)
    from ngsolve import x, y, z
    from kelvin_source import (biot_savart_A_at_points,
                                B_s_at_obs_with_kelvin,
                                kelvin_map_3d, kelvin_pullback_vector,
                                is_in_kelvin_exterior_domain,
                                line_integral_A_filaments,
                                compute_back_reaction)

    fil_paths = [make_coil_polyline()]
    fil_currents = [complex(I_TOTAL)]

    mesh, R_K, offset_z = make_kelvin_mesh(
        maxh_wp, maxh_air, R_kelvin, fes_order)
    ne = mesh.GetNE(VOL)
    nv = mesh.nv
    kc = np.array([0.0, 0.0, offset_z])

    # --- A_s for line integral (HCurl, needed later) ---
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
    gf_h1_A = GridFunction(fes_h1)
    fv = gf_h1_A.vec.FV()
    for vi in range(nv):
        for k in range(3):
            fv[vi * 3 + k] = float(A_all[vi, k])
    fes_As = HCurl(mesh, order=fes_order, complex=False)
    gf_As = GridFunction(fes_As)
    gf_As.Set(gf_h1_A)

    # --- B_s = curl(A_s) directly projected (NOT via curl of HCurl) ---
    # This avoids the piecewise-constant curl error that corrupts the
    # Kelvin domain RHS.  B_s is evaluated analytically at vertices
    # (Biot-Savart + 2-form pullback) and projected into H1(dim=3).
    B_all = B_s_at_obs_with_kelvin(
        all_pts, fil_paths, fil_currents,
        kelvin_center=kc, R_kelvin=R_K)
    gf_h1_B = GridFunction(fes_h1)
    fv_B = gf_h1_B.vec.FV()
    for vi in range(nv):
        for k in range(3):
            fv_B[vi * 3 + k] = float(B_all[vi, k])

    # nu with Kelvin modulation
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
            nu_dict[m] = NU_0 / MU_R
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # FEM: Periodic HCurl with GND Dirichlet
    base_fes = HCurl(mesh, order=fes_order, complex=False,
                      nograds=True, dirichlet="GND")
    fes = Periodic(base_fes)
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=2 * fes_order)
    non_kelvin = [m for m in materials if "kelvin" not in m.lower()]
    if non_kelvin:
        a += 1e-6 * NU_0 * u * v * dx("|".join(non_kelvin))
    a.Assemble()

    # RHS: -int nu_cf * B_s . curl(v) dV   (correct, full domain)
    # B_s is directly projected — no curl(HCurl) degradation.
    f = LinearForm(fes)
    f += -nu_cf * gf_h1_B * curl(v) * dx(bonus_intorder=2 * fes_order + 2)
    f.Assemble()

    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(
            fes.FreeDofs(), inverse="pardiso") * f.vec

    # Line integral (filaments are in the physical inner domain)
    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=6)
    back = compute_back_reaction(fil_currents, delta_phi, I_TOTAL)

    return {
        'Delta_L': back['Delta_L'],
        'ne': ne,
        'ndof': fes.ndof,
    }


def main():
    dL_ana = delta_L_analytical()
    print(f"Analytical (dipole, l=1 only): {dL_ana*1e12:.1f} pH")
    print(f"  True value > dipole (higher multipoles add positively)")
    print()

    R_K = 0.06  # Kelvin radius (must enclose coil + workpiece)

    # ---- h-refinement ----
    print(f"h-refinement (Kelvin R={R_K}, order=1)")
    print(f"  {'maxh_wp':>8s} {'ne':>8s} {'ndof':>8s} "
          f"{'Delta_L[pH]':>12s} {'vs dipole':>10s}")
    print(f"  {'-'*55}")

    prev_dL = None
    for maxh_wp in [0.005, 0.003, 0.002, 0.0015, 0.001]:
        maxh_air = max(maxh_wp * 4, 0.012)
        t0 = time.perf_counter()
        r = run_kelvin(maxh_wp, maxh_air, R_K, fes_order=1)
        t = time.perf_counter() - t0
        err = (r['Delta_L'] - dL_ana) / dL_ana * 100
        rate = ""
        if prev_dL is not None and r['Delta_L'] != prev_dL:
            rate = f"  d={abs(r['Delta_L'] - prev_dL)*1e12:.0f}pH"
        prev_dL = r['Delta_L']
        print(f"  {maxh_wp:8.4f} {r['ne']:8d} {r['ndof']:8d} "
              f"{r['Delta_L']*1e12:12.1f} {err:+9.1f}%  "
              f"({t:.1f}s){rate}")

    # ---- p-refinement ----
    print(f"\np-refinement (Kelvin R={R_K}, maxh_wp=0.002)")
    print(f"  {'order':>8s} {'ne':>8s} {'ndof':>8s} "
          f"{'Delta_L[pH]':>12s} {'vs dipole':>10s}")
    print(f"  {'-'*55}")

    for order in [1, 2]:
        t0 = time.perf_counter()
        r = run_kelvin(0.002, 0.012, R_K, fes_order=order)
        t = time.perf_counter() - t0
        err = (r['Delta_L'] - dL_ana) / dL_ana * 100
        print(f"  {order:8d} {r['ne']:8d} {r['ndof']:8d} "
              f"{r['Delta_L']*1e12:12.1f} {err:+9.1f}%  ({t:.1f}s)")

    print(f"\n  Dipole is l=1 ONLY. True Delta_L includes l>=2.")
    print(f"  Converged FEM exceeding dipole = higher multipoles captured.")


if __name__ == "__main__":
    main()
