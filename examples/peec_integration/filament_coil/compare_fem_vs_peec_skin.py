"""compare_fem_vs_peec_skin.py

2-way IH inductance + heating comparison, same mesh and same SIBC workpiece:

  (1) A-V: full 3D FEM, coil volume meshed with coil_sigma
      - Coil skin/proximity eddies captured by the FEM itself
      - L_av = 2 * int nu |curl A|^2 / I^2
      - P_av = 0.5 * Re(Z_s) * H_t_rms^2 * A_wp
      - Heavy (ndof ~ 280k for this sample), but 3D reference truth

  (2) PEEC + back-reaction: PEEC circuit (L_peec) + FEM scattered A_r
      - PEEC filaments from STEP give Biot-Savart excitation A_s
      - FEM solves scattered A_r with SIBC Robin on workpiece surface
      - L_total = L_peec + Delta_L (line-integral of A_r along filaments)
      - P_wp    = 0.5 * Delta_R * |I|^2

Both paths use SIBC for the workpiece.  3D without a 2D-axisym shortcut
needs SIBC to avoid meshing the skin layer.
"""
import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
for p in (SRC, os.path.join(SRC, 'radia'), os.path.join(SRC, 'radia', 'panels'), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

SAMPLES = os.path.join(SRC, 'radia', 'panels', 'samples')

MU_0 = 4 * math.pi * 1e-7
NU_0 = 1.0 / MU_0


def fmt_pct(val, ref):
    if ref == 0 or ref is None:
        return "n/a"
    return f"{(val - ref) / abs(ref) * 100:+.2f}%"


def solve_peec_back_reaction(vol_file, step, frequency, mat,
                             half_thickness, nwinc, nhinc,
                             peec_sigma=5.8e7, I_total=1.0):
    """PEEC L_peec + FEM-SIBC scattered A_r for back-reaction Delta_L.

    Returns dict with L_peec, Delta_L, L_total, Delta_R, P_wp.
    """
    from ngsolve import (Mesh, HCurl, BilinearForm, LinearForm,
                         GridFunction, Periodic, curl, dx, ds,
                         VOL, x, y, z, TaskManager)
    from coil_from_cad import filaments_from_step
    from calc_fem_kelvin import detect_kelvin_offset
    from kelvin_source import (project_A_s_to_hcurl,
                               line_integral_A_filaments,
                               compute_back_reaction)

    # --- 1. PEEC topology + port impedance + per-filament currents ---
    t0 = time.perf_counter()
    topo = filaments_from_step(step, sigma=peec_sigma,
                               nwinc=nwinc, nhinc=nhinc)
    solver = topo['solver']
    Z_port = solver.compute_port_impedance(frequency)
    L_peec = float(Z_port.imag) / (2 * math.pi * frequency)
    R_peec = float(Z_port.real)

    fil_paths = topo['filament_paths']
    seg_of_fil = topo['seg_of_filament']
    I_branch = solver.compute_branch_currents(frequency, [I_total])
    fil_currents = np.array([
        complex(np.mean([I_branch[s] for s in segs]))
        for segs in seg_of_fil
    ])
    print(f"  PEEC: n_fil={len(fil_paths)}, L_peec={L_peec*1e9:.3f} nH, "
          f"sum|I|={abs(fil_currents.sum()):.4f} A, "
          f"t={time.perf_counter()-t0:.1f}s")

    # --- 2. Mesh + Kelvin detection + SIBC boundary ---
    t0 = time.perf_counter()
    mesh = Mesh(vol_file)
    materials = mesh.GetMaterials()
    boundaries = set(mesh.GetBoundaries())

    has_kelvin = 'kelvin' in materials
    has_sibc = 'sibc' in boundaries
    if not has_sibc:
        raise RuntimeError("Mesh has no 'sibc' boundary — need SIBC workpiece")

    if has_kelvin:
        kelvin_center = np.array(detect_kelvin_offset(mesh))
        kelvin_verts = set()
        for el in mesh.Elements(VOL):
            if el.mat == 'kelvin':
                for v in el.vertices:
                    kelvin_verts.add(v.nr)
        coords = np.array([mesh.vertices[v].point for v in kelvin_verts])
        dists = np.linalg.norm(coords - kelvin_center[None, :], axis=1)
        a_kelvin = float(np.max(dists))
        n_ident = mesh.ngmesh.GetNrIdentifications()
        has_kelvin_periodic = n_ident > 0
    else:
        kelvin_center = np.array([0.0, 0.0, 0.0])
        a_kelvin = 0.0
        has_kelvin_periodic = False
    print(f"  MESH: ne={mesh.GetNE(VOL)}, has_kelvin={has_kelvin}, "
          f"periodic={has_kelvin_periodic}, t={time.perf_counter()-t0:.1f}s")

    # --- 3. nu_cf: coil treated as air (PEEC owns the coil self-L) ---
    kx, ky, kz = kelvin_center
    nu_dict = {}
    for m in materials:
        if 'kelvin' in m.lower():
            dx_k = x - kx
            dy_k = y - ky
            dz_k = z - kz
            rp_sq = dx_k * dx_k + dy_k * dy_k + dz_k * dz_k + 1e-20
            nu_dict[m] = NU_0 * rp_sq / a_kelvin ** 2
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # --- 4. Project A_s (filament Biot-Savart, Kelvin pullback) to HCurl ---
    t0 = time.perf_counter()
    gf_As = project_A_s_to_hcurl(
        mesh, fil_paths, fil_currents,
        kelvin_center=kelvin_center if has_kelvin else None,
        R_kelvin=a_kelvin if has_kelvin else None,
        factor_mode='pullback', is_complex=True)
    print(f"  A_s projection: t={time.perf_counter()-t0:.1f}s")

    # --- 5. FES + SIBC Robin system for scattered A_r ---
    Z_s = mat.dowell_Zs(frequency, half_thickness)
    omega = 2 * math.pi * frequency
    robin = 1j * omega / Z_s
    print(f"  Z_s={Z_s:.4e}, Robin={robin:.4e}")

    dirichlet_bnd = 'GND' if 'GND' in boundaries else ''
    base_fes = HCurl(mesh, order=1, nograds=True, complex=True,
                     dirichlet=dirichlet_bnd)
    fes = Periodic(base_fes) if has_kelvin_periodic else base_fes

    u, v_ = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v_) * dx(bonus_intorder=4)

    non_kelvin_mats = [m for m in materials if 'kelvin' not in m.lower()]
    if non_kelvin_mats:
        a += 1e-6 * NU_0 * u * v_ * dx('|'.join(non_kelvin_mats))

    a += robin * u.Trace() * v_.Trace() * ds('sibc')

    f = LinearForm(fes)
    f += -robin * gf_As * v_.Trace() * ds('sibc')

    t0 = time.perf_counter()
    a.Assemble()
    f.Assemble()
    t_asm = time.perf_counter() - t0

    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse='pardiso') * f.vec
    t_solve = time.perf_counter() - t0
    print(f"  FEM: ndof={fes.ndof}, t_asm={t_asm:.1f}s, t_solve={t_solve:.1f}s")

    # --- 6. Back-reaction ΔL, ΔR via line integral of A_r along filaments ---
    delta_phi = line_integral_A_filaments(gfu, mesh, fil_paths, n_quad=4)
    back = compute_back_reaction(fil_currents, delta_phi, I_total)
    Delta_L = back['Delta_L']
    Delta_R = omega * back['Delta_R_over_omega']
    P_wp = 0.5 * Delta_R * abs(I_total) ** 2

    return {
        'L_peec': L_peec,
        'R_peec': R_peec,
        'Delta_L': Delta_L,
        'L_total': L_peec + Delta_L,
        'Delta_R': Delta_R,
        'P_wp': P_wp,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vol", default=os.path.join(SAMPLES, "ih_fem_kelvin_skin.vol"))
    ap.add_argument("--step", default=os.path.join(SAMPLES, "ih_fem_kelvin_skin_coil.step"))
    ap.add_argument("--frequency", type=float, default=7000)
    ap.add_argument("--peec-nwinc", type=int, default=3)
    ap.add_argument("--peec-nhinc", type=int, default=3)
    ap.add_argument("--material", default="steel")
    ap.add_argument("--half-thickness", type=float, default=0.0125)
    args = ap.parse_args()

    from calc_fem_kelvin import solve_fem
    from calc_common import EMMaterial
    mat = EMMaterial.from_name(args.material)
    f = args.frequency

    runs = {}

    # --- Run 1: A-V (full FEM with coil_sigma + SIBC) ---
    print("=" * 60)
    print(f"Run 1: A-V (coil volume meshed + coil_sigma=5.8e7, SIBC wp)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["av"] = solve_fem(vol_file=args.vol, frequency=f, mat=mat,
                           impedance_model="sibc",
                           half_thickness=args.half_thickness,
                           coil_sigma=5.8e7, solver="pardiso")
    runs["t_av"] = time.perf_counter() - t0

    # --- Run 2: PEEC + back-reaction (scattered A_r via SIBC Robin) ---
    print("=" * 60)
    print(f"Run 2: PEEC nwinc={args.peec_nwinc} + back-reaction (scattered A_r)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["peec"] = solve_peec_back_reaction(
        vol_file=args.vol, step=args.step, frequency=f, mat=mat,
        half_thickness=args.half_thickness,
        nwinc=args.peec_nwinc, nhinc=args.peec_nhinc)
    runs["t_peec"] = time.perf_counter() - t0

    # --- Results ---
    av = runs["av"]
    pe = runs["peec"]

    L_av = av.get("L", 0) or 0
    L_pe = pe["L_total"]
    P_av = av.get("P_total", 0) or 0
    P_pe = pe["P_wp"]

    print()
    print("=" * 60)
    print("COMPARISON")
    print("=" * 60)
    print(f"{'':18} {'A-V':>14} {'PEEC+BR':>14} {'PE/AV':>9}")
    print("-" * 60)
    print(f"{'L [nH]':18} {L_av*1e9:14.3f} {L_pe*1e9:14.3f} "
          f"{fmt_pct(L_pe, L_av):>9}")
    print(f"{'P [W]':18} {P_av:14.4e} {P_pe:14.4e} "
          f"{fmt_pct(P_pe, P_av):>9}")
    print()
    print(f"  PEEC decomposition:")
    print(f"    L_peec (circuit): {pe['L_peec']*1e9:.3f} nH")
    print(f"    Delta_L (FEM BR): {pe['Delta_L']*1e9:+.3f} nH  "
          f"{'(Lenz, PASS)' if pe['Delta_L'] < 0 else '(positive - flux conc.?)'}")
    print(f"    L_total         : {pe['L_total']*1e9:.3f} nH")
    print()
    print(f"  Times: A-V {runs['t_av']:.1f}s, PEEC+BR {runs['t_peec']:.1f}s")

    if "error" in av:
        print(f"  A-V ERROR: {av['error']}")


if __name__ == "__main__":
    main()
