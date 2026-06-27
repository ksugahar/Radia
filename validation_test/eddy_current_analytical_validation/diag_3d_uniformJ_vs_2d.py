"""Diagnose 3D SIBC vs 2D axisym SIBC discrepancy.

Bypasses the PEEC filament Biot-Savart A_s projection by imposing a
uniform volume J in the 3D "coil" material (= same as the 2D axisym
reference). Historical reference: the 3D "scattered" FEM path
(``solve_fem_biot_savart``) that motivated this diagnostic was retired
2026-04-24 due to an unresolved 3.4x P_wp under-prediction.  Production
now uses ``calc_fem_kelvin --formulation total`` or ``calc_fem_coilmesh``.

2D axisym SIBC (Cu, 7 kHz, closed torus):
  L = 57.65 nH, P = 6.62e-5 W   (validated against full-eddy 0.9% / 14%)

3D scattered-BS (Cu, 7 kHz, closed torus):
  L = 52.18 nH, P = 1.65e-5 W   (-9.5% / -75% vs 2D SIBC)

If uniform-J 3D reproduces ~58 / 6.6e-5, the A_s HCurl projection is
the culprit. Otherwise the Kelvin/mesh/formulation is off.
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src', 'radia'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'src', 'radia', 'panels'))

from em_material import EMMaterial, MU_0
NU_0 = 1.0 / MU_0


def solve_3d_uniformJ(vol_file, mat, frequency=7000, I_total=1.0,
                       a_coil=0.003, half_thickness=0.0125,
                       fes_order=1, solver='pardiso'):
    from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                         GridFunction, Integrate, Conj, curl, dx, ds, CF,
                         BND, VOL, TaskManager, sqrt, IfPos, x, y, z)
    from calc_fem_kelvin import detect_kelvin_offset

    mesh = Mesh(vol_file)
    materials = mesh.GetMaterials()
    boundaries = set(mesh.GetBoundaries())

    has_kelvin = 'kelvin' in materials
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

    # nu_cf: Kelvin weight, everything else nu_0
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

    # Uniform azimuthal J in the coil: J = I/(pi a^2) * e_phi
    # e_phi = (-y/r, x/r, 0) in 3D, with r = sqrt(x^2+y^2)
    J0 = I_total / (math.pi * a_coil ** 2)
    r_xy = sqrt(x * x + y * y)
    r_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_cf = J0 * CF((-y / r_safe, x / r_safe, 0))

    # SIBC Robin
    Z_s = mat.dowell_Zs(frequency, half_thickness)
    omega = 2 * math.pi * frequency
    robin = 1j * omega / Z_s
    print(f'  Z_s = {Z_s:.4e}, Robin = {robin:.4e}')

    dirichlet_bnd = 'GND' if 'GND' in boundaries else ''
    base_fes = HCurl(mesh, order=fes_order, nograds=True, complex=True,
                     dirichlet=dirichlet_bnd)
    fes = Periodic(base_fes) if has_kelvin_periodic else base_fes
    u, v = fes.TnT()
    print(f'  ndof = {fes.ndof}')

    a_bf = BilinearForm(fes)
    a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
    # Tiny mass term for gauge (same as scattered path)
    non_kelvin_mats = [m for m in materials if 'kelvin' not in m.lower()]
    if non_kelvin_mats:
        a_bf += 1e-6 * NU_0 * u * v * dx('|'.join(non_kelvin_mats))
    a_bf += robin * u.Trace() * v.Trace() * ds('sibc')

    f_lf = LinearForm(fes)
    f_lf += J_cf * v * dx('coil')

    t0 = time.perf_counter()
    a_bf.Assemble()
    f_lf.Assemble()
    t_asm = time.perf_counter() - t0

    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    with TaskManager():
        gfu.vec.data = a_bf.mat.Inverse(fes.FreeDofs(), inverse=solver) * f_lf.vec
    t_solve = time.perf_counter() - t0
    print(f'  t_asm = {t_asm:.1f}s, t_solve = {t_solve:.1f}s')

    # L from energy (volume) + skin-layer surface term
    # L_vol = 2 * W_vol / I^2  where W_vol = 0.5 * int nu |curl A|^2 dV
    W_vol = Integrate(0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)),
                      mesh, order=10).real
    L_vol = 2 * W_vol / I_total ** 2

    # Surface |A_t|^2 integral for H_t, P, L_skin
    from ngsolve import specialcf
    n_bnd = specialcf.normal(3)
    A_sq = sum(gfu[i].real * gfu[i].real + gfu[i].imag * gfu[i].imag
               for i in range(3))
    Adn_re = sum(gfu[i].real * n_bnd[i] for i in range(3))
    Adn_im = sum(gfu[i].imag * n_bnd[i] for i in range(3))
    An_sq = Adn_re * Adn_re + Adn_im * Adn_im
    At_sq = A_sq - An_sq

    wp_region = mesh.Boundaries('sibc')
    A_wp = Integrate(CF(1), mesh, BND, definedon=wp_region).real
    At_int = Integrate(At_sq, mesh, BND, definedon=wp_region).real
    H_t_rms = abs(1j * omega / Z_s) * math.sqrt(max(At_int, 0) / max(A_wp, 1e-30))
    P_surf = 0.5 * Z_s.real * H_t_rms ** 2 * A_wp
    L_skin = omega * Z_s.imag / (abs(Z_s) ** 2) * At_int / I_total ** 2
    L_total = L_vol + L_skin

    return {
        'L': L_total,
        'L_vol': L_vol,
        'L_skin': L_skin,
        'P_surf': P_surf,
        'H_t_rms': H_t_rms,
        'At_int': At_int,
        'A_wp': A_wp,
        'Z_s': Z_s,
    }


def main():
    print('=' * 70)
    print('3D uniform-J SIBC  vs  2D axisym SIBC (Cu, 7 kHz, closed torus)')
    print('=' * 70)
    vol = os.path.join(HERE, '..', '..', 'src', 'radia', 'panels', 'samples',
                        'ih_closed_torus.vol')
    mat = EMMaterial.from_name('copper')
    r = solve_3d_uniformJ(vol, mat, frequency=7000, I_total=1.0,
                           a_coil=0.003, half_thickness=0.0125,
                           fes_order=1)
    print()
    print(f'3D uniform-J SIBC (Cu, 7kHz):')
    print(f'  L        = {r["L"]*1e9:.3f} nH  (L_vol {r["L_vol"]*1e9:.3f} + L_skin {r["L_skin"]*1e9:+.3f})')
    print(f'  P_surf   = {r["P_surf"]:.4e} W')
    print(f'  H_t_rms  = {r["H_t_rms"]:.2f} A/m')
    print(f'  At_int   = {r["At_int"]:.4e}')
    print()
    print('2D axisym SIBC reference (Cu, 7kHz):')
    print(f'  L        = 57.647 nH')
    print(f'  P        = 6.6211e-05 W')


if __name__ == '__main__':
    main()
