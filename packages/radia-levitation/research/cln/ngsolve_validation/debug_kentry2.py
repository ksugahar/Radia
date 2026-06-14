"""Debug 2: compare full quadrature point set Phi values + manually sum K[26,24]."""
import sys
from pathlib import Path
import numpy as np
import math

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))

from hex_vim_cupy import evaluate_basis
import radia_vim


def hand_compute_K_entry(basis, i_dof, j_dof, a, b, c,
                         n_th=8, n_ph=12, n_rh=8, n_c=4):
    """Manual (slow) Python translation of compute_K_entry_spherical_duffy.
    Should match C++ exactly within FP."""
    pi = math.pi
    th_n, th_w = np.polynomial.legendre.leggauss(n_th)
    ph_n, ph_w = np.polynomial.legendre.leggauss(n_ph)
    rh_n, rh_w = np.polynomial.legendre.leggauss(n_rh)
    c_n, c_w = np.polynomial.legendre.leggauss(n_c)

    result = 0.0
    for it in range(n_th):
        theta = (th_n[it] + 1) * pi / 2
        jac_th = pi / 2
        s_th = math.sin(theta)
        c_th = math.cos(theta)

        for ip in range(n_ph):
            phi_ang = (ph_n[ip] + 1) * pi
            jac_ph = pi
            c_ph = math.cos(phi_ang)
            s_ph = math.sin(phi_ang)

            omega_x = s_th * c_ph
            omega_y = s_th * s_ph
            omega_z = c_th
            omega_max = max(abs(omega_x), abs(omega_y), abs(omega_z))
            if omega_max < 1e-30:
                continue
            rho_max = 1.0 / omega_max
            Lomega = math.sqrt(a*a*omega_x**2 + b*b*omega_y**2 + c*c*omega_z**2)
            sphere_w = th_w[it] * ph_w[ip] * jac_th * jac_ph * s_th / Lomega

            for ir in range(n_rh):
                rho_val = (rh_n[ir] + 1) * rho_max / 2
                jac_rh = rho_max / 2
                u_x = rho_val * omega_x
                u_y = rho_val * omega_y
                u_z = rho_val * omega_z
                cx_lo = abs(u_x) / 2
                cy_lo = abs(u_y) / 2
                cz_lo = abs(u_z) / 2
                cx_jac = (1.0 - abs(u_x)) / 2
                cy_jac = (1.0 - abs(u_y)) / 2
                cz_jac = (1.0 - abs(u_z)) / 2
                if cx_jac <= 0 or cy_jac <= 0 or cz_jac <= 0:
                    continue

                c_inner = 0.0
                for icx in range(n_c):
                    cx = cx_lo + (c_n[icx] + 1) * cx_jac
                    wcx = c_w[icx]
                    for icy in range(n_c):
                        cy = cy_lo + (c_n[icy] + 1) * cy_jac
                        wcy = c_w[icy]
                        for icz in range(n_c):
                            cz = cz_lo + (c_n[icz] + 1) * cz_jac
                            wcz = c_w[icz]

                            rx = cx + u_x / 2
                            ry = cy + u_y / 2
                            rz = cz + u_z / 2
                            rpx = cx - u_x / 2
                            rpy = cy - u_y / 2
                            rpz = cz - u_z / 2

                            phi_i = basis.evaluate(i_dof, rx, ry, rz)
                            phi_j = basis.evaluate(j_dof, rpx, rpy, rpz)

                            dot = (phi_i[0] * phi_j[0] * a*a +
                                   phi_i[1] * phi_j[1] * b*b +
                                   phi_i[2] * phi_j[2] * c*c)
                            c_inner += wcx * wcy * wcz * dot
                c_inner *= cx_jac * cy_jac * cz_jac

                result += sphere_w * rh_w[ir] * jac_rh * rho_val * c_inner

    return result


def main():
    p = 2
    a, b, c = 5e-3, 2e-3, 1e-3
    basis = radia_vim.HDivDivFreeHexBasis(p)

    print("Hand-computing K[24, 26] via Python loop (slow)...")
    val_py = hand_compute_K_entry(basis, 24, 26, a, b, c)
    print(f"  Python hand K[24, 26] = {val_py:.6e}")

    val_cpp = radia_vim.compute_K_entry_spherical_duffy(
        basis, 24, 26, a, b, c, 8, 12, 8, 4)
    print(f"  C++ K[24, 26]         = {val_cpp:.6e}")

    print()
    print("Hand-computing K[26, 24] via Python loop (different basis args)...")
    val_py_swap = hand_compute_K_entry(basis, 26, 24, a, b, c)
    print(f"  Python hand K[26, 24] = {val_py_swap:.6e}")

    val_cpp_swap = radia_vim.compute_K_entry_spherical_duffy(
        basis, 26, 24, a, b, c, 8, 12, 8, 4)
    print(f"  C++ K[26, 24]         = {val_cpp_swap:.6e}")


if __name__ == "__main__":
    main()
