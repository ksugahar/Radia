"""Debug script: re-implement spherical Duffy quadrature in pure Python and
compare to the C++ implementation. If both agree but disagree with Mathematica,
the math is wrong. If C++ disagrees with Python, the C++ has a bug.
"""
import numpy as np
import math
from radia_vim import HDivDivFreeHexBasis, compute_K_entry_spherical_duffy


def basis_block_A_p3_idx0_py(x, y, z):
    """Block A, order 3, (i,j,k)=(0,0,0).
    L_2(s) = (s²-1)/2,  L_2'(s) = s
    P_1(η) = η
    S.x =  8 L_2(ξ) * P_1(η) * L_2'(ζ)
    S.z = -8 L_2'(ξ) * P_1(η) * L_2(ζ)
    """
    xi = 2*x - 1; eta = 2*y - 1; zeta = 2*z - 1
    Sx =  8 * (xi*xi - 1)/2 * eta * zeta
    Sy = 0.0
    Sz = -8 * xi * eta * (zeta*zeta - 1)/2
    return Sx, Sy, Sz


def K00_spherical_duffy_py(a, b, c, n_t, n_p, n_r, n_c):
    pi = math.pi
    # GL nodes on [-1, 1]
    gl_t = np.polynomial.legendre.leggauss(n_t)
    gl_p = np.polynomial.legendre.leggauss(n_p)
    gl_r = np.polynomial.legendre.leggauss(n_r)
    gl_c = np.polynomial.legendre.leggauss(n_c)

    result = 0.0
    for it in range(n_t):
        theta = (gl_t[0][it] + 1) * pi / 2
        jac_t = pi / 2
        sth = math.sin(theta)
        cth = math.cos(theta)
        for ip in range(n_p):
            phi = (gl_p[0][ip] + 1) * pi
            jac_p = pi
            cph = math.cos(phi)
            sph = math.sin(phi)
            wx = sth * cph; wy = sth * sph; wz = cth
            wmax = max(abs(wx), abs(wy), abs(wz))
            if wmax < 1e-30:
                continue
            rho_max = 1.0 / wmax
            Lwn = math.sqrt(a*a*wx*wx + b*b*wy*wy + c*c*wz*wz)

            sphere_w = gl_t[1][it] * gl_p[1][ip] * jac_t * jac_p * sth / Lwn

            for ir in range(n_r):
                rho = (gl_r[0][ir] + 1) * rho_max / 2
                jac_r = rho_max / 2
                ux = rho * wx; uy = rho * wy; uz = rho * wz

                cx_lo = abs(ux)/2; cx_hi = 1 - abs(ux)/2; cx_jac = (cx_hi - cx_lo)/2
                cy_lo = abs(uy)/2; cy_hi = 1 - abs(uy)/2; cy_jac = (cy_hi - cy_lo)/2
                cz_lo = abs(uz)/2; cz_hi = 1 - abs(uz)/2; cz_jac = (cz_hi - cz_lo)/2
                if cx_jac <= 0 or cy_jac <= 0 or cz_jac <= 0:
                    continue

                c_inner = 0.0
                for icx in range(n_c):
                    cx = cx_lo + (gl_c[0][icx] + 1) * cx_jac
                    for icy in range(n_c):
                        cy = cy_lo + (gl_c[0][icy] + 1) * cy_jac
                        for icz in range(n_c):
                            cz = cz_lo + (gl_c[0][icz] + 1) * cz_jac

                            rx = cx + ux/2; ry = cy + uy/2; rz = cz + uz/2
                            rpx = cx - ux/2; rpy = cy - uy/2; rpz = cz - uz/2

                            Sx_a, Sy_a, Sz_a = basis_block_A_p3_idx0_py(rx, ry, rz)
                            Sx_b, Sy_b, Sz_b = basis_block_A_p3_idx0_py(rpx, rpy, rpz)

                            dot = Sx_a*Sx_b*a*a + Sy_a*Sy_b*b*b + Sz_a*Sz_b*c*c

                            w = gl_c[1][icx] * gl_c[1][icy] * gl_c[1][icz]
                            c_inner += w * dot
                c_inner *= cx_jac * cy_jac * cz_jac

                result += sphere_w * gl_r[1][ir] * jac_r * rho * c_inner
    return result


def main():
    basis = HDivDivFreeHexBasis(order=3)

    # Verify that Python basis and C++ basis agree at a few points
    for pt in [(0.3, 0.7, 0.2), (0.5, 0.5, 0.5), (0.1, 0.9, 0.8)]:
        cpp = basis.evaluate(0, *pt)
        py  = basis_block_A_p3_idx0_py(*pt)
        diff = max(abs(cpp[c] - py[c]) for c in range(3))
        print(f"basis check pt={pt} cpp={cpp} py={py} max|diff|={diff:.3e}")
    print()

    # Run Python and C++ Spherical Duffy at the same N
    for nt, np_, nr, nc in [(8, 16, 8, 4), (12, 24, 12, 6)]:
        py_val  = K00_spherical_duffy_py(5e-3, 2e-3, 1e-3, nt, np_, nr, nc)
        cpp_val = compute_K_entry_spherical_duffy(basis, 0, 0, 5e-3, 2e-3, 1e-3,
                                                   nt, np_, nr, nc)
        print(f"N=({nt},{np_},{nr},{nc}):  py={py_val:.10e}  cpp={cpp_val:.10e}  "
              f"py/cpp={py_val/cpp_val:.4f}")


if __name__ == "__main__":
    main()
