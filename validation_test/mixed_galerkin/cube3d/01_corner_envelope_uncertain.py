"""
3D cube mixed Galerkin: rank-N Foster bulk + corner-aware envelope.

Phase 6 of the 2026-05-28 -> 2026-06-12 research sprint.

## STATUS: PROVISIONAL — task #183 (3D cube Foster reference audit) is open.

This script's result depends on `_references/cube3d_foster.Y_exact_cube3d`,
which is itself provisional pending the audit.  The Phase 6 conclusion
reported a 20% wall-band peak; given the 2D-square experience (Phase
8c, see square2d/02_foster_reference_audit.py), most of that 20% is
probably Foster truncation artifact and the true mixed-Galerkin error
is in the 0.1% - 1% range.  The audit will confirm or refute this.

## Setup

The natural 3-D extension of the 2D square corner envelope.  Volume
element  dV = dx dy dz  (no curvilinear factor since cube is
rectangular).  Bulk Foster modes are orthogonal sine modes
phi_{npq}(x, y, z) = sin(n pi x/L) sin(p pi y/L) sin(q pi z/L) for
n, p, q all odd.

Surface envelope (1-DOF, corner-aware):

    psi(x, y, z; s) = f(x) f(y) f(z)
    f(xi) = cosh((xi - L/2) t) / cosh(L t / 2) - 1

psi vanishes on all 6 faces (any face has one f = 0).  Near each face
(away from edges and corners), psi ~ exp(-d_face t) - 1, the standard
slab SIBC envelope.  Near each 90deg edge (where 2 faces meet), psi
~ r^2 sin(2 theta), the 2D Wiener-Hopf wedge asymptote.  Near each
octant corner (where 3 faces meet), psi ~ x y z ~ r^3 Y_{3,0}, the
spherical-harmonic eigenfunction of the octant Dirichlet problem.

These qualitative asymptotes match the true solution, but a single
tensor-product DOF cannot simultaneously give the correct
quantitative weight to all 3 boundary-element classes (face, edge,
corner).  Whether this insufficiency or the Foster truncation
dominates the apparent 20% error is the subject of task #183.

## Provisional result (vs Foster N=99, which is under-converged)

    wall band max ~ 20% (Phase 6)

Pending Phase 8c-style audit, treat as ROUGH UPPER BOUND only.
"""

from __future__ import annotations

import math
import cmath
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import (  # noqa: E402
    K_SIBC_cube3d,
    Y_DC_cube3d,
    Y_exact_cube3d,
    Y_foster_cube3d,
    Y_mellin_cube3d,
)

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
K_SIBC = K_SIBC_cube3d(L, SIGMA, MU)
V_CUBE = L**3

N_PSI = 49  # sine modes for psi expansion per axis (3D = 49^3 ~ 1.2e5 terms)
N_BULK = 1


def f_n(n: int, t):
    return -4 * L**2 * t**2 / (n * math.pi * (n**2 * math.pi**2 + L**2 * t**2))


def F_0(t):
    return (2.0 / t) * cmath.tanh(L * t / 2) - L


def Y_mixed_galerkin(s):
    """rank-1 Foster bulk + 1-DOF corner envelope psi = f(x) f(y) f(z)."""
    t = cmath.sqrt(s * MS)

    odd_n = np.arange(1, 2 * N_PSI, 2)
    f_arr = np.array([f_n(int(n), t) for n in odd_n], dtype=complex)

    # Surface block via 3D outer product (n, p, q) all odd.
    nn, pp, qq = np.meshgrid(odd_n, odd_n, odd_n, indexing="ij")
    lam_3d = (nn**2 + pp**2 + qq**2) * math.pi**2 / L**2
    f_3d = np.einsum("i,j,k->ijk", f_arr, f_arr, f_arr)
    K_ss = (L / 2)**3 * np.sum((lam_3d + s * MS) * f_3d**2)
    b_psi = F_0(t)**3

    # Bulk: phi_{1,1,1} = sin(pi x/L) sin(pi y/L) sin(pi z/L)
    lambda_111 = 3 * math.pi**2 / L**2
    K_b = (lambda_111 + s * MS) * (L / 2)**3
    K_bs = (lambda_111 + s * MS) * (L / 2)**3 * f_arr[0]**3
    b_b = 8 * L**3 / math.pi**3

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


def main():
    print("=== 3D cube mixed Galerkin: rank-1 bulk + corner envelope (PROVISIONAL) ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e} S/m")
    print(f"Y_DC   = {Y_DC:.4e} S*m^2")
    print(f"K_SIBC = {K_SIBC:.4e}")
    print()
    print("*** REMINDER: this result is PROVISIONAL pending task #183.")
    print("*** Foster ref is N=99/199/399 Aitken at low-mid f, Mellin at high f.")
    print()

    print(f"{'f (Hz)':>10}  {'|Y_ref|':>12}  {'|Y_mixed|':>12}  {'rel err vs ref':>14}")
    for f in [1.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_ref = Y_exact_cube3d(s, L, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        err = abs(Y_ref - Y_m) / abs(Y_ref)
        print(f"{f:10.2e}  {abs(Y_ref):12.4e}  {abs(Y_m):12.4e}  {err*100:13.3f}%")

    print()
    print("Also compare to Mellin asymptote directly at high f:")
    print(f"{'f (Hz)':>10}  {'|Y_Mellin|':>12}  {'|Y_mixed|':>12}  {'rel err vs Mellin':>16}")
    for f in [1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_M = Y_mellin_cube3d(s, L, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        err = abs(Y_M - Y_m) / abs(Y_M)
        print(f"{f:10.2e}  {abs(Y_M):12.4e}  {abs(Y_m):12.4e}  {err*100:15.3f}%")


if __name__ == "__main__":
    main()
