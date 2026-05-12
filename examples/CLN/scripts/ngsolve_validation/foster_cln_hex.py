"""Foster-CLN library for arbitrary 8-node hexahedron (Phase 12).

Strategy:
  1. Hex defined by 8 vertices x_i, y_i, z_i (i=0..7) in standard FE ordering
  2. Trilinear shape functions N_i(xi, eta, zeta) on reference cube [-1, 1]^3
  3. Jacobian J(xi, eta, zeta) = det(d(x,y,z)/d(xi,eta,zeta))
  4. K element integrated over reference cube via mp.quad

Vertex ordering (standard FE 8-node hex):
  0: (-1,-1,-1)  4: (-1,-1,+1)
  1: (+1,-1,-1)  5: (+1,-1,+1)
  2: (+1,+1,-1)  6: (+1,+1,+1)
  3: (-1,+1,-1)  7: (-1,+1,+1)

Trilinear shape functions:
  N_i(xi, eta, zeta) = (1/8)(1 + xi_i*xi)(1 + eta_i*eta)(1 + zeta_i*zeta)

Test case 1: cuboid (degenerate hex) - should match foster_cln_python_pure result
Test case 2: deformed hex (parallelepiped) - new ground

Date: 2026-05-02
"""

import time
import json

import numpy as np
import mpmath as mp

mp.mp.dps = 50

# === Reference vertices (8-node hex in canonical ordering) ===
ref_vertices = np.array([
    [-1, -1, -1],  # 0
    [+1, -1, -1],  # 1
    [+1, +1, -1],  # 2
    [-1, +1, -1],  # 3
    [-1, -1, +1],  # 4
    [+1, -1, +1],  # 5
    [+1, +1, +1],  # 6
    [-1, +1, +1],  # 7
])


def shape_funcs(xi, eta, zeta):
    """Trilinear shape functions N_i(xi, eta, zeta) for i=0..7."""
    N = []
    for v in ref_vertices:
        xi_i, eta_i, zeta_i = v
        N.append((1 + xi_i*xi) * (1 + eta_i*eta) * (1 + zeta_i*zeta) / 8)
    return N


def shape_funcs_derivatives(xi, eta, zeta):
    """dN_i/dxi, dN_i/deta, dN_i/dzeta for i=0..7."""
    dN_dxi = []
    dN_deta = []
    dN_dzeta = []
    for v in ref_vertices:
        xi_i, eta_i, zeta_i = v
        dN_dxi.append(xi_i * (1 + eta_i*eta) * (1 + zeta_i*zeta) / 8)
        dN_deta.append((1 + xi_i*xi) * eta_i * (1 + zeta_i*zeta) / 8)
        dN_dzeta.append((1 + xi_i*xi) * (1 + eta_i*eta) * zeta_i / 8)
    return dN_dxi, dN_deta, dN_dzeta


def map_point(verts_phys, xi, eta, zeta):
    """Map reference point (xi, eta, zeta) to physical (x, y, z) using vertex coords."""
    N = shape_funcs(xi, eta, zeta)
    x = mp.fsum(N[i] * verts_phys[i, 0] for i in range(8))
    y = mp.fsum(N[i] * verts_phys[i, 1] for i in range(8))
    z = mp.fsum(N[i] * verts_phys[i, 2] for i in range(8))
    return x, y, z


def jacobian_det(verts_phys, xi, eta, zeta):
    """Determinant of Jacobian matrix |d(x,y,z)/d(xi,eta,zeta)|."""
    dN_dxi, dN_deta, dN_dzeta = shape_funcs_derivatives(xi, eta, zeta)
    J = mp.matrix(3, 3)
    for k in range(3):  # k: xi, eta, zeta
        dN_dk = [dN_dxi, dN_deta, dN_dzeta][k]
        for axis in range(3):  # axis: x, y, z
            J[k, axis] = mp.fsum(dN_dk[i] * verts_phys[i, axis] for i in range(8))
    # det(J)
    detJ = (J[0,0]*(J[1,1]*J[2,2] - J[1,2]*J[2,1])
            - J[0,1]*(J[1,0]*J[2,2] - J[1,2]*J[2,0])
            + J[0,2]*(J[1,0]*J[2,1] - J[1,1]*J[2,0]))
    return detJ


def hex_self_newton(verts_phys, p_idx=(0,0,0), p2_idx=(0,0,0), max_recurse=10):
    """Compute K[(p,q,r), (p',q',r')] = ∫_hex ∫_hex x^p y^q z^r x'^p' y'^q' z'^r' / |r-r'| dV dV'.

    For now, use 6D mp.quad. For 8-node hex with possibly 17-digit precision.
    """
    p, q, r = p_idx
    p2, q2, r2 = p2_idx
    verts = mp.matrix(verts_phys.tolist())

    def integrand(xi1, eta1, zeta1, xi2, eta2, zeta2):
        x1, y1, z1 = map_point(verts_phys, xi1, eta1, zeta1)
        x2, y2, z2 = map_point(verts_phys, xi2, eta2, zeta2)
        d = mp.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        if d == 0:
            return mp.mpf(0)  # singular point, 0 measure
        J1 = jacobian_det(verts_phys, xi1, eta1, zeta1)
        J2 = jacobian_det(verts_phys, xi2, eta2, zeta2)
        return (x1**p * y1**q * z1**r) * (x2**p2 * y2**q2 * z2**r2) / d * J1 * J2

    # 6D integration (slow!)
    val = mp.quad(integrand,
                  [-1, 1], [-1, 1], [-1, 1],
                  [-1, 1], [-1, 1], [-1, 1])
    return val


# === Test 1: cuboid (degenerate hex) ===
print("=" * 64)
print("Test 1: cuboid as degenerate hex - verify against foster_cln_python_pure")
print("=" * 64)

a, b, c = mp.mpf("5e-3"), mp.mpf("2e-3"), mp.mpf("1e-3")

# Cuboid vertices: [-a/2, a/2] x [-b/2, b/2] x [-c/2, c/2]
cuboid_verts = np.array([
    [-a/2, -b/2, -c/2],
    [+a/2, -b/2, -c/2],
    [+a/2, +b/2, -c/2],
    [-a/2, +b/2, -c/2],
    [-a/2, -b/2, +c/2],
    [+a/2, -b/2, +c/2],
    [+a/2, +b/2, +c/2],
    [-a/2, +b/2, +c/2],
], dtype=object)

print(f"Computing K[(0,1,0),(0,1,0)] = D_yy for cuboid via hex 6D integration...")
print(f"(Warning: 6D mp.quad is SLOW. May take 30-60 min for first pass.)")
t0 = time.time()
# Just sanity check with reduced integration depth
mp.mp.dps = 30  # lower for speed during initial test
K_test = hex_self_newton(cuboid_verts, (0, 1, 0), (0, 1, 0))
mp.mp.dps = 50
print(f"  K_hex[(0,1,0),(0,1,0)] = {mp.nstr(K_test, 20)}")
print(f"  Time: {time.time()-t0:.1f} s")

ref = mp.mpf("5.24918270035279297821091928151251136119220629893627828114e-21")
print(f"  Reference (cuboid)     = {mp.nstr(ref, 20)}")
print(f"  Relative diff = {mp.nstr(abs((K_test - ref)/ref), 6)}")
