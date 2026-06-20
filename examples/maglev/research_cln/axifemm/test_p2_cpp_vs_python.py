"""test_p2_cpp_vs_python.py — Phase B2 build smoke test.

Verify that the C++ AxiHenrotteFE_P2_Triangle (Phase B2) matches the Python
AxifemmP2Triangle (Phase B1c) at machine precision:
  - shape function values at random interior points
  - shape function gradients at random interior points

Does NOT test the integrator matrices (P2TriangleStiffness/SigmaMass) since
those require an NGSolve mesh + GetFE dispatch; they will be tested next
session via a disk eigenvalue benchmark.
"""
from __future__ import annotations

import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

from radia.axifem import AxiHenrotteFE_P2_Triangle as P2_CPP

from axifemm_p2_triangle import AxifemmP2Triangle


def _make_p2_triangle_nodes(r0, z0, r1, z1, r2, z2):
    rn = [r0, r1, r2,
          0.5 * (r0 + r1), 0.5 * (r1 + r2), 0.5 * (r2 + r0)]
    zn = [z0, z1, z2,
          0.5 * (z0 + z1), 0.5 * (z1 + z2), 0.5 * (z2 + z0)]
    return rn, zn


def test_one(label, r0, z0, r1, z1, r2, z2):
    print(f"\n[{label}]  vertices: "
          f"({r0:.4g}, {z0:.4g}), ({r1:.4g}, {z1:.4g}), ({r2:.4g}, {z2:.4g})")
    rn, zn = _make_p2_triangle_nodes(r0, z0, r1, z1, r2, z2)
    py_tri = AxifemmP2Triangle(np.array(rn), np.array(zn))
    cpp_tri = P2_CPP(rn, zn)

    # Sample at 5 random interior points in reference triangle T_0
    rng = np.random.default_rng(seed=42)
    n_samples = 5
    max_shape_err = 0.0
    max_grad_err = 0.0

    for k in range(n_samples):
        # Random (xi, eta) inside reference triangle
        xi = rng.uniform(0, 1)
        eta = rng.uniform(0, 1 - xi)

        # Physical (r, z) via affine map
        rp = rn[0] + (rn[1] - rn[0]) * xi + (rn[2] - rn[0]) * eta
        zp = zn[0] + (zn[1] - zn[0]) * xi + (zn[2] - zn[0]) * eta

        # Python evaluation
        phi_py = py_tri.shape(rp, zp)
        grad_py = py_tri.shape_grad(rp, zp)  # (6, 2) d phi / d(r, z)
        # Convert to d/d(xi, eta) like C++ does:
        # d phi / dxi = d phi/dr * dr/dxi + d phi/dz * dz/dxi
        dr_dxi = rn[1] - rn[0]
        dr_deta = rn[2] - rn[0]
        dz_dxi = zn[1] - zn[0]
        dz_deta = zn[2] - zn[0]
        grad_py_xi = np.column_stack([
            grad_py[:, 0] * dr_dxi + grad_py[:, 1] * dz_dxi,
            grad_py[:, 0] * dr_deta + grad_py[:, 1] * dz_deta,
        ])

        # C++ evaluation (via the IntegrationPoint interface)
        ip = (xi, eta)
        phi_cpp = np.array(cpp_tri.CalcShape(ip))
        grad_cpp = np.array(cpp_tri.CalcDShape(ip))  # (6, 2) d phi / d(xi, eta)

        shape_err = np.linalg.norm(phi_py - phi_cpp) / max(
            np.linalg.norm(phi_py), 1e-30)
        grad_err = np.linalg.norm(grad_py_xi - grad_cpp) / max(
            np.linalg.norm(grad_py_xi), 1e-30)
        max_shape_err = max(max_shape_err, shape_err)
        max_grad_err = max(max_grad_err, grad_err)

    tag_shape = "OK" if max_shape_err < 1e-10 else "FAIL"
    tag_grad = "OK" if max_grad_err < 1e-8 else "FAIL"
    print(f"  shape max rel err: {max_shape_err:.2e} [{tag_shape}]")
    print(f"  grad max rel err: {max_grad_err:.2e} [{tag_grad}]")
    return max_shape_err < 1e-10 and max_grad_err < 1e-8


def main():
    print("=" * 64)
    print("Phase B2 build smoke test: C++ P2 vs Python P2 reference")
    print("=" * 64)
    triangles = [
        ("general", 1e-3, 0.0, 5e-3, 0.0, 3e-3, 2e-3),
        ("axis-touch-vertex0", 0.0, 0.0, 4e-3, 0.0, 2e-3, 3e-3),
        ("axis-both-vertex01", 0.0, 0.0, 0.0, 2e-3, 5e-3, 1e-3),
        ("skew-general", 2e-3, 1e-3, 6e-3, 0.5e-3, 4e-3, 4e-3),
    ]
    all_ok = True
    for tri_args in triangles:
        ok = test_one(*tri_args)
        if not ok:
            all_ok = False
    print()
    print("PASS — C++ P2 matches Python P2 at machine precision"
          if all_ok else "FAIL — see above")


if __name__ == "__main__":
    main()
