"""Validation helpers for stream-function REGCOIL/cohomology checks.

This module is the validation-test copy of the reusable helper layer from the
public REGCOIL fusion demo.  Executable checks should import this regression
surface, while public demos stay in the docs/examples tier.
"""
from __future__ import annotations

import math

import numpy as np


R_MAJOR = 0.30
A_PLASMA = 0.06
A_WIND = 0.12
_MU0_4PI = 1.0e-7


def _torus_surface_vol(a, maxh, path):
    """Save a torus surface mesh for winding/plasma-boundary tests."""
    from netgen.occ import WorkPlane, Axes, Pnt, Y, Z, Axis, OCCGeometry

    wire = WorkPlane(Axes(Pnt(R_MAJOR, 0, 0), n=Y)).Circle(a).Wire()
    surf = wire.Revolve(Axis(Pnt(0, 0, 0), Z), 360)
    OCCGeometry(surf).GenerateMesh(maxh=maxh).Save(path)
    return path


def _plasma_points_normals(plasma, eval_max):
    """Circular-plasma sample points plus analytic outward torus normals."""
    pts = np.array([list(v.point) for v in plasma.vertices])
    if len(pts) > eval_max:
        pts = pts[np.linspace(0, len(pts) - 1, eval_max).astype(int)]
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    nrm = np.column_stack([
        pts[:, 0] - R_MAJOR * np.cos(phi),
        pts[:, 1] - R_MAJOR * np.sin(phi),
        pts[:, 2],
    ])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1.0e-30
    rho = np.hypot(pts[:, 0], pts[:, 1]) - R_MAJOR
    theta = np.arctan2(pts[:, 2], rho)
    return pts, nrm, theta, phi


def _A_normal(C, fes, n_cf, pts, nrm):
    """B.n design matrix for winding-surface stream-function DOFs."""
    A3 = C._assemble_biot_savart(fes, n_cf, pts, [0, 1, 2]).reshape(
        len(pts), 3, fes.ndof
    )
    return np.einsum("mc,mcj->mj", nrm, A3)


def _design(C, fes, coil, A_n, Bn, regularize, alpha_rel):
    """Solve the small REGCOIL-style normal-field fit used by validation."""
    from ngsolve import GridFunction

    S = C._seminorm(fes, regularize).toarray()
    AtA = A_n.T @ A_n
    alpha = alpha_rel * np.trace(AtA) / fes.ndof
    psi = np.linalg.solve(AtA + alpha * S, A_n.T @ Bn)
    res = float(np.linalg.norm(A_n @ psi - Bn) / (np.linalg.norm(Bn) + 1.0e-30))
    gfu = GridFunction(fes)
    gfu.vec.FV().NumPy()[:] = psi
    return psi, res, float(C._peak_current_density(fes, coil, gfu))


def _secular_currents(n_cf):
    """Analytic torus cohomology generators as surface currents."""
    from ngsolve import Cross, CoefficientFunction as CF, x, y, z, sqrt

    rho = sqrt(x * x + y * y)
    s = rho - R_MAJOR
    den = s * s + z * z
    grad_zeta = CF((-y, x, 0)) / (x * x + y * y)
    grad_theta = CF((-z * x / (rho * den), -z * y / (rho * den), s / den))
    return [
        ("net_poloidal_TF", Cross(n_cf, grad_zeta)),
        ("net_toroidal", Cross(n_cf, grad_theta)),
    ]


def _assemble_secular_normal(coil, pts, nrm, K_list):
    """Normal field at plasma points from each explicit secular current."""
    from ngsolve import x, y, z, sqrt, Integrate, ds

    out = np.zeros((len(pts), len(K_list)))
    for k, (_name, K) in enumerate(K_list):
        for m, p in enumerate(pts):
            dxt, dyt, dzt = p[0] - x, p[1] - y, p[2] - z
            r2 = dxt * dxt + dyt * dyt + dzt * dzt
            r3 = r2 * sqrt(r2)
            cr = (
                K[1] * dzt - K[2] * dyt,
                K[2] * dxt - K[0] * dzt,
                K[0] * dyt - K[1] * dxt,
            )
            out[m, k] = sum(
                nrm[m, c] * _MU0_4PI * Integrate(cr[c] / r3 * ds, coil)
                for c in range(3)
            )
    return out


def _tf_field_1overR(coil, n_cf, radii):
    """Toroidal field of the unit net-poloidal-current secular term."""
    from ngsolve import x, y, z, sqrt, Integrate, ds

    _name, K = _secular_currents(n_cf)[0]
    by = np.zeros(len(radii))
    for i, r in enumerate(radii):
        dxt, dyt, dzt = r - x, -y, -z
        r3 = (dxt * dxt + dyt * dyt + dzt * dzt) * sqrt(
            dxt * dxt + dyt * dyt + dzt * dzt
        )
        by[i] = _MU0_4PI * Integrate((K[2] * dxt - K[0] * dzt) / r3 * ds, coil)
    return by


def _betti1_winding_surface(a, maxh):
    """First Betti number b1 of the torus winding surface, gmsh-free."""
    nu = max(6, int(round(2.0 * math.pi * R_MAJOR / maxh)))
    nv = max(6, int(round(2.0 * math.pi * a / maxh)))

    def vid(i, j):
        return (i % nu) * nv + (j % nv)

    edges = set()
    nfaces = 0
    for i in range(nu):
        for j in range(nv):
            v00, v10 = vid(i, j), vid(i + 1, j)
            v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
            for tri in ((v00, v10, v11), (v00, v11, v01)):
                nfaces += 1
                for e in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                    edges.add(frozenset(e))
    chi = nu * nv - len(edges) + nfaces
    return 2 - chi
