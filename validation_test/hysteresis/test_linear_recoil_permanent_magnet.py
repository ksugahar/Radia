"""Analytic sphere gate for the level-2 linear-recoil permanent magnet."""

import numpy as np
import ngsolve as ng
from netgen.occ import OCCGeometry, Pnt, Sphere

from radia import vim


MU0 = 4.0e-7 * np.pi


def test_curved_sphere_matches_the_analytic_linear_recoil_load_line():
    mu_rec = 1.05
    chi = mu_rec - 1.0
    B_r = np.array([0.0, 0.0, 1.2])
    H_ext = np.array([2.0e4, -1.0e4, 3.0e4])
    mesh = ng.Mesh(
        OCCGeometry(Sphere(Pnt(0.0, 0.0, 0.0), 1.0)).GenerateMesh(maxh=0.7)
    )
    mesh.Curve(2)

    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            mu_r=mu_rec,
            B_r=B_r,
            H_ext=ng.CoefficientFunction(tuple(H_ext)),
            curve_order=2,
            gram_eps=1.0e-8,
            tol=1.0e-9,
        )

    # A sphere has N=I/3.  With M=B_r/mu0+chi*H and H=H_ext-N*M,
    # every Cartesian component has this exact closed-form load line.
    expected = (B_r / MU0 + chi * H_ext) / (1.0 + chi / 3.0)
    relative_error = np.linalg.norm(result["M_avg"] - expected) / np.linalg.norm(expected)

    assert result["permanent_magnet_level"] == 2
    assert abs(result["demag"] - 1.0 / 3.0) < 5.0e-4
    assert relative_error < 2.0e-5
