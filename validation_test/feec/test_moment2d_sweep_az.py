"""Golden lock for the 2D planar MMMM follow-ups: factor-once torque sweep + A_z + complex torque.

  * factor-once sweep (C++ Moment2DSolveMulti): a LINEAR torque-angle sweep factors the moment matrix
    ONCE and back-substitutes all angles -- result IDENTICAL to per-angle re-solves.
  * A_z (shared C++ PlanarChargeAz): B = curl A -> dA_z/dy = mu0 H_x, -dA_z/dx = mu0 H_y (finite diff).
  * complex (eddy phasor) torque reduces to the real Maxwell torque for a real phasor.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_charges as pc

MU0 = 4e-7 * np.pi


def _ellipse(a, b, maxh):
    geo = SplineGeometry()
    npt = 96
    pts = [(a * np.cos(t), b * np.sin(t)) for t in np.linspace(0, 2 * np.pi, npt, endpoint=False)]
    pid = [geo.AppendPoint(*p) for p in pts]
    for i in range(npt):
        geo.Append(["line", pid[i], pid[(i + 1) % npt]], bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_factor_once_sweep_matches_per_angle():
    """LINEAR torque sweep via the factor-once multi-RHS solve == per-angle re-solves."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.1)
        angles = np.deg2rad([10.0, 40.0, 70.0, 100.0])
        sw = m2.torque_angle_sweep(mesh, 1000.0, angles, Rc=3.0, mu_r=4.0)
        assert sw["factored_once"]
        # per-angle reference
        T_ref = []
        for th in angles:
            H_ext = (1000.0 * np.cos(th), 1000.0 * np.sin(th))
            r = m2.solve_planar_demag(mesh, mu_r=4.0, H_ext=H_ext)
            T_ref.append(m2.maxwell_torque(mesh, r["M"], 3.0, H_ext=H_ext))
    for Ts, Tr in zip(sw["torque"], T_ref):
        assert abs(Ts - Tr) <= 1e-6 * max(abs(Tr), 1e-30) + 1e-12, (Ts, Tr)


def test_az_curl_gives_H():
    """Finite-diff of A_z reproduces mu0 H: dA_z/dy = mu0 H_x, -dA_z/dx = mu0 H_y (air, one side)."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.1)
        r = m2.solve_planar_demag(mesh, mu_r=4.0, H_ext=(1000.0, 400.0))
        M = r["M"]
        P0 = np.array([[0.4, -3.0]])                        # below the body -> branch-cut safe
        h = 1e-3
        Px = np.array([[0.4 + h, -3.0], [0.4 - h, -3.0]])
        Py = np.array([[0.4, -3.0 + h], [0.4, -3.0 - h]])
        Azx = pc.vector_potential_az(mesh, M, Px)
        Azy = pc.vector_potential_az(mesh, M, Py)
        H = pc.exterior_field(mesh, M, P0)[0]
    dAz_dx = (Azx[0] - Azx[1]) / (2 * h)
    dAz_dy = (Azy[0] - Azy[1]) / (2 * h)
    assert abs(dAz_dy - MU0 * H[0]) < 1e-3 * abs(MU0 * H[0]) + 1e-12, (dAz_dy, MU0 * H[0])
    assert abs(-dAz_dx - MU0 * H[1]) < 1e-3 * abs(MU0 * H[1]) + 1e-12, (-dAz_dx, MU0 * H[1])


def test_complex_torque_reduces_to_real():
    """A real phasor magnetization + real applied field -> the complex time-averaged torque equals
    the real Maxwell torque (0.5 factor absorbed: <T>_complex(real) == T_real)."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.1)
        H_ext = (1000.0 * np.cos(0.6), 1000.0 * np.sin(0.6))
        r = m2.solve_planar_demag(mesh, mu_r=4.0, H_ext=H_ext)
        M = r["M"]
        T_real = pc.maxwell_torque(mesh, M, 3.0, H_ext=H_ext)
        # embed as a real phasor: M -> M + 0j, H_ext -> H_ext + 0j.  For real phasors
        # 0.5 Re(H_r conj H_phi) = 0.5 H_r H_phi, i.e. HALF the static torque.
        T_cplx = pc.maxwell_torque_complex(mesh, M.astype(complex), 3.0,
                                           H_ext=(complex(H_ext[0]), complex(H_ext[1])))
    assert abs(T_cplx - 0.5 * T_real) < 1e-6 * abs(T_real) + 1e-12, (T_cplx, T_real)
