"""Golden lock for the SHARED 2D planar exterior field + Maxwell torque (radia.planar_charges,
C++ rad_planar_charges), driven by the MMMM demag solve (radia.mmmm2d).

The exterior field of a solved planar body is the field of its magnetization (M.n equivalent bound
charge, log kernel); the Maxwell-stress torque on a circle in air reproduces the reluctance torque
T = mu0 A (M_avg x H0).  This is the SAME shared routine the HDiv-VIM uses (both feed a per-element
M), so it locks the common motor-postprocessing layer.

Checks (linear elliptic cylinder a=2 (x), b=1 (y), mu_r):
  * Maxwell-circle torque == mu0 A (M_avg x H0) to a few percent, at a generic field angle
  * torque ~ 0 when H0 is along a principal axis (M parallel H0)
  * torque sign flips between +45 deg and -45 deg
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_charges as pc

MU0 = 4e-7 * np.pi


def _ellipse(a, b, maxh, quad=False):
    geo = SplineGeometry()
    npt = 128
    pts = [(a * np.cos(t), b * np.sin(t)) for t in np.linspace(0, 2 * np.pi, npt, endpoint=False)]
    pid = [geo.AppendPoint(*p) for p in pts]
    for i in range(npt):
        geo.Append(["line", pid[i], pid[(i + 1) % npt]], bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh, quad_dominated=quad))


def _torque_and_analytic(mesh, mu_r, theta, Rc=3.0, H0=1000.0):
    H_ext = (H0 * np.cos(theta), H0 * np.sin(theta))
    r = m2.solve_planar_demag(mesh, mu_r=mu_r, H_ext=H_ext)
    A = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))
    Mx, My = r["M_avg"]
    T_analytic = MU0 * A * (Mx * H_ext[1] - My * H_ext[0])       # mu0 A (M x H0)
    T_maxwell = m2.maxwell_torque(mesh, r["M"], Rc, H_ext=H_ext)
    return T_maxwell, T_analytic


def test_maxwell_torque_matches_closed_form():
    """Maxwell-circle torque == mu0 A (M_avg x H0) at a generic angle (3-way-style check)."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.08)
        Tm, Ta = _torque_and_analytic(mesh, mu_r=4.0, theta=np.deg2rad(45.0))
    assert abs(Ta) > 0, "analytic torque must be nonzero at 45 deg"
    assert abs(Tm - Ta) / abs(Ta) < 3e-2, (Tm, Ta)


def test_torque_zero_on_principal_axis():
    """H0 along x (a principal axis) -> M parallel H0 -> torque ~ 0."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.08)
        Tm, Ta = _torque_and_analytic(mesh, mu_r=4.0, theta=0.0)
        # normalize by a 45-deg torque scale
        Tref, _ = _torque_and_analytic(mesh, mu_r=4.0, theta=np.deg2rad(45.0))
    assert abs(Tm) < 1e-2 * abs(Tref), (Tm, Tref)


def test_torque_sign_flips():
    """Torque at +45 deg and -45 deg have opposite sign (restoring toward the easy axis)."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.08)
        Tp, _ = _torque_and_analytic(mesh, mu_r=4.0, theta=np.deg2rad(45.0))
        Tn, _ = _torque_and_analytic(mesh, mu_r=4.0, theta=np.deg2rad(-45.0))
    assert Tp * Tn < 0, (Tp, Tn)


def test_shared_layer_is_method_agnostic():
    """planar_charges.maxwell_torque takes (mesh, M_elem) -- it does not know which solver produced
    M.  Feeding a hand-set uniform M reproduces mu0 A (M x H0) (the HDiv-VIM would feed its own M)."""
    with ng.TaskManager():
        mesh = _ellipse(2.0, 1.0, 0.08)
        nEl = mesh.ne
        M = np.tile(np.array([3.0e5, 1.0e5]), (nEl, 1))          # arbitrary uniform magnetization
        A = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))
        H0 = (0.0, 0.0)                                          # self-torque of a uniform M is ~0
        T = pc.maxwell_torque(mesh, M, Rc=3.0, H_ext=H0)
    # a uniform magnetization in NO applied field exerts no net self-torque about its centroid
    assert abs(T) < 1e-6 * MU0 * A * 3.0e5, T
