"""Golden lock for the SHARED 2D planar Maxwell FORCE (radia.planar_charges / rad_planar_charges).

Force per unit length on a circle in air: F_i = mu0 Rc oint [H_r H_i - 1/2 |H|^2 n_i] dphi.  The
maglev / actuator inter-body force -- shared by MMMM (radia.mmmm2d) and HDiv-VIM (both feed M.n clouds).

  * a UNIFORM applied field exerts ZERO net force on a body (only torque)
  * Newton's 3rd law: F(on A) == -F(on B) for two bodies
  * two x-magnetized bodies separated along x ATTRACT (moments along the line)
  * complex (eddy) force reduces to the real force for a real phasor
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_charges as pc

MU0 = 4e-7 * np.pi


def _disk(cx, a=1.0, maxh=0.12):
    geo = SplineGeometry(); geo.AddCircle((cx, 0.0), r=a, bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_uniform_field_zero_force():
    """A magnetized body in a UNIFORM applied field feels no net force (only torque)."""
    with ng.TaskManager():
        mesh = _disk(0.0)
        r = m2.solve_planar_demag(mesh, mu_r=4.0, H_ext=(1000.0, 500.0))
        F = m2.maxwell_force(mesh, r["M"], Rc=3.0, H_ext=(1000.0, 500.0))
        scale = MU0 * 3.0 * np.linalg.norm(r["M_avg"]) ** 2
    assert np.linalg.norm(F) < 1e-3 * scale, (F, scale)


def test_two_body_newton_and_attraction():
    """Two x-magnetized disks separated along x: Newton's 3rd law + attraction (moments on the line)."""
    with ng.TaskManager():
        mA = _disk(-3.0); mB = _disk(+3.0)
        MA = np.tile([5.0e5, 0.0], (mA.ne, 1))
        MB = np.tile([5.0e5, 0.0], (mB.ne, 1))
        F_A = m2.force_between([(mA, MA), (mB, MB)], Rc=1.6, center=(-3.0, 0.0))
        F_B = m2.force_between([(mA, MA), (mB, MB)], Rc=1.6, center=(+3.0, 0.0))
    # Newton's 3rd law
    assert np.allclose(F_A, -F_B, rtol=2e-2, atol=1e-9), (F_A, F_B)
    # A is at x=-3, B at x=+3; aligned moments along the separation -> ATTRACT -> F_A points +x
    assert F_A[0] > 0 and abs(F_A[1]) < 0.05 * abs(F_A[0]), F_A


def test_complex_force_reduces_to_real():
    """A real phasor magnetization + real applied field -> the complex time-averaged force is HALF
    the real Maxwell force (the 0.5 Re factor), component-wise."""
    with ng.TaskManager():
        mesh = _disk(0.0)
        # non-uniform 'applied' field is needed for a nonzero force; use a second body as the source
        mB = _disk(4.0)
        MB = np.tile([6.0e5, 0.0], (mB.ne, 1))
        r = m2.solve_planar_demag(mesh, mu_r=4.0, H_ext=(1000.0, 0.0))
        F_real = m2.force_between([(mesh, r["M"]), (mB, MB)], Rc=1.6, center=(0.0, 0.0))
        # embed body A as a real phasor in the same combined cloud
        XqA, QA = pc.mn_edge_cloud(mesh, r["M"].astype(complex))
        XqB, QB = pc.mn_edge_cloud(mB, MB.astype(complex))
        import radia.planar_charges as _pc
        phi = np.linspace(0, 2 * np.pi, 1440, endpoint=False)
        P = np.stack([1.6 * np.cos(phi), 1.6 * np.sin(phi)], axis=1)
        Hc = _pc.field_complex(np.vstack([XqA, XqB]), np.concatenate([QA, QB]), P)
        c, s = np.cos(phi), np.sin(phi)
        Hr = Hc[:, 0] * c + Hc[:, 1] * s
        H2 = Hc[:, 0] * np.conj(Hc[:, 0]) + Hc[:, 1] * np.conj(Hc[:, 1])
        Fx = np.real(Hr * np.conj(Hc[:, 0]) - 0.5 * H2 * c).sum()
        Fy = np.real(Hr * np.conj(Hc[:, 1]) - 0.5 * H2 * s).sum()
        F_cplx = MU0 * 1.6 * (2 * np.pi / 1440) * 0.5 * np.array([Fx, Fy])
    assert np.allclose(F_cplx, 0.5 * F_real, rtol=1e-6, atol=1e-9), (F_cplx, F_real)
