"""Golden lock for PERMANENT-MAGNET (fixed-M) source regions in the 2D planar MMMM (radia.mmmm2d).

A hard PM has a RIGID magnetization; its field sources the soft-iron demag solve.  The PM field at
the iron centroids is the SHARED planar_charges.magnet_field (built on the shared exterior_field),
so the same one-way source pattern serves BOTH MMMM and the HDiv-VIM (add magnet_field to H_ext).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_charges as pc


def _disk(cx, a=1.0, maxh=0.13):
    geo = SplineGeometry(); geo.AddCircle((cx, 0.0), r=a, bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_magnet_source_equals_explicit_field():
    """solve(magnets=[(pm, M)]) == solve(H_ext = magnet_field at the iron centroids) -- the PM source
    injection is exactly an applied field (rigorous consistency)."""
    with ng.TaskManager():
        iron = _disk(2.0); mag = _disk(-2.0)
        Mm = np.tile([8.0e5, 0.0], (mag.ne, 1))
        rA = m2.solve_planar_demag(iron, mu_r=50.0, H_ext=(0.0, 0.0), magnets=[(mag, Mm)])
        _, _, centroids, _ = m2._extract_geometry(iron)
        Hmag = pc.magnet_field([(mag, Mm)], centroids)
        rB = m2.solve_planar_demag(iron, mu_r=50.0, H_ext=Hmag)
    assert np.allclose(rA["M"], rB["M"], rtol=1e-10, atol=1e-6), (rA["M_avg"], rB["M_avg"])


def test_pm_magnetizes_iron():
    """A magnet (M along +x) at x=-2 magnetizes a soft-iron disk at x=+2 -- the iron M is nonzero and
    points along the local magnet field (+x on the axis)."""
    with ng.TaskManager():
        iron = _disk(2.0); mag = _disk(-2.0)
        Mm = np.tile([8.0e5, 0.0], (mag.ne, 1))
        r = m2.solve_planar_demag(iron, mu_r=200.0, H_ext=(0.0, 0.0), magnets=[(mag, Mm)])
    Mx, My = r["M_avg"]
    assert Mx > 0.0 and abs(My) < 0.1 * abs(Mx), (Mx, My)


def test_pm_attracts_iron():
    """The magnet attracts the magnetized iron: the force on the iron points TOWARD the magnet (-x)."""
    with ng.TaskManager():
        iron = _disk(2.0); mag = _disk(-2.0)
        Mm = np.tile([8.0e5, 0.0], (mag.ne, 1))
        r = m2.solve_planar_demag(iron, mu_r=200.0, H_ext=(0.0, 0.0), magnets=[(mag, Mm)])
        F = m2.force_between([(iron, r["M"]), (mag, Mm)], Rc=1.6, center=(2.0, 0.0))
    assert F[0] < 0.0 and abs(F[1]) < 0.1 * abs(F[0]), F      # pulled toward the magnet at x=-2
