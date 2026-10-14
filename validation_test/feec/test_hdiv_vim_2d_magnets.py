"""Golden: the HDiv-VIM 2D planar layer adopts the SHARED separate-body permanent-magnet source
(radia.vim.PlanarSolve / vim.Solve magnets=, via the shared planar_charges.magnet_field_cf).

A rigid PM disk magnetises a nearby soft-iron disk. Also: no magnet + no applied field -> M ~ 0.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

from radia.vim import Solve

MU0 = 4e-7 * np.pi
A = 1.0
MREM = 8.0e5


def _disk(cx, maxh=0.16):
    geo = SplineGeometry(); geo.AddCircle((cx, 0.0), r=A, bc="e")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_hdiv_vim_magnets_magnetize_nearby_iron():
    """HDiv-VIM magnets= magnetizes the iron toward a nearby permanent magnet."""
    iron_c, mag_c = 3.0 * A, -3.0 * A
    with ng.TaskManager():
        iron_h = _disk(iron_c); mag_h = _disk(mag_c)
        Mm_h = np.tile([MREM, 0.0], (mag_h.ne, 1))
        rH = Solve(iron_h, 200.0, ng.CoefficientFunction((0.0, 0.0)),
                              magnets=[(mag_h, Mm_h)])
    MH = rH["M_avg"]
    assert MH[0] > 0 and abs(MH[1]) < 0.1 * abs(MH[0]), MH        # magnetised toward the magnet (+x)
    assert 1.0e3 < MH[0] < 5.0e5, MH


def test_hdiv_vim_no_source_gives_zero():
    """No applied field and no magnet -> the soft iron is unmagnetised."""
    with ng.TaskManager():
        iron = _disk(0.0)
        r = Solve(iron, 200.0, ng.CoefficientFunction((0.0, 0.0)))
    assert np.linalg.norm(r["M_avg"]) < 1e-3 * MREM, r["M_avg"]


def test_hdiv_vim_magnets_3d_rejected():
    """magnets= is 2D-only in vim.Solve -- a 3D mesh must fail loud."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, OCCGeometry
    with ng.TaskManager():
        m3 = ng.Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=0.5))
        with pytest.raises(NotImplementedError):
            Solve(m3, 100.0, ng.CoefficientFunction((1.0, 0.0, 0.0)),
                             magnets=[(m3, np.zeros((m3.ne, 2)))])
