"""Production curved-geometry HDiv-VIM validation.

The old entry-by-entry test called public pybind probes backed by the same C++
kernel as the Gram entry, so it was not an independent reference.  The durable
contracts are end-to-end solve convergence and magnetic-moment accuracy against
the analytic linear sphere.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Sphere, Pnt, OCCGeometry  # noqa: E402
from radia.vim import Solve  # noqa: E402


def test_curved_demag_solve_runs_and_converges():
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=1.0))
    with ng.TaskManager():
        result = Solve(
            mesh, mu_r=100.0,
            H_ext=ng.CoefficientFunction((0, 0, 1.0)),
            order=1, curve_order=2)
    assert result["curve_order"] == 2
    assert result["iters"] < 400
    assert 0.20 < result["demag"] < 0.45
    assert abs(result["M_avg"][2]) > 1.0
    assert abs(result["M_avg"][0]) < 0.05
    assert abs(result["M_avg"][1]) < 0.05


def test_curved_moment_beats_flat():
    """P2 geometry recovers the analytic sphere moment at least 10x better."""
    radius, mu_r, applied = 1.0, 100.0, 1.0e4
    magnetization = 3.0 * (mu_r - 1.0) / (mu_r + 2.0) * applied
    exact_moment = magnetization * (4.0 / 3.0 * np.pi * radius ** 3)

    def moment_error(curved):
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), radius)).GenerateMesh(maxh=0.6))
        with ng.TaskManager():
            if curved:
                mesh.Curve(2)
            result = Solve(
                mesh, mu_r=mu_r,
                H_ext=ng.CoefficientFunction((0, 0, applied)),
                order=1, curve_order=2 if curved else None)
            volume = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))
        moment = np.asarray(result["M_avg"], float) * volume
        return float(np.linalg.norm(
            [moment[0], moment[1], moment[2] - exact_moment]) / exact_moment)

    flat_error = moment_error(False)
    curved_error = moment_error(True)
    assert curved_error < 5.0e-3
    assert curved_error < flat_error / 10.0
