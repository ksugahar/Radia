"""Device-built Newton tangents of matching-trace mixed Omega against the CPU route.

Guards: the device tangent and its preconditioner stay consistent with the CPU linearization, through
the host-evaluated point tensor and its update on every Newton step, for the production entry
``solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin(..., Backend="gpu")``.
Assumes: the nonlinear residual stays the same host quadrature in both routes, so identical stopping
tests hold on identical inputs. SOLUTION_REL is a PROPOSED detector; the CPU direct solve is the
reference, not validation of the formulation. The linear matching-trace case keeps the tolerances of
the CPU iterative solve in ``test_kelvin_mixed_omega.py``.
"""
import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

from radia.mixed_omega_newton import solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin as solve

SOLUTION_REL = 1e-6
ITERATION_LIMIT = 8

pytestmark = pytest.mark.skipif(not os.environ.get("RADIA_FEM_NATIVE_TEST_DIR"), reason="the FEM device native lane is not selected")


@pytest.fixture(scope="module")
def case():
    spec = importlib.util.spec_from_file_location("mixed_helpers", Path(__file__).with_name("test_kelvin_mixed_omega.py"))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    from netgen.meshing import Element0D

    mesh, source, _, table = helper._picard_case(maxh=0.65)
    material = mesh.GetMaterials().index("total") + 1
    vertex = next(v for e in mesh.ngmesh.Elements3D() if e.index == material
                  for v in e.vertices if mesh.ngmesh.Points()[v].p[0] > 0.0)
    gauge_index = len(mesh.GetBBBoundaries()) + 1
    mesh.ngmesh.Add(Element0D(vertex, index=gauge_index))
    mesh.ngmesh.SetCD3Name(gauge_index, "GND")
    with ng.TaskManager():
        mesh = ng.Mesh(mesh.ngmesh)
    return mesh, source, table


def _run(case, order, **options):
    from radia.kelvin_solver import project_source_interface_potential

    mesh, source, table = case
    with ng.TaskManager():
        trace = project_source_interface_potential(mesh, source, "source_total_interface", order=order)["potential"]
    arguments = dict(bh_table=table, nonlinear_materials=("total",), reduced_materials=("reduced",), total_materials=("total",),
                     interface_boundary="source_total_interface", dirichlet_bbbnd="GND", kelvin_mats=(), tolerance=1e-6,
                     residual_tolerance=1e-9, order=order, condense_matching_trace=True)
    arguments.update(options)
    with ng.TaskManager():
        return solve(mesh, source, trace, 1.0, (3.0, 0.0, 0.0), **arguments)


@pytest.mark.parametrize("order", [1, 2])
def test_matching_trace_device_matches_the_direct_solve(order):
    """The linear entry of ``test_matching_trace_condensation_matches_multiplier_and_cg``, device against direct."""
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X
    from radia.kelvin_solver import solve_magnetostatic_matching_trace_total_reduced_omega

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "natural"
    reduced.faces.Min(X).name = "norm_boundary"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "natural"
    total.faces.Max(X).name = "norm_boundary"
    total.faces.Min(X).name = "source_total_interface"
    mesh = ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=0.65))
    source_h = ng.CoefficientFunction((ng.x * ng.x + 0.1 * ng.y, 0.2 * ng.y, -0.15 * ng.z))
    trace_space = ng.H1(mesh, order=order, definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    trace.Set(-ng.z, definedon=mesh.Boundaries("source_total_interface"))
    common = dict(mu_r_by_material={"reduced": 1.0, "total": 2.0}, reduced_materials=("reduced",), total_materials=("total",),
                  interface_boundary="source_total_interface", order=order, dirichlet_boundary="norm_boundary")
    with ng.TaskManager():
        direct = solve_magnetostatic_matching_trace_total_reduced_omega(mesh, source_h, trace, solver="direct", **common)
        device = solve_magnetostatic_matching_trace_total_reduced_omega(mesh, source_h, trace, solver="cg", cg_preconditioner="bddc_h1amg",
                                                                        cg_tolerance=1.0e-12, Backend="gpu", **common)
    assert device["backend"] == "gpu", f"backend {device['backend']}, expected gpu"
    assert device["linear_residual_relative"] < 1.0e-9, f"true residual {device['linear_residual_relative']}"
    for point in ((-0.65, 0.1, 0.15), (-0.2, -0.2, 0.25), (0.2, 0.15, -0.1), (0.65, -0.1, 0.2)):
        actual = np.asarray(device["H_cf"](mesh(*point)), dtype=float)
        expected = np.asarray(direct["H_cf"](mesh(*point)), dtype=float)
        assert np.allclose(actual, expected, rtol=1.0e-9, atol=1.0e-10), f"H at {point}: device {actual}, direct {expected}"


@pytest.mark.parametrize("order", [1, 2])
def test_device_tangent_reaches_the_cpu_newton_solution(case, order):
    reference = _run(case, order, linear_solver="direct")
    device = _run(case, order, linear_solver="cg", Backend="gpu")
    assert device["backend"] == "gpu", f"backend {device['backend']}, expected gpu"
    stats = device["nonlinear_stats"]
    assert stats["converged"], f"device Newton stopped: {stats}"
    assert len(stats["history"]) <= reference["nonlinear_stats"]["iterations"] + ITERATION_LIMIT, (
        f"device Newton took {len(stats['history'])} steps, CPU {reference['nonlinear_stats']['iterations']}")
    expected = reference["solution"].vec.FV().NumPy()
    error = np.linalg.norm(device["solution"].vec.FV().NumPy() - expected) / np.linalg.norm(expected)
    assert error <= SOLUTION_REL, f"device Newton solution off by {error}, limit {SOLUTION_REL}"
    assert stats["backend"] == "gpu", f"tangent backend {stats['backend']}, expected gpu"
    families = {row["payload"].split(".")[0] for row in stats["device_report"]["ledger"]}
    assert families <= {"mesh", "table", "coef", "rhs", "lift", "result", "scalar", "probe"}, f"ledger families {families}"
    assert all(row["linear_iterations"] > 0 for row in stats["history"]), f"linear iterations {stats['history']}"
