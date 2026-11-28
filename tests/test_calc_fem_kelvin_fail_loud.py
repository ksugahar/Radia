"""Fail-loud contracts for FEM-Kelvin solver mode and Karl acceptance."""

import sys
from pathlib import Path

import pytest

PANELS = Path(__file__).parents[1] / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))
import calc_fem_kelvin as solver  # noqa: E402


def _validate(**overrides):
    arguments = dict(
        solver="sparsecholesky", formulation="total", max_iter=4,
        esim_per_panel=False, periodic_kelvin=False)
    arguments.update(overrides)
    solver._validate_mode_contracts(**arguments)


@pytest.mark.parametrize("overrides, message", [
    ({"max_iter": 0}, "positive"),
    ({"esim_per_panel": True}, "local-Z surface loss"),
    ({"solver": "iccg", "periodic_kelvin": True}, "periodic Kelvin"),
    ({"solver": "bddc", "formulation": "scattered"}, "only for.*sparsecholesky"),
])
def test_unsupported_modes_fail_before_solving(overrides, message):
    with pytest.raises(ValueError, match=message):
        _validate(**overrides)


def test_supported_total_and_direct_scattered_modes_remain_available():
    _validate(solver="bddc", formulation="total")
    _validate(solver="sparsecholesky", formulation="scattered",
              periodic_kelvin=True)


def test_karl_acceptance_publishes_impedance_used_by_final_field():
    converged, published = solver._accepted_karl_impedance(
        iteration=2, max_iter=5, relative_change=1.0e-4,
        tolerance=1.0e-3, impedance_used_for_solve=2.0 + 3.0j)
    assert converged
    assert published == 2.0 + 3.0j

    converged, published = solver._accepted_karl_impedance(
        iteration=4, max_iter=5, relative_change=2.0e-3,
        tolerance=1.0e-3, impedance_used_for_solve=2.0 + 3.0j)
    assert not converged
    assert published is None


def test_fem_cell_convergence_gate_keeps_diagnostics():
    with pytest.raises(
            RuntimeError,
            match=r"FEM scalar ESIM.*iterations=9.*relative_change=0.02"):
        solver.require_esim_converged(
            {"converged": False, "iterations": 9,
             "relative_change": 0.02, "Z": 1.0 + 1.0j},
            "FEM scalar ESIM cell solve")


def test_unknown_backend_fails_before_any_mesh_work():
    with pytest.raises(ValueError, match="Backend must be"):
        solver.solve_fem(vol_file="not_read.vol", Backend="tpu")


def test_periodic_kelvin_dc_peec_checks_all_free_rows(tmp_path, monkeypatch):
    """Guards: the edge forest and DC PEEC inductance both preserve the solve."""
    import math

    import numpy as np
    from netgen.occ import Glue, IdentificationType, OCCGeometry, Pnt, Sphere, Vertex
    from ngsolve import Mesh, TaskManager

    from radia import coil_from_cad
    from radia.peec_bundle import build_bundle_solver

    inner = Sphere(Pnt(0, 0, 0), 0.06)
    inner.mat("air")
    inner.faces.name = "kelvin_int"
    center = Pnt(0.18, 0, 0)
    outer = Sphere(center, 0.06)
    outer.mat("kelvin")
    outer.faces.name = "kelvin_ext"
    ground = Vertex(center)
    ground.name = "GND"
    geo = Glue([inner, outer, ground])
    faces = {face.name: face for face in geo.faces}
    faces["kelvin_int"].Identify(faces["kelvin_ext"], "kelvin", IdentificationType.PERIODIC)
    vol = tmp_path / "periodic.vol"
    with TaskManager():
        mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.018))
        assert mesh.ngmesh.GetNrIdentifications() == 1
        mesh.ngmesh.Save(str(vol))

    # A single open polygonal turn has one ampere in each series branch.
    # Substitute only the CAD extraction so this test needs no Cubit.
    theta = np.linspace(0, 2 * math.pi, 21, endpoint=False)
    points = np.column_stack((0.025 * np.cos(theta), 0.025 * np.sin(theta), np.zeros_like(theta)))
    paths = [list(zip(points[:-1], points[1:]))]
    circuit, segments, _, _ = build_bundle_solver(paths, 0.001, 0.001, 5.8e7)
    step = tmp_path / "synthetic.step"
    step.touch()
    monkeypatch.setattr(coil_from_cad, "filaments_from_step", lambda *_args, **_kwargs: {
        "filament_paths": paths,
        "seg_of_filament": segments,
        "solver": circuit,
    })

    result = solver.solve_fem(vol_file=str(vol), fes_order=1, frequency=0,
                              peec_step=str(step), reg=1e-6, nthreads=1)
    assert result["kelvin_gradient_gauge_dofs"] > 0
    residual = result["linear_true_relative_residual"]
    assert residual <= solver.LINEAR_TRUE_RESIDUAL_LIMIT, residual
    expected = float(np.sum(circuit.L))
    assert result["L_peec"] == pytest.approx(expected), (result["L_peec"], expected)
