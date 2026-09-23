"""Contracts of the first-order finite-box lane (validation_test/c_type_p1_ams_box).

The lane's reduced-A engine is exercised on a small synthetic mesh with the
same label contract as the C-type box mesh: an iron cube and a dummy coil
block inside an air box, ``outer`` faces, ``iron_air_interface`` faces and a
``GND`` vertex.  The tests pin what the benchmark relies on:

* the tabulated flux-side law equals the production inverse of the PCHIP
  B(H) law and is positive and monotone;
* the compiled AMS + CG linear solve reaches the same field as the direct
  solve and honours its true-residual contract;
* Newton on a linear law is one exact step, and on the sample B-H table it
  converges with full steps;
* the mesh contract rejects a mesh that does not carry the lane's labels.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

LANE = Path(__file__).resolve().parents[1] / "validation_test" / "c_type_p1_ams_box"
MU0 = 4.0e-7 * math.pi
_MODULE = None


def lane():
    global _MODULE
    if _MODULE is None:
        spec = importlib.util.spec_from_file_location("c_type_p1_box_lane", LANE / "run_p1_box.py")
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _MODULE = module
    return _MODULE


def _bh_table():
    return [[0.0, 0.0], [100.0, 0.6], [400.0, 1.2], [2000.0, 1.6], [20000.0, 1.9]]


@pytest.fixture(scope="module")
def box_mesh():
    from netgen.occ import Box, Glue, OCCGeometry, Pnt

    iron = Box(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2))
    iron.mat("iron")
    iron.faces.name = "iron_air_interface"
    coil = Box(Pnt(0.5, -0.1, -0.1), Pnt(0.7, 0.1, 0.1))
    coil.mat("coil")
    coil.faces.name = "coil_air_interface"
    outer = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))
    outer.faces.name = "outer"
    outer.mat("air")
    air = outer - iron - coil
    for vertex in air.vertices:
        if all(abs(vertex.p[k] + 1.0) < 1e-9 for k in range(3)):
            vertex.name = "GND"
    mesh = ng.Mesh(OCCGeometry(Glue([iron, coil, air])).GenerateMesh(maxh=0.35))
    assert "GND" in set(mesh.GetBBBoundaries())
    return mesh


def _engine(mesh, solver, source=(0.0, 0.0, 0.8), **overrides):
    module = lane()
    settings = dict(linear_solver=solver, cg_tolerance=1.0e-9, cg_max_iterations=1000,
                    ams_num_smooth=1, source_projection_order=2)
    settings.update(overrides)
    return module.ReducedAP1Box(mesh, ng.CoefficientFunction(tuple(source)), **settings)


def _points():
    return np.asarray([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.0, 0.0, 0.3], [0.5, 0.5, 0.5]])


def test_the_flux_side_law_matches_the_production_inverse():
    from radia.vector_potential_solver import _build_nu_of_b_interpolator

    law = lane().SoftIronLaw(_bh_table())
    nu_of_b, b_sat = _build_nu_of_b_interpolator(_bh_table())
    assert law.b_saturation == b_sat
    samples = np.array([0.05, 0.3, 0.9, 1.3, 1.7, 2.5])
    expected = np.array([nu_of_b(b) for b in samples])
    np.testing.assert_allclose(law.reluctivity(samples), expected, rtol=1.0e-4)
    assert np.all(law.differential_reluctivity(samples) > 0.0)
    # Beyond the table the law continues with the vacuum slope.
    assert law.differential_reluctivity(np.array([5.0]))[0] == pytest.approx(1.0 / MU0, rel=0.05)
    with pytest.raises(ValueError):
        lane().SoftIronLaw([[0.0, 0.0]])


@pytest.mark.parametrize("smoothing_sweeps", [1, 2])
def test_ams_and_direct_linear_solves_agree(box_mesh, smoothing_sweeps):
    ams = _engine(box_mesh, "ams", ams_num_smooth=smoothing_sweeps)
    direct = _engine(box_mesh, "direct")
    with ng.TaskManager():
        pass
    field_ams, stats_ams, _ = ams.run_linear(1000.0, _points())
    field_direct, _, _ = direct.run_linear(1000.0, _points())
    assert np.all(np.isfinite(field_ams))
    scale = float(np.max(np.linalg.norm(field_direct, axis=1)))
    assert np.max(np.linalg.norm(field_ams - field_direct, axis=1)) / scale < 1.0e-6
    row = stats_ams["history"][0]
    assert row["relative_residual"] <= 1.0e-9
    assert row["cg_iterations"] > 0
    # Iron with mu_r = 1000 concentrates the flux: the centre field exceeds the source.
    assert field_direct[0, 2] > 0.8


def test_ungauged_shifted_iccg_matches_the_gauged_ams_field(box_mesh):
    """The singular curl-curl system with a consistent right-hand side, solved by
    Radia's shifted incomplete-Cholesky CG, gives the same B as the gauged AMS
    solve (A differs by a gradient, B does not)."""
    module = lane()
    with pytest.raises(ValueError, match="requires linear_solver='iccg'"):
        _engine(box_mesh, "ams", gauge_epsilon=0.0)
    iccg = _engine(box_mesh, "iccg", gauge_epsilon=0.0, cg_tolerance=1.0e-9)
    field_iccg, stats, _ = iccg.run_linear(1000.0, _points())
    field_ams, _, _ = _engine(box_mesh, "ams").run_linear(1000.0, _points())
    row = stats["history"][0]
    assert row["relative_residual"] <= 1.0e-9
    assert row["cg_iterations"] > 0
    scale = float(np.max(np.linalg.norm(field_ams, axis=1)))
    assert np.max(np.linalg.norm(field_iccg - field_ams, axis=1)) / scale < 1.0e-5
    assert iccg.describe()["gauge"].startswith("nograds space only")
    assert module.GAUGE_EPSILON > 0.0


def test_newton_on_a_linear_law_is_one_exact_step(box_mesh):
    module = lane()
    engine = _engine(box_mesh, "ams")
    mu_r = 500.0
    # A "table" that is exactly linear up to far beyond the fields reached here.
    table = [[0.0, 0.0], [10.0, 10.0 * MU0 * mu_r], [1.0e4, 1.0e4 * MU0 * mu_r]]
    law = module.SoftIronLaw(table)
    field_newton, stats, _ = engine.run_newton(
        law, newton_tolerance=1.0e-8, tolerance=1.0e-6, max_iterations=10,
        max_halvings=4, observation=_points())
    assert stats["converged"]
    assert stats["iterations"] <= 2
    assert all(row["step_length"] == 1.0 for row in stats["history"])
    field_linear, _, _ = _engine(box_mesh, "ams").run_linear(mu_r, _points())
    scale = float(np.max(np.linalg.norm(field_linear, axis=1)))
    assert np.max(np.linalg.norm(field_newton - field_linear, axis=1)) / scale < 1.0e-5


@pytest.mark.parametrize("shift", [1e-2, 1e-6])
@pytest.mark.parametrize("project", [False, True])
def test_shifted_ams_preserves_the_ungauged_operator_and_field(box_mesh, shift, project):
    ams = _engine(box_mesh, "ams", gauge_epsilon=0., ams_preconditioner_shift=shift,
                  cg_tolerance=1e-7, ams_project_gradients=project)
    reference = _engine(box_mesh, "iccg", gauge_epsilon=0., cg_tolerance=1e-7)
    field, stats, _ = ams.run_linear(1000., _points())
    expected, _, _ = reference.run_linear(1000., _points())
    np.testing.assert_allclose(field, expected, rtol=1e-5, atol=1e-7)
    row = stats['history'][0]
    assert row['relative_residual'] <= 1e-7
    assert row['cg_restarts'] == 0
    # Matrix setup must not overwrite the operator that CG actually solves.
    a, _ = ams._picard_forms()
    with ng.TaskManager():
        a.Assemble()
    original = a.mat.AsVector().FV().NumPy().copy()
    ams._ams_prepare(a.mat)
    np.testing.assert_array_equal(a.mat.AsVector().FV().NumPy(), original)
    assert not np.array_equal(ams._ams_shift_form.mat.AsVector().FV().NumPy(), original)
    if project:
        pre = ams._ams_prepare(a.mat)
        gradient, _ = ams._gradient_projection
        x = np.random.default_rng(18).normal(size=ams.fes.ndof)
        projected = pre.project(x)
        assert np.linalg.norm(gradient.T @ projected) < 1e-10 * np.linalg.norm(x)
        np.testing.assert_allclose(pre.project(projected), projected, atol=1e-10)


def test_shifted_ams_fails_closed_on_invalid_shift_and_iteration_limit(box_mesh):
    for shift in (-1., float('nan'), float('inf')):
        with pytest.raises(ValueError, match='ams_preconditioner_shift'):
            _engine(box_mesh, 'ams', gauge_epsilon=0., ams_preconditioner_shift=shift)
    with pytest.raises(ValueError, match='requires AMS'):
        _engine(box_mesh, 'iccg', gauge_epsilon=0., ams_preconditioner_shift=.01)
    with pytest.raises(RuntimeError, match='true relative residual'):
        _engine(box_mesh, 'ams', gauge_epsilon=0., ams_preconditioner_shift=.01,
                cg_max_iterations=1).run_linear(1000., _points())


def test_newton_accepts_an_initially_zero_residual(box_mesh):
    engine = _engine(box_mesh, "ams", source=(0.0, 0.0, 0.0))
    field, stats, _ = engine.run_newton(
        lane().SoftIronLaw(_bh_table()), newton_tolerance=1e-8,
        tolerance=1e-6, max_iterations=10, max_halvings=4, observation=_points())
    assert stats["converged"] and stats["iterations"] == 0
    assert stats["final_residual_relative"] == 0.0
    np.testing.assert_array_equal(field, 0.0)


def test_beta_zero_ams_preserves_the_ungauged_field(box_mesh):
    ams = _engine(box_mesh, "ams", gauge_epsilon=0., ams_beta_zero=True,
                  cg_tolerance=1e-7)
    reference = _engine(box_mesh, "iccg", gauge_epsilon=0., cg_tolerance=1e-7)
    field, stats, _ = ams.run_linear(1000., _points())
    expected, _, _ = reference.run_linear(1000., _points())
    assert ams._ams.beta_zero is True
    assert ams._ams_shift_form is None and ams._gradient_projection is None
    assert stats['history'][0]['relative_residual'] <= 1e-7
    np.testing.assert_allclose(field, expected, rtol=1e-5, atol=1e-7)
    for kwargs in ({"gauge_epsilon": 1e-6}, {"ams_preconditioner_shift": .01}):
        with pytest.raises(ValueError, match="ams_beta_zero"):
            _engine(box_mesh, "ams", ams_beta_zero=True,
                    **{"gauge_epsilon": 0., **kwargs})


def test_cli_preserves_runtime_failure_as_json(tmp_path, monkeypatch):
    module = lane()
    output = tmp_path / "failure.json"
    monkeypatch.setattr(sys, "argv", ["run_p1_box.py", "--vol", str(tmp_path / "missing.vol"),
                                     "--output", str(output)])
    def fail(options):
        raise RuntimeError("native backend unavailable")
    monkeypatch.setattr(module, "run", fail)
    with pytest.raises(RuntimeError, match="native backend unavailable"):
        module.main()
    report = json.loads(output.read_text())
    assert not report["completed"] and not report["passed"]
    assert report["error"]["message"] == "native backend unavailable"


def test_cli_returns_failure_for_nonconvergence(tmp_path, monkeypatch):
    module = lane()
    monkeypatch.setattr(sys, "argv", ["run_p1_box.py", "--vol", "unused.vol",
                                     "--output", str(tmp_path / "result.json")])
    monkeypatch.setattr(module, "run", lambda options: {"completed": True, "passed": False})
    with pytest.raises(SystemExit) as exc:
        module.main()
    assert exc.value.code == 1


def test_newton_converges_on_the_sample_table_with_full_steps(box_mesh):
    module = lane()
    engine = _engine(box_mesh, "ams")
    law = module.SoftIronLaw(_bh_table())
    field, stats, _ = engine.run_newton(
        law, newton_tolerance=1.0e-8, tolerance=1.0e-5, max_iterations=30,
        max_halvings=6, observation=_points())
    assert stats["converged"]
    assert stats["final_residual_relative"] <= 1.0e-8
    assert stats["history"][-1]["step_length"] == 1.0
    assert not any(row.get("line_search_failed") for row in stats["history"])
    assert np.all(np.isfinite(field))
    # The 0.8 T source drives the cube past the knee: part of the iron carries a
    # reluctivity well above the initial one, so the converged material state is
    # genuinely nonlinear rather than the initial linear guess.
    assert stats["history"][-1]["nu_max"] > 5.0 * law.nu_initial
    mu_r_initial = 1.0 / (MU0 * law.nu_initial)
    field_linear, _, _ = _engine(box_mesh, "ams").run_linear(mu_r_initial, _points())
    assert not np.allclose(field, field_linear, rtol=1.0e-4, atol=0.0)


def test_picard_and_newton_agree_on_a_mildly_saturated_cube(box_mesh):
    """A 0.1 T source keeps the cube on the shoulder of the table, where the
    damped Picard map contracts; both loops must then land on one field."""
    module = lane()
    law = module.SoftIronLaw(_bh_table())
    source = (0.0, 0.0, 0.1)
    field_newton, newton_stats, _ = _engine(box_mesh, "ams", source).run_newton(
        law, newton_tolerance=1.0e-9, tolerance=1.0e-6, max_iterations=30,
        max_halvings=6, observation=_points())
    assert newton_stats["converged"]
    field_picard, picard_stats, _ = _engine(box_mesh, "ams", source).run_picard(
        law, relax=0.3, anderson_depth=3, tolerance=1.0e-6, max_iterations=120,
        observation=_points())
    assert picard_stats["converged"]
    assert picard_stats["iterations"] > newton_stats["iterations"]
    scale = float(np.max(np.linalg.norm(field_newton, axis=1)))
    assert np.max(np.linalg.norm(field_picard - field_newton, axis=1)) / scale < 1.0e-4


def test_the_mesh_contract_rejects_a_mesh_without_the_lane_labels(tmp_path, box_mesh):
    module = lane()
    from netgen.occ import Box, OCCGeometry, Pnt

    plain = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.5))
    vol = tmp_path / "plain.vol"
    plain.ngmesh.Save(str(vol))
    with pytest.raises(FileNotFoundError):
        module.check_mesh_contract(vol, plain)
    vol.with_suffix(".json").write_text(json.dumps({"vol_sha256": module.sha256(vol)}))
    with pytest.raises(RuntimeError, match="materials"):
        module.check_mesh_contract(vol, plain)
    good = tmp_path / "box.vol"
    box_mesh.ngmesh.Save(str(good))
    good.with_suffix(".json").write_text(json.dumps({"vol_sha256": "0" * 64}))
    with pytest.raises(RuntimeError, match="differs"):
        module.check_mesh_contract(good, box_mesh)
    good.with_suffix(".json").write_text(json.dumps({"vol_sha256": module.sha256(good)}))
    contract = module.check_mesh_contract(good, box_mesh)
    assert contract["vol_sha256"] == module.sha256(good)


def test_inexact_newton_preserves_field_and_outer_gates(box_mesh):
    law = lane().SoftIronLaw(_bh_table())
    options = dict(newton_tolerance=1e-8, tolerance=1e-5, max_iterations=30,
                   max_halvings=6, observation=_points())
    reference, fixed_stats, _ = _engine(box_mesh, "ams").run_newton(law, **options)
    engine = _engine(box_mesh, "ams", ams_update_every=2)
    field, stats, _ = engine.run_newton(law, inexact_linear=True, **options)
    assert fixed_stats["converged"] and stats["converged"]
    assert stats["final_residual_relative"] <= options["newton_tolerance"]
    assert stats["final_relative_change"] <= options["tolerance"]
    assert engine.cg_tolerance == 1e-9
    assert stats["inexact_linear"] is True
    for row in stats["history"]:
        assert 1e-9 <= row["linear_tolerance"] <= 0.01
        assert row["relative_residual"] <= row["linear_tolerance"]
    assert stats["history"][-1]["linear_tolerance"] < stats["history"][0]["linear_tolerance"]
    np.testing.assert_allclose(field, reference, rtol=1e-6, atol=1e-8)


def test_inexact_newton_restores_tolerance_after_linear_failure(box_mesh, monkeypatch):
    engine = _engine(box_mesh, "ams")
    def fail(*args, **kwargs):
        assert engine.cg_tolerance == 0.01
        raise RuntimeError("injected linear failure")
    monkeypatch.setattr(engine, "_solve_linear", fail)
    with pytest.raises(RuntimeError, match="injected linear failure"):
        engine.run_newton(lane().SoftIronLaw(_bh_table()), inexact_linear=True,
                          newton_tolerance=1e-8, tolerance=1e-5,
                          max_iterations=30, max_halvings=6, observation=_points())
    assert engine.cg_tolerance == 1e-9


def test_periodic_true_residual_checks_preserve_linear_field(box_mesh):
    engine = _engine(box_mesh, "ams", cg_check_interval=4)
    field, stats, _ = engine.run_linear(1000., _points())
    expected, _, _ = _engine(box_mesh, "direct").run_linear(1000., _points())
    np.testing.assert_allclose(field, expected, rtol=1e-6, atol=1e-8)
    row = stats["history"][0]
    assert row["relative_residual"] <= engine.cg_tolerance
    assert row["true_residual_checks"] < row["cg_iterations"]
    for interval in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            _engine(box_mesh, "ams", cg_check_interval=interval)
    with pytest.raises(ValueError, match="requires AMS"):
        _engine(box_mesh, "direct", cg_check_interval=4)
    with pytest.raises(RuntimeError, match="true relative residual"):
        _engine(box_mesh, "ams", cg_check_interval=100,
                cg_max_iterations=2).run_linear(1000., _points())


def test_native_ams_diagnostics_preserve_the_linear_solution(box_mesh):
    engine = _engine(box_mesh, "ams", ams_print_level=1)
    field, stats, _ = engine.run_linear(1000., _points())
    expected, _, _ = _engine(box_mesh, "ams").run_linear(1000., _points())
    np.testing.assert_allclose(field, expected, rtol=1e-8, atol=1e-9)
    assert stats["history"][0]["relative_residual"] <= engine.cg_tolerance
    assert engine.describe()["ams_print_level"] == 1
    with pytest.raises(ValueError, match="ams_print_level"):
        _engine(box_mesh, "ams", ams_print_level=2)
