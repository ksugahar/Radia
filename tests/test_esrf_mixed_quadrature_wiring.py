"""Exercise actual adapter/caller AST with doubles, without native execution."""
from __future__ import annotations

import argparse
import ast
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
ADAPTER = ROOT / "validation_test/c_type_three_engine/run_three_engine.py"
RUNNER = ROOT / "validation_test/esrf_three_engine/run_coil_yoke_three_engine.py"
HYBRID_RUNNER = ROOT / "validation_test/esrf_three_engine/run_hybrid_undulator_three_engine.py"


def compile_functions(nodes, namespace):
    tree = ast.Module(body=[ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0), *nodes], type_ignores=[])
    exec(compile(ast.fix_missing_locations(tree), "<reviewed-local-caller>", "exec"), namespace)  # noqa: S102 - local reviewed AST


@pytest.mark.parametrize("method,order", [("picard", 1), ("picard", 2), ("newton", 2)])
def test_hybrid_material_controls_reach_solver_and_checkpoint(method, order):
    calls = {}
    def solve(*args, **kwargs):
        calls.update(kwargs)
        return {"B_cf": object(), "fes": SimpleNamespace(ndof=10),
                "nonlinear_stats": {"converged": True},
                "static_electromagnet_contract": {"source_trace": {}}}
    tree = ast.parse(HYBRID_RUNNER.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                and n.name == "_solve_mixed_omega")
    ns = {"np": np, "ng": SimpleNamespace(TaskManager=nullcontext),
          "time": SimpleNamespace(perf_counter=lambda: 0),
          "MIXED_DOMAIN": object(), "MIXED_DOMAIN_LABEL": "mixed Omega",
          "solve_static_electromagnet_mixed_total_reduced_omega": solve,
          "_evaluate_cf": lambda *args: np.ones((2, 3))}
    compile_functions([node], ns)
    bonus = 8 if method == "newton" else None
    _, diag = ns["_solve_mixed_omega"](
        SimpleNamespace(ne=2, nv=8), SimpleNamespace(field_cf=object()), "BH", None,
        order=order, kelvin_radius=.18, kelvin_center=(0, 0, 0),
        nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80,
        nonlinear_relaxation=.2, source_potential_tolerance=.05,
        nonlinear_method=method, bonus_intorder=12, material_bonus_intorder=bonus)
    assert calls["nonlinear_method"] == method
    assert calls["nonlinear_material_bonus_intorder"] == bonus
    assert calls["bonus_intorder"] == 12
    assert calls["source_potential_contract"] == "global_physical"
    assert calls["nonlinear_material_sampling"] == (
        "integration_point" if method == "newton" else "element_centroid")
    assert calls["nonlinear_material_update_order"] == (
        1 if method == "picard" and order > 1 else None)
    assert diag["nonlinear_method"] == method
    # A changed material rule must invalidate previous field checkpoints.
    contracts = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == "_checkpoint_contract"
                 and any(k.arg == "assets" for k in n.keywords)]
    assert len(contracts) == 1
    settings = {k.arg: ast.unparse(k.value) for k in contracts[0].keywords}
    for key in ("mixed_method", "mixed_bonus", "mixed_material_bonus"):
        assert settings[key] == f"options.{key}"


@pytest.mark.parametrize("exact,source_order,bonus,missing", [(False, None, 4, False), (True, 3, 8, False), (True, 4, 12, False), (True, 3, 8, True)])
@pytest.mark.parametrize("method", ["picard", "newton"])
def test_actual_adapter_arguments_and_diagnostics(exact, source_order, bonus, missing, method):
    calls = {}
    token = object()
    def kelvin(*args):
        calls["kelvin"] = args
        return None if missing else token
    def solve(*args, **kwargs):
        calls["solver"] = kwargs
        source = {"projection_order": source_order or 2, "iron_relative_harmonic_norm": .01}
        if not exact:
            source["kelvin_relative_tangential_residual"] = .001
        return {"B_cf": object(), "static_electromagnet_contract": {"source_trace": source},
                "fes": SimpleNamespace(ndof=10), "nonlinear_stats": {"converged": True}}
    ns = {"np": np, "ng": SimpleNamespace(TaskManager=nullcontext),
          "time": SimpleNamespace(perf_counter=lambda: 0),
          "rad": SimpleNamespace(RadiaField=lambda *args: "source", KelvinRadiaFieldStrength=kelvin),
          "MIXED_DOMAIN": object(), "solve_static_electromagnet_mixed_total_reduced_omega": solve,
          "evaluate_cf": lambda *args: np.ones((2, 3))}
    node = next(n for n in ast.parse(ADAPTER.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == "solve_omega")
    compile_functions([node], ns)
    if missing:
        with pytest.raises(RuntimeError, match="not constructed"):
            ns["solve_omega"](None, 7, None, nonlinear=True, order=2,
                nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80,
                nonlinear_verbose=False, kelvin_center=(0, 0, 0), kelvin_radius=.16,
                points=None, source_trace_tolerance=.05, exact_exterior_source=True)
        assert "solver" not in calls
        return
    _, diag = ns["solve_omega"](SimpleNamespace(ne=2, nv=8), 7, "BH", nonlinear=True,
        order=2, nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80,
        nonlinear_verbose=False, kelvin_center=(1., 2., 3.), kelvin_radius=.16,
        points=np.ones((2, 3)), source_trace_tolerance=.05,
        source_projection_order=source_order, bonus_intorder=bonus, exact_exterior_source=exact,
        nonlinear_method=method, material_bonus_intorder=8 if method == "newton" else None,
        source_load="surface_flux")
    assert calls["solver"]["source_projection_order"] == source_order
    assert calls["solver"]["bonus_intorder"] == bonus
    assert calls["solver"]["nonlinear_method"] == method
    assert calls["solver"]["nonlinear_material_bonus_intorder"] == (8 if method == "newton" else None)
    assert calls["solver"]["reduced_source_load"] == "surface_flux"
    assert calls["solver"]["total_source_load"] == "surface_flux"
    assert calls["solver"]["nonlinear_material_update_order"] == (1 if method == "picard" else None)
    assert calls["solver"]["kelvin_source_h"] is (token if exact else None)
    assert diag["bonus_intorder"] == bonus
    if exact:
        assert calls["kelvin"] == (7, (1., 2., 3.), .16, (0., 0., 0.))
        assert diag["source_trace"]["kelvin_relative_tangential_residual"] is None
    else:
        assert "kelvin" not in calls
        assert diag["source_trace"]["kelvin_relative_tangential_residual"] == .001


@pytest.mark.parametrize("source,bonus,exact", [(0, 4, False), (True, 4, False),
                                              (None, -1, False), (None, 4, 1)])
def test_adapter_invalid_options_fail_before_native(source, bonus, exact):
    node = next(n for n in ast.parse(ADAPTER.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == "solve_omega")
    ns = {"time": SimpleNamespace(perf_counter=lambda: 0)}
    compile_functions([node], ns)
    with pytest.raises(ValueError):
        ns["solve_omega"](None, 1, None, nonlinear=True, order=2, nonlinear_tolerance=2e-5,
                         nonlinear_maximum_iterations=80, nonlinear_verbose=False,
                         kelvin_center=(0, 0, 0), kelvin_radius=.16, points=None,
                         source_trace_tolerance=.05, source_projection_order=source,
                         bonus_intorder=bonus, exact_exterior_source=exact)


def test_actual_runner_call_and_checkpoint_settings():
    tree = ast.parse(RUNNER.read_text(encoding="utf-8"))
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    caller = next(n for n in main.body if isinstance(n, ast.FunctionDef) and n.name == "run_mixed")
    options = SimpleNamespace(fem_order=2, nonlinear_tolerance=2e-5, source_trace_tolerance=.05,
                              mixed_relaxation=.3, mixed_anderson_depth=0, mixed_source_order=3,
                              mixed_bonus=12, mixed_exact_exterior_source=True,
                              mixed_method="newton", mixed_material_bonus=8,
                              mixed_source_load="surface_flux")
    calls = {}
    def solve(*args, **kwargs):
        calls.update(kwargs)
        return "result"
    ns = {"engines": SimpleNamespace(solve_omega=solve), "options": options, "fem_mesh": None,
          "coil": 7, "bh_table": "BH", "mixed_nonlinear_maximum_iterations": 80,
          "kelvin_center": (0, 0, 0), "case": SimpleNamespace(kelvin_radius_m=.16),
          "field_points": np.ones((2, 3)), "np": np}
    compile_functions([caller], ns)
    assert ns["run_mixed"](None) == "result"
    assert calls["source_projection_order"] == 3
    assert calls["nonlinear_method"] == "newton"
    assert calls["material_bonus_intorder"] == 8
    assert calls["source_load"] == "surface_flux"
    assert calls["bonus_intorder"] == 12 and calls["exact_exterior_source"] is True
    assignment = next(n for n in main.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "engine_settings" for t in n.targets))
    ns["reduced_a_settings"] = {}
    compile_functions([assignment], ns)
    settings = ns["engine_settings"]["mixed_total_reduced_omega"]
    assert settings["nonlinear_method"] == "newton"
    assert settings["material_bonus_intorder"] == 8
    assert settings["source_load"] == "surface_flux"
    assert settings["source_projection_order"] == 3
    assert settings["bonus_intorder"] == 12 and settings["exact_exterior_source"] is True


@pytest.mark.parametrize("extra,expected", [([], (None, 4, False)),
    (["--mixed-source-order", "3", "--mixed-bonus", "12", "--mixed-exact-exterior-source"], (3, 12, True))])
def test_actual_cli_defaults_and_explicit_settings(monkeypatch, extra, expected):
    tree = ast.parse(RUNNER.read_text(encoding="utf-8"))
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    validator = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_validated_tolerance")
    prefix = []
    for node in main.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "case" for t in node.targets):
            break
        prefix.append(node)
    ns = {"argparse": argparse, "Path": Path, "np": np, "argv": [
        "--case", "6", "--assets-dir", "unused", "--fem-mesh", "unused.vol",
        "--fem-mesh-report", "unused.json", "--output", "unused-out.json", *extra],
        "__doc__": "test"}
    compile_functions([validator, *prefix], ns)
    options = ns["options"]
    assert (options.mixed_source_order, options.mixed_bonus, options.mixed_exact_exterior_source) == expected
