"""Bound point-location work independently of mesh or solver dependencies."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest


def quadrature():
    source = Path(__file__).parents[1] / "src/radia/panels/calc_fem_kelvin.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                and n.name == "_element_line_quadrature")
    namespace = {}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), "exec"), namespace)
    return namespace[node.name]


def test_ambiguous_element_faces_stop_at_a_global_budget():
    calls = []
    def mesh(*point, **kwargs):
        calls.append(point)
        return SimpleNamespace(nr=len(calls) % 2)
    with pytest.raises(RuntimeError, match="point-location budget"):
        list(quadrature()(mesh, np.zeros(3), np.ones(3), [0.5], [1.0], max_evaluations=30))
    assert len(calls) == 30


def test_single_element_keeps_the_quadrature_weights():
    def mesh(*point, **kwargs):
        return SimpleNamespace(nr=1)
    values = list(quadrature()(mesh, np.zeros(3), np.ones(3), [.25, .75], [.5, .5]))
    assert [w for _, w in values] == [.5, .5]
