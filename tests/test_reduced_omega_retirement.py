"""Retired Kelvin API must fail before touching a mesh or assembling forms."""
import ast
from pathlib import Path

import pytest


def test_plain_reduced_omega_fails_without_numerical_dependencies():
    source = Path(__file__).resolve().parents[1] / "src/radia/kelvin_solver.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == "solve_magnetostatic_reduced_omega_kelvin")
    # Execute the actual entry point in isolation: no numerical globals may run.
    namespace = {}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"),
         namespace)
    with pytest.raises(NotImplementedError, match="Migrate to .*mixed_total_reduced"):
        namespace[function.name](None, None, None, None, mu_r_by_material={})
