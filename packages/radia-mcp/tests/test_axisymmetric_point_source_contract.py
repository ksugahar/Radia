from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOLVE = (
    ROOT
    / "packages"
    / "radia-mcp"
    / "src"
    / "radia_mcp"
    / "radia_ngsolve"
    / "solve.py"
)


def _function(tree: ast.Module, name: str) -> ast.FunctionDef:
    return next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    )


def test_nonlinear_axisymmetric_solver_preserves_legacy_positional_prefix():
    tree = ast.parse(SOLVE.read_text(encoding="utf-8"), filename=str(SOLVE))
    function = _function(tree, "solve_axi_magnetostatic_nonlinear")
    names = [argument.arg for argument in function.args.args]

    assert names[:9] == [
        "mesh",
        "nu_of_B",
        "Jr",
        "order",
        "dirichlet",
        "relax",
        "max_iter",
        "tol",
        "min_iter",
    ]
    assert names[9:] == [
        "magnets",
        "ring_currents",
        "point_potentials",
        "dof_constraints",
        "mixed_boundaries",
    ]


def test_axisymmetric_source_rows_are_normalized_without_array_truth_tests():
    source = SOLVE.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(SOLVE))

    assert "ring_currents or []" not in source
    assert "point_potentials or []" not in source
    for name in (
        "axisymmetric_ring_current_load_contract",
        "axisymmetric_point_potential_constraint_contract",
        "solve_axi_magnetostatic",
        "solve_axi_magnetostatic_nonlinear",
    ):
        segment = ast.get_source_segment(source, _function(tree, name))
        assert segment is not None
        assert "_axisymmetric_triplet_rows" in segment
