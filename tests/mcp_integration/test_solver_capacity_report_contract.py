"""The MCP report checklist must name real application result fields."""
import ast
from pathlib import Path


def test_ih_report_fields_exist_in_application_result():
    root = Path(__file__).resolve().parents[2]
    source = root / "src/radia/panels/calc_fem_kelvin.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    result_keys = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
            if any(isinstance(t, ast.Name) and t.id == "result" for t in node.targets):
                result_keys.update(k.value for k in node.value.keys
                                   if isinstance(k, ast.Constant) and isinstance(k.value, str))
    from radia_mcp.matrix_solvers.direct_solvers_knowledge import SOLVER_CAPACITY
    fields = {"linear_solver_requested", "linear_solver", "bddc_ams_coarse_cycles",
              "ndof", "ne", "linear_krylov_iterations", "linear_true_relative_residual",
              "linear_true_residual_limit", "t_solve_s", "t_total_s"}
    assert fields <= result_keys
    assert all(field in SOLVER_CAPACITY for field in fields)
