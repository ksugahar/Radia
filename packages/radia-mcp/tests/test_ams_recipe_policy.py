"""Check the executable AMS recipe without importing the MCP server stack."""
import ast
from pathlib import Path
import runpy


def test_ams_recipe_keeps_shift_out_of_physical_system():
    path = Path(__file__).parents[1] / "src/radia_mcp/matrix_solvers/preconditioners_knowledge.py"
    knowledge = runpy.run_path(str(path))
    recipe = next(value for value in knowledge.values()
                  if isinstance(value, str) and "## Code recipe" in value)
    code = recipe.split("## Code recipe", 1)[1].split("```python", 1)[1].split("```", 1)[0]
    tree = ast.parse(code)
    additions = [node for node in ast.walk(tree) if isinstance(node, ast.AugAssign)
                 and isinstance(node.target, ast.Name)]
    system = [node for node in additions if node.target.id == "a"]
    surrogate = [node for node in additions if node.target.id == "ar"]
    assert system and surrogate
    assert all("eps" not in {n.id for n in ast.walk(node.value)
                              if isinstance(n, ast.Name)} for node in system)
    assert any("eps" in {n.id for n in ast.walk(node.value)
                         if isinstance(n, ast.Name)} for node in surrogate)
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "ComplexCompactAMSPreconditioner")
    assert not call.args
    assert {kw.arg for kw in call.keywords} >= {
        "a_real_mat", "grad_mat", "freedofs", "coord_x", "coord_y", "coord_z"}
