"""Executable-manual contracts without loading native solver libraries."""
import ast
from pathlib import Path
import re


SOURCE = Path(__file__).resolve().parents[1] / 'src/radia_mcp/matrix_solvers/preconditioners_knowledge.py'


def ams_text():
    tree = ast.parse(SOURCE.read_text(encoding='utf-8'))
    return next(ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'AMS_HIPTMAIR_XU'
                        for t in node.targets))


def test_ams_recipe_imports_its_ngsolve_operations():
    recipe = re.findall(r'```python\n(.*?)```', ams_text(), flags=re.S)[-1]
    tree = ast.parse(recipe)
    imports = {alias.name for node in ast.walk(tree)
               if isinstance(node, ast.ImportFrom) and node.module == 'ngsolve'
               for alias in node.names}
    assert {'curl', 'dx', 'CF', 'HCurl', 'TaskManager'} <= imports
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == 'ComplexCompactAMSPreconditioner')
    assert not call.args
    assert {'a_real_mat', 'grad_mat', 'freedofs', 'coord_x', 'coord_y', 'coord_z'} <= {
        keyword.arg for keyword in call.keywords}
    for region in (node for node in ast.walk(tree) if isinstance(node, ast.With)):
        assert call not in list(ast.walk(region))


def test_ams_knowledge_bounds_claims_and_does_not_link_private_memory():
    text = ams_text()
    assert 'stable up to p=10' not in text
    assert 'HYPRE supports high-order' in text
    assert 'not just once' in text
    assert 'not independently certified release evidence' in text
    assert '`validation_test/`' in text
    assert 'memory/' not in text
