"""Guards: workflow examples reference callable post-processing tools correctly.

Assumes: examples describe the public API; this does not execute a solver/render.
"""
import ast
import inspect
import re

from radia_mcp.gmsh import server
from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


def test_workflow_calls_match_public_signatures_and_inspect_before_render():
    body = get_gmsh_documentation("workflow")
    calls = []
    for code in re.findall(r"```python\n(.*?)```", body, re.DOTALL):
        tree = ast.parse(code)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id.startswith("gmsh_"):
                name = node.func.id
                tool = getattr(server, name, None)
                assert callable(tool), f"workflow references missing tool {name}"
                inspect.signature(tool).bind(*[object() for _ in node.args],
                                             **{arg.arg: object() for arg in node.keywords})
                calls.append(name)
    assert calls, "workflow lost its executable tool examples"
    assert calls.index("gmsh_inspect_msh") < calls.index("gmsh_write_post_launch_artifact"), calls
    assert calls.index("gmsh_validate_msh") < calls.index("gmsh_render"), calls
