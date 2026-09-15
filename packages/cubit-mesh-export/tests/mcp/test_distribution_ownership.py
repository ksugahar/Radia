"""Cubit-only installation owns its API, MCP entry point and runtime dependency."""
from importlib import metadata
import ast
from pathlib import Path

from packaging.requirements import Requirement


def test_exporter_owns_mcp_without_radia_dependency():
    distribution = metadata.distribution("cubit-mesh-export")
    requires = {Requirement(item).name for item in distribution.requires or []}
    assert "mcp" in requires
    assert not {"radia", "radia-mcp", "cae-mcp-core"}.intersection(requires)
    entry = next(e for e in distribution.entry_points if e.name == "mcp-server-cubit")
    assert entry.value == "cubit_mesh_export.mcp.server:main"
    from cubit_mesh_export.mcp.api_reference import get_api_reference
    assert get_api_reference("all")


def test_cubit_runtime_has_no_other_product_imports():
    import cubit_mesh_export.mcp
    root = Path(cubit_mesh_export.mcp.__file__).parent
    violations = []
    for path in root.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
            modules = ([node.module or ""] if isinstance(node, ast.ImportFrom)
                       else [alias.name for alias in node.names] if isinstance(node, ast.Import)
                       else [])
            if any(name.split('.')[0] in {"radia", "radia_mcp", "cae_mcp_core"} for name in modules):
                violations.append(f"{path.relative_to(root)}:{node.lineno}")
    assert not violations, violations


def test_cubit_support_contains_only_cubit_example_providers():
    from cubit_mesh_export.mcp._support import examples
    assert set(examples.FAMILIES) == {'cubit'}
    assert set(examples.REFRESH_FUNCS) == set(examples.FAMILIES['cubit'])
    assert not hasattr(examples, 'refresh_build123d_examples')
