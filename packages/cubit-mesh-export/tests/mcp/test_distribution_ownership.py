"""Cubit-only installation owns its API, MCP entry point and runtime dependency."""
from importlib import metadata

from packaging.requirements import Requirement


def test_exporter_owns_mcp_without_radia_dependency():
    distribution = metadata.distribution("cubit-mesh-export")
    requires = {Requirement(item).name for item in distribution.requires or []}
    assert "cae-mcp-core" in requires
    assert not {"radia", "radia-mcp"}.intersection(requires)
    entry = next(e for e in distribution.entry_points if e.name == "mcp-server-cubit")
    assert entry.value == "cubit_mesh_export.mcp.server:main"
    from cubit_mesh_export.mcp.api_reference import get_api_reference
    assert get_api_reference("all")
