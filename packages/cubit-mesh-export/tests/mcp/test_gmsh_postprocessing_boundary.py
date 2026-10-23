"""Guards: Gmsh inspection cannot claim solver-input or metadata authority.
Assumes: MSH 4.2 is visualization-only; checked .vol and contract own solver labels.
"""
import ast
from pathlib import Path
from cubit_mesh_export.mcp import server
from cubit_mesh_export.mcp.gmsh_v42 import gmsh_v42_mixed_order_series_gate


def test_gmsh_inventory_cannot_certify_a_solver_mesh():
    result = gmsh_v42_mixed_order_series_gate([])
    assert result["artifact_role"] == "post-processing-only"
    assert result["solver_mesh"] is False
    assert result["solver_metadata_authority"] == "checked .vol and its versioned contract"
    assert result["success"] is False


def test_removed_gmsh_to_vol_workflow_is_unavailable():
    result = server.generate_cubit_script("tet_gmsh")
    assert "Unknown workflow" in result, result
    assert "export netgen" in server.generate_cubit_script("tet_netgen")


def test_exporter_has_no_gmsh_meshing_or_solver_conversion_calls():
    # Inspect shipped production code: stored fixtures never run a Gmsh mesher.
    root = Path(server.__file__).resolve().parent.parent
    violations = []
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                name = ast.unparse(node.func)
                if name == "gmsh.model.mesh.generate" or name.endswith("ReadGmshMesh"):
                    violations.append((str(path.relative_to(root)), node.lineno, name))
    assert not violations, violations
