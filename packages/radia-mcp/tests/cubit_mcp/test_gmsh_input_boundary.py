"""Guards: Gmsh MSH is never meshed, converted or staged as solver input.
Assumes: checked .vol and its contract own solver labels.
"""
import ast
from pathlib import Path
from radia_mcp.cubit import server


def test_removed_gmsh_to_vol_workflow_is_unavailable():
    result = server.generate_cubit_script("tet_gmsh")
    assert "Unknown workflow" in result, result
    script = server.generate_cubit_script("tet_netgen")
    assert "cubit_export_vol" in script and "cubit.cmd('export vol" not in script, script


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


def test_stage_rejects_msh_before_starting_cubit(tmp_path, monkeypatch):
    import json
    def unexpected_session():
        raise AssertionError("invalid solver input must not start Cubit")
    monkeypatch.setattr(server, "_cubit_session_or_error", unexpected_session)
    mesh = tmp_path / "input.msh"
    mesh.write_text("$MeshFormat\n4.2 0 8\n$EndMeshFormat\n", encoding="ascii")
    for kwargs in [{"path": str(mesh)}, {"extra_commands": [f'import mesh "{mesh}"', 'export vol "bad.vol"']}]:
        result = json.loads(server.cubit_stage(**kwargs))
        assert result["status"] == "error", result
        assert result["stage"] == "input", result
        assert "not a supported input" in result["error"], result
