"""Native Gmsh readback of post-view MSH files: what the viewer actually loads.

Guards: `radia_mcp.gmsh.readback` reports node, element and group counts from
native Gmsh 5 for the MSH bytes alone (never an adjacent companion) and
refuses files that Gmsh cannot open or loads with warnings or unconnected
nodes; the viewer preflight runs it. Real VOL conversions are validated in
validation_test/radia_mcp/test_gmsh_post_view_native_readback.py.
Assumes: the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE); without it the
tests skip visibly.
"""
import json
import os

import pytest

from radia_mcp.gmsh.readback import native_readback, require_native_readback

pytestmark = pytest.mark.skipif(not os.environ.get("RADIA_GMSH_PROVENANCE"),
                                reason="needs the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE)")

HEADER = "$MeshFormat\n4.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n2 1 \"skin\"\n$EndPhysicalNames\n"
# Entity and element layout as native Gmsh 5 itself writes one discrete triangle.
ENTITIES = "$Entities\n0 0 1 0\n1 0 0 0 1 1 0 1 1 0 0\n$EndEntities\n"


def _triangle(extra_node: bool) -> str:
    count = 4 if extra_node else 3
    coords = ["0 0 0", "1 0 0", "0 1 0", "0.3 0.3 0"][:count]
    nodes = (f"$Nodes\n1 {count} 1 {count}\n2 1 0 {count}\n" + "".join(f"{i}\n" for i in range(1, count + 1))
             + "".join(c + "\n" for c in coords) + "$EndNodes\n")
    elements = "$Elements\n1 1 1 1\n2 1 2 1\n1 1 2 3\n$EndElements\n"
    return HEADER + ENTITIES + nodes + elements


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="ascii", newline="\n")
    return path


def test_clean_file_reports_what_native_gmsh_loaded(tmp_path):
    result = native_readback(_write(tmp_path, "clean.msh", _triangle(False)))
    assert result["ok"] is True, result
    assert (result["nodes"], result["elements"], result["physical_groups"]) == (3, {"2": 1}, [[2, "skin"]]), result
    assert result["unreferenced_nodes"] == 0 and result["warnings"] == [] and result["errors"] == []


def test_unconnected_node_is_refused(tmp_path):
    path = _write(tmp_path, "orphan.msh", _triangle(True))
    result = native_readback(path)
    assert (result["ok"], result["nodes"], result["unreferenced_nodes"]) == (False, 4, 1), result
    with pytest.raises(RuntimeError, match="'unreferenced_nodes': 1"):
        require_native_readback(path)


def test_file_loaded_with_a_native_warning_is_refused(tmp_path):
    # Native Gmsh still loads it but warns that the node-tag header is wrong.
    result = native_readback(_write(tmp_path, "header.msh", _triangle(False).replace("1 3 1 3\n", "1 3 1 4\n")))
    assert result["ok"] is False and result["nodes"] == 3, result
    assert any("node tags reported in section header" in w for w in result["warnings"]), result


def test_file_native_gmsh_cannot_open_is_refused(tmp_path):
    truncated = _triangle(False).split("$Elements")[0].replace("0 1 0\n$EndNodes", "0 1\n$EndNodes")
    result = native_readback(_write(tmp_path, "broken.msh", truncated))
    assert result["ok"] is False and "native gmsh.open failed" in result.get("error", ""), result


def test_expectations_must_match_the_native_load(tmp_path):
    path = _write(tmp_path, "clean.msh", _triangle(False))
    assert require_native_readback(path, expected_groups={(2, "skin")}, expected_elements={2: 1})["nodes"] == 3
    with pytest.raises(RuntimeError, match="physical groups"):
        require_native_readback(path, expected_groups={(2, "wall")})
    with pytest.raises(RuntimeError, match="element counts"):
        require_native_readback(path, expected_elements={2: 2})


def test_viewer_preflight_runs_the_native_readback(tmp_path, monkeypatch, capsys):
    from radia_mcp.mesh_viewer import launcher
    monkeypatch.setattr(launcher, "require_interactive_desktop", lambda: {})
    assert launcher.main([str(_write(tmp_path, "clean.msh", _triangle(False))), "--check"]) == 0
    identity = json.loads(capsys.readouterr().out)
    assert identity["native_readback"]["nodes"] == 3, identity
    assert launcher.main([str(_write(tmp_path, "orphan.msh", _triangle(True))), "--check"]) == 1
    assert "'unreferenced_nodes': 1" in capsys.readouterr().err


def test_adjacent_companion_is_not_read(tmp_path):
    # Native Gmsh reads FILE.msh.opt beside the file it opens; this one is invalid
    # and would log an error. The readback must see the MSH bytes only.
    path = _write(tmp_path, "view.msh", _triangle(False))
    (tmp_path / "view.msh.opt").write_text("this is not a Gmsh option script;\n", encoding="ascii")
    result = native_readback(path)
    assert result["ok"] is True and result["errors"] == [], result
    assert result["path"] == str(path.resolve()) and len(result["sha256"]) == 64
