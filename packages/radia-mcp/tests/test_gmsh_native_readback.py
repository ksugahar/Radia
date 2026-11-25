"""Native Gmsh readback of post-view MSH files: what the viewer actually loads.

Guards: `radia_mcp.gmsh.readback` reports node, element and group counts from
native Gmsh 5 and refuses files that Gmsh cannot open or loads with warnings or
unconnected nodes; `convert_post_view` publishes a VOL post view only after
that readback matches the selected groups and element count.
Assumes: the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE); the end-to-end
cases also need NGSolve/Netgen OCC and Radia's production MSH exporter. Without
them the tests skip visibly with the missing piece named.
"""
import hashlib
import json
import os
import tempfile
from pathlib import Path

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


@pytest.fixture
def checked_vol(tmp_path, monkeypatch):
    occ = pytest.importorskip("netgen.occ", reason="needs Netgen OCC to build the VOL fixture")
    ng = pytest.importorskip("ngsolve", reason="needs NGSolve to load the VOL")
    exporter = pytest.importorskip("radia.gmsh_post_export", reason="needs Radia's production MSH exporter")
    if not hasattr(exporter, "PayloadV42"):
        pytest.skip("installed radia lacks the production MSH 4.2 exporter")
    source_dir = tmp_path / "in"; source_dir.mkdir()
    scratch = tmp_path / "out"; scratch.mkdir()
    iron = occ.Box((0, 0, 0), (0.01, 0.01, 0.01)); iron.mat("iron"); iron.faces.name = "iron_wall"
    air = occ.Box((0.01, 0, 0), (0.02, 0.01, 0.01)); air.mat("air"); air.faces.name = "air_wall"
    source = source_dir / "two_region.vol"
    occ.OCCGeometry(occ.Glue([iron, air])).GenerateMesh(maxh=0.005).Save(str(source))
    mesh = ng.Mesh(str(source))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    # The handoff is caller-asserted (documented contract): counts come from NGSolve itself.
    report = {"schema": "cubit-mesh-export.vol-check.v1", "passed": True,
              "mesh": {"dimension": mesh.dim, "curve_order": 1, "n_elements": mesh.ne, "n_points": mesh.nv},
              "labels": {"schema": "radia.vol-label-contract.v1", "passed": True,
                         "materials": sorted(set(mesh.GetMaterials())),
                         "boundaries": sorted(set(mesh.GetBoundaries()))}}
    handoff = tmp_path / "handoff.json"
    handoff.write_text(json.dumps({"schema": "radia.checked-vol-handoff.v1", "source_sha256": digest,
                                   "report": report}), encoding="utf-8")
    monkeypatch.setenv("RADIA_MESH_VIEWER_SCRATCH", str(scratch))
    return source, scratch, handoff, digest, mesh, ng


def _convert(checked_vol, scope):
    from radia_mcp.mesh_viewer import post_view
    source, scratch, handoff, digest, _, _ = checked_vol
    operation = Path(tempfile.mkdtemp(prefix=".radia-post-receipt-", dir=scratch))
    return post_view.convert_post_view(source, str(scratch / f"view_{scope}.msh"), str(handoff), digest, scope,
                                       staging_dir=operation / "staging")


@pytest.mark.parametrize("scope", ["materials", "boundaries"])
def test_checked_vol_post_view_loads_cleanly_in_native_gmsh(checked_vol, scope, monkeypatch):
    from radia_mcp.gmsh import readback as module
    _, _, _, _, mesh, ng = checked_vol
    boundary = scope == "boundaries"
    selected = sum(1 for _ in mesh.Elements(ng.BND if boundary else ng.VOL))
    calls = []
    original = module.require_native_readback
    monkeypatch.setattr(module, "require_native_readback",
                        lambda path, **kwargs: calls.append(kwargs) or original(path, **kwargs))
    try:
        result = _convert(checked_vol, scope)
    except RuntimeError as error:
        if boundary and "'unreferenced_nodes': 1" in str(error):
            pytest.skip("installed radia exporter predates boundary-node compaction; "
                        f"the conversion correctly refused its output: {error}")
        raise
    readback = result["native_readback"]
    names = mesh.GetBoundaries() if boundary else mesh.GetMaterials()
    # The conversion binds the native load to its own selection, not just to "it opened".
    assert calls == [{"expected_groups": {(mesh.dim - int(boundary), n) for n in names},
                      "expected_elements": {mesh.dim - int(boundary): selected}}], calls
    assert readback["elements"] == {str(mesh.dim - int(boundary)): selected}, readback
    assert {tuple(g) for g in readback["physical_groups"]} == {(mesh.dim - int(boundary), n) for n in names}
    assert readback["unreferenced_nodes"] == 0
    if not boundary:
        assert readback["nodes"] == mesh.nv, readback
    # The published bytes are the bytes that were read back.
    assert native_readback(result["output"])["nodes"] == readback["nodes"]
