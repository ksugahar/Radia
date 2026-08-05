"""Tests for radia_mcp.gmsh.msh_inspect (MSH v4.1 inspect/validate tools).

Structural tests are pure Python (run in the minimal-dep matrix).  The
Jacobian tests need the gmsh Python package and skip when it is absent;
the check itself runs gmsh in a subprocess, never in this process.
"""

import importlib.util

import pytest

from radia_mcp.gmsh.msh_inspect import (
    inspect_msh,
    validate_geo,
    validate_msh,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

# Unit tet (volume 1/6) + one boundary tri + a 2-step scalar NodeData view.
_BASE_MSH = """$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 2 "surface"
3 1 "block"
$EndPhysicalNames
$Entities
0 0 1 1
1 0 0 0 1 1 0 1 2 0
1 0 0 0 1 1 1 1 1 1 1
$EndEntities
$Nodes
2 4 1 4
2 1 0 3
1
2
3
0 0 0
1 0 0
0 1 0
3 1 0 1
4
0 0 1
$EndNodes
$Elements
2 2 1 2
2 1 2 1
1 1 2 3
3 1 4 1
2 1 2 3 4
$EndElements
$NodeData
1
"T"
1
0.0
3
0
1
4
1 1.0
2 2.0
3 3.0
4 4.0
$EndNodeData
$NodeData
1
"T"
1
0.5
3
1
1
4
1 1.5
2 2.5
3 3.5
4 4.5
$EndNodeData
"""

_TRI6_MSH = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 6 1 6
2 1 0 6
1
2
3
4
5
6
0 0 0
1 0 0
0 1 0
0.5 0 0
0.5 0.5 0
0 0.5 0
$EndNodes
$Elements
1 1 1 1
2 1 9 1
1 1 2 3 4 5 6
$EndElements
"""


def _write(tmp_path, text, name="case.msh"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# ======================================================================
# inspect_msh
# ======================================================================

def test_inspect_reports_structure(tmp_path):
    info = inspect_msh(_write(tmp_path, _BASE_MSH))

    assert info["ok"] is True
    assert info["version"] == "4.1"
    assert info["ascii"] is True
    assert info["parse_errors"] == []
    assert {p["name"] for p in info["physical_names"]} == {"surface", "block"}
    assert info["entities"] == {"points": 0, "curves": 0,
                                "surfaces": 1, "volumes": 1}
    assert info["nodes"]["count"] == 4
    assert info["elements"]["count"] == 2
    by_type = {t["name"]: t for t in info["elements"]["by_type"]}
    assert by_type["tri3"]["count"] == 1
    assert by_type["tet4"]["count"] == 1
    assert info["elements"]["max_dim"] == 3
    assert info["bbox"] == {"min": [0.0, 0.0, 0.0], "max": [1.0, 1.0, 1.0]}
    assert info["high_order"] is False

    assert len(info["views"]) == 1
    view = info["views"][0]
    assert view["name"] == "T"
    assert view["section"] == "NodeData"
    assert view["steps"] == 2
    assert view["components"] == 1
    assert view["time_range"] == [0.0, 0.5]
    assert view["entries_per_step"] == 4


def test_inspect_flags_high_order_display_hint(tmp_path):
    info = inspect_msh(_write(tmp_path, _TRI6_MSH))

    assert info["ok"] is True
    assert info["high_order"] is True
    assert any("NumSubEdges" in hint for hint in info["hints"])


def test_inspect_missing_file(tmp_path):
    info = inspect_msh(tmp_path / "nope.msh")
    assert info["ok"] is False
    assert "not found" in info["error"]


# ======================================================================
# validate_msh: structural checks
# ======================================================================

def test_validate_clean_mesh_passes(tmp_path):
    result = validate_msh(_write(tmp_path, _BASE_MSH))

    assert result["ok"] is True
    assert result["status"] == "ok"
    assert result["errors"] == []
    assert all(result["checks"].values())
    for name in ("format_is_v41", "is_ascii", "sections_balanced",
                 "node_tags_unique", "element_tags_unique",
                 "element_node_refs_exist", "data_declared_counts_match",
                 "data_tags_exist", "data_components_valid",
                 "view_components_consistent"):
        assert result["checks"][name] is True, name


def test_validate_detects_nodedata_count_mismatch(tmp_path):
    broken = _BASE_MSH.replace("1\n4\n1 1.0", "1\n5\n1 1.0", 1)
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_declared_counts_match"] is False
    assert any("declares 5 entries" in e for e in result["errors"])


def test_validate_detects_undefined_node_reference(tmp_path):
    broken = _BASE_MSH.replace("2 1 2 3 4\n", "2 1 2 3 99\n")
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["element_node_refs_exist"] is False
    assert any("99" in e for e in result["errors"])


def test_validate_detects_duplicate_element_tags(tmp_path):
    broken = _BASE_MSH.replace("2 1 2 3 4\n", "1 1 2 3 4\n")
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["element_tags_unique"] is False


def test_validate_detects_invalid_component_count(tmp_path):
    broken = _BASE_MSH.replace("0\n1\n4\n1 1.0", "0\n2\n4\n1 1.0", 1)
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_components_valid"] is False
    assert result["checks"]["view_components_consistent"] is False


def test_validate_rejects_v22(tmp_path):
    v22 = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    result = validate_msh(_write(tmp_path, v22))

    assert result["ok"] is False
    assert result["checks"]["format_is_v41"] is False


def test_validate_detects_unbalanced_sections(tmp_path):
    broken = _BASE_MSH.replace("$EndElements\n", "")
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["sections_balanced"] is False


# ======================================================================
# validate_geo
# ======================================================================

def test_validate_geo_clean(tmp_path):
    _write(tmp_path, _BASE_MSH)
    geo = tmp_path / "case.geo"
    geo.write_text('Merge "case.msh";\nMesh.NumSubEdges = 4;\n'
                   "View[0].IntervalsType = 2;\n", encoding="utf-8")

    result = validate_geo(geo)
    assert result["ok"] is True
    assert result["checks"]["merge_targets_exist"] is True
    assert result["checks"]["no_invalid_options"] is True
    assert result["merge_targets"][0]["exists"] is True
    assert result["numsubedges"] == 4
    assert result["max_view_index"] == 0


def test_validate_geo_detects_missing_merge_and_invalid_options(tmp_path):
    geo = tmp_path / "display.geo"
    geo.write_text('Merge "missing.msh";\n'
                   "Mesh.Volumes = 1;\n"
                   "General.GraphicsSizeX = 800;\n", encoding="utf-8")

    result = validate_geo(geo)
    assert result["ok"] is False
    assert result["checks"]["merge_targets_exist"] is False
    assert result["checks"]["no_invalid_options"] is False
    options = {o["option"] for o in result["invalid_options"]}
    assert options == {"Mesh.Volumes", "General.GraphicsSizeX"}
    assert any("Mesh.VolumeEdges" in o["fix"] for o in result["invalid_options"])


def test_validate_geo_missing_file(tmp_path):
    result = validate_geo(tmp_path / "nope.geo")
    assert result["ok"] is False


# ======================================================================
# Jacobian subprocess check (needs the gmsh Python package)
# ======================================================================

@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_validate_jacobians_positive_unit_tet(tmp_path):
    result = validate_msh(_write(tmp_path, _BASE_MSH), check_jacobians=True)

    jac = result["jacobian"]
    assert jac["ran"] is True, jac.get("error")
    assert jac["total_negative"] == 0
    assert result["checks"]["jacobians_positive"] is True
    assert result["ok"] is True
    assert jac["total_volume_dim3"] == pytest.approx(1.0 / 6.0, rel=1e-9)


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_validate_jacobians_detect_inverted_tet(tmp_path):
    inverted = _BASE_MSH.replace("2 1 2 3 4\n", "2 1 2 4 3\n")
    result = validate_msh(_write(tmp_path, inverted), check_jacobians=True)

    jac = result["jacobian"]
    assert jac["ran"] is True, jac.get("error")
    assert jac["total_negative"] > 0
    assert result["checks"]["jacobians_positive"] is False
    assert result["ok"] is False
