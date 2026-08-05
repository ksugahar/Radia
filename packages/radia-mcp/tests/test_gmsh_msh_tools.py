"""Tests for radia_mcp.gmsh.msh_inspect (MSH v4.1 inspect/validate tools).

Structural tests are pure Python (run in the minimal-dep matrix).  The
Jacobian tests need the gmsh Python package and skip when it is absent;
the check itself runs gmsh in a subprocess, never in this process.
"""

import importlib.util

import pytest

from radia_mcp.gmsh.msh_inspect import (
    field_stats,
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


# ======================================================================
# Field value statistics + NaN/Inf gating
# ======================================================================

def test_field_stats_scalar_two_steps(tmp_path):
    result = field_stats(_write(tmp_path, _BASE_MSH))

    assert result["ok"] is True
    assert len(result["views"]) == 1
    view = result["views"][0]
    assert view["name"] == "T"
    assert view["steps"] == 2
    step0 = view["per_step"][0]
    assert step0["metric"] == "value"
    assert step0["min"] == 1.0
    assert step0["max"] == 4.0
    assert step0["mean"] == pytest.approx(2.5)
    assert step0["rms"] == pytest.approx((30.0 / 4.0) ** 0.5)
    assert step0["nan"] == 0
    step1 = view["per_step"][1]
    assert step1["time"] == 0.5
    assert step1["max"] == 4.5
    assert view["overall"]["min"] == 1.0
    assert view["overall"]["max"] == 4.5


def test_field_stats_vector_magnitude(tmp_path):
    vec = _BASE_MSH.split("$NodeData")[0] + (
        "$NodeData\n1\n\"D\"\n1\n0.0\n3\n0\n3\n4\n"
        "1 3.0 4.0 0.0\n2 0.0 0.0 1.0\n3 1.0 0.0 0.0\n4 0.0 2.0 0.0\n"
        "$EndNodeData\n")
    result = field_stats(_write(tmp_path, vec))

    view = result["views"][0]
    assert view["components"] == 3
    step = view["per_step"][0]
    assert step["metric"] == "magnitude"
    assert step["max"] == 5.0        # |(3,4,0)|
    assert step["min"] == 1.0
    assert step["comp_min"] == 0.0
    assert step["comp_max"] == 4.0


def test_field_stats_view_name_filter_lists_available(tmp_path):
    msh = _write(tmp_path, _BASE_MSH)
    ok = field_stats(msh, view_name="T")
    assert ok["ok"] is True and len(ok["views"]) == 1

    missing = field_stats(msh, view_name="nope")
    assert missing["ok"] is False
    assert "'T'" in missing["error"] or "T" in missing["error"]


def test_validate_detects_nan_values(tmp_path):
    broken = _BASE_MSH.replace("3 3.0\n", "3 nan\n", 1)
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_values_finite"] is False
    assert any("NaN" in e for e in result["errors"])

    stats = field_stats(tmp_path / "case.msh")
    assert stats["views"][0]["overall"]["nan"] == 1


def test_validate_detects_wrong_row_width(tmp_path):
    broken = _BASE_MSH.replace("3 3.0\n", "3 3.0 9.9\n", 1)
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_row_width_matches"] is False


# ======================================================================
# Headless rendering (needs gmsh + an FLTK graphics context)
# ======================================================================

def _skip_if_no_graphics(result):
    if not result.get("ran"):
        pytest.skip(f"no gmsh graphics context: {result.get('error')}")


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_render_png_writes_image(tmp_path):
    from radia_mcp.gmsh.render import render_png

    msh = _write(tmp_path, _BASE_MSH)
    out = tmp_path / "case.png"
    result = render_png(msh, out, width=500, height=400)
    _skip_if_no_graphics(result)

    assert result["ok"] is True, result
    assert out.is_file()
    assert result["png_size"] is not None
    assert result["png_size"][1] == 400
    assert result["n_views"] == 1


@pytest.mark.skipif(not _GMSH_AVAILABLE or
                    importlib.util.find_spec("PIL") is None,
                    reason="gmsh or Pillow not installed")
def test_export_animation_two_step_gif(tmp_path):
    from radia_mcp.gmsh.render import export_animation

    msh = _write(tmp_path, _BASE_MSH)
    gif = tmp_path / "case.gif"
    result = export_animation(msh, gif, keep_frames=True,
                              width=400, height=300, delay_ms=100)
    _skip_if_no_graphics(result)

    assert result["ok"] is True, result
    assert result["num_steps"] == 2
    assert gif.is_file()
    frames = sorted((tmp_path / "case_frames").glob("frame_*.png"))
    assert len(frames) == 2


# ======================================================================
# Lint rule: invalid GMSH option names
# ======================================================================

def test_lint_flags_invalid_gmsh_option_names(tmp_path):
    from radia_mcp.gmsh.rules import check_invalid_gmsh_option_names

    script = tmp_path / "viewer.py"
    script.write_text(
        'import gmsh\n'
        'gmsh.option.setNumber("Mesh.Volumes", 1)\n'
        'gmsh.option.setNumber("General.GraphicsSizeX", 800)\n'
        'gmsh.option.setNumber("Mesh.VolumeEdges", 1)  # valid\n'
        'geo = "Mesh.SurfaceFaces = 1;"  # valid\n',
        encoding="utf-8")
    lines = script.read_text(encoding="utf-8").splitlines(keepends=True)
    findings = check_invalid_gmsh_option_names(str(script), lines)

    assert len(findings) == 2
    assert {f["line"] for f in findings} == {2, 3}
    assert all(f["rule"] == "invalid-gmsh-option" for f in findings)


def test_knowledge_records_mcp_tooling():
    from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation

    workflow = get_gmsh_documentation("workflow")
    assert "gmsh_inspect_msh" in workflow
    assert "gmsh_validate_msh" in workflow
    assert "gmsh_render" in workflow

    high_order = get_gmsh_documentation("high_order")
    assert "check_jacobians=True" in high_order
    assert "AdaptVisualizationGrid" in high_order

    pitfalls = get_gmsh_documentation("pitfalls")
    assert "1b." in pitfalls
    assert "gmsh_validate_msh" in pitfalls
