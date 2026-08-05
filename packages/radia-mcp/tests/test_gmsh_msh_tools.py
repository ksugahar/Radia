"""Tests for radia_mcp.gmsh.msh_inspect (MSH v4.1 inspect/validate tools).

Structural tests are pure Python (run in the minimal-dep matrix).  The
Jacobian tests need the gmsh Python package and skip when it is absent;
the check itself runs gmsh in a subprocess, never in this process.
"""

import importlib.util
from pathlib import Path

import pytest

from radia_mcp.gmsh.msh_inspect import (
    audit_msh_directory,
    diff_msh,
    field_stats,
    inspect_msh,
    main as msh_inspect_main,
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


def test_validate_geo_deep_counts_views_and_flags_out_of_range(tmp_path):
    _write(tmp_path, _BASE_MSH)  # contributes ONE view ("T", 2 steps)
    geo = tmp_path / "case.geo"
    geo.write_text('Merge "case.msh";\n'
                   "View[0].Visible = 1;\n"
                   "View[3].Visible = 1;\n", encoding="utf-8")

    result = validate_geo(geo)
    assert result["merged_views_total"] == 1
    assert result["merge_targets"][0]["views"] == 1
    assert result["checks"]["view_indices_in_range"] is False
    assert result["ok"] is False
    assert any("View[3]" in e for e in result["errors"])
    assert result["sidecars"]["geo_opt"] is False

    # In-range references pass, and the sidecar report notices .geo.opt
    geo.write_text('Merge "case.msh";\nView[0].Visible = 1;\n',
                   encoding="utf-8")
    (tmp_path / "case.geo.opt").write_text("// sidecar\n", encoding="utf-8")
    result = validate_geo(geo)
    assert result["ok"] is True
    assert result["checks"]["view_indices_in_range"] is True
    assert result["sidecars"]["geo_opt"] is True


# ======================================================================
# diff_msh
# ======================================================================

def test_diff_msh_identical(tmp_path):
    a = _write(tmp_path, _BASE_MSH, "a.msh")
    b = _write(tmp_path, _BASE_MSH, "b.msh")
    result = diff_msh(a, b)

    assert result["ok"] is True
    assert result["identical_structure"] is True
    assert result["fields_match"] is True
    assert result["differences"] == []


def test_diff_msh_detects_structure_change(tmp_path):
    a = _write(tmp_path, _BASE_MSH, "a.msh")
    b = _write(tmp_path, _TRI6_MSH, "b.msh")
    result = diff_msh(a, b)

    assert result["identical_structure"] is False
    assert any("node count" in d for d in result["differences"])
    assert any("tri6" in d for d in result["differences"])
    assert result["views"]["only_a"] == ["T"]


def test_diff_msh_detects_field_drift(tmp_path):
    a = _write(tmp_path, _BASE_MSH, "a.msh")
    drifted = _BASE_MSH.replace("4 4.5\n", "4 4.6\n")
    b = _write(tmp_path, drifted, "b.msh")
    result = diff_msh(a, b)

    assert result["identical_structure"] is True
    assert result["fields_match"] is False
    common = result["views"]["common"][0]
    assert common["max_rel_delta"] == pytest.approx(0.1 / 4.6, rel=1e-6)
    assert any("drift" in d for d in result["differences"])


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
# ElementData / ElementNodeData coverage
# ======================================================================

_ELEMENT_DATA = (
    "$ElementData\n1\n\"quality\"\n1\n0.0\n3\n0\n1\n2\n"
    "1 0.9\n2 0.8\n$EndElementData\n"
    "$ElementNodeData\n1\n\"nodal_flux\"\n1\n0.0\n3\n0\n1\n1\n"
    "2 4 1.0 2.0 3.0 4.0\n$EndElementNodeData\n")


def test_element_data_sections_parse_validate_and_stats(tmp_path):
    msh = _write(tmp_path, _BASE_MSH + _ELEMENT_DATA)

    result = validate_msh(msh)
    assert result["ok"] is True, result["errors"]

    info = inspect_msh(msh)
    sections = {(v["section"], v["name"]) for v in info["views"]}
    assert ("ElementData", "quality") in sections
    assert ("ElementNodeData", "nodal_flux") in sections

    stats = field_stats(msh, view_name="quality")
    step = stats["views"][0]["per_step"][0]
    assert step["min"] == 0.8
    assert step["max"] == 0.9


def test_element_data_undefined_element_tag_fails(tmp_path):
    broken = _BASE_MSH + _ELEMENT_DATA.replace("1 0.9\n", "99 0.9\n")
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_tags_exist"] is False
    assert any("99" in e for e in result["errors"])


def test_elementnodedata_width_mismatch_fails(tmp_path):
    broken = _BASE_MSH + _ELEMENT_DATA.replace(
        "2 4 1.0 2.0 3.0 4.0\n", "2 4 1.0 2.0 3.0\n")
    result = validate_msh(_write(tmp_path, broken))

    assert result["ok"] is False
    assert result["checks"]["data_row_width_matches"] is False


# ======================================================================
# Directory audit + CLI
# ======================================================================

def test_audit_msh_directory_reports_issues(tmp_path):
    _write(tmp_path, _BASE_MSH, "good.msh")
    _write(tmp_path, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n", "legacy.msh")
    sub = tmp_path / "sub"
    sub.mkdir()
    _write(sub, _TRI6_MSH, "curved.msh")

    audit = audit_msh_directory(tmp_path)
    assert audit["ok"] is True
    assert audit["files_scanned"] == 3
    assert audit["clean"] is False
    assert audit["by_status"] == {"ok": 2, "needs_attention": 1}
    issue = audit["issues"][0]
    assert issue["path"] == "legacy.msh"
    assert "format_is_v41" in issue["failed_checks"]
    highorder = [f for f in audit["files"] if f["path"].endswith("curved.msh")]
    assert highorder[0]["high_order"] is True


def test_cli_exit_codes_and_modes(tmp_path, capsys):
    good = _write(tmp_path, _BASE_MSH, "good.msh")
    bad = _write(tmp_path, _BASE_MSH.replace("1\n4\n1 1.0", "1\n5\n1 1.0", 1),
                 "bad.msh")

    assert msh_inspect_main([str(good), "--validate"]) == 0
    assert msh_inspect_main([str(bad), "--validate"]) == 1
    assert msh_inspect_main([str(tmp_path)]) == 1        # audit: bad inside
    assert msh_inspect_main([str(good), "--diff", str(good)]) == 0
    assert msh_inspect_main([str(good), "--diff", str(bad)]) == 0  # same data
    out = capsys.readouterr().out
    assert "[ok]" in out and "[needs_attention]" in out

    geo = tmp_path / "case.geo"
    geo.write_text('Merge "good.msh";\nView[0].Visible = 1;\n',
                   encoding="utf-8")
    assert msh_inspect_main([str(geo)]) == 0
    assert msh_inspect_main([str(good), "--stats", "--json"]) == 0


# ======================================================================
# Mesh quality gate (needs gmsh)
# ======================================================================

_TET10_TEMPLATE = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 10 1 10
3 1 0 10
1
2
3
4
5
6
7
8
9
10
0 0 0
1 0 0
0 1 0
0 0 1
{e01}
0.5 0.5 0
0 0.5 0
0 0 0.5
0 0.5 0.5
0.5 0 0.5
$EndNodes
$Elements
1 1 1 1
3 1 11 1
1 1 2 3 4 5 6 7 8 9 10
$EndElements
"""


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_mesh_quality_affine_tet10_is_perfect(tmp_path):
    from radia_mcp.gmsh.msh_inspect import mesh_quality

    msh = _write(tmp_path, _TET10_TEMPLATE.format(e01="0.5 0 0"))
    q = mesh_quality(msh, threshold=0.5)
    assert q["ran"] is True, q.get("error")
    assert q["ok"] is True
    bt = q["by_type"][0]
    assert bt["min_scaled"] == pytest.approx(1.0)
    assert bt["negative"] == 0
    assert bt["below_threshold"] == 0


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_mesh_quality_flags_degrading_curved_element(tmp_path):
    from radia_mcp.gmsh.msh_inspect import mesh_quality

    # Mid-edge node pushed sideways: NOT inverted (sign gate passes)
    # but the scaled Jacobian collapses to ~0.49.
    msh = _write(tmp_path, _TET10_TEMPLATE.format(e01="0.5 0.12 0.05"))
    v = validate_msh(msh, check_jacobians=True)
    assert v["ok"] is True, "sign gate must PASS for this element"

    q = mesh_quality(msh, threshold=0.5)
    assert q["ok"] is False
    bt = q["by_type"][0]
    assert bt["negative"] == 0
    assert bt["below_threshold"] == 1
    assert bt["min_scaled"] == pytest.approx(0.4895, abs=0.01)
    assert bt["worst"][0]["tag"] == 1


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_mesh_quality_counts_inverted_separately(tmp_path):
    from radia_mcp.gmsh.msh_inspect import mesh_quality

    msh = _write(tmp_path, _TET10_TEMPLATE.format(e01="0.5 0.3 0.15"))
    q = mesh_quality(msh, threshold=0.5)
    assert q["ok"] is False
    bt = q["by_type"][0]
    assert bt["negative"] == 1
    assert bt["below_threshold"] == 0  # inverted, not double-counted


# ======================================================================
# Dynamic option probing (needs gmsh)
# ======================================================================

@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_probe_options_flags_missing_and_reports_kind():
    from radia_mcp.gmsh.msh_inspect import probe_options

    result = probe_options(["Mesh.NumSubEdges", "View[0].Visible",
                            "General.Color.Background", "Mesh.Volumes"])
    assert result["ran"] is True
    assert result["ok"] is False
    assert result["missing"] == ["Mesh.Volumes"]
    opts = result["options"]
    assert opts["Mesh.NumSubEdges"]["kind"] == "number"
    assert opts["View[0].Visible"]["exists"] is True
    assert opts["View[0].Visible"]["normalized"] == "View.Visible"
    assert opts["General.Color.Background"]["kind"] == "color"


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_validate_geo_check_options_catches_typo(tmp_path):
    _write(tmp_path, _BASE_MSH)
    geo = tmp_path / "case.geo"
    geo.write_text('Merge "case.msh";\n'
                   "Mesh.NumSubEdgs = 4;\n"     # typo
                   "View[0].IntervalsType = 2;\n", encoding="utf-8")

    result = validate_geo(geo, check_options=True)
    assert result["checks"]["option_names_exist"] is False
    assert result["ok"] is False
    assert any("NumSubEdgs" in e for e in result["errors"])

    geo.write_text('Merge "case.msh";\nMesh.NumSubEdges = 4;\n'
                   "View[0].IntervalsType = 2;\n", encoding="utf-8")
    result = validate_geo(geo, check_options=True)
    assert result["checks"]["option_names_exist"] is True
    assert result["ok"] is True


# ======================================================================
# verify_artifact (one-call gate runner)
# ======================================================================

@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_verify_artifact_geo_runs_all_gates(tmp_path):
    from radia_mcp.gmsh.verify import verify_artifact

    _write(tmp_path, _BASE_MSH)
    geo = tmp_path / "case.geo"
    geo.write_text('Merge "case.msh";\nView[0].Visible = 1;\n',
                   encoding="utf-8")

    result = verify_artifact(geo)
    assert result["ok"] is True
    assert set(result["passed"]) == {"geo:case.geo", "msh:case.msh"}
    assert result["failed"] == []
    assert result["jacobians_checked"] is True


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_verify_artifact_msh_reports_failed_gate(tmp_path):
    from radia_mcp.gmsh.verify import verify_artifact

    inverted = _BASE_MSH.replace("2 1 2 3 4\n", "2 1 2 4 3\n")
    msh = _write(tmp_path, inverted)

    result = verify_artifact(msh)
    assert result["ok"] is False
    assert result["failed"] == ["msh:case.msh"]
    gate = result["gates"][0]
    assert "jacobians_positive" in gate["failed_checks"]


def test_verify_artifact_rejects_unknown_type(tmp_path):
    from radia_mcp.gmsh.verify import verify_artifact

    other = tmp_path / "case.step"
    other.write_text("dummy", encoding="utf-8")
    result = verify_artifact(other)
    assert result["ok"] is False
    assert "unsupported artifact type" in result["error"]


# ======================================================================
# Lint fixtures lock (selftest companions)
# ======================================================================

def test_lint_fixture_bad_script_trips_every_rule_class():
    from radia_mcp.gmsh.server import _lint_file

    fixture = (Path(__file__).parent / "mcp_server" / "fixtures"
               / "bad_gmsh_script.py")
    findings = _lint_file(str(fixture))
    rules = {f["rule"] for f in findings}
    assert {"pip-gmsh-import", "gmsh-mesh-generation", "meshio-removed",
            "gmsh-builder-removed", "invalid-gmsh-option",
            "readgmsh-deprecated"} <= rules


def test_lint_fixture_clean_script_has_no_gmsh_findings():
    from radia_mcp.gmsh.server import _lint_file

    fixture = (Path(__file__).parent / "mcp_server" / "fixtures"
               / "clean_gmsh_script.py")
    findings = _lint_file(str(fixture))
    gmsh_findings = [f for f in findings
                     if f["rule"].startswith(("gmsh-", "pip-gmsh", "meshio-",
                                              "msh-", "numsubedges",
                                              "readgmsh", "invalid-gmsh"))]
    assert gmsh_findings == []


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
    blank = result.get("blank_check", {})
    if blank.get("ran"):
        assert blank["looks_blank"] is False, result


@pytest.mark.skipif(not _GMSH_AVAILABLE, reason="gmsh package not installed")
def test_render_png_with_structured_cut_plane(tmp_path):
    from radia_mcp.gmsh.render import render_png

    msh = _write(tmp_path, _BASE_MSH)
    out = tmp_path / "cut.png"
    result = render_png(msh, out, width=400, height=300,
                        cut_plane={"enabled": True, "normal": [0, -1, 0],
                                   "offset": 0.5})
    _skip_if_no_graphics(result)
    assert result["ok"] is True, result
    assert out.is_file()


@pytest.mark.skipif(not _GMSH_AVAILABLE or
                    importlib.util.find_spec("PIL") is None,
                    reason="gmsh or Pillow not installed")
def test_render_png_flags_blank_image(tmp_path):
    from radia_mcp.gmsh.render import render_png

    msh = _write(tmp_path, _BASE_MSH)
    out = tmp_path / "blank.png"
    # hide the view AND all mesh display: nothing gets drawn.
    # (gmsh.open() on a view-bearing .msh flips Mesh.SurfaceFaces to 1
    # by itself -- observed on gmsh 4.15.2 -- so it must be overridden.)
    result = render_png(msh, out, width=500, height=400,
                        auto_mesh_display=False,
                        options={"View[0].Visible": 0,
                                 "Mesh.SurfaceFaces": 0,
                                 "Mesh.VolumeEdges": 0,
                                 "Mesh.SurfaceEdges": 0,
                                 "Mesh.Lines": 0})
    _skip_if_no_graphics(result)

    assert result["ok"] is True
    blank = result["blank_check"]
    assert blank["ran"] is True
    assert blank["looks_blank"] is True
    assert "blank" in result.get("note", "")


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
