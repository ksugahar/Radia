"""Fail-fast checks for the dependency-free MSH data reader."""

import csv
import re

import pytest
from radia_mcp.gmsh.msh_inspect import read_msh_data
from radia_mcp.gmsh.post_process import _PLOT_SCRIPT, export_view_csv
from radia_mcp.gmsh.render import render_montage

_VALID = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
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
"""


def _write(tmp_path, text):
    path = tmp_path / "case.msh"
    path.write_text(text, encoding="utf-8")
    return path


def test_reader_and_csv_accept_consistent_counts_without_gmsh(tmp_path):
    msh = _write(tmp_path, _VALID)
    data = read_msh_data(msh, include_elements=True)

    assert sorted(data["nodes"]) == [1, 2, 3, 4]
    assert data["elements"][1]["nodes"] == [1, 2, 3, 4]
    assert data["views"][0]["rows"][4] == [4.0]

    csv_path = tmp_path / "view.csv"
    result = export_view_csv(msh, csv_path)
    assert result["ok"] is True, result.get("error")
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    assert len(rows) == 5


@pytest.mark.parametrize(
    ("broken", "message"),
    [
        (_VALID.replace("1 4 1 4\n3 1 0 4", "1 5 1 4\n3 1 0 4"),
         "$Nodes declares 5 nodes"),
        (_VALID.replace("1 1 1 1\n3 1 4 1", "1 2 1 1\n3 1 4 1"),
         "$Elements declares 2 elements"),
        (_VALID.replace("0\n1\n4\n1 1.0", "0\n1\n5\n1 1.0"),
         "$NodeData"),
        (_VALID.replace("4 4.0", "4"),
         "has 0 values; expected 1"),
        (_VALID.replace("1\n2\n3\n4\n0 0 0", "1\n2\n3\n3\n0 0 0"),
         "duplicate node tag 3"),
        (_VALID.replace("3 3.0\n4 4.0", "3 3.0\n3 4.0"),
         "duplicate data tag 3"),
    ],
)
def test_reader_rejects_malformed_counts_widths_and_duplicates(
        tmp_path, broken, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        read_msh_data(_write(tmp_path, broken))


def test_reader_rejects_bad_element_node_data_width(tmp_path):
    element_node_data = """$ElementNodeData
1
"E"
1
0.0
3
0
1
1
1 4 10 20 30
$EndElementNodeData
"""
    msh = _write(tmp_path, _VALID + element_node_data)

    with pytest.raises(ValueError, match="has 3 values; expected 4"):
        read_msh_data(msh)


def test_reader_matches_element_node_data_count_to_connectivity(tmp_path):
    element_node_data = """$ElementNodeData
1
"E"
1
0.0
3
0
1
1
1 3 10 20 30
$EndElementNodeData
"""
    msh = _write(tmp_path, _VALID + element_node_data)

    with pytest.raises(ValueError, match="mesh element has 4"):
        read_msh_data(msh)


def test_reader_rejects_a_truncated_mesh_format_header(tmp_path):
    msh = _write(tmp_path, _VALID.replace("4.1 0 8", "4.1"))

    with pytest.raises(ValueError, match="malformed \\$MeshFormat header"):
        read_msh_data(msh)


def test_csv_returns_a_structured_error_for_a_malformed_reader_input(tmp_path):
    msh = _write(tmp_path, _VALID.replace("4 4.0", "4"))

    result = export_view_csv(msh, tmp_path / "must_not_exist.csv")

    assert result["ok"] is False
    assert "$NodeData" in result["error"]
    assert not (tmp_path / "must_not_exist.csv").exists()


def test_generated_post_plots_do_not_embed_a_figure_title():
    assert ".set_title(" not in _PLOT_SCRIPT


@pytest.mark.parametrize("cols", [0, -1, "not-an-integer"])
def test_montage_rejects_invalid_column_counts(cols, tmp_path):
    result = render_montage([], tmp_path / "grid.png", cols=cols)

    assert result["ok"] is False
    assert "positive integer" in result["error"]
