import importlib.util
import sys
from pathlib import Path

import pytest

from radia.peec_mesh_import import GMSHCenterlineReader, read_gmsh_centerline


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "validation_test" / "peec_integration" / "verification" / "cubit_mesh_generation"
DOC_READER = (ROOT / "docs" / "peec_integration" / "demos" /
              "gmsh_models" / "gmsh_ascii_reader.py")


def test_gmsh_centerline_reader_loads_straight_wire_fixture():
    mesh_path = FIXTURES / "circular_wire_centerline.msh"

    reader = GMSHCenterlineReader(mesh_path)
    nodes, edges = reader.read()

    assert len(nodes) == 11
    assert len(edges) == 10
    assert reader.get_total_length() == pytest.approx(1.0)
    assert reader.summary()["physical_groups"] == {1: "wire_centerline"}


def test_read_gmsh_centerline_helper_matches_reader_api():
    mesh_path = FIXTURES / "circular_coil_centerline.msh"

    nodes, edges = read_gmsh_centerline(mesh_path)
    reader = GMSHCenterlineReader(mesh_path)
    reader.read()

    assert nodes == reader.nodes
    assert edges == reader.edge_elements
    assert len(edges) == 36
    assert reader.get_total_length() == pytest.approx(2.0 * 3.14159265359 * 0.05, rel=2e-3)


def test_gmsh_centerline_reader_rejects_other_v4_versions(tmp_path):
    mesh_path = tmp_path / "v40.msh"
    mesh_path.write_text(
        "$MeshFormat\n4.0 0 8\n$EndMeshFormat\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="standard is v4.1"):
        GMSHCenterlineReader(mesh_path).read()


def test_gmsh_centerline_reader_requires_mesh_format(tmp_path):
    mesh_path = tmp_path / "missing_format.msh"
    mesh_path.write_text("$Nodes\n0 0 0 0\n$EndNodes\n", encoding="utf-8")

    with pytest.raises(ValueError, match=r"\$MeshFormat section missing"):
        GMSHCenterlineReader(mesh_path).read()


def test_docs_gmsh_reader_rejects_multiple_physical_tags():
    spec = importlib.util.spec_from_file_location("docs_gmsh_ascii_reader",
                                                  DOC_READER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.modules.pop(spec.name, None)
    lines = [
        "$Entities",
        "1 0 0 0",
        "1 0 0 0 2 10 11",
        "$EndEntities",
    ]

    with pytest.raises(ValueError, match="multiple physical tags"):
        module._parse_entity_physicals(lines)
