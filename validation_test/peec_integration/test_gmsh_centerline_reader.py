from pathlib import Path

import pytest

from radia.peec_mesh_import import GMSHCenterlineReader, read_gmsh_centerline


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "validation_test" / "peec_integration" / "verification" / "cubit_mesh_generation"


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
