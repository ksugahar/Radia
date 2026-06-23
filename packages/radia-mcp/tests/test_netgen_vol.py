import pytest

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol


TET_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
1
1 1 0 1 1
surfaceelements
4
1 1 1 0 3 1 2 3
1 1 1 0 3 1 4 2
1 1 1 0 3 2 4 3
1 1 1 0 3 3 4 1
volumeelements
1
1 4 1 2 3 4
points
4
0 0 0
1 0 0
0 1 0
0 0 1
pointelements
0
materials
1
1 air
bcnames
1
1 outer
endmesh
"""


def test_parse_tri_tet_vol_summary():
    mesh = parse_netgen_tri_tet_vol(TET_VOL)

    assert mesh.summary() == {
        "points": 4,
        "surface_triangles": 4,
        "tetrahedra": 1,
        "materials": 1,
        "boundary_names": 1,
    }
    assert mesh.trace_node_ids() == (1, 2, 3, 4)


def test_fem_bem_trace_view_preserves_one_based_connectivity():
    view = parse_netgen_tri_tet_vol(TET_VOL).fem_bem_trace_view()

    assert view["tetrahedra"] == [[1, 2, 3, 4]]
    assert view["tetrahedron_material_numbers"] == [1]
    assert view["surface_triangles"][0] == [1, 2, 3]
    assert view["surface_boundary_numbers"] == [1, 1, 1, 1]
    assert view["trace_node_ids"] == [1, 2, 3, 4]
    assert view["total_volume"] == pytest.approx(1.0 / 6.0)
    assert view["total_surface_area"] == pytest.approx(1.5 + 0.5 * 3**0.5)
    assert view["policy"] == "netgen_vol_tri_tet_only_shared_one_based_nodes"


def test_geometry_metrics_for_single_tetrahedron():
    mesh = parse_netgen_tri_tet_vol(TET_VOL)

    assert mesh.bounding_box() == {"x": (0.0, 1.0), "y": (0.0, 1.0), "z": (0.0, 1.0)}
    assert mesh.tetrahedron_signed_volumes() == pytest.approx((1.0 / 6.0,))
    assert mesh.tetrahedron_volumes() == pytest.approx((1.0 / 6.0,))
    assert mesh.total_volume() == pytest.approx(1.0 / 6.0)
    assert sorted(mesh.surface_triangle_areas()) == pytest.approx([0.5, 0.5, 0.5, 0.5 * 3**0.5])
    assert mesh.surface_area_by_boundary_number() == {1: pytest.approx(1.5 + 0.5 * 3**0.5)}


def test_quad_surface_is_rejected_not_split():
    bad = TET_VOL.replace("1 1 1 0 3 1 2 3", "1 1 1 0 4 1 2 3 4", 1)

    with pytest.raises(ValueError, match="tri/tet-only policy rejected surface element"):
        parse_netgen_tri_tet_vol(bad)


def test_hex_volume_is_rejected_not_tetized():
    bad = TET_VOL.replace("1 4 1 2 3 4", "1 8 1 2 3 4 1 2 3 4", 1)

    with pytest.raises(ValueError, match="tri/tet-only policy rejected volume element"):
        parse_netgen_tri_tet_vol(bad)


def test_out_of_range_node_is_rejected():
    bad = TET_VOL.replace("1 4 1 2 3 4", "1 4 1 2 3 5", 1)

    with pytest.raises(ValueError, match="references node 5"):
        parse_netgen_tri_tet_vol(bad)
