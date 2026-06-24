import math

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


TET_VOL_SURFACE_UV = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
1
1 1 0 1 1
surfaceelementsuv
4
1 1 1 0 3 1 2 3 0 0 1 0 0 1
1 1 1 0 3 1 4 2 0 0 0 1 1 0
1 1 1 0 3 2 4 3 1 0 0 1 0 0
1 1 1 0 3 3 4 1 0 1 0 0 1 0
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


TET_VOL_FOUR_BOUNDARIES = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
4
1 1 0 1 1
2 1 0 1 1
3 1 0 1 1
4 1 0 1 1
surfaceelements
4
1 1 1 0 3 1 2 3
2 2 1 0 3 1 4 2
3 3 1 0 3 2 4 3
4 4 1 0 3 3 4 1
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
4
1 face123
2 face142
3 face243
4 face341
endmesh
"""


FOUR_TET_WITH_INTERIOR_NODE_VOL = """\
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
1 1 1 0 3 1 3 4
1 1 1 0 3 2 4 3
volumeelements
4
1 4 1 2 3 5
1 4 1 4 2 5
1 4 1 3 4 5
1 4 2 4 3 5
points
5
0 0 0
1 0 0
0 1 0
0 0 1
0.25 0.25 0.25
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


OPEN_SURFACE_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
1
1 1 0 1 1
surfaceelements
1
1 1 1 0 3 1 2 3
volumeelements
0
points
3
0 0 0
1 0 0
0 1 0
pointelements
0
materials
0
bcnames
1
1 patch
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


def test_surfaceelementsuv_is_parsed_as_triangles():
    mesh = parse_netgen_tri_tet_vol(TET_VOL_SURFACE_UV)

    assert mesh.summary()["surface_triangles"] == 4
    assert mesh.surface_triangles[0].nodes == (1, 2, 3)
    assert mesh.total_volume() == pytest.approx(1.0 / 6.0)
    assert mesh.total_surface_area() == pytest.approx(1.5 + 0.5 * 3**0.5)


def test_fem_bem_trace_view_preserves_one_based_connectivity():
    view = parse_netgen_tri_tet_vol(TET_VOL).fem_bem_trace_view()

    assert view["tetrahedra"] == [[1, 2, 3, 4]]
    assert view["tetrahedron_material_numbers"] == [1]
    assert view["surface_triangles"][0] == [1, 2, 3]
    assert view["surface_boundary_numbers"] == [1, 1, 1, 1]
    assert view["trace_node_ids"] == [1, 2, 3, 4]
    assert view["trace_node_ids_by_boundary_number"] == {1: [1, 2, 3, 4]}
    assert view["total_volume"] == pytest.approx(1.0 / 6.0)
    assert view["total_surface_area"] == pytest.approx(1.5 + 0.5 * 3**0.5)
    assert view["policy"] == "netgen_vol_tri_tet_only_shared_one_based_nodes"


def test_trace_node_ids_by_boundary_number():
    mesh = parse_netgen_tri_tet_vol(TET_VOL_FOUR_BOUNDARIES)

    assert mesh.trace_node_ids_by_boundary_number() == {
        1: (1, 2, 3),
        2: (1, 2, 4),
        3: (2, 3, 4),
        4: (1, 3, 4),
    }
    assert mesh.fem_bem_trace_view()["trace_node_ids_by_boundary_number"] == {
        1: [1, 2, 3],
        2: [1, 2, 4],
        3: [2, 3, 4],
        4: [1, 3, 4],
    }


def test_first_order_fem_bem_topology_for_unit_tetrahedron():
    topology = parse_netgen_tri_tet_vol(TET_VOL).first_order_fem_bem_topology()

    assert topology["policy"] == "first_order_h1_p1_hcurl_nedelec0_bem_p1_rwg_only"
    assert topology["h1"]["node_ids"] == [1, 2, 3, 4]
    assert topology["h1"]["trace_node_ids"] == [1, 2, 3, 4]
    assert topology["hcurl"]["edges"] == [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
    assert topology["hcurl"]["tet_edges"] == [[1, 2, 3, 4, 5, 6]]
    assert topology["hcurl"]["tet_edge_signs"] == [[1, 1, 1, 1, 1, 1]]
    assert topology["scalar_bem"]["global_node_ids"] == [1, 2, 3, 4]
    assert topology["rwg"]["dof_edge_ids"] == [1, 2, 3, 4, 5, 6]
    assert topology["rwg"]["hcurl_edge_ids"] == [1, 2, 3, 4, 5, 6]
    assert topology["trace"]["h1_to_scalar_bem_rows"] == [1, 2, 3, 4]
    assert topology["trace"]["h1_to_scalar_bem_cols"] == [1, 2, 3, 4]
    assert topology["trace"]["rwg_to_hcurl_edge_ids"] == [1, 2, 3, 4, 5, 6]


def test_first_order_topology_compacts_boundary_nodes_with_interior_node():
    topology = parse_netgen_tri_tet_vol(FOUR_TET_WITH_INTERIOR_NODE_VOL).first_order_fem_bem_topology()

    assert topology["h1"]["node_ids"] == [1, 2, 3, 4, 5]
    assert topology["h1"]["trace_node_ids"] == [1, 2, 3, 4]
    assert topology["scalar_bem"]["node_ids"] == [1, 2, 3, 4]
    assert topology["scalar_bem"]["global_node_ids"] == [1, 2, 3, 4]
    assert topology["trace"]["h1_to_scalar_bem_rows"] == [1, 2, 3, 4]
    assert topology["trace"]["h1_to_scalar_bem_cols"] == [1, 2, 3, 4]
    assert len(topology["hcurl"]["edges"]) == 10
    assert topology["rwg"]["hcurl_edge_ids"] == [1, 2, 3, 5, 6, 8]


def test_first_order_topology_balances_rwg_orientation_with_hcurl_trace():
    topology = parse_netgen_tri_tet_vol(FOUR_TET_WITH_INTERIOR_NODE_VOL).first_order_fem_bem_topology()

    hcurl_edges = {tuple(edge) for edge in topology["hcurl"]["edges"]}
    rwg_edges = {tuple(edge) for edge in topology["rwg"]["dof_edges_global"]}
    edge_signs = {edge_id: [] for edge_id in topology["rwg"]["dof_edge_ids"]}
    for tri_edges, tri_signs in zip(topology["rwg"]["tri_edges"], topology["rwg"]["tri_edge_signs"]):
        for edge_id, sign in zip(tri_edges, tri_signs):
            edge_signs[edge_id].append(sign)

    assert sorted(hcurl_edges - rwg_edges) == [(1, 5), (2, 5), (3, 5), (4, 5)]
    assert rwg_edges.issubset(hcurl_edges)
    assert all(sorted(signs) == [-1, 1] for signs in edge_signs.values())
    assert -1 in {sign for row in topology["hcurl"]["tet_edge_signs"] for sign in row}
    assert topology["trace"]["rwg_to_hcurl_edge_ids"] == [1, 2, 3, 5, 6, 8]


def test_geometry_metrics_for_single_tetrahedron():
    mesh = parse_netgen_tri_tet_vol(TET_VOL)

    assert mesh.bounding_box() == {"x": (0.0, 1.0), "y": (0.0, 1.0), "z": (0.0, 1.0)}
    assert mesh.tetrahedron_signed_volumes() == pytest.approx((1.0 / 6.0,))
    assert mesh.tetrahedron_volumes() == pytest.approx((1.0 / 6.0,))
    edge_lengths = mesh.tetrahedron_edge_lengths()
    assert len(edge_lengths) == 1
    assert edge_lengths[0] == pytest.approx((1.0, 1.0, 1.0, 2**0.5, 2**0.5, 2**0.5))
    assert mesh.tetrahedron_edge_length_ratios() == pytest.approx((2**0.5,))
    assert mesh.tetrahedron_edge_length_summary() == {
        "tetrahedra": 1,
        "min_edge": pytest.approx(1.0),
        "max_edge": pytest.approx(2**0.5),
        "mean_edge": pytest.approx((3.0 + 3.0 * 2**0.5) / 6.0),
        "max_edge_ratio": pytest.approx(2**0.5),
    }
    assert mesh.total_volume() == pytest.approx(1.0 / 6.0)
    assert sorted(mesh.surface_triangle_areas()) == pytest.approx([0.5, 0.5, 0.5, 0.5 * 3**0.5])
    assert mesh.surface_area_by_boundary_number() == {1: pytest.approx(1.5 + 0.5 * 3**0.5)}


def test_tetrahedron_quality_rows_for_right_and_equilateral_tets():
    right = parse_netgen_tri_tet_vol(TET_VOL)
    row = right.tetrahedron_quality_rows()[0]
    expected_surface_area = 1.5 + 0.5 * math.sqrt(3.0)
    expected_inradius = 0.5 / expected_surface_area
    expected_circumradius = math.sqrt(3.0) / 2.0
    expected_quality = 3.0 * expected_inradius / expected_circumradius

    assert row["volume"] == pytest.approx(1.0 / 6.0)
    assert row["surface_area"] == pytest.approx(expected_surface_area)
    assert row["inradius"] == pytest.approx(expected_inradius)
    assert row["circumradius"] == pytest.approx(expected_circumradius)
    assert row["radius_ratio_quality"] == pytest.approx(expected_quality)
    assert row["edge_ratio"] == pytest.approx(math.sqrt(2.0))
    assert right.tetrahedron_quality_summary()["min_radius_ratio_quality"] == pytest.approx(expected_quality)

    h = math.sqrt(3.0) / 2.0
    z = math.sqrt(2.0 / 3.0)
    equilateral_vol = TET_VOL.replace(
        "0 0 0\n1 0 0\n0 1 0\n0 0 1",
        f"0 0 0\n1 0 0\n0.5 {h} 0\n0.5 {math.sqrt(3.0) / 6.0} {z}",
    )
    equilateral = parse_netgen_tri_tet_vol(equilateral_vol)
    eq_summary = equilateral.tetrahedron_quality_summary()
    assert eq_summary["min_radius_ratio_quality"] == pytest.approx(1.0)
    assert eq_summary["max_radius_ratio_quality"] == pytest.approx(1.0)
    assert eq_summary["max_edge_ratio"] == pytest.approx(1.0)


def test_surface_closure_summary_for_single_tetrahedron():
    mesh = parse_netgen_tri_tet_vol(TET_VOL)

    assert mesh.surface_vector_area() == pytest.approx((0.0, 0.0, 0.0))
    assert mesh.surface_signed_volume_from_triangles() == pytest.approx(-1.0 / 6.0)
    summary = mesh.surface_closure_summary()
    assert summary["surface_triangles"] == 4
    assert summary["tetrahedra"] == 1
    assert summary["surface_vector_area_norm_over_area"] == pytest.approx(0.0)
    assert summary["surface_signed_volume"] == pytest.approx(-1.0 / 6.0)
    assert summary["surface_abs_volume"] == pytest.approx(mesh.total_volume())
    assert summary["surface_abs_volume_rel_error"] == pytest.approx(0.0)
    assert summary["boundary_orientation"] == "inward"


def test_surface_edge_manifold_summary_for_closed_and_open_surfaces():
    closed = parse_netgen_tri_tet_vol(TET_VOL).surface_edge_manifold_summary()

    assert closed == {
        "trace_nodes": 4,
        "surface_edges": 6,
        "surface_triangles": 4,
        "closed_edges": 6,
        "open_edges": 0,
        "is_closed_manifold": True,
        "euler_characteristic": 2,
    }

    open_patch = parse_netgen_tri_tet_vol(OPEN_SURFACE_VOL).surface_edge_manifold_summary()
    assert open_patch == {
        "trace_nodes": 3,
        "surface_edges": 3,
        "surface_triangles": 1,
        "closed_edges": 0,
        "open_edges": 3,
        "is_closed_manifold": False,
        "euler_characteristic": 1,
    }


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
