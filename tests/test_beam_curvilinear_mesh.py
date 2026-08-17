import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")
orbit_module = pytest.importorskip("radia.accelerator_magnet_topopt")
beam_module = pytest.importorskip("radia.beam")
_native = pytest.importorskip("radia._radia_pybind")

Box = occ.Box
Pnt = occ.Pnt
PlanarDesignOrbit = orbit_module.PlanarDesignOrbit
bishop_rmf_frame = beam_module.bishop_rmf_frame
build_curvilinear_beam_mesh = beam_module.build_curvilinear_beam_mesh
project_design_orbit_gauge = beam_module.project_design_orbit_gauge


def _arc_orbit(n_station=7):
    radius = 0.25
    angle = np.linspace(0.0, 0.4, n_station)
    positions = np.column_stack(
        (radius * np.sin(angle), np.zeros_like(angle), radius * np.cos(angle))
    )
    tangents = np.column_stack(
        (np.cos(angle), np.zeros_like(angle), -np.sin(angle))
    )
    return PlanarDesignOrbit(
        positions=positions,
        tangents=tangents,
        magnetic_rigidity=1.5,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )


def test_planar_bishop_frame_keeps_y_normal_to_median_plane():
    orbit = _arc_orbit()
    frame = bishop_rmf_frame(
        orbit.positions,
        orbit.tangents,
        initial_horizontal=np.cross(orbit.bend_axis, orbit.tangents[0]),
    )

    np.testing.assert_allclose(
        frame.vertical,
        np.broadcast_to(orbit.bend_axis, frame.vertical.shape),
        atol=2.0e-15,
    )
    np.testing.assert_allclose(
        np.einsum("ij,ij->i", frame.tangent, frame.horizontal), 0.0, atol=2.0e-15
    )
    np.testing.assert_allclose(
        np.cross(frame.tangent, frame.horizontal), frame.vertical, atol=2.0e-15
    )
    np.testing.assert_allclose(
        np.linalg.norm(frame.horizontal, axis=1), 1.0, atol=2.0e-15
    )


def test_double_reflection_survives_straight_to_curved_transition():
    positions = np.array(
        [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.10, 0.0, 0.0], [0.15, 0.0, 0.01]]
    )
    tangents = np.array(
        [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.98, 0.0, 0.20]]
    )
    frame = bishop_rmf_frame(
        positions, tangents, initial_horizontal=[0.0, 0.0, -1.0]
    )

    assert np.all(np.isfinite(frame.horizontal))
    np.testing.assert_allclose(
        np.linalg.norm(frame.vertical, axis=1), 1.0, atol=2.0e-15
    )
    np.testing.assert_allclose(
        np.einsum("ij,ij->i", frame.tangent, frame.vertical), 0.0, atol=2.0e-15
    )


def test_curvilinear_beam_mesh_lofts_centered_rmf_sections():
    orbit = _arc_orbit()
    result = build_curvilinear_beam_mesh(
        orbit,
        half_width_m=0.008,
        half_height_m=0.004,
        maxh_m=0.006,
        curve_order=2,
    )

    assert result.mesh.dim == 3
    assert result.topology == "bishop-rmf-four-quad-loft-chain"
    assert result.mesh.ne == 4 * (len(result.longitudinal_stations_m) - 1)
    assert result.hex_connectivity.shape == (result.mesh.ne, 8)
    assert {str(element.type) for element in result.mesh.Elements()} == {"ET.HEX"}
    np.testing.assert_allclose(
        result.x_nodes_m, [-0.008, -0.004, 0.0, 0.004, 0.008]
    )
    np.testing.assert_allclose(result.y_nodes_m, [-0.004, 0.004])
    assert result.mesh.GetMaterials() == ("beam_tube",)
    assert set(result.mesh.GetBoundaries()) == {"beam_tube_boundary"}
    assert result.mesh.GetCurveOrder() == 2
    assert float(ng.Integrate(1.0, result.mesh)) > 0.0
    np.testing.assert_allclose(
        result.frame.vertical,
        np.broadcast_to(orbit.bend_axis, result.frame.vertical.shape),
        atol=2.0e-15,
    )
    np.testing.assert_allclose(
        result.frame.positions_m,
        orbit.position_at(result.longitudinal_stations_m),
        atol=0.0,
    )
    expected_vertices = (
        result.frame.positions_m[:, None, None, :]
        + result.x_nodes_m[None, None, :, None]
        * result.frame.horizontal[:, None, None, :]
        + result.y_nodes_m[None, :, None, None]
        * result.frame.vertical[:, None, None, :]
    )
    np.testing.assert_allclose(result.vertices_m, expected_vertices, atol=2.0e-15)
    # x=0 is the shared face between the two central strips; y=0 is inside it.
    np.testing.assert_allclose(
        np.mean(result.vertices_m[:, :, 2, :], axis=1),
        result.frame.positions_m,
        atol=2.0e-15,
    )


def test_curvilinear_mesh_accepts_p6_hcurl_vector_potential():
    result = build_curvilinear_beam_mesh(
        _arc_orbit(n_station=4),
        half_width_m=0.008,
        half_height_m=0.004,
        maxh_m=0.010,
        horizontal_subdivisions_per_macro_strip=2,
        vertical_layers=3,
        curve_order=2,
    )

    space = ng.HCurl(result.mesh, order=6)
    vector_potential = ng.GridFunction(space, name="A_p6")
    assert result.topology == (
        "bishop-rmf-four-macro-strip-8-x-3-y-loft-chain"
    )
    assert result.hexes_per_longitudinal_cell == 24
    assert result.mesh.ne == 24 * (len(result.longitudinal_stations_m) - 1)
    np.testing.assert_allclose(
        result.x_nodes_m,
        [-0.008, -0.006, -0.004, -0.002, 0.0, 0.002, 0.004, 0.006, 0.008],
    )
    np.testing.assert_allclose(
        result.y_nodes_m,
        [-0.004, -0.004 / 3.0, 0.004 / 3.0, 0.004],
    )
    assert not np.any(np.isclose(result.y_nodes_m, 0.0))
    orbit_x_index = int(np.argmin(np.abs(result.x_nodes_m)))
    np.testing.assert_allclose(
        np.mean(result.vertices_m[:, 1:3, orbit_x_index, :], axis=1),
        result.frame.positions_m,
        atol=2.0e-15,
    )
    assert space.ndof > 0
    assert vector_potential.space is space


@pytest.mark.parametrize("vertical_layers", [0, 2, 4])
def test_curvilinear_mesh_rejects_nonpositive_or_even_vertical_layers(
    vertical_layers,
):
    with pytest.raises(ValueError, match="positive odd integer"):
        build_curvilinear_beam_mesh(
            _arc_orbit(n_station=4),
            half_width_m=0.008,
            half_height_m=0.004,
            maxh_m=0.010,
            vertical_layers=vertical_layers,
            curve_order=2,
        )


@pytest.mark.parametrize("subdivisions", [0, -1])
def test_curvilinear_mesh_rejects_invalid_horizontal_subdivisions(subdivisions):
    with pytest.raises(ValueError, match="positive integer"):
        build_curvilinear_beam_mesh(
            _arc_orbit(n_station=4),
            half_width_m=0.008,
            half_height_m=0.004,
            maxh_m=0.010,
            horizontal_subdivisions_per_macro_strip=subdivisions,
            curve_order=2,
        )


def test_p6_hp_mesh_covers_twenty_millimetre_transverse_aperture():
    result = build_curvilinear_beam_mesh(
        _arc_orbit(n_station=4),
        half_width_m=0.020,
        half_height_m=0.020,
        maxh_m=0.050,
        horizontal_subdivisions_per_macro_strip=2,
        vertical_layers=9,
        curve_order=2,
    )
    space = ng.HCurl(result.mesh, order=6)

    assert result.hexes_per_longitudinal_cell == 72
    assert result.mesh.ne == 72 * (len(result.longitudinal_stations_m) - 1)
    assert result.x_nodes_m[0] == pytest.approx(-0.020)
    assert result.x_nodes_m[-1] == pytest.approx(0.020)
    assert result.y_nodes_m[0] == pytest.approx(-0.020)
    assert result.y_nodes_m[-1] == pytest.approx(0.020)
    assert not np.any(np.isclose(result.y_nodes_m, 0.0))
    assert np.max(np.diff(result.x_nodes_m)) <= 0.005 + 1.0e-15
    assert np.max(np.diff(result.y_nodes_m)) <= 0.0045
    assert space.ndof > 0


def test_p5_design_orbit_gauge_zeroes_As_Ay_without_changing_curl():
    result = build_curvilinear_beam_mesh(
        _arc_orbit(n_station=5),
        half_width_m=0.008,
        half_height_m=0.004,
        maxh_m=0.008,
        curve_order=2,
    )
    vector_potential = ng.CoefficientFunction(
        (
            0.2 + 0.1 * ng.x,
            0.3 + 0.2 * ng.z,
            -0.1 + 0.05 * ng.y,
        )
    )
    with ng.TaskManager():
        gauged = project_design_orbit_gauge(
            vector_potential,
            result,
            order=5,
            gauge_tolerance=2.0e-6,
        )

    assert np.max(np.abs(gauged.As_before_t_m)) > 0.1
    assert np.max(np.abs(gauged.Ay_before_t_m)) > 0.1
    segments = len(result.frame.positions_m) - 1
    assert len(gauged.frame.positions_m) == 12 * segments + 1
    # The gate is the independent audit frame, not the constrained stations.
    assert len(gauged.audit_frame.positions_m) == 31 * segments + 1
    assert gauged.maximum_orbit_gauge_residual_t_m < 2.0e-6
    assert gauged.maximum_constraint_residual_t_m < 2.0e-6
    assert gauged.schur_retained_rank <= gauged.constraint_count
    assert gauged.curl_change_l2_t_m32 < 1.0e-9


def test_python_rmf_matches_earlytimes_native_frame():
    positions = np.array(
        [[0.00, 0.0, 0.00], [0.05, 0.0, 0.01], [0.10, 0.0, 0.03]]
    )
    tangents = np.array(
        [[1.0, 0.0, 0.1], [1.0, 0.0, 0.2], [1.0, 0.0, 0.3]]
    )
    seed = np.array([0.0, 1.0, 0.0])
    frame = bishop_rmf_frame(
        positions, tangents, initial_horizontal=seed
    )

    mesh = Box(Pnt(-0.05, -0.05, -0.05), Pnt(0.15, 0.05, 0.08)).GenerateMesh(
        maxh=0.03
    )
    has_hcurl_input = "field_representation" in (
        _native._beam_grid_function_linear_map.__doc__ or ""
    )
    space = ng.HCurl(mesh, order=1) if has_hcurl_input else ng.VectorH1(mesh, order=1)
    field = ng.GridFunction(space)
    with ng.TaskManager():
        field.Set(ng.CoefficientFunction((0.0, 0.0, 0.0)))
        options = {
            "sample_radius_m": 1.0e-3,
            "periodic_frame": False,
        }
        if has_hcurl_input:
            options["field_representation"] = "hcurl_vector_potential"
        report = _native._beam_grid_function_linear_map(
            field,
            np.full(len(positions), 0.05),
            positions,
            tangents,
            1.5,
            seed,
            **options,
        )

    np.testing.assert_allclose(report["frame_tangent"], frame.tangent, atol=2.0e-15)
    np.testing.assert_allclose(
        report["frame_horizontal"], frame.horizontal, atol=2.0e-15
    )
    np.testing.assert_allclose(report["frame_vertical"], frame.vertical, atol=2.0e-15)
