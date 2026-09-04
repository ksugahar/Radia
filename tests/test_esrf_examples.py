"""Geometry/source contracts for the seven original ESRF Radia examples."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import radia as rad
from radia.vim._nonlinear import _bh_table_funcs

from radia.esrf_examples import (
    MM,
    MU0,
    build_esrf_coils,
    build_esrf_cubit_hdiv_iron,
    build_esrf_fixed_magnetization_blocks,
    build_esrf_fixed_magnetization_source,
    build_esrf_hdiv_iron,
    build_esrf_occ,
    esrf_fixed_magnetization_coefficient,
    export_esrf_cubit_assets,
    get_esrf_fixed_magnetization_by_material,
    get_esrf_cubit_mesh_policy,
    get_esrf_bh_table,
    get_esrf_example_spec,
    list_esrf_example_specs,
    validate_esrf_radia_reference,
)


def test_all_seven_specs_are_source_identified_and_si_scaled():
    specs = list_esrf_example_specs()
    assert [s.number for s in specs] == list(range(1, 8))
    assert [s.source_notebook for s in specs] == [f"Example#{i}.nb" for i in range(1, 8)]
    assert get_esrf_example_spec(1).parameters_si["size_m"] == (MM, MM, MM)
    assert get_esrf_example_spec(7).parameters_si["iron_length_m"] == 400 * MM


@pytest.mark.parametrize("number", range(1, 8))
def test_every_example_builds_nonempty_occ_source(number):
    pytest.importorskip("netgen.occ")
    shapes = build_esrf_occ(number)
    assert shapes
    assert all(len(list(shape.solids)) > 0 for shape in shapes.values())
    assert all(solid.mass > 0.0 for shape in shapes.values() for solid in shape.solids)


@pytest.mark.parametrize("number", (2, 5, 6, 7))
def test_current_source_paths_are_closed_and_finite(number):
    coils = build_esrf_coils(number)
    assert coils
    for coil in coils:
        assert coil.is_closed
        wires, current = coil.to_wire_segments(n_arc=12)
        assert wires
        assert np.isfinite(current)
        assert np.isfinite(np.asarray(wires)).all()


def test_expected_coil_multiplicity_after_symmetry():
    assert len(build_esrf_coils(2)) == 10
    assert len(build_esrf_coils(5)) == 1
    assert len(build_esrf_coils(6)) == 8
    assert len(build_esrf_coils(7)) == 8


def test_example3_fixed_magnetization_blocks_preserve_legacy_field_symmetry():
    blocks = build_esrf_fixed_magnetization_blocks(3)
    assert len(blocks) == 24
    assert [block.material_name for block in blocks] == [
        f"pm_{index:03d}" for index in range(24)
    ]
    assert {block.symmetry_path for block in blocks} == {
        ("base",),
        ("base", "mirror_x"),
        ("base", "mirror_z"),
        ("base", "mirror_x", "mirror_z"),
        ("base", "mirror_y"),
        ("base", "mirror_x", "mirror_y"),
        ("base", "mirror_z", "mirror_y"),
        ("base", "mirror_x", "mirror_z", "mirror_y"),
    }
    m = np.asarray([block.magnetization_A_m for block in blocks])
    remanence = get_esrf_example_spec(3).parameters_si["remanence_T"]
    np.testing.assert_allclose(np.linalg.norm(m, axis=1), remanence / MU0)
    np.testing.assert_allclose(m[:, (0, 2)], 0.0, atol=1.0e-12)
    # The legacy source uses TrfZerPerp at x/y and TrfZerPara at z.  Thus the
    # x image preserves My, the z image reverses My, and the y image reverses
    # My.  This is a B-parity prescription, not one generic CAD reflection.
    signs = np.sign(m[:, 1]).astype(int).tolist()
    assert signs == [
        -1, 1, -1, -1, 1, -1,
        1, -1, 1, 1, -1, 1,
        1, -1, 1, 1, -1, 1,
        -1, 1, -1, -1, 1, -1,
    ]
    np.testing.assert_allclose(m.sum(axis=0), 0.0, atol=1.0e-7)


def test_example3_compound_keeps_every_fixed_source_solid():
    shapes = build_esrf_occ(3, include_coils=False)
    blocks = build_esrf_fixed_magnetization_blocks(3)
    assert len(list(shapes["magnet"].solids)) == len(blocks)
    assert sum(block.shape.mass for block in blocks) == pytest.approx(
        shapes["magnet"].mass, rel=1.0e-12
    )


def test_example3_fixed_source_material_map_is_complete_and_unambiguous():
    by_material = get_esrf_fixed_magnetization_by_material(3)
    assert set(by_material) == {f"pm_{index:03d}" for index in range(24)}
    vectors = np.asarray(list(by_material.values()))
    assert np.count_nonzero(vectors[:, 1] > 0.0) == 12
    assert np.count_nonzero(vectors[:, 1] < 0.0) == 12
    np.testing.assert_allclose(vectors[:, (0, 2)], 0.0, atol=1.0e-12)


def test_fixed_pm_source_requires_exact_source_mesh_labels():
    class MeshLabels:
        def GetMaterials(self):
            return ("pm_000",)

    with pytest.raises(ValueError, match="must exactly match"):
        esrf_fixed_magnetization_coefficient(MeshLabels(), 3)


def test_example1_source_constructor_is_native_hdiv_path():
    import ngsolve as ng
    from netgen.csg import CSGeometry, OrthoBrick, Pnt

    geometry = CSGeometry()
    geometry.Add(OrthoBrick(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2)).mat(
        "pm_000"
    ))
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.3))
    with ng.TaskManager():
        coefficient = esrf_fixed_magnetization_coefficient(mesh, 1)
        source = build_esrf_fixed_magnetization_source(
            mesh, 1, order=1, field_cf_algorithm="direct"
        )
    sampled = np.asarray(coefficient(mesh(0.0, 0.0, 0.0)), dtype=float)
    expected = np.asarray(get_esrf_fixed_magnetization_by_material(1)["pm_000"])
    np.testing.assert_allclose(sampled, expected, rtol=0.0, atol=1.0e-7)
    assert source.permanent_magnet_model == "fixed-given"
    assert source.field_cf_algorithm == "direct"


def test_example7_coils_preserve_notebook_minus_45_degree_phase():
    coils = build_esrf_coils(7)
    assert [coil.current for coil in coils] == pytest.approx([
        8 * 533.3, -8 * 533.3, 8 * 533.3, -8 * 533.3,
        7 * 533.3, -7 * 533.3, 7 * 533.3, -7 * 533.3,
    ])
    centres = []
    for coil in coils:
        endpoints = np.asarray([
            point
            for segment in coil.segments
            for point in (segment.start_pos, segment.end_pos)
        ])
        centres.append(endpoints.mean(axis=0))
    centres = np.asarray(centres)
    expected = 178 * MM / np.sqrt(2.0)
    assert np.abs(centres[:, 0]) == pytest.approx(np.full(8, expected))
    assert np.abs(centres[:, 1]) == pytest.approx(np.full(8, expected))
    assert centres[:, 2] == pytest.approx(np.zeros(8), abs=1.0e-14)
    # The unrotated source centre is on +y (90 degrees); the notebook's
    # -45-degree transform therefore places the first winding at +45 degrees.
    expected_angles = np.array([45.0, 135.0, -135.0, -45.0] * 2)
    actual_angles = np.degrees(np.arctan2(centres[:, 1], centres[:, 0]))
    assert actual_angles == pytest.approx(expected_angles, abs=1.0e-12)
    radii = [39.1 * MM] * 4 + [53.3 * MM] * 4
    heights = [80 * MM] * 4 + [70 * MM] * 4
    for coil, radius, height in zip(coils, radii, heights):
        assert coil._width == pytest.approx(14.2 * MM)
        assert coil._height == pytest.approx(height)
        assert [segment.length for segment in coil.segments[::2]] == \
            pytest.approx([402 * MM, 402 * MM])
        assert [segment.radius for segment in coil.segments[1::2]] == \
            pytest.approx([radius, radius])


def test_example6_coils_follow_diagonal_alternating_quadrupole_symmetry():
    coils = build_esrf_coils(6)
    assert [coil.current for coil in coils] == pytest.approx([
        -1560.0, 1560.0, -1560.0, 1560.0,
        -960.0, 960.0, -960.0, 960.0,
    ])
    assert [coil._width for coil in coils] == pytest.approx(
        [13 * MM] * 4 + [16 * MM] * 4)
    assert [coil._height for coil in coils] == pytest.approx(
        [40 * MM] * 4 + [20 * MM] * 4)
    centres = []
    for coil in coils:
        endpoints = np.asarray([
            point
            for segment in coil.segments
            for point in (segment.start_pos, segment.end_pos)
        ])
        centres.append(endpoints.mean(axis=0))
    centres = np.asarray(centres)
    assert centres[:, 0] == pytest.approx(np.zeros(8), abs=1.0e-14)
    for pack, radius in ((centres[:4], 50 * MM), (centres[4:], 60 * MM)):
        expected = radius / np.sqrt(2.0)
        assert np.abs(pack[:, 1]) == pytest.approx(np.full(4, expected))
        assert np.abs(pack[:, 2]) == pytest.approx(np.full(4, expected))
        assert np.sign(pack[:, 1]).tolist() == [-1.0, -1.0, 1.0, 1.0]
        assert np.sign(pack[:, 2]).tolist() == [1.0, -1.0, -1.0, 1.0]


def test_example6_iron_contains_both_halves_of_all_four_poles():
    iron = build_esrf_occ(6, include_coils=False)["iron"]
    assert len(list(iron.solids)) == 40
    # The x-z-2 mm ObjCutMag end plane is a physical chamfer.  Radia's
    # transformed surface polygon, rather than ObjGeoVol on the mutated parent
    # container, is the source of truth for this cut volume.
    assert iron.mass == pytest.approx(1024780.794801163e-9, rel=2.0e-10)
    end_vertices = np.asarray([
        [vertex.p.x, vertex.p.y, vertex.p.z]
        for vertex in iron.vertices
        if abs(abs(vertex.p.x) - 30 * MM) < 1.0e-7
    ])
    assert np.min(np.linalg.norm(end_vertices[:, 1:], axis=1)) == pytest.approx(
        28 * MM, abs=2.0e-10)


def test_example7_iron_preserves_original_symmetric_end_chamfer():
    spec = get_esrf_example_spec(7)
    assert spec.parameters_si["end_chamfer_m"] == 7 * MM
    assert "symmetric 7 mm 45-degree pole-end chamfer" in spec.cad_fidelity
    iron = build_esrf_occ(7, include_coils=False)["iron"]
    assert len(list(iron.solids)) == 4
    assert iron.mass == pytest.approx(0.044355134334166775, rel=2.0e-12)
    vertices = np.asarray([
        [vertex.p.x, vertex.p.y, vertex.p.z]
        for vertex in iron.vertices
    ])
    assert np.max(np.abs(vertices[:, 2])) == pytest.approx(200 * MM)
    chamfered = np.abs(vertices[:, 2]) < 200 * MM - 1.0e-12
    assert int(np.count_nonzero(chamfered)) > 0
    expected_end = (
        193 * MM
        + (np.abs(vertices[chamfered, 0])
           + np.abs(vertices[chamfered, 1])
           - 36 * np.sqrt(2.0) * MM) / np.sqrt(2.0)
    )
    assert np.abs(vertices[chamfered, 2]) == pytest.approx(
        expected_end, abs=5.0e-15)


def test_invalid_example_number_is_rejected():
    with pytest.raises(ValueError, match="1..7"):
        get_esrf_example_spec(8)


def test_hdiv_iron_removes_same_material_internal_interfaces():
    raw = build_esrf_occ(5, include_coils=False)["iron"]
    normalized = build_esrf_hdiv_iron(5)
    assert len(list(raw.solids)) == 12
    assert len(list(normalized.solids)) == 1
    assert normalized.mass > 0.0


def test_example5_hdiv_image_domain_is_exact_quarter():
    full = build_esrf_hdiv_iron(5)
    quarter = build_esrf_hdiv_iron(5, image="+x-z")
    assert len(list(quarter.solids)) == 1
    assert quarter.mass == pytest.approx(0.25 * full.mass, rel=1.0e-9)
    vertices = np.asarray(
        [[vertex.p.x, vertex.p.y, vertex.p.z] for vertex in quarter.vertices]
    )
    assert vertices[:, 0].min() == pytest.approx(0.0, abs=1.0e-12)
    assert vertices[:, 2].min() == pytest.approx(0.0, abs=1.0e-12)
    assert vertices[:, 0].max() > 0.0 and vertices[:, 2].max() > 0.0


def test_example5_cubit_hdiv_parts_preserve_exact_quarter_volume():
    partitioned = build_esrf_cubit_hdiv_iron(5, image="+x-z")
    normalized = build_esrf_hdiv_iron(5, image="+x-z")
    assert len(list(partitioned.solids)) == 6
    assert len(list(normalized.solids)) == 1
    assert partitioned.mass == pytest.approx(normalized.mass, rel=1.0e-12)


def test_hdiv_image_contract_rejects_unregistered_reductions():
    with pytest.raises(ValueError, match="example 5 only"):
        build_esrf_hdiv_iron(6, image="+x-z")


@pytest.mark.parametrize("number", [3, 5, 6, 7])
def test_original_iron_law_is_converted_to_monotone_total_bh(number):
    options = {"sample_count": 31} if number in (5, 6) else {}
    table = np.asarray(get_esrf_bh_table(number, **options), dtype=float)
    assert table.ndim == 2 and table.shape[1] == 2
    assert np.isfinite(table).all()
    assert np.all(np.diff(table[:, 0]) > 0.0)
    assert np.all(np.diff(table[:, 1]) > 0.0)
    assert table[0] == pytest.approx([0.0, 0.0])


def test_example7_bdm_pchip_reproduces_tabulated_knots():
    table = np.asarray(get_esrf_bh_table(7), dtype=float)
    _, b_interpolator, _, _, _ = _bh_table_funcs(
        table[:, 0], table[:, 1])

    # HDiv uses a monotone PCHIP B(H), whereas MatSatIsoTab evaluates a
    # piecewise-linear law between rows.  The common physical contract is the
    # source table itself, so test every selected knot instead of requiring a
    # false inter-knot identity between distinct interpolation policies.
    selected = table[[1, 5, 10, 15, 20, 25, 30, -1]]
    np.testing.assert_allclose(
        b_interpolator(selected[:, 0]), selected[:, 1], rtol=2.0e-12, atol=1.0e-12
    )


def test_example1_published_radia_field_vector():
    result = validate_esrf_radia_reference(1)
    assert result["passed"]
    assert result["max_abs_error_T"] < 1.0e-5


def test_example2_coilbuilder_matches_native_racetrack_profile():
    result = validate_esrf_radia_reference(2, n_points=31)
    assert result["passed"]
    assert result["profile_relative_l2"] < 0.01


@pytest.mark.parametrize("number", [2, 5, 6, 7])
def test_manifest_preserves_physical_coil_source_without_meshing(number, tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(number, tmp_path)
    assert manifest["coil_sources"]
    assert all(source["closed"] for source in manifest["coil_sources"])
    assert all(source["closure_gap_m"] < 1.0e-12
               for source in manifest["coil_sources"])
    assert "coil" not in manifest["solver_step_files"]
    journal = Path(manifest["journal"]).read_text(encoding="utf-8")
    if number == 2:
        assert "mesh volume" not in journal
        assert "export netgen" not in journal
    else:
        assert "coil.step" not in journal


def test_example5_coilbuilder_source_metadata(tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(5, tmp_path)
    assert len(manifest["coil_sources"]) == 1
    source = manifest["coil_sources"][0]
    assert source["current_A"] == -2000.0
    assert source["segment_count"] == 8


def test_example3_cubit_export_separates_response_and_fixed_pm_source_meshes(tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(3, tmp_path)
    assert set(manifest["solver_step_files"]) == {"iron"}
    assert manifest["solver_vol"].endswith("model.vol")
    assert len(manifest["fixed_magnetization_sources"]) == 24
    assert len(manifest["fixed_magnetization_source_step_files"]) == 24
    assert manifest["fixed_magnetization_source_vol"].endswith(
        "magnet_source.vol"
    )
    assert manifest["fixed_magnetization_source_journal"]
    assert {source["material"] for source in manifest[
        "fixed_magnetization_sources"
    ]} == set(manifest["fixed_magnetization_source_step_files"])
    assert {source["model"] for source in manifest[
        "fixed_magnetization_sources"
    ]} == {"fixed-given MagnetizationSource"}
    source_journal = Path(
        manifest["fixed_magnetization_source_journal"]
    ).read_text(encoding="utf-8")
    assert 'block 1 name "pm_000"' in source_journal
    assert 'block 24 name "pm_023"' in source_journal
    assert 'export netgen "' in source_journal
    assert "magnet_source.vol" in source_journal


def test_example3_fixed_pm_source_mesh_size_is_independent_and_recorded(tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(
        3,
        tmp_path,
        mesh_size_m=0.01,
        fixed_magnetization_source_mesh_size_m=0.05,
        order=2,
    )
    assert manifest["mesh_size_m"] == pytest.approx(0.01)
    assert manifest["fixed_magnetization_source_mesh_size_m"] == pytest.approx(0.05)
    journal = Path(manifest["fixed_magnetization_source_journal"]).read_text(
        encoding="utf-8"
    )
    assert "volume all size 0.05" in journal


def test_source_only_example2_has_no_solver_volume_mesh(tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(2, tmp_path)
    journal = Path(manifest["journal"]).read_text(encoding="utf-8")
    assert manifest["solver_step_files"] == {}
    assert manifest["solver_vol"] is None
    assert "coil.step" in journal
    assert "mesh volume" not in journal
    assert "export netgen" not in journal


@pytest.mark.parametrize("number", range(1, 8))
def test_all_esrf_examples_have_cubit_hex_preferred_mesh_policy(number):
    policy = get_esrf_cubit_mesh_policy(number)
    assert policy["all_examples_use_cubit"]
    assert policy["preferred_volume_family"] == "HEX"
    assert policy["regions"]["air"]["fem_preferred_family"] == "HEX"
    assert policy["regions"]["air"]["fem_allowed_families"] == ["HEX", "TET"]
    assert "WEDGE" in policy["regions"]["transition"]["allowed_families"]
    assert policy["regions"]["coil"]["solver_mesh"] is False
    assert policy["regions"]["coil"]["source"] == "CoilBuilder"
    assert policy["regions"]["permanent_magnet"]["source_mesh"] is True
    assert policy["regions"]["permanent_magnet"]["response_unknown"] is False


@pytest.mark.parametrize("number", [3, 5, 6])
def test_every_esrf_yoke_requires_cubit_hex(number):
    policy = get_esrf_cubit_mesh_policy(number)
    assert policy["regions"]["iron"] == {
        "present": True,
        "required_family": "HEX",
        "scheme": "auto",
    }


def test_example7_explicitly_uses_curved_tet_until_cad_partitioning_exists():
    policy = get_esrf_cubit_mesh_policy(7)
    assert policy["regions"]["iron"] == {
        "present": True,
        "required_family": "TET",
        "scheme": "tetmesh",
    }
    assert policy["solver_routes"]["bdm2_ima"]["iron_family"] == "TET"
    assert policy["solver_routes"]["nonlinear_hcurl_fem"]["iron_family"] == "TET"
    assert "cannot map/submap" in policy["cad_volume_fallback_reason"]


def test_cubit_journal_uses_hex_capable_auto_for_c_yoke(tmp_path):
    pytest.importorskip("netgen.occ")
    manifest = export_esrf_cubit_assets(5, tmp_path)
    journal = Path(manifest["journal"]).read_text(encoding="utf-8")
    assert "volume all scheme auto" in journal
    assert "volume all scheme tetmesh" not in journal
