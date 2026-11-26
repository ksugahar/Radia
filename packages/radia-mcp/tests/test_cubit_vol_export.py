"""Guards: the Cubit-independent parts of the curving worker for Cubit .vol export.

Seam-consistent UV shifts, the geometric periodicity and pole tests, boundary
naming, the nine-field periodic contract, identification-section stripping,
point pairing on a translational fixture and the no-replace bundle publication.

Assumes: analytic stand-ins replace Cubit surfaces (a cylinder, a disc and a
plane with their closed-form parametrizations); no Cubit is imported.  The
named periodic fixture needs a Netgen build with
``Mesh.AddNamedPointIdentification``; without it that test skips visibly.
"""
from __future__ import annotations

import math

import pytest

from radia_mcp.cubit import periodic
from radia_mcp.cubit import vol_export as ve


# -- analytic Cubit stand-ins --------------------------------------------------

class _Surface:
    def __init__(self, evaluate, u_range, v_range, box, planar):
        self._evaluate, self._u, self._v, self._box, self._planar = evaluate, u_range, v_range, box, planar

    def get_param_range_U(self):
        return self._u

    def get_param_range_V(self):
        return self._v

    def position_from_u_v(self, u, v):
        return self._evaluate(u, v)

    def bounding_box(self):
        return self._box

    def is_planar(self):
        return self._planar


class _Cubit:
    def __init__(self, surfaces, periodic_ids, names=None):
        self._surfaces, self._periodic, self._names = surfaces, periodic_ids, names or {}

    def surface(self, sid):
        return self._surfaces[sid]

    def is_periodic(self, kind, sid):
        assert kind == "surface"
        return sid in self._periodic

    def get_entity_name(self, kind, sid):
        return self._names.get((kind, sid), f"Surface {sid}")


def _cylinder(u, v):
    return (math.cos(u), math.sin(u), v)


def _disc(u, v):  # polar parametrization: u = radius, v = angle
    return (u * math.cos(v), u * math.sin(v), 0.0)


def _plane(u, v):
    return (u, v, 0.0)


CUBIT = _Cubit(
    {1: _Surface(_cylinder, (0.0, 2 * math.pi), (0.0, 1.0), (-1, -1, 0, 1, 1, 1), False),
     2: _Surface(_disc, (0.0, 1.0), (0.0, 2 * math.pi), (-1, -1, 0, 1, 1, 0), True),
     3: _Surface(_plane, (0.0, 1.0), (0.0, 1.0), (0, 0, 0, 1, 1, 0), True),
     4: _Surface(_plane, (0.0, 1.0), (0.0, 1.0), (0, 0, 0, 1, 1, 0), True)},
    periodic_ids={1, 2, 4},
    names={("surface", 3): "inlet face"})


def test_seam_shift_moves_parameters_to_the_reference_side():
    period = 2 * math.pi
    cases = [(0.1, 6.2, 6.2 - period), (6.2, 0.1, 0.1 + period), (1.0, 1.5, 1.5), (0.1, 0.2, 0.2)]
    for u_ref, u, expected in cases:
        actual = ve._seam(u_ref, u, period)
        assert actual == pytest.approx(expected, abs=1e-15), (u_ref, u, actual, expected)
    assert ve._seam(0.1, 6.2, 0.0) == 6.2, "non-periodic directions must stay unchanged"


def test_periods_are_detected_geometrically():
    actual = {sid: ve._periods(CUBIT, sid) for sid in (1, 2, 3)}
    expected = {1: (2 * math.pi, 0.0), 2: (0.0, 2 * math.pi), 3: (0.0, 0.0)}
    for sid in expected:
        assert actual[sid] == pytest.approx(expected[sid]), (sid, actual[sid], expected[sid])


def test_reported_periodicity_without_a_closing_direction_fails_loudly():
    with pytest.raises(ve.ExportError, match="no parameter direction closes"):
        ve._periods(CUBIT, 4)


def test_pole_detection_finds_the_disc_centre_only():
    actual = {sid: ve._has_pole(CUBIT, sid) for sid in (1, 2, 3)}
    assert actual == {1: False, 2: True, 3: False}, actual


def test_boundary_names_follow_sideset_then_surface_name_without_spaces():
    sideset_label = {5: "outer wall", -7: "free set"}
    actual = [ve._boundary_name(CUBIT, bc, sideset_label) for bc in (5, -7, 3, 1)]
    expected = ["outer_wall", "free_set", "inlet_face", "Surface_1"]
    assert actual == expected, (actual, expected)
    assert ve._explicit_surface_label(CUBIT, 1) == "", "default Cubit names are not explicit labels"
    assert ve._explicit_surface_label(CUBIT, 3) == "inlet face"


def test_mesh_area_of_triangle_and_bilinear_quad():
    tri = [(0, 0, 0), (2, 0, 0), (0, 1, 0)]
    quad = [(0, 0, 0), (2, 0, 0), (2, 3, 0), (0, 3, 0)]
    assert ve._mesh_area(tri) == pytest.approx(1.0, rel=1e-15)
    assert ve._mesh_area(quad) == pytest.approx(6.0, rel=1e-13)


# -- periodic contract -----------------------------------------------------------

FIELDS = {"periodic_name": "kt", "periodic_master": "source_patch",
          "periodic_slave": "target_patch", "periodic_unit": "m",
          "periodic_patches": "4:1,2:3", "periodic_tx": -2.0, "periodic_ty": 0.0,
          "periodic_tz": 0.0, "periodic_match": 0.125}


def test_contract_needs_all_nine_fields_or_none():
    assert periodic.parse_contract({}) is None
    contract = periodic.parse_contract(FIELDS)
    assert contract.boundary_pairs == [(4, 1), (2, 3)], contract.boundary_pairs
    partial = dict(FIELDS)
    del partial["periodic_match"]
    with pytest.raises(periodic.PeriodicError, match="all nine"):
        periodic.parse_contract(partial)


@pytest.mark.parametrize("text", ["", "4:1,", "4-1", "4:1:2", "a:1", ":1", "4:"])
def test_patch_list_syntax_is_strict(text):
    with pytest.raises(periodic.PeriodicError):
        periodic.parse_patches(text)


@pytest.mark.parametrize("change, message", [
    ({"periodic_name": "other"}, "incomplete"),
    ({"periodic_unit": "mm"}, "incomplete"),
    ({"periodic_slave": "source_patch"}, "incomplete"),
    ({"periodic_match": 0.0}, "incomplete"),
    ({"periodic_tx": float("nan")}, "Nonfinite"),
    ({"periodic_patches": "4:1,4:3"}, "one-to-one"),
    ({"periodic_patches": "4:1,1:3"}, "overlap"),
])
def test_contract_validation_rejects(change, message):
    contract = periodic.parse_contract({**FIELDS, **change})
    with pytest.raises(periodic.PeriodicError, match=message):
        periodic.validate_contract(contract)


def test_geometry_serialization_drops_only_identification_sections():
    body = "mesh3d\npoints\n2\n0 0 0\n1 0 0\n"
    ident = ("identifications\n2\n1 7 1\n2 8 1\nidentificationtypes\n1\n2\n"
             "identificationnames\n1\nkt\n")
    tail = "endmesh\n"
    assert periodic.geometry_serialization(body + ident + tail) == periodic.geometry_serialization(body + tail)
    assert "identifications" not in periodic.geometry_serialization(body + ident + tail)
    with pytest.raises(periodic.PeriodicError, match="one named periodic group"):
        periodic.geometry_serialization(body + "identificationtypes\n2\n2 2\n" + tail)


# -- named periodic fixture (surface-only translational mesh) ---------------------

def _named_fixture(ngm):
    """Points 1-6 at x=2, 7-12 at x=0; FD1..4 carry BC 4, 2, 1, 3."""
    mesh = ngm.Mesh(dim=3)
    for side in range(2):
        for y in range(3):
            for z in range(2):
                mesh.Add(ngm.MeshPoint(ngm.Pnt(2.0 if side == 0 else 0.0, y, z)))
    for index, bc in enumerate((4, 2, 1, 3)):
        fd = mesh.Add(ngm.FaceDescriptor(surfnr=index + 1, domin=1, domout=0, bc=bc))
        mesh.SetBCName(index, "source_patch" if bc in (4, 2) else "target_patch")
        base = 1 + 2 * index if index < 2 else 7 + 2 * (index - 2)
        for tri in ((base, base + 2, base + 1), (base + 1, base + 2, base + 3)):
            mesh.Add(ngm.Element2D(fd, [ngm.PointId(p) for p in tri]))
    return mesh


def _save_text(tmp_path):
    def save(mesh):
        path = tmp_path / "serialized.vol"
        mesh.Save(str(path))
        return path.read_bytes().replace(b"\r\n", b"\n").decode()
    return save


def _netgen():
    ngm = pytest.importorskip("netgen.meshing", reason="Netgen is not installed")
    if not hasattr(ngm.Mesh, "AddNamedPointIdentification"):
        pytest.skip("this Netgen build lacks Mesh.AddNamedPointIdentification")
    return ngm


def test_named_periodic_pairs_points_by_translation(tmp_path):
    ngm = _netgen()
    mesh = _named_fixture(ngm)
    result = periodic.apply_identification(mesh, periodic.parse_contract(FIELDS), _save_text(tmp_path))
    actual = sorted(result["pairs"])
    expected = [(i, i + 6) for i in range(1, 7)]
    assert actual == expected, (actual, expected)
    assert result["descriptor_to_bc"] == {1: 4, 2: 2, 3: 1, 4: 3}, result["descriptor_to_bc"]
    assert mesh.GetNrIdentifications() == 1
    assert "PERIODIC" in str(mesh.GetIdentificationType(1)).upper()
    stored = sorted(tuple(int(getattr(x, "nr", x)) for x in row) for row in mesh.GetIdentifications())
    assert stored == [(i, i + 6, 1) for i in range(1, 7)], stored


@pytest.mark.parametrize("damage, message", [
    ("partial", "omits a labeled boundary"),
    ("unmatched", "candidates"),
    ("ambiguous", "candidates"),
    ("occupied", "Pre-existing"),
])
def test_named_periodic_rejections_leave_no_identification(tmp_path, damage, message):
    ngm = _netgen()
    mesh = _named_fixture(ngm)
    fields = dict(FIELDS)
    if damage == "partial":
        fields["periodic_patches"] = "4:1"
    elif damage == "unmatched":
        mesh.Points()[ngm.PointId(7)] = ngm.MeshPoint(ngm.Pnt(-1.0, 0.0, 0.0))
    elif damage == "ambiguous":
        mesh.Points()[ngm.PointId(8)] = ngm.MeshPoint(ngm.Pnt(0.0, 0.0, 0.0))
    elif damage == "occupied":
        mesh.AddPointIdentification(ngm.PointId(1), ngm.PointId(7), 4, 2)
    expected_groups = mesh.GetNrIdentifications()
    with pytest.raises(periodic.PeriodicError, match=message):
        periodic.apply_identification(mesh, periodic.parse_contract(fields), _save_text(tmp_path))
    assert mesh.GetNrIdentifications() == expected_groups, (mesh.GetNrIdentifications(), expected_groups)


def test_bundle_publication_never_replaces_and_needs_fresh_paths(tmp_path):
    final = tmp_path / "mesh.vol"
    periodic.publish_bundle(final, b"vol\n", b"{}\n")
    assert final.read_bytes() == b"vol\n"
    assert (tmp_path / "mesh.vol.json").read_bytes() == b"{}\n"
    assert not (tmp_path / "mesh.vol.partial.vol").exists()
    with pytest.raises(periodic.PeriodicError, match="already exists"):
        periodic.publish_bundle(final, b"other\n", b"{}\n")
    assert final.read_bytes() == b"vol\n", "an existing bundle must stay untouched"
    with pytest.raises(periodic.PeriodicError, match=".vol suffix"):
        periodic.publish_bundle(tmp_path / "mesh.txt", b"vol\n", b"{}\n")
