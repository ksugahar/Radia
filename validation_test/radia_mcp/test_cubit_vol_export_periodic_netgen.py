"""Guards: named periodic point pairing of the Cubit .vol worker on a real Netgen mesh.

Assumes: a Netgen build with ``Mesh.AddNamedPointIdentification``; without it
the tests skip visibly. No Cubit is imported.
"""
from __future__ import annotations

import pytest

from radia_mcp.cubit import periodic

FIELDS = {"periodic_name": "kt", "periodic_master": "source_patch",
          "periodic_slave": "target_patch", "periodic_unit": "m",
          "periodic_patches": "4:1,2:3", "periodic_tx": -2.0, "periodic_ty": 0.0,
          "periodic_tz": 0.0, "periodic_match": 0.125}


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
