"""Contract for the Kelvin point identifications: reading them across Netgen
6.2.2606/2607, and using the .vol pairs or pairing a flat mesh after load.

6.2.2607 appends the identification number to each pair.  The reader must accept
both shapes, and a changed shape must raise instead of reading as "no Kelvin
identification" (an earlier catch-all hid exactly that change).
"""
from types import SimpleNamespace

import pytest

from radia.kelvin_identify_ngsolve import _identified_point_pairs, has_kelvin_identification


def _pid(nr):
    return SimpleNamespace(nr=nr)


def test_pairs_and_triples_read_the_same_point_numbers():
    pairs = [(_pid(3), _pid(7)), (_pid(4), _pid(8))]
    triples = [(_pid(3), _pid(7), 1), (_pid(4), _pid(8), 1)]
    assert _identified_point_pairs(pairs) == [(3, 7), (4, 8)]
    assert _identified_point_pairs(triples) == [(3, 7), (4, 8)]


@pytest.mark.parametrize("entry", [(_pid(1),), (_pid(1), _pid(2), 1, "x")])
def test_other_entry_shapes_raise(entry):
    with pytest.raises(ValueError, match="identification entry"):
        _identified_point_pairs([entry])


def test_mesh_errors_propagate_instead_of_reading_as_unidentified():
    class BrokenMesh:
        def GetBoundaries(self):
            raise RuntimeError("mesh access failed")

    with pytest.raises(RuntimeError, match="mesh access failed"):
        has_kelvin_identification(BrokenMesh())


def _kelvin_box_mesh(tmp_path, strip):
    """Air box + translated Kelvin box, meshed with the pairs, optionally stripped."""
    occ = pytest.importorskip("netgen.occ")
    ngsolve = pytest.importorskip("ngsolve")
    air = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    kel = occ.Box(occ.Pnt(3, 0, 0), occ.Pnt(4, 1, 1))
    air.mat("air")
    kel.mat("kelvin")
    air.faces.name = "kelvin_int"
    kel.faces.name = "kelvin_ext"
    trafo = occ.gp_Trsf.Translation(occ.Vec(3, 0, 0))
    for a, b in zip(air.faces, kel.faces):
        a.Identify(b, "kelvin", occ.IdentificationType.PERIODIC, trafo)
    path = tmp_path / "k.vol"
    occ.OCCGeometry(occ.Glue([air, kel])).GenerateMesh(maxh=0.5).Save(str(path))
    if strip:
        lines, kept, i = path.read_text().splitlines(), [], 0
        while i < len(lines):
            if lines[i].strip() in ("identifications", "identificationtypes",
                                    "identificationnames"):
                i += 2 + int(lines[i + 1].split()[0])
                continue
            kept.append(lines[i])
            i += 1
        path.write_text("\n".join(kept) + "\n")
    return ngsolve.Mesh(str(path))


def _pairs(mesh):
    return {frozenset(p.nr for p in row[:2]) for row in mesh.ngmesh.GetIdentifications()}


@pytest.mark.usefixtures("ngsolve_taskmanager")  # caller-owned TaskManager
def test_ensure_uses_vol_pairs_or_pairs_a_flat_mesh_after_load(tmp_path):
    from radia.kelvin_identify_ngsolve import ensure_kelvin_identification
    (tmp_path / "a").mkdir()
    paired = _kelvin_box_mesh(tmp_path / "a", strip=False)
    expected = _pairs(paired)
    assert has_kelvin_identification(paired) and expected
    assert ensure_kelvin_identification(paired)["source"] == "vol"
    assert _pairs(paired) == expected

    (tmp_path / "b").mkdir()
    bare = _kelvin_box_mesh(tmp_path / "b", strip=True)
    assert not tuple(bare.ngmesh.GetIdentifications())
    info = ensure_kelvin_identification(bare)
    assert info["source"] == "post_hoc" and info["n_unmatched"] == 0, info
    assert max(abs(c - e) for c, e in zip(info["kelvin_offset"], (3, 0, 0))) < 1e-12
    assert _pairs(bare) == expected


@pytest.mark.usefixtures("ngsolve_taskmanager")  # caller-owned TaskManager
def test_ensure_refuses_partial_and_curved_meshes(tmp_path):
    from radia.kelvin_identify_ngsolve import ensure_kelvin_identification
    (tmp_path / "p").mkdir()
    partial = _kelvin_box_mesh(tmp_path / "p", strip=True)
    partial.ngmesh.AddPointIdentification(1, 2, identnr=1, type=2)
    with pytest.raises(ValueError, match="do not pair every"):
        ensure_kelvin_identification(partial)
    (tmp_path / "c").mkdir()
    curved = _kelvin_box_mesh(tmp_path / "c", strip=True)
    curved.Curve(2)
    with pytest.raises(ValueError, match="curve order 2"):
        ensure_kelvin_identification(curved)
    assert not tuple(curved.ngmesh.GetIdentifications())


def _true_pairs(tmp_path, name, inner):
    """The exported (inner, outer) point pairs, oriented by ``inner``."""
    (tmp_path / name).mkdir()
    paired = _kelvin_box_mesh(tmp_path / name, strip=False)
    rows = [(r[0].nr, r[1].nr) for r in paired.ngmesh.GetIdentifications()]
    return sorted((a, b) if a in inner else (b, a) for a, b in rows)


@pytest.mark.usefixtures("ngsolve_taskmanager")  # caller-owned TaskManager
@pytest.mark.parametrize("split_groups", [True, False])
def test_ensure_refuses_a_covering_but_wrong_kelvin_mapping(tmp_path, split_groups):
    from radia.kelvin_identify_ngsolve import ensure_kelvin_identification
    from radia.kelvin_identify_ngsolve import _boundary_vertex_ids
    (tmp_path / "w").mkdir()
    wrong = _kelvin_box_mesh(tmp_path / "w", strip=True)
    pairs = _true_pairs(tmp_path, "t", _boundary_vertex_ids(wrong, "kelvin_int"))
    if split_groups:   # true partners, but split over identification numbers 1 and 2
        rows = [(a, b, 1 + (i >= len(pairs) // 2)) for i, (a, b) in enumerate(pairs)]
    else:              # one group, every vertex covered, partners rotated
        outers = [b for _, b in pairs]
        rows = [(a, b, 1) for (a, _), b in zip(pairs, outers[1:] + outers[:1])]
    for a, b, nr in rows:
        wrong.ngmesh.AddPointIdentification(a, b, identnr=nr, type=2)
    assert has_kelvin_identification(wrong)  # coverage alone would accept it
    with pytest.raises(ValueError, match="identification numbers" if split_groups else "one translation"):
        ensure_kelvin_identification(wrong)


@pytest.mark.usefixtures("ngsolve_taskmanager")  # caller-owned TaskManager
def test_incomplete_pairing_registers_nothing(tmp_path):
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    (tmp_path / "s").mkdir()
    bare = _kelvin_box_mesh(tmp_path / "s", strip=True)
    with pytest.raises(ValueError, match="no pair was registered"):
        add_kelvin_identification(bare, kelvin_offset=(3.0, 0.0, 0.4), point_tolerance=1e-3,
                                  skip_if_existing=False, require_complete=True)
    assert not tuple(bare.ngmesh.GetIdentifications())
