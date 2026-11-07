"""Static solver input (.jcf) export: writer and named-settings contracts, without Cubit.

Guards: the writer stays standalone (stdlib + NumPy), magnets off +z get their own coordinate
system, flow-in faces from a sideset go to the coil element the current leaves, pyramid faces map to
local face numbers, and the settings resolver fails loudly on unknown keys, unmaterialled
blocks and empty sidesets.
Assumes: the local face table and the magnet coordinate-system form, both checked against the
solver outside this repository; this file does not validate numerics.
"""
from __future__ import annotations

import ast
import importlib
import sys
import types
import zipfile
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "cubit_mesh_export"


@pytest.fixture(scope="module")
def mods():
    # Load the two modules from source without running the package __init__ (native runtime).
    pkg = types.ModuleType("_cme_jcf")
    pkg.__path__ = [str(SRC)]
    sys.modules["_cme_jcf"] = pkg
    try:
        yield importlib.import_module("_cme_jcf.jcf_writer"), importlib.import_module("_cme_jcf.jmag")
    finally:
        for name in [n for n in sys.modules if n == "_cme_jcf" or n.startswith("_cme_jcf.")]:
            del sys.modules[name]


def _stack(jw, top=None):
    """Element 1: hex z 0..1 (part 1). Element 2: hex z 1..2 (part 2) or, with top='pyramid',
    a pyramid on element 1's top face (apex node 9)."""
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    nodes = {i + 1: (x, y, z) for i, (x, y, z) in enumerate([(x, y, z) for z in (0, 1, 2) for (x, y) in sq])}
    els = [(1, jw.HEX, 1, [1, 2, 3, 4, 5, 6, 7, 8])]
    if top == "pyramid":
        nodes = {k: v for k, v in nodes.items() if k <= 8}
        nodes[9] = (0.5, 0.5, 1.5)
        els.append((2, jw.PYRAMID, 2, [5, 6, 7, 8, 9]))
    else:
        els.append((2, jw.HEX, 2, [5, 6, 7, 8, 9, 10, 11, 12]))
    return jw.mesh_from_data(nodes, els, "m")


def _root(path):
    return zipfile.ZipFile(path).read("root.xml").decode()


def test_writer_is_standalone():
    for name in ("jcf_writer.py", "jmag.py"):
        tree = ast.parse((SRC / name).read_text(encoding="utf-8"))
        imported = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        imported |= {n.module.split(".")[0] for n in ast.walk(tree)
                     if isinstance(n, ast.ImportFrom) and n.module and n.level == 0}
        assert imported <= {"__future__", "collections", "os", "struct", "tempfile", "zipfile",
                             "dataclasses", "pathlib", "xml", "numpy", "cubit"}, (name, imported)


NEAR_MINUS_Z = [(1e-6, 0, -1), (2e-6, 0, -1), (3e-7, -4e-7, -1), (1e-9, 1e-9, -1), (0.3, 0.2, -1)]
FRAME_DIRS = [(1, 0, 0), (0, 1, 0), (0, 0, -1), (1, 1, 1), (-1, 0, 0.2), *NEAR_MINUS_Z,
              (1e308, 0, 0), (1e308, 1e308, 0), (1e-320, 0, 0), (0, 0, 1e-300)]


@pytest.mark.parametrize("d", FRAME_DIRS)
def test_magnet_frame_has_direction_as_z(mods, d):
    jw, _ = mods
    x, y = jw.magnet_frame(d)
    z = np.cross(x, y)
    u = np.asarray(d, float) / np.abs(np.asarray(d, float)).max()
    assert np.allclose([x @ x, y @ y, x @ y], [1, 1, 0], rtol=0, atol=1e-14)
    assert np.allclose(z, u / np.linalg.norm(u), rtol=0, atol=1e-14), (x, y, z)


def test_unit_vectors_fail_loudly(mods):
    jw, _ = mods
    for bad in [(0, 0, 0), (np.inf, 0, 0), (np.nan, 0, 1), (1, 0)]:
        with pytest.raises(jw.JcfError, match="finite nonzero"):
            jw.magnet_frame(bad)


def test_plus_x_frame_matches_the_reference_form(mods):
    jw, _ = mods
    x, y = jw.magnet_frame((1, 0, 0))
    assert x.tolist() == [0.0, 0.0, -1.0] and y.tolist() == [0.0, 1.0, 0.0]


def test_magnet_off_z_gets_its_own_coordinate_system(mods, tmp_path):
    jw, _ = mods
    m = _stack(jw)
    mats = [{"part": 1, "name": "pm", "magnet": {"br": 1.4, "hcb": 1061000}}, {"part": 2, "name": "air", "mu_r": 1}]
    jw.write_jcf(m, tmp_path / "z.jcf", {"unit": "m", "materials": mats, "magnets": [{"part": 1, "direction": [0, 0, 1]}]})
    jw.write_jcf(m, tmp_path / "x.jcf", {"unit": "m", "materials": mats, "magnets": [{"part": 1, "direction": [1, 0, 0]}]})
    z, x = _root(tmp_path / "z.jcf"), _root(tmp_path / "x.jcf")
    assert "local_coordinate_id" not in z and z.count("<coordinate_system ") == 1
    assert '<parameter name="local_coordinate_id">2</parameter>' in x
    assert '<coordinate_system index="2" type="1">' in x
    # The solver rejects an unnormalized axis2 (error 105030 for 1e308); it must be written unit.
    jw.write_jcf(m, tmp_path / "big.jcf", {"unit": "m", "materials": mats, "magnets": [{"part": 1, "direction": [1e308, 0, 0]}]})
    big = _root(tmp_path / "big.jcf")
    i = big.index('<condition name="PermanentMagnet">')
    assert "<x>1</x>\n\t\t\t\t<y>0</y>\n\t\t\t\t<z>0</z>" in big[i:big.index("</condition>", i)]   # axis2 normalized


def test_flowin_by_sideset_nodes_equals_flowin_by_part(mods, tmp_path):
    jw, _ = mods
    m = _stack(jw)
    cu = {"mu_r": 1, "resistivity": 1.673e-8}
    mats = [dict(cu, part=1, name="a"), dict(cu, part=2, name="b")]
    by_part = jw._flowin_faces_by_part(m, {1, 2}, 2, [0, 0, 1])
    by_nodes = jw._flowin_faces_by_nodes(m, {1, 2}, [(5, 6, 7, 8)], [0, 0, 1], "cut")
    assert by_part.tolist() == by_nodes.tolist() == [[1, 2]]        # top face of element 1
    assert jw._flowin_faces_by_nodes(m, {1, 2}, [(5, 6, 7, 8)], [0, 0, -1], "cut").tolist() == [[2, 1]]
    with pytest.raises(jw.JcfError, match="coil element faces"):
        jw._flowin_faces_by_nodes(m, {1, 2}, [(1, 2, 3, 4)], [0, 0, 1], "bottom")
    r = jw.write_jcf(m, tmp_path / "c.jcf", {"unit": "m", "materials": mats, "coils": [
        {"parts": [1, 2], "current": 3, "flowin": {"face_nodes": [(5, 6, 7, 8)], "normal": [0, 0, 1], "label": "cut"}}]})
    assert ("face", 1) in r["groups"] and '<parameter name="flowin_type">0</parameter>' in _root(tmp_path / "c.jcf")


def test_pyramid_faces_map_to_local_face_numbers(mods):
    jw, _ = mods
    m = _stack(jw, "pyramid")
    lateral = jw._exterior_faces(m, [(5, 6, 9)], "side")
    assert lateral.tolist() == [[2, 3]]                       # pyramid local face 3 = nodes 0, 1, 4
    assert jw._exterior_faces(m, [(1, 2, 3, 4)], "bottom").tolist() == [[1, 1]]
    with pytest.raises(jw.JcfError, match="not an exterior"):
        jw._exterior_faces(m, [(5, 6, 7, 8)], "inner")


class FakeCubit:
    """Just enough of the Cubit Python API for name resolution."""

    def __init__(self, blocks, sidesets):
        self.blocks, self.sidesets = blocks, sidesets       # {id: name}, {id: (name, [tri node tuples])}

    def get_block_id_list(self):
        return list(self.blocks)

    def get_sideset_id_list(self):
        return list(self.sidesets)

    def get_exodus_entity_name(self, kind, i):
        return self.blocks[i] if kind == "block" else self.sidesets[i][0]

    def get_sideset_curves(self, sid):
        return []

    get_sideset_edges = get_sideset_curves
    get_sideset_surfaces = get_sideset_curves
    get_sideset_quads = get_sideset_curves

    def get_sideset_tris(self, sid):
        return list(range(len(self.sidesets[sid][1])))

    def get_expanded_connectivity(self, kind, i):
        return next(s[1][i] for s in self.sidesets.values() if i < len(s[1]))


LIB = {"pm": {"magnet": {"br": 1.4, "hcb": 1061000, "direction": [1, 0, 0]}}, "air": {"mu_r": 1}}


def test_settings_resolve_names(mods):
    _, jm = mods
    cb = FakeCubit({1: "magnet", 2: "air"}, {7: ("sym", [(5, 6, 9)])})
    spec = jm.settings_to_spec(cb, {"material_library": LIB, "block_materials": {"magnet": "pm", "air": "air"},
                                    "boundaries": {"sym": "symmetry"}})
    assert spec["magnets"] == [{"part": 1, "direction": [1.0, 0.0, 0.0]}]
    assert spec["boundaries"] == [{"type": "symmetry", "face_nodes": [(5, 6, 9)], "label": "sym"}]
    assert spec["materials"][0] == {"magnet": {"br": 1.4, "hcb": 1061000}, "part": 1, "name": "magnet"}


@pytest.mark.parametrize("settings, match", [
    ({"symmetry_sidesets": ["sym"]}, "unknown analysis settings"),
    ({"block_materials": {"magnet": "pm"}}, r"blocks \[2\] have no material"),
    ({"boundaries": {"empty": "symmetry"}}, "free mesh"),
    ({"boundaries": {"nope": "symmetry"}}, "expected one Cubit sideset"),
    ({"coils": [{"blocks": ["magnet"], "current": 1, "flowin": {"sideset": "sym"}}]}, "exactly 'sideset' and 'normal'"),
])
def test_settings_fail_loudly(mods, settings, match):
    _, jm = mods
    cb = FakeCubit({1: "magnet", 2: "air"}, {7: ("sym", [(5, 6, 9)]), 8: ("empty", [])})
    base = {"material_library": LIB, "block_materials": {"magnet": "pm", "air": "air"}}
    with pytest.raises(ValueError, match=match):
        jm.settings_to_spec(cb, {**base, **settings})
