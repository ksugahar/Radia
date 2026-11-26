"""Validation: curved .vol export through the native Cubit plugin, end to end.

Guards: ``cubit_export_vol`` (isolated worker + ``cubit_mesh_export.dll``)
on a real Cubit unit cylinder (radius 1, height 1, exact volume pi) with a
named vertex nodeset: the curved volume error must fall as the geometric order
rises from 1 to 2 and 3, the nodeset must arrive as CD3 name 1, and the result
carries its companion JSON and sha256.  Two more cases add the Kelvin
open-boundary sphere, to a full air + coil model with a z mesh seam and to a
quarter model reduced on x and z: ``kelvin_ext`` must be an exact translation
of ``kelvin_int`` (same node and triangle counts, a one-to-one point and
triangle map, every point within 1e-12 R of its translated position).

Assumes: a licensed Cubit (``RADIA_CUBIT_BIN``), the plugin directory
(``RADIA_CUBIT_PLUGIN_DIR``, built from ``packages/radia-mcp/native/cubit_vol_export``)
and a worker Python with the fork NGSolve the plugin was built against
(``RADIA_CUBIT_WORKER_PYTHON``).  Without them the tests skip visibly.  No
acceptance threshold is introduced: the measured errors are reported.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import pytest

from radia_mcp.cubit import server

CYLINDER = """reset
create cylinder height 1 radius 1
volume 1 scheme tetmesh
volume 1 size 0.35
mesh volume 1
block 1 volume 1
block 1 name "body"
nodeset 1 vertex all
nodeset 1 name "corners"
"""

AIR_COIL = """reset
create sphere radius 0.1
create cylinder height 0.02 radius 0.02
subtract volume 2 from volume 1 keep_tool
webcut volume 3 with plane zplane
imprint all
merge all
volume all size 0.015
volume all scheme tetmesh
mesh volume all
block 1 volume 2
block 1 name "coil"
block 2 volume all except 2
block 2 name "air"
"""


AIR_COIL_QUARTER = """reset
create sphere radius 0.1
create cylinder height 0.02 radius 0.02
subtract volume 2 from volume 1 keep_tool
webcut volume all with plane zplane
delete volume with z_coord < 0
webcut volume all with plane xplane
delete volume with x_coord < 0
imprint volume all
merge volume all
volume all size 0.015
volume all scheme tetmesh
mesh volume all
block 1 add volume with z_coord < 0.0101 and x_coord < 0.0101
block 1 name "coil"
block 2 add volume all except volume in block 1
block 2 name "air"
"""


def _runtime():
    names = ("RADIA_CUBIT_BIN", "RADIA_CUBIT_PLUGIN_DIR", "RADIA_CUBIT_WORKER_PYTHON")
    missing = [name for name in names if not os.environ.get(name)]
    if missing:
        pytest.skip(f"native Cubit export runtime not configured: {missing}")
    return {name: os.environ[name] for name in names}


def _export(tmp_path, journal_text, name, order, **kwargs):
    rt = _runtime()
    journal = tmp_path / f"{name}.jou"
    journal.write_text(journal_text, encoding="utf-8")
    vol = tmp_path / f"{name}_p{order}.vol"
    result = json.loads(server.cubit_export_vol(str(journal), str(vol), order,
                                                cubit_bin=rt["RADIA_CUBIT_BIN"],
                                                worker_python=rt["RADIA_CUBIT_WORKER_PYTHON"],
                                                plugin_dir=rt["RADIA_CUBIT_PLUGIN_DIR"], **kwargs))
    assert result["status"] == "ok", result
    assert result["sha256"] == hashlib.sha256(vol.read_bytes()).hexdigest()
    assert Path(result["companion_json"]).is_file(), result
    return vol, result


def test_cylinder_volume_error_falls_with_order_and_nodeset_survives(tmp_path):
    errors = {}
    for order in (1, 2, 3):
        vol, result = _export(tmp_path, CYLINDER, "cyl", order)
        text = vol.read_text(encoding="utf-8")
        assert "cd3names\n1\n1\tcorners\n" in text
        import ngsolve as ng
        volume = ng.Integrate(ng.CF(1), ng.Mesh(str(vol)))
        errors[order] = abs(volume - math.pi) / math.pi
    print("relative volume error by order:", errors)
    assert errors[1] > errors[2] > errors[3], errors


def _surface(mesh, label):
    from netgen.meshing import PointId
    names = mesh.GetBoundaries()
    points = mesh.ngmesh.Points()
    tris = [tuple(int(v.nr) for v in el.vertices[:3])
            for el in mesh.ngmesh.Elements2D() if names[el.index - 1] == label]
    nodes = sorted({n for t in tris for n in t})
    return tris, nodes, {n: tuple(points[PointId(n)].p) for n in nodes}


@pytest.mark.parametrize("name,journal,options", [
    ("seam", AIR_COIL, {}),
    ("quarter", AIR_COIL_QUARTER, {"kelvin_sym_x": "ht", "kelvin_sym_z": "bn"}),
])
def test_kelvin_exterior_is_exact_translation_of_interior(tmp_path, name, journal, options):
    vol, result = _export(tmp_path, journal, name, 2, add_kelvin=True, **options)
    import ngsolve as ng
    mesh = ng.Mesh(str(vol))
    offset, R = result["kelvin"]["center"], result["kelvin"]["R"]
    tri_i, nodes_i, x_i = _surface(mesh, "kelvin_int")
    tri_e, nodes_e, x_e = _surface(mesh, "kelvin_ext")
    assert (len(nodes_i), len(tri_i)) == (len(nodes_e), len(tri_e)), (
        len(nodes_i), len(tri_i), len(nodes_e), len(tri_e))
    mapping, worst = {}, 0.0
    for n in nodes_i:
        q = [a + b for a, b in zip(x_i[n], offset)]
        dist, m = min((max(abs(a - b) for a, b in zip(x_e[k], q)), k) for k in nodes_e)
        mapping[n] = m
        worst = max(worst, dist)
    print(name, "kelvin nodes", len(nodes_i), "triangles", len(tri_i), "max deviation / R", worst / R)
    assert worst <= 1e-12 * R, (worst, 1e-12 * R)
    assert len(set(mapping.values())) == len(nodes_e)
    ext = {frozenset(t) for t in tri_e}
    unmapped = [t for t in tri_i if frozenset(mapping[n] for n in t) not in ext]
    assert not unmapped, unmapped[:5]
