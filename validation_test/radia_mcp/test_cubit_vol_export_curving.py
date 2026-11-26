"""Validation: curved .vol export through the native Cubit plugin, end to end.

Guards: ``cubit_export_vol`` (isolated worker + ``cubit_mesh_export.dll``)
on a real Cubit unit cylinder (radius 1, height 1, exact volume pi) with a
named vertex nodeset: the curved volume error must fall as the geometric order
rises from 1 to 2 and 3, the nodeset must arrive as CD3 name 1, and the result
carries its companion JSON and sha256.  A second case adds the Kelvin
open-boundary sphere to an air + coil model and requires ``kelvin_int`` and
``kelvin_ext`` boundaries in the ``.vol``.

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


def test_kelvin_shell_is_added_and_labelled(tmp_path):
    vol, result = _export(tmp_path, AIR_COIL, "kelvin", 2, add_kelvin=True)
    assert result["kelvin"] is not None, result
    text = vol.read_text(encoding="utf-8")
    assert "kelvin_int" in text and "kelvin_ext" in text
