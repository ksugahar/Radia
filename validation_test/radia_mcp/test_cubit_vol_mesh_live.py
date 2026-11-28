"""Validation: the in-memory Cubit mesh equals the export vol file at every order.

Guards: ``radia_mcp.cubit.vol_mesh`` (``_cubit_to_vol_mesh``).  On a real Cubit
unit cylinder and on the Kelvin quarter air + coil model, the mesh built in
memory at order 2 and curved again at orders 3, 1 and 2 maps every
tetrahedron's points where the ``export vol`` file of the same order, written
in the same Cubit session, maps them (difference <= 1e-12 of the coordinate
scale), with the same counts, labels and Kelvin periodic pairs (the same
vertex pairs and groups).  Curving back to order 2 reproduces the first
order-2 mesh, the cylinder volume error falls from order 1 to 3, and curving
again is refused, before the mesh is touched, after a move, a flipped surface
normal, a reversed curve, a renamed block, sideset or nodeset, or a new volume.

The full-sphere Kelvin model has a same-material webcut seam at z = 0.  The
export curves with the seam triangles and then removes them, so the file's
seam faces carry the projection noise of a planar face; a re-curving has no
seam triangles left.  There the 1e-12 bound holds for the elements off the
seam, and the seam elements' difference is reported, not bounded.

Assumes: a licensed Cubit (``RADIA_CUBIT_BIN``), the plugin directory with
``cubit_to_vol_export.dll`` and ``_cubit_to_vol_mesh`` (``RADIA_CUBIT_PLUGIN_DIR``)
and the worker Python they were built for (``RADIA_CUBIT_WORKER_PYTHON``).
Without them the tests skip visibly.  The export vol file is the reference;
its own accuracy is guarded by test_cubit_vol_export_curving.py.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import pytest

from test_cubit_vol_export_curving import AIR_COIL, AIR_COIL_QUARTER, CYLINDER

HERE = Path(__file__).resolve().parent
SOURCE_ROOT = HERE.parents[1] / "packages" / "radia-mcp" / "src"


def _session(tmp_path, case, journal):
    names = ("RADIA_CUBIT_BIN", "RADIA_CUBIT_PLUGIN_DIR", "RADIA_CUBIT_WORKER_PYTHON")
    missing = [name for name in names if not os.environ.get(name)]
    if missing:
        pytest.skip(f"native Cubit export runtime not configured: {missing}")
    (tmp_path / f"{case}.jou").write_text(journal, encoding="utf-8")
    out = tmp_path / f"{case}.json"
    run = subprocess.run(
        [os.environ["RADIA_CUBIT_WORKER_PYTHON"], "-I", "-B", str(HERE / "cubit_vol_mesh_session.py"),
         str(out), os.environ["RADIA_CUBIT_PLUGIN_DIR"], os.environ["RADIA_CUBIT_BIN"],
         str(SOURCE_ROOT), str(tmp_path), case],
        capture_output=True, encoding="utf-8", errors="replace", timeout=600)
    assert out.is_file(), run.stdout[-2000:] + run.stderr[-2000:]
    result = json.loads(out.read_text(encoding="utf-8"))
    assert result["status"] == "ok", result.get("traceback")
    return result


def _report(result):
    print(result["case"], json.dumps([{k: row.get(k) for k in ("order", "max_off_z0", "max_on_z0", "volume")}
                                      for row in result["orders"]]),
          "back to 2:", result["recurve_back_to_2_off_z0"], result["recurve_back_to_2_on_z0"])


@pytest.mark.parametrize("case,journal", [("cylinder", CYLINDER), ("kelvin_quarter", AIR_COIL_QUARTER)])
def test_live_mesh_matches_export_vol_at_every_order(tmp_path, case, journal):
    result = _session(tmp_path, case, journal)
    _report(result)
    tol = 1e-12 * result["coordinate_scale"]
    assert result["is_live"]
    for row in result["orders"]:
        assert row["same_shape"], row
        assert row["live"] == row["file"], row
        assert row["max_point_difference"] <= tol, (row["order"], row["max_point_difference"], tol)
    assert max(result["recurve_back_to_2_off_z0"], result["recurve_back_to_2_on_z0"]) <= tol, result
    if case == "kelvin_quarter":
        assert all(row["live"]["periodic_pairs"] > 0 for row in result["orders"]), result["orders"]
        assert "sideset_name" in result["edits"], result["edits"]
    else:
        assert "nodeset_name" in result["edits"], result["edits"]
        error = {row["order"]: abs(row["volume"] - math.pi) / math.pi for row in result["orders"]}
        assert error[1] > error[2] > error[3], error
    assert all(edit["refused"] for edit in result["edits"].values()), result["edits"]
    assert {"move", "surface_normal", "curve_sense", "block_name", "new_volume"} <= set(result["edits"])
    assert result["changed_model_refused"], result.get("changed_model_message")
    assert result["bad_order"] == {**result["bad_order"], "refused": True, "curved_after": 3}, result["bad_order"]


def test_seam_model_matches_export_vol_off_the_removed_seam(tmp_path):
    result = _session(tmp_path, "kelvin_seam", AIR_COIL)
    _report(result)
    tol = 1e-12 * result["coordinate_scale"]
    for row in result["orders"]:
        assert row["same_shape"], row
        assert row["live"] == row["file"], row
        assert row["live"]["periodic_pairs"] > 0, row
        assert row["max_off_z0"] <= tol, (row["order"], row["max_off_z0"], tol)
    assert result["recurve_back_to_2_off_z0"] <= tol, result
    assert all(edit["refused"] for edit in result["edits"].values()), result["edits"]
    assert {"move", "surface_normal", "curve_sense", "block_name", "new_volume"} <= set(result["edits"])
    assert result["changed_model_refused"], result.get("changed_model_message")
    assert result["bad_order"] == {**result["bad_order"], "refused": True, "curved_after": 3}, result["bad_order"]
