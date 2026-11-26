"""Validation: the external-Python Cubit curving worker reproduces exact geometry.

Guards: ``radia_mcp.cubit.vol_export`` run end to end through its CLI on a
real Cubit unit cylinder (radius 1, height 1, exact volume pi); the curved
volume error must fall as the geometric order rises from 1 to 2 and 3.

Assumes: a licensed Cubit whose bin directory is in ``RADIA_CUBIT_BIN`` and a
Netgen/NGSolve build with the callback-geometry bindings.  Without them the
test skips visibly.  No acceptance threshold is introduced here: the check
is the order-monotone error decrease, and the measured errors are reported.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

JOURNAL = """reset
create cylinder height 1 radius 1
volume 1 scheme tetmesh
volume 1 size 0.35
mesh volume 1
block 1 volume 1
block 1 name "body"
"""


def _cubit_bin():
    path = os.environ.get("RADIA_CUBIT_BIN")
    if not path or not Path(path).is_dir():
        pytest.skip("RADIA_CUBIT_BIN does not name a Cubit bin directory")
    return path


def _volume(vol_path, order):
    from ngsolve import CF, Integrate, Mesh
    mesh = Mesh(str(vol_path))  # curved elements come from the file's curvedelements section
    return Integrate(CF(1.0), mesh, order=2 * order + 2)


def test_curved_cylinder_volume_error_falls_with_order(tmp_path):
    cubit_bin = _cubit_bin()
    pytest.importorskip("ngsolve")
    journal = tmp_path / "cyl.jou"
    journal.write_text(JOURNAL, encoding="utf-8")
    errors = {}
    for order in (1, 2, 3):
        out = tmp_path / f"cyl_p{order}.vol"
        status = tmp_path / f"cyl_p{order}.json"
        run = subprocess.run(
            [sys.executable, "-B", "-m", "radia_mcp.cubit.vol_export",
             "--cubit-bin", cubit_bin, "--journal", str(journal), "--out", str(out),
             "--order", str(order), "--result", str(status)],
            capture_output=True, cwd=tmp_path, timeout=1800)
        assert run.returncode == 0, (order, run.returncode, run.stderr[-2000:].decode(errors="replace"))
        assert json.loads(status.read_text(encoding="utf-8"))["status"] == "ok"
        errors[order] = abs(_volume(out, order) - math.pi) / math.pi
    print("relative volume error by order:", errors)
    assert errors[1] > errors[2] > errors[3], errors
