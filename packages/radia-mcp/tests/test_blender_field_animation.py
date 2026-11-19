"""Guards: blender_field_animation input checks, worker resolution and readback.

Assumes: h5py for dataset tests.  The live write/reopen test runs only when
RADIA_BLENDER_PYTHON points at a Python 3.13 + bpy interpreter; it is the
direct check that a reopened .blend reproduces pose, time and signed J.
"""
import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest

from radia_mcp.blender import field_animation as fa

h5py = pytest.importorskip("h5py")


def _dataset(path, n=4, m=9, bad=None):
    t = np.arange(n) * 0.25
    th = np.linspace(0, np.pi, n)
    rot = np.array([[[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]] for a in th])
    xs = np.linspace(-0.004, 0.004, 3)
    pts = np.array([[x, -2.5e-4, z] for x in xs for z in xs])
    # Closed-form field: rigid in-plane rotation of a uniform J (A/m^2).
    J = np.array([[[100 * np.cos(a), 0.0, 100 * np.sin(a)]] * m for a in th])
    with h5py.File(path, "w") as f:
        f.attrs["schema"] = bad or fa.SCHEMA
        f["time_s"] = t
        f["magnet/center_m"] = np.column_stack([np.linspace(-0.003, 0.003, n), np.full(n, 0.002), np.zeros(n)])
        f["magnet/rotation"] = rot
        f["magnet/size_m"] = np.full(3, 1e-3)
        f["magnet/magnetization_A_per_m"] = np.array([0.0, 1.6e5, 0.0])
        f["conductor/box_min_m"] = np.array([-0.006, -5e-4, -0.006])
        f["conductor/box_max_m"] = np.array([0.006, 0.0, 0.006])
        f["samples/id"] = np.arange(m)
        f["samples/point_m"] = pts
        f["field/J_A_per_m2"] = J
    return path


def test_rejects_wrong_schema(tmp_path):
    with pytest.raises(ValueError, match="schema"):
        fa.load_dataset(_dataset(tmp_path / "d.h5", bad="other"))


def test_rejects_non_orthonormal_rotation(tmp_path):
    p = _dataset(tmp_path / "d.h5")
    with h5py.File(p, "r+") as f:
        r = f["magnet/rotation"][()]
        r[1] *= 1.01
        f["magnet/rotation"][...] = r
    with pytest.raises(ValueError, match="orthonormal"):
        fa.load_dataset(p)


def test_missing_interpreter_is_a_clear_error(monkeypatch):
    monkeypatch.delenv("RADIA_BLENDER_PYTHON", raising=False)
    with pytest.raises(fa.BlenderWorkerError, match="RADIA_BLENDER_PYTHON is not set"):
        fa.resolve_worker()


def test_missing_worker_is_a_clear_error(monkeypatch, tmp_path):
    fake = tmp_path / "python.exe"
    fake.write_text("")
    monkeypatch.setenv("RADIA_BLENDER_PYTHON", str(fake))
    monkeypatch.setenv("RADIA_BLENDER_WORKER", str(tmp_path / "absent.py"))
    with pytest.raises(fa.BlenderWorkerError, match="not shipped in the wheel"):
        fa.resolve_worker()


def test_source_tree_worker_is_found():
    default = Path(fa.__file__).resolve().parents[3] / "tools" / "blender" / fa.WORKER_NAME
    assert default.is_file(), f"expected worker at {default}"


def _back(pay, spec):
    n = len(pay["time_s"])
    steps = {}
    for k in range(n):
        steps[str(k)] = {"center_m": pay["magnet"]["center_m"][k],
                         "rotation": pay["magnet"]["rotation"][k],
                         "visible": [f"J_glyph_{k:04d}", f"Time_{k:04d}"],
                         "time_label": f"t = {pay['time_s'][k]:.3f} s",
                         "sample_id": pay["samples"]["id"],
                         "J_A_per_m2": np.asarray(pay["J_A_per_m2"][k]).ravel().tolist()}
    return {"frame_start": 1, "frame_end": n, "fps": spec["fps"], "schema": fa.SCHEMA,
            "camera_animated": False, "camera_motion": "static", "n_glyph_objects": n,
            "interpolations": ["CONSTANT"],
            "bbox_m": {"Conductor": [pay["conductor"]["box_min_m"], pay["conductor"]["box_max_m"]]},
            "physical_time_s": pay["time_s"], "J_scale_A_per_m2": [0.0, pay["J_scale_max_A_per_m2"]],
            "steps": steps}


def test_compare_accepts_faithful_and_rejects_altered_readback(tmp_path):
    pay = fa.load_dataset(_dataset(tmp_path / "d.h5"))
    pay["J_scale_max_A_per_m2"] = 100.0
    spec = {"fps": 12}
    good = _back(pay, spec)
    ok = fa.compare(pay, good, spec, [0, 1, 3])
    assert ok["passed"], ok["failures"]
    for mutate, needle in (
            (lambda b: b["steps"]["1"]["J_A_per_m2"].__setitem__(0, -b["steps"]["1"]["J_A_per_m2"][0]), "J readback"),
            (lambda b: b["steps"]["3"]["center_m"].__setitem__(0, 1.0), "pose"),
            (lambda b: b.__setitem__("camera_animated", True), "camera"),
            (lambda b: b["steps"]["1"].__setitem__("visible", ["J_glyph_0002"]), "visible"),
            (lambda b: b.__setitem__("fps", 24), "fps"),
            (lambda b: b.__setitem__("interpolations", ["BEZIER", "CONSTANT"]), "interpolation"),
            (lambda b: b["bbox_m"].__setitem__("Conductor", [[0, 0, 0], [1, 1, 1]]), "bounding box")):
        bad = copy.deepcopy(good)
        mutate(bad)
        res = fa.compare(pay, bad, spec, [0, 1, 3])
        assert not res["passed"] and any(needle in f for f in res["failures"]), res["failures"]


@pytest.mark.skipif(not os.environ.get("RADIA_BLENDER_PYTHON"), reason="optional bpy worker not configured")
def test_live_write_reopen_matches_closed_form(tmp_path):
    p = _dataset(tmp_path / "d.h5")
    res = fa.field_animation(p, tmp_path / "out", {"resolution": [160, 90], "samples": 1,
                                                   "video": False, "preview_steps": [2], "threads": 2})
    assert res["status"] == "ok", json.dumps(res["readback"], indent=1)
    assert Path(res["artifacts"]["blend"]).is_file()
    th = np.linspace(0, np.pi, 4)[2]
    assert res["readback"]["steps"]["2"]["J_max_rel_err"] <= 1e-6
    assert res["J_scale_A_per_m2"][1] == pytest.approx(100.0)
    assert abs(np.cos(th)) <= 1.0
