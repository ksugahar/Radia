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


def _b_dataset(path, both=False):
    """Rotating uniform B in a sphere, no magnet (field key B_T, unit T)."""
    n, th = 3, np.linspace(0, np.pi, 3)
    pts = np.array([[x, 0.0, z] for x in (-2e-3, 0.0, 2e-3) for z in (-2e-3, 0.0, 2e-3)])
    B = np.array([[[0.1 * np.cos(a), 0.0, 0.1 * np.sin(a)]] * len(pts) for a in th])
    with h5py.File(path, "w") as f:
        f.attrs["schema"] = fa.SCHEMA
        f["time_s"] = np.arange(n) * 0.5
        f["body/box_min_m"] = np.full(3, -5e-3)
        f["body/box_max_m"] = np.full(3, 5e-3)
        f["body"].attrs["shape"] = "sphere"
        f["samples/id"] = np.arange(len(pts))
        f["samples/point_m"] = pts
        f["field/B_T"] = B
        if both:
            f["field/J_A_per_m2"] = B
    return path


def test_b_field_keeps_its_name_and_unit(tmp_path):
    pay = fa.load_dataset(_b_dataset(tmp_path / "b.h5"))
    assert pay["field"] == {"key": "B_T", "symbol": "B", "unit": "T"}
    assert pay["magnet"] is None and pay["region"]["name"] == "Body" and pay["region"]["shape"] == "sphere"


def test_rejects_two_fields(tmp_path):
    with pytest.raises(ValueError, match="exactly one"):
        fa.load_dataset(_b_dataset(tmp_path / "b.h5", both=True))


OFFSET = 1e-4


def _glyphs(pay, spec, k):
    """Arrows as a faithful renderer would draw them (closed-form, not the worker)."""
    out = []
    jmax = pay["scale_max"]
    for sid, p, v in zip(pay["samples"]["id"], pay["samples"]["point_m"], pay["vectors"][k]):
        mag = float(np.linalg.norm(v))
        d = np.array([v[0], 0.0, v[2]]) / np.hypot(v[0], v[2])
        c = np.array([p[0], p[1] + OFFSET, p[2]])
        L = spec["arrow_length_m"] * min(mag / jmax, 1.0)
        out.append({"id": int(sid), "base": list(c - 0.5 * L * d), "tip": list(c + 0.5 * L * d),
                    "color": list(fa.expected_color(mag / jmax, spec.get("colormap", "viridis")))})
    return out


def _back(pay, spec):
    n = len(pay["time_s"])
    steps = {}
    for k in range(n):
        steps[str(k)] = {"glyphs": _glyphs(pay, spec, k), "center_m": pay["magnet"]["center_m"][k],
                         "rotation": pay["magnet"]["rotation"][k],
                         "visible": [f"J_glyph_{k:04d}", f"Time_{k:04d}"],
                         "time_label": f"t = {pay['time_s'][k]:.3f} s",
                         "sample_id": pay["samples"]["id"],
                         "vectors": np.asarray(pay["vectors"][k]).ravel().tolist()}
    return {"frame_start": 1, "frame_end": n, "fps": spec["fps"], "schema": fa.SCHEMA,
            "arrow_length_m": spec["arrow_length_m"], "glyph_display_offset_m": OFFSET,
            "camera_animated": False, "camera_motion": "static", "n_glyph_objects": n,
            "style": spec.get("style", "standard"), "camera": spec.get("camera", "isometric"),
            "colormap": spec.get("colormap", "viridis"),
            "interpolations": ["CONSTANT"],
            "bbox_m": {"Conductor": [pay["region"]["box_min_m"], pay["region"]["box_max_m"]]},
            "field_key": "J_A_per_m2", "field_unit": "A/m^2", "magnet_present": True,
            "physical_time_s": pay["time_s"], "field_scale": [0.0, pay["scale_max"]],
            "steps": steps}


def test_compare_accepts_faithful_and_rejects_altered_readback(tmp_path):
    pay = fa.load_dataset(_dataset(tmp_path / "d.h5"))
    pay["scale_max"] = 100.0
    spec = {"fps": 12, "arrow_length_m": 1e-3}
    good = _back(pay, spec)
    ok = fa.compare(pay, good, spec, [0, 1, 3])
    assert ok["passed"], ok["failures"]
    for mutate, needle in (
            (lambda b: b["steps"]["1"]["vectors"].__setitem__(0, -b["steps"]["1"]["vectors"][0]), "J readback"),
            (lambda b: b.__setitem__("field_unit", "T"), "field"),
            (lambda b: b["steps"]["3"]["center_m"].__setitem__(0, 1.0), "pose"),
            (lambda b: b.__setitem__("camera_animated", True), "camera"),
            (lambda b: b.__setitem__("colormap", "turbo"), "colormap"),
            (lambda b: b.__setitem__("camera", "top"), "camera in file"),
            (lambda b: b["steps"]["1"].__setitem__("visible", ["J_glyph_0002"]), "visible"),
            (lambda b: b.__setitem__("fps", 24), "fps"),
            # Visible glyphs alone altered; hidden J data untouched.
            (lambda b: b["steps"]["1"]["glyphs"][2].update(base=b["steps"]["1"]["glyphs"][2]["tip"],
                                                           tip=b["steps"]["1"]["glyphs"][2]["base"]), "direction"),
            (lambda b: b["steps"]["3"]["glyphs"][0]["base"].__setitem__(0, 0.5), "position"),
            (lambda b: b["steps"]["0"]["glyphs"][4].__setitem__("color", [1.0, 0.0, 0.0]), "color"),
            (lambda b: b["steps"]["1"]["glyphs"].pop(), "arrow ids"),
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
    assert res["readback"]["steps"]["2"]["field_max_rel_err"] <= 1e-6
    assert res["field_scale"][1] == pytest.approx(100.0)
    assert abs(np.cos(th)) <= 1.0


_FLIP = """
import sys, bpy
bpy.ops.wm.open_mainfile(filepath=sys.argv[1])
me = bpy.data.objects["J_glyph_0002"].data
for a in range(len(me.vertices) // 7):        # reverse every arrow of step 2 only
    v = [me.vertices[7 * a + i].co.copy() for i in range(7)]
    base, tip = (v[0] + v[6]) / 2, v[3]
    for i in range(7):
        me.vertices[7 * a + i].co = base + tip - v[i]
bpy.ops.wm.save_mainfile(filepath=sys.argv[1])
"""


@pytest.mark.skipif(not os.environ.get("RADIA_BLENDER_PYTHON"), reason="optional bpy worker not configured")
def test_live_rendered_glyph_edit_is_rejected(tmp_path):
    import subprocess
    p = _dataset(tmp_path / "d.h5")
    res = fa.field_animation(p, tmp_path / "out", {"resolution": [160, 90], "samples": 1,
                                                   "video": False, "preview_steps": [2], "threads": 2})
    assert res["status"] == "ok", json.dumps(res["readback"], indent=1)
    py, worker, _ = fa.resolve_worker()
    script = tmp_path / "flip.py"
    script.write_text(_FLIP, encoding="utf-8")
    blend = res["artifacts"]["blend"]
    subprocess.run([py, "-I", "-B", str(script), blend], check=True, capture_output=True,
                   stdin=subprocess.DEVNULL, timeout=300)
    back = fa._run(py, worker, "inspect", {"blend_path": blend, "steps": [0, 2]}, tmp_path, 300)
    pay = fa.load_dataset(p)
    pay["scale_max"] = res["field_scale"][1]
    v = fa.compare(pay, back, res["spec"], [0, 2])
    # Hidden J data is untouched, so only the rendered-arrow check can catch this.
    assert v["steps"]["2"]["field_max_rel_err"] <= 1e-6
    assert not v["passed"] and any("step 2: rendered arrow direction" in f for f in v["failures"]), v["failures"]
    assert not any(f.startswith("step 0") for f in v["failures"]), v["failures"]


def test_worker_rejects_unknown_spec_options_and_distinguishes_colormaps(monkeypatch):
    import importlib.util
    import sys
    import types

    for mod in ("bpy", "bmesh", "mathutils"):
        if mod not in sys.modules:
            stub = types.ModuleType(mod)
            if mod == "mathutils":
                stub.Matrix = object
                stub.Vector = object
            monkeypatch.setitem(sys.modules, mod, stub)

    worker_path = Path(__file__).resolve().parents[1] / "tools" / "blender" / "field_animation_worker.py"
    spec = importlib.util.spec_from_file_location("field_animation_worker_test", worker_path)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)

    for bad in ({"camera": "Top"}, {"colormap": "jet"}, {"style": "cinematic"}):
        with pytest.raises(ValueError):
            worker._validate_spec_options(bad)

    with pytest.raises(RuntimeError):  # no table from the request: fail, no built-in palette
        worker._cmap(0.5, "viridis")
    assert set(worker._COLORMAP_NAMES) == set(fa.COLORMAPS)
    # The worker draws exactly the table the wrapper sends, with the same lookup.
    for name, rgb in fa.COLORMAPS.items():
        worker._LUT = (name, [tuple(c) for c in rgb])
        for t in (0.0, 0.1, 0.37, 0.5, 0.8, 1.0):
            assert np.allclose(worker._cmap(t, name), fa.expected_color(t, name), atol=0), (name, t)
    assert worker._VALID_STYLES == fa.STYLES and worker._VALID_CAMERAS == fa.CAMERAS


def test_compare_checks_colours_against_the_requested_colormap(tmp_path):
    """Guards: a turbo render is judged with turbo colours, not the viridis default."""
    pay = fa.load_dataset(_dataset(tmp_path / "d.h5"))
    pay["scale_max"] = 100.0
    turbo = {"fps": 12, "arrow_length_m": 1e-3, "colormap": "turbo"}
    ok = fa.compare(pay, _back(pay, turbo), turbo, [0, 1, 3])
    assert ok["passed"], ok["failures"]
    wrong = _back(pay, {"fps": 12, "arrow_length_m": 1e-3})  # viridis colours
    wrong["colormap"] = "turbo"
    bad = fa.compare(pay, wrong, turbo, [0, 1, 3])
    assert not bad["passed"] and any("color" in f for f in bad["failures"]), bad["failures"]


@pytest.mark.parametrize("opt,value", [("style", "cinematic"), ("camera", "Top"), ("colormap", "jet")])
def test_unknown_render_options_fail_before_blender_starts(tmp_path, monkeypatch, opt, value):
    called = []
    monkeypatch.setattr(fa, "resolve_worker", lambda: called.append(1))
    with pytest.raises(ValueError, match=opt):
        fa.field_animation(_dataset(tmp_path / "d.h5"), tmp_path / "out", {opt: value})
    assert not called, "Blender worker was resolved for an invalid option"


def test_colormap_tables_match_upstream_samples():
    """Guards: viridis/turbo/inferno equal the pinned upstream 256-entry tables and
    the declared lookup (256 bins, no interpolation) at bin edges and interiors.
    Assumes: matplotlib's ListedColormap(N=256) is the independent oracle."""
    mpl = pytest.importorskip("matplotlib")
    lut = fa.COLORMAP_LUT
    assert lut["schema"] == "radia.colormap-lut.v1"
    assert {k: v["license"] for k, v in lut["maps"].items()} == {
        "viridis": "CC0-1.0", "inferno": "CC0-1.0", "turbo": "Apache-2.0"}
    ts = np.concatenate([np.linspace(0.0, 1.0, 1001), np.arange(257) / 256.0,
                         np.arange(256) / 256.0 - 1e-9, [-0.2, 1.3]])
    for name in fa.COLORMAPS:
        assert len(fa.COLORMAPS[name]) == 256
        ref = mpl.colormaps[name]
        assert ref.N == 256
        for t in ts:
            want = np.asarray(ref(float(np.clip(t, 0.0, 1.0)))[:3])
            got = fa.expected_color(t, name)
            assert np.array_equal(got, want), (name, t, got, want)
