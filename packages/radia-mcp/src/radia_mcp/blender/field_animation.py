"""Solver-faithful 3-D field animation through an optional headless Blender.

The adapter reads a ``radia.field-animation.v1`` HDF5 dataset (solved
conductor J(x, t) plus the magnet pose of the same solver steps), drives the
explicit bpy worker in one process to write the ``.blend``, previews and MP4,
then reopens the ``.blend`` in a SECOND process and compares pose, timing,
colour scale and signed J per sample identifier with the dataset.

The worker lives outside the wheel.  Resolution (fail loud, no fallback):
    RADIA_BLENDER_PYTHON  explicit Python 3.13 interpreter with bpy installed
    RADIA_BLENDER_WORKER  path to field_animation_worker.py; defaults to the
                          editable source tree (packages/radia-mcp/tools/blender)
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

SCHEMA = "radia.field-animation.v1"
WORKER_NAME = "field_animation_worker.py"
DEFAULT_SPEC = {
    "fps": 12, "samples": 16, "threads": 8, "resolution": [1280, 720],
    "title": "Moving magnet over a copper plate: eddy current J",
    "arrow_length_m": None, "arrow_width_m": None, "glyph_stride": 1,
    "preview_steps": None, "video": True, "J_scale": "global",
}


class BlenderWorkerError(RuntimeError):
    """The optional worker is missing, misconfigured or failed."""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_worker():
    """Return (interpreter, worker_script) or raise with the setting to fix."""
    py = os.environ.get("RADIA_BLENDER_PYTHON", "").strip()
    if not py:
        raise BlenderWorkerError(
            "RADIA_BLENDER_PYTHON is not set: point it to a Python 3.13 interpreter "
            "with bpy installed (pip install -r packages/radia-mcp/tools/blender/requirement.txt).")
    if not Path(py).is_file():
        raise BlenderWorkerError(f"RADIA_BLENDER_PYTHON does not exist: {py}")
    worker = os.environ.get("RADIA_BLENDER_WORKER", "").strip()
    if worker:
        wpath = Path(worker)
    else:
        wpath = Path(__file__).resolve().parents[3] / "tools" / "blender" / WORKER_NAME
    if not wpath.is_file():
        raise BlenderWorkerError(
            f"Blender worker not found at {wpath}; the worker is not shipped in the wheel. "
            "Set RADIA_BLENDER_WORKER to tools/blender/field_animation_worker.py of a source checkout.")
    probe = subprocess.run([py, "-I", "-B", "-c",
                            "import sys, bpy; print(sys.version_info[:2] == (3, 13), bpy.app.version_string)"],
                           capture_output=True, text=True, timeout=300,
                           stdin=subprocess.DEVNULL)  # never inherit the MCP stdio pipe
    if probe.returncode != 0 or not probe.stdout.startswith("True"):
        raise BlenderWorkerError(
            f"RADIA_BLENDER_PYTHON is not Python 3.13 with bpy: rc={probe.returncode} "
            f"{probe.stdout.strip()} {probe.stderr.strip()[-400:]}")
    return py, wpath, probe.stdout.split(None, 1)[1].strip()


def load_dataset(path):
    """Read and check a field-animation dataset; return a JSON-able payload."""
    import h5py
    with h5py.File(path, "r") as f:
        schema = f.attrs.get("schema")
        if isinstance(schema, bytes):
            schema = schema.decode()
        if schema != SCHEMA:
            raise ValueError(f"dataset schema {schema!r} is not {SCHEMA!r}")
        t = f["time_s"][()]
        J = f["field/J_A_per_m2"][()]
        ids = f["samples/id"][()]
        pts = f["samples/point_m"][()]
        center = f["magnet/center_m"][()]
        rot = f["magnet/rotation"][()]
        n, m = len(t), len(ids)
        checks = [(J.shape == (n, m, 3), f"J shape {J.shape} != ({n}, {m}, 3)"),
                  (pts.shape == (m, 3), f"sample points {pts.shape}"),
                  (center.shape == (n, 3), f"magnet centre {center.shape}"),
                  (rot.shape == (n, 3, 3), f"magnet rotation {rot.shape}"),
                  (len(np.unique(ids)) == m, "sample ids are not unique"),
                  (n >= 2 and np.all(np.diff(t) > 0), "time_s must increase"),
                  (np.all(np.isfinite(J)), "J contains non-finite values")]
        for ok, msg in checks:
            if not ok:
                raise ValueError(f"invalid dataset: {msg}")
        orth = np.abs(np.einsum("nji,njk->nik", rot, rot) - np.eye(3)).max()
        if orth > 1e-9:
            raise ValueError(f"magnet rotation is not orthonormal (max error {orth:.2e})")
        prov = {k: (v.decode() if isinstance(v, bytes) else
                    v.item() if hasattr(v, "item") else v) for k, v in f.attrs.items()}
        return {
            "schema": SCHEMA, "provenance": prov, "time_s": t.tolist(),
            "magnet": {"center_m": center.tolist(), "rotation": rot.tolist(),
                       "size_m": f["magnet/size_m"][()].tolist(),
                       "magnetization_A_per_m": f["magnet/magnetization_A_per_m"][()].tolist()},
            "conductor": {"box_min_m": f["conductor/box_min_m"][()].tolist(),
                          "box_max_m": f["conductor/box_max_m"][()].tolist()},
            "samples": {"id": ids.tolist(), "point_m": pts.tolist()},
            "J_A_per_m2": J.tolist(),
        }


def _run(py, worker, mode, req, workdir, timeout):
    rq = Path(workdir) / f"{mode}_request.json"
    req["result_path"] = str(Path(workdir) / f"{mode}_result.json")
    rq.write_text(json.dumps(req), encoding="utf-8")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    try:
        p = subprocess.run([py, "-I", "-B", str(worker), mode, str(rq)], capture_output=True,
                           text=True, timeout=timeout, env=env, cwd=workdir,
                           stdin=subprocess.DEVNULL)  # never inherit the MCP stdio pipe
    except subprocess.TimeoutExpired as exc:
        raise BlenderWorkerError(f"worker {mode} exceeded {timeout} s") from exc
    if p.returncode != 0:
        raise BlenderWorkerError(f"worker {mode} failed rc={p.returncode}: {p.stderr.strip()[-1500:]}")
    return json.loads(Path(req["result_path"]).read_text(encoding="utf-8"))


def field_animation(dataset_path, out_dir, spec=None, timeout_s=3600):
    """Write, reopen and verify one animation; return the result record."""
    dataset_path = Path(dataset_path).resolve()
    out = Path(out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    s = dict(DEFAULT_SPEC)
    unknown = set(spec or {}) - set(DEFAULT_SPEC)
    if unknown:
        raise ValueError(f"unknown spec keys: {sorted(unknown)}")
    s.update(spec or {})
    py, worker, blender_version = resolve_worker()
    pay = load_dataset(dataset_path)
    J = np.asarray(pay["J_A_per_m2"])
    n = J.shape[0]
    mag = np.linalg.norm(J, axis=-1)
    if s["J_scale"] == "global":
        jmax = float(mag.max())
    else:
        jmax = float(s["J_scale"])
    if not jmax > 0:
        raise ValueError(f"|J| scale must be positive, got {jmax}")
    pay["J_scale_max_A_per_m2"] = jmax
    stride = int(s["glyph_stride"])
    if stride > 1:
        keep = slice(None, None, stride)
        pay["samples"] = {k: v[keep] for k, v in pay["samples"].items()}
        pay["J_A_per_m2"] = J[:, keep].tolist()
    lo = np.array(pay["conductor"]["box_min_m"]); hi = np.array(pay["conductor"]["box_max_m"])
    pitch = (hi[0] - lo[0]) / max(np.sqrt(len(pay["samples"]["id"])), 1.0)
    s["arrow_length_m"] = s["arrow_length_m"] or 1.5 * pitch
    s["arrow_width_m"] = s["arrow_width_m"] or 0.5 * pitch
    if s["preview_steps"] is None:
        s["preview_steps"] = sorted({1, n // 4, n // 2, (3 * n) // 4, n - 1})
    stem = dataset_path.stem
    with tempfile.TemporaryDirectory(prefix="radia_blender_") as wd:
        pp = Path(wd) / "payload.json"
        pp.write_text(json.dumps(pay), encoding="utf-8")
        wreq = {"payload_path": str(pp), "spec": s,
                "blend_path": str(out / f"{stem}.blend"),
                "preview_pattern": str(out / f"{stem}_step{{step:04d}}.png"),
                "video_path": str(out / f"{stem}.mp4") if s["video"] else None}
        written = _run(py, worker, "write", wreq, wd, timeout_s)
        check = sorted(set([0] + list(s["preview_steps"])))
        back = _run(py, worker, "inspect", {"blend_path": written["blend_path"], "steps": check},
                    wd, 600)
    verification = compare(pay, back, s, check)
    artifacts = {"blend": written["blend_path"], "previews": written["previews"],
                 "video": written["video_path"]}
    hashes = {p: _sha256(p) for p in [artifacts["blend"], *artifacts["previews"],
                                      *( [artifacts["video"]] if artifacts["video"] else [])]}
    return {
        "status": "ok" if verification["passed"] else "mismatch",
        "artifacts": artifacts, "sha256": hashes,
        "dataset": {"path": str(dataset_path), "sha256": _sha256(dataset_path),
                    "provenance": pay["provenance"], "steps": n,
                    "samples_drawn": len(pay["samples"]["id"])},
        "blender": blender_version, "worker_sha256": _sha256(worker),
        "spec": s, "J_scale_A_per_m2": [0.0, jmax],
        "scene": {k: back.get(k) for k in ("render", "length_unit", "bbox_m", "interpolations", "objects")},
        "timing": {"fps": s["fps"], "frame_start": back["frame_start"],
                   "frame_end": back["frame_end"],
                   "interpolation": "none: pose and J held constant between solver steps",
                   "camera": back["camera_motion"]},
        "readback": verification,
        "note": "Rendering readability is not numerical validation; the J values are the solver's.",
    }


def compare(pay, back, spec, steps):
    """Independent comparison of the reopened .blend with the dataset."""
    n = len(pay["time_s"])
    fails = []
    if back["frame_end"] - back["frame_start"] + 1 != n:
        fails.append(f"frame count {back['frame_end'] - back['frame_start'] + 1} != steps {n}")
    if back["fps"] != spec["fps"]:
        fails.append(f"fps {back['fps']} != {spec['fps']}")
    if back["schema"] != SCHEMA:
        fails.append(f"schema {back['schema']!r}")
    if back["camera_animated"] or back["camera_motion"] != "static":
        fails.append("camera is animated")
    if back.get("interpolations") != ["CONSTANT"]:
        fails.append(f"keyframe interpolation {back.get('interpolations')} is not CONSTANT only")
    lo = np.array(pay["conductor"]["box_min_m"]); hi = np.array(pay["conductor"]["box_max_m"])
    cb = back.get("bbox_m", {}).get("Conductor")
    if cb is None or np.abs(np.array(cb) - np.array([lo, hi])).max() > 1e-8:
        fails.append(f"conductor bounding box {cb} != dataset {[lo.tolist(), hi.tolist()]}")
    if back["n_glyph_objects"] != n:
        fails.append(f"glyph objects {back['n_glyph_objects']} != {n}")
    dt = np.abs(np.array(back["physical_time_s"]) - np.array(pay["time_s"])).max()
    if dt > 1e-12:
        fails.append(f"physical time differs by {dt:.3e} s")
    if abs(back["J_scale_A_per_m2"][1] - pay["J_scale_max_A_per_m2"]) > 1e-9 * pay["J_scale_max_A_per_m2"]:
        fails.append("|J| scale differs")
    per_step = {}
    ids_ref = np.asarray(pay["samples"]["id"])
    for k in steps:
        b = back["steps"][str(k)]
        ids = np.asarray(b["sample_id"])
        order = np.argsort(ids)
        if not np.array_equal(ids[order], np.sort(ids_ref)):
            fails.append(f"step {k}: sample identifiers differ")
            continue
        Jb = np.asarray(b["J_A_per_m2"]).reshape(-1, 3)
        pos = {int(i): j for j, i in enumerate(ids_ref)}
        Jr = np.asarray(pay["J_A_per_m2"][k])[[pos[int(i)] for i in ids]]
        scale = max(float(np.abs(Jr).max()), 1e-300)
        jerr = float(np.abs(Jb - Jr).max() / scale) if np.abs(Jr).max() > 0 else float(np.abs(Jb).max())
        perr = float(np.abs(np.asarray(b["center_m"]) - np.asarray(pay["magnet"]["center_m"][k])).max())
        rerr = float(np.abs(np.asarray(b["rotation"]) - np.asarray(pay["magnet"]["rotation"][k])).max())
        vis = [v for v in b["visible"] if v.startswith("J_glyph_")]
        ok_vis = vis == [f"J_glyph_{k:04d}"]
        label_ok = f"{pay['time_s'][k]:.3f} s" in b["time_label"]
        per_step[str(k)] = {"J_max_rel_err": jerr, "pose_center_max_err_m": perr,
                       "pose_rotation_max_err": rerr, "visible_glyphs": vis,
                       "time_label": b["time_label"]}
        if jerr > 1e-6:
            fails.append(f"step {k}: J readback rel err {jerr:.2e}")
        if perr > 1e-8 or rerr > 1e-6:
            fails.append(f"step {k}: pose readback err {perr:.2e} m / {rerr:.2e}")
        if not ok_vis:
            fails.append(f"step {k}: visible glyph objects {vis}")
        if not label_ok:
            fails.append(f"step {k}: time label {b['time_label']!r}")
    return {"passed": not fails, "failures": fails, "steps": per_step,
            "tolerances": {"J_rel": 1e-6, "center_m": 1e-8, "rotation": 1e-6,
                           "note": "float32 mesh coordinates; J stored as float attribute"}}
