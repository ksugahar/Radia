"""Numpy rasterizers over gmsh-probed fields: ray-cast volume + LIC.

gmsh's renderer has neither a ray-caster nor line integral convolution.
Instead of approximating them inside gmsh (the slice-stack and the
dense-streamline texture remain available for that), these two verbs do
the real thing OUTSIDE gmsh's renderer:

* the FIELD comes from gmsh -- every sample is a ``gmsh.view.probe``
  evaluation on the actual mesh/view, with the measured "found" gating
  (vector list views return a nearest value at any distance; the
  distance is what decides inside/outside);
* the PICTURE is computed in numpy -- ``volume_raycast`` marches rays
  front-to-back through an emission-absorption medium on a regular
  resample grid (the same resample-to-image approach ParaView's GPU
  volume mode uses), and ``lic`` convolves white noise along the vector
  field with RK2 advection so every pixel carries flow direction.

The trade is stated, not hidden: output is a standalone PNG with
labelled axes (axis-equal per lab policy, no in-figure title) and no
gmsh interactivity.  Optional STEP/BREP inputs add depth-composited CAD
to the ray-cast view or a section outline to the LIC view.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess

_MAX_GRID = 128
_MAX_IMAGE_SIZE = 2048
_MAX_DEPTH_STEPS = 2048
_MAX_LIC_RESOLUTION = 2048
_MAX_LIC_KERNEL = 256
_MIN_STEP_REL_SIZE = 0.004
_MAX_CAD_NODES = 500_000
_MAX_CAD_TRIANGLES = 500_000


def _bounded_int(name: str, value: Any, minimum: int, maximum: int) -> int:
    """Return an exact integer inside a resource-safe public range."""
    import operator

    try:
        parsed = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer, got {value!r}") from exc
    if not minimum <= parsed <= maximum:
        raise ValueError(
            f"{name} must be in [{minimum}, {maximum}], got {parsed}")
    return parsed


def _finite_float(name: str, value: Any, *, positive: bool = False,
                  nonnegative: bool = False) -> float:
    """Convert a scalar and reject NaN/Inf before allocating or spawning."""
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number, got {value!r}") \
            from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if positive and parsed <= 0.0:
        raise ValueError(f"{name} must be positive, got {parsed}")
    if nonnegative and parsed < 0.0:
        raise ValueError(f"{name} must be nonnegative, got {parsed}")
    return parsed

_PROBE_GRID_SCRIPT = r"""
import json
import sys

import numpy as np

pts_path, vals_path, found_path, cfg_path, out_path = sys.argv[1:6]
result = {"ok": False}
try:
    with open(cfg_path, encoding="utf-8") as f:
        cfg = json.load(f)
    pts = np.load(pts_path)
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(cfg["path"])
        tags = list(gmsh.view.getTags())
        if not tags:
            raise RuntimeError("no post-processing views in the file")
        sel = cfg.get("view")
        if sel is None:
            tag = tags[0]
        elif isinstance(sel, int):
            if not 0 <= sel < len(tags):
                raise RuntimeError(
                    f"view index {sel} out of range ({len(tags)} views)")
            tag = tags[sel]
        else:
            names = [gmsh.option.getString(
                f"View[{gmsh.view.getIndex(t)}].Name") for t in tags]
            if sel not in names:
                raise RuntimeError(
                    f"no view named {sel!r} (views: {names})")
            tag = tags[names.index(sel)]
        step = int(cfg.get("step", 0))
        dmax = float(cfg.get("distance_max", 0.0))
        n = pts.shape[0]
        vals = np.zeros((n, 9), dtype=np.float64)
        found = np.zeros(n, dtype=bool)
        ncomp = 0
        probe = gmsh.view.probe
        for i in range(n):
            v, dist = probe(tag, float(pts[i, 0]), float(pts[i, 1]),
                            float(pts[i, 2]), step=step, distanceMax=dmax)
            # measured gating: vector LIST views return the nearest
            # value at ANY distance -- found is decided by the distance
            # (len() not truthiness: gmsh may hand back a numpy array)
            if len(v) and float(dist) <= max(dmax, 0.0):
                found[i] = True
                m = min(len(v), 9)
                vals[i, :m] = v[:m]
                if m > ncomp:
                    ncomp = m
        np.save(vals_path, vals[:, :max(ncomp, 1)])
        np.save(found_path, found)
        result = {"ok": True, "ncomp": int(max(ncomp, 1)),
                  "n_found": int(found.sum()), "n_points": int(n)}
    finally:
        gmsh.finalize()
except Exception as exc:                                   # noqa: BLE001
    result = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


_TESSELLATE_SCRIPT = r"""
import json
import sys

import numpy as np

nodes_path, tris_path, cfg_path, out_path = sys.argv[1:5]
result = {"ok": False}
try:
    with open(cfg_path, encoding="utf-8") as f:
        cfg = json.load(f)
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cad")
        for f in cfg["files"]:
            gmsh.model.occ.importShapes(f)
        gmsh.model.occ.synchronize()
        bb = gmsh.model.getBoundingBox(-1, -1)
        diag = max(((bb[3] - bb[0]) ** 2 + (bb[4] - bb[1]) ** 2
                    + (bb[5] - bb[2]) ** 2) ** 0.5, 1e-30)
        # display tessellation only -- never written as a mesh file and
        # never a solver mesh (GMSH mesh generation stays banned for
        # solving; this is what the GUI itself does to shade a face)
        gmsh.option.setNumber("Mesh.MeshSizeMax",
                              diag * float(cfg.get("rel_size", 0.04)))
        gmsh.option.setNumber("Mesh.MeshSizeMin", diag * 0.004)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        xyz = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
        remap = {int(t): i for i, t in enumerate(tags)}
        tris = []
        etypes, _etags, enodes = gmsh.model.mesh.getElements(2)
        for et, conn in zip(etypes, enodes):
            arr = np.asarray(conn, dtype=np.int64)
            if int(et) == 2:                      # 3-node triangles
                tris.append(arr.reshape(-1, 3))
            elif int(et) == 3:                    # quads -> two tris
                q = arr.reshape(-1, 4)
                tris.append(np.concatenate([q[:, [0, 1, 2]],
                                            q[:, [0, 2, 3]]]))
        if not tris:
            raise RuntimeError("tessellation produced no surface elements")
        tri = np.concatenate(tris)
        tri = np.vectorize(remap.__getitem__)(tri)
        np.save(nodes_path, xyz)
        np.save(tris_path, tri.astype(np.int64))
        result = {"ok": True, "n_nodes": int(len(xyz)),
                  "n_triangles": int(len(tri))}
    finally:
        gmsh.finalize()
except Exception as exc:                                   # noqa: BLE001
    result = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _tessellate_step(step_files, *, rel_size: float = 0.04,
                     timeout_s: float = 600.0):
    """Display tessellation of STEP/BREP files -> (nodes, triangles).

    Runs gmsh's OCC import + 2D mesh in a subprocess purely to obtain
    shaded-display triangles (the same thing the gmsh GUI does when it
    shades a geometry face).  Coordinates come through unchanged, so a
    radia/netgen STEP (meters) lands 1:1 on the field data; an external
    mm STEP will be 1000x off -- same caveat as gmsh_render.
    """
    import numpy as np

    rel_size = _finite_float("step_rel_size", rel_size, positive=True)
    if not _MIN_STEP_REL_SIZE <= rel_size <= 1.0:
        raise ValueError(
            f"step_rel_size must be in [{_MIN_STEP_REL_SIZE}, 1], got "
            f"{rel_size}")
    timeout_s = _finite_float("timeout_s", timeout_s, positive=True)
    files = [str(Path(f)) for f in step_files]
    missing = [f for f in files if not Path(f).is_file()]
    if missing:
        return None, None, {"ok": False,
                            "error": f"STEP file(s) not found: {missing}"}
    with tempfile.TemporaryDirectory(prefix="radia_mcp_tess_") as work:
        w = Path(work)
        (w / "cfg.json").write_text(
            json.dumps({"files": files, "rel_size": rel_size}),
            encoding="utf-8")
        res = run_gmsh_json_subprocess(
            _TESSELLATE_SCRIPT,
            [str(w / "nodes.npy"), str(w / "tris.npy"),
             str(w / "cfg.json")],
            timeout_s=timeout_s, prefix="radia_mcp_raster_")
        if not res.get("ok"):
            return None, None, res
        nodes = np.load(w / "nodes.npy", allow_pickle=False)
        tris = np.load(w / "tris.npy", allow_pickle=False)
    if (nodes.ndim != 2 or nodes.shape[1] != 3
            or not np.isfinite(nodes).all()):
        return None, None, {"ok": False,
                            "error": "STEP tessellation returned invalid nodes"}
    if (tris.ndim != 2 or tris.shape[1] != 3
            or not np.issubdtype(tris.dtype, np.integer)
            or (tris.size and (tris.min() < 0 or tris.max() >= len(nodes)))):
        return None, None, {"ok": False,
                            "error": "STEP tessellation returned invalid triangles"}
    if len(nodes) > _MAX_CAD_NODES or len(tris) > _MAX_CAD_TRIANGLES:
        return None, None, {
            "ok": False,
            "error": "STEP display tessellation exceeds the raster safety "
                     "limit; increase step_rel_size",
        }
    return nodes, tris, res


def _cad_depth_buffer(nodes, tris, centre, right, up, d,
                      r0, r1, u0, u1, W, H):
    """Orthographic z-buffer of the CAD triangles in screen space.

    Returns (depth, shade): per-pixel nearest surface depth in
    d-projection units (-inf where no CAD) and the Lambert factor
    |n . d| of that surface.  Pure numpy scanline over each triangle's
    pixel bbox -- a few thousand display triangles rasterize in
    seconds, and per-ray occlusion against the volume then follows
    from comparing depths during the march.
    """
    import numpy as np

    rel = nodes - centre[None, :]
    sx = rel @ right
    sy = rel @ up
    sd = rel @ d
    px = (sx - r0) / (r1 - r0) * W - 0.5
    py = (sy - u0) / (u1 - u0) * H - 0.5

    depth = np.full((H, W), -np.inf)
    shade = np.zeros((H, W))
    v0, v1, v2 = tris[:, 0], tris[:, 1], tris[:, 2]
    e1 = nodes[v1] - nodes[v0]
    e2 = nodes[v2] - nodes[v0]
    n = np.cross(e1, e2)
    nn = np.linalg.norm(n, axis=1)
    ok = nn > 1e-30
    lambert = np.zeros(len(tris))
    lambert[ok] = np.abs((n[ok] / nn[ok, None]) @ d)

    for t in range(len(tris)):
        ia, ib, ic = tris[t]
        xs = (px[ia], px[ib], px[ic])
        ys = (py[ia], py[ib], py[ic])
        x_lo = max(int(np.floor(min(xs))), 0)
        x_hi = min(int(np.ceil(max(xs))), W - 1)
        y_lo = max(int(np.floor(min(ys))), 0)
        y_hi = min(int(np.ceil(max(ys))), H - 1)
        if x_hi < x_lo or y_hi < y_lo:
            continue
        gx, gy = np.meshgrid(np.arange(x_lo, x_hi + 1),
                             np.arange(y_lo, y_hi + 1))
        det = ((ys[1] - ys[2]) * (xs[0] - xs[2])
               + (xs[2] - xs[1]) * (ys[0] - ys[2]))
        if abs(det) < 1e-12:
            continue
        w0 = ((ys[1] - ys[2]) * (gx - xs[2])
              + (xs[2] - xs[1]) * (gy - ys[2])) / det
        w1 = ((ys[2] - ys[0]) * (gx - xs[2])
              + (xs[0] - xs[2]) * (gy - ys[2])) / det
        w2 = 1.0 - w0 - w1
        inside = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9)
        if not inside.any():
            continue
        z = w0 * sd[ia] + w1 * sd[ib] + w2 * sd[ic]
        sub = depth[y_lo:y_hi + 1, x_lo:x_hi + 1]
        upd = inside & (z > sub)
        sub[upd] = z[upd]
        shade[y_lo:y_hi + 1, x_lo:x_hi + 1][upd] = lambert[t]
    return depth, shade


def _probe_grid(path: Path, points, *, view, step: int = 0,
                timeout_s: float = 900.0):
    """Probe a flat (N, 3) point array; returns (values, found, ncomp).

    Values travel through .npy files, not JSON -- a 64^3 grid is two
    million floats and would bloat a JSON round-trip pointlessly.
    """
    import numpy as np

    timeout_s = _finite_float("timeout_s", timeout_s, positive=True)
    pts = np.ascontiguousarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or not np.isfinite(pts).all():
        raise ValueError("probe points must be a finite (N, 3) array")
    with tempfile.TemporaryDirectory(prefix="radia_mcp_probe_grid_") as work:
        w = Path(work)
        np.save(w / "pts.npy", pts)
        cfg = {"path": str(path), "view": view, "step": step,
               "distance_max": 0.0}
        (w / "cfg.json").write_text(json.dumps(cfg), encoding="utf-8")
        res = run_gmsh_json_subprocess(
            _PROBE_GRID_SCRIPT,
            [str(w / "pts.npy"), str(w / "vals.npy"), str(w / "found.npy"),
             str(w / "cfg.json")],
            timeout_s=timeout_s, prefix="radia_mcp_raster_")
        if not res.get("ok"):
            return None, None, res
        vals = np.load(w / "vals.npy", allow_pickle=False)
        found = np.load(w / "found.npy", allow_pickle=False)
    if vals.ndim != 2 or vals.shape[0] != len(pts):
        return None, None, {"ok": False,
                            "error": "gmsh probe returned an invalid value array"}
    if found.shape != (len(pts),) or found.dtype.kind != "b":
        return None, None, {"ok": False,
                            "error": "gmsh probe returned an invalid found mask"}
    if not np.isfinite(vals[found]).all():
        return None, None, {"ok": False,
                            "error": "gmsh probe returned non-finite field values"}
    return vals, found, res


def _magnitude(vals, ncomp: int):
    import numpy as np

    if ncomp == 1:
        return vals[:, 0].copy()
    return np.sqrt((vals[:, :ncomp] ** 2).sum(axis=1))


_AXIS_DIRS = {
    "+x": (1.0, 0.0, 0.0), "-x": (-1.0, 0.0, 0.0),
    "+y": (0.0, 1.0, 0.0), "-y": (0.0, -1.0, 0.0),
    "+z": (0.0, 0.0, 1.0), "-z": (0.0, 0.0, -1.0),
    "iso": (1.0, 1.0, 1.0),
}


def _camera_frame(view_dir):
    """Right-handed screen frame (right, up, d) with d pointing AT the
    camera (out of the screen): right x up = d, so the picture is a
    physically correct view from the d side, never a mirror image."""
    import numpy as np

    if isinstance(view_dir, str):
        if view_dir not in _AXIS_DIRS:
            raise ValueError(
                f"unknown view_dir {view_dir!r} (available: "
                f"{', '.join(sorted(_AXIS_DIRS))}, or a 3-vector)")
        d = np.array(_AXIS_DIRS[view_dir], dtype=float)
    else:
        try:
            d = np.asarray(view_dir, dtype=float).reshape(3)
        except (TypeError, ValueError) as exc:
            raise ValueError("view_dir must be a finite 3-vector") from exc
    if not np.isfinite(d).all():
        raise ValueError("view_dir must be a finite 3-vector")
    nd = float(np.linalg.norm(d))
    if nd < 1e-30:
        raise ValueError("view_dir must be a non-zero vector")
    d = d / nd
    u0 = (np.array([0.0, 1.0, 0.0]) if abs(d[2]) > 0.9
          else np.array([0.0, 0.0, 1.0]))
    right = np.cross(u0, d)
    right /= np.linalg.norm(right)
    up = np.cross(d, right)
    return right, up, d


def _axis_label(vec) -> str:
    """"x [m]" when the screen axis is a world axis, signed; else generic."""
    import numpy as np

    v = np.asarray(vec, dtype=float)
    k = int(np.argmax(np.abs(v)))
    if abs(abs(v[k]) - 1.0) < 1e-9 and np.abs(np.delete(v, k)).max() < 1e-9:
        name = "xyz"[k]
        return f"{name} [m]" if v[k] > 0 else f"-{name} [m]"
    return "screen axis [m]"


def _bbox_from_msh(path: Path):
    from .msh_inspect import read_msh_data

    pts = list(read_msh_data(path)["nodes"].values())
    if not pts:
        raise ValueError(f"{path.name} holds no nodes")
    if any(len(p) != 3 or not all(math.isfinite(float(v)) for v in p)
           for p in pts):
        raise ValueError(f"{path.name} holds non-finite or non-3D nodes")
    lo = [min(p[i] for p in pts) for i in range(3)]
    hi = [max(p[i] for p in pts) for i in range(3)]
    if not any(hi[i] > lo[i] for i in range(3)):
        raise ValueError(f"{path.name} has a degenerate bounding box")
    return lo, hi


def volume_raycast(path: str | Path,
                   png_out: str | Path | None = None, *,
                   view: str | int | None = None,
                   grid: int = 64,
                   view_dir: str | list[float] = "iso",
                   image_size: int = 560,
                   n_steps: int | None = None,
                   value_range: list[float] | None = None,
                   alpha: float = 0.05,
                   alpha_power: float = 2.0,
                   cmap: str = "jet",
                   colorbar: bool = True,
                   step_files: list[str | Path] | None = None,
                   step_color: tuple[float, float, float] = (0.55, 0.57,
                                                             0.62),
                   step_rel_size: float = 0.04,
                   timeout_s: float = 900.0) -> dict[str, Any]:
    """TRUE ray-cast volume rendering (emission-absorption, front-to-back).

    The field is resampled onto a ``grid^3`` regular grid by probing
    gmsh (real field evaluations on the mesh, outside = transparent),
    then orthographic rays march near-to-far through the grid:

        C += T * a_k * colour(v_k);   T *= (1 - a_k)

    with per-sample opacity ``a = alpha * t^alpha_power`` for the
    normalized value t.  This is per-RAY compositing -- occlusion is
    physical, unlike the slice-stack of ``volume_render`` -- and the
    same resample-to-image approach ParaView's GPU volume mode uses.

    Honest limits: the medium lives on the resample grid (``grid``
    controls fidelity, cost is one probe per grid point), opacity is
    defined per depth SAMPLE so the look depends on ``n_steps`` (the
    returned ``transmittance_min`` makes that dependence checkable:
    for a uniform field it equals ``(1-alpha)^n_inside`` exactly), and
    the output is a standalone labelled PNG with no gmsh interactivity.

    Args:
        path: .msh/.pos holding the field.
        png_out: output PNG (default: alongside the input).
        view: view name or index (default: the first view).
        grid: resample resolution per axis (8..128; 64 -> 262k probes).
        view_dir: "+x".."-z", "iso", or a world-space 3-vector pointing
            from the scene TOWARD the camera.
        image_size: image width in pixels (64..2048; height follows aspect).
        n_steps: depth samples per ray (2..2048; default: 1.5 * grid).
        value_range: [lo, hi] normalization (default: the probed
            min/max; pass it explicitly to compare figures).
        alpha: opacity per depth sample at t = 1.
        alpha_power: opacity exponent (2 fades low values out).
        cmap: matplotlib colormap name.
        step_files: STEP/BREP files rendered as OPAQUE shaded surfaces
            INSIDE the volume compositing -- each ray stops where it
            meets the CAD, so the geometry occludes the field behind it
            and the field in front glows over it (the mixed
            geometry+volume scene ParaView builds with its depth
            buffer).  Radia/netgen STEP is meters and lands 1:1;
            external mm CAD is 1000x off (same caveat as gmsh_render).
        step_color: CAD base colour (Lambert-shaded).
        step_rel_size: display-tessellation size vs the CAD diagonal.
    """
    import numpy as np

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    grid = _bounded_int("grid", grid, 8, _MAX_GRID)
    image_size = _bounded_int(
        "image_size", image_size, 64, _MAX_IMAGE_SIZE)
    if n_steps is not None:
        n_steps = _bounded_int(
            "n_steps", n_steps, 2, _MAX_DEPTH_STEPS)
    alpha = _finite_float("alpha", alpha, positive=True)
    if alpha > 1.0:
        raise ValueError(f"alpha must be <= 1, got {alpha}")
    alpha_power = _finite_float(
        "alpha_power", alpha_power, nonnegative=True)
    timeout_s = _finite_float("timeout_s", timeout_s, positive=True)
    step_rel_size = _finite_float(
        "step_rel_size", step_rel_size, positive=True)
    if not _MIN_STEP_REL_SIZE <= step_rel_size <= 1.0:
        raise ValueError(
            f"step_rel_size must be in [{_MIN_STEP_REL_SIZE}, 1], got "
            f"{step_rel_size}")
    try:
        cad_rgb = np.asarray(step_color, dtype=float).reshape(3)
    except (TypeError, ValueError) as exc:
        raise ValueError("step_color must contain three finite RGB values") \
            from exc
    if not np.isfinite(cad_rgb).all() or np.any((cad_rgb < 0.0)
                                                | (cad_rgb > 1.0)):
        raise ValueError("step_color values must be finite and in [0, 1]")
    camera_frame = _camera_frame(view_dir)
    explicit_range: tuple[float, float] | None = None
    if value_range is not None:
        if not isinstance(value_range, (list, tuple)) or len(value_range) != 2:
            raise ValueError("value_range must contain exactly [lo, hi]")
        explicit_range = (
            _finite_float("value_range[0]", value_range[0]),
            _finite_float("value_range[1]", value_range[1]),
        )
        if not explicit_range[1] > explicit_range[0]:
            raise ValueError("value_range must satisfy finite lo < hi")

    lo, hi = _bbox_from_msh(src)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    span = hi - lo
    if not np.all(span > 0.0):
        raise ValueError(f"{src.name} must have a nondegenerate 3D bounding box")

    # --- resample the field onto the regular grid (cell centres) -----
    axes = [lo[i] + span[i] * (np.arange(grid) + 0.5) / grid
            for i in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    vals, found, info = _probe_grid(src, pts, view=view,
                                    timeout_s=timeout_s)
    if vals is None:
        return {"ok": False, "error": f"grid probe failed: {info}"}
    mag = _magnitude(vals, info["ncomp"])
    if not np.isfinite(mag[found]).all():
        return {"ok": False,
                "error": "grid probe returned non-finite field magnitudes"}
    V = mag.reshape(grid, grid, grid)
    M = found.reshape(grid, grid, grid).astype(np.float64)
    if not found.any():
        return {"ok": False,
                "error": "no grid point landed inside the mesh"}

    constant_range = False
    if explicit_range is None:
        v_lo = float(mag[found].min())
        v_hi = float(mag[found].max())
        constant_range = bool(np.isclose(
            v_hi, v_lo, rtol=1.0e-12, atol=1.0e-15))
    else:
        v_lo, v_hi = explicit_range
    if constant_range:
        # A uniform auto-scaled field is the strongest visible value, not
        # an all-zero image.  The found mask still keeps the exterior clear.
        T_norm = M.copy()
        colour_span = max(abs(v_lo), 1.0)
        colorbar_lo = v_lo - colour_span
        colorbar_hi = v_hi
    else:
        T_norm = np.clip((V - v_lo) / (v_hi - v_lo), 0.0, 1.0)
        colorbar_lo, colorbar_hi = v_lo, v_hi

    # --- optional CAD overlay ----------------------------------------
    cad_nodes = cad_tris = None
    cad_info: dict[str, Any] | None = None
    if step_files:
        cad_nodes, cad_tris, cad_info = _tessellate_step(
            step_files, rel_size=step_rel_size, timeout_s=timeout_s)
        if cad_nodes is None:
            return {"ok": False,
                    "error": f"STEP tessellation failed: {cad_info}"}

    # --- camera and image plane --------------------------------------
    right, up, d = camera_frame
    centre = 0.5 * (lo + hi)
    corners = np.array([[lo[0] if i & 1 else hi[0],
                         lo[1] if i & 2 else hi[1],
                         lo[2] if i & 4 else hi[2]] for i in range(8)])
    proj_pts = (corners if cad_nodes is None
                else np.vstack([corners, cad_nodes]))
    rel = proj_pts - centre
    pr, pu, pd = rel @ right, rel @ up, rel @ d
    margin = 0.02 * max(np.ptp(pr), np.ptp(pu))
    r0, r1 = pr.min() - margin, pr.max() + margin
    u0, u1 = pu.min() - margin, pu.max() + margin
    W = image_size
    H = max(round(W * (u1 - u0) / (r1 - r0)), 16)
    if H > _MAX_IMAGE_SIZE:
        raise ValueError(
            f"projected image height {H} exceeds {_MAX_IMAGE_SIZE}; "
            "reduce image_size or choose another view_dir")
    xs = r0 + (r1 - r0) * (np.arange(W) + 0.5) / W
    ys = u0 + (u1 - u0) * (np.arange(H) + 0.5) / H
    SX, SY = np.meshgrid(xs, ys)                      # [iy, ix], y up
    base = (centre[None, None, :] + SX[..., None] * right[None, None, :]
            + SY[..., None] * up[None, None, :])

    n_depth = n_steps if n_steps is not None else int(1.5 * grid)
    ds = float(np.ptp(pd)) / n_depth
    s_near = pd.max() - 0.5 * ds                      # near-to-far

    from matplotlib import colormaps

    cmap_f = colormaps[cmap]
    lut = np.asarray(cmap_f(np.linspace(0.0, 1.0, 256)))[:, :3]

    from scipy.ndimage import map_coordinates

    cad_depth = cad_shade = None
    cad_hit = np.zeros((H, W), dtype=bool)
    if cad_nodes is not None:
        cad_depth, cad_shade = _cad_depth_buffer(
            cad_nodes, cad_tris, centre, right, up, d,
            r0, r1, u0, u1, W, H)
    def _composite_cad(mask, C, T):
        if not mask.any():
            return
        lam = 0.35 + 0.65 * cad_shade[mask]
        C[mask] += T[mask, None] * cad_rgb[None, :] * lam[:, None]
        T[mask] = 0.0

    C = np.zeros((H, W, 3))
    T = np.ones((H, W))
    inside_hits = np.zeros((H, W), dtype=np.int64)
    for k in range(n_depth):
        s_k = s_near - k * ds
        if cad_depth is not None:
            # the ray has reached the CAD surface within this step:
            # composite the opaque surface FIRST, then T = 0 blocks
            # everything behind it (per-ray occlusion, not a paste-on)
            hit_now = (~cad_hit) & (cad_depth >= s_k)
            _composite_cad(hit_now, C, T)
            cad_hit |= hit_now
        P = base + s_k * d[None, None, :]
        idx = ((P - lo[None, None, :]) / span[None, None, :]) * grid - 0.5
        coords = np.stack([idx[..., 0], idx[..., 1], idx[..., 2]])
        t_s = map_coordinates(T_norm, coords, order=1, mode="constant",
                              cval=0.0)
        m_s = map_coordinates(M, coords, order=1, mode="constant",
                              cval=0.0)
        m_s = (m_s > 0.5)
        inside_hits += m_s
        a = alpha * np.power(t_s, alpha_power) * m_s
        col = lut[np.clip((t_s * 255).astype(np.intp), 0, 255)]
        C += (T * a)[..., None] * col
        T *= (1.0 - a)
    if cad_depth is not None:
        # CAD beyond the far marching plane (outside the field bbox)
        _composite_cad((~cad_hit) & (cad_depth > -np.inf), C, T)
        cad_hit |= cad_depth > -np.inf

    img = C + T[..., None]                            # white background

    out = (Path(png_out) if png_out is not None
           else src.with_name(src.stem + "_raycast.png"))
    out.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    fig, ax = plt.subplots(figsize=(6.4, 6.4 * H / W + 0.4), dpi=140)
    ax.imshow(np.clip(img, 0.0, 1.0), origin="lower",
              extent=[r0, r1, u0, u1], aspect="equal",
              interpolation="bilinear")
    ax.set_xlabel(_axis_label(right))
    ax.set_ylabel(_axis_label(up))
    if colorbar:
        sm = ScalarMappable(norm=Normalize(colorbar_lo, colorbar_hi),
                            cmap=cmap_f)
        fig.colorbar(sm, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    result = {"ok": True, "png": str(out), "grid": grid,
              "n_depth_samples": n_depth,
              "value_range": [v_lo, v_hi],
              "constant_auto_range": constant_range,
              "n_probes": int(info["n_points"]),
              "found_fraction": float(found.mean()),
              "transmittance_min": float(T.min()),
              "max_inside_samples": int(inside_hits.max()),
              "view_dir": (view_dir if isinstance(view_dir, str)
                           else [float(v) for v in view_dir]),
              "method": ("ray-cast volume rendering (emission-absorption, "
                         "front-to-back, per-ray occlusion) on a regular "
                         "resample grid probed from gmsh; standalone PNG")}
    if cad_info is not None:
        result["step_files"] = [str(Path(f)) for f in step_files]
        result["cad_triangles"] = cad_info["n_triangles"]
        result["cad_covered_fraction"] = float(cad_hit.mean())
    return result


_LIC_PLANES = {"xy": (0, 1, 2), "yz": (1, 2, 0), "xz": (0, 2, 1)}


def lic(path: str | Path,
        png_out: str | Path | None = None, *,
        view: str | int | None = None,
        plane: str = "xy", offset: float = 0.0,
        resolution: int = 420,
        kernel: int = 18,
        cmap: str = "jet",
        color_by_magnitude: bool = True,
        seed: int = 0,
        step_files: list[str | Path] | None = None,
        step_rel_size: float = 0.03,
        timeout_s: float = 900.0) -> dict[str, Any]:
    """TRUE line integral convolution on a section plane.

    White noise is convolved along the in-plane vector field (RK2
    advection, box kernel of half-length ``kernel`` pixels, forward and
    backward), so EVERY pixel carries the local flow direction -- the
    property the dense-streamline ``flow_texture`` substitute cannot
    provide.  The field samples are gmsh probes on the actual mesh;
    ``color_by_magnitude`` modulates the streaks with |v| through the
    colormap (the "Surface LIC" look).

    Honest limits: the picture is direction texture, not trajectories
    -- individual curves cannot be probed the way ``flow_texture``
    lines can -- and it lives on a regular resample of the plane
    (``resolution`` pixels across the larger side).  Standalone
    labelled PNG, axis-equal, with an optional CAD section outline.

    Args:
        path: .msh/.pos holding a vector view.
        png_out: output PNG (default: alongside the input).
        view: view name or index (default: the first view).
        plane: "xy" | "yz" | "xz" section plane.
        offset: signed offset of the plane from the bbox centre along
            its normal.
        resolution: pixels across the larger in-plane span (64..2048).
        kernel: convolution half-length in pixels (2..256).
        cmap: colormap for the |v| modulation.
        color_by_magnitude: False gives the plain grey LIC texture.
        seed: noise seed (deterministic output).
        step_files: STEP/BREP files whose SECTION OUTLINE (triangle-
            plane intersection of a display tessellation) is drawn in
            black over the texture -- the conductor cross-section on a
            field-line figure.  Meter coordinates, as everywhere.
        step_rel_size: display-tessellation size vs the CAD diagonal.
    """
    import numpy as np

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if plane not in _LIC_PLANES:
        raise ValueError(f"unknown plane {plane!r} "
                         f"(available: {', '.join(sorted(_LIC_PLANES))})")
    resolution = _bounded_int(
        "resolution", resolution, 64, _MAX_LIC_RESOLUTION)
    kernel = _bounded_int("kernel", kernel, 2, _MAX_LIC_KERNEL)
    seed = _bounded_int("seed", seed, 0, 2 ** 64 - 1)
    offset = _finite_float("offset", offset)
    timeout_s = _finite_float("timeout_s", timeout_s, positive=True)
    step_rel_size = _finite_float(
        "step_rel_size", step_rel_size, positive=True)
    if not _MIN_STEP_REL_SIZE <= step_rel_size <= 1.0:
        raise ValueError(
            f"step_rel_size must be in [{_MIN_STEP_REL_SIZE}, 1], got "
            f"{step_rel_size}")

    a0, a1, an = _LIC_PLANES[plane]
    lo, hi = _bbox_from_msh(src)
    span = [hi[i] - lo[i] for i in range(3)]
    if span[a0] <= 0 or span[a1] <= 0:
        return {"ok": False,
                "error": f"{src.name} is degenerate in the {plane} plane"}
    mid = 0.5 * (lo[an] + hi[an]) + offset
    if not lo[an] - 1e-12 <= mid <= hi[an] + 1e-12:
        return {"ok": False,
                "error": f"plane offset {offset} puts the section at "
                         f"{'xyz'[an]} = {mid:.6g}, outside the mesh "
                         f"({lo[an]:.6g} .. {hi[an]:.6g})"}

    # square pixels (axis-equal at the data level)
    h = max(span[a0], span[a1]) / (resolution - 1)
    W = max(round(span[a0] / h) + 1, 8)
    H = max(round(span[a1] / h) + 1, 8)
    us = lo[a0] + h * np.arange(W)
    vs = lo[a1] + h * np.arange(H)
    UU, VV = np.meshgrid(us, vs)                      # [iy, ix]
    pts = np.zeros((H * W, 3))
    pts[:, a0] = UU.ravel()
    pts[:, a1] = VV.ravel()
    pts[:, an] = mid
    vals, found, info = _probe_grid(src, pts, view=view,
                                    timeout_s=timeout_s)
    if vals is None:
        return {"ok": False, "error": f"plane probe failed: {info}"}
    if info["ncomp"] < 3:
        return {"ok": False,
                "error": f"view is scalar ({info['ncomp']} component); "
                         "LIC needs a vector view"}
    if not found.any():
        return {"ok": False, "error": "the section plane misses the mesh"}

    U = vals[:, a0].reshape(H, W)
    Vc = vals[:, a1].reshape(H, W)
    Min = found.reshape(H, W)
    mag = np.sqrt(U * U + Vc * Vc)
    if not np.isfinite(mag[Min]).all():
        return {"ok": False,
                "error": "plane probe returned non-finite field magnitudes"}
    nz = mag > (mag[Min].max() * 1e-9 if Min.any() else 1e-30)
    du = np.where(nz, U / np.where(nz, mag, 1.0), 0.0)
    dv = np.where(nz, Vc / np.where(nz, mag, 1.0), 0.0)

    from scipy.ndimage import map_coordinates

    rng = np.random.default_rng(seed)
    noise = rng.random((H, W))

    def _samp(field, py, px):
        return map_coordinates(field, [py, px], order=1, mode="nearest")

    iy0, ix0 = np.mgrid[0:H, 0:W].astype(np.float64)
    acc = noise.copy()
    wsum = np.ones((H, W))
    for sgn in (1.0, -1.0):
        py, px = iy0.copy(), ix0.copy()
        alive = Min.astype(np.float64)
        for _ in range(kernel):
            # RK2 advection, one pixel per step
            k1u = _samp(du, py, px)
            k1v = _samp(dv, py, px)
            pym = py + sgn * 0.5 * k1v
            pxm = px + sgn * 0.5 * k1u
            k2u = _samp(du, pym, pxm)
            k2v = _samp(dv, pym, pxm)
            px = px + sgn * k2u
            py = py + sgn * k2v
            inside = ((px >= 0) & (px <= W - 1) & (py >= 0) & (py <= H - 1))
            alive = alive * inside * (_samp(Min.astype(np.float64),
                                            py, px) > 0.5)
            acc += alive * _samp(noise, py, px)
            wsum += alive
    tex = acc / wsum

    # contrast stretch (LIC output is intrinsically low-contrast)
    inside_tex = tex[Min]
    p2, p98 = np.percentile(inside_tex, [2.0, 98.0])
    tex = np.clip((tex - p2) / max(p98 - p2, 1e-12), 0.0, 1.0)

    from matplotlib import colormaps

    if color_by_magnitude:
        m_lo = float(mag[Min].min())
        m_hi = float(mag[Min].max())
        m_norm = np.clip((mag - m_lo) / max(m_hi - m_lo, 1e-30), 0.0, 1.0)
        rgb = np.asarray(colormaps[cmap](m_norm))[:, :, :3]
        img = rgb * (0.25 + 0.75 * tex[..., None])
    else:
        m_lo = m_hi = 0.0
        img = np.repeat(tex[..., None], 3, axis=2)
    img = np.where(Min[..., None], img, 1.0)          # outside: white

    out = (Path(png_out) if png_out is not None
           else src.with_name(src.stem + "_lic.png"))
    out.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    # optional CAD section outline (triangle-plane intersection)
    outline_segments: list = []
    if step_files:
        cad_nodes, cad_tris, cad_info = _tessellate_step(
            step_files, rel_size=step_rel_size, timeout_s=timeout_s)
        if cad_nodes is None:
            return {"ok": False,
                    "error": f"STEP tessellation failed: {cad_info}"}
        sd = cad_nodes[:, an] - mid
        sgn = np.sign(sd)
        tri_s = sgn[cad_tris]
        crossing = np.nonzero((tri_s.max(axis=1) > 0)
                              & (tri_s.min(axis=1) < 0))[0]
        for t in crossing:
            pts_uv = []
            for e0, e1 in ((0, 1), (1, 2), (2, 0)):
                i0, i1 = cad_tris[t, e0], cad_tris[t, e1]
                d0, d1 = sd[i0], sd[i1]
                if d0 * d1 < 0.0:
                    w = d0 / (d0 - d1)
                    p = cad_nodes[i0] + w * (cad_nodes[i1] - cad_nodes[i0])
                    pts_uv.append((p[a0], p[a1]))
            if len(pts_uv) == 2:
                outline_segments.append(pts_uv)

    fig, ax = plt.subplots(figsize=(6.4, 6.4 * H / W + 0.4), dpi=140)
    ax.imshow(img, origin="lower",
              extent=[us[0], us[-1], vs[0], vs[-1]], aspect="equal",
              interpolation="bilinear")
    if outline_segments:
        from matplotlib.collections import LineCollection

        ax.add_collection(LineCollection(outline_segments, colors="black",
                                         linewidths=1.1))
    ax.set_xlabel(f"{'xyz'[a0]} [m]")
    ax.set_ylabel(f"{'xyz'[a1]} [m]")
    if color_by_magnitude:
        sm = ScalarMappable(norm=Normalize(m_lo, m_hi),
                            cmap=colormaps[cmap])
        fig.colorbar(sm, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    result = {"ok": True, "png": str(out),
              "size": [int(W), int(H)], "kernel": int(kernel),
              "plane": plane, "offset": float(offset),
              "found_fraction": float(found.mean()),
              "magnitude_range": [m_lo, m_hi],
              "method": ("line integral convolution (box kernel, RK2 "
                         "advection) of white noise along the gmsh-probed "
                         "vector field; every pixel filled; standalone "
                         "PNG")}
    if step_files:
        result["step_files"] = [str(Path(f)) for f in step_files]
        result["step_outline_segments"] = len(outline_segments)
    return result
