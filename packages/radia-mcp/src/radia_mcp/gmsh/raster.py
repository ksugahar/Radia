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
labelled axes (axis-equal per lab policy, no in-figure title) -- there
is no CAD overlay and no gmsh interactivity in these two figures.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess

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


def _probe_grid(path: Path, points, *, view, step: int = 0,
                timeout_s: float = 900.0):
    """Probe a flat (N, 3) point array; returns (values, found, ncomp).

    Values travel through .npy files, not JSON -- a 64^3 grid is two
    million floats and would bloat a JSON round-trip pointlessly.
    """
    import numpy as np

    pts = np.ascontiguousarray(points, dtype=np.float64)
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
        vals = np.load(w / "vals.npy")
        found = np.load(w / "found.npy")
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
        d = np.asarray(view_dir, dtype=float).reshape(3)
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
    lo = [min(p[i] for p in pts) for i in range(3)]
    hi = [max(p[i] for p in pts) for i in range(3)]
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
    the output is a standalone labelled PNG -- no CAD overlay, no gmsh
    interactivity.

    Args:
        path: .msh/.pos holding the field.
        png_out: output PNG (default: alongside the input).
        view: view name or index (default: the first view).
        grid: resample resolution per axis (64 -> 262k probes).
        view_dir: "+x".."-z", "iso", or a world-space 3-vector pointing
            from the scene TOWARD the camera.
        image_size: image width in pixels (height follows the aspect).
        n_steps: depth samples per ray (default: 1.5 * grid).
        value_range: [lo, hi] normalization (default: the probed
            min/max; pass it explicitly to compare figures).
        alpha: opacity per depth sample at t = 1.
        alpha_power: opacity exponent (2 fades low values out).
        cmap: matplotlib colormap name.
    """
    import numpy as np

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if grid < 8:
        raise ValueError(f"grid must be >= 8, got {grid}")
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    if image_size < 64:
        raise ValueError(f"image_size must be >= 64, got {image_size}")

    lo, hi = _bbox_from_msh(src)
    lo = np.asarray(lo)
    hi = np.asarray(hi)
    span = np.maximum(hi - lo, 1e-30)

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
    V = mag.reshape(grid, grid, grid)
    M = found.reshape(grid, grid, grid).astype(np.float64)
    if not found.any():
        return {"ok": False,
                "error": "no grid point landed inside the mesh"}

    if value_range is None:
        v_lo = float(mag[found].min())
        v_hi = float(mag[found].max())
    else:
        v_lo, v_hi = (float(v) for v in value_range)
    if not v_hi > v_lo:
        v_hi = v_lo + 1.0
    T_norm = np.clip((V - v_lo) / (v_hi - v_lo), 0.0, 1.0)

    # --- camera and image plane --------------------------------------
    right, up, d = _camera_frame(view_dir)
    centre = 0.5 * (lo + hi)
    corners = np.array([[lo[0] if i & 1 else hi[0],
                         lo[1] if i & 2 else hi[1],
                         lo[2] if i & 4 else hi[2]] for i in range(8)])
    rel = corners - centre
    pr, pu, pd = rel @ right, rel @ up, rel @ d
    margin = 0.02 * max(np.ptp(pr), np.ptp(pu))
    r0, r1 = pr.min() - margin, pr.max() + margin
    u0, u1 = pu.min() - margin, pu.max() + margin
    W = int(image_size)
    H = max(int(round(W * (u1 - u0) / (r1 - r0))), 16)
    xs = r0 + (r1 - r0) * (np.arange(W) + 0.5) / W
    ys = u0 + (u1 - u0) * (np.arange(H) + 0.5) / H
    SX, SY = np.meshgrid(xs, ys)                      # [iy, ix], y up
    base = (centre[None, None, :] + SX[..., None] * right[None, None, :]
            + SY[..., None] * up[None, None, :])

    n_depth = int(n_steps) if n_steps is not None else int(1.5 * grid)
    if n_depth < 2:
        raise ValueError(f"n_steps must be >= 2, got {n_steps}")
    ds = float(np.ptp(pd)) / n_depth
    s_near = pd.max() - 0.5 * ds                      # near-to-far

    from matplotlib import colormaps

    cmap_f = colormaps[cmap]
    lut = np.asarray(cmap_f(np.linspace(0.0, 1.0, 256)))[:, :3]

    from scipy.ndimage import map_coordinates

    C = np.zeros((H, W, 3))
    T = np.ones((H, W))
    inside_hits = np.zeros((H, W), dtype=np.int64)
    for k in range(n_depth):
        P = base + (s_near - k * ds) * d[None, None, :]
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
        sm = ScalarMappable(norm=Normalize(v_lo, v_hi), cmap=cmap_f)
        fig.colorbar(sm, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    return {"ok": True, "png": str(out), "grid": grid,
            "n_depth_samples": n_depth,
            "value_range": [v_lo, v_hi],
            "n_probes": int(info["n_points"]),
            "found_fraction": float(found.mean()),
            "transmittance_min": float(T.min()),
            "max_inside_samples": int(inside_hits.max()),
            "view_dir": (view_dir if isinstance(view_dir, str)
                         else [float(v) for v in view_dir]),
            "method": ("ray-cast volume rendering (emission-absorption, "
                       "front-to-back, per-ray occlusion) on a regular "
                       "resample grid probed from gmsh; standalone PNG, "
                       "no CAD overlay")}


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
    labelled PNG, axis-equal; no CAD overlay.

    Args:
        path: .msh/.pos holding a vector view.
        png_out: output PNG (default: alongside the input).
        view: view name or index (default: the first view).
        plane: "xy" | "yz" | "xz" section plane.
        offset: signed offset of the plane from the bbox centre along
            its normal.
        resolution: pixels across the larger in-plane span (>= 64).
        kernel: convolution half-length in pixels (>= 2).
        cmap: colormap for the |v| modulation.
        color_by_magnitude: False gives the plain grey LIC texture.
        seed: noise seed (deterministic output).
    """
    import numpy as np

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if plane not in _LIC_PLANES:
        raise ValueError(f"unknown plane {plane!r} "
                         f"(available: {', '.join(sorted(_LIC_PLANES))})")
    if resolution < 64:
        raise ValueError(f"resolution must be >= 64, got {resolution}")
    if kernel < 2:
        raise ValueError(f"kernel must be >= 2, got {kernel}")

    a0, a1, an = _LIC_PLANES[plane]
    lo, hi = _bbox_from_msh(src)
    span = [hi[i] - lo[i] for i in range(3)]
    if span[a0] <= 0 or span[a1] <= 0:
        return {"ok": False,
                "error": f"{src.name} is degenerate in the {plane} plane"}
    mid = 0.5 * (lo[an] + hi[an]) + float(offset)
    if not lo[an] - 1e-12 <= mid <= hi[an] + 1e-12:
        return {"ok": False,
                "error": f"plane offset {offset} puts the section at "
                         f"{'xyz'[an]} = {mid:.6g}, outside the mesh "
                         f"({lo[an]:.6g} .. {hi[an]:.6g})"}

    # square pixels (axis-equal at the data level)
    h = max(span[a0], span[a1]) / (resolution - 1)
    W = max(int(round(span[a0] / h)) + 1, 8)
    H = max(int(round(span[a1] / h)) + 1, 8)
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

    fig, ax = plt.subplots(figsize=(6.4, 6.4 * H / W + 0.4), dpi=140)
    ax.imshow(img, origin="lower",
              extent=[us[0], us[-1], vs[0], vs[-1]], aspect="equal",
              interpolation="bilinear")
    ax.set_xlabel(f"{'xyz'[a0]} [m]")
    ax.set_ylabel(f"{'xyz'[a1]} [m]")
    if color_by_magnitude:
        sm = ScalarMappable(norm=Normalize(m_lo, m_hi),
                            cmap=colormaps[cmap])
        fig.colorbar(sm, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    return {"ok": True, "png": str(out),
            "size": [int(W), int(H)], "kernel": int(kernel),
            "plane": plane, "offset": float(offset),
            "found_fraction": float(found.mean()),
            "magnitude_range": [m_lo, m_hi],
            "method": ("line integral convolution (box kernel, RK2 "
                       "advection) of white noise along the gmsh-probed "
                       "vector field; every pixel filled; standalone "
                       "PNG, no CAD overlay")}
