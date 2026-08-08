"""Headless gmsh rendering: PNG screenshots and animation export.

Codifies the lab animation/display knowledge as executable tools:

- Every render runs the gmsh Python API in a SUBPROCESS with an FLTK
  graphics context (``gmsh.write`` of image formats needs OpenGL; a
  crashing gmsh must not kill the MCP server).
- ``-noconfig`` + explicit window geometry prevents the stale off-screen
  window-position pitfall (GMSH restores last-used coordinates from
  %APPDATA%/gmsh-options, which may point at a detached monitor).
- ``Mesh.NumSubEdges = 4`` by default so curved high-order MESH edges
  render curved (the whole point of high-order export).
- ``View[i].AdaptVisualizationGrid = 1`` by default on post views:
  elements with more than 8 nodes (TET10, HEX20, ...) are otherwise
  SILENTLY skipped by the view renderer, and view display ignores
  ``Mesh.NumSubEdges`` entirely.
- Multi-view animations set ``PostProcessing.Link = 1`` and
  ``AnimationCycle = 0`` and drive every target view's ``TimeStep``
  explicitly (Cycle=1 cycles views instead of synchronizing them).
"""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import run_gmsh_json_subprocess
from .post_display import CAMERA_PRESETS, cut_plane_option_values

_RENDER_SCRIPT = r"""
import json
import sys

cfg_path, out_path = sys.argv[1], sys.argv[2]
with open(cfg_path, encoding="utf-8") as f:
    cfg = json.load(f)
result = {"ok": False, "ran": False}
try:
    import gmsh
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.GraphicsWidth", cfg["width"])
        gmsh.option.setNumber("General.GraphicsHeight", cfg["height"])
        gmsh.option.setNumber("General.GraphicsPositionX", 100)
        gmsh.option.setNumber("General.GraphicsPositionY", 100)
        gmsh.option.setNumber("General.MenuPositionX", 100)
        gmsh.option.setNumber("General.MenuPositionY", 100)
        gmsh.open(cfg["path"])
        for extra in cfg.get("merge_files") or []:
            gmsh.merge(extra)
        gmsh.option.setNumber("Mesh.NumSubEdges", cfg["numsubedges"])

        if cfg.get("geometry_display"):
            # merged CAD (STEP/BREP) shows as shaded solid faces; points
            # off (vertex dots clutter overlay figures)
            gmsh.option.setNumber("Geometry.Surfaces", 1)
            gmsh.option.setNumber("Geometry.SurfaceType", 2)
            gmsh.option.setNumber("Geometry.Points", 0)

        view_tags = list(gmsh.view.getTags())
        n_views = len(view_tags)
        result["n_views"] = n_views

        if cfg.get("auto_mesh_display") and n_views == 0:
            gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
            gmsh.option.setNumber("Mesh.ColorCarousel", 2)

        rot = cfg.get("rotation")
        if rot is not None:
            gmsh.option.setNumber("General.Trackball", 0)
            gmsh.option.setNumber("General.RotationX", rot[0])
            gmsh.option.setNumber("General.RotationY", rot[1])
            gmsh.option.setNumber("General.RotationZ", rot[2])
            gmsh.option.setNumber("General.RotationCenterGravity", 1)

        if cfg.get("adapt_views"):
            for i in range(n_views):
                gmsh.option.setNumber(f"View[{i}].AdaptVisualizationGrid", 1)
                gmsh.option.setNumber(f"View[{i}].MaxRecursionLevel", 2)
                gmsh.option.setNumber(f"View[{i}].TargetError", 1e-4)

        if cfg.get("smooth_normals", True):
            # gmsh default is OFF, which leaves isosurfaces and skins
            # visibly faceted; per-vertex normal averaging is the
            # single biggest surface-quality lever.
            for i in range(n_views):
                gmsh.option.setNumber(f"View[{i}].SmoothNormals", 1)

        if cfg.get("clip_views"):
            for i in range(n_views):
                gmsh.option.setNumber(f"View[{i}].Clip", 1)

        # ---- structured figure controls -----------------------------
        # Applied BEFORE the raw option passthrough, so an explicit
        # options={} entry always wins over the structured form.
        def _targets(spec):
            if spec is None:
                return list(range(n_views))
            return [int(i) for i in spec]

        col = cfg.get("color")
        if col:
            idxs = _targets(col.get("views"))
            rng = col.get("range")
            if rng == "shared" and idxs:
                lo = min(gmsh.option.getNumber(f"View[{i}].Min") for i in idxs)
                hi = max(gmsh.option.getNumber(f"View[{i}].Max") for i in idxs)
                rng = [lo, hi]
            if isinstance(rng, (list, tuple)):
                result["color_range"] = [float(rng[0]), float(rng[1])]
            for i in idxs:
                if isinstance(rng, (list, tuple)):
                    gmsh.option.setNumber(f"View[{i}].RangeType", 2)
                    gmsh.option.setNumber(f"View[{i}].CustomMin", float(rng[0]))
                    gmsh.option.setNumber(f"View[{i}].CustomMax", float(rng[1]))
                for key, opt in (("scale_type", "ScaleType"),
                                 ("intervals", "NbIso"),
                                 ("intervals_type", "IntervalsType"),
                                 ("colormap", "ColormapNumber"),
                                 ("alpha", "ColormapAlpha"),
                                 ("show_scale", "ShowScale"),
                                 ("saturate", "SaturateValues"),
                                 ("invert", "ColormapInvert")):
                    if col.get(key) is not None:
                        gmsh.option.setNumber(f"View[{i}].{opt}",
                                              float(col[key]))
                if col.get("format"):
                    gmsh.option.setString(f"View[{i}].Format", col["format"])

        gly = cfg.get("glyphs")
        if gly:
            for i in _targets(gly.get("views")):
                for key, opt in (("vector_type", "VectorType"),
                                 ("sampling", "Sampling"),
                                 ("size_max", "ArrowSizeMax"),
                                 ("size_min", "ArrowSizeMin"),
                                 ("center", "CenterGlyphs"),
                                 ("location", "GlyphLocation"),
                                 ("line_width", "LineWidth"),
                                 ("point_size", "PointSize")):
                    if gly.get(key) is not None:
                        gmsh.option.setNumber(f"View[{i}].{opt}",
                                              float(gly[key]))

        planes = cfg.get("clip_planes") or []
        for k, pl in enumerate(planes):
            for letter, val in zip("ABCD", pl["plane"]):
                gmsh.option.setNumber(f"General.Clip{k}{letter}", float(val))
            mask = 1 << k
            for target in pl["targets"]:
                cur = gmsh.option.getNumber(target)
                gmsh.option.setNumber(target, float(int(cur) | mask))
            if pl.get("whole_elements"):
                gmsh.option.setNumber("General.ClipWholeElements", 1)

        ax = cfg.get("axes")
        if ax:
            gmsh.option.setNumber("General.Axes", float(ax["mode"]))
            for j, letter in enumerate("XYZ"):
                if ax.get("labels"):
                    gmsh.option.setString(f"General.AxesLabel{letter}",
                                          ax["labels"][j])
                if ax.get("format"):
                    gmsh.option.setString(f"General.AxesFormat{letter}",
                                          ax["format"][j])
                if ax.get("tics"):
                    gmsh.option.setNumber(f"General.AxesTics{letter}",
                                          float(ax["tics"][j]))

        annots = cfg.get("annotations") or []
        if annots:
            atag = gmsh.view.add("annotations")
            for a in annots:
                gmsh.view.addListDataString(
                    atag, [a["x"], a["y"]], [a["text"]],
                    ["Align", a.get("align", "Left"),
                     "Font", a.get("font", "Helvetica"),
                     "FontSize", str(int(a.get("size", 18)))])
            aidx = gmsh.view.getIndex(atag)
            gmsh.option.setNumber(f"View[{aidx}].ShowScale", 0)
            result["annotation_view"] = aidx

        for name, val in (cfg.get("options") or {}).items():
            gmsh.option.setNumber(name, float(val))
        for name, val in (cfg.get("string_options") or {}).items():
            gmsh.option.setString(name, str(val))

        def _revert_empty_adaptive_views():
            # Adaptive visualization skips LINE-element views entirely
            # (PEEC filament exports render as an empty 1e+200 range).
            # After the first draw, any view whose computed range is
            # empty gets its adaptation reverted and the scene redrawn.
            if not cfg.get("adapt_views"):
                return []
            reverted = []
            for i in range(n_views):
                vmin = gmsh.option.getNumber(f"View[{i}].Min")
                vmax = gmsh.option.getNumber(f"View[{i}].Max")
                if vmin > vmax:
                    gmsh.option.setNumber(
                        f"View[{i}].AdaptVisualizationGrid", 0)
                    reverted.append(i)
            if reverted:
                gmsh.fltk.update()
            return reverted

        if cfg["mode"] == "png":
            ts = cfg.get("time_step")
            if ts is not None:
                for i in range(n_views):
                    gmsh.option.setNumber(f"View[{i}].TimeStep", int(ts))
            gmsh.fltk.initialize()
            gmsh.fltk.update()
            result["adapt_reverted_views"] = _revert_empty_adaptive_views()
            gmsh.write(cfg["png_out"])
            gmsh.fltk.finalize()
            result.update({"ok": True, "ran": True, "png": cfg["png_out"]})
            # Self-check: an image that is one flat color (plus the tiny
            # axes triad) means nothing was actually drawn.
            try:
                from PIL import Image
                img = Image.open(cfg["png_out"]).convert("RGB")
                w, h = img.size
                colors = img.getcolors(maxcolors=1 << 22)
                if colors:
                    background = max(count for count, _ in colors)
                    content = 1.0 - background / float(w * h)
                else:
                    content = 1.0
                result["blank_check"] = {
                    "ran": True,
                    "content_fraction": round(content, 5),
                    "looks_blank": content < 0.005,
                }
            except ImportError:
                result["blank_check"] = {
                    "ran": False, "error": "Pillow not installed"}
        else:
            frames_dir = cfg["frames_dir"]
            frames = []
            orbit = cfg.get("orbit")
            if orbit:
                # Camera-orbit animation: the data stays fixed (one time
                # step); the camera sweeps `degrees` around one axis.
                gmsh.option.setNumber("General.Trackball", 0)
                gmsh.option.setNumber("General.RotationCenterGravity", 1)
                base = cfg.get("rotation") or [0.0, 0.0, 0.0]
                axis_opt = {"x": "General.RotationX",
                            "y": "General.RotationY",
                            "z": "General.RotationZ"}[orbit["axis"]]
                axis_idx = {"x": 0, "y": 1, "z": 2}[orbit["axis"]]
                ts = cfg.get("time_step")
                if ts is not None:
                    for i in range(n_views):
                        gmsh.option.setNumber(f"View[{i}].TimeStep",
                                              int(ts))
                n_frames = int(orbit["n_frames"])
                gmsh.fltk.initialize()
                gmsh.fltk.update()
                result["adapt_reverted_views"] = \
                    _revert_empty_adaptive_views()
                for k in range(n_frames):
                    gmsh.option.setNumber(
                        axis_opt,
                        base[axis_idx] + orbit["degrees"] * k / n_frames)
                    gmsh.fltk.update()
                    frame = frames_dir + f"/frame_{k:04d}.png"
                    gmsh.write(frame)
                    frames.append(frame)
                gmsh.fltk.finalize()
                result.update({"ok": True, "ran": True,
                               "num_steps": n_frames, "frames": frames,
                               "orbit": orbit})
            else:
                targets = cfg.get("view_indices")
                if targets is None:
                    targets = list(range(n_views))
                if not targets:
                    raise RuntimeError(
                        "no post-processing views to animate (the file has "
                        "no NodeData/ElementData sections); for a mesh-only "
                        "fly-around pass orbit_axis instead")
                steps = cfg.get("num_steps")
                if steps is None:
                    steps = max(int(gmsh.option.getNumber(
                        f"View[{i}].NbTimeStep")) for i in targets)
                if cfg.get("link_views", True):
                    gmsh.option.setNumber("PostProcessing.Link", 1)
                    gmsh.option.setNumber("PostProcessing.AnimationCycle", 0)
                gmsh.fltk.initialize()
                gmsh.fltk.update()
                result["adapt_reverted_views"] = \
                    _revert_empty_adaptive_views()
                for step in range(int(steps)):
                    for i in targets:
                        gmsh.option.setNumber(f"View[{i}].TimeStep", step)
                    gmsh.fltk.update()
                    frame = frames_dir + f"/frame_{step:04d}.png"
                    gmsh.write(frame)
                    frames.append(frame)
                gmsh.fltk.finalize()
                result.update({"ok": True, "ran": True,
                               "num_steps": int(steps),
                               "frames": frames, "view_indices": targets})
            gif_out = cfg.get("gif_out")
            if gif_out:
                try:
                    from PIL import Image
                except ImportError:
                    result["ok"] = False
                    result["error"] = ("Pillow not installed (pip install "
                                       "pillow); PNG frames were written but "
                                       "the GIF was not assembled")
                else:
                    images = [Image.open(p).convert(
                        "P", palette=Image.Palette.ADAPTIVE) for p in frames]
                    images[0].save(gif_out, save_all=True,
                                   append_images=images[1:],
                                   duration=int(cfg.get("delay_ms", 40)),
                                   loop=0, disposal=2)
                    result["gif"] = gif_out
    finally:
        gmsh.finalize()
except Exception as exc:
    result["error"] = f"{type(exc).__name__}: {exc}"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(result, f)
"""


def _png_size(path: Path) -> list[int] | None:
    try:
        with open(path, "rb") as f:
            header = f.read(24)
    except OSError:
        return None
    if header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        return None
    return [int.from_bytes(header[16:20], "big"),
            int.from_bytes(header[20:24], "big")]


def _resolve_rotation(camera_preset: str | None,
                      rotation: list[float] | tuple[float, ...] | None
                      ) -> list[float] | None:
    if rotation is not None:
        rot = [float(x) for x in rotation]
        if len(rot) != 3:
            raise ValueError("rotation must contain three values")
        return rot
    if camera_preset is None:
        return None
    if camera_preset not in CAMERA_PRESETS:
        raise ValueError(
            f"unknown camera_preset: {camera_preset} "
            f"(available: {', '.join(sorted(CAMERA_PRESETS))})")
    return [float(x) for x in CAMERA_PRESETS[camera_preset]]


_INTERVAL_STYLES = {"iso": 1, "continuous": 2, "discrete": 3, "numeric": 4}
_SCALE_TYPES = {"linear": 1, "log": 2, "double_log": 3}
_VECTOR_TYPES = {"segment": 1, "arrow": 2, "pyramid": 3, "arrow3d": 4,
                 "displacement": 5, "comet": 6}
_GLYPH_LOCATIONS = {"cog": 1, "vertex": 2}
_CLIP_TARGETS = {"views": "View[0].Clip", "mesh": "Mesh.Clip",
                 "geometry": "Geometry.Clip"}
_AXES_MODES = {"none": 0, "box": 1, "frame": 2, "open": 3, "full": 4,
               "open_grid": 5}


def _pick(mapping: dict[str, int], value, what: str):
    """Named choice -> gmsh code.  Ints pass through (escape hatch for a
    value gmsh gained after this table was written); unknown names raise
    with the valid list rather than silently defaulting."""
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    key = str(value).lower()
    if key not in mapping:
        raise ValueError(
            f"unknown {what}: {value!r} (available: "
            f"{', '.join(sorted(mapping))}, or a raw gmsh code)")
    return float(mapping[key])


def _build_color(color: dict[str, Any] | None) -> dict[str, Any] | None:
    """Colour-scale control.  ``range`` is the important one: an explicit
    [min, max] (or "shared" across the views of THIS render) is what makes
    two figures comparable -- gmsh otherwise autoscales each view to its
    own extrema, so the same colour means different values per panel."""
    if not color:
        return None
    out: dict[str, Any] = {}
    rng = color.get("range")
    if rng is not None and rng != "auto":
        if rng == "shared":
            out["range"] = "shared"
        else:
            lo, hi = (float(v) for v in rng)
            if not hi > lo:
                raise ValueError(f"color range must be increasing, got {rng}")
            out["range"] = [lo, hi]
    if color.get("log"):
        out["scale_type"] = _SCALE_TYPES["log"]
    if color.get("scale") is not None:
        out["scale_type"] = _pick(_SCALE_TYPES, color["scale"], "color scale")
    if color.get("intervals") is not None:
        out["intervals"] = int(color["intervals"])
    if color.get("style") is not None:
        out["intervals_type"] = _pick(_INTERVAL_STYLES, color["style"],
                                      "color style")
    for key in ("colormap", "alpha"):
        if color.get(key) is not None:
            out[key] = float(color[key])
    for key in ("show_scale", "saturate", "invert"):
        if color.get(key) is not None:
            out[key] = 1.0 if color[key] else 0.0
    if color.get("format"):
        out["format"] = str(color["format"])
    if color.get("views") is not None:
        out["views"] = [int(i) for i in color["views"]]
    return out or None


def _build_glyphs(glyphs: dict[str, Any] | None) -> dict[str, Any] | None:
    """Vector-glyph control.  ``sampling=N`` draws every Nth element --
    the difference between a readable arrow field and a solid blue mat."""
    if not glyphs:
        return None
    out: dict[str, Any] = {}
    if glyphs.get("type") is not None:
        out["vector_type"] = _pick(_VECTOR_TYPES, glyphs["type"],
                                   "glyph type")
    if glyphs.get("location") is not None:
        out["location"] = _pick(_GLYPH_LOCATIONS, glyphs["location"],
                                "glyph location")
    if glyphs.get("sampling") is not None:
        n = int(glyphs["sampling"])
        if n < 1:
            raise ValueError(f"glyph sampling must be >= 1, got {n}")
        out["sampling"] = n
    for key in ("size_max", "size_min", "line_width", "point_size"):
        if glyphs.get(key) is not None:
            out[key] = float(glyphs[key])
    if glyphs.get("center") is not None:
        out["center"] = 1.0 if glyphs["center"] else 0.0
    if glyphs.get("views") is not None:
        out["views"] = [int(i) for i in glyphs["views"]]
    return out or None


def _build_clip(clip: list[dict[str, Any]] | None) -> list[dict[str, Any]]:
    """Up to 6 clipping planes, given as an outward normal + offset:
    the kept half-space is ``n . x + d >= 0``."""
    if not clip:
        return []
    if len(clip) > 6:
        raise ValueError(f"gmsh supports at most 6 clip planes, got {len(clip)}")
    out = []
    for k, spec in enumerate(clip):
        if "plane" in spec:
            plane = [float(v) for v in spec["plane"]]
            if len(plane) != 4:
                raise ValueError("clip plane must be [a, b, c, d]")
        else:
            n = [float(v) for v in spec["normal"]]
            if len(n) != 3:
                raise ValueError("clip normal must be [nx, ny, nz]")
            plane = n + [float(spec.get("offset", 0.0))]
        names = spec.get("apply_to") or ["views", "mesh", "geometry"]
        targets = []
        for name in names:
            if name not in _CLIP_TARGETS:
                raise ValueError(
                    f"unknown clip target {name!r} "
                    f"(available: {', '.join(sorted(_CLIP_TARGETS))})")
            targets.append(_CLIP_TARGETS[name])
        out.append({"plane": plane, "targets": targets,
                    "whole_elements": bool(spec.get("whole_elements"))})
    return out


def _build_axes(axes: dict[str, Any] | bool | None) -> dict[str, Any] | None:
    """Labelled axis box.  Publication figures need the axis LABELS to
    carry units; gmsh keeps label and format per axis."""
    if axes is None or axes is False:
        return None
    if axes is True:
        axes = {"mode": "box"}
    out: dict[str, Any] = {"mode": _pick(_AXES_MODES, axes.get("mode", "box"),
                                         "axes mode")}
    for key in ("labels", "format"):
        if axes.get(key) is not None:
            vals = list(axes[key])
            if len(vals) != 3:
                raise ValueError(f"axes {key} needs three entries (x, y, z)")
            out[key] = [str(v) for v in vals]
    if axes.get("tics") is not None:
        tics = list(axes["tics"])
        if len(tics) != 3:
            raise ValueError("axes tics needs three entries (x, y, z)")
        out["tics"] = [int(v) for v in tics]
    return out


def _build_annotations(annotations) -> list[dict[str, Any]]:
    """2D text overlays (window pixel coordinates; negative counts from
    the opposite edge, which is gmsh's own convention)."""
    if not annotations:
        return []
    out = []
    for a in annotations:
        if isinstance(a, str):
            a = {"text": a}
        if not a.get("text"):
            raise ValueError("annotation needs a 'text'")
        out.append({"text": str(a["text"]),
                    "x": float(a.get("x", 20.0)),
                    "y": float(a.get("y", 30.0)),
                    "align": str(a.get("align", "Left")),
                    "font": str(a.get("font", "Helvetica")),
                    "size": int(a.get("size", 18))})
    return out


def _figure_cfg(color, glyphs, clip, axes, annotations) -> dict[str, Any]:
    return {"color": _build_color(color),
            "glyphs": _build_glyphs(glyphs),
            "clip_planes": _build_clip(clip),
            "axes": _build_axes(axes),
            "annotations": _build_annotations(annotations)}


def _run_render(cfg: dict[str, Any], timeout_s: float) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="radia_mcp_gmsh_cfg_") as work:
        cfg_path = Path(work) / "render.json"
        cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
        return run_gmsh_json_subprocess(
            _RENDER_SCRIPT, [str(cfg_path)],
            timeout_s=timeout_s, prefix="radia_mcp_gmsh_render_")


def _merge_cut_plane(cut_plane: dict[str, Any] | None,
                     options: dict[str, float] | None
                     ) -> tuple[dict[str, float], bool]:
    """Fold a structured cut_plane into raw options (user options win)."""
    merged: dict[str, float] = {}
    clip_views = False
    if cut_plane is not None:
        clip_options = cut_plane_option_values(cut_plane)
        clip_views = clip_options.get("Mesh.Clip") == 1.0
        merged.update(clip_options)
    merged.update(options or {})
    return merged, clip_views


def _axis_equal_note(options: dict[str, float]) -> str | None:
    """Warn when a General.Scale* override breaks the 1:1:1 aspect.

    Policy (gmsh_usage topic "policy"): spatial figures -- contour /
    flux-line / streamline / section renders -- must be axis equal.
    Overriding a Scale axis is allowed only as a LABELED exaggeration.
    """
    scaled = {name: float(val) for name, val in options.items()
              if name in ("General.ScaleX", "General.ScaleY",
                          "General.ScaleZ") and float(val) != 1.0}
    if not scaled:
        return None
    parts = ", ".join(f"{k}={v:g}" for k, v in sorted(scaled.items()))
    return (f"axis-equal policy: {parts} overrides the 1:1:1 aspect -- "
            f"contour/flux-line/streamline/section figures must stay "
            f"axis equal; keep an exaggerated-scale figure only with "
            f"the factor stated in its caption")


_CAD_SUFFIXES = {".step", ".stp", ".brep", ".iges", ".igs"}


def _resolve_merge_files(merge_files) -> tuple[list[str], bool] | dict:
    """Validate overlay files; returns (paths, any_cad) or an error dict."""
    paths: list[str] = []
    any_cad = False
    for m in merge_files or []:
        p = Path(m)
        if not p.is_file():
            return {"ok": False, "error": f"merge file not found: {p}"}
        if p.suffix.lower() in _CAD_SUFFIXES:
            any_cad = True
        paths.append(str(p))
    return paths, any_cad


def render_png(path: str | Path,
               png_out: str | Path | None = None, *,
               width: int = 1000, height: int = 800,
               numsubedges: int = 4,
               camera_preset: str | None = None,
               rotation: list[float] | None = None,
               time_step: int | None = None,
               cut_plane: dict[str, Any] | None = None,
               merge_files: list[str | Path] | None = None,
               geometry_display: bool | None = None,
               color: dict[str, Any] | None = None,
               glyphs: dict[str, Any] | None = None,
               clip: list[dict[str, Any]] | None = None,
               axes: dict[str, Any] | bool | None = None,
               annotations: list[Any] | None = None,
               options: dict[str, float] | None = None,
               string_options: dict[str, str] | None = None,
               auto_mesh_display: bool = True,
               adapt_views: bool = True,
               smooth_normals: bool = True,
               timeout_s: float = 300.0) -> dict[str, Any]:
    """Render a .msh/.geo file to PNG in a gmsh subprocess.

    Figure controls (named values, no gmsh option strings needed):

    ``camera_preset`` -- ``"+x" "-x" "+y" "-y" "+z" "-z" "iso"``: the
    named axis points AT the camera, so ``"+y"`` shows the x-z plane
    face-on.  (Guessing raw ``rotation`` angles is how a plane ends up
    rendered edge-on as a single line.)

    ``color`` -- ``{"range": [lo, hi] | "shared", "log": bool,
    "intervals": n, "style": "continuous|iso|discrete|numeric",
    "format": "%.3g", "colormap": n, "alpha": a, "show_scale": bool,
    "saturate": bool, "views": [i]}``.  ``range`` is the one that
    matters for publication: gmsh autoscales EVERY view to its own
    extrema, so two panels of the same quantity are not comparable
    until they share a range.  ``"shared"`` unifies the views of this
    render; pass an explicit ``[lo, hi]`` to unify across files (get it
    from ``gmsh_field_stats``).

    ``glyphs`` -- ``{"type": "arrow3d", "sampling": n, "size_max": px,
    "center": bool, "location": "cog|vertex", "views": [i]}``.
    ``sampling`` draws every n-th element: the difference between a
    readable arrow field and a solid mat of arrowheads.

    ``clip`` -- up to 6 planes ``{"normal": [nx, ny, nz], "offset": d,
    "apply_to": ["views", "mesh", "geometry"], "whole_elements": bool}``;
    the kept half-space is ``n . x + d >= 0``.

    ``axes`` -- ``True`` or ``{"mode": "box|frame|open|full|open_grid",
    "labels": ["x [m]", ...], "format": ["%.3g", ...], "tics": [5,5,5]}``.

    ``annotations`` -- ``["text"]`` or ``[{"text":, "x":, "y":, "align":,
    "size":}]`` drawn in window pixels (negative counts from the far
    edge).

    Opening a ``.geo`` auto-loads its exact ``.geo.opt`` sidecar, so the
    launch artifact renders exactly as a user double-click would show it.
    ``cut_plane`` takes the same structured dict as the post-display
    contract ({enabled, normal, offset, ...}); views get Clip=1.
    ``smooth_normals`` (default on) averages surface normals per vertex
    -- the gmsh default leaves isosurfaces visibly faceted.

    ``merge_files`` overlays additional files into the same scene (the
    lab's Merge workflow: coil STEP + filament .msh + field .msh in one
    picture).  When a CAD file (.step/.brep) is merged,
    ``geometry_display`` defaults on: shaded solid faces, vertex dots
    off.  Radia/netgen STEP files carry meter coordinates and overlay
    the field data 1:1 (measured).
    """
    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    resolved = _resolve_merge_files(merge_files)
    if isinstance(resolved, dict):
        return resolved
    merge_paths, any_cad = resolved
    if geometry_display is None:
        geometry_display = any_cad
    out = Path(png_out) if png_out is not None else src.with_suffix(".png")
    out.parent.mkdir(parents=True, exist_ok=True)

    merged_options, clip_views = _merge_cut_plane(cut_plane, options)
    cfg = {
        "mode": "png",
        "path": str(src),
        "merge_files": merge_paths,
        "geometry_display": bool(geometry_display),
        "png_out": str(out),
        "width": int(width),
        "height": int(height),
        "numsubedges": int(numsubedges),
        "rotation": _resolve_rotation(camera_preset, rotation),
        "time_step": time_step,
        **_figure_cfg(color, glyphs, clip, axes, annotations),
        "options": merged_options,
        "string_options": string_options or {},
        "auto_mesh_display": bool(auto_mesh_display),
        "adapt_views": bool(adapt_views),
        "smooth_normals": bool(smooth_normals),
        "clip_views": clip_views,
    }
    result = _run_render(cfg, timeout_s)
    result["input"] = str(src)
    if result.get("ok"):
        size = _png_size(out)
        result["png_size"] = size
        notes = []
        scale_note = _axis_equal_note(merged_options)
        if scale_note:
            notes.append(scale_note)
        if size and size[0] < int(width):
            notes.append(
                f"exported width {size[0]} < requested {width}: the FLTK "
                f"sidebar consumes part of the window on Windows builds; "
                f"request a larger General.GraphicsWidth for exact sizes")
        blank = result.get("blank_check") or {}
        if blank.get("looks_blank"):
            notes.append(
                "image looks blank (content fraction "
                f"{blank.get('content_fraction')}): nothing was drawn -- "
                "check view visibility, mesh display options, and that the "
                "file actually contains elements/views")
        if notes:
            result["note"] = " | ".join(notes)
    return result


def export_animation(path: str | Path,
                     gif_out: str | Path | None = None, *,
                     keep_frames: bool = False,
                     view_indices: list[int] | None = None,
                     num_steps: int | None = None,
                     delay_ms: int = 40,
                     width: int = 1000, height: int = 800,
                     numsubedges: int = 4,
                     camera_preset: str | None = None,
                     rotation: list[float] | None = None,
                     cut_plane: dict[str, Any] | None = None,
                     merge_files: list[str | Path] | None = None,
                     geometry_display: bool | None = None,
                     color: dict[str, Any] | None = None,
                     glyphs: dict[str, Any] | None = None,
                     clip: list[dict[str, Any]] | None = None,
                     axes: dict[str, Any] | bool | None = None,
                     annotations: list[Any] | None = None,
                     options: dict[str, float] | None = None,
                     string_options: dict[str, str] | None = None,
                     adapt_views: bool = True,
                     smooth_normals: bool = True,
                     link_views: bool = True,
                     orbit_axis: str | None = None,
                     orbit_degrees: float = 360.0,
                     orbit_frames: int = 36,
                     time_step: int | None = None,
                     timeout_s: float = 900.0) -> dict[str, Any]:
    """Export a time-stepped post view animation as GIF (+ PNG frames).

    Steps every target view's ``TimeStep`` explicitly with linked views,
    writes one PNG per step, and assembles the GIF with Pillow inside
    the subprocess.  ``keep_frames=True`` retains the per-step PNGs in a
    ``<gif stem>_frames`` directory next to the GIF.

    ``orbit_axis`` switches to a CAMERA-ORBIT animation instead: the
    data stays at ``time_step`` while the camera sweeps
    ``orbit_degrees`` around the axis in ``orbit_frames`` frames (the
    ParaView camera-path fly-around; works on mesh-only files too).
    """
    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if orbit_axis is not None and orbit_axis not in ("x", "y", "z"):
        return {"ok": False,
                "error": f"orbit_axis must be x|y|z, got {orbit_axis!r}"}
    resolved = _resolve_merge_files(merge_files)
    if isinstance(resolved, dict):
        return resolved
    merge_paths, any_cad = resolved
    if geometry_display is None:
        geometry_display = any_cad
    gif = Path(gif_out) if gif_out is not None else src.with_suffix(".gif")
    gif.parent.mkdir(parents=True, exist_ok=True)

    if keep_frames:
        frames_dir = gif.parent / f"{gif.stem}_frames"
        frames_dir.mkdir(parents=True, exist_ok=True)
        temp_frames = None
    else:
        temp_frames = tempfile.mkdtemp(prefix="radia_mcp_gmsh_frames_")
        frames_dir = Path(temp_frames)

    merged_options, clip_views = _merge_cut_plane(cut_plane, options)
    cfg = {
        "mode": "animation",
        "path": str(src),
        "merge_files": merge_paths,
        "geometry_display": bool(geometry_display),
        "gif_out": str(gif),
        "frames_dir": str(frames_dir).replace("\\", "/"),
        "width": int(width),
        "height": int(height),
        "numsubedges": int(numsubedges),
        "rotation": _resolve_rotation(camera_preset, rotation),
        "view_indices": view_indices,
        "num_steps": num_steps,
        "delay_ms": int(delay_ms),
        **_figure_cfg(color, glyphs, clip, axes, annotations),
        "options": merged_options,
        "string_options": string_options or {},
        "auto_mesh_display": orbit_axis is not None,
        "adapt_views": bool(adapt_views),
        "smooth_normals": bool(smooth_normals),
        "link_views": bool(link_views),
        "clip_views": clip_views,
        "orbit": ({"axis": orbit_axis, "degrees": float(orbit_degrees),
                   "n_frames": max(2, int(orbit_frames))}
                  if orbit_axis is not None else None),
        "time_step": time_step,
    }
    try:
        result = _run_render(cfg, timeout_s)
    finally:
        if temp_frames is not None:
            shutil.rmtree(temp_frames, ignore_errors=True)

    result["input"] = str(src)
    if temp_frames is not None:
        result.pop("frames", None)  # already deleted with the temp dir
    else:
        result["frames_dir"] = str(frames_dir)
    if result.get("ok"):
        first_frame_size = None
        if keep_frames:
            frames = sorted(frames_dir.glob("frame_*.png"))
            if frames:
                first_frame_size = _png_size(frames[0])
        result["frame_size"] = first_frame_size
        result["gif_size_bytes"] = gif.stat().st_size if gif.is_file() else None
        scale_note = _axis_equal_note(merged_options)
        if scale_note:
            result["note"] = scale_note
    return result


def render_montage(images: list[str | Path],
                   out_png: str | Path, *,
                   cols: int | None = None,
                   labels: list[str] | None = None,
                   background: str = "white") -> dict[str, Any]:
    """Compose rendered PNGs into one comparison grid (side-by-side).

    The ParaView comparative-views analog for static output: before/
    after, per-frequency, or per-design renders lined up in a single
    image.  Cells take the largest input size; optional labels are
    drawn into the top-left corner of each cell.
    """
    if cols is not None:
        try:
            n_cols = int(cols)
        except (TypeError, ValueError):
            return {"ok": False,
                    "error": f"cols must be a positive integer, got {cols!r}"}
        if n_cols < 1:
            return {"ok": False,
                    "error": f"cols must be a positive integer, got {cols!r}"}
    else:
        n_cols = 0
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return {"ok": False,
                "error": "Pillow not installed (pip install pillow)"}
    paths = [Path(p) for p in images]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        return {"ok": False, "error": f"missing images: {missing}"}
    if len(paths) < 2:
        return {"ok": False, "error": "montage needs at least 2 images"}
    if labels is not None and len(labels) != len(paths):
        return {"ok": False,
                "error": f"got {len(labels)} labels for {len(paths)} images"}

    imgs = []
    for p in paths:
        with Image.open(p) as im:
            imgs.append(im.convert("RGB"))
    cell_w = max(im.width for im in imgs)
    cell_h = max(im.height for im in imgs)
    n = len(imgs)
    if not n_cols:
        n_cols = min(n, 3 if n != 4 else 2)
    n_rows = (n + n_cols - 1) // n_cols
    canvas = Image.new("RGB", (cell_w * n_cols, cell_h * n_rows),
                       background)
    draw = ImageDraw.Draw(canvas)
    for k, im in enumerate(imgs):
        r, c = divmod(k, n_cols)
        x = c * cell_w + (cell_w - im.width) // 2
        y = r * cell_h + (cell_h - im.height) // 2
        canvas.paste(im, (x, y))
        if labels is not None:
            draw.text((c * cell_w + 8, r * cell_h + 6), labels[k],
                      fill="black")
    out = Path(out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out)
    return {"ok": True, "png": str(out),
            "grid": [n_rows, n_cols], "cell_size": [cell_w, cell_h],
            "n_images": n}
