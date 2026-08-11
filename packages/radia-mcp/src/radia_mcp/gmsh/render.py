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
import math
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

        # ---- view selector (hide-others) ----------------------------
        # Resolved HERE, not in the wrapper: views contributed by
        # merge_files exist only after the merges above, so a
        # pre-read of the primary file (read_msh_data) cannot see them.
        sel = cfg.get("view_select")
        if sel is not None:
            names = [gmsh.option.getString(f"View[{i}].Name")
                     for i in range(n_views)]
            result["view_names"] = names
            chosen = []
            for want in sel:
                if isinstance(want, str):
                    hits = [i for i, nm in enumerate(names) if nm == want]
                    if not hits:
                        raise ValueError(
                            f"unknown view name {want!r} "
                            f"(available: {names})")
                    chosen.extend(hits)
                else:
                    if want < 0 or want >= n_views:
                        raise ValueError(
                            f"view index {want} outside [0, {n_views}) "
                            f"(available: {names})")
                    chosen.append(want)
            chosen = sorted(set(chosen))
            result["view_selected"] = chosen
            for i in range(n_views):
                # applied BEFORE the raw option passthrough below, so an
                # explicit options={"View[2].Visible": 1} still wins
                gmsh.option.setNumber(f"View[{i}].Visible",
                                      1.0 if i in chosen else 0.0)

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
            indices = [int(i) for i in spec]
            if any(i < 0 or i >= n_views for i in indices):
                raise IndexError(
                    f"view index outside [0, {n_views}): {indices}")
            return indices

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
                option_names = ([f"View[{i}].Clip" for i in range(n_views)]
                                if target == "__all_views__" else [target])
                for option_name in option_names:
                    cur = gmsh.option.getNumber(option_name)
                    gmsh.option.setNumber(option_name,
                                          float(int(cur) | mask))
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

        if cfg["mode"] == "volume":
            # Pseudo-volume rendering: stack semi-transparent cut planes.
            # gmsh has no ray-caster; CutPlane per offset + a
            # value-dependent alpha (ColormapAlphaPower) is the
            # classical substitute.
            src_tag = (gmsh.view.getTags()[int(cfg["view"])]
                       if isinstance(cfg["view"], int)
                       else next(t for t in gmsh.view.getTags()
                                 if gmsh.option.getString(
                                     f"View[{gmsh.view.getIndex(t)}].Name")
                                 == cfg["view"]))
            src_idx = gmsh.view.getIndex(src_tag)
            a, b, c = cfg["normal"]
            made = []
            for off in cfg["offsets"]:
                gmsh.plugin.setNumber("CutPlane", "A", a)
                gmsh.plugin.setNumber("CutPlane", "B", b)
                gmsh.plugin.setNumber("CutPlane", "C", c)
                gmsh.plugin.setNumber("CutPlane", "D", -off)
                gmsh.plugin.setNumber("CutPlane", "View", src_idx)
                made.append(gmsh.plugin.run("CutPlane"))
            # hide every pre-existing view, not just the source: a
            # second field view would keep drawing over the stack and
            # add a meaningless second colour bar
            for t in view_tags:
                gmsh.option.setNumber(
                    f"View[{gmsh.view.getIndex(t)}].Visible", 0)
            for t in made:
                i = gmsh.view.getIndex(t)
                # Plugin-created views do not reliably inherit display
                # options.  Copy the selected source's colour contract so
                # every slice uses one transfer function and one range.
                for opt in ("RangeType", "CustomMin", "CustomMax",
                            "ScaleType", "NbIso", "IntervalsType",
                            "ColormapNumber", "ColormapInvert",
                            "SaturateValues"):
                    gmsh.option.setNumber(
                        f"View[{i}].{opt}",
                        gmsh.option.getNumber(f"View[{src_idx}].{opt}"))
                gmsh.option.setString(
                    f"View[{i}].Format",
                    gmsh.option.getString(f"View[{src_idx}].Format"))
                gmsh.option.setNumber(f"View[{i}].ColormapAlpha",
                                      cfg["alpha"])
                gmsh.option.setNumber(f"View[{i}].ColormapAlphaPower",
                                      cfg["alpha_power"])
                gmsh.option.setNumber(f"View[{i}].ShowScale", 0)
                gmsh.option.setNumber(f"View[{i}].SmoothNormals", 1)
            if made:
                last = gmsh.view.getIndex(made[-1])
                gmsh.option.setNumber(f"View[{last}].ShowScale", 1)
            if cfg.get("slices_pos"):
                for t in made:
                    gmsh.view.write(t, cfg["slices_pos"], append=True)
            gmsh.fltk.initialize()
            gmsh.fltk.update()
            gmsh.write(cfg["png_out"])
            gmsh.fltk.finalize()
            result.update({"ok": True, "ran": True, "png": cfg["png_out"],
                           "n_slices": len(made),
                           "slices_pos": cfg.get("slices_pos")})
        elif cfg["mode"] == "png":
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
                    # a `view` selector doubles as the animation target:
                    # stepping a hidden view only costs time
                    targets = result.get("view_selected")
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
                if not frames:
                    # Pillow would die on images[0] with a bare
                    # IndexError; say what actually went wrong instead.
                    raise RuntimeError(
                        "no frames were written, so there is nothing to "
                        "assemble into a GIF")
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
    # Every mode sets ok=True BEFORE its last work item (GIF assembly,
    # the blank check, gmsh.finalize), so an exception raised after that
    # point used to be reported as a SUCCESS carrying an error string.
    # MEASURED on gmsh 4.15.2 before this line existed:
    # export_animation(num_steps=0) returned
    # {"ok": True, "error": "IndexError: list index out of range"} and
    # wrote no GIF.  ok is the caller's only gate -- clear it here.
    result["ok"] = False
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
_CLIP_TARGETS = {"views": "__all_views__", "mesh": "Mesh.Clip",
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
            if not isinstance(rng, (list, tuple)) or len(rng) != 2:
                raise ValueError("color range must be [min, max], 'shared', "
                                 "or 'auto'")
            lo, hi = (float(v) for v in rng)
            if not math.isfinite(lo) or not math.isfinite(hi) or not hi > lo:
                raise ValueError(f"color range must be increasing, got {rng}")
            out["range"] = [lo, hi]
    if color.get("log"):
        out["scale_type"] = _SCALE_TYPES["log"]
    if color.get("scale") is not None:
        out["scale_type"] = _pick(_SCALE_TYPES, color["scale"], "color scale")
    if color.get("intervals") is not None:
        intervals = int(color["intervals"])
        if intervals < 1:
            raise ValueError("color intervals must be >= 1")
        out["intervals"] = intervals
    if color.get("style") is not None:
        out["intervals_type"] = _pick(_INTERVAL_STYLES, color["style"],
                                      "color style")
    if color.get("colormap") is not None:
        out["colormap"] = float(color["colormap"])
    if color.get("alpha") is not None:
        alpha = float(color["alpha"])
        if not math.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
            raise ValueError("color alpha must be in [0, 1]")
        out["alpha"] = alpha
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
            if "normal" not in spec:
                raise ValueError("clip needs a 'plane' or 'normal'")
            n = [float(v) for v in spec["normal"]]
            if len(n) != 3:
                raise ValueError("clip normal must be [nx, ny, nz]")
            plane = n + [float(spec.get("offset", 0.0))]
        if any(not math.isfinite(v) for v in plane):
            raise ValueError("clip plane values must be finite")
        if not any(plane[i] for i in range(3)):
            raise ValueError("clip plane normal must be nonzero")
        names = spec.get("apply_to") or ["views", "mesh", "geometry"]
        if isinstance(names, str):
            names = [names]
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


def _build_view_select(view) -> list[Any] | None:
    """Normalize a ``view`` selector into a list of names / indices.

    Only the SHAPE is checked here.  Name and index resolution happens
    inside the render worker, after every ``merge_files`` entry has been
    merged: a merged file contributes views that a pre-read of the
    primary file cannot see, which is why ``render_panels``' pre-read
    (``read_msh_data``) selector could not be reused for this.
    """
    if view is None:
        return None
    items = list(view) if isinstance(view, (list, tuple)) else [view]
    if not items:
        raise ValueError(
            "view selector is empty; pass view=None to show every view")
    out: list[Any] = []
    for v in items:
        if isinstance(v, bool):
            raise TypeError(
                f"view must be a view name or a non-negative index, "
                f"got {v!r}")
        if isinstance(v, str):
            if not v:
                raise ValueError("view name must not be empty")
            out.append(v)
        elif isinstance(v, int):
            if v < 0:
                raise ValueError(f"view index must be >= 0, got {v}")
            out.append(int(v))
        else:
            raise TypeError(
                f"view must be a view name or a non-negative index, "
                f"got {v!r}")
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
               view: str | int | list[str | int] | None = None,
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

    ``view`` -- isolate one (or a few) post views by NAME or index:
    ``view="B_magnitude"``, ``view=2``, ``view=["B", "coil"]``.  Every
    other view is hidden, which is what a figure of one quantity needs
    -- a second visible view keeps drawing over the picture and adds a
    second colour bar.  It replaces the hand-written
    ``options={"View[0].Visible": 1, "View[1].Visible": 0, ...}`` block:
    an explicit ``View[i].Visible`` in ``options`` still overrides the
    selector.  An unknown name fails with the available names listed.

    ``color`` -- ``{"range": [lo, hi] | "shared", "log": bool,
    "intervals": n, "style": "continuous|iso|discrete|numeric",
    "format": "%.3g", "colormap": n, "alpha": a, "show_scale": bool,
    "saturate": bool, "views": [i]}``.  ``range`` is the one that
    matters for publication: gmsh autoscales EVERY view to its own
    extrema, so two panels of the same quantity are not comparable
    until they share a range.  ``"shared"`` unifies the views of this
    render; pass an explicit ``[lo, hi]`` to unify across files (get it
    from ``gmsh_field_range``).

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
        "view_select": _build_view_select(view),
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
                     view: str | int | list[str | int] | None = None,
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

    ``view`` isolates post views by NAME or index (``view="B"``,
    ``view=2``, ``view=["B", "coil"]``); every other view is hidden, and
    an explicit ``options={"View[i].Visible": ...}`` still wins.  When
    ``view_indices`` is not given it defaults to the views ``view``
    resolved to, so one argument replaces the usual
    ``view_indices=[4]`` + a wall of ``View[i].Visible`` entries.

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
    if num_steps is not None:
        # A non-positive count produced ZERO frames, and the empty frame
        # list then died inside GIF assembly -- reported (MEASURED) as
        # ok=True with "IndexError: list index out of range" and no GIF.
        # Reject it here, where the offending value is still in hand.
        if isinstance(num_steps, bool) or not isinstance(num_steps, int):
            return {"ok": False,
                    "error": f"num_steps must be an integer >= 1, or None "
                             f"to animate every time step, got "
                             f"{num_steps!r}"}
        if num_steps < 1:
            return {"ok": False,
                    "error": f"num_steps must be >= 1, or None to animate "
                             f"every time step in the target views, got "
                             f"{num_steps}"}
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
        "view_select": _build_view_select(view),
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


def _union_bbox(paths: list[Path]) -> tuple[list[float], list[float]]:
    """Union node bounding box across .msh inputs (pure-Python reader)."""
    from .msh_inspect import read_msh_data

    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    for p in paths:
        if p.suffix.lower() not in (".msh", ".pos"):
            continue
        for pt in read_msh_data(p)["nodes"].values():
            for i in range(3):
                lo[i] = min(lo[i], pt[i])
                hi[i] = max(hi[i], pt[i])
    if any(v == float("inf") for v in lo):
        raise ValueError("no .msh input carried nodes for the shared frame")
    return lo, hi


def _write_frame_geo(lo, hi, path: Path, pad: float = 0.02) -> Path:
    """8 corner points spanning the union box.

    gmsh refits the scene on every draw and ignores General.Min*/Max*
    and ZoomFactor (measured), so the ONLY way to give several renders
    one common zoom is to give them one common bounding box.  The points
    are hidden with Geometry.Points=0 and still count toward the fit.
    """
    span = [hi[i] - lo[i] for i in range(3)]
    diag = max(sum(s * s for s in span) ** 0.5, 1e-12)
    m = pad * diag
    lines = []
    n = 1
    for x in (lo[0] - m, hi[0] + m):
        for y in (lo[1] - m, hi[1] + m):
            for z in (lo[2] - m, hi[2] + m):
                lines.append(f"Point({n}) = {{{x:.9g}, {y:.9g}, {z:.9g}, 1}};")
                n += 1
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def render_panels(items: list[dict[str, Any] | str | Path],
                  out_png: str | Path, *,
                  cols: int | None = None,
                  camera_preset: str | None = "iso",
                  rotation: list[float] | None = None,
                  share_camera: bool = True,
                  share_color: bool = True,
                  view: str | int | None = None,
                  color: dict[str, Any] | None = None,
                  width: int = 620, height: int = 560,
                  background: str = "white",
                  keep_panels: bool = True,
                  **render_kwargs) -> dict[str, Any]:
    """Multi-panel figure with ONE camera, ONE zoom and ONE colour scale.

    ``render_montage`` pastes independently rendered images together:
    every panel gets its own auto-fit and its own autoscaled colour bar,
    so the panels look comparable while encoding different scales.  This
    renders the panels as a set:

    * ``share_camera`` gives every panel the same view AND the same zoom,
      by merging a hidden 8-point frame spanning the union bounding box
      (gmsh refits on each draw and ignores General.Min*/Max* and
      ZoomFactor -- measured -- so a common box is the mechanism).
    * ``share_color`` pins every panel to the union colour range
      (``field_range`` over the inputs), so one colour means one value
      across the whole figure.

    Args:
        items: at least two paths, or dicts
            ``{"path":, "label":, "merge_files":, ...}``
            whose extra keys override the shared render arguments.
        out_png: montage output path.
        cols: montage columns (default: one row).
        camera_preset / rotation: the shared view.
        share_camera: merge the common frame so the zoom matches too.
        share_color: compute and apply the union colour range.
        color: extra colour options merged into the shared range.
        keep_panels: keep the individual panel PNGs next to the montage.

    Returns:
        montage result plus ``shared_range``, ``panels`` and ``frame``.
    """
    specs: list[dict[str, Any]] = []
    for it in items:
        spec = dict(it) if isinstance(it, dict) else {"path": str(it)}
        if "path" not in spec:
            raise ValueError(f"panel item without a 'path': {spec}")
        specs.append(spec)
    if len(specs) < 2:
        raise ValueError("render_panels needs at least two items")
    if isinstance(view, bool):
        raise TypeError("view must be a name or a non-negative index")
    srcs = [Path(s["path"]) for s in specs]
    missing = [str(p) for p in srcs if not p.is_file()]
    if missing:
        return {"ok": False, "error": f"file(s) not found: {missing}"}

    out = Path(out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    panel_dir = out.with_name(out.stem + "_panels")
    panel_dir.mkdir(parents=True, exist_ok=True)

    effective_view = view
    shared_color = dict(color or {})
    range_info = None
    if share_color and "range" not in shared_color:
        from .msh_inspect import read_msh_data
        from .post_process import field_range

        msh = [p for p in srcs if p.suffix.lower() == ".msh"]
        if len(msh) != len(srcs):
            raise ValueError(
                "automatic share_color requires .msh inputs; pass an "
                "explicit color={'range': [...]}, or share_color=False")
        if msh:
            if effective_view is None and len(msh) > 1:
                # Sharing a range across DIFFERENT quantities produces a
                # colour bar that means nothing (T and A/m^2 on one
                # scale).  Refuse rather than emit a plausible figure.
                names = [{v["name"] for v in read_msh_data(p)["views"]}
                         for p in msh]
                common = set.intersection(*names) if names else set()
                if not common:
                    raise ValueError(
                        "share_color=True but the panels have no view in "
                        f"common ({[sorted(n) for n in names]}); pass "
                        "view=<name>, an explicit color={'range': [...]}, "
                        "or share_color=False for a mixed-quantity figure")
                if len(common) > 1:
                    raise ValueError(
                        "share_color=True but the panels have multiple views "
                        f"in common ({sorted(common)}); pass view=<name> so "
                        "one physical quantity owns the shared colour bar")
                effective_view = next(iter(common))
            range_info = field_range(msh, view=effective_view)
            if not range_info.get("ok"):
                raise ValueError(
                    f"could not compute a shared colour range: "
                    f"{range_info.get('error', range_info)}")
            shared_color["range"] = range_info["range"]

    frame = None
    if share_camera:
        try:
            lo, hi = _union_bbox(srcs)
        except ValueError as exc:
            raise ValueError(
                "share_camera=True requires inputs with readable mesh nodes; "
                "pass share_camera=False for unsupported formats") from exc
        frame = _write_frame_geo(lo, hi, panel_dir / "_shared_frame.geo")

    shared_merges = [str(m) for m in
                     (render_kwargs.pop("merge_files", None) or [])]
    shared_options = dict(render_kwargs.pop("options", None) or {})
    pngs: list[str] = []
    labels: list[str] = []
    for k, spec in enumerate(specs):
        label = spec.pop("label", Path(spec["path"]).stem)
        merges = shared_merges + [
            str(m) for m in (spec.pop("merge_files", None) or [])]
        panel_color = dict(shared_color)
        panel_color.update(spec.pop("color", None) or {})
        opts = dict(shared_options)
        opts.update(spec.pop("options", None) or {})
        if isinstance(effective_view, (str, int)):
            # one quantity per comparison figure: a second visible view
            # would add a second colour bar that the shared range does
            # not describe
            from .msh_inspect import read_msh_data

            if Path(spec["path"]).suffix.lower() == ".msh":
                for i, v in enumerate(read_msh_data(spec["path"])["views"]):
                    visible = (v["name"] == effective_view
                               if isinstance(effective_view, str)
                               else i == effective_view)
                    opts.setdefault(f"View[{i}].Visible",
                                    1 if visible else 0)
        if frame is not None:
            merges.append(str(frame))
            opts.setdefault("Geometry.Points", 0)
            opts.setdefault("Geometry.PointNumbers", 0)
            opts.setdefault("Geometry.Curves", 0)
        panel_args = {
            "width": width, "height": height,
            "camera_preset": camera_preset, "rotation": rotation,
        }
        panel_args.update(render_kwargs)
        panel_args.update({k2: v for k2, v in spec.items() if k2 != "path"})
        panel_args.update({"color": panel_color or None,
                           "merge_files": merges or None,
                           "options": opts or None})
        png = panel_dir / f"panel{k}_{Path(spec['path']).stem}.png"
        res = render_png(spec["path"], png, **panel_args)
        if not res.get("ok"):
            return {"ok": False, "error": f"panel {k} failed: {res}",
                    "panel": str(png)}
        pngs.append(str(png))
        labels.append(label)

    montage = render_montage(pngs, out, cols=cols or len(pngs),
                             labels=labels, background=background)
    montage.update({
        "panels": pngs if keep_panels else [],
        "shared_range": shared_color.get("range"),
        "range_info": range_info,
        "frame": str(frame) if frame else None,
        "shared_camera": bool(frame),
        "view": effective_view,
    })
    if not keep_panels:
        for p in pngs:
            Path(p).unlink(missing_ok=True)
    return montage


def volume_render(path: str | Path,
                  png_out: str | Path | None = None, *,
                  view: str | int = 0,
                  n_slices: int = 24,
                  axis: str = "z",
                  alpha: float = 0.35,
                  alpha_power: float = 2.0,
                  camera_preset: str | None = "iso",
                  rotation: list[float] | None = None,
                  color: dict[str, Any] | None = None,
                  width: int = 900, height: int = 800,
                  keep_slices: bool = False,
                  timeout_s: float = 900.0) -> dict[str, Any]:
    """Pseudo-volume rendering by compositing semi-transparent slices.

    gmsh has NO volume renderer -- there is no ray-caster and no 3D
    texture path.  This is the classical substitute and is named for
    what it does: ``n_slices`` cut planes perpendicular to ``axis`` are
    extracted and drawn semi-transparently, so the eye integrates them
    into a volume.  A value-dependent opacity
    (``ColormapAlphaPower``: alpha grows as value**power) acts as the
    transfer function, so low values fade out instead of fogging the
    picture.

    Honest limits vs a real volume renderer: the compositing is
    per-slice, not per-ray, so opacity depends on how many slices the
    eye line crosses; slices perpendicular to a strongly oblique view
    read as stripes (keep the axis roughly along the view direction);
    and it costs one CutPlane pass per slice.

    Args:
        path: ASCII MSH v4.x file holding a scalar field.
        png_out: output PNG (default: alongside the input).
        view: source view name or index.
        n_slices: number of cut planes (24 is a good default; >64 is slow).
        axis: "x" | "y" | "z" -- the stacking direction.
        alpha: base opacity of each slice.
        alpha_power: opacity exponent (0 = uniform, 2 = low values fade).
        keep_slices: also write the slice stack as a .pos file.
    """
    from .msh_inspect import read_msh_data

    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    if src.suffix.lower() != ".msh":
        return {"ok": False,
                "error": "volume_render requires an ASCII MSH v4.x file"}
    if axis not in ("x", "y", "z"):
        raise ValueError(f"axis must be x|y|z, got {axis!r}")
    if isinstance(n_slices, bool) or int(n_slices) != n_slices:
        raise ValueError(f"n_slices must be an integer, got {n_slices}")
    n_slices = int(n_slices)
    if n_slices < 2:
        raise ValueError(f"n_slices must be >= 2, got {n_slices}")
    alpha = float(alpha)
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    alpha_power = float(alpha_power)
    if not math.isfinite(alpha_power) or alpha_power < 0.0:
        raise ValueError(f"alpha_power must be finite and >= 0, got {alpha_power}")

    data = read_msh_data(src)
    if not data["views"]:
        return {"ok": False, "error": f"{src.name} carries no view"}
    if isinstance(view, bool):
        return {"ok": False, "error": "view must be a name or an index"}
    if isinstance(view, int):
        if view < 0 or view >= len(data["views"]):
            return {"ok": False, "error": f"view index {view} out of range"}
        chosen = data["views"][view]
    else:
        matches = [entry for entry in data["views"]
                   if entry["name"] == view]
        if len(matches) != 1:
            return {"ok": False,
                    "error": f"expected one view {view!r}, got {len(matches)}"}
        chosen = matches[0]
    if chosen["components"] != 1:
        return {"ok": False,
                "error": f"volume_render needs a scalar view; "
                         f"{chosen['name']!r} has {chosen['components']} "
                         "components"}
    pts = list(data["nodes"].values())
    if not pts:
        return {"ok": False, "error": f"{src.name} holds no nodes"}
    k = "xyz".index(axis)
    lo = min(p[k] for p in pts)
    hi = max(p[k] for p in pts)
    if hi <= lo:
        return {"ok": False,
                "error": f"{src.name} is degenerate along {axis}"}
    # interior offsets: the bounding faces carry no interior information
    offsets = [lo + (hi - lo) * (i + 0.5) / n_slices for i in range(n_slices)]

    out = (Path(png_out) if png_out is not None
           else src.with_name(src.stem + f"_volume_{axis}.png"))
    out.parent.mkdir(parents=True, exist_ok=True)
    slices_pos = (str(out.with_suffix(".pos")) if keep_slices else None)

    cfg = {
        "mode": "volume",
        "path": str(src),
        "png_out": str(out),
        "view": view,
        "normal": [1.0 if i == k else 0.0 for i in range(3)],
        "offsets": offsets,
        "alpha": float(alpha),
        "alpha_power": float(alpha_power),
        "slices_pos": slices_pos,
        "width": int(width),
        "height": int(height),
        "numsubedges": 4,
        "rotation": _resolve_rotation(camera_preset, rotation),
        **_figure_cfg(color, None, None, None, None),
        "options": {},
        "string_options": {},
        "auto_mesh_display": False,
        "adapt_views": False,
        "smooth_normals": True,
        "geometry_display": False,
        "merge_files": [],
    }
    result = _run_render(cfg, timeout_s)
    result.setdefault("n_slices", n_slices)
    result.setdefault("axis", axis)
    result.setdefault("method",
                      "slice-stack compositing; NOT ray-cast volume "
                      "rendering (gmsh has no volume renderer)")
    return result
