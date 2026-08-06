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
        gmsh.option.setNumber("Mesh.NumSubEdges", cfg["numsubedges"])

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

        if cfg.get("clip_views"):
            for i in range(n_views):
                gmsh.option.setNumber(f"View[{i}].Clip", 1)

        for name, val in (cfg.get("options") or {}).items():
            gmsh.option.setNumber(name, float(val))
        for name, val in (cfg.get("string_options") or {}).items():
            gmsh.option.setString(name, str(val))

        if cfg["mode"] == "png":
            ts = cfg.get("time_step")
            if ts is not None:
                for i in range(n_views):
                    gmsh.option.setNumber(f"View[{i}].TimeStep", int(ts))
            gmsh.fltk.initialize()
            gmsh.fltk.update()
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


def render_png(path: str | Path,
               png_out: str | Path | None = None, *,
               width: int = 1000, height: int = 800,
               numsubedges: int = 4,
               camera_preset: str | None = None,
               rotation: list[float] | None = None,
               time_step: int | None = None,
               cut_plane: dict[str, Any] | None = None,
               options: dict[str, float] | None = None,
               string_options: dict[str, str] | None = None,
               auto_mesh_display: bool = True,
               adapt_views: bool = True,
               timeout_s: float = 300.0) -> dict[str, Any]:
    """Render a .msh/.geo file to PNG in a gmsh subprocess.

    Opening a ``.geo`` auto-loads its exact ``.geo.opt`` sidecar, so the
    launch artifact renders exactly as a user double-click would show it.
    ``cut_plane`` takes the same structured dict as the post-display
    contract ({enabled, normal, offset, ...}); views get Clip=1.
    """
    src = Path(path)
    if not src.is_file():
        return {"ok": False, "error": f"file not found: {src}"}
    out = Path(png_out) if png_out is not None else src.with_suffix(".png")
    out.parent.mkdir(parents=True, exist_ok=True)

    merged_options, clip_views = _merge_cut_plane(cut_plane, options)
    cfg = {
        "mode": "png",
        "path": str(src),
        "png_out": str(out),
        "width": int(width),
        "height": int(height),
        "numsubedges": int(numsubedges),
        "rotation": _resolve_rotation(camera_preset, rotation),
        "time_step": time_step,
        "options": merged_options,
        "string_options": string_options or {},
        "auto_mesh_display": bool(auto_mesh_display),
        "adapt_views": bool(adapt_views),
        "clip_views": clip_views,
    }
    result = _run_render(cfg, timeout_s)
    result["input"] = str(src)
    if result.get("ok"):
        size = _png_size(out)
        result["png_size"] = size
        notes = []
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
                     options: dict[str, float] | None = None,
                     string_options: dict[str, str] | None = None,
                     adapt_views: bool = True,
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
        "gif_out": str(gif),
        "frames_dir": str(frames_dir).replace("\\", "/"),
        "width": int(width),
        "height": int(height),
        "numsubedges": int(numsubedges),
        "rotation": _resolve_rotation(camera_preset, rotation),
        "view_indices": view_indices,
        "num_steps": num_steps,
        "delay_ms": int(delay_ms),
        "options": merged_options,
        "string_options": string_options or {},
        "auto_mesh_display": orbit_axis is not None,
        "adapt_views": bool(adapt_views),
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
    n_cols = int(cols) if cols else min(n, 3 if n != 4 else 2)
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
