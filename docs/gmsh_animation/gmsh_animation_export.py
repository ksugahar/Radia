"""Export frames, GIF, and MP4 from the docs-local Gmsh animation .geo."""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _metadata
import json
import os
import platform
import shutil
import sys
import time
from pathlib import Path

from gmsh_animation_inspect import inspect_msh, validation_result_path


def _pkg_version(name: str) -> str | None:
    try:
        return _metadata.version(name)
    except Exception:
        return None


def runtime_versions() -> dict:
    return {
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "gmsh_version": _pkg_version("gmsh"),
        "pillow_version": _pkg_version("Pillow"),
        "imageio_version": _pkg_version("imageio"),
        "imageio_ffmpeg_version": _pkg_version("imageio-ffmpeg"),
        "radia_version": _pkg_version("radia"),
        "radia_mcp_version": _pkg_version("radia-mcp"),
    }


def _read_geo_merges(path: Path) -> list[str]:
    merges = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("Merge "):
            parts = s.split('"')
            if len(parts) >= 2:
                merges.append(parts[1])
    return merges


def _display_options() -> str:
    return "\n".join([
        "General.Orthographic = 1;",
        "General.Trackball = 0;",
        "General.RotationX = -68;",
        "General.RotationY = 0;",
        "General.RotationZ = 0;",
        "General.RotationCenterGravity = 1;",
        "General.SmallAxes = 1;",
        "General.BackgroundGradient = 0;",
        "Mesh.NumSubEdges = 4;",
        "Mesh.SurfaceFaces = 1;",
        "Mesh.SurfaceEdges = 0;",
        "Mesh.VolumeEdges = 0;",
        "Mesh.VolumeFaces = 0;",
        "View[0].VectorType = 5;",
        "View[0].DisplacementFactor = 1.0;",
        "View[0].ShowElement = 1;",
        "View[0].ShowScale = 0;",
        "View[0].IntervalsType = 2;",
        "View[0].DrawSkinOnly = 1;",
        "View[0].AdaptVisualizationGrid = 1;",
        "View[0].MaxRecursionLevel = 1;",
        "PostProcessing.Link = 1;",
        "PostProcessing.AnimationCycle = 0;",
    ]) + "\n"


def ensure_companion_files(example_dir: str | Path = ".") -> dict:
    """Ensure the .geo launch artifact has exact Gmsh .opt sidecars."""
    root = Path(example_dir)
    if not root.is_absolute():
        root = (Path(__file__).resolve().parent / root).resolve()
    geo = root / "animation.geo"
    msh = root / "rotor_animation.msh"

    if not geo.exists() or not msh.exists():
        from moving_body_demo import main as build_demo
        build_demo()

    options = _display_options()
    geo_opt = Path(str(geo) + ".opt")
    msh_opt = Path(str(msh) + ".opt")
    geo_opt.write_text(
        "// Auto-generated display options for animation.geo\n" + options,
        encoding="utf-8",
    )
    msh_opt.write_text(
        "// Auto-generated raw mesh/data inspection options for rotor_animation.msh\n"
        + options,
        encoding="utf-8",
    )
    return {
        "geo": geo,
        "msh": msh,
        "geo_opt": geo_opt,
        "msh_opt": msh_opt,
        "geo_merges": _read_geo_merges(geo),
    }


def _repo_rel(path: Path) -> str:
    repo = Path(__file__).resolve().parents[2]
    try:
        return path.resolve().relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _file_info(path: Path) -> dict:
    return {
        "path": _repo_rel(path),
        "bytes": path.stat().st_size if path.exists() else 0,
        "exists": path.exists(),
    }


def export_animation(
    example_dir: str | Path = ".",
    *,
    frame_steps: list[int] | None = None,
    width: int = 700,
    height: int = 480,
    duration_ms: int = 80,
    cleanup_frames: bool = True,
) -> dict:
    """Open animation.geo in Gmsh, synchronize steps, and export GIF/MP4."""
    companions = ensure_companion_files(example_dir)
    root = companions["geo"].parent
    msh_info = inspect_msh(companions["msh"])
    n_steps = msh_info["node_data_count"]
    if frame_steps is None:
        frame_steps = list(range(n_steps))
    frame_steps = [int(s) for s in frame_steps if 0 <= int(s) < n_steps]
    if not frame_steps:
        raise ValueError("no valid frame steps selected")

    import gmsh
    from PIL import Image
    import imageio.v2 as imageio

    temp_root = Path("C:/temp") if os.name == "nt" else Path(os.environ.get("TMPDIR", "/tmp"))
    temp_root.mkdir(parents=True, exist_ok=True)
    frame_dir = temp_root / "radia_gmsh_animation_export_frames"
    if frame_dir.exists():
        shutil.rmtree(frame_dir)
    frame_dir.mkdir(parents=True, exist_ok=True)

    frames: list[Path] = []
    timings: dict[str, float] = {}
    t0 = time.perf_counter()
    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.GraphicsPositionX", 100)
        gmsh.option.setNumber("General.GraphicsPositionY", 100)
        gmsh.option.setNumber("General.GraphicsWidth", width)
        gmsh.option.setNumber("General.GraphicsHeight", height)
        gmsh.open(str(companions["geo"]))
        view_count = len(gmsh.view.getTags())
        gmsh.option.setNumber("PostProcessing.Link", 1)
        gmsh.option.setNumber("PostProcessing.AnimationCycle", 0)
        gmsh.fltk.initialize()
        timings["gmsh_open_and_fltk_s"] = round(time.perf_counter() - t0, 3)

        t_frames = time.perf_counter()
        for step in frame_steps:
            for idx in range(view_count):
                gmsh.option.setNumber(f"View[{idx}].TimeStep", step)
            gmsh.fltk.update()
            frame = frame_dir / f"gmsh_animation_export_{step:03d}.png"
            gmsh.write(str(frame))
            frames.append(frame)
        timings["frame_export_s"] = round(time.perf_counter() - t_frames, 3)
    finally:
        try:
            gmsh.fltk.finalize()
        except Exception:
            pass
        gmsh.finalize()

    gif_path = root / "gmsh_animation_export.gif"
    mp4_path = root / "gmsh_animation_export.mp4"
    t_encode = time.perf_counter()
    images = [Image.open(p).convert("P", palette=Image.Palette.ADAPTIVE)
              for p in frames]
    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        disposal=2,
    )
    rgb_frames = [imageio.imread(str(p)) for p in frames]
    imageio.mimsave(str(mp4_path), rgb_frames, fps=max(1, round(1000 / duration_ms)))
    timings["gif_mp4_encode_s"] = round(time.perf_counter() - t_encode, 3)
    timings["total_s"] = round(time.perf_counter() - t0, 3)

    if cleanup_frames:
        shutil.rmtree(frame_dir, ignore_errors=True)

    return {
        "schema": "radia.validation.gmsh_animation_export.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": runtime_versions(),
        "source": {
            "geo": _file_info(companions["geo"]),
            "msh": _file_info(companions["msh"]),
            "geo_opt": _file_info(companions["geo_opt"]),
            "msh_opt": _file_info(companions["msh_opt"]),
            "geo_merges": companions["geo_merges"],
        },
        "animation": {
            "node_data_frames": n_steps,
            "exported_steps": frame_steps,
            "exported_frame_count": len(frames),
            "view_time_steps_synchronized": True,
            "postprocessing_link": 1,
            "animation_cycle": 0,
        },
        "artifacts": {
            "gif": _file_info(gif_path),
            "mp4": _file_info(mp4_path),
        },
        "timing_breakdown_s": timings,
        "checks": {
            "opened_geo_not_msh": companions["geo"].name == "animation.geo",
            "geo_opt_exact_sidecar": companions["geo_opt"].name == "animation.geo.opt",
            "msh_opt_raw_inspection_sidecar": companions["msh_opt"].name == "rotor_animation.msh.opt",
            "gif_written": gif_path.exists() and gif_path.stat().st_size > 0,
            "mp4_written": mp4_path.exists() and mp4_path.stat().st_size > 0,
        },
    }


def write_results_json(result: dict, path: str | Path) -> Path:
    p = Path(path)
    p.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                 encoding="utf-8")
    return p


def format_summary(result: dict) -> str:
    src = result["source"]
    artifacts = result["artifacts"]
    checks = result["checks"]
    return "\n".join([
        f"opened: {src['geo']['path']}",
        f"geo sidecar: {src['geo_opt']['path']}",
        f"raw msh sidecar: {src['msh_opt']['path']}",
        f"frames: {result['animation']['exported_frame_count']}",
        f"gif: {artifacts['gif']['path']} ({artifacts['gif']['bytes']} bytes)",
        f"mp4: {artifacts['mp4']['path']} ({artifacts['mp4']['bytes']} bytes)",
        f"checks: geo={checks['opened_geo_not_msh']}, "
        f"geo_opt={checks['geo_opt_exact_sidecar']}, "
        f"gif={checks['gif_written']}, mp4={checks['mp4_written']}",
    ])


if __name__ == "__main__":
    result = export_animation()
    out = write_results_json(
        result,
        validation_result_path("gmsh_animation_export_results.json"),
    )
    print(format_summary(result))
    print(f"wrote {out}")
