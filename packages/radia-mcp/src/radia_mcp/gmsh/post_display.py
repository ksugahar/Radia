"""Gmsh post-processing launch artifact helpers.

These helpers keep Radia and MATLAB/Gypsilab acoustic artifacts on the same
display contract: the raw data live in a single Gmsh MSH v4.1 file, while the
human launch target is a sibling ``case.geo`` plus exact sidecars
``case.geo.opt`` and ``case.msh.opt``.  The contract is intentionally about
post-processing only; it does not create geometry or generate meshes.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


CAMERA_PRESETS: dict[str, tuple[float, float, float]] = {
    "z_up_xz_from_positive_y": (-68.0, 0.0, 0.0),
    "positive_y_oblique": (75.0, 0.0, -20.0),
    "front_xz": (-90.0, 0.0, 0.0),
    "custom": (0.0, 0.0, 0.0),
}


def build_gmsh_post_display_contract(
    msh_path: str | Path,
    *,
    output_base: str | Path | None = None,
    title: str = "Gmsh post launch artifact",
    camera_preset: str = "z_up_xz_from_positive_y",
    rotation: tuple[float, float, float] | list[float] | None = None,
    scale: tuple[float, float, float] | list[float] = (1.25, 1.25, 1.25),
    cut_plane: dict[str, Any] | None = None,
    views: list[dict[str, Any]] | None = None,
    mesh: dict[str, Any] | None = None,
    animation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a public-safe manifest for a Gmsh .geo/.opt launch artifact."""

    msh = Path(msh_path)
    if output_base is None:
        base = msh.with_suffix("")
    else:
        base = Path(output_base)
    geo = Path(str(base) + ".geo")
    geo_opt = Path(str(geo) + ".opt")
    opt = Path(str(base) + ".opt")
    msh_opt = msh.with_name(msh.name + ".opt")
    display_json = Path(str(base) + ".display.json")

    if camera_preset not in CAMERA_PRESETS:
        raise ValueError(f"unknown camera_preset: {camera_preset}")
    rot = tuple(float(x) for x in (rotation if rotation is not None else CAMERA_PRESETS[camera_preset]))
    if len(rot) != 3:
        raise ValueError("rotation must contain three values")
    scale_tuple = tuple(float(x) for x in scale)
    if len(scale_tuple) != 3 or any(x <= 0 for x in scale_tuple):
        raise ValueError("scale must contain three positive values")

    norm_views = [_normalize_view(i, view) for i, view in enumerate(views or [])]
    norm_cut = _normalize_cut_plane(cut_plane)
    norm_mesh = {
        "surface_faces": False,
        "surface_edges": False,
        "volume_faces": False,
        "volume_edges": False,
        "num_sub_edges": 4,
    }
    norm_mesh.update(mesh or {})
    norm_animation = {
        "delay": 0.04,
        "cycle": False,
        "step": 1,
        "link_time_steps": True,
    }
    norm_animation.update(animation or {})

    manifest = {
        "schema": "cae-ai-lab.gmsh-post-launch.v1",
        "title": str(title),
        "gmsh_msh": str(msh),
        "gmsh_geo": str(geo),
        "gmsh_geo_opt": str(geo_opt),
        "gmsh_opt": str(opt),
        "gmsh_msh_opt": str(msh_opt),
        "display_json": str(display_json),
        "launch_target": str(geo),
        "mesh_inspection_target": str(msh),
        "gmsh_msh_version": "4.1",
        "camera": {
            "preset": camera_preset,
            "rotation": list(rot),
            "scale": list(scale_tuple),
            "axis_up": "z",
        },
        "cut_plane": norm_cut,
        "mesh": norm_mesh,
        "animation": norm_animation,
        "views": norm_views,
    }
    manifest["checks"] = gmsh_post_display_manifest_gate(manifest)["checks"]
    manifest["status"] = "ok" if all(manifest["checks"].values()) else "needs_attention"
    manifest["pass"] = manifest["status"] == "ok"
    return manifest


def write_gmsh_post_launch_artifact(
    msh_path: str | Path,
    *,
    output_base: str | Path | None = None,
    title: str = "Gmsh post launch artifact",
    camera_preset: str = "z_up_xz_from_positive_y",
    rotation: tuple[float, float, float] | list[float] | None = None,
    scale: tuple[float, float, float] | list[float] = (1.25, 1.25, 1.25),
    cut_plane: dict[str, Any] | None = None,
    views: list[dict[str, Any]] | None = None,
    mesh: dict[str, Any] | None = None,
    animation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write ``case.geo``, ``case.geo.opt``, ``case.msh.opt`` and JSON."""

    manifest = build_gmsh_post_display_contract(
        msh_path,
        output_base=output_base,
        title=title,
        camera_preset=camera_preset,
        rotation=rotation,
        scale=scale,
        cut_plane=cut_plane,
        views=views,
        mesh=mesh,
        animation=animation,
    )
    geo = Path(manifest["gmsh_geo"])
    geo_opt = Path(manifest["gmsh_geo_opt"])
    opt = Path(manifest["gmsh_opt"])
    msh_opt = Path(manifest["gmsh_msh_opt"])
    display_json = Path(manifest["display_json"])
    for path in (geo, geo_opt, opt, msh_opt, display_json):
        path.parent.mkdir(parents=True, exist_ok=True)

    merge_target = _merge_target(Path(manifest["gmsh_msh"]), geo)
    body = _options_body(manifest)
    geo.write_text(
        f"// {manifest['title']}\n"
        "// Open this .geo for post-processing; open .msh for raw mesh inspection.\n"
        f'Merge "{merge_target}";\n'
        f"{body}",
        encoding="utf-8",
    )
    geo_opt.write_text("// Exact .geo.opt post-display sidecar.\n" + body, encoding="utf-8")
    opt.write_text(geo_opt.read_text(encoding="utf-8"), encoding="utf-8")
    msh_opt.write_text(_mesh_options_body(manifest), encoding="utf-8")

    manifest["checks"] = gmsh_post_display_manifest_gate(manifest, require_files=True)["checks"]
    manifest["status"] = "ok" if all(manifest["checks"].values()) else "needs_attention"
    manifest["pass"] = manifest["status"] == "ok"
    display_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def gmsh_post_display_manifest_gate(
    manifest: dict[str, Any],
    *,
    require_files: bool = False,
) -> dict[str, Any]:
    """Validate the shared Radia/Gypsilab Gmsh post-display contract."""

    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a dictionary")

    camera = manifest.get("camera") or {}
    cut_plane = manifest.get("cut_plane") or {}
    views = manifest.get("views") or []
    checks = {
        "schema_recorded": manifest.get("schema") == "cae-ai-lab.gmsh-post-launch.v1",
        "msh_version_is_41": str(manifest.get("gmsh_msh_version", "")) == "4.1",
        "launch_target_is_geo": str(manifest.get("launch_target", "")).endswith(".geo"),
        "geo_opt_exact_autoload": str(manifest.get("gmsh_geo_opt", "")).endswith(".geo.opt"),
        "msh_opt_exact_autoload": str(manifest.get("gmsh_msh_opt", "")).endswith(".msh.opt"),
        "plain_opt_only_compatibility": str(manifest.get("gmsh_opt", "")).endswith(".opt")
        and not str(manifest.get("gmsh_opt", "")).endswith((".geo.opt", ".msh.opt")),
        "camera_z_up_recorded": camera.get("axis_up") == "z",
        "camera_rotation_recorded": _has_three_numbers(camera.get("rotation")),
        "view_metadata_recorded": isinstance(views, list) and len(views) > 0,
        "views_have_names": isinstance(views, list) and all(str(v.get("name", "")).strip() for v in views),
        "cut_plane_metadata_recorded": (not cut_plane.get("enabled"))
        or (_has_three_numbers(cut_plane.get("normal")) and "offset" in cut_plane),
    }
    if require_files:
        checks.update(
            {
                "geo_file_exists": Path(str(manifest.get("gmsh_geo", ""))).is_file(),
                "geo_opt_file_exists": Path(str(manifest.get("gmsh_geo_opt", ""))).is_file(),
                "msh_opt_file_exists": Path(str(manifest.get("gmsh_msh_opt", ""))).is_file(),
                "display_json_path_recorded": str(manifest.get("display_json", "")).endswith(".display.json"),
            }
        )
    return {
        "schema": "cae-ai-lab.gmsh-post-launch-gate.v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
    }


def _normalize_cut_plane(cut_plane: dict[str, Any] | None) -> dict[str, Any]:
    base = {
        "enabled": False,
        "normal": [0.0, -1.0, 0.0],
        "offset": 0.0,
        "whole_elements": False,
        "only_volume": False,
        "only_draw_intersecting_volume": False,
    }
    if cut_plane:
        base.update(cut_plane)
    if base["enabled"] and not _has_three_numbers(base.get("normal")):
        raise ValueError("enabled cut_plane requires a three-number normal")
    base["normal"] = [float(x) for x in base["normal"]]
    base["offset"] = float(base["offset"])
    return base


def _normalize_view(index: int, view: dict[str, Any]) -> dict[str, Any]:
    base = {
        "index": index,
        "name": f"view_{index}",
        "kind": "scalar",
        "visible": True,
        "time_step": 0,
        "intervals_type": 2,
        "nb_iso": 32,
        "colormap_number": 2,
        "show_element": False,
        "show_scale": True,
        "light": False,
    }
    base.update(view)
    if not str(base["name"]).strip():
        raise ValueError("view name must be nonempty")
    return base


def _options_body(manifest: dict[str, Any]) -> str:
    camera = manifest["camera"]
    mesh = manifest["mesh"]
    cut_plane = manifest["cut_plane"]
    animation = manifest["animation"]
    lines = [
        "General.InitialModule = 5;",
        "General.Trackball = 0;",
        f"General.RotationX = {camera['rotation'][0]:.17g};",
        f"General.RotationY = {camera['rotation'][1]:.17g};",
        f"General.RotationZ = {camera['rotation'][2]:.17g};",
        "General.RotationCenterGravity = 1;",
        f"General.ScaleX = {camera['scale'][0]:.17g};",
        f"General.ScaleY = {camera['scale'][1]:.17g};",
        f"General.ScaleZ = {camera['scale'][2]:.17g};",
        "General.Color.Background = {255,255,255};",
        "General.Color.Foreground = {0,0,0};",
        "General.Color.Text = {0,0,0};",
        f"Mesh.NumSubEdges = {int(mesh['num_sub_edges'])};",
        f"Mesh.SurfaceFaces = {_bool(mesh['surface_faces'])};",
        f"Mesh.SurfaceEdges = {_bool(mesh['surface_edges'])};",
        f"Mesh.VolumeFaces = {_bool(mesh['volume_faces'])};",
        f"Mesh.VolumeEdges = {_bool(mesh['volume_edges'])};",
        "Mesh.Light = 1;",
        "Mesh.LightTwoSide = 1;",
    ]
    lines.extend(_cut_plane_lines(cut_plane))
    lines.append(f"PostProcessing.Link = {_bool(animation['link_time_steps'])};")
    for view in manifest["views"]:
        lines.extend(_view_lines(view, bool(cut_plane.get("enabled"))))
    lines.extend(
        [
            f"PostProcessing.AnimationDelay = {float(animation['delay']):.17g};",
            f"PostProcessing.AnimationCycle = {_bool(animation['cycle'])};",
            f"PostProcessing.AnimationStep = {int(animation['step'])};",
        ]
    )
    return "\n".join(lines) + "\n"


def _mesh_options_body(manifest: dict[str, Any]) -> str:
    camera = manifest["camera"]
    max_view = max((int(view.get("index", i)) for i, view in enumerate(manifest["views"])), default=-1)
    lines = [
        "// Mesh/data inspection sidecar. Open .geo for the post display.",
        "General.InitialModule = 1;",
        "General.Trackball = 0;",
        f"General.RotationX = {camera['rotation'][0]:.17g};",
        f"General.RotationY = {camera['rotation'][1]:.17g};",
        f"General.RotationZ = {camera['rotation'][2]:.17g};",
        f"Mesh.NumSubEdges = {int(manifest['mesh']['num_sub_edges'])};",
        "Mesh.SurfaceFaces = 1;",
        "Mesh.SurfaceEdges = 1;",
        "Mesh.VolumeFaces = 0;",
        "Mesh.VolumeEdges = 0;",
        "Mesh.ColorCarousel = 2;",
    ]
    for idx in range(max_view + 1):
        lines.append(f"View[{idx}].Visible = 0;")
    return "\n".join(lines) + "\n"


def cut_plane_option_values(cut_plane: dict[str, Any] | None) -> dict[str, float]:
    """Normalized cut-plane state as {gmsh option name: numeric value}.

    Shared with the renderer: the same cut_plane dict accepted by the
    launch-artifact contract can be applied through option.setNumber.
    """
    cut = _normalize_cut_plane(cut_plane)
    if not cut["enabled"]:
        return {"Mesh.Clip": 0.0}
    normal = cut["normal"]
    norm = sum(float(x) ** 2 for x in normal) ** 0.5
    if norm <= 0:
        raise ValueError("cut plane normal must be nonzero")
    n = [float(x) / norm for x in normal]
    return {
        "General.Clip0A": n[0],
        "General.Clip0B": n[1],
        "General.Clip0C": n[2],
        "General.Clip0D": float(cut["offset"]),
        "General.ClipFactor": 5.0,
        "General.ClipWholeElements": float(_bool(cut["whole_elements"])),
        "General.ClipOnlyVolume": float(_bool(cut["only_volume"])),
        "General.ClipOnlyDrawIntersectingVolume": float(
            _bool(cut["only_draw_intersecting_volume"])),
        "Mesh.Clip": 1.0,
    }


def _cut_plane_lines(cut_plane: dict[str, Any]) -> list[str]:
    if not cut_plane.get("enabled"):
        return ["Mesh.Clip = 0;"]
    normal = cut_plane["normal"]
    norm = sum(float(x) ** 2 for x in normal) ** 0.5
    if norm <= 0:
        raise ValueError("cut plane normal must be nonzero")
    n = [float(x) / norm for x in normal]
    return [
        f"General.Clip0A = {n[0]:.17g};",
        f"General.Clip0B = {n[1]:.17g};",
        f"General.Clip0C = {n[2]:.17g};",
        f"General.Clip0D = {float(cut_plane['offset']):.17g};",
        "General.ClipFactor = 5;",
        f"General.ClipWholeElements = {_bool(cut_plane['whole_elements'])};",
        f"General.ClipOnlyVolume = {_bool(cut_plane['only_volume'])};",
        f"General.ClipOnlyDrawIntersectingVolume = {_bool(cut_plane['only_draw_intersecting_volume'])};",
        "Mesh.Clip = 1;",
    ]


def _view_lines(view: dict[str, Any], cut_plane_enabled: bool) -> list[str]:
    idx = int(view.get("index", 0))
    lines = [
        f"View[{idx}].Visible = {_bool(view.get('visible', True))};",
        f'View[{idx}].Name = "{_escape(str(view["name"]))}";',
        f"View[{idx}].TimeStep = {int(view.get('time_step', 0))};",
        f"View[{idx}].IntervalsType = {int(view.get('intervals_type', 2))};",
        f"View[{idx}].NbIso = {int(view.get('nb_iso', 32))};",
        f"View[{idx}].ColormapNumber = {int(view.get('colormap_number', 2))};",
        f"View[{idx}].ShowElement = {_bool(view.get('show_element', False))};",
        f"View[{idx}].ShowScale = {_bool(view.get('show_scale', True))};",
        f"View[{idx}].Light = {_bool(view.get('light', False))};",
        f"View[{idx}].Clip = {_bool(view.get('clip', cut_plane_enabled))};",
    ]
    if "range" in view and view["range"] is not None:
        lo, hi = view["range"]
        lines.extend(
            [
                f"View[{idx}].RangeType = 2;",
                f"View[{idx}].CustomMin = {float(lo):.17g};",
                f"View[{idx}].CustomMax = {float(hi):.17g};",
            ]
        )
    kind = str(view.get("kind", "scalar")).lower()
    if kind == "vector":
        lines.append(f"View[{idx}].VectorType = {int(view.get('vector_type', 4))};")
        pixels = float(view.get("arrow_size_pixels", 0))
        if pixels > 0:
            lines.append(f"View[{idx}].ArrowSizeMin = {pixels:.17g};")
            lines.append(f"View[{idx}].ArrowSizeMax = {pixels:.17g};")
    elif kind == "displacement":
        lines.extend(
            [
                f"View[{idx}].VectorType = 5;",
                f"View[{idx}].DisplacementFactor = {float(view.get('displacement_factor', 1.0)):.17g};",
                f"View[{idx}].DrawSkinOnly = {_bool(view.get('draw_skin_only', True))};",
                f"View[{idx}].SmoothNormals = {_bool(view.get('smooth_normals', True))};",
            ]
        )
    return lines


def _merge_target(msh: Path, geo: Path) -> str:
    if msh.parent == geo.parent:
        target = msh.name
    else:
        target = str(msh)
    return _escape(target.replace("\\", "/"))


def _has_three_numbers(value: Any) -> bool:
    try:
        if value is None or len(value) != 3:
            return False
        [float(x) for x in value]
        return True
    except (TypeError, ValueError):
        return False


def _bool(value: Any) -> int:
    return 1 if bool(value) else 0


def _escape(text: str) -> str:
    return text.replace('"', '\\"')
