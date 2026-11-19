"""Optional headless Blender worker for ``blender_field_animation``.

Runs ONLY under an explicit Python 3.13 interpreter with ``bpy`` installed
(see requirement.txt), started by the radia-mcp server as
``python -I -B field_animation_worker.py <mode> <request.json>``.  It is not
part of the radia-mcp wheel and imports nothing from radia or radia_mcp.

Modes
  write    build the scene from a solver payload, save the .blend, render
           preview frames and an MP4 (Cycles on the CPU, no GPU).
  inspect  reopen a saved .blend in a fresh process and report what it holds:
           scene summary, timing, magnet pose and the stored J per sample.

Fidelity rules: every animated quantity comes from the payload.  Magnet pose
and the J glyphs change only at solver steps (CONSTANT keyframes, no
interpolation).  The camera is static.  The J glyph objects are drawn a small
display offset above the plate top; their stored attributes keep the solved
sample point and signed J vector.
"""
from __future__ import annotations

import json
import math
import sys

import bpy  # noqa: I001  (bpy must load before bmesh)
import bmesh
import numpy as np
from mathutils import Matrix, Vector

# Perceptually ordered colour stops (viridis subset) for |J| in [0, J_max].
_STOPS = [(0.0, (0.267, 0.005, 0.329)), (0.25, (0.229, 0.322, 0.546)),
          (0.5, (0.128, 0.567, 0.551)), (0.75, (0.369, 0.789, 0.383)),
          (1.0, (0.993, 0.906, 0.144))]


def _cmap(t):
    t = min(max(t, 0.0), 1.0)
    for (a, ca), (b, cb) in zip(_STOPS, _STOPS[1:]):
        if t <= b:
            w = (t - a) / (b - a)
            return tuple(ca[i] + w * (cb[i] - ca[i]) for i in range(3))
    return _STOPS[-1][1]


def _material(name, color, alpha=1.0, emission=0.0, metallic=0.0, rough=0.5, attr=None):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    bsdf = nt.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Metallic"].default_value = metallic
    bsdf.inputs["Roughness"].default_value = rough
    bsdf.inputs["Alpha"].default_value = alpha
    if attr:
        node = nt.nodes.new("ShaderNodeAttribute")
        node.attribute_name = attr
        nt.links.new(node.outputs["Color"], bsdf.inputs["Base Color"])
        nt.links.new(node.outputs["Color"], bsdf.inputs["Emission Color"])
    if emission:
        bsdf.inputs["Emission Strength"].default_value = emission
        if not attr:
            bsdf.inputs["Emission Color"].default_value = (*color, 1.0)
    return m


def _box(name, lo, hi, mat):
    me = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bmesh.ops.create_cube(bm, size=1.0)
    for v in bm.verts:
        v.co = Vector([lo[i] + (v.co[i] + 0.5) * (hi[i] - lo[i]) for i in range(3)])
    bm.to_mesh(me)
    bm.free()
    ob = bpy.data.objects.new(name, me)
    ob.data.materials.append(mat)
    bpy.context.scene.collection.objects.link(ob)
    return ob


def _text(name, body, loc, size, mat, rot=(math.radians(90), 0, 0)):
    cu = bpy.data.curves.new(name, "FONT")
    cu.body = body
    cu.size = size
    ob = bpy.data.objects.new(name, cu)
    ob.location = loc
    ob.rotation_euler = rot
    ob.data.materials.append(mat)
    bpy.context.scene.collection.objects.link(ob)
    return ob


def _look_at(direction, up=(0.0, 1.0, 0.0)):
    """Rotation whose -Z looks along direction with +Y as close to up as possible."""
    z = -Vector(direction).normalized()
    x = Vector(up).cross(z).normalized()
    y = z.cross(x)
    return Matrix((x, y, z)).transposed().to_euler()


def _show_only_at(ob, frame, first, last):
    """Visible exactly on one frame, CONSTANT keys so nothing is interpolated."""
    for f, hide in ((first, True), (frame, False), (frame + 1, True)):
        if f < first or f > last + 1:
            continue
        ob.hide_render = hide
        ob.hide_viewport = hide
        ob.keyframe_insert("hide_render", frame=f)
        ob.keyframe_insert("hide_viewport", frame=f)
    if frame == first:
        ob.hide_render = ob.hide_viewport = False
        ob.keyframe_insert("hide_render", frame=first)
        ob.keyframe_insert("hide_viewport", frame=first)


def _constant(ob):
    ad = ob.animation_data
    if ad and ad.action:
        for layer in getattr(ad.action, "layers", []):
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for fc in bag.fcurves:
                        for k in fc.keyframe_points:
                            k.interpolation = "CONSTANT"
        for fc in getattr(ad.action, "fcurves", []):
            for k in fc.keyframe_points:
                k.interpolation = "CONSTANT"


def _arrows(name, pts, vecs, jmax, length, width, y_disp):
    """Flat arrows in the display plane y = y_disp, coloured by |J|."""
    me = bpy.data.meshes.new(name)
    verts, faces, cols = [], [], []
    for p, v in zip(pts, vecs):
        mag = float(np.linalg.norm(v))
        if mag <= 0.0:
            continue
        d = np.array([v[0], 0.0, v[2]])
        dn = np.linalg.norm(d)
        if dn <= 0.0:
            continue
        d /= dn
        n = np.array([-d[2], 0.0, d[0]])
        L = length * min(mag / jmax, 1.0)
        base = np.array([p[0], y_disp, p[2]]) - 0.5 * L * d
        tip = base + L * d
        neck = base + 0.65 * L * d
        w = width * (0.35 + 0.65 * min(mag / jmax, 1.0))
        q = [base + 0.3 * w * n, neck + 0.3 * w * n, neck + w * n, tip,
             neck - w * n, neck - 0.3 * w * n, base - 0.3 * w * n]
        i0 = len(verts)
        verts.extend(tuple(x) for x in q)
        faces.extend([(i0, i0 + 1, i0 + 5, i0 + 6), (i0 + 2, i0 + 3, i0 + 4)])
        c = _cmap(mag / jmax)
        cols.extend([c] * 7)
    me.from_pydata(verts, [], faces)
    attr = me.color_attributes.new("J_color", "FLOAT_COLOR", "POINT")
    for i, c in enumerate(cols):
        attr.data[i].color = (*c, 1.0)
    ob = bpy.data.objects.new(name, me)
    bpy.context.scene.collection.objects.link(ob)
    return ob


def _data_points(name, ids, pts, vecs):
    """Non-rendered carrier of the solved samples (id, point, signed J)."""
    me = bpy.data.meshes.new(name)
    me.from_pydata([tuple(p) for p in pts], [], [])
    a_id = me.attributes.new("sample_id", "INT", "POINT")
    a_j = me.attributes.new("J_A_per_m2", "FLOAT_VECTOR", "POINT")
    a_id.data.foreach_set("value", [int(i) for i in ids])
    a_j.data.foreach_set("vector", np.asarray(vecs, dtype=np.float64).ravel().tolist())
    ob = bpy.data.objects.new(name, me)
    ob.hide_render = True
    bpy.context.scene.collection.objects.link(ob)
    return ob


def write(req):
    pay = json.load(open(req["payload_path"], encoding="utf-8"))
    spec = req["spec"]
    bpy.ops.wm.read_factory_settings(use_empty=True)
    sc = bpy.context.scene
    n = len(pay["time_s"])
    first = 1
    sc.frame_start, sc.frame_end = first, first + n - 1
    sc.render.fps = int(spec["fps"])
    sc.unit_settings.system = "METRIC"
    sc.unit_settings.length_unit = "MILLIMETERS"
    lo, hi = pay["conductor"]["box_min_m"], pay["conductor"]["box_max_m"]
    span = max(hi[i] - lo[i] for i in range(3))
    jmax = float(pay["J_scale_max_A_per_m2"])
    y_disp = hi[1] + 0.02 * span * 0.1

    world = bpy.data.worlds.new("World")
    sc.world = world
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.012, 0.015, 0.022, 1)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.6

    copper = _material("Conductor", (0.80, 0.42, 0.25), metallic=0.9, rough=0.28)
    _box("Conductor", lo, hi, copper)
    north = _material("Magnet_N", (0.85, 0.08, 0.12), metallic=0.3, rough=0.18)
    south = _material("Magnet_S", (0.08, 0.22, 0.85), metallic=0.3, rough=0.18)
    s = pay["magnet"]["size_m"]
    mag = bpy.data.objects.new("Magnet", None)
    sc.collection.objects.link(mag)
    mdir = np.asarray(pay["magnet"]["magnetization_A_per_m"], float)
    ax = int(np.argmax(np.abs(mdir)))
    sign = 1.0 if mdir[ax] >= 0 else -1.0
    for half, mat in (("N", north), ("S", south)):
        hlo = [-0.5 * s[i] for i in range(3)]
        hhi = [0.5 * s[i] for i in range(3)]
        if (half == "N") == (sign > 0):
            hlo[ax] = 0.0
        else:
            hhi[ax] = 0.0
        part = _box(f"Magnet_{half}", hlo, hhi, mat)
        part.parent = mag
    mag.rotation_mode = "QUATERNION"
    for k in range(n):
        R = Matrix([list(r) for r in pay["magnet"]["rotation"][k]])
        mag.location = Vector(pay["magnet"]["center_m"][k])
        mag.keyframe_insert("location", frame=first + k)
        mag.rotation_quaternion = R.to_quaternion()
        mag.keyframe_insert("rotation_quaternion", frame=first + k)
    _constant(mag)

    glyph_mat = _material("J_glyph", (1, 1, 1), emission=3.0, rough=0.6, attr="J_color")
    text_mat = _material("Label", (0.92, 0.93, 0.95), emission=2.0)
    ids = pay["samples"]["id"]
    pts = pay["samples"]["point_m"]

    # Static camera; HUD texts and the legend are children of the camera, so
    # they stay fixed on screen and never imply physical motion.
    cam = bpy.data.objects.new("Camera", bpy.data.cameras.new("Camera"))
    sc.collection.objects.link(cam)
    sc.camera = cam
    c = Vector(((lo[0] + hi[0]) / 2, hi[1], (lo[2] + hi[2]) / 2))
    cam.location = c + Vector((0.0, 0.95 * span, 1.55 * span))
    cam.rotation_euler = _look_at(c - cam.location)
    cam.data.lens = 40
    cam.data.clip_start = span * 1e-3
    depth = 0.5 * span
    hw = depth * cam.data.sensor_width / (2 * cam.data.lens)
    aspect = spec["resolution"][1] / spec["resolution"][0]

    def hud(ob, x, y):
        ob.parent = cam
        ob.location = (x * hw, y * hw * aspect, -depth)
        ob.rotation_euler = (0.0, 0.0, 0.0)
        return ob

    for k in range(n):
        J = pay["J_A_per_m2"][k]
        g = _arrows(f"J_glyph_{k:04d}", pts, J, jmax, float(spec["arrow_length_m"]),
                    float(spec["arrow_width_m"]), y_disp)
        g.data.materials.append(glyph_mat)
        dpt = _data_points(f"J_data_{k:04d}", ids, pts, J)
        lab = hud(_text(f"Time_{k:04d}", f"t = {pay['time_s'][k]:.3f} s   (solver step {k} of {n - 1})",
                        (0, 0, 0), 0.045 * hw, text_mat), -0.95, -0.92)
        for ob in (g, dpt, lab):
            _show_only_at(ob, first + k, first, first + n - 1)
            _constant(ob)
        dpt.hide_render = True

    hud(_text("Title", spec["title"], (0, 0, 0), 0.055 * hw, text_mat), -0.95, 0.86)
    tt = pay["time_s"]
    phys = tt[-1] - tt[0] + (tt[1] - tt[0])
    rate = phys / (n / sc.render.fps)
    sc["playback_physical_s_per_video_s"] = rate
    hud(_text("Subtitle", f"Magnet pose and J from the same solver steps; camera static; no interpolation; "
              f"playback {rate:.3g} s physical per 1 s video",
              (0, 0, 0), 0.032 * hw, text_mat), -0.95, 0.76)
    # Fixed legend: vertical colour bar + bounds in SI, on screen.
    for i in range(32):
        t0, t1 = i / 32, (i + 1) / 32
        seg = _box(f"Legend_{i:02d}", (0.0, 0.0, 0.0), (0.025 * hw, (t1 - t0) * 0.9 * hw * aspect, 1e-6 * hw),
                   _material(f"Legend_{i:02d}", _cmap(0.5 * (t0 + t1)), emission=1.0))
        hud(seg, 0.78, -0.45 + 0.9 * t0)
    hud(_text("Legend_min", "0", (0, 0, 0), 0.032 * hw, text_mat), 0.82, -0.46)
    hud(_text("Legend_max", f"{jmax:.3g}", (0, 0, 0), 0.032 * hw, text_mat), 0.82, 0.42)
    hud(_text("Legend_unit", "|J| (A/m²)", (0, 0, 0), 0.032 * hw, text_mat), 0.74, 0.52)

    sun = bpy.data.objects.new("Key", bpy.data.lights.new("Key", "SUN"))
    sun.data.energy = 3.0
    sun.data.angle = math.radians(8)
    sun.rotation_euler = _look_at((0.4, -1.0, -0.6))
    sc.collection.objects.link(sun)
    rim = bpy.data.objects.new("Rim", bpy.data.lights.new("Rim", "SUN"))
    rim.data.energy = 1.2
    rim.rotation_euler = _look_at((-0.5, -0.6, 0.9))
    sc.collection.objects.link(rim)
    sc.view_settings.view_transform = "AgX"


    sc["radia_schema"] = pay["schema"]
    sc["radia_source"] = json.dumps(pay["provenance"])
    sc["physical_time_s"] = pay["time_s"]
    sc["frame_map"] = "frame = frame_start + solver step index; values held constant between steps"
    sc["J_scale_A_per_m2"] = [0.0, jmax]
    sc["camera_motion"] = "static"
    sc["glyph_display_offset_m"] = y_disp - pay["samples"]["point_m"][0][1]

    r = sc.render
    r.engine = "CYCLES"
    sc.cycles.device = "CPU"
    sc.cycles.samples = int(spec["samples"])
    sc.cycles.seed = 0
    r.threads_mode = "FIXED"
    r.threads = int(spec["threads"])
    r.resolution_x, r.resolution_y = spec["resolution"]
    r.film_transparent = False
    bpy.ops.wm.save_as_mainfile(filepath=req["blend_path"])
    previews = []
    for k in spec["preview_steps"]:
        sc.frame_set(first + k)
        if hasattr(r.image_settings, "media_type"):
            r.image_settings.media_type = "IMAGE"
        r.image_settings.file_format = "PNG"
        r.filepath = req["preview_pattern"].format(step=k)
        bpy.ops.render.render(write_still=True)
        previews.append(r.filepath)
    if req.get("video_path"):
        if hasattr(r.image_settings, "media_type"):  # Blender >= 5: video is a media type
            r.image_settings.media_type = "VIDEO"
        r.image_settings.file_format = "FFMPEG"
        r.ffmpeg.format = "MPEG4"
        r.ffmpeg.codec = "H264"
        r.ffmpeg.constant_rate_factor = "HIGH"
        r.filepath = req["video_path"]
        bpy.ops.render.render(animation=True)
    return {"blend_path": req["blend_path"], "previews": previews,
            "video_path": req.get("video_path"), "frames": n, "blender": bpy.app.version_string}


def _fcurves(action):
    """F-curves of an action on Blender 4.4+ slotted actions and older layouts."""
    for layer in getattr(action, "layers", []):
        for strip in layer.strips:
            for bag in strip.channelbags:
                yield from bag.fcurves
    yield from getattr(action, "fcurves", [])


def _interpolations():
    """Every keyframe interpolation used in the file (fidelity: CONSTANT only)."""
    seen = set()
    for action in bpy.data.actions:
        for fc in _fcurves(action):
            seen.update(k.interpolation for k in fc.keyframe_points)
    return seen


def _world_bbox(ob):
    pts = [ob.matrix_world @ Vector(c) for c in ob.bound_box]
    return [[min(p[i] for p in pts) for i in range(3)], [max(p[i] for p in pts) for i in range(3)]]


def inspect(req):
    bpy.ops.wm.open_mainfile(filepath=req["blend_path"])
    sc = bpy.context.scene
    out = {"blender": bpy.app.version_string, "frame_start": sc.frame_start,
           "frame_end": sc.frame_end, "fps": sc.render.fps,
           "schema": sc.get("radia_schema"), "camera_motion": sc.get("camera_motion"),
           "J_scale_A_per_m2": list(sc.get("J_scale_A_per_m2", [])),
           "physical_time_s": list(sc.get("physical_time_s", [])),
           "camera_animated": bool(sc.camera and sc.camera.animation_data),
           "objects": sorted(o.name for o in bpy.data.objects
                             if not o.name.startswith(("J_glyph_", "J_data_", "Time_", "Legend_"))),
           "n_glyph_objects": sum(o.name.startswith("J_glyph_") for o in bpy.data.objects),
           "render": {"engine": sc.render.engine, "device": sc.cycles.device,
                      "resolution": [sc.render.resolution_x, sc.render.resolution_y],
                      "samples": sc.cycles.samples, "view_transform": sc.view_settings.view_transform},
           "length_unit": sc.unit_settings.length_unit,
           "bbox_m": {o.name: _world_bbox(o) for o in bpy.data.objects
                      if o.type == "MESH" and o.name in ("Conductor", "Magnet_N", "Magnet_S")},
           "interpolations": sorted(_interpolations()),
           "steps": {}}
    for k in req["steps"]:
        sc.frame_set(sc.frame_start + k)
        mag = bpy.data.objects["Magnet"]
        mw = mag.matrix_world
        visible = sorted(o.name for o in bpy.data.objects
                         if o.name.startswith(("J_glyph_", "Time_")) and not o.hide_render)
        d = bpy.data.objects[f"J_data_{k:04d}"].data
        ids = [0] * len(d.vertices)
        d.attributes["sample_id"].data.foreach_get("value", ids)
        J = [0.0] * (3 * len(d.vertices))
        d.attributes["J_A_per_m2"].data.foreach_get("vector", J)
        out["steps"][str(k)] = {
            "center_m": list(mw.translation),
            "rotation": [list(mw.to_3x3()[i]) for i in range(3)],
            "visible": visible,
            "time_label": bpy.data.objects[f"Time_{k:04d}"].data.body,
            "sample_id": ids, "J_A_per_m2": J,
        }
    return out


def main():
    mode, path = sys.argv[1], sys.argv[2]
    req = json.load(open(path, encoding="utf-8"))
    res = {"write": write, "inspect": inspect}[mode](req)
    with open(req["result_path"], "w", encoding="utf-8") as f:
        json.dump(res, f)


if __name__ == "__main__":
    main()
