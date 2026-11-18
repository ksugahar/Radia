"""Export a Radia Gmsh B section to Blender using headless Python processes.

Requires Gmsh 5 in the calling Python and a Blender executable.
Example: python tools/export_blender_b_section.py sphere.msh sphere.blend --blender PATH
"""
import argparse
from array import array
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


FIELD = "B_magnitude_T"
COLOR = "B_magnitude_color"


def read_section(source, axis, offset):
    import gmsh

    if gmsh.__version__.split(".")[0] != "5":
        raise ValueError("Gmsh 5 is required")
    if Path(f"{source}.opt").exists():
        raise ValueError("Remove the adjacent Gmsh .opt from the input directory")
    with source.open("rb") as file:
        if file.readline().strip() != b"$MeshFormat" or file.readline().split() != [b"4.2", b"0", b"8"]:
            raise ValueError("Expected an ASCII MSH 4.2 file")

    gmsh.initialize([], readConfigFiles=False)
    try:
        gmsh.open(str(source))
        matches = [tag for tag in gmsh.view.getTags()
                   if gmsh.view.option.getString(tag, "Name") == "B"]
        if len(matches) != 1:
            raise ValueError("Expected exactly one B view")
        view = matches[0]
        association, _, _, _, components = gmsh.view.getModelData(view, 0)
        if association != "ElementNodeData" or components != 3:
            raise ValueError("Expected three-component B ElementNodeData")
        if gmsh.view.option.getNumber(view, "NbTimeStep") != 1:
            raise ValueError("Expected one B time step")
        settings = dict(zip("ABC", (int(c == axis) for c in "xyz")))
        settings.update(D=-offset, View=gmsh.view.getIndex(view))
        for key, value in settings.items():
            gmsh.plugin.setNumber("CutPlane", key, value)
        cut = gmsh.plugin.run("CutPlane")
        types, counts, blocks = gmsh.view.getListData(cut)
    finally:
        gmsh.finalize()

    vertices, faces, magnitudes = [], [], []
    for kind, count, block in zip(types, counts, blocks, strict=True):
        size = {"VT": 3, "VQ": 4}.get(kind)
        if size is None or len(block) != count * size * 6:
            raise ValueError(f"Unexpected CutPlane data block: {kind}")
        for start in range(0, len(block), size * 6):
            entry = block[start:start + size * 6]
            # Keep every polygon's own vertices so material-side values stay separate.
            faces.append(list(range(len(vertices), len(vertices) + size)))
            vertices.extend(zip(entry[:size], entry[size:2 * size], entry[2 * size:3 * size]))
            magnitudes.extend(math.hypot(*entry[i:i + 3]) for i in range(3 * size, 6 * size, 3))
    if not faces or not all(math.isfinite(v) for v in magnitudes):
        raise ValueError("The section is empty or includes non-finite field values")
    if not all(math.isfinite(c) for point in vertices for c in point):
        raise ValueError("The section includes non-finite coordinates")
    with source.open("rb") as file:
        source_sha256 = hashlib.file_digest(file, "sha256").hexdigest()
    return dict(vertices=vertices, faces=faces, magnitudes=magnitudes,
                axis=axis, offset_m=offset, gmsh_version=gmsh.__version__,
                source_sha256=source_sha256)


def write_blend(bpy, section, output):
    from mathutils import Vector

    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    mesh = bpy.data.meshes.new("Magnetic flux density section")
    mesh.from_pydata(section["vertices"], [], section["faces"])
    mesh.update()
    obj = bpy.data.objects.new("|B| section", mesh)
    scene.collection.objects.link(obj)
    for key in ("axis", "offset_m", "source_sha256", "gmsh_version"):
        obj[key] = section[key]
    obj["field_units"] = "T"
    values = section["magnitudes"]
    lo, hi = min(values), max(values)
    obj["color_min_T"], obj["color_max_T"] = lo, hi
    attr = mesh.attributes.new(FIELD, "FLOAT", "POINT")
    attr.data.foreach_set("value", array("f", values))
    colors = mesh.color_attributes.new(COLOR, "FLOAT_COLOR", "POINT")
    mesh.color_attributes.active_color = colors
    rgba = array("f")
    for value in values:
        t = (value - lo) / (hi - lo) if hi > lo else 0.5
        rgba.extend((0.02 + 0.93 * t, 0.12 + 0.68 * t, 0.55 - 0.50 * t, 1.0))
    colors.data.foreach_set("color", rgba)

    material = bpy.data.materials.new("|B| color")
    material.use_nodes = True
    nodes = material.node_tree.nodes
    nodes.clear()
    color_node = nodes.new("ShaderNodeVertexColor")
    color_node.layer_name = COLOR
    emission = nodes.new("ShaderNodeEmission")
    surface = nodes.new("ShaderNodeOutputMaterial")
    material.node_tree.links.new(color_node.outputs["Color"], emission.inputs["Color"])
    material.node_tree.links.new(emission.outputs[0], surface.inputs["Surface"])
    mesh.materials.append(material)

    bounds = list(zip(*section["vertices"]))
    center = Vector([(min(v) + max(v)) / 2 for v in bounds])
    span = max(max(v) - min(v) for v in bounds)
    normal = Vector([int(c == section["axis"]) for c in "xyz"])
    rotation = (-normal).to_track_quat("-Z", "Y")
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "VIEW_3D":
                space = area.spaces.active
                space.shading.color_type = "VERTEX"
                space.shading.light = "FLAT"
                space.overlay.show_overlays = False
                view = space.region_3d
                view.view_rotation = rotation
                view.view_perspective = "ORTHO"
                view.view_location = center
                view.view_distance = 2 * span
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.wm.save_as_mainfile(filepath=str(output))


def verify_blend(bpy, section, output):
    bpy.ops.wm.open_mainfile(filepath=str(output), use_scripts=False)
    obj = bpy.data.objects["|B| section"]
    mesh = obj.data
    attr = mesh.attributes[FIELD]
    colors = mesh.color_attributes[COLOR]
    coords = array("f", [0]) * (3 * len(mesh.vertices))
    actual = array("f", [0]) * len(attr.data)
    mesh.vertices.foreach_get("co", coords)
    attr.data.foreach_get("value", actual)
    expected_coords = array("f", (c for vertex in section["vertices"] for c in vertex))
    expected_values = array("f", section["magnitudes"])
    if (coords != expected_coords or actual != expected_values
            or [list(p.vertices) for p in mesh.polygons] != section["faces"]):
        raise RuntimeError("Saved mesh or |B| values differ from the Gmsh section")
    if (attr.domain != "POINT" or attr.data_type != "FLOAT"
            or colors.domain != "POINT" or len(colors.data) != len(actual)
            or mesh.color_attributes.active_color.name != COLOR
            or obj["field_units"] != "T"
            or any(obj[key] != section[key] for key in
                   ("axis", "offset_m", "source_sha256", "gmsh_version"))):
        raise RuntimeError("Saved field attributes or provenance are invalid")
    print(json.dumps(dict(blender=bpy.app.version_string, vertices=len(mesh.vertices),
                          faces=len(mesh.polygons), range_T=[min(actual), max(actual)],
                          source_sha256=obj["source_sha256"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--blender", required=True, type=Path)
    parser.add_argument("--axis", choices=tuple("xyz"), default="y")
    parser.add_argument("--offset", type=float, default=0.0, help="Cut coordinate in meters")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.suffix.lower() != ".blend" or output.exists():
        parser.error("Output must be a new .blend file")
    if not math.isfinite(args.offset):
        parser.error("--offset must be finite")
    section = read_section(args.source.resolve(), args.axis, args.offset)
    with TemporaryDirectory(dir=output.parent) as tmp:
        staged = Path(tmp) / output.name
        data = Path(tmp) / "section.json"
        data.write_text(json.dumps(section), encoding="utf-8")
        command = [str(args.blender.resolve()), "--background", "--factory-startup",
                   "--disable-autoexec", "--threads", "1", "--python-exit-code", "1", "--python",
                   str(Path(__file__).resolve()), "--"]
        for operation in ("write", "verify"):
            subprocess.run([*command, operation, str(data), str(staged)], check=True)
        os.link(staged, output)


if __name__ == "__main__":
    if "--" in sys.argv:
        import bpy

        operation, data_path, blend_path = sys.argv[sys.argv.index("--") + 1:]
        section = json.loads(Path(data_path).read_text(encoding="utf-8"))
        {"write": write_blend, "verify": verify_blend}[operation](bpy, section, Path(blend_path))
    else:
        main()
