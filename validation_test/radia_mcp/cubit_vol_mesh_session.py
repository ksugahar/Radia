"""One Cubit session for test_cubit_vol_mesh_live.py; run by the worker Python.

argv: result.json plugin_dir cubit_bin radia_mcp_source_root work_dir case
For ``case`` (cylinder, kelvin_quarter or kelvin_seam) it builds the in-memory
mesh, curves it again at several orders and, in the same session, writes
``export vol`` files at the same orders.  For each order it records the largest
distance between the mapped points of the in-memory mesh and of the reloaded
file, separately for the elements with a vertex on the z = 0 plane: in
kelvin_seam that plane is a same-material webcut seam, whose triangles the
export removes after its first curving, so a re-curving no longer projects
onto it.
"""
import json
import os
import sys
import traceback
from pathlib import Path

out_json, plugin_dir, cubit_bin, source_root, work, case = sys.argv[1:7]
sys.path.insert(0, source_root)
result = {"case": case}
try:
    import numpy as np
    import ngsolve as ng            # the fork nglib first, as the export worker does
    import netgen
    os.add_dll_directory(os.path.dirname(netgen.__file__))
    sys.path.append(cubit_bin)
    os.add_dll_directory(cubit_bin)
    os.environ["CUBIT_PLUGIN_DIR"] = plugin_dir   # export vol, for the file side
    import cubit
    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    from radia_mcp.cubit import kelvin as _kelvin
    from radia_mcp.cubit import vol_mesh
    from radia_mcp.cubit.vol_export import export_command

    journal = Path(work) / f"{case}.jou"
    if not cubit.cmd(f'playback "{journal.as_posix()}"'):
        raise RuntimeError(f"playback failed: {journal}")
    options = {"kelvin_seam": {"add_kelvin": True},
               "kelvin_quarter": {"add_kelvin": True, "kelvin_sym_x": "ht", "kelvin_sym_z": "bn"}}
    kelvin_info = _kelvin.apply_options(options[case]) if case in options else None
    kelvin_args = _kelvin.export_arguments(kelvin_info) if kelvin_info else ""

    REF = [(0.25, 0.25, 0.25), (0.5, 0.0, 0.0), (0.0, 0.5, 0.5), (0.1, 0.2, 0.3), (0.5, 0.5, 0.0)]

    def mapped(mesh):
        ir = ng.IntegrationRule(REF, [0] * len(REF))
        values = ng.CF((ng.x, ng.y, ng.z))(mesh.MapToAllElements({ng.ET.TET: ir}, ng.VOL))
        return np.asarray(values).reshape(mesh.ne, len(REF), 3)

    def summary(mesh):
        pairs = mesh.GetPeriodicNodePairs(ng.VERTEX) if kelvin_info else []
        pair_set = sorted([int(a), int(b), int(group)] for (a, b), group in pairs)
        return {"np": mesh.nv, "ne": mesh.ne, "nse": len(mesh.ngmesh.Elements2D()),
                "boundaries": list(mesh.GetBoundaries()), "materials": list(mesh.GetMaterials()),
                "periodic_pairs": len(pairs), "periodic_pair_set": pair_set}

    live = vol_mesh.mesh_from_cubit(2, kelvin=kelvin_info, plugin_dir=plugin_dir)
    mesh = ng.Mesh(live)
    result["is_live"] = vol_mesh.is_live(mesh)
    points = mesh.ngmesh.Points()
    on_z0 = np.array([any(abs(points[v].p[2]) < 1e-12 for v in el.vertices)
                      for el in mesh.ngmesh.Elements3D()])
    first_two = mapped(mesh)
    rows = []
    for order in (2, 3, 1, 2):
        if rows:                       # the mesh was built at order 2
            vol_mesh.curve(mesh, order)
        path = Path(work) / f"{case}_p{order}_{len(rows)}.vol"
        if not cubit.cmd(export_command(path, order, None, kelvin_args)):
            raise RuntimeError(f"export vol failed at order {order}")
        filed = ng.Mesh(str(path))
        row = {"order": order, "live": summary(mesh), "file": summary(filed)}
        a, b = mapped(mesh), mapped(filed)
        row["same_shape"] = a.shape == b.shape
        if row["same_shape"]:
            per_element = abs(a - b).max(axis=(1, 2))
            row["max_point_difference"] = float(per_element.max())
            row["max_off_z0"] = float(per_element[~on_z0].max(initial=0.0))
            row["max_on_z0"] = float(per_element[on_z0].max(initial=0.0))
        with ng.TaskManager():
            row["volume"] = ng.Integrate(ng.CF(1), mesh)
        rows.append(row)
    result["orders"] = rows
    back = abs(mapped(mesh) - first_two).max(axis=(1, 2))
    result["recurve_back_to_2_off_z0"] = float(back[~on_z0].max(initial=0.0))
    result["recurve_back_to_2_on_z0"] = float(back[on_z0].max(initial=0.0))
    result["coordinate_scale"] = float(abs(first_two).max())

    # Each edit of the Cubit model after a build must be refused, not curved on
    # stale geometry or labels: a fresh mesh is built, the model edited, and
    # curving that mesh again must raise.
    def refused(*commands):
        probe = ng.Mesh(vol_mesh.mesh_from_cubit(2, kelvin=kelvin_info, plugin_dir=plugin_dir))
        for command in commands:
            if not cubit.cmd(command):
                raise RuntimeError(f"Cubit refused {command!r}")
        # The refusal comes before any curving: the order stays at 2.
        try:
            vol_mesh.curve(probe, 3)
            return {"refused": False, "message": "", "order_after": probe.GetCurveOrder()}
        except RuntimeError as error:
            return {"refused": "changed" in str(error) and probe.GetCurveOrder() == 2,
                    "message": str(error), "order_after": probe.GetCurveOrder()}
    # A bad order is refused before the mesh is touched and leaves it usable.
    probe = ng.Mesh(vol_mesh.mesh_from_cubit(2, kelvin=kelvin_info, plugin_dir=plugin_dir))
    bad = {}
    try:
        vol_mesh.curve(probe, 6)
        bad["refused"] = False
    except Exception as error:
        bad["refused"] = "1-5" in str(error) and probe.GetCurveOrder() == 2
        bad["message"] = str(error)
    vol_mesh.curve(probe, 3)
    bad["curved_after"] = probe.GetCurveOrder()
    result["bad_order"] = bad
    surface = min(cubit.get_entities("surface"))
    curve = min(cubit.get_entities("curve"))
    block = min(cubit.get_block_id_list())
    # (edit, undo): a few edits make the model unexportable (a reversed arc, a
    # Kelvin sideset under another name), so each is undone before the next build.
    edits = {
        "move": ("volume all move x 0.001", None),
        "surface_normal": (f"surface {surface} normal opposite", None),
        "curve_sense": (f"curve {curve} tangent opposite", f"curve {curve} tangent opposite"),
        "block_name": (f'block {block} name "renamed_material"', None),
    }
    if cubit.get_sideset_id_list():
        sideset = min(cubit.get_sideset_id_list())
        name = cubit.get_exodus_entity_name("sideset", sideset)
        edits["sideset_name"] = (f'sideset {sideset} name "renamed_side"', f'sideset {sideset} name "{name}"')
    if cubit.get_nodeset_id_list():
        nodeset = min(cubit.get_nodeset_id_list())
        name = cubit.get_exodus_entity_name("nodeset", nodeset)
        edits["nodeset_name"] = (f'nodeset {nodeset} name "renamed_points"', f'nodeset {nodeset} name "{name}"')
    edits["new_volume"] = ("create brick x 0.01", None)
    result["edits"] = {}
    for label, (edit, undo) in edits.items():
        result["edits"][label] = refused(edit)
        if undo and not cubit.cmd(undo):
            raise RuntimeError(f"Cubit refused {undo!r}")
    # The first mesh was built before all of these edits.
    try:
        vol_mesh.curve(mesh, 3)
        result["changed_model_refused"] = False
    except RuntimeError as error:
        result["changed_model_refused"] = "changed" in str(error)
        result["changed_model_message"] = str(error)
    result["status"] = "ok"
except Exception:
    result["status"] = "error"
    result["traceback"] = traceback.format_exc()
Path(out_json).write_text(json.dumps(result, indent=1), encoding="utf-8")
