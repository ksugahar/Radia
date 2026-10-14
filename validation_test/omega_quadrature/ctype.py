"""C-type coil source; uses a checked Cubit Kelvin mesh and sibling report."""
import json
from pathlib import Path
import numpy as np
import ngsolve as ng
import radia as rad
from radia.kelvin_identify_ngsolve import detect_kelvin_offset

def create_case(mesh_path):
    from radia.coil_builder import CoilBuilder

    level = Path(mesh_path).parent
    report = json.loads((level / "mesh_result.json").read_text(encoding="utf-8"))
    mesh_path = Path(mesh_path)
    mesh = ng.Mesh(str(mesh_path))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    physical = np.asarray(report["kelvin_physical_center_m"], dtype=float)
    radius, straight_x, straight_y = 0.0225, 0.050, 0.0625
    centre = np.array([0.0, 0.13125, 0.0])
    start = centre + np.array([0.5 * straight_x + radius, -0.5 * straight_y, 0.0])
    builder = (CoilBuilder(-2000.0).set_start(start).set_cross_section(0.035, 0.105)
               .add_straight(straight_y).add_arc(radius, 90.0)
               .add_straight(straight_x).add_arc(radius, 90.0)
               .add_straight(straight_y).add_arc(radius, 90.0)
               .add_straight(straight_x).add_arc(radius, 90.0))
    rad.UtiDelAll()
    coil = rad.ObjCnt(builder.to_radia(arc_max_segment_length=0.004))
    center = tuple((physical + offset).tolist())
    radius = float(report["kelvin_radius_m"])
    return {"mesh": mesh, "H_s": rad.RadiaField(coil, "h"),
            "H_ext": rad.KelvinRadiaFieldStrength(coil, center, radius, tuple(physical)),
            "radius": radius, "center": center, "mu_r": 1000.0,
            "controls": {"case": "ctype", "current_A": -2000.0, "mu_r": 1000.0,
                         "physical_center": physical.tolist(), "kelvin_center": center,
                         "kelvin_radius": radius, "arc_segment_length": 0.004,
                         "mesh_report": report}}
