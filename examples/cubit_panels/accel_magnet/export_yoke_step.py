#!/usr/bin/env python
"""Export C-type dipole yoke quarter model as STEP.

Uses original Trelis.jou geometry (matches ELF coil model).
NO coordinate transformation — Cubit coords = beam coords.

Coordinate system:
  x: beam direction (62.5mm pole, return yoke extends to 151.25mm)
  y: transverse (105mm, symmetric about y=0)  -- NOT a symmetry plane for FEM
  z: field/vertical (25mm, gap at z=0)

Quarter model: x > 0, z > 0 (IMA '+x-z')
  x=0: sym_tangential (B parallel, natural BC for Omega)
  z=0: sym_normal (B perpendicular, Dirichlet Omega=0)

Gap center at origin (0, 0, 0).

Output: yoke_dipole.step

Usage:
    python export_yoke_step.py
"""

import sys
import os

cubit_path = os.environ.get("CUBIT_PATH",
                            r"C:\Program Files\Coreform Cubit 2025.3\bin")
if cubit_path not in sys.path:
    sys.path.insert(0, cubit_path)

import contextlib, io
with contextlib.redirect_stderr(io.StringIO()):
    import cubit
    cubit.init(['cubit', '-nojournal', '-nographics'])

cubit.cmd("reset")
cubit.cmd("set journal off")

# === Original Trelis.jou geometry (7 volumes, quarter C-dipole) ===
cubit.cmd("brick x 62.5 y 105 z 25")
cubit.cmd("brick x 80 y 50 z 25")
cubit.cmd("move Volume 2 location 71.25 -27.5 0")
cubit.cmd("brick x 40 y 100 z 25")
cubit.cmd("move Volume 3 location 131.25 -2.5 0")

cubit.cmd("modify curve 35 26 36 chamfer radius 8")
cubit.cmd("webcut volume 3 with plane yplane offset 39.5")
cubit.cmd("webcut volume 3 with plane yplane offset 27.5")
cubit.cmd("webcut volume 1 3 with plane yplane offset -2.5")

cubit.cmd("imprint volume 7 3 2 1 6")
cubit.cmd("merge volume 7 3 2 1 6")
cubit.cmd("imprint volume 5 4")
cubit.cmd("merge volume 5 4")

# NO reflects, NO rotation
# Cut at symmetry planes: keep x > 0 and z > 0 (IMA '+x-z')
cubit.cmd("webcut volume all with plane xplane offset 0")
cubit.cmd("webcut volume all with plane zplane offset 0")

# Remove x < 0 and z < 0 parts
vol_ids = cubit.parse_cubit_list("volume", "all")
remove_ids = []
for vid in vol_ids:
    bbox = cubit.get_bounding_box("volume", vid)
    x_center = 0.5 * (bbox[0] + bbox[1])
    z_center = 0.5 * (bbox[6] + bbox[7])
    if x_center < -1e-6 or z_center < -1e-6:
        remove_ids.append(vid)

if remove_ids:
    ids_str = " ".join(str(v) for v in remove_ids)
    cubit.cmd(f"delete volume {ids_str}")
    print(f"Removed {len(remove_ids)} volumes (x<0 or z<0)")

# Unite, scale to meters
cubit.cmd("delete mesh volume all propagate")
cubit.cmd("unite volume all")
cubit.cmd("volume all scale 0.001")

# Export STEP
out_dir = os.path.dirname(os.path.abspath(__file__))
step_file = os.path.join(out_dir, "yoke_dipole.step")
cubit.cmd(f'export step "{step_file}" overwrite')

# Print info
vol_ids = cubit.parse_cubit_list("volume", "all")
vid = vol_ids[0] if vol_ids else 1
bbox = cubit.get_bounding_box("volume", vid)
print(f"Exported: {step_file}")
print(f"Quarter model (x>0, z>0), gap center at origin")
print(f"BBox (m): x=[{bbox[0]:.6f}, {bbox[1]:.6f}] "
      f"y=[{bbox[3]:.6f}, {bbox[4]:.6f}] "
      f"z=[{bbox[6]:.6f}, {bbox[7]:.6f}]")
print(f"Symmetry: x=0 (sym_tangential), z=0 (sym_normal)")
