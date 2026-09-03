"""
Visualize the PEEC two-layer spiral model in Cubit.

Reads point data from point.xlsx (provided by Sato-sensei) and creates
curves in Cubit. Overlays the original Parasolid CAD model for comparison.

Data:
  point.xlsx -- 3D waypoints (x, y, z) in mm
  3Dmodel_14turns.x_t -- Parasolid CAD model (meters)

Units:
  point.xlsx is in mm; Cubit/Parasolid uses meters.
  Coordinates are converted mm -> m before creating geometry.
"""

import sys
import os
import numpy as np
import pandas as pd

# --- Cubit setup ---
# Auto-detect Cubit installation
from radia.install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# --- Paths ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
POINT_FILE = os.path.join(SCRIPT_DIR, "provided_by_sato", "point.xlsx")
CAD_FILE = os.path.join(SCRIPT_DIR, "provided_by_sato", "3Dmodel_14turns.x_t")
OUT_CUB5 = os.path.join(SCRIPT_DIR, "peec_visualization.cub5")

# =============================================================
#  1. Read point data from xlsx
# =============================================================
df = pd.read_excel(POINT_FILE, header=0)
df3 = df.iloc[:, :3].apply(pd.to_numeric, errors="coerce")
df3 = df3.replace([np.inf, -np.inf], np.nan).dropna()

# Convert mm -> meters
x_mm = df3.iloc[:, 0].to_numpy()
y_mm = df3.iloc[:, 1].to_numpy()
z_mm = df3.iloc[:, 2].to_numpy()
x = x_mm * 1e-3
y = y_mm * 1e-3
z = z_mm * 1e-3

n_pts = len(x)
print(f"Read {n_pts} points from {POINT_FILE}")
print(f"  X: [{x_mm.min():.3f}, {x_mm.max():.3f}] mm")
print(f"  Y: [{y_mm.min():.3f}, {y_mm.max():.3f}] mm")
print(f"  Z: [{z_mm.min():.4f}, {z_mm.max():.4f}] mm")

# =============================================================
#  2. Import original Parasolid CAD model
# =============================================================
cubit.cmd(f'import parasolid "{CAD_FILE}"')
vol_ids = cubit.parse_cubit_list("volume", "all")
print(f"\nImported Parasolid model ({len(vol_ids)} volumes)")

for vid in vol_ids:
	cubit.cmd(f'color volume {vid} gray')
	cubit.cmd(f'volume {vid} visibility on')

# =============================================================
#  3. Create curves from point data
# =============================================================
# Detect layer by z coordinate
z_mid = (z.max() + z.min()) / 2  # midpoint between top and bottom

# Create Cubit vertices
vert_ids = []
for i in range(n_pts):
	cubit.cmd(f'create vertex {x[i]:.10g} {y[i]:.10g} {z[i]:.10g}')
	vid = cubit.get_last_id("vertex")
	vert_ids.append(vid)

print(f"Created {len(vert_ids)} vertices")

# Create curves between consecutive points, colored by layer
curve_ids_top = []
curve_ids_bottom = []
curve_ids_via = []

for i in range(n_pts - 1):
	cubit.cmd(f'create curve vertex {vert_ids[i]} {vert_ids[i+1]}')
	cid = cubit.get_last_id("curve")

	z0, z1 = z[i], z[i + 1]
	# Via: z changes sign or crosses midpoint
	if (z0 > z_mid) != (z1 > z_mid):
		curve_ids_via.append(cid)
	elif z0 > z_mid:
		curve_ids_top.append(cid)
	else:
		curve_ids_bottom.append(cid)

print(f"Created curves: {len(curve_ids_top)} top, "
      f"{len(curve_ids_via)} via, {len(curve_ids_bottom)} bottom")

# Color all curves black (top, via, bottom)
for cid in curve_ids_top + curve_ids_via + curve_ids_bottom:
	cubit.cmd(f'color curve {cid} black')

# =============================================================
#  4. Export visualization
# =============================================================
cubit.cmd(f'save cub5 "{OUT_CUB5}" overwrite')
print(f"\nSaved CUB5: {OUT_CUB5}")

# =============================================================
#  5. Summary
# =============================================================
print(f"\n{'='*60}")
print(f"  PEEC Visualization Summary (from point.xlsx)")
print(f"{'='*60}")
print(f"  Points: {n_pts}")
print(f"  Curves: {len(curve_ids_top) + len(curve_ids_via) + len(curve_ids_bottom)}")
print(f"  Top layer:    {len(curve_ids_top)} curves, z = {z_mm[z_mm > z_mid*1e3].mean():.4f} mm")
print(f"  Bottom layer: {len(curve_ids_bottom)} curves, z = {z_mm[z_mm < z_mid*1e3].mean():.4f} mm")
print(f"  Via:          {len(curve_ids_via)} curves")
print(f"\n  CAD model: gray ({len(vol_ids)} volumes from Parasolid)")

# CAD bounding box
for vid in [vol_ids[0], vol_ids[-1]]:
	bbox = cubit.get_bounding_box("volume", vid)
	print(f"  CAD Vol {vid}: x=[{bbox[0]*1e3:.1f}, {bbox[1]*1e3:.1f}] "
	      f"y=[{bbox[3]*1e3:.1f}, {bbox[4]*1e3:.1f}] mm")

print("\nDone.")
