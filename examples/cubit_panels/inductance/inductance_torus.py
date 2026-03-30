"""
Torus with gap - inductance extraction test model.

Creates a torus (ring conductor) with a small section subtracted
to define source and sink faces for BEM inductance extraction.

Supports both tet (tri surface) and hex (quad surface) mesh types.

Usage:
    python inductance_torus.py              # tet mesh (default)
    python inductance_torus.py --mesh hex   # hex sweep mesh

Run in Cubit: play "inductance_torus.py"  (uses tet mesh)
"""

import argparse
import sys
import os

# Parse args (only when run standalone, not from Cubit playback)
mesh_type = "tet"
if __name__ == "__main__" or "__main__" not in dir():
	try:
		parser = argparse.ArgumentParser(description="Create torus model")
		parser.add_argument("--mesh", choices=["tet", "hex"], default="tet")
		args, _ = parser.parse_known_args()
		mesh_type = args.mesh
	except:
		pass  # Cubit playback: no argparse

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

import cubit

try:
	cubit.init(['cubit', '-nojournal', '-batch'])
except:
	pass  # Already initialized in GUI

cubit.cmd('reset')

# ============================================================
# Parameters
# ============================================================
R = 0.05    # Major radius [m]
a = 0.005   # Minor radius [m]
gap = 0.005 # Gap width [m]

# ============================================================
# Step 1: Create torus with gap
# ============================================================
cubit.cmd(f'create torus major radius {R} minor radius {a}')

# Cutting brick at x=R, spanning the gap in y-direction
cubit.cmd(f'brick x {3*a} y {gap} z {3*a}')
cubit.cmd(f'move Volume 2 x {R} y 0 z 0 include_merged')

# Subtract to create gap
cubit.cmd('subtract volume 2 from volume 1')

print(f"Torus: R={R}, a={a}, gap={gap}, mesh={mesh_type}")

# ============================================================
# Step 2: Mesh (tet or hex)
# ============================================================
vol_ids = list(cubit.get_entities("volume"))
vid = vol_ids[0]

if mesh_type == "hex":
	cubit.cmd(f'volume {vid} scheme sweep')
	cubit.cmd('curve all interval 8')
else:
	cubit.cmd(f'volume {vid} scheme tetmesh')
	cubit.cmd(f'volume {vid} size {a/3}')
cubit.cmd(f'mesh volume {vid}')

ne_hex = len(cubit.get_volume_hexes(vid))
ne_tet = len(cubit.get_volume_tets(vid))
print(f"  Mesh: {ne_hex} hexes, {ne_tet} tets")

# ============================================================
# Step 3: Find source and sink surfaces (cut faces at gap)
# ============================================================
surfaces = cubit.get_relatives("volume", vid, "surface")
source_surfs = []
sink_surfs = []

for sid in surfaces:
	s = cubit.surface(sid)
	if s.is_planar():
		cx, cy, cz = s.center_point()
		if cy > 0:
			source_surfs.append(sid)
			print(f"  Source surface: {sid} (centroid y={cy:.6f})")
		else:
			sink_surfs.append(sid)
			print(f"  Sink surface:   {sid} (centroid y={cy:.6f})")

if not source_surfs or not sink_surfs:
	print("ERROR: Could not find source/sink surfaces.")

# ============================================================
# Step 4: Register blocks (tri AND quad for universal support)
# ============================================================
cubit.cmd('set duplicate block elements on')

cubit.cmd(f'block 1 add volume {vid}')
cubit.cmd('block 1 name "conductor"')

# Boundary block: both tri and quad
has_tri = len(cubit.get_entities("tri")) > 0
has_quad = len(cubit.get_entities("quad")) > 0
if has_tri:
	cubit.cmd('block 2 add tri all')
if has_quad:
	cubit.cmd('block 2 add quad all')
cubit.cmd('block 2 name "boundary"')

# Source block
for sid in source_surfs:
	if has_tri:
		cubit.cmd(f'block 3 add tri in surface {sid}')
	if has_quad:
		cubit.cmd(f'block 3 add quad in surface {sid}')
cubit.cmd('block 3 name "source"')

# Sink block
for sid in sink_surfs:
	if has_tri:
		cubit.cmd(f'block 4 add tri in surface {sid}')
	if has_quad:
		cubit.cmd(f'block 4 add quad in surface {sid}')
cubit.cmd('block 4 name "sink"')

# Count surface elements
n_tri = len(cubit.get_entities("tri"))
n_quad = len(cubit.get_entities("quad"))
print(f"  Surface elements: {n_tri} tris, {n_quad} quads")

# ============================================================
# Step 5: Save
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
suffix = "_hex" if mesh_type == "hex" else ""
cub5_file = os.path.join(script_dir, f'inductance_torus{suffix}.cub5')
cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

print(f"\nSaved: {cub5_file}")
print(f"\nTo extract inductance:")
print(f"  Tools > Radia-NGSolve > Inductance...")
print(f"  Source: source, Sink: sink")
print(f"  Expected L ~ {4e-7 * 3.14159 * R * (7/4 + 0.5*3.14159 - 2) * 1e9:.1f} nH (Neumann approx)")
