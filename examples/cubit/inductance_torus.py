"""
Torus with gap - inductance extraction test model.

Creates a torus (ring conductor) with a small section subtracted
to define source and sink faces for BEM inductance extraction.

Run in Cubit: play "inductance_torus.py"
Or standalone: python inductance_torus.py
"""

import sys
import os

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

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

print(f"Torus: R={R}, a={a}, gap={gap}")
print(f"  Volume ID: 1")

# ============================================================
# Step 2: Mesh
# ============================================================
cubit.cmd('volume 1 scheme tetmesh')
cubit.cmd(f'volume 1 size {a/3}')
cubit.cmd('mesh volume 1')

ne_tet = len(cubit.get_volume_tets(1))
print(f"  Mesh: {ne_tet} tets")

# ============================================================
# Step 3: Find source and sink surfaces (cut faces at gap)
# ============================================================
# The subtract creates two planar surfaces at y ~ +gap/2 and y ~ -gap/2.
# All other surfaces on the torus are curved (not planar).
# Find them by checking which surfaces are planar.

surfaces = cubit.get_relatives("volume", 1, "surface")
source_surfs = []
sink_surfs = []
tol = gap * 0.1

for sid in surfaces:
	s = cubit.surface(sid)
	if s.is_planar():
		# Get surface centroid y-coordinate
		cx, cy, cz = s.center_point()
		if cy > 0:
			source_surfs.append(sid)
			print(f"  Source surface: {sid} (centroid y={cy:.6f})")
		else:
			sink_surfs.append(sid)
			print(f"  Sink surface:   {sid} (centroid y={cy:.6f})")

if not source_surfs or not sink_surfs:
	print("ERROR: Could not find source/sink surfaces.")
	print("  Planar surfaces found:")
	for sid in surfaces:
		s = cubit.surface(sid)
		if s.is_planar():
			cx, cy, cz = s.center_point()
			print(f"    Surface {sid}: centroid=({cx:.4f}, {cy:.4f}, {cz:.4f})")

# ============================================================
# Step 4: Register blocks
# ============================================================
cubit.cmd('block 1 add volume 1')
cubit.cmd('block 1 name "conductor"')

cubit.cmd('block 2 add tri all')
cubit.cmd('block 2 name "boundary"')

# Source block: triangles on source surfaces
for sid in source_surfs:
	cubit.cmd(f'block 3 add tri in surface {sid}')
cubit.cmd('block 3 name "source"')

# Sink block: triangles on sink surfaces
for sid in sink_surfs:
	cubit.cmd(f'block 4 add tri in surface {sid}')
cubit.cmd('block 4 name "sink"')

n_src = sum(len(cubit.get_surface_tris(sid)) for sid in source_surfs)
n_snk = sum(len(cubit.get_surface_tris(sid)) for sid in sink_surfs)
print(f"  Source block: {n_src} tris")
print(f"  Sink block:   {n_snk} tris")

# ============================================================
# Step 5: Save
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
cub5_file = os.path.join(script_dir, 'inductance_torus.cub5')
cubit.cmd(f'save cub5 "{cub5_file}" overwrite')

print(f"\nSaved: {cub5_file}")
print(f"\nTo extract inductance:")
print(f"  Tools > Radia-NGSolve > Inductance...")
print(f"  Source: source, Sink: sink")
print(f"  Expected L ~ {4e-7 * 3.14159 * R * (7/4 + 0.5*3.14159 - 2) * 1e9:.1f} nH (Neumann approx)")
