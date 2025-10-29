#!/usr/bin/env python
"""
Test case_01 with actual Cubit mesh generation
Modified to work without OpenGL and my_module dependencies
"""

import os, sys
from numpy import *

# Add Cubit to path
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import cubit

print("="*60)
print("Case 01: Brick webcut with Cubit + Radia")
print("="*60)

cubit.init(['cubit','-nojournal','-batch','-nographics','-noecho'])
cubit.reset()

# Create brick
lx = 10
cubit.silent_cmd(f'create brick x {lx} y {lx} z {lx}')
volume_id = cubit.get_last_id('volume')
print(f"Created brick volume {volume_id}")

# Webcut at 20 degree angle
angle = 20
cubit.silent_cmd(f'webcut volume {volume_id} with plane yplane offset 0.4 rotate {angle} about x center 0 -5 0 group_results')
webcut_group_id = cubit.get_last_id('group')
volume_ids = cubit.get_group_volumes(webcut_group_id)
print(f"Webcut created volumes: {volume_ids}")

# Set mesh intervals
cubit.silent_cmd(f'volume {volume_ids[0]} interval 2')
cubit.silent_cmd(f'volume {volume_ids[1]} interval 2')

# Mesh
cubit.silent_cmd(f'imprint all')
cubit.silent_cmd(f'merge all')
cubit.silent_cmd(f'compress')
cubit.silent_cmd(f'mesh vol all')
cubit.silent_cmd(f'compress')
print("Meshing completed")

# Create block
cubit.silent_cmd(f'block 1 add volume {volume_ids[0]} {volume_ids[1]}')
cubit.silent_cmd(f'block 1 name "york"')

# Import Radia
import radia as rad

# Convert Cubit hexes to Radia polyhedra
hex = []
total_hexes = 0

for block_id in cubit.get_block_id_list():
	name = cubit.get_exodus_entity_name("block", block_id)
	print(f'Block name = "{name}"')
	volume_list = cubit.get_block_volumes(block_id)

	for volume_id in volume_list:
		hex_list = cubit.get_volume_hexes(volume_id)
		if len(hex_list) > 0:
			print(f'Volume {volume_id}: {len(hex_list)} hexahedra')
			total_hexes += len(hex_list)

			for hex_id in hex_list:
				connectivity_list = cubit.get_connectivity("hex", hex_id)
				coords = []
				for node_id in connectivity_list:
					coord = cubit.get_nodal_coordinates(node_id)
					coords.append(list(coord))

				face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]
				M = [0, 0, 0.5]

				try:
					hex_obj = rad.ObjPolyhdr(coords, face, M)
					hex.append(hex_obj)
				except Exception as e:
					print(f'  ERROR creating hex {hex_id}: {e}')

print()
print(f"Total hexes from Cubit: {total_hexes}")
print(f"Successfully created Radia objects: {len(hex)}")

if len(hex) > 0:
	# Create container
	g = rad.ObjCnt(hex)
	rad.ObjDrwAtr(g, [0,0,1], 0.1)
	print(f"Container object ID: {g}")

	# Calculate field at origin
	field = rad.Fld(g, 'b', [0, 0, 0])
	print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

	# Calculate field at other points
	test_points = [[0, 0, 5], [0, 0, -5], [5, 0, 0]]
	for pt in test_points:
		field = rad.Fld(g, 'b', pt)
		print(f"Field at {pt}: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

	print()
	print("="*60)
	print("SUCCESS: All hexahedra created and field calculated!")
	print("="*60)
else:
	print()
	print("="*60)
	print("ERROR: No hexahedra were successfully created")
	print("="*60)
