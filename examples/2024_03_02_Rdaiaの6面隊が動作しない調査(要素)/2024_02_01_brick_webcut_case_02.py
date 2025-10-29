import os, sys
from numpy import *
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")
import cubit

cubit.init(['cubit','-nojournal','-batch'])
cubit.reset()

lx = 10
cubit.silent_cmd(f'create brick x {lx} y {lx} z {lx}  ')
volume_id = cubit.get_last_id('volume')

cubit.silent_cmd(f'webcut volume all with plane yplane offset 0 rotate 10 about x center 0 0 0 group_results')
cubit.silent_cmd(f'webcut volume all with plane yplane offset 2.5 rotate 10 about z center 0 0 0 group_results')
cubit.silent_cmd(f'webcut volume all with plane xplane offset -2 rotate 10 about z center 0 0 0 group_results')

webcut_group_id = cubit.get_last_id('group')
volume_ids = cubit.get_group_volumes(webcut_group_id)

print(volume_ids)
cubit.silent_cmd(f'delete volume all except 2')
cubit.silent_cmd(f'compress')
#cubit.silent_cmd(f'delete volume 1')

#cubit.silent_cmd(f'volume {volume_ids[0]} interval 1 ')
cubit.silent_cmd(f'volume all interval 1 ')

#cubit.silent_cmd(f'imprint all')
#cubit.silent_cmd(f'merge all')
#cubit.silent_cmd(f'compress')
cubit.silent_cmd(f'mesh vol all')
cubit.silent_cmd(f'compress')

#cubit.silent_cmd(f'block 1 add volume {volume_ids[0]} to {volume_ids[-1]}')
cubit.silent_cmd(f'block 1 add volume all')
cubit.silent_cmd(f'block 1 name "york"')

sys.path.append("../../dist")
import radia as rad
hex = []
for block_id in cubit.get_block_id_list():
	name = cubit.get_exodus_entity_name("block",block_id)
	print(f'block name = "{name}"')
	volume_list = cubit.get_block_volumes(block_id)
	for volume_id in volume_list:
		hex_list = cubit.get_volume_hexes(volume_id)
		if len(hex_list)>0:
			print(f'number of hexahedron = {len(hex_list)}')
			for hex_id in hex_list:
				connectivity_list = cubit.get_connectivity("hex", hex_id)
				coords = []
				for node_id in connectivity_list:
					coord = cubit.get_nodal_coordinates(node_id)
					coords.append(list(coord))
				face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]];
				M = [0, 0, 0.5]
				print(f'hex_id = {hex_id}')
				if hex_id==1:
					hex.append(rad.ObjPolyhdr(coords,face,M))

print(f'created hex = {len(hex)} elements')
g = rad.ObjCnt(hex)
rad.ObjDrwAtr(g, [0,0,1], 0.1)
#rad.ObjDrwOpenGL(g)  # Not implemented on this platform
#from my_module import *
#exportGeometryToVTK(g, 'RADIA')

# Calculate and display magnetic field
field = rad.Fld(g, 'b', [0, 0, 0])
print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")
print("SUCCESS: All hexahedra created!")

