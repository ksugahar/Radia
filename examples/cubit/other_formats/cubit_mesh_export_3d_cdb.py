########################################################################
###	3D CDB file (ANSYS format)
###	Status: Experimental - may have issues
########################################################################

from datetime import datetime

def export_3D_CDB(cubit, FileName):
	"""Export mesh to ANSYS CDB format.

	Note: This is experimental code in other_formats/.
	"""
	# First pass: count nodes and elements
	node_set = set()
	elem_count = 0
	for block_id in cubit.get_block_id_list():
		elem_types = ["hex", "tet", "wedge", "pyramid"]
		for elem_type in elem_types:
			if elem_type == "hex":
				func = getattr(cubit, f"get_block_{elem_type}es")
			else:
				func = getattr(cubit, f"get_block_{elem_type}s")
			for element_id in func(block_id):
				elem_count += 1
				node_ids = cubit.get_connectivity(elem_type, element_id)
				node_set.update(node_ids)
	node_count = len(node_set)

	# Write file
	fid = open(FileName,'w',encoding='UTF-8')
	fid.write(f'/COM,ANSYS RELEASE 15.0\n')
	fid.write(f'! Exported from Coreform Cubit 2024.3\n')
	fid.write(f'! {datetime.now().strftime("%Y/%m/%d %I:%M:%S %p")}\n')
	fid.write(f'/PREP7\n')
	fid.write(f'/TITLE,\n')
	fid.write(f'*IF,_CDRDOFF,EQ,1,THEN\n')
	fid.write(f'_CDRDOFF= \n')
	fid.write(f'*ELSE\n')
	fid.write(f'NUMOFF,NODE, {node_count:8d}\n')
	fid.write(f'NUMOFF,ELEM, {elem_count:8d}\n')
	fid.write(f'NUMOFF,TYPE,        1\n')
	fid.write(f'*ENDIF\n')
	fid.write(f'DOF,DELETE\n')
	fid.write(f'ET,       1,185\n')
	fid.write(f'NBLOCK,6,SOLID,  {node_count:8d},  {node_count:8d}\n')
	fid.write(f'(3i9,6e20.13)\n')

	# Write nodes
	node_list = set()
	for block_id in cubit.get_block_id_list():
		elem_types = ["hex", "tet", "wedge", "pyramid"]
		for elem_type in elem_types:
			if elem_type == "hex":
				func = getattr(cubit, f"get_block_{elem_type}es")
			else:
				func = getattr(cubit, f"get_block_{elem_type}s")
			for element_id in func(block_id):
				node_ids = cubit.get_connectivity(elem_type, element_id)
				node_list.update(node_ids)
	fid.write(f'{len(node_list)}\n')
	for node_id in node_list:
		coord = cubit.get_nodal_coordinates(node_id)
		fid.write(f'{node_id:9d}{0:9d}{0:9d}{coord[0]:20.13e}{coord[1]:20.13e}{coord[2]:20.13e}\n')

	fid.write(f'N,R5.3,LOC, -1\n')

	# Write elements
	elem_count = 0
	for block_id in cubit.get_block_id_list():
		fid.write(f'EBLOCK, 19, SOLID\n')
		fid.write(f'(19i9)\n')
		for volume_id in cubit.get_block_volumes(block_id):
			hex_list = cubit.get_volume_hexes(volume_id)
			if len(hex_list)>0:
				for hex_id in hex_list:
					elem_count += 1
					node_list = cubit.get_connectivity("hex", hex_id)
					fid.write(f'{block_id:9d}{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{8:9d}{0:9d}{elem_count:9d}{node_list[0]:9d}{node_list[1]:9d}{node_list[2]:9d}{node_list[3]:9d}{node_list[4]:9d}{node_list[5]:9d}{node_list[6]:9d}{node_list[7]:9d}\n')
			wedge_list = cubit.get_volume_wedges(volume_id)
			if len(wedge_list)>0:
				for wedge_id in wedge_list:
					elem_count += 1
					node_list = cubit.get_connectivity("wedge", wedge_id)
					fid.write(f'{block_id:9d}{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{8:9d}{0:9d}{elem_count:9d}{node_list[0]:9d}{node_list[1]:9d}{node_list[2]:9d}{node_list[2]:9d}{node_list[3]:9d}{node_list[4]:9d}{node_list[5]:9d}{node_list[5]:9d}\n')
			tet_list = cubit.get_volume_tets(volume_id)
			if len(tet_list)>0:
				for tet_id in tet_list:
					elem_count += 1
					node_list = cubit.get_connectivity("tet", tet_id)
					fid.write(f'{block_id:9d}{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{8:9d}{0:9d}{elem_count:9d}{node_list[0]:9d}{node_list[1]:9d}{node_list[2]:9d}{node_list[2]:9d}{node_list[3]:9d}{node_list[3]:9d}{node_list[3]:9d}{node_list[3]:9d}\n')
			pyramid_list = cubit.get_volume_pyramids(volume_id)
			if len(pyramid_list)>0:
				for pyramid_id in pyramid_list:
					elem_count += 1
					node_list = cubit.get_connectivity("pyramid", pyramid_id)
					fid.write(f'{block_id:9d}{1:9d}{1:9d}{1:9d}{0:9d}{0:9d}{0:9d}{0:9d}{8:9d}{0:9d}{elem_count:9d}{node_list[0]:9d}{node_list[1]:9d}{node_list[2]:9d}{node_list[3]:9d}{node_list[4]:9d}{node_list[4]:9d}{node_list[4]:9d}{node_list[4]:9d}\n')
		fid.write(f'       -1\n')

	for block_id in cubit.get_block_id_list():
		name = cubit.get_exodus_entity_name("block",block_id)
		elem_list = []
		volume_list = cubit.get_block_volumes(block_id)
		for volume_id in volume_list:
			hex_list = cubit.get_volume_hexes(volume_id)
			elem_list.extend(hex_list)
			wedge_list = cubit.get_volume_wedges(volume_id)
			elem_list.extend(wedge_list)
			tet_list = cubit.get_volume_tets(volume_id)
			elem_list.extend(tet_list)
			pyramid_list = cubit.get_volume_pyramids(volume_id)
			elem_list.extend(pyramid_list)
		fid.write(f'CMBLOCK,{name:<8},ELEM,{len(elem_list):8d}\n')
		fid.write(f'(8i10)\n')
		for n in range(0, len(elem_list), 8):
			strLine = ""
			for m in range(n, min(n+8, len(elem_list))):
				strLine += f'{elem_list[m]:10d}'
			fid.write(f'{strLine}\n')

	for nodeset_id in cubit.get_nodeset_id_list():
		name = cubit.get_exodus_entity_name("nodeset",nodeset_id)

		node_list = set()
		for block_id in cubit.get_block_id_list():
			elem_types = ["tri", "quad"]
			for elem_type in elem_types:
				func = getattr(cubit, f"get_block_{elem_type}s")
				for element_id in func(block_id):
					node_ids = cubit.get_connectivity(elem_type, element_id)
					node_list.update(node_ids)

		fid.write(f'CMBLOCK,{name:<8},NODE,{len(node_list):8d}\n')
		fid.write(f'(8i10)\n')
		node_list = list(node_list)
		for n in range(0, len(node_list), 8):
			strLine = ""
			for m in range(n, min(n+8, len(node_list))):
				strLine += f'{node_list[m]:10d}'
			fid.write(f'{strLine}\n')

	fid.close()
	return cubit
