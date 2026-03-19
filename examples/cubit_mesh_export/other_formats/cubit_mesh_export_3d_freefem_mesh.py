########################################################################
###	FreeFEM mesh format
########################################################################

def export_3D_mesh(cubit, FileName):

	block_count = cubit.get_block_count()
	block_list = cubit.get_block_id_list()
	node_count = cubit.get_node_count()

	with open(FileName, 'w') as fid:
		fid.write("MeshVersionFormatted 2\n")
		fid.write("\n")
		fid.write("Dimension 3\n")
		fid.write("\n")
		fid.write("Vertices\n")
		fid.write(f'{node_count}\n')
		node_list = set()
		for block_id in cubit.get_block_id_list():
			elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face"]
			for elem_type in elem_types:
				if elem_type == "hex":
					func = getattr(cubit, f"get_block_{elem_type}es")
				else:
					func = getattr(cubit, f"get_block_{elem_type}s")
				for element_id in func(block_id):
					node_ids = cubit.get_connectivity(elem_type, element_id)
					node_list.update(node_ids)
		for node_id in node_list:
			coord = cubit.get_nodal_coordinates(node_id)
			fid.write(f'{coord[0]} {coord[1]} {coord[2]} {0}\n')
		fid.write("\n")
		fid.write("Tetrahedra\n")
		tet_list = set()
		for block_id in cubit.get_block_id_list():
			tet_list.update(cubit.get_block_tets(block_id))
		fid.write(f'{len(tet_list)}\n')
		for block_id in cubit.get_block_id_list():
			tet_list = cubit.get_block_tets(block_id)
			if len(tet_list)>0:
				for tet_id in tet_list:
					node_list = cubit.get_connectivity("tet", tet_id)
					fid.write(f'{node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {block_id}\n')

		fid.write("Triangles\n")
		tri_list = set()
		for block_id in cubit.get_block_id_list():
			tri_list.update(cubit.get_block_tris(block_id))
		fid.write(f'{len(tri_list)}\n')
		for block_id in cubit.get_block_id_list():
			tet_list = cubit.get_block_tets(block_id)
			if len(tri_list)>0:
				for tri_id in tri_list:
					node_list = cubit.get_connectivity("tri", tri_id)
					fid.write(f'{node_list[0]} {node_list[1]} {node_list[2]} {block_id}\n')
		fid.write("\n")
		fid.write("End\n")
		fid.close()
	return cubit
