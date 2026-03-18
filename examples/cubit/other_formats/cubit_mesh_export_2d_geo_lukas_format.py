########################################################################
###	Lukas FEM 2D file
########################################################################

def export_2D_geo_mesh(cubit, FileName):

	import numpy
	import scipy.io

	nodes = []
	node_list = []

	N = cubit.get_node_count()
	M = cubit.get_tri_count()

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
	for node_id in node_list:
		coord = cubit.get_nodal_coordinates(node_id)
		nodes.append([coord[0],coord[1]])

	for nodeset_id in cubit.get_nodeset_id_list():
		name = cubit.get_exodus_entity_name("nodeset",nodeset_id)
		curve_list = cubit.get_nodeset_curves(nodeset_id)
		node_list = []
		for curve_id in curve_list:
			node_list += cubit.get_curve_nodes(curve_id)
		nodeset = numpy.array([(name, node_list)], dtype=[('name', 'U20'), ('DBCnodes', 'O')])
		try:
			nodesets = append(nodesets,nodeset)
		except:
			nodesets = nodeset

	conn_matrix = numpy.zeros((M,3))
	center_x = numpy.zeros((M))
	center_y = numpy.zeros((M))
	block_count = cubit.get_block_count()
	regions = numpy.rec.array([("", [], [])]*(block_count), dtype=[('name', 'U20'), ('Elements', 'O'), ('Nodes', 'O')])

	for block_id in cubit.get_block_id_list():
		Elements = []
		name = cubit.get_exodus_entity_name("block",block_id)
		surface_list = cubit.get_block_surfaces(block_id)
		Nodes = []
		Elements = []
		for surface_id in surface_list:
			tri_list = cubit.get_surface_tris(surface_id)
			Elements	+= tri_list
			for tri_id in tri_list:
				x = []
				y = []
				node_list = cubit.get_connectivity('tri',tri_id)
				Nodes += node_list
				conn_matrix[tri_id-1,:] = node_list
				for node_id in node_list:
					coord = cubit.get_nodal_coordinates(node_id)
					x.append(coord[0])
					y.append(coord[1])
				center_x[tri_id-1] = numpy.mean(x)
				center_y[tri_id-1] = numpy.mean(y)
		regions[block_id-1][0] = name
		regions[block_id-1][1] = Elements
		regions[block_id-1][2] = Nodes

	geo = {'conn_matrix':conn_matrix, 'nodes':nodes, 'M':M, 'N':N, 'nodesets':nodesets , 'center_x':center_x, 'center_y':center_y, 'regions':regions }
	scipy.io.savemat(FileName, {'geo': geo}, format='5', long_field_names=True)
	return cubit

