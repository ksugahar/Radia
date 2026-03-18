import math
from typing import Any, List, Optional, Tuple

########################################################################
###	Common utility functions
########################################################################

def _block_contains_geometry(cubit: Any, block_id: int) -> bool:
	"""Check if a block contains geometry (volume, surface, curve, vertex) instead of mesh elements."""
	try:
		if len(cubit.get_block_volumes(block_id)) > 0:
			return True
	except:
		pass
	try:
		if len(cubit.get_block_surfaces(block_id)) > 0:
			return True
	except:
		pass
	try:
		if len(cubit.get_block_curves(block_id)) > 0:
			return True
	except:
		pass
	try:
		if len(cubit.get_block_vertices(block_id)) > 0:
			return True
	except:
		pass
	return False


def _get_block_elements(cubit: Any, block_id: int, elem_type: str) -> Tuple[int, ...]:
	"""Get mesh elements from a block, supporting both geometry and mesh element blocks.

	This function mimics Cubit's 'draw <elem_type> in block X' behavior.
	When a block contains geometry (volume, surface, curve, vertex),
	it uses parse_cubit_list to get the associated mesh elements.
	When a block contains direct mesh elements, it uses the standard get_block_* functions.

	Args:
		cubit: Cubit Python interface object
		block_id: Block ID
		elem_type: Element type ('tet', 'hex', 'wedge', 'pyramid', 'tri', 'face',
		           'quad', 'edge', 'node')

	Returns:
		Tuple of element IDs
	"""
	# Check if block contains geometry (volume, surface, curve, vertex)
	if _block_contains_geometry(cubit, block_id):
		# Use parse_cubit_list for geometry-based blocks
		# This mimics Cubit's 'draw <elem_type> in block X' behavior
		try:
			elements = cubit.parse_cubit_list(elem_type, f"in block {block_id}")
			return tuple(elements)
		except:
			return ()

	# For mesh element blocks, use the standard get_block_* functions
	try:
		if elem_type == "hex":
			return cubit.get_block_hexes(block_id)
		elif elem_type == "tet":
			return cubit.get_block_tets(block_id)
		elif elem_type == "wedge":
			return cubit.get_block_wedges(block_id)
		elif elem_type == "pyramid":
			return cubit.get_block_pyramids(block_id)
		elif elem_type == "tri":
			return cubit.get_block_tris(block_id)
		elif elem_type == "face":
			return cubit.get_block_faces(block_id)
		elif elem_type == "quad":
			return cubit.get_block_quads(block_id)
		elif elem_type == "edge":
			return cubit.get_block_edges(block_id)
		elif elem_type == "node":
			return cubit.get_block_nodes(block_id)
		else:
			return ()
	except:
		return ()




def _warn_mixed_element_types_in_blocks(cubit: Any) -> None:
	"""Warn if any block contains multiple 3D or 2D element types.

	When a block contains multiple element types (e.g., tet + hex + pyramid),
	Cubit's 'block X element type' command only reports the last set type via
	get_block_element_type(). While Cubit does convert all elements to 2nd order
	when you issue commands like 'block 1 element type tetra10', this warning
	helps users understand they need to issue multiple commands for mixed blocks.

	This function prints a warning for each block with mixed element types.
	"""
	# 3D element types to check
	volume_elem_types = ["hex", "tet", "wedge", "pyramid"]
	# 2D element types to check
	surface_elem_types = ["tri", "face"]

	for block_id in cubit.get_block_id_list():
		# Check 3D elements
		found_3d = []
		for elem_type in volume_elem_types:
			elements = _get_block_elements(cubit, block_id, elem_type)
			if len(elements) > 0:
				found_3d.append(elem_type)

		if len(found_3d) > 1:
			block_name = cubit.get_exodus_entity_name("block", block_id)
			type_str = ", ".join(found_3d)
			print(f"WARNING: Block {block_id} ('{block_name}') contains multiple 3D element types: {type_str}")
			print(f"  To convert all elements to 2nd order, you need to issue separate commands:")
			for elem_type in found_3d:
				if elem_type == "tet":
					print(f"    block {block_id} element type tetra10")
				elif elem_type == "hex":
					print(f"    block {block_id} element type hex20")
				elif elem_type == "wedge":
					print(f"    block {block_id} element type wedge15")
				elif elem_type == "pyramid":
					print(f"    block {block_id} element type pyramid13")

		# Check 2D elements
		found_2d = []
		for elem_type in surface_elem_types:
			elements = _get_block_elements(cubit, block_id, elem_type)
			if len(elements) > 0:
				found_2d.append(elem_type)

		if len(found_2d) > 1:
			block_name = cubit.get_exodus_entity_name("block", block_id)
			type_str = ", ".join(["tri" if t == "tri" else "quad" for t in found_2d])
			print(f"WARNING: Block {block_id} ('{block_name}') contains multiple 2D element types: {type_str}")
			print(f"  To convert all elements to 2nd order, you need to issue separate commands:")
			for elem_type in found_2d:
				if elem_type == "tri":
					print(f"    block {block_id} element type tri6")
				elif elem_type == "face":
					print(f"    block {block_id} element type quad8")

########################################################################
###	Gmsh format (unified interface)
########################################################################

def export_Gmesh(cubit: Any, FileName: str, version: str = "2.2", DIM: str = "auto") -> Any:
	"""Export mesh to Gmsh format.

	Unified interface for Gmsh export supporting v2.2 and v4.1 formats.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .msh file
		version: Gmsh format version ("2.2" or "4.1")
		DIM: Dimension mode (v4.1 only)
			- "auto": Auto-detect (3D if volume elements exist, else 2D)
			- "2D": 2D mode (orient surface normals to +z, z-coords set to 0)
			- "3D": 3D mode (no normal orientation)

	Returns:
		cubit: The cubit object (for method chaining)
	"""
	if version in ("2", "2.2"):
		return _export_gmsh_v2(cubit, FileName)
	elif version in ("4", "4.1"):
		return _export_gmsh_v4(cubit, FileName, DIM=DIM)
	else:
		raise ValueError(
			f"Unsupported Gmsh version '{version}'. Use '2.2' or '4.1'.")

########################################################################
###	Gmsh v2.2 implementation (internal)
########################################################################

def _export_gmsh_v2(cubit: Any, FileName: str) -> Any:
	"""Export mesh to Gmsh format version 2.2.

	Exports 1D, 2D, and 3D mesh elements from Cubit to Gmsh v2.2 format.
	Supports both first-order and second-order elements.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .msh file

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 1st order: Point, Line, Triangle, Quad, Tetrahedron, Hexahedron, Wedge, Pyramid
		- 2nd order: Line3, Triangle6, Quad8/9, Tetrahedron10/11, Hexahedron20, Wedge15, Pyramid13
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	with open(FileName, 'w') as fid:

		fid.write("$MeshFormat\n")
		fid.write("2.2 0 8\n")
		fid.write("$EndMeshFormat\n")

		fid.write("$PhysicalNames\n")
		fid.write(f'{cubit.get_block_count()}\n')
		for block_id in cubit.get_block_id_list():
			name = cubit.get_exodus_entity_name("block", block_id)
			if len(_get_block_elements(cubit, block_id, "node")) > 0:
				fid.write(f'0 {block_id} "{name}"\n')
			elif len(_get_block_elements(cubit, block_id, "edge")) > 0:
				fid.write(f'1 {block_id} "{name}"\n')
			elif len(_get_block_elements(cubit, block_id, "tri")) + len(_get_block_elements(cubit, block_id, "face")) > 0:
				fid.write(f'2 {block_id} "{name}"\n')
			else:
				fid.write(f'3 {block_id} "{name}"\n')
		fid.write('$EndPhysicalNames\n')

		fid.write("$Nodes\n")
		node_list = set()
		for block_id in cubit.get_block_id_list():
			elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
			for elem_type in elem_types:
				for element_id in _get_block_elements(cubit, block_id, elem_type):
					node_ids = cubit.get_expanded_connectivity(elem_type, element_id)
					node_list.update(node_ids)

		fid.write(f'{len(node_list)}\n')
		for node_id in node_list:
			coord = cubit.get_nodal_coordinates(node_id)
			fid.write(f'{node_id} {coord[0]} {coord[1]} {coord[2]}\n')
		fid.write('$EndNodes\n')

		hex_list = set()
		tet_list = set()
		wedge_list = set()
		pyramid_list = set()
		tri_list = set()
		quad_list = set()
		edge_list = set()
		node_list = set()

		for block_id in cubit.get_block_id_list():
			tet_list.update(_get_block_elements(cubit, block_id, "tet"))
			hex_list.update(_get_block_elements(cubit, block_id, "hex"))
			wedge_list.update(_get_block_elements(cubit, block_id, "wedge"))
			pyramid_list.update(_get_block_elements(cubit, block_id, "pyramid"))
			tri_list.update(_get_block_elements(cubit, block_id, "tri"))
			quad_list.update(_get_block_elements(cubit, block_id, "face"))
			edge_list.update(_get_block_elements(cubit, block_id, "edge"))
			node_list.update(_get_block_elements(cubit, block_id, "node"))

		element_id = 0
		fid.write('$Elements\n')
		fid.write(f'{len(hex_list) + len(tet_list) + len(wedge_list) + len(pyramid_list) + len(tri_list) + len(quad_list) + len(edge_list) + len(node_list)}\n')

		for block_id in cubit.get_block_id_list():

			tet_list = _get_block_elements(cubit, block_id, "tet")
			for tet_id in tet_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("tet", tet_id)
				if len(node_list)==4:
					fid.write(f'{element_id} { 4} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]}\n')
				elif len(node_list)==10:
					fid.write(f'{element_id} {11} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]} {node_list[9]} {node_list[8]}\n')
				elif len(node_list)==11:
					fid.write(f'{element_id} {35} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]} {node_list[9]} {node_list[8]} {node_list[10]}\n')

			hex_list = _get_block_elements(cubit, block_id, "hex")
			for hex_id in hex_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("hex", hex_id)
				if len(node_list)==8:
					fid.write(f'{element_id} { 5} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]}\n')
				elif len(node_list)==20:
					fid.write(f'{element_id} {17} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]} {node_list[8]} {node_list[11]} {node_list[12]} {node_list[9]} {node_list[13]} {node_list[10]} {node_list[14]} {node_list[15]} {node_list[16]} {node_list[19]} {node_list[17]} {node_list[18]}\n')

			wedge_list = _get_block_elements(cubit, block_id, "wedge")
			for wedge_id in wedge_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("wedge", wedge_id)
				if len(node_list)==6:
					fid.write(f'{element_id} { 6} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]}\n')
				elif len(node_list)==15:
					fid.write(f'{element_id} {18} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[8]} {node_list[9]} {node_list[7]} {node_list[10]} {node_list[11]} {node_list[12]} {node_list[14]} {node_list[13]}\n')

			pyramid_list = _get_block_elements(cubit, block_id, "pyramid")
			for pyramid_id in pyramid_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("pyramid", pyramid_id)
				if len(node_list)==5:
					fid.write(f'{element_id} { 7} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]}\n')
				elif len(node_list)==13:
					fid.write(f'{element_id} {19} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[8]} {node_list[9]} {node_list[6]} {node_list[10]} {node_list[7]} {node_list[11]} {node_list[12]} \n')

			tri_list = _get_block_elements(cubit, block_id, "tri")
			for tri_id in tri_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("tri", tri_id)
				if len(node_list)==3:
					fid.write(f'{element_id} { 2} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]}\n')
				elif len(node_list)==6:
					fid.write(f'{element_id} { 9} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]}\n')
				elif len(node_list)==7:
					fid.write(f'{element_id} {42} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]}\n')

			quad_list = _get_block_elements(cubit, block_id, "face")
			for quad_id in quad_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("quad", quad_id)
				if len(node_list)==4:
					fid.write(f'{element_id} { 3} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]}\n')
				elif len(node_list)==8:
					fid.write(f'{element_id} {16} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]}\n')
				elif len(node_list)==9:
					fid.write(f'{element_id} {10} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]} {node_list[8]}\n')

			edge_list = _get_block_elements(cubit, block_id, "edge")
			for edge_id in edge_list:
				element_id += 1
				node_list = cubit.get_expanded_connectivity("edge", edge_id)
				if len(node_list)==2:
					fid.write(f'{element_id} {1} {2} {block_id} {block_id} {node_list[0]} {node_list[1]}\n')
				elif len(node_list)==3:
					fid.write(f'{element_id} {8} {2} {block_id} {block_id} {node_list[0]} {node_list[1]} {node_list[2]}\n')

			node_list = _get_block_elements(cubit, block_id, "node")
			for node_id in node_list:
				element_id += 1
				fid.write(f'{element_id} {15} {2} {block_id} {block_id} {node_id}\n')

		fid.write('$EndElements\n')
	return cubit

########################################################################
###	Gmsh v4.1 implementation (internal)
########################################################################

def _export_gmsh_v4(cubit: Any, FileName: str, DIM: str = "auto") -> Any:
	"""Export mesh to Gmsh format version 4.1.

	Exports 1D, 2D, and 3D mesh elements from Cubit to Gmsh v4.1 format.
	Supports both first-order and second-order elements.
	Includes full $Entities section with geometry topology.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .msh file
		DIM: Dimension mode
			- "auto": Auto-detect (3D if volume elements exist, else 2D)
			- "2D": 2D mode (orient surface normals to +z direction, z-coordinates set to 0)
			- "3D": 3D mode (no normal orientation)

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 1st order: Point, Line, Triangle, Quad, Tetrahedron, Hexahedron, Wedge, Pyramid
		- 2nd order: Line3, Triangle6, Quad8/9, Tetrahedron10/11, Hexahedron20, Wedge15, Pyramid13
	"""
	_warn_mixed_element_types_in_blocks(cubit)


	# ============================================================
	# Gmsh element type codes (same as v2)
	# ============================================================
	GMSH_POINT = 15
	GMSH_LINE2 = 1
	GMSH_LINE3 = 8
	GMSH_TRI3 = 2
	GMSH_TRI6 = 9
	GMSH_TRI7 = 42
	GMSH_QUAD4 = 3
	GMSH_QUAD8 = 16
	GMSH_QUAD9 = 10
	GMSH_TET4 = 4
	GMSH_TET10 = 11
	GMSH_TET11 = 35
	GMSH_HEX8 = 5
	GMSH_HEX20 = 17
	GMSH_WEDGE6 = 6
	GMSH_WEDGE15 = 18
	GMSH_PYRAMID5 = 7
	GMSH_PYRAMID13 = 19

	# ============================================================
	# Collect block information and determine dimension
	# ============================================================

	# Collect elements by type from all blocks
	all_tets = set()
	all_hexes = set()
	all_wedges = set()
	all_pyramids = set()
	all_tris = set()
	all_quads = set()
	all_edges = set()
	all_points = set()

	# Map element to block_id
	elem_to_block = {}

	for block_id in cubit.get_block_id_list():
		for tet_id in _get_block_elements(cubit, block_id, "tet"):
			all_tets.add(tet_id)
			elem_to_block[('tet', tet_id)] = block_id
		for hex_id in _get_block_elements(cubit, block_id, "hex"):
			all_hexes.add(hex_id)
			elem_to_block[('hex', hex_id)] = block_id
		for wedge_id in _get_block_elements(cubit, block_id, "wedge"):
			all_wedges.add(wedge_id)
			elem_to_block[('wedge', wedge_id)] = block_id
		for pyramid_id in _get_block_elements(cubit, block_id, "pyramid"):
			all_pyramids.add(pyramid_id)
			elem_to_block[('pyramid', pyramid_id)] = block_id
		for tri_id in _get_block_elements(cubit, block_id, "tri"):
			all_tris.add(tri_id)
			elem_to_block[('tri', tri_id)] = block_id
		for quad_id in _get_block_elements(cubit, block_id, "face"):
			all_quads.add(quad_id)
			elem_to_block[('quad', quad_id)] = block_id
		for edge_id in _get_block_elements(cubit, block_id, "edge"):
			all_edges.add(edge_id)
			elem_to_block[('edge', edge_id)] = block_id
		for node_id in _get_block_elements(cubit, block_id, "node"):
			all_points.add(node_id)
			elem_to_block[('node', node_id)] = block_id

	has_3d_elements = len(all_tets) + len(all_hexes) + len(all_wedges) + len(all_pyramids) > 0
	has_2d_elements = len(all_tris) + len(all_quads) > 0

	# Determine dimension mode
	if DIM == "auto":
		actual_dim = "3D" if has_3d_elements else "2D"
	else:
		actual_dim = DIM

	# ============================================================
	# Collect nodes
	# ============================================================

	node_set = set()
	for block_id in cubit.get_block_id_list():
		elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
		for elem_type in elem_types:
			for element_id in _get_block_elements(cubit, block_id, elem_type):
				node_ids = cubit.get_expanded_connectivity(elem_type, element_id)
				node_set.update(node_ids)

	sorted_nodes = sorted(node_set)
	num_nodes = len(sorted_nodes)
	min_node_tag = min(sorted_nodes) if sorted_nodes else 0
	max_node_tag = max(sorted_nodes) if sorted_nodes else 0

	# ============================================================
	# Collect geometry entities for $Entities section
	# ============================================================

	# Get volumes and surfaces associated with blocks
	volume_set = set()
	surface_set = set()
	curve_set = set()
	vertex_set = set()

	# Physical tags: entity_id -> block_id
	volume_physical = {}
	surface_physical = {}
	curve_physical = {}
	vertex_physical = {}

	for block_id in cubit.get_block_id_list():
		# 3D blocks -> volumes
		for vol_id in cubit.get_block_volumes(block_id):
			volume_set.add(vol_id)
			volume_physical[vol_id] = block_id
			# Get surfaces of this volume
			try:
				vol = cubit.volume(vol_id)
				for surf in vol.surfaces():
					surface_set.add(surf.id())
					# Get curves of this surface
					for curve in surf.curves():
						curve_set.add(curve.id())
						# Get vertices of this curve
						for vert in curve.vertices():
							vertex_set.add(vert.id())
			except:
				# Fallback: use parse_cubit_list
				surf_ids = cubit.parse_cubit_list("surface", f"in volume {vol_id}")
				surface_set.update(surf_ids)

		# 2D blocks -> surfaces
		tri_list = _get_block_elements(cubit, block_id, "tri")
		quad_list = _get_block_elements(cubit, block_id, "face")
		if len(tri_list) + len(quad_list) > 0:
			# Find surfaces containing these elements
			for tri_id in tri_list:
				try:
					owner = cubit.get_geometry_owner("tri", tri_id)
					if owner and "surface" in owner.lower():
						surf_id = int(owner.split()[1])
						surface_set.add(surf_id)
						surface_physical[surf_id] = block_id
				except:
					pass
			for quad_id in quad_list:
				try:
					owner = cubit.get_geometry_owner("quad", quad_id)
					if owner and "surface" in owner.lower():
						surf_id = int(owner.split()[1])
						surface_set.add(surf_id)
						surface_physical[surf_id] = block_id
				except:
					pass

	# If no geometry found, create minimal entities based on blocks
	if not volume_set and not surface_set:
		# Use block IDs as entity tags
		for block_id in cubit.get_block_id_list():
			if len(_get_block_elements(cubit, block_id, "tet")) + len(_get_block_elements(cubit, block_id, "hex")) + \
			   len(_get_block_elements(cubit, block_id, "wedge")) + len(_get_block_elements(cubit, block_id, "pyramid")) > 0:
				volume_set.add(block_id)
				volume_physical[block_id] = block_id
			elif len(_get_block_elements(cubit, block_id, "tri")) + len(_get_block_elements(cubit, block_id, "face")) > 0:
				surface_set.add(block_id)
				surface_physical[block_id] = block_id

	# Always create entities for 1D/0D blocks (geometry collection doesn't cover these)
	for block_id in cubit.get_block_id_list():
		if len(_get_block_elements(cubit, block_id, "edge")) > 0 and block_id not in curve_set:
			curve_set.add(block_id)
			curve_physical[block_id] = block_id
		if len(_get_block_elements(cubit, block_id, "node")) > 0 and block_id not in vertex_set:
			vertex_set.add(block_id)
			vertex_physical[block_id] = block_id

	# ============================================================
	# Write Gmsh v4 file
	# ============================================================

	with open(FileName, 'w') as fid:

		# ------------------------------------------------------------
		# $MeshFormat
		# ------------------------------------------------------------
		fid.write('$MeshFormat\n')
		fid.write('4.1 0 8\n')
		fid.write('$EndMeshFormat\n')

		# ------------------------------------------------------------
		# $PhysicalNames
		# ------------------------------------------------------------
		fid.write('$PhysicalNames\n')
		fid.write(f'{cubit.get_block_count()}\n')
		for block_id in cubit.get_block_id_list():
			name = cubit.get_exodus_entity_name("block", block_id)
			if len(_get_block_elements(cubit, block_id, "node")) > 0:
				dim = 0
			elif len(_get_block_elements(cubit, block_id, "edge")) > 0:
				dim = 1
			elif len(_get_block_elements(cubit, block_id, "tri")) + len(_get_block_elements(cubit, block_id, "face")) > 0:
				dim = 2
			else:
				dim = 3
			fid.write(f'{dim} {block_id} "{name}"\n')
		fid.write('$EndPhysicalNames\n')

		# ------------------------------------------------------------
		# $Entities
		# ------------------------------------------------------------
		num_points_ent = len(vertex_set)
		num_curves_ent = len(curve_set)
		num_surfaces_ent = len(surface_set) if surface_set else len([b for b in cubit.get_block_id_list() if len(_get_block_elements(cubit, b, "tri")) + len(_get_block_elements(cubit, b, "face")) > 0])
		num_volumes_ent = len(volume_set) if volume_set else len([b for b in cubit.get_block_id_list() if len(_get_block_elements(cubit, b, "tet")) + len(_get_block_elements(cubit, b, "hex")) + len(_get_block_elements(cubit, b, "wedge")) + len(_get_block_elements(cubit, b, "pyramid")) > 0])

		fid.write('$Entities\n')
		fid.write(f'{num_points_ent} {num_curves_ent} {num_surfaces_ent} {num_volumes_ent}\n')

		# Points (vertices)
		for vert_id in sorted(vertex_set):
			phys_tag = vertex_physical.get(vert_id, 0)
			num_phys = 1 if phys_tag else 0
			phys_str = f' {phys_tag}' if phys_tag else ''
			if phys_tag and vert_id == phys_tag:
				# Synthetic entity from node block: compute coords from block nodes
				node_ids = list(_get_block_elements(cubit, vert_id, "node"))
				if node_ids:
					coord = cubit.get_nodal_coordinates(node_ids[0])
				else:
					coord = (0, 0, 0)
			else:
				# Geometric vertex
				try:
					coord = cubit.vertex(vert_id).coordinates()
				except:
					coord = (0, 0, 0)
			fid.write(f'{vert_id} {coord[0]} {coord[1]} {coord[2]} {num_phys}{phys_str}\n')

		# Curves
		for curve_id in sorted(curve_set):
			phys_tag = curve_physical.get(curve_id, 0)
			num_phys = 1 if phys_tag else 0
			phys_str = f' {phys_tag}' if phys_tag else ''
			if phys_tag and curve_id == phys_tag:
				# Synthetic entity from edge block: compute bbox from edge nodes
				edge_nodes = set()
				for edge_id in _get_block_elements(cubit, curve_id, "edge"):
					nodes = cubit.get_expanded_connectivity("edge", edge_id)
					edge_nodes.update(nodes)
				if edge_nodes:
					coords = [cubit.get_nodal_coordinates(n) for n in edge_nodes]
					minx = min(c[0] for c in coords)
					maxx = max(c[0] for c in coords)
					miny = min(c[1] for c in coords)
					maxy = max(c[1] for c in coords)
					minz = min(c[2] for c in coords)
					maxz = max(c[2] for c in coords)
				else:
					minx = miny = minz = maxx = maxy = maxz = 0
				fid.write(f'{curve_id} {minx} {miny} {minz} {maxx} {maxy} {maxz} {num_phys}{phys_str} 0\n')
			else:
				# Geometric curve
				try:
					bbox = cubit.get_bounding_box("curve", curve_id)
					minx, maxx = bbox[0], bbox[1]
					miny, maxy = bbox[3], bbox[4]
					minz, maxz = bbox[6], bbox[7]
					curve = cubit.curve(curve_id)
					vert_ids = [v.id() for v in curve.vertices()]
					num_bnd = len(vert_ids)
					vert_str = ' '.join(str(v) for v in vert_ids)
					fid.write(f'{curve_id} {minx} {miny} {minz} {maxx} {maxy} {maxz} {num_phys}{phys_str} {num_bnd} {vert_str}\n')
				except:
					fid.write(f'{curve_id} 0 0 0 0 0 0 {num_phys}{phys_str} 0\n')

		# Surfaces
		if surface_set:
			for surf_id in sorted(surface_set):
				try:
					bbox = cubit.get_bounding_box("surface", surf_id)
					minx, maxx = bbox[0], bbox[1]
					miny, maxy = bbox[3], bbox[4]
					minz, maxz = bbox[6], bbox[7]
					# Physical tag
					phys_tag = surface_physical.get(surf_id, 0)
					num_phys = 1 if phys_tag else 0
					phys_str = f' {phys_tag}' if phys_tag else ''
					# Get bounding curves
					try:
						surf = cubit.surface(surf_id)
						curve_ids = [c.id() for c in surf.curves()]
						num_bnd = len(curve_ids)
						curve_str = ' '.join(str(c) for c in curve_ids)
					except:
						num_bnd = 0
						curve_str = ''
					fid.write(f'{surf_id} {minx} {miny} {minz} {maxx} {maxy} {maxz} {num_phys}{phys_str} {num_bnd} {curve_str}\n')
				except:
					fid.write(f'{surf_id} 0 0 0 0 0 0 0 0\n')
		else:
			# Fallback: use blocks as surfaces
			for block_id in cubit.get_block_id_list():
				if len(_get_block_elements(cubit, block_id, "tri")) + len(_get_block_elements(cubit, block_id, "face")) > 0:
					fid.write(f'{block_id} 0 0 0 0 0 0 1 {block_id} 0\n')

		# Volumes
		if volume_set:
			for vol_id in sorted(volume_set):
				try:
					bbox = cubit.get_bounding_box("volume", vol_id)
					minx, maxx = bbox[0], bbox[1]
					miny, maxy = bbox[3], bbox[4]
					minz, maxz = bbox[6], bbox[7]
					# Physical tag
					phys_tag = volume_physical.get(vol_id, 0)
					num_phys = 1 if phys_tag else 0
					phys_str = f' {phys_tag}' if phys_tag else ''
					# Get bounding surfaces
					try:
						vol = cubit.volume(vol_id)
						surf_ids = [s.id() for s in vol.surfaces()]
						num_bnd = len(surf_ids)
						surf_str = ' '.join(str(s) for s in surf_ids)
					except:
						num_bnd = 0
						surf_str = ''
					fid.write(f'{vol_id} {minx} {miny} {minz} {maxx} {maxy} {maxz} {num_phys}{phys_str} {num_bnd} {surf_str}\n')
				except:
					fid.write(f'{vol_id} 0 0 0 0 0 0 0 0\n')
		else:
			# Fallback: use blocks as volumes
			for block_id in cubit.get_block_id_list():
				if len(_get_block_elements(cubit, block_id, "tet")) + len(_get_block_elements(cubit, block_id, "hex")) + \
				   len(_get_block_elements(cubit, block_id, "wedge")) + len(_get_block_elements(cubit, block_id, "pyramid")) > 0:
					fid.write(f'{block_id} 0 0 0 0 0 0 1 {block_id} 0\n')

		fid.write('$EndEntities\n')

		# ------------------------------------------------------------
		# $Nodes (v4 format: entity blocks)
		# ------------------------------------------------------------
		# For simplicity, put all nodes in one entity block
		fid.write('$Nodes\n')
		# numEntityBlocks numNodes minNodeTag maxNodeTag
		fid.write(f'1 {num_nodes} {min_node_tag} {max_node_tag}\n')
		# entityDim entityTag parametric numNodesInBlock
		entity_dim = 3 if has_3d_elements else 2
		entity_tag = 1
		fid.write(f'{entity_dim} {entity_tag} 0 {num_nodes}\n')
		# Node tags
		for node_id in sorted_nodes:
			fid.write(f'{node_id}\n')
		# Node coordinates
		for node_id in sorted_nodes:
			coord = cubit.get_nodal_coordinates(node_id)
			if actual_dim == "2D":
				fid.write(f'{coord[0]} {coord[1]} 0\n')
			else:
				fid.write(f'{coord[0]} {coord[1]} {coord[2]}\n')
		fid.write('$EndNodes\n')

		# ------------------------------------------------------------
		# $Elements (v4 format: entity blocks)
		# ------------------------------------------------------------

		# Helper function: get element nodes with optional normal flip for 2D
		def get_element_nodes(elem_type, elem_id, flip_normal=False):
			nodes = list(cubit.get_expanded_connectivity(elem_type, elem_id))
			if flip_normal:
				if elem_type == "tri":
					if len(nodes) == 3:
						nodes = [nodes[0], nodes[2], nodes[1]]
					elif len(nodes) == 6:
						nodes = [nodes[0], nodes[2], nodes[1], nodes[5], nodes[4], nodes[3]]
				elif elem_type in ["quad", "face"]:
					if len(nodes) == 4:
						nodes = [nodes[0], nodes[3], nodes[2], nodes[1]]
					elif len(nodes) == 8:
						nodes = [nodes[0], nodes[3], nodes[2], nodes[1], nodes[7], nodes[6], nodes[5], nodes[4]]
			return nodes

		# Check if normal flip is needed for 2D elements
		def needs_normal_flip(elem_type, elem_id):
			if actual_dim != "2D":
				return False
			try:
				owner = cubit.get_geometry_owner(elem_type, elem_id)
				if owner and "surface" in owner.lower():
					surf_id = int(owner.split()[1])
					normal = cubit.get_surface_normal(surf_id)
					return normal[2] < 0
			except:
				pass
			return False

		# Count elements per block for entity blocks
		block_elements = {}  # block_id -> {'tet': [...], 'hex': [...], ...}
		for block_id in cubit.get_block_id_list():
			block_elements[block_id] = {
				'tet': list(_get_block_elements(cubit, block_id, "tet")),
				'hex': list(_get_block_elements(cubit, block_id, "hex")),
				'wedge': list(_get_block_elements(cubit, block_id, "wedge")),
				'pyramid': list(_get_block_elements(cubit, block_id, "pyramid")),
				'tri': list(_get_block_elements(cubit, block_id, "tri")),
				'quad': list(_get_block_elements(cubit, block_id, "face")),
				'edge': list(_get_block_elements(cubit, block_id, "edge")),
				'node': list(_get_block_elements(cubit, block_id, "node")),
			}

		# Count total elements and entity blocks
		total_elements = len(all_tets) + len(all_hexes) + len(all_wedges) + len(all_pyramids) + \
						 len(all_tris) + len(all_quads) + len(all_edges) + len(all_points)

		# Count entity blocks (one per element type per block)
		num_entity_blocks = 0
		for block_id, elems in block_elements.items():
			for elem_type, elem_list in elems.items():
				if len(elem_list) > 0:
					num_entity_blocks += 1

		if total_elements == 0:
			min_elem_tag = 0
			max_elem_tag = 0
		else:
			min_elem_tag = 1
			max_elem_tag = total_elements

		fid.write('$Elements\n')
		fid.write(f'{num_entity_blocks} {total_elements} {min_elem_tag} {max_elem_tag}\n')

		element_tag = 0

		for block_id in cubit.get_block_id_list():
			elems = block_elements[block_id]

			# Tetrahedra
			if len(elems['tet']) > 0:
				tet_list = elems['tet']
				# Determine element type from first element
				sample_nodes = cubit.get_expanded_connectivity("tet", tet_list[0])
				if len(sample_nodes) == 4:
					gmsh_type = GMSH_TET4
				elif len(sample_nodes) == 10:
					gmsh_type = GMSH_TET10
				else:
					gmsh_type = GMSH_TET11

				fid.write(f'3 {block_id} {gmsh_type} {len(tet_list)}\n')
				for tet_id in tet_list:
					element_tag += 1
					nodes = cubit.get_expanded_connectivity("tet", tet_id)
					if len(nodes) == 4:
						node_str = ' '.join(str(n) for n in nodes)
					elif len(nodes) == 10:
						# Cubit -> Gmsh node order for TET10
						reordered = [nodes[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]]
						node_str = ' '.join(str(n) for n in reordered)
					else:  # TET11
						reordered = [nodes[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 10]]
						node_str = ' '.join(str(n) for n in reordered)
					fid.write(f'{element_tag} {node_str}\n')

			# Hexahedra
			if len(elems['hex']) > 0:
				hex_list = elems['hex']
				sample_nodes = cubit.get_expanded_connectivity("hex", hex_list[0])
				gmsh_type = GMSH_HEX20 if len(sample_nodes) == 20 else GMSH_HEX8

				fid.write(f'3 {block_id} {gmsh_type} {len(hex_list)}\n')
				for hex_id in hex_list:
					element_tag += 1
					nodes = cubit.get_expanded_connectivity("hex", hex_id)
					if len(nodes) == 8:
						node_str = ' '.join(str(n) for n in nodes)
					else:  # HEX20
						# Cubit -> Gmsh node order for HEX20
						reordered = [nodes[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 9, 13, 10, 14, 15, 16, 19, 17, 18]]
						node_str = ' '.join(str(n) for n in reordered)
					fid.write(f'{element_tag} {node_str}\n')

			# Wedges
			if len(elems['wedge']) > 0:
				wedge_list = elems['wedge']
				sample_nodes = cubit.get_expanded_connectivity("wedge", wedge_list[0])
				gmsh_type = GMSH_WEDGE15 if len(sample_nodes) == 15 else GMSH_WEDGE6

				fid.write(f'3 {block_id} {gmsh_type} {len(wedge_list)}\n')
				for wedge_id in wedge_list:
					element_tag += 1
					nodes = cubit.get_expanded_connectivity("wedge", wedge_id)
					if len(nodes) == 6:
						node_str = ' '.join(str(n) for n in nodes)
					else:  # WEDGE15
						reordered = [nodes[i] for i in [0, 1, 2, 3, 4, 5, 6, 8, 9, 7, 10, 11, 12, 14, 13]]
						node_str = ' '.join(str(n) for n in reordered)
					fid.write(f'{element_tag} {node_str}\n')

			# Pyramids
			if len(elems['pyramid']) > 0:
				pyramid_list = elems['pyramid']
				sample_nodes = cubit.get_expanded_connectivity("pyramid", pyramid_list[0])
				gmsh_type = GMSH_PYRAMID13 if len(sample_nodes) == 13 else GMSH_PYRAMID5

				fid.write(f'3 {block_id} {gmsh_type} {len(pyramid_list)}\n')
				for pyramid_id in pyramid_list:
					element_tag += 1
					nodes = cubit.get_expanded_connectivity("pyramid", pyramid_id)
					if len(nodes) == 5:
						node_str = ' '.join(str(n) for n in nodes)
					else:  # PYRAMID13
						reordered = [nodes[i] for i in [0, 1, 2, 3, 4, 5, 8, 9, 6, 10, 7, 11, 12]]
						node_str = ' '.join(str(n) for n in reordered)
					fid.write(f'{element_tag} {node_str}\n')

			# Triangles
			if len(elems['tri']) > 0:
				tri_list = elems['tri']
				sample_nodes = cubit.get_expanded_connectivity("tri", tri_list[0])
				if len(sample_nodes) == 3:
					gmsh_type = GMSH_TRI3
				elif len(sample_nodes) == 6:
					gmsh_type = GMSH_TRI6
				else:
					gmsh_type = GMSH_TRI7

				fid.write(f'2 {block_id} {gmsh_type} {len(tri_list)}\n')
				for tri_id in tri_list:
					element_tag += 1
					flip = needs_normal_flip("tri", tri_id)
					nodes = get_element_nodes("tri", tri_id, flip)
					node_str = ' '.join(str(n) for n in nodes)
					fid.write(f'{element_tag} {node_str}\n')

			# Quads
			if len(elems['quad']) > 0:
				quad_list = elems['quad']
				sample_nodes = cubit.get_expanded_connectivity("quad", quad_list[0])
				if len(sample_nodes) == 4:
					gmsh_type = GMSH_QUAD4
				elif len(sample_nodes) == 8:
					gmsh_type = GMSH_QUAD8
				else:
					gmsh_type = GMSH_QUAD9

				fid.write(f'2 {block_id} {gmsh_type} {len(quad_list)}\n')
				for quad_id in quad_list:
					element_tag += 1
					flip = needs_normal_flip("quad", quad_id)
					nodes = get_element_nodes("quad", quad_id, flip)
					node_str = ' '.join(str(n) for n in nodes)
					fid.write(f'{element_tag} {node_str}\n')

			# Edges
			if len(elems['edge']) > 0:
				edge_list = elems['edge']
				sample_nodes = cubit.get_expanded_connectivity("edge", edge_list[0])
				gmsh_type = GMSH_LINE3 if len(sample_nodes) == 3 else GMSH_LINE2

				fid.write(f'1 {block_id} {gmsh_type} {len(edge_list)}\n')
				for edge_id in edge_list:
					element_tag += 1
					nodes = cubit.get_expanded_connectivity("edge", edge_id)
					node_str = ' '.join(str(n) for n in nodes)
					fid.write(f'{element_tag} {node_str}\n')

			# Points (nodes as elements)
			if len(elems['node']) > 0:
				node_list = elems['node']
				fid.write(f'0 {block_id} {GMSH_POINT} {len(node_list)}\n')
				for node_id in node_list:
					element_tag += 1
					fid.write(f'{element_tag} {node_id}\n')

		fid.write('$EndElements\n')

	return cubit

########################################################################
###	Nastran file
########################################################################

def export_Nastran(cubit: Any, FileName: str, DIM: str = "3D", PYRAM: bool = True) -> Any:
	"""Export mesh to Nastran format.

	Exports mesh elements from Cubit to NX Nastran bulk data format.
	Supports 2D and 3D meshes with first-order elements only.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .nas or .bdf file
		DIM: Dimension mode - "2D" for 2D meshes, "3D" for 3D meshes (default: "3D")
		PYRAM: If True, export pyramids as CPYRAM; if False, convert to degenerate hex (default: True)

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 3D: CTETRA (tet4), CHEXA (hex8), CPENTA (wedge6), CPYRAM (pyramid5)
		- 2D: CTRIA3 (tri3), CQUAD4 (quad4)
		- 1D: CROD (edge/bar)
		- 0D: CMASS (point mass)

	Note:
		Only first-order elements are supported. Second-order elements are not exported.
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	import datetime
	formatted_date_time = datetime.datetime.now().strftime("%d-%b-%y at %H:%M:%S")
	with open(FileName,'w',encoding='UTF-8') as fid:
		fid.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
		fid.write("$\n")
		fid.write("$                         CUBIT NX Nastran Translator\n")
		fid.write("$\n")
		fid.write(f"$            File: {FileName}\n")
		fid.write(f"$      Time Stamp: {formatted_date_time}\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$                        PLEASE CHECK YOUR MODEL FOR UNITS CONSISTENCY.\n")
		fid.write("$\n")
		fid.write("$       It should be noted that load ID's from CUBIT may NOT correspond to Nastran SID's\n")
		fid.write("$ The SID's for the load and restraint sets start at one and increment by one:i.e.,1,2,3,4...\n")
		fid.write("$\n")
		fid.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$ -------------------------\n")
		fid.write("$ Executive Control Section\n")
		fid.write("$ -------------------------\n")
		fid.write("$\n")
		fid.write("SOL 101\n")
		fid.write("CEND\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$ --------------------\n")
		fid.write("$ Case Control Section\n")
		fid.write("$ --------------------\n")
		fid.write("$\n")
		fid.write("ECHO = SORT\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$ Name: Initial\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$ Name: Default Set\n")
		fid.write("$\n")
		fid.write("SUBCASE = 1\n")
		fid.write("$\n")
		fid.write("LABEL = Default Set\n")
		fid.write("$\n")
		fid.write("$ -----------------\n")
		fid.write("$ Bulk Data Section\n")
		fid.write("$ -----------------\n")
		fid.write("$\n")
		fid.write("BEGIN BULK\n")
		fid.write("$\n")
		fid.write("$ Params\n")
		fid.write("$\n")
		fid.write("$\n")
		fid.write("$ Node cards\n")
		fid.write("$\n")

		node_list = set()
		for block_id in cubit.get_block_id_list():
			elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
			for elem_type in elem_types:
				for element_id in _get_block_elements(cubit, block_id, elem_type):
					node_ids = cubit.get_connectivity(elem_type, element_id)
					node_list.update(node_ids)
		for node_id in node_list:
			coord = cubit.get_nodal_coordinates(node_id)
			if DIM == "3D":
				fid.write(f"GRID*   {node_id:>16}{0:>16}{coord[0]:>16.5f}{coord[1]:>16.5f}\n*       {coord[2]:>16.5f}\n")
			else:
				fid.write(f"GRID*   {node_id:>16}{0:>16}{coord[0]:>16.5f}{coord[1]:>16.5f}\n*       {0}\n")

		element_id = 0
		fid.write("$\n")
		fid.write("$ Element cards\n")
		fid.write("$\n")
		for block_id in cubit.get_block_id_list():
			name = cubit.get_exodus_entity_name("block",block_id)
			fid.write("$\n")
			fid.write(f"$ Name: {name}\n")
			fid.write("$\n")
			tet_list = _get_block_elements(cubit, block_id, "tet")

			if DIM=="3D":
				for tet_id in tet_list:
					node_list = cubit.get_connectivity('tet',tet_id)
					element_id += 1
					fid.write(f"CTETRA  {element_id:>8}{block_id:>8}{node_list[0]:>8}{node_list[1]:>8}{node_list[2]:>8}{node_list[3]:>8}\n")
				hex_list = _get_block_elements(cubit, block_id, "hex")
				for hex_id in hex_list:
					node_list = cubit.get_connectivity('hex',hex_id)
					element_id += 1
					fid.write(f"CHEXA   {element_id:>8}{block_id:>8}{node_list[0]:>8}{node_list[1]:>8}{node_list[2]:>8}{node_list[3]:>8}{node_list[4]:>8}{node_list[5]:>8}+\n+       {node_list[6]:>8}{node_list[7]:>8}\n")
				wedge_list = _get_block_elements(cubit, block_id, "wedge")
				for wedge_id in wedge_list:
					node_list = cubit.get_connectivity('wedge',wedge_id)
					element_id += 1
					fid.write(f"CPENTA  {element_id:>8}{block_id:>8}{node_list[0]:>8}{node_list[1]:>8}{node_list[2]:>8}{node_list[3]:>8}{node_list[4]:>8}{node_list[5]:>8}\n")
				pyramid_list = _get_block_elements(cubit, block_id, "pyramid")
				for pyramid_id in pyramid_list:
					node_list = cubit.get_connectivity('pyramid',pyramid_id)
					if PYRAM:
						element_id += 1
						fid.write(f"CPYRAM  {element_id:>8}{block_id:>8}{node_list[0]:>8}{node_list[1]:>8}{node_list[2]:>8}{node_list[3]:>8}{node_list[4]:>8}\n")
					else:
						element_id += 1
						fid.write(f"CHEXA   {element_id:>8}{block_id:>8}{node_list[0]:>8}{node_list[1]:>8}{node_list[2]:>8}{node_list[3]:>8}{node_list[4]:>8}{node_list[4]:>8}+\n+       {node_list[4]:>8}{node_list[4]:>8}\n")

			tri_list = _get_block_elements(cubit, block_id, "tri")
			for tri_id in tri_list:
				node_list = cubit.get_connectivity('tri',tri_id)
				element_id += 1
				if DIM=="3D":
					fid.write(f"CTRIA3  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[1]:<8}{node_list[2]:<8}\n")
				else:
					surface_id = int(cubit.get_geometry_owner("tri", tri_id).split()[1])
					normal = cubit.get_surface_normal(surface_id)
					if normal[2] > 0:
						fid.write(f"CTRIA3  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[1]:<8}{node_list[2]:<8}\n")
					else:
						fid.write(f"CTRIA3  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[2]:<8}{node_list[1]:<8}\n")
			quad_list = _get_block_elements(cubit, block_id, "face")
			for quad_id in quad_list:
				node_list = cubit.get_connectivity('quad',quad_id)
				element_id += 1
				if DIM=="3D":
					fid.write(f"CQUAD4  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[1]:<8}{node_list[2]:<8}{node_list[3]:<8}\n")
				else:
					surface_id = int(cubit.get_geometry_owner("quad", quad_id).split()[1])
					normal = cubit.get_surface_normal(surface_id)
					node_list = cubit.get_connectivity('quad',quad_id)
					if normal[2] > 0:
						fid.write(f"CQUAD4  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[1]:<8}{node_list[2]:<8}{node_list[3]:<8}\n")
					else:
						fid.write(f"CQUAD4  {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[3]:<8}{node_list[2]:<8}{node_list[1]:<8}\n")
			edge_list = _get_block_elements(cubit, block_id, "edge")
			for edge_id in edge_list:
				element_id += 1
				node_list = cubit.get_connectivity('edge', edge_id)
				fid.write(f"CROD    {element_id:<8}{block_id:<8}{node_list[0]:<8}{node_list[1]:<8}\n")
			node_list = _get_block_elements(cubit, block_id, "node")
			for node_id in node_list:
				element_id += 1
				fid.write(f"CMASS   {element_id:<8}{block_id:<8}{node_id:<8}\n")
		fid.write("$\n")
		fid.write("$ Property cards\n")
		fid.write("$\n")

		for block_id in cubit.get_block_id_list():
			name = cubit.get_exodus_entity_name("block",block_id)
			fid.write("$\n")
			fid.write(f"$ Name: {name}\n")
			if len(_get_block_elements(cubit, block_id, "node")) > 0:
				fid.write(f"PMASS    {block_id:< 8}{block_id:< 8}\n")
			elif len(_get_block_elements(cubit, block_id, "edge")) > 0:
				fid.write(f"PROD     {block_id:< 8}{block_id:< 8}\n")
			elif len(_get_block_elements(cubit, block_id, "tri")) + len(_get_block_elements(cubit, block_id, "face")) > 0:
				fid.write(f"PSHELL   {block_id:< 8}{block_id:< 8}\n")
			else:
				fid.write(f"PSOLID   {block_id:< 8}{block_id:< 8}\n")
			fid.write("$\n")

		fid.write("ENDDATA\n")
	return cubit

########################################################################
###	ELF meg file
########################################################################

def export_meg(cubit: Any, FileName: str, DIM: str = 'T', MGR2: Optional[List[List[float]]] = None) -> Any:
	"""Export mesh to ELF/MESH MEG format.

	Exports mesh elements from Cubit to ELF/MAGIC MEG format for finite element analysis.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .meg file
		DIM: Dimension mode - 'T' for 3D, 'R' for axisymmetric, 'K' for 2D (default: 'T')
		MGR2: Optional list of spatial nodes [[x1,y1,z1], [x2,y2,z2], ...] (default: None)

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 3D: Tetrahedron, Hexahedron, Wedge, Pyramid
		- 2D: Triangle, Quadrilateral
		- 1D: Edge
		- 0D: Node
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	if MGR2 is None:
		MGR2 = []

	with open(FileName,'w',encoding='UTF-8') as fid:
		fid.write("BOOK  MEP  3.50\n")
		fid.write("* ELF/MESH VERSION 7.3.0\n")
		fid.write("* SOLVER = ELF/MAGIC\n")
		fid.write("MGSC 0.001\n")
		fid.write("* NODE\n")

		node_list = set()
		for block_id in cubit.get_block_id_list():
			elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge"]
			for elem_type in elem_types:
				for element_id in _get_block_elements(cubit, block_id, elem_type):
					node_ids = cubit.get_connectivity(elem_type, element_id)
					node_list.update(node_ids)
		for node_id in node_list:
			coord = cubit.get_nodal_coordinates(node_id)
			if DIM=='T':
				fid.write(f"MGR1 {node_id} 0 {coord[0]} {coord[1]} {coord[2]}\n")
			if DIM=='K':
				fid.write(f"MGR1 {node_id} 0 {coord[0]} {coord[1]} {0}\n")
			if DIM=='R':
				fid.write(f"MGR1 {node_id} 0 {coord[0]} {0} {coord[2]}\n")

		element_id = 0
		fid.write("* ELEMENT K\n")
		for block_id in cubit.get_block_id_list():
			name = cubit.get_exodus_entity_name("block",block_id)

			if DIM=='T':
				tet_list = _get_block_elements(cubit, block_id, "tet")
				for tet_id in tet_list:
					node_list = cubit.get_connectivity('tet',tet_id)
					element_id += 1
					fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]}\n")

				hex_list = _get_block_elements(cubit, block_id, "hex")
				for hex_id in hex_list:
					node_list = cubit.get_connectivity('hex',hex_id)
					element_id += 1
					fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]} {node_list[6]} {node_list[7]}\n")

				wedge_list = _get_block_elements(cubit, block_id, "wedge")
				for wedge_id in wedge_list:
					node_list = cubit.get_connectivity('wedge',wedge_id)
					element_id += 1
					fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[5]}\n")

				pyramid_list = _get_block_elements(cubit, block_id, "pyramid")
				for pyramid_id in pyramid_list:
					node_list = cubit.get_connectivity('pyramid',pyramid_id)
					element_id += 1
					fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]} {node_list[4]} {node_list[4]} {node_list[4]} {node_list[4]}\n")

			tri_list = _get_block_elements(cubit, block_id, "tri")
			for tri_id in tri_list:
				node_list = cubit.get_connectivity('tri',tri_id)
				element_id += 1
				fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]}\n")

			quad_list = _get_block_elements(cubit, block_id, "face")
			for quad_id in quad_list:
				node_list = cubit.get_connectivity('quad',quad_id)
				element_id += 1
				fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]} {node_list[2]} {node_list[3]}\n")

			edge_list = _get_block_elements(cubit, block_id, "edge")
			for edge_id in edge_list:
				node_list = cubit.get_connectivity('edge',edge_id)
				element_id += 1
				fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_list[0]} {node_list[1]}\n")

			node_list = _get_block_elements(cubit, block_id, "node")
			for node_id in node_list:
				element_id += 1
				fid.write(f"{name[0:4]}{DIM} {element_id} 0 {block_id} {node_id}\n")

		fid.write("* NODE\n")
		for node_id in range(len(MGR2)):
			fid.write(f"MGR2 {node_id+1} 0 {MGR2[node_id][0]} {MGR2[node_id][1]} {MGR2[node_id][2]}\n")
		fid.write("BOOK  END\n")
	return cubit

########################################################################
###	vtk format
########################################################################

def export_vtk(cubit: Any, FileName: str) -> Any:
	"""Export mesh to Legacy VTK format.

	Exports mesh elements from Cubit to VTK (Visualization Toolkit) Legacy format.
	Automatically detects element order (1st or 2nd) based on node count.
	Supports mixed-order meshes (1st and 2nd order elements in same file).

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .vtk file

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 1st order: Tet4, Hex8, Wedge6, Pyramid5, Triangle3, Quad4, Line2, Point
		- 2nd order: Tet10, Hex20, Wedge15, Pyramid13, Triangle6, Quad8, Line3

	Note:
		VTK cell data includes scalar values to distinguish element types.
		Element order is automatically detected from node count per element.
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	with open(FileName,'w') as fid:
		fid.write('# vtk DataFile Version 3.0\n')
		fid.write(f'Unstructured Grid {FileName}\n')
		fid.write('ASCII\n')
		fid.write('DATASET UNSTRUCTURED_GRID\n')
		# First, collect all unique node IDs from all elements
		node_list = set()
		hex_list = set()
		tet_list = set()
		wedge_list = set()
		pyramid_list = set()
		tri_list = set()
		quad_list = set()
		edge_list = set()
		nodes_list = set()

		for block_id in cubit.get_block_id_list():
			tet_list.update(_get_block_elements(cubit, block_id, "tet"))
			hex_list.update(_get_block_elements(cubit, block_id, "hex"))
			wedge_list.update(_get_block_elements(cubit, block_id, "wedge"))
			pyramid_list.update(_get_block_elements(cubit, block_id, "pyramid"))
			tri_list.update(_get_block_elements(cubit, block_id, "tri"))
			quad_list.update(_get_block_elements(cubit, block_id, "face"))
			edge_list.update(_get_block_elements(cubit, block_id, "edge"))
			nodes_list.update(_get_block_elements(cubit, block_id, "node"))

		# Collect all node IDs from all element types
		elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
		for block_id in cubit.get_block_id_list():
			for elem_type in elem_types:
				for element_id in _get_block_elements(cubit, block_id, elem_type):
					node_ids = cubit.get_expanded_connectivity(elem_type, element_id)
					node_list.update(node_ids)

		# Write POINTS header
		fid.write(f'POINTS {len(node_list)} float\n')

		# Create mapping from Cubit node ID to VTK index (0-indexed)
		# Write coordinates in sorted order for consistency
		node_id_to_vtk_index = {}
		vtk_index = 0
		for node_id in sorted(node_list):
			coord = cubit.get_nodal_coordinates(node_id)
			fid.write(f'{coord[0]} {coord[1]} {coord[2]}\n')
			node_id_to_vtk_index[node_id] = vtk_index
			vtk_index += 1

		# Pre-scan all elements to get node counts for each element
		# This enables auto-detection of element order and mixed-order support
		tet_node_counts = {tet_id: len(cubit.get_expanded_connectivity("tet", tet_id)) for tet_id in tet_list}
		hex_node_counts = {hex_id: len(cubit.get_expanded_connectivity("hex", hex_id)) for hex_id in hex_list}
		wedge_node_counts = {wedge_id: len(cubit.get_expanded_connectivity("wedge", wedge_id)) for wedge_id in wedge_list}
		pyramid_node_counts = {pyramid_id: len(cubit.get_expanded_connectivity("pyramid", pyramid_id)) for pyramid_id in pyramid_list}
		tri_node_counts = {tri_id: len(cubit.get_expanded_connectivity("tri", tri_id)) for tri_id in tri_list}
		quad_node_counts = {quad_id: len(cubit.get_expanded_connectivity("quad", quad_id)) for quad_id in quad_list}
		edge_node_counts = {edge_id: len(cubit.get_expanded_connectivity("edge", edge_id)) for edge_id in edge_list}

		# Calculate total CELLS size: sum of (1 + node_count) for each element
		num_cells = len(tet_list) + len(hex_list) + len(wedge_list) + len(pyramid_list) + len(tri_list) + len(quad_list) + len(edge_list) + len(nodes_list)
		cells_size = (
			sum(1 + n for n in tet_node_counts.values()) +
			sum(1 + n for n in hex_node_counts.values()) +
			sum(1 + n for n in wedge_node_counts.values()) +
			sum(1 + n for n in pyramid_node_counts.values()) +
			sum(1 + n for n in tri_node_counts.values()) +
			sum(1 + n for n in quad_node_counts.values()) +
			sum(1 + n for n in edge_node_counts.values()) +
			2 * len(nodes_list)  # Each node: 1 (count) + 1 (node_id)
		)
		fid.write(f'CELLS {num_cells} {cells_size}\n')

		for tet_id in tet_list:
			node_list = cubit.get_expanded_connectivity("tet", tet_id)
			if len(node_list)==4:
				fid.write(f'4 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]}\n')
			elif len(node_list)==10:
				fid.write(f'10 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]} {node_id_to_vtk_index[node_list[8]]} {node_id_to_vtk_index[node_list[9]]}\n')
		for hex_id in hex_list:
			node_list = cubit.get_expanded_connectivity("hex", hex_id)
			if len(node_list)==8:
				fid.write(f'8 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]}\n')
			elif len(node_list)==20:
				fid.write(f'20 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]} {node_id_to_vtk_index[node_list[8]]} {node_id_to_vtk_index[node_list[9]]} {node_id_to_vtk_index[node_list[10]]} {node_id_to_vtk_index[node_list[11]]} {node_id_to_vtk_index[node_list[16]]} {node_id_to_vtk_index[node_list[17]]} {node_id_to_vtk_index[node_list[18]]} {node_id_to_vtk_index[node_list[19]]} {node_id_to_vtk_index[node_list[12]]} {node_id_to_vtk_index[node_list[13]]} {node_id_to_vtk_index[node_list[14]]} {node_id_to_vtk_index[node_list[15]]}\n')
		for wedge_id in wedge_list:
			node_list = cubit.get_expanded_connectivity("wedge", wedge_id)
			if len(node_list)==6:
				fid.write(f'6 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} \n')
			elif len(node_list)==15:
				fid.write(f'15 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]} {node_id_to_vtk_index[node_list[8]]} {node_id_to_vtk_index[node_list[12]]} {node_id_to_vtk_index[node_list[13]]} {node_id_to_vtk_index[node_list[14]]} {node_id_to_vtk_index[node_list[9]]} {node_id_to_vtk_index[node_list[10]]} {node_id_to_vtk_index[node_list[11]]} \n')

		for pyramid_id in pyramid_list:
			node_list = cubit.get_expanded_connectivity("pyramid", pyramid_id)
			if len(node_list)==5:
				fid.write(f'5 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} \n')
			elif len(node_list)==13:
				fid.write(f'13 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]} {node_id_to_vtk_index[node_list[8]]} {node_id_to_vtk_index[node_list[9]]} {node_id_to_vtk_index[node_list[10]]} {node_id_to_vtk_index[node_list[11]]} {node_id_to_vtk_index[node_list[12]]} \n')
		for tri_id in tri_list:
			node_list = cubit.get_expanded_connectivity("tri", tri_id)
			if len(node_list)==3:
				fid.write(f'3 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} \n')
			elif len(node_list)==6:
				fid.write(f'6 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} \n')
		for quad_id in quad_list:
			node_list = cubit.get_expanded_connectivity("quad", quad_id)
			if len(node_list)==4:
				fid.write(f'4 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} \n')
			elif len(node_list)==8:
				fid.write(f'8 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} {node_id_to_vtk_index[node_list[3]]} {node_id_to_vtk_index[node_list[4]]} {node_id_to_vtk_index[node_list[5]]} {node_id_to_vtk_index[node_list[6]]} {node_id_to_vtk_index[node_list[7]]}\n')
		for edge_id in edge_list:
			node_list = cubit.get_expanded_connectivity("edge", edge_id)
			if len(node_list)==2:
				fid.write(f'2 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} \n')
			elif len(node_list)==3:
				fid.write(f'3 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} {node_id_to_vtk_index[node_list[2]]} \n')
		for node_id in nodes_list:
			fid.write(f'1 {node_id_to_vtk_index[node_id]} \n')

		fid.write(f'CELL_TYPES {num_cells}\n')
		# VTK cell types based on node count (auto-detection)
		# 1st order: TET=10, HEX=12, WEDGE=13, PYRAMID=14, TRI=5, QUAD=9, LINE=3, POINT=1
		# 2nd order: TET=24, HEX=25, WEDGE=26, PYRAMID=27, TRI=22, QUAD=23, LINE=21
		for tet_id in tet_list:
			fid.write('24\n' if tet_node_counts[tet_id] == 10 else '10\n')
		for hex_id in hex_list:
			fid.write('25\n' if hex_node_counts[hex_id] == 20 else '12\n')
		for wedge_id in wedge_list:
			fid.write('26\n' if wedge_node_counts[wedge_id] == 15 else '13\n')
		for pyramid_id in pyramid_list:
			fid.write('27\n' if pyramid_node_counts[pyramid_id] == 13 else '14\n')
		for tri_id in tri_list:
			fid.write('22\n' if tri_node_counts[tri_id] == 6 else '5\n')
		for quad_id in quad_list:
			fid.write('23\n' if quad_node_counts[quad_id] == 8 else '9\n')
		for edge_id in edge_list:
			fid.write('21\n' if edge_node_counts[edge_id] == 3 else '3\n')
		for node_id in nodes_list:
			fid.write('1\n')
		fid.write(f'CELL_DATA {num_cells}\n')
		fid.write('SCALARS scalars float\n')
		fid.write('LOOKUP_TABLE default\n')
		for tet_id in tet_list:
			fid.write('1\n')
		for hex_id in hex_list:
			fid.write('2\n')
		for wedge_id in wedge_list:
			fid.write('3\n')
		for pyramid_id in pyramid_list:
			fid.write('4\n')
		for tri_id in tri_list:
			fid.write('5\n')
		for quad_id in quad_list:
			fid.write('6\n')
		for edge_id in edge_list:
			fid.write('0\n')
		for node_id in nodes_list:
			fid.write('-1\n')
	return cubit

########################################################################
###	NGSolve/Netgen format
########################################################################

def export_NetgenMesh(cubit: Any, geometry_file: Optional[str] = None, geometry: Any = None) -> Any:
	"""Export Cubit mesh to Netgen mesh format.

	Creates a netgen.meshing.Mesh object directly from Cubit mesh data.
	When geometry is provided (file or OCC object), the mesh can be curved
	using ngsolve.Mesh.Curve(order) for high-order geometry approximation.

	Note: netgen.meshing.Mesh is the core mesh data structure.
	      ngsolve.Mesh is a wrapper/view for FEM analysis.

	Args:
		cubit: Cubit Python interface object
		geometry_file: Path to geometry file (.step, .stp, .brep, .iges) for
		               mesh.Curve() support. If None, mesh is created without
		               geometry reference (Curve() will not work).
		geometry: An netgen.occ.OCCGeometry object. If provided, this takes
		          precedence over geometry_file. Useful for avoiding seam
		          issues with cylindrical surfaces (see note below).

	Returns:
		netgen.meshing.Mesh: Netgen mesh object ready for use with NGSolve

	Supported elements:
		- 3D: Tetrahedron (4-node), Hexahedron (8-node), Wedge (6-node), Pyramid (5-node)
		- 2D boundary: Triangle (3-node), Quadrilateral (4-node)
		- 1D boundary: Edge (2-node)

	Usage:
		# Basic usage with STEP file
		ngmesh = export_NetgenMesh(cubit, geometry_file="geometry.step")
		mesh = ngsolve.Mesh(ngmesh)
		mesh.Curve(3)

		# Using OCC geometry directly (recommended for cylinders)
		from netgen import occ
		from netgen.occ import OCCGeometry
		cyl = occ.Cylinder(occ.Pnt(0,0,-1), occ.Vec(0,0,1), radius=0.5, height=2)
		geo = OCCGeometry(cyl)
		ngmesh = export_NetgenMesh(cubit, geometry=geo)
		mesh = ngsolve.Mesh(ngmesh)
		mesh.Curve(3)

		# Without geometry (no Curve support)
		ngmesh = export_NetgenMesh(cubit)
		mesh = ngsolve.Mesh(ngmesh)

	Note:
		- Only 1st order elements are transferred; use mesh.Curve(order) for
		  high-order geometry approximation
		- The geometry should match the geometry used to create the mesh
		- Block names are used as Materials/Boundaries in NGSolve

	Cylindrical Surface Note:
		When using STEP files exported from Cubit with cylindrical surfaces,
		OCC may split the surface at the seam line (y=0), creating 2 faces
		where Cubit has 1 surface. Mesh elements crossing this seam may not
		curve correctly. To avoid this:
		1. Use OCC primitives (occ.Cylinder) via the geometry parameter, or
		2. Ensure mesh elements don't cross y=0 on cylindrical surfaces
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	from netgen.meshing import Mesh, MeshPoint, Element3D, Element2D, Element1D, FaceDescriptor
	from netgen.csg import Pnt

	# ============================================================
	# Node ordering conversion tables: Cubit -> Netgen
	# ============================================================

	# Tetrahedron: Cubit [0,1,2,3] -> Netgen [0,1,2,3]
	TET_ORDERING = [0, 1, 2, 3]

	# Hexahedron: Cubit [0,1,2,3,4,5,6,7] -> Netgen [0,1,5,4,3,2,6,7]
	HEX_ORDERING = [0, 1, 5, 4, 3, 2, 6, 7]

	# Wedge/Prism: Cubit [0,1,2,3,4,5] -> Netgen [0,2,1,3,5,4]
	WEDGE_ORDERING = [0, 2, 1, 3, 5, 4]

	# Pyramid: Cubit [0,1,2,3,4] -> Netgen [3,2,1,0,4]
	PYRAMID_ORDERING = [3, 2, 1, 0, 4]

	# Triangle: Cubit [0,1,2] -> Netgen [0,1,2]
	TRI_ORDERING = [0, 1, 2]

	# Quadrilateral: Cubit [0,1,2,3] -> Netgen [0,1,2,3]
	QUAD_ORDERING = [0, 1, 2, 3]

	# ============================================================
	# Create netgen mesh
	# ============================================================

	ngmesh = Mesh(dim=3)

	# Load and attach geometry if provided
	# Priority: geometry parameter > geometry_file parameter
	geo = None
	if geometry is not None:
		# Use provided OCC geometry object directly
		geo = geometry
		ngmesh.SetGeometry(geo)
	elif geometry_file is not None:
		from netgen.occ import OCCGeometry
		geo = OCCGeometry(geometry_file)
		ngmesh.SetGeometry(geo)

	# ============================================================
	# Collect all nodes from blocks
	# ============================================================

	node_map = {}  # Cubit node ID -> netgen point index
	all_nodes = set()

	for block_id in cubit.get_block_id_list():
		elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
		for elem_type in elem_types:
			for element_id in _get_block_elements(cubit, block_id, elem_type):
				# Use get_connectivity for 1st order nodes only
				node_ids = cubit.get_connectivity(elem_type, element_id)
				all_nodes.update(node_ids)

	# Add nodes to netgen mesh
	for node_id in sorted(all_nodes):
		coord = cubit.get_nodal_coordinates(node_id)
		pnt_idx = ngmesh.Add(MeshPoint(Pnt(coord[0], coord[1], coord[2])))
		node_map[node_id] = pnt_idx

	# ============================================================
	# Add 3D volume elements
	# ============================================================

	material_index = 0

	for block_id in cubit.get_block_id_list():
		block_name = cubit.get_exodus_entity_name("block", block_id)

		tet_list = _get_block_elements(cubit, block_id, "tet")
		hex_list = _get_block_elements(cubit, block_id, "hex")
		wedge_list = _get_block_elements(cubit, block_id, "wedge")
		pyramid_list = _get_block_elements(cubit, block_id, "pyramid")

		# Only process blocks with 3D elements
		if len(tet_list) + len(hex_list) + len(wedge_list) + len(pyramid_list) > 0:
			material_index += 1
			ngmesh.SetMaterial(material_index, block_name)

			# Add tetrahedra
			for tet_id in tet_list:
				nodes = cubit.get_connectivity("tet", tet_id)
				ng_nodes = [node_map[nodes[i]] for i in TET_ORDERING]
				ngmesh.Add(Element3D(material_index, ng_nodes))

			# Add hexahedra
			for hex_id in hex_list:
				nodes = cubit.get_connectivity("hex", hex_id)
				ng_nodes = [node_map[nodes[i]] for i in HEX_ORDERING]
				ngmesh.Add(Element3D(material_index, ng_nodes))

			# Add wedges/prisms
			for wedge_id in wedge_list:
				nodes = cubit.get_connectivity("wedge", wedge_id)
				ng_nodes = [node_map[nodes[i]] for i in WEDGE_ORDERING]
				ngmesh.Add(Element3D(material_index, ng_nodes))

			# Add pyramids
			for pyramid_id in pyramid_list:
				nodes = cubit.get_connectivity("pyramid", pyramid_id)
				ng_nodes = [node_map[nodes[i]] for i in PYRAMID_ORDERING]
				ngmesh.Add(Element3D(material_index, ng_nodes))

	# ============================================================
	# Add 2D surface elements (boundary faces)
	# ============================================================

	fd_index = 0

	# Build OCC face mapping when geometry is provided (for Curve() support)
	occ_face_info = None
	elem_to_occ_face = {}  # (elem_type, elem_id) -> OCC face index

	if geo is not None:
		# Analyze OCC faces to build mapping
		occ_faces = list(geo.shape.faces)
		occ_face_info = []
		for i, face in enumerate(occ_faces):
			center = face.center
			bb = face.bounding_box
			# bb is ((xmin, ymin, zmin), (xmax, ymax, zmax))
			z_range = bb[1][2] - bb[0][2]
			y_min, y_max = bb[0][1], bb[1][1]
			occ_face_info.append({
				'index': i,
				'center': (center.x, center.y, center.z),
				'bbox': bb,
				'z_range': z_range,
				'y_min': y_min,
				'y_max': y_max,
				'is_planar': z_range < 0.01,
			})

		# Build tri_id -> Cubit surface_id mapping
		tri_to_surface = {}
		surface_ids = cubit.get_entities("surface")
		for sid in surface_ids:
			for tri_id in cubit.get_surface_tris(sid):
				tri_to_surface[tri_id] = sid

		# Build quad_id -> Cubit surface_id mapping
		quad_to_surface = {}
		for sid in surface_ids:
			for quad_id in cubit.get_surface_quads(sid):
				quad_to_surface[quad_id] = sid

		# Helper function to get element centroid
		def get_elem_centroid(elem_type, elem_id):
			nodes = cubit.get_connectivity(elem_type, elem_id)
			coords = [cubit.get_nodal_coordinates(n) for n in nodes]
			return [sum(c[i] for c in coords) / len(coords) for i in range(3)]

		# Helper function to check if element crosses y=0 seam
		def crosses_seam(elem_type, elem_id, tolerance=0.001):
			nodes = cubit.get_connectivity(elem_type, elem_id)
			coords = [cubit.get_nodal_coordinates(n) for n in nodes]
			y_vals = [c[1] for c in coords]
			y_min, y_max = min(y_vals), max(y_vals)
			# Element crosses seam if vertices are on both sides of y=0
			return y_min < -tolerance and y_max > tolerance

		# Map Cubit surfaces to OCC faces based on centroid
		surface_to_occ = {}
		for sid in surface_ids:
			surf_centroid = cubit.get_surface_centroid(sid)
			surf_type = cubit.get_surface_type(sid).lower()

			# For planar surfaces, match by z-coordinate
			if "plane" in surf_type:
				for info in occ_face_info:
					if info['is_planar']:
						if abs(info['center'][2] - surf_centroid[2]) < 0.1:
							surface_to_occ[sid] = [info['index']]
							break
			else:
				# For curved surfaces (cone, cylinder, etc.), find all matching OCC faces
				# OCC may split curved surfaces into multiple faces
				matching = []
				for info in occ_face_info:
					if not info['is_planar']:
						# Check if centroid z-range overlaps
						if (info['bbox'][0][2] <= surf_centroid[2] <= info['bbox'][1][2]):
							matching.append(info['index'])
				if matching:
					surface_to_occ[sid] = matching

		# Map each triangle to its OCC face
		seam_crossing_count = 0
		for tri_id, sid in tri_to_surface.items():
			if sid not in surface_to_occ:
				continue
			occ_faces_for_surface = surface_to_occ[sid]
			if len(occ_faces_for_surface) == 1:
				elem_to_occ_face[('tri', tri_id)] = occ_faces_for_surface[0]
			else:
				# Multiple OCC faces for this surface (e.g., cylinder split)
				# Check if element crosses the seam (y=0)
				if crosses_seam("tri", tri_id):
					# Element crosses seam - mark as None to exclude from Curve()
					elem_to_occ_face[('tri', tri_id)] = None
					seam_crossing_count += 1
				else:
					# Determine by y-coordinate of triangle centroid
					centroid = get_elem_centroid("tri", tri_id)
					for occ_idx in occ_faces_for_surface:
						info = occ_face_info[occ_idx]
						# Check if centroid y is within OCC face bbox
						if info['y_min'] - 0.01 <= centroid[1] <= info['y_max'] + 0.01:
							elem_to_occ_face[('tri', tri_id)] = occ_idx
							break

		# Map each quad to its OCC face
		for quad_id, sid in quad_to_surface.items():
			if sid not in surface_to_occ:
				continue
			occ_faces_for_surface = surface_to_occ[sid]
			if len(occ_faces_for_surface) == 1:
				elem_to_occ_face[('quad', quad_id)] = occ_faces_for_surface[0]
			else:
				# Check if element crosses the seam (y=0)
				if crosses_seam("face", quad_id):
					# Element crosses seam - mark as None to exclude from Curve()
					elem_to_occ_face[('quad', quad_id)] = None
					seam_crossing_count += 1
				else:
					centroid = get_elem_centroid("face", quad_id)
					for occ_idx in occ_faces_for_surface:
						info = occ_face_info[occ_idx]
						if info['y_min'] - 0.01 <= centroid[1] <= info['y_max'] + 0.01:
							elem_to_occ_face[('quad', quad_id)] = occ_idx
							break

		# Print warning if seam-crossing elements found
		if seam_crossing_count > 0:
			print(f"  Note: {seam_crossing_count} boundary element(s) cross y=0 seam - excluded from Curve()")
			print(f"  For better results, use OCC geometry primitives (geometry= parameter)")

	# Create FaceDescriptors and add boundary elements
	if occ_face_info is not None:
		# Create one FaceDescriptor per OCC face
		# surfnr is 1-indexed in Netgen (OCC Face i -> surfnr=i+1)
		occ_to_fd_index = {}
		for info in occ_face_info:
			fd_index += 1
			fd = FaceDescriptor(bc=fd_index, surfnr=info['index'] + 1)
			fd.bcname = f"face_{info['index']}"
			ngmesh.Add(fd)
			ngmesh.SetBCName(fd_index - 1, f"face_{info['index']}")
			occ_to_fd_index[info['index']] = fd_index

		# Create a special FaceDescriptor for seam-crossing elements (surfnr=0)
		fd_index += 1
		fd_seam = FaceDescriptor(bc=fd_index, surfnr=0)  # surfnr=0 -> no geometry curving
		fd_seam.bcname = "seam_crossing"
		ngmesh.Add(fd_seam)
		ngmesh.SetBCName(fd_index - 1, "seam_crossing")
		fd_seam_index = fd_index

		# Add boundary elements with correct FaceDescriptor
		for block_id in cubit.get_block_id_list():
			tri_list = _get_block_elements(cubit, block_id, "tri")
			quad_list = _get_block_elements(cubit, block_id, "face")

			for tri_id in tri_list:
				nodes = cubit.get_connectivity("tri", tri_id)
				if all(n in node_map for n in nodes):
					ng_nodes = [node_map[nodes[i]] for i in TRI_ORDERING]
					key = ('tri', tri_id)
					if key in elem_to_occ_face:
						occ_idx = elem_to_occ_face[key]
						if occ_idx is None:
							bc_idx = fd_seam_index  # Seam-crossing element
						else:
							bc_idx = occ_to_fd_index[occ_idx]
					else:
						bc_idx = 1  # Fallback
					ngmesh.Add(Element2D(bc_idx, ng_nodes))

			for quad_id in quad_list:
				nodes = cubit.get_connectivity("quad", quad_id)
				if all(n in node_map for n in nodes):
					ng_nodes = [node_map[nodes[i]] for i in QUAD_ORDERING]
					key = ('quad', quad_id)
					if key in elem_to_occ_face:
						occ_idx = elem_to_occ_face[key]
						if occ_idx is None:
							bc_idx = fd_seam_index  # Seam-crossing element
						else:
							bc_idx = occ_to_fd_index[occ_idx]
					else:
						bc_idx = 1  # Fallback
					ngmesh.Add(Element2D(bc_idx, ng_nodes))
	else:
		# No geometry file - use simple block-based FaceDescriptors
		for block_id in cubit.get_block_id_list():
			block_name = cubit.get_exodus_entity_name("block", block_id)

			tri_list = _get_block_elements(cubit, block_id, "tri")
			quad_list = _get_block_elements(cubit, block_id, "face")

			# Only process blocks with 2D elements
			if len(tri_list) + len(quad_list) > 0:
				fd_index += 1
				fd = FaceDescriptor(bc=fd_index, surfnr=fd_index)
				fd.bcname = block_name
				ngmesh.Add(fd)
				ngmesh.SetBCName(fd_index - 1, block_name)

				# Add triangles
				for tri_id in tri_list:
					nodes = cubit.get_connectivity("tri", tri_id)
					if all(n in node_map for n in nodes):
						ng_nodes = [node_map[nodes[i]] for i in TRI_ORDERING]
						ngmesh.Add(Element2D(fd_index, ng_nodes))

				# Add quadrilaterals
				for quad_id in quad_list:
					nodes = cubit.get_connectivity("quad", quad_id)
					if all(n in node_map for n in nodes):
						ng_nodes = [node_map[nodes[i]] for i in QUAD_ORDERING]
						ngmesh.Add(Element2D(fd_index, ng_nodes))

	# ============================================================
	# Add 1D edge elements (if any)
	# ============================================================

	for block_id in cubit.get_block_id_list():
		block_name = cubit.get_exodus_entity_name("block", block_id)
		edge_list = _get_block_elements(cubit, block_id, "edge")

		if len(edge_list) > 0:
			fd_index += 1
			fd = FaceDescriptor(bc=fd_index, surfnr=fd_index)
			fd.bcname = block_name
			ngmesh.Add(fd)

			for edge_id in edge_list:
				nodes = cubit.get_connectivity("edge", edge_id)
				if all(n in node_map for n in nodes):
					ng_nodes = [node_map[n] for n in nodes]
					ngmesh.Add(Element1D(ng_nodes, index=fd_index))

	return ngmesh

########################################################################
###	VTK XML format (VTU)
########################################################################

def export_vtu(cubit: Any, FileName: str, binary: bool = False) -> Any:
	"""Export mesh to VTK XML format (VTU - Unstructured Grid).

	Exports mesh elements from Cubit to VTK XML format (.vtu).
	This is the modern VTK format recommended for ParaView and other tools.
	Automatically detects element order (1st or 2nd) based on node count.

	Args:
		cubit: Cubit Python interface object
		FileName: Output file path for the .vtu file
		binary: If True, write binary data (appended format). Default: False (ASCII)

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 1st order: Tet4, Hex8, Wedge6, Pyramid5, Triangle3, Quad4, Line2, Point
		- 2nd order: Tet10, Hex20, Wedge15, Pyramid13, Triangle6, Quad8, Line3

	VTK XML vs Legacy:
		- XML format is more efficient and extensible
		- Supports compression (when binary=True)
		- Better metadata support
		- Recommended format for modern workflows

	Example:
		cubit_mesh_export.export_vtu(cubit, "mesh.vtu")
		cubit_mesh_export.export_vtu(cubit, "mesh.vtu", binary=True)  # Binary mode
	"""
	_warn_mixed_element_types_in_blocks(cubit)

	import base64
	import struct

	# VTK cell type codes
	VTK_VERTEX = 1
	VTK_LINE = 3
	VTK_QUADRATIC_EDGE = 21
	VTK_TRIANGLE = 5
	VTK_QUADRATIC_TRIANGLE = 22
	VTK_QUAD = 9
	VTK_QUADRATIC_QUAD = 23
	VTK_TETRA = 10
	VTK_QUADRATIC_TETRA = 24
	VTK_HEXAHEDRON = 12
	VTK_QUADRATIC_HEXAHEDRON = 25
	VTK_WEDGE = 13
	VTK_QUADRATIC_WEDGE = 26
	VTK_PYRAMID = 14
	VTK_QUADRATIC_PYRAMID = 27

	# Collect all elements from blocks
	hex_list = set()
	tet_list = set()
	wedge_list = set()
	pyramid_list = set()
	tri_list = set()
	quad_list = set()
	edge_list = set()
	nodes_list = set()
	node_set = set()

	for block_id in cubit.get_block_id_list():
		tet_list.update(_get_block_elements(cubit, block_id, "tet"))
		hex_list.update(_get_block_elements(cubit, block_id, "hex"))
		wedge_list.update(_get_block_elements(cubit, block_id, "wedge"))
		pyramid_list.update(_get_block_elements(cubit, block_id, "pyramid"))
		tri_list.update(_get_block_elements(cubit, block_id, "tri"))
		quad_list.update(_get_block_elements(cubit, block_id, "face"))
		edge_list.update(_get_block_elements(cubit, block_id, "edge"))
		nodes_list.update(_get_block_elements(cubit, block_id, "node"))

	# Collect all node IDs
	elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
	for block_id in cubit.get_block_id_list():
		for elem_type in elem_types:
			for element_id in _get_block_elements(cubit, block_id, elem_type):
				node_ids = cubit.get_expanded_connectivity(elem_type, element_id)
				node_set.update(node_ids)

	# Create node ID to VTK index mapping
	node_id_to_vtk_index = {}
	vtk_index = 0
	sorted_nodes = sorted(node_set)
	for node_id in sorted_nodes:
		node_id_to_vtk_index[node_id] = vtk_index
		vtk_index += 1

	# Pre-scan elements for node counts (for auto order detection)
	tet_node_counts = {tet_id: len(cubit.get_expanded_connectivity("tet", tet_id)) for tet_id in tet_list}
	hex_node_counts = {hex_id: len(cubit.get_expanded_connectivity("hex", hex_id)) for hex_id in hex_list}
	wedge_node_counts = {wedge_id: len(cubit.get_expanded_connectivity("wedge", wedge_id)) for wedge_id in wedge_list}
	pyramid_node_counts = {pyramid_id: len(cubit.get_expanded_connectivity("pyramid", pyramid_id)) for pyramid_id in pyramid_list}
	tri_node_counts = {tri_id: len(cubit.get_expanded_connectivity("tri", tri_id)) for tri_id in tri_list}
	quad_node_counts = {quad_id: len(cubit.get_expanded_connectivity("quad", quad_id)) for quad_id in quad_list}
	edge_node_counts = {edge_id: len(cubit.get_expanded_connectivity("edge", edge_id)) for edge_id in edge_list}

	num_points = len(sorted_nodes)
	num_cells = len(tet_list) + len(hex_list) + len(wedge_list) + len(pyramid_list) + len(tri_list) + len(quad_list) + len(edge_list) + len(nodes_list)

	# Build connectivity and offsets arrays
	connectivity = []
	offsets = []
	cell_types = []
	block_ids = []
	current_offset = 0

	# Helper to get VTK type based on element type and node count
	def get_vtk_type(elem_type, node_count):
		if elem_type == "tet":
			return VTK_QUADRATIC_TETRA if node_count == 10 else VTK_TETRA
		elif elem_type == "hex":
			return VTK_QUADRATIC_HEXAHEDRON if node_count == 20 else VTK_HEXAHEDRON
		elif elem_type == "wedge":
			return VTK_QUADRATIC_WEDGE if node_count == 15 else VTK_WEDGE
		elif elem_type == "pyramid":
			return VTK_QUADRATIC_PYRAMID if node_count == 13 else VTK_PYRAMID
		elif elem_type == "tri":
			return VTK_QUADRATIC_TRIANGLE if node_count == 6 else VTK_TRIANGLE
		elif elem_type == "quad":
			return VTK_QUADRATIC_QUAD if node_count == 8 else VTK_QUAD
		elif elem_type == "edge":
			return VTK_QUADRATIC_EDGE if node_count == 3 else VTK_LINE
		else:
			return VTK_VERTEX

	# Process tetrahedra
	for tet_id in tet_list:
		nodes = cubit.get_expanded_connectivity("tet", tet_id)
		vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("tet", len(nodes)))
		# Find block ID for this element
		for block_id in cubit.get_block_id_list():
			if tet_id in _get_block_elements(cubit, block_id, "tet"):
				block_ids.append(block_id)
				break

	# Process hexahedra
	for hex_id in hex_list:
		nodes = cubit.get_expanded_connectivity("hex", hex_id)
		# HEX20 node reordering: Cubit -> VTK
		if len(nodes) == 20:
			vtk_nodes = [node_id_to_vtk_index[nodes[i]] for i in [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]]
		else:
			vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("hex", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if hex_id in _get_block_elements(cubit, block_id, "hex"):
				block_ids.append(block_id)
				break

	# Process wedges
	for wedge_id in wedge_list:
		nodes = cubit.get_expanded_connectivity("wedge", wedge_id)
		# WEDGE15 node reordering: Cubit -> VTK
		if len(nodes) == 15:
			vtk_nodes = [node_id_to_vtk_index[nodes[i]] for i in [0,1,2,3,4,5,6,7,8,12,13,14,9,10,11]]
		else:
			vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("wedge", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if wedge_id in _get_block_elements(cubit, block_id, "wedge"):
				block_ids.append(block_id)
				break

	# Process pyramids
	for pyramid_id in pyramid_list:
		nodes = cubit.get_expanded_connectivity("pyramid", pyramid_id)
		vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("pyramid", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if pyramid_id in _get_block_elements(cubit, block_id, "pyramid"):
				block_ids.append(block_id)
				break

	# Process triangles
	for tri_id in tri_list:
		nodes = cubit.get_expanded_connectivity("tri", tri_id)
		vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("tri", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if tri_id in _get_block_elements(cubit, block_id, "tri"):
				block_ids.append(block_id)
				break

	# Process quads
	for quad_id in quad_list:
		nodes = cubit.get_expanded_connectivity("quad", quad_id)
		vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("quad", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if quad_id in _get_block_elements(cubit, block_id, "face"):
				block_ids.append(block_id)
				break

	# Process edges
	for edge_id in edge_list:
		nodes = cubit.get_expanded_connectivity("edge", edge_id)
		vtk_nodes = [node_id_to_vtk_index[n] for n in nodes]
		connectivity.extend(vtk_nodes)
		current_offset += len(nodes)
		offsets.append(current_offset)
		cell_types.append(get_vtk_type("edge", len(nodes)))
		for block_id in cubit.get_block_id_list():
			if edge_id in _get_block_elements(cubit, block_id, "edge"):
				block_ids.append(block_id)
				break

	# Process point elements
	for node_id in nodes_list:
		connectivity.append(node_id_to_vtk_index[node_id])
		current_offset += 1
		offsets.append(current_offset)
		cell_types.append(VTK_VERTEX)
		for block_id in cubit.get_block_id_list():
			if node_id in _get_block_elements(cubit, block_id, "node"):
				block_ids.append(block_id)
				break

	# Write VTU file
	with open(FileName, 'w', encoding='UTF-8') as fid:
		fid.write('<?xml version="1.0"?>\n')
		fid.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">\n')
		fid.write('  <UnstructuredGrid>\n')
		fid.write(f'    <Piece NumberOfPoints="{num_points}" NumberOfCells="{num_cells}">\n')

		# Points
		fid.write('      <Points>\n')
		fid.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
		for node_id in sorted_nodes:
			coord = cubit.get_nodal_coordinates(node_id)
			fid.write(f'          {coord[0]} {coord[1]} {coord[2]}\n')
		fid.write('        </DataArray>\n')
		fid.write('      </Points>\n')

		# Cells
		fid.write('      <Cells>\n')

		# Connectivity
		fid.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
		fid.write('          ')
		for i, c in enumerate(connectivity):
			fid.write(f'{c} ')
			if (i + 1) % 20 == 0:
				fid.write('\n          ')
		fid.write('\n')
		fid.write('        </DataArray>\n')

		# Offsets
		fid.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
		fid.write('          ')
		for i, o in enumerate(offsets):
			fid.write(f'{o} ')
			if (i + 1) % 20 == 0:
				fid.write('\n          ')
		fid.write('\n')
		fid.write('        </DataArray>\n')

		# Cell types
		fid.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
		fid.write('          ')
		for i, t in enumerate(cell_types):
			fid.write(f'{t} ')
			if (i + 1) % 20 == 0:
				fid.write('\n          ')
		fid.write('\n')
		fid.write('        </DataArray>\n')

		fid.write('      </Cells>\n')

		# Cell Data (Block IDs)
		fid.write('      <CellData Scalars="BlockID">\n')
		fid.write('        <DataArray type="Int32" Name="BlockID" format="ascii">\n')
		fid.write('          ')
		for i, b in enumerate(block_ids):
			fid.write(f'{b} ')
			if (i + 1) % 20 == 0:
				fid.write('\n          ')
		fid.write('\n')
		fid.write('        </DataArray>\n')
		fid.write('      </CellData>\n')

		# Point Data (Node IDs)
		fid.write('      <PointData Scalars="NodeID">\n')
		fid.write('        <DataArray type="Int32" Name="NodeID" format="ascii">\n')
		fid.write('          ')
		for i, node_id in enumerate(sorted_nodes):
			fid.write(f'{node_id} ')
			if (i + 1) % 20 == 0:
				fid.write('\n          ')
		fid.write('\n')
		fid.write('        </DataArray>\n')
		fid.write('      </PointData>\n')

		fid.write('    </Piece>\n')
		fid.write('  </UnstructuredGrid>\n')
		fid.write('</VTKFile>\n')

	return cubit


########################################################################
###	Exodus II format (Cubit's native format)
########################################################################

def export_exodus(
	cubit: Any,
	filename: str,
	overwrite: bool = True,
	large_model: bool = False
) -> Any:
	"""Export mesh to Exodus II format.

	Exports the current Cubit mesh to Exodus II format (.exo, .e, .g).
	This is Cubit's native format and supports all element types and orders.

	Args:
		cubit: Cubit Python interface object
		filename: Output file path (typically .exo, .e, or .g extension)
		overwrite: Whether to overwrite existing file (default: True)
		large_model: Use 64-bit integers for large models (default: False)

	Returns:
		cubit: The cubit object (for method chaining)

	Supported elements:
		- 0D: NODE
		- 1D: BAR, BAR2, BAR3
		- 2D: TRI, TRI6, TRI7, QUAD, QUAD8, QUAD9
		- 3D: TET, TET10, TET11, HEX, HEX20, HEX27, WEDGE, WEDGE15, PYRAMID, PYRAMID13

	Example:
		>>> cubit.cmd("create brick x 1 y 1 z 1")
		>>> cubit.cmd("volume 1 scheme tetmesh")
		>>> cubit.cmd("mesh volume 1")
		>>> cubit.cmd("block 1 add tet all")
		>>> cubit.cmd("block 1 name 'solid'")
		>>> export_exodus(cubit, "cube.exo")

	Note:
		Exodus format natively supports nodesets and sidesets.
		For element order control (e.g., converting to TET10),
		blocks must contain mesh elements rather than geometry.
	"""
	_warn_mixed_element_types_in_blocks(cubit)


	# Build export command
	cmd_parts = ['export mesh']
	cmd_parts.append(f'"{filename}"')

	# Add options
	if overwrite:
		cmd_parts.append('overwrite')

	if large_model:
		cmd_parts.append('large')

	# Execute export
	cmd = ' '.join(cmd_parts)
	cubit.cmd(cmd)

	# Print summary
	_print_exodus_summary(cubit, filename)

	return cubit


def _print_exodus_summary(cubit: Any, filename: str) -> None:
	"""Print summary of exported Exodus mesh."""
	print(f"\nExodus export: {filename}")
	print("-" * 50)

	# Count elements by type
	element_counts = {}

	# 3D elements
	try:
		tet_count = cubit.get_tet_count()
		if tet_count > 0:
			element_counts['Tetrahedra'] = tet_count
	except:
		pass

	try:
		hex_count = cubit.get_hex_count()
		if hex_count > 0:
			element_counts['Hexahedra'] = hex_count
	except:
		pass

	try:
		wedge_count = cubit.get_wedge_count()
		if wedge_count > 0:
			element_counts['Wedges'] = wedge_count
	except:
		pass

	try:
		pyramid_count = cubit.get_pyramid_count()
		if pyramid_count > 0:
			element_counts['Pyramids'] = pyramid_count
	except:
		pass

	# 2D elements
	try:
		tri_count = cubit.get_tri_count()
		if tri_count > 0:
			element_counts['Triangles'] = tri_count
	except:
		pass

	try:
		quad_count = cubit.get_quad_count()
		if quad_count > 0:
			element_counts['Quads'] = quad_count
	except:
		pass

	# Print counts
	print(f"Nodes: {cubit.get_node_count()}")
	for elem_type, count in element_counts.items():
		print(f"{elem_type}: {count}")

	# Blocks
	block_ids = cubit.get_block_id_list()
	print(f"\nBlocks: {len(block_ids)}")
	for block_id in block_ids:
		name = cubit.get_exodus_entity_name("block", block_id)
		if name:
			print(f"  Block {block_id}: \"{name}\"")
		else:
			print(f"  Block {block_id}")

	# Nodesets
	nodeset_ids = cubit.get_nodeset_id_list()
	if len(nodeset_ids) > 0:
		print(f"\nNodesets: {len(nodeset_ids)}")
		for ns_id in nodeset_ids:
			name = cubit.get_exodus_entity_name("nodeset", ns_id)
			if name:
				print(f"  Nodeset {ns_id}: \"{name}\"")
			else:
				print(f"  Nodeset {ns_id}")

	# Sidesets
	sideset_ids = cubit.get_sideset_id_list()
	if len(sideset_ids) > 0:
		print(f"\nSidesets: {len(sideset_ids)}")
		for ss_id in sideset_ids:
			name = cubit.get_exodus_entity_name("sideset", ss_id)
			if name:
				print(f"  Sideset {ss_id}: \"{name}\"")
			else:
				print(f"  Sideset {ss_id}")

	print("-" * 50)


########################################################################
###	Function Aliases (snake_case naming convention)
###	For backward compatibility, original names are preserved.
###	New code should prefer snake_case names.
########################################################################

# Gmsh format alias
export_gmesh = export_Gmesh

# Nastran format alias
export_nastran = export_Nastran

# Netgen format alias
export_netgen = export_NetgenMesh


# ============================================================
# SetGeomInfo UV Utilities (for mesh.Curve() support)
# ============================================================

def compute_cylinder_uv(x: float, y: float, z: float, radius: float, height: float,
                        center: tuple = (0, 0, 0), axis: str = 'z') -> tuple:
    """
    Compute UV parameters for a cylindrical surface.

    Parameters:
    - x, y, z: Point coordinates
    - radius: Cylinder radius
    - height: Cylinder height
    - center: Center of cylinder (default origin)
    - axis: Cylinder axis ('x', 'y', or 'z')

    Returns:
    - (u, v): UV parameters where u is angle [0,1] and v is height [0,1]
    """
    cx, cy, cz = center

    if axis == 'z':
        dx, dy, dz = x - cx, y - cy, z - cz
        theta = math.atan2(dy, dx)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)
        v = (dz + height/2) / height
    elif axis == 'y':
        dx, dy, dz = x - cx, y - cy, z - cz
        theta = math.atan2(dz, dx)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)
        v = (dy + height/2) / height
    else:  # x-axis
        dx, dy, dz = x - cx, y - cy, z - cz
        theta = math.atan2(dz, dy)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)
        v = (dx + height/2) / height

    return u, max(0, min(1, v))


def compute_sphere_uv(x: float, y: float, z: float, radius: float,
                      center: tuple = (0, 0, 0)) -> tuple:
    """
    Compute UV parameters for a spherical surface.

    Parameters:
    - x, y, z: Point coordinates
    - radius: Sphere radius
    - center: Center of sphere (default origin)

    Returns:
    - (u, v): UV parameters (azimuth [0,1], polar [0,1])
    """
    cx, cy, cz = center
    dx, dy, dz = x - cx, y - cy, z - cz

    # Azimuthal angle (u)
    theta = math.atan2(dy, dx)
    if theta < 0:
        theta += 2 * math.pi
    u = theta / (2 * math.pi)

    # Polar angle (v)
    r = math.sqrt(dx*dx + dy*dy + dz*dz)
    if r > 1e-10:
        phi = math.acos(max(-1, min(1, dz / r)))
        v = phi / math.pi
    else:
        v = 0.5

    return u, v


def compute_torus_uv(x: float, y: float, z: float, major_radius: float, minor_radius: float,
                     center: tuple = (0, 0, 0), axis: str = 'z') -> tuple:
    """
    Compute UV parameters for a toroidal surface.

    Parameters:
    - x, y, z: Point coordinates
    - major_radius: Distance from center to tube center
    - minor_radius: Tube radius
    - center: Center of torus (default origin)
    - axis: Torus axis ('x', 'y', or 'z')

    Returns:
    - (u, v): UV parameters (major angle [0,1], minor angle [0,1])
    """
    cx, cy, cz = center
    dx, dy, dz = x - cx, y - cy, z - cz

    if axis == 'z':
        # Major angle (around the torus)
        theta = math.atan2(dy, dx)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)

        # Minor angle (around the tube)
        r_major = math.sqrt(dx*dx + dy*dy)
        phi = math.atan2(dz, r_major - major_radius)
        if phi < 0:
            phi += 2 * math.pi
        v = phi / (2 * math.pi)
    elif axis == 'y':
        theta = math.atan2(dz, dx)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)

        r_major = math.sqrt(dx*dx + dz*dz)
        phi = math.atan2(dy, r_major - major_radius)
        if phi < 0:
            phi += 2 * math.pi
        v = phi / (2 * math.pi)
    else:  # x-axis
        theta = math.atan2(dz, dy)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)

        r_major = math.sqrt(dy*dy + dz*dz)
        phi = math.atan2(dx, r_major - major_radius)
        if phi < 0:
            phi += 2 * math.pi
        v = phi / (2 * math.pi)

    return u, v


def set_cylinder_geominfo(ngmesh: Any, radius: float, height: float,
                          center: tuple = (0, 0, 0), axis: str = 'z',
                          tol: float = 0.01) -> int:
    """
    Set geominfo (UV parameters) for cylindrical surface elements.

    This enables mesh.Curve() to work correctly on external meshes.

    Parameters:
    - ngmesh: Netgen mesh
    - radius: Cylinder radius
    - height: Cylinder height
    - center: Center of cylinder
    - axis: Cylinder axis ('x', 'y', or 'z')
    - tol: Tolerance for detecting cylinder surface

    Returns:
    - Number of elements modified
    """
    points = [(p.p[0], p.p[1], p.p[2]) for p in ngmesh.Points()]
    cx, cy, cz = center
    modified = 0

    for el in ngmesh.Elements2D():
        verts = list(el.vertices)

        # Check if all vertices are on cylindrical surface
        on_cylinder = True
        for v in verts:
            x, y, z = points[v.nr - 1]
            dx, dy, dz = x - cx, y - cy, z - cz

            if axis == 'z':
                r = math.sqrt(dx*dx + dy*dy)
            elif axis == 'y':
                r = math.sqrt(dx*dx + dz*dz)
            else:
                r = math.sqrt(dy*dy + dz*dz)

            if abs(r - radius) > tol:
                on_cylinder = False
                break

        if on_cylinder:
            for i, v in enumerate(verts):
                x, y, z = points[v.nr - 1]
                u, uv_v = compute_cylinder_uv(x, y, z, radius, height, center, axis)
                el.SetGeomInfo(i, u, uv_v)
                modified += 1

    return modified


def set_sphere_geominfo(ngmesh: Any, radius: float, center: tuple = (0, 0, 0),
                        tol: float = 0.01) -> int:
    """
    Set geominfo (UV parameters) for spherical surface elements.

    Parameters:
    - ngmesh: Netgen mesh
    - radius: Sphere radius
    - center: Center of sphere
    - tol: Tolerance for detecting sphere surface

    Returns:
    - Number of elements modified
    """
    points = [(p.p[0], p.p[1], p.p[2]) for p in ngmesh.Points()]
    cx, cy, cz = center
    modified = 0

    for el in ngmesh.Elements2D():
        verts = list(el.vertices)

        # Check if all vertices are on sphere surface
        on_sphere = True
        for v in verts:
            x, y, z = points[v.nr - 1]
            dx, dy, dz = x - cx, y - cy, z - cz
            r = math.sqrt(dx*dx + dy*dy + dz*dz)

            if abs(r - radius) > tol:
                on_sphere = False
                break

        if on_sphere:
            for i, v in enumerate(verts):
                x, y, z = points[v.nr - 1]
                u, uv_v = compute_sphere_uv(x, y, z, radius, center)
                el.SetGeomInfo(i, u, uv_v)
                modified += 1

    return modified


def set_torus_geominfo(ngmesh: Any, major_radius: float, minor_radius: float,
                       center: tuple = (0, 0, 0), axis: str = 'z',
                       tol: float = 0.01) -> int:
    """
    Set geominfo (UV parameters) for toroidal surface elements.

    Parameters:
    - ngmesh: Netgen mesh
    - major_radius: Distance from center to tube center
    - minor_radius: Tube radius
    - center: Center of torus
    - axis: Torus axis ('x', 'y', or 'z')
    - tol: Tolerance for detecting torus surface

    Returns:
    - Number of elements modified
    """
    points = [(p.p[0], p.p[1], p.p[2]) for p in ngmesh.Points()]
    cx, cy, cz = center
    modified = 0

    for el in ngmesh.Elements2D():
        verts = list(el.vertices)

        # Check if all vertices are on torus surface
        on_torus = True
        for v in verts:
            x, y, z = points[v.nr - 1]
            dx, dy, dz = x - cx, y - cy, z - cz

            if axis == 'z':
                r_major = math.sqrt(dx*dx + dy*dy)
                dist = math.sqrt((r_major - major_radius)**2 + dz*dz)
            elif axis == 'y':
                r_major = math.sqrt(dx*dx + dz*dz)
                dist = math.sqrt((r_major - major_radius)**2 + dy*dy)
            else:
                r_major = math.sqrt(dy*dy + dz*dz)
                dist = math.sqrt((r_major - major_radius)**2 + dx*dx)

            if abs(dist - minor_radius) > tol:
                on_torus = False
                break

        if on_torus:
            for i, v in enumerate(verts):
                x, y, z = points[v.nr - 1]
                u, uv_v = compute_torus_uv(x, y, z, major_radius, minor_radius, center, axis)
                el.SetGeomInfo(i, u, uv_v)
                modified += 1

    return modified


def compute_cone_uv(x: float, y: float, z: float, base_radius: float, height: float,
                    center: tuple = (0, 0, 0), axis: str = 'z') -> tuple:
    """
    Compute UV parameters for a point on a cone surface.

    The cone has its apex at center + (0, 0, height/2) and base at center + (0, 0, -height/2)
    for axis='z'.

    Parameters:
    - x, y, z: Point coordinates
    - base_radius: Radius at the base of the cone
    - height: Height of the cone
    - center: Center of the cone (midpoint of axis)
    - axis: Cone axis ('x', 'y', or 'z')

    Returns:
    - (u, v) UV parameters where u is angle [0, 2*pi] and v is height [0, 1]
    """
    cx, cy, cz = center
    dx, dy, dz = x - cx, y - cy, z - cz

    if axis == 'z':
        u = math.atan2(dy, dx)
        v = (dz + height/2) / height  # 0 at base, 1 at apex
    elif axis == 'y':
        u = math.atan2(dx, dz)
        v = (dy + height/2) / height
    else:  # axis == 'x'
        u = math.atan2(dz, dy)
        v = (dx + height/2) / height

    # Normalize u to [0, 2*pi]
    if u < 0:
        u += 2 * math.pi

    return u, v


def set_cone_geominfo(ngmesh: Any, base_radius: float, height: float,
                      center: tuple = (0, 0, 0), axis: str = 'z',
                      tol: float = 0.01) -> int:
    """
    Set geominfo (UV parameters) for conical surface elements.

    Parameters:
    - ngmesh: Netgen mesh
    - base_radius: Radius at the base of the cone
    - height: Height of the cone
    - center: Center of cone (midpoint of axis)
    - axis: Cone axis ('x', 'y', or 'z')
    - tol: Tolerance for detecting cone surface

    Returns:
    - Number of vertex geominfo entries modified
    """
    points = [(p.p[0], p.p[1], p.p[2]) for p in ngmesh.Points()]
    cx, cy, cz = center
    modified = 0

    for el in ngmesh.Elements2D():
        verts = list(el.vertices)

        # Check if all vertices are on cone surface (not base)
        on_cone = True
        for v in verts:
            x, y, z = points[v.nr - 1]
            dx, dy, dz = x - cx, y - cy, z - cz

            if axis == 'z':
                local_h = dz + height / 2  # 0 at base, height at apex
                expected_r = base_radius * (1 - local_h / height)
                actual_r = math.sqrt(dx*dx + dy*dy)
            elif axis == 'y':
                local_h = dy + height / 2
                expected_r = base_radius * (1 - local_h / height)
                actual_r = math.sqrt(dx*dx + dz*dz)
            else:
                local_h = dx + height / 2
                expected_r = base_radius * (1 - local_h / height)
                actual_r = math.sqrt(dy*dy + dz*dz)

            # Check if point is on cone surface (not at apex or base)
            if local_h < 0 or local_h > height:
                on_cone = False
                break
            if expected_r < tol and actual_r < tol:
                # At apex, skip
                on_cone = False
                break
            if abs(actual_r - expected_r) > tol:
                on_cone = False
                break

        if on_cone:
            for i, v in enumerate(verts):
                x, y, z = points[v.nr - 1]
                u, uv_v = compute_cone_uv(x, y, z, base_radius, height, center, axis)
                el.SetGeomInfo(i, u, uv_v)
                modified += 1

    return modified


########################################################################
###	Name-based face mapping for OCC-Cubit interoperability
########################################################################

def name_occ_faces(shape: Any, prefix: str = "occ_face_") -> int:
    """
    Assign unique names to all faces of an OCC shape.

    This enables reliable face mapping between OCC and Cubit through STEP.
    The names are preserved when exporting/importing STEP files.

    Parameters:
    - shape: OCC shape (from netgen.occ)
    - prefix: Name prefix (default: "occ_face_")

    Returns:
    - Number of faces named

    Example:
        from netgen.occ import Box, Cylinder, gp_Pnt, gp_Ax2, gp_Dir

        # Create geometry
        brick = Box(gp_Pnt(-1,-1,-1), gp_Pnt(1,1,1))
        cyl = Cylinder(gp_Ax2(gp_Pnt(0,0,-2), gp_Dir(0,0,1)), 0.3, 4)
        shape = brick - cyl

        # Name faces before STEP export
        cubit_mesh_export.name_occ_faces(shape)
        shape.WriteStep("geometry.step")
    """
    count = 0
    for i, face in enumerate(shape.faces):
        face.name = f"{prefix}{i}"
        count += 1
    return count


def export_netgen_with_names(cubit: Any, geometry: Any) -> Any:
    """
    Export Cubit mesh to Netgen format using name-based face mapping.

    This function uses face names in STEP to establish exact correspondence
    between OCC faces and Cubit surfaces. Use this after:
    1. Creating geometry in OCC with name_occ_faces()
    2. Exporting STEP from OCC
    3. Importing STEP into Cubit
    4. Meshing in Cubit

    Parameters:
    - cubit: Cubit Python interface
    - geometry: OCCGeometry object (loaded from the same STEP file)

    Returns:
    - Netgen mesh with correct face indices

    Example:
        from netgen.occ import OCCGeometry
        import cubit_mesh_export

        # After OCC created and exported STEP with named faces
        geo = OCCGeometry("geometry.step")

        # After Cubit imported STEP and meshed
        ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)

        # Apply SetGeomInfo for curved surfaces
        cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=0.3, height=2.0)

        # Now mesh.Curve() works correctly
        mesh = Mesh(ngmesh)
        mesh.Curve(2)
    """
    _warn_mixed_element_types_in_blocks(cubit)

    import re
    from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Element3D, Element2D, FaceDescriptor
    from netgen.csg import Pnt

    # Node ordering conversion tables: Cubit -> Netgen
    TET_ORDERING = [0, 1, 2, 3]
    HEX_ORDERING = [0, 1, 5, 4, 3, 2, 6, 7]
    WEDGE_ORDERING = [0, 2, 1, 3, 5, 4]
    PYRAMID_ORDERING = [3, 2, 1, 0, 4]
    TRI_ORDERING = [0, 1, 2]
    QUAD_ORDERING = [0, 1, 2, 3]

    # Build Cubit surface ID -> OCC face index mapping using names
    surface_to_occ = {}
    surface_ids = cubit.get_entities('surface')
    for sid in surface_ids:
        name = cubit.get_entity_name('surface', sid)
        match = re.match(r'occ_face_(\d+)', name)
        if match:
            surface_to_occ[sid] = int(match.group(1))

    if not surface_to_occ:
        raise ValueError(
            "No OCC face names found in Cubit surfaces. "
            "Make sure to use name_occ_faces() before exporting STEP from OCC."
        )

    # Create Netgen mesh
    ngmesh = NetgenMesh()
    ngmesh.dim = 3

    # Collect all nodes from all blocks
    node_map = {}
    all_nodes = set()

    for block_id in cubit.get_block_id_list():
        elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face"]
        for elem_type in elem_types:
            for element_id in _get_block_elements(cubit, block_id, elem_type):
                node_ids = cubit.get_connectivity(elem_type, element_id)
                all_nodes.update(node_ids)

    # Add nodes to mesh
    for nid in sorted(all_nodes):
        coords = cubit.get_nodal_coordinates(nid)
        pi = ngmesh.Add(MeshPoint(Pnt(*coords)))
        node_map[nid] = pi

    # Add 3D volume elements (iterate all blocks)
    material_index = 0

    for block_id in cubit.get_block_id_list():
        block_name = cubit.get_exodus_entity_name("block", block_id)

        tet_list = _get_block_elements(cubit, block_id, "tet")
        hex_list = _get_block_elements(cubit, block_id, "hex")
        wedge_list = _get_block_elements(cubit, block_id, "wedge")
        pyramid_list = _get_block_elements(cubit, block_id, "pyramid")

        if len(tet_list) + len(hex_list) + len(wedge_list) + len(pyramid_list) > 0:
            material_index += 1
            ngmesh.SetMaterial(material_index, block_name)

            for tet_id in tet_list:
                nodes = cubit.get_connectivity("tet", tet_id)
                ng_nodes = [node_map[nodes[i]] for i in TET_ORDERING]
                ngmesh.Add(Element3D(material_index, ng_nodes))

            for hex_id in hex_list:
                nodes = cubit.get_connectivity("hex", hex_id)
                ng_nodes = [node_map[nodes[i]] for i in HEX_ORDERING]
                ngmesh.Add(Element3D(material_index, ng_nodes))

            for wedge_id in wedge_list:
                nodes = cubit.get_connectivity("wedge", wedge_id)
                ng_nodes = [node_map[nodes[i]] for i in WEDGE_ORDERING]
                ngmesh.Add(Element3D(material_index, ng_nodes))

            for pyramid_id in pyramid_list:
                nodes = cubit.get_connectivity("pyramid", pyramid_id)
                ng_nodes = [node_map[nodes[i]] for i in PYRAMID_ORDERING]
                ngmesh.Add(Element3D(material_index, ng_nodes))

    # Create FaceDescriptors for each OCC face
    num_occ_faces = len(geometry.shape.faces)
    fd_map = {}
    for occ_idx in range(num_occ_faces):
        fd = FaceDescriptor(bc=occ_idx+1, surfnr=occ_idx+1)
        fd.bcname = f'face_{occ_idx}'
        fd_idx = ngmesh.Add(fd)
        ngmesh.SetBCName(fd_idx - 1, f'face_{occ_idx}')
        fd_map[occ_idx] = fd_idx

    # Build element -> surface mapping
    tri_to_surface = {}
    quad_to_surface = {}
    for sid in surface_ids:
        for tri_id in cubit.get_surface_tris(sid):
            tri_to_surface[tri_id] = sid
        for quad_id in cubit.get_surface_quads(sid):
            quad_to_surface[quad_id] = sid

    # Add 2D boundary elements with correct face indices (iterate all blocks)
    for block_id in cubit.get_block_id_list():
        tri_list = _get_block_elements(cubit, block_id, "tri")
        quad_list = _get_block_elements(cubit, block_id, "face")

        for tri_id in tri_list:
            nodes = cubit.get_connectivity('tri', tri_id)
            if all(n in node_map for n in nodes):
                ng_nodes = [node_map[nodes[i]] for i in TRI_ORDERING]
                sid = tri_to_surface.get(tri_id)
                if sid and sid in surface_to_occ:
                    occ_idx = surface_to_occ[sid]
                    ngmesh.Add(Element2D(fd_map[occ_idx], ng_nodes))

        for quad_id in quad_list:
            nodes = cubit.get_connectivity('face', quad_id)
            if all(n in node_map for n in nodes):
                ng_nodes = [node_map[nodes[i]] for i in QUAD_ORDERING]
                sid = quad_to_surface.get(quad_id)
                if sid and sid in surface_to_occ:
                    occ_idx = surface_to_occ[sid]
                    ngmesh.Add(Element2D(fd_map[occ_idx], ng_nodes))

    # Set geometry reference
    ngmesh.SetGeometry(geometry)

    return ngmesh


########################################################################
###	Auto-register Cubit toolbar panels on first import
########################################################################

def _auto_register_panels():
    """Register Cubit toolbar panels if not already registered.

    This runs once when cubit_mesh_export is first imported inside Cubit GUI.
    It checks ~/.cubit for the marker block and adds it if missing.
    Only takes effect when Cubit is restarted (modifies startup file).
    Skipped in batch/subprocess mode (no GUI).
    """
    import os

    # Only run in Cubit GUI environment (PySide6 or PyQt5 must be available)
    try:
        try:
            import PySide6
        except ImportError:
            import PyQt5
    except ImportError:
        return  # No Qt available -- not in Cubit GUI

    marker = "## BEGIN cubit_mesh_export toolbar"
    cubit_file = os.path.join(os.path.expanduser("~"), ".cubit")

    # Check if already registered
    if os.path.isfile(cubit_file):
        try:
            with open(cubit_file, "r", encoding="utf-8", errors="ignore") as f:
                if marker in f.read():
                    return  # Already registered
        except Exception:
            return

    # Register panels silently
    try:
        from radia.install_panels import install_panels
        install_panels()
    except Exception:
        try:
            # Fallback: direct import when running from src/radia/
            from install_panels import install_panels
            install_panels()
        except Exception:
            pass  # Not critical -- user can run manually

_auto_register_panels()
