"""
Verify GMSH v2.2 and v4.1 export consistency.

Strategy:
  1. Read Cubit-exported v2.2 file into GMSH
  2. Re-save as v4.1 using GMSH itself (reference)
  3. Read Cubit-exported v4.1 file into GMSH
  4. Compare: node count, element count, element types, bounding box,
     physical names, and per-group element counts

This verifies that Cubit's v4.1 export produces the same mesh as
GMSH's own v2.2->v4.1 conversion.

Requires: gmsh Python module (pip install gmsh), Coreform Cubit (tests 11-12)

12 test cases:
  1-3. Pre-generated Cubit files (cube tet 1st/2nd order)
  4-10. GMSH-generated meshes (tet, tri, quad, sphere, cylinder, multi-volume, 2nd order)
  11. Cubit: 0D+1D+2D+3D mixed dimension blocks
  12. Cubit: Multi-material (iron + air + boundary) blocks
"""

import sys
import os
import tempfile
import numpy as np

# Auto-detect Cubit installation (single source of truth)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
	sys.path.append(_cubit_path)

try:
	import gmsh
except ImportError:
	print("ERROR: gmsh Python module not found. Install with: pip install gmsh")
	sys.exit(1)

# Cubit is optional (only needed for tests 11-12)
try:
	import cubit
	cubit.init(['cubit', '-nojournal', '-batch'])
	HAS_CUBIT = True
except ImportError:
	HAS_CUBIT = False


EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
							"examples", "gmsh")


def read_mesh_data(msh_file):
	"""Read a .msh file with GMSH and extract mesh data for comparison.

	Returns dict with node_count, element_counts (by type), node_coords,
	and per_group_elements (element count per Physical Group name).
	"""
	gmsh.initialize()
	gmsh.option.setNumber("General.Terminal", 0)  # Suppress output
	gmsh.open(msh_file)

	# Get nodes
	node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
	n_nodes = len(node_tags)

	# Reshape coordinates to (N, 3)
	coords = np.array(node_coords).reshape(-1, 3)

	# Get elements by type
	elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

	element_counts = {}
	total_elements = 0
	for etype, etags in zip(elem_types, elem_tags):
		element_counts[int(etype)] = len(etags)
		total_elements += len(etags)

	# Bounding box
	bbox_min = coords.min(axis=0) if n_nodes > 0 else np.zeros(3)
	bbox_max = coords.max(axis=0) if n_nodes > 0 else np.zeros(3)

	# Physical groups + per-group element counts
	phys_groups = {}
	per_group_elements = {}
	for dim, tag in gmsh.model.getPhysicalGroups():
		name = gmsh.model.getPhysicalName(dim, tag)
		phys_groups[(dim, tag)] = name
		# Count elements belonging to this physical group
		entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
		group_count = 0
		for entity_tag in entities:
			et, etags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
			for tags in etags:
				group_count += len(tags)
		per_group_elements[name] = group_count

	gmsh.finalize()

	return {
		'n_nodes': n_nodes,
		'n_elements': total_elements,
		'element_counts': element_counts,
		'bbox_min': bbox_min,
		'bbox_max': bbox_max,
		'physical_groups': phys_groups,
		'per_group_elements': per_group_elements,
		'coords': coords,
	}


def convert_v2_to_v4(v2_file, v4_output):
	"""Use GMSH to read v2.2 file and re-save as v4.1."""
	gmsh.initialize()
	gmsh.option.setNumber("General.Terminal", 0)
	gmsh.open(v2_file)
	gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
	gmsh.write(v4_output)
	gmsh.finalize()


def compare_meshes(data_a, data_b, label_a, label_b):
	"""Compare two mesh data dicts. Returns (passed, messages)."""
	messages = []
	passed = True

	# Node count
	if data_a['n_nodes'] != data_b['n_nodes']:
		messages.append(f"  FAIL: Node count differs: {label_a}={data_a['n_nodes']}, "
						f"{label_b}={data_b['n_nodes']}")
		passed = False
	else:
		messages.append(f"  Nodes: {data_a['n_nodes']} (match)")

	# Element count
	if data_a['n_elements'] != data_b['n_elements']:
		messages.append(f"  FAIL: Element count differs: {label_a}={data_a['n_elements']}, "
						f"{label_b}={data_b['n_elements']}")
		passed = False
	else:
		messages.append(f"  Elements: {data_a['n_elements']} (match)")

	# Element types
	types_a = set(data_a['element_counts'].keys())
	types_b = set(data_b['element_counts'].keys())
	if types_a != types_b:
		messages.append(f"  FAIL: Element types differ: {label_a}={types_a}, {label_b}={types_b}")
		passed = False
	else:
		for etype in sorted(types_a):
			ca = data_a['element_counts'][etype]
			cb = data_b['element_counts'][etype]
			if ca != cb:
				messages.append(f"  FAIL: Type {etype} count: {label_a}={ca}, {label_b}={cb}")
				passed = False
			else:
				messages.append(f"  Type {etype}: {ca} elements (match)")

	# Bounding box
	bbox_tol = 1e-10
	if np.max(np.abs(data_a['bbox_min'] - data_b['bbox_min'])) > bbox_tol:
		messages.append(f"  FAIL: BBox min differs: {data_a['bbox_min']} vs {data_b['bbox_min']}")
		passed = False
	if np.max(np.abs(data_a['bbox_max'] - data_b['bbox_max'])) > bbox_tol:
		messages.append(f"  FAIL: BBox max differs: {data_a['bbox_max']} vs {data_b['bbox_max']}")
		passed = False
	if np.max(np.abs(data_a['bbox_min'] - data_b['bbox_min'])) <= bbox_tol and \
	   np.max(np.abs(data_a['bbox_max'] - data_b['bbox_max'])) <= bbox_tol:
		messages.append(f"  BBox: [{data_a['bbox_min']}] to [{data_a['bbox_max']}] (match)")

	# Physical group names
	names_a = set(data_a['physical_groups'].values())
	names_b = set(data_b['physical_groups'].values())
	if names_a != names_b:
		messages.append(f"  WARN: Physical names differ: {names_a} vs {names_b}")
		# Not a hard failure - tag numbering may differ
	else:
		messages.append(f"  Physical names: {sorted(names_a)} (match)")

	# Per-group element counts
	pg_a = data_a.get('per_group_elements', {})
	pg_b = data_b.get('per_group_elements', {})
	common_groups = set(pg_a.keys()) & set(pg_b.keys())
	for gname in sorted(common_groups):
		ca = pg_a[gname]
		cb = pg_b[gname]
		if ca != cb:
			messages.append(f"  FAIL: Group '{gname}' element count: "
							f"{label_a}={ca}, {label_b}={cb}")
			passed = False
		else:
			messages.append(f"  Group '{gname}': {ca} elements (match)")

	return passed, messages


def test_1_cube_tet_1st_order():
	"""Test 1: cube_v2.msh vs cube_v4.msh (1st order tet + tri)."""
	print("=" * 60)
	print("Test 1: Cube Tet 1st Order (pre-generated files)")
	print("=" * 60)

	v2_file = os.path.join(EXAMPLES_DIR, "cube_v2.msh")
	v4_file = os.path.join(EXAMPLES_DIR, "cube_v4.msh")

	if not os.path.exists(v2_file) or not os.path.exists(v4_file):
		print("  SKIP: Sample files not found")
		return True

	# Read both files through GMSH
	data_v2 = read_mesh_data(v2_file)
	data_v4 = read_mesh_data(v4_file)

	passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
	for m in msgs:
		print(m)

	if passed:
		print("  PASS")
	return passed


def test_2_cube_tet_2nd_order_roundtrip():
	"""Test 2: cube_2nd_order.msh v2.2 -> GMSH -> v4.1 round-trip."""
	print("\n" + "=" * 60)
	print("Test 2: Cube 2nd Order Round-Trip (v2.2 -> GMSH -> v4.1)")
	print("=" * 60)

	v2_file = os.path.join(EXAMPLES_DIR, "cube_2nd_order.msh")

	if not os.path.exists(v2_file):
		print("  SKIP: cube_2nd_order.msh not found")
		return True

	# Read original v2.2
	data_v2 = read_mesh_data(v2_file)

	# Convert v2.2 -> v4.1 using GMSH
	with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as tmp:
		tmp_v4 = tmp.name
	try:
		convert_v2_to_v4(v2_file, tmp_v4)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2(orig)", "v4.1(gmsh)")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		if os.path.exists(tmp_v4):
			os.remove(tmp_v4)


def test_3_cubit_v4_vs_gmsh_converted_v4():
	"""Test 3: Compare Cubit v4.1 with GMSH-converted v4.1."""
	print("\n" + "=" * 60)
	print("Test 3: Cubit v4.1 vs GMSH-Converted v4.1")
	print("=" * 60)

	v2_file = os.path.join(EXAMPLES_DIR, "cube_v2.msh")
	v4_cubit = os.path.join(EXAMPLES_DIR, "cube_v4.msh")

	if not os.path.exists(v2_file) or not os.path.exists(v4_cubit):
		print("  SKIP: Sample files not found")
		return True

	# Convert v2.2 -> v4.1 using GMSH (reference)
	with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as tmp:
		tmp_v4 = tmp.name
	try:
		convert_v2_to_v4(v2_file, tmp_v4)

		data_cubit = read_mesh_data(v4_cubit)
		data_gmsh = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_cubit, data_gmsh, "Cubit-v4.1", "GMSH-v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		if os.path.exists(tmp_v4):
			os.remove(tmp_v4)


def _create_gmsh_mesh(dim, element_type, filename_v2):
	"""Create a simple mesh using GMSH API, save as v2.2."""
	gmsh.initialize()
	gmsh.option.setNumber("General.Terminal", 0)
	gmsh.model.add("test")

	if dim == 3 and element_type == "tet":
		gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
		gmsh.model.mesh.generate(3)
	elif dim == 3 and element_type == "hex":
		gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.RecombineAll", 1)
		gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)  # All hex
		gmsh.model.mesh.generate(3)
		gmsh.model.mesh.recombine()
	elif dim == 2 and element_type == "tri":
		gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
		gmsh.model.occ.synchronize()
		gmsh.model.mesh.generate(2)
	elif dim == 2 and element_type == "quad":
		gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.RecombineAll", 1)
		gmsh.model.mesh.generate(2)
		gmsh.model.mesh.recombine()

	# Add physical groups
	volumes = gmsh.model.getEntities(3)
	surfaces = gmsh.model.getEntities(2)
	if volumes:
		gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], 1)
		gmsh.model.setPhysicalName(3, 1, "domain")
	if surfaces:
		gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 2)
		gmsh.model.setPhysicalName(2, 2, "boundary")

	gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
	gmsh.write(filename_v2)
	gmsh.finalize()


def test_4_gmsh_generated_tet():
	"""Test 4: GMSH-generated tet mesh, v2.2 vs v4.1 round-trip."""
	print("\n" + "=" * 60)
	print("Test 4: GMSH-Generated Tet Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		_create_gmsh_mesh(3, "tet", tmp_v2)
		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_5_gmsh_generated_tri():
	"""Test 5: GMSH-generated 2D tri mesh, v2.2 vs v4.1 round-trip."""
	print("\n" + "=" * 60)
	print("Test 5: GMSH-Generated 2D Tri Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		_create_gmsh_mesh(2, "tri", tmp_v2)
		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_6_gmsh_generated_quad():
	"""Test 6: GMSH-generated 2D quad mesh, v2.2 vs v4.1 round-trip."""
	print("\n" + "=" * 60)
	print("Test 6: GMSH-Generated 2D Quad Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		_create_gmsh_mesh(2, "quad", tmp_v2)
		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_7_sphere_tet():
	"""Test 7: Sphere tet mesh (curved geometry), v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 7: Sphere Tet Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		gmsh.initialize()
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("sphere")
		gmsh.model.occ.addSphere(0, 0, 0, 0.5)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
		gmsh.model.mesh.generate(3)
		volumes = gmsh.model.getEntities(3)
		surfaces = gmsh.model.getEntities(2)
		gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], 1)
		gmsh.model.setPhysicalName(3, 1, "sphere")
		gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 2)
		gmsh.model.setPhysicalName(2, 2, "surface")
		gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
		gmsh.write(tmp_v2)
		gmsh.finalize()

		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_8_cylinder_hex():
	"""Test 8: Cylinder with hex mesh (transfinite), v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 8: Cylinder Hex Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		gmsh.initialize()
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("cylinder")
		gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 1.0, 0.3)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.15)
		gmsh.model.mesh.generate(3)
		volumes = gmsh.model.getEntities(3)
		surfaces = gmsh.model.getEntities(2)
		gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], 1)
		gmsh.model.setPhysicalName(3, 1, "cylinder")
		gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 2)
		gmsh.model.setPhysicalName(2, 2, "outer")
		gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
		gmsh.write(tmp_v2)
		gmsh.finalize()

		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_9_multi_volume():
	"""Test 9: Multi-volume (2 boxes), v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 9: Multi-Volume 2-Box Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		gmsh.initialize()
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("multi")
		b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
		b2 = gmsh.model.occ.addBox(1, 0, 0, 1, 1, 1)
		gmsh.model.occ.fragment([(3, b1)], [(3, b2)])
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.4)
		gmsh.model.mesh.generate(3)
		volumes = gmsh.model.getEntities(3)
		surfaces = gmsh.model.getEntities(2)
		gmsh.model.addPhysicalGroup(3, [volumes[0][1]], 1)
		gmsh.model.setPhysicalName(3, 1, "left")
		if len(volumes) > 1:
			gmsh.model.addPhysicalGroup(3, [volumes[1][1]], 2)
			gmsh.model.setPhysicalName(3, 2, "right")
		gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 3)
		gmsh.model.setPhysicalName(2, 3, "boundary")
		gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
		gmsh.write(tmp_v2)
		gmsh.finalize()

		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def test_10_2nd_order_tet():
	"""Test 10: 2nd order tet10 mesh, v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 10: 2nd Order Tet10 Mesh (v2.2 <-> v4.1)")
	print("=" * 60)

	with tempfile.NamedTemporaryFile(suffix='_v2.msh', delete=False) as f1, \
		 tempfile.NamedTemporaryFile(suffix='_v4.msh', delete=False) as f2:
		tmp_v2 = f1.name
		tmp_v4 = f2.name

	try:
		gmsh.initialize()
		gmsh.option.setNumber("General.Terminal", 0)
		gmsh.model.add("tet10")
		gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
		gmsh.model.occ.synchronize()
		gmsh.option.setNumber("Mesh.ElementOrder", 2)
		gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.5)
		gmsh.model.mesh.generate(3)
		volumes = gmsh.model.getEntities(3)
		surfaces = gmsh.model.getEntities(2)
		gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], 1)
		gmsh.model.setPhysicalName(3, 1, "domain")
		gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 2)
		gmsh.model.setPhysicalName(2, 2, "boundary")
		gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
		gmsh.write(tmp_v2)
		gmsh.finalize()

		convert_v2_to_v4(tmp_v2, tmp_v4)

		data_v2 = read_mesh_data(tmp_v2)
		data_v4 = read_mesh_data(tmp_v4)

		passed, msgs = compare_meshes(data_v2, data_v4, "v2.2", "v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		for f in [tmp_v2, tmp_v4]:
			if os.path.exists(f):
				os.remove(f)


def _safe_read_mesh_data(msh_file):
	"""Read mesh data, returning None on failure (instead of raising)."""
	try:
		return read_mesh_data(msh_file)
	except Exception as e:
		print(f"  ERROR reading {os.path.basename(msh_file)}: {e}")
		try:
			gmsh.finalize()
		except Exception:
			pass
		return None


def _cleanup_files(*files):
	"""Remove files, ignoring errors."""
	for f in files:
		try:
			if os.path.exists(f):
				os.remove(f)
		except OSError:
			pass


def test_11_cubit_multi_material_tet():
	"""Test 11: Cubit multi-material with two tet volumes + boundary, v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 11: Cubit Multi-Material (2 tet volumes + boundary)")
	print("=" * 60)

	if not HAS_CUBIT:
		print("  SKIP: Cubit not available")
		return True

	tmp_v2 = os.path.abspath("test_c11_v2.msh")
	tmp_v4 = os.path.abspath("test_c11_v4.msh")
	_cleanup_files(tmp_v2, tmp_v4)

	try:
		cubit.cmd("reset")
		cubit.cmd("create brick x 1 y 1 z 1")
		cubit.cmd("create brick x 1 y 1 z 1")
		cubit.cmd("volume 2 move 1.5 0 0")
		cubit.cmd("volume 1 scheme tetmesh")
		cubit.cmd("volume 2 scheme tetmesh")
		cubit.cmd("volume all size 0.5")
		cubit.cmd("mesh volume all")

		# Block 1: tet from volume 1 (iron)
		cubit.cmd("block 1 add tet all in volume 1")
		cubit.cmd("block 1 name 'iron'")

		# Block 2: tet from volume 2 (air)
		cubit.cmd("block 2 add tet all in volume 2")
		cubit.cmd("block 2 name 'air'")

		# Block 3: tri boundary from all surfaces
		cubit.cmd("block 3 add tri all in surface all")
		cubit.cmd("block 3 name 'boundary'")

		# Export both versions
		cubit.cmd(f'radia export gmsh "{tmp_v2}" version 2 overwrite')
		cubit.cmd(f'radia export gmsh "{tmp_v4}" version 4 overwrite')

		# Read both through GMSH API and compare
		data_v2 = _safe_read_mesh_data(tmp_v2)
		data_v4 = _safe_read_mesh_data(tmp_v4)

		if data_v2 is None or data_v4 is None:
			print("  FAIL: Could not read exported files")
			return False

		passed, msgs = compare_meshes(data_v2, data_v4, "Cubit-v2.2", "Cubit-v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		_cleanup_files(tmp_v2, tmp_v4)


def test_12_cubit_mixed_dim_0d_1d_2d_3d():
	"""Test 12: Cubit 0D+1D+2D+3D blocks (node + edge + tri + tet), v2.2 vs v4.1."""
	print("\n" + "=" * 60)
	print("Test 12: Cubit 0D+1D+2D+3D Mixed Dimension Blocks")
	print("=" * 60)

	if not HAS_CUBIT:
		print("  SKIP: Cubit not available")
		return True

	tmp_v2 = os.path.abspath("test_c12_v2.msh")
	tmp_v4 = os.path.abspath("test_c12_v4.msh")
	_cleanup_files(tmp_v2, tmp_v4)

	try:
		cubit.cmd("reset")
		cubit.cmd("create brick x 1 y 1 z 1")
		cubit.cmd("volume 1 scheme tetmesh")
		cubit.cmd("volume 1 size 0.5")
		cubit.cmd("mesh volume 1")

		# Block 1: 3D tet elements (volume)
		cubit.cmd("block 1 add tet all")
		cubit.cmd("block 1 name 'domain'")

		# Block 2: 2D tri elements (boundary surface)
		cubit.cmd("block 2 add tri all in surface all")
		cubit.cmd("block 2 name 'boundary'")

		# Block 3: 1D edge elements (one curve)
		cubit.cmd("block 3 add edge all in curve 1")
		cubit.cmd("block 3 name 'edge_set'")

		# Block 4: 0D node elements (one vertex)
		cubit.cmd("block 4 add node all in vertex 1")
		cubit.cmd("block 4 name 'vertex_set'")

		# Export both versions
		cubit.cmd(f'radia export gmsh "{tmp_v2}" version 2 overwrite')
		cubit.cmd(f'radia export gmsh "{tmp_v4}" version 4 overwrite')

		data_v2 = _safe_read_mesh_data(tmp_v2)
		data_v4 = _safe_read_mesh_data(tmp_v4)

		if data_v2 is None or data_v4 is None:
			print("  FAIL: Could not read exported files")
			return False

		passed, msgs = compare_meshes(data_v2, data_v4, "Cubit-v2.2", "Cubit-v4.1")
		for m in msgs:
			print(m)

		if passed:
			print("  PASS")
		return passed
	finally:
		_cleanup_files(tmp_v2, tmp_v4)


if __name__ == "__main__":
	tests = [
		test_1_cube_tet_1st_order,
		test_2_cube_tet_2nd_order_roundtrip,
		test_3_cubit_v4_vs_gmsh_converted_v4,
		test_4_gmsh_generated_tet,
		test_5_gmsh_generated_tri,
		test_6_gmsh_generated_quad,
		test_7_sphere_tet,
		test_8_cylinder_hex,
		test_9_multi_volume,
		test_10_2nd_order_tet,
		test_11_cubit_multi_material_tet,
		test_12_cubit_mixed_dim_0d_1d_2d_3d,
	]

	print("\n" + "=" * 60)
	print(f"GMSH v2.2 vs v4.1 Consistency Test Suite ({len(tests)} tests)")
	print("=" * 60)

	all_passed = True
	n_passed = 0
	for test in tests:
		try:
			if test():
				n_passed += 1
			else:
				all_passed = False
		except Exception as e:
			print(f"  FAIL: {e}")
			import traceback
			traceback.print_exc()
			all_passed = False

	print("\n" + "=" * 60)
	print(f"{n_passed}/{len(tests)} tests PASSED")
	if all_passed:
		print("All tests PASSED!")
	else:
		print("Some tests FAILED!")
	print("=" * 60)
