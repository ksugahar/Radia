"""
Lint rules for Cubit mesh export Python scripts.

Each rule function receives (filepath, lines) and returns a list of findings.
A finding is a dict with keys: line, severity, rule, message.
"""

import re
import os
from typing import List, Dict


def check_missing_block_registration(filepath: str, lines: List[str]) -> List[Dict]:
	"""CRITICAL: Export called without any block registration."""
	findings = []
	export_pattern = re.compile(
		r'(?:cubit_mesh_export\.|cme\.)?export_'
		r'(?:Gmsh_ver[24]|gmsh_v[24]|Nastran|nastran|'
		r'curved|exodus)\s*\('
		r'|export\s+(?:netgen|gmsh|nastran|vtk)\s'
	)
	block_pattern = re.compile(r'cubit\.cmd\s*\(\s*["\'].*\bblock\b', re.IGNORECASE)

	export_line = None
	has_block = False

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		if export_pattern.search(stripped) and export_line is None:
			export_line = i
		if block_pattern.search(stripped):
			has_block = True

	if export_line is not None and not has_block:
		findings.append({
			'line': export_line,
			'severity': 'CRITICAL',
			'rule': 'missing-block-registration',
			'message': (
				'Export function called but no block registration found. '
				'Add: cubit.cmd("block 1 add tet all") before export.'
			),
		})
	return findings


def check_geometry_block_for_2nd_order(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Geometry block used with element type conversion."""
	findings = []

	# Find blocks that add geometry (volume/surface)
	geom_block_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+(?:add\s+)?(?:volume|surface)\b',
		re.IGNORECASE
	)
	# Find blocks that set element type
	elem_type_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+element\s+type\b',
		re.IGNORECASE
	)

	geom_blocks = {}  # block_id -> line_number
	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		match = geom_block_pattern.search(stripped)
		if match:
			geom_blocks[match.group(1)] = i

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		match = elem_type_pattern.search(stripped)
		if match and match.group(1) in geom_blocks:
			findings.append({
				'line': i,
				'severity': 'HIGH',
				'rule': 'geometry-block-2nd-order',
				'message': (
					f'Block {match.group(1)} contains geometry (volume/surface), '
					f'not mesh elements. Element type conversion has no effect on '
					f'geometry blocks. Use "block {match.group(1)} add tet all" instead.'
				),
			})
	return findings


def check_missing_cubit_init(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Script imports cubit but does not call cubit.init()."""
	findings = []
	has_cubit_import = any(
		re.search(r'^\s*import\s+cubit\b', line) or
		re.search(r'^\s*from\s+cubit\b', line)
		for line in lines
	)
	if not has_cubit_import:
		return findings

	has_init = any('cubit.init' in line for line in lines)
	if not has_init:
		findings.append({
			'line': 1,
			'severity': 'HIGH',
			'rule': 'missing-cubit-init',
			'message': (
				'Script imports cubit but does not call cubit.init(). '
				"Add: cubit.init(['cubit', '-nojournal', '-batch'])"
			),
		})
	return findings


def check_wrong_connectivity_for_2nd_order(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: get_connectivity used when 2nd order elements are present."""
	findings = []
	has_2nd_order = any(
		re.search(r'element\s+type\s+(?:tetra10|hex20|hex27|wedge15|tri6|quad[89])', line, re.IGNORECASE)
		for line in lines
	)
	if not has_2nd_order:
		return findings

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		if 'get_connectivity(' in stripped and 'get_expanded_connectivity(' not in stripped:
			findings.append({
				'line': i,
				'severity': 'HIGH',
				'rule': 'wrong-connectivity-2nd-order',
				'message': (
					'get_connectivity() returns only corner nodes (1st order). '
					'For 2nd order elements, use get_expanded_connectivity() instead.'
				),
			})
	return findings


def check_element_type_before_add(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Element type set before elements added to block."""
	findings = []
	elem_type_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+element\s+type\b',
		re.IGNORECASE
	)
	block_add_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+add\b',
		re.IGNORECASE
	)

	elem_type_lines = {}  # block_id -> first line with element type
	block_add_lines = {}  # block_id -> first line with add

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue

		match = elem_type_pattern.search(stripped)
		if match:
			bid = match.group(1)
			if bid not in elem_type_lines:
				elem_type_lines[bid] = i

		match = block_add_pattern.search(stripped)
		if match:
			bid = match.group(1)
			if bid not in block_add_lines:
				block_add_lines[bid] = i

	for bid, type_line in elem_type_lines.items():
		if bid in block_add_lines and type_line < block_add_lines[bid]:
			findings.append({
				'line': type_line,
				'severity': 'HIGH',
				'rule': 'element-type-before-add',
				'message': (
					f'Block {bid} element type set at line {type_line} before '
					f'elements added at line {block_add_lines[bid]}. '
					f'Element type must be set AFTER adding elements to the block.'
				),
			})
	return findings


def check_missing_step_reimport(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Deprecated API usage detection."""
	findings = []
	# Detect usage of deleted APIs
	deleted_apis = [
		('export_netgen', 'export_NGSolveCurvedMesh'),
		('export_NetgenMesh', 'export_NGSolveCurvedMesh'),
		('export_netgen_with_names', 'export_NGSolveCurvedMesh'),
		('name_occ_faces', 'export_NGSolveCurvedMesh (no OCC needed)'),
		('set_cylinder_geominfo', 'export_NGSolveCurvedMesh'),
		('set_sphere_geominfo', 'export_NGSolveCurvedMesh'),
		('set_torus_geominfo', 'export_NGSolveCurvedMesh'),
		('set_cone_geominfo', 'export_NGSolveCurvedMesh'),
	]
	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		for old_api, replacement in deleted_apis:
			if old_api in stripped:
				findings.append({
					'line': i,
					'severity': 'CRITICAL',
					'rule': 'deleted-api-usage',
					'message': (
						f'{old_api}() has been DELETED. '
						f'Use {replacement}() instead. '
						'See: netgen_workflow_guide("deleted_apis")'
					),
				})
				break
	return findings


def check_hardcoded_absolute_paths(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Hardcoded absolute paths in sys.path or path assignments."""
	findings = []
	patterns = [
		re.compile(r'sys\.path\.(?:insert|append)\s*\(\s*(?:\d+\s*,\s*)?r?["\'][A-Za-z]:\\'),
		re.compile(r'sys\.path\.(?:insert|append)\s*\(\s*(?:\d+\s*,\s*)?r?["\'][A-Za-z]:/'),
	]
	# Allow known Cubit and NGSolve paths (common in examples)
	allow_patterns = [
		re.compile(r'Coreform Cubit', re.IGNORECASE),
		re.compile(r'NGSolve', re.IGNORECASE),
		re.compile(r'os\.environ\.get\s*\(\s*["\']CUBIT_PATH["\']'),  # CUBIT_PATH env var pattern
	]

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		for pattern in patterns:
			if pattern.search(stripped):
				# Skip allowed paths
				if any(ap.search(stripped) for ap in allow_patterns):
					continue
				findings.append({
					'line': i,
					'severity': 'MODERATE',
					'rule': 'hardcoded-absolute-path',
					'message': (
						'Hardcoded absolute path detected. '
						'Use os.environ.get("CUBIT_PATH", "...") for Cubit paths '
						'or os.path.dirname(__file__) for relative paths.'
					),
				})
				break
	return findings


def check_missing_boundary_block(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Volume element blocks but no surface element blocks."""
	findings = []
	has_3d_block = any(
		re.search(r'block\s+\d+\s+add\s+(?:tet|hex|wedge|pyramid)\b', line, re.IGNORECASE)
		for line in lines
	)
	has_2d_block = any(
		re.search(r'block\s+\d+\s+add\s+(?:tri|quad|face)\b', line, re.IGNORECASE)
		for line in lines
	)
	# Also check geometry-based surface blocks
	has_surface_block = any(
		re.search(r'block\s+\d+\s+(?:add\s+)?surface\b', line, re.IGNORECASE)
		for line in lines
	)
	# Only flag for Netgen/curved exports (boundary elements most important)
	has_netgen_export = any(
		re.search(r'export_NGSolveCurvedMesh', line)
		for line in lines
	)

	if has_3d_block and not has_2d_block and not has_surface_block and has_netgen_export:
		findings.append({
			'line': 1,
			'severity': 'MODERATE',
			'rule': 'missing-boundary-block',
			'message': (
				'Volume element block found but no surface element block. '
				'For Netgen export with boundary conditions, add: '
				'cubit.cmd("block 2 add tri all")'
			),
		})
	return findings


def check_missing_mesh_command(filepath: str, lines: List[str]) -> List[Dict]:
	"""CRITICAL: Export called but no mesh command found."""
	findings = []
	export_pattern = re.compile(
		r'(?:cubit_mesh_export\.|cme\.)?export_'
		r'(?:Gmsh_ver[24]|gmsh_v[24]|Nastran|nastran|'
		r'curved|exodus)\s*\('
	)
	mesh_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\'].*\b(?:mesh\s+(?:volume|surface|curve)|tetmesh|trimesh)\b',
		re.IGNORECASE
	)

	export_line = None
	has_mesh = False

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		if export_pattern.search(stripped) and export_line is None:
			export_line = i
		if mesh_pattern.search(stripped):
			has_mesh = True

	if export_line is not None and not has_mesh:
		findings.append({
			'line': export_line,
			'severity': 'CRITICAL',
			'rule': 'missing-mesh-command',
			'message': (
				'Export function called but no mesh command found. '
				'Add: cubit.cmd("mesh volume all") before export.'
			),
		})
	return findings


def check_export_file_extension(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Export file extension doesn't match format."""
	findings = []
	# Map export functions to expected extensions
	extension_map = {
		r'export_(?:Gmsh_ver[24]|gmsh_v[24])': ['.msh'],
		r'export_(?:Nastran|nastran)': ['.bdf', '.nas', '.dat'],
				r'export_exodus': ['.exo', '.e', '.ex2'],
	}

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		for func_pattern, valid_exts in extension_map.items():
			match = re.search(
				rf'(?:cubit_mesh_export\.|cme\.)?{func_pattern}\s*\(\s*cubit\s*,\s*["\']([^"\']+)["\']',
				stripped
			)
			if match:
				filename = match.group(1)
				_, ext = os.path.splitext(filename)
				ext = ext.lower()
				if ext and ext not in valid_exts:
					findings.append({
						'line': i,
						'severity': 'MODERATE',
						'rule': 'wrong-file-extension',
						'message': (
							f'File "{filename}" has extension "{ext}" but '
							f'expected {" or ".join(valid_exts)} for this export format.'
						),
					})
	return findings


def check_noheal_for_named_workflow(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Deleted name-based workflow detection."""
	# name_occ_faces and export_netgen_with_names are deleted.
	# Detection is handled by check_missing_step_reimport.
	return []


def check_setgeominfo_without_geometry(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Deleted SetGeomInfo API usage detected."""
	# This rule now detects usage of deleted SetGeomInfo functions
	# The actual detection is handled by check_missing_step_reimport
	# which catches all deleted APIs. This function is kept for
	# backward compatibility but returns empty.
	return []


def check_curve_without_setgeominfo(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: Manual mesh.Curve() without export_NGSolveCurvedMesh (legacy pattern)."""
	findings = []
	has_curve = False
	curve_line = None
	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		if re.search(r'\.Curve\s*\(\s*[2-9]', stripped):
			has_curve = True
			if curve_line is None:
				curve_line = i

	if not has_curve:
		return findings

	# Check for export_NGSolveCurvedMesh (which handles Curve internally)
	has_export_NGSolveCurvedMesh = any('export_NGSolveCurvedMesh' in line for line in lines)
	# Check for OCC native mesh (which has built-in geometry)
	has_occ_mesh = any('GenerateMesh' in line for line in lines)

	if not has_export_NGSolveCurvedMesh and not has_occ_mesh:
		findings.append({
			'line': curve_line,
			'severity': 'MODERATE',
			'rule': 'curve-without-export-curved',
			'message': (
				'mesh.Curve() called but no export_NGSolveCurvedMesh() found. '
				'For Cubit meshes, use export_NGSolveCurvedMesh(cubit, order=N) which '
				'handles curving automatically via ACIS CallbackGeometry. '
				'See: netgen_workflow_guide("export_NGSolveCurvedMesh")'
			),
		})
	return findings


def check_nodeset_sideset_usage(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: nodeset/sideset used with non-Exodus export functions."""
	findings = []
	nodeset_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\'].*\b(?:nodeset|sideset)\b', re.IGNORECASE
	)
	non_exodus_export = re.compile(
		r'(?:cubit_mesh_export\.|cme\.)?export_'
		r'(?:Gmsh_ver[24]|gmsh_v[24]|Nastran|nastran|'
		r'NetgenMesh|netgen|netgen_with_names)\s*\('
	)

	nodeset_lines = []
	has_non_exodus_export = False

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		if nodeset_pattern.search(stripped):
			nodeset_lines.append(i)
		if non_exodus_export.search(stripped):
			has_non_exodus_export = True

	if nodeset_lines and has_non_exodus_export:
		findings.append({
			'line': nodeset_lines[0],
			'severity': 'HIGH',
			'rule': 'nodeset-sideset-usage',
			'message': (
				'nodeset/sideset commands found but export uses '
				'blocks-only policy. Nodesets and sidesets are ignored by '
				'non-Exodus export functions. Use blocks for boundary conditions: '
				'cubit.cmd("block N add tri all in surface X"). '
				'See: cubit_scripting_guide("blocks_only_policy")'
			),
		})
	return findings


def check_missing_name_occ_faces(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Deleted API detection (handled elsewhere)."""
	# export_netgen_with_names and name_occ_faces are deleted.
	# Detection is handled by check_missing_step_reimport.
	return []


def check_missing_block_names(filepath: str, lines: List[str]) -> List[Dict]:
	"""LOW: Blocks registered without names (best practice)."""
	findings = []
	# Only check files that have export calls
	export_pattern = re.compile(
		r'(?:cubit_mesh_export\.|cme\.)?export_'
		r'(?:Gmsh_ver[24]|gmsh_v[24]|Nastran|nastran|'
		r'curved|exodus)\s*\('
	)
	has_export = any(
		export_pattern.search(line) for line in lines
		if not line.strip().startswith('#')
	)
	if not has_export:
		return findings

	block_add_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+add\b', re.IGNORECASE
	)
	block_name_pattern = re.compile(
		r'cubit\.cmd\s*\(\s*["\']block\s+(\d+)\s+name\b', re.IGNORECASE
	)

	added_blocks = set()
	named_blocks = set()

	for line in lines:
		stripped = line.strip()
		if stripped.startswith('#'):
			continue
		match = block_add_pattern.search(stripped)
		if match:
			added_blocks.add(match.group(1))
		match = block_name_pattern.search(stripped)
		if match:
			named_blocks.add(match.group(1))

	unnamed = added_blocks - named_blocks
	if unnamed:
		unnamed_sorted = sorted(unnamed, key=int)
		findings.append({
			'line': 1,
			'severity': 'LOW',
			'rule': 'missing-block-names',
			'message': (
				f'Block(s) {", ".join(unnamed_sorted)} registered without names. '
				f'Best practice: cubit.cmd(\'block {unnamed_sorted[0]} name "domain"\')'
			),
		})
	return findings


def check_non_ascii_bytes(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Non-ASCII bytes in Python file loaded by Cubit (cp932 crash).

	Cubit's embedded Python uses cp932 on Japanese Windows.
	Any non-ASCII byte (em dash, smart quotes, etc.) causes
	UnicodeDecodeError when startup.py reads register_toolbar.py.
	"""
	findings = []
	try:
		with open(filepath, 'rb') as f:
			raw = f.read()
	except (OSError, IOError):
		return findings

	for i, b in enumerate(raw):
		if b > 127:
			# Find line number
			line_num = raw[:i].count(b'\n') + 1
			context = raw[max(0, i - 10):i + 10]
			findings.append({
				'line': line_num,
				'severity': 'HIGH',
				'rule': 'non-ascii-byte',
				'message': (
					f'Non-ASCII byte 0x{b:02x} found. '
					f'Cubit cp932 will crash. Use ASCII only. '
					f'Context: {context!r}'
				),
			})
			break  # One finding is enough
	return findings


def check_qt_imports(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: Qt class used but not imported (NameError at runtime).

	Checks that all Qt classes referenced in the code are present
	in PySide6/PyQt5 import statements.
	"""
	findings = []
	qt_classes_used = set()
	qt_classes_imported = set()

	# Common Qt classes used in Cubit panels
	known_qt = {
		'QAction', 'QApplication', 'QMainWindow', 'QMenu', 'QMenuBar',
		'QMessageBox', 'QToolBar', 'QWidget', 'QDialog', 'QFileDialog',
		'QLabel', 'QPushButton', 'QVBoxLayout', 'QHBoxLayout',
		'QFormLayout', 'QLineEdit', 'QComboBox', 'QSpinBox',
		'QTableWidget', 'QDialogButtonBox', 'QProcess', 'QTimer',
	}

	import_pattern = re.compile(r'from\s+(?:PySide6|PyQt5)\.\w+\s+import\s+(.+)')

	for i, line in enumerate(lines, 1):
		stripped = line.strip()
		if stripped.startswith('#'):
			continue

		# Collect imports
		m = import_pattern.match(stripped)
		if m:
			names = m.group(1)
			# Handle multi-line imports (parenthesized)
			names = names.replace('(', '').replace(')', '').replace('\\', '')
			for name in names.split(','):
				name = name.strip()
				if name:
					qt_classes_imported.add(name)

		# Collect usage
		for cls in known_qt:
			if cls in stripped and cls not in qt_classes_imported:
				# Check it's actually used as a class (not in a string/comment)
				if re.search(rf'\b{cls}\s*\(', stripped):
					qt_classes_used.add((cls, i))

	for cls, line_num in sorted(qt_classes_used):
		if cls not in qt_classes_imported:
			findings.append({
				'line': line_num,
				'severity': 'HIGH',
				'rule': 'missing-qt-import',
				'message': f'{cls} used but not imported from PySide6/PyQt5.',
			})
	return findings


# All rules in execution order
ALL_RULES = [
	check_missing_block_registration,
	check_missing_mesh_command,
	check_geometry_block_for_2nd_order,
	check_missing_cubit_init,
	check_wrong_connectivity_for_2nd_order,
	check_element_type_before_add,
	check_setgeominfo_without_geometry,
	check_nodeset_sideset_usage,
	check_missing_name_occ_faces,
	check_missing_step_reimport,
	check_noheal_for_named_workflow,
	check_hardcoded_absolute_paths,
	check_missing_boundary_block,
	check_export_file_extension,
	check_curve_without_setgeominfo,
	check_missing_block_names,
	check_non_ascii_bytes,
	check_qt_imports,
]
