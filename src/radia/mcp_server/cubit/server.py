"""
Cubit Mesh Export MCP Server

Provides tools for:
- Export function documentation (Gmsh, Nastran, Netgen, Exodus)
- Netgen/NGSolve high-order curving workflows (SetGeomInfo, name-based)
- Cubit scripting patterns and best practices
- Cubit Python API reference (600+ functions, entity classes)
- Linting Python scripts for Cubit export convention violations

Usage:
    python server.py              # Start MCP server (stdio transport)
    python server.py --selftest   # Run self-test on examples/
"""

import os
import sys
from pathlib import Path

from mcp.server.fastmcp import FastMCP

# Import knowledge bases and rules (relative imports for pip package)
from .rules import ALL_RULES
from .export_knowledge import get_export_documentation
from .netgen_workflow_knowledge import get_netgen_documentation
from .cubit_scripting_knowledge import get_cubit_documentation
from .cubit_forum_tips import get_forum_tips
from .cubit_api_reference import get_api_reference
from .panel_conventions_knowledge import PANEL_CONVENTIONS, LABEL_GUIDE

# Create MCP server
mcp = FastMCP("cubit-export")

# Project root (current working directory for pip-installed package)
PROJECT_ROOT = Path.cwd()


# ============================================================
# Lint helpers
# ============================================================

def _lint_file(filepath: str) -> list[dict]:
	"""Run all lint rules on a single file."""
	try:
		with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
			lines = f.readlines()
	except (OSError, IOError) as e:
		return [{'line': 0, 'severity': 'ERROR', 'rule': 'read-error',
		         'message': f'Cannot read file: {e}'}]

	findings = []
	for rule_fn in ALL_RULES:
		findings.extend(rule_fn(filepath, lines))

	# Sort by severity, then line number
	severity_order = {'CRITICAL': 0, 'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'ERROR': -1}
	findings.sort(key=lambda f: (severity_order.get(f['severity'], 9), f['line']))
	return findings


def _format_findings(filepath: str, findings: list[dict]) -> str:
	"""Format findings for display."""
	if not findings:
		return f"[OK] {filepath}: No issues found."

	lines = [f"[{len(findings)} issue(s)] {filepath}:"]
	for f in findings:
		lines.append(
			f"  L{f['line']:>4d} [{f['severity']}] {f['rule']}: {f['message']}"
		)
	return '\n'.join(lines)


# ============================================================
# Netgen code examples
# ============================================================

EXAMPLES = {}

EXAMPLES['cylinder'] = '''# Cylinder with radia_export netgen (high-order curving)
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF, BND
import cubit

R, H, ORDER = 0.5, 2.0, 3

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

# Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Export with curving via C++ APREPRO command
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['sphere'] = '''# Sphere with radia_export netgen
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

R, ORDER = 0.5, 3

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {R}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 4/3 * math.pi * R**3
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['torus'] = '''# Torus with radia_export netgen
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

R_MAJOR, R_MINOR, ORDER = 1.0, 0.3, 3

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create torus major {R_MAJOR} minor {R_MINOR}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.08")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 2 * math.pi**2 * R_MAJOR * R_MINOR**2
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['cone'] = '''# Cone with radia_export netgen
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

R_BASE, H, ORDER = 0.5, 2.0, 3

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cone height {H} radius {R_BASE} top 0")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 1/3 * math.pi * R_BASE**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['hex_cylinder'] = '''# Hex-meshed Cylinder with radia_export netgen
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

R, H = 0.5, 2.0

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

# Hex meshing requires webcut for cylinder
cubit.cmd("webcut volume all with plane xplane")
cubit.cmd("webcut volume all with plane yplane")
cubit.cmd("merge all")
cubit.cmd("volume all scheme auto")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")

cubit.cmd("block 1 add volume all")

vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['complex_named'] = '''# Complex Geometry with Boolean Operations (radia_export netgen)
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

BRICK_SIZE, R_HOLE, ORDER = 2.0, 0.3, 3

# Create geometry directly in Cubit (no OCC needed!)
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create brick x {BRICK_SIZE} y {BRICK_SIZE} z {BRICK_SIZE}")
cubit.cmd(f"create cylinder height {BRICK_SIZE*2} radius {R_HOLE}")
cubit.cmd("subtract volume 2 from volume 1")

# Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Export with curving via C++ APREPRO command
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = BRICK_SIZE**3 - math.pi * R_HOLE**2 * BRICK_SIZE
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['vol_2nd_order'] = '''# Netgen .vol 2nd Order Export (Recommended Workflow)
import sys, os, math
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

R = 1.0

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {R}")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Register blocks
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')

# Export to Netgen .vol (order 2 with ACIS curving)
cubit.cmd('radia_export netgen "sphere.vol" order 2 overwrite')

# Read into NGSolve
mesh = Mesh("sphere.vol")

expected_vol = 4/3 * math.pi * R**3
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['poisson_sphere'] = '''# Complete FEM Example: Poisson on Sphere
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import *
import cubit

R, ORDER = 1.0, 3

# --- Mesh generation ---
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create sphere radius {R}")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.2")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# --- Export with curving via C++ APREPRO command ---
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

# --- Solve Poisson equation ---
fes = H1(mesh, order=ORDER, dirichlet="boundary")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += 1 * v * dx
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

print(f"DOFs: {fes.ndof}")
print(f"Integral of solution: {Integrate(gfu, mesh):.6f}")
'''


# ============================================================
# Tools
# ============================================================

@mcp.tool()
def cubit_docs(topic: str = "all") -> str:
	"""
	Get Cubit documentation: export formats, scripting guide, and API reference.

	Unified documentation tool covering export functions, Python scripting
	patterns, and the Cubit API (600+ functions). Use topic prefixes to
	navigate: export_*, scripting_*, api_*.

	Args:
	    topic: Documentation topic. Options:
	        "all"                    - Overview of all categories
	        --- Export formats ---
	        "export_overview"        - All export commands overview
	        "gmsh_v2"               - Gmsh v4.1 format (radia_export gmsh)
	        "gmsh_v4"               - Gmsh v4.1 format policy
	        "netgen"                - Netgen .vol export (radia_export netgen)
	        "nastran"               - Nastran BDF (radia_export nastran)
	        "exodus"                - Exodus II (Cubit built-in)
	        "export_comparison"     - Format comparison and decision matrix
	        "export_decision_guide" - Decision tree for format selection
	        --- Scripting guide ---
	        "scripting_overview"     - Why Cubit, element types, workflow
	        "scripting_lab_policy"   - Role split vs build123d (hex vs tet), translation guidance
	        "scripting_build123d_crossref" - Cubit `.jou` ↔ build123d mapping table
	        "scripting_blocks"       - Block registration
	        "scripting_blocks_only_policy" - Why blocks only (no nodesets/sidesets)
	        "scripting_element_order"  - 1st/2nd order, connectivity APIs
	        "scripting_mesh_schemes"   - tetmesh, map, sweep, trimesh
	        "scripting_step_exchange"  - STEP import/export, heal vs noheal
	        "scripting_initialization" - Cubit Python init boilerplate
	        "scripting_mixed_elements" - Mixed element types and warnings
	        "scripting_common_mistakes" - Frequent Cubit API pitfalls
	        "scripting_hex_workflow"    - Hexahedral mesh generation
	        "scripting_tet_workflow"    - Tetrahedral mesh generation
	        "scripting_2d_mesh"         - 2D surface meshing and export
	        "scripting_troubleshooting" - Common problems and debugging
	        "scripting_design_philosophy" - Why module reads Python API
	        "scripting_aprepro_journal"  - Running scripts via play command
	        --- API reference ---
	        "api_core"               - cmd, parse_cubit_list, init, naming
	        "api_geometry_queries"   - center_point, bounding_box, normals
	        "api_mesh_access"        - connectivity, nodal_coordinates
	        "api_blocks_sets"        - block/nodeset/sideset/group functions
	        "api_mesh_settings"      - scheme, size, intervals, quality
	        "api_entity_classes"     - Volume, Surface, Curve, Vertex
	        "api_graphics_selection" - Graphics control, selection
	        "api_advanced"           - Merge detection, geometry analysis
	        --- Radia-NGSolve panels ---
	        "panel_conventions"      - Analysis window conventions (TITLE, LABELS, etc.)
	        "panel_labels"           - Label guide (blocks/sidesets -> .vol -> NGSolve)
	"""
	topic = topic.lower().strip()

	if topic == "all":
		return (
			"# Cubit Documentation\n\n"
			"Use topic prefixes to navigate:\n"
			"- `export_*` - Export format documentation (gmsh_v2, netgen, nastran, etc.)\n"
			"- `scripting_*` - Python scripting guide (blocks, mesh_schemes, etc.)\n"
			"- `api_*` - API reference (core, geometry_queries, mesh_access, etc.)\n"
			"- `panel_*` - Radia-NGSolve analysis window conventions\n\n"
			+ get_export_documentation("overview")
		)

	# Export topics: strip prefix
	if topic.startswith("export_"):
		return get_export_documentation(topic[7:])  # remove "export_"

	# Scripting topics: strip prefix
	if topic.startswith("scripting_"):
		return get_cubit_documentation(topic[10:])  # remove "scripting_"

	# API topics: strip prefix
	if topic.startswith("api_"):
		return get_api_reference(topic[4:])  # remove "api_"

	# Panel topics
	if topic == "panel_conventions":
		return PANEL_CONVENTIONS
	if topic == "panel_labels":
		return LABEL_GUIDE

	# Try without prefix (backward compat for simple names)
	result = get_export_documentation(topic)
	if not result.startswith("Unknown"):
		return result
	result = get_cubit_documentation(topic)
	if not result.startswith("Unknown"):
		return result
	result = get_api_reference(topic)
	if not result.startswith("Unknown"):
		return result

	return (
		f"Unknown topic: '{topic}'. Use prefixes: export_*, scripting_*, api_*. "
		f"Examples: export_overview, scripting_blocks, api_core"
	)


@mcp.tool()
def netgen_workflow_guide(workflow: str = "overview") -> str:
	"""
	Get step-by-step Netgen/NGSolve workflow documentation.

	The Netgen high-order curving workflow is complex, involving Cubit meshing,
	OCC geometry, STEP exchange, and SetGeomInfo UV parameters. This tool
	provides detailed guidance for each workflow variant.

	Args:
	    workflow: Workflow to document. Options:
	        "overview"          - Decision tree: which workflow to use
	        "radia_export netgen"     - radia_export netgen() API reference
	        "simple_cylinder"   - Cylinder example
	        "simple_sphere"     - Sphere example
	        "simple_torus"      - Torus example
	        "complex"           - Complex geometry (Boolean operations, any shape)
	        "accuracy"          - Accuracy guide: order selection and verification
	        "kelvin_auto"       - Kelvin auto-add (auto R, symmetry, mesh copy)
	        "gmsh_2nd_order"    - Alternative: Gmsh 2nd order (simplest)
	        "troubleshooting"   - Common errors and fixes
	        "deleted_apis"      - Migration guide from old APIs
	"""
	return get_netgen_documentation(workflow)


@mcp.tool()
def netgen_code_example(shape: str = "cylinder") -> str:
	"""
	Get a ready-to-run Netgen export code example.

	Returns a complete Python script for exporting Cubit mesh to Netgen
	with high-order curving. Each example follows the recommended workflow.

	Args:
	    shape: Geometry shape for the example. Options:
	        "cylinder"       - Cylinder with radia_export netgen
	        "sphere"         - Sphere with radia_export netgen
	        "torus"          - Torus with radia_export netgen
	        "cone"           - Cone with radia_export netgen
	        "hex_cylinder"   - Hex-meshed cylinder with radia_export netgen
	        "complex_named"  - Boolean geometry with radia_export netgen
	        "gmsh_2nd_order" - Gmsh 2nd order alternative (simplest)
	        "poisson_sphere" - Complete FEM: Poisson equation on sphere
	"""
	shape = shape.lower().strip()
	if shape in EXAMPLES:
		return EXAMPLES[shape]
	else:
		return (
			f"Unknown shape: '{shape}'. "
			f"Available: {', '.join(EXAMPLES.keys())}"
		)


@mcp.tool()
def cubit_forum_tips(topic: str = "all") -> str:
	"""
	Get practical Cubit meshing tips sourced from the Coreform forum.

	Real-world solutions to common meshing problems, collected from
	forum.coreform.com discussions with high view counts.

	Args:
	    topic: Tip category. Options:
	        "all"                  - All tips
	        "mesh_quality"         - Smoothing, Jacobian fixes, quality diagnostics
	        "biased_graded"        - Curve bias, circle scheme, auto factor, refinement
	        "hex_decomposition"    - Webcut patterns, torus meshing, embedded spheres
	        "sweep_tips"           - Transform translate, copy mesh, interval matching
	        "parallel_meshing"     - HPC tetmesher, multiprocessing pattern
	        "python_api"           - Entity names, connectivity, batch mode, STL engine
	        "export_tips"          - Abaqus parts, Exodus sizing, LS-DYNA
	        "geometry_healing"     - STEP/STL import, heal, stitch, bounding voids
	        "sculpt_workflow"      - Sculpt setup, command-line, capture, HEX27
	        "boundary_layer"       - BL setup, reversal nodes fix, composite surface
	        "python_performance"   - Performance settings, batch queries, entity tracking
	        "quality_diagnostics"  - Quality extraction, bad element selection, normals
	        "solver_workflows"     - FEniCS, VTK, CalculiX, OpenFOAM integration
	        "advanced_meshing"     - Void, crack, dome, hemisphere, topography, THex

	    For troubleshooting, use cubit_docs(topic="scripting_troubleshooting").
	"""
	return get_forum_tips(topic)


@mcp.tool()
def lint_cubit_script(filepath: str) -> str:
	"""
	Lint a Python script for Cubit mesh export convention violations.

	Checks for (16 rules):
	- Missing block registration before export (CRITICAL)
	- Missing mesh command before export (CRITICAL)
	- Geometry block with 2nd order conversion (HIGH)
	- Missing cubit.init() (HIGH)
	- get_connectivity instead of get_expanded_connectivity for 2nd order (HIGH)
	- Element type set before adding to block (HIGH)
	- Deleted SetGeomInfo API usage (CRITICAL)
	- nodeset/sideset usage with blocks-only export formats (HIGH)
	- Deleted API usage (export_netgen, set_*_geominfo, name_occ_faces) (CRITICAL)
	- Deleted STEP reimport / name-based workflow APIs (CRITICAL)
	- Hardcoded absolute paths (MODERATE)
	- Missing boundary element block for Netgen (MODERATE)
	- Wrong file extension for export format (MODERATE)
	- mesh.Curve() without radia_export netgen (MODERATE)
	- Blocks registered without names (LOW)

	Args:
	    filepath: Absolute or relative path to the Python file to check.
	"""
	p = Path(filepath)
	if not p.is_absolute():
		p = PROJECT_ROOT / p

	if not p.exists():
		return f"Error: File not found: {p}"
	if not p.suffix == '.py':
		return f"Error: Not a Python file: {p}"

	findings = _lint_file(str(p))
	return _format_findings(str(p), findings)


@mcp.tool()
def lint_cubit_directory(directory: str = "examples") -> str:
	"""
	Lint all Python scripts in a directory for Cubit export convention violations.

	Recursively scans .py files and reports findings grouped by file.

	Args:
	    directory: Directory path relative to project root (default: "examples").
	"""
	d = Path(directory)
	if not d.is_absolute():
		d = PROJECT_ROOT / d

	if not d.exists():
		return f"Error: Directory not found: {d}"

	# Collect all .py files
	py_files = sorted(d.rglob("*.py"))

	if not py_files:
		return f"No Python files found in {directory}."

	# Lint all files
	total_findings = 0
	file_results = []
	summary_by_severity = {'CRITICAL': 0, 'HIGH': 0, 'MODERATE': 0, 'LOW': 0}

	for py_file in py_files:
		findings = _lint_file(str(py_file))
		if findings:
			total_findings += len(findings)
			rel_path = py_file.relative_to(PROJECT_ROOT) if py_file.is_relative_to(PROJECT_ROOT) else py_file
			file_results.append(_format_findings(str(rel_path), findings))
			for f in findings:
				sev = f['severity']
				if sev in summary_by_severity:
					summary_by_severity[sev] += 1

	# Build output
	output_parts = [
		f"Cubit Export Lint Report: {len(py_files)} files scanned, {total_findings} issues found.",
		"",
		f"Summary: {summary_by_severity['CRITICAL']} CRITICAL, "
		f"{summary_by_severity['HIGH']} HIGH, "
		f"{summary_by_severity['MODERATE']} MODERATE, "
		f"{summary_by_severity['LOW']} LOW",
		"",
	]

	if file_results:
		output_parts.append("=" * 70)
		output_parts.extend(file_results)
	else:
		output_parts.append("All files passed!")

	return '\n'.join(output_parts)


# ============================================================
# New tools: script generation and mesh quality
# ============================================================

SCRIPT_TEMPLATES = {}

SCRIPT_TEMPLATES['tet_netgen'] = '''# Cubit -> NGSolve: Tet mesh with high-order curving (radia_export netgen)
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF, BND
import cubit

# === PARAMETERS ===
MESH_SIZE = 0.15
CURVE_ORDER = 3

# === GEOMETRY ===
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
# TODO: Create your geometry here
# cubit.cmd("create cylinder height 2 radius 0.5")

# === MESH ===
cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")

# === BLOCKS ===
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# === EXPORT WITH CURVING (C++ APREPRO command) ===
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
mesh = Mesh(vol_path)

# === VERIFY ===
vol = Integrate(CF(1), mesh)
area = Integrate(CF(1), mesh, VOL_or_BND=BND)
print(f"Volume: {vol:.6f}, Surface area: {area:.6f}")
'''

SCRIPT_TEMPLATES['tet_netgen_named'] = '''# Cubit -> NGSolve: Complex geometry with Boolean operations (radia_export netgen)
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF, BND
import cubit

MESH_SIZE = 0.15
CURVE_ORDER = 3

# === GEOMETRY (create directly in Cubit) ===
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
# TODO: Build your geometry using Cubit Boolean operations
# cubit.cmd("create brick x 2 y 2 z 2")
# cubit.cmd("create cylinder height 4 radius 0.3")
# cubit.cmd("subtract volume 2 from volume 1")

# === MESH ===
cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# === EXPORT WITH CURVING (C++ APREPRO command) ===
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
mesh = Mesh(vol_path)

vol = Integrate(CF(1), mesh)
print(f"Volume: {vol:.6f}")
'''

SCRIPT_TEMPLATES['tet_vol'] = '''# Cubit -> Netgen .vol 2nd order (recommended workflow)
import sys, os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

MESH_SIZE = 0.15

# === GEOMETRY AND MESH ===
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
# TODO: Create your geometry
# cubit.cmd("create sphere radius 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")

# === BLOCKS ===
cubit.cmd("block 1 add volume all")

# === EXPORT (.vol with order 2 curving) ===
cubit.cmd('radia_export netgen "mesh.vol" order 2 overwrite')

# === IMPORT TO NGSOLVE ===
mesh = Mesh("mesh.vol")
vol = Integrate(CF(1), mesh)
print(f"Volume: {vol:.6f}")
'''

SCRIPT_TEMPLATES['hex_netgen'] = '''# Cubit -> NGSolve: Hex mesh with high-order curving (radia_export netgen)
import sys, os, math, tempfile
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ngsolve import Mesh, Integrate, CF
import cubit

MESH_SIZE = 0.15
CURVE_ORDER = 3

# === GEOMETRY ===
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
# TODO: Create your geometry
# cubit.cmd("create cylinder height 2 radius 0.5")

# === HEX MESH (may need webcut for non-sweepable geometries) ===
# cubit.cmd("webcut volume all with plane xplane")
# cubit.cmd("webcut volume all with plane yplane")
# cubit.cmd("merge all")
cubit.cmd("volume all scheme auto")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")

# === BLOCKS ===
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# === EXPORT WITH CURVING (C++ APREPRO command) ===
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'radia_export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
mesh = Mesh(vol_path)

vol = Integrate(CF(1), mesh)
print(f"Volume: {vol:.6f}")
'''

MESH_QUALITY_KNOWLEDGE = """
# Mesh Quality Guide

## Quality Metrics in Cubit

```python
# Check overall mesh quality
cubit.cmd("quality volume all shape")        # Shape metric (0-1, higher=better)
cubit.cmd("quality volume all aspect ratio") # Aspect ratio (1=ideal)
cubit.cmd("quality volume all jacobian")     # Jacobian (must be >0)
cubit.cmd("quality volume all condition")    # Condition number
```

## Minimum Requirements by Application

| Application | Shape | Aspect Ratio | Jacobian |
|-------------|-------|-------------|----------|
| General FEM | > 0.2 | < 10 | > 0 |
| High-order curving | > 0.3 | < 5 | > 0 |
| Radia BEM (hex) | > 0.5 | < 3 | > 0 (convex) |
| CFD | > 0.3 | < 5 | > 0 |

## Improving Mesh Quality

### Tetrahedral Meshes
```python
# Reduce element size
cubit.cmd("volume 1 size 0.05")

# Use quality-based smoothing
cubit.cmd("volume 1 smooth scheme condition number beta 2 cpu 10")
cubit.cmd("smooth volume 1")

# Webcut for better element shapes
cubit.cmd("webcut volume 1 with plane xplane")
```

### Hexahedral Meshes
```python
# Interval control for uniformity
cubit.cmd("curve 1 interval 10")

# Auto scheme selection
cubit.cmd("volume 1 scheme auto")

# Smoothing after meshing
cubit.cmd("smooth volume 1")
```

## Element Count Guidelines

```python
# Check counts
n_tet = cubit.get_tet_count()
n_hex = cubit.get_hex_count()
n_tri = cubit.get_tri_count()
n_quad = cubit.get_quad_count()
print(f"3D: {n_tet} tet + {n_hex} hex, 2D: {n_tri} tri + {n_quad} quad")
```

| Element Count | Suitability |
|---------------|-------------|
| < 1,000 | Quick prototyping, coarse analysis |
| 1,000 - 10,000 | Standard analysis |
| 10,000 - 100,000 | Fine analysis, convergence studies |
| > 100,000 | Very fine, may need parallel solver |

## Common Quality Issues

### 1. Negative Jacobian (Invalid Elements)
- **Cause**: Collapsed or inverted elements
- **Fix**: Refine mesh, webcut complex regions, use tet instead of hex

### 2. High Aspect Ratio
- **Cause**: Non-uniform sizing, thin geometry features
- **Fix**: Local size control, bias on curves

### 3. Poor Shape Quality Near Curved Surfaces
- **Cause**: Mesh doesn't follow curvature well
- **Fix**: Smaller element size near curved surfaces
  ```python
  cubit.cmd("surface 1 size 0.02")  # Finer on specific surface
  ```

### 4. Transition Element Issues
- **Cause**: Size transition too abrupt
- **Fix**: Gradual size change
  ```python
  cubit.cmd("volume 1 sizing function type skeleton scale 1.5")
  ```
"""


@mcp.tool()
def generate_cubit_script(workflow: str = "tet_netgen") -> str:
	"""
	Generate a template Cubit Python script for common workflows.

	Returns a ready-to-customize script with TODO markers for
	geometry-specific sections. Follows all project conventions.

	Args:
	    workflow: Script template to generate. Options:
	        "tet_netgen"       - Tet mesh -> radia_export netgen with ACIS curving
	        "tet_netgen_named" - Complex geometry with radia_export netgen
	        "tet_gmsh"         - Tet mesh -> Netgen .vol 2nd order (simplest)
	        "hex_netgen"       - Hex mesh -> radia_export netgen with curving
	"""
	workflow = workflow.lower().strip()
	if workflow in SCRIPT_TEMPLATES:
		return SCRIPT_TEMPLATES[workflow]
	else:
		return (
			f"Unknown workflow: '{workflow}'. "
			f"Available: {', '.join(SCRIPT_TEMPLATES.keys())}"
		)


@mcp.tool()
def get_lint_rules() -> str:
	"""
	List all available Cubit export lint rules with descriptions.

	Returns a summary of each rule including its name, severity,
	description, and fix guidance. Useful for understanding what
	the linter checks before running lint_cubit_script.
	"""
	rules_info = [
		{
			'rule': 'missing-block-registration',
			'severity': 'CRITICAL',
			'description': 'Export function called but no block registration found.',
			'trigger': 'radia_export gmsh without any block/mesh',
			'fix': 'cubit.cmd("block 1 add tet all")',
		},
		{
			'rule': 'missing-mesh-command',
			'severity': 'CRITICAL',
			'description': 'Export function called but no mesh command found.',
			'trigger': 'radia_export without mesh volume/surface command',
			'fix': 'cubit.cmd("mesh volume all")',
		},
		{
			'rule': 'geometry-block-2nd-order',
			'severity': 'HIGH',
			'description': 'Geometry block (volume/surface) used with element type conversion. Has no effect.',
			'trigger': 'block 1 add volume 1 + block 1 element type tetra10',
			'fix': 'Use mesh elements: block 1 add tet all in volume 1',
		},
		{
			'rule': 'missing-cubit-init',
			'severity': 'HIGH',
			'description': 'Script imports cubit but does not call cubit.init().',
			'trigger': 'import cubit without cubit.init()',
			'fix': "cubit.init(['cubit', '-nojournal', '-batch'])",
		},
		{
			'rule': 'wrong-connectivity-2nd-order',
			'severity': 'HIGH',
			'description': 'get_connectivity() used with 2nd order elements. Returns only corner nodes.',
			'trigger': 'block 1 element type tetra10 + get_connectivity("tet", tid)',
			'fix': 'Use get_expanded_connectivity("tet", tid)',
		},
		{
			'rule': 'element-type-before-add',
			'severity': 'HIGH',
			'description': 'Element type set before elements added to block. Type gets reset.',
			'trigger': 'block 1 element type tetra10 before block 1 add tet all',
			'fix': 'Add elements first, then set element type',
		},
		{
			'rule': 'deleted-api-usage',
			'severity': 'CRITICAL',
			'description': 'Deleted API (export_netgen, set_*_geominfo, name_occ_faces) detected.',
			'trigger': 'Any usage of export_netgen, export_NetgenMesh, set_*_geominfo, name_occ_faces',
			'fix': 'Use radia_export netgen "f.vol" order N',
		},
		{
			'rule': 'nodeset-sideset-usage',
			'severity': 'HIGH',
			'description': 'nodeset/sideset commands found with non-Exodus export. Blocks-only policy.',
			'trigger': 'nodeset/sideset with non-Exodus export',
			'fix': 'Use blocks: cubit.cmd("block 2 add tri all in surface 1")',
		},
		{
			'rule': 'deleted-api-usage',
			'severity': 'CRITICAL',
			'description': '(See above) Catches all deleted API usage.',
			'trigger': 'export_netgen_with_names, name_occ_faces',
			'fix': 'Use radia_export netgen "f.vol" order N',
		},
		{
			'rule': 'deleted-api-usage',
			'severity': 'CRITICAL',
			'description': '(See above) Catches export_netgen and other deleted APIs.',
			'trigger': 'export_netgen(cubit, geometry=geo)',
			'fix': 'Use radia_export netgen "f.vol" order N',
		},
		{
			'rule': 'deleted-api-noheal',
			'severity': 'INFO',
			'description': 'Name-based workflow (name_occ_faces + noheal) is deleted. Use radia_export netgen.',
			'trigger': 'name_occ_faces(shape) + STEP import',
			'fix': 'Use radia_export netgen(cubit, order=N) - no STEP needed',
		},
		{
			'rule': 'hardcoded-absolute-path',
			'severity': 'MODERATE',
			'description': 'Hardcoded absolute paths in sys.path (except Coreform Cubit/NGSolve).',
			'trigger': 'sys.path.insert(0, "S:/Projects/mymodule")',
			'fix': 'sys.path.insert(0, os.path.dirname(__file__))',
		},
		{
			'rule': 'missing-boundary-block',
			'severity': 'MODERATE',
			'description': 'Volume element block found but no surface element block for Netgen export.',
			'trigger': 'block 1 add tet all + radia_export netgen(...) without tri block',
			'fix': 'cubit.cmd("block 2 add tri all")',
		},
		{
			'rule': 'wrong-file-extension',
			'severity': 'MODERATE',
			'description': 'Export file extension does not match the format.',
			'trigger': 'radia_export gmsh "mesh.vtk" (wrong extension)',
			'fix': 'Use correct extension: .msh, .bdf, .exo, .vol',
		},
		{
			'rule': 'curve-without-export-curved',
			'severity': 'MODERATE',
			'description': 'Manual mesh.Curve() is not needed. radia_export netgen embeds curving in .vol.',
			'trigger': 'mesh.Curve(3) without radia_export netgen',
			'fix': 'Use cubit.cmd(\'radia_export netgen "mesh.vol" order 3 overwrite\') then mesh = Mesh("mesh.vol")',
		},
		{
			'rule': 'missing-block-names',
			'severity': 'LOW',
			'description': 'Blocks registered without name command. Named blocks improve readability.',
			'trigger': 'block 1 add tet all without block 1 name "..."',
			'fix': 'cubit.cmd(\'block 1 name "domain"\')',
		},
		{
			'rule': 'non-ascii-byte',
			'severity': 'HIGH',
			'description': 'Non-ASCII byte in Python file. Cubit cp932 (Japanese Windows) will crash.',
			'trigger': 'Em dash, smart quotes, or any non-ASCII character in .py loaded by Cubit',
			'fix': 'Replace with ASCII equivalents: -- instead of em dash, straight quotes, etc.',
		},
		{
			'rule': 'missing-qt-import',
			'severity': 'HIGH',
			'description': 'Qt class used but not imported. Causes NameError at runtime in Cubit.',
			'trigger': 'QMenu(...) without "from PySide6.QtWidgets import QMenu"',
			'fix': 'Add the missing class to PySide6/PyQt5 import statement.',
		},
	]

	lines = ["# Cubit Export Lint Rules (16 rules)", ""]
	for r in rules_info:
		lines.append(f"## [{r['severity']}] {r['rule']}")
		lines.append(f"{r['description']}")
		lines.append(f"- **Trigger**: {r['trigger']}")
		lines.append(f"- **Fix**: `{r['fix']}`")
		lines.append("")

	return '\n'.join(lines)


# ============================================================
# Preview tool (OCP CAD Viewer bridge via STEP)
# ============================================================

_DEFAULT_CUBIT_EXE = r"C:\Program Files\Coreform Cubit 2025.3\bin\coreform_cubit.exe"


@mcp.tool()
def preview_jou(jou_path: str,
                cubit_exe: str = "",
                keep_step: bool = False,
                timeout_s: int = 120) -> str:
	"""
	Run a Cubit journal (.jou) in batch mode, export the resulting
	geometry as STEP, load it via build123d, and show it in the OCP
	CAD Viewer (VSCode panel).

	Lab policy (2026-04-19): Cubit runs out-of-process (separate exe
	with multi-second startup), so STEP is the mandatory bridge to
	the OCP viewer. The symmetric build123d-side preview tool
	(`preview_shape` on mcp-server-build123d) goes **direct** with
	no STEP roundtrip.

	Prerequisites:
	  - Coreform Cubit installed; binary at `cubit_exe` (default:
	    Coreform Cubit 2025.3 in Program Files).
	  - `ocp_vscode` (OCP CAD Viewer) panel open in VSCode.

	Pipeline:
	  1. Write a wrapper .jou: `playback "<jou_path>"` + `export step
	     "<tmp>.step" overwrite`.
	  2. Invoke `coreform_cubit -batch -nographics -nojournal -input
	     wrapper.jou`. Capture stdout/stderr.
	  3. Verify STEP exists; `from build123d import import_step`.
	  4. `ocp_vscode.show(part)`.
	  5. Return summary JSON (Cubit exit + shape stats).

	Args:
	    jou_path: absolute or relative path to a .jou file.
	    cubit_exe: Path to coreform_cubit.exe (default: Coreform Cubit
	        2025.3 bin location). Override for other versions.
	    keep_step: keep the intermediate .step file after preview.
	    timeout_s: abort Cubit subprocess after this many seconds.
	"""
	import json as _json
	import subprocess as _sp
	import tempfile as _tf
	import traceback as _tb
	from pathlib import Path

	p = Path(jou_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	if not p.exists():
		return _json.dumps({"status": "error", "stage": "input",
		                    "error": f"Journal not found: {p}"})

	exe = cubit_exe or _DEFAULT_CUBIT_EXE
	if not Path(exe).exists():
		return _json.dumps({
			"status": "error", "stage": "cubit_binary",
			"error": f"Cubit exe not found: {exe}. Pass cubit_exe="
		             f"<path-to-coreform_cubit.exe>.",
		})

	work = Path(_tf.mkdtemp(prefix="cubit_preview_"))
	step_path = work / "preview.step"
	wrapper = work / "wrapper.jou"
	# Cubit journal expects forward slashes in paths.
	abs_jou = str(p.resolve()).replace("\\", "/")
	abs_step = str(step_path.resolve()).replace("\\", "/")
	wrapper.write_text(
		f'playback "{abs_jou}"\n'
		f'export step "{abs_step}" overwrite\n'
	)

	cmd = [exe, "-batch", "-nographics", "-nojournal",
	       "-input", str(wrapper)]
	try:
		proc = _sp.run(cmd, capture_output=True, text=True,
		               timeout=timeout_s)
	except _sp.TimeoutExpired:
		return _json.dumps({
			"status": "error", "stage": "cubit_run",
			"error": f"Cubit timed out after {timeout_s}s",
			"cmd": cmd,
		})

	if not step_path.exists():
		return _json.dumps({
			"status": "error", "stage": "cubit_export",
			"error": "Cubit did not produce STEP file",
			"cubit_exit": proc.returncode,
			"cubit_stdout_tail": proc.stdout[-800:],
			"cubit_stderr_tail": proc.stderr[-800:],
		})

	# Load STEP and hand to OCP viewer.
	try:
		from build123d import import_step
	except ImportError as e:
		return _json.dumps({
			"status": "error", "stage": "build123d_import",
			"error": f"build123d not installed: {e}",
		})
	try:
		shape = import_step(str(step_path))
	except Exception:
		return _json.dumps({
			"status": "error", "stage": "step_load",
			"error": _tb.format_exc(),
			"step_path": str(step_path),
		})

	try:
		from ocp_vscode import show
	except ImportError as e:
		return _json.dumps({
			"status": "error", "stage": "ocp_import",
			"error": f"ocp_vscode not installed: {e}. "
			         f"pip install ocp-vscode",
		})
	try:
		show(shape, names=[p.stem])
	except Exception:
		return _json.dumps({
			"status": "error", "stage": "show",
			"error": _tb.format_exc(),
			"hint": "Is the OCP CAD Viewer panel open in VSCode?",
			"step_path": str(step_path),
		})

	info = {
		"status": "ok",
		"stage": "shown",
		"viewer": "ocp_vscode",
		"jou": str(p),
		"cubit_exit": proc.returncode,
		"step_bytes": step_path.stat().st_size,
	}
	try:
		info["volume"] = round(shape.volume, 6)
	except Exception:
		pass
	try:
		info["face_count"] = len(shape.faces())
		info["edge_count"] = len(shape.edges())
	except Exception:
		pass
	try:
		bb = shape.bounding_box()
		info["bounding_box"] = {
			"min": [round(bb.min.X, 4), round(bb.min.Y, 4),
			        round(bb.min.Z, 4)],
			"max": [round(bb.max.X, 4), round(bb.max.Y, 4),
			        round(bb.max.Z, 4)],
		}
	except Exception:
		pass

	if keep_step:
		info["step_path"] = str(step_path)
	else:
		# Clean up temp work dir.
		try:
			step_path.unlink(missing_ok=True)
			wrapper.unlink(missing_ok=True)
			work.rmdir()
		except OSError:
			pass

	return _json.dumps(info, indent=2)


# ============================================================
# Cubit GUI launcher (Cheap phase of LLM-driven Cubit loop)
# ============================================================
#
# Motivation (2026-04-19): in the collaborative workflow (student +
# Claude Code in VSCode), the student wants to drive mesh creation via
# LLM but watch the result in Cubit's own native window (not OCP).
# Phase split:
#   - Cheap (this tool): each call launches a fresh Cubit GUI process
#     with a wrapper .jou that either imports a file or executes a
#     list of commands. Cubit stays open; student inspects, closes,
#     or triggers another round.  Slow loop but zero dev cost.
#   - Medium (TODO): use Cubit's Python binding
#     (C:/Program Files/Coreform Cubit 2025.3/bin/python3/) to drive a
#     persistent Cubit process in-process — sub-second command cycles.
#     Requires binding mode that keeps GUI alive; experimental.


@mcp.tool()
def open_in_cubit(path: str = "",
                  commands: list = None,
                  cubit_exe: str = "",
                  detach: bool = True) -> str:
	"""
	Open a file, or execute a list of commands, in **Cubit GUI**.
	This is the "Cheap phase" of the LLM→Cubit loop: each invocation
	launches a fresh Cubit GUI process so the student can watch the
	result in Cubit's native window. Use for inspection / measurement
	after Claude generates a .jou.

	File-type dispatch (by extension):
	  .jou  → `playback "<path>"`
	  .step .stp → `import step "<path>"`
	  .brep .brp → `import acis "<path>"`  (Cubit uses ACIS for BREP)
	  .msh         → `import mesh "<path>"`
	  .vol         → `import mesh "<path>"`  (Netgen .vol passes through)
	  .g .e .exo   → `import mesh "<path>"`
	Any `commands` given are appended AFTER the file load (useful for
	"import this STEP then mesh it").

	Args:
	    path: file to open in Cubit. Leave empty to run only `commands`.
	    commands: extra Cubit commands to run after loading. Each item
	        is one line (no trailing newline).
	    cubit_exe: override Cubit binary path.
	    detach: if True (default), return immediately after launching
	        Cubit; student interacts with the GUI window. If False,
	        wait for Cubit to exit and return its output (useful for
	        smoke tests or headless CI).

	Returns:
	    JSON with launch status, pid (when detach=True), wrapper .jou
	    path (kept for the student to inspect/edit), exit code (when
	    detach=False), and a copy of the commands that were issued.
	"""
	import json as _json
	import subprocess as _sp
	import tempfile as _tf
	from pathlib import Path

	exe = cubit_exe or _DEFAULT_CUBIT_EXE
	if not Path(exe).exists():
		return _json.dumps({"status": "error", "stage": "cubit_binary",
		                    "error": f"Cubit exe not found: {exe}"})

	work = Path(_tf.mkdtemp(prefix="cubit_gui_"))
	wrapper = work / "wrapper.jou"
	lines = []

	if path:
		p = Path(path)
		if not p.is_absolute():
			p = PROJECT_ROOT / p
		if not p.exists():
			return _json.dumps({"status": "error", "stage": "input",
			                    "error": f"File not found: {p}"})
		abs_path = str(p.resolve()).replace("\\", "/")
		suffix = p.suffix.lower()
		if suffix == ".jou":
			lines.append(f'playback "{abs_path}"')
		elif suffix in (".step", ".stp"):
			lines.append(f'import step "{abs_path}"')
		elif suffix in (".brep", ".brp"):
			# Cubit imports BREP via the ACIS kernel.
			lines.append(f'import acis "{abs_path}"')
		elif suffix in (".msh", ".vol", ".g", ".e", ".exo"):
			lines.append(f'import mesh "{abs_path}"')
		else:
			return _json.dumps({"status": "error", "stage": "input",
			                    "error": f"Unsupported extension: {suffix}"})

	if commands:
		lines.extend(str(c).rstrip() for c in commands)

	if not lines:
		return _json.dumps({"status": "error", "stage": "input",
		                    "error": "Either `path` or `commands` is required."})

	wrapper.write_text("\n".join(lines) + "\n")

	# GUI mode: no -batch, no -nographics. -nojournal prevents auto-
	# journaling of the wrapper playback (keeps the directory clean).
	# -input <wrapper> plays our commands after GUI initialization.
	cmd = [exe, "-nojournal", "-input", str(wrapper)]

	try:
		if detach:
			proc = _sp.Popen(
				cmd,
				stdout=_sp.DEVNULL,
				stderr=_sp.DEVNULL,
				creationflags=getattr(_sp, "DETACHED_PROCESS", 0)
				             | getattr(_sp, "CREATE_NEW_PROCESS_GROUP", 0),
			)
			info = {
				"status": "ok",
				"mode": "gui_detached",
				"pid": proc.pid,
				"wrapper": str(wrapper),
				"commands": lines,
				"note": ("Cubit GUI launched. Student: inspect in the "
				         "Cubit window. Close the Cubit window when done."),
			}
		else:
			proc = _sp.run(cmd, capture_output=True, text=True, timeout=180)
			info = {
				"status": "ok",
				"mode": "gui_blocking",
				"exit": proc.returncode,
				"wrapper": str(wrapper),
				"commands": lines,
				"stdout_tail": proc.stdout[-400:] if proc.stdout else "",
				"stderr_tail": proc.stderr[-400:] if proc.stderr else "",
			}
	except _sp.TimeoutExpired:
		return _json.dumps({"status": "error", "stage": "timeout",
		                    "wrapper": str(wrapper)})
	except OSError as e:
		return _json.dumps({"status": "error", "stage": "launch",
		                    "error": str(e)})

	return _json.dumps(info, indent=2)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_cubit_mesh(geometry: str, element_type: str = "tet") -> str:
	"""Create a Cubit meshing script for the given geometry."""
	return (
		f"Create a Cubit meshing script for: {geometry}\n"
		f"Element type: {element_type}\n\n"
		"Follow these conventions:\n"
		"1. cubit.init(['cubit', '-nojournal', '-batch'])\n"
		"2. Register blocks with names: block N add {tet|hex} all; block N name '...'\n"
		"3. Add boundary block: block N add {tri|quad} all\n"
		"4. Use relative paths for imports\n"
		"5. For curved mesh: cubit.cmd('radia_export netgen \"mesh.vol\" order N overwrite')\n"
		"6. Verify volume after curving\n"
		"7. For Gmsh visualization: cubit.cmd('radia_export gmsh \"mesh.msh\" order N overwrite')\n"
	)


@mcp.prompt()
def cubit_to_ngsolve(mesh_file: str) -> str:
	"""Set up Cubit mesh -> NGSolve high-order FEM workflow."""
	return (
		f"Set up a Cubit -> NGSolve high-order workflow for: {mesh_file}\n\n"
		"Workflow:\n"
		"1. Create geometry and mesh in Cubit\n"
		"2. cubit.cmd('radia_export netgen \"mesh.vol\" order N overwrite')\n"
		"3. mesh = Mesh('mesh.vol')  # already curved, no mesh.Curve() needed\n"
	)


# ============================================================
# MCP Resources
# ============================================================

@mcp.resource("cubit://export-decision")
def cubit_export_decision_guide() -> str:
	"""Quick decision guide for choosing Cubit export format."""
	return (
		"# Cubit Export Format Decision Guide\n\n"
		"| Need | Command | Max Order |\n"
		"|------|---------|----------|\n"
		"| NGSolve FEM (recommended) | radia_export netgen \"f.vol\" order N | 1-5 |\n"
		"| GMSH visualization | radia_export gmsh \"f.msh\" order N | 1-3 |\n"
		"| Nastran / JMAG | radia_export nastran \"f.bdf\" order N | 1-2 |\n"
		"| ParaView | radia_export vtk \"f.vtk\" order N | 1-2 |\n"
		"| ELF/MAGIC | export meg \"f.meg\" | 1 |\n"
		"| Cubit-native archival | export mesh \"f.exo\" | all |\n"
	)


@mcp.resource("cubit://element-types")
def cubit_element_types_reference() -> str:
	"""Cubit element type reference for block registration."""
	return (
		"# Cubit Element Types\n\n"
		"## 3D Elements\n"
		"| Type | Nodes | 2nd Order | Cubit Name |\n"
		"|------|-------|-----------|------------|\n"
		"| Tetrahedron | 4 | 10 (tetra10) | tet |\n"
		"| Hexahedron | 8 | 20 (hex20) / 27 (hex27) | hex |\n"
		"| Wedge | 6 | 15 (wedge15) | wedge |\n"
		"| Pyramid | 5 | 13 (pyramid13) | pyramid |\n\n"
		"## 2D Elements\n"
		"| Type | Nodes | 2nd Order | Cubit Name |\n"
		"|------|-------|-----------|------------|\n"
		"| Triangle | 3 | 6 (tri6) | tri |\n"
		"| Quadrilateral | 4 | 8 (quad8) / 9 (quad9) | quad |\n\n"
		"## Block Registration Pattern\n"
		"```python\n"
		"cubit.cmd('block 1 add tet all')     # Add elements first\n"
		"cubit.cmd('block 1 element type tetra10')  # Then set type\n"
		"cubit.cmd('block 1 name \"domain\"')  # Then name\n"
		"```\n"
	)


# ============================================================
# Self-test
# ============================================================

def _selftest():
	"""Run lint on fixtures and optionally examples/."""
	print("=" * 70)
	print("Cubit Export Lint Self-Test")
	print("=" * 70)
	print()

	# --- Fixtures validation ---
	fixtures_dir = (
		Path(__file__).parent.parent.parent.parent.parent / "tests"
		/ "mcp_server" / "fixtures"
	)
	if not fixtures_dir.exists():
		fixtures_dir = Path(__file__).parent / "fixtures"

	if fixtures_dir.exists():
		bad_file = fixtures_dir / "bad_cubit_script.py"
		clean_file = fixtures_dir / "clean_cubit_script.py"
		if bad_file.exists():
			findings = _lint_file(str(bad_file))
			print(f"  bad_cubit_script.py: {len(findings)} finding(s)")
			if not findings:
				print("  WARNING: bad_cubit_script.py has no findings")
		if clean_file.exists():
			findings = _lint_file(str(clean_file))
			print(f"  clean_cubit_script.py: {len(findings)} finding(s)")
			if findings:
				for f in findings:
					print(f"    L{f['line']} [{f['severity']}] {f['rule']}: {f['message']}")
				print("  FAIL: clean script should have zero findings")
				sys.exit(1)
		print("  fixture validation: PASSED")
		print()

	# --- Examples scan ---
	examples_dir = PROJECT_ROOT / "examples"
	if not examples_dir.exists():
		if not fixtures_dir.exists():
			print(f"SKIP: No fixtures or examples/ found")
		return

	result = lint_cubit_directory("examples")
	print(result)


def main():
	"""Entry point for mcp-server-cubit command."""
	if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
		_selftest()
	else:
		mcp.run(transport="stdio")


if __name__ == '__main__':
	main()
