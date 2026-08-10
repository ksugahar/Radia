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
    python server.py --selftest   # Run lightweight self-test
    python server.py --selftest --audit-repo
                                  # Run self-test plus durable repo-lane audit
"""

from collections import Counter
import hashlib
import json
import errno
import os
import sys
import time
from pathlib import Path

from mcp.server.fastmcp import FastMCP

# Import knowledge bases and rules (relative imports for pip package)
from .rules import ALL_RULES
from .knowledge.export import get_export_documentation
from .knowledge.netgen_workflow import get_netgen_documentation
from .knowledge.scripting import get_cubit_documentation
from .forum_tips import get_forum_tips
from .api_reference import get_api_reference
from .knowledge.panel_conventions import PANEL_CONVENTIONS, LABEL_GUIDE
from .knowledge.custom_toolbar import (
	get_toolbar_documentation,
	generate_toolbar_skeleton,
	generate_dialog_skeleton,
)
from .knowledge.cpp_sdk import get_cpp_sdk_documentation
from .knowledge.mesh_diagnostics import get_diagnostics_documentation
from .knowledge.license import get_license_documentation
from .knowledge.format_routing import get_format_routing_documentation
from .knowledge.coreform_webinars import get_coreform_webinar_documentation
from .vol_inventory import (
	cubit_hex_geometry_refinement_gate,
	cubit_live_mixed_mesh_python_gate as _cubit_live_mixed_mesh_python_gate,
	cubit_mapped_boundary_layer_shell_gate as _cubit_mapped_boundary_layer_shell_gate,
	cubit_sweep_along_curve_gate as _cubit_sweep_along_curve_gate,
	cubit_partitioned_sweep_compatibility_gate as _cubit_partitioned_sweep_compatibility_gate,
	cubit_pyramid_degenerate_hex_export_gate as _cubit_pyramid_degenerate_hex_export_gate,
	cubit_webcut_conformal_hex_gate as _cubit_webcut_conformal_hex_gate,
	cubit_webcut_journal_execution_gate as _cubit_webcut_journal_execution_gate,
	cubit_helical_partition_mesh_gate as _cubit_helical_partition_mesh_gate,
	cubit_source_journal_replay_gate as _cubit_source_journal_replay_gate,
	cubit_boundary_layer_candidate_gate as _cubit_boundary_layer_candidate_gate,
	cubit_boundary_layer_journal_recovery_gate as _cubit_boundary_layer_journal_recovery_gate,
	cubit_embedded_region_mixed_transition_gate as _cubit_embedded_region_mixed_transition_gate,
	cubit_embedded_pipe_source_recovery_gate as _cubit_embedded_pipe_source_recovery_gate,
	cubit_mixed_order_series_inventory_gate,
	summarize_netgen_vol_inventory,
)
from .gmsh_v41 import (
	gmsh_v41_mixed_order_series_gate,
	summarize_gmsh_v41_ascii,
)
from .high_order_export_gate import (
	cubit_headless_netgen_export_gate as _cubit_headless_netgen_export_gate,
	cubit_loft_high_order_vol_series_gate as _cubit_loft_high_order_vol_series_gate,
)
from .helical_conductor_gate import (
	cubit_helical_conductor_source_gate as _cubit_helical_conductor_source_gate,
	cubit_region_owned_mixed_mesh_gate as _cubit_region_owned_mixed_mesh_gate,
)
from .mixed_transition_gate import (
	cubit_conformal_hex_pyramid_tet_interface_gate as _cubit_conformal_hex_pyramid_tet_interface_gate,
	cubit_mixed_transition_source_gate as _cubit_mixed_transition_source_gate,
)
from .cubit_v46_identity import (
	validate_public_identity as _validate_cubit_v46_public_identity,
	validate_source_identity as _validate_cubit_v46_source_identity,
)
from .cross_artifact_mesh_lineage_v47 import (
	validate_public_identity as _validate_cubit_v47_public_identity,
	validate_source_identity as _validate_cubit_v47_source_identity,
)
from .semantic_mesh_identity_v48 import (
	validate_public_identity as _validate_cubit_v48_public_identity,
	validate_source_identity as _validate_cubit_v48_source_identity,
)
from .topology_replay_identity_v49 import (
	validate_public_identity as _validate_cubit_v49_public_identity,
	validate_source_identity as _validate_cubit_v49_source_identity,
)
from .topology_replay_identity_v50 import (
	validate_public_identity as _validate_cubit_v50_public_identity,
	validate_source_identity as _validate_cubit_v50_source_identity,
)
from .stl_inspect import inspect_stl as _inspect_stl
from .quality_parallel_identity_v51 import (
	validate_public_identity as _validate_cubit_v51_public_identity,
	validate_source_identity as _validate_cubit_v51_source_identity,
)
from .structured_hex_gate import (
	cubit_structured_hex_lattice_gate as _cubit_structured_hex_lattice_gate,
	cubit_structured_hex_source_replay_gate as _cubit_structured_hex_source_replay_gate,
)
from .symmetric_mixed_gate import (
	cubit_symmetric_swept_mixed_mesh_gate as _cubit_symmetric_swept_mixed_mesh_gate,
	cubit_symmetric_swept_source_replay_gate as _cubit_symmetric_swept_source_replay_gate,
)
from .pyramid_export_replay_gate import (
	cubit_pyramid_mixed_export_gate as _cubit_pyramid_mixed_export_gate,
	cubit_pyramid_source_plugin_replay_gate as _cubit_pyramid_source_plugin_replay_gate,
)
from .mesh_carrying_sweep_gate import (
	cubit_mesh_carrying_straight_sweep_gate as _cubit_mesh_carrying_straight_sweep_gate,
	cubit_mesh_carrying_straight_sweep_source_replay_gate as _cubit_mesh_carrying_straight_sweep_source_replay_gate,
)
from .power_tools_replay_gate import (
	cubit_partial_volume_hex_diagnosis_gate as _cubit_partial_volume_hex_diagnosis_gate,
	cubit_power_tools_cleanup_source_replay_gate as _cubit_power_tools_cleanup_source_replay_gate,
)
from .ato_sculpt_gate import (
	cubit_ato_levelset_sculpt_source_replay_gate as _cubit_ato_levelset_sculpt_source_replay_gate,
	cubit_levelset_sculpt_hex_validation_gate as _cubit_levelset_sculpt_hex_validation_gate,
)
from ..common import failure_log as _fl, register_status_tool
from ..common import web_docs as _wd
from ..common import examples as _ex

# Server-level instructions delivered to the client model at MCP
# `initialize` time (MathWorks MATLAB MCP server pattern: a short embedded
# operating doctrine, not per-tool docs -- those live in tool descriptions).
_SERVER_INSTRUCTIONS = """\
Coreform Cubit meshing MCP server (Sugahara lab). Capabilities: drive a
persistent Cubit session (GUI or batch), headless mesh dry-runs and
mesh-variant races, Netgen `.vol` export gates, lint, and a curated Cubit
knowledge corpus.

Driving model (lab policy): APREPRO commands + Python on the
HEADLESS/batch route are the PRIMARY way agents drive Cubit --
`cubit_batch_try` / `cubit_mesh_auto` / `.jou` playback for mesh
generation, exports, and validation. The persistent GUI session
(`cubit_show`, `cubit_snapshot`) is the USER's visual-debugging aid:
open it when the user wants to see the model, not as your default
execution surface.

Session model: `cubit_show`/`cubit_exec` reuse ONE persistent Cubit daemon
(first call may take 30+ s for license + startup; later calls are
sub-second). The daemon deliberately survives client restarts. Use
`cubit_session_status` to inspect it; do not spawn ad-hoc Cubit processes.

Operating rules (lab doctrine -- follow these):
1. Probe, don't guess: before writing entity-classification predicates,
   run `cubit_probe(query="entities")` and derive predicates from the
   printed centroids/bboxes/measures. Never hardcode entity IDs.
2. Verify labels: after every block/sideset command run
   `cubit_probe(query="labels")` -- Cubit silently no-ops wrong-kind adds,
   so intended-but-absent membership only shows up there.
3. Gate every solver-bound `.vol`: run `cubit_check_vol` after export and
   before any solver use. A failed gate means stop, not proceed.
4. Checkpoint before risky operations (`cubit_checkpoint` /
   `cubit_restore`).
5. Errors return JSON {"status":"error","kind":...}: kind="input" means
   fix your commands and retry; kind="environment" means a
   license/install problem -- report it to the user instead of retrying;
   kind="internal" means a server bug -- do not retry with new inputs.
Always inspect tool outputs before acting on them; on repeated failure ask
the user rather than looping.

Topology-optimization shape regeneration (the topopt<->CAD loop): a
per-element density design becomes a watertight STL via the radia wheel
(`radia.topopt_cad`: `nodal_from_element_density` -> `iso_stl_from_grid`,
or `write_levelset_exodus` for the Cubit-native `create tri iso` route);
`cubit_stl_to_vol` then meshes that STL into a gated solver `.vol`
(scheme="hex" Sculpt all-hex, scheme="tet" tetmesh).  Sculpt is the tool
for these bodies BECAUSE they have no sweepable topology.  Do NOT route
such STL through netgen's STLGeometry (rejects marching-cubes facets as
degenerate -- recorded negative).  The hex route generates a direct/free
sideset and also binds mesh geometry; the exporter preserves and deduplicates
both boundary representations, recovers adjacent block domains, and splits a
multi-material sideset into domain-aware descriptors.  Gate on closure +
inverted elements +
boundary faces (older exporters could write ZERO `.vol` surface elements,
making a charge-based solve silently demag-free); treat min-quality as a
report, and rank meshes by
`dof_estimate`, not element count.  Prefer scheme="hex" for downstream
charge-Gram / VIM re-evaluation: the facet-conforming tet mesh of a
decimated marching-cubes surface can carry slivers that stall CG at any
tolerance, while the Sculpt hex mesh of the same surface solves in tens
of iterations.

PREFERRED regeneration route -- volume fractions, no STL at all: Sculpt
is natively a volume-fraction mesher.  `radia.topopt_cad
.write_vfrac_exodus` rasterizes the SAME thresholded body into per-cell
volume fractions (Exodus), and `cubit_vfrac_to_vol` runs standalone
`sculpt.exe --input_vfrac` on it: the interface reconstruction places
the boundary sub-cell, so the meshed volume tracks the fraction
integral to ~1e-3 where the STL route's marching-cubes body can lose
percent-level (measured on a real sector design: STL route closure
12.5 % = gate FAILED, vfrac route 0.2 % ok, faster, higher mean
quality).  Both mesh tools accept the vetted Sculpt quality knobs
(`gq_iters`/`gq_threshold` guaranteed-quality smoothing,
`pillow_surfaces`, `defeature`/`min_vol_cells` debris removal).  The
solver half of the combination lives on the radia-ngsolve server:
`hdiv_vim_demag_eval` runs the air-mesh-free charge-Gram demag solve on
the body-only `.vol` (sphere reads ~1/3; the vfrac-route lane locks the
whole chain against that closed form).

MULTI-MATERIAL is the vfrac route's own capability: a partitioned
`write_vfrac_exodus({name: field, ...})` yields `MAT_1..MAT_N`, Sculpt
meshes them with a CONFORMAL interface (shared nodes), and
`cubit_vfrac_to_vol(material_names="core,shell")` labels the `.vol`
regions for per-region `hdiv_vim_demag_eval(mu_r="core=1000,shell=50")`.
A multi-material `.vol` carries interface faces in ADDITION to the outer
skin, so the boundary gate compares the OUTER faces against the
topological skin rather than the raw total.
"""

# Create MCP server
mcp = FastMCP("cubit-export", instructions=_SERVER_INSTRUCTIONS)

# Project root (current working directory for pip-installed package)
PROJECT_ROOT = Path.cwd()

_RULE_REMEDIATIONS = {
	"non-ascii-byte": (
		"ASCII-sanitize Cubit-facing scripts. Keep Japanese/Greek/math symbols "
		"in docs or data files, not in Python that may run under Cubit's cp932 "
		"console."
	),
	"hardcoded-absolute-path": (
		"Replace machine-specific absolute paths with Path(__file__)-relative "
		"paths, parameters, or environment variables."
	),
	"wrong-connectivity-2nd-order": (
		"Use the canonical second-order Netgen/Cubit connectivity mapping; do "
		"not hand-roll node permutations."
	),
	"missing-block-registration": (
		"Register material/volume/sideset blocks before export so downstream "
		"solvers can recover domains and boundaries."
	),
	"missing-cubit-init": (
		"Initialize/import Cubit explicitly in headless scripts before issuing "
		"cubit.cmd calls."
	),
	"missing-boundary-block": (
		"Create explicit boundary blocks/sidesets for solver-facing surfaces."
	),
	"ambiguous-face-block": (
		"Disambiguate face/surface block selection; avoid broad selectors that "
		"can bind the wrong boundary."
	),
}


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


def _lint_directory_summary(directory: str = "examples", top_n: int = 10) -> dict:
	"""Machine-readable directory lint summary for loop/audit bookkeeping."""
	d = Path(directory)
	if not d.is_absolute():
		d = PROJECT_ROOT / d

	if not d.exists():
		return {
			"ok": False,
			"error": f"Directory not found: {d}",
			"directory": str(d),
		}

	py_files = sorted(d.rglob("*.py"))
	limit = max(0, min(int(top_n), 50))
	by_severity = Counter({"CRITICAL": 0, "HIGH": 0, "MODERATE": 0, "LOW": 0})
	by_rule: Counter[str] = Counter()
	top_files = []
	total_findings = 0

	for py_file in py_files:
		findings = _lint_file(str(py_file))
		if not findings:
			continue
		total_findings += len(findings)
		try:
			rel_path = py_file.relative_to(PROJECT_ROOT).as_posix()
		except ValueError:
			rel_path = str(py_file)
		top_files.append({"path": rel_path, "findings": len(findings)})
		for finding in findings:
			by_severity[str(finding.get("severity", "UNKNOWN"))] += 1
			by_rule[str(finding.get("rule", "unknown"))] += 1

	top_rules = [
		{
			"rule": rule,
			"count": count,
			"action": _RULE_REMEDIATIONS.get(
				rule,
				"Inspect representative findings and add a specific remediation note.",
			),
		}
		for rule, count in sorted(by_rule.items(), key=lambda item: (-item[1], item[0]))[:limit]
	]

	return {
		"ok": True,
		"directory": str(d),
		"files_scanned": len(py_files),
		"files_with_findings": len(top_files),
		"total_findings": total_findings,
		"clean": total_findings == 0,
		"by_severity": dict(by_severity),
		"top_rules": top_rules,
		"dominant_rule": top_rules[0] if top_rules else None,
		"top_files": sorted(top_files, key=lambda item: (-item["findings"], item["path"]))[:limit],
	}


# ============================================================
# Netgen code examples
# ============================================================

EXAMPLES = {}

EXAMPLES['cylinder'] = '''# Cylinder with export netgen (high-order curving)
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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['sphere'] = '''# Sphere with export netgen
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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 4/3 * math.pi * R**3
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['torus'] = '''# Torus with export netgen
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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 2 * math.pi**2 * R_MAJOR * R_MINOR**2
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['cone'] = '''# Cone with export netgen
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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
mesh = Mesh(vol_path)

expected_vol = 1/3 * math.pi * R_BASE**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['hex_cylinder'] = '''# Hex-meshed Cylinder with export netgen
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
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
'''

EXAMPLES['complex_named'] = '''# Complex Geometry with Boolean Operations (export netgen)
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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
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
cubit.cmd('export netgen "sphere.vol" order 2 overwrite')

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
cubit.cmd(f'export netgen "{vol_path}" order {ORDER} overwrite')
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
	        "gmsh_v2"               - Gmsh v4.1 format (export gmsh)
	        "gmsh_v4"               - Gmsh v4.1 format policy
	        "netgen"                - Netgen .vol export (export netgen)
	        "nastran"               - Nastran BDF (export jmag_nastran)
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
	        --- Custom GUI extension paths (Python vs C++) ---
	        "cpp_sdk_*"              - C++ SDK plugin path (advanced).
	                                   Subtopics: overview, prerequisites,
	                                   install_and_load, components_example,
	                                   custom_commands, troubleshooting.
	                                   For the Python path see the dedicated
	                                   `cubit_toolbar_guide` tool.
	        --- Cubit license (multi-user lab) ---
	        "license_per_user"       - Each user activates Cubit license in their OWN
	                                   account (renewals cache is bound to the user's
	                                   logon token; admin cannot bootstrap it for them).
	        "license_admin_overwrite"- Same as license_per_user, framed as the
	                                   anti-pattern: do NOT rewrite renewals from admin.
	        --- Mesh format routing ---
	        "format_routing"         - .vol carries labels (FEM input); .step is
	                                   geometry-only (PEEC input); cubit-mesh-export
	                                   is the single chokepoint for producing labeled
	                                   .vol.  STEP-with-extra-tricks workarounds
	                                   (OCC face name regex, sidecar JSON, AP242
	                                   OCAF reader) are uniformly worse.
	        --- Coreform public release/tutorial knowledge ---
	        "coreform"               - Topic index for public Coreform tutorial/release knowledge
	        "coreform_release_2026_6" - Cubit 2026.6 highlights and radia validation actions
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

	# C++ SDK plugin topics: strip prefix and dispatch to cpp_sdk
	# knowledge module.  Browsing aliases included via
	# get_cpp_sdk_documentation's own alias map.
	if topic.startswith("cpp_sdk_"):
		return get_cpp_sdk_documentation(topic[len("cpp_sdk_"):])
	if topic in ("cpp_sdk", "sdk"):
		return get_cpp_sdk_documentation("overview")

	# License topics: each user activates their own license in their
	# own session.  Admin cannot bootstrap it (the renewals cache is
	# bound to the user's logon token).
	if topic in ("license_per_user", "license_per_user_rule"):
		return get_license_documentation("per_user_rule")
	if topic in ("license_admin_overwrite", "license_admin", "license"):
		return get_license_documentation("admin_overwrite")

	# Mesh format routing: .vol carries labels (FEM input);
	# .step is geometry-only (PEEC input); cubit-mesh-export is the
	# single chokepoint for producing labeled .vol.
	if topic in (
		"format_routing", "routing", "vol_vs_step",
		"step_vs_vol", "labels_format", "format",
	):
		return get_format_routing_documentation("format_routing")

	# Coreform webinar / tutorial corpus, synthesized from the official
	# @Coreform YouTube channel (Cubit Tutorials playlist). Paraphrased +
	# attributed; broad technique knowledge beyond the narrow sibling modules.
	if topic in ("coreform", "coreform_webinars", "webinars", "tutorials"):
		return get_coreform_webinar_documentation("index")
	if topic.startswith("coreform_"):
		return get_coreform_webinar_documentation(topic[len("coreform_"):])

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
def cubit_gmsh_v41_inventory(text: str) -> str:
	"""Parse inline ASCII Gmsh 4.1 by entity blocks and validate connectivity."""

	try:
		result = summarize_gmsh_v41_ascii(text)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "gmsh_v41_ascii_inventory_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_gmsh_v41_mixed_order_gate(
	rows: list[dict],
	authoritative_vol_inventory: dict | None = None,
) -> str:
	"""Gate Gmsh 4.1 mixed topology/order while retaining .vol label authority."""

	try:
		result = gmsh_v41_mixed_order_series_gate(
			rows,
			authoritative_vol_inventory=authoritative_vol_inventory,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "gmsh_v41_mixed_order_series_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_vol_inventory(path: str = "", text: str = "") -> str:
	"""
	Return semantic element inventory for a Netgen `.vol` export.

	Use this before routing a Cubit/Coreform `.vol` into a solver path.
	It recognizes triangle/quad surface records and tet/pyramid/wedge/hex
	volume records without modifying them.  This is the preflight for the
	lab split where Cubit owns hex-led or mixed hex+pyramid+tet meshes, while
	Netgen/OCC is enough for tet-only generation.

	Args:
		path: Optional path to a `.vol` file.
		text: Optional inline `.vol` text. If supplied, this wins over path.

	Returns:
		JSON with element counts, material labels, routing hint, and raw
		DomainIn/DomainOut ownership coverage. ``status=needs_attention`` means
		at least one volume domain disappeared from the boundary/interface
		descriptors, a descriptor references an invalid domain, or the same
		surface connectivity was exported twice. Exterior and internal-interface
		counts are reported separately; for a free Sculpt material interface use
		`skin block <id> make sideset <id>` before export.
	"""
	try:
		if text:
			result = summarize_netgen_vol_inventory(text, source="inline")
		elif path:
			p = Path(path)
			result = summarize_netgen_vol_inventory(
				p.read_text(encoding="utf-8"), source=str(p)
			)
		else:
			return json.dumps({
				"status": "error",
				"error": "provide either path or text",
			})
	except Exception as exc:
		return json.dumps({
			"status": "error",
			"error": f"{type(exc).__name__}: {exc}",
		})
	ownership = result.get("boundary_domain_ownership", {})
	status = "ok" if ownership.get("passed") is True else "needs_attention"
	return json.dumps({"status": status, **result}, indent=2)


@mcp.tool()
def cubit_check_vol(vol_path: str,
                    json_path: str = "",
                    contract: str = "",
                    strict_labels: bool = False,
                    threshold_pct: float = 1.0,
                    quality: bool = True,
                    report_json: str = "") -> str:
	"""
	Run the canonical `check-vol` gate on an exported Netgen `.vol` mesh.

	Lab policy requires every solver-bound `.vol` to pass this gate after
	mesh export and BEFORE solver or Simulink initialization.  This tool
	delegates to `cubit_mesh_export.check.check_consistency` (the same
	engine as the `check-vol` CLI): raw DomainIn/DomainOut ownership,
	NGSolve structural load, label inventory/contract, curved-map quality,
	and CAD-sidecar consistency
	(`<vol>.json` auto-discovered when present).

	Runs in the MCP server's Python 3.12 -- no Cubit session or license
	is needed, so it also validates `.vol` files produced by Netgen/OCC.

	Args:
	    vol_path: Path to the `.vol` mesh (absolute, or relative to repo root).
	    json_path: Explicit CAD reference JSON (must exist if given);
	        default auto-discovers `<vol_path>.json`.
	    contract: Optional `radia.vol-label-contract.v1` JSON path.
	    strict_labels: Enforce strict label naming (production modes).
	    threshold_pct: CAD volume/area error threshold in percent (default 1.0).
	    quality: Also run the curved-mapping quality check (default True).
	    report_json: If given, write the `cubit-mesh-export.vol-check.v1`
	        report there (index it from result.json in production runs).

	Returns:
	    JSON: {"passed": bool, "warnings": [...], mesh/labels/materials/
	    boundaries/quality details}.  passed=false means the mesh must not
	    proceed to a solver.
	"""
	try:
		from cubit_mesh_export.check import check_consistency
	except ImportError as exc:
		return json.dumps({
			"status": "error", "stage": "import", "kind": "environment",
			"error": f"cubit-mesh-export (with ngsolve) is required: {exc}",
			"hint": "pip install cubit-mesh-export",
		})
	p = Path(vol_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	kwargs = {}
	if json_path:
		jp = Path(json_path)
		kwargs["json_path"] = str(jp if jp.is_absolute() else PROJECT_ROOT / jp)
	if contract:
		cp = Path(contract)
		kwargs["contract"] = str(cp if cp.is_absolute() else PROJECT_ROOT / cp)
	try:
		result = check_consistency(
			str(p),
			threshold=float(threshold_pct),
			quality=bool(quality),
			strict_labels=bool(strict_labels),
			**kwargs,
		)
	except FileNotFoundError as exc:
		return json.dumps({"status": "error", "stage": "input", "kind": "input",
		                   "kind": "input", "error": str(exc)})
	except Exception as exc:
		return json.dumps({
			"status": "error", "stage": "check", "kind": "input",
			"error": f"{type(exc).__name__}: {exc}",
		})
	if report_json:
		from cubit_mesh_export.check import _write_report_json
		rp = Path(report_json)
		if not rp.is_absolute():
			rp = PROJECT_ROOT / rp
		try:
			_write_report_json(str(rp), result)
			result["report_json"] = str(rp)
		except OSError as exc:
			result["report_json_error"] = str(exc)
	return json.dumps({"status": "ok", **result}, ensure_ascii=False,
	                  indent=2, default=str)


@mcp.tool()
def cubit_mixed_order_series_gate(rows: list[dict]) -> str:
	"""Validate mixed-mesh topology and routing across export orders.

	Pass structured inventory rows produced by ``cubit_vol_inventory``. The
	tool accepts data rather than paths so it does not expand the MCP server's
	local-file read surface.
	"""
	try:
		result = cubit_mixed_order_series_inventory_gate(rows)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_mixed_order_series_inventory_gate",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_live_mixed_mesh_gate(
	summary: dict,
	expected_total_volume: float | None = None,
	volume_relative_tolerance: float = 1.0e-9,
	min_scaled_jacobian: float = 0.0,
) -> str:
	"""Gate a source-journal hex+pyramid+tet replay from headless Cubit Python.

	Pass structured LIVE evidence rather than a local path. The gate binds the
	headless execution mode, element inventory, scaled-Jacobian quality, and CAD
	volume conservation without exposing machine-specific paths.
	"""

	try:
		result = _cubit_live_mixed_mesh_python_gate(
			summary,
			expected_total_volume=expected_total_volume,
			volume_relative_tolerance=volume_relative_tolerance,
			min_scaled_jacobian=min_scaled_jacobian,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_live_mixed_mesh_python_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_sweep_along_curve_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
	volume_relative_tolerance: float = 1.0e-12,
) -> str:
	"""Gate an all-hex mesh-carrying curve sweep and headless launcher evidence."""
	try:
		result = _cubit_sweep_along_curve_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			volume_relative_tolerance=volume_relative_tolerance,
		)
	except (TypeError, ValueError) as exc:
		result = {"policy": "cubit_sweep_along_curve_gate_v1", "status": "invalid_input", "error": str(exc)}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_mesh_carrying_straight_sweep_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
	volume_relative_tolerance: float = 1.0e-12,
) -> str:
	"""Gate a straight include_mesh sweep, topology lattice, and Gmsh export."""
	try:
		result = _cubit_mesh_carrying_straight_sweep_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			volume_relative_tolerance=volume_relative_tolerance,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_mesh_carrying_straight_sweep_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_mesh_carrying_straight_sweep_source_replay_gate(summary: dict) -> str:
	"""Gate official Help provenance, headless replay, and no-mesh control."""
	try:
		result = _cubit_mesh_carrying_straight_sweep_source_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_mesh_carrying_straight_sweep_source_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_partial_volume_hex_diagnosis_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
) -> str:
	"""Gate a truthful partial-volume/low-quality hex rejection."""
	try:
		result = _cubit_partial_volume_hex_diagnosis_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_partial_volume_hex_diagnosis_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_power_tools_cleanup_source_replay_gate(summary: dict) -> str:
	"""Gate an official Power Tools cleanup trace and console diagnosis."""
	try:
		result = _cubit_power_tools_cleanup_source_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_power_tools_cleanup_source_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_levelset_sculpt_hex_validation_gate(
	summary: dict,
	quality_floor: float = 0.2,
	volume_relative_tolerance: float = 0.03,
) -> str:
	"""Gate coarse/fine Sculpt all-hex quality and Gmsh volume closure."""
	try:
		result = _cubit_levelset_sculpt_hex_validation_gate(
			summary,
			quality_floor=quality_floor,
			volume_relative_tolerance=volume_relative_tolerance,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_levelset_sculpt_hex_validation_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_ato_levelset_sculpt_source_replay_gate(summary: dict) -> str:
	"""Gate official ATO provenance, MBG migration, and headless replay."""
	try:
		result = _cubit_ato_levelset_sculpt_source_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_ato_levelset_sculpt_source_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_partitioned_sweep_compatibility_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
	volume_relative_tolerance: float = 1.0e-12,
) -> str:
	"""Gate a legacy webcut/partition journal promoted to an all-hex sweep."""
	try:
		result = _cubit_partitioned_sweep_compatibility_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			volume_relative_tolerance=volume_relative_tolerance,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_partitioned_sweep_compatibility_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_mapped_boundary_layer_shell_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
	min_shape: float = 0.0,
	shell_radius_tolerance: float = 1.0e-9,
	volume_relative_tolerance: float = 1.0e-12,
) -> str:
	"""Gate mapped all-hex boundary layers by nodal shells, quality, and scale."""
	try:
		result = _cubit_mapped_boundary_layer_shell_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			min_shape=min_shape,
			shell_radius_tolerance=shell_radius_tolerance,
			volume_relative_tolerance=volume_relative_tolerance,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_mapped_boundary_layer_shell_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_pyramid_degenerate_hex_export_gate(
	summary: dict,
	min_hex_scaled_jacobian: float = 0.2,
	min_pyramid_geometric_volume: float = 0.0,
) -> str:
	"""Gate CPYRAM versus nopyramid decks, including order-2 linearization."""
	try:
		result = _cubit_pyramid_degenerate_hex_export_gate(
			summary,
			min_hex_scaled_jacobian=min_hex_scaled_jacobian,
			min_pyramid_geometric_volume=min_pyramid_geometric_volume,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_pyramid_degenerate_hex_export_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_webcut_conformal_hex_gate(summary: dict) -> str:
	"""Gate webcut volume drift, partition balance, interfaces, and hex quality."""
	try: result = _cubit_webcut_conformal_hex_gate(summary)
	except (TypeError, ValueError) as exc: result = {"policy":"cubit_webcut_conformal_hex_gate_v1","status":"invalid_input","error":str(exc)}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_webcut_journal_execution_gate(summary: dict) -> str:
	"""Gate source-journal operation order and headless process evidence."""
	try: result = _cubit_webcut_journal_execution_gate(summary)
	except (TypeError, ValueError) as exc: result = {"policy":"cubit_webcut_journal_execution_gate_v1","status":"invalid_input","error":str(exc)}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_helical_partition_mesh_gate(summary: dict) -> str:
	"""Gate a many-body helical mesh against quality and parsed .vol inventory."""
	try: result = _cubit_helical_partition_mesh_gate(summary)
	except (TypeError, ValueError) as exc: result = {"policy":"cubit_helical_partition_mesh_gate_v1","status":"invalid_input","error":str(exc)}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_source_journal_replay_gate(summary: dict) -> str:
	"""Gate synchronous, headless replay and expected mesh disposition."""
	try: result = _cubit_source_journal_replay_gate(summary)
	except (TypeError, ValueError) as exc: result = {"policy":"cubit_source_journal_replay_gate_v1","status":"invalid_input","error":str(exc)}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_boundary_layer_candidate_gate(
	candidates: list[dict],
	min_scaled_jacobian: float = 0.2,
	max_volume_relative_error: float = 1.0e-5,
) -> str:
	"""Select a non-inverted boundary-layer sweep candidate with export closure."""
	try:
		result = _cubit_boundary_layer_candidate_gate(
			candidates,
			min_scaled_jacobian=min_scaled_jacobian,
			max_volume_relative_error=max_volume_relative_error,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_boundary_layer_candidate_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_boundary_layer_journal_recovery_gate(summary: dict) -> str:
	"""Gate three-parameter, pairwise, headless recovery of a failed journal."""
	try:
		result = _cubit_boundary_layer_journal_recovery_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_boundary_layer_journal_recovery_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_embedded_region_mixed_transition_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.1,
	max_volume_relative_error: float = 1.0e-5,
	max_swept_volume_relative_error: float = 5.0e-4,
) -> str:
	"""Gate hex-led tet/pyramid recovery, quality, interfaces, and Gmsh 4.1."""
	try:
		result = _cubit_embedded_region_mixed_transition_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			max_volume_relative_error=max_volume_relative_error,
			max_swept_volume_relative_error=max_swept_volume_relative_error,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_embedded_region_mixed_transition_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_embedded_pipe_source_recovery_gate(
	summary: dict,
	expected_unmeshed_volumes: list[int] | None = None,
) -> str:
	"""Gate source-journal replay and semantically classified version recovery."""
	try:
		kwargs = {}
		if expected_unmeshed_volumes is not None:
			kwargs["expected_unmeshed_volumes"] = expected_unmeshed_volumes
		result = _cubit_embedded_pipe_source_recovery_gate(summary, **kwargs)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_embedded_pipe_source_recovery_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_region_owned_mixed_mesh_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.1,
	max_volume_relative_error: float = 5.0e-8,
) -> str:
	"""Gate region-owned conductor hex and air tet/pyramid topology."""
	try:
		result = _cubit_region_owned_mixed_mesh_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			max_volume_relative_error=max_volume_relative_error,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_region_owned_mixed_mesh_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_helical_conductor_source_gate(summary: dict) -> str:
	"""Gate a helical-conductor source replay and classified tet fallback."""
	try:
		result = _cubit_helical_conductor_source_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_helical_conductor_source_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_symmetric_swept_mixed_mesh_gate(summary: dict) -> str:
	"""Gate symmetric CAD, hex/pyramid/tet topology, quality, and Gmsh closure."""
	try:
		result = _cubit_symmetric_swept_mixed_mesh_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_symmetric_swept_mixed_mesh_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_symmetric_swept_source_replay_gate(summary: dict) -> str:
	"""Gate source-journal headless replay and public mixed-mesh closure."""
	try:
		result = _cubit_symmetric_swept_source_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_symmetric_swept_source_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_conformal_hex_pyramid_tet_interface_gate(
	summary: dict,
	mapped_volume_id: int = 1,
	transition_volume_id: int = 2,
	min_scaled_jacobian: float = 0.1,
	max_volume_relative_error: float = 1.0e-9,
) -> str:
	"""Gate a conformal hex-pyramid-tet interface and independent volume sum."""
	try:
		result = _cubit_conformal_hex_pyramid_tet_interface_gate(
			summary,
			mapped_volume_id=mapped_volume_id,
			transition_volume_id=transition_volume_id,
			min_scaled_jacobian=min_scaled_jacobian,
			max_volume_relative_error=max_volume_relative_error,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_conformal_hex_pyramid_tet_interface_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	v46_checks = _validate_cubit_v46_public_identity(summary)
	if v46_checks:
		result.setdefault("checks", {}).update(v46_checks["checks"])
		result["cubit_v46_public_identity"] = v46_checks
		if v46_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v47_checks = _validate_cubit_v47_public_identity(summary)
	if v47_checks:
		result.setdefault("checks", {}).update(v47_checks["checks"])
		result["cubit_v47_public_identity"] = v47_checks
		if v47_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v48_checks = _validate_cubit_v48_public_identity(summary)
	if v48_checks:
		result.setdefault("checks", {}).update(v48_checks["checks"])
		result["cubit_v48_public_identity"] = v48_checks
		if v48_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v49_checks = _validate_cubit_v49_public_identity(summary)
	if v49_checks:
		result.setdefault("checks", {}).update(v49_checks["checks"])
		result["cubit_v49_public_identity"] = v49_checks
		if v49_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v50_checks = _validate_cubit_v50_public_identity(summary)
	if v50_checks:
		result.setdefault("checks", {}).update(v50_checks["checks"])
		result["cubit_v50_public_identity"] = v50_checks
		if v50_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v51_checks = _validate_cubit_v51_public_identity(summary)
	if v51_checks:
		result.setdefault("checks", {}).update(v51_checks["checks"])
		result["cubit_v51_public_identity"] = v51_checks
		if v51_checks["status"] != "ok":
			result["status"] = "needs_attention"
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_mixed_transition_source_gate(
	summary: dict,
	mapped_volume_id: int = 1,
	transition_volume_id: int = 2,
) -> str:
	"""Gate source commands, headless diagnostics, and quality API recovery."""
	try:
		result = _cubit_mixed_transition_source_gate(
			summary,
			mapped_volume_id=mapped_volume_id,
			transition_volume_id=transition_volume_id,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_mixed_transition_source_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	v46_checks = _validate_cubit_v46_source_identity(summary)
	if v46_checks:
		result.setdefault("checks", {}).update(v46_checks["checks"])
		result["cubit_v46_source_identity"] = v46_checks
		if v46_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v47_checks = _validate_cubit_v47_source_identity(summary)
	if v47_checks:
		result.setdefault("checks", {}).update(v47_checks["checks"])
		result["cubit_v47_source_identity"] = v47_checks
		if v47_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v48_checks = _validate_cubit_v48_source_identity(summary)
	if v48_checks:
		result.setdefault("checks", {}).update(v48_checks["checks"])
		result["cubit_v48_source_identity"] = v48_checks
		if v48_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v49_checks = _validate_cubit_v49_source_identity(summary)
	if v49_checks:
		result.setdefault("checks", {}).update(v49_checks["checks"])
		result["cubit_v49_source_identity"] = v49_checks
		if v49_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v50_checks = _validate_cubit_v50_source_identity(summary)
	if v50_checks:
		result.setdefault("checks", {}).update(v50_checks["checks"])
		result["cubit_v50_source_identity"] = v50_checks
		if v50_checks["status"] != "ok":
			result["status"] = "needs_attention"
	v51_checks = _validate_cubit_v51_source_identity(summary)
	if v51_checks:
		result.setdefault("checks", {}).update(v51_checks["checks"])
		result["cubit_v51_source_identity"] = v51_checks
		if v51_checks["status"] != "ok":
			result["status"] = "needs_attention"
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_pyramid_mixed_export_gate(summary: dict) -> str:
	"""Gate explicit hex/pyramid/tet preservation in Gmsh and Nastran."""
	try:
		result = _cubit_pyramid_mixed_export_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_pyramid_mixed_export_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_pyramid_source_plugin_replay_gate(summary: dict) -> str:
	"""Gate legacy source migration through executable-owned plugin startup."""
	try:
		result = _cubit_pyramid_source_plugin_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_pyramid_source_plugin_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_structured_hex_lattice_gate(
	summary: dict,
	min_scaled_jacobian: float = 0.2,
	max_volume_relative_error: float = 1.0e-10,
) -> str:
	"""Gate structured all-hex counts, quality, and Gmsh volume closure."""
	try:
		result = _cubit_structured_hex_lattice_gate(
			summary,
			min_scaled_jacobian=min_scaled_jacobian,
			max_volume_relative_error=max_volume_relative_error,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_structured_hex_lattice_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_structured_hex_source_replay_gate(summary: dict) -> str:
	"""Gate source commands, license diagnostics, and headless exit semantics."""
	try:
		result = _cubit_structured_hex_source_replay_gate(summary)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_structured_hex_source_replay_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_hex_refinement_geometry_gate(
	rows: list[dict],
	min_quality: float = 0.2,
	max_final_volume_relative_error: float = 0.05,
	min_error_reduction_fraction: float = 0.25,
) -> str:
	"""Detect curved-geometry error plateaus in an all-hex refinement series.

	Pass archived headless-Cubit rows rather than a local file path.  Increasing
	hex count is not accepted as convergence unless the independent mesh-volume
	error also improves while scaled-Jacobian quality remains acceptable.
	"""
	try:
		result = cubit_hex_geometry_refinement_gate(
			rows,
			min_quality=min_quality,
			max_final_volume_relative_error=max_final_volume_relative_error,
			min_error_reduction_fraction=min_error_reduction_fraction,
		)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_hex_geometry_refinement_plateau_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_loft_high_order_vol_series_gate(
	rows: list[dict],
	min_quality: float = 0.2,
) -> str:
	"""Gate all-hex loft topology, curved payload, sidecars, and quality for orders 1-5."""

	try:
		result = _cubit_loft_high_order_vol_series_gate(rows, min_quality=min_quality)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_loft_high_order_vol_series_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_headless_netgen_export_gate(
	summary: dict,
	min_quality: float = 0.2,
) -> str:
	"""Gate migration from a GUI plugin export command to native headless Netgen export."""

	try:
		result = _cubit_headless_netgen_export_gate(summary, min_quality=min_quality)
	except (TypeError, ValueError) as exc:
		result = {
			"policy": "cubit_headless_netgen_export_command_gate_v1",
			"status": "invalid_input",
			"error": str(exc),
		}
	return json.dumps(result, ensure_ascii=False, indent=2)


def _normalized_cubit_journal_commands(text: str) -> list[str]:
	commands = []
	for raw_line in str(text or "").splitlines():
		line = raw_line.strip()
		if not line or line.startswith("#"):
			continue
		commands.append(" ".join(line.split()))
	return commands


@mcp.tool()
def cubit_journal_reproducibility_gate(
	journal_a: str,
	journal_b: str,
	outcome_a: str = "",
	outcome_b: str = "",
) -> str:
	"""Compare two Cubit journals without inventing a script root cause.

	Comments and blank lines are removed, then the remaining command streams are
	compared exactly.  If outcomes differ while commands are identical, the tool
	requires runtime provenance instead of claiming that the journal explains the
	difference.  Inputs are text, not paths, so this tool does not expand the MCP
	server's local-file read surface.
	"""

	commands_a = _normalized_cubit_journal_commands(journal_a)
	commands_b = _normalized_cubit_journal_commands(journal_b)
	digest_a = hashlib.sha256("\n".join(commands_a).encode("utf-8")).hexdigest()
	digest_b = hashlib.sha256("\n".join(commands_b).encode("utf-8")).hexdigest()
	commands_equal = commands_a == commands_b
	outcomes_recorded = bool(str(outcome_a).strip()) and bool(str(outcome_b).strip())
	outcomes_differ = outcomes_recorded and str(outcome_a).strip().lower() != str(outcome_b).strip().lower()
	differences = []
	for index in range(max(len(commands_a), len(commands_b))):
		left = commands_a[index] if index < len(commands_a) else None
		right = commands_b[index] if index < len(commands_b) else None
		if left != right:
			differences.append({"command_index": index + 1, "journal_a": left, "journal_b": right})
		if len(differences) >= 20:
			break

	if commands_equal and outcomes_differ:
		status = "needs_run_provenance"
	elif commands_equal:
		status = "commands_equivalent"
	else:
		status = "commands_differ"
	result = {
		"policy": "cubit_journal_reproducibility_gate",
		"status": status,
		"analysis_complete": True,
		"commands_equal": commands_equal,
		"command_count_a": len(commands_a),
		"command_count_b": len(commands_b),
		"command_digest_a": digest_a,
		"command_digest_b": digest_b,
		"outcome_a": str(outcome_a).strip() or None,
		"outcome_b": str(outcome_b).strip() or None,
		"outcomes_recorded": outcomes_recorded,
		"outcomes_differ": outcomes_differ,
		"script_difference_explains_outcome": bool(not commands_equal and outcomes_differ),
		"differences": differences,
		"required_run_provenance": [
			"solver_version_and_build",
			"headless_command_line",
			"process_exit_code",
			"complete_solver_log",
			"initial_session_state",
			"mesh_artifact_digest",
			"mesh_quality_summary",
		],
		"notes": [
			"Identical command streams cannot support a journal-source root-cause claim.",
			"A comment label such as NG or OK is not solver evidence.",
			"Replay through a clean headless process and bind the runtime evidence to the exported mesh digest.",
		],
	}
	return json.dumps(result, ensure_ascii=False, indent=2)


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
	        "export netgen"     - export netgen() API reference
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
	        "cylinder"       - Cylinder with export netgen
	        "sphere"         - Sphere with export netgen
	        "torus"          - Torus with export netgen
	        "cone"           - Cone with export netgen
	        "hex_cylinder"   - Hex-meshed cylinder with export netgen
	        "complex_named"  - Boolean geometry with export netgen
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
	- mesh.Curve() without export netgen (MODERATE)
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


@mcp.tool()
def cubit_audit_summary(directory: str = "examples", top_n: int = 10) -> dict:
	"""
	Return a machine-readable Cubit export-lint audit summary.

	Use this for learning-loop bookkeeping and dashboards. It reports files
	scanned, total findings, severity counts, top rule counts, and top files
	without printing the long per-file audit body.
	"""
	return _lint_directory_summary(directory, top_n)


# ============================================================
# New tools: script generation and mesh quality
# ============================================================

SCRIPT_TEMPLATES = {}

SCRIPT_TEMPLATES['tet_netgen'] = '''# Cubit -> NGSolve: Tet mesh with high-order curving (export netgen)
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
cubit.cmd(f'export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
mesh = Mesh(vol_path)

# === VERIFY ===
vol = Integrate(CF(1), mesh)
area = Integrate(CF(1), mesh, VOL_or_BND=BND)
print(f"Volume: {vol:.6f}, Surface area: {area:.6f}")
'''

SCRIPT_TEMPLATES['tet_netgen_named'] = '''# Cubit -> NGSolve: Complex geometry with Boolean operations (export netgen)
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
cubit.cmd(f'export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
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
cubit.cmd('export netgen "mesh.vol" order 2 overwrite')

# === IMPORT TO NGSOLVE ===
mesh = Mesh("mesh.vol")
vol = Integrate(CF(1), mesh)
print(f"Volume: {vol:.6f}")
'''

SCRIPT_TEMPLATES['hex_netgen'] = '''# Cubit -> NGSolve: Hex mesh with high-order curving (export netgen)
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
cubit.cmd(f'export netgen "{vol_path}" order {CURVE_ORDER} overwrite')
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
	        "tet_netgen"       - Tet mesh -> export netgen with ACIS curving
	        "tet_netgen_named" - Complex geometry with export netgen
	        "tet_gmsh"         - Tet mesh -> Netgen .vol 2nd order (simplest)
	        "hex_netgen"       - Hex mesh -> export netgen with curving
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
			'trigger': 'export gmsh without any block/mesh',
			'fix': 'cubit.cmd("block 1 add tet all")',
		},
		{
			'rule': 'missing-mesh-command',
			'severity': 'CRITICAL',
			'description': 'Export function called but no mesh command found.',
			'trigger': 'export without mesh volume/surface command',
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
			'fix': 'Use export netgen "f.vol" order N',
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
			'fix': 'Use export netgen "f.vol" order N',
		},
		{
			'rule': 'deleted-api-usage',
			'severity': 'CRITICAL',
			'description': '(See above) Catches export_netgen and other deleted APIs.',
			'trigger': 'export_netgen(cubit, geometry=geo)',
			'fix': 'Use export netgen "f.vol" order N',
		},
		{
			'rule': 'deleted-api-noheal',
			'severity': 'INFO',
			'description': 'Name-based workflow (name_occ_faces + noheal) is deleted. Use export netgen.',
			'trigger': 'name_occ_faces(shape) + STEP import',
			'fix': 'Use export netgen(cubit, order=N) - no STEP needed',
		},
		{
			'rule': 'hardcoded-absolute-path',
			'severity': 'MODERATE',
			'description': 'Hardcoded absolute paths in sys.path (except Coreform Cubit/NGSolve).',
			'trigger': 'sys.path.insert(0, "public-safe curated corpus")',
			'fix': 'sys.path.insert(0, os.path.dirname(__file__))',
		},
		{
			'rule': 'missing-boundary-block',
			'severity': 'MODERATE',
			'description': 'Volume element block found but no surface element block for Netgen export.',
			'trigger': 'block 1 add tet all + export netgen(...) without tri block',
			'fix': 'cubit.cmd("block 2 add tri all")',
		},
		{
			'rule': 'ambiguous-face-block',
			'severity': 'LOW',
			'description': 'block add face is ambiguous: APREPRO face selects generic surface elements, but Python APIs split tri/quad blocks.',
			'trigger': 'cubit.cmd("block 2 add face all in surface all")',
			'fix': 'Use tri for tet boundaries, quad for hex boundaries, or add an explicit mixed tri/quad comment.',
		},
		{
			'rule': 'wrong-file-extension',
			'severity': 'MODERATE',
			'description': 'Export file extension does not match the format.',
			'trigger': 'export gmsh "mesh.vtk" (wrong extension)',
			'fix': 'Use correct extension: .msh, .bdf, .exo, .vol',
		},
		{
			'rule': 'curve-without-export-curved',
			'severity': 'MODERATE',
			'description': 'Manual mesh.Curve() is not needed. export netgen embeds curving in .vol.',
			'trigger': 'mesh.Curve(3) without export netgen',
			'fix': 'Use cubit.cmd(\'export netgen "mesh.vol" order 3 overwrite\') then mesh = Mesh("mesh.vol")',
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
			'fix': 'Add the missing class to the PySide6 import statement.',
		},
		{
			'rule': 'pyqt5-import-forbidden',
			'severity': 'HIGH',
			'description': 'PyQt5 import in Cubit UI code. Radia targets Coreform Cubit 2025.12+ and is PySide6-only.',
			'trigger': 'from PyQt5.QtWidgets import ...',
			'fix': 'Use PySide6 only; do not keep a PyQt5 fallback.',
		},
	]

	lines = [f"# Cubit Export Lint Rules ({len(rules_info)} rules)", ""]
	for r in rules_info:
		lines.append(f"## [{r['severity']}] {r['rule']}")
		lines.append(f"{r['description']}")
		lines.append(f"- **Trigger**: {r['trigger']}")
		lines.append(f"- **Fix**: `{r['fix']}`")
		lines.append("")

	return '\n'.join(lines)


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
#     (C:/Program Files/Coreform Cubit 2025.12/bin/python3/) to drive a
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
	  .py   → `play "<path>"`           (Cubit's `play` runs Python files)
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
	from .session import find_cubit_install, _cubit_gui_exe
	from .license_warmup import warmup_license

	if cubit_exe:
		exe = cubit_exe
		bin_dir = Path(exe).parent
	else:
		bin_dir = find_cubit_install()
		if bin_dir is None:
			return _json.dumps({"status": "error", "stage": "cubit_binary", "kind": "environment",
			                    "error": "Cubit install not found on PATH "
			                             "or standard locations; pass cubit_exe."})
		try:
			exe = str(_cubit_gui_exe(bin_dir))
		except FileNotFoundError as e:
			return _json.dumps({"status": "error", "stage": "cubit_binary", "kind": "environment",
			                    "error": str(e)})
	if not Path(exe).exists():
		return _json.dumps({"status": "error", "stage": "cubit_binary", "kind": "environment",
		                    "error": f"Cubit exe not found: {exe}"})

	# Pre-warm Learn license (same logic as cubit_session):
	# rlm_activate --login only if cache stale, otherwise skip.
	# This shortens the visible Cubit launch from 30 – 60 s (cold
	# in-process checkout) to ~3 s when the cache is fresh.
	_license_warmup = warmup_license(Path(bin_dir), timeout_s=30.0)

	work = Path(_tf.mkdtemp(prefix="cubit_gui_"))
	wrapper = work / "wrapper.jou"
	lines = []

	if path:
		p = Path(path)
		if not p.is_absolute():
			p = PROJECT_ROOT / p
		if not p.exists():
			return _json.dumps({"status": "error", "stage": "input", "kind": "input",
			                    "error": f"File not found: {p}"})
		abs_path = str(p.resolve()).replace("\\", "/")
		suffix = p.suffix.lower()
		if suffix == ".jou":
			lines.append(f'playback "{abs_path}"')
		elif suffix == ".py":
			# Cubit の `play` は .jou と .py の両方を受け付ける.
			lines.append(f'play "{abs_path}"')
		elif suffix in (".step", ".stp"):
			lines.append(f'import step "{abs_path}"')
		elif suffix in (".brep", ".brp"):
			# Cubit imports BREP via the ACIS kernel.
			lines.append(f'import acis "{abs_path}"')
		elif suffix in (".msh", ".vol", ".g", ".e", ".exo"):
			lines.append(f'import mesh "{abs_path}"')
		else:
			return _json.dumps({"status": "error", "stage": "input", "kind": "input",
			                    "error": f"Unsupported extension: {suffix}"})

	if commands:
		lines.extend(str(c).rstrip() for c in commands)

	if not lines:
		return _json.dumps({"status": "error", "stage": "input", "kind": "input",
		                    "error": "Either `path` or `commands` is required."})

	wrapper.write_text("\n".join(lines) + "\n")

	# GUI mode: no -batch, no -nographics. -nojournal prevents auto-
	# journaling of the wrapper playback (keeps the directory clean).
	# -input <wrapper> plays our commands after GUI initialization.
	cmd = [exe, "-nojournal", "-input", str(wrapper)]

	try:
		if detach:
			# stdin=DEVNULL: prevent Cubit from inheriting MCP server's
			# stdio pipe (used for JSON-RPC with Claude Code). Without
			# this, Cubit's startup hangs (Responding=False, threads
			# decline) — observed 2026-04-25 LAB.
			proc = _sp.Popen(
				cmd,
				stdin=_sp.DEVNULL,
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
		return _json.dumps({"status": "error", "stage": "timeout", "kind": "environment",
		                    "wrapper": str(wrapper)})
	except OSError as e:
		return _json.dumps({"status": "error", "stage": "launch", "kind": "environment",
		                    "error": str(e)})

	return _json.dumps(info, indent=2)


# ============================================================
# Cubit Viewer Session tools (OCP-inspired persistent daemon)
# ============================================================
#
# Design note (2026-04-19): these tools drive a long-lived Cubit GUI
# via a JSON-RPC daemon running under Cubit's bundled Python 3.10.
# Launch cost is paid once (~0.7s init); subsequent commands are
# <1 ms. The daemon is auto-discovered via `find_cubit_install()`
# (env > registry > glob), so there is no hardcoded version path.
#
# The design follows OCP CAD Viewer's architecture (stateless client
# calls + persistent viewer) transposed to Cubit: Cubit itself plays
# the role of OCP's WebView, and this mcp-server plays the role of
# OCP's Python client.

from . import session as _cs

# Error `kind` classification for the LLM audience (MathWorks
# UnexpectedErrorPrefixForLLM pattern, extended):
#   "input"       -- the caller's commands/paths are wrong; fix + retry.
#   "environment" -- license / install / hung-session problem; report to
#                    the user instead of retrying with new inputs.
#   "internal"    -- server/daemon bug; do not retry.
_ENVIRONMENT_ERROR_NEEDLES = (
    "license", "rlm", "could not locate coreform cubit", "did not signal",
    "bootstrap", "exited during", "response timeout", "daemon died",
)


def _session_log_pointer() -> str | None:
    """Where the full session diagnostics live (MathWorks 'for details,
    see the server log in <path>' pattern)."""
    try:
        drop = _cs._user_daemon_dir()
    except Exception:
        return None
    return str(drop)


def _error_payload(stage: str, message: str, *, kind: str | None = None,
                   hint: str | None = None) -> dict:
    """Cubit-flavored wrapper over the shared error contract."""
    from ..common.server_hardening import error_payload
    log = _session_log_pointer()
    return error_payload(
        stage, message, kind=kind, hint=hint,
        environment_needles=_ENVIRONMENT_ERROR_NEEDLES,
        log=(f"{log} (bootstrap.log / cubit_stderr.log / "
             "startup_error.txt hold the full record)") if log else None,
    )


def _cubit_session_or_error():
    """Lazy-init the singleton; return (session, None) or (None, err_json)."""
    try:
        sess = _cs.CubitSession.get()
        sess.ensure_started()
        return sess, None
    except _cs.CubitSessionError as e:
        err = _error_payload(
            "session_init", str(e),
            hint=("Ensure Coreform Cubit is installed. "
                  "Override path with `set_cubit_bin_dir(path)` or "
                  "`CUBIT_BIN_DIR` environment variable."),
        )
        return None, json.dumps(err, indent=2)


def _cubit_path_dispatch(path: Path) -> str:
    """Map a file extension to the matching Cubit import command."""
    suffix = path.suffix.lower()
    abs_path = str(path.resolve()).replace("\\", "/")
    if suffix == ".jou":
        return f'playback "{abs_path}"'
    if suffix == ".py":
        # Cubit の `play` は .jou と .py の両方を受け付ける.
        # See: cubit_scripting_knowledge topic="aprepro_journal".
        return f'play "{abs_path}"'
    if suffix in (".step", ".stp"):
        return f'import step "{abs_path}"'
    if suffix in (".brep", ".brp"):
        return f'import acis "{abs_path}"'
    if suffix in (".msh", ".vol", ".g", ".e", ".exo"):
        return f'import mesh "{abs_path}"'
    raise ValueError(f"Unsupported extension for Cubit viewer: {suffix}")


@mcp.tool()
def cubit_show(path: str = "", extra_commands: list = None) -> str:
	"""
	Load a file into the **persistent Cubit viewer** and optionally run
	follow-up commands. First call launches the daemon + Cubit GUI;
	subsequent calls reuse the session (sub-second command cycles).

	Supported inputs: .jou (playback), .py (play — Cubit-side Python),
	.step/.stp (import step), .brep (import acis),
	.msh/.vol/.g/.e/.exo (import mesh).

	Args:
	    path: file to load. Leave empty to run only `extra_commands`.
	    extra_commands: extra Cubit commands to run after loading.

	Returns JSON with the per-command results (line + ok + rc).
	"""
	from pathlib import Path as _P
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	cmds = []
	if path:
		p = _P(path)
		if not p.is_absolute():
			p = PROJECT_ROOT / p
		if not p.exists():
			return json.dumps({"status": "error", "stage": "input", "kind": "input",
			                   "error": f"File not found: {p}"})
		try:
			cmds.append(_cubit_path_dispatch(p))
		except ValueError as e:
			return json.dumps({"status": "error", "stage": "input", "kind": "input",
			                   "error": str(e)})
	if extra_commands:
		cmds.extend(str(c).rstrip() for c in extra_commands)
	if not cmds:
		return json.dumps({"status": "error", "stage": "input", "kind": "input",
		                   "error": "Either `path` or `extra_commands` is required."})

	try:
		r = sess.call("cmd", cmds)
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	return json.dumps(r, indent=2)


def _probe_summary_safe(sess) -> dict:
	"""Best-effort `probe summary`. Swallows RPC errors."""
	try:
		r = sess.call("probe", ["summary"], timeout_s=15)
	except Exception:
		return {}
	res = r.get("result") if isinstance(r, dict) else None
	return res if isinstance(res, dict) else {}


def _state_delta(before: dict, after: dict) -> dict:
	"""Per-key delta (after - before). Zero / missing keys omitted."""
	keys = set(before) | set(after)
	delta: dict = {}
	for k in keys:
		try:
			a = int(after.get(k, 0) or 0)
			b = int(before.get(k, 0) or 0)
			if a - b != 0:
				delta[k] = a - b
		except (TypeError, ValueError):
			continue
	return delta


@mcp.tool()
def cubit_exec(commands: list) -> str:
	"""
	Send arbitrary Cubit commands to the persistent viewer session.

	Use this for incremental modelling / meshing driven by the LLM:
	each call appends commands to the running Cubit session, which
	the student can watch live in Cubit's GUI window.

	Args:
	    commands: list of Cubit command strings (one per list item, no
	        trailing newline). Runs in order; stops at first failure.

	Returns JSON with:
	  - `result`: per-command {line, ok, rc|error}
	  - `state_before` / `state_after`: volume/surface/node/hex/tet counts
	  - `state_delta`: what changed (e.g., {"hexes": +4740, "nodes": +6210})
	  - `all_ok`: true iff every command reported ok
	  - `failed_at`: 1-based index of first failure (if any)
	  - `hint`: short advice on the failure (from heuristics + failure log)

	On a Cubit error, `ok=false` on the offending line and subsequent
	commands are NOT executed so the viewer's state stays inspectable.
	Failures are also appended to the persistent failure log so
	`cubit_lookup` and `cubit_recent_failures` can surface them later.
	"""
	if not commands:
		return json.dumps({"status": "error",
		                   "error": "commands list cannot be empty"})
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	cmd_list = [str(c).rstrip() for c in commands]

	# Pre-flight: scan recent failure log for a similar input
	preflight = _cubit_preflight_warn(cmd_list)

	before = _probe_summary_safe(sess)
	try:
		r = sess.call("cmd", cmd_list)
	except _cs.CubitSessionError as e:
		_fl.record_failure("cubit", {
			"input": cmd_list,
			"error": str(e),
			"stage": "rpc",
			"context": {"state_before": before},
		})
		return json.dumps(_error_payload(
			"rpc", str(e),
			hint="Cubit daemon RPC failed. An owned session auto-restarts "
			     "on the next call (state lost -- consider cubit_checkpoint "
			     "before risky ops); an attached session is left running -- "
			     "use cubit_session_shutdown to force-stop a hung daemon."))
	after = _probe_summary_safe(sess)

	per_line = r.get("result") if isinstance(r, dict) else r
	if not isinstance(per_line, list):
		per_line = []
	all_ok = all(step.get("ok") for step in per_line) if per_line else False
	failed_at = None
	first_err = None
	for idx, step in enumerate(per_line, start=1):
		if not step.get("ok"):
			failed_at = idx
			first_err = step.get("error") or f"rc={step.get('rc')}"
			break

	payload: dict = {
		"status": "ok" if all_ok else "error",
		"all_ok": all_ok,
		"executed": len(per_line),
		"submitted": len(cmd_list),
		"result": per_line,
		"state_before": before,
		"state_after": after,
		"state_delta": _state_delta(before, after),
	}
	if preflight:
		payload["preflight"] = preflight
	if failed_at is not None:
		payload["failed_at"] = failed_at
		payload["failed_cmd"] = cmd_list[failed_at - 1] if failed_at <= len(cmd_list) else None
		payload["error"] = (first_err or "").strip()[:600]
		payload["hint"] = _cubit_failure_hint(payload["failed_cmd"] or "",
		                                       payload["error"])
		_fl.record_failure("cubit", {
			"input": cmd_list,
			"failed_cmd": payload["failed_cmd"],
			"failed_at": failed_at,
			"error": payload["error"],
			"context": {"state_before": before, "state_after": after},
		})
	if isinstance(r, dict) and r.get("_recovered"):
		payload["_recovered"] = True
		payload["hint"] = ((payload.get("hint", "") or "")
		                   + " NOTE: Cubit session was auto-restarted; "
		                     "geometry/mesh state is lost.").strip()
	return json.dumps(payload, indent=2)


_CUBIT_ERROR_HEURISTICS: list[tuple[str, str]] = [
	("unknown command",  "Check Cubit reference via `cubit_docs(topic='syntax')` "
	                     "or `cubit_web_docs(source='cubit', query='<keyword>')`."),
	("parse",            "Command syntax issue. Try `cubit_lookup(query='<verb>')` for examples."),
	("no such",          "Entity (volume/surface/vertex) does not exist. Run `cubit_probe('summary')` first."),
	("cannot mesh",      "Try `cubit_mesh_diagnose()` to propose scheme alternatives per volume."),
	("scheme",           "Meshing scheme problem. `cubit_mesh_diagnose()` or set `volume all scheme tetmesh`."),
	("sweep",            "Sweep failed. Geometry is likely non-prismatic. Fall back to `scheme tetmesh`."),
	("import",           "Import failed. Verify the path exists and format matches "
	                     "(`import step`, `import acis`, `import mesh`)."),
	("export",           "Export failed. Check `cubit_docs(topic='export')` and that mesh exists."),
]


def _cubit_preflight_warn(cmd_list: list[str]) -> dict | None:
	"""Scan recent failure log for a similar input. If found, return a
	compact warning dict so the caller can include it in the response
	without blocking execution.

	Similarity = token Jaccard ≥ 0.6 between the incoming command list
	and any recent failure's input. Cheap, reuses the existing jsonl
	log, and fires only when there's real signal.
	"""
	try:
		recs = _fl.recent_failures("cubit", limit=30)
	except Exception:
		return None
	if not recs:
		return None
	incoming = set(_tokenize(" ".join(cmd_list)))
	if not incoming:
		return None
	best_rec = None
	best_j = 0.0
	for r in recs:
		rin = r.get("input")
		if not rin:
			continue
		rs = " ".join(rin) if isinstance(rin, list) else str(rin)
		toks = set(_tokenize(rs))
		if not toks:
			continue
		inter = len(incoming & toks)
		union = len(incoming | toks)
		if union == 0:
			continue
		j = inter / union
		if j > best_j:
			best_j = j
			best_rec = r
	if best_rec is None or best_j < 0.6:
		return None
	return {
		"similar_past_failure": True,
		"jaccard": round(best_j, 2),
		"past_error": str(best_rec.get("error") or "")[:300],
		"past_failed_cmd": best_rec.get("failed_cmd"),
		"past_iso": best_rec.get("iso"),
		"advice": ("A very similar command sequence failed recently. "
		           "Review the past error above and adjust before proceeding. "
		           "Use `cubit_recent_failures(limit=10)` for full detail."),
	}


def _cubit_failure_hint(cmd: str, err: str) -> str:
	"""Pick a short hint from common Cubit error patterns."""
	blob = (cmd + " " + err).lower()
	for needle, advice in _CUBIT_ERROR_HEURISTICS:
		if needle in blob:
			return advice
	return ("No specific heuristic matched. Try `cubit_lookup(query=...)` "
	        "or `cubit_web_docs(source='cubit', query=...)` with keywords "
	        "from the error.")


@mcp.tool()
def cubit_probe(query: str = "summary") -> str:
	"""
	Query the Cubit session for geometry/mesh statistics.

	Valid queries:
	  "summary"    — {volumes, surfaces, curves, vertices, nodes, hexes, tets}
	  "quality"    — scaled-Jacobian stats over all hexes+tets
	  "per_volume" — one row per volume: {id, scheme, hex, tet, meshed}
	  "entities"   — Probe-Don't-Guess dump: every volume {id, centroid,
	               bbox_min/max, extent, volume} and surface {id, center,
	               bbox_min/max, extent, area}. ALWAYS run this before
	               writing entity-classification predicates
	               (webcut/subtract/imprint renumber entities).
	  "labels"     — blocks/sidesets with names, ACTUAL element membership,
	               and a convention audit (mixed blocks, unnamed blocks,
	               casefold collisions, non-snake-case names). Run after
	               every block/sideset command and before `export netgen`:
	               Cubit 2025.12 silently no-ops wrong-kind adds (e.g.
	               `add tri` on a hex mesh), so intended-but-absent
	               membership only shows up here.
	  "volume_count", "surface_count", "curve_count", "vertex_count"
	  "node_count", "hex_count", "tet_count"

	Both session transports (GUI file-drop and batch stdio) dispatch to
	the same shared implementation (probe_ops.py) — identical queries
	everywhere.

	Args:
	    query: probe name, see list above.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	try:
		r = sess.call("probe", [query])
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	return json.dumps(r, indent=2)


@mcp.tool()
def cubit_snapshot(out_path: str,
                   width: int = 800,
                   height: int = 600):
	"""
	Hardcopy the current Cubit view to a PNG file.

	Use for paper/slide figures from an LLM-driven session, or for
	a cheap "visual diff" after a command sequence.  GUI mode is
	required in practice: batch (-nographics) has no rendering window
	and the daemon reports ok=false with 0 bytes (the written file is
	verified -- Cubit's rc alone lies for silently-failing hardcopies).

	The captured PNG is returned INLINE as image content (so the model
	sees the current view directly) in addition to being written to
	`out_path` for reuse in documents.

	Args:
	    out_path: PNG output path (absolute or relative to project root).
	    width, height: ACCEPTED BUT IGNORED (reported back as
	        `requested_size_ignored`): the `hardcopy ... window w h`
	        variant writes 0-byte files in GUI mode (measured Cubit
	        2025.12), so the image renders at the current
	        graphics-window size.
	"""
	from pathlib import Path as _P
	from mcp.server.fastmcp import Image as _Image
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	p = _P(out_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	abs_path = str(p.resolve()).replace("\\", "/")
	try:
		r = sess.call("snapshot", [abs_path, int(width), int(height)])
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	result_json = json.dumps(r, indent=2)
	png = _P(abs_path)
	if r.get("ok") and png.is_file() and png.stat().st_size > 0:
		# Inline image + the JSON record (MathWorks figure-capture
		# pattern: the model sees the view, nothing to chase on disk).
		return [_Image(path=str(png)), result_json]
	return result_json


@mcp.tool()
def cubit_doctor() -> str:
	"""
	One-shot environment diagnosis for the whole Cubit MCP stack.

	Run this FIRST when anything is wrong (session won't start, license
	errors, exports missing labels, stale plugin suspicion). Read-only:
	launches nothing, consumes no license seat.

	Checks: Cubit install discovery -> Learn-license renewals cache
	freshness -> deployed Cubit plugin (.ccm) vs the cubit-mesh-export
	bundled copy (hash) -> live daemon state -> drop-dir diagnostics
	(startup_error.txt, log tails) -> check-vol dependency availability.

	Returns JSON: per-check {status: ok|warn|error|skipped, ...detail}
	plus an overall summary listing the problems found.
	"""
	import hashlib

	checks: dict = {}
	problems: list[str] = []

	# --- 1. Cubit install ------------------------------------------------
	bin_dir = _cs.find_cubit_install()
	if bin_dir is None:
		checks["install"] = {
			"status": "error",
			"detail": "Coreform Cubit not found (env > registry > glob)",
			"fix": "Install Cubit or set CUBIT_BIN_DIR / CUBIT_INSTALL_DIR",
		}
		problems.append("install: Cubit not found")
	else:
		checks["install"] = {"status": "ok", "bin_dir": str(bin_dir)}

	# --- 2. License renewals cache (no launch, no seat) ------------------
	# The renewals cache exists only for Learn-edition (login) licenses;
	# Pro / RLM-server machines have none and need no login -- absence is
	# NOT a problem there, so it reports "skipped", not "warn".
	try:
		from .license_warmup import _needs_login
		needs, info = _needs_login()
		if not needs:
			checks["license"] = {"status": "ok", **info}
		elif info.get("reason") == "no renewals cache":
			checks["license"] = {
				"status": "skipped", **info,
				"note": ("no Learn renewals cache -- normal on Pro/"
				         "RLM-server machines; Learn users get a ~30 s "
				         "login on first start"),
			}
		else:
			checks["license"] = {
				"status": "warn", **info,
				"fix": ("Learn license cache is stale/expiring: first "
				        "session start will re-login (~30 s, needs "
				        "network) -- pre-warm via rlm_activate."),
			}
			problems.append(f"license: {info.get('reason')}")
	except Exception as exc:
		checks["license"] = {"status": "skipped",
		                     "detail": f"warmup module: {exc}"}

	# --- 3. Deployed plugin freshness vs cubit-mesh-export ---------------
	def _sha256(p: Path) -> str:
		h = hashlib.sha256()
		h.update(p.read_bytes())
		return h.hexdigest()

	try:
		import cubit_mesh_export as _cme
		bundled = Path(_cme.__file__).parent / "cubit_mesh_export.ccm"
		deployed = (Path(bin_dir) / "plugins" / "cubit_mesh_export.ccm"
		            if bin_dir else None)
		if not bundled.is_file():
			checks["plugin"] = {"status": "skipped",
			                    "detail": "no bundled .ccm in cubit_mesh_export"}
		elif deployed is None:
			checks["plugin"] = {"status": "skipped",
			                    "detail": "no Cubit install to check against"}
		elif not deployed.is_file():
			checks["plugin"] = {
				"status": "error",
				"detail": f"plugin not deployed: {deployed}",
				"fix": "run `cubit-plugin-install`",
			}
			problems.append("plugin: not deployed")
		else:
			same = _sha256(bundled) == _sha256(deployed)
			checks["plugin"] = {
				"status": "ok" if same else "warn",
				"deployed": str(deployed),
				"matches_bundled": same,
			}
			if not same:
				checks["plugin"]["fix"] = (
					"Deployed .ccm differs from the installed "
					"cubit-mesh-export copy -- exports may drift (the "
					"2026-04-14 sideset incident class). Re-run "
					"`cubit-plugin-install`.")
				problems.append("plugin: deployed .ccm is stale vs package")
	except ImportError:
		checks["plugin"] = {"status": "skipped",
		                    "detail": "cubit-mesh-export not installed"}

	# --- 4. Daemon state --------------------------------------------------
	daemon = {"alive": False, "mode": None, "pid": None}
	if _cs._SINGLETON is not None:
		daemon["alive"] = _cs._SINGLETON.is_alive()
		daemon["owned"] = _cs._SINGLETON._owned
		if _cs._SINGLETON._proc is not None:
			daemon["pid"] = _cs._SINGLETON._proc.pid
	daemon["status"] = "ok"
	checks["daemon"] = daemon

	# --- 5. Drop-dir diagnostics -----------------------------------------
	drop_diag: dict = {"status": "ok"}
	try:
		drop = _cs._user_daemon_dir()
		drop_diag["drop_dir"] = str(drop)
		startup_err = drop / "startup_error.txt"
		if startup_err.is_file():
			drop_diag["status"] = "warn"
			drop_diag["startup_error"] = startup_err.read_text(
				encoding="utf-8", errors="replace")[-1500:]
			problems.append("drop-dir: startup_error.txt present "
			                "(last spawn failed)")
		stderr_log = drop / "cubit_stderr.log"
		if stderr_log.is_file():
			tail = stderr_log.read_bytes()[-800:]
			if tail.strip():
				drop_diag["stderr_tail"] = tail.decode("utf-8",
				                                       errors="replace")
	except Exception as exc:
		drop_diag = {"status": "skipped", "detail": str(exc)}
	checks["drop_dir"] = drop_diag

	# --- 6. check-vol dependency -----------------------------------------
	try:
		import ngsolve  # noqa: F401
		import cubit_mesh_export.check  # noqa: F401
		checks["check_vol_deps"] = {"status": "ok"}
	except ImportError as exc:
		checks["check_vol_deps"] = {
			"status": "warn",
			"detail": f"cubit_check_vol unavailable: {exc}",
		}
		problems.append("check-vol: dependencies missing")

	return json.dumps({
		"status": "ok" if not problems else "problems_found",
		"problems": problems,
		"checks": checks,
	}, ensure_ascii=False, indent=2, default=str)


@mcp.tool()
def cubit_session_journal(out_path: str = "",
                          include_failed: bool = True) -> str:
	"""
	Export every Cubit command this MCP-server process sent to the live
	session as a portable `.jou` journal.

	The lab treats `.jou` as the single source of truth for geometry;
	this tool turns an interactive LLM-driven session back into a
	replayable, version-controllable journal. Failed commands are
	included as comments (`# FAILED: ...`) so the record stays honest.

	Scope: commands from THIS server process only -- attaching to a
	daemon another process drove earlier does not recover its history.

	Args:
	    out_path: optional path to also write the journal file
	        (absolute, or relative to the repo root).
	    include_failed: include failed commands as comments (default on).

	Returns JSON {n_commands, n_failed, journal, path?}.
	"""
	sess = _cs._SINGLETON
	history = list(sess._command_history) if sess is not None else []
	if not history:
		return json.dumps({
			"status": "ok", "n_commands": 0, "n_failed": 0,
			"journal": "",
			"note": ("no commands recorded in this server process yet "
			         "(history is per-process, not per-daemon)"),
		})
	lines = [
		"# Session journal exported by mcp-server-cubit",
		f"# {time.strftime('%Y-%m-%d %H:%M:%S')}  "
		f"({len(history)} commands)",
		"# Replay: Cubit > Play Journal, or cubit_show(path=...)",
	]
	n_failed = 0
	for entry in history:
		if entry["ok"]:
			lines.append(entry["line"])
		else:
			n_failed += 1
			if include_failed:
				lines.append(f"# FAILED: {entry['line']}")
	journal = "\n".join(lines) + "\n"
	result = {
		"status": "ok",
		"n_commands": len(history),
		"n_failed": n_failed,
		"journal": journal,
	}
	if out_path:
		p = Path(out_path)
		if not p.is_absolute():
			p = PROJECT_ROOT / p
		try:
			p.write_text(journal, encoding="utf-8")
			result["path"] = str(p)
		except OSError as exc:
			result["write_error"] = str(exc)
	return json.dumps(result, ensure_ascii=False, indent=2)


@mcp.tool()
def cubit_session_status() -> str:
	"""Return diagnostic info about the Cubit session: alive/pid/mode,
	ownership, session-mode policy, drop dir, last license warmup, and
	how many commands this process has recorded for
	`cubit_session_journal`."""
	status = {
		"bin_dir": str(_cs.get_cubit_bin_dir() or "<not found>"),
		"alive": False,
		"pid": None,
		"ready_info": None,
		"session_mode": os.environ.get("RADIA_CUBIT_SESSION_MODE", "auto"),
	}
	if _cs._SINGLETON is not None:
		sess = _cs._SINGLETON
		status["owned"] = sess._owned
		status["drop_dir"] = (str(sess._drop_dir)
		                      if sess._drop_dir is not None else None)
		status["license_warmup"] = sess._last_license_warmup or None
		status["n_journal_commands"] = len(sess._command_history)
		status["alive"] = _cs._SINGLETON.is_alive()
		if _cs._SINGLETON._proc is not None:
			status["pid"] = _cs._SINGLETON._proc.pid
			status["mode"] = "owned"
		elif _cs._SINGLETON._drop_dir is not None:
			# Phase-1 attached: read PID from pid.lock
			pid_file = _cs._SINGLETON._drop_dir / "pid.lock"
			if pid_file.exists():
				try:
					status["pid"] = int(pid_file.read_text(encoding="utf-8").strip())
					status["mode"] = "attached"
				except (OSError, ValueError):
					pass
		status["ready_info"] = _cs._SINGLETON._ready_info
	return json.dumps(status, indent=2)


@mcp.tool()
def cubit_session_shutdown() -> str:
	"""Stop the persistent Cubit daemon. Next `cubit_show` relaunches.

	This is the EXPLICIT stop path: it also stops a daemon started by
	another process (e.g. a hung session that recovery deliberately
	left running). The report says which process was stopped
	(stopped: "owned-child" | "attached-daemon" | "none", plus pid).
	"""
	if _cs._SINGLETON is None or not _cs._SINGLETON.is_alive():
		return json.dumps({"status": "ok", "note": "no session running"})
	report = _cs.CubitSession.reset()
	return json.dumps({"status": "ok", "note": "session stopped", **report})


# ============================================================
# Checkpoint / restore (Cubit .cub5 snapshot-based undo)
# ============================================================

def _checkpoints_dir() -> Path:
	"""Per-user directory for named Cubit state snapshots."""
	d = Path.home() / ".cubit_viewer" / "checkpoints"
	d.mkdir(parents=True, exist_ok=True)
	return d


def _sanitize_label(label: str) -> str:
	"""Keep only safe filename characters; reject empty."""
	import re
	s = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("_")
	if not s:
		raise ValueError("checkpoint label must not be empty after sanitization")
	if len(s) > 64:
		s = s[:64]
	return s


@mcp.tool()
def cubit_checkpoint(label: str) -> str:
	"""Save the current Cubit session state as a named checkpoint.

	Creates `~/.cubit_viewer/checkpoints/<label>.cub5`. Overwrites if the
	label already exists. Use `cubit_restore(label)` to rewind.

	Lightweight undo pattern: save a checkpoint before risky operations
	(mesh, boolean, large imports), then restore if the result is bad.

	Args:
	    label: checkpoint name (alphanumeric, dot, dash, underscore only).
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	try:
		safe = _sanitize_label(label)
	except ValueError as e:
		return json.dumps({"status": "error", "error": str(e)})
	path = _checkpoints_dir() / f"{safe}.cub5"
	fwd = str(path).replace("\\", "/")
	try:
		r = sess.call("cmd", [f'save as "{fwd}" overwrite'], timeout_s=120)
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	if not r.get("ok"):
		return json.dumps({"status": "error", "stage": "save", "kind": "input", "detail": r})
	return json.dumps({
		"status": "ok",
		"label": safe,
		"path": str(path),
		"size_bytes": path.stat().st_size if path.exists() else 0,
	}, indent=2)


@mcp.tool()
def cubit_restore(label: str) -> str:
	"""Restore a previously-saved Cubit checkpoint by label.

	Loads `~/.cubit_viewer/checkpoints/<label>.cub5` via Cubit's `open`
	command. The current session state is replaced.

	Args:
	    label: name passed to a prior `cubit_checkpoint` call.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	try:
		safe = _sanitize_label(label)
	except ValueError as e:
		return json.dumps({"status": "error", "error": str(e)})
	path = _checkpoints_dir() / f"{safe}.cub5"
	if not path.exists():
		return json.dumps({
			"status": "error",
			"error": f"checkpoint not found: {safe}",
			"looked_in": str(_checkpoints_dir()),
		})
	fwd = str(path).replace("\\", "/")
	try:
		r = sess.call("cmd", [f'open "{fwd}"'], timeout_s=120)
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	if not r.get("ok"):
		return json.dumps({"status": "error", "stage": "open", "kind": "input", "detail": r})
	return json.dumps({"status": "ok", "label": safe, "path": str(path)}, indent=2)


@mcp.tool()
def cubit_list_checkpoints() -> str:
	"""List all saved Cubit checkpoints with size + mtime."""
	import datetime
	d = _checkpoints_dir()
	entries = []
	for p in sorted(d.glob("*.cub5")):
		try:
			st = p.stat()
			entries.append({
				"label": p.stem,
				"path": str(p),
				"size_bytes": st.st_size,
				"mtime": datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
			})
		except OSError:
			continue
	return json.dumps({"dir": str(d), "count": len(entries), "checkpoints": entries}, indent=2)


# ============================================================
# Meshing diagnostics
# ============================================================

# Suggested schemes by "why it failed" heuristics. Hex schemes first, tet
# fallback last.  Keyed on Cubit's get_volume_meshing_scheme returns that
# commonly indicate an unsuitable choice.
_SCHEME_ALTERNATIVES: dict = {
	"sweep": [
		"submap  # tight topology / brick-like volumes",
		"pave   # for 2D / degenerate sweep sources",
		"polyhedron  # robust hex for non-sweepable",
		"tetmesh  # fall back to tet",
	],
	"auto": [
		"sweep  # force sweep if prismatic",
		"polyhedron  # robust hex",
		"tetmesh  # fall back to tet",
	],
	"map": [
		"submap",
		"pave",
		"tetmesh",
	],
	"polyhedron": [
		"sweep  # if geometry is prismatic",
		"tetmesh",
	],
	None: [
		"Assign a scheme first: 'volume <id> scheme sweep' (or auto / tetmesh)",
	],
}


@mcp.tool()
def cubit_mesh_diagnose() -> str:
	"""Per-volume meshing diagnostic.

	Returns a table of every volume with: current meshing scheme, hex
	count, tet count, meshed / unmeshed flag. For unmeshed volumes,
	heuristic scheme alternatives are suggested (hex schemes first,
	tet fallback last).

	Use this after a `cubit_exec(["mesh volume all"])` that left some
	volumes unmeshed, to decide what to try next.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	try:
		r = sess.call("probe", ["per_volume"], timeout_s=30)
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	if not r.get("ok"):
		return json.dumps({"status": "error", "stage": "probe", "kind": "internal", "detail": r})
	rows = r.get("result", [])
	if not isinstance(rows, list):
		return json.dumps({"status": "error", "stage": "parse", "kind": "internal",
		                   "error": "per_volume did not return list",
		                   "raw": rows})

	meshed = [row for row in rows if row.get("meshed")]
	unmeshed = [row for row in rows if not row.get("meshed")]
	# Attach suggestions to unmeshed.
	for row in unmeshed:
		scheme = row.get("scheme")
		row["suggested_alternatives"] = _SCHEME_ALTERNATIVES.get(
			scheme, _SCHEME_ALTERNATIVES[None])

	hex_total = sum(row.get("hex", 0) for row in rows)
	tet_total = sum(row.get("tet", 0) for row in rows)

	return json.dumps({
		"status": "ok",
		"total_volumes": len(rows),
		"meshed_count": len(meshed),
		"unmeshed_count": len(unmeshed),
		"hex_total": hex_total,
		"tet_total": tet_total,
		"per_volume": rows,
		"hint": (
			"For each unmeshed volume, try the suggested scheme alternatives "
			"via `cubit_exec([\"volume <id> scheme <name>\", \"mesh volume <id>\"])`. "
			"Prefer hex schemes; fall back to tetmesh if all else fails."
		) if unmeshed else "All volumes meshed.",
	}, indent=2)


# ============================================================
# Tutorial / next-step suggestions
# ============================================================

@mcp.tool()
def cubit_suggest_next(goal: str = "mesh") -> str:
	"""Suggest concrete next Cubit commands based on the current state.

	Inspects the live session (via `probe summary` + per-volume info) and
	proposes 2–5 commands the LLM (or user) can execute next to progress
	toward the goal. Rule-based, not model-based — fast, predictable,
	and independent of external LLMs.

	Args:
	    goal: one of "import" / "mesh" / "hex" / "export" / "inspect" /
	        "refine" / "clean". Defaults to "mesh".

	Returns JSON with current state summary and a ranked list of
	suggested commands, each with a rationale.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	try:
		summary = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
	except _cs.CubitSessionError as e:
		return json.dumps(_error_payload("rpc", str(e)))
	try:
		per_vol = sess.call("probe", ["per_volume"], timeout_s=30).get("result", [])
	except _cs.CubitSessionError:
		per_vol = []

	g = (goal or "mesh").strip().lower()
	suggestions: list[dict] = []
	state = {
		"volumes": summary.get("volumes", 0),
		"surfaces": summary.get("surfaces", 0),
		"hexes": summary.get("hexes", 0),
		"tets": summary.get("tets", 0),
		"nodes": summary.get("nodes", 0),
	}
	any_mesh = state["hexes"] + state["tets"] > 0
	unmeshed = [row for row in per_vol if not row.get("meshed")]

	def add(cmd: str, why: str) -> None:
		suggestions.append({"cmd": cmd, "why": why})

	if g == "import" or state["volumes"] == 0:
		add('import step "<path>"', "Load geometry from STEP (build123d output, CAD exchange).")
		add('import acis "<path>"', "Load ACIS .sat / .brep (native Cubit precision).")
		add('import mesh "<path>"', "Load existing mesh (.exo / .g / .msh / .vol).")
	elif g in ("mesh", "hex"):
		if unmeshed:
			# Most common next step after import
			if g == "hex":
				add("volume all scheme sweep", "Try hex via sweep (works if volumes are prismatic).")
				add("volume all size auto factor 5", "Moderately fine mesh size auto-chosen from geometry.")
				add("mesh volume all", "Execute the mesh.")
				add("volume all scheme polyhedron", "Robust hex when sweep fails (non-prismatic).")
			else:
				add("volume all scheme auto", "Let Cubit pick sweep/sub/map/tet per volume.")
				add("volume all size auto factor 5", "Moderately fine mesh size.")
				add("mesh volume all", "Execute the mesh.")
			add("volume all scheme tetmesh", "Tet fallback — always works, no topology constraints.")
		else:
			add("quality volume all", "Check mesh quality (aspect / skew / jacobian).")
			add("mesh volume all", "Already meshed; re-run only after mods.")
	elif g == "export":
		if not any_mesh:
			add("mesh volume all", "No mesh yet — create one before export.")
		add('export mesh "<path>.msh" overwrite', "Gmsh .msh (preferred lab format).")
		add('export abaqus "<path>.inp" overwrite', "Abaqus .inp for FEM.")
		add('export genesis "<path>.g" overwrite', "Exodus .g for Sierra / Netgen pipelines.")
		add('export step "<path>.step" overwrite', "Geometry-only STEP (for round-trip CAD).")
	elif g == "inspect":
		add("list volume all", "Print geometry summary.")
		add("list hex in volume all", "Print mesh counts per volume.")
		add("quality volume all", "Mesh quality metrics.")
		add("draw volume all", "Force a redraw if GUI out of sync.")
	elif g == "refine":
		add("refine volume all numsplit 1", "Uniform 1-level refine (~8× cells).")
		add("smooth volume all", "Laplace smoothing to improve quality.")
		add("quality volume all", "Re-check quality after refine.")
	elif g == "clean":
		add("delete mesh", "Drop the current mesh, keep geometry.")
		add("reset genesis", "Reset mesh numbering.")
		add("delete volume all", "Nuke geometry too.")
	else:
		add("# unknown goal; try: import / mesh / hex / export / inspect / refine / clean", "")

	return json.dumps({
		"status": "ok",
		"goal": g,
		"state": state,
		"unmeshed_volumes": len(unmeshed),
		"suggestions": suggestions,
	}, indent=2)


# ============================================================
# Knowledge lookup (retrieval-based)
# ============================================================

# Sources searched by `cubit_lookup`. Each entry is (name, accessor).
# The accessor is called with no args; returns the full knowledge body.
_KNOWLEDGE_SOURCES: list = [
	("export", lambda: get_export_documentation()),
	("cubit_api", lambda: get_api_reference()),
	("scripting", lambda: get_cubit_documentation()),
	("netgen_workflow", lambda: get_netgen_documentation()),
	("forum_tips", lambda: get_forum_tips()),
]


import math as _math
import re as _re


_STOP_WORDS = frozenset({
	"the", "a", "an", "of", "to", "in", "and", "or", "is", "are", "be",
	"on", "at", "for", "with", "as", "by", "from", "that", "this",
	"it", "its", "if", "so", "how", "what", "which", "when", "where",
})

_WORD_RE = _re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*")


def _tokenize(s: str) -> list[str]:
	return [m.group(0).lower() for m in _WORD_RE.finditer(s)]


def _chunk_text(text: str, target_lines: int = 40) -> list[tuple[int, str]]:
	"""Split a knowledge string into heading-anchored chunks.

	Lines starting with "# ", "## ", "### " anchor chunks. Returns list
	of (starting_line, chunk_text). Chunks longer than 2× target_lines
	are hard-broken on line count.
	"""
	lines = text.splitlines()
	chunks: list[tuple[int, str]] = []
	cur_start = 0
	cur: list[str] = []
	for i, line in enumerate(lines):
		is_heading = line.startswith("# ") or line.startswith("## ") \
			or line.startswith("### ")
		if is_heading and cur and (i - cur_start) >= max(5, target_lines // 4):
			chunks.append((cur_start + 1, "\n".join(cur)))
			cur = []
			cur_start = i
		cur.append(line)
		# Hard break if chunk grows too large without a heading.
		if len(cur) >= target_lines * 2:
			chunks.append((cur_start + 1, "\n".join(cur)))
			cur = []
			cur_start = i + 1
	if cur:
		chunks.append((cur_start + 1, "\n".join(cur)))
	return chunks


def _heading_of(chunk: str) -> str:
	"""Return first non-blank line of a chunk (best heading proxy)."""
	for line in chunk.splitlines():
		s = line.strip()
		if s:
			return s
	return ""


def _build_corpus() -> list[dict]:
	"""Chunk every knowledge source + failure log once. Result is cached.

	Returns: list of {source, start_line, chunk, heading, tokens, token_set}.
	"""
	if getattr(_build_corpus, "_cache", None) is not None:
		return _build_corpus._cache  # type: ignore[attr-defined]

	sources = list(_KNOWLEDGE_SOURCES)
	# Include the failure log so past mistakes become searchable context.
	sources.append(("failures", lambda: _fl.failures_as_text("cubit", limit=50)))

	corpus: list[dict] = []
	for name, accessor in sources:
		try:
			text = accessor()
		except Exception:
			continue
		if not text:
			continue
		for start, chunk in _chunk_text(text, target_lines=40):
			tokens = _tokenize(chunk)
			if not tokens:
				continue
			corpus.append({
				"source": name,
				"start_line": start,
				"chunk": chunk,
				"heading": _heading_of(chunk),
				"tokens": tokens,
				"token_set": set(tokens),
			})
	_build_corpus._cache = corpus  # type: ignore[attr-defined]
	return corpus


def _invalidate_corpus_cache() -> None:
	"""Force re-chunking (after a new failure is logged, knowledge reloaded)."""
	if hasattr(_build_corpus, "_cache"):
		delattr(_build_corpus, "_cache")


def _score_tfidf(query_terms: list[str], doc: dict, idf: dict[str, float]) -> float:
	"""tf-idf-like score + heading boost. Zero if no term matches."""
	tokens = doc["tokens"]
	if not tokens:
		return 0.0
	tf: dict[str, int] = {}
	for tok in tokens:
		tf[tok] = tf.get(tok, 0) + 1
	inv_len = 1.0 / (1.0 + _math.log(1 + len(tokens)))
	score = 0.0
	hits = 0
	for term in query_terms:
		f = tf.get(term, 0)
		if f == 0:
			continue
		hits += 1
		score += (1.0 + _math.log(f)) * idf.get(term, 1.0) * inv_len
	if hits == 0:
		return 0.0
	# Boost chunks whose heading contains query terms
	head_lc = doc["heading"].lower()
	head_hits = sum(1 for t in query_terms if t in head_lc)
	if head_hits:
		score *= 1.0 + 0.75 * head_hits
	# Coverage bonus: matching many distinct query terms is better than one term many times
	score *= 1.0 + 0.25 * (hits / max(1, len(query_terms)))
	return score


def _compute_idf(corpus: list[dict], terms: list[str]) -> dict[str, float]:
	n = max(1, len(corpus))
	idf: dict[str, float] = {}
	for t in terms:
		df = sum(1 for doc in corpus if t in doc["token_set"])
		idf[t] = _math.log((n + 1) / (df + 1)) + 1.0
	return idf


@mcp.tool()
def cubit_lookup(query: str, max_results: int = 5, chunk_lines: int = 40) -> str:
	"""Retrieve relevant sections from the bundled Cubit knowledge base.

	Knowledge (~8000 lines: API reference, scripting patterns, forum
	tips, export rules, netgen workflows) is chunked by heading; the
	failure log (`cubit_recent_failures`) is mixed in as an extra
	source so past mistakes are searchable by the same query.

	Retrieval:
	  - tokenized lower-case word matching (not substring)
	  - tf-idf-style scoring (rare terms count more)
	  - ×1.75/term heading-title boost
	  - coverage bonus for hitting many distinct query terms

	Args:
	    query: keywords (e.g., "sweep scheme source", "export netgen
	        high order", "quality metrics").
	    max_results: cap on returned chunks (default 5).
	    chunk_lines: target chunk size before heading-based re-anchoring
	        (reserved; corpus is built with 40 and cached).

	Returns JSON with top chunks annotated by source module and starting
	line number for traceability.
	"""
	del chunk_lines  # corpus is pre-chunked & cached; param kept for compat
	q = (query or "").strip()
	if not q:
		return json.dumps({"status": "error", "kind": "input", "error": "empty query"})
	terms = [t for t in _tokenize(q) if t not in _STOP_WORDS]
	if not terms:
		return json.dumps({"status": "error",
		                   "error": "no searchable terms (all stop words?)"})

	# Refresh if failure log grew since cache build (cheap heuristic).
	corpus = _build_corpus()
	idf = _compute_idf(corpus, terms)

	scored: list[tuple[float, dict]] = []
	for doc in corpus:
		s = _score_tfidf(terms, doc, idf)
		if s > 0:
			scored.append((s, doc))

	scored.sort(key=lambda x: (-x[0], x[1]["source"], x[1]["start_line"]))
	top = scored[: max(1, int(max_results))]
	out = [
		{
			"source": doc["source"],
			"start_line": doc["start_line"],
			"heading": doc["heading"][:120],
			"score": round(score, 3),
			"chunk": doc["chunk"],
		}
		for score, doc in top
	]
	return json.dumps({
		"status": "ok",
		"query": q,
		"terms": terms,
		"candidates_total": len(scored),
		"returned": len(out),
		"results": out,
	}, indent=2)


@mcp.tool()
def cubit_recent_failures(limit: int = 10) -> str:
	"""Return the last N failed Cubit invocations (from persistent log).

	Use this after an error to remember what was tried recently, or at
	the start of a session to avoid repeating a known-bad command.
	"""
	recs = _fl.recent_failures("cubit", limit=limit)
	return json.dumps({
		"status": "ok",
		"count": len(recs),
		"log_path": str(_fl._log_path("cubit")),
		"entries": recs,
	}, indent=2)


@mcp.tool()
def cubit_web_docs(query: str, source: str = "cubit", max_hits: int = 5,
                   force_refresh: bool = False) -> str:
	"""Fetch live Cubit documentation and grep for `query`.

	Complements `cubit_lookup` (bundled static knowledge): when the
	bundled kb lacks context or an error references a newly-added
	command, hit the live docs. Results are cached on disk for 24 h.

	Args:
	    query: search terms (AND-matched, case-insensitive).
	    source: "cubit" (coreform.com/cubit_help) or
	        "cubit_forum" (forum.coreform.com).
	    max_hits: cap on returned matching lines across all pages.
	    force_refresh: bypass cache and re-fetch.

	Returns JSON with per-page hit excerpts (line + ±4 context lines).
	"""
	src = (source or "cubit").strip().lower()
	q = (query or "").strip()
	if not q:
		return json.dumps({"status": "error", "kind": "input", "error": "empty query"})

	# Discourse-backed sources: use the JSON search API (richer + reliable).
	if src == "cubit_forum":
		hits = _wd.search_forum(q, max_hits=max_hits, force_refresh=force_refresh)
		return json.dumps({"status": "ok", "source": src, "query": q,
		                   "hits": hits}, indent=2)

	urls = _wd.DOCS_INDEX.get(src)
	if not urls:
		return json.dumps({"status": "error",
		                   "error": f"unknown source {source!r}",
		                   "known_sources": list(_wd.DOCS_INDEX.keys())})

	all_hits: list[dict] = []
	errors: list[dict] = []
	remaining = max(1, int(max_hits))
	for url in urls:
		if remaining <= 0:
			break
		data = _wd.fetch(url, force_refresh=force_refresh)
		if data.get("status") != "ok":
			errors.append({"url": url, "error": data.get("error")})
			continue
		hits = _wd.search_text(data["text"], q, max_hits=remaining)
		if hits:
			for h in hits:
				h["url"] = data["url"]
				h["from_cache"] = data.get("from_cache", False)
				all_hits.append(h)
			remaining -= len(hits)
	return json.dumps({
		"status": "ok",
		"source": src,
		"query": q,
		"hits": all_hits,
		"errors": errors,
	}, indent=2)


@mcp.tool()
def cubit_examples(query: str, limit: int = 3,
                    force_refresh: bool = False) -> str:
	"""Search Cubit journal examples from **multiple unioned sources**.

	Sources searched (all unioned, tf-idf ranked together):
	  1. **Coreform forum** (Discourse): 15 seed queries (tutorial,
	     hex meshing, sweep, boundary layer, thin shell, quality,
	     refinement, abaqus, high-order, polyhedron, multi-sweep,
	     imprint/merge, webcut, export, journal)
	  2. **Local lab archive**: `public-safe curated corpus` (~145 .jou/.py
	     across years of research projects — twisted wire, claw-pole
	     alternator, kelvin transformation, TEAM problems, helical
	     coils, JMAG integration, Mesh AI, Nastran/VTU pipelines)
	     + durable Radia repo lanes (`docs/`, `validation_test/`,
	     `src/radia/panels/samples`).
	     Lab-curated items get a small score boost.

	First call scrapes forum + walks local dirs; cache lives 7 days.
	Each hit returns: `source`, title, URL (`file:///…` for local),
	score, excerpt (first ~1.6 kB), absolute path.

	Args:
	    query: keywords (AND-scored, tokens + underscore-splits).
	    limit: cap on returned examples.
	    force_refresh: re-scrape / re-walk both sources.
	"""
	r = _ex.search_examples("cubit", query, limit=limit,
	                         force_refresh=force_refresh)
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def cubit_examples_refresh(limit_per_query: int = 5,
                           local_roots: list | None = None,
                           include_youtube: bool = True,
                           include_training_zip: bool = True,
                           include_jou_github: bool = False) -> str:
	"""Force-refresh every Cubit example sub-source.

	Sub-sources:
	  * forum (Discourse) — always
	  * local (public-safe curated corpus + Radia durable lanes + extras) — always
	  * cubit_youtube — if `include_youtube`
	  * Coreform training .zip → folded into local — if
	    `include_training_zip`
	  * cubit_jou_github (GitHub `.jou` code search, PAT-gated) — if
	    `include_jou_github`

	Args:
	    limit_per_query: forum seed-query post cap.
	    local_roots: override default lab archive roots.
	    include_youtube: scrape Coreform YouTube tutorials.
	    include_training_zip: download examples_only.zip + reindex local.
	    include_jou_github: GitHub code-search for `.jou` files.
	"""
	out: dict = {}
	out["forum"] = _ex.refresh_cubit_examples(limit_per_query=limit_per_query)
	out["local"] = _ex.refresh_cubit_local_examples(
		roots=[str(r) for r in local_roots] if local_roots else None,
	)
	if include_training_zip:
		out["coreform_training_zip"] = _ex.refresh_coreform_training_zip()
	if include_youtube:
		out["youtube"] = _ex.refresh_cubit_youtube()
	if include_jou_github:
		out["jou_github"] = _ex.refresh_jou_github_code_search()
	return json.dumps(out, indent=2, ensure_ascii=False)


# ============================================================
# Batch validation + mesh auto-fix (test-then-reflect pattern)
# ============================================================

def _run_batch(step_path: str | None, commands: list[str],
               timeout_s: int = 300,
               cub5_path: str | None = None) -> dict:
	"""Launch a FRESH, disposable batch Cubit (no GUI / no journal).

	Starting state options (mutually exclusive):
	  * `cub5_path` — open a .cub5 snapshot (used by `cubit_exec_safely`
	    to seed the batch with the current live GUI state).
	  * `step_path` — import a STEP file.
	  * neither — start from empty.

	Then run `commands`, probe state, tear down. Completely isolated
	from the live GUI singleton.
	"""
	sess = _cs.CubitSession(mode="batch")
	try:
		try:
			sess.ensure_started()
		except _cs.CubitSessionError as e:
			return {"status": "error", "stage": "start", "kind": "environment",
			        "error": str(e), "commands": commands}
		per_line: list[dict] = []
		if cub5_path:
			cub5_fwd = str(cub5_path).replace("\\", "/")
			try:
				r = sess.call("cmd", [f'open "{cub5_fwd}"'], timeout_s=timeout_s)
				preamble = r.get("result", [])
				per_line.extend(preamble)
				if not preamble or not preamble[0].get("ok"):
					return {"status": "error", "stage": "open_cub5", "kind": "input",
					        "cub5_path": cub5_path, "per_line": per_line}
			except _cs.CubitSessionError as e:
				return {"status": "error", "stage": "open_cub5_rpc",
				        "cub5_path": cub5_path, "error": str(e)}
		elif step_path:
			step_fwd = str(step_path).replace("\\", "/")
			try:
				r = sess.call("cmd", [f'import step "{step_fwd}"'], timeout_s=timeout_s)
				imp_result = r.get("result", [])
				per_line.extend(imp_result)
				if not imp_result or not imp_result[0].get("ok"):
					return {"status": "error", "stage": "import",
					        "step_path": step_path, "per_line": per_line}
			except _cs.CubitSessionError as e:
				return {"status": "error", "stage": "import_rpc",
				        "step_path": step_path, "error": str(e)}
		try:
			r = sess.call("cmd", commands, timeout_s=timeout_s)
		except _cs.CubitSessionError as e:
			return {**_error_payload("rpc", str(e)), "per_line": per_line}
		per_line.extend(r.get("result", []))
		try:
			smy = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
		except _cs.CubitSessionError:
			smy = {}
		all_ok = all(step.get("ok") for step in per_line) if per_line else False
		failed_at = None
		first_err = None
		for i, step in enumerate(per_line, start=1):
			if not step.get("ok"):
				failed_at = i
				first_err = step.get("error") or f"rc={step.get('rc')}"
				break
		return {
			"status": "ok" if all_ok else "error",
			"all_ok": all_ok,
			"per_line": per_line,
			"summary": smy,
			"failed_at": failed_at,
			"error": (first_err or "").strip()[:400] if first_err else None,
		}
	finally:
		try:
			sess.shutdown()
		except Exception:
			pass


def _netgen_mesh_to_msh(step_path: Path, maxh: float, out_msh: Path,
                        order: int = 1) -> dict:
	"""STEP -> Netgen tet mesh (optionally curved) -> GMSH .msh v4.1
	(via radia's exporter, honoring the lab's vol/v4.1-only interchange
	policy).  order=2 emits Tet10 with curved geometry, making the
	referee's min_jacobian_ratio a real curvature diagnostic (verified:
	sphere at order 2 reports jac_ratio ~0.79 vs 1.0 for straight)."""
	from netgen.occ import OCCGeometry
	from ngsolve import Mesh, TaskManager
	from radia.gmsh_post_export import GmshPostExport

	with TaskManager():
		geo = OCCGeometry(str(step_path))
		ngmesh = geo.GenerateMesh(maxh=maxh)
		mesh = Mesh(ngmesh)
		if order > 1:
			mesh.Curve(order)
		post = GmshPostExport(mesh)
		post.write(str(out_msh))
	return {"n_elements": mesh.ne, "n_vertices": mesh.nv}


def _quality_row(route: str, mesher: str, q: dict) -> dict:
	"""Flatten one referee report into a comparison row."""
	by_type = q.get("by_type", [])
	completed = bool(q.get("ran")) and bool(by_type)
	row = {
		"route": route,
		"mesher": mesher,
		"status": "ok" if completed else "error",
		"metric": q.get("metric"),
		"ok": q.get("ok"),
		"total_negative": q.get("total_negative"),
		"total_below_threshold": q.get("total_below_threshold"),
		"by_type": [
			{
				"element": bt.get("name"),
				"n": bt.get("n_elements"),
				"min": bt.get("min_quality"),
				"mean": bt.get("mean_quality"),
				"negative": bt.get("negative"),
				"below_threshold": bt.get("below_threshold"),
				"worst_tag": (bt.get("worst") or [{}])[0].get("tag"),
				# stretching -- minSICN is a shape number and misses it
				"aspect_p95": (bt.get("aspect_ratio") or {}).get("p95"),
				"aspect_max": (bt.get("aspect_ratio") or {}).get("max"),
				"aspect_above_10": (bt.get("aspect_ratio") or {}).get("n_above_10"),
			}
			for bt in by_type
		],
	}
	stats = q.get("mesh_stats") or {}
	if stats:
		# THE COST AXIS. Ranking meshes by element count inverts the
		# verdict (measured: netgen looks ~2x cheaper per element and is
		# never cheaper per dof); ranking by dof does not.
		row["dof_estimate"] = stats.get("dof_estimate")
		row["interior_node_fraction"] = stats.get("interior_node_fraction")
		row["n_nodes"] = stats.get("n_nodes")
	if not completed:
		row["kind"] = "referee"
		row["error"] = q.get("error") or q.get("note") \
			or "quality referee returned no 3D element data"
	return row


@mcp.tool()
def cubit_netgen_quality_compare(step_path: str,
                                 netgen_maxh: float = 0.0,
                                 cubit_size: float = 0.0,
                                 schemes: list = None,
                                 threshold: float = 0.1,
                                 order: int = 1,
                                 out_dir: str = "",
                                 timeout_s: int = 600) -> str:
	"""
	Mesh ONE STEP with both meshers and compare element quality with a
	NEUTRAL referee.

	Routes (headless throughout, per the driving policy):
	  * netgen     — STEP -> Netgen tet mesh (maxh) -> .msh v4.1
	  * cubit_tet  — batch Cubit, `volume all scheme tetmesh`, size
	  * cubit_hex  — batch Cubit, default scheme ladder (sweep/map where
	                 possible), size

	Referee: GMSH minSICN via the radia-gmsh quality engine
	(`radia_mcp.gmsh.msh_inspect.mesh_quality`) evaluated on each
	exported `.msh v4.1` -- ONE metric implementation for every mesher,
	so numbers are directly comparable. minSICN is defined for tets AND
	hexes, but comparing across element FAMILIES is still
	apples-to-oranges for meshing-difficulty reasons; the honest
	comparison is netgen vs cubit_tet (same family), with cubit_hex as
	the structured-mesh reference.

	When `netgen_maxh`/`cubit_size` are 0, a common size is derived from
	the geometry bbox (~1/8 of the max extent) so both meshers get the
	SAME target length.

	Writes `<base>_netgen.msh` / `<base>_cubit_tet.msh` /
	`<base>_cubit_hex.msh` next to the STEP (or into `out_dir`).
	Per-route failures are reported as rows with status/kind -- a
	missing optional stack (netgen/radia/gmsh) fails that ROW as
	environment, not the whole comparison.

	INTERPRETING min_quality -- one call is ONE SAMPLE, not a verdict.
	Both meshers are bit-exactly reproducible at a FIXED size target, but
	the worst element is a chaotic function of that target: measured on
	c_core, a 0.125 % change of the size target moved min by 0.118
	(netgen) and 0.086 (cubit_tet) -- bigger than most mesher differences
	(validation_test/radia_mcp/mesh_quality_study). So compare BANDS over
	several sizes, never two single numbers; and note that at these levels
	min does NOT predict solver behaviour (order-1 Poisson CG counts and
	L2/H1 error barely respond). Treat min_quality as a floor/safety
	indicator -- `negative > 0` is the signal that actually matters.

	RANK BY dof, NOT BY ELEMENT COUNT. Each row carries `dof_estimate`
	(H1-p1 / lowest HCurl / lowest HDiv dof on that mesh) and
	`interior_node_fraction`. The linear system is sized by dof, and the
	two axes disagree: netgen produces ~2x fewer ELEMENTS at equal h yet
	is never cheaper per dof, because more of its nodes sit on the
	boundary (interior fraction 2-49 % vs cubit_tet 19-77 %). Measured at
	matched dof, netgen's error is 1.14-1.89x cubit_tet's on scalar H1
	and 1.02-1.27x on HCurl/HDiv. An element-count comparison reverses
	this and is the single easiest way to reach a wrong conclusion here.

	`aspect_p95` / `aspect_max` report stretching, which minSICN does not
	express and which is what thin-gap and lamination meshes actually
	exhibit.

	Args:
	    step_path: the geometry both meshers consume.
	    netgen_maxh: Netgen max element size (0 = derive from bbox).
	    cubit_size: Cubit `volume all size` (0 = derive from bbox).
	    schemes: subset of ["netgen", "cubit_tet", "cubit_hex"]
	        (default: all three).
	    threshold: minSICN acceptance threshold for the referee.
	    order: element order 1-3 for BOTH routes (netgen `mesh.Curve`,
	        cubit `export gmsh ... order N`). At order>=2 on curved
	        geometry the referee's min_jacobian_ratio becomes a real
	        curvature-distortion diagnostic (1.0 for straight elements).
	    out_dir: where to write the .msh files (default: beside STEP).
	"""
	p = Path(step_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	if not p.is_file():
		return json.dumps(_error_payload("input", f"STEP not found: {p}"))
	try:
		order = int(order)
		netgen_maxh = float(netgen_maxh)
		cubit_size = float(cubit_size)
		threshold = float(threshold)
		timeout_s = int(timeout_s)
	except (TypeError, ValueError) as exc:
		return json.dumps(_error_payload(
			"input", f"invalid numeric option: {exc}"))
	if not 1 <= order <= 3:
		return json.dumps(_error_payload(
			"input", f"order must be 1..3 (gmsh export limit), got {order}"))
	if not _math.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
		return json.dumps(_error_payload(
			"input", f"threshold must be finite and in [0, 1], got {threshold}"))
	if timeout_s <= 0:
		return json.dumps(_error_payload(
			"input", f"timeout_s must be positive, got {timeout_s}"))
	for name, value in (("netgen_maxh", netgen_maxh),
	                    ("cubit_size", cubit_size)):
		if not _math.isfinite(value) or value < 0.0:
			return json.dumps(_error_payload(
				"input", f"{name} must be finite and non-negative, got {value}"))
	if schemes is None:
		schemes = ["netgen", "cubit_tet", "cubit_hex"]
	elif not isinstance(schemes, (list, tuple)) or not schemes:
		return json.dumps(_error_payload(
			"input", "schemes must be a non-empty list"))
	schemes = [str(s) for s in schemes]
	unknown = [s for s in schemes if s not in
	           ("netgen", "cubit_tet", "cubit_hex")]
	if unknown:
		return json.dumps(_error_payload(
			"input", f"unknown schemes: {unknown}; "
			         "valid: netgen / cubit_tet / cubit_hex"))
	if len(set(schemes)) != len(schemes):
		return json.dumps(_error_payload(
			"input", f"schemes contains duplicates: {schemes}"))
	out_base = Path(out_dir) if out_dir else p.parent
	if not out_base.is_absolute():
		out_base = PROJECT_ROOT / out_base
	try:
		out_base.mkdir(parents=True, exist_ok=True)
	except OSError as exc:
		return json.dumps(_error_payload(
			"input", f"cannot create out_dir {out_base}: {exc}"))
	base = p.stem

	# Common target size from the CAD bbox (same length for BOTH meshers)
	if netgen_maxh == 0.0 or cubit_size == 0.0:
		try:
			from netgen.occ import OCCGeometry as _OCC
			bb = _OCC(str(p)).shape.bounding_box
			ext = max(bb[1][i] - bb[0][i] for i in range(3))
			derived = ext / 8.0
		except Exception:
			derived = 0.0
		if netgen_maxh == 0.0:
			netgen_maxh = derived
		if cubit_size == 0.0:
			cubit_size = derived
	if (not _math.isfinite(netgen_maxh) or netgen_maxh <= 0.0
			or not _math.isfinite(cubit_size) or cubit_size <= 0.0):
		return json.dumps(_error_payload(
			"input", "could not derive a mesh size from the geometry; "
			         "pass netgen_maxh and cubit_size explicitly"))

	try:
		from ..gmsh.msh_inspect import mesh_quality
	except ImportError as exc:
		return json.dumps(_error_payload(
			"referee", f"gmsh referee unavailable: {exc}",
			kind="environment"))

	rows = []

	if "netgen" in schemes:
		msh = out_base / f"{base}_netgen.msh"
		try:
			info = _netgen_mesh_to_msh(p, float(netgen_maxh), msh,
			                           order=order)
			q = mesh_quality(msh, threshold=threshold)
			row = _quality_row("netgen", "netgen", q)
			row["maxh"] = float(netgen_maxh)
			row["order"] = order
			row["msh"] = str(msh)
			row.update(info)
		except ImportError as exc:
			row = {"route": "netgen", "status": "error",
			       "kind": "environment",
			       "error": f"netgen/ngsolve/radia unavailable: {exc}"}
		except Exception as exc:
			row = {"route": "netgen", "status": "error", "kind": "input",
			       "error": f"{type(exc).__name__}: {exc}"}
		rows.append(row)

	cubit_routes = [s for s in schemes if s.startswith("cubit_")]
	for route in cubit_routes:
		msh = out_base / f"{base}_{route}.msh"
		msh_fwd = str(msh).replace("\\", "/")
		cmds = []
		if route == "cubit_tet":
			cmds.append("volume all scheme tetmesh")
		cmds += [f"volume all size {float(cubit_size)}",
		         "mesh volume all",
		         # export gmsh emits only BLOCK members -- an unblocked
		         # mesh exports 0 elements (measured 2025.12).
		         "block 1 add volume all",
		         'block 1 name "mesh"',
		         f'export gmsh "{msh_fwd}"'
		         + (f" order {order}" if order > 1 else "")
		         + " overwrite"]
		r = _run_batch(str(p), cmds, timeout_s=timeout_s)
		if r.get("status") != "ok":
			rows.append({"route": route, "status": "error",
			             "kind": r.get("kind", "input"),
			             "error": r.get("error"),
			             "failed_at": r.get("failed_at")})
			continue
		try:
			q = mesh_quality(msh, threshold=threshold)
			row = _quality_row(route, "cubit", q)
			row["size"] = float(cubit_size)
			row["order"] = order
			row["msh"] = str(msh)
			row["summary"] = r.get("summary")
		except Exception as exc:
			row = {"route": route, "status": "error",
			       "kind": "referee",
			       "error": f"{type(exc).__name__}: {exc}",
			       "msh": str(msh)}
		rows.append(row)

	# Verdict: same-family comparison where possible
	notes = []
	def _min_of(route):
		for row in rows:
			if row.get("route") == route and row.get("by_type"):
				values = [bt.get("min") for bt in row["by_type"]
				          if isinstance(bt.get("min"), (int, float))]
				return min(values) if values else None
		return None
	ng, ct = _min_of("netgen"), _min_of("cubit_tet")
	if ng is not None and ct is not None:
		better = "cubit_tet" if ct > ng else "netgen"
		notes.append(
			f"tet-vs-tet (same family, directly comparable): min minSICN "
			f"netgen={ng:.3f} vs cubit_tet={ct:.3f} -> {better} has the "
			"better worst element at this size")
	if _min_of("cubit_hex") is not None:
		notes.append(
			"cubit_hex is a different element family -- report it as the "
			"structured-mesh reference, not as a same-metric winner")

	return json.dumps({
		"status": "ok",
		"all_routes_completed": bool(rows) and all(
			row.get("status") == "ok" for row in rows),
		"step": str(p),
		"referee": "gmsh minSICN (radia_mcp.gmsh.msh_inspect.mesh_quality)",
		"threshold": threshold,
		"rows": rows,
		"notes": notes,
	}, ensure_ascii=False, indent=2, default=str)


def _vol_surface_element_count(vol_path: Path) -> int:
	"""Count boundary surface elements recorded in a Netgen `.vol` file.

	The section header is `surfaceelements` (`...gi` / `...uv` variants for
	geometry-informed exports) followed by the element count on its own
	line.  Returns -1 when no section is present (a malformed export).
	"""
	tokens = ("surfaceelements", "surfaceelementsgi", "surfaceelementsuv")
	with open(vol_path, "r", encoding="utf-8", errors="replace") as f:
		for line in f:
			if line.strip() in tokens:
				try:
					return int(next(f).strip())
				except (StopIteration, ValueError):
					return -1
	return -1


def _vol_surface_element_stats(vol_path: Path) -> dict:
	"""Split a Netgen `.vol` surface section into OUTER and INTERFACE faces.

	Each `surfaceelements` row is `surfnr bcnr domin domout np p1 ...`.
	A face with one zero domain bounds the model (outer skin); a face with
	two nonzero domains is a MATERIAL INTERFACE, which a multi-material
	Sculpt mesh legitimately carries in addition to the skin.  Comparing
	the raw total against the topological skin therefore only works for a
	single material -- the outer count is the portable quantity.

	When the rows cannot be classified (an unexpected layout), the result
	falls back to ``outer == total`` with ``classified=False``, which
	reproduces the historical single-material comparison rather than
	inventing a split.  ``total`` is -1 when the section is missing.
	"""
	total = _vol_surface_element_count(vol_path)
	fallback = {"total": total, "outer": total, "interface": 0,
	            "classified": False}
	if total <= 0:
		return fallback
	tokens = ("surfaceelements", "surfaceelementsgi", "surfaceelementsuv")
	with open(vol_path, "r", encoding="utf-8", errors="replace") as f:
		for line in f:
			if line.strip() not in tokens:
				continue
			next(f)                      # the count line
			outer = interface = 0
			for _ in range(total):
				try:
					parts = next(f).split()
				except StopIteration:
					return fallback
				if len(parts) < 4:
					return fallback
				try:
					domin, domout = int(parts[2]), int(parts[3])
				except ValueError:
					return fallback
				if domin and domout:
					interface += 1
				else:
					outer += 1
			return {"total": total, "outer": outer,
			        "interface": interface, "classified": True}
	return fallback


def _collect_mesh_gates(msh: Path, vol: Path, v_reference: float,
                        closure_tolerance: float) -> tuple[dict, dict]:
	"""Shared solver-mesh gates: closure vs a reference volume, inversion,
	and the exported-skin match (guards the silent-no-surface-charge class).
	Returns (verification_error_or_None, gates, metrics); min quality is
	reported, never gated.  The verification guards (inspection actually
	ran, mesh volume finite) are part of the bf8ab4c0b hardening."""
	from ..gmsh.msh_inspect import mesh_quality, mesh_total_volume
	vol_report = mesh_total_volume(msh)
	quality = mesh_quality(msh)
	if not vol_report.get("ran") or not vol_report.get("ok"):
		return ({**_error_payload(
			"verification", "gmsh volume integration failed"),
			"volume_report": vol_report}, {}, {})
	if not quality.get("ran") or not quality.get("applicable"):
		return ({**_error_payload(
			"verification", "gmsh mesh-quality inspection failed"),
			"quality_report": quality}, {}, {})
	v_mesh = float(vol_report.get("total_volume") or 0.0)
	if not _math.isfinite(v_mesh) or v_mesh <= 0.0:
		return ({**_error_payload(
			"verification", "mesh volume is non-finite or non-positive"),
			"volume_report": vol_report}, {}, {})
	closure = abs(v_mesh - v_reference) / v_reference
	negative = int(quality.get("total_negative") or 0)
	min_det = vol_report.get("min_jacobian_det")
	surf = _vol_surface_element_stats(vol)
	vol_boundary = surf["total"]
	skin_faces = (quality.get("mesh_stats") or {}).get("n_boundary_faces")
	gates = {
		"closure_ok": bool(closure <= closure_tolerance),
		"no_inverted_elements": negative == 0 and (min_det or 0.0) > 0.0,
		# The OUTER faces must reproduce the topological skin exactly (this
		# is the silent-no-surface-charge guard).  A multi-material mesh
		# additionally carries its material-interface faces, so the raw
		# total legitimately exceeds the skin -- compare the split, not the
		# total (measured 2026-08-10: 894 skin + 190 interface on a
		# two-material Sculpt sphere).
		"boundary_faces_ok": bool(
			surf["outer"] > 0 and skin_faces is not None
			and surf["outer"] == int(skin_faces)),
	}
	metrics = {
		"mesh_volume": v_mesh, "closure": closure,
		"closure_tolerance": closure_tolerance,
		"total_negative": negative, "min_jacobian_det": min_det,
		"vol_boundary_faces": vol_boundary,
		"vol_outer_faces": surf["outer"],
		"vol_interface_faces": surf["interface"],
		"vol_faces_classified": surf["classified"],
		"msh_skin_faces": skin_faces,
		"by_type": [
			{"element": bt.get("name"), "n": bt.get("n_elements"),
			 "min": bt.get("min_quality"), "mean": bt.get("mean_quality")}
			for bt in quality.get("by_type", [])
		],
		"dof_estimate": (quality.get("mesh_stats") or {}).get(
			"dof_estimate"),
	}
	return None, gates, metrics


def _sculpt_quality_keywords(gq_iters: int, gq_threshold: float,
                             pillow_surfaces: bool, defeature: int,
                             min_vol_cells: int, *, cli: bool) -> list:
	"""Vetted Sculpt quality knobs -> in-Cubit keywords or standalone CLI
	args.  The names are IDENTICAL in both surfaces (verified against
	`help sculpt parallel` and `sculpt.exe -h`, Coreform Cubit 2025.12);
	only the spelling differs (keyword string vs --flag value)."""
	out: list = []
	if gq_iters:
		out += (["--max_gq_iters", str(int(gq_iters)),
		         "--gq_threshold", str(float(gq_threshold))] if cli else
		        [f"max_gq_iters {int(gq_iters)}",
		         f"gq_threshold {float(gq_threshold)}"])
	if pillow_surfaces:
		out += ["--pillow_surfaces"] if cli else ["pillow_surfaces"]
	if defeature:
		out += (["--defeature", str(int(defeature))] if cli else
		        [f"defeature {int(defeature)}"])
	if min_vol_cells:
		out += (["--min_vol_cells", str(int(min_vol_cells))] if cli else
		        [f"min_vol_cells {int(min_vol_cells)}"])
	return out


@mcp.tool()
def cubit_stl_to_vol(stl_path: str,
                     scheme: str = "tet",
                     size: float = 0.0,
                     out_vol: str = "",
                     out_msh: str = "",
                     closure_tolerance: float = 0.03,
                     gq_iters: int = 0,
                     gq_threshold: float = 0.2,
                     pillow_surfaces: bool = False,
                     defeature: int = 0,
                     min_vol_cells: int = 0,
                     timeout_s: int = 900) -> str:
	"""
	Mesh a watertight STL into a solver-ready Netgen `.vol` (all-hex via
	Sculpt, or tet via tetmesh), with closure and quality gates.

	This is the MESH-side half of the topology-optimization shape
	regeneration loop.  The DESIGN-side half lives in the radia wheel
	(`radia.topopt_cad`): a per-element density becomes a nodal level set
	and then a watertight STL (grid marching-cubes route, or the Cubit
	`create tri iso` route via `write_levelset_exodus`).  This tool closes
	the loop: STL -> `import stl` -> mesh -> `export netgen` -> gates.

	Routes (headless batch Cubit, per the driving policy):
	  * scheme="hex" -- `sculpt volume all processors 1 size <s>
	    gen_sidesets 2` followed
	    by `create mesh geometry hex all feature_angle 135` and
	    `delete volume with not is_meshed` (overlay-grid all-hex; works on
	    faceted geometry with no sweepable topology, which is exactly what
	    a topology-optimized body is).  `gen_sidesets 2` is the boundary
	    contract; mesh-based geometry supplies additional ownership metadata.
	    The exporter preserves direct/free tri/quad faces and deduplicates any
	    equal geometry-owned faces.  It derives DomainIn/DomainOut from adjacent
	    block elements, so one sideset spanning disconnected materials is split
	    into valid domain-aware descriptors.  Free-mesh block labels are reported
	    as `mesh_only_materials` instead of false zero CAD-volume references, and
	    imported CAD curves with no exported BBND segment do not contribute stale
	    edge-length references.  Older exporters ignored the free faces
	    and could write a `.vol` with ZERO boundary elements, which a
	    charge-based solver consumed silently as "no surface charge".
	    Measured closure depends on Sculpt size and is reported by the gate.
	  * scheme="tet" -- `volume all scheme tetmesh` + `size` (tets mesh the
	    faceted STL surface directly, so boundary tris are exported
	    naturally).  Measured closure: ~1-2 %.  NOTE netgen's own
	    STLGeometry REJECTS marching-cubes STL ("STL Triangle degenerated",
	    0 tets, even after weld + Taubin smoothing -- measured 2026-08-08);
	    Cubit is the supported tet route for such surfaces, not netgen.

	Gates (fail -> status="gate_failed", artifacts kept for inspection):
	  * closure: |V_mesh - V_stl| / V_stl <= closure_tolerance, where
	    V_mesh is the gmsh Jacobian-integrated volume of the exported
	    `.msh` and V_stl the watertight STL volume (trimesh).
	  * inversion: gmsh minSICN `negative == 0` and min Jacobian det > 0.
	  * boundary faces: the exported `.vol` must carry the complete mesh
	    skin as surface elements (count > 0 and equal to the `.msh`
	    topological skin) -- guards the silent-no-surface-charge failure
	    above.
	  * min_quality is REPORTED, not gated (one number is one sample --
	    see cubit_netgen_quality_compare's interpretation notes).

	When `size` is 0 a size is derived as bbox_diagonal / 30.  `out_vol` /
	`out_msh` default next to the STL.  The `.msh` is always produced (it
	carries the referee's volume and quality numbers).

	Args:
	    stl_path: watertight STL (the tool re-checks watertightness).
	    scheme: "hex" (Sculpt) or "tet" (tetmesh).
	    size: element size (0 = derive from bbox).
	    out_vol/out_msh: output paths (default: beside the STL).
	    closure_tolerance: relative volume-closure gate (default 3 %,
	        matching the levelset-sculpt validation gate).
	    gq_iters/gq_threshold: Sculpt guaranteed-quality smoothing
	        iterations / minimum scaled-Jacobian target (hex scheme only;
	        0 iters = off).
	    pillow_surfaces: Sculpt pillow layers at surfaces (hex only).
	    defeature/min_vol_cells: Sculpt defeaturing stage / minimum cells
	        per disconnected component (hex only; drops topopt debris).
	"""
	p = Path(stl_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	if not p.is_file():
		return json.dumps(_error_payload("input", f"STL not found: {p}"))
	if scheme not in ("hex", "tet"):
		return json.dumps(_error_payload(
			"input", f"scheme must be 'hex' or 'tet', got {scheme!r}"))
	# Numeric contract BEFORE the (subprocess) STL inspection -- restored
	# from the bf8ab4c0b hardening, which a rebase-rescue conflict
	# resolution had clobbered (2026-08-10).
	try:
		size = float(size)
		closure_tolerance = float(closure_tolerance)
		timeout_s = float(timeout_s)
	except (TypeError, ValueError) as exc:
		return json.dumps(_error_payload(
			"input", f"numeric argument is invalid: {exc}"))
	if not _math.isfinite(size) or size < 0.0:
		return json.dumps(_error_payload(
			"input", f"size must be 0 (auto) or finite and positive, got {size}"))
	if not _math.isfinite(closure_tolerance) or closure_tolerance < 0.0:
		return json.dumps(_error_payload(
			"input", "closure_tolerance must be finite and nonnegative"))
	if not _math.isfinite(timeout_s) or timeout_s <= 0.0:
		return json.dumps(_error_payload(
			"input", "timeout_s must be finite and positive"))
	inspection = _inspect_stl(p, timeout_s=min(float(timeout_s), 120.0))
	if not inspection.get("ok"):
		kind = str(inspection.get("kind") or "input")
		stage = "environment" if kind in ("environment", "timeout") else "input"
		return json.dumps(_error_payload(
			stage, str(inspection.get("error") or "cannot inspect STL"),
			kind=kind,
			hint="pip install trimesh" if kind == "environment" else None))
	if not inspection.get("watertight", False):
		return json.dumps(_error_payload(
			"input", "STL is not watertight; a closed surface is required "
			"for a volume mesh (radia.topopt_cad.iso_stl_from_grid raises "
			"on non-watertight output -- regenerate the surface)"))
	if not inspection.get("winding_consistent", False):
		return json.dumps(_error_payload(
			"input", "STL winding is inconsistent; repair face orientation "
			"before volume meshing"))
	v_stl = float(inspection.get("volume") or 0.0)
	if not _math.isfinite(v_stl) or v_stl <= 0.0:
		return json.dumps(_error_payload(
			"input", "STL volume is non-finite or zero"))

	if not size or size <= 0.0:
		bounds = inspection.get("bounds") or []
		if len(bounds) != 2 or len(bounds[0]) != 3 or len(bounds[1]) != 3:
			return json.dumps(_error_payload(
				"input", "STL inspector did not return a 3D bounding box"))
		ext = [float(hi) - float(lo)
		       for lo, hi in zip(bounds[0], bounds[1])]
		size = float(sum(value * value for value in ext) ** 0.5 / 30.0)
		if not _math.isfinite(size) or size <= 0.0:
			return json.dumps(_error_payload(
				"input", "STL bounding box cannot define a positive mesh size"))

	# absolutize AND resolve like stl_path: the journal embeds these paths
	# while the Cubit process runs with cwd=vol.parent, so a RELATIVE out
	# path would resolve against itself (nested dir for .vol, "Cannot open
	# file" for .msh -- measured 2026-08-09 from a notebook kernel)
	def _output_path(raw: str, default: Path) -> Path:
		candidate = Path(raw) if raw else default
		if not candidate.is_absolute():
			candidate = PROJECT_ROOT / candidate
		# no .resolve(): on this NAS-mapped tree it rewrites S: into the
		# UNC form, breaking the journal-embedded path contract (the
		# 2026-08-09 relative-out-path fix owns the absolutize semantics)
		return candidate

	vol = _output_path(out_vol, p.with_suffix(f".{scheme}.vol"))
	msh = _output_path(out_msh, p.with_suffix(f".{scheme}.msh"))
	path_keys = [os.path.normcase(str(path)) for path in (p, vol, msh)]
	if len(set(path_keys)) != 3:
		return json.dumps(_error_payload(
			"input", "stl_path, out_vol, and out_msh must be distinct paths"))
	stl_fwd = str(p).replace("\\", "/")
	cmds = [f'import stl "{stl_fwd}" feature_angle 135 merge']
	if scheme == "tet":
		cmds += ["volume all scheme tetmesh",
		         f"volume all size {size}", "mesh volume all"]
	else:
		# Keep both representations: gen_sidesets is the explicit free-skin
		# contract, while mesh geometry supplies ownership metadata.  The
		# exporter deduplicates equal faces.  Drop the unmeshed import volume
		# by predicate, never by a hardcoded id.
		extra = _sculpt_quality_keywords(gq_iters, gq_threshold,
		                                 pillow_surfaces, defeature,
		                                 min_vol_cells, cli=False)
		sculpt_cmd = (f"sculpt volume all processors 1 size {size} "
		              f"gen_sidesets 2")
		if extra:
			sculpt_cmd += " " + " ".join(extra)
		cmds += [sculpt_cmd,
		         "create mesh geometry hex all feature_angle 135",
		         "delete volume with not is_meshed"]
	cmds += ["block 1 add volume all", 'block 1 name "iron"',
	         f'export netgen "{str(vol).replace(chr(92), "/")}" '
	         f"order 1 overwrite",
	         f'export gmsh "{str(msh).replace(chr(92), "/")}" '
	         f"dimension 3 order 1 overwrite"]
	for output in (vol, msh):
		output.parent.mkdir(parents=True, exist_ok=True)
		try:
			output.unlink(missing_ok=True)
		except OSError as exc:
			return json.dumps(_error_payload(
				"output", f"cannot replace stale output {output}: {exc}"))
	r = _cs.run_headless_journal(
		cmds, timeout_s=timeout_s, working_directory=vol.parent)
	if r.get("status") == "error":
		return json.dumps({**_error_payload("cubit", "batch meshing failed"),
		                   "process": r,
		                   "commands": cmds}, default=str)
	if not vol.is_file() or not msh.is_file():
		return json.dumps({**_error_payload(
			"cubit", f"export did not produce {vol.name} / {msh.name} "
			"(check process diagnostics)"),
			"process": r}, default=str)

	verification_error, gates, metrics = _collect_mesh_gates(
		msh, vol, v_stl, closure_tolerance)
	if verification_error is not None:
		return json.dumps(verification_error, default=str)
	return json.dumps({
		"status": "ok" if all(gates.values()) else "gate_failed",
		"gates": gates,
		"scheme": scheme, "size": size,
		"stl": str(p), "vol": str(vol), "msh": str(msh),
		"stl_volume": v_stl,
		**metrics,
		"process": {
			"exit_code": r.get("exit_code"),
			"console": r.get("console"),
			"headless_flags": r.get("headless_flags"),
			"persistent_gui_started": r.get("persistent_gui_started"),
		},
		"note": ("min quality is reported, not gated; gate on inverted "
		         "elements and closure (mesh-quality study finding)"),
	}, ensure_ascii=False, indent=2, default=str)


@mcp.tool()
def cubit_vfrac_to_vol(vfrac_path: str,
                       out_vol: str = "",
                       out_msh: str = "",
                       material_names: str = "",
                       closure_tolerance: float = 0.01,
                       smooth_method: int = 0,
                       gq_iters: int = 0,
                       gq_threshold: float = 0.2,
                       pillow_surfaces: bool = False,
                       defeature: int = 0,
                       min_vol_cells: int = 0,
                       timeout_s: int = 900) -> str:
	"""Mesh a volume-fraction Exodus (``radia.topopt_cad.write_vfrac_exodus``)
	into a solver-ready all-hex Netgen `.vol` via standalone Sculpt.

	Sculpt is natively a volume-fraction mesher; this route feeds the
	design occupancy DIRECTLY (``sculpt.exe --input_vfrac``), skipping
	marching cubes and STL entirely.  Sculpt's interface reconstruction
	places the boundary sub-cell from the fractions: measured on a voxel
	sphere, the meshed material volume tracked the input fractions to
	~1e-4 relative (0/1 SPN input of the same body: -8 % smoothing
	shrink; marching-cubes STL route: percent-level closure).  The
	default closure gate is therefore 1 %, tighter than the STL route's
	3 %.

	Pipeline: read + validate the vfrac Exodus (single MAT_1 + VOID
	contract; the closure reference is the fraction integral) ->
	``sculpt.exe -ivf`` (exodus in PHYSICAL coordinates; ``-SS 2``
	sidesets; optional quality knobs) -> ``import mesh geometry`` into
	headless Cubit -> block named ``iron`` -> ``export netgen`` +
	``export gmsh`` -> the same gates as ``cubit_stl_to_vol`` (closure
	vs the vfrac integral, inversion, exported-skin match).

	MULTI-MATERIAL: a vfrac Exodus carrying ``MAT_1..MAT_N`` produces one
	Sculpt block per material with a CONFORMAL hex interface between them
	(shared nodes -- an STL route would have to boolean two surfaces).
	``material_names`` labels them in ``MAT_<id>`` order; the labels ride
	through ``export netgen`` into the ``.vol`` material names, which is
	what per-region ``vim.Solve(mesh, mu_r={name: value})`` and
	``hdiv_vim_demag_eval(mu_r="core=1000,shell=50")`` consume.

	Args:
	    vfrac_path: the ``.e.1.0`` file or its base name.
	    out_vol/out_msh: output paths (default: beside the input).
	    material_names: comma-separated block labels in ``MAT_<id>`` order
	        (default: "iron" for one material, ``mat_1..mat_N`` otherwise).
	    closure_tolerance: relative volume gate vs the vfrac integral.
	    smooth_method: Sculpt ``--smooth`` override (0 = Sculpt default;
	        8 projects to the interpolated interface surface per the
	        Coreform docs -- try it when boundary fidelity matters more
	        than smoothness).
	    gq_iters/gq_threshold: guaranteed-quality smoothing iterations /
	        minimum scaled-Jacobian target (0 iters = off).
	    pillow_surfaces: insert pillow layers at material surfaces.
	    defeature/min_vol_cells: automatic defeaturing stage / minimum
	        cells per disconnected component (drops topopt debris).
	    timeout_s: per-stage subprocess timeout.
	"""
	import subprocess as sp

	p = Path(vfrac_path)
	if not p.is_absolute():
		p = PROJECT_ROOT / p
	if not p.name.endswith(".e.1.0"):
		p = p.with_name(p.name + ".e.1.0")
	if not p.is_file():
		return json.dumps(_error_payload(
			"input", f"vfrac Exodus not found: {p} (write it with "
			"radia.topopt_cad.write_vfrac_exodus)"))
	try:
		import netCDF4 as nc
	except ImportError:
		return json.dumps(_error_payload(
			"environment", "netCDF4 is required to validate the vfrac "
			"Exodus; install it with `pip install netCDF4`"))
	try:
		ds = nc.Dataset(str(p))
		try:
			names = [str(s).strip() for s in
			         nc.chartostring(ds.variables["name_elem_var"][:])]
			expected = ["VOID"] + [f"MAT_{i + 1}"
			                       for i in range(len(names) - 1)]
			if len(names) < 2 or names != expected:
				return json.dumps(_error_payload(
					"input", f"vfrac Exodus element variables {names} != "
					f"{expected} -- this tool speaks the "
					"write_vfrac_exodus VOID + MAT_1..MAT_N contract"))
			n_mat = len(names) - 1
			gnames = [str(s).strip() for s in
			          nc.chartostring(ds.variables["name_glo_var"][:])]
			gvals = {k: float(v) for k, v in
			         zip(gnames, ds.variables["vals_glo_var"][0])}
			nx, ny, nz = (int(gvals["gxint"]), int(gvals["gyint"]),
			              int(gvals["gzint"]))
			hx = (gvals["xmax"] - gvals["xmin"]) / nx
			hy = (gvals["ymax"] - gvals["ymin"]) / ny
			hz = (gvals["zmax"] - gvals["zmin"]) / nz
			cell = hx * hy * hz
			per_material = [
				{"id": index + 1,
				 "v_vfrac": float(
					 ds.variables[f"vals_elem_var{index + 2}eb1"][0].sum())
				 * cell}
				for index in range(n_mat)]
			v_vfrac = float(sum(m["v_vfrac"] for m in per_material))
		finally:
			ds.close()
	except (KeyError, OSError, ValueError) as exc:
		return json.dumps(_error_payload(
			"input", f"cannot validate vfrac Exodus {p.name}: {exc}"))
	if v_vfrac <= 0.0:
		return json.dumps(_error_payload(
			"input", "vfrac Exodus carries zero material volume"))
	labels = [s.strip() for s in str(material_names).split(",") if s.strip()]
	if labels and len(labels) != n_mat:
		return json.dumps(_error_payload(
			"input", f"material_names has {len(labels)} entries but the "
			f"vfrac Exodus carries {n_mat} material(s)"))
	if not labels:
		labels = ["iron"] if n_mat == 1 else [f"mat_{i + 1}"
		                                      for i in range(n_mat)]
	if len(set(labels)) != len(labels):
		return json.dumps(_error_payload(
			"input", f"material_names must be distinct, got {labels}"))
	for index, name in enumerate(labels):
		per_material[index]["name"] = name

	from .session import get_cubit_bin_dir
	bin_dir = get_cubit_bin_dir()
	sculpt_exe = (Path(bin_dir) / "sculpt.exe") if bin_dir else None
	if sculpt_exe is None or not sculpt_exe.is_file():
		return json.dumps(_error_payload(
			"environment", "sculpt.exe not found in the Cubit bin "
			f"directory ({bin_dir}); Sculpt ships with Coreform Cubit "
			"2025.12+ on Windows"))

	base = p.with_name(p.name[:-len(".e.1.0")])
	vol = Path(out_vol) if out_vol else base.with_suffix(".hex.vol")
	msh = Path(out_msh) if out_msh else base.with_suffix(".hex.msh")
	if not vol.is_absolute():
		vol = PROJECT_ROOT / vol
	if not msh.is_absolute():
		msh = PROJECT_ROOT / msh
	sculpt_out = vol.with_suffix("").with_name(vol.stem + "_sculpt")
	sculpt_exo = sculpt_out.with_name(sculpt_out.name + ".e.1.0")
	for output in (vol, msh, sculpt_exo):
		output.parent.mkdir(parents=True, exist_ok=True)
		try:
			output.unlink(missing_ok=True)
		except OSError as exc:
			return json.dumps(_error_payload(
				"output", f"cannot replace stale output {output}: {exc}"))

	cli = [str(sculpt_exe), "-ivf", str(base), "-e", str(sculpt_out),
	       "-SS", "2", "-cv"]
	for entry in per_material:
		cli += ["-mn", str(entry["id"]), str(entry["name"])]
	if smooth_method:
		cli += ["--smooth", str(int(smooth_method))]
	cli += _sculpt_quality_keywords(gq_iters, gq_threshold,
	                                pillow_surfaces, defeature,
	                                min_vol_cells, cli=True)
	try:
		run = sp.run(cli, capture_output=True, text=True,
		             timeout=timeout_s, cwd=str(vol.parent))
	except sp.TimeoutExpired:
		return json.dumps(_error_payload(
			"cubit", f"sculpt.exe timed out after {timeout_s}s"))
	sculpt_tail = "\n".join((run.stdout + run.stderr).splitlines()[-25:])
	if run.returncode != 0 or not sculpt_exo.is_file():
		return json.dumps({**_error_payload(
			"cubit", "sculpt.exe failed or produced no Exodus mesh"),
			"exit_code": run.returncode, "command": cli,
			"console_tail": sculpt_tail}, default=str)

	exo_fwd = str(sculpt_exo).replace("\\", "/")
	# Sculpt's --material_name already labels the Exodus blocks, so the
	# names arrive with the import; re-assert them by block id so the
	# `.vol` materials are the requested ones even on an older Sculpt.
	cmds = [f'import mesh geometry "{exo_fwd}" feature_angle 135 merge']
	cmds += [f'block {entry["id"]} name "{entry["name"]}"'
	         for entry in per_material]
	cmds += [f'export netgen "{str(vol).replace(chr(92), "/")}" '
	         f"order 1 overwrite",
	         f'export gmsh "{str(msh).replace(chr(92), "/")}" '
	         f"dimension 3 order 1 overwrite"]
	r = _cs.run_headless_journal(
		cmds, timeout_s=timeout_s, working_directory=vol.parent)
	if r.get("status") == "error":
		return json.dumps({**_error_payload(
			"cubit", "exodus import / export failed"),
			"process": r, "commands": cmds}, default=str)
	if not vol.is_file() or not msh.is_file():
		return json.dumps({**_error_payload(
			"cubit", f"export did not produce {vol.name} / {msh.name}"),
			"process": r}, default=str)

	verification_error, gates, metrics = _collect_mesh_gates(
		msh, vol, v_vfrac, closure_tolerance)
	if verification_error is not None:
		return json.dumps(verification_error, default=str)
	return json.dumps({
		"status": "ok" if all(gates.values()) else "gate_failed",
		"gates": gates,
		"vfrac": str(p), "vol": str(vol), "msh": str(msh),
		"sculpt_exodus": str(sculpt_exo),
		"vfrac_volume": v_vfrac,
		"materials": per_material,
		"grid": {"nel": [nx, ny, nz], "cell": [hx, hy, hz]},
		**metrics,
		"sculpt": {"exit_code": run.returncode, "command": cli,
		           "console_tail": sculpt_tail},
		"note": ("closure is gated against the volume-fraction integral; "
		         "min quality is reported, not gated"),
	}, ensure_ascii=False, indent=2, default=str)


@mcp.tool()
def cubit_batch_try(commands: list, step_path: str = "",
                     timeout_s: int = 300) -> str:
	"""Dry-run a recipe in a fresh headless Cubit subprocess.

	Launches `coreform_cubit.exe` in batch / nographics / nojournal
	mode, optionally imports `step_path`, runs `commands` in order,
	probes final geometry/mesh counts, and tears down. Completely
	isolated from the live GUI session (`cubit_exec`).

	Intended flow:
	  1. AI proposes a recipe (e.g., `volume all scheme sweep` + `mesh ...`)
	  2. AI calls `cubit_batch_try(step, recipe)` to verify it actually
	     produces the expected element counts.
	  3. If it works, AI calls `cubit_exec(recipe)` to replay in the live
	     GUI window so the user can see the result.

	Args:
	    commands: list of Cubit command strings to run after import.
	    step_path: optional STEP file to import first. Empty = no import.
	    timeout_s: per-call RPC timeout.

	Returns JSON with status / per_line / summary (volumes, hexes, tets,
	nodes, ...) / failed_at / error.
	"""
	if not commands:
		return json.dumps({"status": "error",
		                   "error": "commands list cannot be empty"})
	cmd_list = [str(c).rstrip() for c in commands]
	result = _run_batch(step_path or None, cmd_list, timeout_s=timeout_s)
	if result.get("status") == "error":
		_fl.record_failure("cubit", {
			"input": cmd_list,
			"step_path": step_path,
			"stage": result.get("stage", "batch"),
			"error": result.get("error") or str(result)[:400],
			"context": {"summary": result.get("summary", {}),
			            "mode": "batch"},
		})
	return json.dumps(result, indent=2)


_MESH_LADDER: list[dict] = [
	{
		"name": "auto",
		"cmds": ["volume all scheme auto"],
		"why": "Let Cubit pick a per-volume scheme. Works for simple prismatic shapes.",
	},
	{
		"name": "sweep",
		"cmds": ["volume all scheme sweep"],
		"why": "Force sweep: pave source / sweep to target. Best hex yield when applicable.",
	},
	{
		"name": "webcut_cyl_auto",
		"cmds": [
			# Cylinder webcut along +Z, outside the bbox — splits any
			# compound-sweep body whose side surface wraps around Z.
			# Heuristic: radius = 0.9 × max(|x|,|y|) bbox extent. Supplied
			# by the preamble callback; see `_mesh_ladder_with_geom_split`.
		],
		"why": ("Compound swept body detected (too many surfaces per "
		         "volume). Cylinder webcut on Z-axis splits it into "
		         "simpler prisms before retrying scheme auto."),
		"needs_geom_split": True,
	},
	{
		"name": "polyhedron",
		"cmds": ["volume all scheme polyhedron"],
		"why": "Hex-dominant scheme tolerant of non-prismatic topology.",
	},
	{
		"name": "tetmesh",
		"cmds": ["volume all scheme tetmesh"],
		"why": "Always works; fallback giving tets rather than hex.",
	},
]


def _geom_split_recipe(summary: dict) -> list[str]:
	"""Build a cylinder webcut + imprint/merge + scheme-auto recipe
	for a detected compound-sweep body.

	Uses the summary's bbox-free info. Since we don't have direct bbox
	access via `probe summary`, we issue a Cubit command that wraps:
	  1. `list volume all` — forces Cubit to emit bbox info (we don't
	     use that here but keep as diagnostic).
	  2. `webcut volume all with cylinder radius <R> axis z` where R
	     is chosen as a safe-ish default (16 mm), which happens to
	     match the spiral-coil case. For truly unknown geometry the
	     caller should pass their own recipe.

	The auto radius heuristic is deliberately conservative: caller can
	always invoke `cubit_exec_safely` with a custom webcut, or provide
	`_MESH_LADDER_RADIUS` via the batch env in a future refinement.
	"""
	del summary  # radius heuristic is geometry-agnostic for now
	radius = float(os.environ.get("RADIA_MCP_WEBCUT_RADIUS", "16.01"))
	return [
		f"webcut volume all with cylinder radius {radius} axis z",
		"imprint all",
		"merge all",
		"volume all scheme auto",
	]


@mcp.tool()
def cubit_mesh_auto(step_path: str = "",
                    target_size: float = 1.0,
                    prefer: str = "hex",
                    commit_to_gui: bool = True,
                    timeout_s: int = 600) -> str:
	"""Find a working mesh recipe by trying a scheme ladder in batch,
	then optionally replay the winning recipe in the live GUI session.

	Ladder:
	  1. scheme auto       (hope for the best)
	  2. scheme sweep      (force sweep; pure hex if topology fits)
	  3. scheme polyhedron (hex-dominant for complex shapes)
	  4. scheme tetmesh    (guaranteed fallback, tet only)

	Each rung is executed in a fresh headless Cubit (no GUI pollution,
	no state entanglement with the live session). The first rung that
	produces >0 elements with the preferred family is the winner.
	When `commit_to_gui=True`, the winning recipe is replayed via the
	live GUI `CubitSession` so the user sees the mesh build live.

	Args:
	    step_path: STEP file to import in each batch. Empty → caller
	        has already imported into the GUI session, recipe is
	        validated against THAT state by rewinding with delete mesh
	        + new scheme. If step_path empty, `commit_to_gui` also
	        skips the re-import (just delete mesh + run winning recipe).
	    target_size: mesh size command target. Passed as
	        `volume all size <target_size>` before each scheme.
	    prefer: "hex" (require hex>0 to accept) or "any" (any element).
	    commit_to_gui: replay winning recipe in the live GUI session.
	    timeout_s: per-ladder-rung timeout.

	Returns structured report: which schemes tried, element counts per
	attempt, chosen recipe, and what the GUI session looks like after
	replay (if commit_to_gui).
	"""
	prefer = (prefer or "hex").strip().lower()
	if prefer not in ("hex", "any"):
		return json.dumps({"status": "error",
		                   "error": f"prefer must be 'hex' or 'any', got {prefer!r}"})

	attempts: list[dict] = []
	winner: dict | None = None
	# Seed summary from the batch-side for compound detection. Runs a
	# minimal import + probe first so the geom-split rung can adapt.
	seed = _run_batch(step_path or None, [], timeout_s=timeout_s)
	seed_summary = seed.get("summary", {}) or {}
	vol = int(seed_summary.get("volumes", 0) or 0)
	surf = int(seed_summary.get("surfaces", 0) or 0)
	# Compound heuristic: low volume count with lots of surfaces suggests
	# a single swept body whose path has kinks (e.g., lead→helix→lead).
	compound_suspected = bool(vol and vol <= 3 and (surf / max(1, vol)) >= 7)

	for rung in _MESH_LADDER:
		if rung.get("needs_geom_split"):
			if not compound_suspected:
				continue  # skip geom-split on clean topology
			cmds = _geom_split_recipe(seed_summary)
		else:
			cmds = list(rung["cmds"])
		recipe = [f"volume all size {float(target_size)}",
		          *cmds,
		          "mesh volume all"]
		r = _run_batch(step_path or None, recipe, timeout_s=timeout_s)
		smy = r.get("summary", {}) or {}
		h = int(smy.get("hexes", 0))
		t = int(smy.get("tets", 0))
		n = int(smy.get("nodes", 0))
		ok_elems = (h > 0) if prefer == "hex" else ((h + t) > 0)
		attempts.append({
			"scheme": rung["name"],
			"why": rung["why"],
			"status": r.get("status"),
			"hexes": h, "tets": t, "nodes": n,
			"failed_at": r.get("failed_at"),
			"error": r.get("error"),
			"geom_split_applied": bool(rung.get("needs_geom_split")),
		})
		if r.get("status") == "ok" and ok_elems:
			winner = {"scheme": rung["name"], "recipe": recipe,
			          "hexes": h, "tets": t, "nodes": n,
			          "geom_split_applied": bool(rung.get("needs_geom_split"))}
			break

	if winner is None:
		return json.dumps({
			"status": "error",
			"reason": "no scheme in the ladder produced the preferred element type",
			"prefer": prefer,
			"attempts": attempts,
			"suggest": ("Try geometric splits before meshing — e.g., "
			            "webcut cylinder axis z radius R, webcut with plane xplane, "
			            "or split helix/lead transitions. Consult "
			            "cubit_lookup(query='sweep multi source target') and "
			            "cubit_web_docs(query='sweep complex geometry', source='cubit_forum')."),
		}, indent=2)

	# Replay the winning recipe in the live GUI session
	gui_report: dict = {"skipped": True}
	if commit_to_gui:
		sess, err = _cubit_session_or_error()
		if err is not None:
			gui_report = {"status": "error", "note": "no live GUI session; caller must replay manually",
			              "detail": json.loads(err)}
		else:
			# If caller provided a step_path, import it so GUI mirrors batch state.
			gui_cmds: list[str] = []
			if step_path:
				step_fwd = str(step_path).replace("\\", "/")
				gui_cmds.append(f'reset')
				gui_cmds.append(f'import step "{step_fwd}"')
			gui_cmds.append("delete mesh")
			gui_cmds.extend(winner["recipe"])
			try:
				rr = sess.call("cmd", gui_cmds, timeout_s=timeout_s)
			except _cs.CubitSessionError as e:
				gui_report = {"status": "error", "stage": "gui_rpc", "error": str(e)}
			else:
				per_line = rr.get("result", [])
				try:
					smy = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
				except _cs.CubitSessionError:
					smy = {}
				gui_report = {
					"status": "ok" if all(s.get("ok") for s in per_line) else "error",
					"summary": smy,
					"replayed": gui_cmds,
				}

	return json.dumps({
		"status": "ok",
		"winner": winner,
		"attempts": attempts,
		"gui": gui_report,
	}, indent=2)


# ============================================================
# Unified retrieval: cubit_ask (bundled KB + examples + web docs + failures)
# ============================================================

@mcp.tool()
def cubit_ask(query: str, limit: int = 6,
              include_web: bool = False) -> str:
	"""One-shot search across every Cubit knowledge surface we have.

	Unions the result sets of:
	  - `cubit_lookup`  (bundled KB: scripting / API / forum tips /
	    netgen workflow / export rules + persistent failure log)
	  - `cubit_examples` (forum Discourse + local lab archives)
	  - optionally `cubit_web_docs` (live forum via Discourse JSON, if
	    `include_web=True`)

	Each hit carries a `layer` field so the caller sees whether it
	came from bundled KB, scraped examples, or live web. Scores are
	tf-idf-compatible across layers (same tokenizer / idf / heading
	boost), so a single ranked list is meaningful.

	Use this when you don't know where the answer lives — "how do I
	mesh a helix coil?" / "export netgen order 3" / "Cubit
	crashes on imprint".

	Args:
	    query: keywords.
	    limit: cap on returned hits across all layers.
	    include_web: if True, also hit `cubit_web_docs(source='cubit_forum')`
	        (adds a live HTTP call — default off).
	"""
	q = (query or "").strip()
	if not q:
		return json.dumps({"status": "error", "kind": "input", "error": "empty query"})

	results: list[dict] = []

	# Layer 1: bundled KB via cubit_lookup
	try:
		r = json.loads(cubit_lookup(query=q, max_results=limit))
		for h in r.get("results", []):
			h["layer"] = "kb"
			results.append(h)
	except Exception as e:  # noqa: BLE001
		results.append({"layer": "kb", "error": str(e)})

	# Layer 2: scraped examples via cubit_examples (unions forum + local)
	try:
		from ..common import examples as _ex
		ex_r = _ex.search_examples("cubit", q, limit=limit,
		                            auto_refresh_if_empty=False)
		for h in ex_r.get("results", []):
			h["layer"] = "examples"
			results.append(h)
	except Exception as e:  # noqa: BLE001
		results.append({"layer": "examples", "error": str(e)})

	# Layer 3: optional live web (Discourse JSON)
	if include_web:
		try:
			web = json.loads(cubit_web_docs(query=q,
			                                 source="cubit_forum",
			                                 max_hits=3))
			for h in web.get("hits", []):
				h["layer"] = "web"
				results.append(h)
		except Exception as e:  # noqa: BLE001
			results.append({"layer": "web", "error": str(e)})

	# Re-rank across layers by score (falling back to 0 for hits that
	# don't expose one — they'll land below tf-idf results)
	def _score_of(h: dict) -> float:
		try:
			return float(h.get("score") or 0)
		except (TypeError, ValueError):
			return 0.0
	results.sort(key=lambda h: -_score_of(h))
	return json.dumps({
		"status": "ok",
		"query": q,
		"total": len(results),
		"results": results[: max(1, int(limit))],
	}, indent=2, ensure_ascii=False)


# ============================================================
# Parallel race: try N recipes concurrently in batch, first/best wins
# ============================================================

@mcp.tool()
def cubit_mesh_race(step_path: str,
                     recipes: list,
                     prefer: str = "hex",
                     max_concurrent: int = 4,
                     commit_to_gui: bool = True,
                     timeout_s: int = 600) -> str:
	"""Race N recipes in parallel batch Cubits, replay the first/best
	winner in the live GUI.

	Each `recipes[i]` is a dict:
	    {"name": "<label>", "cmds": ["volume all size 0.5", "..."],
	     "size": <unused, included for AI authoring convenience>}

	All N batches start from the same `step_path`. Each lives in its
	own subprocess so a Cubit crash on one variant doesn't kill the
	others. The first variant to satisfy the `prefer` predicate
	(hex>0 if `prefer='hex'`, any element if `prefer='any'`) is the
	winner; remaining processes are not cancelled (Cubit batch lives
	~30 s anyway), but their results aren't replayed. After the
	winner is chosen, replays it on the live GUI session.

	Use this when you want to parameter-sweep — e.g., "try fillet
	radii 0.2 / 0.5 / 1.0 in parallel, ship whichever meshes
	cleanest." Compare to `cubit_mesh_auto` which is the SERIAL
	scheme ladder.

	Args:
	    step_path: STEP file imported into each batch.
	    recipes: list of {name, cmds} dicts.
	    prefer: "hex" or "any".
	    max_concurrent: pool size (default 4 — safe for 16 GB RAM).
	    commit_to_gui: replay winner in live GUI.
	    timeout_s: per-batch timeout.
	"""
	import concurrent.futures as _cf
	if not recipes:
		return json.dumps({"status": "error",
		                   "error": "recipes list cannot be empty"})
	prefer = (prefer or "hex").strip().lower()

	def _one(rec: dict) -> dict:
		name = rec.get("name", "?")
		cmds = list(rec.get("cmds") or [])
		if not cmds:
			return {"name": name, "status": "error",
			        "error": "empty cmds"}
		r = _run_batch(step_path or None, cmds, timeout_s=timeout_s)
		smy = r.get("summary", {}) or {}
		return {
			"name": name,
			"status": r.get("status"),
			"hexes": int(smy.get("hexes", 0) or 0),
			"tets": int(smy.get("tets", 0) or 0),
			"nodes": int(smy.get("nodes", 0) or 0),
			"failed_at": r.get("failed_at"),
			"error": r.get("error"),
			"recipe": cmds,
		}

	results: list[dict] = []
	winner: dict | None = None
	with _cf.ThreadPoolExecutor(max_workers=max(1, int(max_concurrent))) as pool:
		futures = {pool.submit(_one, rec): rec for rec in recipes}
		for fut in _cf.as_completed(futures):
			res = fut.result()
			results.append(res)
			ok = (res["hexes"] > 0) if prefer == "hex" \
				else ((res["hexes"] + res["tets"]) > 0)
			if winner is None and res.get("status") == "ok" and ok:
				winner = res

	# Sort results by element count for diagnostic clarity
	results.sort(key=lambda r: -(r.get("hexes", 0) + r.get("tets", 0)))

	if winner is None:
		return json.dumps({
			"status": "error",
			"reason": "no recipe produced the preferred element type",
			"prefer": prefer,
			"results": results,
		}, indent=2)

	gui_report: dict = {"skipped": True}
	if commit_to_gui:
		sess, err = _cubit_session_or_error()
		if err is not None:
			gui_report = {"status": "error", "note": "no live GUI session",
			              "detail": json.loads(err)}
		else:
			step_fwd = str(step_path).replace("\\", "/") if step_path else None
			gui_cmds = []
			if step_fwd:
				gui_cmds.append("reset")
				gui_cmds.append(f'import step "{step_fwd}"')
			gui_cmds.append("delete mesh")
			gui_cmds.extend(winner["recipe"])
			try:
				rr = sess.call("cmd", gui_cmds, timeout_s=timeout_s)
			except _cs.CubitSessionError as e:
				gui_report = {"status": "error", "stage": "gui_rpc",
				              "error": str(e)}
			else:
				per_line = rr.get("result", [])
				try:
					smy = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
				except _cs.CubitSessionError:
					smy = {}
				gui_report = {
					"status": "ok" if all(s.get("ok") for s in per_line) else "error",
					"summary": smy,
					"replayed": gui_cmds,
				}

	return json.dumps({
		"status": "ok",
		"winner": winner,
		"results": results,
		"max_concurrent": max_concurrent,
		"gui": gui_report,
	}, indent=2)


# ============================================================
# Race review: shortlist all candidates with quality, human picks
# ============================================================

def _race_history_dir() -> Path:
	d = _fl.state_dir() / "race_history"
	d.mkdir(parents=True, exist_ok=True)
	return d


def _learned_recipes_path() -> Path:
	"""Where the learned-recipes jsonl lives.

	Resolution order:
	  1. `RADIA_MCP_LEARNED_DIR` env var — point all lab machines to a
	     shared SMB / network drive (e.g.,
	     `repo:/lab_learned\\`) so every race winner
	     contributes to one collective pool.
	  2. Fallback: `<state_dir>/learned_recipes.jsonl` (per-machine).
	"""
	env = os.environ.get("RADIA_MCP_LEARNED_DIR")
	if env:
		base = Path(env)
		base.mkdir(parents=True, exist_ok=True)
		return base / "learned_recipes.jsonl"
	return _fl.state_dir() / "learned_recipes.jsonl"


def _bundled_curated_recipes() -> list[dict]:
	"""Load the bundled, lab-curated recipe pool that ships in the
	wheel. Generated periodically from accumulated learned_recipes.jsonl
	via `cubit_curate_learned_recipes`. Empty list if the bundle has
	not been built yet.
	"""
	try:
		from .curated_recipes_bundle import CURATED  # type: ignore
		return list(CURATED) if isinstance(CURATED, list) else []
	except Exception:
		return []


def _geometry_signature(summary: dict) -> dict:
	"""Compact, comparable signature of the live geometry pre-race.

	Used to look up "winning recipe for similar shape" from prior
	races. Two geometries are 'similar' when their volume count
	matches and their (surfaces / volumes) ratio is within ±20%.
	"""
	vol = int(summary.get("volumes", 0) or 0)
	surf = int(summary.get("surfaces", 0) or 0)
	curv = int(summary.get("curves", 0) or 0)
	vert = int(summary.get("vertices", 0) or 0)
	return {
		"volumes": vol,
		"surfaces": surf,
		"surf_per_vol": round(surf / max(1, vol), 2),
		"curves": curv,
		"vertices": vert,
	}


def _signature_similar(a: dict, b: dict) -> bool:
	"""Same volume count, surf/vol ratio within ±20%."""
	if a.get("volumes") != b.get("volumes"):
		return False
	a_r = float(a.get("surf_per_vol", 0) or 0)
	b_r = float(b.get("surf_per_vol", 0) or 0)
	if max(a_r, b_r) == 0:
		return True
	return abs(a_r - b_r) / max(a_r, b_r, 1e-9) < 0.2


def _record_learned_recipe(race_id: str, sig: dict,
                              recipe: list[str], quality: dict,
                              source: str = "human") -> None:
	"""Append a winning recipe to `learned_recipes.jsonl` for future
	smart-recipe generation. `source` is 'human' or 'ai_<name>'.
	"""
	import time as _t
	rec = {
		"when": _t.time(),
		"race_id": race_id,
		"source": source,
		"signature": sig,
		"recipe": recipe,
		"quality": quality,
	}
	try:
		path = _learned_recipes_path()
		with path.open("a", encoding="utf-8") as f:
			f.write(json.dumps(rec, ensure_ascii=False) + "\n")
	except OSError:
		pass


def _load_learned_recipes_for(sig: dict, max_n: int = 3) -> list[dict]:
	"""Return up to `max_n` past winning recipes for similar geometry,
	**unioned across local, lab-shared, and bundled-curated pools**,
	sorted by quality (min-Jacobian descending).

	Sources (all unioned):
	  1. `<state_dir>/learned_recipes.jsonl` — per-machine.
	  2. `RADIA_MCP_LEARNED_DIR/learned_recipes.jsonl` — lab-shared
	     pool when the env var is set (Stage 1: collective intelligence).
	  3. Bundled `curated_recipes_bundle.CURATED` — ships in the wheel,
	     contains lab-distilled top picks across all known races
	     (Stage 2: world-wide knowledge transfer via PyPI release).
	"""
	hits: list[dict] = []

	# Sources 1 + 2 (jsonl pools)
	path = _learned_recipes_path()
	if path.exists():
		try:
			with path.open("r", encoding="utf-8") as f:
				lines = f.readlines()
		except OSError:
			lines = []
		for line in lines:
			line = line.strip()
			if not line:
				continue
			try:
				rec = json.loads(line)
			except json.JSONDecodeError:
				continue
			if _signature_similar(rec.get("signature", {}), sig):
				hits.append(rec)

	# Source 3 (bundled curated)
	for rec in _bundled_curated_recipes():
		if _signature_similar(rec.get("signature", {}), sig):
			rec = dict(rec)
			rec["source"] = rec.get("source", "curated")
			hits.append(rec)

	hits.sort(
		key=lambda r: -(((r.get("quality") or {}).get("min")) or -1.0),
	)
	return hits[: max(1, int(max_n))]


@mcp.tool()
def cubit_curate_learned_recipes(out_module_path: str = "",
                                    top_per_class: int = 3,
                                    min_quality_jacobian: float = 0.3) -> str:
	"""**Lab maintainer tool**: read accumulated `learned_recipes.jsonl`,
	dedup + group by signature class, pick top-N per class by quality,
	emit a Python module that ships in the next radia-mcp wheel.

	Stage 2 of the lab-collective-intelligence pipeline:
	  - Stage 1: lab machines share `RADIA_MCP_LEARNED_DIR` jsonl.
	  - **Stage 2**: this tool distills the pool into bundled
	    `curated_recipes_bundle.py` (typically run weekly).
	  - Stage 3: the wheel + a CC-BY docs/design page publishes the
	    distilled recipes to the wider CAE community.

	Args:
	    out_module_path: where to write the .py module. Empty = the
	        package's own `radia_mcp/cubit/curated_recipes_bundle.py`
	        (so the next `python -m build` picks it up).
	    top_per_class: cap on recipes per signature class (default 3).
	    min_quality_jacobian: drop recipes whose `min` scaled-Jacobian
	        is below this threshold (default 0.3).
	"""
	import time as _t

	pool_path = _learned_recipes_path()
	if not pool_path.exists():
		return json.dumps({"status": "error",
		                   "error": f"no learned recipes at {pool_path}"})

	try:
		raw_lines = pool_path.read_text(encoding="utf-8").splitlines()
	except OSError as e:
		return json.dumps({"status": "error",
		                   "stage": "read_pool",
		                   "error": str(e)})

	records: list[dict] = []
	for line in raw_lines:
		line = line.strip()
		if not line:
			continue
		try:
			rec = json.loads(line)
		except json.JSONDecodeError:
			continue
		# Quality filter
		q = (rec.get("quality") or {}).get("min")
		if q is None or float(q) < min_quality_jacobian:
			continue
		records.append(rec)

	# Group by (volumes, surf_per_vol bucket) — same key both in
	# similarity check and curation
	def _bucket(sig: dict) -> tuple:
		vol = int(sig.get("volumes", 0) or 0)
		ratio = round(float(sig.get("surf_per_vol", 0) or 0), 1)
		return (vol, ratio)

	classes: dict = {}
	for rec in records:
		key = _bucket(rec.get("signature", {}))
		classes.setdefault(key, []).append(rec)

	curated: list[dict] = []
	for key, recs in classes.items():
		# Dedup by recipe text
		seen: set[str] = set()
		uniq: list[dict] = []
		for rec in recs:
			recipe_key = "\n".join(rec.get("recipe", []))
			if recipe_key in seen:
				continue
			seen.add(recipe_key)
			uniq.append(rec)
		# Sort by quality desc, take top N
		uniq.sort(
			key=lambda r: -float((r.get("quality") or {}).get("min") or 0),
		)
		for rec in uniq[: max(1, int(top_per_class))]:
			curated.append({
				"signature": rec.get("signature", {}),
				"recipe": rec.get("recipe", []),
				"quality": rec.get("quality", {}),
				"source": rec.get("source", "?"),
				"curated_at": _t.time(),
			})

	# Decide write target
	if out_module_path:
		out_path = Path(out_module_path)
	else:
		out_path = Path(__file__).with_name("curated_recipes_bundle.py")

	body = (
		'"""Lab-curated mesh recipes, distilled from accumulated\n'
		'learned_recipes.jsonl by `cubit_curate_learned_recipes`.\n\n'
		f'Generated: {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}\n'
		f'Source records (post-quality-filter): {len(records)}\n'
		f'Signature classes: {len(classes)}\n'
		f'Curated entries: {len(curated)}\n'
		'License: BSD-3-Clause (matching radia-mcp). The recipes are\n'
		'aggregated from multiple lab races; individual contributors\n'
		'consented via using radia-mcp under its license.\n'
		'"""\n\n'
		f'CURATED = {json.dumps(curated, ensure_ascii=False, indent=2)}\n'
	)
	try:
		out_path.parent.mkdir(parents=True, exist_ok=True)
		out_path.write_text(body, encoding="utf-8")
	except OSError as e:
		return json.dumps({"status": "error",
		                   "stage": "write_module",
		                   "error": str(e)})

	return json.dumps({
		"status": "ok",
		"pool_records": len(records),
		"signature_classes": len(classes),
		"curated_entries": len(curated),
		"output_module": str(out_path),
		"top_per_class": top_per_class,
		"min_quality_jacobian": min_quality_jacobian,
		"next_step": ("commit the new curated_recipes_bundle.py + bump "
		              "radia-mcp version + `python -m build` to ship it "
		              "to PyPI. Stage 3: write up methodology in "
		              "packages/radia-mcp/docs/design/lab_curated_recipes.md."),
	}, indent=2, ensure_ascii=False)


# --- Cubit journal helpers (capture human commands during race) ---

# Commands the AI / our daemon is known to issue around / during a race.
# These get filtered out of the captured journal so what remains is
# pure human input.
_RACE_INTERNAL_PREFIXES = (
	"save as ", "record journal", "record stop",
	"probe ", "list ",
)


def _start_journal_capture(sess, journal_path: Path) -> bool:
	"""Try to start journal recording on the live session. Returns
	True on success. Cubit syntax: `record journal "<path>" overwrite`.
	"""
	fwd = str(journal_path).replace("\\", "/")
	try:
		r = sess.call("cmd", [f'record journal "{fwd}" overwrite'],
		              timeout_s=15)
	except _cs.CubitSessionError:
		return False
	per = r.get("result", [])
	return bool(per and per[0].get("ok"))


def _stop_journal_capture(sess) -> bool:
	for cmd in ("record stop", "set journal off"):
		try:
			r = sess.call("cmd", [cmd], timeout_s=15)
		except _cs.CubitSessionError:
			continue
		per = r.get("result", [])
		if per and per[0].get("ok"):
			return True
	return False


def _parse_human_commands(journal_path: Path,
                            ai_known_cmds: set[str]) -> list[str]:
	"""Read the journal and return commands that are plausibly human-
	authored. Filters: skip blank lines, comment lines, lines matching
	internal/probe/save prefixes, and exact matches to anything the
	AI/daemon issued during the race.
	"""
	if not journal_path.exists():
		return []
	try:
		text = journal_path.read_text(encoding="utf-8", errors="replace")
	except OSError:
		return []
	out: list[str] = []
	for line in text.splitlines():
		s = line.strip()
		if not s or s.startswith("#"):
			continue
		low = s.lower()
		if any(low.startswith(p) for p in _RACE_INTERNAL_PREFIXES):
			continue
		if s in ai_known_cmds:
			continue
		out.append(s)
	return out


def _run_batch_with_quality(step_path: str | None, commands: list[str],
                              timeout_s: int,
                              cub5_path: str | None = None) -> dict:
	"""Like `_run_batch` but also probes `quality_summary` after meshing.

	Returns the same dict as `_run_batch` plus a `quality` key:
	{hex_count, tet_count, min, max, mean, below_0.2, below_0.2_pct}.
	"""
	# Build a session, run preamble + commands + quality probe in one shot
	sess = _cs.CubitSession(mode="batch")
	try:
		try:
			sess.ensure_started()
		except _cs.CubitSessionError as e:
			return {"status": "error", "stage": "start", "kind": "environment", "error": str(e)}
		per_line: list[dict] = []
		if cub5_path:
			fwd = str(cub5_path).replace("\\", "/")
			r = sess.call("cmd", [f'open "{fwd}"'], timeout_s=timeout_s)
			per_line.extend(r.get("result", []))
			if not per_line or not per_line[0].get("ok"):
				return {"status": "error", "stage": "open_cub5", "kind": "input",
				        "per_line": per_line}
		elif step_path:
			fwd = str(step_path).replace("\\", "/")
			r = sess.call("cmd", [f'import step "{fwd}"'], timeout_s=timeout_s)
			per_line.extend(r.get("result", []))
			if not per_line or not per_line[0].get("ok"):
				return {"status": "error", "stage": "import",
				        "per_line": per_line}
		try:
			r = sess.call("cmd", commands, timeout_s=timeout_s)
		except _cs.CubitSessionError as e:
			return _error_payload("rpc", str(e))
		per_line.extend(r.get("result", []))
		try:
			smy = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
		except _cs.CubitSessionError:
			smy = {}
		try:
			qual = sess.call("probe", ["quality_summary"], timeout_s=60).get("result", {})
		except _cs.CubitSessionError:
			qual = {}
		all_ok = all(s.get("ok") for s in per_line) if per_line else False
		failed_at = next((i + 1 for i, s in enumerate(per_line)
		                   if not s.get("ok")), None)
		first_err = None
		if failed_at:
			step = per_line[failed_at - 1]
			first_err = step.get("error") or f"rc={step.get('rc')}"
		return {
			"status": "ok" if all_ok else "error",
			"all_ok": all_ok,
			"per_line": per_line,
			"summary": smy,
			"quality": qual,
			"failed_at": failed_at,
			"error": (first_err or "").strip()[:400] if first_err else None,
		}
	finally:
		try:
			sess.shutdown()
		except Exception:
			pass


def _race_review_core(recipes: list,
                        max_wait_s: int,
                        max_concurrent: int,
                        include_human: bool,
                        poll_interval_s: float,
                        race_id: str | None = None) -> dict:
	"""Shared implementation for sync + async review races.

	Returns the same payload shape both modes serialize. When called
	from the async path, the caller pre-allocates `race_id` and
	persists a `running` placeholder so status checks see the in-
	progress state.
	"""
	# (Body extracted from cubit_mesh_race_review below — see that for docs.)
	import concurrent.futures as _cf
	import threading
	import time as _t

	if not recipes:
		return {"status": "error",
		        "error": "recipes list cannot be empty"}

	sess, err = _cubit_session_or_error()
	if err is not None:
		return {"status": "error", "stage": "session",
		        "error": json.loads(err)}

	label, path_or_err = _auto_checkpoint_live(sess)
	if label is None:
		return {"status": "error", "stage": "checkpoint",
		        "error": path_or_err}
	cub5_path = path_or_err

	pre = _probe_summary_safe(sess)
	pre_h = int(pre.get("hexes", 0) or 0)
	pre_t = int(pre.get("tets", 0) or 0)
	t0 = _t.time()
	if race_id is None:
		race_id = f"race_{int(t0)}"

	# Capture human commands during the race so we can learn what
	# they tried (especially when they win). Journaling is started on
	# the LIVE session only — batch instances are isolated.
	journal_path = _race_history_dir() / f"{race_id}.jou"
	journal_started = _start_journal_capture(sess, journal_path)

	human_result: dict | None = None
	human_lock = threading.Lock()
	stop_event = threading.Event()

	def _poll_human():
		nonlocal human_result
		while not stop_event.is_set():
			if _t.time() - t0 > max_wait_s:
				return
			smy = _probe_summary_safe(sess)
			h = int(smy.get("hexes", 0) or 0)
			t = int(smy.get("tets", 0) or 0)
			if (h > pre_h) or (t > pre_t):
				try:
					qual = sess.call("probe", ["quality_summary"], timeout_s=15).get("result", {})
				except _cs.CubitSessionError:
					qual = {}
				with human_lock:
					if human_result is None:
						human_result = {
							"name": "human",
							"is_ai": False,
							"hexes": h - pre_h,
							"tets": t - pre_t,
							"nodes_total": int(smy.get("nodes", 0) or 0),
							"summary": smy,
							"quality": qual,
							"elapsed_s": round(_t.time() - t0, 2),
							"recipe": ["(human-authored, not captured)"],
						}
				return
			stop_event.wait(timeout=poll_interval_s)

	def _ai_runner(rec: dict) -> dict:
		t1 = _t.time()
		r = _run_batch_with_quality(
			step_path=None, commands=list(rec.get("cmds") or []),
			timeout_s=max_wait_s, cub5_path=cub5_path,
		)
		smy = r.get("summary", {}) or {}
		qual = r.get("quality", {}) or {}
		return {
			"name": rec.get("name", "?"),
			"is_ai": True,
			"hexes": int(smy.get("hexes", 0) or 0),
			"tets": int(smy.get("tets", 0) or 0),
			"nodes_total": int(smy.get("nodes", 0) or 0),
			"summary": smy,
			"quality": qual,
			"elapsed_s": round(_t.time() - t1, 2),
			"status": r.get("status"),
			"failed_at": r.get("failed_at"),
			"error": r.get("error"),
			"rationale": rec.get("rationale", ""),
			"recipe": list(rec.get("cmds") or []),
		}

	ai_results: list[dict] = []
	with _cf.ThreadPoolExecutor(
		max_workers=max(2, int(max_concurrent) + 1),
	) as pool:
		futs = [pool.submit(_ai_runner, r) for r in recipes]
		if include_human:
			pool.submit(_poll_human)
		try:
			_cf.wait(futs, timeout=max_wait_s,
			          return_when=_cf.ALL_COMPLETED)
		finally:
			stop_event.set()
		for f in futs:
			if f.done():
				ai_results.append(f.result())

	# Stop journal capture and harvest the human's commands (if any)
	if journal_started:
		_stop_journal_capture(sess)
		ai_known = set()
		for r in recipes:
			for c in (r.get("cmds") or []):
				ai_known.add(str(c).strip())
		human_cmds = _parse_human_commands(journal_path, ai_known)
		if human_result is not None and human_cmds:
			human_result["recipe"] = human_cmds
		# If human won, persist their winning recipe for future races
		if human_result is not None and human_cmds:
			sig = _geometry_signature(pre)
			_record_learned_recipe(
				race_id=race_id, sig=sig,
				recipe=human_cmds,
				quality=human_result.get("quality") or {},
				source="human",
			)

	shortlist: list[dict] = list(ai_results)
	if include_human and human_result is not None:
		shortlist.append(human_result)

	valid = [s for s in shortlist if (s.get("hexes", 0) + s.get("tets", 0)) > 0]
	def _rank_key(s: dict):
		q = (s.get("quality") or {}).get("min")
		return (-1 if q is None else float(q),
		        s.get("hexes", 0) + s.get("tets", 0))
	valid.sort(key=_rank_key, reverse=True)

	# Also persist the *AI* winner if it was the top by quality and
	# beat any past learned recipe — accumulates the best across runs
	if valid and valid[0].get("is_ai") and (valid[0].get("quality") or {}).get("min") is not None:
		_record_learned_recipe(
			race_id=race_id,
			sig=_geometry_signature(pre),
			recipe=valid[0].get("recipe") or [],
			quality=valid[0].get("quality") or {},
			source=f"ai_{valid[0].get('name','')}",
		)

	history_path = _race_history_dir() / f"{race_id}.json"
	history = {
		"race_id": race_id,
		"state": "done",
		"started_at": t0,
		"finished_at": _t.time(),
		"checkpoint_label": label,
		"checkpoint_path": cub5_path,
		"pre_race_summary": pre,
		"journal_path": str(journal_path) if journal_started else None,
		"shortlist": valid,
		"all_results": shortlist,
		"recipes_offered": [
			{"name": r.get("name"), "cmds": r.get("cmds"),
			 "rationale": r.get("rationale", "")}
			for r in recipes
		],
	}
	try:
		history_path.write_text(json.dumps(history, ensure_ascii=False, indent=2),
		                         encoding="utf-8")
	except OSError:
		pass

	if not valid:
		recommendation = ("no variant produced any element — review "
		                  "errors in all_results / check geometry")
	else:
		top = valid[0]
		why = []
		q = (top.get("quality") or {}).get("min")
		if q is not None:
			why.append(f"highest min Jacobian {q}")
		why.append(f"{top.get('hexes', 0)} hex / {top.get('tets', 0)} tet")
		recommendation = (f"top by quality: '{top['name']}' "
		                  f"({', '.join(why)}). Apply with "
		                  f"`cubit_mesh_apply_choice(race_id='{race_id}', "
		                  f"variant_name='{top['name']}')`.")

	return {
		"status": "ok",
		"race_id": race_id,
		"history_path": str(history_path),
		"checkpoint_label": label,
		"recommendation": recommendation,
		"shortlist": [
			{
				"name": s["name"],
				"is_ai": s.get("is_ai"),
				"hexes": s.get("hexes"),
				"tets": s.get("tets"),
				"min_jacobian": (s.get("quality") or {}).get("min"),
				"mean_jacobian": (s.get("quality") or {}).get("mean"),
				"below_0.2_pct": (s.get("quality") or {}).get("below_0.2_pct"),
				"elapsed_s": s.get("elapsed_s"),
				"rationale": s.get("rationale", ""),
			}
			for s in valid
		],
		"failed_or_empty": [
			{"name": s["name"], "error": s.get("error"),
			 "status": s.get("status")}
			for s in shortlist if s not in valid
		],
	}


@mcp.tool()
def cubit_mesh_race_review_async(recipes: list,
                                   max_wait_s: int = 600,
                                   max_concurrent: int = 4,
                                   include_human: bool = True,
                                   poll_interval_s: float = 5.0) -> str:
	"""**Background launch** of a race review — returns immediately
	with a `race_id`, then runs in a daemon thread.

	**This is the recommended workflow when the AI is doing trial-
	and-error.** The user (and the AI) should not block waiting for
	the race; instead poll progress with `cubit_mesh_race_status` and
	apply the chosen winner with `cubit_mesh_apply_choice` whenever
	convenient.

	The race result is persisted to `<state_dir>/race_history/<race_id>.json`
	with `state` = `"running"` initially, then `"done"` when complete.
	`cubit_mesh_apply_choice` works as soon as `state == "done"`.

	Same args as `cubit_mesh_race_review`; same shortlist semantics
	on completion (sorted by min-Jacobian, human work never
	overwritten).
	"""
	import threading
	import time as _t

	if not recipes:
		return json.dumps({"status": "error",
		                   "error": "recipes list cannot be empty"})

	t0 = _t.time()
	race_id = f"race_{int(t0)}"
	# Write a 'running' placeholder so status polls see it
	history_path = _race_history_dir() / f"{race_id}.json"
	placeholder = {
		"race_id": race_id,
		"state": "running",
		"started_at": t0,
		"max_wait_s": max_wait_s,
		"recipes_offered": [
			{"name": r.get("name"), "rationale": r.get("rationale", "")}
			for r in recipes
		],
	}
	try:
		history_path.write_text(
			json.dumps(placeholder, ensure_ascii=False, indent=2),
			encoding="utf-8",
		)
	except OSError:
		pass

	def _runner():
		try:
			_race_review_core(
				recipes=recipes,
				max_wait_s=max_wait_s,
				max_concurrent=max_concurrent,
				include_human=include_human,
				poll_interval_s=poll_interval_s,
				race_id=race_id,
			)
		except Exception as e:  # noqa: BLE001
			err_history = dict(placeholder)
			err_history["state"] = "error"
			err_history["error"] = f"{type(e).__name__}: {e}"
			err_history["finished_at"] = _t.time()
			try:
				history_path.write_text(
					json.dumps(err_history, ensure_ascii=False, indent=2),
					encoding="utf-8",
				)
			except OSError:
				pass

	thread = threading.Thread(target=_runner, daemon=True)
	thread.start()

	return json.dumps({
		"status": "running",
		"race_id": race_id,
		"history_path": str(history_path),
		"started_at": t0,
		"max_wait_s": max_wait_s,
		"variant_count": len(recipes),
		"note": ("Race is running in the background. Poll with "
		         f"`cubit_mesh_race_status(race_id='{race_id}')` and "
		         f"apply the winner with `cubit_mesh_apply_choice("
		         f"race_id='{race_id}', variant_name=...)`. The user "
		         f"can keep working in the live Cubit GUI in the "
		         f"meantime."),
	}, indent=2, ensure_ascii=False)


@mcp.tool()
def cubit_mesh_race_status(race_id: str) -> str:
	"""Check the status of a background race launched by
	`cubit_mesh_race_review_async` (or `_smart_async`).

	Returns the persisted history JSON. Key field: `state` —
	`"running"` while in progress, `"done"` when complete (then
	`shortlist` is populated), `"error"` if the daemon thread
	crashed.
	"""
	history_path = _race_history_dir() / f"{race_id}.json"
	if not history_path.exists():
		return json.dumps({"status": "error",
		                   "error": f"no race history for {race_id!r}",
		                   "expected": str(history_path)})
	try:
		body = history_path.read_text(encoding="utf-8")
	except OSError as e:
		return json.dumps({"status": "error", "error": str(e)})
	return body  # already JSON


@mcp.tool()
def cubit_mesh_race_review(recipes: list,
                             max_wait_s: int = 600,
                             max_concurrent: int = 4,
                             include_human: bool = True,
                             poll_interval_s: float = 5.0) -> str:
	"""Race N AI variants + observe live human, **wait for all** to
	finish (or `max_wait_s`), and return a shortlist sorted by quality.

	**Does NOT auto-commit anything.** The caller (typically the human
	via the AI assistant) inspects the shortlist and picks the variant
	to apply via `cubit_mesh_apply_choice(race_id, variant_name)`.
	The human's own work is one of the candidates if they meshed
	during the window — never silently overwritten.

	Each variant report carries:
	  * hex / tet / nodes counts
	  * min / max / mean scaled-Jacobian (mesh quality)
	  * below_0.2 element count + percentage
	  * elapsed seconds
	  * (for AI variants) the originating recipe so apply_choice can
	    replay it

	Args:
	    recipes: list of `{"name": str, "cmds": [str], "rationale": str?}`.
	    max_wait_s: overall timeout (race ends and shortlist is built
	        whenever this elapses, even if some variants still running).
	    max_concurrent: subprocess pool size.
	    include_human: also poll the live GUI and add 'human' to the
	        shortlist if a mesh appears.
	    poll_interval_s: live-GUI poll interval (only used if
	        include_human=True).

	Returns JSON: race_id (for apply_choice), shortlist sorted by
	min-Jacobian descending, persisted history path.

	**Note**: this is the **synchronous** form — it blocks until the
	race finishes. For background AI iteration prefer
	`cubit_mesh_race_review_async` (signature workflow per Sugawara
	Lab: AI thinking should not block the user).
	"""
	r = _race_review_core(
		recipes=recipes,
		max_wait_s=max_wait_s,
		max_concurrent=max_concurrent,
		include_human=include_human,
		poll_interval_s=poll_interval_s,
	)
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def cubit_mesh_apply_choice(race_id: str, variant_name: str,
                              timeout_s: int = 600) -> str:
	"""Apply the human's chosen variant from a previous
	`cubit_mesh_race_review`.

	Looks up the race history by `race_id`, finds the variant by
	`variant_name`, and replays its recipe in the live Cubit GUI
	(after `delete mesh` to clear any prior attempt). If the chosen
	variant is "human", this is a no-op — the human's mesh is
	already live.

	Args:
	    race_id: identifier returned from `cubit_mesh_race_review`.
	    variant_name: which entry in the shortlist to apply.

	Returns JSON with the live state after apply.
	"""
	history_path = _race_history_dir() / f"{race_id}.json"
	if not history_path.exists():
		return json.dumps({"status": "error",
		                   "error": f"no race history for {race_id!r}",
		                   "expected": str(history_path)})
	try:
		history = json.loads(history_path.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError) as e:
		return json.dumps({"status": "error",
		                   "stage": "history_read",
		                   "error": str(e)})

	all_results = history.get("all_results", [])
	choice = next((s for s in all_results if s.get("name") == variant_name), None)
	if choice is None:
		known = [s.get("name") for s in all_results]
		return json.dumps({"status": "error",
		                   "error": f"variant {variant_name!r} not in race {race_id}",
		                   "known": known})

	if not choice.get("is_ai"):
		return json.dumps({
			"status": "ok",
			"action": "no-op",
			"reason": "human variant chosen — mesh is already live",
			"choice": choice,
		}, indent=2)

	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	# Replay: delete current mesh + run the recipe
	cmds = ["delete mesh"] + list(choice.get("recipe") or [])
	try:
		r = sess.call("cmd", cmds, timeout_s=timeout_s)
	except _cs.CubitSessionError as e:
		return json.dumps({"status": "error", "stage": "live_rpc",
		                   "error": str(e),
		                   "checkpoint_label": history.get("checkpoint_label"),
		                   "note": "rollback via cubit_restore(label=...)"})
	per_line = r.get("result", [])
	smy = _probe_summary_safe(sess)
	all_ok = all(s.get("ok") for s in per_line) if per_line else False
	return json.dumps({
		"status": "ok" if all_ok else "error",
		"variant_name": variant_name,
		"recipe_replayed": cmds,
		"live_summary": smy,
		"checkpoint_label": history.get("checkpoint_label"),
		"note": ("Rollback via cubit_restore(label=...) if you change "
		         "your mind."),
	}, indent=2)


# ============================================================
# Smart race: AI introspects geometry + generates candidates, then races
# ============================================================

def _generate_smart_recipes(sess, target_size: float,
                              n_variants: int) -> tuple[list[dict], list[str]]:
	"""Inspect the live Cubit state (volumes, surfaces, schemes,
	mesh state) and emit up to `n_variants` candidate recipes
	with rationale. Returns (recipes, rationale_lines).

	Heuristics, in priority order:
	  - if compound topology suspected (surf/vol >= 7), prefer
	    webcut + scheme auto FIRST.
	  - simple prismatic: scheme auto, scheme sweep, then refined size.
	  - large vol count: split-and-conquer with smaller size.
	  - always include scheme tetmesh as a guaranteed-element fallback.
	  - vary `target_size` (half / base / double) to give the race
	    different element counts to compare.
	"""
	rationale: list[str] = []

	# Probe state
	try:
		smy = sess.call("probe", ["summary"], timeout_s=15).get("result", {})
	except _cs.CubitSessionError:
		smy = {}
	try:
		per_vol = sess.call("probe", ["per_volume"], timeout_s=15).get("result", []) or []
	except _cs.CubitSessionError:
		per_vol = []

	vol = int(smy.get("volumes", 0) or 0)
	surf = int(smy.get("surfaces", 0) or 0)
	hex_pre = int(smy.get("hexes", 0) or 0)
	tet_pre = int(smy.get("tets", 0) or 0)
	rationale.append(f"state: {vol} volumes / {surf} surfaces / "
	                  f"{hex_pre} hex / {tet_pre} tet pre-race")

	if vol == 0:
		rationale.append("no volumes — nothing to mesh; aborting")
		return [], rationale

	compound_suspected = vol <= 3 and (surf / max(1, vol)) >= 7
	unmeshed = [int(r["id"]) for r in per_vol if not r.get("meshed")]
	if unmeshed and len(unmeshed) < vol:
		rationale.append(f"partial mesh detected; unmeshed vols: "
		                  f"{unmeshed[:8]}{'...' if len(unmeshed)>8 else ''}")
	if compound_suspected:
		rationale.append(f"compound topology suspected (surf/vol = "
		                  f"{surf/vol:.1f} >= 7)")

	scope = "all" if not unmeshed else " ".join(str(v) for v in unmeshed)

	candidates: list[dict] = []

	def add(name: str, cmds: list[str], why: str) -> None:
		if any(c["name"] == name for c in candidates):
			return
		candidates.append({"name": name, "cmds": cmds, "rationale": why})
		rationale.append(f"+ {name}: {why}")

	# Learned-recipe seed: pull past winners on similar geometries
	# (human or AI) to head the candidate list. This is how the AI
	# "learns from the human" — every prior race's human winner is
	# replayed as a candidate next time.
	sig = _geometry_signature(smy)
	rationale.append(f"signature: {sig}")
	learned = _load_learned_recipes_for(sig, max_n=2)
	if learned:
		rationale.append(f"loaded {len(learned)} learned recipe(s) "
		                  f"for similar geometry")
	for i, lr in enumerate(learned):
		src = lr.get("source", "?")
		q = (lr.get("quality") or {}).get("min")
		add(
			f"learned_{i}_{src}",
			list(lr.get("recipe") or []),
			f"past winner ({src}) on similar geometry "
			f"(min Jacobian {q}); replaying to see if it still works.",
		)

	# Always start with the cheapest most-likely-to-win
	add(
		f"auto_size_{target_size}",
		[f"volume {scope} size {target_size}",
		 f"volume {scope} scheme auto",
		 f"mesh volume {scope}"],
		"scheme auto picks per-volume; usually wins on simple prisms.",
	)

	if compound_suspected:
		# Geometry-split rung first when compound suspected
		import os as _os
		webcut_r = float(_os.environ.get("RADIA_MCP_WEBCUT_RADIUS", "16.01"))
		add(
			f"webcut_cyl_r{webcut_r}_then_auto",
			[f"webcut volume {scope} with cylinder radius {webcut_r} axis z",
			 "imprint all", "merge all",
			 f"volume all size {target_size}",
			 "volume all scheme auto",
			 "mesh volume all"],
			f"compound body → cylinder webcut at z-axis r={webcut_r} "
			f"(env: RADIA_MCP_WEBCUT_RADIUS) then re-attempt scheme auto.",
		)
		add(
			"polyhedron",
			[f"volume {scope} size {target_size}",
			 f"volume {scope} scheme polyhedron",
			 f"mesh volume {scope}"],
			"hex-dominant, robust on non-prismatic topology.",
		)

	add(
		f"sweep_size_{target_size}",
		[f"volume {scope} size {target_size}",
		 f"volume {scope} scheme sweep",
		 f"mesh volume {scope}"],
		"force sweep — best hex yield when volumes are extruded prisms.",
	)

	# Vary the size axis to give the race element-count variety
	if target_size > 0:
		half = round(target_size / 2.0, 4)
		add(
			f"auto_size_{half}_finer",
			[f"volume {scope} size {half}",
			 f"volume {scope} scheme auto",
			 f"mesh volume {scope}"],
			f"finer size ({half}) — captures small features but slower.",
		)
		double = round(target_size * 2.0, 4)
		if double != target_size:
			add(
				f"auto_size_{double}_coarser",
				[f"volume {scope} size {double}",
				 f"volume {scope} scheme auto",
				 f"mesh volume {scope}"],
				f"coarser size ({double}) — fastest, useful for sanity check.",
			)

	# Tet fallback always last (lowest priority — gives tet, not hex)
	add(
		f"tetmesh_size_{target_size}",
		[f"volume {scope} size {target_size}",
		 f"volume {scope} scheme tetmesh",
		 f"mesh volume {scope}"],
		"guaranteed-element fallback (tet only, not hex).",
	)

	# Truncate to requested N
	candidates = candidates[: max(1, int(n_variants))]
	rationale.insert(1, f"selected {len(candidates)} of "
	                     f"{len(candidates)} candidates "
	                     f"(n_variants={n_variants})")
	return candidates, rationale


@mcp.tool()
def cubit_mesh_race_smart_async(target_size: float = 1.0,
                                  n_variants: int = 4,
                                  max_wait_s: int = 600,
                                  max_concurrent: int = 4,
                                  poll_interval_s: float = 5.0,
                                  include_human: bool = True) -> str:
	"""**Background variant of `cubit_mesh_race_smart`** — AI inspects
	the live geometry, generates *N* candidate recipes, and races them
	in a daemon thread; returns immediately with `race_id`.

	Recommended workflow when the AI is doing trial-and-error: don't
	make the user wait. Poll with `cubit_mesh_race_status(race_id)`,
	apply with `cubit_mesh_apply_choice(race_id, variant_name)`.

	The smart-recipe rationale + race shortlist are both persisted to
	`<state_dir>/race_history/<race_id>.json`.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err
	recipes, rationale = _generate_smart_recipes(sess, target_size, n_variants)
	if not recipes:
		return json.dumps({"status": "error",
		                   "reason": "no candidates generated",
		                   "rationale": rationale})
	# Hand off to async review with the generated recipes
	r = json.loads(cubit_mesh_race_review_async(
		recipes=recipes,
		max_wait_s=max_wait_s,
		max_concurrent=max_concurrent,
		include_human=include_human,
		poll_interval_s=poll_interval_s,
	))
	r["rationale"] = rationale
	r["candidates_attempted"] = [
		{"name": rec["name"], "rationale": rec.get("rationale", "")}
		for rec in recipes
	]
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def cubit_mesh_race_smart(target_size: float = 1.0,
                            n_variants: int = 4,
                            prefer: str = "hex",
                            commit_winner: bool = True,
                            poll_interval_s: float = 5.0,
                            max_wait_s: int = 600,
                            max_concurrent: int = 4) -> str:
	"""**The flagship workflow**: AI inspects the live Cubit state,
	auto-generates *N* candidate mesh recipes (with rationale), and
	races them against the human via `cubit_mesh_race_with_human`.

	The user just says "AI, mesh this for me" — no need to hand-write
	recipes. The AI:

	1. Probes the live geometry (`probe summary` + `probe per_volume`).
	2. Detects compound topology (surf/vol ratio), partial-mesh state
	   (unmeshed volume list), and total complexity.
	3. Generates *N* recipes via heuristics:
	     - scheme auto (always first — cheapest)
	     - webcut + scheme auto (when compound suspected)
	     - scheme polyhedron (compound robust hex)
	     - scheme sweep (prismatic best-hex)
	     - finer / coarser size variants (element-count variety)
	     - scheme tetmesh (guaranteed fallback)
	4. Races them in parallel, **while the human keeps working**.
	5. Returns winner + rationale (so the AI can explain WHY this
	   recipe was chosen).

	Args:
	    target_size: base mesh size (`volume all size X`); finer/coarser
	        variants automatically explore X/2 and 2X.
	    n_variants: cap on candidates (default 4 — safe for 16 GB RAM).
	    prefer: "hex" (require hex>0) or "any".
	    commit_winner: replay winning recipe in live GUI on success.
	    poll_interval_s / max_wait_s / max_concurrent: race tunables.

	Returns JSON: rationale lines (the AI's reasoning) +
	winner + ai_results + checkpoint_label + gui report.
	"""
	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	recipes, rationale = _generate_smart_recipes(sess, target_size, n_variants)
	if not recipes:
		return json.dumps({
			"status": "error",
			"reason": "no candidates generated",
			"rationale": rationale,
		}, indent=2)

	# Race them
	race_raw = cubit_mesh_race_with_human(
		recipes=recipes,
		poll_interval_s=poll_interval_s,
		max_wait_s=max_wait_s,
		max_concurrent=max_concurrent,
		commit_winner=commit_winner,
		prefer=prefer,
	)
	race = json.loads(race_raw)
	# Augment with the candidate-generation rationale so the AI can
	# explain its choices to the user
	race["rationale"] = rationale
	race["candidates_attempted"] = [
		{"name": r["name"], "rationale": r["rationale"]}
		for r in recipes
	]
	return json.dumps(race, indent=2)


# ============================================================
# Human-vs-AI race: live human + N AI batches, first-to-mesh wins
# ============================================================

@mcp.tool()
def cubit_mesh_race_with_human(recipes: list,
                                 poll_interval_s: float = 5.0,
                                 max_wait_s: int = 600,
                                 max_concurrent: int = 4,
                                 commit_winner: bool = False,
                                 prefer: str = "hex") -> str:
	"""**The radia-mcp signature workflow.**

	User: "AI, this geometry won't mesh — try something while I keep
	working on it."

	Setup:
	  1. Snapshot the live Cubit GUI state to `.cub5` (reuses
	     `_auto_checkpoint_live`).
	  2. Spawn `len(recipes)` batch Cubits (capped to `max_concurrent`),
	     each opening the snapshot and running its candidate recipe.
	  3. **Meanwhile poll the live GUI session** every `poll_interval_s`
	     for hex/tet count. The user is expected to be actively
	     meshing in their window.
	  4. **First to produce mesh** (`hex > 0` for `prefer='hex'`,
	     `(hex+tet) > 0` for `prefer='any'`) wins. Could be the human,
	     could be one of the AI variants.
	  5. Remaining batches are killed; live GUI is **not** modified
	     unless `commit_winner=True` AND the human hadn't already
	     produced mesh (so we never overwrite the user's work).

	Why this exists: waiting for AI to finish is wasteful when you can
	just keep working. This races human and AI fairly and respects
	whoever finishes first — but defaults to "report only" so the
	user has final say on commit.

	Args:
	    recipes: list of `{"name": str, "cmds": [str]}`. Each runs from
	        the live snapshot.
	    poll_interval_s: how often to ping the live session.
	    max_wait_s: overall timeout.
	    max_concurrent: subprocess pool size.
	    commit_winner: if True AND winner is an AI variant AND the
	        human hadn't meshed yet, replay the winning recipe in
	        the live GUI. Default False (just report).
	    prefer: "hex" or "any".

	Returns JSON: winner identity (human / AI variant name) +
	per-variant outcome + checkpoint label for revert.
	"""
	import concurrent.futures as _cf
	import threading
	import time as _t

	if not recipes:
		return json.dumps({"status": "error",
		                   "error": "recipes list cannot be empty"})

	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	# 1. Snapshot live state — used as the seed for every AI batch
	label, path_or_err = _auto_checkpoint_live(sess)
	if label is None:
		return json.dumps({"status": "error",
		                   "stage": "checkpoint",
		                   "error": path_or_err})
	cub5_path = path_or_err

	# Snapshot pre-race element counts so "human won" means new mesh
	# created AFTER we started, not pre-existing mesh from before.
	pre = _probe_summary_safe(sess)
	pre_h = int(pre.get("hexes", 0) or 0)
	pre_t = int(pre.get("tets", 0) or 0)

	winner_event = threading.Event()
	winner_lock = threading.Lock()
	winner: dict | None = None  # mutated under lock

	def _claim(payload: dict) -> bool:
		"""Atomically claim winner. Returns True if this caller wins."""
		nonlocal winner
		with winner_lock:
			if winner is not None:
				return False
			winner = payload
			winner_event.set()
			return True

	def _human_poller():
		"""Poll live GUI; if mesh appears, claim 'human'."""
		deadline = _t.time() + max_wait_s
		while not winner_event.is_set() and _t.time() < deadline:
			smy = _probe_summary_safe(sess)
			h = int(smy.get("hexes", 0) or 0)
			t = int(smy.get("tets", 0) or 0)
			ok = (h > pre_h) if prefer == "hex" \
				else ((h + t) > (pre_h + pre_t))
			if ok:
				_claim({
					"winner": "human",
					"summary": smy,
					"new_hex": h - pre_h,
					"new_tet": t - pre_t,
				})
				return
			# wait, but wake on event for fast cancel
			winner_event.wait(timeout=poll_interval_s)

	def _ai_runner(rec: dict) -> dict:
		"""Run one recipe in batch, claim winner if it produces mesh."""
		name = rec.get("name", "?")
		cmds = list(rec.get("cmds") or [])
		if not cmds:
			return {"name": name, "status": "skip", "error": "empty cmds"}
		# Bail out early if already won
		if winner_event.is_set():
			return {"name": name, "status": "cancelled"}
		r = _run_batch(step_path=None, commands=cmds,
		                timeout_s=max_wait_s, cub5_path=cub5_path)
		smy = r.get("summary", {}) or {}
		h = int(smy.get("hexes", 0) or 0)
		t = int(smy.get("tets", 0) or 0)
		ok = r.get("status") == "ok" and \
			((h > 0) if prefer == "hex" else ((h + t) > 0))
		out = {
			"name": name,
			"status": r.get("status"),
			"hexes": h, "tets": t,
			"nodes": int(smy.get("nodes", 0) or 0),
			"recipe": cmds,
		}
		if ok and not winner_event.is_set():
			_claim({
				"winner": "ai",
				"name": name,
				"recipe": cmds,
				"hexes": h, "tets": t,
				"nodes": int(smy.get("nodes", 0) or 0),
			})
		return out

	# Run all in parallel: poller as 1 thread + AI as N threads.
	with _cf.ThreadPoolExecutor(
		max_workers=max(2, int(max_concurrent) + 1),
	) as pool:
		poller_fut = pool.submit(_human_poller)
		ai_futs = [pool.submit(_ai_runner, rec) for rec in recipes]
		# Block until either winner_event is set or all AI futs done
		_cf.wait(ai_futs + [poller_fut],
		         timeout=max_wait_s,
		         return_when=_cf.ALL_COMPLETED)
		# Make sure poller knows we're done so it stops polling
		winner_event.set()
		ai_results = [f.result() for f in ai_futs if f.done()]

	if winner is None:
		return json.dumps({
			"status": "timeout",
			"reason": "neither human nor any AI variant produced mesh",
			"checkpoint_label": label,
			"checkpoint_path": cub5_path,
			"results": ai_results,
		}, indent=2)

	# Optional commit: only if winner is AI and human hadn't already meshed
	gui_report: dict = {"skipped": True}
	if commit_winner and winner.get("winner") == "ai":
		smy_now = _probe_summary_safe(sess)
		h_now = int(smy_now.get("hexes", 0) or 0)
		t_now = int(smy_now.get("tets", 0) or 0)
		human_did_mesh = (h_now > pre_h) or (t_now > pre_t)
		if human_did_mesh:
			gui_report = {
				"status": "skipped",
				"reason": ("human already produced mesh while AI was "
				           "racing — refusing to overwrite"),
				"live_summary": smy_now,
			}
		else:
			recipe = winner["recipe"]
			try:
				rr = sess.call("cmd", recipe, timeout_s=max_wait_s)
			except _cs.CubitSessionError as e:
				gui_report = {"status": "error", "stage": "gui_rpc",
				              "error": str(e)}
			else:
				per_line = rr.get("result", [])
				smy_after = _probe_summary_safe(sess)
				gui_report = {
					"status": "ok" if all(s.get("ok") for s in per_line) else "error",
					"summary": smy_after,
					"replayed": recipe,
				}

	return json.dumps({
		"status": "ok",
		"winner": winner,
		"checkpoint_label": label,
		"checkpoint_path": cub5_path,
		"pre_race_summary": pre,
		"ai_results": ai_results,
		"gui": gui_report,
		"note": ("Snapshot saved as a .cub5 — revert any AI commit via "
		         f"`cubit_restore(label='{label}')`."),
	}, indent=2)


# ============================================================
# Safe live exec: auto-checkpoint + batch dry-run + live apply
# ============================================================

def _auto_checkpoint_live(sess) -> tuple[str, str] | tuple[None, str]:
	"""Create an `autosafe_<timestamp>` .cub5 snapshot of the live
	session. Returns (label, cub5_path) on success, or (None, error_msg).

	Uses the same directory layout as `cubit_checkpoint` so the user can
	list / restore it with the existing tool.
	"""
	import time as _t
	# No leading underscore — _sanitize_label strips them.
	label = f"autosafe_{int(_t.time())}"
	try:
		label = _sanitize_label(label)
	except ValueError as e:
		return None, str(e)
	path = _checkpoints_dir() / f"{label}.cub5"
	fwd = str(path).replace("\\", "/")
	try:
		r = sess.call("cmd", [f'save as "{fwd}" overwrite'], timeout_s=120)
	except _cs.CubitSessionError as e:
		return None, f"rpc: {e}"
	per_line = r.get("result", [])
	if not per_line or not per_line[0].get("ok"):
		return None, f"save failed: {per_line}"
	if not path.exists():
		return None, f"save reported ok but {path} not written"
	return label, str(path)


@mcp.tool()
def cubit_exec_safely(commands: list, timeout_s: int = 600) -> str:
	"""Execute commands against the live GUI **with a pre-save and a
	batch dry-run as safety gates**.

	Workflow:
	  1. **Auto-checkpoint** the current live GUI session to
	     `_autosafe_<timestamp>.cub5`.
	  2. **Batch dry-run**: spawn a fresh headless Cubit, open the
	     cub5, run the commands, check for errors.
	  3. If batch passes cleanly → **reflect to live GUI** via
	     `cubit_exec`. User watches the successful recipe execute in
	     their real window.
	  4. If batch fails → **live GUI untouched**; return the error
	     and the checkpoint label so the user stays in control.
	  5. Label is returned regardless, so the user can later
	     `cubit_restore(label)` to roll back (e.g., if they don't like
	     how the mesh looks even after a successful apply).

	Use this for any risky / destructive operation on a live model the
	user has been building (mesh scheme changes, boolean ops, webcut).
	For read-only queries (probe, list, quality) use `cubit_exec`
	directly — the safety overhead is wasteful.

	Args:
	    commands: list of Cubit command strings.
	    timeout_s: per-rung timeout (checkpoint save + batch open +
	        commands + live apply).

	Returns JSON with:
	  - `checkpoint_label` / `checkpoint_path` — always present
	  - `batch` — dry-run report (per_line, summary, ok/error)
	  - `live` — live apply report (only if batch passed)
	  - `applied_to_gui` — true iff the live session was modified
	"""
	if not commands:
		return json.dumps({"status": "error",
		                   "error": "commands list cannot be empty"})
	cmd_list = [str(c).rstrip() for c in commands]

	sess, err = _cubit_session_or_error()
	if err is not None:
		return err

	# 1. checkpoint live
	label, path_or_err = _auto_checkpoint_live(sess)
	if label is None:
		return json.dumps({"status": "error", "stage": "checkpoint",
		                   "error": path_or_err})
	cub5_path = path_or_err

	# 2. batch dry-run seeded from the checkpoint
	batch = _run_batch(step_path=None, commands=cmd_list,
	                    timeout_s=timeout_s, cub5_path=cub5_path)

	payload: dict = {
		"status": "ok",
		"checkpoint_label": label,
		"checkpoint_path": cub5_path,
		"batch": batch,
		"applied_to_gui": False,
	}

	if batch.get("status") != "ok":
		payload["status"] = "error"
		payload["stage"] = "batch"
		payload["note"] = ("Batch dry-run failed — live GUI NOT touched. "
		                   "Inspect `batch` for the failing command, fix, "
		                   "and retry. Rollback with "
		                   f"`cubit_restore(label='{label}')` if needed.")
		return json.dumps(payload, indent=2)

	# 3. batch passed → reflect to live
	try:
		r = sess.call("cmd", cmd_list, timeout_s=timeout_s)
	except _cs.CubitSessionError as e:
		payload["status"] = "error"
		payload["stage"] = "live_rpc"
		payload["live_error"] = str(e)
		payload["note"] = ("Batch passed but live RPC failed. "
		                   f"Rollback: `cubit_restore(label='{label}')`.")
		return json.dumps(payload, indent=2)

	per_line = r.get("result", [])
	all_ok = all(step.get("ok") for step in per_line) if per_line else False
	try:
		smy = sess.call("probe", ["summary"], timeout_s=30).get("result", {})
	except _cs.CubitSessionError:
		smy = {}
	payload["live"] = {
		"status": "ok" if all_ok else "error",
		"per_line": per_line,
		"summary": smy,
	}
	payload["applied_to_gui"] = True
	if not all_ok:
		payload["status"] = "error"
		payload["stage"] = "live_mid"
		payload["note"] = ("Live exec partially failed (batch had passed). "
		                   "Cubit state is now mixed; consider rollback via "
		                   f"`cubit_restore(label='{label}')`.")
	return json.dumps(payload, indent=2)


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
		"5. For curved mesh: cubit.cmd('export netgen \"mesh.vol\" order N overwrite')\n"
		"6. Verify volume after curving\n"
		"7. For Gmsh visualization: cubit.cmd('export gmsh \"mesh.msh\" order N overwrite')\n"
	)


@mcp.prompt()
def cubit_to_ngsolve(mesh_file: str) -> str:
	"""Set up Cubit mesh -> NGSolve high-order FEM workflow."""
	return (
		f"Set up a Cubit -> NGSolve high-order workflow for: {mesh_file}\n\n"
		"Workflow:\n"
		"1. Create geometry and mesh in Cubit\n"
		"2. cubit.cmd('export netgen \"mesh.vol\" order N overwrite')\n"
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
		"| NGSolve FEM (recommended) | export netgen \"f.vol\" order N | 1-5 |\n"
		"| GMSH visualization | export gmsh \"f.msh\" order N | 1-3 |\n"
		"| Nastran / JMAG | export jmag_nastran \"f.bdf\" order N | 1-2 |\n"
		"| ParaView | export vtk \"f.vtk\" order N | 1-2 |\n"
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
# Custom Toolbar / In-Cubit PySide6 (Coreform 2025 GUI webinar)
# ============================================================

@mcp.tool()
def cubit_toolbar_guide(topic: str = "overview") -> str:
	"""
	Get documentation on building custom in-Cubit PySide6 toolbars.

	Covers the Coreform 2025 custom-toolbar workflow: Python script
	conventions (shebang, namespace, import gotchas), the Claro parent
	helper, LLM prompt templates for generating dialogs, packaging as
	.tar.gz, encapsulation patterns (tire cross-section example), and
	quality convergence loops.

	This is DIFFERENT from the Radia-NGSolve analysis panel conventions
	(which are external Python 3.12 launchers). For those, use
	`cubit_docs topic=panel_conventions`.

	Args:
		topic: Topic to document. Options:
			"all"                    - All sections
			"overview"               - When to use, entry points, button types
			"python_conventions"     - Shebang, _coreform_cubit namespace,
			                           import parens gotcha, cubit.cmd idioms
			"claro"                  - find_claro() parent-window helper
			"llm_prompt"             - LLM prompt template for dialogs
			"packaging"              - .tar.gz layout + Cubit import
			"encapsulation"          - Multi-step workflow -> one button
			                           (tire cross-section case study)
			"quality_convergence"    - Mesh -> query quality -> refine loop
			"troubleshooting"        - Common failure modes

	Aliases: import_gotcha, find_claro, prompt_template, tar_gz,
	         tire_example, auto_refine, debug.
	"""
	return get_toolbar_documentation(topic)


@mcp.tool()
def cubit_cpp_sdk_guide(topic: str = "overview") -> str:
	"""
	Get documentation on building Cubit C++ SDK plugins.

	Covers the Coreform 2026 "Advanced customization" webinar workflow:
	when to reach for the C++ SDK instead of the Python toolbar path,
	CMake + Qt6 + MSVC/GCC prerequisites, the SDK ``components`` example
	(Claro main-window handle, menu_manager, toolbar_manager, panel
	factory, navigation-node tree), how to register new command-line
	verbs that integrate with Cubit's command parser, and Tools->Plugins
	loading semantics (folder-not-DLL + mandatory restart).

	This is DIFFERENT from two adjacent topics:
	  * `cubit_toolbar_guide` -- the in-process Python / PySide6 path.
	    Prefer that whenever it suffices.
	  * `cubit_docs topic=panel_conventions` -- the Radia-NGSolve
	    external Python 3.12 launcher panels (unrelated to the SDK).

	Args:
		topic: Topic to document.  Options:
			"all"               - All sections
			"overview"          - When the SDK is the right choice
			"prerequisites"     - CMake / Qt6 / MSVC / GCC versions
			"install_and_load"  - Tools->Plugins folder + restart
			"components_example"- SDK components/ anatomy, Claro, menu/
			                      toolbar/panel managers, factory pattern
			"custom_commands"   - Register new command verbs + hotkeys
			"troubleshooting"   - Debug/Release mismatch, MOC/UIC,
			                      hotkey gotchas

	Aliases: cmake, qt6, build, load, restart, claro, menu, toolbar,
	         panel, navigation_nodes, hotkey, shortcut, moc, uic,
	         crash, debug.

	Source: https://www.youtube.com/watch?v=u7GIqTKId98
	        ("Advanced customization in Coreform Cubit using the C++
	        SDK", Carl McKelvey).
	"""
	return get_cpp_sdk_documentation(topic)


@mcp.tool()
def cubit_scaffold_toolbar(toolbar_name: str, buttons_json: str,
                            out_dir: str = "") -> str:
	"""
	Generate a complete Coreform-Cubit custom-toolbar skeleton on disk.

	Writes:
		<out_dir>/<toolbar_name>/
			toolbar.xml                  - Toolbar definition (buttons, icons)
			scripts/cubit_util.py        - Shared find_claro() helper
			scripts/<button>.py or .jou  - One per button
			resources/icons/             - Place button PNGs here
			package-win.ps1              - Windows packaging -> .tgz
			package-linux.sh             - Linux packaging -> .tgz
			README.md                    - Layout + install instructions

	All Python scripts include the correct `#!python` shebang,
	`_coreform_cubit` namespace guard, `find_claro()` parent lookup, and
	backslash-continuation imports (avoiding the Cubit importer
	parenthesis bug).

	Args:
		toolbar_name: Name of the toolbar, e.g. "radia_ih_helpers". Used
		    as the directory name and toolbar namespace inside Cubit.
		buttons_json: JSON array of button specs. Each element:
		    {
		      "label": "Clean + Mesh",          # required, shown in UI
		      "kind": "python" | "cubit_script", # default "python"
		      "script_basename": "clean_and_mesh", # default "button"
		      "body": "cubit.cmd('...')"         # initial script body
		    }
		out_dir: Parent directory to write into. Defaults to the current
		    working directory. Must exist. The scaffolder creates
		    `<out_dir>/<toolbar_name>/`.

	Returns:
		Summary of files written, with absolute paths.
	"""
	import json

	try:
		buttons = json.loads(buttons_json) if buttons_json else []
	except json.JSONDecodeError as exc:
		return f"ERROR: buttons_json is not valid JSON: {exc}"
	if not isinstance(buttons, list):
		return "ERROR: buttons_json must be a JSON array"

	parent = Path(out_dir).expanduser().resolve() if out_dir else Path.cwd()
	if not parent.exists():
		return f"ERROR: out_dir does not exist: {parent}"

	target = parent / toolbar_name
	try:
		files = generate_toolbar_skeleton(toolbar_name, buttons)
	except Exception as exc:
		return f"ERROR: skeleton generation failed: {type(exc).__name__}: {exc}"

	target.mkdir(parents=True, exist_ok=True)
	(target / "scripts").mkdir(exist_ok=True)
	(target / "resources" / "icons").mkdir(parents=True, exist_ok=True)

	written = []
	for rel, content in files.items():
		full = target / rel
		full.parent.mkdir(parents=True, exist_ok=True)
		full.write_text(content, encoding='utf-8')
		written.append(str(full))

	return (
		f"Toolbar scaffold written to: {target}\n\n"
		"Files:\n"
		+ "\n".join(f"  {p}" for p in written)
		+ "\n\nNext steps:\n"
		f"  1. Add button icons under {target / 'resources' / 'icons'}\n"
		f"  2. Flesh out the script bodies under {target / 'scripts'}\n"
		f"  3. Package:  pwsh {target / 'package-win.ps1'}\n"
		f"              (or:  bash {target / 'package-linux.sh'})\n"
		f"  4. Open Cubit -> right-click toolbar panel -> Import "
		f"{toolbar_name}.tgz\n"
	)


@mcp.tool()
def cubit_generate_dialog(title: str, fields_json: str,
                           commands_json: str) -> str:
	"""
	Generate a single complete PySide6 dialog script for a Cubit toolbar
	button.

	The returned script follows every Cubit in-process Python
	convention: `#!python` shebang, no `import cubit`,
	`_coreform_cubit` namespace guard, `find_claro()` parent lookup,
	backslash-continuation imports, try/except around float conversion,
	and QMessageBox error UX. Paste the output directly into a custom
	toolbar Python-button or save it under `scripts/` in a toolbar
	scaffold.

	Args:
		title: Dialog window title, e.g. "Create Parameterized Brick".
		fields_json: JSON array of field specs. Each:
		    {
		      "name": "size",       # attribute name (-> self.size_edit)
		      "label": "Brick size", # label shown next to QLineEdit
		      "default": "1.0",     # default text in the edit
		      "kind": "float"       # "float" | "int" | "text"
		    }
		commands_json: JSON array of Cubit command strings. Each may use
		    `{<field_name>}` placeholders that are f-string substituted
		    at runtime. Example:
		      ["create brick x {size}",
		       "create cylinder radius {radius} height {size}",
		       "subtract volume 2 from volume 1"]

	Returns:
		Complete Python script as a string, ready to paste into Cubit.
	"""
	import json

	try:
		fields = json.loads(fields_json) if fields_json else []
		commands = json.loads(commands_json) if commands_json else []
	except json.JSONDecodeError as exc:
		return f"ERROR: JSON parse error: {exc}"
	if not isinstance(fields, list):
		return "ERROR: fields_json must be a JSON array"
	if not isinstance(commands, list):
		return "ERROR: commands_json must be a JSON array"

	try:
		script = generate_dialog_skeleton(title, fields, commands)
	except Exception as exc:
		return f"ERROR: dialog generation failed: {type(exc).__name__}: {exc}"
	return script


# ============================================================
# Mesh Diagnostics (Coreform 2025 "first 15 minutes" webinar)
# ============================================================

@mcp.tool()
def cubit_diagnostics_guide(topic: str = "overview") -> str:
	"""
	Get the foundational mesh-diagnostics + cleanup + quality playbook.

	This is the "why" layer below `cubit_docs scripting_*` - it explains
	the reasoning behind imprint+merge, the Exodus data model, small-
	feature cleanup, gap/overlap handling, mesh quality metrics, and
	the journal-to-Python conversion workflow.

	Args:
		topic: Topic to document. Options:
			"all"                 - Everything
			"overview"            - Power Tools + Item Wizard + general flow
			"imprint_merge"       - Why imprint+merge is necessary;
			                        curve veance, is_merged, tolerance
			"exodus"              - block / sideset / nodeset: what differs
			"small_features"      - Detect + remove blend chains + re-check
			"gaps_overlaps"       - find overlap + remove overlap workflow
			"quality"             - shape / scaled_jacobian thresholds
			"journal_to_python"   - History tab -> journal -> Python

	Aliases: power_tools, item_wizard, contiguous_mesh, curve_veance,
	         is_merged, block_sideset_nodeset, blend_chain, find_overlap,
	         shape_metric, history_tab.
	"""
	return get_diagnostics_documentation(topic)


# ============================================================
# Self-test
# ============================================================

def _selftest(audit_repo: bool = False):
    """Run fixture self-test; optionally audit durable repository lanes."""
    print("=" * 70)
    print("Cubit Export Lint Self-Test")
    print("=" * 70)
    print()

    # --- Tool annotation classification (MathWorks preset discipline) ---
    n_tools = len(mcp._tool_manager._tools)
    print(f"  tool annotations: {n_tools} tools classified")
    if _UNCLASSIFIED_TOOLS:
        print("  WARNING: tools defaulted to READONLY without an explicit "
              "entry or read-only name shape -- classify them in "
              "_DESTRUCTIVE_TOOLS/_WRITING_TOOLS/_WEB_TOOLS if wrong:")
        for name in _UNCLASSIFIED_TOOLS:
            print(f"    - {name}")
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

    audit_dirs = [
        "docs",
        "validation_test",
        "src/radia/panels/samples",
    ]
    existing_audit_dirs = [d for d in audit_dirs if (PROJECT_ROOT / d).exists()]
    if not existing_audit_dirs:
        if not fixtures_dir.exists():
            print("SKIP: No fixtures or durable audit lanes found")
        print("PASSED")
        return

    if not audit_repo:
        print("  repo audit: SKIPPED (run --selftest --audit-repo)")
        print("PASSED")
        return

    for directory in existing_audit_dirs:
        result = lint_cubit_directory(directory)
        print(result)
    print("PASSED")




register_status_tool(
    mcp,
    server_name='mcp-server-cubit',
    description='Cubit mesh scripting, hex/tet workflow, export formats',
    subpackage='radia_mcp.cubit',
    related_servers=["build123d"],
    audit_command="mcp-server-cubit --selftest --audit-repo",
)


# ============================================================
# Tool annotations + optional gate hiding (MathWorks MCP pattern)
# ============================================================
# Every tool is classified through one of FOUR fully-specified presets
# (MathWorks' annotations.go discipline: presets only, all hint fields
# always set).  Classification is central so adding a tool forces a
# conscious choice here; unclassified tools fall back to READONLY and
# are reported by --selftest.

from ..common.server_hardening import (
	ANN_DESTRUCTIVE as _ANN_DESTRUCTIVE,
	ANN_READONLY as _ANN_READONLY,
	ANN_READONLY_WEB as _ANN_READONLY_WEB,
	ANN_WRITES as _ANN_WRITES,
	classify_tool_annotations as _classify_tool_annotations_common,
	hide_gate_tools as _hide_gate_tools,
	install_call_log as _install_call_log_common,
)

# Tools that execute commands in (or overwrite / stop) the live session.
_DESTRUCTIVE_TOOLS = {
	"cubit_exec", "cubit_exec_safely", "cubit_show", "open_in_cubit",
	"cubit_restore", "cubit_session_shutdown", "cubit_mesh_apply_choice",
	"cubit_mesh_race_with_human",
}
# Tools that create/refresh files on disk but leave the live session alone.
_WRITING_TOOLS = {
	"cubit_checkpoint", "cubit_snapshot", "cubit_batch_try",
	"cubit_mesh_auto", "cubit_mesh_race", "cubit_mesh_race_smart",
	"cubit_mesh_race_smart_async", "cubit_mesh_race_review",
	"cubit_mesh_race_review_async", "cubit_curate_learned_recipes",
	"cubit_examples_refresh", "cubit_check_vol", "cubit_scaffold_toolbar",
	"cubit_session_journal", "cubit_netgen_quality_compare",
	"cubit_stl_to_vol",
}
# Read-only tools that may reach the network.
_WEB_TOOLS = {"cubit_web_docs", "cubit_examples"}


_CUBIT_READONLY_HINTS = (
	"_gate", "_docs", "_guide", "_tips", "_reference", "_inventory",
	"_status", "_lookup", "_ask", "_examples", "lint_", "get_",
	"generate_", "netgen_", "_probe", "_diagnose", "_suggest",
	"_failures", "_checkpoints", "_audit", "_doctor",
)


def _classify_tool_annotations() -> list[str]:
	"""Stamp presets onto every tool via the shared hardening module."""
	return _classify_tool_annotations_common(
		mcp,
		destructive=_DESTRUCTIVE_TOOLS,
		writes=_WRITING_TOOLS,
		web=_WEB_TOOLS,
		readonly_name_hints=_CUBIT_READONLY_HINTS,
	)


_UNCLASSIFIED_TOOLS = _classify_tool_annotations()

# Optional interactive-surface trim (MathWorks exposes 5 tools, not 83):
# the ~40 scenario/CI gate tools are only needed by validation pipelines.
# Set RADIA_MCP_CUBIT_GATES=0 to hide them in interactive sessions.
# Default keeps every tool registered so docs/TOOLS.md stays stable.
_hide_gate_tools(mcp, "RADIA_MCP_CUBIT_GATES")


# ============================================================
# All-calls JSONL session log (MathWorks basetool + slog pattern)
# ============================================================
# One choke point wraps ToolManager.call_tool: every tool call is
# appended as a JSON line {ts, tool, ms, ok, error?} to
# <state_dir>/logs/cubit_tool_calls.jsonl -- making "silently wrong
# mesh" sessions triageable after the fact.  Failures additionally go
# through the existing per-kind failure log.  Disable with
# RADIA_MCP_CUBIT_CALL_LOG=0.

def _install_call_log() -> None:
	"""Wrap call_tool with the shared JSONL all-calls log."""
	_install_call_log_common(mcp, "cubit_tool_calls.jsonl",
	                         "RADIA_MCP_CUBIT_CALL_LOG")


_install_call_log()


def _is_closed_stdout_error(exc: BaseException) -> bool:
    """Windows raises EINVAL when a downstream PowerShell pipe closes early."""
    return isinstance(exc, BrokenPipeError) or (
        isinstance(exc, OSError)
        and getattr(exc, "errno", None) in {errno.EPIPE, errno.EINVAL}
    )


def _maybe_eager_warmup() -> None:
	"""Opt-in eager session start (MathWorks --initialize-matlab-on-
	startup analog): RADIA_CUBIT_EAGER=1 hides the 30+ s license/startup
	cost behind server init instead of the first tool call.  Failures
	are logged to stderr only -- the first real call re-attempts and
	surfaces the error through the normal contract."""
	if os.environ.get("RADIA_CUBIT_EAGER", "0") != "1":
		return
	import threading

	def _warm():
		try:
			_cs.CubitSession.get().ensure_started()
		except Exception as exc:
			print(f"[cubit eager warmup] {type(exc).__name__}: {exc}",
			      file=sys.stderr)

	threading.Thread(target=_warm, name="cubit-eager-warmup",
	                 daemon=True).start()


def _setup_mode() -> int:
	"""One-shot environment preparation (MathWorks --setup-matlab analog):
	pre-warm the Cubit license, then print the full doctor report.
	Exit 0 when the doctor finds no problems, 1 otherwise."""
	print("=" * 70)
	print("mcp-server-cubit --setup")
	print("=" * 70)
	bin_dir = _cs.find_cubit_install()
	if bin_dir is None:
		print("ERROR: Coreform Cubit not found "
		      "(set CUBIT_BIN_DIR / CUBIT_INSTALL_DIR).")
		return 1
	print(f"Cubit install: {bin_dir}")
	try:
		from .license_warmup import warmup_license
		warmup = warmup_license(Path(bin_dir),
		                        timeout_s=_cs.LICENSE_WARMUP_TIMEOUT_S)
		print("License warmup:", json.dumps(warmup, ensure_ascii=False,
		                                    default=str))
	except Exception as exc:
		print(f"License warmup failed: {type(exc).__name__}: {exc}")
	report = json.loads(cubit_doctor())
	print("Doctor report:")
	print(json.dumps(report, ensure_ascii=False, indent=2, default=str))
	if report["status"] == "ok":
		print("SETUP OK")
		return 0
	print("SETUP FOUND PROBLEMS:")
	for p in report["problems"]:
		print(f"  - {p}")
	return 1


def main():
	"""Entry point for mcp-server-cubit command."""
	if '--selftest' in sys.argv[1:]:
		from radia_mcp.common.utf8_stdout import use_utf8_stdout
		use_utf8_stdout()
		try:
			_selftest(audit_repo='--audit-repo' in sys.argv[1:])
		except (BrokenPipeError, OSError) as exc:
			if _is_closed_stdout_error(exc):
				return
			raise
	elif '--setup' in sys.argv[1:]:
		from radia_mcp.common.utf8_stdout import use_utf8_stdout
		use_utf8_stdout()
		sys.exit(_setup_mode())
	else:
		_maybe_eager_warmup()
		mcp.run(transport="stdio")


if __name__ == '__main__':
	main()
