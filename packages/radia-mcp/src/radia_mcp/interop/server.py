"""mcp-server-radia-interop — CAD-MCP mesh dispatcher.

**Sugawara Lab (菅原研) stance**: the official primary pair is
**build123d + Cubit**. All other CAD sources (CadQuery, OpenSCAD,
FreeCAD, Blender, Onshape, AutoCAD, KiCad, …) are treated as
compatibility layers routed through this interop module.

New lab work should be authored in **build123d** and meshed via
**Cubit** (`build123d_to_cubit_hex` + `cubit_mesh_auto`). This
module exists only to let legacy scripts and other CAD MCPs reuse
the same backend.

Tools:
  any_step_to_cubit_hex(step_path)          — universal STEP → hex
  openscad_to_cubit_hex(code)               — compat: OpenSCAD CLI
  freecad_to_cubit_hex(script)              — compat: FreeCADCmd
  list_cad_mcp_interop()                    — lists adapters with roles
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .comsol_livelink_knowledge import (
	get_comsol_livelink,
	TOPICS as _COMSOL_LIVELINK_TOPICS,
)
from .comsol_lab_tips_knowledge import (
	get_comsol_lab_tips_documentation,
	TOPICS as _COMSOL_LAB_TIPS_TOPICS,
)

mcp = FastMCP("mcp-server-radia-interop")


def _which(*candidates: str) -> str | None:
	for c in candidates:
		p = shutil.which(c)
		if p:
			return p
	return None


# Common Windows install locations — many CAD tools are installed but not
# added to PATH. Glob-expand and prefer highest-version-looking hit.
_WIN_OPENSCAD_GLOBS = [
	r"C:\Program Files\OpenSCAD\openscad.exe",
	r"C:\Program Files (x86)\OpenSCAD\openscad.exe",
]
_WIN_FREECAD_GLOBS = [
	r"C:\Program Files\FreeCAD*\bin\FreeCADCmd.exe",
	r"C:\Program Files\FreeCAD*\bin\FreeCAD.exe",
	r"C:\Program Files (x86)\FreeCAD*\bin\FreeCADCmd.exe",
]
_POSIX_OPENSCAD_GLOBS = [
	"/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD",
	"/usr/local/bin/openscad",
	"/usr/bin/openscad",
]
_POSIX_FREECAD_GLOBS = [
	"/Applications/FreeCAD.app/Contents/Resources/bin/FreeCADCmd",
	"/usr/local/bin/FreeCADCmd",
	"/usr/bin/FreeCADCmd",
]


def _glob_first(patterns: list[str]) -> str | None:
	import glob as _glob
	for pat in patterns:
		hits = sorted(_glob.glob(pat), reverse=True)  # newest version first
		for h in hits:
			if os.path.exists(h):
				return h
	return None


def _find_openscad() -> str | None:
	"""Locate the OpenSCAD binary. Priority:
	  1. $OPENSCAD env var
	  2. PATH via `shutil.which`
	  3. Platform-specific install globs (Windows Program Files, macOS)
	"""
	env = os.environ.get("OPENSCAD")
	if env and os.path.exists(env):
		return env
	w = _which("openscad", "OpenSCAD")
	if w:
		return w
	if os.name == "nt":
		return _glob_first(_WIN_OPENSCAD_GLOBS)
	return _glob_first(_POSIX_OPENSCAD_GLOBS)


def _find_freecad() -> str | None:
	"""Locate the FreeCADCmd / FreeCAD binary. Priority:
	  1. $FREECAD_CMD env var
	  2. PATH
	  3. Windows `Program Files/FreeCAD*` globs (latest-version wins)
	  4. macOS / Linux default install paths
	"""
	env = os.environ.get("FREECAD_CMD")
	if env and os.path.exists(env):
		return env
	w = _which("FreeCADCmd", "freecadcmd", "freecad-cmd",
	           "FreeCAD", "freecad")
	if w:
		return w
	if os.name == "nt":
		return _glob_first(_WIN_FREECAD_GLOBS)
	return _glob_first(_POSIX_FREECAD_GLOBS)


# ---------------------------------------------------------------------------
# Universal: any STEP → Cubit hex mesh
# ---------------------------------------------------------------------------

@mcp.tool()
def any_step_to_cubit_hex(step_path: str,
                           target_size: float = 1.0,
                           prefer: str = "hex",
                           commit_to_gui: bool = True,
                           timeout_s: int = 600) -> str:
	"""Universal CAD-MCP mesh backend: accept ANY STEP file and run it
	through `cubit_mesh_auto` (batch ladder → winner replayed in live
	Cubit GUI).

	Upstream sources that have been tested:
	  * build123d (via `build123d_to_cubit_hex`)
	  * cadquery (via `cadquery_to_cubit_hex`)
	  * OpenSCAD (via `openscad_to_cubit_hex`)
	  * FreeCAD / FreeCAD-MCP (via `freecad_to_cubit_hex` or any
	    FreeCAD STEP export)
	  * Blender-MCP, Onshape-MCP, AutoCAD-MCP, KiCad-MCP — any of
	    them; feed the exported STEP directly.

	This is the "advertised" backend door for the CAD-MCP ecosystem:
	other MCP authors can call this from their server knowing it
	will give a reliable hex mesh with GUI preview.

	Args:
	    step_path: path to a STEP file from any source.
	    target_size: `volume all size` for the mesh.
	    prefer: "hex" (require hex>0) or "any".
	    commit_to_gui: replay winning recipe in live Cubit GUI.
	    timeout_s: per-rung timeout.
	"""
	if not Path(step_path).exists():
		return json.dumps({"status": "error",
		                   "error": f"STEP not found: {step_path}"})
	try:
		from ..cubit.server import cubit_mesh_auto
	except ImportError as e:
		return json.dumps({"status": "error",
		                   "stage": "cubit_import",
		                   "error": str(e)})
	return cubit_mesh_auto(step_path=step_path.replace("\\", "/"),
	                        target_size=target_size,
	                        prefer=prefer,
	                        commit_to_gui=commit_to_gui,
	                        timeout_s=timeout_s)


# ---------------------------------------------------------------------------
# OpenSCAD adapter (CLI: openscad -o out.step in.scad)
# ---------------------------------------------------------------------------

@mcp.tool()
def openscad_to_cubit_hex(code: str,
                           target_size: float = 1.0,
                           prefer: str = "hex",
                           commit_to_gui: bool = True,
                           timeout_s: int = 600) -> str:
	"""Execute OpenSCAD code, export STEP, run through `cubit_mesh_auto`.

	Requires the `openscad` CLI on PATH (or `OPENSCAD` env var).
	OpenSCAD ≥ 2021.01 supports STEP export via `-o out.step`.

	Args:
	    code: OpenSCAD source (self-contained).
	    target_size / prefer / commit_to_gui / timeout_s: forwarded
	        to `cubit_mesh_auto`.
	"""
	exe = _find_openscad()
	if not exe:
		return json.dumps({
			"status": "error",
			"stage": "openscad_missing",
			"error": ("OpenSCAD CLI not found. Install OpenSCAD ≥ 2021.01 "
			          "(PATH, Windows Program Files, or $OPENSCAD env var)."),
		})

	tmp = Path(tempfile.mkdtemp(prefix="openscad_to_cubit_"))
	scad_path = tmp / "model.scad"
	step_path = tmp / "model.step"
	scad_path.write_text(code, encoding="utf-8")
	try:
		r = subprocess.run(
			[exe, "-o", str(step_path), str(scad_path)],
			capture_output=True, text=True, timeout=timeout_s,
		)
	except subprocess.TimeoutExpired:
		return json.dumps({"status": "error",
		                   "stage": "openscad_timeout",
		                   "error": f"exceeded {timeout_s}s"})
	except OSError as e:
		return json.dumps({"status": "error",
		                   "stage": "openscad_spawn",
		                   "error": str(e)})
	if r.returncode != 0 or not step_path.exists():
		return json.dumps({"status": "error",
		                   "stage": "openscad_export",
		                   "returncode": r.returncode,
		                   "stderr": (r.stderr or "")[-800:],
		                   "stdout": (r.stdout or "")[-400:]})

	# Hand off
	return any_step_to_cubit_hex(step_path=str(step_path),
	                              target_size=target_size,
	                              prefer=prefer,
	                              commit_to_gui=commit_to_gui,
	                              timeout_s=timeout_s)


# ---------------------------------------------------------------------------
# FreeCAD adapter (FreeCADCmd subprocess)
# ---------------------------------------------------------------------------

_FREECAD_WRAPPER = """\
# Auto-generated by radia-mcp freecad_to_cubit_hex.
# User script runs here; whatever ends up in the variable `shape` (or
# the module-level `__result__`) is exported to `{out_step}`.
import sys, traceback, Part
__result__ = None

try:
{indented}
except Exception as e:
    traceback.print_exc()
    sys.exit(1)

# Pick up result: prefer __result__, else `shape`, else last Part.Shape var
result = None
if __result__ is not None:
    result = __result__
elif 'shape' in dir() and hasattr(shape, 'exportStep'):
    result = shape
else:
    for _n in reversed(list(dir())):
        _o = eval(_n)
        if hasattr(_o, 'exportStep'):
            result = _o
            break

if result is None:
    print('NO_SHAPE')
    sys.exit(2)
try:
    result.exportStep(r"{out_step}")
except Exception:
    traceback.print_exc()
    sys.exit(3)
print('WROTE_STEP')
"""


@mcp.tool()
def freecad_to_cubit_hex(script: str,
                          target_size: float = 1.0,
                          prefer: str = "hex",
                          commit_to_gui: bool = True,
                          timeout_s: int = 600) -> str:
	"""Execute a FreeCAD script in a FreeCADCmd subprocess, export the
	resulting shape to STEP, run through `cubit_mesh_auto`.

	Requires `FreeCADCmd` (headless FreeCAD) on PATH (or `FREECAD_CMD`
	env var).  The script runs with `Part` pre-imported.  To indicate
	the shape to export, set `__result__ = your_shape` or leave a
	variable named `shape` in the namespace.

	Example (minimal):
	    ```python
	    box = Part.makeBox(10, 10, 10)
	    __result__ = box
	    ```

	Args:
	    script: FreeCAD Python snippet (runs with `Part` imported).
	    target_size / prefer / commit_to_gui / timeout_s: forwarded
	        to `cubit_mesh_auto`.
	"""
	exe = _find_freecad()
	if not exe:
		return json.dumps({
			"status": "error",
			"stage": "freecad_missing",
			"error": ("FreeCADCmd not found. Install FreeCAD (PATH, "
			          "Windows Program Files, /Applications, or "
			          "$FREECAD_CMD env var)."),
		})

	tmp = Path(tempfile.mkdtemp(prefix="freecad_to_cubit_"))
	step_path = tmp / "model.step"
	wrapper = tmp / "_wrapper.py"
	indented = "\n".join("    " + ln for ln in script.splitlines())
	wrapper.write_text(
		_FREECAD_WRAPPER.format(
			indented=indented,
			out_step=str(step_path).replace("\\", "/"),
		),
		encoding="utf-8",
	)
	try:
		r = subprocess.run(
			[exe, str(wrapper)],
			capture_output=True, text=True, timeout=timeout_s,
		)
	except subprocess.TimeoutExpired:
		return json.dumps({"status": "error",
		                   "stage": "freecad_timeout",
		                   "error": f"exceeded {timeout_s}s"})
	except OSError as e:
		return json.dumps({"status": "error",
		                   "stage": "freecad_spawn",
		                   "error": str(e)})
	if r.returncode != 0 or not step_path.exists():
		return json.dumps({"status": "error",
		                   "stage": "freecad_export",
		                   "returncode": r.returncode,
		                   "stderr": (r.stderr or "")[-800:],
		                   "stdout": (r.stdout or "")[-400:]})

	return any_step_to_cubit_hex(step_path=str(step_path),
	                              target_size=target_size,
	                              prefer=prefer,
	                              commit_to_gui=commit_to_gui,
	                              timeout_s=timeout_s)


# ---------------------------------------------------------------------------
# FreeCAD safety pattern: checkpoint-then-batch-then-commit
# ---------------------------------------------------------------------------

_FREECAD_SAVE_AS = """\
# Save the active document to a checkpoint .FCStd
import sys, FreeCAD
doc = FreeCAD.ActiveDocument
if doc is None:
    print('NO_ACTIVE_DOC')
    sys.exit(2)
doc.saveAs(r"{out}")
print('SAVED_FCSTD')
"""


_FREECAD_TRY_WRAPPER = """\
# Open checkpoint, run user commands, report success
import sys, traceback, FreeCAD, Part
try:
    doc = FreeCAD.openDocument(r"{checkpoint}")
except Exception:
    traceback.print_exc()
    sys.exit(2)

try:
{indented}
    print('__FCMD_TRY_OK__')
except Exception:
    print('__FCMD_TRY_ERR__')
    traceback.print_exc()
    sys.exit(3)
"""


@mcp.tool()
def freecad_exec_safely(commands: str,
                          live_doc_path: str = "",
                          out_apply_path: str = "",
                          timeout_s: int = 300) -> str:
	"""Cubit-style safety pattern for FreeCAD: snapshot → batch dry-run
	→ commit-on-success.  Mirror of `cubit_exec_safely`.

	**Sugawara Lab note**: this tool exists for the FreeCAD community
	and external users — the lab itself authors in **build123d**, not
	FreeCAD.  Lab-internal iteration uses `build123d_try_race` (N
	subprocess variants in parallel, valid+largest-volume winner) for
	parameter sweeps.  `freecad_exec_safely` is the friendly-interop
	port of the safety pattern, not a lab-recommended workflow.

	Workflow:
	  1. **Snapshot**: take the supplied `live_doc_path` (or, if you
	     have neka-nat/freecad-mcp running and `FREECAD_LIVE` env set,
	     ask the live FreeCAD to save its ActiveDocument to a temp
	     `.FCStd` checkpoint via FreeCADCmd → impossible without
	     direct live access, so for now we require `live_doc_path` to
	     point at a saved-on-disk FreeCAD document).
	  2. **Batch dry-run**: spawn FreeCADCmd, open the checkpoint,
	     run `commands` (Python with `FreeCAD` and `Part` pre-imported),
	     check for clean exit.
	  3. **Commit**: on success, save the modified doc to
	     `out_apply_path` (or overwrite `live_doc_path` if not given).
	     On failure: leave the live doc untouched and return the
	     checkpoint path for revert.

	Args:
	    commands: Python source. Use `doc = FreeCAD.ActiveDocument`
	        to access the opened checkpoint; calls like
	        `doc.addObject('Part::Box','b')` mutate it; the wrapper
	        calls `doc.recompute()` and `doc.save()` automatically
	        when commands finish without exception.
	    live_doc_path: existing `.FCStd` representing the "live" state
	        (we treat the file on disk as the live state).
	    out_apply_path: where the modified doc is written on success.
	        Empty = overwrite `live_doc_path`.
	    timeout_s: subprocess timeout.

	Returns JSON with: checkpoint path / batch_passed / live_applied /
	stdout / stderr.
	"""
	exe = _find_freecad()
	if not exe:
		return json.dumps({"status": "error",
		                   "stage": "freecad_missing",
		                   "error": "FreeCADCmd not found"})
	if not live_doc_path:
		return json.dumps({"status": "error",
		                   "error": "live_doc_path required (path to .FCStd)"})
	src = Path(live_doc_path)
	if not src.exists():
		return json.dumps({"status": "error",
		                   "error": f"live_doc not found: {live_doc_path}"})

	tmp = Path(tempfile.mkdtemp(prefix="freecad_safe_"))
	checkpoint = tmp / f"checkpoint_{int(__import__('time').time())}.FCStd"
	# Snapshot: copy the live doc to checkpoint (immutable rollback target)
	shutil.copy2(src, checkpoint)

	# Batch dry-run on the checkpoint
	indented = "\n".join("    " + ln for ln in commands.splitlines())
	wrapper_src = _FREECAD_TRY_WRAPPER.format(
		checkpoint=str(checkpoint).replace("\\", "/"),
		indented=indented,
	) + "\n    doc.recompute()\n    doc.save()\n"
	# Re-build with recompute/save inside the try block
	wrapper_src = (
		"import sys, traceback, FreeCAD, Part\n"
		f"doc = FreeCAD.openDocument(r'{str(checkpoint).replace(chr(92), '/')}')\n"
		"try:\n"
		f"{indented}\n"
		"    doc.recompute()\n"
		"    doc.save()\n"
		"    print('__FCMD_TRY_OK__')\n"
		"except Exception:\n"
		"    print('__FCMD_TRY_ERR__')\n"
		"    traceback.print_exc()\n"
		"    sys.exit(3)\n"
	)
	wrapper_path = tmp / "_try.py"
	wrapper_path.write_text(wrapper_src, encoding="utf-8")

	try:
		r = subprocess.run([exe, str(wrapper_path)],
		                    capture_output=True, text=True, timeout=timeout_s)
	except subprocess.TimeoutExpired:
		return json.dumps({"status": "error",
		                   "stage": "freecad_timeout",
		                   "checkpoint_path": str(checkpoint),
		                   "error": f"exceeded {timeout_s}s"})

	out = r.stdout or ""
	err = r.stderr or ""
	batch_passed = (r.returncode == 0) and ("__FCMD_TRY_OK__" in out)
	payload: dict = {
		"status": "ok",
		"checkpoint_path": str(checkpoint),
		"batch_passed": batch_passed,
		"batch_returncode": r.returncode,
		"batch_stdout": out[-2000:],
		"batch_stderr": err[-2000:],
		"live_applied": False,
	}

	if not batch_passed:
		payload["status"] = "error"
		payload["stage"] = "batch"
		payload["note"] = (f"Batch dry-run failed — live doc untouched. "
		                   f"Inspect batch_stderr; revert via the checkpoint "
		                   f"at {checkpoint}.")
		return json.dumps(payload, indent=2, ensure_ascii=False)

	# Commit: the batch saved the modified file to `checkpoint`. Now
	# move/copy it to `out_apply_path` (or overwrite live_doc_path).
	dest = Path(out_apply_path) if out_apply_path else src
	try:
		shutil.copy2(checkpoint, dest)
	except OSError as e:
		payload["status"] = "error"
		payload["stage"] = "commit"
		payload["error"] = f"copy to {dest}: {e}"
		return json.dumps(payload, indent=2, ensure_ascii=False)

	payload["live_applied"] = True
	payload["live_doc_path"] = str(dest)
	payload["note"] = (f"Batch passed and committed to {dest}. "
	                   f"Rollback by copying {checkpoint} back over it.")
	return json.dumps(payload, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# COMSOL LiveLink knowledge (Java + MATLAB + MPh Python)
# ---------------------------------------------------------------------------

@mcp.tool()
def interop_comsol_livelink(topic: str = "overview") -> str:
	"""COMSOL LiveLink (Java + MATLAB + MPh Python) knowledge.

	The lab uses COMSOL as a reference comparison code; .java exports
	from the GUI are the canonical scripting artifact. This tool returns
	distilled knowledge drawn from S:/COMSOL/ -- the lab's model1.java,
	model2.java, ForJavaScript.java practice files, the COMSOL_compile.ps1
	helper, and the symmetry / cheat-sheet notes in Comsol_Tips.docx.

	Args:
	    topic: one of
	      "overview"          -- when to use Java vs MATLAB vs MPh Python
	      "java_api"          -- com.comsol.model class hierarchy + idioms
	      "matlab_livelink"   -- mli usage from MATLAB
	      "mph_python"        -- MPh Python (lab default for sweeps)
	      "lab_examples"      -- S:/COMSOL/ lab files walk-through
	      "compile_workflow"  -- comsolcompile + COMSOL_compile.ps1 recipe
	      "when_to_use_which" -- decision table + perf / license notes
	      "pitfalls"          -- license, encoding, escaping, version drift
	      "servers"           -- comsolserver / comsolmphserver / comsolbatch
	      "references"        -- lab file paths + COMSOL doc paths
	      "all"               -- everything (concatenated)
	"""
	return get_comsol_livelink(topic)


# ---------------------------------------------------------------------------
# Sugahara Lab COMSOL practical tips (5 yrs of NAS-stashed notes)
# ---------------------------------------------------------------------------

@mcp.tool()
def interop_comsol_lab_tips(topic: str = "overview") -> str:
	"""Sugahara Lab (Kindai University) COMSOL practical tips compendium.

	Five years (2020-2025) of lab-discovered COMSOL Multiphysics
	know-how distilled into ~20 topics. Sources are NAS-stashed plain
	text seed notes, vendor (KESCO Japan) inquiry transcripts,
	bilingual .docx Q&A, .pptx accelerator-magnet setup decks, and a
	2024-03 mesh-type investigation folder.

	Pairs with the COMSOL fork's `docs/LAB_TIPS.md` (kept in sync
	intentionally) -- this MCP surface ships with radia-mcp on PyPI
	so the tips show up in any Sugahara Lab MCP client without
	requiring the COMSOL fork to be installed.

	The biggest single trip-wire covered here: symmetry-plane current
	scaling (topic "coil_current_scaling"). A 1/8 model with three
	symmetry planes can need an 8x scaling on the post-processed
	quantity -- COMSOL does NOT warn -- and the wrong answer looks
	"almost right".

	Args:
		topic: One of (call with "overview" first if unsure):
			"overview"                    -- topic index + meta-stance
			"solver_choice"               -- MUMPS up to 1M DOF
			"iteration_limit"             -- Fully Coupled 25 -> 250 bump
			"acdc_iterative_solvers"      -- GMRES OK, BiCGStab/CG erratic
			"convergence_curve"           -- residual plot reading
			"symmetry_plane_bc"           -- Mag Insul vs PMC decision
			"version_55_unsupported"      -- 5.5 sunset; re-validate old
			"dof_count_check"             -- numberofdofs Help shortcut
			"tokyo_tech_hpc"              -- FNL license + Batch Run
			"coil_current_scaling"        -- CRITICAL symmetry scaling
			"homogenized_multiturn_coil"  -- N*I setup; (3)x(4) trap
			"coil_geometry_analysis"      -- when the analysis step pays
			"infinite_element_domain"     -- accelerator open-bdy recipe
			"element_order"               -- 1/2/3-ji yousou cost vs acc
			"mesh_statistics"             -- Free Tetrahedral Info 1/2
			"stripping_mph_for_support"   -- Save Solutions=OFF for KESCO
			"heat_transfer_bc_dimcheck"   -- ht BC dimensional bug
			"square_coil_mesh_study"      -- 2024-03 ObjectLayer vs hex
			"cubit_to_comsol_handoff"     -- Nastran .bdf gotchas
			"lab_file_index"              -- every source path on S:/COMSOL/
			"all"                         -- everything concatenated
	"""
	return get_comsol_lab_tips_documentation(topic)


# ---------------------------------------------------------------------------
# Introspection
# ---------------------------------------------------------------------------

@mcp.tool()
def list_cad_mcp_interop() -> str:
	"""List registered CAD-MCP interop adapters + their availability.

	Reports which upstream CAD tools are detected on this machine
	(OpenSCAD, FreeCAD) so callers know which one-shot tools will
	work end-to-end vs which will return a structured "not installed"
	error.
	"""
	adapters = [
		{
			"name": "build123d_to_cubit_hex",
			"tool": "gumyr/build123d",
			"requires": "pip install radia-mcp[build123d]",
			"detected": _try_import("build123d"),
			"where": "radia_mcp.build123d.server",
			"role": "PREFERRED — lab standard upstream CAD",
		},
		{
			"name": "any_step_to_cubit_hex",
			"tool": "any source with STEP export",
			"requires": "any STEP file on disk",
			"detected": True,
			"role": "universal fallback for any upstream CAD MCP",
		},
		{
			"name": "cadquery_to_cubit_hex",
			"tool": "CadQuery",
			"requires": "pip install radia-mcp[cadquery]",
			"detected": _try_import("cadquery"),
			"where": "radia_mcp.build123d.server",
			"role": "compat — OCCT sibling, cadquery-mcp interop",
		},
		{
			"name": "openscad_to_cubit_hex",
			"tool": "OpenSCAD ≥ 2021.01",
			"requires": "PATH / Windows Program Files / $OPENSCAD",
			"detected": _find_openscad() or False,
			"role": "compat — legacy OpenSCAD scripts",
		},
		{
			"name": "freecad_to_cubit_hex",
			"tool": "FreeCADCmd (headless FreeCAD)",
			"requires": "PATH / Windows Program Files/FreeCAD* / $FREECAD_CMD",
			"detected": _find_freecad() or False,
			"role": ("friendly — Sugawara Lab respects the FreeCAD "
			         "community; first-class interop path, not first-"
			         "class authoring"),
		},
	]
	# Convert detected-path strings to bool-ish for readability
	for a in adapters:
		d = a["detected"]
		if isinstance(d, str):
			a["detected_path"] = d
			a["detected"] = True
	return json.dumps({
		"status": "ok",
		"backend": "cubit_mesh_auto (batch ladder → live GUI replay)",
		"family": "radia-mcp universal CAD-MCP mesh dispatcher",
		"lab": "Sugawara Lab (菅原研), Kindai University",
		"primary_pair": "build123d (CAD authoring) + Cubit (hex mesh)",
		"preferred_upstream": "build123d (lab standard) — use "
		                       "build123d_to_cubit_hex / build123d_try / "
		                       "generate_build123d_script(pattern='...')",
		"note": ("OpenSCAD / FreeCAD / CadQuery adapters exist for "
		         "interop with legacy scripts or other CAD MCPs, but "
		         "new lab work should be authored in build123d and "
		         "meshed in Cubit."),
		"adapters": adapters,
	}, indent=2, ensure_ascii=False)


def _try_import(module: str) -> bool:
	try:
		__import__(module)
		return True
	except ImportError:
		return False


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------



register_status_tool(
    mcp,
    server_name='mcp-server-radia-interop',
    description='Cross-CAD interop (STEP/IGES/CadQuery <-> Cubit/Netgen)',
    subpackage='radia_mcp.interop',
    related_servers=["cubit", "build123d"],
)


def main():
	if "--selftest" in sys.argv:
		from radia_mcp.common.utf8_stdout import use_utf8_stdout
		use_utf8_stdout()
		print("radia-interop MCP server self-test:")
		print(list_cad_mcp_interop())
		return
	mcp.run(transport="stdio")


if __name__ == "__main__":
	main()
