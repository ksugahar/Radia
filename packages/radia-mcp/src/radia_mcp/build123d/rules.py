"""Lint rules for build123d scripts.

Static analysis (no execution) catching common pitfalls that would
otherwise produce cryptic OCCT errors or silent no-ops at runtime.

Each rule is a function `(filepath: str, lines: list[str]) -> list[dict]`
returning zero-or-more findings. Finding dict:
    {"line": int, "severity": "CRITICAL|HIGH|MODERATE|LOW", "rule": str,
     "message": str}

Wired into `lint_build123d_script` / `lint_build123d_directory` MCP
tools in build123d/server.py.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List


# ---------------------------------------------------------------------------
# Rule: context-manager presence
# ---------------------------------------------------------------------------

_BP_CONTEXT_RE = re.compile(r"\bwith\s+BuildPart\s*\(")
_BS_CONTEXT_RE = re.compile(r"\bwith\s+BuildSketch\s*\(")
_BL_CONTEXT_RE = re.compile(r"\bwith\s+BuildLine\s*\(")
_EXTRUDE_RE = re.compile(r"\bextrude\s*\(")
_REVOLVE_RE = re.compile(r"\brevolve\s*\(")
_SWEEP_RE = re.compile(r"\bsweep\s*\(")


def check_missing_buildpart_context(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: `extrude()` / `sweep()` / `revolve()` called outside a
	`with BuildPart() as ...:` block in the Builder API will either
	raise or silently produce an unbound shape. The Algebra API
	doesn't need BuildPart but then uses different syntax entirely.
	"""
	findings: List[Dict] = []
	in_bp = False
	bp_indent: int | None = None
	for i, raw in enumerate(lines, 1):
		line = raw.rstrip("\n")
		stripped = line.strip()
		if stripped.startswith("#"):
			continue
		indent = len(line) - len(line.lstrip())
		if _BP_CONTEXT_RE.search(stripped):
			in_bp = True
			bp_indent = indent
			continue
		if in_bp and bp_indent is not None:
			# Exit BP context when indent drops back to or below BP's line
			if stripped and indent <= bp_indent and not stripped.startswith(("with ", "\t", " ")):
				in_bp = False
				bp_indent = None
		if (_EXTRUDE_RE.search(stripped) or _SWEEP_RE.search(stripped) or
		    _REVOLVE_RE.search(stripped)):
			# Allow algebra form: `part = extrude(sketch, amount)` outside ctx
			# but only when result is assigned (algebra API).
			is_algebra_call = re.match(r"^\s*\w+\s*=\s*(extrude|sweep|revolve)\s*\(",
			                           line) is not None
			if not in_bp and not is_algebra_call:
				findings.append({
					"line": i,
					"severity": "HIGH",
					"rule": "missing-buildpart-context",
					"message": (
						"`extrude/sweep/revolve` called outside "
						"`with BuildPart() as _bp:` context. Either wrap "
						"in a BuildPart context (Builder API) or assign "
						"the result to a variable (Algebra API)."),
				})
	return findings


# ---------------------------------------------------------------------------
# Rule: sweep without path
# ---------------------------------------------------------------------------

_SWEEP_ARGS_RE = re.compile(r"\bsweep\s*\(([^)]*)\)")


def check_sweep_path_arg(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: `sweep(sketch)` without `path=...` relies on a
	BuildPart's pending curve, which is a frequent cause of
	`BRep_API: Make Shape failed` errors when no pending curve exists.
	"""
	findings: List[Dict] = []
	in_bp = False
	for i, raw in enumerate(lines, 1):
		line = raw.rstrip("\n")
		stripped = line.strip()
		if stripped.startswith("#"):
			continue
		if _BP_CONTEXT_RE.search(stripped):
			in_bp = True
			continue
		if stripped and not line.startswith((" ", "\t")):
			in_bp = False
		m = _SWEEP_ARGS_RE.search(stripped)
		if not m:
			continue
		args = m.group(1)
		if "path" not in args and "=" not in args and len(args.strip()) < 80:
			# Positional sweep(sketch) only — probably missing path
			findings.append({
				"line": i,
				"severity": "MODERATE",
				"rule": "sweep-no-path",
				"message": (
					"`sweep(sketch)` with a single positional arg relies "
					"on the BuildPart's pending curve. Prefer explicit "
					"`sweep(sketch, path=...)` to avoid `BRep_API: Make "
					"Shape failed` when no pending curve exists."),
			})
	return findings


# ---------------------------------------------------------------------------
# Rule: missing export
# ---------------------------------------------------------------------------

_EXPORT_RE = re.compile(r"\b(export_step|export_brep|export_stl|export_3mf)\s*\(")


def check_missing_export(filepath: str, lines: List[str]) -> List[Dict]:
	"""LOW: script creates geometry but never exports. Fine for
	interactive exploration; worth flagging for pipeline scripts.
	"""
	findings: List[Dict] = []
	has_geometry = False
	has_export = False
	for raw in lines:
		stripped = raw.strip()
		if stripped.startswith("#"):
			continue
		if (_BP_CONTEXT_RE.search(stripped) or _BS_CONTEXT_RE.search(stripped)
		        or "Part(" in stripped):
			has_geometry = True
		if _EXPORT_RE.search(stripped):
			has_export = True
	if has_geometry and not has_export:
		findings.append({
			"line": 1,
			"severity": "LOW",
			"rule": "missing-export",
			"message": (
				"Script creates geometry but never calls `export_step` "
				"/ `export_brep` / `export_stl`. Add an export if this "
				"is a pipeline script."),
		})
	return findings


# ---------------------------------------------------------------------------
# Rule: Plane argument ambiguity
# ---------------------------------------------------------------------------

_PLANE_ARG_RE = re.compile(r"\bBuildSketch\s*\(\s*([^)]*)\)")


def check_buildsketch_plane_kwarg(filepath: str, lines: List[str]) -> List[Dict]:
	"""MODERATE: `BuildSketch(some_plane_expr)` positional can match
	`workplanes` or `mode` depending on type. Prefer
	`BuildSketch(plane)` when the plane is literal, or
	`BuildSketch(workplanes=...)` when explicit.
	"""
	findings: List[Dict] = []
	for i, raw in enumerate(lines, 1):
		stripped = raw.strip()
		if stripped.startswith("#"):
			continue
		m = _PLANE_ARG_RE.search(stripped)
		if not m:
			continue
		args = m.group(1)
		if not args.strip():
			continue
		# Heuristic: positional that isn't a Plane / Location but is a
		# variable name ≥ 2 chars → might be ambiguous
		if "=" not in args and args.strip().isidentifier() and \
		   args.strip() not in ("Plane.XY", "Plane.XZ", "Plane.YZ"):
			findings.append({
				"line": i,
				"severity": "MODERATE",
				"rule": "buildsketch-ambiguous-arg",
				"message": (
					f"`BuildSketch({args.strip()})` positional arg — "
					"prefer `BuildSketch(workplanes=...)` or explicit "
					"`Plane(...)` to avoid `workplanes` vs `mode` "
					"resolution surprises."),
			})
	return findings


# ---------------------------------------------------------------------------
# Rule: unclosed BuildLine polyline
# ---------------------------------------------------------------------------

_POLYLINE_RE = re.compile(r"\b(Polyline|polyline)\s*\(")
_CLOSE_RE = re.compile(r"\.close\(\)|\bclose\(\)")


def check_unclosed_polyline_for_extrude(filepath: str, lines: List[str]) -> List[Dict]:
	"""HIGH: a Polyline / BuildLine in a BuildSketch is almost always
	meant to be closed before extrude. Missing `close()` produces
	`OCC::BRepPrimAPI: Make Shape failed`.
	"""
	findings: List[Dict] = []
	in_bs = False
	bs_start_line = 0
	saw_polyline = False
	saw_close = False
	saw_extrude = False
	for i, raw in enumerate(lines, 1):
		stripped = raw.strip()
		if stripped.startswith("#"):
			continue
		if _BS_CONTEXT_RE.search(stripped):
			in_bs = True
			bs_start_line = i
			saw_polyline = saw_close = False
			continue
		if in_bs:
			if _POLYLINE_RE.search(stripped):
				saw_polyline = True
			if _CLOSE_RE.search(stripped):
				saw_close = True
			# Heuristic: blank line or dedent ends the with-block
			if not raw.startswith((" ", "\t")) and stripped:
				in_bs = False
				# At end of BS block check
				if saw_polyline and not saw_close:
					findings.append({
						"line": bs_start_line,
						"severity": "HIGH",
						"rule": "polyline-not-closed",
						"message": (
							"Polyline declared in BuildSketch but never "
							"`close()`-d. `extrude()` on an open curve "
							"will fail with BRep_API errors."),
					})
		if _EXTRUDE_RE.search(stripped):
			saw_extrude = True
	return findings


# ---------------------------------------------------------------------------
# Rule: deprecated / renamed API names
# ---------------------------------------------------------------------------

_DEPRECATED = {
	# (bad → good) — kept small, only well-known cases
	"Loc": "Location",
	"make_face": "make_face",  # kept as sentinel — ok
}


def check_deprecated_api(filepath: str, lines: List[str]) -> List[Dict]:
	"""LOW: flag obviously stale API names."""
	findings: List[Dict] = []
	for i, raw in enumerate(lines, 1):
		if raw.strip().startswith("#"):
			continue
		# Old cadquery-ism: `cq.Workplane(` in a build123d file is a smell
		if re.search(r"\bcq\.Workplane\s*\(", raw):
			findings.append({
				"line": i,
				"severity": "LOW",
				"rule": "cadquery-in-build123d",
				"message": (
					"`cq.Workplane(...)` inside a build123d script — "
					"this is CadQuery, not build123d. Use "
					"`execute_cadquery` tool instead, or port to "
					"`BuildPart` / `Part` primitives."),
			})
	return findings


# ---------------------------------------------------------------------------
# Rule: micro-edge radius (CAE quality)
# ---------------------------------------------------------------------------

_FILLET_RE = re.compile(r"\b(fillet|chamfer)\s*\([^)]*\bradius\s*=\s*([0-9.eE+-]+)")


def check_tiny_fillet_radius(filepath: str, lines: List[str]) -> List[Dict]:
	"""LOW: fillet/chamfer with radius < 0.001 mm produces micro-edges
	that break mesh quality. CAE workflows typically want
	`min(edge_length) / char_length ≥ 0.001`.
	"""
	findings: List[Dict] = []
	for i, raw in enumerate(lines, 1):
		if raw.strip().startswith("#"):
			continue
		m = _FILLET_RE.search(raw)
		if not m:
			continue
		try:
			r = float(m.group(2))
		except ValueError:
			continue
		if r < 1e-3 and r > 0:
			findings.append({
				"line": i,
				"severity": "LOW",
				"rule": "micro-fillet-radius",
				"message": (
					f"`{m.group(1)}(radius={r})` — below 1e-3 mm typically "
					"creates micro-edges that wreck downstream mesh "
					"quality. Confirm units or raise the radius."),
			})
	return findings


ALL_RULES = [
	check_missing_buildpart_context,
	check_sweep_path_arg,
	check_missing_export,
	check_buildsketch_plane_kwarg,
	check_unclosed_polyline_for_extrude,
	check_deprecated_api,
	check_tiny_fillet_radius,
]
