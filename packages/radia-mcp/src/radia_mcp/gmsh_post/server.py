"""mcp-server-gmsh-post — MSH v4.1 post-processing tools.

Policy: lab standardises on MSH v4.1. Anything older gets lifted here.

Tools:
  gmsh_post_inspect(msh_path)
  gmsh_post_validate(msh_path)
  gmsh_post_convert(in_path, out_path, binary=False)    # → v4.1
  gmsh_post_write_node_data(msh_path, values, out_path, view_name,
                             time=0.0, time_step=0)
  gmsh_post_write_element_data(msh_path, values, out_path, view_name,
                                time=0.0, time_step=0)
  gmsh_post_spec(topic="overview")
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP

from .spec import get_spec as _get_spec, SPEC as _SPEC

# Late-bind the auto-generated API reference into the spec dict. Kept in
# a separate module so the wheel ships a static snapshot we don't have
# to regenerate at server start.
try:
	from .gmsh_api import get_api_reference as _get_api_reference
	_SPEC["api_reference"] = _get_api_reference()
except Exception:
	# Fresh install before _gen_api_reference has been run — leave the
	# placeholder None so get_spec reports "unknown topic" cleanly.
	pass

mcp = FastMCP("mcp-server-gmsh-post")


# ---------------------------------------------------------------------------
# gmsh import — lazy so the server loads even if the user hasn't installed
# the gmsh extra. Tools that need gmsh return a structured error instead.
# ---------------------------------------------------------------------------

def _ensure_gmsh():
	try:
		import gmsh  # noqa: F401
		return gmsh, None
	except ImportError:
		return None, ("gmsh Python package not installed. Install with "
		              "`pip install radia-mcp[gmsh]` or `pip install gmsh`.")


class _GmshSession:
	"""Context manager that initializes / finalizes gmsh cleanly.

	Gmsh's global state makes nested use tricky — we init at entry,
	finalize at exit, and bubble errors as Python exceptions.
	"""

	def __init__(self):
		self.gmsh = None

	def __enter__(self):
		import gmsh
		gmsh.initialize([], False)
		gmsh.option.setNumber("General.Terminal", 0)
		self.gmsh = gmsh
		return gmsh

	def __exit__(self, exc_type, exc, tb):
		if self.gmsh is not None:
			try:
				self.gmsh.finalize()
			except Exception:
				pass
		return False  # propagate


# ---------------------------------------------------------------------------
# msh header sniffer (stdlib only — works without gmsh installed)
# ---------------------------------------------------------------------------

def _read_msh_header(path: Path, max_lines: int = 300) -> dict:
	"""Scan the first ~300 lines of a .msh to infer version + sections.

	Intentionally stdlib-only: callers can ask "is this v4.1?" without
	paying for gmsh initialization, and the inspect tool uses this as a
	fast-path summary before (optionally) delegating to the gmsh API.
	"""
	info: dict[str, Any] = {
		"version": None,
		"file_type": None,          # 0 ASCII, 1 binary
		"data_size": None,
		"sections": [],
		"first_error": None,
	}
	try:
		with path.open("r", encoding="utf-8", errors="replace") as f:
			seen_format = False
			for i, line in enumerate(f):
				if i >= max_lines and seen_format:
					break
				s = line.strip()
				if s.startswith("$End"):
					continue
				if s.startswith("$"):
					info["sections"].append(s[1:])
					continue
				if not seen_format and info["sections"] == ["MeshFormat"]:
					parts = s.split()
					if len(parts) >= 3:
						info["version"] = parts[0]
						try:
							info["file_type"] = int(parts[1])
							info["data_size"] = int(parts[2])
						except ValueError:
							pass
					seen_format = True
	except OSError as e:
		info["first_error"] = f"read error: {e}"
	return info


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@mcp.tool()
def gmsh_post_inspect(msh_path: str) -> str:
	"""Structured summary of a .msh file.

	Reports: format version / binary-or-ASCII, sections present,
	node / element counts, entity block counts, physical names, and
	any post-processing view blocks ($NodeData / $ElementData).

	Uses a stdlib header scan first (fast, always works) and then,
	if gmsh is available, adds topology-level counts via `gmsh.open`.

	Args:
	    msh_path: path to the .msh file (absolute or project-relative).
	"""
	path = Path(msh_path)
	if not path.exists():
		return json.dumps({"status": "error",
		                   "error": f"file not found: {msh_path}"})

	hdr = _read_msh_header(path)
	result: dict[str, Any] = {
		"status": "ok",
		"path": str(path.resolve()),
		"size_bytes": path.stat().st_size,
		"format_version": hdr["version"],
		"binary": (hdr["file_type"] == 1) if hdr["file_type"] is not None else None,
		"sections": hdr["sections"],
	}

	gmsh, err = _ensure_gmsh()
	if err:
		result["gmsh_detail"] = err
		return json.dumps(result, indent=2)

	try:
		with _GmshSession() as g:
			g.open(str(path))
			# Physical groups
			pgs = []
			for dim, tag in g.model.getPhysicalGroups():
				try:
					name = g.model.getPhysicalName(dim, tag)
				except Exception:
					name = ""
				try:
					ents = g.model.getEntitiesForPhysicalGroup(dim, tag)
				except Exception:
					ents = []
				pgs.append({"dim": int(dim), "tag": int(tag),
				            "name": name,
				            "entities": [int(e) for e in ents]})
			# Entities
			ents = g.model.getEntities()
			ent_by_dim: dict[int, int] = {0: 0, 1: 0, 2: 0, 3: 0}
			for dim, _ in ents:
				ent_by_dim[dim] = ent_by_dim.get(dim, 0) + 1
			# Mesh counts
			node_tags, _, _ = g.model.mesh.getNodes()
			elem_types, elem_tags, _ = g.model.mesh.getElements()
			elem_counts: dict[str, int] = {}
			for et, etags in zip(elem_types, elem_tags):
				try:
					name, _dim, _order, _nn, _coords, _ = g.model.mesh.getElementProperties(et)
				except Exception:
					name = f"type{et}"
				elem_counts[f"{name} (type {et})"] = len(etags)

			# View blocks (post-processing)
			view_tags = list(g.view.getTags()) if hasattr(g, "view") else []
			views = []
			for vt in view_tags:
				try:
					vname = g.view.option.getString(int(vt), "Name")
				except Exception:
					vname = f"view{int(vt)}"
				views.append({"tag": int(vt), "name": vname})

			result["entities"] = {"by_dim": {int(k): int(v)
			                                  for k, v in ent_by_dim.items()},
			                      "total": int(len(ents))}
			result["physical_groups"] = pgs
			result["node_count"] = len(node_tags)
			result["element_counts"] = elem_counts
			result["total_elements"] = sum(elem_counts.values())
			result["post_views"] = views
	except Exception as e:  # noqa: BLE001
		result["gmsh_error"] = f"{type(e).__name__}: {e}"

	return json.dumps(result, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_validate(msh_path: str) -> str:
	"""Check a .msh against the lab v4.1 compliance checklist.

	Applies rules from the bundled spec (see `gmsh_post_spec
	topic='compat_checks'`):
	  - $MeshFormat is 4.1
	  - $Entities present
	  - $PhysicalNames covers every element dim
	  - No duplicate nodes (spatial check)
	  - Consistent element order across $Elements blocks

	Returns structured JSON: severity levels (error / warn / ok),
	one entry per rule, plus a `passed_all` flag.
	"""
	path = Path(msh_path)
	if not path.exists():
		return json.dumps({"status": "error",
		                   "error": f"file not found: {msh_path}"})

	findings: list[dict[str, Any]] = []

	def add(severity: str, rule: str, msg: str, **extra):
		f = {"severity": severity, "rule": rule, "message": msg}
		f.update(extra)
		findings.append(f)

	hdr = _read_msh_header(path)
	# 1. version
	v = hdr.get("version")
	if v != "4.1":
		add("error", "version",
		    f"$MeshFormat is {v!r}, expected '4.1'. Lab policy is v4.1-only.",
		    actual=v, expected="4.1")
	else:
		add("ok", "version", "v4.1")

	sections = hdr.get("sections", [])
	# 2. $MeshFormat first
	if sections and sections[0] != "MeshFormat":
		add("error", "meshformat_first",
		    f"First section is {sections[0]!r}; $MeshFormat must be first.")
	# 3. $Entities present
	if "Entities" not in sections:
		add("error", "entities_present",
		    "$Entities missing — ngsolve / radia will flatten topology.")
	else:
		add("ok", "entities_present", "$Entities section found")
	# 4. $PhysicalNames recommended
	if "PhysicalNames" not in sections:
		add("warn", "physical_names",
		    "$PhysicalNames absent; downstream regions will be anonymous.")

	gmsh, err = _ensure_gmsh()
	if err:
		add("warn", "gmsh_dep", err)
	else:
		try:
			with _GmshSession() as g:
				g.open(str(path))
				# 5. $PhysicalNames covers every element dim
				present_dims: set[int] = set()
				for et, _tags in zip(*g.model.mesh.getElements()[0:2]):
					try:
						_, d, *_ = g.model.mesh.getElementProperties(et)
						present_dims.add(d)
					except Exception:
						pass
				pg_dims = {dim for dim, _ in g.model.getPhysicalGroups()}
				missing_dims = present_dims - pg_dims
				if missing_dims:
					add("warn", "physicalnames_coverage",
					    f"Element dims without PhysicalNames: "
					    f"{sorted(missing_dims)}; named-group lookup breaks.")
				else:
					add("ok", "physicalnames_coverage",
					    f"PhysicalNames cover all element dims ({sorted(present_dims)})")
				# 6. duplicate nodes (spatial dedup)
				node_tags, coords_flat, _ = g.model.mesh.getNodes()
				if len(node_tags) > 0:
					import itertools
					rounded = [
						(round(coords_flat[3 * i], 10),
						 round(coords_flat[3 * i + 1], 10),
						 round(coords_flat[3 * i + 2], 10))
						for i in range(len(node_tags))
					]
					dup_count = len(rounded) - len(set(rounded))
					if dup_count > 0:
						add("warn", "duplicate_nodes",
						    f"{dup_count} coincident nodes detected "
						    f"(likely from imprint/merge). Fix with "
						    f"`gmsh.model.mesh.removeDuplicateNodes()`.")
					else:
						add("ok", "duplicate_nodes", f"no duplicates "
						    f"({len(node_tags)} unique)")
				# 7. element order consistency
				orders: set[int] = set()
				for et in g.model.mesh.getElements()[0]:
					try:
						_, _, order, *_ = g.model.mesh.getElementProperties(et)
						orders.add(order)
					except Exception:
						pass
				if len(orders) > 1:
					add("warn", "element_order_mixed",
					    f"Multiple element orders in $Elements: {sorted(orders)}. "
					    f"radia_ngsolve high-order assembly expects a single order.")
				elif orders:
					add("ok", "element_order_uniform",
					    f"all elements order {orders.pop()}")
		except Exception as e:  # noqa: BLE001
			add("error", "gmsh_open", f"{type(e).__name__}: {e}")

	errors = [f for f in findings if f["severity"] == "error"]
	warns = [f for f in findings if f["severity"] == "warn"]
	return json.dumps({
		"status": "ok" if not errors else "fail",
		"path": str(path.resolve()),
		"passed_all": not errors and not warns,
		"errors": len(errors),
		"warnings": len(warns),
		"findings": findings,
	}, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_convert(in_path: str, out_path: str,
                      binary: bool = False) -> str:
	"""Rewrite a .msh as MSH v4.1 (ASCII by default).

	Intended as a *boundary lift*: anything older than v4.1 on the
	input becomes v4.1 on the output, idempotent on already-4.1 files.
	Uses gmsh internally — preserves $PhysicalNames, lifts flat v2.2
	into entity-grouped v4.1 (gmsh synthesises entities per physical
	group when source has none).

	Args:
	    in_path: path to source .msh (any gmsh-readable version).
	    out_path: where to write the v4.1 result.
	    binary: True for binary encoding (~2-3× smaller for huge meshes).
	"""
	gmsh, err = _ensure_gmsh()
	if err:
		return json.dumps({"status": "error", "error": err})

	inp = Path(in_path)
	out = Path(out_path)
	if not inp.exists():
		return json.dumps({"status": "error",
		                   "error": f"input not found: {in_path}"})
	out.parent.mkdir(parents=True, exist_ok=True)

	try:
		with _GmshSession() as g:
			g.open(str(inp))
			g.option.setNumber("Mesh.MshFileVersion", 4.1)
			g.option.setNumber("Mesh.Binary", 1 if binary else 0)
			g.write(str(out))
	except Exception as e:  # noqa: BLE001
		return json.dumps({"status": "error",
		                   "error": f"{type(e).__name__}: {e}"})

	hdr = _read_msh_header(out)
	return json.dumps({
		"status": "ok",
		"input": str(inp.resolve()),
		"output": str(out.resolve()),
		"input_size": inp.stat().st_size,
		"output_size": out.stat().st_size,
		"output_version": hdr.get("version"),
		"binary": binary,
	}, indent=2, ensure_ascii=False)


def _append_data_block(target: Path, view_msh: Path, data_type: str) -> None:
	"""Extract `$<data_type>...$End<data_type>` from `view_msh` (a
	gmsh-written view file) and append it to `target` (an already-written
	v4.1 mesh).

	Why: `gmsh.view.write` writes a standalone .msh with its own
	$MeshFormat header. Appending the file verbatim produces a doubled
	header that breaks v4.1 spec compliance. We only want the data
	block.
	"""
	start_tag = f"$${data_type}".replace("$$", "$")  # "$NodeData"
	end_tag = f"$End{data_type}"
	with view_msh.open("r", encoding="utf-8", errors="replace") as f:
		lines = f.readlines()
	block: list[str] = []
	in_block = False
	for line in lines:
		s = line.rstrip("\n")
		if not in_block and s == start_tag:
			in_block = True
		if in_block:
			block.append(line)
			if s == end_tag:
				break
	if not block:
		return
	# Ensure target ends with a newline before appending
	with target.open("rb") as f:
		f.seek(-1, 2)
		last = f.read(1)
	with target.open("a", encoding="utf-8") as f:
		if last != b"\n":
			f.write("\n")
		f.writelines(block)


def _parse_values_json(values_json: Any) -> dict[int, list[float]]:
	"""Normalize `values_json` into {int tag: [float, ...]}.

	Accepted inputs:
	  - dict mapping tag → scalar or list
	  - JSON string of such a dict
	"""
	if isinstance(values_json, str):
		values_json = json.loads(values_json)
	if not isinstance(values_json, dict):
		raise ValueError("values must be a dict {tag: scalar|list}")
	out: dict[int, list[float]] = {}
	for k, v in values_json.items():
		tag = int(k)
		if isinstance(v, (int, float)):
			out[tag] = [float(v)]
		elif isinstance(v, (list, tuple)):
			out[tag] = [float(x) for x in v]
		else:
			raise ValueError(f"value for tag {tag} not numeric: {v!r}")
	return out


def _write_model_data(msh_in: str, msh_out: str, view_name: str,
                      values: dict[int, list[float]], data_type: str,
                      time: float, time_step: int) -> dict:
	"""Shared helper for gmsh_post_write_{node,element}_data."""
	gmsh, err = _ensure_gmsh()
	if err:
		return {"status": "error", "error": err}
	inp = Path(msh_in)
	out = Path(msh_out)
	if not inp.exists():
		return {"status": "error", "error": f"input not found: {msh_in}"}
	if not values:
		return {"status": "error", "error": "values is empty"}
	out.parent.mkdir(parents=True, exist_ok=True)

	try:
		with _GmshSession() as g:
			g.open(str(inp))
			model = g.model.getCurrent()
			vtag = g.view.add(view_name)
			tags = list(values.keys())
			data = list(values.values())
			g.view.addModelData(
				vtag, time_step, model, data_type,
				tags, data, time=time,
				numComponents=len(next(iter(values.values()))),
			)
			g.option.setNumber("Mesh.MshFileVersion", 4.1)
			g.option.setNumber("Mesh.Binary", 0)
			# 1. write the mesh (no view) to the output path
			g.write(str(out))
			# 2. write the view to a temp file (standalone .msh with
			#    $MeshFormat + $NodeData/$ElementData only)
			import tempfile as _tmp
			tmp_view = Path(_tmp.mkdtemp(prefix="gmsh_view_")) / "view.msh"
			g.view.write(int(vtag), str(tmp_view))
			# 3. extract the data block from the temp and append to out
			_append_data_block(Path(out), tmp_view, data_type)
			try:
				import shutil as _sh
				_sh.rmtree(tmp_view.parent, ignore_errors=True)
			except Exception:
				pass
	except Exception as e:  # noqa: BLE001
		return {"status": "error", "error": f"{type(e).__name__}: {e}"}

	return {
		"status": "ok",
		"input": str(inp.resolve()),
		"output": str(out.resolve()),
		"view_name": view_name,
		"data_type": data_type,
		"values_written": len(values),
		"num_components": len(next(iter(values.values()))),
		"time": time,
		"time_step": time_step,
	}


@mcp.tool()
def gmsh_post_write_node_data(msh_path: str, out_path: str,
                              view_name: str,
                              values: Any,
                              time: float = 0.0,
                              time_step: int = 0) -> str:
	"""Append a `$NodeData` view to a v4.1 mesh and write the result.

	Args:
	    msh_path: existing .msh (v4.1 or will be lifted).
	    out_path: output path (gets mesh + $NodeData block).
	    view_name: view label that appears in Gmsh GUI.
	    values: dict `{node_tag: scalar | [vx, vy, vz]}` or a JSON
	        string of the same. Component count auto-detected from
	        the first value.
	    time: view time (default 0.0).
	    time_step: view time step (default 0).

	Example — scalar B magnitude at nodes:
	    gmsh_post_write_node_data(
	        msh_path="mesh.msh", out_path="mesh_with_B.msh",
	        view_name="B_magnitude",
	        values={1: 0.0, 2: 0.1, 3: 0.2, ...}
	    )
	"""
	try:
		vals = _parse_values_json(values)
	except (ValueError, json.JSONDecodeError) as e:
		return json.dumps({"status": "error", "error": str(e)})
	r = _write_model_data(msh_path, out_path, view_name, vals,
	                      data_type="NodeData", time=time,
	                      time_step=time_step)
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_write_element_data(msh_path: str, out_path: str,
                                 view_name: str,
                                 values: Any,
                                 time: float = 0.0,
                                 time_step: int = 0) -> str:
	"""Append a `$ElementData` view to a v4.1 mesh and write the result.

	Args:
	    msh_path: existing .msh.
	    out_path: output.
	    view_name: view label.
	    values: dict `{element_tag: scalar | vector}` or JSON string.
	    time / time_step: view coordinates.
	"""
	try:
		vals = _parse_values_json(values)
	except (ValueError, json.JSONDecodeError) as e:
		return json.dumps({"status": "error", "error": str(e)})
	r = _write_model_data(msh_path, out_path, view_name, vals,
	                      data_type="ElementData", time=time,
	                      time_step=time_step)
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_spec(topic: str = "overview") -> str:
	"""Return a section of the bundled MSH v4.1 spec.

	Topics: overview / meshformat / entities / nodes / elements /
	physicalnames / nodedata / write_node_data / compat_checks /
	cubit_to_v41 / view_data_cookbook / physical_groups_cookbook /
	api_reference. Use `topic="all"` for the whole thing.
	"""
	return _get_spec(topic)


# ---------------------------------------------------------------------------
# gmsh_post_api — focused tf-idf search over the auto-generated API ref
# ---------------------------------------------------------------------------

import math as _math
import re as _re

_WORD_RE = _re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*")
_STOP = frozenset({"the", "a", "an", "of", "to", "in", "and", "or", "is",
                    "are", "be", "on", "at", "for", "with", "as"})


def _tokenize(s: str) -> list[str]:
	out: list[str] = []
	for m in _WORD_RE.finditer(s):
		w = m.group(0).lower()
		out.append(w)
		if "_" in w:
			out.extend(p for p in w.split("_") if p)
	return out


def _chunk_by_heading(text: str, target_lines: int = 20) -> list[tuple[int, str]]:
	lines = text.splitlines()
	chunks: list[tuple[int, str]] = []
	cur_start = 0
	cur: list[str] = []
	for i, line in enumerate(lines):
		is_h = (line.startswith("# ") or line.startswith("## ")
		        or line.startswith("### "))
		if is_h and cur and (i - cur_start) >= max(3, target_lines // 4):
			chunks.append((cur_start + 1, "\n".join(cur)))
			cur = []
			cur_start = i
		cur.append(line)
		if len(cur) >= target_lines * 2:
			chunks.append((cur_start + 1, "\n".join(cur)))
			cur = []
			cur_start = i + 1
	if cur:
		chunks.append((cur_start + 1, "\n".join(cur)))
	return chunks


@mcp.tool()
def gmsh_post_api(query: str, max_results: int = 5) -> str:
	"""Search the bundled gmsh API reference (auto-generated via
	`inspect.getmembers(gmsh)` + submodules — ~650 function entries
	across `gmsh.model.mesh`, `gmsh.view`, `gmsh.option`, etc.).

	Use when you know the user is asking about a specific gmsh
	function (`gmsh.view.addModelData`, `gmsh.model.mesh.getNodes`,
	…) and want the exact signature plus its one-line docstring.

	Args:
	    query: function name fragment or keywords.
	    max_results: cap on chunks returned.
	"""
	q = (query or "").strip()
	if not q:
		return json.dumps({"status": "error", "error": "empty query"})
	api = _SPEC.get("api_reference") or ""
	if not api:
		return json.dumps({"status": "error",
		                   "error": ("api_reference topic unavailable — "
		                              "run `python -m radia_mcp.gmsh_post._gen_api_reference` "
		                              "to regenerate.")})
	terms = [t for t in _tokenize(q) if t not in _STOP]
	if not terms:
		return json.dumps({"status": "error",
		                   "error": "no searchable terms"})
	chunks = _chunk_by_heading(api, target_lines=12)
	docs = []
	for start, text in chunks:
		tokens = _tokenize(text)
		if not tokens:
			continue
		docs.append({"start_line": start, "chunk": text,
		             "heading": text.splitlines()[0][:120] if text else "",
		             "tokens": tokens, "token_set": set(tokens)})
	n = max(1, len(docs))
	idf = {
		t: _math.log((n + 1) / (sum(1 for d in docs if t in d["token_set"]) + 1)) + 1.0
		for t in terms
	}
	scored = []
	for doc in docs:
		tf = {}
		for tok in doc["tokens"]:
			tf[tok] = tf.get(tok, 0) + 1
		inv_len = 1.0 / (1.0 + _math.log(1 + len(doc["tokens"])))
		score = 0.0
		hits = 0
		for term in terms:
			f = tf.get(term, 0)
			if f == 0:
				continue
			hits += 1
			score += (1.0 + _math.log(f)) * idf[term] * inv_len
		if hits == 0:
			continue
		# Heading boost
		head_lc = doc["heading"].lower()
		head_hits = sum(1 for t in terms if t in head_lc)
		if head_hits:
			score *= 1.0 + 0.75 * head_hits
		score *= 1.0 + 0.25 * (hits / max(1, len(terms)))
		scored.append((score, doc))
	scored.sort(key=lambda x: (-x[0], x[1]["start_line"]))
	top = scored[: max(1, int(max_results))]
	return json.dumps({
		"status": "ok",
		"query": q,
		"terms": terms,
		"total_scored": len(scored),
		"returned": len(top),
		"results": [
			{"start_line": d["start_line"], "heading": d["heading"],
			 "score": round(s, 3), "chunk": d["chunk"]}
			for s, d in top
		],
	}, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Post-processing tools: quality / extract_physical / boundary / add_view_csv
# ---------------------------------------------------------------------------

@mcp.tool()
def gmsh_post_quality(msh_path: str,
                       quality_measure: str = "minSJ") -> str:
	"""Compute per-element mesh quality + summary statistics.

	Wraps `gmsh.model.mesh.getElementQualities`. Returns a histogram
	(8 bins) + min / max / median / count below threshold so the
	caller can gate meshes with obvious degeneracy before sending
	them to a solver.

	Args:
	    msh_path: v4.1 mesh file.
	    quality_measure: `minSJ` (minimum scaled Jacobian, default),
	        `maxAngle`, `innerRadius`, `volume`, `minIsotropy`,
	        `gamma`, `minEdge`, `maxEdge`, `minSICN`, `minSIGE`.
	        See gmsh docs for full list.

	Returns JSON: {count, min, max, median, below_0.2_pct, histogram}.
	"""
	path = Path(msh_path)
	if not path.exists():
		return json.dumps({"status": "error",
		                   "error": f"not found: {msh_path}"})
	gmsh, err = _ensure_gmsh()
	if err:
		return json.dumps({"status": "error", "error": err})
	try:
		with _GmshSession() as g:
			g.open(str(path))
			# Collect all element tags
			elem_types, elem_tags, _ = g.model.mesh.getElements()
			all_tags: list[int] = []
			for tags in elem_tags:
				all_tags.extend(int(t) for t in tags)
			if not all_tags:
				return json.dumps({"status": "ok",
				                   "path": str(path.resolve()),
				                   "count": 0,
				                   "note": "no elements in mesh"})
			try:
				q_raw = g.model.mesh.getElementQualities(
					elementTags=all_tags, qualityName=quality_measure,
				)
			except Exception as e:
				return json.dumps({"status": "error",
				                   "stage": "getElementQualities",
				                   "error": f"{type(e).__name__}: {e}"})
			q_vals = [float(v) for v in q_raw]
	except Exception as e:  # noqa: BLE001
		return json.dumps({"status": "error",
		                   "error": f"{type(e).__name__}: {e}"})

	if not q_vals:
		return json.dumps({"status": "ok", "count": 0,
		                   "note": "no quality values returned"})
	qmin, qmax = min(q_vals), max(q_vals)
	srt = sorted(q_vals)
	median = srt[len(srt) // 2]
	below = sum(1 for v in q_vals if v < 0.2)
	# 8-bin histogram between min and max (linear)
	nbins = 8
	rng = qmax - qmin
	hist = [0] * nbins
	if rng > 0:
		for v in q_vals:
			b = min(int((v - qmin) / rng * nbins), nbins - 1)
			hist[b] += 1
	else:
		hist[0] = len(q_vals)
	return json.dumps({
		"status": "ok",
		"path": str(path.resolve()),
		"quality_measure": quality_measure,
		"count": len(q_vals),
		"min": round(qmin, 6),
		"max": round(qmax, 6),
		"median": round(median, 6),
		"mean": round(sum(q_vals) / len(q_vals), 6),
		"below_0.2_count": below,
		"below_0.2_pct": round(100 * below / len(q_vals), 2),
		"histogram": hist,
		"histogram_bins": [
			round(qmin + (qmax - qmin) * i / nbins, 4)
			for i in range(nbins + 1)
		],
		"recommendation": (
			"ready for solver" if below == 0
			else f"{below} elements below q=0.2 — consider remesh / "
			     "refine / quality-fix pass"
		),
	}, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_extract_physical(msh_path: str,
                                 group_name: str,
                                 out_path: str) -> str:
	"""Extract the mesh subset belonging to a named physical group.

	Loads `msh_path`, finds the group matching `group_name` (by
	`$PhysicalNames`), identifies all entities in that group + their
	elements, and writes a fresh v4.1 mesh containing only those
	elements / their nodes. Useful for per-region debug or feeding
	a single-material solver.

	Args:
	    msh_path: source v4.1 mesh.
	    group_name: name as declared in `$PhysicalNames`.
	    out_path: destination.
	"""
	src = Path(msh_path)
	if not src.exists():
		return json.dumps({"status": "error",
		                   "error": f"not found: {msh_path}"})
	out = Path(out_path)
	out.parent.mkdir(parents=True, exist_ok=True)
	gmsh, err = _ensure_gmsh()
	if err:
		return json.dumps({"status": "error", "error": err})
	try:
		with _GmshSession() as g:
			g.open(str(src))
			# Find the group by name across all dims
			target = None
			for dim, tag in g.model.getPhysicalGroups():
				try:
					name = g.model.getPhysicalName(dim, tag)
				except Exception:
					name = ""
				if name == group_name:
					target = (int(dim), int(tag))
					break
			if target is None:
				# Enumerate for diagnostic
				known = []
				for dim, tag in g.model.getPhysicalGroups():
					try:
						known.append(g.model.getPhysicalName(dim, tag))
					except Exception:
						pass
				return json.dumps({"status": "error",
				                   "error": f"no physical group named "
				                             f"{group_name!r}",
				                   "known_groups": known})
			dim, tag = target
			ents = list(g.model.getEntitiesForPhysicalGroup(dim, tag))
			# Delete every element not belonging to these entities, then write
			all_ents = list(g.model.getEntities(dim))
			keep = {int(e) for e in ents}
			to_remove = [(dim, int(e)) for _, e in all_ents
			             if int(e) not in keep]
			if to_remove:
				g.model.removeEntities(to_remove, recursive=False)
			g.model.occ.synchronize()
			g.option.setNumber("Mesh.MshFileVersion", 4.1)
			g.option.setNumber("Mesh.Binary", 0)
			g.write(str(out))
	except Exception as e:  # noqa: BLE001
		return json.dumps({"status": "error",
		                   "error": f"{type(e).__name__}: {e}"})
	return json.dumps({
		"status": "ok",
		"input": str(src.resolve()),
		"output": str(out.resolve()),
		"group_name": group_name,
		"dim": dim,
		"tag": tag,
		"entities_kept": len(ents),
		"input_size": src.stat().st_size,
		"output_size": out.stat().st_size,
	}, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_boundary(msh_path: str, out_path: str) -> str:
	"""Extract the boundary (dim-1 facets) of a volume mesh and write
	it as a standalone v4.1 mesh.

	Handy for solvers that want surface meshes separately (BEM,
	radiation exchange, radia ngbem), or for visual sanity: does my
	exterior look right?

	Args:
	    msh_path: source v4.1 mesh (must contain dim-3 elements).
	    out_path: destination v4.1 surface mesh.
	"""
	src = Path(msh_path)
	if not src.exists():
		return json.dumps({"status": "error",
		                   "error": f"not found: {msh_path}"})
	out = Path(out_path)
	out.parent.mkdir(parents=True, exist_ok=True)
	gmsh, err = _ensure_gmsh()
	if err:
		return json.dumps({"status": "error", "error": err})
	try:
		with _GmshSession() as g:
			g.open(str(src))
			# Classify boundary: sum of getBoundary over all dim-3 entities
			vol_ents = g.model.getEntities(3)
			if not vol_ents:
				return json.dumps({"status": "error",
				                   "error": "no dim-3 entities in mesh"})
			boundary = g.model.getBoundary(
				dimTags=vol_ents, combined=True, oriented=False, recursive=False,
			)
			surface_tags = [int(t) for d, t in boundary if d == 2]
			# Remove everything that isn't on the boundary surfaces
			# (this is the blunt-hammer approach; gmsh does fine here)
			all_s2 = {int(t) for _, t in g.model.getEntities(2)}
			drop2 = [(2, t) for t in (all_s2 - set(surface_tags))]
			# Drop volumes too — we want a 2D-only output
			drop3 = [(3, int(t)) for _, t in vol_ents]
			if drop2 or drop3:
				g.model.removeEntities(drop2 + drop3, recursive=False)
			g.option.setNumber("Mesh.MshFileVersion", 4.1)
			g.option.setNumber("Mesh.Binary", 0)
			g.write(str(out))
	except Exception as e:  # noqa: BLE001
		return json.dumps({"status": "error",
		                   "error": f"{type(e).__name__}: {e}"})
	return json.dumps({
		"status": "ok",
		"input": str(src.resolve()),
		"output": str(out.resolve()),
		"boundary_surfaces": len(surface_tags),
		"input_size": src.stat().st_size,
		"output_size": out.stat().st_size,
	}, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_post_add_view_from_csv(msh_path: str,
                                  csv_path: str,
                                  out_path: str,
                                  tag_col: int = 0,
                                  value_cols: Any = 1,
                                  data_type: str = "NodeData",
                                  view_name: str = "csv_view",
                                  has_header: bool = True,
                                  time: float = 0.0,
                                  time_step: int = 0) -> str:
	"""Load a CSV of `tag, value[, value, ...]` and write an
	`$NodeData` / `$ElementData` block onto an existing v4.1 mesh.

	This is the most frequent post workflow: solver dumps a flat CSV
	with "node_id, B_magnitude" (or x,y,z + Bx,By,Bz), and we want it
	overlaid on the mesh for Gmsh GUI visualisation.

	Args:
	    msh_path: input v4.1 mesh (unchanged on disk).
	    csv_path: CSV file. Columns are 0-indexed.
	    out_path: output mesh + view.
	    tag_col: column holding the entity tag (node or element).
	    value_cols: single int (scalar) or list of int (vector /
	        tensor). Default: column 1 as a 1-component scalar.
	    data_type: "NodeData" (default) or "ElementData".
	    view_name: label shown in Gmsh GUI.
	    has_header: skip the first line if True.
	    time / time_step: view coordinates.
	"""
	src = Path(msh_path)
	csv = Path(csv_path)
	if not src.exists():
		return json.dumps({"status": "error",
		                   "error": f"mesh not found: {msh_path}"})
	if not csv.exists():
		return json.dumps({"status": "error",
		                   "error": f"csv not found: {csv_path}"})
	if isinstance(value_cols, int):
		value_cols = [value_cols]
	try:
		value_cols = [int(c) for c in value_cols]
	except (TypeError, ValueError):
		return json.dumps({"status": "error",
		                   "error": f"value_cols must be int or list of int"})

	# Parse CSV (stdlib only — no pandas dep for MCP).
	import csv as _csv
	values: dict[int, list[float]] = {}
	try:
		with csv.open("r", encoding="utf-8", errors="replace", newline="") as f:
			reader = _csv.reader(f)
			rows = list(reader)
	except OSError as e:
		return json.dumps({"status": "error", "error": f"csv read: {e}"})
	rows = rows[1:] if has_header else rows
	for lineno, row in enumerate(rows, start=2 if has_header else 1):
		if len(row) <= max(tag_col, max(value_cols)):
			continue
		try:
			tag = int(row[tag_col])
			v = [float(row[c]) for c in value_cols]
		except (ValueError, IndexError):
			continue
		values[tag] = v
	if not values:
		return json.dumps({"status": "error",
		                   "error": "no valid (tag, value) rows parsed"})

	# Delegate to the existing _write_model_data helper
	r = _write_model_data(str(src), str(out_path), view_name, values,
	                       data_type=data_type, time=time,
	                       time_step=time_step)
	r["csv_rows_used"] = len(values)
	r["csv_path"] = str(csv.resolve())
	return json.dumps(r, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Community retrieval: gmsh GitLab issues + StackOverflow [gmsh]
# ---------------------------------------------------------------------------

from ..common import examples as _ex


@mcp.tool()
def gmsh_examples(query: str, limit: int = 5,
                   force_refresh: bool = False) -> str:
	"""Search scraped gmsh community threads (unioned).

	Two sources:
	  1. **gmsh/gmsh GitLab issues** (`gitlab.onelab.info`, ~3000
	     threads) — bugs, feature requests, author Q&A; often the only
	     place a nuanced usage detail is explained.
	  2. **StackOverflow / SciComp.SE `[gmsh]`** — top-voted Q&A with
	     accepted answers folded in.

	Both fetched anonymously; 7-day cache. Hits carry `source` so the
	caller sees GitLab issue vs SO answer provenance.

	Args:
	    query: keywords (e.g., "Physical surface volume", "boundary
	        layer", "transfinite", "3D Delaunay refine").
	    limit: cap on returned threads.
	    force_refresh: re-scrape sources, ignoring cache TTL.
	"""
	r = _ex.search_examples("gmsh", query, limit=limit,
	                         force_refresh=force_refresh)
	return json.dumps(r, indent=2, ensure_ascii=False)


@mcp.tool()
def gmsh_examples_refresh(max_issues: int = 60,
                           max_so_questions: int = 30,
                           include_youtube: bool = True) -> str:
	"""Force-refresh all gmsh community sources.

	Args:
	    max_issues: cap for GitLab issue walk.
	    max_so_questions: cap for StackExchange questions per site
	        (stackoverflow + scicomp).
	    include_youtube: also scrape gmsh YouTube tutorials.
	"""
	out = {
		"gitlab_issues": _ex.refresh_gmsh_issues(max_issues=max_issues),
		"stackoverflow": _ex.refresh_gmsh_stackoverflow(
			max_questions=max_so_questions,
		),
	}
	if include_youtube:
		out["youtube"] = _ex.refresh_gmsh_youtube()
	return json.dumps(out, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
	if "--selftest" in sys.argv:
		print("gmsh-post MCP server self-test:")
		for t in list(_SPEC.keys())[:3]:
			body = _get_spec(t)
			print(f"  spec '{t}': {len(body)} chars")
		print("OK")
		return
	mcp.run(transport="stdio")


if __name__ == "__main__":
	main()
