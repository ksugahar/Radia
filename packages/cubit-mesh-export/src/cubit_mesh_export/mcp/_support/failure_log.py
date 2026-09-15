"""Append-only JSONL log of failed tool executions.

Purpose: let LLM-driven sessions learn from past mistakes across calls
(and across conversations). Each failure is recorded with the offending
input, the error, a short context snapshot, and timestamp. Lookup tools
can query recent entries to surface "we tried that; it failed with X"
hints without forcing the user to paste error messages again.

Location: `<state_dir>/logs/{kind}_failures.jsonl`

`state_dir` resolves to (in order):
  1. $RADIA_MCP_STATE_DIR if set
  2. Platform-appropriate dir under the user profile
     - Windows: %LOCALAPPDATA%/radia-mcp/
     - Other:   ~/.local/state/radia-mcp/
  3. Fallback: ~/.radia-mcp/
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Iterable


def state_dir() -> Path:
	env = os.environ.get("RADIA_MCP_STATE_DIR")
	if env:
		return Path(env)
	if os.name == "nt":
		base = os.environ.get("LOCALAPPDATA")
		if base:
			return Path(base) / "radia-mcp"
	xdg = os.environ.get("XDG_STATE_HOME")
	if xdg:
		return Path(xdg) / "radia-mcp"
	home = Path.home()
	linux_default = home / ".local" / "state" / "radia-mcp"
	if linux_default.exists() or os.name != "nt":
		return linux_default if os.name != "nt" else home / ".radia-mcp"
	return home / ".radia-mcp"


def _log_path(kind: str) -> Path:
	d = state_dir() / "logs"
	d.mkdir(parents=True, exist_ok=True)
	return d / f"{kind}_failures.jsonl"


def record_failure(kind: str, payload: dict[str, Any]) -> Path:
	"""Append a failure record to `<kind>_failures.jsonl`.

	Args:
	    kind: "cubit" / "build123d" / etc. Becomes the log file prefix.
	    payload: free-form dict. Callers typically include:
	        - "input":   the tool's input (commands / script)
	        - "error":   short error message
	        - "traceback": (optional) full traceback string
	        - "context": (optional) state snapshot (volumes, nodes, ...)

	Returns the log path so callers can cite it.
	"""
	path = _log_path(kind)
	rec = dict(payload)
	rec.setdefault("ts", time.time())
	rec.setdefault("iso", time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(rec["ts"])))
	rec["kind"] = kind
	try:
		with path.open("a", encoding="utf-8") as f:
			f.write(json.dumps(rec, ensure_ascii=False) + "\n")
	except OSError:
		pass
	return path


def recent_failures(kind: str, limit: int = 10) -> list[dict[str, Any]]:
	"""Read the last `limit` records from the failure log."""
	path = _log_path(kind)
	if not path.exists():
		return []
	try:
		with path.open("r", encoding="utf-8") as f:
			lines = f.readlines()
	except OSError:
		return []
	out: list[dict[str, Any]] = []
	for line in lines[-max(1, int(limit)):]:
		line = line.strip()
		if not line:
			continue
		try:
			out.append(json.loads(line))
		except json.JSONDecodeError:
			continue
	return out


def iter_failures(kind: str) -> Iterable[dict[str, Any]]:
	"""Iterate over all failure records for `kind`."""
	path = _log_path(kind)
	if not path.exists():
		return
	try:
		with path.open("r", encoding="utf-8") as f:
			for line in f:
				line = line.strip()
				if not line:
					continue
				try:
					yield json.loads(line)
				except json.JSONDecodeError:
					continue
	except OSError:
		return


def failures_as_text(kind: str, limit: int = 50) -> str:
	"""Render recent failures as a markdown-ish block for retrieval/chunking.

	Used by `*_lookup` tools: the failure log is fed in as an additional
	knowledge source, so past mistakes are searchable with the same query
	interface as the static knowledge base.
	"""
	recs = recent_failures(kind, limit=limit)
	if not recs:
		return ""
	lines: list[str] = [f"# {kind} failure log (last {len(recs)})"]
	for r in recs:
		lines.append("")
		lines.append(f"## {r.get('iso', '?')} — {str(r.get('error',''))[:80]}")
		inp = r.get("input")
		if inp:
			s = inp if isinstance(inp, str) else json.dumps(inp, ensure_ascii=False)
			lines.append("**input:**")
			lines.append("```")
			lines.append(s[:1200])
			lines.append("```")
		err = r.get("error")
		if err:
			lines.append(f"**error:** {err}")
		ctx = r.get("context")
		if ctx:
			lines.append(f"**context:** {json.dumps(ctx, ensure_ascii=False)[:400]}")
	return "\n".join(lines)
