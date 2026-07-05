"""radia_mcp.mathematica — Wolfram Mathematica subprocess bridge for MCP.

Exposes `mathematica_evaluate` and 8 high-level wrappers (simplify, to_tex,
check_identity, vector_calc, unit_convert, solve, integrate, differentiate,
status) as MCP tools.  The intent is "AI as Mathematica copilot": Claude
generates Wolfram Language code from natural-language intent, runs it via
this tool, interprets the result, and explains it back to the user.

Requires `wolframscript.exe` on PATH (Mathematica install) or the
WOLFRAMSCRIPT_PATH environment variable.  Without it, all tools return
{"ok": False, "error": "..."} cleanly — the server itself loads.

Promoted on 2026-05-20 from
  legacy private source tree
into radia-mcp as part of consolidating the LAB MCP knowledge under a single
public PyPI package.  Companion knowledge for symbolic differential-forms
manipulation: pairs well with `radia_mcp.differential_forms` for verifying
Maxwell equations, vector-calculus identities, and unit conversions.

Domain notes (read before answering FEM / NGSolve questions):
- `notes_fem_hcurl.md` — hierarchical H(curl) shape functions on tetrahedra,
  the meaning of NGSolve `nograds=True`, tree-cotree gauge fixing and its
  generalizations, and Mathematica verification recipes.
"""
