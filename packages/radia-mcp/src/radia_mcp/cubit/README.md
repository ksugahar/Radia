# `radia_mcp.cubit` — Cubit mesh scripting MCP server

**83 MCP tools** — the largest subpackage in the radia-mcp wheel.
Production-grade workflow for Coreform Cubit hex/tet meshing via
Python (`cubit.cmd`), STEP import, NGSolve `.vol` export with the
canonical `check-vol` gate, a persistent Cubit session (daemon), and
the Sugahara lab's curated Cubit knowledge corpus.

The authoritative, auto-generated tool list lives in
[`docs/TOOLS.md`](../../../docs/TOOLS.md#mcp-server-cubit) — regenerate
with `scripts/gen_tools_doc.py` after adding/renaming tools.

## Quick start

```bash
pip install radia-mcp
mcp-server-cubit                    # stdio server
mcp-server-cubit --selftest         # lightweight self-test
# in Claude Code: register `cubit` in .mcp.json
```

Then in a session:

```
> cubit_status()                            # live capability list
> cubit_show(path="coil.step")              # open in persistent Cubit session
> cubit_probe(query="entities")             # Probe-Don't-Guess dump
> cubit_mesh_auto(step_path="coil.step")    # scheme-ladder auto mesh
> cubit_check_vol(vol_path="coil.vol")      # canonical check-vol gate
```

## Tool families (83 total)

| Family | Examples |
|---|---|
| **Live session** | `cubit_show`, `cubit_exec`, `cubit_exec_safely`, `cubit_probe`, `cubit_snapshot`, `cubit_session_status`, `cubit_session_shutdown` |
| **Checkpoint / restore** | `cubit_checkpoint`, `cubit_restore`, `cubit_list_checkpoints` |
| **Headless batch** | `cubit_batch_try`, `cubit_mesh_auto`, `open_in_cubit` |
| **Mesh race (variant exploration)** | `cubit_mesh_race`, `cubit_mesh_race_smart[_async]`, `cubit_mesh_race_review[_async]`, `cubit_mesh_race_status`, `cubit_mesh_apply_choice`, `cubit_mesh_race_with_human`, `cubit_curate_learned_recipes` |
| **Export / .vol gates** | `cubit_check_vol` (canonical check-vol), `cubit_vol_inventory`, `cubit_gmsh_v41_inventory`, `cubit_headless_netgen_export_gate`, `cubit_mixed_order_series_gate`, ~30 further scenario gates (`*_gate`) |
| **Diagnostics** | `cubit_mesh_diagnose`, `cubit_suggest_next`, `cubit_recent_failures`, `cubit_diagnostics_guide` |
| **Lint** | `lint_cubit_script`, `lint_cubit_directory`, `cubit_audit_summary`, `get_lint_rules`, `generate_cubit_script` |
| **Knowledge** | `cubit_ask`, `cubit_lookup`, `cubit_docs`, `cubit_forum_tips`, `cubit_web_docs`, `cubit_examples[_refresh]`, `netgen_workflow_guide`, `netgen_code_example` |
| **Toolbar / SDK** | `cubit_toolbar_guide`, `cubit_scaffold_toolbar`, `cubit_generate_dialog`, `cubit_cpp_sdk_guide` |
| **Status / meta** | `cubit_status` |

Run `cubit_status()` for the live, definitive list.

## Persistent session architecture

```
Claude Code (MCP client)
    → mcp-server-cubit (system Python 3.12)
        → daemon.py (Cubit's bundled Python 3.10, subprocess, JSON-RPC)
            → cubit.cmd(...) → Cubit GUI or batch session
```

`cubit_show` / `cubit_exec` reuse one persistent Cubit process
(license-friendly, <1 ms per command after init); `cubit_batch_try` /
`cubit_mesh_auto` spawn fresh headless subprocesses for dry-runs.
`cubit_check_vol` needs **no Cubit at all** — it runs the
`cubit_mesh_export.check` engine (NGSolve) in the server's Python 3.12.

## Lab-specific workflows

This server encodes the Sugahara lab's hard-won Cubit usage rules
(see `CLAUDE.md` "Cubit Block/Sideset Label Convention",
"Journal File Portability Policy", "AI-Driven Cubit: Probe, Don't
Guess", "Mesh Export Consistency Check Policy"):

1. **Never hardcode entity IDs** — `lint_cubit_script` flags them.
   Use `get_last_id()` + geometric predicates.
2. **Probe before classifying** — `cubit_probe(query="entities")`
   returns every volume/surface with centroid + bbox extent +
   measure in one call; derive classification predicates from the
   printed values, never from a-priori reasoning about the `.jou`.
3. **Separate blocks for material vs boundary** —
   `cubit_probe(query="labels")` audits the live session: mixed
   volume+surface blocks (label LOST on export), unnamed blocks,
   casefold collisions, non-snake-case names.
4. **Every solver-bound `.vol` passes check-vol** —
   `cubit_check_vol(vol_path, contract=..., strict_labels=True,
   report_json=...)` before any solver / Simulink initialization.

## Knowledge corpus shipped

- **Coreform Cubit official help** (curated index + `cubit_web_docs`)
- **Cubit Discourse forum tips** (`cubit_forum_tips`)
- **Coreform webinars** (`knowledge/coreform_webinars.py`)
- **Cubit Python API reference** (600+ functions, `api_reference.py`)
- **In-tree examples** at `src/radia/panels/samples/*.jou`

## When to use this vs Cubit GUI

| Task | This server | Cubit GUI |
|---|---|---|
| Reproducible mesh from STEP | ✅ | (manual) |
| Mesh-quality probing | ✅ | (manual) |
| Cross-version `.jou` | ✅ (lint) | (no help) |
| `.vol` gate before solver | ✅ (`cubit_check_vol`) | ❌ |
| Single-shot mesh you'll never re-do | ❌ overkill | ✅ |
| Visual debugging | `cubit_snapshot` | ✅ |
| License-saving batch runs | ✅ (`-batch -nographics`) | ❌ |

## Cross-references

- `mcp-server-gmsh` — mesh post/visualization (`.msh v4.1`)
- `mcp-server-build123d` — Pythonic OCCT front-end → STEP → Cubit
- `mcp-server-fem` — FEM formulations consuming the resulting `.vol`
- `cubit-mesh-export` (PyPI) — `check-vol` CLI + Cubit plugin binaries

## Source

- `src/radia_mcp/cubit/server.py` — tool registration
- `src/radia_mcp/cubit/session.py` + `daemon.py` — persistent session
  (Python 3.12 client ↔ Cubit Python 3.10 JSON-RPC daemon)
- `src/radia_mcp/cubit/label_audit.py` — cubit-free block/sideset
  convention audit (shared by daemon `probe("labels")` and tests)
- `src/radia_mcp/cubit/knowledge/` — sub-modules (export,
  netgen_workflow, scripting, cpp_sdk, custom_toolbar, ...)
- `src/radia_mcp/cubit/vol_inventory.py`, `gmsh_v41.py`, `*_gate.py`
  — export gates
