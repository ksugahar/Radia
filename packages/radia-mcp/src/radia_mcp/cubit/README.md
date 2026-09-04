# `radia_mcp.cubit` — Cubit mesh scripting MCP server

The production `core` profile exposes the operating workflows directly and
groups fine-grained scenario/identity checks behind `cubit_validation_catalog`
and `cubit_validation_run`. The underlying Python checks remain independently
testable.
Production-grade workflow for Coreform Cubit hex/tet meshing via
Python (`cubit.cmd`), STEP import, NGSolve `.vol` export with the
canonical `check-vol` gate, a persistent Cubit session (daemon), and
the Sugahara lab's curated Cubit knowledge corpus.

The live `tools/list` response is authoritative. Use `cubit_status`,
`cubit_topics`, and the `mcp-server-radia-meta` catalog for discovery. For an
offline local snapshot only, run `scripts/gen_tools_doc.py`; its output is not
version-controlled.

## Quick start

```bash
pip install radia-mcp
mcp-server-cubit                    # stdio server
mcp-server-cubit --selftest         # lightweight self-test
```

Client setup (stdio transport; the server needs no network port):

```bash
# Claude Code
claude mcp add cubit -- mcp-server-cubit
# remove again with: claude mcp remove cubit
```

```json
// Claude Desktop / VS Code (mcp.json style)
{ "mcpServers": { "cubit": { "command": "mcp-server-cubit" } } }
```

Environment knobs:

| Variable | Effect |
|---|---|
| `CUBIT_BIN_DIR` / `CUBIT_INSTALL_DIR` | Override Coreform Cubit install discovery |
| `RADIA_CUBIT_SESSION_MODE` | `auto` (default): attach to the live shared daemon, else spawn. `new`: always spawn a fresh daemon in a private per-process drop dir (hermetic CI; removed on shutdown). `existing`: attach only — fail loud when no shared daemon is running |
| `RADIA_MCP_TOOL_PROFILE=full` | Restore legacy individual validation tools for migration/debugging; production defaults to `core` |
| `RADIA_MCP_CUBIT_GATES=0` | In the `full` compatibility profile, additionally hide direct `*_gate` tools |
| `RADIA_CUBIT_EAGER=1` | Start the Cubit session in the background at server startup (hides the 30+ s first-call cost) |
| `RADIA_MCP_CUBIT_CALL_LOG=0` | Disable the all-calls JSONL log (`<state_dir>/logs/cubit_tool_calls.jsonl`) |

One-shot environment preparation (license warmup + full doctor report,
exit 1 when problems are found):

```bash
mcp-server-cubit --setup
```

First stop when anything misbehaves: `cubit_doctor()` — a read-only
one-shot diagnosis of install discovery, license cache, deployed-plugin
freshness (hash vs the cubit-mesh-export copy), daemon state, drop-dir
startup diagnostics, and check-vol dependencies.

Then in a session:

```
> cubit_status()                            # live capability list
> cubit_show(path="coil.step")              # open in persistent Cubit session
> cubit_probe(query="entities")             # Probe-Don't-Guess dump
> cubit_mesh_auto(step_path="coil.step")    # scheme-ladder auto mesh
> cubit_check_vol(vol_path="coil.vol")      # canonical check-vol gate
```

## Tool families (85 total)

| Family | Examples |
|---|---|
| **Live session** | `cubit_show`, `cubit_exec`, `cubit_exec_safely`, `cubit_probe`, `cubit_snapshot`, `cubit_session_status`, `cubit_session_shutdown`, `cubit_session_journal` (export the session as a replayable `.jou`) |
| **Environment** | `cubit_doctor` (one-shot install/license/plugin/daemon diagnosis) |
| **Checkpoint / restore** | `cubit_checkpoint`, `cubit_restore`, `cubit_list_checkpoints` |
| **Headless batch** | `cubit_batch_try`, `cubit_mesh_auto`, `open_in_cubit` |
| **Mesh race (variant exploration)** | `cubit_mesh_race`, `cubit_mesh_race_smart[_async]`, `cubit_mesh_race_review[_async]`, `cubit_mesh_race_status`, `cubit_mesh_apply_choice`, `cubit_mesh_race_with_human`, `cubit_curate_learned_recipes` |
| **Export / .vol gates** | `cubit_check_vol` (canonical check-vol), `cubit_vol_inventory`, `cubit_gmsh_v41_inventory`, `cubit_headless_netgen_export_gate`, `cubit_mixed_order_series_gate`, ~30 further scenario gates (`*_gate`) |
| **Cross-mesher quality** | `cubit_netgen_quality_compare` — one STEP through Netgen tet + Cubit tet + Cubit hex, judged by ONE gmsh minSICN referee (same metric implementation for every route; tet-vs-tet is the directly comparable pair, hex reported as the structured reference) |
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
        → GUI: bootstrap.py (file-drop JSON-RPC inside Cubit's Qt/Python)
          batch: daemon.py (Cubit's bundled Python 3.10, stdio JSON-RPC)
            → probe_ops.py (SHARED probe queries — no transport drift)
                → cubit.cmd(...) / cubit API
```

`cubit_show` / `cubit_exec` reuse one persistent Cubit process
(license-friendly, <1 ms per command after init); `cubit_batch_try` /
`cubit_mesh_auto` spawn fresh headless subprocesses for dry-runs.
`cubit_check_vol` needs **no Cubit at all** — it runs the
`cubit_mesh_export.check` engine (NGSolve) in the server's Python 3.12.

Session robustness (MathWorks MATLAB-MCP patterns, 2026-08-05):

- **Startup failures report the real error**: the GUI bootstrap writes
  `startup_error.txt` on any in-process exception, Cubit's own console
  goes to `cubit_stdout.log` / `cubit_stderr.log` in the per-user drop
  dir, and the ready-poll surfaces those instead of a bare timeout.
- **Ownership-tagged cleanup**: recovery paths only ever kill a Cubit
  THIS process spawned; a live daemon another window started is detached
  from, never terminated. The explicit `cubit_session_shutdown` tool is
  the one way to stop a foreign/hung daemon, and it reports which
  process it stopped.
- **Errors carry `kind`**: `"input"` (fix your commands and retry),
  `"environment"` (license/install/hung — tell the user), `"internal"`
  (server bug — do not retry), plus a `log` pointer to the drop-dir
  diagnostics. Server-level MCP `instructions` teach connecting models
  the same contract.
- **Every tool is annotation-classified** (read-only / read-only+web /
  file-writing / session-destructive presets) so MCP clients can gate
  permissions correctly.
- **`cubit_snapshot` returns the PNG inline** as MCP image content —
  the model sees the current view directly.

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

## Driving policy: headless is primary, GUI is the user's debugger

Lab policy (2026-08-05): agents drive Cubit through **APREPRO + Python on
the headless/batch route** — `.jou` playback, `cubit_batch_try`,
`cubit_mesh_auto`, the batch stdio daemon. That is the primary path for
mesh generation, exports, gates, and validation, and it must never
require a GUI window. The **persistent GUI session** (`cubit_show`,
`cubit_snapshot`) is maintained as the **user's visual-debugging aid** —
open it when a human wants to watch the model or capture a figure
(`cubit_snapshot` needs the rendering window; batch reports an honest
`ok=false` there).

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

## Licensing boundary and data collection

- This server drives **your own licensed Coreform Cubit install**: it
  ships no Cubit binaries and no license/credentials, and must not be
  used to share one Cubit license between users (Coreform RLM licensing
  is per-user; each user activates their own seat).
- The persistent session and headless batch runs each consume a license
  seat while alive.
- **This server sends no usage data anywhere.** The only records are
  local diagnostics under the per-user drop dir and the local failure
  log.

## Cross-server API compatibility (cubit ↔ build123d ↔ external CAD)

`mcp-server-cubit` and `mcp-server-build123d` share one hardening layer
(`radia_mcp.common.server_hardening`: annotation presets, error-kind
contract, gate hiding, all-calls JSONL log) and one **probe contract**:

| Concept | build123d (CAD side) | Cubit (mesh side) | History-based CAD (e.g. CST) |
|---|---|---|---|
| Replayable history | the Python script | `.jou` journal (`cubit_session_journal`) | history list |
| Named bodies | `part.label` → STEP names | entity names → blocks | component/solid names |
| Per-body probe | `build123d_probe(path, "entities")` | `cubit_probe("entities")` | volume-evidence rows |
| Naming audit | `build123d_probe(path, "labels")` | `cubit_probe("labels")` | — |

Both `entities` probes emit the SAME core keys per body —
`{id, centroid, bbox_min, bbox_max, extent, volume}` (faces:
`{id, center, bbox_min, bbox_max, extent, area}`) — locked by
`PROBE_SOLID_CORE_KEYS` / `PROBE_FACE_CORE_KEYS` in
`radia_mcp.common.server_hardening` and the
`test_b3d_cubit_probe_compat.py` contract test, so an agent can author a
labeled STEP, mesh it, and compare per-body volumes/centroids directly
(verified end-to-end: identical to 2e-16 relative on a 2-solid
assembly). External-CAD evidence rows (Cubit/CST) connect through the
build123d volume-crosscheck tools with mandatory units.

## Cross-references

- `mcp-server-gmsh` — mesh post/visualization (`.msh v4.1`)
- `mcp-server-build123d` — Pythonic OCCT front-end → STEP → Cubit
- `mcp-server-fem` — FEM formulations consuming the resulting `.vol`
- `cubit-mesh-export` (PyPI) — `check-vol` CLI + Cubit plugin binaries

## Source

- `src/radia_mcp/cubit/server.py` — tool registration, instructions,
  annotation presets
- `src/radia_mcp/cubit/session.py` + `bootstrap.py` (GUI file-drop) +
  `daemon.py` (batch stdio) — persistent session
- `src/radia_mcp/cubit/probe_ops.py` — SHARED probe queries for both
  runners (summary/quality/per_volume/entities/labels)
- `src/radia_mcp/cubit/label_audit.py` — cubit-free block/sideset
  convention audit (shared by `probe("labels")` and tests)
- `src/radia_mcp/cubit/knowledge/` — sub-modules (export,
  netgen_workflow, scripting, cpp_sdk, custom_toolbar, ...)
- `src/radia_mcp/cubit/vol_inventory.py`, `gmsh_v41.py`, `*_gate.py`
  — export gates
