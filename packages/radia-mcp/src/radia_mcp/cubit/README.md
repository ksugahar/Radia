# `radia_mcp.cubit` — Cubit mesh scripting MCP server

**43 MCP tools** — the largest subpackage in the radia-mcp wheel.
Production-grade workflow for Coreform Cubit hex/tet meshing via
Python (`cubit.cmd`), STEP import, NGSolve `.vol` export, and the
Sugahara lab's curated Cubit-2025.12 scrape index (787-file local
archive + Coreform Discourse + YouTube transcripts).

## Quick start

```bash
pip install radia-mcp
mcp-server-cubit                    # stdio server
# in Claude Code: register `cubit` in .mcp.json
```

Then in a session:

```
> cubit_status()                    # what tools are available?
> cubit_topics()                    # list the 43 topic areas
> cubit_mesh_auto(step="coil.step", scheme="hex_sweep")
```

## Tool families (43 total)

| Family | Count | Examples |
|---|---:|---|
| **Execution** | 5 | `cubit_exec`, `cubit_exec_safely`, `cubit_batch`, `cubit_mesh_auto`, `cubit_run_journal` |
| **Mesh / scheme** | 9 | `cubit_scheme_ladder`, `cubit_split_geometry_at_imprint`, `cubit_imprint_merge`, `cubit_hex_sweep`, `cubit_tet_fallback`, `cubit_get_mesh_quality`, `cubit_get_volume_count`, `cubit_get_element_count`, `cubit_size_function` |
| **Checkpoint / restore** | 4 | `cubit_save_cub5`, `cubit_restore_cub5`, `cubit_save_journal`, `cubit_diff_cub5` |
| **Knowledge** | 13 | `cubit_ask`, `cubit_search_scrape`, `cubit_docs`, `cubit_youtube_search`, `cubit_youtube_transcript`, `cubit_discourse_search`, `cubit_local_archive_grep`, `cubit_training_recipe`, `cubit_lab_jou_snippet`, `cubit_block_sideset_convention`, `cubit_journal_portability_lint`, `cubit_high_order_curving`, `cubit_export_format_reference` |
| **Export** | 7 | `cubit_export_netgen`, `cubit_export_gmsh`, `cubit_export_nastran`, `cubit_export_vtk`, `cubit_export_femeem`, `cubit_export_curved_high_order`, `cubit_export_companion_json` |
| **Geometry inspect** | 4 | `cubit_get_volume_info`, `cubit_get_surface_info`, `cubit_get_curve_info`, `cubit_probe_entity` |
| **Status / meta** | 2 | `cubit_status`, (built-in MCP `tools/list`) |
| **Bibliography** | 1 | `cubit_bibliography_index` |

Run `cubit_status()` for the live, definitive list.

## Lab-specific workflows

This server encodes the Sugahara lab's hard-won Cubit usage rules
(see `CLAUDE.md` "Cubit Block/Sideset Label Convention",
"Journal File Portability Policy", "AI-Driven Cubit: Probe, Don't
Guess"):

1. **Never hardcode entity IDs** — `cubit_journal_portability_lint`
   flags them. Use `get_last_id()` + geometric predicates.
2. **Separate blocks for material vs boundary** — `cubit_block_
   sideset_convention` documents the priority ladder.
3. **Probe before classifying** — `cubit_probe_entity` returns
   centroid + bbox + area / volume; classify in code, not by
   visual inspection.

## Knowledge corpus shipped

- **Coreform Cubit official help** (offline copy + scrape index)
- **Cubit Discourse forum** (search + cached threads)
- **Coreform training videos** (YouTube transcripts, search indexed)
- **Sugahara lab Cubit archive** at `public-safe curated corpus` (787 files,
  proven .jou snippets, lab-specific recipes)
- **In-tree examples** at `src/radia/panels/samples/*.jou`

## When to use this vs Cubit GUI

| Task | This server | Cubit GUI |
|---|---|---|
| Reproducible mesh from STEP | ✅ | (manual) |
| Mesh-quality probing | ✅ | (manual) |
| Cross-version `.jou` | ✅ (lint) | (no help) |
| Single-shot mesh you'll never re-do | ❌ overkill | ✅ |
| Visual debugging | ❌ | ✅ |
| Cubit license-saving (no GUI per run) | ✅ (`-batch -nographics`) | ❌ |

## Cross-references

- `mcp-server-gmsh` — alternative mesh post (`.msh v4.1`)
- `mcp-server-build123d` — Pythonic OCCT front-end → STEP → Cubit
- `mcp-server-fem` — what FEM formulations consume the resulting
  `.vol`

## Source

- `src/radia_mcp/cubit/server.py` — tool registration
- `src/radia_mcp/cubit/cubit_knowledge.py` — main knowledge text
- `src/radia_mcp/cubit/knowledge/` — sub-modules (cpp_sdk,
  custom_toolbar, scrape_index, etc.)
- `src/radia_mcp/cubit/bibliography_index_knowledge.py` —
  auto-generated bibliography
