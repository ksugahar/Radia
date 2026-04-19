# radia-mcp

[![PyPI](https://img.shields.io/pypi/v/radia-mcp.svg)](https://pypi.org/project/radia-mcp/)
[![Python](https://img.shields.io/pypi/pyversions/radia-mcp.svg)](https://pypi.org/project/radia-mcp/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

> **First-and-only public Model Context Protocol (MCP) server suite for
> Coreform Cubit, Gmsh, build123d, and the Radia CAE ecosystem.**
> Pioneering MCP territory for mesh generators worldwide.

Authored by the **Sugawara Lab (菅原研究室)**, Kindai University —
where the lab-standard primary pair is **build123d (CAD authoring) +
Cubit (hex meshing)**, with Gmsh as the post-processing workhorse.

---

## Why this exists

`radia-mcp` lets an AI assistant drive a real CAE pipeline end-to-end:

```
build123d script  →  STEP  →  Cubit hex mesh  →  .msh v4.1 / .vol  →  ngsolve / radia
        ↑                            ↓
   AI authors                   AI iterates safely
   (Builder API)                 (auto-checkpoint
                                  + batch dry-run
                                  + live GUI replay)
```

What makes it different from typical CAD-MCP servers (FreeCAD, OpenSCAD,
Blender):

- **Test-then-reflect safety pattern**: every risky operation
  auto-checkpoints to `.cub5`, runs in a disposable batch Cubit first,
  reflects to the live GUI only when batch passes cleanly. The user
  watches success paths, not failures.
- **Scheme ladder with geometry split**: `cubit_mesh_auto` walks
  `auto → sweep → polyhedron → tetmesh`, auto-detects compound bodies
  (volume/surface ratio heuristic) and auto-`webcut`s before retrying.
- **Universal mesh backend** for the CAD-MCP ecosystem:
  `any_step_to_cubit_hex` accepts STEP from any upstream MCP
  (FreeCAD-MCP / OpenSCAD-MCP / cadquery-mcp / Blender-MCP / Onshape-MCP
  / KiCad-MCP) and dispatches to the same hex-mesh backend.
- **Deeply scraped knowledge** (≈800 lab files + 65 build123d examples
  + 15 bd_warehouse + GitHub issues + GitLab issues + StackOverflow +
  YouTube tutorial transcripts + Coreform training pack), all
  searchable via tf-idf with heading boost.

---

## Demo (placeholder — recordings to be added)

`build123d_to_cubit_hex(script=generate_build123d_script("helix_coil")["script"], target_size=1.0)`
in one call:

```
1. build123d Builder API: 4-turn helix + 2 radial leads (3 prismatic bodies)
2. STEP export → /tmp/coil.step
3. cubit_batch_try (headless): scheme ladder
   ├─ auto → 1668 hex / 0 tet ✓ WIN
   └─ (sweep / polyhedron / tetmesh skipped)
4. live Cubit GUI replay → user sees the winning recipe execute
5. .cub5 checkpoint preserved for rollback
```

Total time: ~30 s on a 4-core laptop. AI watches state delta, applies
hint heuristics on failure, and never leaves a half-broken mesh in the
user's window.

---

## Install

```bash
pip install radia-mcp                 # core (Cubit + build123d + gmsh-post + interop)
pip install radia-mcp[build123d]      # adds build123d itself
pip install radia-mcp[cadquery]       # adds CadQuery (interop with cadquery-mcp)
pip install radia-mcp[gmsh]           # adds gmsh Python bindings (gmsh_post operations)
pip install radia-mcp[youtube]        # adds youtube-transcript-api (tutorial scrape)
pip install radia-mcp[radia]          # adds Radia core (radia-coupled servers)
pip install radia-mcp[full]           # everything above
```

Requires Python ≥ 3.10. Coreform Cubit is auto-discovered from
standard install paths; set `CUBIT_BIN_DIR` env var if installed
non-standard.

---

## The four MCP servers

| Server | Entry point | Tools | Highlights |
|---|---|---|---|
| **Cubit** | `mcp-server-cubit` | 52 | `cubit_mesh_auto`, `cubit_exec_safely`, `cubit_ask`, scheme ladder + geometry split, .cub5 checkpoint/restore, scrape index over Coreform forum + S:\\CoreformCubit lab archive (787 files) + YouTube + Coreform training |
| **build123d** | `mcp-server-build123d` | 37 | `build123d_to_cubit_hex`, `lint_build123d_script` (7 rules), `build123d_try` (subprocess isolation), `build123d_inspect_step`, `build123d_heal`, `build123d_api`, 13 Radia/general templates (`magnet_ring`, `c_core`, `racetrack_coil`, …), CadQuery + bd_warehouse interop |
| **gmsh-post** | `mcp-server-gmsh-post` | 21 | `gmsh_post_inspect/validate/convert` (MSH v4.1 only), `gmsh_post_quality`, `gmsh_post_extract_physical`, `gmsh_post_boundary`, `gmsh_post_add_view_from_csv`, scrape over GitLab issues + StackOverflow + YouTube |
| **radia-interop** | `mcp-server-radia-interop` | 4 | `any_step_to_cubit_hex` (any source) + `freecad_to_cubit_hex` (FreeCADCmd subprocess) + `openscad_to_cubit_hex` (OpenSCAD CLI) + `list_cad_mcp_interop` |

`pip install radia-mcp[radia]` additionally enables the radia-coupled
servers (radia-ngsolve / ih / electromagnet / peec / gmsh authoring /
elf), which depend on the Radia core.

---

## Claude Code / Desktop configuration

Add to your MCP client config (Claude Desktop, Claude Code, Cursor,
Continue, …):

```json
{
  "mcpServers": {
    "cubit":         {"command": "mcp-server-cubit"},
    "build123d":     {"command": "mcp-server-build123d"},
    "gmsh-post":     {"command": "mcp-server-gmsh-post"},
    "radia-interop": {"command": "mcp-server-radia-interop"}
  }
}
```

For local development from a checkout (no install needed):

```json
{
  "mcpServers": {
    "cubit": {
      "command": "python",
      "args": ["-m", "radia_mcp.cubit.server"]
    }
  }
}
```

---

## Knowledge bases shipped (offline)

| Server | Bundled knowledge | Lines |
|---|---|---|
| Cubit | `cubit_api_reference.py` (600+ functions) + scripting + forum tips + netgen workflow + export rules + panels | ≈ 29 000 |
| build123d | Auto-generated API reference (`inspect.getmembers`, 142 classes / 65 functions) + 18 curated topics (Plane/Axis/Location cookbook, Builder ↔ Algebra rosetta, joints, assemblies, CAE workflow, …) | 1 673 |
| gmsh-post | Auto-generated API reference (651 entries across `gmsh.model`, `.view`, `.option`, `.fltk`, …) + MSH v4.1 spec + view_data + physical_groups cookbooks | 2 008 |

Plus persistent **failure log** per kind, fed into every `*_lookup` /
`*_ask` retrieval so past mistakes are searchable next session.

## Live-scraped knowledge (cached 7 days)

| Source | Kind | Volume |
|---|---|---|
| Coreform forum (Discourse) | `cubit_examples` | full archive walk, ≈ 60 code-bearing posts |
| `S:\CoreformCubit` + `Radia/01_GitHub/examples` | `cubit_examples` | 787 `.jou` / `.py` files (lab archive + Coreform training pack) |
| Cubit YouTube tutorials | `cubit_examples` | 5 transcripts |
| GitHub `.jou` code search (PAT) | `cubit_examples` | 15 files |
| `gumyr/build123d/examples` | `build123d_examples` | 65 scripts |
| `gumyr/bd_warehouse` | `build123d_examples` | 15 modules |
| `gumyr/build123d` Issues + GraphQL Discussions | `build123d_examples` | 60 + 50 |
| build123d YouTube | `build123d_examples` | 4 transcripts |
| `gitlab.onelab.info/gmsh/gmsh` issues | `gmsh_examples` | walking 3000+ |
| StackOverflow / SciComp.SE `[gmsh]` | `gmsh_examples` | 30 votes-sorted |
| gmsh YouTube | `gmsh_examples` | 5 transcripts |

`{cubit,build123d,gmsh}_examples_refresh()` re-scrapes anytime;
otherwise cache TTL is 7 days.

---

## Lab stance (Sugawara Lab, 菅原研究室)

| Tool | Position |
|---|---|
| **build123d** | **主力 (push)** — new lab work is authored here |
| **Cubit** | **主力 (push)** — hex mesh + visualization |
| **FreeCAD** | **応援 (friendly)** — first-class interop, not first-class authoring; lab respects the FreeCAD community |
| CadQuery | interop / compat — OCCT sibling |
| OpenSCAD | interop / compat — legacy scripts |

The MCP-tool layer encodes this stance: `list_cad_mcp_interop` reports
`primary_pair: "build123d (CAD authoring) + Cubit (hex mesh)"` and
flags FreeCAD as `friendly`, others as `compat`.

---

## Design pattern: Checkpoint-Batch-Commit-Race (CBCR)

The signature workflow `cubit_mesh_race_with_human` realizes a
publishable design pattern we call **CBCR**: the user and *N* AI
agents race the same starting state; first-to-finish wins; the
human's in-progress work is never silently overwritten.

Full design document (defensive prior-art publication, BSD-3-Clause):
[`docs/design/checkpoint_batch_commit_race.md`](docs/design/checkpoint_batch_commit_race.md).

## Architecture (Plan A)

The Cubit server uses **Plan A**: launch
`coreform_cubit.exe -nojournal cubit_bootstrap.py`, which installs a
PyQt5 `QTimer` (200 ms poll) inside Cubit's Qt event loop. The MCP
side drops `*.req.json` into a temp directory; the bootstrap polls,
runs `cubit.cmd()` on the Qt main thread, and writes `out/*.resp.json`
back. **No sockets, no pipes** — just atomically-renamed JSON files.

This sidesteps every ABI / event-loop / GIL issue the COMSOL-style
"external Python + library load" approach hits.

`cubit_exec_safely` adds a second safety layer on top:

1. Auto-`save as` the live GUI state to `~/.cubit_viewer/checkpoints/autosafe_<ts>.cub5`.
2. Spawn a fresh batch Cubit, `open` the checkpoint, run candidate commands.
3. If batch passes (and `cubit.get_error_count()` didn't tick up — silent-error guard), replay on the live GUI.
4. Otherwise: live GUI untouched, checkpoint label returned for `cubit_restore`.

---

## Quickstart

### build123d → Cubit hex mesh, one call

```python
build123d_to_cubit_hex(
    script=generate_build123d_script("helix_coil")["script"],
    target_size=1.0,
    prefer="hex",
    commit_to_gui=True,
)
# ⇒ STEP exported → batch ladder picks scheme auto → live Cubit GUI
#   shows 1668 hex / 0 tet / 3780 nodes
```

### Safe Cubit live edit

```python
cubit_exec_safely(commands=[
    "volume all size 1.0",
    "volume all scheme tetmesh",
    "mesh volume all",
])
# ⇒ auto-checkpoint → batch dry-run → only on pass: live GUI mesh
#   if dry-run fails, GUI is untouched + rollback label returned
```

### Search any layer of knowledge

```python
cubit_ask("hex meshing tutorial")
# ⇒ unioned ranked hits across bundled KB + lab archive (787 files)
#   + Coreform forum + YouTube transcripts
```

### Convert + post-process a Cubit-exported mesh

```python
gmsh_post_convert("/tmp/cubit_export.msh", "/tmp/v41.msh")  # any → v4.1
gmsh_post_validate("/tmp/v41.msh")                          # spec compliance
gmsh_post_quality("/tmp/v41.msh")                           # min Jacobian etc.
gmsh_post_add_view_from_csv(
    msh_path="/tmp/v41.msh",
    csv_path="/tmp/B_field.csv",
    out_path="/tmp/v41_with_B.msh",
    view_name="B_magnitude",
)
```

---

## Contributing

Bug reports + PRs welcome — particularly for:

- **New CAD-MCP adapters** in `radia_mcp.interop.server` (Blender,
  Onshape, KiCad, AutoCAD via `any_step_to_cubit_hex` plumbing).
- **Additional scrape sub-sources** under `radia_mcp.common.examples`
  (mailing list archives, more YouTube channels, blog posts).
- **Cookbook topics** for `build123d_usage` / `gmsh_post_spec` —
  worked-example knowledge fragments are always welcome.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the lightweight workflow.
History of every release is in [CHANGELOG.md](CHANGELOG.md).

## License

BSD-3-Clause — see [LICENSE](LICENSE).

## Acknowledgments

- **Coreform** for Cubit and the open Discourse forum that we scrape
  with appreciation.
- **Roger Maitland (gumyr)** for build123d and bd_warehouse — the
  lab-standard upstream CAD.
- **Christophe Geuzaine** for Gmsh and the open issue tracker on
  `gitlab.onelab.info`.
- **MCP / Anthropic** for the protocol that made AI-driven CAE
  practical.
