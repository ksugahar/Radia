# radia-mcp

[![PyPI](https://img.shields.io/pypi/v/radia-mcp.svg)](https://pypi.org/project/radia-mcp/)
[![Python](https://img.shields.io/pypi/pyversions/radia-mcp.svg)](https://pypi.org/project/radia-mcp/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

> **First-and-only public Model Context Protocol (MCP) server suite for
> Coreform Cubit, Gmsh, build123d, and the Radia CAE ecosystem —
> including differential geometry and Mathematica integration.**
> Pioneering MCP territory for mesh generators worldwide.

**Killer demo (30 seconds)**: ask Claude to derive the Kelvin transform
factor by hand — it gets stuck on a 3×3 Jacobian + 27-term Laplacian.
Then ask it to use `differential-forms` + `mathematica` together:

```
> verify_with_mathematica(identity="kelvin")
```

Claude pulls the recipe, sends it to Wolfram, reports back

```
k=1:  Laplacian(psi) = 0          [harmonic — factor R/|y| is correct]
k=2:  Laplacian(psi) ≠ 0          [factor wrong]
k=3:  Laplacian(psi) ≠ 0          [factor wrong]
```

— in 8 seconds.  The `kelvin_factor = R/|y|` that takes half a day
in vector calculus appears as the conformal weight λ^((n−2k)/2) of a
k-form, and Mathematica verifies it symbolically.

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
- **Symbolic verification of FEM formulations** (no other MCP server
  does this).  Pair `differential-forms` (theory: Bossavit Whitney
  complex, Arnold-Falk-Winther FEEC, Kameari edge elements) with
  `mathematica` (Wolfram subprocess bridge) to:
  - Verify d² = 0 on any specific function
  - Compute element mass / stiffness matrices on a unit tetrahedron
    in closed form
  - Derive `kelvin_factor = R/|y|` for the Kelvin transformation from
    the conformal-weight formula
  - Output paper-quality TeX for hand-tuned identities
  Backed by 21-PDF cohomology bibliography (Bossavit 1998, FEEC 2006,
  Whitney 1957, Kameari 2011, Codecasa 2010, 新しい計算電磁気学 2003 ...).

---

## Validated EM & multiphysics engineering (`radia-ngsolve`)

Most CAE-MCP servers stop at *generating* a script. `radia-ngsolve` goes further: an
AI assistant can **solve a real electromagnetics / multiphysics engineering problem in
open-source [NGSolve](https://ngsolve.org), and every capability is cross-checked against
a closed-form analytic solution**. The result is a library you can trust an AI to drive.

### 20+ "COMSOL-class" problems, each validated and baked in

Each model ships as a **reusable helper + runnable example + regression test + queryable
knowledge** (`examples/comsol_class/`, `ngsolve_usage(...)`), so the server gets smarter,
not just a pile of scripts. Agreement with the closed form is typically **well under 1 %**:

| Domain | Examples (validation vs closed form) |
|--------|--------------------------------------|
| **Electrostatics** | 3-D capacitance (<0.2 %), layered dielectrics (0.1 %), capacitance matrix (0.03 %), electrostatic force (0.4 %), two-wire line (0.1 %) |
| **Magnetostatics** | μ-metal shielding sphere/cylinder (0.9–1.3 %), Halbach PM dipole (0.3 %), cylinder magnet on-axis (centre 0.01 %), finite solenoid + Nagaoka inductance (0.5 %), Helmholtz uniformity (<1 %), iron-yoke dipole gap field (0.5 %) |
| **Field quality** | accelerator-magnet **multipoles** b_n / a_n — normal-quad main term **0.02 %**, allowed/forbidden harmonics resolved |
| **Eddy / AC** | round-wire skin effect R_ac/R_dc vs Kelvin ber/bei (0.07 %), induction heating |
| **Force / torque** | weighted Maxwell-stress ("eggshell") force & torque, busbar Lorentz force (~1 %) |
| **Multiphysics** | **electro-thermal** Joule heating (exact), **electro-thermo-mechanical** thermal-stress chain (exact), **magneto-mechanical** Lorentz→beam deflection vs Euler-Bernoulli (0.02 %) |

Many cases are additionally cross-checked, internally, against a reference commercial FEM
solver (three-way agreement analytic = open-source = commercial); the published numbers
above are all against the **closed-form analytic**, the unimpeachable reference.

### Designer-facing, not textbook

The helpers take what an engineer actually has — material, geometry, excitation — and
return what they want: capacitance, shielding factor, multipole spectrum, inductance,
force/torque, temperature rise, thermal stress, deflection. Ask the MCP server
`ngsolve_usage("field_quality")`, `("solenoid")`, `("c_magnet")`, `("elasticity")`, … for
the validated recipe (70+ topics), or `lint_radia_script` to catch the known FEM traps
before they cost a debug session.

### How it's built

Pure-Python NGSolve (H1 / HCurl / axisymmetric H1Henrotte / VectorH1 elasticity), small
reusable solvers chained for couplings, fast pytest contracts for package health, and
explicit validation scripts under `validation/` that lock every number above. Open
source (BSD-3); validated against closed-form analytics, with commercial solvers used
only as an internal benchmark.

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

## ★ Discovery — start here

Before calling any specific server, ask **`mcp-server-radia-meta`** which
server has the knowledge you need. It is the authoritative cross-server
catalog of all 36 radia-mcp servers and answers "which tool covers
concept X?" without trial-and-error.

```python
# 4 catalog tools + 1 health probe
radia_mcp_overview()                  # all 36 servers + 43 tags
radia_mcp_get("optuna")               # full info for one server
radia_mcp_by_tag("optimization")      # filter (returns 6 matches)
radia_mcp_related("optuna")           # cross-link map (optuna → bayesian-opt, ...)
radia_mcp_health()                    # importability probe of all 36
```

Then drill into a specific server with its `<short>_status()` (auto-
introspected tool list + dep probe) and `<short>_topics()` (for
dispatcher-style servers: the topic enum) tools.

**Pattern**: 3-call discovery instead of guess-and-error —

```
radia_mcp_by_tag("optimization")
  → [optuna, bayesian-opt, evolutionary, topology-optimization,
     data-assimilation, gnn, pinn]
  → optuna_status()              # confirm tools available
  → optuna_topics()              # ['usage', 'algorithm',
                                  'lab_applications', 'all']
  → optuna(topic='algorithm')    # get the content
```

## MCP servers

The catalog is the **source of truth** — call `radia_mcp_overview()` for
the full 36-server live list. The historically primary servers are
shown below for reference; everything else is discoverable via meta.

### Standalone (no Radia core dependency — `pip install radia-mcp`)

| Server | Entry point | Tools | Highlights |
|---|---|---|---|
| **★ meta** | `mcp-server-radia-meta` | 6 | Cross-server catalog + health probe — RECOMMENDED FIRST CALL |
| **literature-index** | `mcp-server-literature-index` | 9 | Full-text search across 2,339 lab literature files in W:/03_文献・論文 (ChromaDB + semantic search) |
| **Cubit** | `mcp-server-cubit` | 52 | `cubit_mesh_auto`, `cubit_exec_safely`, `cubit_ask`, scheme ladder + geometry split, .cub5 checkpoint/restore, scrape index over Coreform forum + S:\\CoreformCubit lab archive (787 files) + YouTube + Coreform training |
| **build123d** | `mcp-server-build123d` | 37 | `build123d_to_cubit_hex`, `lint_build123d_script` (7 rules), `build123d_try` (subprocess isolation), `build123d_inspect_step`, `build123d_heal`, `build123d_api`, 13 Radia/general templates, CadQuery + bd_warehouse interop |
| **gmsh-post** | `mcp-server-gmsh-post` | 21 | `gmsh_post_inspect/validate/convert` (MSH v4.1 only), `gmsh_post_quality`, scrape over GitLab issues + StackOverflow + YouTube |
| **differential-forms** | `mcp-server-differential-forms` | 8 | Differential geometry for computational EM: tangent spaces, k-forms, wedge product, exterior derivative, **Hodge star, Whitney complex, de Rham, tree-cotree, FEEC (Arnold-Falk-Winther 2006)**, Mathematica recipes for symbolic verification. Distilled from Bossavit 1998 + Whitney 1957 + Kameari 2011 + 新しい計算電磁気学 2003 + Codecasa 2010. |
| **mathematica** | `mcp-server-mathematica` | 10 | Wolfram Mathematica subprocess bridge: `mathematica_evaluate` + 9 high-level helpers (simplify, to_tex, vector_calc, unit_convert, solve, integrate, differentiate, check_identity, status). Pairs with `differential-forms` for symbolic verification of d²=0, Stokes, Whitney elements, Kelvin transform, Maxwell identities. Requires `wolframscript` on PATH. |

### Radia-coupled (`pip install radia-mcp[radia]`)

| Server | Entry point | Highlights |
|---|---|---|
| **radia-ngsolve** | `mcp-server-radia-ngsolve` | NGSolve FEM/BEM, Whitney elements via H1/HCurl/HDiv, PEEC inductance, **closed-form formulas** (Wakao-Igarashi Part 1-9, cuboid average B, Bessel impedance, etc.) |
| **ih** | `mcp-server-ih` | Induction heating workflow: workpiece SIBC, ESIM nonlinear cell problem, Karl iteration, screening physics |
| **peec** | `mcp-server-peec` | PEEC: Loop-Star, FastHenry, PyPEECBuilder, Bessel/Dowell/ESIM SIBC, PRIMA model-order reduction, SPICE extraction |
| **electromagnet** | `mcp-server-electromagnet` | Accelerator magnet: CoilBuilder, Hantila polarization, B-input Play/Energy hysteresis, IMA sign selection, multipole harmonics |

---

## Claude Code / Desktop configuration

Add to your MCP client config (Claude Desktop, Claude Code, Cursor,
Continue, …):

```json
{
  "mcpServers": {
    "radia-meta":          {"command": "mcp-server-radia-meta"},
    "literature-index":    {"command": "mcp-server-literature-index"},
    "cubit":               {"command": "mcp-server-cubit"},
    "build123d":           {"command": "mcp-server-build123d"},
    "gmsh-post":           {"command": "mcp-server-gmsh-post"},
    "differential-forms":  {"command": "mcp-server-differential-forms"},
    "mathematica":         {"command": "mcp-server-mathematica"}
  }
}
```

Registering more than the 8 above is rarely necessary — once you have
**meta**, `radia_mcp_get(name)` returns the entry point for any of the
36 servers and you can register them on demand. The full list of
catalog-driven server names is in
`radia_mcp.meta.catalog.CATALOG`.

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
PySide6 `QTimer` (200 ms poll) inside Cubit's Qt event loop. The MCP
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
