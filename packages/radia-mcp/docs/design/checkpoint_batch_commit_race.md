# Checkpoint-Batch-Commit-Race: a design pattern for AI-driven interactive CAE

| | |
|---|---|
| **Author** | Kengo Sugahara (菅原 賢悟) — Sugawara Lab, Kindai University |
| **Contact** | ksugahar@ele.kindai.ac.jp |
| **First public disclosure** | 2026-04-20, this repository, BSD-3-Clause |
| **Reference implementation** | [`radia-mcp`](https://pypi.org/project/radia-mcp/) ≥ 0.25.0 |
| **Status** | Proposed pattern, in production use at the Sugawara Lab |
| **License of this document** | BSD-3-Clause (matching the project) |

> **Defensive publication.** This document, together with the
> reference implementation in `radia-mcp` (publicly released on PyPI
> on 2026-04-20), constitutes prior-art disclosure of the pattern
> described below. The Sugawara Lab does not claim a patent and
> intends the pattern to remain freely usable.

---

## Abstract

Interactive computer-aided engineering (CAE) tools — Coreform Cubit
for hex meshing, FreeCAD for assembly, ngsolve for FEM — have a
universal user-interface pain point: when an operation fails or
stalls, the user must wait for the AI assistant to suggest a fix
before they can do anything else, or the AI overwrites the user's
in-progress work. This document proposes
**Checkpoint-Batch-Commit-Race (CBCR)**, a four-stage cooperative
pattern that lets a human operator and *N* AI agents work
concurrently on the same starting state, with first-to-finish
semantics and a hard guarantee that the human's work is never
silently overwritten. CBCR is implemented in the open-source
`radia-mcp` Model Context Protocol server suite for Cubit, gmsh,
and build123d, and is portable to any CAE tool with a
snapshot/restore primitive (`.cub5`, `.FCStd`, git, …).

---

## 1. Problem statement

The dominant interaction loop with AI assistants in CAE software is
serial:

```
user: "this geometry won't mesh, can you try?"
AI:   ...thinking...
AI:   "try this scheme"
user: ...waits, then runs it...
result: failure or success
```

There are two well-known failure modes:

1. **Idle-waiting**: the user stops working while the AI thinks. For
   a meshing scheme search that takes 30 s × 4 candidates = 2 min,
   this is wasteful — the user could have meshed the geometry
   manually in the same time.
2. **Silent overwrite**: when the AI eventually returns, naïvely
   replaying its commands on the live session destroys whatever the
   user did during those 2 minutes (rerun the same `mesh volume all`
   that just succeeded by hand, lose mesh-size tweaks, etc.).

Existing AI-pair-programming tools (GitHub Copilot, Cursor) sidestep
this by working at the text-buffer level where a CRDT or simple line
diff can merge edits. **CAE tools have no such structure**: the
"document" is a multi-megabyte solid-model database with tightly
coupled topology, mesh element ids, and physical-group tags.
Naïve merge is impossible.

CBCR is the lightweight alternative: snapshot, race, commit only on
clear win, never overwrite ongoing human work.

---

## 2. The pattern

CBCR has four stages. Each stage corresponds to an MCP tool in the
reference implementation; the names below match `radia-mcp` 0.25.0.

### 2.1 Stage 1 — Checkpoint

Before any risky operation, snapshot the current live CAE-tool state
to an immutable file using the tool's own native save mechanism.

| Tool | Native snapshot | radia-mcp helper |
|---|---|---|
| Coreform Cubit | `save as "<path>.cub5" overwrite` | `cubit_checkpoint(label)` |
| FreeCAD | `Document.saveAs(path).FCStd` | `freecad_exec_safely`'s internal copy |
| build123d | git commit (script is the state) | not needed; stateless |

The snapshot label is returned to the caller so a manual revert is
always available (`cubit_restore(label)` in the reference impl).

### 2.2 Stage 2 — Batch (parallel race)

Spawn *N* disposable batch instances of the same CAE tool from the
snapshot. Each instance runs one *recipe* (a sequence of commands)
in isolation. Subprocess isolation matters for two reasons:

- **Crash containment**: a degenerate-geometry crash on one variant
  doesn't kill the others.
- **State isolation**: each variant starts from the identical
  snapshot, so per-recipe outcomes are independently comparable.

In the reference implementation, batch Cubit instances are launched
with `coreform_cubit.exe -batch -nographics -nojournal` running a
small Python daemon (`cubit_daemon.py`) that exchanges JSON-RPC
messages over stdio. Concurrency is a `ThreadPoolExecutor` whose
size defaults to `max(1, ncpu // 2)` — the same machine usually
hosts both the live GUI and the batch race, so saturating the CPU
hurts the human's responsiveness.

### 2.3 Stage 3 — Live-state polling (the human side of the race)

Concurrently with the batch race, a polling thread interrogates the
live CAE session at a configurable interval (default 5 s). The
poll asks for a small, side-effect-free metric — for Cubit it is
`probe summary` (volumes / surfaces / nodes / hexes / tets). The
thread compares each sample against the pre-race baseline.

The polling is what makes CBCR fundamentally different from plain
parallel hyperparameter search: **the human is one of the racers**.
If the human meshes the geometry by hand before any AI variant
finishes, the polling thread detects the new mesh and claims
victory.

### 2.4 Stage 4 — Commit (with non-overwrite guarantee)

A shared `winner` variable is protected by a mutex. The first
batch-thread or polling-thread to satisfy the success predicate
(`hex > pre_hex` for `prefer="hex"`) atomically claims it; any
later claim is ignored.

Commit policy:

| Winner | Default action | `commit_winner=True` action |
|---|---|---|
| Human | report only | report only |
| AI variant | report only | replay recipe in live, **iff** human hasn't meshed since |

The "iff human hasn't meshed since" guard re-probes the live state
right before the commit and refuses if the live mesh has already
grown beyond the pre-race baseline. **The user's work is sacred**
— even an explicit `commit_winner=True` request will not overwrite
human progress.

---

## 3. Reference implementation

`radia-mcp` 0.25.0 ships
[`cubit_mesh_race_with_human`](https://github.com/ksugahar/Radia/blob/main/packages/radia-mcp/src/radia_mcp/cubit/server.py)
as the canonical CBCR tool. The signature:

```python
@mcp.tool()
def cubit_mesh_race_with_human(recipes: list,
                               poll_interval_s: float = 5.0,
                               max_wait_s: int = 600,
                               max_concurrent: int = 4,
                               commit_winner: bool = False,
                               prefer: str = "hex") -> str: ...
```

### 3.1 Plan A — the live-Cubit transport

The polling thread can only ask the live Cubit because we control
the live session via **Plan A**: a custom `cubit_bootstrap.py`
runs inside Cubit's GUI process, installed by
`coreform_cubit.exe -nojournal cubit_bootstrap.py`. The bootstrap
registers a `PySide6.QtCore.QTimer` (200 ms) that polls a temp drop
directory for `*.req.json` files, runs `cubit.cmd()` on the Qt
main thread, and writes responses to `out/*.resp.json`. Atomically
renamed JSON files; no sockets, no pipes, no ABI conflicts.

This in-GUI-Python design is itself non-novel (FreeCAD-MCP uses
XML-RPC inside the FreeCAD process; Blender add-ons routinely
register polling timers), but it is *necessary* for CBCR because
without a live-state interrogation channel the polling thread has
no way to detect human progress.

### 3.2 Algorithm pseudocode

```python
def cbcr(live, recipes, prefer="hex", poll_s=5):
    snapshot_path, label = take_snapshot(live)
    pre = probe(live)
    winner_event = Event()
    winner = [None]  # mutable for closure
    lock = Lock()

    def claim(payload):
        with lock:
            if winner[0] is None:
                winner[0] = payload
                winner_event.set()
                return True
        return False

    def poll_human():
        while not winner_event.is_set():
            if probe(live) > pre by `prefer` predicate:
                claim({"who": "human", "summary": probe(live)})
                return
            winner_event.wait(poll_s)

    def run_ai(recipe):
        result = batch_run(snapshot_path, recipe)
        if result satisfies prefer:
            claim({"who": "ai", "name": recipe.name,
                   "recipe": recipe.cmds, **result.counts})
        return result

    with ThreadPoolExecutor() as pool:
        pool.submit(poll_human)
        for r in recipes:
            pool.submit(run_ai, r)
        winner_event.wait(timeout=max_wait_s)

    return winner[0], snapshot_path, label
```

The MCP tool wraps this with optional commit-on-AI-win, the
non-overwrite re-probe guard, and structured JSON output for
the LLM caller.

---

## 4. Worked example — 4-turn spiral coil hex mesh

### 4.1 Setup

A 4-turn spiral coil with two radial leads (parameters: R=15 mm,
pitch=8 mm, lead length=20 mm, square cross-section 2×2 mm)
exported as STEP from build123d (3 separate prismatic bodies
joined topologically by `imprint all; merge all`).

Imported into a live Cubit GUI via `import step "spiral.step"`
(3 volumes, 18 surfaces post-imprint).

### 4.2 The race

```python
cubit_mesh_race_with_human(
    recipes=[
        {"name": "auto_size_1.0",       "cmds": [
            "volume all size 1.0",
            "volume all scheme auto",
            "mesh volume all"]},
        {"name": "sweep_size_1.0",      "cmds": [
            "volume all size 1.0",
            "volume all scheme sweep",
            "mesh volume all"]},
        {"name": "polyhedron_size_1.5", "cmds": [
            "volume all size 1.5",
            "volume all scheme polyhedron",
            "mesh volume all"]},
    ],
    poll_interval_s=3.0,
    commit_winner=True,
)
```

Meanwhile the user, in the live Cubit GUI window, can keep
typing: experimenting with `volume 2 scheme tetmesh`, adjusting
mesh size, undoing, redoing — completely unaware of the race in
the background.

### 4.3 Outcome (measured 2026-04-19)

The `auto_size_1.0` variant finishes in ~12 s with 1668 hex / 0 tet
/ 3780 nodes — *before* the user has time to type a single Cubit
command. Because `commit_winner=True` AND the user has not yet
meshed (live `hexes == 0`), the recipe is replayed in the live GUI;
the user sees the mesh appear in their window. Total time from
"AI, try meshing" to mesh-on-screen: **~14 s**, two of which were
the snapshot save.

If the user had been faster — say, manually hitting the toolbar
"Mesh" button while the race ran — the polling thread would have
detected the live mesh first; the AI batches would still complete
but their results would be discarded; the user's manual mesh
remains; the snapshot label is returned in case they want to undo
later.

---

## 5. Variants

The same pattern lifts cleanly to other CAE tools:

### 5.1 `build123d_try_race`

build123d is stateless (each `execute_build123d` is a fresh Python
namespace), so the snapshot stage is degenerate (the script *is* the
state). Variants are different scripts (e.g., a parameter sweep
over fillet radii); they run in subprocess pool; the winner is
chosen by validity + largest volume (`valid_largest_volume`),
first-valid (`first_valid`), or first-no-error (`first_ok`). No
human race because build123d's natural editing surface is the
source code, version-controlled by git.

### 5.2 `freecad_exec_safely`

Implements the snapshot + batch-dry-run + commit half of the
pattern (no live-poll race, since the lab does not author in
FreeCAD — see [README §lab stance](../../README.md)). Snapshot
is `.FCStd` copy; batch is `FreeCADCmd`; commit is file replace.

### 5.3 Other CAE / CAD tools

The pattern requires only:

- A snapshot/restore primitive (filesystem-level is fine).
- A way to launch a disposable instance from the snapshot.
- A small read-only metric to poll the live instance with.

Tools that satisfy these natively today: Coreform Cubit, FreeCAD,
ngsolve (via .vol round-trip), Salome, Abaqus/CAE (via .cae journal
replay), Ansys Workbench (via Mechanical APDL save). Tools that
would need an in-process polling adapter: Blender (add-on),
SolidWorks (COM API), Onshape (REST API; trivially polled).

---

## 6. Prior art

CBCR composes well-known building blocks:

- **Parallel hyperparameter racing** (Hyperband [Li et al. 2017],
  ASHA [Li et al. 2020], BOHB [Falkner et al. 2018]). These race
  *N* configurations and pick the empirical winner; they have no
  human in the loop, no live-state polling, no commit guard.
- **AI pair programming** (GitHub Copilot, Cursor, Continue, Cline).
  The AI suggests at the cursor; the human accepts or rejects.
  Concurrency is left-to-right text editing with simple line
  merge; no race, no batched alternative attempts, no
  non-overwrite primitive.
- **Computer-supported cooperative work (CSCW)** literature on
  CRDTs (Shapiro 2011) and operational transformation (Ellis &
  Gibbs 1989) handles concurrent text edits with eventual
  consistency. CBCR is much simpler because we *refuse* to merge
  — the winner takes all, the loser's work is preserved
  out-of-band by the snapshot.
- **Speculative execution** in compilers / CPUs / databases.
  Logically related (commit only on clear win) but operates on
  bounded operations, not multi-second user-visible state.
- **Live-coding environments** (LightTable, Pharo) and notebook
  parallelism (Jupyter). Different scope: serial, single-actor.
- **MCP ecosystem** (Anthropic 2024). CBCR is implemented as an
  MCP tool but is not specific to MCP; the same protocol could
  expose it via JSON-RPC, gRPC, REST, or direct Python import.

To the author's knowledge, **no existing CAE tool implements the
human-vs-AI race with snapshot-protected non-overwrite semantics
described here**. (The closest is FreeCAD's transaction system,
which handles undo/redo for a single actor.)

---

## 7. Why defensive publication, not patent

The Sugawara Lab evaluated three options — defensive publication,
JP provisional patent, or trade secret — and chose option 1 for
the following reasons:

1. **Cost**: defensive publication is free; a patent costs
   ¥100k–¥300k just to file in JP, and orders of magnitude more
   for international coverage.
2. **Mission fit**: a research lab's job is to publish, not to
   monitor and litigate infringement.
3. **Ecosystem health**: AI + interactive CAE is a young field;
   establishing the pattern as freely usable accelerates everyone
   (including the lab's own follow-on work).
4. **Already disclosed**: `radia-mcp` 0.25.0 was published on PyPI
   on 2026-04-20. The 6-month JP grace period (特許法 30 条) is
   the only patent option still alive; the EU and US public-
   disclosure clocks have effectively closed.
5. **Citation > royalty**: a paper or design-doc citation is the
   academic currency that matters for a university lab.

The license is BSD-3-Clause; the pattern is intended to be
freely adopted by FreeCAD-MCP, OpenSCAD-MCP, Blender-MCP,
Ansys/SolidWorks/Onshape MCPs, and any future CAE-MCP author.

---

## 8. Citation

If you use or describe this pattern, the suggested citation is:

> Sugahara, K. (2026). *Checkpoint-Batch-Commit-Race: a design
> pattern for AI-driven interactive CAE*. Sugawara Lab, Kindai
> University. radia-mcp documentation, 2026-04-20.
> https://github.com/ksugahar/Radia
> https://pypi.org/project/radia-mcp/0.25.0/

For BibTeX:

```bibtex
@misc{sugahara2026cbcr,
  author       = {Sugahara, Kengo},
  title        = {{Checkpoint-Batch-Commit-Race: a design pattern
                   for AI-driven interactive CAE}},
  year         = {2026},
  month        = {April},
  howpublished = {radia-mcp documentation, Sugawara Lab,
                   Kindai University},
  url          = {https://github.com/ksugahar/Radia},
  note         = {radia-mcp v0.25.0 reference implementation,
                   PyPI, 2026-04-20.},
}
```

---

## 9. Future work

- **Concurrent merge**: when both human and AI produce different
  meshes within the race window, present a side-by-side comparison
  + quality metric (`gmsh_post_quality`) for the user to choose.
  Genuinely novel; today CBCR refuses to merge.
- **Cross-tool race**: race a Cubit recipe against a Gmsh recipe
  against a NETGEN recipe for the same STEP. The pattern
  generalizes to "first valid mesh from any backend wins".
- **Persistent race history**: log every race outcome to the
  failure-log infrastructure so future invocations can prefer
  variants that historically won faster on similar geometries.
- **Web-UI for the race**: a small dashboard showing live progress
  of *N* AI variants + the human's editor, with a "claim" button.

---

## Appendix A — Full cubit_mesh_race_with_human source

See [`packages/radia-mcp/src/radia_mcp/cubit/server.py`](../../src/radia_mcp/cubit/server.py),
function `cubit_mesh_race_with_human` (≈ 130 lines). License:
BSD-3-Clause.
