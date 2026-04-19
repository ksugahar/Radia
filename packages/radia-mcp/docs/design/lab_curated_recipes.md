# Lab-curated mesh recipes — collective intelligence pipeline

| | |
|---|---|
| **Author** | Kengo Sugahara (菅原 賢悟) — Sugawara Lab, Kindai University |
| **Status** | Operating, distillation runs ad-hoc (target: weekly) |
| **License of this document** | CC-BY-4.0 (recipes themselves are BSD-3-Clause via the wheel) |
| **First public disclosure** | 2026-04-20, this repository |
| **Reference implementation** | [`radia-mcp`](https://pypi.org/project/radia-mcp/) ≥ 0.30.0 |

---

## 1. Why this exists

The Sugawara Lab's `radia-mcp` 0.25 introduced
[CBCR](checkpoint_batch_commit_race.md): a human and *N* AI agents
race to mesh a stuck geometry; first-to-finish wins; the human's
work is never overwritten. By 0.29, every winning recipe (human or
AI) is **persisted** to a per-machine `learned_recipes.jsonl`, and
`cubit_mesh_race_smart_async` automatically replays similar past
winners as candidates next time.

That gives one machine a private learning loop. But the lab has
many machines (LAB / 100号機 / student laptops) and many years of
mesh problems. **A truly smart MCP needs to share knowledge across
machines and across time** — and ideally publish the distilled
knowledge to the wider CAE community.

This document describes the three-stage pipeline that does that.

---

## 2. The pipeline

```
                 (per-machine)
race winners → state_dir/learned_recipes.jsonl    ← Stage 0: 0.29 baseline
                     │
                     ▼  set RADIA_MCP_LEARNED_DIR=//share/lab_learned/
        lab-shared learned_recipes.jsonl           ← Stage 1: 0.30, env var
                     │
                     ▼  cubit_curate_learned_recipes (weekly)
            curated_recipes_bundle.py              ← Stage 2: ships in wheel
                     │  python -m build && twine upload
                     ▼
              radia-mcp on PyPI                    ← Stage 3: world-wide
                     │
                     ▼
       cubit_mesh_race_smart_async (anywhere)
```

### Stage 0 — local race-winner persistence (already shipped, 0.29)

After every race in `cubit_mesh_race_review_async`, the module
extracts the human's command sequence from a per-race Cubit
journal, computes a geometry signature
`{volumes, surfaces, surf_per_vol, curves, vertices}`, and appends
`{when, race_id, source, signature, recipe, quality}` as a JSON
line to `<state_dir>/learned_recipes.jsonl`. The next
`_generate_smart_recipes` call queries similar signatures
(volume count exact + surf/vol ratio within ±20%) and seeds the
race with past winners.

### Stage 1 — lab-shared pool (`RADIA_MCP_LEARNED_DIR`, 0.30+)

```
# All lab machines: set this in shell rc / cubit_session.py env
export RADIA_MCP_LEARNED_DIR='S:\Radia\01_GitHub\lab_learned'
```

`_learned_recipes_path()` checks the env var first; when set, the
jsonl lives on the lab-shared SMB drive instead of the per-machine
state dir. **Every race on every lab machine contributes to one
common pool**, and every machine's smart-recipe generator sees the
union.

Concretely: a student wins a tough mesh on their laptop at 14:00;
by 14:30 the AI on the LAB workstation already lists their winning
recipe as the top candidate for the same geometry class.

The shared file is append-only JSONL (no locking conflicts); concurrent
writes from multiple processes are safe because each line is a
single atomic write.

### Stage 2 — distillation into the wheel (`cubit_curate_learned_recipes`, 0.30+)

```python
cubit_curate_learned_recipes(top_per_class=3, min_quality_jacobian=0.3)
```

Reads `RADIA_MCP_LEARNED_DIR/learned_recipes.jsonl`, drops anything
below the quality threshold, groups by signature class
`(volumes, round(surf_per_vol, 1))`, dedups by recipe text, takes
the top *N* per class by quality, and writes
`packages/radia-mcp/src/radia_mcp/cubit/curated_recipes_bundle.py`.
The next `python -m build` ships it inside the wheel.

The bundled module is loaded on import alongside the local jsonl.
**Fresh `pip install radia-mcp` anywhere in the world picks up the
lab's accumulated wisdom on first use** — no shared filesystem
required.

The curation is intentionally lightweight (no ML), so the
maintainer can audit every entry before the wheel ships. Recommended
cadence: weekly during heavy lab use, monthly otherwise.

### Stage 3 — public release as a CC-BY dataset

The recipes themselves continue to ship under the wheel's
BSD-3-Clause license, but **this methodology document and the
distilled summary statistics** are published under CC-BY-4.0 so
they can be:

- Cited in the Sugawara Lab's papers without complication.
- Reused by other research labs / CAE software vendors for their
  own AI-mesh-assistant designs.
- Compared against future lab-curated datasets from elsewhere.

Suggested first publication target: a short paper / dataset note
to ASME IDETC/CIE or a similar CAE+AI venue, after the bundle has
~100 distilled recipes (estimated 3-6 months of routine lab use).

---

## 3. Operational notes

### Setting up a fresh lab machine

```bash
pip install radia-mcp[full]
# Optional but recommended for lab-shared learning:
setx RADIA_MCP_LEARNED_DIR "S:\Radia\01_GitHub\lab_learned"
# or, on POSIX:
export RADIA_MCP_LEARNED_DIR=/mnt/lab/learned
```

Restart the MCP host (Claude Code, Claude Desktop, …) so the env
var takes effect.

### Curation workflow (lab maintainer)

```bash
# 1. From any machine attached to the shared pool:
python -c "
import radia_mcp.cubit.server as cs
import json
print(json.loads(cs.cubit_curate_learned_recipes(top_per_class=3))['curated_entries'])
"
# 2. Review the diff:
git -C s:/Radia/01_GitHub diff packages/radia-mcp/src/radia_mcp/cubit/curated_recipes_bundle.py
# 3. Commit + bump version + ship:
cd s:/Radia/01_GitHub/packages/radia-mcp
sed -i 's/version = "0.X.Y"/version = "0.X.(Y+1)"/' pyproject.toml
sed -i 's/__version__ = "0.X.Y"/__version__ = "0.X.(Y+1)"/' src/radia_mcp/__init__.py
rm -rf build/ dist/radia_mcp-0.X.*
python -m build && python -m build --wheel
PYTHONIOENCODING=utf-8 TWINE_USERNAME=__token__ TWINE_PASSWORD=$PYPI_TOKEN \
  twine upload --disable-progress-bar dist/radia_mcp-0.X.*
```

### Privacy / contributor consent

Recipe entries record `source` (`human` / `ai_<variant_name>`),
`signature` (geometry shape only — no project-identifying data),
`recipe` (Cubit text commands), and `quality` (numbers).  No
filenames, project names, user names, or identifying metadata are
stored. Contributors using radia-mcp under BSD-3-Clause are
implicitly consenting to recipe collection; lab policy is to keep
the pool internal until distillation, then ship only the distilled
top-N (no raw history) in the public wheel.

If a contributor wants their machine excluded:
`export RADIA_MCP_LEARNED_DIR=$HOME/.private_radia_learned`
(point at a private dir) — local racing still works, recipes just
don't enter the shared pool.

---

## 4. Citation

If you use the curated recipe bundle or describe this pipeline:

```bibtex
@misc{sugahara2026labcurated,
  author       = {Sugahara, Kengo},
  title        = {{Lab-curated mesh recipes: a collective intelligence
                   pipeline for AI-driven CAE meshing}},
  year         = {2026},
  howpublished = {radia-mcp documentation, Sugawara Lab,
                   Kindai University},
  url          = {https://github.com/ksugahar/Radia},
  note         = {radia-mcp v0.30.0+; companion to the CBCR pattern
                   (Sugahara 2026, checkpoint\_batch\_commit\_race.md).},
}
```

---

## 5. Future work

- **Auto-curation in CI**: GitHub Action that runs
  `cubit_curate_learned_recipes` weekly, opens a PR, ships when
  merged. Removes the manual maintainer step.
- **Cross-lab pool**: optional opt-in to share with other research
  groups via a hosted aggregator (Discourse / GitHub LFS).
- **Recipe embedding**: replace surf/vol ratio with a richer
  geometric embedding (BREP feature vector) so similarity matches
  finer-grained shape classes.
- **Per-solver curation**: split bundles by downstream solver
  (radia / ngsolve / Abaqus / JMAG) since recipe quality depends
  on what the mesh feeds.
