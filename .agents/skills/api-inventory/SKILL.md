---
name: api-inventory
description: Take inventory (stock-take) of the radia package's PROPRIETARY Python+pybind API surface and categorize every API FAMILY per the "Reduce Proprietary API Surface" policy -- plumbing (delete/delegate to netgen/ngsolve/MKL/OCC/GMSH/Cubit) vs method (keep; maybe demote to internal) vs user-intent (keep/promote) vs deprecated (drop) -- producing a committed audit doc + a phased reduction plan. Use after adding/removing public APIs, before a deprecation / un-pybind pass, or periodically to drive the surface toward the 2-layer (intent / internal) target. READ-ONLY analysis: safe to run while another agent edits source; do the actual removals SEPARATELY, avoiding active-dev areas (e.g. HACApK).
---

# api-inventory

The radia API surface is large (~200 pybind `.def` in `src/lib/radia_pybind.cpp`
+ ~95 top-level python modules under `src/radia/` + subpackages
`vim/bem/open_boundary/maglev/analytical_formulas/tools` + the
`axifem/cln_core/peec_matrices/sparsesolv_ngsolve` `.pyd`s). The
**"Reduce Proprietary API Surface"** policy (CLAUDE.md / AGENTS.md, 2026-06-19)
wants it driven toward a **2-layer shape**:

- a lean **intent-based USER layer** -- `SoftIron("yoke.vol", mu_r=)` / `Magnet` /
  `CoilBuilder` / `rad.Fld` / `rad.Solve` / materials;
- over an **internal layer** of the genuine METHODS NGSolve lacks -- `rad.Fld`
  analytic open-boundary field, MMM/MSC, yano-MSC, HDiv-VIM, axifem (Henrotte),
  DtN/FEM-Kelvin, PEEC, BEM, sparsesolv, HACApK, analytical_formulas,
  coil_builder, levitation, stream-function;
- with all **PLUMBING delegated** to netgen/ngsolve (mesh gen/IO, geometry/CAD,
  visualization/mesh export, generic linalg) and the geometry PRIMITIVES
  (`ObjHexahedron`/`ObjRecMag`/...) demoted from the user's hand-built-mesh API
  to an internal representation behind `.vol` -> `soft_iron_from_mesh` + intent
  objects.

This skill takes stock so the reduction proceeds deliberately, not by guesswork.

## When to use
- After adding or removing public APIs (keep the surface honest).
- Before a deprecation / un-pybind pass (know what is plumbing vs method first).
- Periodically, to measure progress toward the 2-layer target.

## The five buckets (decision rule per API family)
For each API FAMILY (group by name prefix / role -- do NOT enumerate 200
functions) assign exactly one bucket:

| Bucket | Meaning | Action |
|---|---|---|
| **plumbing-delete** | netgen/ngsolve/MKL/OCC/GMSH/Cubit already provides it | delete; delegate; name the replacement |
| **method-keep** | a genuine method NGSolve lacks (Radia's reason to exist) | keep as-is |
| **method-demote** | a kept method that should be un-pybind'd behind the intent layer over time | keep != expose; demote gradually |
| **user-intent** | the intended lean user layer | keep / promote |
| **deprecated-drop** | already-removed shells, back-compat shims, dead/legacy, duplicates | retire (per "No Development Cruft in SOURCE") |

The discriminator is plumbing-vs-method, NOT uses-ngsolve-vs-not: a method built
ON ngsolve (axifem, HDiv-VIM) is still a method -- KEEP.

## Step 1 -- scout the surface (read-only; refresh the workflow's AREAS lists)
```bash
ls src/radia/*.pyd                                   # compiled modules
ls src/radia/*.py                                    # python module surface
ls -d src/radia/*/                                   # subpackages
grep -c '\.def' src/lib/radia_pybind.cpp             # pybind family count
sed -n '1,80p' src/radia/__init__.py                 # re-export hub + 2-layer wrappers
```
(PowerShell: `Get-ChildItem 'src/radia/*.pyd'`, `Select-String -Path src/lib/radia_pybind.cpp -Pattern '\.def'`.)
If the module layout changed since the bundled workflow was written, update the
`AREAS` file lists in `inventory_workflow.js` to match.

## Step 2 -- categorize (workflow for the full surface)
Run the bundled fan-out workflow (5 read-only area readers -> synthesis). This
requires the multi-agent opt-in (ultracode or an explicit user request):
```
Workflow({scriptPath: ".agents/skills/api-inventory/inventory_workflow.js"})
```
It returns: `counts_by_bucket`, `top_deletion_candidates`,
`top_demotion_candidates`, `recommended_user_layer`, `safe_now_vs_blocked`, and a
ready-to-commit `markdown_report`. For a small/targeted check (one module), a
single inline agent with the five-bucket framework above is enough -- no workflow.

## Step 3 -- the deliverable
Write the synthesis `markdown_report` to a dated, committed audit doc:
`docs/api_inventory/API_SURFACE_<YYYY-MM-DD>.md`. Dating it makes the next run a
diff ("what shrank since last time"). Commit it (a tracked doc, not in an
active-dev area).

## Step 4 -- act (separately, carefully)
- **plumbing-delete / deprecated-drop**: remove and delegate to the named
  netgen/ngsolve replacement. Per "No Development Cruft in SOURCE": keep one
  canonical version, distill any lesson to `memory/` first, recover from git if
  needed.
- **method-demote**: un-pybind gradually behind the intent layer -- panels /
  examples / CoilBuilder depend on the primitives today, so DEMOTE first, remove
  after migration.
- **AVOID active-dev areas.** Categorize HACApK (`src/ext/HACApK/`,
  `src/core/rad_hacapk.*`) but never edit it while another agent works it; same
  for any file in the current uncommitted `M` set. The inventory is read-only and
  always safe; only the Step-4 edits carry collision risk.

## Related
- "Reduce Proprietary API Surface" + "No Development Cruft in SOURCE" policies (CLAUDE.md / AGENTS.md).
- `inventory_workflow.js` -- the bundled, proven fan-out implementation.
- `ipynb-gui-health`, `panel-cli-diff`, `panel-review` -- adjacent surface-health skills.
