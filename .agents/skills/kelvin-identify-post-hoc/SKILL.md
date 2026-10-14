---
name: kelvin-identify-post-hoc
description: Add Kelvin Periodic Identifications to an existing NGSolve mesh AFTER load (no Cubit launcher / OCC Identify needed). Use this when a `.vol` file has `kelvin_int` / `kelvin_ext` boundary labels but no `Identifications` section -- typical when the Cubit C++ exporter's all-or-nothing vertex matching skipped due to a single tolerance miss, or when the mesh came from outside the Cubit panel pipeline.
---

# kelvin-identify-post-hoc

Add Kelvin point-pair `Identifications` to an already-loaded NGSolve
mesh without going through `add_kelvin_cubit` (Cubit launcher) or OCC
`Identify(IdentificationType.PERIODIC)` (which has to be called
BEFORE `GenerateMesh`).

## When to invoke this skill

- User has a `.vol` with `kelvin_int` and `kelvin_ext` labels but
  loading the mesh + building `Periodic(H1(...))` gives `slaved=0`
- `calc_fem_kelvin.py --require-kelvin` fails with "kelvin material
  missing" but the .vol DOES have kelvin labels — likely the
  Identifications section is missing
- User asks to "add Identify to an existing .vol"
- User pastes an `add_kelvin_cubit` script that skipped C++
  identification on a 1/8-octant geometry

## Phase 1: Diagnose the mesh

Before adding identifications, confirm the mesh has the prerequisites.

```python
from ngsolve import Mesh
from radia import (has_kelvin_identification,
                   detect_kelvin_offset)

mesh = Mesh(vol_path)
print(f"materials:  {mesh.GetMaterials()}")
print(f"boundaries: {mesh.GetBoundaries()}")
print(f"already has identifications: {has_kelvin_identification(mesh)}")
print(f"detected kelvin offset:      {detect_kelvin_offset(mesh)}")
```

**Pass conditions** (all required):

| Check | Required value | Failure means |
|---|---|---|
| `kelvin_int` in boundaries | True | Rebuild mesh with `add_kelvin_cubit` or OCC Identify |
| `kelvin_ext` in boundaries | True | Same |
| `has_kelvin_identification` | False | Skip Phase 2 — the mesh already works |
| `detect_kelvin_offset != (0,0,0)` | Non-trivial vector | Two spheres at same centre — degenerate |

If `has_kelvin_identification` returns True, the work is already done.
Verify with Phase 3 directly.

## Phase 2: Add identifications

```python
from radia import add_kelvin_identification

info = add_kelvin_identification(mesh)   # auto-detects offset
print(f"[kelvin] {info['n_pairs']} pairs, "
      f"max_dist={info['max_dist']:.2e} m, "
      f"unmatched={info['n_unmatched']}")
```

**Pass conditions**:

| Diagnostic | Pass | Fail means |
|---|---|---|
| `info['n_pairs'] > 0` | required | No vertices matched — wrong offset or labels |
| `info['n_unmatched'] == 0` | preferred | Some vertices outside tolerance — mesh is not strictly 1:1 |
| `info['max_dist'] < 1e-3` | preferred | Worse-than-mm jitter; check tolerance / offset |

**If `n_unmatched > 0`**: the two surfaces are NOT a clean rigid
translate of each other.  This is the helper's load-bearing
pre-condition (see memory `feedback_kelvin_1_8_blocker.md`).
Options:
- Pass an explicit `kelvin_offset` if auto-detect is wrong
- Increase `point_tolerance` (default = `max(5%*|offset|, 1mm)`)
- Rebuild the mesh via `add_kelvin_cubit(symmetry=[...])` or OCC
  `Identify(...)` before `GenerateMesh` — these guarantee 1:1

**If `info["was_pre_existing"] is True`**: the mesh already has
identifications.  Don't re-add unless you explicitly know they're
broken (then pass `skip_if_existing=False`).

## Phase 3: Verify FES coupling

Two sub-second checks, both must pass.  Run them BEFORE any expensive
physics solve (AGENTS.md "Verify-First Policy").

### Check 1: slaved DOF count

```python
from ngsolve import H1, Periodic

fb = H1(mesh, order=1, dirichlet="GND")
fp = Periodic(fb)
slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
assert slaved > 0, f"Periodic slaved 0 DOFs"
```

`slaved` should be roughly equal to the number of vertices on the
kelvin boundary (minus the GND vertex).

### Check 2: functional boundary test

```python
from ngsolve import GridFunction, Integrate

fp = Periodic(H1(mesh, order=1, dirichlet="GND"))
gfu = GridFunction(fp)
gfu.vec[:] = 0
gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))

a_int = float(Integrate(gfu*gfu, mesh,
              definedon=mesh.Boundaries("kelvin_int")))
a_ext = float(Integrate(gfu*gfu, mesh,
              definedon=mesh.Boundaries("kelvin_ext")))
ratio = a_ext / a_int
assert abs(ratio - 1.0) < 1e-3, f"Periodic ratio not 1.0: {ratio}"
```

If `ratio` is not 1.0, the identification did not slave the boundary
correctly.  Likely cause: the vertex pairing was structurally wrong
(see "If n_unmatched > 0" above).

## Phase 4: Persist (optional)

NGSolve does NOT save identifications back to `.vol` automatically
when you read + modify a mesh.  If you want the identifications to
travel with the file, save explicitly:

```python
mesh.ngmesh.Save("model_with_kelvin_idents.vol")
```

Caveat: NGSolve `Identifications` cannot be cleared from Python (the
list returned by `GetIdentifications()` is a copy).  If the saved
.vol gets re-loaded and re-identified by code that does NOT pass
`skip_if_existing=True`, duplicate pairs accumulate.  Default
`skip_if_existing=True` in `add_kelvin_identification` avoids this.

## Tools available

| Tool | Source |
|---|---|
| `radia.add_kelvin_identification` | `src/radia/kelvin_identify_ngsolve.py` (public) |
| `radia.detect_kelvin_offset` | same module |
| `radia.has_kelvin_identification` | same module |
| MCP: `kelvin_identify_post_hoc(topic, vol_path)` | `radia-mcp` package |
| Test fixtures | `tests/panels/test_kelvin_identify_post_hoc.py` |

The MCP tool's `topic="verify"` mode with a `vol_path` argument runs
Phase 1 + Phase 2 in a single call and returns a JSON diagnostic.

## Related knowledge

- Theory: `kelvin_transformation` MCP tool (Kelvin inversion, material
  scaling, formulations).  This skill is about the IDENTIFY step;
  the transformation tool covers the material side.
- 1:1 mesh blocker: memory `feedback_kelvin_1_8_blocker.md` (deterministic
  copy-mesh anchor selection).
- Detect-offset bug history: memory + AGENTS.md "Verify-First Policy"
  (pre-2026-04-25 implementation returned kelvin REGION centroid
  instead of `mean(kelvin_ext) - mean(kelvin_int)` — wrong by ~14 cm
  on 1/4-octant geometry).
