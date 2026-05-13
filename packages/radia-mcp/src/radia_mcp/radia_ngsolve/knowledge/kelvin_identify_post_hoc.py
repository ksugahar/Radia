"""Knowledge: post-hoc Kelvin Periodic identification on an existing mesh.

Companion to the `kelvin_transformation` knowledge module (theory of the
Kelvin inversion).  This file documents the PRACTICAL helper that lets
you add the Kelvin point-pair `Identifications` AFTER the mesh is
already loaded, WITHOUT going through the Cubit panel launcher or
re-running OCC `Identify(IdentificationType.PERIODIC)`.

Public API: `radia.kelvin_identify_ngsolve` (also re-exported at
`radia.add_kelvin_identification`, `radia.detect_kelvin_offset`,
`radia.has_kelvin_identification`).
"""

OVERVIEW = """
# Post-hoc Kelvin Identification (launcher-free)

## When to use this

You have an NGSolve mesh with `kelvin_int` / `kelvin_ext` BND labels
but **no `Identifications` section yet**, and you need to enable
`Periodic(H1(mesh, ...))` to couple the two surfaces for the Kelvin
open-boundary truncation.  Reaching this state typically means:

  - The Cubit C++ exporter (`radia_export netgen`) ran but its
    all-or-nothing vertex matching skipped due to a single tolerance
    miss (typical 1/8-octant geometry).
  - The mesh came from somewhere outside the Cubit panel pipeline --
    a colleague's `.vol`, a hand-built OCC repro, or a Cubit Python
    script that did NOT call `add_kelvin_cubit()` (just exported
    geometry + labels).
  - You're scripting from a Jupyter notebook with NGSolve only and
    no Cubit license available.

## The 3-line solution

```python
from ngsolve import Mesh
from radia import add_kelvin_identification

mesh = Mesh("my_model.vol")              # has kelvin_int / kelvin_ext labels
info = add_kelvin_identification(mesh)   # auto-detects offset + pairs vertices
print(f"Added {info['n_pairs']} pairs, max_dist={info['max_dist']:.2e} m")
```

After this call, `mesh.GetIdentifications()` is populated and
`Periodic(H1(mesh, ...))` correctly slaves DOFs across the two
Kelvin surfaces (verified by `Set(1)|kelvin_int -> kelvin_ext`
boundary-integral ratio = 1.0).

## Critical pre-condition: 1:1 vertex correspondence

The helper requires that `kelvin_int` and `kelvin_ext` surfaces have
a **1:1 vertex correspondence** under the Kelvin rigid translation.
Two paths produce such a mesh:

  - Cubit `copy mesh source vertex ... target vertex` (used by
    `add_kelvin_cubit`)
  - NGSolve OCC `Identify(IdentificationType.PERIODIC)` called
    BEFORE `GenerateMesh()`

A mesh where the two surfaces are meshed INDEPENDENTLY (e.g. two
separately defined spheres glued without `Identify`) is **NOT
supported** -- Netgen produces different vertex distributions on
each surface, the Hungarian pairing gives garbage pairs, and the
resulting Periodic FES does not couple.

To diagnose: after calling `add_kelvin_identification`, inspect
`info["n_unmatched"]`:
  - `0`  -> full 1:1 coverage (good)
  - `>0` -> partial coverage; the open-boundary truncation will be
            incomplete; rebuild the mesh with one of the two 1:1
            paths above
"""

API = """
# API: `radia.add_kelvin_identification(mesh, ...)`

## Signature

```python
def add_kelvin_identification(
    mesh,                                    # ngsolve.Mesh
    kelvin_offset=None,                      # auto-detected if None
    *,
    inner_bnd="kelvin_int",                  # interior sphere BND label
    outer_bnd="kelvin_ext",                  # exterior shell BND label
    idnr=1,                                  # netgen identification number
    point_tolerance=None,                    # auto: max(5%*|offset|, 1mm)
    skip_if_existing=True,                   # no-op if mesh already has idents
) -> dict
```

## Return value

A diagnostic dict with keys:

| Key | Type | Meaning |
|---|---|---|
| `n_pairs` | int | Number of identification pairs registered |
| `n_unmatched` | int | Vertices that exceeded `point_tolerance` |
| `kelvin_offset` | tuple | Offset used (auto-detected or provided) |
| `max_dist` | float | Largest residual distance after translation [m] |
| `point_tolerance` | float | Tolerance used [m] |
| `was_pre_existing` | bool | True if `skip_if_existing` returned early |

## Auto-detection

When `kelvin_offset=None`, the helper calls `detect_kelvin_offset(mesh)`
which returns `mean(kelvin_ext_vertices) - mean(kelvin_int_vertices)`.
This is EXACT for any 1:1 copy-mesh (point-mean is translation-equivariant).

Pre-2026-04-25 bug: `detect_kelvin_offset` once returned the kelvin
*region centroid* (mean of all vertices in the kelvin block), which
gave a translation off by ~14 cm on 1/4-octant geometries.  The
current implementation uses `kelvin_ext` / `kelvin_int` BND vertices
ONLY -- exact regardless of partial-octant geometry.

## Fail-fast on missing labels

If `inner_bnd` or `outer_bnd` is not in `mesh.GetBoundaries()`, the
helper raises `ValueError` listing the missing labels and the
available alternatives.  It does NOT silently no-op (per CLAUDE.md
"No Fallbacks" policy).
"""

WORKFLOW = """
# Recommended workflow

## Step 0: Pre-flight check (optional but recommended)

```python
from radia import has_kelvin_identification, detect_kelvin_offset

if has_kelvin_identification(mesh):
    print("Mesh already has identifications -- nothing to do.")
else:
    offset = detect_kelvin_offset(mesh)
    print(f"Will pair across offset = {offset}")
```

## Step 1: Add identifications

```python
from radia import add_kelvin_identification
info = add_kelvin_identification(mesh)
assert info["n_unmatched"] == 0, (
    f"Partial 1:1 coverage ({info['n_unmatched']} unmatched).  "
    f"The two Kelvin surfaces are NOT a rigid translate of each "
    f"other -- rebuild the mesh via add_kelvin_cubit or OCC Identify.")
```

## Step 2: Verify FES coupling (sub-second sanity check)

```python
from ngsolve import H1, Periodic, GridFunction, Integrate

# FreeDofs count check
fb = H1(mesh, order=1, dirichlet="GND")
fp = Periodic(fb)
slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
assert slaved > 0, f"Periodic slaved 0 DOFs after identification"

# Functional boundary test
gfu = GridFunction(fp)
gfu.vec[:] = 0
gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
ratio = (Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_ext"))
       / Integrate(gfu*gfu, mesh, definedon=mesh.Boundaries("kelvin_int")))
assert abs(ratio - 1.0) < 1e-3, f"Periodic ratio not 1.0: {ratio}"
```

## Step 3: Run the physics solve

`Periodic(H1(...))` or `Periodic(HCurl(...))` is now coupled across
the Kelvin surfaces.  Combine with material scaling
(`nu' = (rho'/R)^2 * nu_0` etc -- see the `kelvin_transformation`
tool) to get the open-boundary FEM result.
"""

SNIPPET = """
# Minimal copy-paste snippet

```python
from ngsolve import Mesh, H1, Periodic
from radia import add_kelvin_identification

mesh = Mesh("model_with_kelvin_int_kelvin_ext.vol")
info = add_kelvin_identification(mesh)
print(f"[kelvin] {info['n_pairs']} pairs, max_dist={info['max_dist']:.2e}, "
      f"unmatched={info['n_unmatched']}")
assert info["n_unmatched"] == 0, "Mesh is not 1:1; rebuild via Cubit/OCC."

# Now Periodic FES couples the Kelvin surfaces.
fp = Periodic(H1(mesh, order=1, dirichlet="GND"))
# ... assemble + solve as usual ...
```
"""

CAVEATS = """
# Caveats

## 1. Helper requires 1:1 mesh (see OVERVIEW section)

The helper is a vertex-pairing operation.  It cannot create a 1:1
correspondence where none exists.  If your mesh has different vertex
counts / positions on the two surfaces, the Hungarian pairing will
either fail (n_unmatched > 0) or produce mathematically meaningless
pairs that do NOT slave DOFs in `Periodic` FES.

## 2. NGSolve Identifications cannot be cleared from Python

`ngmesh.GetIdentifications()` returns a Python LIST that is a COPY
of the internal storage; calling `.clear()` on it has no effect on
the mesh.  If you need to re-pair vertices on a mesh that already
has identifications, use `skip_if_existing=False` -- but be aware
this APPENDS new pairs without removing the existing ones (a small
inefficiency; the FES result is still correct because identifications
are stored as sets of equivalence classes).

## 3. C++ Cubit-export path is preferred when available

When `radia_export netgen` succeeds in writing identifications during
the .vol export (typical 1/2 / 1/4 octant geometry), use those.  This
post-hoc helper is the FALLBACK for cases where the C++ all-or-nothing
matching skipped due to a single tolerance miss -- it accepts the same
mesh but tolerates a few unmatched vertices.

## 4. Tolerance default

Default `point_tolerance` = `max(5% * |kelvin_offset|, 1 mm)`.  For
typical IH / accelerator-magnet geometries with offset ~0.1-1 m, this
is 5-50 mm -- much larger than the sub-mm Netgen vertex jitter on a
true 1:1 mesh, but tight enough to reject the totally-wrong-offset
case.  Pass `point_tolerance` explicitly if your mesh has unusual
scale or you want to verify the residual is sub-mm.

## 5. Identification number `idnr`

Default `idnr=1` matches the Cubit C++ exporter convention.  You only
need to vary this if you plan to add multiple independent identification
sets to the same mesh (rare; e.g. Kelvin + a separate Periodic BC).
"""


def get_post_hoc_documentation(topic: str = "all") -> str:
    """Return documentation for the post-hoc Kelvin identification helper.

    Args:
        topic: which section to return.  See `TOPICS` for valid choices.
    """
    sections = {
        "overview": OVERVIEW,
        "api": API,
        "workflow": WORKFLOW,
        "snippet": SNIPPET,
        "caveats": CAVEATS,
    }
    if topic == "all":
        return "\n\n".join(sections.values())
    if topic not in sections:
        return (f"Unknown topic '{topic}'. Available: "
                f"{sorted(sections.keys())} (or 'all').")
    return sections[topic]


TOPICS = sorted(["overview", "api", "workflow", "snippet", "caveats", "all"])
