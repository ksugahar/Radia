"""Bundled MSH v4.1 spec reference (excerpted from gmsh.info docs).

**Lab policy (2026-04-19, strict)**:
  - MSH is v4.1 **only**. v2.2 is **completely retired** from the Radia
    repository. No radia-emitted pipeline produces v2.2.
  - The **one exception** is high-order curved meshes, which bypass MSH
    entirely and use NETGEN's native `.vol` format via `radia_export
    netgen "file.vol" order N` (read directly by ngsolve).
  - Incoming v2.2 from external tools is lifted to v4.1 at the
    gmsh_post boundary (`gmsh_post_convert`) before any lab code
    consumes it.

Kept in-source (no external dep) so `gmsh_post_spec` works offline.
"""

SPEC: dict[str, str] = {
	"overview": """# MSH v4.1 — the single mesh interchange format (lab standard)

**Policy (2026-04-19, strict)**: Radia lab uses MSH **v4.1 only**. v2.2
is **completely retired from the radia repo** — no radia-emitted pipeline
produces v2.2. If a v2.2 file enters from an external tool, it is lifted
to v4.1 at the gmsh_post boundary before any lab code touches it.

**Exception**: high-order curved meshes bypass MSH entirely and use
NETGEN's native `.vol` format via `radia_export netgen "file.vol"
order N`; ngsolve reads `.vol` directly.

MSH v4.1 preserves CAD topology: nodes and elements are grouped per
geometric entity (point / curve / surface / volume). ngsolve, radia,
and gmsh GUI consume this format cleanly.

## Section order (ASCII)

```
$MeshFormat                 # required, first
4.1 <file-type> <data-size>
$EndMeshFormat
$PhysicalNames              # strongly recommended
$Entities                   # required in v4+
$PartitionedEntities        # optional (parallel)
$Nodes                      # required
$Elements                   # required
$Periodic                   # optional
$GhostElements              # optional (parallel)
$Parametrizations           # optional (CAD parametrization)
$NodeData / $ElementData / $ElementNodeData   # post-processing views
$InterpolationScheme        # optional
```

Binary vs ASCII: `file-type` = 0 for ASCII, 1 for binary. Default to
ASCII when you plan to author $NodeData by hand; binary for 10M+ element
shipping meshes.
""",

	"meshformat": """# $MeshFormat

```
$MeshFormat
4.1 <file-type> <data-size>
$EndMeshFormat
```

- `4.1` is the only version we write.
- `file-type`: 0 = ASCII, 1 = binary.
- `data-size`: 8 (bytes of `size_t` on 64-bit).

For binary files, the single 4-byte int `1` (endianness marker) follows
the version line before the binary section begins.
""",

	"entities": """# $Entities — geometric entity list (required in v4.1)

Every $Nodes / $Elements block is keyed to a geometric entity.
$Entities declares them up-front with bounding boxes, physical-tag
links, and boundary topology.

```
$Entities
<numPoints> <numCurves> <numSurfaces> <numVolumes>
# points:
<pointTag> <x> <y> <z> <numPhysicalTags> <physTag1> ...
# curves:
<curveTag> <xMin> <yMin> <zMin> <xMax> <yMax> <zMax>
<numPhysicalTags> <physTag1> ...
<numBoundingPoints> <pointTag1> ...
# surfaces, volumes: same "bbox + physical tags + boundary" pattern
$EndEntities
```

**Radia / ngsolve requirement**: this block must be complete. An empty
or incomplete `$Entities` forces ngsolve / radia_ngsolve to flatten the
mesh into one anonymous region — physical groups stop meaning anything.

If the source tool cannot emit `$Entities` natively, gmsh_post rebuilds
it by synthesising one entity per physical group on the fly.
""",

	"nodes": """# $Nodes — entity-grouped

```
$Nodes
<numEntityBlocks> <totalNodes> <minNodeTag> <maxNodeTag>
# per block:
<entityDim> <entityTag> <parametric> <numNodesInBlock>
<nodeTag1>            # tags first
<nodeTag2>
...
<x1> <y1> <z1>        # then coordinates
<x2> <y2> <z2>
...
# if parametric == 1, parametric coords follow
$EndNodes
```

- Tags are global (not per-block) and need not be contiguous.
- `minNodeTag` / `maxNodeTag` bracket the range so readers can
  pre-allocate arrays.
- Entity hierarchy: dim ∈ {0,1,2,3} = point/curve/surface/volume.
""",

	"elements": """# $Elements — entity-grouped

```
$Elements
<numEntityBlocks> <totalElements> <minElemTag> <maxElemTag>
# per block:
<entityDim> <entityTag> <elementType> <numElementsInBlock>
<elemTag> <nodeTag1> <nodeTag2> ...
<elemTag> <nodeTag1> <nodeTag2> ...
...
$EndElements
```

One block = one element type in one entity.  A surface with both quads
and tris becomes two blocks, even though they share the entity.

## Element types (common)

| type | name | nodes |
|------|------|-------|
| 1  | line (2)                 | 2  |
| 2  | tri (3)                  | 3  |
| 3  | quad (4)                 | 4  |
| 4  | tet (4)                  | 4  |
| 5  | hex (8)                  | 8  |
| 6  | prism / wedge (6)        | 6  |
| 7  | pyramid (5)              | 5  |
| 8  | line2 (3)                | 3  |
| 9  | tri2 (6)                 | 6  |
| 10 | quad2 (9)                | 9  |
| 11 | tet2 (10)                | 10 |
| 12 | hex2 (27)                | 27 |
| 15 | vertex (point)           | 1  |

Higher-order types 16+ encode order-3 etc. — see the gmsh header for
the full table.
""",

	"physicalnames": """# $PhysicalNames

```
$PhysicalNames
<numNames>
<dim> <tag> "<name>"
...
$EndPhysicalNames
```

- `dim` ∈ {0..3}.
- `tag` is the integer physical tag referenced from `$Entities`.
- `name` is double-quoted.

Named blocks on the source side (Cubit's `block N name "magnet"`)
become one PhysicalName entry with matching (dim, tag) after v4.1
export.

**Best practice**: name every block you care about downstream. An
unnamed block still meshes but becomes an anonymous region in ngsolve.
""",

	"nodedata": """# $NodeData — scalar / vector / tensor field at nodes

```
$NodeData
<numStringTags>
"<viewName>"
... other string tags
<numRealTags>
<time>
... other real tags
<numIntTags>
<timeStep>
<numComponents>      # 1 = scalar, 3 = vector, 9 = tensor
<numNodes>
<partitionIndex>     # optional
<nodeTag1> <val1_1> <val1_2> ... <val1_numComponents>
<nodeTag2> <val2_1> ...
...
$EndNodeData
```

Typical minimal invocation (scalar at nodes):
```
$NodeData
1
"B_magnitude"
1
0.0
3
0
1
<numNodes>
<nodeTag> <value>
...
$EndNodeData
```

`$ElementData` has the same structure but keyed by element tag.
`$ElementNodeData` gives per-node-per-element values (discontinuous
fields, DG bases).
""",

	"write_node_data": """# Appending $NodeData to an existing v4.1 mesh

Typical flow:

1. Mesh is already written as v4.1 (Cubit → gmsh → file.msh).
2. You compute a field per node (from ngsolve / radia / FEM solver).
3. Write the new v4.1 file that is the mesh + an additional $NodeData
   block at the end.

`gmsh_post_write_node_data(msh_path, values, out_path, view_name)`
implements this: it reads the input mesh via the gmsh API
(`gmsh.open`), calls `gmsh.view.add` + `gmsh.view.addModelData`, and
writes the result as v4.1 ASCII. `values` is `{nodeTag: value}` or
`{nodeTag: [vx, vy, vz]}` for vectors.

The resulting `.msh` is still a valid v4.1 mesh file; opening it in
the Gmsh GUI shows the view under the mesh.
""",

	"compat_checks": """# v4.1 compliance / lab-compat checklist

`gmsh_post_validate` inspects a .msh against this list:

1. **$MeshFormat is 4.1** — reject anything else. v2.2 is deprecated;
   the converter lifts it but validation flags lingering v2.2.
2. **$Entities present and non-empty** — required. If absent, ngsolve
   will flatten the topology.
3. **$PhysicalNames covers every $Elements dim** — otherwise named
   regions disappear downstream.
4. **No duplicate nodes** — Cubit after `imprint all` can emit nodes
   at identical coords on shared faces. Validator runs an O(N log N)
   spatial check; fix with `gmsh.model.mesh.removeDuplicateNodes()`.
5. **Element type set matches downstream solver** — ngsolve wants
   linear or order-2 types; radia_ngsolve wants matching order across
   all $Elements blocks.
6. **$MeshFormat is the first section** — a leading BOM / blank line
   breaks most parsers.
""",

	"cubit_to_v41": """# Lab export workflow (v4.1 everywhere, .vol for HO curved)

Two mesh-export paths out of Cubit, zero v2.2:

## Path A — linear mesh → MSH v4.1

If Cubit's `export mesh` still emits v2.2 on your build, lift at the
gmsh_post boundary:

```python
# Cubit side (linear hex/tet mesh):
cubit.cmd('export mesh "/tmp/raw.msh" overwrite')

# Lift to v4.1 once — any lab consumer assumes v4.1:
import radia_mcp.gmsh_post.server as g
g.gmsh_post_convert('/tmp/raw.msh', '/tmp/final_v41.msh')
```

`cubit_mesh_auto` (radia-mcp ≥ 0.14.0) auto-invokes the conversion on
its output, so the GUI-replay artefact is v4.1 regardless of the Cubit
build's default.

## Path B — high-order curved mesh → NETGEN `.vol`

MSH is **not** used for high-order curved meshes in this lab. Use
radia's NETGEN exporter, which writes the native format ngsolve
consumes directly:

```python
cubit.cmd('radia_export netgen "/tmp/mesh.vol" order 2 overwrite')
```

Then in ngsolve:
```python
from ngsolve import Mesh
mesh = Mesh('/tmp/mesh.vol')
```

No MSH, no conversion step, no round-trip loss. HO curved elements
(radia_export's SetGeomInfo-annotated curving) survive the `.vol`
boundary; MSH would mangle them.
""",
}


SPEC["view_data_cookbook"] = """# Gmsh view data cookbook — which `$*Data` block to use

Post-processing data in MSH v4.1 lives in three similar-but-distinct
block types. Picking the right one is the #1 stumbling point.

## Decision tree

```
Is the value defined at nodes?            → $NodeData
Is it a single value per element?          → $ElementData
Is it a different value at each node of each element
  (e.g., discontinuous DG, per-face data)? → $ElementNodeData
```

## Shared structure

```
$<DataType>
<numStringTags>
"<viewName>"
... other string tags ...
<numRealTags>
<time>
... other real tags ...
<numIntTags>
<timeStep>          # for transient data, increments per frame
<numComponents>     # 1 scalar / 3 vector / 9 tensor
<numEntries>
<partitionIndex>    # optional (parallel meshes)
<tag> <v1> <v2> ... <v_numComponents>
<tag> <v1> <v2> ...
...
$End<DataType>
```

## Minimal scalar $NodeData

```python
import gmsh
gmsh.initialize()
gmsh.open("mesh.msh")
vtag = gmsh.view.add("B_magnitude")
# node_tags and values are arrays of the same length
gmsh.view.addModelData(vtag, step=0, modelName=gmsh.model.getCurrent(),
                       dataType="NodeData",
                       tags=node_tags, data=[[v] for v in values],
                       time=0.0, numComponents=1)
gmsh.write("mesh_with_B.msh")
gmsh.finalize()
```

The radia-mcp tool `gmsh_post_write_node_data` wraps this exactly.

## Vector field: magnetic B at elements

```python
vtag = gmsh.view.add("B_vector")
gmsh.view.addModelData(vtag, 0, model, "ElementData",
                       elem_tags, [[Bx, By, Bz] for ...],
                       time=0.0, numComponents=3)
```

## Transient (multi-timestep) data

Loop over time steps, calling `addModelData` with incrementing `step`:

```python
for k, t_now in enumerate(times):
    gmsh.view.addModelData(vtag, k, model, "NodeData",
                           tags, data_at_t[k],
                           time=t_now, numComponents=1)
```

The GUI `Tools › Animate` slider then scrubs through steps.

## View options worth knowing

```python
gmsh.view.option.setNumber(vtag, "Visible", 1)         # show
gmsh.view.option.setNumber(vtag, "IntervalsType", 3)    # 1 iso, 2 filled, 3 num
gmsh.view.option.setNumber(vtag, "NbIso", 15)           # isoline count
gmsh.view.option.setNumber(vtag, "RangeType", 2)        # 1 default, 2 custom
gmsh.view.option.setNumber(vtag, "CustomMin", 0.0)
gmsh.view.option.setNumber(vtag, "CustomMax", 1.5)
gmsh.view.option.setString(vtag, "ColormapName", "Viridis")
```

## InterpolationScheme (high-order / DG)

For HO meshes the $ElementNodeData values need a matching
`$InterpolationScheme` so Gmsh knows how to evaluate between nodes:

```
$InterpolationScheme
"myScheme"
<numElementTopologies>
<elementType>
<numBases>
<order>
<basis matrix: n rows of numMonomials floats>
... and a coefficient matrix ...
$EndInterpolationScheme
```

Most lab users don't write this by hand — `radia_export netgen ...`
+ ngsolve HO assembly skip MSH for HO curved fields entirely and go
straight to NETGEN `.vol`.

## Pitfalls

1. **Wrong key** — passing `node_tags` to an `ElementData` view gives
   a nonsense (but not always erroring) view.
2. **Length mismatch** — `tags` and `data` must be the same length.
   Gmsh silently truncates to the shorter; symptoms: missing colors
   on half the mesh.
3. **Units / scaling** — the raw numbers go into the view as-is. Set
   `CustomMin/Max` before sharing, or use `Tools › Options › View` UI.
4. **Appending to an existing file** — `gmsh.write(path)` rewrites
   the whole mesh + loaded views. `gmsh.view.write(tag, path,
   append=True)` appends ONLY the view but emits a fresh
   `$MeshFormat` header → breaks v4.1 spec.  radia-mcp's
   `gmsh_post_write_node_data` strips the redundant header after
   append (see `_append_data_block`).
"""


SPEC["physical_groups_cookbook"] = """# Physical groups + entities cookbook

Physical groups are the **named-region** mechanism in Gmsh and the
only reason downstream FEM solvers can distinguish "iron" from "coil"
from "air". They live at the intersection of geometry (`$Entities`)
and semantics (`$PhysicalNames`).

## Mental model

```
$Entities           Physical group
   dim, tag    ←→   dim, tag  ←→  name
(geometric)        (PhysicalGroup)     (human-readable)
```

A **physical group** is a tuple `(dim, tag)` with a name. It groups
one or more geometric entities of the same dimension. Element blocks
inherit the group of the entity they live on.

## Creating physical groups from a script

```python
import gmsh
gmsh.initialize()
# ... build geometry with gmsh.model.occ / .geo ...
gmsh.model.occ.synchronize()

# Mark volume tags 1-3 as "iron", 4 as "coil", 5 as "air"
iron_tag = gmsh.model.addPhysicalGroup(dim=3, tags=[1, 2, 3])
gmsh.model.setPhysicalName(dim=3, tag=iron_tag, name="iron")
coil_tag = gmsh.model.addPhysicalGroup(3, [4])
gmsh.model.setPhysicalName(3, coil_tag, "coil")
air_tag = gmsh.model.addPhysicalGroup(3, [5])
gmsh.model.setPhysicalName(3, air_tag, "air")

# Surface groups for boundary conditions
gmsh.model.addPhysicalGroup(2, [10, 11], name="inlet")
gmsh.model.addPhysicalGroup(2, [20], name="wall")

gmsh.model.mesh.generate(3)
gmsh.write("with_groups.msh")
```

## Reading physical groups from a .msh

```python
for dim, tag in gmsh.model.getPhysicalGroups():
    name = gmsh.model.getPhysicalName(dim, tag)
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    print(f"dim={dim} tag={tag} name={name!r} entities={entities}")
```

## Conventions that downstream tools rely on

| consumer | group dimension expected | notes |
|---|---|---|
| ngsolve `Mesh("x.msh")` | 3 for materials, 2 for boundaries | `.Materials()` / `.Boundaries()` return names |
| radia_ngsolve HO | 3 for materials; 2 for Neumann; 1 for Dirichlet edges | matching names drive ngbem assembly |
| Abaqus export | 2 for surfaces (block / contact) | `export abaqus "..." overwrite` reads names for `*Nset` |
| Nastran export | 3 for `*Material`, 2 for `*Pset` | block names feed `PSOLID`/`PSHELL` |
| Cubit | 3 for block, 2 for sideset, 1 for nodeset | `block N add volume ...; block N name "..."` |

## Tag collisions across dimensions

Gmsh allows tag `1` to exist in dim 2 and dim 3 simultaneously —
identity = (dim, tag), not tag alone.  Some older downstream tools
(legacy Nastran, radia < 4.6) can't cope.  Keep tags globally unique
when you know you'll ship to such a tool:

```python
# Reserve ranges per dimension: volumes 100-199, surfaces 200-299, curves 300-399
gmsh.model.addPhysicalGroup(3, [1, 2], tag=100, name="iron")
gmsh.model.addPhysicalGroup(2, [10],   tag=200, name="inlet")
gmsh.model.addPhysicalGroup(1, [20],   tag=300, name="fixed_edge")
```

## Cubit → Gmsh v4.1 round-trip

Named Cubit blocks emit as `$PhysicalNames` entries:

| Cubit side | MSH v4.1 |
|---|---|
| `block 1 add volume 1 2`<br>`block 1 name "iron"` | `$PhysicalNames: 3 1 "iron"` + entity link |
| `sideset 10 add surface 20` | goes to `$Elements` but NOT `$PhysicalNames` unless you also `block N add tri`... |

After `gmsh_post_convert` lifts v2.2 to v4.1, re-verify with
`gmsh_post_inspect` — if a block dim doesn't show in
`physical_groups`, the source didn't name it.

## Extract one physical group into its own mesh

`gmsh_post_extract_physical(msh, "iron", "iron_only.msh")` (new in
0.20.0) reads the source, identifies elements in the named group,
and writes them as a fresh v4.1 file with their own entities /
physicals preserved.  Useful for per-region mesh debug or handing
one sub-mesh to a single-material solver.
"""


SPEC["api_reference"] = None  # filled in below from the generated file


def get_spec(topic: str = "overview") -> str:
	t = (topic or "overview").strip().lower()
	if t in ("all", "full"):
		return "\n\n---\n\n".join(f"{v}" for v in SPEC.values() if v)
	if t not in SPEC or SPEC[t] is None:
		return (f"# Unknown topic {topic!r}\n\nAvailable topics: "
		        + ", ".join(k for k, v in SPEC.items() if v))
	return SPEC[t]
