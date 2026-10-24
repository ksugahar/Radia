# Gmsh development runtime

Use the newest official Gmsh **5 development** wheel. MSH output and input are
ASCII **4.2** only. A checked `.vol` and its versioned label contract remain the
solver mesh authority; `.msh` is exclusively a derived visualization or field
file. Conductor display surfaces/centerlines and streamfunction wires retain
their discretization; their display artifacts use strict MSH 4.2.

The `gmsh` and `shape-regeneration` extras allow `gmsh>=5.0.0.dev0`.
PyPI may only have the stable Gmsh 4 release, so obtain the development wheel
from its public official development index. Install Gmsh alone from that index,
then install the package and its other dependencies from PyPI:

```text
python -m pip install --pre -i https://gmsh.info/python-packages-dev gmsh
python -m pip install "radia-mcp[gmsh]"
```

For reproducible private operations, download the selected wheel separately:

```text
python -m pip download --pre --no-deps --only-binary=:all: -i https://gmsh.info/python-packages-dev gmsh -d wheelhouse
python -m pip hash wheelhouse/SELECTED_GMSH_WHEEL.whl
python -m pip install --no-index --no-deps wheelhouse/SELECTED_GMSH_WHEEL.whl
python -m radia_mcp.gmsh.runtime --wheel wheelhouse/SELECTED_GMSH_WHEEL.whl --record gmsh-provenance.json
```

Archive the wheel with its SHA256, official source URL, fetch timestamp and
index listing before use. Retain the GPL licence with the privately archived
artifact. Set `RADIA_GMSH_PROVENANCE` to the absolute path of the resulting
provenance JSON. Recording checks that the installed Python API and native
library match the wheel; it does not install or qualify numerical semantics.
Native tools fail explicitly if runtime identity, Gmsh major 5, or MSH 4.2 is
missing. Pure Python structural checks require no native Gmsh installation.

For a hash-locked private CI install, create a requirements file with the exact
resolved development version and its actual archive digest:

```text
gmsh==RESOLVED_VERSION --hash=sha256:VERIFIED_DIGEST
```

Then install with `--require-hashes --no-deps --no-index --find-links wheelhouse
-r gmsh-requirements.txt`. If resolving dependencies in the same command, all
requirements and their transitive dependencies must also be pinned and hashed.
No private wheelhouse location belongs in public package metadata.

Each native result/session response includes the exact distribution/native
version, wheel SHA256 and verified runtime payload hashes. The version range
does not automatically approve a newly downloaded nightly: before adoption,
rerun per-writer Gmsh load/field roundtrips, opposite-direction API-generated
fixtures, strict format negatives and the measured probe/plugin/rendering
semantics. Preserve the previous evidence as dated provenance without calling
it a measurement of the new build. Record failures; never relax goldens to
adopt a nightly. A Windows artifact does not establish Linux qualification.

`gmsh.write` writes the mesh, while `gmsh.view.write` writes views and can also
write a mesh. When appending multiple model-based views to one file, set
`PostProcessing.SaveMesh=0` after the initial mesh write. List-based `.pos`
views remain explicit view artifacts; a plugin writing an unsupported MSH format is rejected.
Changing only the version header does not produce 4.2: field/ghost records are
blocked, and non-point entity records include embedded entity lists.

The pure Python field/structural readers support ordinary unpartitioned element
meshes, embedded entities and NodeData/ElementData/ElementNodeData blocks.
PartitionedEntities, GhostElements, Polytopes, Edges and Faces sections currently
raise explicit unsupported-section errors; they are never parsed as older layouts.
