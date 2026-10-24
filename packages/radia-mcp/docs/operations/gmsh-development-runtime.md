# Gmsh development runtime

Use the newest qualified official Gmsh **5** wheel. MSH output and input are
ASCII **4.2** only. A checked `.vol` and its versioned label contract remain the
solver mesh authority; `.msh` is exclusively a derived visualization or field
file. Conductor display surfaces/centerlines and streamfunction wires retain
their discretization; their display artifacts use strict MSH 4.2.

The Gmsh dependency requirements live only in the package pyproject extras.
Install Gmsh from its dedicated official development index first, then install
`radia-mcp[gmsh]` from PyPI. Restrict the development index and pre-release
selection to Gmsh; other package dependencies use the normal PyPI index.
README, MCP knowledge and viewer deployment plans link to this canonical recipe:

```text
pip install --pre --upgrade -i https://gmsh.info/python-packages-dev gmsh
pip install "radia-mcp[gmsh]"
```

Gmsh is imported lazily on native post-processing paths. Missing-runtime
errors name the declared extra. MCP tools never install packages or select
an alternative runtime automatically.

For reproducible operations, capture pip's JSON installation report using
`--report` on the first, Gmsh-only command above. Use the exact wheel URL
and SHA256 from that report to retain the selected Gmsh artifact. Download that artifact
without resolving another version, verify its SHA256 against the report,
and bind it to the actually installed runtime:

```text
python -m radia_mcp.gmsh.runtime --wheel wheelhouse/SELECTED_GMSH_WHEEL.whl --record gmsh-provenance.json
```

Archive the wheel with its SHA256, official source URL, fetch timestamp and
index listing before use. Retain the GPL licence with the privately archived
artifact. Set `RADIA_GMSH_PROVENANCE` to the absolute path of the resulting
provenance JSON. Recording checks that the installed Python API and native
library match the wheel; it does not install or qualify numerical semantics.
Native tools fail explicitly if runtime identity, Gmsh major 5, or MSH 4.2 is
missing. Pure Python structural checks require no native Gmsh installation.

Private CI installs the retained artifact with hash verification. Generate
its resolved version and hash from that artifact rather than adding another
policy pin to a script or document. No private wheelhouse path belongs in
public metadata.

When adopting the released Gmsh on PyPI, update the extra requirement and
this canonical install recipe together. Runtime guards check major 5,
provenance and real MSH 4.2 I/O; they do not impose a second development-only
version policy. Released and development builds need the same qualification.

Each native result/session response includes the exact distribution/native
version, wheel SHA256 and verified runtime payload hashes. The version range
does not automatically approve a newly downloaded nightly: before adoption,
rerun per-writer Gmsh load/field roundtrips, opposite-direction API-generated
fixtures, strict format negatives and the measured probe/plugin/rendering
semantics. Preserve the previous evidence as dated provenance without calling
it a measurement of the new build. Record failures; never relax goldens to
adopt a build. A Windows artifact does not establish Linux qualification.

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
