Synthetic post-processing-only fixtures for the core's native field-tag
preflight. They are never solver meshes. `generate_native.py` uses the direct
Gmsh API with no MCP dependency: it writes each vector ElementNodeData file,
opens it again, and records `getModelData` in `inputs/native-oracle.json`.
The point fixture has one sample per element, the line fixture two; both have
three components per sample. Only CRLF-to-LF normalization is applied, and both
hashes are retained. The same point fixture reproduced the original core
failure before the component-factor fix.

Regeneration requires the qualified Gmsh development runtime with
`RADIA_GMSH_PROVENANCE` pointing at its official-wheel provenance record.
Normal core tests read the committed files and require neither Gmsh nor MCP.
