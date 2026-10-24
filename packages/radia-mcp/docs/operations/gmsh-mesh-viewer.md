# Human Gmsh post views

The same `python -B -m radia_mcp.mesh_viewer FILE.msh` launcher and
`open_mesh_viewer` MCP tool accept ASCII MSH 4.2. Other formats fail before
native loading. Adjacent `.opt` files are rejected. Gmsh is imported only for
the MSH route; installed payload hashes must match `RADIA_GMSH_PROVENANCE`.
Configuration loading and automatic option/session saving are disabled.
The tool returns launch acceptance, not proof of successful rendering.

Use the declared extra and the canonical
[two-step runtime installation](gmsh-development-runtime.md):

```text
pip install --pre --upgrade -i https://gmsh.info/python-packages-dev gmsh
pip install "radia-mcp[gmsh]"
```

## Explicit VOL post-processing view

```text
python -B -m radia_mcp.mesh_viewer inputs/model.vol --post-output views/model.msh --checked-report checks/model.vol-check.json --source-sha256 CHECKED_VOL_SHA256 --group-scope materials
```

The caller supplies a `radia.checked-vol-handoff.v1` record containing
`source_sha256` and the original passing `check-vol` receipt under `report`.
The producer records that digest when checking the source; a historical report
without a checked-input digest cannot be upgraded by guessing or silently
attaching the current file hash. The receipt includes its versioned label
contract. The same checked digest is supplied as `--source-sha256`. This is a handoff contract,
not a new independent mesh certification. Loaded dimension, element/point counts
and full material/boundary name inventories must match the receipt. The source
hash is checked before and after conversion. Only linear 2D/3D maps are supported;
curved maps need a retained-geometry protocol and fail explicitly here.

The existing production `radia.gmsh_post_export.GmshPostExport.write_mesh` writes
into a new staging directory. `RADIA_MESH_VIEWER_SCRATCH` selects an existing
ASCII job directory outside the source directory; an ASCII output directory is
the default. Japanese output directories require that explicit scratch setting
because the native VOL reader does not reliably handle Unicode filenames.
A hash-identical VOL copy is loaded from staging. The exporter's `.geo`/`.opt`
companions remain there and are never executed, opened or merged. After format
and physical-group/selected-element-count checks, only byte-identical MSH content is copied to a new
destination. Existing outputs, source-directory outputs and adjacent option
scripts are rejected. Staging is retained on both success and failure for
caller-owned evidence and cleanup. Receipt directories are also retained; a
failed conversion creates no success receipt JSON. The operator removes these
explicitly after consumers finish and any required evidence is recovered.
The launcher performs no automatic deletion.
No other MSH writer, reverse conversion or solver initialization is used.

`materials` exports all volume material groups (or 2D material groups).
`boundaries` exports all 3D boundary groups as a surface view. These are distinct
scopes; one output does not promise both. Unsupported dimensions, missing groups
or an installed Radia without the MSH 4.2 exporter fail explicitly.

MCP supplies the same explicit `post_output`, `checked_report`, `source_sha256`
and `group_scope` arguments. The detached viewer delegates conversion to a
separate worker so NGSolve and the Gmsh GUI do not share a process;
preflight performs no writes and probes the selected interpreter's installed
exporter in a separate process. A missing production MSH 4.2 exporter returns
an MCP preflight error before the detached GUI is launched. The server does not wait for conversion
or the GUI. The output is a post-processing view, never an input to a solver.

All-user installation and Windows associations require their separate reviewed
installation decision. This launcher does not modify either.
