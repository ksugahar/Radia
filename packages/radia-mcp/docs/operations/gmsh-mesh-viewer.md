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
and physical-group/selected-element-count checks, byte-identical MSH content is
copied to a temporary name beside the destination, hash-verified, and published
by an exclusive hard link. The output filesystem must support hard links;
otherwise publication fails explicitly. The temporary copy is removed on success
or failure. Existing outputs, source-directory outputs and adjacent option
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
and `group_scope` arguments. MCP awaits a separate conversion worker without
blocking its event loop. Only successful conversion, a bound receipt and MSH
preflight can reach GUI launch; every conversion gate failure returns to MCP.
NGSolve and the Gmsh GUI never share a process. The `post_view` response contains
the verified output hash, source hash, group scope, physical groups, receipt and
staging paths. The GUI receives the checked MSH alone, not another conversion
request. A GUI spawn failure reports the retained verified artifact. Launch
acceptance remains distinct from rendering success. The output is a
post-processing view, never an input to a solver.

`RADIA_MESH_VIEWER_TIMEOUT_SECONDS` is an operator setting: default 300 seconds,
finite range 30..900, excluding the human GUI lifetime. Preflight and conversion
share one monotonic deadline and a cleanup reserve. Windows workers are started
suspended, assigned to their explicit kill-on-close job, and resumed; timeout
and cancellation reap the owned descendants, including a venv redirector.
UTF-8 is explicit on outer and inner children. Native non-UTF-8 diagnostics use
replacement decoding, never strict reader-thread decoding. Diagnostic capture
and receipt metadata are capped at 1 MiB each; exceeding a cap fails explicitly.

Native MSH eligibility streams sections and field blocks, keeping a bounded
native-tag bitmap rather than all data tags. It preserves the header/section/
field checks and the native field-tag maximum 2147483647; it is not the full
mesh/reference validator. Detailed `inspect_msh` reporting still materializes
inventories, including the current converter's inspection step. Consequently
this change does not certify arbitrary >2 GiB inspection, VOL conversion or GUI
rendering. A separate scale proof is required; no numerical gates are weakened.

All-user installation and Windows associations require their separate reviewed
installation decision. This launcher does not modify either.
