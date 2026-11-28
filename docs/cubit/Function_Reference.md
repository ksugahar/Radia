# Function Reference

Reference documentation for the `radia_mcp.cubit` APREPRO export command and the `check-vol` checker.

## APREPRO Commands (Recommended)

Native Cubit command registered by the `cubit_to_vol_export.dll` command plugin. No Python import needed.
The command is available in journal files (.jou) and the Cubit command line.

### Mesh Export Commands

| Command | Format | Orders | HO Method | Doc |
|---------|--------|--------|-----------|-----|
| `export vol "f.vol" order N` | Netgen .vol (+ .vol.json) | 1-5 | NetgenCurver + ACIS | [export_NetgenMesh](export_NetgenMesh.md) |

### Build & Installation

The plugin `cubit_to_vol_export.dll` is built from
`packages/radia-mcp/native/cubit_to_vol_export` against the fork Netgen:

```bash
python packages/radia-mcp/native/cubit_to_vol_export/build.py \
  --cubit-sdk "<Cubit>/cmake" --netgen-dir "<venv>/Lib/site-packages/netgen/cmake" --out <plugin_dir>
# The packaged Export Vol menu (radia_mcp.cubit.gui) is the sole GUI entry point.
```

Installation: `python -m radia_mcp.cubit.install --plugin-dir <plugin_dir>`

---

## Command Details

### export vol

```
export vol "filename.vol" [order <1-5>] [overwrite]
```

Exports mesh with high-order curving (CallbackGeometry + ACIS projection).
Produces companion JSON (.vol.json) with CAD reference values for consistency checks.

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 2 | Curve order (1=linear, 2-5=high-order via NetgenCurver) |
| overwrite | off | Overwrite existing file |

## Python API

### In-Memory Curving

Curved meshes come from the Cubit command:

```python
cubit.cmd('export vol "model.vol" order 3 overwrite')
```

The plugin curves the mesh in-process after Cubit has supplied linear
mesh arrays and geometry callbacks.

### check_consistency (radia.vol_check)

```bash
check-vol model.vol                         # Full .vol gate; sidecar optional
check-vol model.vol --strict-labels         # Enforce canonical label names
check-vol model.vol --contract labels.json  # Application/mode label contract
check-vol model.vol --report-json run/vol_check.json
```

```python
from radia.vol_check import check_consistency, check_label_contract
```

The checker runs after export and before solver/Simulink initialization. It
always validates mesh loading, label relations, and the curved NGSolve mapping.
When `model.vol.json` exists it is auto-discovered for CAD volume/area/length
and mesh metadata comparison; `--json` makes an explicit sidecar mandatory.
Material constants remain part of the application DesignSpec/configuration and
are not inferred from `.vol` label strings.

---

## GUI Menu Structure

```
Menu bar: ... Export Vol  Help

Export Vol (PySide6):
  Export .vol...
```

- **Export Vol**: PySide6 dialog with settings persistence (`%APPDATA%/radia-mcp/export_settings.json`)
- Solver orchestration is intentionally outside the Cubit plugin. Downstream
  applications consume the exported mesh through their own process boundary.
- **Mesh p-convergence demo**: a saved notebook under `docs/`; it is not an engineering design panel or Cubit menu action.

---

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| "Interrupt Detected" on AddPoint | ABI mismatch (full Netgen DLL) | Rebuild with compact_netgen |
| cp932 UnicodeDecodeError | Non-ASCII in .py | Use ASCII only + encoding='utf-8' |
