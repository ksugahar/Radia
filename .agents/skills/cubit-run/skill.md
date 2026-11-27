---
name: cubit-run
description: Run Cubit .jou files and export curved .vol meshes via the radia-mcp `export netgen` command plugin
allowed-tools: Bash(python *), Bash(*coreform_cubit*), Bash(ssh *), Bash(cp *), Bash(sed *), Bash(cat *), Read, Glob, Write
---

# Cubit .jou Runner & Tester

> **CURRENT POLICY**: the only export command is `export netgen` from
> `radia_cubit_export.dll` (`radia_mcp.cubit`, deployed by
> `python -m radia_mcp.cubit.install`). Prefer `C:\temp` for temporary
> journals and outputs.

Run Cubit journal (.jou) files and test the `.vol` export orders.

## Arguments

| Argument | Action |
|----------|--------|
| `<path/to/file.jou>` | Run a single .jou file with full export test suite |
| `test` or `tests` | Run all tests in `tests/cubit/export_mixed_test/` |
| `test 100` | Run tests on 100号機 via SSH |

## Full Export Test Suite (per .jou)

For EVERY .jou file tested, run ALL of these export formats:

```bash
CUBIT="C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe"

# Generate test .jou from the user's .jou by appending export commands
# The user's .jou should contain geometry + mesh + blocks only (no export)
# Append these export lines with output to a temp directory:

export netgen "{DIR}/test_o1.vol" order 1 overwrite
export netgen "{DIR}/test_o2.vol" order 2 overwrite
export netgen "{DIR}/test_o3.vol" order 3 overwrite
export netgen "{DIR}/test_o4.vol" order 4 overwrite
export netgen "{DIR}/test_o5.vol" order 5 overwrite
```

### Verification checklist (per file)

1. **File exists** and size > minimum threshold
2. **check-vol**: every .vol passes `check-vol --strict-labels`
3. **p-convergence monotonic**: .vol file sizes should increase with order (order N+1 > order N)
4. **No "Interrupt Detected"** in Cubit output

### Running on 100号機

```bash
# Copy .jou to 100号機 via SMB
cp test.jou "//192.168.11.100/c$/temp/test.jou"

# Run via SSH
ssh 100 '& "C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe" -batch -nographics -nojournal C:\temp\test.jou 2>&1 | Select-String "Exported|succeed|FAIL|Interrupt"'

# Verify: all lines show "Exported" or "succeeded", no "Interrupt"
```

## Existing test suite

```bash
# Run all registered tests (local)
python tests/cubit/export_mixed_test/run_test.py

# Tests included:
# 1. Mixed hex+tet+pyramid: order 1+2, gmsh+nastran+vol
# 2. Hex cylinder: order 2, gmsh+nastran (HEX20 = 20 nodes check)
# 3. Sphere p-convergence: order 1-5, vol
```

## Export commands (APREPRO)

| Command | Format | Orders | Notes |
|---------|--------|--------|-------|
| `export gmsh "f.msh" order N overwrite` | MSH 4.2 post-processing | 1-4 | Solver input is checked `.vol`. |
| `export radia_nastran "f.bdf" order N` | Nastran BDF | 1-2 | NOT `export nastran` (built-in conflict) |
| `export vtk "f.vtk" order N` | VTK Legacy | 1-2 | |
| `export meg "f.meg"` | MEG/ELF | 1 only | |
| `export netgen "f.vol" order N` | Netgen Vol | 1-5 | p-convergence test |
| `coil "script.py"` | CoilBuilder STEP | -- | subprocess Python 3.12 |

## Troubleshooting

- **`Interrupt Detected` on AddPoint**: ABI mismatch. Rebuild ccm with compact_netgen (not full Netgen DLL).
- **HEX20 has 8 nodes**: `edge_ho_nodes_` not registered for volume element edges. Fixed 2026-04-05.
- **Non-zero exit code (0xC0000005)**: Cubit cleanup crash. Check output files, not exit code.
- **cp932 UnicodeDecodeError**: Non-ASCII in .jou filename or .py panel. Use ASCII only.
- **`topLevelWidgets` error**: Batch mode (-nographics). Harmless, panels skip registration.
