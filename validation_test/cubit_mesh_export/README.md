# Cubit Mesh Export Validation Corpus

## Standalone wheel acceptance

`validate_standalone_wheel.py` checks the exact installed wheel without Radia
or radia-mcp in a fresh venv. It verifies installed package bytes and native
payload hashes, runs `standalone_sphere.jou` through Cubit's console with
`-batch -nographics -noinitfile`, and requires order-2 export, strict labels,
CAD measures, NGSolve reload and positive mapped Jacobians. No GUI or system
plugin replacement is involved. This is a licensed-host check, not normal CI.

```powershell
python -m venv C:/temp/cubit-standalone-venv
C:/temp/cubit-standalone-venv/Scripts/python.exe -m pip install <candidate.whl>
C:/temp/cubit-standalone-venv/Scripts/python.exe -I validation_test/cubit_mesh_export/validate_standalone_wheel.py --wheel <candidate.whl> --cubit-bin "C:/Program Files/Coreform Cubit 2025.12/bin" --output C:/temp/cubit-standalone-result
```

The output directory must be new. Retain `result.json`, `vol-check.json`,
the journals, `.vol`/sidecar and `cubit.log` together, including failed runs.
The result records the wheel hash, actual import path, dependency versions,
absence of Radia/MCP, native hashes and exact Cubit command. A passing result
certifies this standalone sphere route; it does not certify the optional
Radia toolbar or replace broader geometry regression evidence.

On a host with an existing exporter, add `--use-installed-plugin`. This
requires the deployed `.ccm` and Netgen DLLs to match the candidate environment
byte-for-byte and explicitly selects that same plugin directory. Selecting a
second directory would register commands twice. An absent or different deployed
binary fails acceptance rather than silently testing a different exporter.

Runnable validation scripts promoted from the former Cubit mesh export examples
tree live here.

Presentation/demo material is kept under `docs/cubit_mesh_export`; these files
are the regression-oriented checks for Netgen `.vol` boundary inventory,
surface quality, pressure/traction resultants, and FEM/BEM topology handling.
`cubit_mesh_export_showcase_results.json` aggregates that checked evidence for
the executed public showcase notebook without placing a duplicate result JSON
under `docs/`.
`p_convergence_demo_results.json` records the Cubit-to-`.vol` p-convergence
run consumed by the executed docs notebook of the same name.

## Geometric curving refit

`geometric_refit_benchmark.py` exports the same TET and HEX unit spheres at
orders 1-5 twice through one candidate `.ccm`, with Netgen's coefficients
(`CUBIT_MESH_EXPORT_GEOMETRIC_REFIT=0`) and with the geometric refit, then
compares boundary radial error, the electrostatic field error of
`paper_sphere_benchmark.py` and check-vol. On a host with an installed
exporter Cubit also loads the deployed plugin, so the commands register
twice; the script therefore requires the candidate's own refit log line
(`geometric refit - ...` or `geometric refit disabled`) before it accepts
the exported files.

```powershell
python validation_test/cubit_mesh_export/geometric_refit_benchmark.py --plugin-dir <dir with candidate cubit_mesh_export.ccm>
```

Pass `--plugin-dir` as an absolute path: Cubit rejects a relative one
("Invalid Plugin Directory") and silently loads the installed plugin.
It writes `geometric_refit_benchmark_results.json` next to the script.
