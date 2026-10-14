# Coreform Cubit and High-Order HEX Validation

The Cubit export tests are opt-in because importing `cubit` can start the
Coreform engine and license checkout during pytest collection.

## Public Contract

Use Cubit's APREPRO export command for Netgen/NGSolve meshes:

```python
cubit.cmd('export netgen "model.vol" order 2 overwrite')
```

The old top-level Python helpers such as `extract_curved_mesh`,
`extract_mesh_data`, and `export_netgen()` are retired from the public API.
The `cubit_mesh_curver` `.pyd` module is still present as a low-level plugin
implementation detail; tests should import it only for module-shape smoke
coverage.

## Running Tests

Set `RADIA_RUN_CUBIT_TESTS=1` explicitly. The helper auto-detects Coreform
Cubit 2025.12+ through the same `install_panels.find_cubit_bin()` path used by
the installer. `CUBIT_PATH` can override discovery when needed.

```powershell
$env:RADIA_RUN_CUBIT_TESTS = "1"
$env:CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.12/bin"
python -m pytest validation_test/cubit -q
```

Coreform Cubit 2025.12 embeds Python 3.10. Real `import cubit` tests skip when
run under an incompatible system Python (for example Python 3.12). To execute
those tests instead of skipping them, run pytest from Cubit's bundled Python
environment:

```powershell
& "C:/Program Files/Coreform Cubit 2025.12/bin/python3/python.exe" -m pytest validation_test/cubit -q
```

Most CI and normal repository test runs intentionally skip this directory.

## NGSolve-Only High-Order HEX Evidence

`test_hex_highorder_fem.py` validates p-refinement, h-refinement, and a smooth
variable-coefficient weak form on genuine 3D HEX elements without requiring a
Cubit license. Regenerate its committed numerical evidence with:

```powershell
python -m pytest validation_test/cubit/test_hex_highorder_fem.py -q
python validation_test/cubit/generate_hex_highorder_fem_results.py
```

The generator owns `hex_highorder_fem_results.json`; do not edit its numerical
values by hand.

## Complex `.vol` volume-accuracy learning data

`vol_accuracy_dataset/` is the durable complex-geometry validation lane. It
builds Boolean-cut flanges and busbars, a stepped spacer, a trimmed half torus,
a mixed tet/wedge boundary-layer solid, and a periodic-seam loft regression;
it also keeps coarse/refined 355-degree sweep twins to distinguish order-3
under-resolution from an exporter defect. It exports each at orders 1 through
3, then records Cubit ACIS volume versus NGSolve-integrated `.vol` volume
together with Jacobian quality. See its README and run
`python validation_test/cubit/vol_accuracy_dataset/generate.py` after closing
interactive Cubit.

## Test Notes

- NGSolve must be imported before Cubit in tests that use both packages, to
  avoid Windows DLL load conflicts.
- Cubit is initialized in batch/no-graphics mode with the 2025.12 plugin
  directory passed through `-commandplugindir`.
- Temporary `.vol` and `.msh` files are written under the system temp directory.
- Tests should prefer `validation_test/cubit/cubit_202512_helpers.py` over hand-written
  Cubit path setup.
