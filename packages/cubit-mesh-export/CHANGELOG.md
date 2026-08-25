# Changelog

All notable changes to `cubit-mesh-export` — the high-order curved
mesh export package for Coreform Cubit (Netgen / GMSH / Nastran /
VTK / MEG / FEMEEM writers + Python bindings for consistency checks).

## Unreleased

## 0.14.10 - Netgen 6.2.2606 ABI migration

Released 2026-08-25.

- Pinned NGSolve and Netgen to 6.2.2606 so the bundled high-order curver,
  deployed Netgen DLLs, `check-vol`, and solver runtime use one C++ ABI.
- Rebuilt the Cubit backend and high-order curver and re-ran the complex-shape
  volume and curved-Jacobian acceptance gates.
- Evaluated curved Jacobians through NGSolve's vectorized mapping route and
  kept periodic Cubit edge interpolation on the local unwrapped branch.

## 0.14.9 - Sculpt/Exodus sideset preservation

Released 2026-08-19.

- Preserve direct/free triangle and quadrilateral faces from Sculpt or
  imported-Exodus sidesets in Netgen exports.  Equal geometry-owned faces are
  deduplicated, free boundaries recover adjacent material domains, and a
  sideset spanning several domain pairs is split into domain-aware descriptors.
  Sideset names/areas survive in `.vol` and companion metadata without area
  double-counting.  `check-vol` now rejects missing or invalid boundary-domain
  ownership before solver use.
- Distinguish free-mesh block labels (`mesh_only_materials`) from materials with
  a real CAD-volume reference, and omit stale CAD curve lengths when no BBND
  mesh segment was exported.  This prevents false zero-volume and edge-length
  comparisons for imported STL/Sculpt models.
- Preserve Skin-generated free material interfaces as one non-duplicated
  `DomainIn -> DomainOut` surface set.  Free-sideset sidecar areas now measure
  only the faces actually exported under each descriptor, and the inventory
  reports exterior/interface counts while rejecting duplicate connectivity.

## 0.14.8 - Reproducible plugin release build

Released 2026-08-06.

- Rebuilt the bundled `cubit_mesh_export.ccm` from the unchanged canonical
  source with the current MSVC release toolchain for the four-machine release
  candidate; public APIs and mesh-format contracts remain unchanged.

## 0.14.7 - Cubit plugin release refresh

Released 2026-08-05.

- Rebuilt the bundled `cubit_mesh_export.ccm` with the current release
  toolchain; the public Python API and mesh-format contracts are unchanged.

## 0.14.6 - Compressed VOL label normalization

Released 2026-08-05.

- Normalize transport-only label whitespace reported by NGSolve so strict
  `check-vol` contracts treat Windows CRLF `.vol.gz` files identically to their
  uncompressed `.vol` sources, including conductor-face adjacency checks.

## 0.14.5 - MEG nonlinear-magnet labels

Released 2026-08-04.

- Corrected the Cubit MEG exporter help contract so `MWL` denotes a nonlinear
  magnet with a fixed local axis and `MWV` denotes a nonlinear magnet with
  direction vectors; the plugin help, Radia Cubit menu, and public docs now
  agree.

## 0.14.4 - Production `.vol` preflight

Released 2026-07-22.

- Finish the standalone `check-vol` gate: a `.vol` now works without a CAD
  sidecar, while an available `.vol.json` is auto-discovered for volume, area,
  total edge-length, element-count, point-count, and curve-order comparisons.
- Add versioned application label contracts, strict canonical-name checks,
  generated-label and case-collision detection, Kelvin/source/symmetry relation
  checks, JSON reports, and stable CLI exit codes.
- Align the CLI and Python API aliases documented by the package (`--quality`,
  `--tet-only`, `conductors=`, and `tet_only=`).
- Treat a consistently negative NGSolve element orientation as valid and use
  `abs(det(J))` for scaled quality; fail only singular maps or sign changes
  within an element. This removes false inversion reports on Cubit `.vol`
  tetrahedra while retaining folded high-order-map detection.
- Sample affine order-1 mappings once per element (their Jacobian is constant)
  while retaining the high-order integration rule for curved elements.

## 0.14.3 - Simulink application handoff documentation

Released 2026-07-21.

- Updated the Cubit export handoff to point users to Radia's production
  Simulink application blocks, with the IH notebook retained only for its
  temporary comparison period.

## 0.14.2 - Curved-mesh and conductor-face quality gates

Released 2026-07-17.

- Add sampled high-order Jacobian checks, required label checks, and
  conductor/SIBC face-adjacency classification to `check-vol` and the Python
  consistency API.
- Require Coreform Cubit 2025.12+ in both `cubit-plugin-install` and
  Radia panel registration; older 2025.3/2025.6 installs are no longer
  selected accidentally.
- Move the generated Radia toolbar startup shim out of the Python package
  tree and into `%ProgramData%/Radia/Cubit/` for `--all-users` installs
  (or `%LOCALAPPDATA%/Radia/Cubit/` for current-user installs), so first
  install no longer rewrites tracked/editable `startup.py`.
- Extend `cubit-plugin-install --verify-only` to verify Radia panel
  startup registration when `radia` is installed, not just Cubit plugin
  binary hashes.

## 0.11.0 — Tier-2 sole-shipper + de-radia rename + `export` command verb

Released 2026-06-01.

Three coupled changes that let `cubit-mesh-export` and `radia` release
fully independently and clean up the radia-prefixed naming:

1. **Sole shipper of the Cubit plugin binary.**  `cubit-mesh-export`
   is now the ONLY package that ships and deploys the Cubit plugin;
   `radia` no longer bundles it.  Previously both wheels carried the
   plugin binary, which forced lockstep releases (a plugin fix meant
   re-releasing `radia` too).  The shared C++ source stays in the
   monorepo (`src/cubit_plugin/`) -- only the ship/release coupling is
   removed.

2. **De-radia file/module rename.**  The deployed plugin files drop
   the `radia_` prefix (the plugin is Cubit-side tooling, not a radia
   runtime component):

   | Old | New |
   |-----|-----|
   | `radia_cubit.ccm` | `cubit_mesh_export.ccm` |
   | `radia_cubit_mesh.cp312-win_amd64.pyd` | `cubit_mesh_curver.cp312-win_amd64.pyd` |
   | `radia_cubit_pybind.cpp` (C++ source) | `cubit_mesh_export_pybind.cpp` |
   | pybind module `radia_cubit_mesh` | `cubit_mesh_curver` |

   `cubit-plugin-install` removes any old `radia_cubit.*` files left by
   a pre-rename deployment, so Cubit does not load both the old and new
   `.ccm` and double-register the export commands.

3. **APREPRO command verb renamed `radia_export <fmt>` -> `export <fmt>`.**
   The mesh-export commands now extend Cubit's native `export` verb
   instead of using a separate `radia_export` verb:

   | Old command | New command |
   |-------------|-------------|
   | `radia_export netgen "f.vol"`  | `export netgen "f.vol"` |
   | `radia_export gmsh "f.msh"`    | `export gmsh "f.msh"` |
   | `radia_export vtk "f.vtk"`     | `export vtk "f.vtk"` |
   | `radia_export femeem "dir"`    | `export femeem "dir"` |
   | `radia_export meg "f.meg"`     | `export meg "f.meg"` |
   | `radia_export nastran "f.bdf"` | `export jmag_nastran "f.bdf"` |

   **Breaking**: existing `.jou` scripts calling `radia_export ...` must
   be updated to `export ...` (the old verb is removed -> Cubit reports
   `Unrecognized Keyword: 'radia_export'`).  Nastran is the one
   exception: Cubit has a built-in `export nastran` (different BDF
   format, no high-order support), so the plugin's BDF writer is exposed
   as `export jmag_nastran` to avoid shadowing the built-in.  The other
   five formats are not built-in Cubit export keywords, so they extend
   `export` cleanly.

## 0.6.0 — Japanese / Unicode path support

Released 2026-04-22.

All 6 mesh exporters now correctly write `.vol` / `.msh` / `.bdf` /
`.vtk` / `.meg` / FEMEEM `in.dat` to paths containing non-ASCII
characters (Japanese, Korean, Greek etc.) on any Windows codepage.

Before: `export netgen "C:/temp/日本語/coil.vol"` raised
  `"No mapping for the Unicode character exists in the target
   multi-byte code page."` and wrote no file.
After: same command writes the file (22,521 bytes on the reference
  sphere test).

### Implementation

New `src/cubit_plugin/utf8_path.hpp` provides one helper:

```cpp
std::filesystem::path u8_string_to_path(const std::string &s);
```

On Windows: `MultiByteToWideChar(CP_UTF8)` → `std::wstring` →
`std::filesystem::path(wstring)`.  Side-steps `std::string` →
`std::filesystem::path` implicit conversion (which uses the
system codepage = cp932 on Japanese Windows) and the previously-
used `CP_ACP` narrow-API pattern.

Applied to all 6 exporters:

* `ExportNetgenCommand.cpp` — `.vol` via `ng_mesh->Save()` + `.vol.json`
* `ExportGmshCommand.cpp` — `.msh` v4.1
* `ExportNastranCommand.cpp` — `.bdf`
* `ExportVtkCommand.cpp` — `.vtk`
* `ExportMegCommand.cpp` — `.meg` (FEMEEM / MAGIC)
* `ExportFemeemCommand.cpp` — `in.dat` + `node.dat` + ... (4 files)

### Smoke test harness

`smoke_test.py` driver wrapper now terminates with explicit
`exit 0` after `export netgen`.  On slower boxes (100号機)
Cubit's headless teardown access-violates (exit code 0xC0000005)
before the mesh DB destructor flushes the `.vol` writer; the
explicit exit forces shutdown through the normal exit handler so
the file is closed first.  Removes a flake that had the smoke
test failing ~50 % on 100号機 while passing on LAB.

### Compatibility

- Requires `radia >= 4.7.0` (matching plugin binary bundle version).
- C++ plugin rebuilt; `radia_cubit.ccl` / `radia_cubit.ccm` /
  `radia_cubit_mesh.pyd` refreshed in the bundle.

## 0.5.x and earlier

See git log.  Key points: 0.5.0 split from `radia` as an independent
package; 0.5.6 moved compact_netgen snapshot in-tree; subsequent
patch releases iterated on NetgenCurver projection robustness.
