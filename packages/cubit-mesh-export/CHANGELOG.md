# Changelog

All notable changes to `cubit-mesh-export` — the high-order curved
mesh export package for Coreform Cubit (Netgen / GMSH / Nastran /
VTK / MEG / FEMEEM writers + Python bindings for consistency checks).

## 1.0.4 - GUI release launch correction

- Pass the installed plugin directory before the GUI probe journal, preventing
  Cubit 2025.12 from treating appended INI plugin arguments as journal files.
- Retain launcher output with GUI evidence; keep nonzero process exits fatal.

## 1.0.3 - Review corrections

- Declare independently maintained runtime support and preserve its BSD notices.
- Use Cubit-owned environment settings, status/log schemas and HTTP identity;
  remove unsupported example providers and Radia installation advice.
- Add `sculpt`, `mesh-quality`, and `youtube` extras with metadata-driven CI.
- Preserve old Radia logs/caches on disk but do not automatically migrate them;
  `CUBIT_MCP_STATE_DIR` selects the new product-owned state.
- Harden the private Netgen quality reference with checked Jacobians and explicit
  MSH entities; cache supported reference-node layouts by element order.

## 1.0.2 - Independent Cubit GUI and MCP release

- Validate coinciding vfrac output paths before checking Sculpt availability;
  run all license-free MCP contracts with their declared test dependencies.
- Own `cubit_mesh_export.mcp` and `mcp-server-cubit`, included by default
  without Radia or radia-mcp. Include the required runtime support inside this distribution; no shared private package is needed.
  Move Cubit API/reference, headless journal workflows and their tests here;
  retain tool contracts but remove the old `radia_mcp.cubit` module path.

- Package the Cubit-private export menu, WorkflowToolbar assets, installer
  and default smoke journal without requiring the Radia distribution.
- Deploy independently of any installed Radia version. Optional combined
  integration is checked explicitly with `cubit-plugin-install --check-radia-compat`.
- Add the LAB/100-only `release_quad.py cubit-dual` exact-artifact gate; never
  reinstall Radia/MCP or deploy Cubit to compute hosts in this lane.
- Release Claro-owned menu actions on application shutdown to prevent
  Windows fast-fail during interpreter teardown.
- Add isolated-wheel registration and explicitly scoped GUI smoke probes.
- Fix all six official toolbar launchers for Cubit's `play` execution, which
  does not define `__file__`; use the declared toolbar working directory.
- Verify real dialog accept/cancel/checker-failure behavior and official package
  import, restart and six-button dispatch using an installed wheel without Radia.
  Tests restore the exact original Cubit preferences after the persistence lane.
- Remove Radia's old GUI, installer, startup and compatibility bridges; migrate
  callers to the exporter rather than retaining legacy entry points.
- Move the installed GUI smoke CLI, embedded probe and focused tests to the
  exporter. Require a clean process exit and terminate only the test-owned
  process on failure. Cross-host GUI hashes now cover installed exporter files.
- Show standalone checker startup failures and open Netgen output in a viewer
  only after successful validation.

## 1.0.1 - Radia 5 compatibility

- Extend the reciprocal Radia compatibility declaration through exactly
  5.0.0. No export implementation, native source, command or label contract
  changes from 1.0.0. Combined installed-wheel and export/check-vol acceptance
  remains a release gate; changing a version bound alone is not validation.

## 1.0.0 - Stable standalone mesh export

- Standalone wheel acceptance passed on LAB without Radia or radia-mcp:
  headless APREPRO sphere, order-2 export, strict labels, CAD measures,
  NGSolve reload and mapped Jacobians. Evidence is retained in
  `validation_test/cubit_mesh_export/standalone_1_0_0_lab_result.json`.
- Radia 4.95.91 declares an exporter upper bound of 0.999.999. Its
  compatibility gate rejects 1.0.0; keep 0.14.17 in that combined environment
  until Radia publishes an updated compatibility declaration.

- Promote the existing standalone export and checker interfaces to 1.0.0;
  no command or label-contract migration is required from 0.14.17.
- Preserve explicitly labelled same-material internal surfaces by rejecting
  unsupported ambiguous boundaries instead of silently deleting their labels.
- Refresh the bundled command backend with the labelled-surface protection.

- Keep mandatory binary hash/size and manifest checks active when building
  from an sdist without the native source tree.
- Replace checkout-time freshness guesses with content-addressed provenance
  for the native source tree and both mandatory `.ccm` / `.pyd` payloads.
- Align package metadata with the shipped `cp312-win_amd64` wheel: CPython
  3.12 and 64-bit Windows are the supported runtime contract.

## 0.14.17 - Refreshed Coreform Cubit backend

Released 2026-09-11.

- Refresh the bundled Coreform Cubit 2025.12 command backend from the current
  reviewed main source while preserving the standalone mesh-export API.

## 0.14.16 - Standalone installation guidance and refreshed backend

Released 2026-09-10.

- Lead the package documentation with the standalone installation path and
  make the Coreform Cubit, Radia, and smoke-test prerequisites explicit.
- Refresh the bundled Coreform Cubit 2025.12 command backend from the current
  reviewed source before publication.

## 0.14.15 - High-order mapping checks and refreshed native payload

Released 2026-09-04.

- Validate high-order `.vol` element mappings before accepting an exported
  mesh, including the Q2-axis geometry exercised by the Radia validation lane.
- Refresh the bundled Cubit 2025.12 command backend from the reviewed source
  and keep the installed toolbar and smoke-test contract aligned with it.

## 0.14.14 - Complete native payload rebuild and Claro export menu

Released 2026-08-30.

- Rebuilt both mandatory Coreform Cubit 2025.12 payloads from the same source:
  the APREPRO command backend (`cubit_mesh_export.ccm`) and the Netgen curver
  (`cubit_mesh_curver.pyd`).
- Made the focused native build worktree-relative and propagate both artifacts
  into the wheel source only after successful builds and SHA-256 verification.
- Extended the pre-wheel freshness gate to reject either stale native payload;
  the installer can no longer receive a current `.ccm` beside an old `.pyd`.
- Emit a native `cp312-cp312-win_amd64` wheel directly with
  `Root-Is-Purelib: false`; release CI no longer retags a pure wheel after the
  fact.
- Bind the gitignored curver to a content-addressed GitHub Release asset via a
  checked manifest; CI verifies its size and SHA-256 before and after wheel
  assembly instead of trusting the mutable historical asset name.
- Initialize Netgen's Windows DLL search path when the package is imported so
  the bundled curver loads in a fresh Python process without an undocumented
  preceding `import netgen`.
- Restored the user-facing **Export** menu through Cubit's official Claro API,
  while keeping PySide6 limited to dialogs and retaining the WorkflowToolbar.
- Refresh any previously imported Coreform WorkflowToolbar in place and detect
  stale toolbar scripts during `--verify-only`.

## 0.14.13 - Reflection-invariant Kelvin validation meshes

Released 2026-08-29.

- Added the exact Cubit/ACIS C-yoke route used by Radia's HDiv/Omega
  cross-validation. The exporter builds a locally refined physical-air sphere
  plus translated Kelvin sphere without a finite outer box and preserves the
  required one-to-one periodic identification.
- Mesh one half and reflect the meshed volumes for exact z-symmetry, remove
  same-material webcut seam descriptors without removing conforming shared
  topology, and fail loudly when automatic Kelvin construction cannot satisfy
  its geometry or identification contract.
- Added a helpers-only installation path so Cubit Python helpers can be
  refreshed without copying or replacing native `.ccm`/`.pyd` binaries.

## 0.14.12 - Validated Nastran mesh interchange

Released 2026-08-27.

- Promoted `export nastran_bdf` to the primary solver-neutral BDF command while
  retaining `export jmag_nastran` as a separately registered compatibility
  alias for existing journals.
- Corrected Nastran interchange semantics: `dimension 2` no longer flattens
  volume meshes or discards off-plane GRID coordinates, surface blocks receive
  PSHELL rather than PSOLID, sideset PIDs cannot collide with block PIDs, and
  nodesets use strict-reader-compatible fixed-field SET1 cards.
- Expanded real-Cubit coverage across TET, HEX, WEDGE, PYRAMID, TRI, QUAD,
  orders 1--2, dimension filtering, groups, properties, and alias behavior;
  added an independent pyNastran referential-integrity/JSON validation gate.

## 0.14.11 - Active Python ABI plugin refresh

Released 2026-08-26.

- Rebuilt the bundled Coreform Cubit backend after pinning the active Python
  ABI and isolating the NGSolve 6.2.2606 build environment. The mesh-format and
  Python API contracts are unchanged.

## 0.14.10 - Netgen 6.2.2606 ABI migration

Released 2026-08-25.

- Fixed high-order edge interpolation at periodic Cubit curve seams.  Cubit
  can return the trimmed-edge endpoint for every interior closest-point query;
  the exporter now detects that stalled projection and evaluates the intended
  point on the locally unwrapped periodic branch.  The circle-to-rectangle
  loft no longer folds one HEX element despite apparently converged volume.
- Run complex-shape CAD-volume and surface-area comparisons in an isolated
  NGSolve process and persist the order 1--5 results together with sampled
  scaled-Jacobian and within-element orientation data.
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
