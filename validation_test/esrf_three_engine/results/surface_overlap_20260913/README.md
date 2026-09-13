# ESRF6 Surface Audit and Known-Bad Mesh Cleanup

## Decision

The selected iron and FEM meshes have no detected positive-area boundary overlap
or strict transverse boundary intersection in the final NGSolve-mapped surface
tessellation audit. This is measured evidence, not a proof of Q2 patch
nonintersection, volume validity, or discretization convergence. ESRF6 numerical
acceptance remains **HOLD**. No field solve was launched during this audit.

The known-bad iron mesh is not the selected input. Six obsolete copies were
found and deleted: two each on hibino, mdx1, and mdx2. The current identity gate
was exercised on both bad hibino files before deletion: both were rejected;
the selected iron/FEM pair was accepted.

## Exact Inputs and Runtime

All geometric runs executed on hibino with NGSolve 6.2.2606 and Shapely 2.1.2,
using the isolated candidate environment at
`C:/temp/esrf6-newton-4d85e72cc/venv`. The audit imports NGSolve, not a different
Radia source overlay. Coordinates and tolerances are in metres.

| Input | SHA-256 |
|---|---|
| Selected iron, 2408 HEX | `fdd1387281e550e27eab435cfac3ec54e0b9d1339438699fdc50c852dfcafc74` |
| Selected FEM, 134686 TET | `dfc12b84f80db1fa17fb5012b6f87c072e8aa1fb3c085b61269f047a879a54ac` |
| Known-bad iron, 4768 HEX | `0d4449962b60757cba78e771a136b2805a4323e26a1ee6d7ae482751bb7e5a5f` |

Selected paths are recorded in each JSON. They are respectively
`C:/temp/radia-ex6-diag/assets_conforming/model.vol` and
`C:/temp/radia-esrf-fem/example6/coil_yoke_kelvin.vol`.

## Method and Limits

The audit enumerates all registered boundary elements, including once-registered
two-owner material interfaces. It separately counts volume-face owners, missing
one-owner boundary registrations, duplicate registrations, and nonmanifold faces.
It does not expand a normal two-owner interface into two duplicate surfaces.

NGSolve `GetTrafo` evaluates a reference grid on each Q2 boundary face. Grids with
2 and 4 subdivisions are triangulated. Radius-binned cKDTree searches identify
possible triangle pairs; centre/radius bounds and axis-aligned bounding boxes
screen them before geometric intersection tests. Same-parent subtriangles are
excluded. Bounds enclose tessellated triangles, not unsampled curved patches.

The final v3 detector uses dominant-axis projection, an independent 2D
separating-axis test (SAT), and fixed-precision polygon intersection. It reports
positive-area coplanar overlaps and strict interior transverse piercings.

- Plane-distance and SAT contact tolerance: `1e-10 m`.
- Coplanar area threshold: `1e-16 m^2`.
- Polygon overlay grid: `1e-12 m`.
- Transverse piercing requires reference barycentric margin greater than `1e-8`.
- SAT excludes edge/vertex contacts and strips no wider than its tolerance.
  Grid rounding can erase tiny true overlaps. These thresholds are detection
  limits, not an accuracy claim or a justification for relaxing solver tolerances.
- Tangencies, subthreshold features, same-parent folding, unregistered internal
  intersections, and intersections between unsampled parts of curved patches
  are not certified absent. No full-volume intersection proof is claimed.

## Results

| Mesh | Subdivisions | Mapped triangles | Candidate pairs | Detected parent pairs |
|---|---:|---:|---:|---:|
| Iron, v1 | 2 | 17600 | 83458 | 0 |
| FEM, v1 | 2 | 48720 | 226829 | 0 |
| Iron, final v3 | 4 | 70400 | 182279 | 0 |
| FEM, final v3 | 4 | 194880 | 529761 | 0 |

Iron has 2200 one-owner boundary faces and 6124 two-owner internal faces.
FEM has 5564 one-owner boundary faces and 6616 once-registered two-owner
material-interface faces. Both have zero missing one-owner registrations,
duplicate boundary face IDs, and faces with more than two volume owners.

For the iron's mapped grid, maximum parent planarity deviation is
`1.185e-15 m` and maximum edge-chord deviation is `1.041e-17 m`.
FEM has 5564 sampled nonplanar boundary faces; its curved geometry was not
silently replaced by corner-only planar elements.

## Detector Regression and Controls

The original v2 detector reported one FEM overlap at boundary IDs 3840 and 3878
(`iron_air_interface`). Its arbitrary rotated 2D projection caused a GEOS overlay
robustness failure near an almost-coincident vertex, returning a spurious
`1.1472689294122711e-6 m^2`. The original JSON, including coordinates, is retained
as `fem_overlap_n4.json`. The actual triangles touch at a point; the result is
**not an observed mesh defect**. Final v3 rejects it independently by SAT.

Controls cover partial overlap, different subdivisions, opposite square
diagonals, normal edge adjacency, separated triangles, transverse intersection,
and that measured grazing case. Full synthetic meshes also verify that a
normal two-tetrahedron shared interface is not reported and a partially
overlapping pair is detected. All control assertions passed on hibino.

A positive control translates the grazing triangle by `1e-6 m` and detects a
thin overlap of `5.098975283407465e-13 m^2`. The two nondegenerate coordinate
projections agree to relative `2.405e-12`; the third projection has zero normal
component and is explicitly unusable for area. See
`surface_overlap_precision_controls.json` rather than treating a degenerate
projection's zero area as negative evidence.

Final script SHA-256:
`8b2cf6f81c3ce11d5e2d01aeaa693ccb03a935ccb4df7e84c3e04eefe533a230`.
The final script and control generators are committed here. Superseded v1/v2
scripts remain only in the durable operational archive for reproducing the
recorded detector failure, not as additional maintained validators.

## Cleanup and Retention

Durable raw archive:
`S:/Radia/validation_artifacts/esrf6_mesh_audit_20260913/`, separated by host.
It holds full before/after inventories, deletion receipts, scripts, original
detector versions, and synthetic mesh inputs. These operational inventories are
not tracked source. Selected numerical JSON and final audit code are tracked here.

| Host | VOL files before | Known-bad deleted | VOL files after | Known-bad remaining |
|---|---:|---:|---:|---:|
| hibino | 188 | 2 | 186 | 0 |
| mdx1 | 158 | 2 | 156 | 0 |
| mdx2 | 162 | 2 | 160 | 0 |

Scope was non-link top-level `C:/temp/radia*`, `esrf*`, and `ctype*` directories
and their `.vol` files. The two removed relative paths on every host were
`radia-validation/assets/example_6_quadrupole/model.vol` and
`radia-xfer6/assets/example_6_quadrupole/model.vol`. Each deletion used a fresh
hash, explicit resolved-path and parent-reparse checks, no Git worktree, no
active compute process, and nonrecursive `Remove-Item -LiteralPath`.

Audit scratch cleanup follows durable recovery, hash verification, and this
evidence commit. The candidate environment and selected physical meshes remain
temporarily because ESRF6 is unfinished. Their owner is the execution task on
branch `codex/hdiv-energy-newton-consistency-20260913`, not the release manager.
Trigger: finish ESRF6, recover required inputs/results with verified hashes, then
remove its candidate/staging/job workspace. Other jobs' or CI-owned workspaces
are outside this cleanup. This report does not claim that all historical
`C:/temp` contents on all hosts have been cleaned.
