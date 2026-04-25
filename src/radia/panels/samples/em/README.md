# EM sample corpus (accelerator magnet, C-type yoke)

## Coil input policy (2026-04-25, user-set)

Analytical coils and PEEC coils use DIFFERENT file formats by
design.  The EM panel (`calc_accel_magnet` / `calc_accel_msc`)
accepts **only** a Python module.  The PEEC path (IH panel,
`calc_peec_inductance` / `calc_peec_bem`) accepts **only** STEP.

| Panel / calc | Coil input | Loader |
|--------------|-----------|--------|
| EM (`calc_accel_{magnet,msc}`) | `.py` module exposing `build_coil() -> CoilBuilder` | `_load_coil_script` (literal import) |
| PEEC (`calc_peec_*`) | `.step` swept-solid coil | `coil_from_cad.filaments_from_step` (walker + nwinc x nhinc cross-section subdivision) |

Why the split:

- CoilBuilder (analytical Biot-Savart) needs exact straight + arc
  segments.  A Python `build_coil()` literal is the highest-
  fidelity source.  STEP's walker-based decomposition is
  approximate (often collapses a racetrack to a single arc).
- PEEC needs the conductor's **cross-section** (to subdivide into
  nwinc x nhinc filaments for skin / proximity effects).  The
  swept-solid STEP carries that cross-section implicitly; a
  centerline `.py` does not.

CoilBuilder -> STEP (swept solid, for PEEC):

```python
# CoilBuilder is the single source of truth for coil geometry.
# For PEEC workflows that require a STEP swept solid, export via:
coil.write_step("coil.step")   # PEEC pipeline consumes this
```

The reverse path (STEP -> CoilBuilder) is intentionally not
supported: STEP's walker-based centerline decomposition is
approximate (often collapses racetracks to a single arc), so
CoilBuilder must remain the authoritative `.py` definition.

---



Known-good (and known-missing) variants of the C-type dipole EM
sample for `calc_accel_magnet.py`, derived from the canonical ELF
CEFC 2020 model at `examples/cubit_panels/accel_magnet/yoke.jou`.

Each variant is a different symmetry reduction of the same physical
dipole.  The `.jou` defines the yoke geometry + air sphere; the
launcher's Auto-Kelvin checkbox (or `auto_add_kelvin_from_current_model`)
adds the Kelvin exterior sphere and the GND nodeset.

| Variant | Status | Golden | Notes |
|---------|--------|--------|-------|
| `em_1-1_full.jou` | **WORKING** (2026-04-25) | [em_full_mu1000.json](../../../../tests/panels/golden/em_full_mu1000.json) | True 1/1 FULL geometry: yoke at full extent in x, y, z; air sphere with z=0 mesh seam (geometry full, seam for Kelvin copy-mesh anchor only).  Auto-Kelvin auto-detects `symmetry=['z']` from the air's volume centroids straddling z=0 (centroid-based heuristic, fixes a vertex-range-based bug that spuriously detected 'x' too).  Largest of the 4 variants in DOFs (~8x the 1/8 model). |
| `em_1-2_half_z.jou` | **WORKING** | [em_sample_mu1000.json](../../../../tests/panels/golden/em_sample_mu1000.json) | Z-plane webcut on air sphere (mesh seam for Kelvin copy-mesh).  Yoke has full x, full z; only coil + coil-mirror about z gives the physical field.  Verified end-to-end at NI=1 and NI=2000, mu_r=1000, PARDISO. |
| `em_1-4_quarter_xz.jou` | **WORKING** (2026-04-25) | [em_quarter_xz_mu1000.json](../../../../tests/panels/golden/em_quarter_xz_mu1000.json) | True 1/4 reduction (x>=0 AND z>=0) via `add_kelvin_cubit(reduction={x: bn=0, z: ht=0})`.  Air quarter-sphere via `intersect` with an octant brick (cleaner than webcut+delete).  `sym_bn=0_x` and `sym_ht=0_z` sidesets propagated to Netgen bcnames; FEM Omega picks Dirichlet on z=0 automatically.  Regression guard only, not a physics-accuracy match to the ELF published -976 mT (coil is not x/z-symmetric). |
| `em_1-8_eighth.jou` | **WORKING** (2026-04-25) | [em_eighth_mu1000.json](../../../../tests/panels/golden/em_eighth_mu1000.json) | True 1/8 reduction (x>=0, y>=0, z>=0) via single-axis Kelvin offset (default along the first reduction axis = x).  The Kelvin sphere's three webcut planes through its centre coincide with the AIR's two non-offset sym planes (sym_<bc>_y and sym_<bc>_z combined with air); the third (perpendicular to offset_dir) becomes a new "kelvin_far" sideset that calc_accel_magnet treats as always-Dirichlet (an "infinity plane" extension of the GND vertex).  Physically-impossible all-bn=0 (B=0) is rejected with a clear ValueError.  C++ exporter uses an all-or-nothing Kelvin identification policy: when copy-mesh leaves a few corner vertices unmatched at large geometry scale, identification is skipped entirely and `add_periodic_kelvin` reconstructs at solve time. |
| `em_elf_quarter_xz.jou` | **WORKING** (2026-04-25) | [em_elf_quarter_xz_mu1000.json](../../../../tests/panels/golden/em_elf_quarter_xz_mu1000.json) | **PHYSICS-ACCURACY** ELF reproduction (NOT just a regression guard).  Yoke geometry is built directly from the ELF .meg file (`em_elf_yoke_builder.py`): 13 CAD volumes including the slanted-face pole-tip bevel, then `unite` -> tetmesh at 5 mm for face-conformality with the air sphere.  Coil from `em_elf_coil.py` (Y_CENTER = 131.25 mm matching ELF's `coil_model.py`).  FEM Omega gives **Bz = -240.0 mT vs ELF -228.1 mT (5.2% match)** -- the closest physics benchmark in this corpus.  Requires the ELF .meg at the LAB path or via `ELF_MEG_PATH` env var. |

## Workflow for each variant

1. Author / copy the `.jou` with the right bricks + air-sphere webcut
   pattern for the symmetry reduction.
2. In the Cubit launcher dialog, check "Add Kelvin open boundary
   (auto)" and select the matching symmetry ("Full", "Half (Z)",
   "Quarter (X,Z)", "Eighth").
3. Export via `radia_export netgen ... order 2 overwrite`.  The
   `.vol` will have `cd3names` (GND) and `pointelements` entries
   from the nodeset propagation added 2026-04-24 (commit 07b15414).
4. Run `calc_accel_magnet.py --coil-script em_sample_coil.{step|py}
   --vol <generated>.vol --formulation omega --material linear
   --mu-r 1000 --solver pardiso`.
5. First-time runs populate the golden JSON under
   `tests/panels/golden/em/<variant>.json`.  Later runs of
   `tests/panels/test_em_golden.py` lock the numbers.

## The full-ELF-CEFC-2020 reference

The published ELF reference (`examples/c_type_electromagnet/mu=1000/
quarter/verify_elf_radia.py`) is the **quarter_xz** variant with
Radia MSC (direct hex-element extraction via `IMA='+x-z'` symmetry),
giving `B_z = -228.1 mT` at `NI=2000, mu_r=1000` (verified by
running the ELF reference 2026-04-25; earlier README revisions
incorrectly cited -976 mT, which was a transcription error).
Matching that number in FEM requires:

- `em_1-4_quarter_xz.jou` (now WORKING) at the ELF-matched yoke
  geometry (NOT the simplified bricks currently shipped -- ELF uses
  17mm pole bevel + 25mm thick C-back at specific positions).
- The ELF racetrack coil (port `examples/c_type_electromagnet/
  mu=1000/coil_model.py` to a CoilBuilder script with Y_CENTER =
  131.25 mm, NOT the 151.25 mm in `em_sample_coil.py`).
- Mesh refinement matching the ELF discretization for FEM
  convergence to ~1% accuracy.

The currently-shipping samples (`em_1-1_full`, `em_1-2_half_z`,
`em_1-4_quarter_xz`, `em_1-8_eighth`) are INTERNAL REGRESSION GUARDS
for the reduction-mode + Auto-Kelvin pipeline plumbing.  None of
them match ELF physics: their yoke uses a simplified brick-only
geometry and their coil (`em_sample_coil.py`) sits at y=151.25mm
instead of the ELF y=131.25mm.  The 4 reductions give 4 different
B numbers because the sym BCs impose different constraints on a
non-symmetric source.

A future ELF-matched sample (`em_elf_quarter_xz.jou` +
`em_elf_coil.py`) would close this loop and target -228.1 mT.
