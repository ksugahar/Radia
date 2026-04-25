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
| `em_1-1_full.jou` | TODO | — | No webcut, full mesh.  Simplest; Auto-Kelvin copies a full sphere 1:1 which should work without any plugin fix. |
| `em_1-2_half_z.jou` | **WORKING** | [em_sample_mu1000.json](../../../../tests/panels/golden/em_sample_mu1000.json) | Z-plane webcut on air sphere (mesh seam for Kelvin copy-mesh).  Yoke has full x, full z; only coil + coil-mirror about z gives the physical field.  Verified end-to-end at NI=1 and NI=2000, mu_r=1000, PARDISO. |
| `em_1-4_quarter_xz.jou` | **WORKING** (2026-04-25) | [em_quarter_xz_mu1000.json](../../../../tests/panels/golden/em_quarter_xz_mu1000.json) | True 1/4 reduction (x>=0 AND z>=0) via `add_kelvin_cubit(reduction={x: bn=0, z: ht=0})`.  Air quarter-sphere via `intersect` with an octant brick (cleaner than webcut+delete).  `sym_bn=0_x` and `sym_ht=0_z` sidesets propagated to Netgen bcnames; FEM Omega picks Dirichlet on z=0 automatically.  Regression guard only, not a physics-accuracy match to the ELF published -976 mT (coil is not x/z-symmetric). |
| `em_1-8_eighth.jou` | TODO (blocker: Kelvin mapping) | — | Offset-then-cut Kelvin strategy geometrically fails for 1/8 (diagonal offset does not intersect the x=0/y=0/z=0 planes so the Kelvin cut faces cannot coincide with the air cut faces).  `add_kelvin_cubit(reduction=...)` raises NotImplementedError for three-axis reduction.  Needs a different Kelvin mapping (e.g. 7 reflected copies, or spherical inversion without translation). |

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
giving `B_z = -976 mT` at `NI=2000, mu_r=1000`.  Matching that
number in FEM requires:

- `em_1-4_quarter_xz.jou` + its Add-Kelvin support (see blockers above).
- `sym_tangential` / `sym_normal` sidesets on the X=0 and Z=0 planes
  so `calc_accel_magnet` swaps the Dirichlet/natural pair correctly
  per formulation.
- Same mesh density order as the ELF reference for convergence.

The currently-shipping `em_1-2_half_z` variant is NOT a direct
substitute for that reference — it's a half-z reduction that
doesn't match the quarter-xz geometry.  The `em_1-2_half_z` golden
is an internal regression guard, not a physics-accuracy check.
