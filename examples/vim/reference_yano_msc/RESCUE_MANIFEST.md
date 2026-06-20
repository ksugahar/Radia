# yano-type past-data rescue (from git history)

Rescued 2026-06-19 from the Radia git history into `C:\temp\yano_rescue\`.
These files were deleted from the working tree during the drop-yano /
HDiv-VIM consolidation. They are recoverable from git anytime, but are
extracted here as the **yano-MSC reference corpus for the CEFC 2026 paper**
(C-yoke: yano-MSC vs HDiv-VIM, loop-mode / iteration headline).

## Provenance (which commit deleted what)

| Item | Deleted in | Date | Recovered from |
|---|---|---|---|
| `examples/c_type_electromagnet/` (whole tree, 35 files) | `20cc1696` "Consolidate in-flight drop-yano working tree" | 2026-06-17 | `20cc1696^` |
| `examples/mmm_eigenvalue_study/null_removed_mmm_msc.py` | `2a0adc3d` "Remove yano-type loop-star / deflation / loop-projection" | 2026-06-09 | `2a0adc3d^` |
| `tests/test_hex_demag_convergence.py` | `9facc2e7` "hdiv-vim: gate removed yano demag path" | 2026-06-19 | `9facc2e7^` |

Model lineage: `examples/c_type_electromagnet` is a direct port of the
**ELF_MAGIC CEFC-2020 C-Type model** (`S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\Cubit`,
Trelis.jou, 1x1x1 / 3x3x3 / 6x6x6 hex). Quarter model + Image symmetry
`+x-z`, racetrack coil, C-type yoke. Verified against ELF_MAGIC.

## KEY: yano-MSC reference numbers (nonlinear, Newton, BH curve)

`nonlinear/quarter/{LU,bicgstab,hacapk}/*DOF.json` — hex MSC, 6 DOF/face.
These are the **baseline the HDiv-VIM is compared against**. Note the large
nonlinear-Newton iteration counts AND the HACApK *linear* Krylov iteration
counts (1500-2700) — the latter is the loop-mode pollution the loop-free
HDiv-VIM collapses.

| solver   | DoF    | nelem | nonl_it | lin_it | Bz (mT)   | t_solve (s) |
|----------|--------|-------|---------|--------|-----------|-------------|
| LU       | 7200   | 1200  | 82      | 0      | -963.783  | 98.3        |
| LU       | 18900  | 3150  | 71      | 0      | -961.505  | 1393.7      |
| BiCGSTAB | 7200   | 1200  | 49      | 0      | -963.868  | 111.9       |
| BiCGSTAB | 18900  | 3150  | 61      | 0      | -961.515  | 635.8       |
| HACApK   | 7200   | 1200  | 182     | 1698   | -963.900  | 32.2        |
| HACApK   | 18900  | 3150  | 174     | 1506   | -961.480  | 99.1        |
| HACApK   | 165600 | 27600 | 214     | 2686   | -954.382  | 2607.1      |

(`165600DOF` = 20x20x20 hex, HACApK eps 1e-4, leaf 10, eta 2.0, compression
0.042, peak 9.6 GB, 8 threads.)

## Inventory (38 files)

- `examples/c_type_electromagnet/generate_hex_mesh.py` (330 L) — **the hex mesh generator** (Cubit -> Netgen .vol, coarse/medium/fine = 1/3/6 intervals)
- `examples/c_type_electromagnet/generate_quarter_mesh.py` — quarter-model hex mesh
- `examples/c_type_electromagnet/coil.step`, `coil_geometry.step` — racetrack coil geometry
- `examples/c_type_electromagnet/dipole_with_coilbuilder.py` — CoilBuilder version of the drive coil
- `examples/c_type_electromagnet/mu=1000/**` — LINEAR (mu_r=1000) ELF-Radia verify (full / quarter / x-mirror / z-mirror) + IMA experiments + wedge experiment
- `examples/c_type_electromagnet/nonlinear/**` — NONLINEAR BH: block-jacobi-Newton benchmark, LU/bicgstab/hacapk verify (full + quarter), the 9 golden DOF JSONs above, summary.json, BH.txt
- `examples/mmm_eigenvalue_study/null_removed_mmm_msc.py` — loop-star null-space removal prototype (loop-mode study)
- `tests/test_hex_demag_convergence.py` — hex demag p/h-convergence golden test (yano path)

## Runnability caveat

Per CLAUDE.md (2026-06-19), the mesh-less hex/wedge yano-MSC `rad.Solve`
path now **raises `Radia::Error203`** (directs to `radia.vim.soft_iron_from_mesh`).
These scripts solved hex soft iron via the OLD mesh-less path, so they will
NOT run as-is on current radia. The durable value here is:
  1. the **golden reference numbers** (JSONs above) — frozen, citeable;
  2. the **hex geometry** (generate_hex_mesh.py + .step) — reusable to
     rebuild the C-yoke hex mesh for an HDiv-VIM re-run;
  3. the **ELF-Radia verification provenance** (mu=1000 + nonlinear verify_elf_radia.py).
