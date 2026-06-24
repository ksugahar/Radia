# six-face surface-charge MSC reference (frozen baseline for the HDiv-VIM C-yoke)

This directory holds the **frozen six-face surface-charge collocation-MSC (distortion-element)
reference numbers** that the live HDiv-VIM C-yoke head-to-head
([`../hdiv_cyoke_headtohead.py`](../hdiv_cyoke_headtohead.py)) is compared against.

The historical collocation six-face surface-charge method has been replaced in live Radia by
the canonical multipole-moment MMM surface-charge formulation; these are **reference data
only**, not a runnable reproduction of the old EIEM2 collocation kernel. Current
mesh-less hex/wedge/pyramid soft iron still solves through `rad.Solve`, but it
uses multipole-moment MMM. Until now these numbers lived only as a comment in
`hdiv_cyoke_headtohead.py` ("the yano JSON is not committed to this repo"); this
directory **commits the source JSONs** so the embedded `YANO_REF` has a tracked
provenance.

## Model

C-type yoke (quarter model + Image symmetry `+x-z`, racetrack coil), a direct
port of the **ELF_MAGIC CEFC-2020 C-Type model**
(`S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type`, Trelis.jou, 1/3/6 hex
intervals). Nonlinear BH (`BH.txt`), 6-DOF/face hexahedral surface charge.
Recovered from git history (commit `20cc1696`, deleted 2026-06-17) — see
[`RESCUE_MANIFEST.md`](RESCUE_MANIFEST.md). Geometry generator:
[`generate_hex_mesh.py`](generate_hex_mesh.py) (requires Cubit).

## golden/ — surface-charge MSC nonlinear (Newton) reference matrix

| file | solver | DoF | nelem | nonl_it | lin_it | Bz (mT) | t_solve (s) |
|---|---|---|---|---|---|---|---|
| `LU_7200DOF.json`        | LU       | 7200   | 1200  | 82  | 0    | -963.783 | 98.3   |
| `LU_18900DOF.json`       | LU       | 18900  | 3150  | 71  | 0    | -961.505 | 1393.7 |
| `BiCGSTAB_7200DOF.json`  | BiCGSTAB | 7200   | 1200  | 49  | 0    | -963.868 | 111.9  |
| `BiCGSTAB_18900DOF.json` | BiCGSTAB | 18900  | 3150  | 61  | 0    | -961.515 | 635.8  |
| `HACApK_7200DOF.json`    | HACApK   | 7200   | 1200  | 182 | 1698 | -963.900 | 32.2   |
| `HACApK_18900DOF.json`   | HACApK   | 18900  | 3150  | 174 | 1506 | -961.480 | 99.1   |
| `HACApK_165600DOF.json`  | HACApK   | 165600 | 27600 | 214 | 2686 | -954.382 | 2607.1 |

`summary_165600.json` is the run record for the 20x20x20 / 165600-DoF case
(HACApK eps 1e-4, leaf 10, eta 2.0, compression 0.042, peak 9.6 GB, 8 threads).

## Honest comparison caveat (Repository-First)

These are the six-face surface-charge **HACApK + Newton** (or LU / BiCGSTAB + Newton) runs.
The large *linear* (Krylov) iteration counts of the HACApK rows (1506-2686) are
the loop-mode pollution the loop-free HDiv-VIM (ker(B) field-null) collapses to
~6 mesh/mu_r-independent Newton iterations. **However**, a strictly fair
same-solver head-to-head (identical Picard + Block-Jacobi for both methods) at
165600 DoF has **not** been run on current radia — the HDiv-VIM inner +N solve
does not yet scale past ~20-40k DoF (see the `status_2026_06_09` note in
`../hdiv_cyoke_headtohead.json`). Read these as the **historical six-face surface-charge
baseline**, not as a head-to-head win at 165600 DoF.
