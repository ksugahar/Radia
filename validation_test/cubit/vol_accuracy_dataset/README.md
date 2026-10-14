# Complex Cubit `.vol` volume-accuracy dataset

This validation lane measures whether moderately complex Cubit CAD survives
`export netgen` as the same geometric volume in NGSolve. It is not an analytic
primitive test. The eight cases deliberately exercise failure-prone CAD and
meshing operations:

- repeated Boolean holes in a perforated flange;
- closed holes plus an open notch in a thin busbar;
- a Boolean-united stepped spacer with a central bore and bolt holes;
- a trimmed circle-to-rectangle loft with a periodic circular-curve seam;
- a webcut, doubly-curved half torus;
- coarse and refined 355-degree circular-profile sweep twins;
- a Boolean cylinder/sphere union with a mixed tet/wedge boundary layer.

Seven are acceptance cases. The loft is a regression for a periodic-curve
seam bug that formerly collapsed one circular mesh edge and reversed 4 of 1536
order-2 Jacobian samples and 9 of 3000 order-3 samples even though the volume
errors were only 0.104% and 0.0633%. After the seam fix, both orders have zero
invalid or orientation-flipped samples; their measured volume errors are
0.0351% and 0.0229%. This case ensures that volume agreement can never hide a
folded high-order map again.

The coarse 355-degree sweep is intentionally diagnostic. At 12 profile and 48
sweep intervals it is valid at order 2, but order 3 has 13 orientation-flipped
samples even though its CAD-volume error is only 0.0349%. The otherwise
identical 16/64 refined twin has zero flips at orders 2 and 3. This records a
mesh-resolution boundary rather than an exporter seam bug: order elevation
alone cannot replace geometric refinement on a tightly curved swept tube.

The collector deliberately disables the CAD-error threshold while measuring so
that poor low-order results remain useful observations. Rows therefore report
the structural/Jacobian gate and the 1% CAD-volume accuracy gate separately.
Acceptance requires every acceptance-role case to pass the structural gate at
all orders and to pass both gates at orders 2 and 3. Diagnostic rows remain in
the dataset specifically so a small volume error cannot erase a known mapping
failure.

For every case, the same Cubit mesh is exported at curve orders 1, 2, and 3.
The exporter sidecar supplies Cubit's ACIS volume for the named material.
Normal Python then runs the production `check_consistency` path, which reloads
the `.vol` with NGSolve and evaluates
`Integrate(CF(1), mesh, definedon=mesh.Materials(name))`. Each learning row
contains both raw volumes, signed and absolute percentage error, mesh size,
mapping sample count, minimum Jacobian metrics, invalid/negative/orientation-
flip counts, and source-command hash.

## Generate

Save and close interactive Cubit first. Cubit's embedded Python lacks NGSolve,
so `generate.py` intentionally launches Cubit as a separate batch process and
performs the NGSolve measurements after it exits.

```powershell
python validation_test/cubit/vol_accuracy_dataset/generate.py
```

Temporary `.vol` artifacts are created only under `C:\temp` and deleted after a
successful run. Use `--keep` to retain them for diagnosis. A failure always
retains the work directory and `cubit.log`.

Tracked outputs:

- `volume_accuracy_dataset.json`: run provenance, aggregate statistics, rows;
- `volume_accuracy_rows.jsonl`: one flat training/analysis row per
  model/order/material.

These values are empirical observations, not handwritten truth. Regenerate
them after changing Cubit, `cubit-mesh-export`, Netgen, NGSolve, CAD projection,
or the case definitions.
