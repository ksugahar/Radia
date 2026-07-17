# ESIM cross-formulation consistency (BIE vs FEM-Kelvin)

Executable validation for the nonlinear per-element ESIM claim:
the scalar-potential **BIE** implementation
(`calc_inductance.py --impedance-model esim`) and the **HCurl A +
Kelvin** implementation (`calc_fem_kelvin.py --impedance esim`) drive
the *same* 1-D ESIM cell solver and the *same* Karl loop, so their
end-to-end workpiece power `P_wp` must agree, and the design-relevant
per-element/uniform ratio must be reproduced across both.

Backs the validation tier **(vi)** of two manuscripts (materials live
outside the repo, under `W:\02_学会資料`):
- SA-26-070 (IEEJ 静止器・回転機合同研究会 @ 八戸, 2026-08)
- IGTE 2026 selected-papers full paper (Kubota et al.)

## Comparison matrix

| Matrix | Geometry | Point | Meshes | Purpose |
|---|---|---|---|---|
| **A** | cylinder dia 20 x 20 mm | I=100 A, f=50 kHz (BH knee) | BIE & FEM share `ih_fem_kelvin_demo.vol` | formulation isolated, **mesh fixed** |
| **B** | cylinder dia 50 x 25 mm (SA-26-070 ch.5) | I=500 A, f=10 kHz (deep saturation) | BIE `ih_bem_sample_p1.vol` vs FEM `ih_fem_kelvin_50mm.vol` (**independent**) | formulation isolated, **mesh not fixed** (harder) |

Material (both sides): steel `em_sample_bh.txt` (CEFC 2020), sigma =
2e6 S/m, cell radius R = 5 mm.

## Golden result (measured 2026-07-18, radia 4.95.17)

| Matrix | variant | BIE P_wp [W] | FEM P_wp [W] | agreement |
|---|---|---|---|---|
| A | per-element | 1.3381 | 1.3117 | +2.0 % |
| A | uniform | 1.6377 | 1.4407 | +13.7 % |
| B | uniform | 162.36 | 161.09 | +0.8 % |
| B | per-element | 81.18 | 74.28 | +9.3 % |

Design-relevant ratio (B, per-element / uniform): **BIE -50.0 % vs
FEM -53.9 %** -> uniform-`Z_s` over-predicts `P_wp` ~2x in deep
saturation, reproduced by both independent implementations.

Interpretation: the matrix-A *uniform* +13.7 % gap (vs per-element
+2.0 %) is not a physics disagreement -- converged `Re(Z_s)` agrees
between the two; the gap is the surface-`|H_t|` **aggregation
convention** feeding the single scalar `Z_s`. Per-element removes that
convention ambiguity. Matrix B's +9 % per-element gap is dominated by
edge-adjacent local `|H_t|` extraction near the cylinder rim (a field
singularity), which the two formulations discretise differently on
their independent meshes; the *ratio* is robust regardless.

## Files

| File | Role | Cubit? |
|---|---|---|
| `make_meshes.py` | regenerate all three `.vol` from tracked generators | **yes** |
| `make_ih_fem_kelvin_50mm.py` | the 50 mm FEM-Kelvin mesh generator (Cubit-embedded) | **yes** |
| `run_cross_formulation.py` | run the 8 calc subprocesses + assert golden bands + write `results.json` | no |
| `results.json` | committed golden record (physics only; timings excluded per Benchmark Policy) | -- |

`.vol` files are gitignored and regenerated on demand. Meshes A and B(BIE)
come from tracked samples (`ih_fem_kelvin_demo.py`, `ih_bem_sample.jou`);
mesh B(FEM) from this lane's `make_ih_fem_kelvin_50mm.py`.

## Run

```powershell
python make_meshes.py            # once; requires Cubit license
python run_cross_formulation.py  # runs + asserts; exit 0 = all bands pass
```

`run_cross_formulation.py` does **not** import cubit (Cubit/NGSolve
separation): it shells out to the two production `calc_*.py`. It errors
with a clear message if a `.vol` is missing (run `make_meshes.py` first).

## Known caveat (excluded from the asserted set)

At the matrix-B deep-saturation point, running the **BIE per-element**
solver on the **FEM's** 50 mm surface mesh (i.e. the *same-mesh* pairing
instead of the independent BIE mesh) enters a shallow limit cycle:
~40 Karl iters, `dZ_max ~ 5e-2` on a few rim DOFs, formally unconverged
(the surface-averaged `<|Z_s|>`, `<H_t>` are stationary, `P_wp ~ 87 W`
quasi-steady). The FEM per-element solver converges on that same mesh in
~19 iters. The lane therefore uses the **independent** BIE mesh
(`ih_bem_sample_p1.vol`, converges in ~22 iters) for matrix B, and this
same-mesh limit cycle is documented but **not** asserted -- it records
that deep-saturation per-element Karl is sensitive to mesh quality /
rim-DOF `|H_t|` extraction, not that the method is wrong.

Per the Benchmark Policy, correctness/agreement assertions may run on any
host; wall-clock timings are **not** part of the committed record.
