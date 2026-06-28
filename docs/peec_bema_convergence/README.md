# PEEC vs BEM-A convergence study on 3turnCoil

Results from the 2026-05-13 convergence sweep on
`src/radia/panels/samples/3turnCoil_work.jou`, comparing the PEEC
n_peri-mesh sweep against the BEM-A surface-mesh sweep.  Recorded in
[memory project_peec_v4_44_0_lead_aware_chain] and
[memory project_3turncoil_baseline_not_analytical].

## Files

| File | Contents |
|---|---|
| `3turnCoil_peec_bema_convergence_2026_05_13.json` | Per-step L_coil values: PEEC at n_peri=8/16/32/64 vs BEM-A at 4028/5938/9848/14972 tris.  `frequency=50000`, Cu wire.  L converges monotonically for PEEC (423.89 → 431.55 nH, L_∞ ≈ 434 nH); BEM-A is non-monotonic (432.57 → 417.56 → 412.18 nH) and 14,972 tris is unreachable on this machine (dense memory + MINRES stall).  Production retains n_peri=16 / L = 426.30 nH, ~1 s runtime. |
| `peec_bema_convergence.ipynb` | Result-bearing public view with the comparison table and convergence plot. |
| `peec_bema_convergence_results.json` | Versioned durable result JSON synchronized with the notebook. |

## Key finding

PEEC and BEM-A agree to 5-7% on this 3-turn coil (peer-validated;
neither is "ground truth" since no analytical solution exists for a
3-turn pancake with lead bars).  The previously-cited
"analytic-equivalent reference 426.25 nH" was inaccurate -- that
number is the legacy `.jou-path PEEC` value (released v4.7.0, dropped
in v4.13.0).  True independent cross-check requires HACApK BEM-A at
>10k tris OR FEM-Kelvin volumetric solve.

See `src/radia/panels/samples/ih_fem_kelvin_demo.py` (canonical
PEEC-FEM-Kelvin demo, v4.45.0) for the FEM-Kelvin path now that the
launcher-free `radia.kelvin_identify_ngsolve` helper is shipped.
