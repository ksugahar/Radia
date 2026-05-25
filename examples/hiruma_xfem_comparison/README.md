# Hiruma 2023 XFEM convergence study — 4-phase comparison

A research-oriented comparison of standard CFEM, Hiruma's PUM/XFEM
enrichment, and Schur-complement augmented CLN model order reduction
on the canonical skin-effect cylinder benchmark of Hiruma 2023.

Reference: Hiruma 2023, IEEE TMag 59(5), DOI 10.1109/TMAG.2023.3246629,
§III.A "Skin Effect in Cylindrical Conductor", Figs. 2-4.

## Problem statement

2D disk Omega, radius r = 0.01 m, sigma = 1e6 S/m, mu_r = 1, uniform
current density J0 = sigma (so E_z = 1 V/m on the boundary).
Governing PDE:

    -div(nu * grad A_z) + j * omega * sigma * A_z = J0      in Omega
                              A_z = 0                      on dOmega

The solution has a `sqrt(s)` corner-singular character (skin-effect
boundary layer) which makes polynomial CFEM converge slowly at
high r/delta.

## Phases

### Phase 1 — `phase1_cfem_cylinder_baseline.py`

Plain CFEM baseline on coarse (~37 DOF, matches Hiruma Fig. 2(b))
and fine reference (~9489 DOF) meshes.  Validates against the
analytical Bessel solution and quantifies the polynomial-FE
convergence floor when the solution has `sqrt(s)` corner
singularity.  Results: `results_phase1.json`.

### Phase 2 — `phase2_xfem_hiruma_enrichment.py`

Hiruma's PUM/XFEM with `psi(x) = exp(-gamma * xi(x))` enrichment
(xi = signed distance from surface).  Compound H1 FE space
`fes_std x fes_enr`.  Recovers near-h^p convergence on the same
coarse mesh, matching Hiruma Fig. 4.

### Phase 3 — `phase3_sqrt_s_schur_comparison.py`

`sqrt(s)` Schur-complement augmented CLN MOR.  Builds the 4-way
comparison table:

| Method                              |  DOF | Error at r/delta=15 |
|-------------------------------------|------|---------------------|
| CFEM fine                           | ~3.5k| <0.2% (reference)   |
| CFEM coarse                         |   44 | -11%                |
| Hiruma XFEM                         |   88 | -2%                 |
| Schur-augmented CLN                 |    5 | (target: -1%)       |

`Y_R(s) = Y_CLN(s) + K_SIBC * sqrt(s) / (s + d)` where
`K_SIBC = 2 pi R sqrt(sigma / mu)`.

### Phase 3b — `phase3b_krylov_galerkin.py`

Krylov-Galerkin (Pade-free) augmented CLN, replacing the
Hankel-ill-conditioned Pade extraction with the canonical CLN
basis `q_k = K_0^{-1} (sigma * M * q_{k-1})`.  Avoids the
high-N float64 breakdown of moment-based Pade.

### Phase 4 — `phase4_xfem_hankel_conditioning.py`

XFEM-CLN Hankel conditioning experiment.  Falsifies the hypothesis
"XFEM enrichment improves Hankel conditioning of the canonical CLN
moment matrix at higher MOR stages N".  Outputs results JSON +
PNG/PDF (`phase4_xfem_hankel_conditioning.{png,pdf}`) showing
`kappa(H_N)` vs N for CFEM vs XFEM bases.

## Running

```bash
cd S:/Radia/01_GitHub/examples/hiruma_xfem_comparison
python phase1_cfem_cylinder_baseline.py
python phase2_xfem_hiruma_enrichment.py
python phase3_sqrt_s_schur_comparison.py
python phase3b_krylov_galerkin.py
python phase4_xfem_hankel_conditioning.py
```

Each writes `results_phaseN.json` next to the script.

Requirements: NGSolve >= 6.2.2603, numpy, scipy, matplotlib (Phase 4).

## Companion knowledge

- `radia_mcp.fem.xfem_em_hiruma_knowledge(topic=...)` — full
  formulation derivation + the Schur sqrt(s) augmentation theory
- `radia_mcp.mor.cln_knowledge(topic=...)` — Cauer/PRIMA MOR
  background (Hiruma's reduced model is a CLN variant)
