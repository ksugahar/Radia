# B-input Play vs Energy: Element-level comparison

This example demonstrates the differences between the two B-input hysteresis
material types in Radia, both built from the same Play shape functions
(Sugahara/Hane B-input identification):

| Material | Class | Inverse algorithm |
|----------|-------|-------------------|
| Type 6 | `radTPlayHysteresisMaterial` (`MatPlayHysteresis`) | Newton on B (state-dependent Jacobian) + bisection fallback |
| Type 5 | `radTEnergyHysteresisMaterial` (`MatEnergyHysteresis`) | K independent per-particle Newtons (Egger framework) |

**Forward (B → H) is algebraically identical between the two** — both compute
`H = sum_k f_k(p_k(B))` from the same shape functions. The two models
differ only in the inverse path and in the energy decomposition (Type 5
exposes per-particle polarisations J_k; Type 6 doesn't). This example
verifies the forward equivalence, then exhibits the inverse differences.

## What's in this directory

- `compare_play_vs_energy.py`     — Radia C++ comparison: drive same H(t)
  through both materials, compare B(t), iteration counts, wallclock per step.
- `compare_play_vs_energy.m`      — MATLAB companion using the forensic
  reference implementations under `w:\02_学会資料\2026年度\2026_09_IGTE_Symposium\菅原\matlab\+bie\`.
- `figures/` — generated plots (BH loops, per-step iter counts, error
  vs MATLAB cross-validation).

## Source data

The Potter-Schmulian B-input.mat (analytical Jiles-Atherton-style model,
20 BMax levels up to 1.9 T) is used so anyone can reproduce the result
deterministically. Path: `W:\999_菅原賢悟\19_磁気ヒステリシス\2024_IGTE_共同研究\2024_03_08_H-input_B-input\Potter_Schmulian\B_input.mat`.

## Background

This work was done in preparation for the IGTE'26 submission on the
B-input energy-based vector hysteresis model. The cross-validation against
the MATLAB implementations exposed and fixed 11 production bugs in
`src/radia/hysteresis_io.py`, `src/radia/energy_play_model.py`, and
`src/core/rad_material_impl.cpp` (commits efd83a7d through ab578758).
