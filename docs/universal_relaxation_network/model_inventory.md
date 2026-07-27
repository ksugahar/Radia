# URN Candidate Model Inventory

This note inventories established impedance, relaxation, resonance, and passive
macromodel families that can inform the Y-domain URN dictionary.  The immediate
motivation is the SA/RM-2026 PCB inductor and NL87 ferrite data: the current
22-basis dictionary can select interpretable mechanisms, but removing frequency
attention degrades the fit, especially around self-resonance.

## Selection Rules

Prefer a candidate basis when it satisfies most of the following:

- It corresponds to a known passive circuit, material relaxation, or transport
  process.
- It can be written as a positive-real admittance or impedance contribution.
- Its parameters have physical meaning: time constant, resonant frequency,
  quality factor, diffusion length, fractional exponent, or asymptotic loss.
- It remains usable for time-domain realization after attention is removed.

## High-Priority Additions

| Candidate | Representative form | Why it matters | URN action |
|---|---|---|---|
| Parallel RLC tank | `Y/G = 1 + Q*(s/w0 + w0/s)` | Produces an impedance maximum at resonance. This is the natural passive shape for inductor self-resonance, unlike a series-RLC admittance peak. | Add `parallel_rlc` to Y-domain URN. Use for PCB self-resonance and high-frequency NL87 roll-off. |
| Multi-resonance Foster branches | `Y = sum_k G_k*(1 + Q_k*(s/w0_k + w0_k/s))` plus RL/RC terms | Captures multiple self-resonances while staying circuit-realizable. | Allow several `parallel_rlc` bases and prune by output ablation. |
| Skin/proximity RL ladder | finite R-L ladder or `Z_series = R0 + K*sqrt(s) + s*L` approximation | Winding resistance rises and internal inductance changes with frequency. Established compact circuits use frequency-independent R-L ladders. | Add either a compact `skin_rl_ladder` macro-basis or a positive-real `sqrt_s_series_impedance` branch. |
| 2-Pi spiral/PCB inductor model | layout-style wideband RLC network with substrate/parasitic capacitance | Known wideband inductor models capture `R(f)`, `L(f)`, parasitic capacitance, and self-resonance with transient-compatible elements. | Add a higher-level "inductor macro" option outside pure basis summation, or compile it into RL/RLC branches. |
| Ferrite magnetic resonance | Lorentz/LLG susceptibility, e.g. `chi(s) = A*w0^2/(s^2 + 2*zeta*w0*s + w0^2)` | Ferrites can show domain-wall resonance and natural ferromagnetic resonance, not only Debye-like relaxation. | Add `magnetic_lorentz` / `magnetic_llg` bases in permeability space, then map to terminal impedance. |
| Magnetoelastic / magnetostrictive resonance | BVD-style motional branch coupled to magnetic port | Magnetostriction couples magnetic excitation to mechanical modes, producing geometry-dependent resonance/antiresonance. | Add `magnetoelastic_bvd` as an optional resonant branch when sharp mechanical resonances are visible. |
| Passive rational macromodel | vector fitting/Loewner/AAA poles with passivity enforcement | Best accuracy and time-domain simulation support, but less physical by default. | Use as benchmark and fallback. Optionally map stable pole pairs into Foster/RLC branches before URN pruning. |

## Relaxation And Material Models

| Family | Representative form | Current status | Notes |
|---|---|---|---|
| Debye | `1/(1+s*tau)` and magnetic duals | Present | Good for a single relaxation time. |
| Cole-Cole | `1/(1+(s*tau)^alpha)` | Present | Broad symmetric distribution of relaxation times. |
| Cole-Davidson | `1/(1+s*tau)^beta` | Missing in Y-domain variant | Useful for asymmetric relaxation peaks; already appears in broader Radia Z-domain URN. |
| Havriliak-Negami | `1/(1+(s*tau)^alpha)^beta` | Missing in Y-domain variant | Generalizes Cole-Cole and Cole-Davidson; strong candidate for ferrite and dielectric broadening. |
| Jonscher universal dielectric response | `Y proportional to s^alpha` or conductivity power law | Partly present via CPE | Keep CPE, but consider explicit low/high asymptotic power-law blends. |
| Lorentz / LLG magnetic resonance | oscillator-like permeability or susceptibility | Missing | Useful when ferrite shows resonant rather than pure relaxation behavior. Needs passivity-constrained oscillator form. |
| Domain-wall resonance | damped oscillator contribution to permeability | Missing | Low-to-mid RF ferrite permeability can include domain-wall motion modes. |
| Magnetoelastic resonance | BVD/motional RLC branch or coupled magnetic-mechanical oscillator | Missing | Captures magnetostriction-driven mechanical resonance. Use only when geometry supports an acoustic mode in band. |

## Electrochemical And Diffusion Models

These are less central for PCB/NL87, but important if URN remains a general
impedance-identification tool.

| Family | Representative form | URN action |
|---|---|---|
| Constant phase element (CPE) | `Y = K*s^beta` or `Z = K/s^beta` | Present as inductive/capacitive CPE. |
| Warburg semi-infinite diffusion | `Z proportional to 1/sqrt(s)` | Present in older Z-domain URN; consider Y-domain finite-length variants only when diffusion data are targeted. |
| Finite-length Warburg open/short | `tanh(sqrt(s)*B)/sqrt(s)` or `coth(sqrt(s)*B)/sqrt(s)` | Candidate for batteries/electrochemistry; not urgent for PCB/NL87. |
| Gerischer | diffusion plus first-order reaction, often `1/sqrt(k+s)` style | Candidate for solid-state electrochemistry; not urgent for PCB/NL87. |
| Randles/ZARC | `Rs + (Rct || CPE) + Warburg` variants | Useful as a macro-template rather than a single basis. |
| DRT | non-parametric distribution over relaxation times | Use as diagnostic to propose how many Debye/Cole/HN bases are needed. |

## What This Means For SA/RM-2026

The present 22-basis Y-domain dictionary is strongest as a mechanism selector,
not as a complete transient-ready equivalent-circuit generator.  The missing
piece for the PCB inductor is not just "more of the same" Debye/CPE bases.  The
dictionary needs at least one basis family whose passive time-domain realization
naturally creates an impedance peak: a parallel RLC tank or a Foster-style
anti-resonance branch.

For NL87, the active-only degradation is less diagnostic of a missing physical
family.  The 22-basis attention-free model recovers much of the fit, suggesting
that pruning may be too aggressive.  Havriliak-Negami or Cole-Davidson magnetic
relaxation bases are still attractive because ferrite permeability spectra often
show asymmetric broadening beyond Debye/Cole-Cole.

Ferrite-specific resonance should be treated as a separate candidate family.
If the measured band includes a genuine magnetic resonance, a damped Lorentz or
LLG-derived permeability term is more physical than approximating the feature
with several Debye/CPE bases.  If the peak position depends strongly on core
geometry or mechanical boundary conditions, a magnetostrictive/magnetoelastic
branch is also plausible: magnetostriction couples magnetic excitation to an
elastic mode, and the terminal impedance can show a motional resonance similar
to a Butterworth-van Dyke resonator.  These resonant bases are especially useful
for time-domain realization because they can be implemented as passive second-
order branches, provided damping and coupling coefficients are constrained
positive.

## Implementation Priority

1. Add `parallel_rlc` to the Y-domain URN and compare attention-free fits.
2. Add `magnetic_lorentz` / `magnetic_llg` and optional `magnetoelastic_bvd`
   bases for ferrite resonance and magnetostriction.
3. Add Y-domain Cole-Davidson and Havriliak-Negami electric/magnetic variants.
4. Add a compact skin/proximity branch for PCB winding loss.
5. Keep passive vector fitting as a reference model and optional route to a
   Foster/Cauer network.
6. Use DRT as a diagnostic layer to propose basis count before training.

## References

- Cole and Cole, "Dispersion and Absorption in Dielectrics I", J. Chem. Phys.,
  1941: https://doi.org/10.1063/1.1750906
- Havriliak and Negami, "A complex plane representation of dielectric and
  mechanical relaxation processes in some polymers", 1967:
  https://doi.org/10.1002/polc.5070140111
- SINTEF Vector Fitting overview:
  https://www.sintef.no/en/software/vector-fitting/
- Wang, "Distribution of relaxation times for impedance analysis in batteries
  and fuel cells", Nature Reviews Clean Technology, 2025:
  https://doi.org/10.1038/s44359-025-00071-z
- Kim and Neikirk, "Compact Skin Effect Circuit":
  https://www.weewave.mer.utexas.edu/MED_files/MED_research/Intrcncts/Skin_Effect_Ldr/MTT_96_skn_ldr.html
- Cao et al., "Frequency-independent equivalent-circuit model for on-chip
  spiral inductors", IEEE JSSC, 2003:
  https://doi.org/10.1109/JSSC.2002.808285
- Critical review on EIS analysis and Warburg/Randles models:
  https://doi.org/10.1063/5.0283768
- Analog Devices, "Ferrite Bead Demystified":
  https://www.analog.com/en/resources/app-notes/an-1368.html
- Perekalina et al., "Resonance of Domain Walls in Cobalt Ferrite", JETP,
  1961: https://jetp.ras.ru/cgi-bin/dn/e_013_02_0303.pdf
- Kutorasinski et al., "Nonlinear magnetic ring model based on impedance
  measurements with DC-bias current", Scientific Reports, 2026:
  https://www.nature.com/articles/s41598-026-39594-1
