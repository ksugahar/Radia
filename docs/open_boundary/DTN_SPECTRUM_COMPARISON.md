# Open-boundary methods, measured on the DtN spectrum

The quantitative companion to [`OPEN_BOUNDARY_MAP.md`](OPEN_BOUNDARY_MAP.md): every
open-boundary closure compared on **one yardstick — the per-multipole Dirichlet-to-Neumann
(DtN) spectral defect**

```
    d_n  =  | λ_h(n) − λ_exact(n) |  /  | λ_exact(n) |
```

across the three regimes, with the measured numbers reproduced in-repo by
[`examples/kelvin_transformation/DtN_spectrum/act7_22_dtn_spectrum_consolidated.py`](../../examples/kelvin_transformation/DtN_spectrum/act7_22_dtn_spectrum_consolidated.py)
(data: `act7_22_dtn_spectrum_consolidated.json`; figure: `…​.pdf`). It consolidates the
scattered per-regime measurements — static [`act7_21`](../../examples/kelvin_transformation/DtN_spectrum/act7_21_lowfreq_openbc_4way.py),
eddy [`act6_09`](../../examples/kelvin_transformation/DtN_spectrum/act6_09_cln_vs_pml.py),
high-frequency `act7_01`/`act7_07` — into one table.

## Why the DtN-spectral defect (the lens)

A **field-error** comparison conflates the *interior FEM error* with the *open-boundary
error*. The **per-multipole DtN defect isolates the boundary operator's accuracy, mode by
mode** — a diagnostic the usual comparisons cannot give. The exact eigenvalue is a property
of the *continuous* exterior operator `Λ_ext`, regime by regime:

| regime | exact `λ_exact(n)` | operator |
|---|---|---|
| static (Laplace) | `−(n+1)/a` (3D), `−n/a` (2D) | real ladder |
| eddy (diffusion) | rational in `q=√s` (reverse-Bessel) | `radia.open_boundary.eddy_dtn` |
| high-freq (Helmholtz) | `z h_n^(1)′(z)/h_n^(1)(z)`, `z=ka` | `radia.open_boundary.wave_dtn` (complex/radiating) |

## The measured comparison (a = 1)

| method | accuracy (per-mode `d_n`) | convergent? | parameter-free? | DC conditioning | cost |
|---|---|---|---|---|---|
| **Kelvin** | **exact** (static `~2e−6`, eddy converges) | **YES** | **YES** | **flat** | sparse |
| **BEM** | **exact** (all regimes) | **YES** | **YES** | — | **DENSE** (`N²`) |
| **PML** | accurate per-mode (`~1e−4` in its home) | no (a tuned **layer**) | no (`σ, L`) | **blows up** (`2.4e4`@DC) | sparse |
| **CFS-PML** | accurate per-mode | no (a tuned layer) | no (`σ, α, L`) | **fixed** (`2.3e3`@DC) | sparse |
| **ballooning** | **fails the LOW modes** (`~(a/R)^{2n+1}`) | no (finite **reach** `R`) | `R` | — | sparse |
| **Robin** (`λ=−1/a`) | **exact `n=0` only**, fails HIGH modes | no (a fixed floor) | YES | — | sparse |

### Per-regime numbers (from `act7_22`)

- **static** (`λ_n=−(n+1)`): Kelvin `≤ 2.1e−6` every mode (converged); **ballooning** `n=0..4 =
  0.33 → 0.024 → 1.6e−3 → 1.1e−4 → 6.9e−6` (fails the *low* slow-decaying modes, shrinks with
  `n` and with `R`); **Robin** `0 → 0.5 → 0.67 → 0.75 → 0.80` (exact `n=0`, fails the *high*
  modes — the opposite failure).
- **eddy** (`s=i·1`): all three closures resolve the mode (`d_n < 5e−2`); the **distinguishing**
  behaviour is **convergence** — Kelvin-built `n=2` defect drops `1.2e−3 → 2.2e−5` under
  `(h, R_mid)` refinement (parameter-free), while PML/CFS-PML are a tuned absorbing layer; and
  **conditioning** — vanilla PML's interior matrix conditioning **blows up toward DC**
  (`cond ≈ 2.4e4` at `s=i·1e−4`) while **CFS-PML removes it** (`cond ≈ 2.3e3`) — *which is why
  CFS-PML exists*.
- **high-freq** (`z=ka=2`, radiating): **PML is accurate in its home** (`d_n ~ 1e−4`); **Kelvin
  exits** — its real-axis DtN cannot carry the radiation `Im` part, and this regime is **outside
  radia's MQS / Laplace-kernel scope** (→ NGSolve).

## The headline (the two classes, measured)

- **Convergent + parameter-free**: **Kelvin** (static / eddy) and **BEM** (all regimes, but
  DENSE) — the discrete operator `S_h → Λ_ext` for every mode under refinement.
- **Fixed-error surrogates**: **PML** (accurate per-mode but **DC-ill-conditioned** + tuned),
  **CFS-PML** (fixes the conditioning at modest accuracy + tuned), **Robin** (exact `n=0` only).
- **Finite-reach**: **ballooning** (cheap, but a finite wall — fails the low modes; shrinks with
  `R`).

> The comparison is **not a single number**: on the *accuracy* axis Kelvin, BEM and PML are all
> good per-mode; the methods separate on **convergence** (Kelvin/BEM vs tuned layer),
> **conditioning** (PML blows up at DC; CFS-PML fixes it; Kelvin flat), **cost** (BEM dense vs
> the rest sparse), and **reach** (ballooning finite). Reporting all axes is the honest result —
> not "method X is best".

## Provenance (cite — this is a measurement / comparison, not a new method)

Kelvin / inversion open boundary: Freeman-Lowther 1988/89; Brunotte-Meunier-Imhoff 1992.
Ballooning / infinite elements: Silvester-Hsieh 1971; Bettess. PML: Bérenger 1994;
Collino-Monk. CFS-PML: Kuzuoğlu-Mittra 1996. BEM-FEM open boundary: standard. The contribution
here is the **unified method × regime × multipole DtN-spectral comparison on one yardstick**.
