# Stream Function Method (SFM) coil design — full documentation

This folder is the **canonical documentation** for the Radia stream-function
coil-design framework: kernel-agnostic (ACA+)+TSVD least-norm solver,
single-stroke chain construction (field_aware / Kuijpers), FE-direct ψ with
H1 / σ-weighted / inductance / L∞ regularisation, Path-A compensated
iteration, Optuna CMA-ES surface-deformation outer loop, **single-current
"sheet-metal" wire distortion** (one bent conductor instead of separate shim
feeds), **FE-direct ψ on arbitrary curved formers** (sphere / conformal /
3D-printed — the case the basis-loop representation cannot do), and
**folded-Tikhonov + ACA-TSVD Pareto optimisation of the (field homogeneity,
peak current density) front** with four stackable levers — Tikhonov α
(L-curve), L∞ minimax, geometry (former size / cylinder length), and
**sheet-metal surface forming** (genuine bending: planar out-of-surface −17 %,
cylinder in-surface −25 %).

The single-file [`../stream_function.md`](../stream_function.md) is the
short *entry point*; this folder is the detailed reference.

## Where to start

| You want to … | Read |
|----------------|------|
| Understand the SFM and ACA+TSVD math | [theory.md](theory.md) |
| **Run the design / pareto / manufacture GUI panel** (`calc_streamfunction.py`) | [**panel.md**](panel.md) |
| Pick the current-confinement BC (`--confine off/on/abe`, Abe edge-equipotential) | [panel.md § confinement](panel.md#current-confinement-boundary-condition---confine-off-on-abe) |
| Draw order-p contours / bubble-system flux lines (`--contour-sub`, `--flux-plot`) | [panel.md § contour=flux-line](panel.md#contour-drawing--flux-line-drawing) |
| Connect contours into one wire | [single_stroke.md](single_stroke.md) |
| Bend that one wire into a manufacturable single-current shim (no extra feeds) | [single_stroke.md § sheet-metal distortion](single_stroke.md#single-current-sheet-metal-coil-distortion-bankin-ho--no-extra-feeds) |
| Design on an arbitrary curved former (sphere / conformal) | [single_stroke.md § arbitrary curved formers](single_stroke.md#arbitrary-curved-formers-sphere--fe-direct-ψ-demo_sphere_fe_directpy) |
| Pick a regularisation for your problem | [regularization.md](regularization.md) |
| Push the **(homogeneity, peak current density)** Pareto front (Tikhonov α / L∞ / geometry / sheet-metal forming) | [regularization.md § Pushing the front](regularization.md#pushing-the-homogeneity-peak-j-pareto-front) |
| Optimise the coil SURFACE geometry (bilevel) | [deformation.md](deformation.md) |
| Look up the Python API | [api.md](api.md) |
| Reproduce a published benchmark | [benchmarks.md](benchmarks.md) |
| Run the public demo gallery | [demo_gallery.ipynb](demo_gallery.ipynb) |
| Audit / migrate the remaining example scripts | [examples_catalog.ipynb](examples_catalog.ipynb) |
| Hook ngsolve.bem H-matrix (2604+) | [ngsbem_integration.md](ngsbem_integration.md) |
| **Design a stellarator coil** (REGCOIL / NESCOIL / FOCUS: winding-surface current potential, net current, coil force/stress, VMEC boundary, winding-shape) | [**fusion.md**](fusion.md) |
| Cite / publish this work | [paper_outline.md](paper_outline.md) |

> **Three different "deformations" — do not confuse them.**
> (1) [deformation.md](deformation.md) reshapes the coil *surface* and
> re-solves ψ for ACCURACY (bilevel geometry optimisation, (RMS, energy)).
> (2) The **sheet-metal surface forming** in
> [regularization.md § Pushing the front](regularization.md#pushing-the-homogeneity-peak-j-pareto-front)
> also reshapes the *surface* + re-solves ψ, but to lower the **peak current
> density** on the (homogeneity, peak) front — with the honest
> standoff-vs-genuine-bending decomposition (`demo_pareto_deform.py` /
> `demo_pareto_cylinder_deform.py`).
> (3) The *sheet-metal wire* distortion in
> [single_stroke.md](single_stroke.md) keeps ψ + the contour levels fixed and
> bends the manufactured *wire* (one current) to cancel the single-stroke
> residual.  (1) and (2) are surface reshapes that re-solve ψ; (3) keeps ψ.

## Quick-start (5 lines)

```python
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel

obs = ...                         # (M, 3) observation points
sources = [radia_obj_1, ...]      # N Radia handles (coils OR magnets)

entry = radia_field_kernel(obs, sources, component=2)   # Bz kernel
res = aca_tsvd(len(obs), len(sources), entry, modes=20)
phi = pseudo_inverse_solve(res, B_target, k_mode=10)
```

For the human-facing path, start from the result-saved public notebooks instead
of running loose example scripts:

| Route | Artifact |
|-------|----------|
| Public demo gallery | [`demo_gallery.ipynb`](demo_gallery.ipynb), synchronized with [`demo_gallery_results.json`](demo_gallery_results.json) |
| Full source/result ledger | [`examples_catalog.ipynb`](examples_catalog.ipynb), synchronized with [`examples_catalog_results.json`](examples_catalog_results.json) |
| Theory and FE-direct psi | [`theory.ipynb`](theory.ipynb) |
| Regularization and Pareto trade-offs | [`regularization.ipynb`](regularization.ipynb) |
| Surface deformation search | [`deformation.ipynb`](deformation.ipynb) |
| Runnable validation/benchmarks | [`validation_test/stream_function/`](../../validation_test/stream_function/) |
| Reusable Stage-2 API/CLI | [`src/radia/panels/calc_streamfunction.py`](../../src/radia/panels/calc_streamfunction.py) |

The transitional demo source names are cataloged in
[`examples.md`](examples.md); new public links should point to the notebooks,
JSON sidecars, `validation_test`, or `src` API rather than the old examples
tree.

## What this is, what it is not

This framework is an **integration of three OSS libraries** —
`NGSolve` / `ngsolve.bem` (BEM operators + high-order FE basis),
`HACApK` (kernel-agnostic ACA+ H-matrix), and `Radia` (MMM/MSC
material kernels + chain construction) — plus SF-coil-design-specific
glue (contour extraction, single-stroke spiral, Path-A iteration,
Optuna CMA-ES bilevel deformation).

The **underlying numerical methods are standard**: ACA+ (Bebendorf
2000), TSVD recompression (Hackbusch 2008 *Hierarchical Matrices*),
Galerkin BEM, H¹ regularisation, CMA-ES.  We do not claim algorithmic
novelty for these.

What may be **new** (subject to literature search before publication):

  - The **Path-A compensated iteration** that folds Kuijpers 2023's
    "deviation = field error" observation back into the SF solve, and
    its **monotone convergence** on FE-direct ψ (vs oscillation on
    grid-sampled ψ).
  - The **integrated open-source pipeline** that lets the same SF
    design loop run on free-space, material (Radia MMM), or surface-
    BEM (ngsolve.bem) kernels by replacing the entry callback.
    See [paper_outline.md](paper_outline.md) "Implementation
    contribution" for the framing.

## Feature comparison (with related tools)

| Capability                       | Our state                   | Related tool / paper       |
| -------------------------------- | --------------------------- | -------------------------- |
| SF inverse design                | ACA+TSVD, kernel-agnostic   | CoilGen (MRI-only OSS)     |
| Single-stroke connection         | 3 methods, Kuijpers Method-1| Kuijpers 2023 (paper only) |
| Single-stroke compensation       | Path-A monotone (FE-direct) | Kuijpers 2023 observes; we iterate |
| Bilevel (inner SF, outer geom)   | Optuna CMA-ES               | Comsol Opt Module (commercial) |
| Regularisation choices           | 5 (L2/H1/σ/L∞/L_diag)       | Liu-Hennig-Korvink 2012 (H²) |
| High-order FE ψ                  | H¹ order p (any)            | NGSolve raw                |
| **Single-current sheet-metal distortion** | plane/cyl/sphere, ONE feed | geometric shim — no separate-feed analogue |
| **Arbitrary curved formers (FE-direct ψ)** | ANY surface (sphere/conformal) | CoilGen / basis-loop need a structured grid |
| Manufacturability gate (inter-turn spacing) | yes (≥ conductor width)     | —                          |
| Material-kernel extension        | callback ready (= TODO demo)| —                          |
| ngsolve.bem 6.2.2604+ alignment  | callback contract matches   | Joachim Schöberl + Pierre Marchand upstream |

## What's solidly proven vs what's pending validation

**SOLID** (numbers reproducible from the demos):

  - Planar uniform Bz (50 cm plane, target at 10 cm, B0=1 mT):
    - Basis-loop baseline:   RMS  2.99 %, p2p 9.59 %
    - + Path-A 100 iter:     RMS  0.58 %, p2p 2.17 %
    - FE-direct H1 baseline: RMS  2.09 %, p2p 6.81 %
    - + Path-A (MONOTONE):   RMS  0.47 %, p2p 1.64 %
    - + bump deform 20 trial: RMS 0.77 %, p2p 3.12 %
    - + order=3 sweet spot:  RMS **0.51 %**, p2p **1.83 %**
  - Cylindrical Gx fingerprint (Hard tier):
    - chain methods: **field_aware (9.3 %)** beats kuijpers (16.24 %)
    - + single-current **sheet-metal distort: 8.5 % → 1.4 %** (one bent wire,
      ~30 mm radial bend) — beats 10 separate electric shims (2.3 %)
    - + a few electric shims on the distorted residual: → **1.0 %** (11 feeds)
    - Path-A oscillates here (tier-bounded; FE-direct does NOT unstick the
      fingerprint topology — the win of FE-direct is generality, below)
  - Cylindrical Gz: smooth helix, trivial single-stroke, done.
  - **SPHERE former** (FE-direct ψ — ANY curved surface, basis-loop cannot):
    - uniform Bz: continuous cres 3e-15 (0 ppm), single-stroke **0.24 %**
    - Z2 shim: continuous cres 5e-14, single-stroke 4.3 % → sheet-metal
      distort **0.36 %** (one current, ~2 mm bend); real inter-turn spacing
      **10.5 mm** ≥ conductor width (manufacturable single wire)
  - HACApK ACA+TSVD on FE-direct matrix: validated equivalent to lstsq.

**PENDING VALIDATION** (industry benchmarks, see `benchmarks.md`):

  - Bilac et al. planar shim — stub
  - Turner cylindrical Gz analytical — stub
  - Lemdiasov-Ludwig 2005 target field — stub
  - CoilGen OSS head-to-head — stub
  - Shielded coil with Radia MMM iron yoke kernel — stub (material-kernel demo via the callback contract)

## Cross-references

  - MCP knowledge: ``streamfunction(topic=session_2026_05_30)``,
    ``streamfunction(topic=single_stroke)`` for the full session-log narrative,
    and ``streamfunction(topic=regularized)`` for the regularisation-folded
    closed form `ψ = S⁻¹V · W⁻¹ · Σ⁻¹ · UᵀB`.
  - Memory entries (LAB-private):
    - ``feedback_single_stroke_chain_orientation_traps``
    - ``feedback_path_a_naive_picard_negative``
    - ``feedback_fmm_vs_aca_distinction``
    - ``feedback_gmsh_gui_invisible_from_background``
    - ``project_fe_direct_psi_path_a_converges`` (FE-direct on plane / cylinder
      / **sphere**; the arbitrary-surface generalisation)
    - ``project_coil_distortion_sheet_metal_2026_05_31`` (single-current
      sheet-metal wire distortion: plane / cylinder / sphere)
