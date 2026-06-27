# `radia.open_boundary` — exact Zs-DtN-CLN open boundary

Production-API examples for the adopted core module
[`radia.open_boundary`](../../src/radia/open_boundary/dtn_cln.py): the **exact
exterior Dirichlet-to-Neumann (DtN) open boundary realised as a Cauer Ladder
Network (CLN)**, for a **separable (spherical) magneto-quasistatic (MQS)**
truncation.

## What it provides

For multipole `n` at a sphere of radius `R0`:

- `eddy_dtn(n, s, R0, mu_sigma)` / `wave_dtn(l, z)` — the **exact** reverse-Bessel
  DtN symbols (diffusion in `q=√s`, wave in `s`; same poles `roots(θ_n)`).
- `cauer_ladder(n)` / `eval_ladder(...)` — the CLN realisation, **exact at `n+1`
  stages**, well-conditioned, passive.
- `companion_poles(n)` — the auxiliary-ODE rates (`roots(θ_n)`, all `Re<0`) for a
  **transient passive Robin** open boundary (Grote-Keller form).
- `sqrt_s_passive_ladder(...)` — a finite passive realisation of the `√s`
  diffusion-memory element.

The **`kelvin_dtn`** companion **BUILDS** the DtN by a Kelvin-FEM (so it carries an
arbitrary / non-separable shape and a **material** exterior, which the closed form
cannot):

- `kelvin_fem_radial_dtn(n, s, ...)` — pure numpy; reproduces the closed-form
  `eddy_dtn` with **no DC floor** (the "Kelvin builds the exact DtN" check).
- `kelvin_dtn_matrix(mesh, order, s, nu=, sigma=)` + `steklov_spectrum(S, Mg)` —
  NGSolve; the arbitrary-shape / iron-exterior DtN matrix and its Steklov ladder
  (point-group split: cube `O_h`, square `C4v`). `band_cln_fit(...)` reduces it.

> **Honest provenance (3 layers):** Kelvin-with-air-exterior = Freeman-Lowther
> 1988/89 (FEMM); transforming `σ` under a conformal map = transformation optics /
> Ward-Pendry 1996 (DC-`σ` cloaks) — both **classical**. Sugahara's own *validated
> fusions* are the contributions: (i) Kelvin inversion as an exact open boundary with
> the **`σ`-conformal** transform so a **conductor crosses the truncation**,
> ECT-validated (IEEE Magnetics 2022 — the basis for the `(a/r)⁴ σ` / `(a/r)² μ`
> weights used here); (ii) the Kelvin material-aware DtN as an inverse-design kernel
> (SF-with-iron). Full record: MCP `kelvin_transformation(topic="material_exterior")`.

## Files

| file | what |
|---|---|
| `demo_dtn_cln_usage.py` | minimal production-API walkthrough; each printed `ok` is asserted |

Run: `python examples/open_boundary/demo_dtn_cln_usage.py` (numpy/scipy only — no
NGSolve / Cubit needed). Expected: the Cauer ladder reproduces the exact DtN to
`~1e-16` at exactly `n+1` stages for `n=1..6`, and all companion poles have `Re<0`.

## Scope (read this — where this is the right tool)

It **wins over a CFS-PML only on its island**: compact / quasi-spherical MQS
problems where you want an **exact, DC-well-conditioned** open boundary **and a
compact passive circuit / ROM**. It does **not** win in general — an elongated /
arbitrary truncation wastes mesh on the spherical shell (a box CFS-PML hugs
better), and genuine wave **radiation** (finite real `kR`) is outside radia's MQS
scope. **Not novel** (Grote-Keller / Hagstrom-Warburton continued-fraction ABCs;
Warburg→Cauer; Kameari CLN) — a verified reusable operator, not a paper claim.

See the selector + the three-sense "superior?" discussion in
[`docs/open_boundary/OPEN_BOUNDARY_MAP.md`](../../docs/open_boundary/OPEN_BOUNDARY_MAP.md),
and the research-stage corpus (non-separable Kelvin-built DtN, reflection FETD,
CLN-vs-PML benchmark) in
[`docs/kelvin/kelvin_dtn_spectrum_archive.ipynb`](../../docs/kelvin/kelvin_dtn_spectrum_archive.ipynb)
(`act6_02`, `act6_09`, `act6_11`, …).
