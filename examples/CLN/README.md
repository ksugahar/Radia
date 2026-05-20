# CLN — Cauer Ladder Network for eddy-current modeling

Working folder for the **Cauer Ladder Network (CLN)** research line: extracting
RL-ladder equivalent circuits of eddy-current decay in 3D conductors via
double-double (DD) verified arithmetic, BEM-Foster decomposition, and the
Hankel-Padé / Stieltjes / Wheeler family of moment-matching methods.

This example folder is the canonical working location (since 2026-05-12) for:

- **IGTE 2026 Symposium digest** (`igte_symposium_2026.tex` + `igtesymp.cls`) —
  sphere-only paper showing how DD GPU pipeline pushes the Cauer-extraction
  precision wall from FP64 stage 4-5 to DD stage 12+.
- **Q&A development log** (`MEMORY.md`, 600+ lines) — running record of design
  decisions during digest preparation, kept in sync with the figures and
  numerical tables.
- **Cross-references** to the radia-mcp tool `cln_sphere_dd_pipeline`
  (under `packages/radia-mcp/.../radia_ngsolve/knowledge/cln_sphere_dd.py`)
  which exposes the same pipeline as an LLM-callable tool.

## Files

| File | Purpose |
|---|---|
| `igte_symposium_2026.tex` | IGTE 2026 digest LaTeX source (2 pages, igtesymp class) |
| `igte_symposium_2026.pdf` | Compiled digest (current state) |
| `igtesymp.cls` | IGTE Symposium 2026 LaTeX class file (provided by organisers) |
| `MEMORY.md` | Q&A development log + design decisions |

## How it relates to the rest of the repo

- **`packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/cln_sphere_dd.py`**
  — markdown doc describing the DD pipeline; surfaced via the `cln_sphere_dd_pipeline` MCP tool.
- **`src/ext/axifemm/`** — C++ Henrotte axisymmetric FE (Q1/Q2/P1/P2) used for
  the FEM cross-check (axifemm_p2_triangle Phase B2; commit `81f6415f`).
- **Other CLN literature**: Cauer 1958, Henrici 1958 QD-Padé, Sugahara TEAM 28
  axisymmetric matlab, Stoll Bessel ground truth, Hiruma 3-term FEM-CLN.

## Status

- IGTE digest: sphere-only pivot complete (2026-05-12 long discussion with
  Prof. Nagamine). Cuboid Outlook section dropped in favour of pure
  sphere + DD-pipeline demonstration.
- Pending: reflect DD 540-cell results (`stage 0-5 reliable`, +1 stage over 270-cell baseline) into the digest tables.

## Schur-F Method (Cauer + SIBC asymptote, 2026-05-18)

A separate paper line in `scripts/` extends CLN by a single non-rational
`sqrt(s)` Schur-complement block, structurally preserving both DC
accuracy and the SIBC asymptote `Y(s) -> K_SIBC / sqrt(s)` at all `s`.
Companion manuscript: Sugahara, Nagamine, Hane,
*"Universal Cauer-SIBC Composition via Schur Complement,"*
submitted to IEEE Transactions on Magnetics, 2026.

The augmentation is
```
K_aug(s) = diag(K_r(s), z(s)),    z(s) = (s + d) / (K_SIBC sqrt(s))
```
with `K_SIBC = S sqrt(sigma/mu)` (`S` = surface area), reading directly
off CAD metadata.  The polarisation-resolved vector extension
(Theorem 1 of the journal manuscript) yields the triple
`K_SIBC^(alpha)` whose sum equals `2 K_SIBC^scalar` (each face is
counted twice, once per parallel polarisation).

| Script | What it verifies |
|---|---|
| `scripts/slab_cpe_schur.wls` | 1D slab: A-axis wall-band peak rel-err |
| `scripts/slab_asymptote_compare.wls` | 1D slab: C-axis ratio r(f) sweep to 1e12 Hz |
| `scripts/rect2D_schur_F.wls` | 2D rectangular prism: Foster + Schur-F, Nagamine cross-check |
| `scripts/cuboid3D_schur_F.wls` | 3D cuboid scalar: A1 plate, full 3-axis evaluation |
| `scripts/cuboid3D_schur_F_both.wls` | 3D cuboid scalar: both A1 + 5x2x1 brick |
| `scripts/cuboid3D_AAA_CFschur.wls` | 3D cuboid: AAA-based CF-Schur revenge (negative result) |
| `scripts/cuboid3D_CFschur_both.wls` | 3D cuboid: continued-fraction SIBC rationalisation (fails on all 3 axes) |
| `scripts/cuboid3D_vector_decomposition.wls` | 3D cuboid vector: Theorem 1 polarisation-resolved K_SIBC trace identity |
| `scripts/cuboid3D_quadrupole.wls` | 2D cuboid multipole (long cuboid): Theorem 3 quadrupole excitation u_ext=(y²-x²)/2, K_SIBC = √(σ/μ) ||u_ext||²_L²(∂Ω) |
| `scripts/cuboid3D_quadrupole_3D.wls` | 3D cuboid multipole: Theorem 3 in genuine 3D with u_ext=xy tesseral harmonic on 5×2×1 mm cuboid, Foster Mmax=99 (~490k modes) |
| `scripts/bem_cln_2cylinder.wls` | 2-cylinder BEM-CLN assembly: Schur-F per-element building blocks + phenomenological coupling, preserves SIBC tail to 5 sig fig at f=10¹² Hz |
| `scripts/bem_cln_Ncylinder.wls` | N-cylinder linear-chain sparse-block BEM-CLN: N=2,3,5,10 scaling demo, DOF = N(N_Cauer+1), Schur-F vs Bessel-exact building blocks (Phase 1/2, phenomenological coupling — superseded) |
| `scripts/bem_cln_2cylinder_rigorous.wls` | **Rigorous BEM-CLN** via polarizability α(s) = a² - Y_cyl(s)/(πσ) and pure 1/D² 2D dipole-dipole coupling. DC limit α(0)=0 + PEC limit α(∞)=a² built in; no saturation factor needed. Schur-F vs Bessel agreement: 0.998-0.999 across N (Phase 2.5) |
| `scripts/bem_cln_2cuboid_3D.wls` | **3D extension** of BEM-CLN: 2-cuboid Mathematica framework with polarizability α(s) = V - Y_cuboid(s)/σ and 3D dipole-dipole coupling μ₀/(4πD³). Sign opposite to 2D (destructive 3D dipole). Coupling magnitude ~3×10⁻¹⁰ at typical separations (Phase 3 B) |
| `scripts/bem_cln_ngsolve_2cuboid.py` | NGSolve A-formulation FEM cross-check (Phase 3 B.2): 2 Cu cuboids + Kelvin sphere, mesh ~83k elements, eddy current solve at 10⁵-10⁷ Hz. Reveals demagnetization-factor gap with scalar-diffusion α (~O(10) discrepancy), confirming need for proper 3D polarizability tensor |
| `scripts/bem_cln_5cylinder_coil.wls` | **Engineering application** (Phase 3 B.3): 5-cylinder coil-like configuration via rigorous 2D Phase 2.5 framework. 25 DOF total (40× reduction vs 1000-mode Foster), +28.8% coupling enhancement at 10⁶ Hz, Schur-F preserves SIBC tail to 10 sig fig at 10¹² Hz |

See [docs/cln/BEM_CLN.md](../../docs/cln/BEM_CLN.md) for the full
BEM-CLN theory and the radia-mcp `bem_cln` tool for an MCP-callable
documentation interface.

The Cauer-SIBC composition is the central method of the IEEE Trans Mag
submission; the scripts above reproduce every numerical claim in the
manuscript Tables 1-5 and Figures 1-2.

## Build the digest

```bash
cd examples/CLN
pdflatex igte_symposium_2026.tex   # or via texcompile skill
```
