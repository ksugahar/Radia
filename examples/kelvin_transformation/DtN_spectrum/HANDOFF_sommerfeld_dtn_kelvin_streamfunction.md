# HANDOFF — Sommerfeld-type DtN-Kelvin stream-function (material-aware coil design)

> **For a COLD Claude Code session in another lab.** You have NO memory of the originating
> conversation. This document is self-contained. Read it, then read the two cross-referenced files
> (`PATHWAY_streamfunction_with_iron.md` here, and the `dtn_coarse_mesh` knowledge module) and you can
> continue. Everything here is **open math** (NGSolve / Kelvin transform / stream-function / open
> boundary) and lives in the **public** Radia repo — no commercial-tool content is involved.

## TL;DR — what you are continuing
Fuse **two modules that already exist in this repo** to design coils **with magnetic material** (iron
yoke / shield / core), which the standard free-space stream-function method cannot do cleanly:
1. **Free-space stream-function coil design** (existing): `radia.streamfunction` / `calc_streamfunction.py`
   / `calc_streamfunction_volume.py`, knowledge `radia_mcp.streamfunction` — ACA-TSVD current-potential
   design (Koiso/Sugahara/Ida). Kernel = **free-space Biot-Savart**.
2. **DtN / FEM-Kelvin core operator** (existing, research-stage): `examples/kelvin_transformation/
   DtN_spectrum/`, knowledge `radia_mcp.radia_ngsolve.dtn_coarse_mesh` (call the MCP tool
   `dtn_coarse_mesh(topic="formulation")` to load the full design log). Generates the exterior
   Dirichlet-to-Neumann / transfer operator **sparsely, Green-function-free, for arbitrary mu(x)**.

**The fusion:** the matrix a stream-function design inverts (psi -> field) becomes, with iron, the
system's *material* Green operator. Build it as the Schur complement of the Kelvin-FEM (which meshes the
iron and carries it as an FE coefficient) -> a **material-aware transfer/DtN matrix M**. Then coil design
is the same clean linear inverse `psi = M^+ B_target` as in free space, but the kernel now contains the
iron. (This is why "the stream-function method WANTS the DtN'd matrix" — the inverse design consumes M.)

## Why free-space fails (the problem you are solving)
Free space: psi->field kernel = Biot-Savart (analytic, easy). With iron: total field = coil field +
**iron reaction**; the kernel becomes the material Green operator. Planar/cylindrical iron -> the (hard)
layered/**Sommerfeld** Green function; **arbitrary** iron -> NO closed-form Green function (a volume
integral equation revives the dense volume unknown). The clean "psi x kernel" structure is lost. The
Kelvin-FEM restores it: it meshes any iron geometry and condenses to M, no Green function.

## What is ALREADY done (verified, committed; run from `examples/kelvin_transformation/DtN_spectrum/`)
See `PATHWAY_streamfunction_with_iron.md` for the full table. The load-bearing demos for THIS task:
- `demo_v_assemble_dtn_matrix.py` — assemble the material-loaded exterior DtN matrix M = Schur complement
  of the sparse Kelvin-FEM; spectrum matches analytic.
- `demo_bb_nonlayered_inclusion_dtn.py` — the SAME for a NON-layered (arbitrary) iron blob, verified by
  the reduced-symmetry (C∞v |m|) splitting — i.e. M is correct for arbitrary iron where no Sommerfeld
  Green function exists.
- `demo_ee_streamfunction_coil_with_iron.py` — coil mode + iron shell: the free-space kernel is off by up
  to ~16x; the material-aware operator matches analytic ~1e-4.
- `demo_ff_streamfunction_design_matrix.py` — **the design = invert M**: with the material-aware M the
  target field is hit (2e-16); the free-space-designed coil **misses by 77%** in the iron system.
- (context) `demo_cc` (it is condensed FEM, not BEM — the line is the Green's function), `demo_dd` (WHEN
  to form M: only because the inverse design consumes it).
Run: `pip install -e packages/radia-mcp` then `python demo_ff_streamfunction_design_matrix.py`
(needs numpy, scipy, ngsolve 6.2.2604, netgen.occ; uses `radia_mcp.radia_ngsolve.fem_bem_coupling`).

## The concrete NEXT task — steps 1-4 DONE (`demo_hh_general_iron_design.py`, streamfunction branch)
Promote `demo_ff` from the concentric/modal toy to a **general coil inverse design with arbitrary
(non-concentric) iron**, wired to the existing stream-function module:
1. **[DONE]** Real winding-surface stream function psi — the order-p H1 nodal trace on the coil surface
   (378 DoFs), not spherical-harmonic modal amplitudes.
2. **[DONE]** Material-aware transfer matrix `M[target, psi-dof]` built directly from the Kelvin-FEM with
   the iron meshed as an arbitrary (non-concentric) blob: ONE sparse factorisation of the Kelvin-FEM +
   one back-substitution per coil DoF (the Dirichlet(coil)->field(target) specialisation of demo_v's Schur
   condensation). Targets sit in the PHYSICAL vacuum region (read directly -> no dependence on the
   inverse-Kelvin map convention).
3. **[DONE]** Inverse design psi = M^+ B_target (folded TSVD; a dense numpy TSVD here — wiring through
   `radia.streamfunction`'s ACA-TSVD is the remaining polish).
4. **[DONE]** Forward check: a fresh full Kelvin-FEM solve of the designed psi (`gf.Set(psi)+solve`, does
   NOT touch M) HITS the target (1e-14) while the free-space-designed psi MISSES by ~43% (stable 31-39%
   across mesh refinement). Physics anchored on the concentric sub-case vs the analytic layered transfer
   (rel 3e-3..2e-2), and M@psi == a fresh solve to 1e-15 (assembly).
5. **[OPEN — the immediate next task]** Benchmark M-build (sparse Kelvin-FEM Schur) vs the dense
   layered-Green / FE-BEM baseline (sparsity, conditioning, FE-coupling) — the selling point is "sparse,
   material-aware, no Green function". Also: a non-spherical coil former (cylinder) and an m≠0 tesseral
   target to drop the residual spherical symmetry.

## CRITICAL Kelvin-FEM gotchas (each is a ~1e7x blow-up or silent error if missed)
1. **Gauge:** the open Kelvin compactification has a constant near-null mode; a single ground POINT has
   ZERO capacity in 3D H1 -> use a **mean-zero NumberSpace constraint** (`int u dx = 0`, solver
   `inverse="umfpack"`) OR a finite ground ball. A bare GND vertex lets the constant blow up.
2. **Source must be EXACTLY discretely neutral** (normalise by the MESHED region volume, net ~1e-15); any
   residual net "charge" excites the near-null mode ~1e7x.
3. **Single-face periodic glue:** keep the truncation sphere ONE face (one `Identify`); carry interfaces
   by `IfPos(z, c1, c2)` coefficients, NOT a geometric hemisphere split (a split = imperfect glue =
   another near-null mode).
4. **Inverse map (periodic-glue convention):** the physical exterior value = the ball value at the
   inverted point `x' = (offset, 0, R^2/r)` with **NO (R/rho) weight** (that weight is the standalone-
   ball convention of `demo_p`).
5. **n=0 monopole** is spurious in the point-grounded ball (scope DtN spectra to n>=1; magnetostatics has
   no monopole anyway).

## Novelty status (VERDICT) — phrase claims carefully
**NOVEL, confidence 0.83** (targeted 9-agent check: 5 search + 4 adversarial refuters, all "not
preempted"). The specific combination — Kelvin/shell-transformed FE open boundary, Schur-condensed to a
DtN/transfer matrix carrying finite-permeability iron in the transformed exterior, used AS the
stream-function coil-design sensitivity kernel — was not found. **It FUSES Sugahara's OWN two threads**
(Kelvin open-boundary FEM: Extended Kelvin IEICE E108-C 2024/25, ECT-with-Kelvin IEEE TMAG 58(9) 2022;
free-space stream-function coil design: Koiso/Sugahara/Ida TSVD+ACA CEFC 2024, ACA+CMA-ES 3D coil IEEJ
2025) — frame the contribution as **the fusion** (defuses self-citation). Elsewhere coil design carries
iron via categorically different kernels (BEM/μ→∞ equipotential = bfieldtools Mäkinen-Zetter 2020;
image/modified-Green = Wang et al. Measurement 2024; saturated dipole = Landreman 2025; magnetization-
response = passive shimming) — never a transformed-FE DtN material operator.
**Recommended wording (NOT a bare "first"):** *"To the best of our knowledge, this is the first method to
use a Kelvin/shell-transformed FE open boundary, condensed by a Schur complement into a DtN (transfer)
matrix that carries the finite-permeability iron in the transformed exterior, directly as the sensitivity
kernel of a stream-function coil-design inverse problem. We are not aware of any prior work fusing these
two ingredients ..."*
**Residual checks before a 'first' claim:** (1) full texts of paywalled shield/coil papers (esp. Wang
2024 Measurement S0263224124008339); (2) Japanese grey-lit (IEEJ 静止器/マグネティックス研究会, J-STAGE,
CEFC/COMPUMAG 2022-26, in Japanese: Sugahara/Koiso/Sato/Ida + ケルビン変換 + 電流ポテンシャル); (3)
2025-26 preprints + the stellarator REGCOIL / current-potential line; (4) confirm no in-press Sugahara
paper already fuses them; (5) match wording to the exact construction; (6) accelerator-magnet (ROXIE/CERN)
field-quality design re-scan.

## Conventions (this repo)
- **Public boundary:** radia-ngsolve / Kelvin / stream-function are OPEN and public. Keep it open math.
  Never bring COMSOL/FEMM/JMAG content into public artifacts.
- **Figures:** all graphs via `radia_mcp.figure` (`paper_figure` -> `emit_paper_figure`), no titles,
  Times, vector PDF — NOT ad-hoc matplotlib.
- **Record insights** into the `radia_mcp.radia_ngsolve.dtn_coarse_mesh` knowledge module (queryable),
  and add an example + commit per discussion.
- **DtN/FEM-Kelvin is a CORE method** (CLAUDE.md core table), promotes into `src/radia/` when stable;
  the stream-function coil design is an APPLICATION DOMAIN (`radia.streamfunction`). This handoff bridges
  the two.
