# (ACA+)+TSVD least-norm solver -- examples

Accelerated least-norm solver for the **stream function method** of coil
design, generalised to **any Radia source family**.

Given M field (observation) points and N basis sources, the field-coupling
matrix is

```
A phi = B            A in R^{M x N},  M < N  (underdetermined)
A(i,j) = (a field component) at observation i produced by basis source j
```

The design problem "find source strengths `phi` that produce a desired field
`B`" is solved by the **TSVD-regularised pseudo-inverse**
`phi = V diag(1/S) U^T B`, truncated to `k` modes.  Building the dense `A` and
doing a full SVD is `O(N M^2)`.  Instead we factor `A ~= C D^T` with **ACA+**
(rank `k_aca << min(M,N)`) and TSVD only the small factors -- about `(M/k)^2`
faster.

Production module: [`src/radia/stream_function.py`](../../src/radia/stream_function.py)
(`aca_tsvd`, `pseudo_inverse_solve`, `solve`, `radia_field_kernel`).
C++ core: `src/core/rad_stream_function.cpp`.  See
[`../../docs/stream_function.md`](../../docs/stream_function.md) for the method,
API, and design notes.

## Kernel-agnostic by design

The solver embeds **no field kernel**.  The matrix entry `A(i,j)` is supplied by
a callback, so the same machinery serves any Radia source using Radia's
*already-implemented* field computation:

| Source family | Radia kernel | Example |
|---------------|--------------|---------|
| coils (thin wires) | Biot-Savart (`ObjFlmCur`, `ObjArcCur`) | `demo_coil_field_synthesis.py` |
| permanent magnets / soft iron | MMM / MSC (`ObjRecMag`, `ObjHexahedron`) | `demo_magnet_array.py` |

`radia_field_kernel(obs, sources, component, field)` builds the callback from a
list of Radia object handles via `radia.Fld` -- no per-application Biot-Savart
code.

ACA+ itself is delegated to the in-repo **HACApK** C library
(`cHACApK_acaplus`), the single source of truth for ACA+ in Radia.

## Files

| File | What it shows |
|------|---------------|
| `demo_coil_field_synthesis.py` | Coil design: N filament loops, solve loop currents for a target axial-region field; ACA+ compression + TSVD L-curve. |
| `demo_magnet_array.py` | Same solver on a permanent-magnet array (MMM/MSC field) -- proves kernel-agnosticism. |
| `bench_aca_vs_dense.py` | `(ACA+)+TSVD` vs naive dense `numpy.linalg.svd`: time / memory / rank, written to `results_aca_vs_dense.json`. |
| `demo_cmaes_magnet_design.py` | The **nonlinear counterpart**: CMA-ES (Optuna `CmaEsSampler`) optimises the magnetization *directions* (angles) of a magnet array for a uniform transverse field. Linear amplitude design -> (ACA+)+TSVD; nonlinear direction design -> CMA-ES (the "+ CMA-ES" half of SA-25-020). Needs `optuna` (optional). |
| `demo_coil_design_gz.py` | **End-to-end coil design**: cylindrical z-gradient (Gz) coil via the stream function method. Target `Bz=Gz*z` -> azimuthal ring currents (ACA+TSVD) -> stream function `psi(z)` -> equal-current wire rings -> verified on-axis gradient linearity. The axisymmetric Gz problem reduces to a full-ring (1D `psi(z)`) basis. |
| `demo_sf_to_peec_gz.py` | **Full workflow, loop closed**: SF design -> **single-stroke** (one continuous wire) smooth helix with blended crossovers -> CAD STEP (build123d Spline + Frenet swept solid) -> PEEC (`L`, `R`) -> exact Biot-Savart field -> verify `Bz` vs the design `Gz*z`. `--with-peec` adds the STEP + PEEC stages (needs build123d, in `radia`). |
| `demo_coil_design_gx.py` | **Transverse gradient (Gx), the 2D case**: a non-axisymmetric target `Bz=Gx*x` gives a genuine 2D surface stream function `psi(phi,z)` (a "fingerprint" pattern) -> marching-squares contour -> wire loops; verified `Bz` matches `Gx*x` to ~0.8% over the DSV. Numpy Biot-Savart kernel (avoids the ObjFlmCur bug); each loop is a SEPARATE closed conductor. |
| `demo_sf_to_peec_gx.py` | **Full workflow for Gx (transverse), loop closed**: SF design -> equal-current contours -> **single-stroke chain** that opens each closed fingerprint loop at its anchor and connects to the next via a cylinder-surface geodesic (= helical arc) -> CAD STEP -> PEEC (`L`, `R`) -> exact Biot-Savart field over the DSV. The Gx companion to `demo_sf_to_peec_gz.py`; the single-stroke chain trades some field accuracy (the connection arcs add stray Bz) for one-piece manufacturability. `--with-peec` adds the STEP + PEEC stages (needs build123d). |
| `demo_planar_uniform_fem_psi_advanced.py` | **FE-direct psi + advanced regularisation + surface deformation**.  Four research knobs in one demo: (a) `--regularize h1_sigma --sigma-cf EXPR` for 1/σ-weighted H1 (true ohmic dissipation with non-uniform conductivity / forbidden regions); (b) `--regularize inductance_diag` for a lumped self-inductance proxy (full inductance form needs ngsolve.bem MaxwellSL, deferred); (c) `--regularize linf --jmax VAL` for L∞ peak-current capping via scipy SLSQP (experimental, slow); (d) `--deform --deform-params {zoff,bump,zoff+bump} --deform-trials N` for an Optuna CMA-ES outer loop over surface deformation parameters.  Measured improvements vs the flat H1 baseline (RMS 2.09 %, p2p 6.81 %): Gaussian 1/σ -> 1.17 %; 1-param zoff deform (10 trials) -> 1.50 %; 3-param bump deform (20 trials, 22 s) -> **0.77 % (-63 %)**, 4-param zoff+bump (50 trials, 33 s) -> 0.85 %. All combinable with Path-A iteration.  Also (e) `--regularize l2_aca`: routes L2 min-norm through `radia.stream_function.aca_tsvd` (validated identical result to numpy lstsq) and (f) `--order N`: arbitrary H1 polynomial order, with **order=3 found to be the sweet spot** (RMS 0.51 %, p2p 1.83 %) at maxh=0.025.  Non-monotone in p due to the FE / contour / single-stroke discretisation interaction.  Deformation freedom + high order can BACKFIRE when baseline is near-optimal (order=3 + bump -> 0.82 % vs order=3 alone 0.51 %); design rule = turn deformation off when single-shot accuracy is already sub-0.5 %. |
| `demo_planar_uniform_fem_psi_aca.py` | **FE-direct psi + Radia HACApK (ACA+)+TSVD** -- same H1 matrix as the `_fem_psi.py` demo, but solved via `radia.stream_function.aca_tsvd` (kernel-agnostic ACA+ + Method-3 TSVD).  Wraps the M×N_free matrix as an `entry(i, j) -> float` callback so the same machinery the basis-loop demo uses carries over to the FE-direct setup.  Best RMS 0.67 % at 100 Path-A iters (slightly worse than direct lstsq because of `aca_eps=1e-10` truncation, but the **factorisation is cached and re-used across iterations** at ~1 ms per back-sub).  Validates the `(A)` path: when material kernels (Radia MMM iron yoke, shielded coil) or large M make full-matrix assembly impractical, the same callback contract works -- just point it at an on-demand integrator.  Ready to swap the matrix assembly for the `ngsolve.bem` H-matrix operator (6.2.2604+) when Joachim ships the ACA-based hierarchical compression (currently 2604 exposes FMM only). |
| `demo_planar_uniform_fem_psi.py` | **psi as DIRECT H1 FE unknown on an NGSolve 2D mesh**, same target as the basis-loop demo below.  ``--regularize h1`` (default) solves ``min psi^T S psi`` s.t. ``A psi = B_target`` for the smoothest psi that hits the target exactly.  Single-shot RMS = 2.09 %.  With ``--compensated-iter 100 --compensated-step 0.05`` the Path-A iteration CONVERGES MONOTONICALLY (iters 40-47 drop 0.62 % -> 0.49 % without backtracking) -> final RMS **0.47 %** (-84 % vs basis-loop), p2p/mean = 1.64 %.  This is the empirical proof that the "naive Picard doesn't converge" finding for the basis-loop case is SPECIFIC to grid-sampled psi -- a continuous FE psi gives a smooth chain-field response and Path A genuinely contracts.  Extends naturally to non-planar OCC surfaces (cylinder, sphere, conformal). |
| `demo_planar_uniform_coil.py` | **Planar uniform-Bz coil (the easy end of single-stroke complexity)**.  Source = square plane at z=0; target = uniform Bz=B0 over a square region at z=h.  SF produces concentric closed contours; single-stroke = spiral (Kuijpers Method-1 with cut line at +x axis).  Baseline RMS = 2.99 %, peak-to-peak / mean = 9.59 %; with `--compensated-iter 100 --compensated-step 0.05` (Path A best-effort) RMS drops to **0.58 %** (-80 %), p2p/mean = 2.17 % (-77 %).  Validates that Path A IS effective on simpler topologies even though it failed for the Gx fingerprint -- a clean demonstration of the "complexity tier" framing in `radia-mcp aca_tsvd(single_stroke)`. |
| `view_sf_coil_gx_gmsh.py` | **GMSH visualiser** for the Gx coil.  Three modes: `--mode contours` (default, **recommended**) writes the SF design's CLOSED CONTOUR FILAMENTS only (each fingerprint loop as its own 1D Physical Group, NO connection arcs) -- the true SF output, with the host cylinder as a 2D reference; `--mode chain` writes the lobe-aware single-stroke chain (= what PEEC sees, with 4-quadrant traversal + 3 inter-quadrant geodesic arcs; the same chain `demo_sf_to_peec_gx.py --chain-method lobe` builds); `--mode step` merges the loft-chain STEP from `--with-peec`.  Includes preventive code (`gmsh.initialize(["-noconfig"])` + explicit `General.GraphicsPositionX/Y` + `Width`/`Height`) so the GMSH window can't restore to an off-screen second-monitor coordinate.  Uses pip-gmsh blocking `fltk.run()` (CLAUDE.md GMSH policy). |

## Running

```bash
python demo_coil_field_synthesis.py
python demo_magnet_array.py
python bench_aca_vs_dense.py
python demo_cmaes_magnet_design.py        # needs optuna (pip install optuna)
python demo_coil_design_gz.py             # end-to-end Gz gradient coil design
python demo_sf_to_peec_gz.py --with-peec  # full SF -> CAD(STEP) -> PEEC -> field
python demo_coil_design_gx.py             # transverse Gx gradient (2D surface SF)
python demo_sf_to_peec_gx.py --with-peec  # full Gx SF -> single-stroke -> CAD -> PEEC
python view_sf_coil_gx_gmsh.py             # SF contour filaments (no connection arcs) -- the real design
python view_sf_coil_gx_gmsh.py --mode chain # single-stroke chain (with connection arcs) -- what PEEC sees
python view_sf_coil_gx_gmsh.py --mode step  # the loft-chain STEP file from --with-peec
```

Each script is standalone (no Cubit, no panel UI).  `matplotlib` is optional:
if installed, the demos save a PNG next to the script; otherwise they print an
ASCII summary only.  `demo_cmaes_magnet_design.py` additionally needs `optuna`
(it prints a friendly message and exits cleanly if optuna is missing).

## Expected results

- **`demo_coil_field_synthesis.py`** (M=25 field points, N=64 loops): smooth
  off-plane field is low rank, so ACA+ stops at `k_aca` well below `min(M,N)=25`.
  The TSVD L-curve shows the residual `||A phi - B|| / ||B||` dropping as modes
  are added, with the solution norm `||phi||` rising -- the usual
  regularisation trade-off.  A few modes already reproduce the target to ~1%.
- **`demo_magnet_array.py`** (N permanent magnets): the `(ACA+)+TSVD`
  factorization reconstructs the Radia MMM/MSC coupling matrix to `< 1e-5`
  relative, identical machinery, zero coil-specific code.
- **`bench_aca_vs_dense.py`**: for a smooth (low-rank) kernel, `(ACA+)+TSVD`
  matches the dense singular values to ~1e-12 while running markedly faster as
  `N` grows and `k_aca` stays small.
- **`demo_cmaes_magnet_design.py`**: a 16-dimensional continuous optimisation
  (one magnetization angle per pixel).  CMA-ES drives the relative field-match
  objective down by ~3x from its first trial, producing an approximately uniform
  transverse field with small cross-components.  The residual is set by the
  finite array (physical limit), not the optimiser.
- **`demo_coil_design_gz.py`**: ACA+ compresses the on-axis ring operator to
  `k_aca ~ 7`; the continuous ring-current solution reproduces `Bz=Gz*z` to
  `~4e-4`; the contoured ~32-wire coil achieves `dBz/dz` within ~0.5% of target
  with ~1.4% on-axis nonlinearity over the DSV -- a textbook generalised
  Maxwell-pair gradient coil.
- **`demo_sf_to_peec_gz.py`** (`--with-peec`, ~16 turns): the single-stroke
  conductor (~15 m, one continuous wire) reproduces the design gradient
  (`dBz/dz ~ 0.99`, ~2.6% nonlinearity); the helix sweeps to a clean STEP solid
  (Frenet frame + auto wire radius so turns don't self-intersect); PEEC returns
  `L ~ 38 uH`, `R ~ 16 mOhm` at 1 kHz.  Confirms the SF design survives the
  single-stroke manufacturing constraint.  (CoilBuilder is for planar
  racetrack/saddle coils; a solenoidal helix uses the smooth-helix + Spline path.)
- **`demo_coil_design_gx.py`** (transverse Gx, the 2D case): the
  non-axisymmetric target `Bz=Gx*x` produces a genuine 2D surface stream
  function `psi(phi,z)` (the classic "fingerprint" pattern).  Unlike the
  low-rank axisymmetric Gz problem, the transverse target fills the operator's
  rank -- ACA+ reaches `k_aca = min(M,N) = 123` (no compression here; the 2D
  problem is intrinsically richer).  Marching-squares contours `psi` into ~68
  saddle-shaped wire loops driven INDEPENDENTLY (each closed loop is a
  separate conductor), and the reconstructed `Bz` matches the design `Gx*x`
  to ~0.8% RMS over the DSV.  Uses a numpy Biot-Savart kernel (the ObjFlmCur
  tilted-loop path is unreliable for these non-planar loops).
- **`demo_sf_to_peec_gx.py`** (`--with-peec`): same Gx SF design, then
  threaded into ONE continuous conductor.  Five single-stroke methods are
  available via `--chain-method`.  **When joining a new coil's contours,
  use the [`single-stroke-chain`](../../.claude/skills/single-stroke-chain/SKILL.md)
  skill** — the connection has no clean closed-form optimum and is a
  reason-and-verify task (build → measure DSV RMS → keep only if it beats
  the baseline → escalate to Path-A).
    - `field_aware` (default, **recommended**, 2026-05-31): the `kuijpers`
      lobe/current-sign visiting ORDER (the dominant factor) + each
      contour's cut chosen by coordinate descent to minimise the AZIMUTHAL
      arc to its chain neighbours (axial `dz` is free).  This lets the rung
      stray fields cancel more symmetrically over the DSV.
    - `kuijpers` (prior best): Kuijpers, Jansen, Lomonova, Compumag 2023
      paper [525] Method-1.  Per-lobe fixed cut phi (`+x` at 0, `-x` at π),
      one straight axial (phi, z) "rung" per contour pair (Fig.4 pattern).
    - `lobe`: 4-quadrant classification + within-quadrant spiral.
    - `greedy`: legacy global nearest-neighbour (visually "wasted" arcs).
    - `nn_blend`: CAUTIONARY negative result — geometric-shortest balanced
      cut.  Shortest rungs but WORST field (geometry-only NN interleaves
      the lobe current signs).  Do not use for field accuracy.
  Field comparison at the same SF design (24 phi x 40 z surface, 12 levels):

  | method       | DSV RMS  | x-axis nonlin |
  |--------------|----------|---------------|
  | **field_aware** | **9.29 %** | **7.20 %**  |
  | kuijpers     | 16.24 %  |  9.73 %       |
  | greedy       | 21.78 %  | 11.40 %       |
  | lobe         | 23.98 %  |  9.67 %       |
  | nn_blend     | 65.08 %  | 38.90 %       |

  `field_aware` gives the lowest DSV RMS — 43 % below `kuijpers` at this
  config, and 30–54 % below across an nlevels 10/12/16 + mesh sweep.  The
  key finding: field impact is NOT predicted by rung length (geometric or
  azimuthal) — `nn_blend` and `field_aware` have near-equal azimuthal rung
  totals but 8× different RMS.  The lever is the current-sign-respecting
  ORDER, then symmetric cut placement so the rung stray fields cancel; no
  single scalar metric captures it (hence the skill).  This SOFTENS the
  HARD-tier "16 % ceiling" — it was a `kuijpers`-method artifact, not a
  fundamental bound.

  **Path A compensated iteration** (research): ``--compensated-iter N
  --compensated-step alpha`` runs N iterations of the fixed-point
  ``phi += alpha * pseudo_inverse(B_target - I_w * Bz_chain_unit)`` to
  bake the chain crossover field back into the SF solve.  Each iteration
  re-uses the cached ACA+TSVD factorisation, so the cost is one
  back-substitution + one chain rebuild per iter.

  Status: the naive Picard form does NOT converge (chain construction
  is too nonlinear in phi -- level-set topology jumps each iter).  The
  iteration OSCILLATES around the baseline but occasionally finds slightly
  better neighbourhoods; the demo tracks the best psi seen and uses it for
  downstream evaluation.  It is step-sensitive (a good ``--compensated-step``
  helps; a bad one finds no gain).  Best observed, stacked on the
  ``field_aware`` chain:

  | chain         | + Path-A (step 0.3, 40 iter) | x-nonlin       |
  |---------------|------------------------------|----------------|
  | field_aware 9.29 % | **8.11 %** (-13 %)      | 7.20 % -> 6.27 % |
  | kuijpers 16.24 %   | 15.28 % (-6 %)          | 9.73 % -> 5.99 % |

  So the two best ideas COMPOSE: ``field_aware`` (better connection) plus
  Path-A (SF design absorbs the residual rung field) reaches 8.11 % DSV RMS,
  roughly half the old ``kuijpers`` baseline.  See
  ``radia-mcp aca_tsvd(single_stroke)`` for the open extensions (Anderson
  acceleration, frozen-topology, B-spline continuous SFD, multivalued-
  potential formulation) that may give true convergence.
  The chain is one physical wire (~22 m); PEEC returns `L ~ 17 uH`,
  `R ~ 80 mOhm` at 1 kHz (port = chain start to chain end, auto-sized round
  filament cross-section so adjacent contour wires don't overlap).  The
  trade-off vs `demo_coil_design_gx.py`: the connection arcs carry the
  series current too and add stray Bz, so the fitted gradient on the x-axis
  is ~0.93 of target with ~11% nonlinearity and ~22% RMS over the full DSV
  -- the price of one-piece manufacturability without symmetric pairing.
  CAD step exports a **multi-piece loft chain**: circular cross-sections
  placed densely along the chain with a parallel-transported (twist-free)
  frame, lofted in short pieces of ~20 sections each (a single OCC loft of
  the full chain fails because the saddle pattern's cumulative twist
  defeats `BRepOffsetAPI_ThruSections`).  This is the build123d analogue
  of "Cubit loft along curve" with a path guide; the result is a Compound
  of solid spools that together cover the wire, ~590 KB STEP.  Further
  accuracy = symmetric pairing of +x/-x lobes + CMA-ES connection
  routing (the "+CMA-ES" half of SA-25-020).

## References

- Sugahara Lab, "ACA-accelerated stream function method + CMA-ES",
  IEEJ Joint Technical Meeting on Static Apparatus / Rotating Machinery,
  SA-25-020 (manuscript Method 2/3).
- HACApK (ppOpen-HPC, MIT): `src/ext/HACApK/`.
- Validation reference: the f2py-wrapped Fortran `coil_solver.f90`
  (`method_aca_tsvd_1/2`), matched bit-for-bit by
  `tests/test_stream_function.py::test_matches_f90_reference`.
