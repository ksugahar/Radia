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

Docs promotion: the maintained narrative lives under
[`docs/stream_function/`](../../docs/stream_function/) as theory,
regularization, deformation, and benchmark pages/notebooks. Historical
source-only archives were pruned from docs; use git history for deleted
development snapshots.

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
| `demo_shim_coil_purity.py` | **Spherical-harmonic shim/gradient coil purity**.  Designs SF coils on ONE cylindrical former for a sequence of NAMED harmonic targets (Gz, Gx, Z2, X2-Y2, ZX, **Z4**) via `calc_streamfunction.py --target-harmonic`, then decomposes the achieved DSV field back into solid harmonics (`--harmonic-lmax 4`) to report the target-harmonic PURITY + the named largest CONTAMINANT -- the canonical MRI gradient/shim quality metrics.  The SAME solid-harmonic table builds the target and analyses the result (exact round-trip).  All targets fit to near-unit purity; the difficulty grows with order (l<=3 residual <0.1 %, the 4th-order Z4 shim ~1.5 % with a named Z3 contaminant -- the field of an l-th harmonic scales as r^l, so high-l shims are intrinsically harder over a fixed DSV).  Mesh-gen + each design run in subprocesses (orchestrator stays NGSolve-free).  Writes `demo_shim_coil_purity.{json,png}`. |
| `demo_active_shield.py` | **Active shielding (primary + shield), Turner 1986**.  Designs a Gz coil two ways -- UNSHIELDED (inner primary only) and SHIELDED (primary + an outer shield cylinder designed JOINTLY via `calc_streamfunction.py --shield-vol`, so primary+shield together hit the DSV target while nulling the field over a LARGE external region) -- and compares the stray field at INDEPENDENT external points.  Genuine **~86x (39 dB)** stray reduction (100-1700x in the well-covered far field), DSV homogeneity preserved.  Documents the CRITICAL lesson (measured): the external null region must COVER the exterior -- a thin mid-plane slice overfits locally (~4x) and worsens the field at larger z; the un-sampled gap between shield and null region stays unshielded.  Writes `demo_active_shield.{json,png}`. |
| `demo_fmm_biot_savart.py` | **FMM Biot-Savart benchmark** -- `ngsolve.bem.BiotSavartCF` (spherical-harmonic multipole expansion, NGSolve 6.2.2604+) vs the LinearForm design-matrix path, evaluating a designed coil's field at EXTERIOR points.  The surface current K = n x grad(psi) is discretised into per-triangle K-dipole filaments (`AddCurrent`); both kernels evaluate the same psi at the same exterior shell.  Honest finding: they agree to ~5 % (the dipole-discretisation error), but BiotSavartCF is **slower at small N** (large multipole setup overhead) and is **EXTERIOR-ONLY** (returns nan inside the source sphere -- so it CANNOT replace the LinearForm path for the DSV-interior design matrix).  Useful only for far-field post-evaluation at very large N.  Writes `demo_fmm_biot_savart.{json,png}`. |
| `view_sf_coil_gx_gmsh.py` | **GMSH visualiser** for the Gx coil.  Three modes: `--mode contours` (default, **recommended**) writes the SF design's CLOSED CONTOUR FILAMENTS only (each fingerprint loop as its own 1D Physical Group, NO connection arcs) -- the true SF output, with the host cylinder as a 2D reference; `--mode chain` writes the lobe-aware single-stroke chain (= what PEEC sees, with 4-quadrant traversal + 3 inter-quadrant geodesic arcs; the same chain `demo_sf_to_peec_gx.py --chain-method lobe` builds); `--mode step` merges the loft-chain STEP from `--with-peec`.  Includes preventive code (`gmsh.initialize(["-noconfig"])` + explicit `General.GraphicsPositionX/Y` + `Width`/`Height`) so the GMSH window can't restore to an off-screen second-monitor coordinate.  Uses pip-gmsh blocking `fltk.run()` (CLAUDE.md GMSH policy). |
| `demo_regcoil_fusion.py` | **Fusion (stellarator Stage-2) coil design -- mini-REGCOIL / NESCOIL**.  The SAME surface stream function the gradient/IH demos use, applied to the fusion winding-surface current-potential problem: given a target normal field `B.n` on a PLASMA boundary, solve for `psi` on a WINDING SURFACE around it whose Biot-Savart `n.B` reproduces it (iso-contours = the coils).  **Four parts**: (1) two PRODUCIBLE targets -- a uniform vertical field (PF / equilibrium coil) and a non-axisymmetric `sin(theta)cos(2 phi)` (stellarator-like) -- reproduced to MACHINE PRECISION (`B.n` residual ~2e-9), proving the SF designer does the Stage-2 forward map exactly; (2) the REGCOIL **L-curve** on a genuinely-hard high-mode target `sin(3 theta)cos(5 phi)` (decays across the plasma-coil gap), sweeping `alpha` to trace the classic `(B.n residual, peak |grad psi|)` trade-off (knee at `alpha_rel ~ 2e-2`); (3) **NET CURRENT (the multivalued / secular term)** -- a single-valued `psi` carries zero net current through each torus hole; the full current potential gains `Psi = psi + (G/2pi)zeta + (I/2pi)theta`, whose TWO extra DOFs are the winding surface's first cohomology generators.  Their COUNT (`b1 = 2`, genus 1) is CONFIRMED gmsh-free via the Euler characteristic of the triangulated torus (`b1 = 2 - chi`); the volume T-Omega cohomology cut engine `src/radia/cohomology_cut.py` is itself gmsh-free (`radia.cohomology`); the net-poloidal-current (TF) secular term is verified to give the textbook `B_tor ~ 1/R` toroidal field (Ampere, `B_tor*R` const to 0.2% inside the tube, ~0 outside).  Key physics: the TF field is TANGENT to the plasma (`B.n` footprint >1000x smaller than the net-toroidal one), so the net poloidal current is a PRESCRIBED engineering parameter (set it for the on-axis `B_tor`: 1 T -> 1.5 MA) -- exactly why REGCOIL takes `net_poloidal_current` as an INPUT; (4) a **VMEC-shaped plasma boundary** -- the circular torus is replaced by a non-axisymmetric **rotating-ellipse** boundary in the VMEC Fourier form (`R = sum RBC cos(m th - n NFP ph)`, `Z = sum ZBS sin(...)`, NFP=3) with analytic surface normals; a real machine's free-boundary equilibrium drops in with `--wout wout_*.nc` (netCDF4 reader, round-trip-verified against the VMEC schema in the golden; stellarator-symmetric wout only -- a non-symmetric `lasym=T` file is rejected, not silently truncated).  **Honest scope**: parts 1-2 single-valued `psi`; part 3 adds the multivalued term (generator count computed, generators analytic on the torus -- general surfaces use `cohomology_cut.py`); part 4's default boundary is an analytic rotating-ellipse MODEL (not a converged equilibrium); coil force / stress and winding-surface optimisation (FOCUS) are the named next steps.  Writes `demo_regcoil_fusion.json` + a 2x2 lab figure (3D coil, L-curve, TF 1/R, VMEC boundary). |
| `demo_regcoil_fusion_advanced.py` | **Advanced fusion: coil force/stress, a REAL VMEC equilibrium, FOCUS standoff + winding-SHAPE** -- the four "named next steps" of `demo_regcoil_fusion.py`, reusing its helpers.  (A) **Coil force/stress**: the Lorentz force per area `f = K x B_avg` (B_avg = the +/-eps average of the coil's own field across the current sheet); for the net-poloidal-current (TF) coil this is verified to equal the magnetic pressure `B_tor^2/(2 mu0)` (ratio ~0.99) and to CONCENTRATE on the INBOARD leg (~5x the outboard) -- why a TF coil is inboard-stress-limited.  (B) **A real equilibrium**: designs against the boundary of the li383 (NCSX-like, NFP=3, quasi-axisymmetric) reference `wout` from simsopt -- `--wout PATH` for any VMEC output, else it fetches li383 (121 kB) to a cache; the Stage-2 forward map reproduces the vertical-field `B.n` on the genuine stellarator boundary to ~4e-8.  (C) **FOCUS standoff study**: sweeping the winding-surface gap shows coil complexity `peak |grad psi|` is MONOTONIC in the gap (closer = simpler, ~50x over the sweep), so the distance optimum is CONSTRAINT-BOUND (push to the minimum engineering standoff `d_min`).  (D) **FOCUS winding-SHAPE** (the core FOCUS contribution): a `_surface_mesh_from_grid` builds an ARBITRARY winding surface (any (theta,phi) point grid) so the winding can be CONFORMAL to the plasma, not just a circular torus.  For an elongated (kappa=2) plasma, blending the winding from CIRCULAR (varying gap) to CONFORMAL (uniform gap, plasma offset along its normals) at the SAME min standoff cuts coil complexity `peak |grad psi|` by **~34%** -- the winding shape, not just distance, is a real design lever.  **Honest scope**: force is the magnetic force per area (the stress DRIVER, N/m^2), not a structural hoop-stress model; the shape study sweeps a 1-parameter circular->conformal blend (a full Fourier-mode winding-surface optimiser is the next step).  Writes `demo_regcoil_fusion_advanced.json` + a 2x2 lab figure (TF stress, li383 boundary, FOCUS standoff, conformal-vs-circular cross-section). |

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
python demo_shim_coil_purity.py           # spherical-harmonic shim/gradient coil purity (l<=4)
python demo_active_shield.py              # active shielding (primary+shield Gz coil, ~86x stray reduction)
python demo_fmm_biot_savart.py            # FMM Biot-Savart (ngsolve.bem) vs LinearForm, exterior eval
python view_sf_coil_gx_gmsh.py             # SF contour filaments (no connection arcs) -- the real design
python view_sf_coil_gx_gmsh.py --mode chain # single-stroke chain (with connection arcs) -- what PEEC sees
python view_sf_coil_gx_gmsh.py --mode step  # the loft-chain STEP file from --with-peec
python demo_regcoil_fusion.py              # fusion Stage-2: winding-surface current potential -> plasma B.n
python demo_regcoil_fusion.py --no-plot    # JSON only (no matplotlib / radia-mcp)
python demo_regcoil_fusion_advanced.py     # coil force/stress + real li383 wout + FOCUS standoff
python demo_regcoil_fusion_advanced.py --no-fetch --no-plot   # offline (skip the li383 download)
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

## End-to-end validation vs an independent codebase

`verify_coil_field_independent.py` closes the design loop on a **real
engineering geometry** (MRI-gradient-coil scale: cylinder r = 0.15 m,
L = 0.5 m, DSV sphere r = 0.05 m) and checks the result against an
**independent field codebase**, not against itself:

1. design a coil with the general FE-direct stream-function stack
   (`calc_streamfunction.py`),
2. discretise it into orientation-consistent equal-current contour turns,
3. compute the designed coil's actual field over the DSV **two ways** --
   the numpy straight-segment Biot-Savart used inside the designer
   (`_segment_field_B`, itself checked vs the circular-loop analytic to
   6 digits) **and** Radia's C++ `rad.ObjFlmCur` + `rad.Fld` (a separate
   codebase).

```
python verify_coil_field_independent.py --order 2 --nlevels 10
```

VERIFIED (order 2, results in `verify_coil_field_independent.json`):

| case | target | turns | design vs **Radia C++** | field over DSV |
|------|--------|-------|-------------------------|----------------|
| uniform  | `Bz = 1` | 10 | **3.5e-11** (vector-diff rel) | uniformity 0.46 % |
| gradient | `Bz = x` | 35 | **1.1e-8** (vector-diff rel)  | G = 0.75 T/m, nonlinearity 4.8 % (35 contours open) |
| gradient_confined | `Bz = x`, `--confine abe` | 21 | **1.2e-8** | G = 0.99 T/m, nonlinearity **1.0 %** (0 open) |

The engines agree to **8-11 digits** -- the designer's field is confirmed by
an independent C++ codebase, so it is not merely self-consistent.  The
*gradient* (unconfined) row shows the honest pre-fix state: on a *finite*
cylinder the equal-current contours run off the ends (35/35 open), closing
them with a rim chord degrades the delivered linearity to 4.8 %.  The
*gradient_confined* row is the **fix**: confining the current to the patch
(`--confine on`, psi=0 on the former edge) closes every contour (0 open) and
the SAME cross-validated coil reaches **1.3 %** nonlinearity on the same short
former.  Locked by `tests/panels/test_streamfunction_golden.py::test_streamfunction_field_cross_codebase`
and `::test_streamfunction_confine_closes_contours`.

### Field-aware single-stroke (no sheet-metal crutch)

The gradient caveat above has a clean fix.  An error-budget sweep
(`calc_streamfunction.py --method manufacture`) shows the ~0.5 single-current
error is **entirely the open-contour rim-chord artifact**, not a fundamental
single-current limit and not the connectors:

| former L | n_open | single-current rms |
|----------|--------|--------------------|
| 0.5 m    | 42/42  | 0.54   |
| 1.0 m    | 10/28  | 0.075  |
| **1.6 m**| **0**  | **0.0073** |

With the former long enough that the contours **close** (`n_open_contours = 0`),
the equal-current discretisation reproduces the target to 0.7 % as SF theory
predicts.  The connectors are then the next error, and `--chain field_aware`
(default) chooses each loop's entry/exit cut to minimise the *full one-current
wire error* `min_I ||I*(loops+connectors) - B||` -- reaching the separate-turns
floor with **no `--distort`**:

| chain | L=1.6 (closed) | L=0.5 (open) |
|-------|----------------|--------------|
| `nn`          | 0.213 | 0.246 |
| `field_aware` | **0.031** | **0.180** |

So the manufacture rule is: **close the contours (`n_open_contours = 0`), then
field-aware chaining gives a single-stroke, single-current wire without the
sheet-metal `--distort` crutch** (which stays as an optional extra).  Two ways
to close them:

- `--confine abe` (**recommended**) -- the canonical Abe edge-equipotential
  current-potential BC (M. Abe, IEEE Trans. Magn., DUCAS; Appendix eq.6
  `T = R T_IN`): each PHYSICAL boundary edge is tied to ONE free constant
  (A-1: no current crosses the edge) with one ground (A-3).  Closes the
  contours on **any** former AND works for gradient *and* solenoid (the two
  cylinder ends take different free constants).  Implemented as a DOF-reduction
  matrix `R`; seam edges are excluded by element adjacency.  The
  gradient_confined row above (short L = 0.5 m, 21 turns, 0 open, **1.0 %**
  separate-turn nonlinearity, cross-validated) uses it.
- `--confine on` -- the simpler special case psi = 0 on every edge.  Fine for
  a single-feed gradient coil; **breaks** solenoid-type targets (the two ends
  are forced to the same psi).
- or simply use a longer former (`L = 1.6 m` above) so the contours close on
  their own.

Caveat: `abe`/`on` make the DESIGN and the SEPARATE-TURN coil accurate and the
contours close; the single-stroke WIRE error additionally depends on the
chaining of the (now boundary-confined) contours -- `abe`'s boundary-hugging
contours can chain into a slightly worse one-wire path than `on`, so the best
*single-stroke* choice is target-specific.  With enough turns the short-former
Gx single-stroke reaches ~1.5-2 % (`--nlevels 30`), no distort.  Locked by
`test_streamfunction_field_aware_chain` and `test_streamfunction_confine_closes_contours`.

## Spherical-harmonic shim/gradient targets + field purity

`--target-harmonic` lets you ask for a field by its MRI gradient/shim NAME
(`Z`, `X`, `Z2`, `ZX`, ..., up to the 4th-order `Z4`, `C4`, ...; or `l=L,m=M` /
`(L,M)`, optionally weighted/summed) instead of a hand-typed polynomial.  In
`design` mode the ACHIEVED Bz over the DSV is then decomposed back into solid
spherical harmonics (`--harmonic-lmax L`), reporting the target's PURITY (its
field fraction) and the largest NAMED contaminant -- the canonical MRI
gradient/shim quality metrics.  The SAME solid-harmonic table builds the target
and analyses the result, so the round-trip is exact.

`demo_shim_coil_purity.py` designs six targets on ONE cylindrical former
(r = 0.15 m, L = 0.5 m; DSV r = 0.05 m; order 2; `--confine abe`):

| target | order l | dominant | purity | impurity (LSQ resid.) | largest contaminant |
|--------|---------|----------|--------|-----------------------|---------------------|
| Z  (Gz)     | 1 | Z  | 1.00000 | 3.1e-05 | Z4  2.1e-05 |
| X  (Gx)     | 1 | X  | 1.00000 | 1.5e-04 | Z2X 6.9e-05 |
| Z2          | 2 | Z2 | 1.00000 | 3.2e-04 | Z3  2.0e-04 |
| C2 (X2-Y2)  | 2 | C2 | 1.00000 | 3.3e-04 | ZC3 2.3e-04 |
| ZX          | 2 | ZX | 1.00000 | 2.7e-04 | Z2X 1.7e-04 |
| **Z4**      | 4 | Z4 | **0.99983** | **1.5e-02** | **Z3 9.5e-03** |

Every target fits to near-unit purity; the difficulty grows with the harmonic
order because over a fixed small DSV the field of an l-th harmonic scales as
`r^l`, so the 4th-order Z4 shim is the hardest -- yet still 0.9998 pure, with
the decomposition naming its ~1 % Z3 contaminant.  The solid-harmonic basis
spans `l <= 4` (1+3+5+7+9 = 25 harmonics); each entry is a Laplacian-zero
polynomial (the precondition for a valid current-free Bz component), golden-
locked in `tests/panels/test_streamfunction_golden.py`
(`test_harmonic_basis_is_harmonic`, `test_harmonic_l4_forms_and_decompose`).

## Active shielding -- primary + shield gradient coil (`--shield-vol`)

An actively-shielded gradient coil (Mansfield & Chapman 1986; Turner 1986)
adds a second, OUTER cylindrical surface (the **shield**) whose current
cancels the **stray field** outside the assembly -- essential in MRI so the
switching gradients do not induce eddy currents in the cryostat.
`calc_streamfunction.py --shield-vol SHIELD.vol --shield-eval-vol EXT.vol`
designs the primary and shield JOINTLY: one stacked least-squares system fits
the Gz target inside the DSV **and** nulls the field at external sample points
(block-diagonal seminorm, the `RegularizedTSVD` machinery unchanged).

`demo_active_shield.py` compares an unshielded (primary-only) and a shielded
(primary+shield) Gz coil at INDEPENDENT external points (so the stray is the
HONEST generalising field, not the circular constraint-point residual):

| design | DSV homogeneity | stray @ external | factor |
|--------|-----------------|------------------|--------|
| unshielded (primary only)  | 4.0e-05 | 1.65   | --   |
| shielded (primary + shield)| 9.5e-05 | 0.019  | **86x (39 dB)** |

The stray-vs-radius profile shows **100-1700x** suppression in the
well-covered far field (r > 0.25 m).  **Critical lesson (measured)**: the
external null region must COVER the exterior, not a thin mid-plane slice -- a
too-small region overfits locally (~4x) and makes the field WORSE at larger z;
the un-sampled gap between the shield and the null region stays unshielded
(~2-3x).  The honest reported metric is `stray_rms` (field at independent
external measure points); the `stray_fit_rms` at the constraint points is the
circular fit residual and is NOT the shielding quality.  An analytic external-
multipole-moment constraint (vs point sampling) is the documented next step.
Golden-locked: `test_streamfunction_active_shielding`.

## FMM Biot-Savart -- `ngsolve.bem.BiotSavartCF` vs LinearForm

NGSolve 6.2.2604+ ships an FMM-style hierarchical Biot-Savart in
`ngsolve.bem` (`BiotSavartCF`, a spherical-harmonic MULTIPOLE expansion).
`demo_fmm_biot_savart.py` benchmarks it against the LinearForm design-matrix
path for evaluating a designed coil's field, with the honest conclusion:

- **Accuracy**: discretising K = n x grad(psi) into per-triangle K-dipole
  filaments and expanding gives ~5 % vs the exact LinearForm surface integral
  (the dipole-discretisation error).
- **Speed**: BiotSavartCF is SLOWER at small N (~1 s multipole setup vs ~0.1 s
  for the LinearForm matvec at N~60); it would only amortise its setup for
  N >> 1000 far-field points.
- **Hard limitation**: BiotSavartCF is a multipole expansion for sources
  ENCLOSED in a sphere -- it returns nan for points INSIDE that sphere.  The
  SF **design matrix** lives at the DSV INTERIOR (inside the coil), where
  BiotSavartCF cannot evaluate at all, so it **cannot replace** the LinearForm
  path for design; it is only a far-field post-evaluation tool.

This is the scope-clarification in the "FMM Removed from Radia core" policy in
practice: FMM math fits free-space far-field Biot-Savart, not the near-field
interior design problem.

## Contour drawing = flux-line drawing (same principle)

The iso-contours of the stream function are drawn by the same rule as magnetic
flux lines: between two adjacent lines flows a FIXED amount (current for psi,
flux for the vector potential) -- Abe's "between nodes i, j flows T_i - T_j".
So equal-`psi`-interval contouring automatically gives wire density
proportional to `|grad psi| = |K|`, the same density rule the flux-line "bubble
system" (Hirahatake/Noguchi/Igarashi/Yamashita) enforces with bubbles of radius
`r ~ 1/sqrt(|B|)`.  Two refinements of `manufacture` follow from this:

- **`--contour-sub N` (order-p contour)** -- by default the wire contours march
  on the vertex (order-1) psi.  For an order-2/3 design that throws away the
  edge/face DOFs.  `--contour-sub 3` subdivides each surface triangle 3x3 and
  evaluates the FULL-order psi at every micro-vertex via `GetTrafo` (the FE
  analogue of the analytical flux-line trace inside an element), so the wire
  follows the curved order-p psi.  LAB Gx order 2: separate-turn
  `loops_homogeneity` 1.3e-4 -> 1.1e-4 and visibly smoother wires.
- **`--flux-plot out.png` (bubble-system flux-line view)** -- renders the
  designed coil's actual B field as flux lines on a cut-plane
  (`--flux-plane {x,y,z}`), seeded by the bubble system so the line density
  reflects `|B|`.  A physical check that the coil produces the intended field
  (e.g. the four-lobe saddle of a Gx gradient coil).

## Fusion (stellarator Stage-2) coil design -- mini-REGCOIL / NESCOIL

The stellarator coil-design community calls the surface stream function a
**winding-surface current potential** and uses it as the unknown of NESCOIL
(Merkel 1987), REGCOIL (Landreman 2017) and the surface part of FOCUS (Zhu
2018).  The "Stage-2" coil problem is mathematically the SAME object this
package solves for MRI gradient / induction-heating coils:

```
given a target normal field  B.n  on the PLASMA boundary,
find a current potential psi on a WINDING SURFACE around it whose
Biot-Savart  n.B  reproduces it          (iso-contours of psi = the coils)
```

`demo_regcoil_fusion.py` runs that forward problem with the existing surface-FE
stream function -- the rows of the design matrix are just the plasma-normal
component `n.B` of the winding-surface Biot-Savart kernel we already assemble.
Four results:

1. **Forward map is exact.**  Two producible targets -- a uniform vertical field
   (a PF / equilibrium / vertical-field coil) and a non-axisymmetric
   `sin(theta) cos(2 phi)` (a stellarator-like shape) -- are reproduced to
   `B.n` residual ~2e-9 (machine precision).  This is honest: the targets are
   smooth and within the winding surface's reach, so the forward map is exact,
   not hidden behind a tolerance.

2. **The REGCOIL trade-off.**  On a genuinely-hard target `sin(3 theta) cos(5
   phi)` (a high mode that decays across the plasma-coil gap, so it is NOT
   cheaply producible), sweeping the regularisation weight `alpha` traces the
   classic L-curve: large `alpha` -> smooth coil, high `B.n` residual; small
   `alpha` -> low residual but the peak surface current density `|grad psi|`
   saturates at the surface's representation limit (~7.9e6, the vertical leg of
   the L-curve), with a knee at `alpha_rel ~ 2e-2`.  `(field error, coil
   complexity)` is exactly the Stage-2 trade-off REGCOIL is built to expose.

3. **Net current = the multivalued / secular term.**  A single-valued `psi`
   carries zero net current through each hole of the winding torus, but a real
   coil set carries NET current.  The full current potential gains a SECULAR
   term `Psi = psi + (G/2pi) zeta + (I/2pi) theta` (`zeta` toroidal, `theta`
   poloidal angle); the TWO extra degrees of freedom are the first cohomology
   generators of the winding surface.  Their COUNT is the surface's first Betti
   number `b1 = 2` (genus 1), which the demo CONFIRMS gmsh-free via the Euler
   characteristic of the triangulated torus (`b1 = 2 - chi = 2g`); the volume
   T-Omega cohomology cut engine `src/radia/cohomology_cut.py` is itself
   gmsh-free (`radia.cohomology`).  On the torus
   the generators themselves are analytic (`grad(zeta)`, `grad(theta)` are
   single-valued vector fields), so the demo uses them directly and verifies the
   net-poloidal-current (TF) secular term gives the textbook toroidal field:
   `B_tor * R` is constant to ~0.2 % inside the tube (Ampere's `1/R`) and ~0
   outside (no enclosed poloidal current).  **The key physics**: the TF field is
   TANGENT to the plasma boundary, so its `B.n` footprint is >1000x smaller than
   the net-toroidal generator's -- i.e. the net poloidal current is NOT fitted
   from `B.n`, it is a PRESCRIBED engineering parameter (set it for the on-axis
   `B_tor`: 1 T at R=0.30 m needs 1.5 MA).  This is exactly why REGCOIL takes
   `net_poloidal_current` as an INPUT, not a fitted output.

4. **A VMEC-shaped plasma boundary.**  The circular plasma torus is replaced by
   a genuinely non-axisymmetric **rotating-ellipse** boundary in the VMEC Fourier
   representation `R = sum RBC(m,n) cos(m th - n NFP ph)`, `Z = sum ZBS(m,n)
   sin(...)` (NFP = 3), with analytic surface normals from the parametric
   tangents.  The vertical-field `B.n` target on this shaped boundary is
   reproduced to ~7e-9.  A real machine's free-boundary equilibrium is dropped in
   with `--wout wout_*.nc` (a standard VMEC output); the netCDF4 reader extracts
   the boundary `rmnc` / `zmns` coefficients and is round-trip-verified against
   the VMEC schema in `tests/panels/test_streamfunction_golden.py`.

**Honest scope** (what this is and is NOT): parts 1-2 use a SINGLE-VALUED `psi`
(correct for PF / RMP / shaping / shim fields); part 3 adds the multivalued
secular term for net-current (TF-type) coils -- the generator COUNT is computed
gmsh-free (Euler characteristic), and on the torus the generators are analytic
(exact).  For a GENERAL winding surface the generators come from the gmsh-free
cohomology cut (`src/radia/cohomology_cut.py`, `radia.cohomology`), not analytic forms.  Part 4's default boundary
is an analytic rotating-ellipse MODEL (the genuine VMEC Fourier shape, not a
converged equilibrium).  Three further capabilities are demonstrated in
`demo_regcoil_fusion_advanced.py`: (a) **coil force/stress** -- the Lorentz force
per area `f = K x B_avg`, verified to equal the magnetic pressure `B^2/(2 mu0)`
and to concentrate on the inboard TF leg; (b) designing against a **real
free-boundary equilibrium** -- the li383 (NCSX-like) `wout` from simsopt, with
`--wout` for any VMEC output; (c) a **FOCUS winding-standoff study** showing coil
complexity is monotonic in the gap, so the distance optimum is the minimum
engineering standoff; and (d) a **winding-SHAPE study** (FOCUS's core
contribution) -- a direct (theta,phi) surface mesher lets the winding be
CONFORMAL to the plasma, and for an elongated plasma a conformal winding (uniform
gap) cuts coil complexity by ~34% vs a circular one (varying gap) at the same
standoff.  The remaining next steps are a full Fourier-mode winding-surface
optimiser (the shape study sweeps a 1-parameter circular->conformal blend) and a
structural hoop-stress model (the demo reports the magnetic force per area, the
stress DRIVER, not the conductor stress itself).

## References

- Sugahara Lab, "ACA-accelerated stream function method + CMA-ES",
  IEEJ Joint Technical Meeting on Static Apparatus / Rotating Machinery,
  SA-25-020 (manuscript Method 2/3).
- P. Merkel, "Solution of stellarator boundary value problems with external
  currents", Nucl. Fusion 27, 867 (1987) -- NESCOIL.
- M. Landreman, "An improved current potential method for fast computation of
  stellarator coil shapes", Nucl. Fusion 57, 046003 (2017) -- REGCOIL.
- C. Zhu et al., "New method to design stellarator coils without the winding
  surface", Nucl. Fusion 58, 016008 (2018) -- FOCUS.
- HACApK (ppOpen-HPC, MIT): `src/ext/HACApK/`.
- Validation reference: the f2py-wrapped Fortran `coil_solver.f90`
  (`method_aca_tsvd_1/2`), matched bit-for-bit by
  `tests/test_stream_function.py::test_matches_f90_reference`.
