# BEM-A for Coil Inductance via Weggler Stabilized EFIE

This example demonstrates **vector-potential surface BEM (BEM-A)** for
the self-inductance of a closed-loop coil with a small gap, where the
gap defines a 2-port (source / sink) terminal pair.

The framework solves the same physical problem as filament PEEC, but
using a **surface mesh of the conductor** with a vector basis
(HDivSurface = RWG triangles).  This makes it the natural fallback
when the centerline of a "deformed coil" cannot be extracted (T-junction,
self-intersecting, Litz, asymmetric cross-section transitions).

## Method: Weggler stabilized EFIE = Lucy decomposition

The classical Maxwell EFIE on `HDivSurface` has a **low-frequency
breakdown** at κ → 0: the irrotational part of the surface current
appears as `V_κ / κ²`, which is unbounded as the frequency drops.

Lucy Weggler's stabilized EFIE [1, 2] solves this elegantly using a
**product-space variational formulation**:

```
HDivSurface(p) × SurfaceL2(p−1)
```

where the divergence operator `div : HDivSurface(p) → SurfaceL2(p−1)`
gives a clean de Rham complex.  The two spaces naturally separate
the solenoidal (closed-current) and irrotational (charge-accumulating)
components of the surface current — **mathematically equivalent to a
Loop-Star decomposition, but achieved at the variational level
without explicit basis construction**.

For coil self-inductance with prescribed port current `I = 1 A`, the
saddle-point system is:

```
[ SL    Dᵀ ] [ J ]   [ 0 ]
[ D     0  ] [ p ] = [ g ]
```

where:

- `J` is the surface current vector (HDivSurface DOFs, edge-based).
- `D : HDivSurface → SurfaceL2` is the divergence operator.
- `SL` is the LaplaceSL operator applied to the vector basis,
  i.e. `SL[i,j] = ∫∫ J_i(r) · J_j(r') / (4π|r − r'|) dS dS'` —
  the mutual-inductance kernel.
- `g` is the source/sink injection: `+1/A_src` on the source cap,
  `−1/A_snk` on the sink cap, zero elsewhere.  This enforces unit
  total current entering through `source` and exiting through `sink`.
- `p` is the Lagrange multiplier (one per SurfaceL2 face) enforcing
  current conservation on each element.

The self-inductance is then a quadratic form on the solved current:

```
L = μ₀ · J^T · SL · J        (in henries; uses I = 1 A)
```

Note that this includes the entire skin distribution captured at the
chosen `fes_order` — order 0 gives RWG (lowest order, piecewise-constant
divergence) and is sufficient for thin-wire-like geometries.

## Why we picked BEM-A over scalar BEM-SIBC

We initially tried two scalar magnetic-potential BEM-SIBC variants for
the same problem; both failed:

1. **Caps as cut surface + integrated jump** (`ψ_src − ψ_snk = −I·A_cap`
   with body Robin everywhere).  Numerically self-consistent but
   physically wrong: cap surfaces having SIBC closure forces `H_n = 0`
   on the cap, which contradicts the loop flux through the gap.
   Result: `L = −2.4 nH` on the gapped torus.
2. **Workpiece flow + centerline filament** (compute `ψ_inc` from a
   1-A filament along the coil's centerline; solve workpiece
   BEM-SIBC).  Computes "induced response on a passive workpiece in
   an external field", **not** "coil with prescribed current via
   topological cut".  Result: `L = −37 nH` (also wrong).

The root issue is that the scalar magnetic potential ψ on a multiply
connected domain is **multi-valued**; for a coil carrying current `I`,
ψ has a jump of `−I` across an Amperian-disk cut.  The Telegen energy
formula `W_mag = (μ₀/4) ∫ |∇ψ|² dV` then has an extra Σ contribution
that depends on `L` itself, making the extraction implicit and
gauge-dependent.

The vector formulation (BEM-A) sidesteps this entirely: the unknown is
the surface current `J`, the port BC is a Dirichlet on `div(J)` (clean
single-valued constraint), and `L = μ₀ J^T SL J` is a direct bilinear form.

## Validation against PEEC golden

`gapped_torus_demo.py` in this folder runs the gapped torus
(`R_maj = 30 mm`, `r_min = 3 mm`, `gap = 5°`) with `fes_order = 0`.
Expected output:

```
L = 86.67 nH  (PEEC golden 85.10 nH, rel err 1.8%)
```

The BEM-A value sits **between** the analytical Maxwell thin-wire
formula (89.80 nH, no gap correction) and the PEEC golden (85.10 nH,
which itself is ~5 % under analytical due to `n_peri = 16` perimeter
discretisation).

## Running

Requires NGSolve ≥ 6.2.2603 (with `ngsolve.bem` for the oracle path)
and `build123d ≥ 0.10` (for the loft demo).

```
# 1. Reference path: ngsolve.bem oracle on the gapped torus.
#    (Slow: ~100 s assembly.)
python gapped_torus_demo.py

# 2. End-to-end demo: build a lofted coil with build123d, then
#    cross-validate PEEC vs BEM-A on the same STEP.  This is the
#    Kubota hand-off entry point.
python loft_coil_peec_vs_bem.py
```

Expect ~100 s runtime at the default mesh density (`ngsolve.bem` is
slow at `N_J ≈ 5000`); see the perf notes below.

The reference solver lives at
`examples/induction_heating/bem_reference/bem_inductance.py
::compute_inductance_source_sink`.  When the **intree HACApK-backed**
assembler (Phase C / D.4) lands, `gapped_torus_demo.py` will switch
to the fast path automatically.

## Status (2026-05-02)

- [x] **Phase C.1-C.3**: intree RWG SS Galerkin assembler
      (`src/radia/bem/efie_rwg.py`) matches ngsolve.bem oracle to
      0.004 % on the gapped torus.
- [x] **Phase C.4a (DC SIBC)**: `compute_dc_resistance_intree()` solves
      the (Z_s · Mass)-based saddle to give DC R for arbitrary-shape
      conductors (R/Z_s = 10.36 vs uniform-J estimate 9.86 on gapped
      torus, ratio 1.05).  Targets PCB trace R extraction.
- [x] **Phase C.5 (cross-validation)**: build123d-generated lofted
      coil (R = 50 mm, 8 × 6 mm rect cross-section, 355 ° arc, 12 lofts)
      handled end-to-end::

          PEEC      L = 139.15 nH
          BEM-A     L = 153.56 nH      (diff +10.35 %)

      ~10 % spread is mesh- and modelling-dependent: PEEC uses
      `n_peri = 16` perimeter filaments at the thin-skin limit, BEM-A
      solves the EFIE on the coarse surface mesh in the PEC limit.
      Both fall inside the PEEC golden hard band [120, 160] nH and
      converge with refinement.  Run via
      ``python loft_coil_peec_vs_bem.py``.
- [ ] **Phase C.4b**: AC SIBC complex saddle
      ([j ω µ₀ SL + Z_s Mass     D^T] [J; p] = [0; g])
      for induction-heating Z_port = R + j ω L.
- [ ] **Phase C.6**: parallel + admissibility cutoff (Phase 1.10
      pattern from scalar SIBC, gives ~50-100 × intree-vs-ngsolve.bem
      assembly speedup once C++ port lands).
- [ ] **Phase C.7**: HACApK MatVec + GMRES (Phase 1.11 pattern).
- [ ] **Phase C.8**: production: `src/radia/bem_inductance_a.py`,
      `calc_coil_bem_a.py` CLI, `tests/panels/test_coil_bem_a_golden.py`.
- [ ] **Phase C++**: C++ port of the assembler (mirror
      `src/core/rad_bem_galerkin.cpp` scalar Phase 1.9).
- [ ] **Higher-order**: HDivSurface(p ≥ 1) for curved meshes.
- [ ] **Joachim forum post** on ngsolve.bem assembly slowness
      (HDivSurface 104 s @ N = 5064; same pattern as the scalar
      `LaplaceSL`/`LaplaceDL` path observed earlier).

## ngsolve.bem performance note

`ngsolve.bem.LaplaceSL` on `HDivSurface` exhibits the same
slow-assembly behaviour we observed in scalar BEM-SIBC (Phase 1.10):
~104 s assembly at `N_J = 5064`, with the matrix internally stored
as a 100 %-fill SparseMatrix and `H-matrix / FMM` flags effectively
ignored.  For one-shot validation this is acceptable; for production
we will either use intree HACApK assembly (the path scalar SIBC took
in Phase 1.10–1.11) or report the slowness upstream and wait for
fixes.

## References

[1] L. Weggler, *Maxwell DtN Stabilized — ngbem demo*,
    https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb
    (accessed 2026-02-23; archived in
    `to_developers/ngsolve/NGSolve-BEM_ stabilized Maxwell.pdf`)

[2] K. Sugahara, *ngbem for Low-Frequency Electromagnetics: Stabilized
    EFIE*, follow-up notebook at
    `to_developers/ngsolve/low_freq_efie_ngbem_applications.ipynb`
    (2026-02-24).

[3] R. Hiptmair, *Boundary Element Methods for Maxwell Transmission
    Problems in Lipschitz Domains*, Numer. Math., 2003 — original
    saddle-point formulation.

[4] G. Vecchi, *Loop-Star decomposition of basis functions in the
    discretization of the EFIE*, IEEE Trans. Antennas Propag., 47(2):
    339–346, 1999 — classical loop-star, mathematically equivalent
    to the Weggler product-space approach.
