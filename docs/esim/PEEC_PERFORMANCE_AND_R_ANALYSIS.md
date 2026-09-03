# PEEC Coil Pipeline: Performance + R Discrepancy Deep Dive

**Audience.**  Authors / readers who want a deeper understanding of:
(1) why PEEC's coil R differs from BEM-A's R, and
(2) why STEP → filament conversion takes ~3 s per call (often the
dominant cost in a fast PEEC analysis).

**Companion docs.**
- [`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) — focused 1-page diagnosis.
- [`IMPLEMENTATION.md`](IMPLEMENTATION.md) — Karl loop + per-DOF mechanics.
- [`CROSS_VALIDATION.md`](CROSS_VALIDATION.md) — benchmark data.

---

## Part 1: Why PEEC R disagrees with BEM-A R

### 1.1 Re-stating the formulae

Both solvers report `R_coil` from the same vacuum coil, but compute it
via fundamentally different formulas (see
[`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) for the
short version):

**PEEC** — full complex impedance solve on a filament-bundle network:
- Per-filament cross-section discretisation: `n_peri` perimeter filaments + `nwinc × nhinc` interior grid.
- Per-filament DC resistance: `R_DC_k = ρ ℓ_k / A_k`.
- Per-filament AC resistance: `R_AC_k = R_DC_k · (R_ac/R_dc)_Dowell` where `(R_ac/R_dc)_Dowell = ξ · (sinh 2ξ + sin 2ξ) / (cosh 2ξ − cos 2ξ)` with `ξ = a_k / δ_coil` (filament radius / coil skin depth).  Implementation: [`esim_cell_problem.py:787-813`](../../src/radia/esim_cell_problem.py#L787-L813).
- Full `K × K` mutual inductance matrix `L_f`.
- Solve `Z_loop = R_f + jω L_f` as complex linear system + port-constraint augmentation.
- `R_coil = Re(V_port / I_port)`.

**BEM-A** — impedance-EFIE (Z_s inside the saddle; sole formulation
since 2026-07-02):
- Solve the COMPLEX saddle `[[jω μ0 SL + Z_s M, D^T], [D, 0]] [J; p] = [0; g]`
  with `Z_s = (1+j)/(σ δ_coil)`; J is the finite-impedance surface
  current.
- `R_coil = Re(Z_s) · (Jᴴ M J)`, `L_coil = μ0 · (Jᴴ SL J)` (external).
  Implementation: [`coil_inductance_ngsolve.py`](../../src/radia/bem/coil_inductance_ngsolve.py).
- The historical perfect-conductor saddle + post-hoc
  `R = Re(Z_s)·JᵀMJ` was REMOVED 2026-07-02: on tightly-wound coils
  the PEC J concentrates singularly at near-contact gaps/edges and
  over-estimated R ~3× (kubota 3-turn 15.14 vs 4.63 mΩ).

### 1.2 What each captures

| Effect | PEEC | BEM-A (impedance-EFIE) |
|---|---|---|
| Ohmic loss in skin layer (`½ Re(Z_s) \|J\|²` integrated) | YES (via Bessel/Dowell on each filament) | YES (Leontovich on the finite-impedance J) |
| Azimuthal J redistribution around the cross-section (proximity) | Partial (perimeter filaments, ~1.2× ceiling) | YES (surface J resolves it) |
| Suppression of the spurious PEC edge/near-contact singularity | n/a | YES (the resistive Z_s M term penalises it) |
| Weak-skin / DC limit R → R_DC | YES | NO (Leontovich SIBC invalid; frequency=0 gives vacuum L, R=0) |

The historical claim that BEM-A's limit is a "lower bound on the true
AC resistance" was **wrong** for the removed PEC path: 2026-07-02
showed it OVER-estimates ~3× on tightly-wound coils (SIBC breakdown at
near-contact J singularities).  The shipped impedance-EFIE agrees with
PEEC / volume PEEC / analytic proximity on those coils (4.63 vs
4.5-4.8 mΩ) and with the Bessel closed form on smooth wires (<1 %).

### 1.3 Empirical data point: gapped-torus benchmark

Production benchmark
([`validation_test/ih_esim_benchmark/results.json`](../../validation_test/ih_esim_benchmark/results.json)),
gapped torus 1 turn + Cu, 50 kHz:

| Source | R_coil_mOhm |
|---|---|
| PEEC (n_peri=16) | 0.233 |
| FEM A-V volumetric (mesh-resolved coil skin) | 0.255 (computed as `2 × (P_coil - P_wp) / I²` from results.json sweep[1].fem_coilmesh: P_total = 1.27e-4, P_coil = 6.08e-5 → R_coil ≈ 0.243 mΩ) |
| BEM-A impedance-EFIE (validation fixture, 7 kHz) | 0.223 mΩ (= Bessel R_ac at 7 kHz within ~7 %) |

PEEC is **within ~10 %** of the volumetric FEM A-V reference for this
geometry.  The historical "BEM-A likely 0.10-0.15 mΩ (lower bound)"
estimate in earlier versions of this table described the removed PEC
post-hoc path and was wrong in both value and direction; the shipped
impedance-EFIE matches the Bessel/PEEC physics (see
[`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) § 2-3).

### 1.4 Complex Leontovich SIBC in the EFIE saddle — SHIPPED 2026-07-02

The proposal below SHIPPED (and is now the sole BEM-A formulation):

```
[jω μ0 SL + Z_s M     D^T]  [J]   [0]
[D                     0 ]  [p] = [g]
```

with `Z_s = (1+j)/(σ δ)` on the coil surface.  The shipped extraction
convention differs from the original proposal in one deliberate way:

```
R_coil = Re(Z_s) · (Jᴴ M J)          # Leontovich dissipation
L_coil = μ0 · (Jᴴ SL J)              # EXTERNAL inductance
```

i.e. the internal surface reactance `Im(Z_s)·(Jᴴ M J)/ω` is NOT folded
into L — it regularizes J inside the solve, but the reported L stays
the pure geometry quantity every downstream consumer and golden
expects (folding it in shifted L by +5.9 % on the 7 kHz gapped torus
for no benefit to the R fix).  Golden:
`validation_test/bem/test_coil_bem_a_impedance_efie.py`.

### 1.5 What this means for the IGTE paper

The paper uses **PEEC's R** as its primary number for L_total + R_total
in the workpiece-coupled case.  BEM-A's R is **not used in the paper's
headline numbers** — it appears only in the `--coil-solver bem-a`
optional path which is for users with non-circular cross-sections
that PEEC's perimeter-filament model cannot resolve.

The honest statement for the paper:

> Coil AC resistance is computed via the PEEC filament bundle's full
> complex impedance solve (Eq. X), giving R_coil within ~10 % of
> volumetric FEM A-V reference for the Cu test case at 50 kHz.  The
> BEM-A coil path uses a perturbative SIBC dissipation integral that
> under-estimates R_coil; this is a known limitation, addressed in
> v4.56 by including the full complex Leontovich SIBC in the saddle
> system.

---

## Part 2: Why STEP → filament conversion is slow (~3 s)

### 2.1 Where the 3 seconds go

The benchmark reports `t_coil_topology_s ≈ 3 s` consistently across
all four test frequencies
([`results.json`](../../validation_test/ih_esim_benchmark/results.json) lines
38, 255, 465, 668).  This is **constant** across frequencies — the
topology extraction is frequency-independent (no warm/cold cache
variance).

Inside [`coil_from_cad.py:filaments_from_step`](../../src/radia/coil_from_cad.py),
the work breaks down approximately as:

| Sub-step | Time | Implementation |
|---|---|---|
| STEP load via build123d / OCC | ~0.5-1.0 s | `import_step(step_path)` |
| Centerline extraction (5 predicates) | ~1.0-1.5 s | `extract_centerline(...)` |
| Cross-section discretisation + filament placement | ~0.3-0.5 s | `coil.to_filaments_peri(...)` |
| PEEC topology + L matrix build (Ruehli kernel) | ~0.1-0.2 s | `build_bundle_solver(...)` |

Total: ~3 s for a typical 16-filament gapped-torus STEP.

### 2.2 Sub-step 1: STEP load

The OCC kernel must parse the B-Rep, build the topology graph, and
construct the solid.  For a simple gapped torus this is ~30 % of the
total `t_coil_topology_s`.

**Optimisation potential**: limited.  STEP parsing is bounded by OCC's
internal data-structure construction; pre-converting to a faster
format (e.g. cached pickle of the topology graph) would help
multi-run scenarios but adds a new caching mechanism to maintain.

### 2.3 Sub-step 2: Centerline extraction

This is the largest single cost.  `extract_centerline` dispatches
across five positive-match predicates (multi-station loft / united
multi-turn / revolution+plane / OPEN / CLOSED) to determine the coil
topology class, then applies the appropriate centerline-recovery
algorithm.

For the canonical IGTE-benchmark coil
([`ih_fem_kelvin_demo_coil.step`](../../src/radia/panels/samples/ih_fem_kelvin_demo_coil.step)),
the dispatch lands on Path 1 (lateral surface UV sampling), which
samples the BSPLINE lateral face at `n_stations × n_peri` (u, v) grid
points.  For 50 stations × 16 perimeter points = 800 OCC surface
evaluations, each ~1-2 ms → ~1-1.5 s.

**Optimisation potential**:
1. **Parallelise UV sampling** — embarrassingly parallel; each station's
   `n_peri` evaluations are independent.  Expected 2-4× speedup on
   multi-core via `concurrent.futures.ThreadPoolExecutor`.
2. **Cache the centerline across frequency sweeps** — see § 2.5.

### 2.4 Sub-step 3: Cross-section discretisation

`CoilBuilder.to_filaments_peri(n_peri)` walks the centerline and
samples the cross-section at each station, placing `n_peri` filament
nodes around the perimeter.  Cost is `O(n_stations × n_peri)`
geometric-primitive evaluations.

Already fast (~0.3-0.5 s) and dominated by Python-side bookkeeping
rather than numeric work.

### 2.5 Caching strategy (recommended for frequency sweeps)

The benchmark's frequency sweep at 4 points spends
`4 × 3 s = 12 s` re-extracting the same coil topology — wasted because
the topology is frequency-independent.

**Proposed cache layer**:

```python
# Pseudo-code, calc_inductance.py refactor
class CoilTopologyCache:
    def __init__(self, step_path, n_peri, ...):
        self.step_path = step_path
        self.n_peri = n_peri
        self._cached_topo = None

    @property
    def topology(self):
        if self._cached_topo is None:
            self._cached_topo = filaments_from_step(
                self.step_path, n_peri=self.n_peri, ...)
        return self._cached_topo

# In the frequency-sweep driver:
cache = CoilTopologyCache(args.coil_step, args.peec_n_peri, ...)
for f in frequencies:
    topo = cache.topology    # ← reused across f, only computed once
    ...solve at f...
```

**Expected impact**: 4-frequency sweep `12 s` → `3 s` topology cost
(75 % reduction).  For 10-point sweeps the savings grow proportionally.

### 2.6 Single-shot optimisation

For one-off `calc_inductance.py` calls (not sweeps), the 3-s overhead
is amortised against the BIE solve (~5 s for typical wp mesh).  Not a
major concern.

The user-facing impact is most visible when:
- The Cubit panel triggers a re-run after a small parameter change
  (e.g. `--current`) — currently re-extracts the topology each time.
- A frequency sweep is being run.
- A multi-config benchmark is being run.

### 2.7 What this means for the IGTE paper

This is **NOT a paper contribution**, it's an implementation-quality
note.  Mention briefly in § 7 reproducibility / wall-time discussion:

> Total runtime per PEEC-BEM evaluation is ~5 s wall, of which ~3 s
> is STEP → filament topology extraction (constant per call,
> independent of frequency).  Frequency sweeps benefit from topology
> caching (planned, v4.56 roadmap).

---

## Part 3: Action items

| Item | Severity | Effort | Tracked |
|---|---|---|---|
| Full complex Leontovich SIBC in BEM-A saddle (align R with PEEC) | Medium | ~1 week | v4.56 roadmap |
| Topology cache for frequency sweeps (75 % cost reduction) | Low | ~3 days | v4.56 roadmap |
| Parallelise UV sampling in centerline extraction (2-4× speedup) | Low | ~2 days | v4.57 roadmap |
| Document the BEM-A R lower-bound nature in user-facing CLI help | Low | ~30 min | next docs pass |

---

**Document version**: 2026-05-30 (radia v4.67.0+ dense-sweep baseline).
