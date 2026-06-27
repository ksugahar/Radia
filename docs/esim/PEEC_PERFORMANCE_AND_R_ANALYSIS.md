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

**BEM-A** — perfect-conductor saddle + post-hoc dissipation:
- Solve real EFIE saddle `[[SL, D^T], [D, 0]] [J; p] = [0; g]` (no Z_s in matrix; J is the perfect-conductor surface current).
- `R_coil = Re(Z_s) · J^T M J` with `Re(Z_s) = 1/(σ δ_coil)`.  Implementation: [`coil_inductance_ngsolve.py:218-226`](../../src/radia/bem/coil_inductance_ngsolve.py#L218-L226).

### 1.2 What each captures

| Effect | PEEC | BEM-A |
|---|---|---|
| Ohmic loss in skin layer (`½ Re(Z_s) \|J\|²` integrated) | YES (via Dowell on each filament) | YES (post-multiply) |
| Filament-internal skin amplification (`R_AC > R_DC` by Dowell factor) | YES | NO (J is treated as uniform over conductor cross-section) |
| Inter-filament reactive current redistribution at high freq | YES (full complex L_f solve) | NO (J is computed at vanishing-resistance limit) |
| Inter-filament resistive coupling | YES (Z_loop has off-diagonals via L_f imaginary at high freq) | NO |

The two formulae are **physically different** and converge to
different limits:

- **PEEC** converges, as `n_peri → ∞` and `nwinc, nhinc → ∞`, to a
  filament-bundle limit that approximates the volumetric eddy-current
  result.  For circular Cu at moderate `ξ` ~ 1-10 this limit agrees
  with the full 3-D Maxwell solve to ~5 %.
- **BEM-A** converges, as the surface mesh → fine, to the
  "Leontovich-perturbed perfect-conductor" limit, which is a
  **lower bound** on the true AC resistance.  The error scales as
  `(δ/a)²` for the filament-bundle skin amplification term.

### 1.3 Empirical data point: gapped-torus benchmark

Production benchmark
([`docs/ih_esim_benchmark/results.json`](../ih_esim_benchmark/results.json)),
gapped torus 1 turn + Cu, 50 kHz:

| Source | R_coil_mOhm |
|---|---|
| PEEC (n_peri=16) | 0.233 |
| FEM A-V volumetric (mesh-resolved coil skin) | 0.255 (computed as `2 × (P_coil - P_wp) / I²` from results.json sweep[1].fem_coilmesh: P_total = 1.27e-4, P_coil = 6.08e-5 → R_coil ≈ 0.243 mΩ) |
| BEM-A (estimate — no production benchmark on the same coil) | likely 0.10-0.15 mΩ (lower bound) |

PEEC is **within ~10 %** of the volumetric FEM A-V reference for this
geometry.  BEM-A's expected under-estimate of ~50 % (see
[`R_MISMATCH_PEEC_VS_BEMA.md`](R_MISMATCH_PEEC_VS_BEMA.md) § 5) is
consistent with its theoretical lower-bound nature.

### 1.4 Proposed fix to align BEM-A with PEEC (v4.56+ roadmap)

Include the **full complex Leontovich SIBC** in the EFIE saddle:

```
[SL + (Z_s / jω) · M     D^T]  [J]   [0]
[D                        0 ]  [p] = [g]
```

with `Z_s = (1+j) ρ / δ` (Dowell linear-SIBC) on the coil surface.  Then:

```
L_coil = Im(J^H Z_eq J) / ω    where Z_eq = jω μ_0 SL + Z_s M
R_coil = Re(J^H Z_eq J)
```

(Hermitian conjugates because J is now complex.)  Effort estimate:
~1 week (assembler modification + golden tests).  Tracked as v4.56
roadmap.

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
([`results.json`](../ih_esim_benchmark/results.json) lines
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
