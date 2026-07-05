# Validation strategy + benchmark catalogue

This page documents the validation strategy for the Radia SF framework
and the shipped benchmark under
[`validation_test/stream_function/benchmarks/`](../../validation_test/stream_function/benchmarks/).
The framework's feature set is described in
[README.md](README.md); this page is about VALIDATING those features
against published references and other tools.

## Shipped benchmarks

### `bench_helmholtz_pair.py` — analytical Maxwell pair baseline ✅

Compares our planar SF coil designing a uniform Bz over a small DSV vs
a classic Helmholtz pair (= two coaxial loops at radius `a` separated
by `a`).  The Helmholtz geometry is the gold-standard *coaxial* layout
for uniform Bz, so this benchmark CONTEXTUALISES our planar (single-
sided) solution.

Sample run:
```bash
cd validation_test/stream_function/benchmarks
python bench_helmholtz_pair.py --a 0.10 --dsv-r 0.025 --B0 0.001 \
    --json result_helmholtz.json
```

Output (sample):

| Metric          | Helmholtz pair  | Our planar SF     | Ratio (us / them) |
|-----------------|-----------------|-------------------|-------------------|
| RMS over DSV    | 0.13 %          | 0.89 %            | 6.7× WORSE        |
| p2p / mean      | 0.52 %          | 3.31 %            | 6.3× WORSE        |
| Wire length     | 1.26 m          | 25.63 m           | 20.4× LONGER      |
| # contours      | 2               | 36                | —                 |

**Interpretation**: a single-sided planar source has fundamentally
*less* geometric leverage than the coaxial Helmholtz pair for centred
uniform-Bz tasks.  The 6× ratio is the inherent geometry-class
penalty, not a flaw in our framework.

**For a fair like-for-like**: use biplanar source (top + bottom).
We can prototype that as `bench_biplanar_uniform.py` (~1 day).

## Future benchmark targets

The following five targets used to exist as TODO-only `.py` stubs. They are
now kept as documentation only until each target has a runnable implementation
and a JSON result. Each is 1-3 days of focused work.

### Bilac Planar Shim

**Reference**: Bilac et al. (Magn Reson Imaging, TBD year).  Planar
shim coil design for MRI; published target spec includes B0 + Z
gradient + ZZ second-order shim.

**Why this benchmark**: closest to our existing planar uniform Bz
demo; provides published numbers from a respected MRI shim coil
design paper.

**Next steps**:
  1. Extract paper PDF (W:\01_paper or W:\02_学会資料\…).
  2. Read target field formula, source plane dimensions, target DSV
     spec.
  3. Reproduce in our pipeline + JSON output.
  4. Cross-tabulate: published numbers vs ours.

### Turner Cylindrical Gz

**Reference**: Turner, R., *J. Phys. D: Appl. Phys.* 19, L147 (1986);
IEEE TMI 5 follow-up.  Cylindrical Gz gradient with analytical SFD.

**Why this benchmark**: Turner provides an ANALYTICAL reference SFD
(Bessel-function expansion).  This validates the SF design INDEPENDENT
of any other tool — we compare our discretised SFD directly to the
analytical formula.

**Next steps**:
  1. Code the Turner Bessel-function SFD analytically (cylinder radius
     `a`, length `L`, target gradient `G`, DSV radius `r`).
  2. Compute our SF design at the same parameters via the existing
     `demo_sf_to_peec_gz.py`.
  3. Compare ψ(φ, z) point-wise + compare resulting Bz on axis.
  4. Report RMS deviation from analytical reference.

### Lemdiasov-Ludwig 2005

**Reference**: Lemdiasov & Ludwig, *Concepts in Magnetic Resonance
Part B* 26B(1), 67-80 (2005). Target field method for MRI gradient
coils with detailed numerics.

**Why this benchmark**: published explicit target specs + reported
performance metrics for multiple coil designs.  Excellent
"second after Turner" benchmark.

### CoilGen Head-To-Head

**Reference**: Schwartz et al., CoilGen
(github.com/Philipp-MR/CoilGen) — the ONLY directly-comparable OSS
SF coil designer.

**Why this benchmark**: head-to-head OSS comparison.  CoilGen has years
of tuning on standard MRI gradient cases; running both on the same
spec and reporting comparable numbers is the most defensible
side-by-side validation we can show in a paper.

**Next steps**:
  1. Install CoilGen + MATLAB / GNU Octave on LAB.
  2. Pick a shipped CoilGen example (e.g., Gx gradient or shim).
  3. Reproduce same spec in our pipeline.
  4. Side-by-side table: RMS, wire length, inductance, compute time.

### Shielded Iron Yoke Material-Kernel Demo

**Reference**: no direct literature equivalent in the OSS space; this
demo extends our SF inverse-design pipeline to magnetic materials via
the kernel-agnostic callback contract.

**Why this benchmark**: the (A) callback contract lets us swap the
free-space Biot-Savart kernel for Radia MMM/MSC (iron yoke, permanent
magnet shield, SIBC workpiece).  Demonstrating this works on a shielded
coil (e.g., coil + iron back plate) makes a useful extension paper /
section because the SF design pipeline + material-kernel evaluation
have not been combined this way in prior open-source work that we are
aware of.  Worth verifying against any commercial-FEM equivalent we
can access before claiming originality.

**Next steps**:
  1. Define source plane + iron back plate (e.g., 5 mm steel plate at
     z = −5 mm below the source).
  2. Build Radia container with the iron yoke as `ObjRecMag` + BH
     curve via `MatSatIsoTab`.
  3. Replace the entry function with one that calls `rad.Solve()` +
     `rad.Fld()` for each (target, basis) pair.
  4. Per entry now ~100× slower than free-space → ACA+ saves the
     `M/k_aca` ratio (~25× at our M).
  5. Compare uniformity with vs without iron back plate.

Expected outcome: iron back plate mirrors the source toward the
target → uniformity improves at the same coil power.

## 2-week validation campaign plan

| Week | Task | Output |
|------|------|--------|
| 1, days 1–3 | Bilac + Turner + Lemdiasov-Ludwig | 3 JSON results + comparison table |
| 1, days 4–5 | CoilGen install + head-to-head | OSS comparison table |
| 2, days 1–4 | Shielded coil + Radia MMM kernel | Material-kernel demo |
| 2, days 5 | Paper draft (IEEE TMag / TMI) | Methods + Results draft |

After this campaign, every claim in the paper is backed by reproducible
numbers in version control and a manuscript draft can be assembled.

## Timing note

NGSolve 6.2.2604 was released 2026-04-30 with the `ngsolve.bem`
H-matrix bridge (Joachim Schöberl + Pierre Marchand).  Other groups
working in the SF + BEM space will likely pick up that infrastructure
over the next 6-12 months and arrive at overlapping capability, so it
is worth running the validation campaign sooner rather than later if
proof-of-priority on any specific contribution is intended.

## Cross-reference

  - Capability matrix: [README.md](README.md)
  - Paper outline: W:\02_学会資料\2025年度\2026_01_JIAM\streamfunction\paper_outline.md (moved out of repo)
  - MCP topic: `streamfunction(topic=session_2026_05_30)` section 11
