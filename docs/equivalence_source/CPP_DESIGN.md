# C++ design — Equivalence-theorem near-field source

**Status**: design review draft, 2026-05-26 (Plan D from chat)
**Owner**: Sugahara Lab (Kindai University)
**Tracks issue**: Python implementation of `radia.equivalence_source` is too slow for production use (O(M · N) per evaluation, single-threaded numpy). Low-frequency static and high-frequency time-harmonic kernels are mathematically different and must both be implemented in C++ via NGSolve's facilities.

This document records the architectural decisions BEFORE Phase A begins so they can be reviewed in one pass. The implementation itself follows in Phase A-E (`src/core/rad_equivalence_source.h/cpp` + NGSolve C++ binding + pybind11).

---

## 1. Goals

1. **Speed**: 100x faster than the v1 Python Stratton-Chu integrator for typical N=10^3 to 10^4 surface panels and M=10^3 to 10^5 observation points.

2. **Two distinct kernels in one module**:
   - **Static / quasi-static (omega == 0)**: `1/R` Laplace Green's function. Real-valued. Used for magnetostatic Radia problems (~99% of lab use cases).
   - **Time-harmonic (omega != 0)**: `exp(-jkR) / (4πR)` with full **dyadic** Green's function `(I + ∇∇/k^2)` for accurate near-field reconstruction. Complex-valued. Used for WPT, EMC, antenna problems.

3. **NGSolve-native integration**: a `NearFieldSourceCoefficientFunction` C++ class that NGSolve's `GridFunction.Set()` can sample directly, enabling proper L2 projection (not the v1 nodal-sampling hack).

4. **HACApK-accelerable**: pluggable kernel evaluator so the v3 acceleration (H-matrix, ACA+ compression, O((M+N) log N)) drops in transparently. FMM is intentionally NOT pursued -- see §10.

5. **Backward-compatible Python API**: existing `NearFieldSource` user code keeps working unchanged. Internal Stratton-Chu loop gets replaced with C++ kernel calls in Phase A; once C++ is verified Phase E removes the Python loop.

## 2. Non-goals (intentionally excluded)

- **IABC / SDI / asymptotic-BC content** -- Radia's Kelvin transformation is the lab default for FEM open boundary. The equivalence-theorem reconstruction is BC-insensitive (Sugahara Lab 2015 empirical result), so any IABC variant is redundant.
- **FMM acceleration** -- not pursued in this design. See §10 for the analysis.
- **Re-implementing the entire FEM solver** -- this is a post-processing layer.
- **Visualisation** -- field probing via the CF interface is sufficient; CST-style 3D rendering is GMSH's job (CLAUDE.md "NGSolve-Through Principle").

## 3. Kernel mathematics

### 3.1 Magnetostatic reduction (omega == 0)

Sources on a closed surface ∂Ω enclosing all primary magnetisation:

```
J_s(r')   = n̂(r') × H(r')                  [A/m]
ρ_m(r')   = μ_0 (n̂(r') · H(r'))            [Wb/m²]
```

Note that even though E = 0, ρ_m is NOT zero (this was the trap in the 2008 Sugahara Lab axisym slide 2 -- see `radia_mcp.fem.fem_equivalence_source('schelkunoff_love')`).

Reconstruction at observation point r outside ∂Ω, scalar-potential + vector-potential form, then differentiate inside:

```
H(r) = (1 / 4π) ∮ [ ∇(1/R) × J_s(r')  -  (n̂(r') · H(r')) ∇(1/R) ] dS(r')

       where R = |r - r'|,  ∇(1/R)|_r = -(r - r') / R³
```

This is exactly what `evaluate_static_H` in the v1 Python module computes, verified to 0.13% error vs analytic dipole on a 6240-panel sphere.

**Cost analysis (direct)**: per observation point r, O(N) operations:
- N face traversals
- 3 cross products + 3 dot products + 1 vector add per face

For M observation points: O(M · N). At M = N = 10^4, that's 10^8 ops — fine in C++ (~50 ms), unacceptable in Python (~5 s).

### 3.2 Time-harmonic, full dyadic GF (omega != 0)

Add the time-harmonic surface sources:

```
M_s(r')   = n̂(r') × E(r')                  [V/m]
```

Reconstruction via the dyadic Green's function `Ḡ_e(r,r') = (I + ∇∇/k²) ψ(r,r')`:

```
ψ(r,r')  = exp(-jkR) / (4πR),   k = ω √(μ_0 ε_0)

E(r) = ∮ [ -jωμ_0 Ḡ_e · J_s  +  ∇ψ × M_s ] dS'
H(r) = ∮ [  jωε_0 Ḡ_e · M_s  +  ∇ψ × J_s ] dS'
```

The dyadic Green-function term already contains the longitudinal
surface-charge contribution; Radia's production C++ kernel does not add
separate `(n̂·E)∇ψ` or `(n̂·H)∇ψ` terms.

The `Ḡ_e · J_s` term expands to:

```
Ḡ_e · J_s = ψ J_s + (1/k²) ∇(∇ · (ψ J_s))
          = ψ J_s + (1/k²) ∇( (∇ψ) · J_s )
```

The second term is what's MISSING from the v1 Python scalar form and what causes the 66% undershoot at R_obs / λ << 1 in Phase 2.

**Working out the gradient** (NGSolve / numpy vectorisable form, per face contribution at obs point r):

```
R_vec    = r - r'
R        = |R_vec|
ψ        = exp(-jkR) / (4π R)
grad_ψ   = -ψ · (jk + 1/R) · R_vec / R                                 // [C³]
H1       = -(jk + 1/R) / R                                              // scalar
H2       = (jk(R) + 1)*(jk(R) + 3) / R² + k²                            // helper
grad_grad_ψ · J_s
         = ψ · [ H2 · (R_vec · J_s / R²) · R_vec / R² 
                 + H1 · J_s  ]                                          // [C³]
```

(The exact expression for `∇∇ψ : v` for a vector v has 4 terms when fully expanded; the C++ implementation will use a compact, numerically-stable form -- see Jackson 3e eq. 9.45, or Balanis Advanced EM Eng eq. 7-23.)

**Important static limit check**: as ω → 0, k → 0, ψ → 1/(4πR), grad_ψ → -R_vec/(4πR³). The `1/k²` factor in front of the dyadic term diverges, but the `∇∇ψ` term has structure that cancels: it becomes `-3 R_vec ⊗ R_vec / R^5 + I / R³`, multiplied by `(1/k²)` it's NOT divergent if the integrand involves J_s appropriately. The continuity relation `jω ρ_e + ∇·J_s = 0` connects ρ_e and J_s so that the divergent parts cancel in the ω → 0 limit and the surviving scalar formula reduces to the magnetostatic case.

The Phase B implementation must verify this limit numerically: at ω = 1 Hz the harmonic kernel must give within 1% of the static kernel.

### 3.3 Why two kernels and not one?

- **Numerical stability**: at ω = 0 the harmonic form has 1/k² · 0/0 cancellations. The static-specific kernel sidesteps this entirely.
- **Speed**: the static kernel is 4-5x faster per face (no complex arithmetic, fewer terms). For Radia's 99% magnetostatic use case, this matters.
- **Code clarity**: two named kernels in C++ are easier to maintain than one templated mess.

## 4. File layout

```
src/core/
  rad_equivalence_source.h               # Public API (NearFieldSource class)
  rad_equivalence_source.cpp             # Constructors, save/load, static kernel
  rad_equivalence_source_dyadic.cpp      # Harmonic dyadic kernel (Phase B)
  rad_equivalence_source_hacapk.h/.cpp   # H-matrix acceleration (Phase D)

src/lib/
  rad_equivalence_source_pybind.cpp      # pybind11 bindings — mirrors v1 Python API

# NGSolve-coupled module (Phase C):
# Uses NGSolve's add_ngsolve_python_module pattern (same as sparsesolv_ngsolve)
src/lib/ngsolve/
  equivalence_source_ngsolve.cpp         # NearFieldSourceCoefficientFunction

# Python (Phase A: shim; Phase E: pure facade):
src/radia/
  equivalence_source.py                   # API surface, delegates to C++
  equivalence_source_ngsolve.pyd          # built by Build.ps1 -> radia wheel
```

### 4.1 `NearFieldSource` C++ class skeleton

```cpp
// rad_equivalence_source.h
namespace radia {

class NearFieldSource {
public:
    // Surface data — flat row-major (CLAUDE.md row-major policy)
    std::vector<double>                m_centroids;  // 3N: [x0,y0,z0, x1,y1,z1, ...]
    std::vector<double>                m_normals;    // 3N: outward unit normals
    std::vector<double>                m_areas;      // N
    std::vector<std::complex<double>>  m_E;          // 3N or empty (static)
    std::vector<std::complex<double>>  m_H;          // 3N
    double                             m_omega;

    NearFieldSource() = default;
    NearFieldSource(int n_faces, double omega, bool with_E);

    // ---- Static kernel (Phase A) -----------------------------------------
    // Compute H at M observation points (real-valued).
    // out_H is M*3, row-major.  Vectorised over (M, N), parallelised over M
    // via NGSolve TaskManager (CLAUDE.md TaskManager policy).
    void EvaluateStaticH(
        const double* obs_points, int M,
        double* out_H) const;

    // ---- Harmonic kernel (Phase B) ---------------------------------------
    void EvaluateHarmonic(
        const double* obs_points, int M, double omega,
        std::complex<double>* out_E,    // M*3 (may be NULL to skip E)
        std::complex<double>* out_H) const;

    // ---- HACApK kernel (Phase D) ------------------------------------------
    // Pre-build the H-matrix; then EvaluateStaticH etc. routes through it
    // when m_hacapk != nullptr.  See rad_equivalence_source_hacapk.h.
    void BuildHMatrix(double aca_tol = 1e-4,
                      int leaf_size = 10,
                      double eta = 2.0);

    // ---- Serialization ----------------------------------------------------
    void SaveJson(const std::string& path) const;
    static NearFieldSource LoadJson(const std::string& path);
};

}  // namespace radia
```

### 4.2 Python API contract (kept stable)

```python
from radia.equivalence_source import NearFieldSource

# Construction (unchanged)
nfs = NearFieldSource.extract_ngsolve(mesh, gf_E=None, gf_H=H_cf,
                                        surface_label="nfs_surface", omega=0)

# Probing — internally delegates to C++ kernels
H_static = nfs.evaluate_static_H(obs_points)           # → C++ static kernel
E, H = nfs.evaluate(obs_points, omega=2*pi*1e6)        # → C++ harmonic dyadic kernel

# NGSolve CF wrap (NEW in Phase C):
cf_H = nfs.as_coefficient_function(component='H')      # native NGSolve CF
gfu.Set(cf_H, definedon=mesh.Materials("outer_air"))   # proper L2 projection

# HACApK acceleration (NEW in Phase D):
nfs.build_hmatrix(aca_tol=1e-4)
H_fast = nfs.evaluate_static_H(obs_points)             # automatic, no API change
```

## 5. NGSolve CoefficientFunction wrapper (Phase C)

The v1 nodal projection (`gfu.vec[v*3+i] = H_at_vertex[v,i]`) is sub-optimal for order >= 2 H1 spaces. The NGSolve-native solution: implement a `CoefficientFunction` subclass that the FE projection machinery calls at each quadrature point.

```cpp
class NearFieldSourceCoefficientFunction
    : public CoefficientFunction
{
    const NearFieldSource& m_nfs;
    int                    m_field;   // 0 = H, 1 = E
public:
    NearFieldSourceCoefficientFunction(const NearFieldSource& nfs, int field)
        : CoefficientFunction(/*dim=*/3, /*complex=*/(nfs.m_omega != 0)),
          m_nfs(nfs), m_field(field) {}

    // Required override -- evaluate at one IntegrationPoint
    void Evaluate(const BaseMappedIntegrationPoint& mip,
                   FlatVector<> values) const override;

    // Vectorised override (NGSolve's batch interface)
    void Evaluate(const BaseMappedIntegrationRule& mir,
                   BareSliceMatrix<> values) const override;
};
```

The vectorised `Evaluate` over an IntegrationRule calls `EvaluateStaticH` (or `EvaluateHarmonic`) on the IR's quadrature-point physical coordinates in one batch -- this is where the C++ TaskManager parallelism pays off.

Python-side: `nfs.as_coefficient_function(component='H')` constructs the CF and returns it. NGSolve's `gfu.Set(cf, ...)` then runs the canonical L2 projection.

## 6. Acceleration via NGSolve.bem ML (Phase D) — see `FMM_DESIGN.md`

**Decision (2026-05-26)**: acceleration of the one-shot evaluator goes
through **NGSolve.bem** Multilevel Expansion (FMM-equivalent), NOT
HACApK and NOT a Radia-vendored FMM library.  Rationale:

1. The one-shot evaluator has **no matrix to recompress** — HACApK's
   strength (ACA+ on a re-used matrix during iterative solve) does
   not apply here.  HACApK remains the right choice for MMM/MSC
   interaction matrices (per CLAUDE.md "Use HACApK Only" policy);
   that is a separate use case.

2. NGSolve.bem 6.2.2603 already ships the full FMM stack:
   `BiotSavartRegularMLCF`, `BiotSavartSingularMLCF`,
   `MaxwellSingleLayer/DoubleLayerPotentialOperator(Curl)`,
   `HelmholtzSingleLayer/DoubleLayerPotentialOperator`, etc.
   These cover all Phase A + Phase B kernels.

3. Adding our own FMM library (ExaFMM-t etc.) into `src/ext/` would
   duplicate NGSolve.bem and violate CLAUDE.md "Complement NGSolve".

**Break-even**: `N_face × N_obs > 10⁹` (≈5 s direct walltime on 8
cores).  Today's examples top out at ~10⁶ pairs — direct C++ wins
below that, so deferring is rational.

**Frequency-dependent FMM path** (per `FMM_DESIGN.md` §3.3): standard
Helmholtz FMM has the well-known **low-frequency breakdown** at
`kR ≪ 1` (Greengard-Huang 2002).  Radia's typical low-frequency
regime (IH 10 kHz, R ≈ 1 m → `kR ~ 10⁻⁵`) falls into the breakdown
territory.  Mitigation: at low frequency, route through the
**Laplace ML path** (no breakdown) and add the small imaginary part
via the direct C++ kernel — the Phase B static-limit test (ω=1 Hz
matches Phase A to 5e-15) justifies this hybrid scheme.
Mid- and high-frequency (`kR ≥ 0.01`, WPT MHz and above) use the
standard Helmholtz / Maxwell ML operators.

**Delivery plan**: see [`FMM_DESIGN.md`](FMM_DESIGN.md) §4 (D1–D5).
Summary: ~1 week of glue Python around NGSolve.bem primitives, no
new Radia C++ kernel.  Implementation deferred until a concrete user
case crosses the break-even.

## 7. Build integration

**`Build.ps1`** changes (Phase A):
- Add `src/core/rad_equivalence_source.cpp` to the radia static-lib target.
- Add `src/lib/rad_equivalence_source_pybind.cpp` to the `_radia_pybind.pyd` target (extends the existing module).

**`Build.ps1`** changes (Phase C):
- Add `add_ngsolve_python_module(equivalence_source_ngsolve ...)` targeting `src/lib/ngsolve/equivalence_source_ngsolve.cpp`. This produces `src/radia/equivalence_source_ngsolve.pyd`.

**Phase D** acceleration: no Build.ps1 / CMake change required — Phase D
adds only Python glue around NGSolve.bem (already a runtime
dependency).  No HACApK link, no FMM library vendor.

**`pyproject.toml` package_data** (Phase C):
- Add `equivalence_source_ngsolve.pyd` to the manifest so the wheel includes it.

## 8. Verification & benchmarks

### 8.1 Regression suite (every phase must keep all green)

- `validation_test/equivalence_source/phase1_static_coil.py` — 0.83% PASS.
- `validation_test/equivalence_source/phase2_wpt_harmonic.py` — 0.12% PASS for E and nonzero H; zero-H observations pass the absolute threshold.
- `validation_test/equivalence_source/phase3_e2e_cubit_to_sol.py` — 18.26% PASS on the 2026-06-28 Cubit `.vol` path (20% band); Phase C should tighten to ~1%.

### 8.2 New benchmarks (introduced with each phase)

- **Phase A**: `validation_test/equivalence_source/bench_static.py` — measure ms/eval for N ∈ {10², 10³, 10⁴, 10⁵} faces, M = 100 obs. Numerical equality and production-scale speed are hard gates.
- **Phase B**: planned `validation_test/equivalence_source/bench_harmonic.py` — same shape, complex arithmetic.
- **Phase C**: `tests/equivalence_source/test_cf_projection.py` — pytest that compares `gfu.Set(nfs.as_coefficient_function())` against the v1 nodal projection on a known dipole; native CF must be within 1% of analytic at obs points.
- **Phase D**: planned `validation_test/equivalence_source/bench_hacapk.py` — scaling study, M=N up to 10⁵. Acceptance: < 10× slower than direct C++ at N = 10³, ≥ 100× faster than direct C++ at N = 10⁵.

### 8.3 Numerical fixtures

Reference data from the 2015_04_12 Femtet directory (Sugahara Lab benchmark grid) lands in `tests/equivalence_source/fixtures/femtet_2015_reference/` (TODO Phase A or B). The grid covers:

- 10 MHz Hertzian dipole, Dirichlet + Neumann + radiation + IABC inner BCs
- 300 MHz Hertzian dipole, same matrix
- WPT parallel plates, 1 MHz, shield on/off

These give us a real-world acceptance band per BC, per geometry, per frequency.

## 9. Migration plan: kill Python Stratton-Chu

| Phase | What | Python kept? |
|---|---|---|
| A | C++ static kernel + Python shim delegates | Python static loop kept as fallback path, gated by `_USE_CPP = True` env var. Default ON. |
| B | C++ harmonic dyadic kernel | Python harmonic loop kept as fallback. |
| C | NGSolve CF wrapper, `as_coefficient_function()` exposed | Python paths still callable. |
| D | HACApK fast path | Direct C++ still callable. |
| **E** | **Delete Python Stratton-Chu integration code** | **Python `evaluate_static_H` / `evaluate()` become pure facades that always call C++.** |

By Phase E the diff is: remove the inner loops in `equivalence_source.py`, keep `NearFieldSource` dataclass + constructors + serialisation as Python (they were always pure-Python convenience anyway). All numerics flow through C++.

## 10. Why FMM is intentionally skipped

CLAUDE.md "FMM (Fast Multipole Method): Removed (2026-03-06)" documents the empirical conclusion that for Radia's MMM/MSC use cases:

1. Dipole approximation accuracy is poor for distributed surface sources
2. Compact geometries have ≥87% near-field pairs → HACApK is 10-100× faster
3. Direct C++ with TaskManager handles N < 10⁴ at production speed
4. HACApK covers all large-scale needs

Equivalence-theorem reconstruction has the same characteristics: closed surface (compact), distributed face sources (not point dipoles), N typically 10³-10⁴. HACApK is the natural acceleration. FMM would add complexity (tree maintenance, M2L translation, P2P near-field) for marginal benefit.

If a future use case demands N > 10⁶ (multi-billion-DOF FEM with massive surface), revisit. Until then, HACApK is the ceiling.

## 11. Open questions (Phase A pre-flight)

- [ ] Does NGSolve's TaskManager work inside a Python-callable C++ extension built via pybind11, or only inside `add_ngsolve_python_module` targets? (If the former is no, the static kernel needs `add_ngsolve_python_module` from Phase A, not just Phase C.)
- [ ] Does the existing `_radia_pybind.pyd` already link the NGSolve TaskManager symbol? If yes, Phase A is straightforward. If no, we may need to split the equivalence-source module out from day one.
- [ ] Serialisation format compatibility: should the new C++ writer produce a binary `.nfs` for speed (~50× smaller files), or keep the JSON format for portability? Default decision: KEEP JSON in Phase A (don't break the existing `.nfs.json` artifacts), add binary `.nfs.bin` as an optional Phase E enhancement.
- [ ] For Phase D, should we use the same HACApK leaf size / ACA tolerance as Radia's MMM (eps=1e-4, leaf=10, eta=2.0), or do equivalence-theorem-specific defaults need exploration?

## 12. References

- **Schelkunoff 1936** — original equivalence theorem
- **Love 1901** — earlier closed-surface equivalence
- **Stratton 1941** — full vector formulation, _Electromagnetic Theory_ §8.14
- **Jackson 1999** — _Classical Electrodynamics_ 3e §10.2 (Kirchhoff integral)
- **Balanis 2012** — _Advanced Engineering Electromagnetics_ §7-3 (dyadic Green's function)
- **Harrington 2001** — _Time-Harmonic Electromagnetic Fields_ §3-5 (Schelkunoff form for radiating fields)
- **Sugahara Lab 2008** axisym foundations (S:\FEMM\等価定理の基礎原理\軸対称\)
- **Sugahara Lab 2015** open-domain solution thesis (S:\FEMM\2015_05_21_等価定理\)
- **Sugahara Lab 2015** FEMTET IABC + reconstruction reports (S:\99_調査済\Femtet\2015_04_12_IABC_recnstruct\)
- **HACApK_LH-Cimplm** library (src/ext/HACApK_LH-Cimplm/) — H-matrix acceleration
- Companion MCP tool: `fem_equivalence_source(topic="overview|schelkunoff_love|...")` (`radia_mcp.fem.equivalence_source_knowledge`)

## 13. Open issue tracker (post-design)

- v4.X.X tag plan once Phase A merges: `equivalence_source-cpp-static`
- v4.Y.Y after Phase B: `equivalence_source-cpp-harmonic`
- v4.Z.Z after Phase D: `equivalence_source-hacapk`
- Python removal lives in the same release as Phase E (no separate tag).
