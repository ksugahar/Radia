"""
Cauer Ladder Network (CLN) — the lab-specialty MOR method for
eddy-current FEM problems. Sugahara Lab is a primary contributor
to the canonical CLN literature.
"""

CLN_OVERVIEW = """
# Cauer Ladder Network (CLN) — MOR for eddy-current FEM

The Cauer Ladder Network method (Kameari-Ebrahimi-Sugahara-Shindo-Matsuo
2018) is a Model Order Reduction technique that replaces a large eddy-
current FEM system by an equivalent electric LADDER NETWORK with a
small number of stages, while preserving the impedance response.

## The big picture

```
Eddy-current FEM (millions of DoFs)
        ↓ (MOR via CLN)
Cauer ladder circuit (N stages, e.g. N=10-60)
        ↓ (can be solved instantly for any input signal)
Time-domain response (transient simulation in real time)
```

A finite-element transient eddy-current simulation that takes ~1-7 HOURS
for a single frequency point becomes a Cauer network simulation that
takes ~MILLISECONDS for any input waveform.

This is essential for:
- Real-time control of converters/inverters
- Online power-loss calculation
- Hardware-in-the-loop simulation of electrical machines
- System-level co-simulation (SPICE coupling)

## The Sugahara Lab's role

Sugahara is a primary author/co-author on the foundational CLN papers:
- 2018a: Kameari-Ebrahimi-**Sugahara**-Shindo-Matsuo
  "Cauer Ladder Network Representation of Eddy-Current Fields for
   Model Order Reduction Using Finite-Element Method"
  IEEE Trans. Mag. 54(3):7201804
- 2018b: **Sugahara**-Kameari-Ebrahimi-Shindo-Matsuo
  "Finite-Element Analysis of Unbounded Eddy-Current Problems Using
   Cauer Ladder Network Method"
  IEEE Trans. Mag. 54(3):7200704
- 2018c: Matsuo-Kameari-**Sugahara**-Shindo
  "Matrix Formulation of the Cauer Ladder Network Method for
   Efficient Eddy-Current Analysis"
  IEEE Trans. Mag. 54(11):7205805
- 2019: Kuriyama-Kameari-Ebrahimi-Fujiwara-**Sugahara**-Shindo-Matsuo
  "Cauer Ladder Network With Multiple Expansion Points for Efficient
   Model Order Reduction of Eddy-Current Field"
  IEEE Trans. Mag. 55(6):7203404
- 2020: Ebrahimi-**Sugahara**-Matsuo-Kaimori-Kameari
  "Modal Decomposition of 3-D Quasi-Static Maxwell Equations by
   Cauer Ladder Network Representation"
  IEEE Trans. Mag. 56(3):7513004

This makes radia-mcp.mor THE canonical MCP knowledge base for CLN.

## Connection to radia core

`radia/lanczos_reduction.py` contains a PRIMA / Lanczos implementation
for PEEC circuits.  PRIMA and CLN are sibling MOR methods (both use
Krylov/Cauer expansion).  The CLN framework is the natural extension
of the PEEC PRIMA work to FEM eddy-current problems.
"""

CLN_BASIC_RECURSION = """
# Basic CLN recursion (A-formulation)

## Setup

Eddy current problem in time-harmonic form (1 conductor + air):
    curl(ν curl A) + jω σ A = j_0   in Ω
with appropriate boundary conditions.

Space-time separation:
    E(t, x) = Σ_n v_n(t) E_n(x)
    H(t, x) = Σ_n i_n(t) H_n(x)

where {E_n} and {H_n} are orthogonal basis fields:
    ∫_{Ω_c} σ E_n · E_m dx = 0   for m ≠ n
    ∫_Ω    μ H_n · H_m dx = 0   for m ≠ n

The time coefficients v_n(t), i_n(t) are the node voltages and
branch currents of a CAUER LADDER CIRCUIT (next section).

## Recursive computation of basis fields

Starting from E_{-1} = H_{-1} = 0 and L_{-1} = 0:

For n = 0, 1, 2, ...:
    curl(E_n − E_{n−1}) = −(L_{n−1})⁻¹ B_{n−1}    (electric step)
    curl(H_n − H_{n−1}) = R_n J_n                 (magnetic step)

    R_n = [∫_Ω σ E_n · E_n dx]⁻¹
    L_n = ∫_Ω μ H_n · H_n dx

Each step is a STATIC FEM problem (no time dependence) → solvable
fast with standard ICCG/PARDISO.

## The resulting Cauer circuit

    vT(t) ──[R_0]──┬──[R_1]──┬─ ... ──[R_n]──
                   |          |
                  [L_0]      [L_1]    [L_n]
                   |          |          |
                   └──────────┴─────────┘

Total impedance Z(jω) computed via continued fraction:
    Z(jω) = R_0 + 1 / (jωL_0 + 1/(R_1 + 1/(jωL_1 + ...)))

After truncating at stage N, the impedance has the EXACT eddy-current
behavior at expansion point + asymptotic correctness.

## Convergence properties

- Resistance R: monotonically convergent
- Inductance L: monotonically convergent
- For N ≥ 30-60, error < 1% over the entire eddy-current frequency
  range (kHz to MHz for typical industrial inductors)

## Where to truncate (Köster-König-Birò 2021)

Use the IMPEDANCE INCREMENT as stopping criterion:
    |Z_n(f) − Z_{n−1}(f)| / |Z_n(f)| < tolerance

This does NOT need a reference FEM impedance value — purely
internal convergence diagnostic.
"""

CLN_MULTIPLE_EXPANSION = """
# Multiple expansion points (Kuriyama et al 2019)

## Limitation of single-expansion CLN

The basic CLN expands around frequency f_0 = 0 (DC limit).  Near f_0,
convergence is fast (few stages).  Far from f_0 (e.g. at hundreds of
kHz for a problem expanded at DC), convergence is slow — may need
N > 100 stages.

## The fix: change expansion point at stage M+1

Instead of always expanding around f_0 = 0, change to expansion point
f_1 (e.g. f_1 = the target operating frequency) at the M+1-th stage:

    Stages 0 to M:    expand around f_0 (typically 0)
    Stages M+1 to N:  expand around f_1 (operating frequency)

The orthogonality conditions need extra terms p_{2i-1} to maintain
the basis structure:
    p_{2i-1} = -1/λ_{2i-1} · a_{2M-1}^T K_0 a_{2i-1}

And one additional inductor L* appears in the circuit:
    L* = -λ_{2M-1} / [(s_1 - s_0) λ_{2M-2} (λ_{2M-1} + p_{2M-3}λ_{2M-3})]

## Results (Kuriyama 2019)

For a wireless power transfer (WPT) system at 85 kHz operating
frequency:
- Single expansion at f_0 = 0: not accurate even after 100 stages
- Single expansion at f_0 = 85 kHz: accurate around 85 kHz with few stages
- Two expansion points (0, 85 kHz), change at stage 2: ACCURATE from
  DC to 1 MHz with only 6 stages

Speedup vs. direct FEM:
- FEM: 31 min for 13 frequency points
- CLN (6 stages): 2.5 min to generate, then trivial cost per frequency

## When to use multiple expansion points

- The signal/load has a known nominal operating frequency
- The impedance Z(f) must be accurate over a wide frequency range
- Computational efficiency matters (real-time / HIL simulation)

## A vs T formulations

Same recursion structure can be applied to:
- A-formulation: vector potential, conductive + non-conductive domain
- T-formulation: current vector potential, conductive domain only

T-formulation gives a simpler circuit but limited applicability
(needs to be extended to T-Ω for general use, per Kuriyama 2019
conclusion).
"""

CLN_NONLINEAR = """
# CLN for nonlinear materials (Sato-Shimotani-Igarashi 2017)

The basic CLN assumes linear materials (μ const, σ const).  For
ferromagnetic media with B-H saturation, the recursion (5) must be
adapted.

## Approach: linearize around an operating point

For a given excitation level, identify the operating point in
B-H space.  Use the LOCAL slope μ_local = dB/dH at that point.
Build the CLN with μ_local.

## Limitations

- Each operating point requires its own CLN
- For varying excitation, multiple CLNs need to be combined
- Hysteresis is not naturally captured

## Extension to hysteresis (Sato-Clemens-Igarashi 2016)

The "adaptive subdomain MOR with discrete empirical interpolation"
method handles nonlinear magneto-quasi-static problems.  Idea:
identify the saturated regions adaptively and use POD / discrete
empirical interpolation on each subdomain.

This is more general than pure CLN but loses the closed-form ladder
circuit interpretation.

## Connection to radia.hantila_solver

`radia/hantila_solver.py` implements Hantila polarization method for
nonlinear MMM.  The Hantila method also splits the constitutive
relation into linear + residual parts and factors the linear
operator once.  This is the same "factor once, iterate on residual"
philosophy as CLN.

A natural future direction: combine CLN-style modal expansion with
Hantila-style nonlinear residual iteration.  This would give a
NONLINEAR CLN that handles full B-H curves without ad-hoc local
linearization.
"""

CLN_APPLICATIONS = """
# CLN applications

## 1. Industrial inductors (Köster-König-Birò 2021)

Geometry: 2 cylindrical ferrite cores with 12 air gaps, aluminum
windings, height 260 mm.  Operating frequency 1 kHz - 1 MHz.

Direct FEM: 1-7 hours per frequency point.
CLN (N=60 stages, mesh m3.2): ~5 min generation, then trivial.
Accuracy: ΔR < 3%, ΔL < 0.3% at 750 kHz.

Different meshes are optimal for different goals:
- Mesh for IMPEDANCE: coarse, just enough to capture conductor topology
- Mesh for FIELD reconstruction: must resolve skin depth at the
  target frequencies

A "winding mesh" (fine in conductor, coarse outside) is the best
trade-off in 2021 author's experience.

## 2. Wireless Power Transfer (Kuriyama 2019)

Two coils (primary 10 turns, secondary open) at 85 kHz operating
frequency.  Multi-expansion-point CLN with f_0=0 and f_1=85 kHz
gives accurate impedance from DC to 1 MHz with only 6 stages.

Useful for system-level WPT design where the magnetic system
appears as a circuit element in a larger SPICE-like simulation.

## 3. Hybrid Twin (digital twin with CLN backend)

Recent work (~2022-2024) combines CLN with machine-learning-
enhanced reduced order models for magnetic bearings.  The CLN
provides the physics-based skeleton; ML corrects the residual.

## 4. Connection to PEEC

PEEC = Partial Element Equivalent Circuit (Radia/Loop-Star method).
Both PEEC and CLN produce circuit-level abstractions of EM problems,
but they differ in:
- PEEC: spatially-distributed circuit elements derived from geometry
  (no MOR initially; MOR via PRIMA-Lanczos comes later)
- CLN: time-evolution-modal expansion, lumped circuit elements

For LARGE problems (e.g. industrial-scale induction heating),
combining PEEC-coil + CLN-workpiece can give the best of both.

## Implementation hooks in radia ecosystem

- `radia.lanczos_reduction`: PRIMA/Lanczos MOR for PEEC
- `radia.sparsesolv_ngsolve`: Compact HX/AMS preconditioner for the
  large H(curl) systems that arise in CLN basis-field computation
- (Future) `radia.cln`: a Python wrapper for the CLN method that
  ties NGSolve FEM to a Cauer circuit output, ready for SPICE
  consumption
"""


def get_cln_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"         - CLN big picture, Sugahara Lab's role
      "recursion"        - Basic A-formulation recursion + Cauer circuit
      "multiple"         - Multiple expansion points (Kuriyama 2019)
      "nonlinear"        - Extension to nonlinear ferromagnetic materials
      "applications"     - Industrial inductors, WPT, hybrid twin
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return CLN_OVERVIEW
    if topic in ("recursion", "basic", "a_formulation"):
        return CLN_BASIC_RECURSION
    if topic in ("multiple", "multi_expansion", "expansion_points"):
        return CLN_MULTIPLE_EXPANSION
    if topic in ("nonlinear", "hysteresis", "saturation"):
        return CLN_NONLINEAR
    if topic in ("applications", "examples", "inductor", "wpt"):
        return CLN_APPLICATIONS
    if topic == "all":
        return "\n\n".join([
            CLN_OVERVIEW, CLN_BASIC_RECURSION, CLN_MULTIPLE_EXPANSION,
            CLN_NONLINEAR, CLN_APPLICATIONS,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, recursion, multiple, nonlinear, applications."
    )
