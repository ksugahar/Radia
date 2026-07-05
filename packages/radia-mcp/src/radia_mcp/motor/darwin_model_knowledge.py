"""Darwin-model time-domain knowledge for radia_mcp.motor.

Reference: Kaimori-Mifune-Kameari-Wakao, "Time-Domain Finite-Element
Analysis of the Darwin Model with Coulomb-type Gauge", in public-safe curated corpus
論文関係/2025年度/136_劉馨遙/モーター/.

The Darwin model is the **A-φ formulation including displacement
current** but discarding electromagnetic-wave propagation — strictly
intermediate between the magnetostatic limit (MQS, used by the
standard machine-magstadyn_a model) and the full Maxwell wave
equations.  It captures:

- inductive coupling (∂A/∂t terms, as in MQS)
- capacitive coupling (∇φ + ∂A/∂t in the E field, as in EQS)
- conductive currents

but neglects radiation (Lorenz wave correction).  This puts the
Darwin model in the same physics regime as the **PEEC** approach used
elsewhere in radia_mcp — both are valid for sizes ≪ λ.

For motors at PWM frequencies (1-20 kHz), the Darwin model is the
right framework when you need to capture **parasitic capacitive
currents** in the stator winding (shaft voltage, common-mode EMI),
which the standard A-formulation MQS cannot describe.
"""

from __future__ import annotations


DARWIN_OVERVIEW = """\
## Darwin model — what it captures and what it discards

The full Maxwell system is:

  ∇ × H  =  J + ∂D/∂t       (Ampere)
  ∇ × E  =  -∂B/∂t          (Faraday)
  ∇ · D  =  ρ
  ∇ · B  =  0

with constitutive laws  B = μH ,  D = εE ,  J = σE + J_s.

Three approximation hierarchies in low-frequency regime:

| Model | Drops | Captures | Regime |
|-------|-------|----------|--------|
| Magnetostatic (MS)  | ∂D/∂t, ∂B/∂t, σ E | static B from J_s, ρ_b | DC only |
| Magneto-Quasi-Static (MQS) | ∂D/∂t | + eddy currents (σ ∂A/∂t) | electric-conductors-dominated |
| Electro-Quasi-Static (EQS) | ∂B/∂t | + displacement, no induction | dielectric-dominated |
| **Darwin** | **none of these** | + inductive + capacitive | both inductive & capacitive |
| Full Maxwell | nothing | radiation | full wave |

The Darwin "drop" is the **Lorenz back-reaction term** that completes
the wave equation.  In a Coulomb gauge, the Darwin equation becomes
**elliptic-in-space, first-order-in-time** — much cheaper than the
full hyperbolic wave equation.

### When you need Darwin (not just MQS)

- PWM common-mode currents through coil-to-frame parasitic
  capacitance.
- Shaft voltage estimation on motor rotor.
- Cable + motor system EMI simulation.
- High-frequency winding analysis where coil-to-coil capacitance
  matters.

For pure low-speed torque waveform analysis, **MQS is sufficient**
and you should stay with `machine_magstadyn_a.pro`.
"""

DARWIN_FORMULATION = """\
## Kaimori-Mifune-Kameari-Wakao Coulomb-gauge Darwin TD formulation

Choose A, φ such that:

  B  =  ∇ × A
  E  = -∇φ - ∂A/∂t

Substitute into Ampere with displacement current:

  ∇ × (ν ∇ × A)  =  σ(-∇φ - ∂A/∂t)  +  J_s  -  ε(∇ ∂φ/∂t + ∂²A/∂t²)

Darwin = **drop ∂²A/∂t²** (the wave second-time derivative on A).
The displacement current ε∇∂φ/∂t **is retained** — it gives the
capacitive coupling.

Apply ∇·  to obtain the φ equation (a charge-conservation form):

  ∇ · (-σ∇φ - σ∂A/∂t + J_s - ε∇∂φ/∂t)  =  0

### Variant 1 — no extra DOF

Solve the coupled (A, φ) system directly.  The matrix is **non-
symmetric** because the cross-terms (σ ∇φ in A-eq, σ ∂A/∂t in φ-eq)
are not transposes of each other.  Solver: GMRES with block
preconditioner (CompactAMS for HCurl block, AMG for H1 block).

### Variant 2 — with redundant variable χ (symmetric form)

Introduce auxiliary `χ` with constraint `χ = ∂φ/∂t` (i.e. χ is the
time-rate of φ).  The system becomes:

  M_σ ∂A/∂t + K_ν A + G_σ φ + G_ε χ  =  f
  G_σ^T A + L_σ φ + L_ε χ            =  0
  ε ∂φ/∂t                            =  ε χ     (definition)

Now the (A, φ, χ) system is **symmetric** — because χ appears
symmetrically as the time-rate of φ.  Solver: PARDISO direct
(real symmetric indefinite) or block-CG with symmetric
preconditioner.

### Numerical verification

The paper verifies the formulation on three test problems:

1. **Cylindrical conductor in uniform external H** at frequency f.
   Compares to the analytical Bessel-function solution.  Both
   variants converge at the expected O(h²) rate; Variant 2 (symmetric)
   is ~30% faster per iteration thanks to the symmetric solve.
2. **LC parallel circuit** modeled as a 3D coil + parallel-plate
   capacitor.  Recovers the analytical resonance frequency
   ω_0 = 1/√(LC) to within 0.5%.
3. **PCB spiral inductor** with substrate capacitance.  Compares to
   a measurement-derived S-parameter file — matches to within 2 dB up
   to 100 MHz (beyond which radiation, which Darwin discards, matters).

### Implications for motor analysis

Variant 2 is the practical choice for motor PWM-EMI studies.  The
symmetric solve allows the same iterative-solver infrastructure as
the MQS A-formulation (just adds a φ block + a χ block).  In the
NGSolve framework:

```python
fes = FESpace([
    HCurl(mesh, order=2),   # A
    H1(mesh, order=2),       # φ
    H1(mesh, order=2),       # χ
])
```

with the corresponding bilinear form blocks.  The symmetric
preconditioner is `CompactAMS` for the (HCurl) A-block, BoomerAMG for
the H1 φ- and χ-blocks, with off-diagonal coupling handled by an
inner GMRES on the Schur complement (well-conditioned because of the
symmetry).
"""

DARWIN_VS_PEEC = """\
## Darwin TD vs PEEC — when to use which

Both Darwin and PEEC capture **inductive + capacitive + resistive**
coupling without radiation.  Comparison:

| Aspect | Darwin (FE) | PEEC (IE) |
|--------|-------------|-----------|
| Discretization | Volumetric tet/hex mesh | Surface mesh + filaments |
| Open boundary | Kelvin / PML / truncation | Natural (Green's function) |
| Nonlinear materials | Easy (ν(B), σ(T)) | Hard (would need iteration) |
| Skin effect | Resolved by mesh (expensive) | SIBC built in (cheap) |
| Coupling to circuit | Per-conductor DOF | Native (port-based) |
| Suitable for | Iron-cored motors, transformers | Air-cored coils, PCB |

For motor PWM-EMI:
- The iron stator core requires Darwin (volumetric, with ν(B)).
- The coil-to-frame parasitic capacitance can be **lumped** into
  external circuit (no need for Darwin in 3D) → MQS + circuit
  is often enough.
- The **shaft voltage** problem genuinely needs Darwin TD if you
  want to predict shaft voltage waveform from CAD geometry without
  measurement-derived capacitances.

For air-cored coils + PCB:
- PEEC is the default (radia_mcp.peec).

### A practical recipe for motor + circuit

1. Use MQS (`machine_magstadyn_a.pro`) for the main torque & loss
   simulation.
2. Extract per-conductor lumped capacitance to frame from a separate
   **electrostatic** solve (no current, no time) — this is just an
   EQS step.
3. Add the lumped capacitances to the external circuit (SPICE-style).
4. Re-run MQS + circuit; the SPICE block now handles the high-frequency
   capacitive coupling.

If step 2 is inadequate (e.g. capacitance depends on B-saturation of
the iron, which is non-linear and time-varying), only then go full
Darwin TD per Kaimori-Wakao 2024.
"""

SECTIONS = {
    "overview": DARWIN_OVERVIEW,
    "formulation": DARWIN_FORMULATION,
    "vs_peec": DARWIN_VS_PEEC,
}


def get_darwin_knowledge(topic: str = "overview") -> str:
    """Return Darwin-model time-domain knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
