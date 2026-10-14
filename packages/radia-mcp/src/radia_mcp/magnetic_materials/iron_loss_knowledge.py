"""Iron loss models: Steinmetz, Bertotti 3-term, Pry-Bean, MSE, iGSE,
Carstensen, anomalous loss.  public-safe curated corpus + adjacent.
"""

STEINMETZ_FAMILY = r"""
# Steinmetz family of iron loss formulas

## Original Steinmetz (1892)
Empirical:  P_iron = k_h * f * B^n,   n ≈ 1.6 - 2.0  (the "Steinmetz
exponent"; ~1.6 for grain-oriented steel, ~2.0 for non-oriented).

VALIDITY: sinusoidal B excitation only.  Non-sinusoidal (e.g., PWM,
short-pulse, switched mode) underestimates loss significantly.

## Modified Steinmetz Equation (MSE, Reinert 2001)

For non-sinusoidal B(t):
```
P_MSE = k_h * f_eq^(n-1) * f * (B_peak)^n
f_eq = (2 / (B_peak² * pi²)) * integral over T of (dB/dt)² dt
```

f_eq is the "equivalent frequency" matching the rms dB/dt of the
actual waveform.  Useful for triangular, trapezoidal, PWM-flat-top.

## Improved Generalized Steinmetz Equation (iGSE, Venkatachalam 2002)

Better than MSE for non-symmetric waveforms (DC bias, asymmetric PWM):
```
P_iGSE = (1/T) * integral over T of  k_i * |dB/dt|^alpha
                                     * (delta_B)^(beta - alpha) dt
k_i = k_h / ((2π)^(alpha-1) * integral_0^(2π) |cos(theta)|^alpha
            * 2^(beta-alpha) dtheta)
```
where (alpha, beta, k_h) are Steinmetz parameters fitted to sinusoidal
data.

Captures: DC bias, minor loops, non-sinusoidal trapezoid.

## When to use which

| Excitation | Use |
|-----------|-----|
| Sinusoidal AC only | Steinmetz original |
| PWM with flat-top, triangular | MSE |
| DC bias + AC, asymmetric | iGSE |
| Minor loops dominant | iGSE or Preisach-based loss |
| Square wave | iGSE with care |

## Reference papers
- Reinert et al., IEEE TIA 37(4):1055, 2001 (MSE)
- Venkatachalam et al., IEEE COMPEL 2002 (iGSE)
- Sullivan et al., Power Electronics Society (review)

## Lab papers (public-safe curated corpus)
- インバータ回路解析と磁気解析の併用計算による...半導体特性の影響まで...
- 電磁鋼板の詳細な積層モデルを用いた高周波用リアクトルの鉄損解析
- 高周波励磁下におけるリアクトル損失の空隙特性
- Iron Loss Comparison between Reactor with Air Gap and Material
- Cauer 等価回路とプレイモデルを用いた...鉄損特性評価
"""


BERTOTTI_3TERM = r"""
# Bertotti 3-term iron loss decomposition

## Formulation (Bertotti 1988)

Iron loss density per cycle:
```
P_iron = P_hysteresis + P_classical_eddy + P_anomalous
       = k_h f B^n + k_c (B f)^2 + k_a (B f)^(1.5)
```

- **P_hysteresis = k_h f B^n**: irreversible domain-wall displacement
  losses (rate-independent, scaled by f because it happens once per cycle)
- **P_classical_eddy = k_c (B f)^2**: macroscopic eddy current loss
  (Faraday-induction in the bulk material)
- **P_anomalous = k_a (B f)^1.5**: excess loss from local eddy currents
  near domain walls (the "anomalous" contribution)

The 1.5 exponent for anomalous loss comes from Bertotti's statistical
domain-wall theory and is well-supported empirically for grain-oriented
silicon steel.

## Statistical foundation (Bertotti 1998 monograph "Hysteresis in
## Magnetism", Ch 12 "Eddy Currents")

The 3-term split is NOT three independent physical mechanisms — it
EMERGES from the statistical average of a random sequence of
**Barkhausen jumps**.  Each jump i induces a local eddy-current
density j_i; the total dissipation ∝ ⟨j²⟩ where j = Σ_i j_i:
```
⟨j²⟩ = N⟨j_i²⟩ + N⟨j_i' j_j⟩        (Bertotti eq 12.6, N = # jumps)
```
- The **diagonal term** N⟨j_i²⟩ (jumps treated independently) →
  **hysteresis loss** P_h, ∝ f (rate-independent loop area × f).
  Comes from the FINE-scale eddy currents around each individual
  jump; lumped into the static B(H) constitutive law.
- The **cross term** N⟨j_i' j_j⟩ (jump correlations) splits further:
  - statistically INDEPENDENT part → **classical loss** P_cl
    (domains averaged away, material = homogeneous medium)
  - CORRELATED part (a jump raises the probability of the next) →
    **excess loss** P_exc (intermediate domain scale)

Loss separation works precisely BECAUSE randomness dominates — it
fails for a single 180° wall (picture-frame crystal) or small-
amplitude oscillation around fixed walls (no statistics).

### Classical loss — exact, material-law-independent (Ch 12.2)

Slab thickness d, conductivity σ.  Maxwell across the slab gives
the instantaneous eddy current and the **classical field**:
```
j(y,t) = σ y (dB/dt)               (parabolic across thickness)
H_cl(t) = (σ d² / 12) (dB/dt)       (eq 12.17, counter-field)
P_cl(t) = (σ d² / 12) (dB/dt)²      (eq 12.13, instantaneous)
```
Time-averaged over a cycle:
```
P_cl = (σ d² / 3) B_max² f²         (triangular B, eq 12.14)
P_cl = (π²/6) σ d² B_max² f²        (sinusoidal B, eq 12.15)
```
This is the LOWEST possible loss for given (B_max, f) — any other
waveform has higher classical loss (⟨(dB/dt)²⟩ ≥ ⟨dB/dt⟩²).  This
is WHY iron-loss standards mandate controlled sinusoidal B.
**Material-law-independent** at low frequency — applies to ANY steel.

Concrete (Bertotti): Si-Fe 50 Hz 1.5 T, σ=2e6 S/m, d=3e-4 m,
dB/dt = 4 B_max f = 300 T/s → H_cl ≈ 5 A/m, already comparable to
the coercive field of grain-oriented steel.

### Skin-effect crossover (Ch 12.2.2, linear B = μH)

The diffusion equation introduces the dimensionless parameter
```
γ ∝ √(σ μ d² f)      (skin depth = d/γ)
```
- γ ≪ 1 (low freq): P_cl ∝ f²  (uniform field, no shielding)
- γ ≫ 1 (high freq): P_cl ∝ f^1.5  (skin-limited, field expelled
  from slab center)
- Crossover at γ ≈ 3, i.e. f ≈ 9/(σ μ d²).  Above this, the
  low-freq formula (12.14) over-predicts and one must solve the
  full diffusion equation.

### Excess loss — statistical magnetic objects (Ch 12.3-12.4)

Excess loss comes from the intermediate scale (magnetic domains).
Bertotti's statistical theory models the field in excess of the
classical + hysteresis background as driving `n(t)` simultaneously
active correlation regions ("magnetic objects", MOs):
```
P_exc = √(σ G V_0 S) ⟨ |dB/dt|^1.5 ⟩          (per cycle)
       = 8 √(σ G V_0 S) (B_max f)^1.5          (sinusoidal)
```
- `G ≈ 0.1356` — dimensionless coefficient (eddy-current damping of
  a moving wall, from the parabolic field profile)
- `V_0` — a phenomenological field (A/m) characterizing the
  statistical distribution of local coercive fields (the spread of
  MO switching strengths)
- `S` — cross-sectional area, `n_0` ∝ √(σ S / V_0) sets how many MOs
  are active simultaneously
- The **1.5 exponent** (√f) is the signature of this √(dB/dt)
  statistical scaling — distinct from classical f² and hysteresis f.

This `V_0`/`n_0`/`G` machinery is what `k_a` in the engineering
3-term formula compresses into a single fitted constant:
`k_a = 8 √(σ G V_0 S)` (per unit volume → divide by S).

## When the statistical picture matters

For routine motor/transformer loss, the fitted 3-term `(k_h, k_c,
k_a)` table above is sufficient.  The Ch 12 statistical foundation
matters when: (a) extrapolating loss outside the measured (B, f)
range — the f / f² / f^1.5 scalings are physically grounded, not
just curve-fits; (b) connecting to the energy-based hysteresis
model (the same Barkhausen-jump pinning picture underlies
`lab_core_production`'s κ pinning field); (c) the skin-effect
crossover invalidates the low-freq P_cl above f ≈ 9/(σμd²).

## Parameter extraction

For commercial NO (non-oriented) and GO (grain-oriented) steel:

| Grade | k_h | n | k_c | k_a |
|-------|-----|---|-----|-----|
| 35JN230 (NO, 0.35 mm) | 0.012 | 1.7 | 4.6e-5 | 0.0001 |
| 35JN300 (NO, 0.35 mm) | 0.018 | 1.7 | 4.6e-5 | 0.00015 |
| 50JN800 (NO, 0.50 mm) | 0.035 | 1.7 | 1.3e-4 | 0.00025 |
| 27ZH95 (GO, 0.27 mm)  | 0.005 | 1.6 | 2.8e-5 | 0.00005 |

(approximate; check manufacturer datasheet for exact values; lab data
in `public-safe curated corpus`)

## Frequency / waveform dependence
- Bertotti 3-term assumes sinusoidal — for PWM, multiply each term by
  appropriate iGSE-style correction
- High frequency (>1 kHz): classical eddy dominates; anomalous may
  saturate
- Low frequency (<100 Hz): hysteresis dominates; classical / anomalous
  small

## Reference
- G. Bertotti, "General properties of power losses in soft
  ferromagnetic materials", IEEE TMAG 24(1):621, 1988
- Lab textbook: `public-safe curated corpus properties of power
  losses in soft ferromagnetic material.pdf` (21 MB)
- Bertotti, "Hysteresis in Magnetism" 1998 (579 MB scan; the canonical
  monograph)

## Pry-Bean classical limit

For thin sheet (thickness d, conductivity σ, sinusoidal B_peak at f):
```
P_classical_per_volume = (pi * σ * d² * f² * B_peak²) / 6
```
This is the analytical limit of Bertotti's k_c when local domain
structure is uniform.  Use as sanity check on k_c calibration.

## Relation to MSFEM-Hollaus effective material

Hollaus's MSFEM (motor_hollaus_eddy) computes the effective complex
μ_eff which has loss tangent matching Bertotti k_c (eddy) and k_h
(hysteresis, if play model embedded).  See motor lamination knowledge.
"""


CARSTENSEN = r"""
# Carstensen AC copper loss + iron loss extension

## Status: PRODUCTION (in radia_mcp.peec.carstensen_ac_loss)

Carstensen 2007 PhD (public-safe curated corpus 2007
PhD; cross-link from peec.carstensen_ac_loss):
- Originally for SRM (switched reluctance motor) winding losses
- Dowell-Kelvin functions for skin + proximity effect in slot windings
- Extended to per-layer iron loss decomposition

## Iron loss method

For laminated stator with windings, per-layer Bertotti applied to:
```
H_t(z) = staircase from window edge to interior
P_layer = k_h * f * B_layer^n + k_c * (B_layer * f)^2
```

This avoids meshing individual laminations.

## Cross-reference

See MCP: `peec_carstensen_ac_loss(topic='motor_application')`.
"""


WAVEFORM_CORRECTIONS = r"""
# Iron loss under non-sinusoidal waveforms

## Why ordinary Steinmetz fails

Steinmetz assumes B(t) = B_peak * sin(2πft).  Real motors have:
- PWM-driven flux ripple
- Field-weakening induced harmonics
- Saturation-induced harmonics

The rms dB/dt under PWM can be 5-10× the sinusoidal equivalent at the
same fundamental frequency, leading to massive iron loss
underestimation (typical 30-100%).

## Decision flowchart

```
Have measured B(t) waveform OR simulated FE B(t)?
├── YES, want fastest estimate:
│   └── MSE with k_h, n from datasheet
│
├── YES, DC bias OR asymmetric PWM:
│   └── iGSE with alpha, beta from datasheet
│
├── NO, only sinusoidal datasheet:
│   └── Steinmetz original + 10-30% safety factor for PWM
│
├── Need PER-LAMINATION accuracy:
│   └── Hollaus MSFEM effective material (motor_hollaus_eddy)
│
├── Need PER-DOMAIN accuracy:
│   └── Vector hysteresis (Energy-Based or E&S) + Bertotti per element
│
└── Need TRUE TRANSIENT:
    └── Time-domain Play / Energy-Based hysteresis solver
        (calc_motor_transient.py + Picard nonlinear FE)
```

## Lab papers (public-safe curated corpus)

- Iron Loss Evaluation of Reactor Core with Air Gaps by Magnetic
  Field (Magnetic-field-based, not assumed-sinusoid)
- Cauer等価回路とプレイモデルを用いた数値解析による
  インバータ励磁下における鉄損特性評価 (PWM-specific)
- 高周波励磁下におけるリアクトル損失の空隙特性
- 電磁鋼板の詳細な積層モデルを用いた高周波用リアクトルの鉄損解析

Cross-reference: see `minor_loops_dc` topic for the Taitoda/
Takahashi/Fujiwara 2015 method that handles DC-biased minor loops
within a non-sinusoidal waveform using bilinear loss coefficients
k_h(B_m, B_max) and an apparent frequency.
"""


MINOR_LOOPS_DC = r"""
# Iron loss with minor loops + DC bias (Taitoda 2015)

## Problem

When the flux density waveform b(t) contains higher harmonics, the
hysteresis loop develops MINOR LOOPS riding on the ascending /
descending branches of the major loop.  Conventional Steinmetz /
Bertotti formulas:
```
P = k_h(B_m) * f + k_e(B_m) * f²
```
calibrate k_h, k_e from SINUSOIDAL excitation only.  When applied
to a minor loop with DC offset, the result can be off by ~17%
(Taitoda 2015 measured) because:
1. The minor loop is on a non-zero DC bias B (DC component of b).
2. The minor loop period ≠ 1/f of the harmonic.

## Proposed method (Taitoda/Takahashi/Fujiwara 2015 IEEE TMAG)

### Step 1: separately measure DC-biased loss coefficients

Replace single-argument k_h(B_m), k_e(B_m) with **bilinear**:
```
k_h(B_m, B_max)   B_max = B_m + ΔB   (max flux = AC amplitude + DC bias)
k_e(B_m, B_max)
```
- Measure on ring core under DC-biased AC excitation:
  b = B_m sin(ωt) + ΔB
- Sweep ΔB ≠ 0 to populate the 2D table
- Linearize P/f vs f to extract slope (k_e) and intercept (k_h)
- Typical measurement frequencies: 50, 150, 750 Hz

### Step 2: decompose b(t) into outline + minor loops

For the actual waveform b_original(t):
```
b_outline(t) = b_original(t) with minor-loop segments removed
b_minor(t)   = each minor-loop trace, attached to its DC bias level
```

### Step 3: apparent frequency for each piece

Each minor loop has its own period T^m_app = t_2 - t_1 where
b(t_1) = b(t_2) AND h(t_1) = h(t_2) on the loop.  The apparent
frequency f^m_app = 1 / T^m_app is generally HIGHER than the harmonic
frequency f_h.  Similarly f^o_app = 1 / T^o_app for the outline.

### Step 4: per-piece loss with the right coefficient set

```
P_minor   = sum over each minor loop of
            [ k_h(B_m, B_max) * f^m_app + k_e(B_m, B_max) * (f^m_app)² ]
            using DC-biased coefficients
P_outline = k_h(B_m_outline) * f^o_app
          + k_e(B_m_outline) * (f^o_app)²
            using AC (DC-bias = 0) coefficients
P_total   = P_minor + P_outline
```

## Accuracy (Taitoda 2015 measurement)

For non-oriented steel 50A470 ring core with 23rd-order harmonic
B_m23 = 0.2-0.8 T superimposed on fundamental:

| Method | Max error vs measurement |
|--------|--------------------------|
| Conventional (single k_h(B_m), single f) | 17 % |
| Proposed (k_h(B_m, B_max), apparent f) | 6 % |

The improvement is largest in the high-flux region where the
DC-biased coefficients deviate most from the sinusoidal coefficients.

## When to use this method

| Scenario | Use this method |
|----------|------------------|
| PWM-driven motor with significant harmonic content | YES |
| Reactor with DC bias on AC excitation | YES |
| Symmetric sinusoidal Epstein-frame measurement only | NO (no minor loops to model) |
| Already running Play / Energy-Based time-domain solver | NO (loss comes out naturally from dissipation integral) |

## Reference

- T. Taitoda, Y. Takahashi, K. Fujiwara (Doshisha University),
  "Iron Loss Estimation Method for a General Hysteresis Loop With
  Minor Loops", IEEE Trans. Magn. 51(11):8112304, Nov 2015.
  doi:10.1109/TMAG.2015.2445930
- File: `public-safe curated corpus
   00_review/Iron_Loss_Estimation_Method_for_a_General_Hysteresis_
   Loop_With_Minor_Loops.pdf`

Cross-ref:
- `iron_loss_knowledge('waveform')` for the general PWM decision tree
- `magnetic_materials_hysteresis('lab_core')` for the time-domain
  alternative (B-input Stop+Energy automatically captures minor-loop
  dissipation)
"""


def get_iron_loss_knowledge(topic: str = "decision") -> str:
    """Dispatch by topic.

    Topics:
        steinmetz_family - Steinmetz + MSE + iGSE family
        bertotti         - Bertotti 3-term decomposition + steel grades
        carstensen       - Carstensen AC copper + iron loss (cross-ref to peec)
        waveform         - Non-sinusoidal (PWM, DC bias) corrections
        minor_loops_dc   - Taitoda 2015 DC-biased k_h(B_m,B_max) + apparent
                            frequency for general hysteresis with minor loops
                            (best accuracy for PWM-driven motors)
        decision         - Decision flowchart (default)
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("steinmetz", "steinmetz_family", "mse", "igse"):
        return STEINMETZ_FAMILY
    if topic in ("bertotti", "bertotti_3term", "three_term", "anomalous"):
        return BERTOTTI_3TERM
    if topic in ("carstensen", "carstensen_loss"):
        return CARSTENSEN
    if topic in ("waveform", "pwm", "non_sinusoidal", "harmonic"):
        return WAVEFORM_CORRECTIONS
    if topic in ("minor_loops_dc", "minor_loops", "taitoda", "dc_biased",
                  "apparent_frequency"):
        return MINOR_LOOPS_DC
    if topic in ("decision", "decision_tree", "flow"):
        return WAVEFORM_CORRECTIONS  # decision flowchart is in waveform topic
    if topic == "all":
        return "\n\n".join([STEINMETZ_FAMILY, BERTOTTI_3TERM, CARSTENSEN,
                              WAVEFORM_CORRECTIONS, MINOR_LOOPS_DC])
    return (f"Unknown topic '{topic}'. Available: steinmetz_family, "
            "bertotti, carstensen, waveform, minor_loops_dc, decision, all.")
