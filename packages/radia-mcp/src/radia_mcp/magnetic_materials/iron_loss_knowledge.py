"""Iron loss models: Steinmetz, Bertotti 3-term, Pry-Bean, MSE, iGSE,
Carstensen, anomalous loss.  W:/磁気特性/鉄損/ + adjacent.
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

## Lab papers (W:/磁気特性/鉄損/)
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

## Parameter extraction

For commercial NO (non-oriented) and GO (grain-oriented) steel:

| Grade | k_h | n | k_c | k_a |
|-------|-----|---|-----|-----|
| 35JN230 (NO, 0.35 mm) | 0.012 | 1.7 | 4.6e-5 | 0.0001 |
| 35JN300 (NO, 0.35 mm) | 0.018 | 1.7 | 4.6e-5 | 0.00015 |
| 50JN800 (NO, 0.50 mm) | 0.035 | 1.7 | 1.3e-4 | 0.00025 |
| 27ZH95 (GO, 0.27 mm)  | 0.005 | 1.6 | 2.8e-5 | 0.00005 |

(approximate; check manufacturer datasheet for exact values; lab data
in `W:/磁気特性/00_教科書/電磁鋼板の磁気特性と取扱方法.pdf`)

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
- Lab textbook: `W:/磁気特性/00_教科書/General properties of power
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

Carstensen 2007 PhD (W:/MOR_モデル縮約/Eddy_Currents/Carstensen 2007
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

## Lab papers (W:/磁気特性/鉄損/)

- Iron Loss Evaluation of Reactor Core with Air Gaps by Magnetic
  Field (Magnetic-field-based, not assumed-sinusoid)
- Cauer等価回路とプレイモデルを用いた数値解析による
  インバータ励磁下における鉄損特性評価 (PWM-specific)
- 高周波励磁下におけるリアクトル損失の空隙特性
- 電磁鋼板の詳細な積層モデルを用いた高周波用リアクトルの鉄損解析
"""


def get_iron_loss_knowledge(topic: str = "decision") -> str:
    """Dispatch by topic.

    Topics:
        steinmetz_family - Steinmetz + MSE + iGSE family
        bertotti         - Bertotti 3-term decomposition + steel grades
        carstensen       - Carstensen AC copper + iron loss (cross-ref to peec)
        waveform         - Non-sinusoidal (PWM, DC bias) corrections
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
    if topic in ("decision", "decision_tree", "flow"):
        return WAVEFORM_CORRECTIONS  # decision flowchart is in waveform topic
    if topic == "all":
        return "\n\n".join([STEINMETZ_FAMILY, BERTOTTI_3TERM, CARSTENSEN,
                              WAVEFORM_CORRECTIONS])
    return (f"Unknown topic '{topic}'. Available: steinmetz_family, "
            "bertotti, carstensen, waveform, decision, all.")
