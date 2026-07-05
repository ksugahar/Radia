"""WPT efficiency analysis (Q, k, kQ) + safety + IEC/SAE standards."""

Q_FACTOR_KQ = r"""
# Q factor, k coupling, kQ product

## Definitions

**Q (quality factor)** per coil:
```
Q_i = ω L_i / R_i    (component Q, R_i = AC resistance)
```

**k (coupling coefficient)**:
```
k = M / √(L_1 L_2),    M = mutual inductance,    0 ≤ k ≤ 1
```

**kQ product** (the key WPT figure of merit):
```
kQ = k √(Q_1 Q_2)
```

## Maximum efficiency

For series resonant SS topology with optimal load matching:
```
η_max = (kQ)² / (1 + √(1 + (kQ)²))²
       ≈ 1 - 2/(kQ)        for kQ >> 1
```

| kQ value | η_max |
|----------|-------|
| 1 | 38% |
| 3 | 75% |
| 10 | 91% |
| 30 | 97% |
| 100 | 99% |

**Target for production WPT**: kQ ≥ 10 (>90% efficiency).

## How to maximize kQ

**Maximize Q**:
- Litz wire (reduces AC resistance via skin/proximity)
- Ferrite plates (reduces eddy loss in chassis)
- Avoid lossy materials near coil
- High-Q capacitors (NP0 ceramic, polypropylene film)

**Maximize k**:
- Close coil spacing
- Aligned axes
- Ferrite flux concentrator
- Larger coil area (within mechanical constraint)

**Trade-off**: Larger coil increases L AND k, but also increases R.
Q is usually proportional to √(coil_diameter). So scaling up helps
both, but ferrite + Litz wire are bigger levers.

## kQ Poincaré distance (Sakai et al.)

Recent unified design metric:
```
kQ-distance = log((1+kQ)/(1-kQ))
```
Maps kQ ∈ [-1,1] to ℝ with constant differential metric, useful for
optimization landscapes.

public-safe curated corpus
- 近傍ダイポールアンテナ間のkQ積とポアンカレ距離 (WPT2018-18)

## Cross-references

- `radia_mcp.peec` — k, L computation via PEEC filaments
- `radia_mcp.peec.carstensen_ac_copper_loss` — AC R for Q
- `radia_mcp.ih` — analogous loss in workpiece (= bad-Q metal load)
"""


COUPLING_K = r"""
# Coupling coefficient k — measurement + computation

## Definition

```
k = M / √(L_1 L_2)
```

where M = mutual inductance, L_i = self inductance.

## Measurement methods

1. **Direct LCR meter** (parasitic-free): measure L_open, L_short on
   primary, then `k = √(1 - L_short/L_open)`
2. **VNA S-parameter**: extract from S11/S21 at multiple frequencies
3. **Time-domain step response**: M = dV/dt fit
4. **Power-balance method** (under load): infer from input/output ratio

## Computation methods

| Method | Tool | Accuracy | When |
|--------|------|----------|------|
| Analytical (Neumann) | hand calc | exact for circular coaxial coils | quick sizing |
| Filament (Grover) | scipy / hand | exact for thin-wire | medium |
| PEEC | `radia.peec_matrices` | engineering | lab production ★ |
| FEM (volume mesh) | NGSolve | benchmark | accuracy check |
| BEM (surface) | `ngsolve.bem` | open boundary | research |

Lab production = PEEC: `from radia.peec_matrices import PyPEECBuilder`

## k as function of geometry

For two coaxial circular planar coils (radius a, separation z):
- k ∝ (a/z)³ for z >> a (far field, dipole-dipole)
- k → ~0.5 max for z ~ a (close coupling)
- Higher with ferrite: k_ferrite ≈ 1.5-2x k_air

## Misalignment sensitivity

dk/d(lateral) is large near nominal alignment. DD/DDQ coil designs are
chosen for **flat k(lateral)** out to ±10 cm at z = 15 cm.

| Coil pair | k(0) | k at 10cm lateral | drop |
|-----------|------|--------------------|------|
| Circular  | 0.30 | 0.10 | -67% |
| Square    | 0.28 | 0.15 | -46% |
| DD        | 0.22 | 0.16 | -27% |
| DDQ       | 0.25 | 0.20 | -20% ★ |

## References

public-safe curated corpus (15 papers on coupling)
public-safe curated corpus (7 papers, dynamic k)
"""


SAFETY_HUMAN_EXPOSURE = r"""
# Safety: human exposure to WPT magnetic fields

WPT at kW level radiates substantial magnetic field. Human exposure
limits from:

- **ICNIRP 2010** (general public, 100 kHz):  H_rms ≤ 21 A/m
- **SAE J2954** (EV-specific):  more restrictive, depends on coil design
- **IEC 62311** (general device safety)
- **Japan: 電波防護指針** (radio protection guideline)

## Compliance methods

1. **Design-side**:
   - Ferrite shielding (block stray flux outside primary)
   - Symmetric coil geometry (cancel far field)
   - DDQ / DD (lower stray field than circular)
   - Reduced power on misalignment

2. **Detection-side**:
   - Living-object detection (LOD) — see `pcb_fod`
   - Auto-shutdown on human proximity

## Special concerns

- **Implanted pacemakers / ICDs**: high field can trigger malfunction
  - HF-band WPT (>1 MHz) particularly concerning
  - IEC 62366-1 + medical device EMC standards
  - public-safe curated corpus 2018 papers (WPT2018-14, 18-53):
    植込み心臓ペースメーカEMI 評価

- **Metallic implants (knee, hip)**: induction heating risk
  - Test in saline phantom + thermocouple measurement
  - Power derating if metallic implant detected

- **Fetal exposure**: special limits in some jurisdictions

## Cross-references

- IEC 61980-1 general EV WPT safety
- public-safe curated corpus
- public-safe curated corpus (Japanese 規制)
"""


STANDARDS = r"""
# WPT standards: IEC, SAE, ISO

## IEC 61980 series (general WPT for EV)

| Part | Scope | Status |
|------|-------|--------|
| **IEC 61980-1** | General requirements | Ed 2.0 (2020) |
| **IEC 61980-2** | Communication / control | Ed 1.0 (2019) |
| **IEC 61980-3** | Specific requirements for MF | FDIS 2024 |

Frequencies: 81.39-90 kHz (EV stationary), with provision for dynamic.

## SAE J2954 (EV stationary WPT, North America)

Released 2018, updated 2020. Power classes:
- WPT1: 3.7 kW (single-phase)
- WPT2: 7.7 kW
- WPT3: 11 kW (three-phase capable)
- WPT4: 22 kW

Operating freq: **85 kHz nominal** (81-90 kHz tuning range).
Coil geometry: standardized to ensure interoperability between
chargers and vehicles from different manufacturers.

## ISO 15118 (V2G communication)

Vehicle-to-grid + V2X communication including WPT alignment messaging.

## QI (low-power)

Wireless Power Consortium standard for phones/wearables:
- Qi 1.x: 5W, 110-205 kHz
- Qi 2.0 (2023): 15W, magnetic-aligned (MagSafe-style)

## FOD compliance

Each standard mandates FOD before/during charging:
- IEC 61980-3 FDIS: specific test items 項目1, 項目2, 項目3
  (public-safe curated corpus)

## Lab compliance status

- 85 kHz development: lab is SAE J2954-aligned
- FOD: multi-method research (5+ approaches in 10_異物検出/)
- Safety: IEC 61980-1 + 電波防護指針 verification

## References

public-safe curated corpus (5 papers, including standards-aligned design)
"""


def get_efficiency_safety_knowledge(topic: str = "q_factor_kQ") -> str:
    """Dispatch efficiency / Q / k / safety / standards topics.

    Topics:
        q_factor_kQ              - Q, k, kQ product, η_max formula (DEFAULT)
        coupling_k               - k computation + measurement + misalignment
        safety_human_exposure    - ICNIRP, pacemaker, implant safety
        standards_IEC_SAE        - IEC 61980, SAE J2954, ISO 15118, Qi
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("q_factor_kQ", "q_factor", "kQ", "q", "kq"):
        return Q_FACTOR_KQ
    if topic in ("coupling_k", "coupling", "k"):
        return COUPLING_K
    if topic in ("safety_human_exposure", "safety", "human", "exposure"):
        return SAFETY_HUMAN_EXPOSURE
    if topic in ("standards_IEC_SAE", "standards", "iec", "sae", "qi"):
        return STANDARDS
    if topic == "all":
        return "\n\n".join([Q_FACTOR_KQ, COUPLING_K, SAFETY_HUMAN_EXPOSURE,
                            STANDARDS])
    return (f"Unknown topic '{topic}'. Available: q_factor_kQ, coupling_k, "
            "safety_human_exposure, standards_IEC_SAE, all.")
