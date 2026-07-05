"""Coil design + compensation topology + resonance matching."""

COIL_DESIGN = r"""
# WPT coil design

## Coil geometry families

| Type | Where | Pros | Cons |
|------|-------|------|------|
| **Planar spiral (circular)** | Phone (Qi), EV pad | Simple, well-modeled | Limited Q |
| **Planar spiral (rectangular)** | EV pad with iron | Better space utilization | Edge effects |
| **Helmholtz pair** | Lab research | Uniform H field | Large size |
| **DD (Double-D)** | EV (newer) | Good lateral tolerance | Complex |
| **DDQ (DD + Quadrature)** | EV (best) | Robust to misalignment | Multi-coil control |
| **Bipolar / TRP** | EV (alt to DD) | Single-coil simpler | Lower coupling |
| **Solenoid** | Cylindrical loads | Strong axial B | Big aspect ratio |

## Design parameters

```
N      = number of turns
d_w    = wire / Litz diameter
r_out  = outer radius
r_in   = inner radius (clear hole)
n_w    = strand count (Litz wire only)
material: copper, optionally with magnetic shield (ferrite plate)
```

## Litz wire (lab core)

★ Litz wire is critical at kHz-MHz for low AC loss. See related MCP:

- `radia_mcp.peec.carstensen_ac_copper_loss` — Carstensen AC copper loss
  with FE-extracted per-layer H
- public-safe curated corpus — 44 papers on homogenization,
  high-freq resistance, magnetic-plated wire, etc.

Key formula (Dowell):

    R_AC / R_DC = ξ * (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
                + (2/3)(m^2 - 1) * ξ * (sinh(ξ) - sin(ξ)) / (cosh(ξ) + cos(ξ))

where ξ = h / δ (layer height / skin depth), m = layer count.

## Ferrite plate (shielding + flux concentrator)

Adds μ_r ≈ 2000-3000 layer behind primary/secondary coils:
- Concentrates flux toward receiver (boosts k)
- Shields from chassis (eddy current loss reduction)
- Cost: ferrite thickness 1-3 mm typical

Lab choice (PEEC+FEM coupled):
- `radia_mcp.peec` — coil filaments
- `radia_mcp.bem.mmm_msc` — ferrite via MMM (Radia MMM)
- Coupled solve for k between primary + secondary with ferrite

## References

public-safe curated corpus (15 papers)
- 受電/送電/送受電 coil design
- Helmholtz, ヘリカル, コイル設計
"""


COMPENSATION = r"""
# Compensation topology — SS / LCC / LCL / CLLLC

For resonant inductive WPT, both primary and secondary coils need
compensation capacitors to resonate at the operating frequency
(typically 85 kHz for EV).

## Topology families

| Topology | Primary | Secondary | Pros | Cons |
|----------|---------|-----------|------|------|
| **SS** | Series C | Series C | Simple, classical | I_p depends on load |
| **SP** | Series C | Parallel C | Constant-V output | I_p depends on k |
| **PS** | Parallel C | Series C | Constant-I source | I_p high at low k |
| **PP** | Parallel C | Parallel C | Symmetric | Sensitive to detuning |
| **LCC** ★ | LCC network | Series C | Constant-I source, low Q sensitivity | Extra inductor |
| **LCL** | LCL network | LCL network | Filter harmonics | More components |
| **CLLLC** | LLC + extra C | LLC + extra C | Wide ZVS range | Complex |
| **LCC-LCC** | LCC | LCC | Constant-V or I selectable | Most components |

## When to use each (lab recommendation)

| Use case | Topology |
|----------|----------|
| Phone Qi @ < 15 W | SS or SP (simple) |
| EV stationary @ 1-22 kW (SAE J2954) | LCC-LCC ★ |
| EV dynamic (running) | LCC-LCC with feed-forward control |
| Bearingless motor + WPT | LCC-S (constant current on rotor) |
| 6.78 MHz IoT | High-Q SS with metamaterial |

## Why LCC dominates EV

LCC compensation creates a **constant-current source** behavior on the
primary side that is **invariant to k variations** (e.g. as EV moves
into/out of charging position). The extra inductor (1-10 μH) provides
this filtering.

Compensation values (LCC, lab-validated):

```
L_f = parasitic + extra (filter inductor)
C_f = 1 / (ω² × L_f)                         (resonance with L_f)
C_p = 1 / (ω² × (L_1 - L_f))                 (resonance with primary L_1)
```

## References

public-safe curated corpus (6 papers)
- ディスクリピーター for arm robot (lab WPT2016-23, 26)
- LCC-LCC topology design
- ZVS / ZCS techniques

## Cross-reference

- `radia_mcp.peec` — primary/secondary coil L via PEEC
- `radia_mcp.mor` — circuit-level reduction for WPT control
"""


RESONANCE_MATCHING = r"""
# Resonance and impedance matching

## Resonant condition

At resonance: `Z_load = Z_source*` (conjugate match).

For SS topology:
- Primary: L_1 series with C_1 → resonant freq f_r = 1/(2π√(L_1 C_1))
- Secondary: L_2 series with C_2 → also f_r
- Operating freq = f_r (match)

## Impedance reflection

From secondary back to primary:
```
Z_reflected = (ω M)² / Z_secondary
            = (ω k √(L_1 L_2))² / (R_2 + j ω L_2 + 1/(j ω C_2) + R_load)
```

At secondary resonance, the j ω L_2 + 1/(jωC_2) terms cancel:
```
Z_reflected = (ω M)² / (R_2 + R_load) = ω² k² L_1 L_2 / (R_2 + R_load)
```

This is the **reflected impedance** the primary inverter sees.

## Matching network design

Goals:
1. Maximize power transfer at fixed f
2. Maintain efficiency across k variation (load movement)
3. Soft switching (ZVS) for inverter MOSFETs

Approaches:
- **Fixed compensation** (SS, simplest)
- **Variable frequency tracking** (PLL on inverter)
- **Variable Cap** (relay-switched or motorized)
- **Variable inductance** (saturable reactor)
- **LCC** (constant-current source, k-insensitive) ★

## Quality factor Q

```
Q = ω L / R    (component Q)
Q_loaded = ω L / (R_coil + R_reflected)    (system Q)
```

Higher Q → better efficiency but narrower bandwidth / harder to control.

Lab target: Q > 100 at 85 kHz for EV WPT (Litz wire + ferrite + careful design).

## References

public-safe curated corpus (3 papers)
"""


def get_coil_compensation_knowledge(topic: str = "coil_design") -> str:
    """Dispatch coil + compensation + resonance topics.

    Topics:
        coil_design              - Coil geometry, Litz wire, ferrite (DEFAULT)
        compensation_topology    - SS / LCC / LCL / CLLLC topology choice
        resonance_matching       - Resonant condition, Z reflection, Q factor
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("coil_design", "coil", "geometry"):
        return COIL_DESIGN
    if topic in ("compensation_topology", "compensation", "topology"):
        return COMPENSATION
    if topic in ("resonance_matching", "resonance", "matching"):
        return RESONANCE_MATCHING
    if topic == "all":
        return "\n\n".join([COIL_DESIGN, COMPENSATION, RESONANCE_MATCHING])
    return (f"Unknown topic '{topic}'. Available: coil_design, "
            "compensation_topology, resonance_matching, all.")
