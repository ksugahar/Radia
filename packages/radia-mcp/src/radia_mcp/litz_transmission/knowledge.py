"""Litz wire + transmission line knowledge."""

OVERVIEW = r"""
# Litz wire + transmission line analysis

For kHz-MHz applications (WPT, IH, induction motors), copper conductors
have severe AC resistance increase due to skin effect + proximity.

## Skin depth

```
δ = sqrt(2 / (ω μ σ))
```

| Frequency | Cu skin depth (mm) |
|-----------|---------------------|
| 50 Hz | 9.3 |
| 1 kHz | 2.1 |
| 10 kHz | 0.66 |
| 85 kHz (EV WPT) | 0.23 |
| 1 MHz | 0.066 |
| 13.56 MHz | 0.018 |

## Solutions for AC loss

| Approach | When | Pro / Con |
|----------|------|-----------|
| **Litz wire** | < 10 MHz | many fine strands, Dowell formula |
| **Hollow conductor** | high power, > 100 kHz | water-cooled tube |
| **Foil winding** | sheet conductor, planar | 1D AC distribution |
| **Plated wire (magnetic)** | NMI | μ-plating reduces flux concentration |
"""


LITZ = r"""
# Litz wire analysis

## Dowell formula (1966)

For multi-layer foil/Litz:
```
R_AC / R_DC = ξ * (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
            + (2/3) (m² - 1) ξ * (sinh(ξ) - sin(ξ)) / (cosh(ξ) + cos(ξ))
```
where ξ = h/δ (layer height / skin depth), m = number of layers.

Implemented in `radia.rad_peec_surface_impedance.cpp` for rectangular
conductors.

## Strand count design

For round Litz wire with N strands of diameter d:
- DC R: same as single-strand of equivalent cross-section
- AC R / DC R: depends on N, d/δ, and twist pitch
- Rule: each strand should be d < δ to avoid skin effect

Lab papers (W:/.../02_リッツ線/04_strand_count_design/):
- リッス線の素線数が交流抵抗に与える影響

## Homogenization (FEM-based)

For COMPUTATIONAL efficiency, homogenize the Litz wire bundle into an
effective material:
- Effective σ_eff (complex, function of ω)
- Effective μ_eff (for magnetic-plated wire)

References:
- Eddy Current Analysis of Litz Wire Using Homogenization-Based FEM
- Fast 3-D Analysis of Eddy Current in Litz Wire Using Integral Equation
- Finite element analysis of multi-turn magnetoplated coils using homogenization
- Voxel-Based FEM Using Homogenization
- 均質化法を用いた磁性メッキ巻線の有限要素解析
- 不整合メッシュを用いた有限要素解析のための均質化法
- 異方性媒質を含む不整合有限要素メッシュ均質化法
- 非適合要素および均質化を用いた電磁界解析

W:/.../02_リッツ線/01_均質化_homogenization/ (8 papers)

## Magnetic-plated wire (lab specialty)

Some lab research focuses on copper wire with thin magnetic plating
(μ-plating) to reduce flux concentration at wire surface and lower
loss. References:
- ワイヤレス給電コイル線材のアルミ化効果
- マクロ透磁率を含む系の等価回路合成
- 渦電流を考慮したマクロ磁気特性
- 渦電流を考慮したマクロ磁気特性細線の解析

## High-frequency resistance models

Beyond Dowell:
- 銅クラッドアルミ線における高周波抵抗解析式の導出
- 複素透磁率を用いた積分方程式

## Cross-references

- `radia_mcp.peec.carstensen_ac_copper_loss` — Carstensen AC copper loss with FE-extracted H
- `radia_mcp.motor.hollaus_eddy` — Hollaus MSFEM for laminated cores
- `radia_mcp.fem.large_scale_special.hollaus_msfem` — MSFEM theory
- `radia_mcp.wpt.coil_compensation` — Litz wire in WPT coils
"""


TRANSMISSION_LINE = r"""
# Multiconductor transmission lines

For PCB / cable / busbar at MHz-GHz, multiconductor transmission line
(MTL) theory governs propagation, coupling, crosstalk.

## Telegrapher's equations

```
d V/d z = -(R + j ω L) I
d I/d z = -(G + j ω C) V
```

For N conductors: V, I are N-vectors; R, L, G, C are N×N matrices.

## Parameter extraction

The L, C, R, G matrices can be computed from cross-section geometry:
- **2D quasi-static FEM** (electrostatic for C, magnetostatic for L)
- **PEEC** (lab method) — extends to 3D
- **Analytical** for simple geometries (microstrip, coaxial)

## Lab approach

- `radia.peec_matrices` for L, M extraction
- `radia.ngsolve` for 2D MTL cross-section
- Cross-link to `radia_mcp.peec`

## Reference

[LOCAL] W:/.../39_部品材料/03_伝送線路/ (22 files)
- The Electromagnetic Theory of Coaxial Transmission Lines and Cylindrical Shields
- Analysis of Multiconductor Transmission Lines
- Continuum Representation of Wound Coils via an Equivalent Foil Approach
- Time-domain homogenization of multiturn windings based on RL Cauer ladder networks

## Cross-references

- `radia_mcp.peec` — PEEC for L, C extraction
- `radia_mcp.mor` — CLN / Cauer ladder for transmission lines
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch Litz / transmission line topics.

    Topics:
        overview         - Skin depth + AC loss summary (DEFAULT)
        litz             - Litz wire (Dowell, homogenization, magnetic-plated)
        transmission     - Multiconductor transmission lines (MTL)
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("litz", "litz_wire", "wire"):
        return LITZ
    if topic in ("transmission", "transmission_line", "mtl"):
        return TRANSMISSION_LINE
    if topic == "all":
        return "\n\n".join([OVERVIEW, LITZ, TRANSMISSION_LINE])
    return f"Unknown topic '{topic}'. Available: overview, litz, transmission, all."
