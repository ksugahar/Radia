"""Alternative WPT: capacitive, microwave/rectenna, metamaterial."""

CAPACITIVE = r"""
# Capacitive WPT (electric field coupling)

Use electric field (between conductor plates) instead of magnetic field
(between coils) to transfer power.

## Principle

```
P_max = (1/2) ε A V²
        ─────────────
              d
```
where A = plate area, d = gap, V = applied voltage. Higher V → higher
power but more breakdown risk.

## Comparison with inductive (magnetic)

| Aspect | Inductive | Capacitive |
|--------|-----------|-----------|
| Frequency | 100 kHz - 6.78 MHz | 1-10 MHz |
| Distance | mm - cm | mm |
| Misalignment tolerance | medium (with DD) | low |
| Magnetic field stray | high | minimal |
| Voltage at coupling | low (current high) | high |
| Best for | EV, phone | tablet, glass-mounted IoT |

## When capacitive WPT wins

- ✓ EMI-sensitive environment (medical, military)
- ✓ Light receivers (no ferrite needed)
- ✓ Glass / transparent applications
- ✓ Phone Qi 2 + glass display

## When inductive wins

- ✓ High power (kW)
- ✓ Large air gap (cm)
- ✓ Misalignment robustness

## Lab status

Reference / niche; lab focus is inductive resonant.

## References

W:/.../05_ワイヤレス給電/03_IEICE/04_dynamic_charging/
- 電界結合型WPTの漏洩電磁界を低減する整流回路配置 (WPT2018-7)
W:/.../05_ワイヤレス給電/99_その他/99_misc/06_論文_電磁界結合/
- 磁界_電界ハイブリッド結合WPTシステム
"""


MICROWAVE_RECTENNA = r"""
# Microwave WPT — rectenna and beam steering

## Principle

Radiate microwave beam from transmit antenna; receive with **rectenna**
(rectifier + antenna). Beam shapes:
- Pencil beam (focused for fixed receiver)
- Steered beam (for moving receiver)
- Sectorized (multiple receivers)

## Frequency bands

| Band | Use |
|------|-----|
| 2.45 GHz (ISM) | low-end consumer, drone |
| 5.8 GHz (ISM) | drone, longer distance |
| 28 GHz (mmWave) | high-gain antenna, drone WPT |
| 24-77 GHz auto radar | very directional, demo only |

## Rectenna design

Critical components:
1. **Antenna** (patch, helix, parabolic) — efficient capture
2. **Schottky diode** — high freq rectification (GaN, GaAs)
3. **Matching circuit** — to maximize RF→DC efficiency
4. **DC filter** — smooth output

Best practical rectenna efficiency: ~70-80% at low power (mW),
~30-50% at higher power (W).

## Drone applications

Lab papers (03_IEICE/08_microwave_WPT/):
- 28GHz レクテナを用いた飛行中ドローン定点給電 (WPT2018-2)
- WPT による Diode 装荷パッチアンテナ
- Rectenna application of quadrature hybrid (WPT2018-11)
- 2 枚の無給電素子を装荷した高利得アンテナのレクテナ
- ドローンへの無線給電のための推定方向誤差を考慮したアンテナ指向性ビーム
- 多関節アームロボットへの電力・信号同時伝送 (microwave variant)

## Lab status

Reference; lab does NOT do microwave WPT (too far from radia core scope).
Cross-reference: ngsolve.bem for full-wave high-freq simulation.

## References

W:/.../05_ワイヤレス給電/03_IEICE/03_2017_技報/ (6 papers includes microwave)
W:/.../05_ワイヤレス給電/03_IEICE/04_2018_技報/ (29 papers, multiple microwave)
"""


METAMATERIAL = r"""
# Metamaterial / negative-index for WPT

Sub-wavelength engineered materials with unusual EM properties
(negative μ, negative ε, near-zero index). Some WPT applications:

## Use cases

1. **High-Q resonator**: SRR (split-ring resonator) array boosts coil Q
2. **Field focusing**: superlens / hyperlens shapes flux to receiver
3. **Cloaking**: reduce stray field outside system
4. **Absorbing surface**: capture wasted RF and convert to DC

## Lab status

Reference only. Cross-link to `radia_mcp.metamaterial` (separate MCP)
for homogenization theory + effective medium approximations.

## References

W:/.../05_ワイヤレス給電/03_IEICE/09_metamaterial_WPT/ (none in this folder
yet; reference is in 39_部品材料/01_メタマテリアル/)
W:/.../99_アプリケーション/05_ワイヤレス給電/A Scalable, Dual-Polarized
Absorber Surface for Electromagnetic Energy Harvesting and Wireless Power
Transfer.pdf

## Cross-reference

- `radia_mcp.metamaterial` — dedicated MCP subpackage
- `radia_mcp.peec` — analyze SRR via PEEC
"""


def get_alternatives_knowledge(topic: str = "capacitive") -> str:
    """Dispatch alternative WPT topics.

    Topics:
        capacitive               - Capacitive (E-field) WPT (DEFAULT)
        microwave_rectenna       - Microwave / rectenna / drone WPT
        metamaterial             - Metamaterial-enhanced WPT
        all                      - Everything
    """
    topic = topic.lower().strip()
    if topic in ("capacitive", "electric_field", "e_field"):
        return CAPACITIVE
    if topic in ("microwave_rectenna", "microwave", "rectenna", "drone"):
        return MICROWAVE_RECTENNA
    if topic in ("metamaterial", "meta"):
        return METAMATERIAL
    if topic == "all":
        return "\n\n".join([CAPACITIVE, MICROWAVE_RECTENNA, METAMATERIAL])
    return (f"Unknown topic '{topic}'. Available: capacitive, "
            "microwave_rectenna, metamaterial, all.")
