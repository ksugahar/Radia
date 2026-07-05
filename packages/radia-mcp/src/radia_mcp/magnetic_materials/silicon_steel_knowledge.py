"""JIS silicon steel grade database + handling notes.

Source: public-safe curated corpus (94 MB,
the canonical JFE/Nippon Steel handling handbook).
"""

JIS_GRADES = r"""
# JIS 電磁鋼板 grade naming convention

## JIS standard format: AABBYYC
- AA: nominal thickness in 0.01 mm (e.g., 35 = 0.35 mm, 50 = 0.50 mm)
- BB: grade letters
  - JN = Japanese Non-oriented (NO)
  - JZ = Japanese Grain-oriented Zero (GO, conventional)
  - JL = Japanese Laser-scribed (high-permeability GO)
  - JF = Japanese High-magnetic-flux (HiB)
- YY: nominal max iron loss at 1.5T, 50Hz in 0.1 W/kg
- C: optional surface treatment code

Examples:
- **35JN230** = 0.35 mm NO, max 2.30 W/kg at 1.5T 50Hz
- **35JN300** = 0.35 mm NO, max 3.00 W/kg at 1.5T 50Hz (cheaper, lossier)
- **50JN400** = 0.50 mm NO, max 4.00 W/kg
- **50JN800** = 0.50 mm NO, max 8.00 W/kg (cheapest NO)
- **23ZH90** = 0.23 mm GO conventional, max 0.90 W/kg
- **27ZH95** = 0.27 mm GO, max 0.95 W/kg
- **35JF180** = 0.35 mm HiB, max 1.80 W/kg

## Typical motor / transformer applications

| Application | Recommended grade | Reason |
|-------------|------------------|--------|
| Mass-market induction motor (400 V) | 50JN800, 50JN1300 | Cheap, low-thickness OK |
| Premium / high-efficiency motor | 35JN230, 35JN270 | Lower iron loss |
| Servo / spindle (high-freq PWM) | 35JN230 + amorphous | Multi-material |
| Distribution transformer (50 Hz) | 27ZH95, 35ZH95 (GO) | GO direction-fixed |
| High-frequency reactor (>1 kHz) | 10JNHF / amorphous / ferrite | Thin lam or amorphous |
| Switching power supply (>50 kHz) | Ferrite (3F3, N87), nanocrystalline | Iron unsuitable |

## Material constants (approximate, sinusoidal at 1.5 T, 50 Hz)

| Grade | k_h | n | k_c | k_a | mu_r_max | B_sat |
|-------|-----|---|-----|-----|----------|-------|
| 35JN230 | 0.012 | 1.7 | 4.6e-5 | 1e-4 | ~10000 | 2.05 T |
| 35JN300 | 0.018 | 1.7 | 4.6e-5 | 1.5e-4 | ~9000 | 2.05 T |
| 50JN800 | 0.035 | 1.7 | 1.3e-4 | 2.5e-4 | ~7000 | 2.05 T |
| 50JN1300 | 0.060 | 1.7 | 1.3e-4 | 4e-4 | ~5000 | 2.05 T |
| 23ZH90 (GO) | 0.005 | 1.6 | 2.0e-5 | 4e-5 | ~50000* | 2.03 T |
| 27ZH95 (GO) | 0.005 | 1.6 | 2.8e-5 | 5e-5 | ~45000* | 2.03 T |

(* GO permeability is HIGHLY anisotropic; ~50000 along rolling direction,
~3000 perpendicular.  Use `rad.MatLin([mu_r_par, mu_r_perp], axis)`.)

NOTE: For PRODUCTION analysis, ALWAYS use manufacturer datasheet
constants rather than these approximate values.

## Datasheet sources

- **JFE Steel** (旧 NKK / 川崎製鉄): https://www.jfe-steel.co.jp
- **Nippon Steel** (旧 新日鉄): https://www.nipponsteel.com
- Both publish detailed PDF datasheets with B(H), P(B, f) curves

The lab's reference text `電磁鋼板の磁気特性と取扱方法.pdf` (94 MB)
covers JFE-specific grades and contains the most complete tabulation
for Japan.
"""


HANDLING_NOTES = r"""
# 電磁鋼板 取扱方法 (handling and processing) notes

Source: `public-safe curated corpus`
(JFE / 川崎製鉄 handbook).

## Critical processing effects on iron loss

| Process step | Iron loss change | Notes |
|--------------|-----------------|-------|
| Slitting / shearing | +5 to +30% | Burrs cause local short circuits |
| Stamping | +10 to +40% | Mechanical stress at die edges |
| Compressive stacking | +5 to +15% | Lamination pressure (axial) |
| Welding (lamination ends) | +3 to +10% | Local short circuit |
| Annealing (post-stamp, 800°C N₂) | -20 to -40% | Stress relief |
| Coating (insulation) | +2 to +5% | Surface oxide layer eats Mn |
| Bake-hardening (paint) | +5 to +10% | Reverses some annealing |

## Practical implications for FE analysis

If you use the DATASHEET iron loss in FE, REAL motor iron loss is
typically **30-50% higher** due to:
- Edge effect from punching (untreated)
- Mechanical clamping in motor stack
- Imperfect insulation (eddy current paths between laminations)

## Building factor (BF)

Industry standard correction:
```
P_iron_real = BF * P_iron_datasheet
BF ≈ 1.4 - 1.7  (typical 1.5 for stamped + glued NO)
BF ≈ 1.2 - 1.3  (for laser-cut GO)
```
Use BF in calc_motor_transient.py iron-loss post-process if datasheet
constants are used.

## Anisotropy notes

NO (Non-Oriented) steel: μ_r isotropic within ±5% in rolling vs
transverse direction.  OK to use `MatLin(mu_r)` scalar.

GO (Grain-Oriented) steel: μ_r anisotropic by factor 10-20.  MUST use
`MatLin([mu_r_par, mu_r_perp], axis=rolling_axis)`.  For motor:
GO usually only in transformer or 1D-flux applications.

## Lamination thickness vs frequency

Skin depth in iron at 50 Hz, μ_r=5000, σ=4e6:
δ = sqrt(2/(ωμσ)) ≈ 0.5 mm

So 0.35 mm and 0.50 mm laminations are BELOW skin depth at 50 Hz — full
flux penetration.  At 500 Hz (10× freq), δ ≈ 0.16 mm: 0.35 mm lam is
2× skin depth → significant eddy effect.

This is why high-frequency motors use thin laminations (0.20 mm or
amorphous foil 25 μm).

## Hollaus MSFEM connection

For 3D analysis of laminated stack:
- Coarse mesh + Hollaus effective material μ_eff(ω) captures
  per-lamination eddy effect without meshing individual sheets
- See MCP: `motor_hollaus_eddy('effective_material')`
- See Radia: `calc_motor_lamination.py` (production)
"""


def get_silicon_steel_knowledge(topic: str = "grades") -> str:
    """Dispatch by topic.

    Topics:
        grades         - JIS grade naming + material constants
        handling       - Processing effects on iron loss + Building Factor
        all            - Both
    """
    topic = topic.lower().strip()
    if topic in ("grades", "jis", "database"):
        return JIS_GRADES
    if topic in ("handling", "processing", "building_factor", "anisotropy"):
        return HANDLING_NOTES
    if topic == "all":
        return JIS_GRADES + "\n\n" + HANDLING_NOTES
    return (f"Unknown topic '{topic}'. Available: grades, handling, all.")
