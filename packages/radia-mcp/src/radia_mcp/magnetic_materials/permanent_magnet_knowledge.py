"""Permanent magnet datasheet handling: NdFeB, SmCo, Ferrite, AlNiCo.

Includes temperature derating, demagnetization curves, BHmax operating point.
"""

PM_FAMILIES = r"""
# Permanent magnet families: comparison table

| Family | Br [T] | Hc [kA/m] | BHmax [kJ/m³] | T_max [°C] | μ_recoil | Cost | Notes |
|--------|--------|-----------|---------------|------------|----------|------|-------|
| NdFeB N52 | 1.4-1.45 | 950-1000 | 400-440 | 80-100 | 1.05 | $$$$ | Strongest; brittle; corrosion |
| NdFeB N42SH | 1.28-1.32 | 920-960 | 320-340 | 150 | 1.05 | $$$$ | High-temp NdFeB |
| NdFeB N40UH | 1.25-1.29 | 900-940 | 300-320 | 180 | 1.05 | $$$$$ | Very high temp |
| NdFeB N35EH | 1.17-1.21 | 870-910 | 260-280 | 200 | 1.05 | $$$$$$ | Extreme temp |
| SmCo 2:17 | 1.05-1.15 | 700-800 | 200-240 | 350 | 1.05 | $$$$$ | Best high-temp; expensive |
| SmCo 1:5 | 0.85-1.0 | 640-720 | 150-180 | 250 | 1.05 | $$$$ | Older SmCo grade |
| Ferrite (sintered Y30) | 0.36-0.40 | 220-260 | 25-30 | 250 | 1.10 | $ | Cheap; weak; T-stable |
| Ferrite (bonded) | 0.20-0.30 | 130-200 | 8-15 | 100 | 1.1-1.2 | $ | Flexible; injection mold |
| AlNiCo 5 | 1.25-1.28 | 50-60 | 40-46 | 525 | 4.5 | $$ | Linear demag; very-high-T |
| AlNiCo 8H | 0.80-0.90 | 130-170 | 38-52 | 525 | 2.0 | $$$ | More coercive |

## Temperature coefficient summary

| Family | αBr [%/K] | αHc [%/K] |
|--------|-----------|-----------|
| NdFeB | -0.10 to -0.13 | -0.50 to -0.65 |
| SmCo | -0.03 to -0.045 | -0.20 to -0.30 |
| Ferrite | -0.20 (POSITIVE Hc αHc = +0.30) | +0.30 |
| AlNiCo | -0.025 | +0.02 to -0.02 |

NdFeB DROPS QUICKLY with temperature.  SmCo MUCH more stable.

## Selection by application

| Application | Magnet | Why |
|-------------|--------|-----|
| Brushless DC servo (-20 to 100°C) | NdFeB N42 / N48 | Highest force density |
| Industrial servo / spindle (high T) | NdFeB N40SH / N35EH | T-margin |
| Aerospace, harsh env | SmCo 2:17 | T_max 350°C; corrosion-resistant |
| Loudspeaker / cheap motor | Ferrite Y30 | Cost; OK for low B |
| Magnetic separator, holding | AlNiCo 5 | Linear B-H; deep field reach |
| Mag levitation, MRI shim | NdFeB N52 + thermal mgmt | Max field per kg |

## Reference
- Coey, "Magnetism and Magnetic Materials", Cambridge 2010 (Ch 13)
- du Trémolet de Lacheisserie, "Magnetism Vol II", Springer 2005

## Lab folder
W:/磁気特性/04_永久磁石/FORC解析と永久磁石材料への適用.pdf — FORC
applied to PM characterisation (Mayergoyz Preisach formalism for
intrinsic vs normal demag curve disentangling).
"""


DEMAG_CURVE = r"""
# Demagnetization curves: normal vs intrinsic, operating point, knee

## Two curves: B-H vs J-H (M-H)

Permanent magnet datasheets show TWO curves in the 2nd quadrant:
- **B-H curve (normal demag)**: B = μ₀(H + M) — the EXTERIOR flux density
- **J-H curve (intrinsic demag)**: J = μ₀M — the MATERIAL magnetization

For NdFeB at room temperature:
- B-H: starts at Br ≈ 1.4 T, intersects H-axis at H_c ≈ 950 kA/m
- J-H: starts at Br ≈ 1.4 T, drops sharply at H_ci > H_c (intrinsic
  coercivity).  The KNEE of J-H is where irreversible demagnetization begins.

KEY POINT: stay above the knee on J-H curve to prevent irreversible
demagnetization.  Working point migration toward knee = thermal /
field overload risk.

## Operating point selection

For a PM in a magnetic circuit:
1. Plot the load line on the B-H plane: slope = -B/H = -μ₀ * permeance/area
2. Intersect with B-H demag curve → operating point (B_op, H_op)
3. Check that H_op > H_knee_on_JH (with safety margin)
4. BHmax point ≈ midpoint of demag curve (gives max field density)

## BHmax design philosophy

To use the PM efficiently:
- Operate AT or NEAR the BHmax point
- BHmax = (B_op * H_op)_max ≈ Br * Hc / 4 (theoretical maximum)
- For NdFeB N52: BHmax ≈ 0.7 T × 700 kA/m = 490 kJ/m³ → actual ~420 kJ/m³

## Irreversible demagnetization

If field swings beyond the knee at high T (e.g., short-circuit current
in motor):
- Some moments lock to wrong orientation
- New (lower) B-H curve replaces original
- Recovery requires re-magnetization

DESIGN RULE: maintain H_at_PM > 0.7 * H_knee under all transients.

## Temperature-derated demag curves

NdFeB datasheet typically shows B-H at:
- 20°C, 60°C, 100°C, 140°C, 180°C
Br DECREASES, Hc DECREASES, knee MOVES into operating region as T
rises.  This is the main reason for upgrading from N42 → N42SH for
hot operation.

## Radia implementation

Permanent magnets are MMMM (surface-charge moment) elements carrying a FIXED magnetization given in
the constructor -- **no Solve needed**: an unsolved element has `sigma == 0`, so `rad.Fld` takes the
analytic magnetization-based open-boundary branch (`B_comp_frM`).  Call `Solve` only when soft iron is
present alongside the magnet (then put both in an `ObjCnt` and Solve the container).

**Rectangular magnet -> use `rad.magnet_box` (the ObjRecMag substitute):**
```python
import radia as rad, math
M = 1.2 / (4 * math.pi * 1e-7)                                  # M = Br / mu_0, A/m (Br=1.2 T -> 954930)
pm = rad.magnet_box([0, 0, 0], [0.02, 0.02, 0.01], [0, 0, M])  # (center, dims=FULL edges, M); meters
B  = rad.Fld(pm, 'b', [0, 0, 0.03])                            # no Solve
```
`magnet_box(center, dimensions, magnetization)` builds an `ObjHexahedron` from the box corners and is
the **drop-in replacement for the retiring `ObjRecMag`** (CLAUDE.md "Reduce Proprietary API Surface":
the surface-current `ObjRecMag` primitive demotes; the surface-charge MMMM block is the user path).  A
uniformly magnetized block has the IDENTICAL external field in both models -- verified to ~5e-8..1e-7
relative vs `ObjRecMag` off-axis, far field, and on the magnetization axis (golden
`tests/test_magnet_box_pm.py`).  Same call shape + units as `ObjRecMag`, so it is a direct swap.

**General shape (non-box) -> build the element directly with M in the constructor:**
```python
# Br = 1.2 T (NdFeB N42); M = Br / mu_0 in A/m = 954930
M_x = 1.2 / (4 * math.pi * 1e-7)
pm_obj = rad.ObjHexahedron(verts, [M_x, 0, 0])   # or ObjTetrahedron / ObjWedge (verts, M)
```
All are MMMM surface-charge moment elements; all evaluate the exact open boundary with no Solve.

UNITS PITFALL (recurring): Radia magnetization is **A/m, NOT Tesla**.  Passing
`[0, 0, 1]` thinking "1 T" gives M = 1 A/m -> field ~0.0002 mT (essentially zero),
a silent wrong result.  Always convert: `M = Br / mu_0` (Br=1.0 T -> 795775 A/m;
Br=1.05 T -> 835559 A/m; Br=1.2 T -> 954930 A/m).  Sanity-check that the probe
field is in the expected mT-to-T range, not micro-Tesla.

### Worked example: SmCo hexagonal array (VERIFIED)

`docs/smco_magnet_array/smco_magnet_array.ipynb` -- 125 cylindrical Sm2Co17
magnets (Br = 1.05 T -> M = 8.36e5 A/m, `ObjCylMag` 16-sided) on a meshed
soft-iron base plate (`MatLin(1000)`, wedge core ring + hex annular sectors).
Self-contained notebook (code + `smco_results.json` + `smco_bz_map.png`).
VERIFIED 2026-06-26 against the built wheel: Solve (LU) converges in 1 iter,
|B| = 177.6 / 68.4 / 202.9 mT at (0,0,20) / (0,0,50) / (30,0,20) mm.  Note the
meshed-disk core ring is a 6-vertex **wedge** -> `ObjWedge` (NOT `ObjHexahedron`);
both element constructors require the magnetization argument.

For demag analysis (full PM demag not yet implemented; the MatMagCurve skeleton was removed 2026-06-26):
- Currently approximate as MatLin(mu_recoil) plus initial M
- Full demag would track local H trajectory and detect knee crossing

## Lab paper: FORC for PM characterisation

`W:/磁気特性/04_永久磁石/FORC解析と永久磁石材料への適用.pdf`:
FORC (First-Order Reversal Curves) measurement disentangles:
- Reversible (recoil-line) component → μ_recoil
- Irreversible (full demag) component → Preisach distribution

For accurate PM modeling under thermal cycling, FORC-based Preisach is
the gold standard.  Not yet in Radia (TODO).
"""


def get_permanent_magnet_knowledge(topic: str = "families") -> str:
    """Dispatch by topic.

    Topics:
        families      - PM family comparison + selection guide
        demag_curve   - B-H vs J-H, BHmax, knee, operating point
        all           - Both
    """
    topic = topic.lower().strip()
    if topic in ("families", "comparison", "selection"):
        return PM_FAMILIES
    if topic in ("demag_curve", "demag", "bhmax", "knee", "operating_point"):
        return DEMAG_CURVE
    if topic == "all":
        return PM_FAMILIES + "\n\n" + DEMAG_CURVE
    return f"Unknown topic '{topic}'. Available: families, demag_curve, all."
