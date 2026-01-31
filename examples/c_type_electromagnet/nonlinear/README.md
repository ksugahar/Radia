# Nonlinear B-H Curve Simulation (20000 AT)

Electromagnet simulation with nonlinear iron (B-H curve) using Image symmetry.

**Updated 2026-01-31**: Uses new Image symmetry API (`image='+x-z'`)

## Reference

- ELF_MAGIC model: `S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1`
- B-H curve: `S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/BHカーブ/iron.bh`

## Model Parameters

- **Geometry**: C-type yoke (quarter model, 13 hexahedral elements)
- **Material**: Nonlinear B-H curve (100 data points)
- **Coil**: 10000 AT (with Image symmetry -> 20000 AT effective)
- **Symmetry**: Image `'+x-z'` (MIMA X symmetric, MIMA -Z antisymmetric)

## Image Symmetry API

```python
import radia as rad

rad.FldUnits('m')

# Create quarter model geometry
yoke = rad.ObjCnt(hex_elements)
mat = rad.MatSatIsoTab(bh_data)  # Nonlinear B-H curve
rad.MatApl(yoke, mat)

# Create coil
coil = rad.ObjRaceTrk(center, radii, lengths, h, nseg, 'man', 'z', j)
model = rad.ObjCnt([yoke, coil])

# Solve with Image symmetry (quarter -> full model)
result = rad.Solve(model, 0.001, 500, 1, image='+x-z')
```

## B-H Curve (Selected Points)

| H (A/m) | B (T) |
|---------|-------|
| 0 | 0 |
| 82 | 1.14 |
| 898 | 1.59 |
| 4582 | 1.81 |
| 17736 | 2.01 |
| 68322 | 2.20 |
| 318000 | 2.56 |

## Comparison with ELF_MAGIC

**Test Date**: 2026-01-31

| Metric | Value |
|--------|-------|
| Radia solver | BiCGSTAB with Image symmetry |
| ELF solver | MMM with MIMA X, MIMA -Z |
| Average |B| difference | ~34% |
| Best match (Elem 9) | 0.7% |

**Notes**:
- Field magnitude differences are due to coil modeling:
  - ELF: Discretized coil (MCL8T elements, 10 segments)
  - Radia: Analytical racetrack coil (ObjRaceTrk)
- Nonlinear convergence is challenging with strong saturation

## Files

| File | Description |
|------|-------------|
| `c_type_electromagnet_nonlinear.py` | Main simulation script |
| `BH.txt` | B-H curve data (100 points) |
| `yoke.vol` | Netgen mesh file |
| `README.md` | This file |

## Run

```bash
python c_type_electromagnet_nonlinear.py
```

## Changelog

### 2026-01-31
- Updated to use new Image symmetry API (`image='+x-z'`)
- Removed deprecated TrfMlt API
- Switched to BiCGSTAB solver for better nonlinear convergence
