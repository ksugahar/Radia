# Nonlinear Electromagnet Verification (20000 AT)

Electromagnet simulation with nonlinear iron (B-H curve).

**Updated 2026-02-05**: Cleaned up folder structure. Full model verified against ELF reference.

## Reference

- ELF_MAGIC model: `S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1`
- B-H curve: `BH.txt` (100 data points extracted from ELF)

## Model Parameters

- **Geometry**: C-type yoke (full model, 52 hexahedral elements)
- **Material**: Nonlinear B-H curve (100 data points)
- **Coil**: 20000 AT racetrack coil

## ELF Reference

- **Bz at (0,0,0)**: -994.12 mT

## Full Model Verification Results (Recommended)

| Solver | Method | Bz (mT) | Error | Status |
|--------|--------|---------|-------|--------|
| **LU** | 0 | -995.01 | **+0.09%** | **PASS** |
| **BiCGSTAB** | 1 | -996.53 | **+0.24%** | **PASS** |
| **HACApK** | 2 | -995.04 | **+0.09%** | **PASS** |

**All three solvers agree within 0.25% of ELF reference.**

## IMA Symmetry Status

### Linear Materials: IMA Works Correctly

| Material | Model | Bz (mT) | Error vs Full |
|----------|-------|---------|---------------|
| **Linear (mu_r=1000)** | Full model | -17.72 | - |
| **Linear (mu_r=1000)** | Quarter IMA (+x-z) | -17.72 | **0.00%** |

**Linear IMA is verified working** - full model and quarter model with IMA produce identical results.

### Nonlinear Materials: KNOWN LIMITATION (~22% Error)

**WARNING**: IMA produces incorrect results for nonlinear materials.

| Material | Model | Bz (mT) | Error vs Full |
|----------|-------|---------|---------------|
| **Nonlinear** | Full (52 elem) | -995.01 | - |
| **Nonlinear** | Quarter IMA (+x-z) | **-770.02** | **~22%** |

### Root Cause

The nonlinear iteration chi update uses `H = M/chi` (constitutive relation) which doesn't account
for demagnetizing field contributions from IMA mirror elements. For linear materials this works
correctly (chi is constant), but for nonlinear materials the H field used for chi lookup is
incomplete.

**Investigation Note (2026-02-05)**: Attempted fix by adding mirror demagnetizing field to H
during chi update made results **worse** (74% error vs 22%). The issue may be deeper in the
nonlinear iteration formulation with IMA, not just the chi update step.

### Workaround

- **Linear materials**: IMA works correctly (0% error)
- **Nonlinear materials**: Use **full model (no IMA)** for accurate results

**Recommendation**: Use **full model (no IMA)** for nonlinear problems.

## Files

| Directory | Description |
|-----------|-------------|
| `full/` | Full model tests (52 elements, recommended) |
| `full/LU/` | LU solver tests |
| `full/bicgstab/` | BiCGSTAB solver tests |
| `full/hacapk/` | HACApK solver tests |
| `quarter/` | Quarter model tests with IMA (+x-z) |
| `quarter/LU/` | LU solver with quarter model |
| `quarter/bicgstab/` | BiCGSTAB solver with quarter model |
| `quarter/hacapk/` | HACApK solver with quarter model |
| `BH.txt` | B-H curve data (100 points) |
| `coil_model.py` | Full racetrack coil model |
| `coil_model_quarter.py` | Quarter coil model for IMA |

## Recommended Usage (Full Model)

```python
import radia as rad
import numpy as np

rad.FldUnits('m')

# Create full model geometry (52 elements)
yoke = rad.ObjCnt(hex_elements)
mat = rad.MatSatIsoTab(bh_data)  # Nonlinear B-H curve
rad.MatApl(yoke, mat)

# Create full racetrack coil
from coil_model import create_racetrack_coil
coil = create_racetrack_coil(20000.0)
model = rad.ObjCnt([yoke, coil])

# Solve (any solver works)
result = rad.Solve(model, 0.001, 100, 0)  # LU solver

# Get field
B = rad.Fld(model, 'b', [0, 0, 0])
print(f"Bz = {B[2]*1000:.2f} mT")  # Expected: ~-995 mT
```

## IMA Symmetry

IMA symmetry (`+x-z`) means:
- **+x**: X-symmetric (same pole)
- **-z**: Z-antisymmetric (opposite pole)

**Linear materials**: IMA works correctly (0% error)
**Nonlinear materials**: Use full model - IMA has ~22% error

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

## Run Verification

```bash
# Recommended: Full model test (all solvers)
cd full/LU && python test_all_solvers.py
```

## Changelog

### 2026-02-05 (Update 2)
- **Linear IMA verified**: 0.00% error (full model vs quarter IMA)
- **Nonlinear IMA**: Documented as known limitation (~22% error)
- Investigated fix attempt - adding mirror H-field made results worse (74% error)
- Updated README with verification results

### 2026-02-05
- Cleaned up debug/temporary scripts
- Fixed typo: quater -> quarter
- Reorganized solver folders

### 2026-02-04
- **Full model verified**: All solvers agree within 0.25% of ELF
- **IMA issue identified**: ~2% error vs full model
- **Recommendation changed**: Use full model, not IMA

### 2026-02-01
- Verified LU solver with IMA (2.87% error, now known to be IMA issue)
- Fixed coil model: Use quarter coil for IMA symmetry

### 2026-01-31
- Updated to use new Image symmetry API (`image='+x-z'`)
- Removed deprecated TrfMlt API
