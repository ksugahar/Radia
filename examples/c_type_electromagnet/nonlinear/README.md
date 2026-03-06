# Nonlinear Electromagnet Verification (20000 AT)

Electromagnet simulation with nonlinear iron (B-H curve).

**Updated 2026-02-07**: IMA verified correct for both linear and nonlinear materials.

## Reference

- B-H curve: `BH.txt` (100 data points)

## Model Parameters

- **Geometry**: C-type yoke (quarter: 13 elements, full: 52 hexahedral elements)
- **Material**: Nonlinear B-H curve (100 data points)
- **Coil**: 20000 AT racetrack coil (FULL coil for all models)

## ELF Reference

- **Bz at (0,0,0)**: -994.12 mT

## Full Model Verification Results

| Solver | Method | Bz (mT) | Error vs ELF | Status |
|--------|--------|---------|-------|--------|
| **LU** | 0 | -995.01 | **+0.09%** | **PASS** |
| **BiCGSTAB** | 1 | -996.53 | **+0.24%** | **PASS** |
| **HACApK** | 2 | -995.04 | **+0.09%** | **PASS** |

**All three solvers agree within 0.25% of ELF reference.**

## IMA Symmetry Verification

### IMA Works Correctly for Both Linear and Nonlinear

| Material | Model | Bz (mT) | Error vs Full |
|----------|-------|---------|---------------|
| **Nonlinear** | Full (52 elem), no IMA | -995.01 | ref |
| **Nonlinear** | Quarter (13 elem) + IMA +x-z | -995.01 | **0.00%** |
| **Nonlinear** | Mirrored quarter (52 elem), no IMA | -995.01 | **0.00%** |
| **Linear (mu_r=1000)** | Full (52 elem) | -231.68 | ref |
| **Linear (mu_r=1000)** | Quarter + IMA +x-z | -231.68 | **0.00%** |

**IMA is verified correct** - full model and quarter model with IMA produce identical results
for both linear and nonlinear materials.

### Previous "22% Error" Explanation

The earlier claim of ~22% IMA error for nonlinear materials was caused by using a **quarter coil**
(`coil_model_quarter.py`) instead of a **full coil** with IMA. The correct setup is:

- **Full coil** (`coil_model.py` with 20000 AT) + IMA `+x-z` for quarter model
- The full coil generates the correct external field; IMA handles the iron symmetry

Using a quarter coil with IMA double-counts the symmetry reduction, leading to incorrect excitation.

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

# Radia always uses meters

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

**Both linear and nonlinear materials**: IMA works correctly (0% error).
Always use **full coil** (not quarter coil) with IMA.

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

### 2026-02-07
- **IMA verified correct for nonlinear**: 0.00% error (full model vs quarter+IMA)
- Previous "22% error" was caused by using quarter coil instead of full coil
- Tested: 1x1x1 (52 elem) and V304 (296 elem) meshes both confirm IMA=Full

### 2026-02-05
- Cleaned up debug/temporary scripts
- Full model verified: All solvers agree within 0.25% of ELF

### 2026-01-31
- Updated to use new Image symmetry API (`image='+x-z'`)
