# Parallel Wires 2D - Adaptive Mesh Refinement with Kelvin Transformation

平行2線問題（2次元静磁場）における Kelvin変換 + 適応メッシュ細分化の実装。

## Problem Setup

- **Problem**: Two parallel wires with opposite currents (+I, -I)
- **Formulation**: 2D A-formulation (scalar potential Az)
- **Infinite domain**: Handled via Kelvin transformation with periodic BC
- **Wire distance**: 1.4 m (wires at x = ±0.7 m)
- **Wire radius**: 0.02 m
- **Kelvin boundary radius**: a = 1.0 m

## Folder Structure

Per-order runner scripts under `order=*/.../*.py` have been promoted to
`docs/kelvin/kelvin_adaptive_mesh_archive.ipynb`; full source text and SHA-256
hashes live in `docs/kelvin/kelvin_adaptive_mesh_archive_results.json`.

```
平行2線/
├── order=2/
│   ├── Refine_with_zz_estimator/   # Doerfler marking (adaptive)
│   ├── Refine_all_elements/         # Uniform refinement
│   └── metric_based/                # Metric-based remeshing
├── order=3/
│   ├── Refine_with_zz_estimator/
│   ├── Refine_all_elements/
│   └── metric_based/
├── order=4/
│   ├── Refine_with_zz_estimator/
│   ├── Refine_all_elements/
│   └── metric_based/
└── readme.md
```

## Refinement Strategies

### 1. Refine_with_zz_estimator (Doerfler Marking)
- Uses ZZ-type error estimator to identify elements with large error
- Doerfler marking strategy: marks elements until cumulative error exceeds threshold
- Refines only marked elements (h-refinement)
- Stop condition: 8 iterations

### 2. Refine_all_elements (Uniform Refinement)
- Refines all elements uniformly at each iteration
- Provides reference convergence rate
- Stop condition: 8 iterations

### 3. metric_based (Metric-Based Remeshing)
- Computes ideal mesh size for each element based on error
- Formula: h_ideal = h_current * (eta_target / eta_element)^(1/(p+1))
- Regenerates entire mesh with local size field
- Stop condition: 8 iterations

## ZZ Error Estimator

The Zienkiewicz-Zhu type error estimator uses flux recovery:

1. Compute B = curl(A) from FEM solution
2. Compute H = nu * B
3. **Recover H* using HCurl interpolation with `recovery_order = solution_order - 1`**
4. Error = ||H* - H||

### Recovery Order Selection

For higher-order elements, using `recovery_order = solution_order - 1` provides better error estimation:

| Solution Order | Recovery Order | HCurl Order |
|----------------|----------------|-------------|
| p = 2          | 1              | HCurl(order=1) |
| p = 3          | 2              | HCurl(order=2) |
| p = 4          | 3              | HCurl(order=3) |

This choice avoids over-fitting in the recovery process and provides more reliable error indicators for adaptive refinement.

## Output

Each script generates:
- **PNG files**: `*_iter_XX.png` for each iteration
  - Top-left: Inner domain ZZ error map + flux lines
  - Top-right: Outer domain ZZ error map + flux lines
  - Bottom-left: DOF vs Error Estimator convergence
  - Bottom-right: DOF vs Magnetic Energy
  - suptitle: "Solution order=X, Iteration Y: ..."
- **MAT files**: `*_iter_XX.mat` containing convergence history

## Expected Convergence Rates

For smooth solutions with polynomial order p:
- Error ~ O(N^(-p/2)) where N = DOFs
- Adaptive refinement achieves optimal rate with fewer DOFs than uniform

## Usage

The historical per-order runner source is now read from
`docs/kelvin/kelvin_adaptive_mesh_archive.ipynb` and its synchronized JSON.

## Dependencies

- NGSolve
- Netgen (OCC geometry)
- NumPy
- Matplotlib
- SciPy (for .mat file export)
