# NGSolve-Radia Integration Examples

C++ CoefficientFunction integration for coupling Radia magnetostatics with NGSolve FEM.

## Quick Start

```bash
# 1. Build modules
cd S:\radia\01_GitHub
.\Build.ps1              # Build radia.pyd
.\build_ngsolve.ps1      # Build rad_ngsolve.pyd

# 2. Run visualization with comparison
cd examples\ngsolve_integration
python visualize_field.py --method both --maxh 0.005
```

## Key Features

### Two Evaluation Methods

**CoefficientFunction (CF)**: Exact Radia evaluation
- Always accurate (0% error)
- Best for: Few points, exact values

**GridFunction (GF)**: FEM interpolation  
- Convergent with mesh refinement
- Best for: Many points, fast evaluation

### Accuracy vs Mesh Size

| Mesh (maxh) | GF Error | Recommendation |
|-------------|----------|----------------|
| 0.020 m | ~60% | Not recommended |
| 0.005 m | ~1% | Production use |
| 0.002 m | <0.01% | High precision |

## Main Script

### visualize_field.py

Comparison tool for CF vs GF evaluation methods.

```bash
# Compare both methods (default)
python visualize_field.py --method both --maxh 0.005

# Exact evaluation only
python visualize_field.py --method cf

# Fast interpolation only
python visualize_field.py --method gf --maxh 0.005
```

Outputs VTK files for Paraview visualization.

## Documentation

- **COMPARISON_GUIDE.md** - Detailed CF vs GF comparison
- **UNIT_CONVERSION.md** - Automatic m/mm conversion
- **COMPLETE_FIX_SUMMARY.md** - Implementation details

## References

- Radia: https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia
- NGSolve: https://ngsolve.org/
