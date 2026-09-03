# Periodic Boundary Validation

This numerical-evidence lane verifies periodic and anti-periodic finite-element
spaces used by sector models. It checks manufactured-solution error, the
anti-periodic sign change, and high-order h convergence.

```powershell
python -m pytest validation_test/periodic_boundary/test_periodic_bc.py -q
python validation_test/periodic_boundary/generate_periodic_bc_results.py
```

The generator owns `periodic_bc_results.json`; do not edit its numerical values
by hand.
