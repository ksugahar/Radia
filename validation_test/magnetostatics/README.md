# Magnetostatics Validation

This directory owns numerical evidence for permanent-magnet, susceptibility,
and nonlinear magnetostatic formulations. These checks are not part of the
per-change implementation-regression lane.

The nonlinear manufactured-solution suite verifies solution error, active
reluctivity variation, and Newton's quadratic residual reduction:

```powershell
python -m pytest validation_test/magnetostatics/test_nonlinear_magnetostatic_newton.py -q
python validation_test/magnetostatics/generate_nonlinear_newton_results.py
```

The generator owns `nonlinear_newton_results.json`; do not edit its numerical
values by hand.
