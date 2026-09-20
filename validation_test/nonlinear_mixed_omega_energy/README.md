# Nonlinear mixed-Omega field/energy refinement

Run the three-level P2-response/P1-material study with:

```powershell
python validation_test\nonlinear_mixed_omega_energy\run_refinement.py C:\temp\nonlinear_mixed_omega_energy.hdf5
```

The script evaluates average and RMS flux density together with magnetic
energy and coenergy from the same converged solution and material-state
identity. It writes the laboratory HDF5 exchange schema and applies the public
refinement gate. The result establishes self-convergence only; it does not by
itself establish agreement with another solver.
