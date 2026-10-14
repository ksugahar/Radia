# Closed-coil RT0 current with SparseCholesky subsolves

The mixed RT0/P0 operator is unchanged. Its pressure Schur complement is applied using the full velocity mass inverse. CG uses a diagonal-mass Schur preconditioner. Both SPD factors use NGSolve SparseCholesky; no pressure penalty or physical mass lumping is introduced. The original free-row mixed residual and normalized current divergence are checked before returning.

Run with a source-matched Radia environment:

```powershell
python validation_test/closed_coil_current/run_schur_validation.py --output private-runtime-path
```

This explicitly uses a synthetic Netgen OCC loop and is a conservation/convergence validation, not a timing comparison. The stored 2607 result covers 30,013 and 30,190 total DOFs at metre and millimetre scales: 25/26 Schur iterations, mixed residual <= 1.82e-14, relative divergence <= 2.14e-12, section current 2 A to rounding precision. Unit tests separately compare all free velocity and pressure unknowns with an independent dense solve of the original saddle system. These results do not establish ESRF or full HDiv-MMM acceptance.
