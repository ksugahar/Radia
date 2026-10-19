# Real thermal BiCGStab regression

The thermal correction changes the BiCGStab shadow-vector choice. It is an
algorithm change, not merely error handling. With the solver's zero initial
iterate, the initial residual is `r0 = b`. Choosing `r0_hat = r0` is a standard
BiCGStab choice; for nonzero real `b`, its initial inner product is `b^T b > 0`.
The zero-RHS path returns zero without dividing by that inner product.
[Netlib Templates, Bi-CGSTAB](https://www.netlib.org/templates/templates.pdf)
gives this initial-residual shadow choice in its algorithm on printed page 24.

For real matrix entries and a real RHS, this choice keeps all arithmetic in
the real subspace. An artificial complex shadow needlessly permits imaginary
truncation error in the thermal increment. A residual accepted at `1e-10`
does not imply that every imaginary component meets the separate `1e-12`
state-admissibility test. Removing that source of imaginary error preserves
the state test rather than discarding imaginary parts or relaxing a limit.

This does not make BiCGStab breakdown impossible. Breakdown and failure of the
original-equation residual still raise; there is no numerical fallback.
Time-domain thermal matrices with complex entries explicitly reject.

The generic `ComplexSolver` default remains the artificial complex shadow.
Only `CoupledState`'s thermal solver selects the RHS shadow. Its eddy solver
still passes its preconditioner and uses the unchanged default. The regression
uses a genuinely complex, nonsymmetric eddy matrix and compares its full
complex solution with an independent NumPy solve.

## Executable gate

Generate the independent fixture with `real_thermal_fixture.py`. Run
`checkRealThermalRegression` on an exact hash-verified MEX. It checks the five
thermal updates against NumPy, independent original thermal/eddy equations,
the real-state constraint, and destruction. No Simulink model is edited.

Also run this same regression against the preserved pre-fix binary as an
explicit negative control. Its candidate-temperature rejection must be
recorded as an expected regression detection, never as acceptance of that
binary. Keep old and corrected source/binary identities separate. Existing
failure evidence is immutable.

The fixture is synthetic fixed coupled algebra. It does not validate physical
PEEC-BEM assembly, dynamic ESIM, or a production Simulink application.
