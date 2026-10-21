# Four-property IH implementation and CPU gate

Start source: `54612dfd70cc729db8cc20734dc902ac8d1a8ea7`.

This preregisters the CPU validation before execution. The validation is a bounded CPU reference run. It runs no GPU computation.

The contract is `radia.ih.four-property.v1`. Temperature is absolute kelvin;
electrical resistivity and constant mass density have distinct units. NGSolve
owns geometry, spaces, integration measures, basis and gradient evaluation.
The new CPU/native operators contract the provided quadrature maps. Maps must
preserve the constant field and its zero gradient. The dense reference route
is bounded to 256 coefficients; it is not a scalability claim or a CPU fallback.

## Acceptance declared before execution

- Every accepted step: input power equals enthalpy increase divided by dt plus
  outward convection/radiation loss, relative defect <=1e-6. Normalize by the
  larger of absolute input, absolute storage plus boundary loss, and 1e-12 W.
- Independently observed peak-phasor port power equals coil resistance loss
  plus integrated body heat, relative defect <=1e-6 using the same floor.
  The voltage projection must be independent of heat weights.
- Linear true residual <=1e-6; nonlinear true residual <=1e-6 and energy gate
  are separate. Normalize the nonlinear residual by the maximum Euclidean
  norm of load, assembled storage, conduction and boundary vectors, with a
  1e-12 W floor. The reference uses the stricter nonlinear tolerance 1e-10.
  Include cooling at zero and small positive load.
- Piecewise-linear cp(T) is integrated exactly for enthalpy, including partial
  intervals. Compare against an independent polynomial primitive.
- Nonuniform kappa(T) conduction and the full Newton Jacobian are checked
  against polynomial weak-form integration and central differences.
- With nonlinear emissivity(T), compare dt, dt/2, dt/4 temperatures against an
  independent scalar RK4 oracle. Expected partitioned temporal order is one;
  require successive error ratios between 1.6 and 2.4. Recheck energy at each dt
  using analytical storage and boundary loss, not the production balance report.
- Reject out-of-range/nonfinite material values and wrong port-power sign.
  Failed candidates must preserve accepted temperature and diagnostics.

`validate_four_property_cpu.py` is an independent CPU reference gate, not a
GPU, full PEEC-BEM, nonlinear BH cell or FEM validation. The native BH gateway
now checks shape, finite values, interval/saturation continuity and positive
differential permeability; native compilation and execution remain required.

## Remaining integration gates

The physical producer must supply checked NGSolve quadrature maps, full EM
operators, loss measures and an independent port-voltage projection. The old
geometry producer and legacy Thermal S-Function remain fixed-property paths;
they are not evidence for this explicit new route. Production SLX wiring and
Toolkit lifecycle acceptance remain a separate stage.

Native build has not been performed and remains unverified.

GPU numerical agreement, independent FEM cross-validation, native artifacts
and release acceptance follow implementation review and new resource slots.

The production thermal energy gate is a constant-field projection of the
assembled residual. It is a consistency guard, not independent evidence.
Independent evidence here comes from analytical enthalpy/boundary expressions
and the separate RK4 reference; port voltage must be independently assembled.
