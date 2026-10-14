# Shared energy and coenergy diagnostics

The current-excited ESRF coil-yoke driver accepts `--energy-quadrature-order N`.
It samples all three formulations on the **same** physical FEM quadrature,
with positive volume weights and the same single-valued isotropic PCHIP B(H)
law and vacuum-slope continuation. Nonlinear energy uses the inverse of that
same law, not a second H(B) interpolant. Linear vacuum contributions use the
usual quadratic density. Iron and air contributions are reported separately.

Stored energy is evaluated from B, coenergy from H. Their sum minus the
measured integral of H dot B is recorded as a constitutive consistency
diagnostic (Fenchel gap), rather than forced to zero algebraically. For HDiv,
B inside iron includes magnetization: B = mu0 (H_source + H_demag + M).
The gap-field shortcut B = mu0 H is invalid inside iron.

The coil-yoke driver includes the Kelvin exterior: it maps quadrature points
to physical exterior points and multiplies volume weights by (R/rho)^6.
Mixed Omega uses the documented H and B form pullbacks. Reduced-A transforms
only curl(A_reaction), then adds the source evaluated at the physical point.
HDiv evaluates its physical field directly. A dipole's analytic exterior
energy tests the volume Jacobian and the Omega field pullbacks.

Calling the observer without an explicit Kelvin map instead returns a clearly
labeled finite-volume partial integral. A Kelvin-only mesh excludes the physical
interior and is also partial. Integration-order and mesh/domain convergence must
still be measured on each real model before these observables replace the
existing accuracy gate. An iron-only integral or a sum over gap observation
points is not an open-boundary energy.
Permanent magnets, hysteresis, anisotropic laws and image-reduced iron meshes
are not supported by this observer and must not be relabeled as this contract.

The comparison reports the relative scalar spread and whether HDiv lies
between reduced-A and mixed Omega. This ordering is a measured diagnostic,
not a guaranteed upper/lower bound. No new numerical tolerance is inferred
from the old field-RMS tolerance. Energy post-processing time is reported
separately from solve/field-observation runtime.

Use `--engines mixed_total_reduced_omega` to finish a selected solve without
waiting for reduced-A. Subset output is explicitly marked as partial, not a
successful three-engine comparison. New energy runs have a distinct checkpoint
contract: old B-observation checkpoints cannot supply energy integrals.

For the direct reduced-A comparison, `--reduced-a-solver direct` and
`--reduced-a-solver sparsecholesky` both select NGSolve SparseCholesky.
PARDISO is not an automatic alternative. Mixed Omega is a remaining migration
task: its multiplier system produced
non-finite values in a small SparseCholesky trial. A guard now rejects such
values immediately; its existing backend defaults are not yet migrated.
Small-system and nonlinear regression checks do not establish large-system readiness: an ESRF case 3 matrix assembled
successfully but its SparseCholesky factorization terminated with an access
violation under NGSolve 6.2.2606 on Windows, both inside and outside TaskManager.
No residual or nonlinear acceptance was obtained for that case.

Repository policy also covers direct subsolves. Existing PARDISO use elsewhere,
including the native HDiv mass factor and older application runners, remains
migration debt until the replacement's residual, matrix-structure support and
representative-size checks pass. This change does not claim those paths have
already been migrated.
