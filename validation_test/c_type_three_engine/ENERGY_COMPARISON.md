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

**These are finite physical-volume integrals, not all-space energies.** The
Kelvin exterior is excluded explicitly. Exterior pullback, exterior-domain
coverage and integration-order convergence must be validated before these
observables can replace the existing accuracy gate. In particular, an iron-only
integral or a sum over gap observation points is not an open-boundary energy.
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

For the direct reduced-A comparison, `--reduced-a-solver sparsecholesky`
selects NGSolve SparseCholesky explicitly. It does not silently substitute
PARDISO or AMS; existing direct/default behavior is preserved. A small-system
test covers both linear and Picard field equivalence and the residual contract.
Large-system memory and factorization behavior still require measurement.
