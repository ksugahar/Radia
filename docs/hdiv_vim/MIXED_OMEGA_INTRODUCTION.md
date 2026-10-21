# Mixed Omega: declared validation scope

The [executed planar-interface introduction](mixed_omega_planar_interface.ipynb)
shows the CPU matching-trace total/reduced scalar-potential route on an
independently specified harmonic problem. Its two materials have relative
permeabilities 100 and 1. The prescribed reduced boundary value is nonzero:
it equals the physical scalar potential minus the declared source potential.
The reconstructed reduced-region field includes that source field.

The example uses a quadratic H1 space and explicit CG with local
preconditioning. NGSolve owns the spaces, orientation, integration, assembly,
and field evaluation. The caller owns the TaskManager region. The notebook
saves the mesh and primary-field WebGUI scenes; viewing them additionally
requires `anywidget`. Notebooks are introductions, not application workbenches.

The true assembled residual and the Euclidean vector-field RMS error are
separate checks. Polynomial reproduction in this example certifies only the
specified geometry, boundary data, order and evaluation procedure.
Overall formulation parity remains **OPEN**; neither this example nor a small
residual establishes an unrestricted whole-interface equivalence theorem.

Kelvin and separated-body coupling have their own executable investigation:
`validation_test/magnetostatics/mixed_omega_kelvin_coupled_probe.py`.
It records the explicitly selected backend, reconstruction route, declared
stencil and numerical gates in JSON. Its matching lift must use the same
operator action or identical quadrature as its operator; matching only the
bonus integration order is insufficient. A fresh MMM residual with voxel
transfer certifies the solve against that voxel-approximated external field.
An explicit harmonic source representation additionally needs its independent
holdout check against direct source evaluation. These certificates have
different scopes and must not be substituted for one another.

The CPU solver remains independent of the optional resident GPU validation
adapter. GPU execution must be selected explicitly and cannot silently replace
a failed CPU route. Solver backend availability is not a performance claim.
