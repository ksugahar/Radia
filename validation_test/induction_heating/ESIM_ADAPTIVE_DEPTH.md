# Adaptive ESIM depth cells

The cell discretization is independent of the outer electromagnetic surface
mesh. A converged panel iteration or a passive impedance alone does not certify
the depth resolution. The finite slab and cylinder use their actual thickness
or radius; the half-space route additionally checks domain-depth doubling.

Both routes solve the same conservative nonuniform diffusion stencil with
row scaling and a checked algebraic residual. The peak-phasor convention is
exp(+i omega t), with `mu = mu_prime - i mu_loss` for passive magnetic loss.
The geometric seed uses a fixed physical grading strength, rather than a
fixed ratio raised to the number of nodes. Increasing the seed node count
therefore does not exponentially collapse the smallest interval.

Each adaptation compares a cell mesh with its nested midpoint refinement.
The relative changes in complex `Z`, integrated active loss, and reactive
power must all be below `mesh_tol`. Local nested Joule-integral differences
mark depth intervals containing 75% of the estimated error for refinement.
This is an empirical convergence estimator, not a rigorous global bound.
The half-space calculation separately doubles the domain until the impedance
and powers agree within the same tolerance. Four depth extensions are allowed.

`ESIMFiniteSlabSolver`, `ESIMCellProblemSolver` and
`EMMaterial.create_esim_solver` accept `mesh_tol=1e-3`, `max_nodes=32769`,
`max_refinements=16`, and `cell_max_iter=200`. `n_nodes` selects the seed,
not the accepted final resolution. `solve(max_iter=...)` explicitly overrides
the nonlinear cap. Exhausted nonlinear, mesh, or depth budgets raise; failed
cells cannot silently populate an impedance table. Panel cache identities
include all these controls. Accepted adaptive fields carry their own mesh;
the solver's declared seed is restored after success or failure.

Run `validate_esim_adaptive_depth.py --output result.json` for synthetic
linear slab/cylinder/half-space analytic comparisons, high-node stability,
and independent seed convergence for a synthetic saturating BH table.
Predeclared limits are 0.3% against analytic solutions and 0.3% across seeds.
Fast tests also check a deliberately insufficient iteration/node budget,
passive complex permeability, and stable cache identity across adaptation.

This validation does not establish accuracy for a supplied laboratory material
or a complete coil. Such cases need matched volumetric FEM, mesh and external
domain convergence, and self-consistent electromagnetic recomputation.
Frozen-field heat rescaling is sensitivity analysis only.

ESIM remains an independent phasor constitutive approximation using the
declared amplitude-dependent law. It is not evidence of a time-resolved
large-amplitude nonlinear BH periodic solution. That solver and its GPU/FEM
acceptance remain separate requirements.
