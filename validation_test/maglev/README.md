# MagLev validation

This directory owns durable numerical evidence for magnetic-levitation and
moving-conductor workflows. The public demonstrations and their figures remain
under `docs/maglev`; checked JSON records belong here.

`demos/` mirrors the public demonstration families:

- root records cover Mixed-Galerkin polarizability, HCurl bulk modes, passive
  ROM fitting, and the moving-magnet magnetic-Reynolds crossover;
- `ellipsoid/` covers analytic and FEM shape-anisotropic polarizability;
- `sphere/` covers induced-dipole force, full eddy-force comparison, and stable
  coil equilibrium;
- `team28/` preserves the full-FEM/CLN force curve and published-height check.

Static evidence tests are fast and do not regenerate the records. Slow solver
reproduction tests remain in this validation lane and use temporary output
directories so a test run never mutates committed evidence.

## ECB plate force references

`ecb_foster_lorentz_reference.py` solves the scalar eddy model of
`radia.maglev.ecb.lorentz` directly on a mirror-symmetric structured mesh and
reconstructs the current as `J = (1/mu) curl(v z)`.  It is a **same-model
self-consistency and Foster-truncation lane**, not an independent physical
reference.  Symmetry, sign, mesh convergence, and the perfect-conductor image
bound catch severe defects but do not validate the scalar ansatz.

`ecb_foster_lorentz_3d_reference.py` is the independent numerical lane.  It
uses a divergence-free three-dimensional HCurl current basis and the
open-boundary Newton-potential VIM interaction, then integrates the resulting
volume current with the shared peak-phasor Lorentz kernel.  Rank and mesh
refinement converge below 0.2%.  On the recorded 100 x 60 x 2 mm plate, the
best lateral-boundary scalar solution still differs from the converged 3-D
lift by about 63%, 87%, and 89% at 50, 500, and 5000 Hz.  The scalar model is
therefore rejected for quantitative three-dimensional plate force.  Use
`compute_lorentz_force_via_hcurl_vim`; retain Foster only as a reduced solution
of its stated scalar model, preferably through the verified direct-solve
fallback.

Until 2026-09-11 the kernel built the current as `-omega sigma Im(v)`: zero
lift for a centred magnet and a horizontal force that mirror symmetry forbids
-- 1906 N at 5 kHz against the 11.7 N image bound.  The summary's `history`
keeps that record; `test_ecb_foster_lorentz_reference_evidence.py` replays the
summary and runs the kernel live on a small mesh.

No experimental data set was supplied for this case.  The 3-D result is an
independent numerical validation, not a measurement-validation claim.
