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

## ECB plate force reference

`ecb_foster_lorentz_reference.py` solves the scalar eddy model of
`radia.maglev.ecb.lorentz` directly on a mirror-symmetric structured mesh and
reconstructs the current as `J = (1/mu) curl(v z)` -- `v` is the reaction-field
z component in tesla, so this is the only reconstruction with the units and
parity of a real eddy current.  It checks what a centred z-dipole over a plate
must show: zero horizontal force, a repulsive lift, mesh convergence, a lift
rising with frequency, and a lift below the infinite perfect-conductor image
bound.  The shipped `compute_lorentz_force_via_foster` is held to the same
symmetry and sign and to the reference lift within the Foster truncation
error.

Until 2026-09-11 the kernel built the current as `-omega sigma Im(v)`: zero
lift for a centred magnet and a horizontal force that mirror symmetry forbids
-- 1906 N at 5 kHz against the 11.7 N image bound.  The summary's `history`
keeps that record; `test_ecb_foster_lorentz_reference_evidence.py` replays the
summary and runs the kernel live on a small mesh.
