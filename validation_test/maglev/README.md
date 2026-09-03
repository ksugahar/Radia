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
