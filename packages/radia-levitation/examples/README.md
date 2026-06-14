# radia-levitation examples

Curated, standalone-runnable demonstrations (README-backed, per the lab
Sample Promotion Ladder).  The broad research corpus is in
`../research/`; the IGTE paper in `../papers/`.

| Folder / script | Topic | Geometry |
|-----------------|-------|----------|
| `cube_alpha_sweep.py` | Mixed-Galerkin `alpha(s)` sweep (CAD-direct edges vs mesh-derived), the headline package API demo | Cu cube |
| `sphere/` | Isotropic levitation force; analytic dipole + eddy FEM + Maxwell-stress cross-check; coil equilibrium | sphere |
| `ellipsoid/` | Shape-anisotropic polarizability tensor `alpha(omega)` (DC + HF anchors + FEM) | triaxial ellipsoid / spheroid |
| `cuboid/` | CLN-SIBC `alpha(s)` core for a brick (modal Foster + CLN + Schur-F) | rectangular cuboid |
| `team28/` | TEAM Problem 28 electrodynamic levitation: full-FEM baseline + CLN force-vs-height + equilibrium | Al disk over coils |

Each subfolder has its own `README.md` with a per-script table (purpose /
run command / headline result).  Start with `cube_alpha_sweep.py` for the
core API, then `sphere/` for the validated force workflow.
