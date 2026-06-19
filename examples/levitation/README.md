# radia.levitation examples

Curated, standalone-runnable demonstrations (README-backed, per the lab
Sample Promotion Ladder).  The broad research corpus is in
`research_cln/`; the IGTE paper in `papers/`.

| Folder / script | Topic | Geometry |
|-----------------|-------|----------|
| `cube_alpha_sweep.py` | Mixed-Galerkin `alpha(s)` sweep (CAD-direct edges vs mesh-derived), the headline package API demo | Cu cube |
| `cube_alpha_sweep_figure.py` | Paper figure: `alpha(s)/V` frequency response (Re exclusion + \|Im\| loss), IEEE single-column; emits `cube_alpha_sweep_results.json` + `.pdf`/`.png` | Cu cube |
| `cube_alpha_tensor.py` | Multi-port (matrix) Mixed-Galerkin admittance `Y_ij(s)`: monopole + 3 dipole ports `{1, x, y, z}` -> matrix-CLN residue `G_n = sigma V b_n b_n^T`, surface-moment `K_mat`, edge-moment `C1_mat`, and a 4-port MIMO LTI export (`build_state_space_mimo`). The monopole port reproduces the validated scalar `alpha(s)`; the dipole block is the (cube-isotropic) multipole structure. NOT the physical vector tensor (that is `ellipsoid/`); emits `cube_alpha_tensor_results.json` | Cu cube |
| `cuboid_vector_bulk.py` | Vector (HCurl) eddy-current Foster bulk -- the de-Rham partner of the scalar H1 bulk (HDiv=demag, HCurl=eddy). Curl-curl GEP `S w = lam M w` (nograds + tree-cotree gauge, shift-invert) -> three distinct leading eddy `tau` (shape split); reproduces the analytic interior-PEC TE-mode `tau` to <2.4% at h=a/28. Interior-PEC model (not the exterior-matched physical tensor; that is `ellipsoid/`); emits `cuboid_vector_bulk_results.json` | Cu cuboid 5x2x1 |
| `sphere/` | Isotropic levitation force; analytic dipole + eddy FEM + Maxwell-stress cross-check; coil equilibrium | sphere |
| `ellipsoid/` | Shape-anisotropic polarizability tensor `alpha(omega)` (DC + HF anchors + FEM) | triaxial ellipsoid / spheroid |
| `cuboid/` | CLN-SIBC `alpha(s)` core for a brick (modal Foster + CLN + Schur-F) | rectangular cuboid |
| `team28/` | TEAM Problem 28 electrodynamic levitation: full-FEM baseline + CLN force-vs-height + equilibrium | Al disk over coils |

Each subfolder has its own `README.md` with a per-script table (purpose /
run command / headline result).  Start with `cube_alpha_sweep.py` for the
core API, then `sphere/` for the validated force workflow.
