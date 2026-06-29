# radia.maglev examples

Curated, standalone-runnable demonstrations (README-backed, per the lab
Sample Promotion Ladder).  The broad research corpus is in
`research_cln/`; the IGTE paper in `papers/`.

| Folder / script | Topic | Geometry |
|-----------------|-------|----------|
| `cube_alpha_sweep.py` | Mixed-Galerkin `alpha(s)` sweep (CAD-direct edges vs mesh-derived), the headline package API demo | Cu cube |
| `cube_alpha_sweep_figure.py` | Paper figure: `alpha(s)/V` frequency response (Re exclusion + \|Im\| loss), IEEE single-column; emits `cube_alpha_sweep_results.json` + `.pdf`/`.png` | Cu cube |
| `cube_alpha_tensor.py` | Multi-port (matrix) Mixed-Galerkin admittance `Y_ij(s)`: monopole + 3 dipole ports `{1, x, y, z}` -> matrix-CLN residue `G_n = sigma V b_n b_n^T`, surface-moment `K_mat`, edge-moment `C1_mat`, and a 4-port MIMO LTI export (`build_state_space_mimo`). The monopole port reproduces the validated scalar `alpha(s)`; the dipole block is the (cube-isotropic) multipole structure. NOT the physical vector tensor (that is `ellipsoid/`); emits `cube_alpha_tensor_results.json` | Cu cube |
| `cuboid_vector_bulk.py` | Vector (HCurl) eddy-current Foster bulk -- the de-Rham partner of the scalar H1 bulk (HDiv=demag, HCurl=eddy). Curl-curl GEP `S w = lam M w` (nograds + tree-cotree gauge, shift-invert) -> three distinct leading eddy `tau` (shape split); reproduces the analytic interior-PEC TE-mode `tau` to <2.4% at h=a/28. Interior-PEC model (not the exterior-matched physical tensor; that is `ellipsoid/`); emits `cuboid_vector_bulk_results.json` | Cu cuboid 5x2x1 |
| `physical_tensor_rom.py` | The EXTERIOR-MATCHED **physical** polarizability tensor as a passive, stable LTI. AAA discovers the dominant real poles (the physical Stoll decay times) + a log-spaced filler + NNLS passive residues (`mixed_galerkin.passive_foster_fit`). Default: analytic sphere Stoll `alpha(s)` -> ROM (pure numpy, ~1 s): AAA poles match `tau_n = mu0 sigma a^2/(n pi)^2` to ~0.00%, LTI matches `alpha(s)` to <0.2% over 1 Hz..1 GHz; emits `physical_tensor_rom_sphere.json`. `--fem`: triaxial ellipsoid 5x3x1.5 mm -> per-axis ROMs -> diagonal MIMO LTI (13 states); samples the verified 3D HCurl tensor (slow, ~tens of min), fits the CONJUGATE (causal/passive convention), band fit ~3.8% (data-limited), dominant decay times shape-split (tau_z=60.2 > tau_x=tau_y=34.6 us); emits `physical_tensor_rom_fem.json`. Sidesteps the Kameari+Kelvin accumulation breakdown by fitting the verified per-frequency solve | Cu sphere / triaxial ellipsoid |
| `rotating_magnet_eddy.py` | **Moving-magnet eddy current + the magnetic-Reynolds crossover** (Yano's rotating-PM-over-plate problem). Three ways to get J_rms / Joule / Lorentz: (1) kinematic source-only `J = -sigma dA_s/dt` (NO eddy FEM), (2) full-FEM A-phi reference (HCurl-nograds A_r + H1 phi, current continuity), (3) constant-basis multiport CLN (Tanimoto-Sugahara-Takahashi-Matsuo TWP28, generalized to a translating+rotating source by SVD-seeded block-Krylov over the transient matrix `A_sys^-1 M`). Headline: the source-only error grows with `Rm = mu0 sigma omega L^2` -- **< 0.5% below Rm~0.1 (Yano's actual case Rm~0.016, error 0.035% -> no per-step FEM needed)**, %-level above Rm~1 where the reaction matters and the CLN reproduces the full-FEM to ~1e-4..1e-6 at ~1000x less per-step cost. Validated J_rms ~ a few 100 A/m^2, \|F\| ~ a few nN (matches Yano). Emits `rotating_magnet_eddy_results.json` + `.png`. See `radia_mcp.maglev` topic `radia_iem_fem` (low-Rm kinematic shortcut) and `cln_mor_control` (high-Rm CLN) | Cu plate + rotating PM |
| `sphere/` | Isotropic levitation force; analytic dipole + eddy FEM + Maxwell-stress cross-check; coil equilibrium | sphere |
| `ellipsoid/` | Shape-anisotropic polarizability tensor `alpha(omega)` (DC + HF anchors + FEM) | triaxial ellipsoid / spheroid |
| `cuboid/` | CLN-SIBC `alpha(s)` core for a brick (modal Foster + CLN + Schur-F) | rectangular cuboid |
| `team28/` | TEAM Problem 28 electrodynamic levitation: full-FEM baseline + CLN force-vs-height + equilibrium | Al disk over coils |

Each subfolder has its own `README.md` with a per-script table (purpose /
run command / headline result).  Start with `cube_alpha_sweep.py` for the
core API, then `sphere/` for the validated force workflow.
