# Commercial-class multiphysics problems, reproduced in NGSolve

Complex problems of the kind a commercial multiphysics package is normally
reached for, solved with the radia-ngsolve stack and validated against closed
forms. Each one is also baked back into the mcp-server (reusable helpers +
queryable knowledge), so the server gets smarter, not just a pile of scripts.

| # | Problem | Status | Where | Validation |
|---|---------|--------|-------|------------|
| 1 | **Induction heating** (EM eddy -> Joule -> thermal, 1-way coupled) | ✅ done | `induction_heating.py` ; helpers `radia_ngsolve.multiphysics` (`joule_loss_density`, `solve_heat_steady`) ; `ngsolve_usage("multiphysics")` | P/L 8.9 %, T_centre 9.6 % vs exact Bessel + 1-D heat |
| 2 | **Rotating-machine torque** (PM rotor, torque-angle) | ✅ done | `motor_torque.py` ; helper `radia_ngsolve.force.eggshell_torque` (3D) ; `ngsolve_usage("multiphysics")` | tau_z(theta) vs m B0 sin(theta), 1-3 % |
| 3 | **Open boundary (Kelvin transform)** | ✅ already in lab | `radia.kelvin_{geometry,material,solver,validate,identify_ngsolve}`, `radia.kelvin_source` ; `kelvin_transformation` MCP tool ; `docs/kelvin/` | proven, high-accuracy (Nagamine CEFC 2026); NOT reinvented here |
| 4 | **TEAM-13 3-D nonlinear** (C-yoke, measured ref) | ⬜ next | (planned) | measured flux-density references |
| 5 | **Magnetic shielding** (mu-metal shell — AC/DC blog) | ✅ done | `magnetic_shielding.py` ; helpers `radia_ngsolve.solve.solve_scattered_uniform_field` + `shell_shielding_factor` ; `ngsolve_usage("multiphysics")` | shielding factor S(μr) at μr=50/100/200 vs exact shell formula, **0.9 %** |
| 6 | **3D capacitance** (electrostatics — AC/DC blog "Computing Capacitance") | ✅ done | `capacitance_3d.py` ; module `radia_ngsolve.electrostatic3d` (`solve_electrostatic_3d`, `capacitance_from_energy`) ; `ngsolve_usage("capacitance")` | spherical capacitor C=4πε₀ab/(b−a) over a gap sweep, **<0.2 %** |
| 7 | **Dielectric capacitance** (layered dielectric — AC/DC "Computing Capacitance") | ✅ done | `capacitance_dielectric.py` ; `electrostatic3d` piecewise-ε (`mesh.MaterialCF`) ; `ngsolve_usage("capacitance")` | series-dielectric spherical capacitor vs exact, **0.1 %** |
| 8 | **Halbach cylinder** (PM dipole — accelerator/undulator magnet) | ✅ done | `halbach_cylinder.py` ; `solve_planar_magnetostatic(magnet_cf=…)` + `halbach_dipole_magnetization` + `halbach_bore_field` ; `ngsolve_usage("multiphysics")` | uniform bore field B=Br·ln(Ro/Ri), **0.3 %** + uniformity |
| 9 | **Capacitance matrix** (multi-conductor — AC/DC "Computing Capacitance") | ✅ done | `capacitance_matrix.py` ; `electrostatic3d.capacitance_matrix` (FEM reaction charge) ; `ngsolve_usage("capacitance")` | closed spherical cap → exact `[[C₀,−C₀],[−C₀,C₀]]`, **0.03 %** + 0 % asymmetry/row-sum |
| 10 | **Electrostatic force** (MEMS actuator — AC/DC electrostatics) | ✅ done | `parallel_plate_force.py` ; `force.electrostatic_eggshell_force` (E-field Maxwell stress) ; `ngsolve_usage("capacitance")` | parallel-plate |F|=ε₀AV²/(2d²) (Neumann sides, exact), **0.4 %** + attractive |
| 11 | **Cylinder magnetic shielding** (2D transverse — twin of #5) | ✅ done | `magnetic_shielding_cylinder.py` ; `solve_planar_magnetostatic(magnet_cf=…)` + `shell_shielding_factor(…,"cylinder")` | S(μr) vs exact `[(μr+1)²−(a/b)²(μr−1)²]/(4μr)`, **1.3 %** |
| 12 | **Two-wire transmission line** (capacitance/length — AC/DC) | ✅ done | `two_wire_line.py` ; `scalar_fem2d.solve_electrostatic` + `capacitance` (2D) | C/length vs exact `πε₀/cosh⁻¹(D/2a)`, **0.1 %** |
| 13 | **Multipole / field-quality** (accelerator-magnet post-processing) | ✅ done | `multipole_field_quality.py` ; module `radia_ngsolve.fieldquality` (`multipole_coefficients`, `line_current_multipoles`, `superpose_multipoles`, `field_errors`) ; `ngsolve_usage("field_quality")` | normal-quad bₙ vs exact line-current superposition: main **b₂ 0.02 %**, allowed 12-pole (n=6) = 256 units = (R_ref/r0)⁴, forbidden < 0.3 units (1e-4) ; **+ independent 3-way reference (b₂ 0.0005 %)** |
| 14 | **Cylinder permanent magnet** (on-axis field — PM design) | ✅ done | `cylinder_magnet.py` ; `solve_axi_magnetostatic(magnets=…)` + `cylinder_magnet_axial_field` ; `ngsolve_usage("permanent_magnet")` | on-axis B_z vs `(Br/2)[…]`: **centre 0.01 %**, profile few-% to z=3L (pole-exit z≈L ~8 %, corner singularity) |
| 15 | **Finite solenoid** (air-core coil — field + inductance) | ✅ done | `solenoid_coil.py` ; `solve_axi_magnetostatic(Jr=…)` + `solenoid_axial_field` / `solenoid_inductance_long` + `force.inductance_axi` ; `ngsolve_usage("solenoid")` | on-axis B_z **<1.7 %** (centre 0.7 %) ; self-L (2W/I²) = Nagaoka k_N vs Wheeler **0.5 %** |
| 16 | **Busbar Lorentz force** (parallel conductors — fault force) | ✅ done | `busbar_force.py` ; `force.lorentz_force_2d` (J×B integral) on planar A_z ; `ngsolve_usage("busbar")` | F/L vs `μ₀I₁I₂/2πd`, attract+repel, **~1 %** + correct signs |
| 17 | **Helmholtz coil** (uniform field + uniformity) | ✅ done | `helmholtz_coil.py` ; `solve_axi_magnetostatic(Jr=…)` + `coil_pair_axial_field` ; `ngsolve_usage("helmholtz")` | on-axis B_z **<0.5 %**, centre vs `(4/5)^{3/2}μ₀NI/a` 0.38 %, bore uniformity **<1 %** |
| 18 | **Iron-yoke dipole** (window-frame, air-gap field) | ✅ done | `c_magnet_gap.py` ; `solve_planar_magnetostatic` (μ_r yoke) + `magnetic_circuit_gap_field` ; `ngsolve_usage("c_magnet")` | gap B vs reluctance `μ₀NI/(g+l_fe/μ_r)`, **0.5 %** (FEM < ideal: fringing) |
| 19 | **Electro-thermal Joule heating** (multiphysics) | ✅ done | `joule_heating.py` ; chain `solve_current_flow`→`joule_heat_source`→`solve_heat_steady` ; `ngsolve_usage("joule_heating")` | bar rise `dT(x)=(q/2k)x(L−x)`, peak `σV²/8k`, **0.00 %** (exact) |
| 20 | **Electro-thermo-mechanical** (3-physics chain) | ✅ done | `electro_thermo_mech.py` ; NEW `radia_ngsolve.elasticity` + chain to `solve_linear_elasticity` ; `ngsolve_usage("elasticity")` | constrained-bar thermal stress `σ_xx=−Eα⟨ΔT⟩` + tension `δ=σ₀L/E`, **0.00 %** |
| 21 | **Magneto-mechanical** (Lorentz force → deflection) | ✅ done | `magneto_mechanical.py` ; `solve_linear_elasticity(body_force=J×B)` + `cantilever_tip_deflection` ; `ngsolve_usage("magneto_mechanical")` | cantilever tip vs Euler-Bernoulli `wL⁴/8EI`, **0.02 %** |
| 22 | **MEMS electro-mechanical** (electrostatic pull → deflection) | ✅ done | `mems_electro_mechanical.py` ; NEW `force.electrostatic_eggshell_force_2d` + chain to `solve_linear_elasticity` ; `ngsolve_usage("electro_mechanical")` | pressure `½ε₀(V₀/d)²` **0.00 %** + cantilever tip vs Euler `wL⁴/8EI` **0.05 %** (electric twin of #21) |
| 23 | **Frozen-permeability superposition** (saturated Ld/Lq engine) | ✅ done | `frozen_permeability.py` ; NEW `solve.frozen_reluctivity` ; `ngsolve_usage("frozen_perm")` | nonlinear λ decomposes into superposable frozen-ν parts, **0.51 %** recombination (iron μ_r 1000→45); self-validating, no external data |

**AC/DC canonical series** (task #24): #5 onward reproduce canonical AC/DC
field-analysis models, each cross-checked against a closed form.

**ALL of #5–#12 additionally confirmed against an independent reference solver**
(2026-06-05) — a complete 8/8 sweep, every one a three-way agreement
reference == analytic == NGSolve:
#5 sphere shielding (S=16.6, 1.5%), #6 spherical cap (0.03%), #7 dielectric cap (0.7%),
#8 Halbach (0.832 T uniform, 0.0%), #9 capacitance matrix (exact), #10 electrostatic
force (2.0%), #11 cylinder shielding (0.9%), #12 two-wire line (0.12%). Spans
electrostatics, dielectrics, capacitance matrices, electrostatic force, permeable
magnetostatics, and permanent magnets. The reference numbers are kept as a stored
regression reference (the reference setup is lab-private). So these reproductions
are validated against both a closed form *and* an independent reference solve.

**Already in the lab** (not re-done here, per "avoid reinventing the wheel"):
*skin effect / AC resistance* of a round wire — `solve_planar_eddy` current-driven,
validated vs the Kelvin ber/bei formula to 0.07 % in `tests/test_planar_eddy.py`.

## Key reusable insight captured (#1)

In a **scattered-field** harmonic eddy solve the source is a background
potential `A0` and the solved `gfA` is the *scattered* potential -- field probes
write `curl(gfA) + B0`. So the total electric field is `E = -j w (A0 + gfA +
grad Phi)`; **omitting A0 overestimates the Joule loss ~10x.** This is baked into
`joule_loss_density(..., A0=...)` so it cannot be forgotten. (The TEAM-21 tests
use a coil source -> `gfA` already total -> `A0=None`; do not copy that for
background-field problems.)

## #3 is intentionally NOT re-implemented

The Sugahara-lab `radia` package already ships a complete, validated Kelvin
open-boundary toolkit (true periodic-BC inversion, not concentric shells). The
right move per the lab's own "avoid reinventing the wheel" policy is to reuse it
(`radia.kelvin_solver`, recipe in the `kelvin` knowledge module) rather than
rebuild a worse copy. The reluctivity rule (HCurl A): `nu_kelvin = (r'/R)^2 nu_0`.

**Companion (why #3 works on coarse meshes), the spectral angle, not a re-impl:**
[`../dtn_spectrum_coarse_mesh_demo.py`](../dtn_spectrum_coarse_mesh_demo.py)
explains Kameari's coarse-mesh accuracy of the Kelvin transform as a *property of
the DtN operator*: every open-BC closure (Kelvin/BEM/PML) discretizes the one
exterior Steklov-Poincaré operator Λ_ext, whose sphere eigenvalues are the
mesh-independent ladder −(n+1)/R. The demo **measures both discretizations**:
(A) the BEM matrix Λ_h reproduces the low (physically-dominant) modes accurately
already on the coarsest mesh (dipole −2/R to 0.07 %), so refining only widens the
accurate band; (B) `kelvin_dtn_eigenvalue` measures the Kelvin closure's effective
DtN by volume FEM — mode n inverts to a degree-n polynomial, exact iff FEM
order≥n, so the dominant dipole (→ linear) is coarse-accurate at order 1. Knowledge
baked into the `dtn_coarse_mesh` MCP tool and gated by
`tests/test_dtn_spectrum_coarse.py`. It does **not** reimplement Kelvin — it
measures the effective DtN of both the BEM and Kelvin realizations of the shared Λ_ext.
