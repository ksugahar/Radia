# sphere/ — isotropic levitation force & dipole validation

AC magnetic levitation force on a conducting sphere, from analytic theory
(Landau-Lifshitz) through validated eddy-current solves.  A sphere is
isotropic, so a single scalar polarizability fully describes its
response — this is the simplest real levitation case and the testing
ground for the Radia (open-boundary external field) + NGSolve (eddy
reaction) workflow.  Three independent force methods (analytic dipole,
full axisymmetric eddy FEM, Maxwell surface stress) converge to pin the
coefficient and the validity window (a/L <~ 0.5).

| Script | Purpose | Run | Headline result |
|--------|---------|-----|-----------------|
| `levitation_sphere_force.py` | Induced-dipole AC levitation force (Landau-Lifshitz `G(x)`), coefficient pinned to the perfect-conductor limit, frequency response reduced by CLN/Cauer | `python levitation_sphere_force.py` (seconds; pure numpy) | `F = (pi a^3/mu0) Re[G] grad(B0^2)`; HF limit matches perfect-conductor energy to 0.09%; 6-stage CLN reproduces the full modal response to <0.0001% |
| `coil_sphere_eddy_force.py` | Full axisymmetric eddy-current Lorentz force vs the induced-dipole force; quantifies the a/L ~ 0.5 approximation error | `python coil_sphere_eddy_force.py` (~1 min; one axisym solve) | a/L=0.5 (5 mm): `F_full/F_dipole = 1.019` (+1.9%); a/L=0.10 (1 mm): 0.998 (point-dipole accurate); coil field matches analytic loop to 0.03% |
| `coil_levitation_equilibrium.py` | Real levitation equilibrium: Radia coil field + verified sphere force -> stable height where lift balances gravity | `python coil_levitation_equilibrium.py` (seconds) | 30 mm coil + 5 mm Cu sphere @ 50 kHz: equilibrium z* = 35.2 mm, lift 46.0 mN = weight; restoring dF/dz = -4.75 N/m (~5 Hz bounce) |
| `maxwell_surface_force_xcheck.py` | Validates the Maxwell surface-stress (SIBC) force against the volume Lorentz force and the L-L analytic result on a real eddy problem | `python maxwell_surface_force_xcheck.py` (~1 min; one axisym solve) | volume Lorentz = surface Maxwell = L-L dipole all agree; surface integral is contour-radius-independent (confirms the Maxwell-stress extraction) |

**Dependencies**: `coil_levitation_equilibrium.py`, `coil_sphere_eddy_force.py`,
and `maxwell_surface_force_xcheck.py` import `levitation_sphere_force`
(same folder); the latter two also use Radia for the open-boundary coil
field.  Run `levitation_sphere_force.py` first if anything looks stale.

**Outputs** (committed next to each script for provenance): `*_results.json`
+ `*.png` (force / Re[G] vs frequency, equilibrium + stability, dipole-error
ratio vs a/L).
