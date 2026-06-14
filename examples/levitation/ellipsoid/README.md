# ellipsoid/ — shape-anisotropic polarizability tensor

The shape-anisotropic polarizability tensor of a conducting ellipsoid,
anchored by two analytic limits — DC (`alpha = 0`) and high frequency
(`alpha_i = -V/(1-N_i)`, with `N_i` the Osborn demagnetizing factor) —
and bridged by full-frequency FEM eddy-current solves.  The anisotropy is
**geometric** (semi-axes a1 != a2 != a3), not material: copper stays a
scalar sigma/mu.  A field along axis i drives eddy currents in the
cross-section perpendicular to i, so unequal axes give three different
eddy time constants and hence a direction-split `alpha = diag(alpha_x,
alpha_y, alpha_z)`.  This is the 3D non-axisymmetric generalization of
the sphere: a body lifts most when the field is along its shortest axis
(largest flux exclusion per unit volume).

| Script | Purpose | Run | Headline result |
|--------|---------|-----|-----------------|
| `ellipsoid_alpha_tensor.py` | Analytic tensor: DC limit, HF perfect-conductor limit `-V/(1-N_i)`, orientation-dependent lift | `python ellipsoid_alpha_tensor.py` (seconds) | sphere `N_i=1/3`; triaxial 5x3x1.5 mm -> kappa = (-109,-130,-228) mm^3, anisotropy `|kappa_z|/|kappa_x| = 2.09`; sum `N_i = 1` to <1e-10 |
| `ellipsoid_alpha_omega_axisym.py` | Full-frequency axial `alpha_c(omega)` of a spheroid by axisymmetric FEM; validates vs analytic sphere `G(x)`, anchors HF at `-V/(1-N_c)` | `python ellipsoid_alpha_omega_axisym.py` (~1-2 min) | FEM sphere reproduces analytic `4 pi a^3 G(x)` to 1.8% over 2-200 kHz; HF error sphere 1.3% / prolate 0.7% / oblate 3.0% |
| `ellipsoid_alpha_tensor_3d.py` | Transverse (m=1) tensor component via 3D HCurl + CompactAMS; completes the full tensor | `python ellipsoid_alpha_tensor_3d.py` (CompactAMS order-1; ~2 min) | 3D `alpha_xx, alpha_zz` match analytic sphere to ~2-3% (a/delta=1.3, 2.0); isotropy <0.1%; **air-shell mesh resolution is critical** (coarse air -> 23% error, fixed by refinement) |

**Dependencies**: `ellipsoid_alpha_omega_axisym.py` and
`ellipsoid_alpha_tensor_3d.py` import `levitation_sphere_force` from
`../sphere/` (for `G_exact`, `delta`) and `ellipsoid_alpha_tensor` (same
folder).  The cross-folder import is wired with one
`sys.path.insert(0, os.path.join(HERE, "..", "sphere"))` line.

**Outputs**: `*_results.json` + `ellipsoid_alpha_tensor_lift.png` /
`ellipsoid_alpha_omega_axisym.png` (demag factors, kappa, lift ratio,
alpha_c(omega) curves with HF anchor checks).
