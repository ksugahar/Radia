# Field reconstruction

Demonstrates that the SAME state variables used to compute port
admittance Y(s) also yield the full spatial field distribution
$H_z(\mathbf{r},t)$, $J(\mathbf{r},t)$ via linear combination of the
underlying mode functions.

| Script | Geometry | Output |
|---|---|---|
| `field_reconstruction_cylinder.py` | 2D Cu cylinder | $H_z(r,t)$, $J_\varphi(r,t)$ inside cylinder vs exact PDE |
| `field_reconstruction_sphere.py` | 3D sphere | Equivalent radial field reconstruction |
| `field_reconstruction_warburg.py` | Warburg surface-mode contribution | Demonstrates surface boundary-layer mode $e^{-(a-r)/\delta_k}$ |
| `field_reconstruction_heatmap.py` | Spatial heatmap over (r, t) | $J_\varphi$ heatmap |
| `fem_krylov_field_recon.py` | 3D cube via FEM Krylov basis | NGSolve Krylov $Q$ + Foster modes → 3D field reconstruction |

These are the key "feature" of hierarchical Cauer over standard ROMs:
typical ROM gives port quantities only, but the Warburg-Schur
construction's bulk Krylov basis + surface Warburg modes carry the
full spatial information.
