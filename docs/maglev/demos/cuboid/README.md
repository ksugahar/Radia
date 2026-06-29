# cuboid/ — CLN-SIBC 3D cuboid core

The non-axisymmetric building block of the 3D levitation route: a cuboid
with a != b != c yields direction-split eddy time constants and hence a
diagonal polarizability tensor.  A single self-contained Python script
ports the lab CLN-SIBC (rank-(1,1) Mixed Galerkin) cuboid model so it is
runnable in-repo with no NGSolve / no Radia dependency.

| Script | Purpose | Run | Headline result |
|--------|---------|-----|-----------------|
| `cln_sibc_cuboid_3d.py` | Foster admittance (modal decomposition) + CLN/Cauer reduction + Schur-F SIBC surface termination -> polarizability `alpha(s) = V - Y(s)/sigma` | `python cln_sibc_cuboid_3d.py` (seconds; pure numpy) | 5x2x1 mm Cu (sigma=5.8e7, mu=mu0): 343 modal modes (m,n,p odd up to 13); `alpha(0)=0`, `alpha(inf)=V`; CLN reduces to a few-stage FOS; Schur-F termination `Y_R(s) = Y_CLN(s) + K_SIBC sqrt(s)/(s+d)` closes the modal series |

**Dependencies**: none — standalone modal/reduction calculation.

**Outputs**: console only (no JSON/PNG).

This is the curated, runnable demo.  The broader 3D cuboid research
(Schur-F / CF / AAA derivations, NGSolve FEM cross-checks, the
`cuboid_521_*` Kameari-Kelvin investigation) lives in
`../../research/cln/`.
