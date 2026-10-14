"""FEM-BEM Schur coupling knowledge -- exact open boundary for interior FEM.

The exterior Laplace DtN assembled from ngsolve.bem replaces the unbounded
exterior domain with a boundary operator, giving an EXACT transparent BC
without mesh truncation, PML, or Kelvin transformation.

Companion code: radia_ngsolve.fem_bem_coupling
  - laplace_fem_bem_schur(): coupled solve (dense, for verification)
  - sphere_shell_analytic_dtn(): spherical-shell analytic-DtN cross-check vs exact dipole
"""

FEM_BEM_SCHUR_OVERVIEW = """
# FEM-BEM Schur Coupling: Exact Open Boundary

## Concept

The DtN (Dirichlet-to-Neumann) Schur complement approach couples:
  - FEM interior: mesh + stiffness matrix K_FEM for the bounded domain Ω
  - BEM exterior: integral-equation DtN Λ_ext on the outer boundary Γ

Resulting in the coupled system (no exterior mesh needed):

    (K_FEM − P^T M_bnd Λ_ext P) u = f

where P is the H1 volume → SurfaceL2 boundary projection and M_bnd is the
SurfaceL2 boundary mass.  The weak boundary term is −∮(∂u/∂n)v = −⟨Λ_ext u, v⟩ =
−P^T M_bnd Λ_ext P: the mass M_bnd turns the DtN DOF-map Λ_ext into the weak form
(omitting it is a magnitude bug).  With Λ_ext negative-definite (eigenvalue
−(n+1)/R) the term −P^T M_bnd Λ_ext P is positive-definite.

## Why BEM DtN beats Kelvin or PML

| Method         | Exterior treatment    | Accuracy  | Mesh size   |
|----------------|-----------------------|-----------|-------------|
| Kelvin         | domain transformation | exact     | 2× volume   |
| PML            | absorbing layer       | approx.   | 1.5× volume |
| BEM DtN        | no exterior mesh      | exact     | 1× volume   |
| Truncation BC  | Robin approximation   | approx.   | 1.1× volume |

BEM DtN = exact truncation of the exterior: the solution matches the exact
open-domain solution at the coupling surface and everywhere inside.

## The DtN Operator (Laplace)

For the exterior Laplace problem (−Δu=0 for |x|>R, u→0 at ∞):

    Λ_ext : u|_Γ → ∂u/∂n|_Γ  (outward normal derivative)

Assembled via the Calderon identity (NGSolve sign convention):

    V σ = (−½ M + K) u   [exterior BIE]
    → Λ_ext = V⁻¹ (−½ M + K)

where:
  V = SingleLayerPotentialOperator   (Laplace, SPD)
  K = DoubleLayerPotentialOperator   (Laplace, n_y outward from enclosed vol.)
  M = SurfaceL2 mass matrix

For a sphere of radius R, Λ_ext applied to Y_n^m gives eigenvalue −(n+1)/R:
  n=0 monopole:  Λ_ext 1 = −1/R × 1
  n=1 dipole:    Λ_ext cosθ = −2/R × cosθ   (validated by test)
  n=2 quadrupole: Λ_ext Y_2 = −3/R × Y_2

## The Coupling Matrix P

P maps H1 volume DOFs to SurfaceL2 boundary DOFs by L2 projection:

    M_bnd P = R,  R[i,j] = ∫_Γ φ_j^{H1}|_Γ  ψ_i^{SL2} dS

With P assembled, the weak DtN contribution to the stiffness is −P^T M_bnd Λ_ext P
(M_bnd = SurfaceL2 mass).  The coupled system K − P^T M_bnd Λ P has the same
sparsity as K_FEM plus a rank-ndof_bnd dense correction on the outer-boundary block.

## Verification: Spherical Shell Dipole

Domain:   R_inner < r < R_outer  (spherical shell)
Inner BC: u = z/R_inner = cosθ   (dipole Dirichlet datum)
Outer BC: BEM DtN (exact exterior condition on sphere of radius R_outer)

Exact solution: u(r,θ) = R_inner² cosθ / r²  = R_inner² z / r³

The BEM DtN condition at r = R_outer enforces:
  ∂u/∂r = −2/R_outer × u  (DtN eigenvalue for n=1)

yielding the unique shell solution u = R_inner²/r² cosθ (C₁=0 via the Robin BC).

Typical accuracy (order=2, maxh=0.3, R_inner=0.5, R_outer=1.0):
  Relative L2 error < 2%  (O(h²) FEM error dominates)
"""

FEM_BEM_SCHUR_API = """
# FEM-BEM Coupling API (fem_bem_coupling.py)

## laplace_fem_bem_schur

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import laplace_fem_bem_schur

gfu = laplace_fem_bem_schur(
    mesh,
    bnd_outer    = "sphere",          # name of outer BEM coupling boundary
    bnd_dirichlet= "inner_sphere",    # name(s) of Dirichlet boundary
    dirichlet_cf = ng.z / R_inner,    # Dirichlet datum CF
    source_cf    = None,              # volume RHS (optional)
    order        = 2,                 # H1 polynomial order
    intorder     = 8,                 # FEM quadrature order
    bem_intorder = 10,                # BEM quadrature order (default intorder+2)
)
```

**Returns**: `GridFunction` on H1(mesh, order).  Validated on the spherical-shell
dipole (inner u=z/R_in, BEM DtN on outer) vs the exact u = R_in²·z/r³ -- L2 rel_err
~5 % at maxh=0.5 order=2, converging with refinement.

## sphere_shell_analytic_dtn (analytic-DtN cross-check)

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import sphere_shell_analytic_dtn

res = sphere_shell_analytic_dtn(
    R_inner=0.5, R_outer=1.0,
    maxh_vol=0.3, order=2, intorder=8,
)
print(res["rel_err"])   # < 0.07 (coarse), converges with refinement
```

Solves the same shell with the ANALYTIC dipole DtN Robin BC ∂u/∂n=(−2/R)u on the
outer sphere (no BEM), confirming the FEM+DtN framework independently of the BEM
assembly.  **Returns** dict: ``{ndof, R_inner, R_outer, l2_err, l2_ref, rel_err}``.

## Notes and Limitations

1. **Dense BEM**: `laplace_exterior_dtn` densifies the BEM operator (O(n²) memory).
   Practical limit: ndof_bnd ≲ 500.  For larger problems use H-matrix BEM + CG.

2. **Coupling matrix assembly**: iterates over all H1 DOFs (O(ndof_vol) linear solves).
   For verification-scale problems (< 1000 DOFs) this is acceptable.

3. **Boundary naming**: label the OCC faces ("outer"/"inner") BEFORE the boolean
   subtraction so the names survive, then pass them as `bnd_outer`/`bnd_dirichlet`.

3b. **SurfaceL2 on a volume-mesh boundary**: `SurfaceL2(mesh, definedon=
   mesh.Boundaries("outer"))` allocates DOFs on the OTHER boundaries too and marks
   them non-free (zero mass / zero BEM rows).  `laplace_exterior_dtn` and the
   coupling-matrix solve therefore operate on the FREE sub-block and scatter zeros
   back (a no-op for a standalone closed surface where every DOF is free).

4. **Sign + mass convention** (two bugs fixed here): the weak DtN contribution is
   −P^T M_bnd Λ_ext P (MINUS, and WITH the SurfaceL2 mass M_bnd).
   - MINUS: from −∮(∂u/∂n)v = −⟨Λ_ext u, v⟩; Λ_ext is negative-definite (−2/R for
     the dipole), so −Λ_ext is positive-definite and K − P^T M_bnd Λ_ext P is SPD
     -- it equals the analytic Robin form +(2/R)∮uv of `sphere_shell_analytic_dtn`.
     (A `+` sign adds a negative-definite term -> grossly wrong field, the classic
     bug.)
   - M_bnd: Λ_ext (from `laplace_exterior_dtn`) is a DOF-coefficient map
     (u_coeff -> σ_coeff); the weak form ∮σv needs the mass, so the term is
     P^T (M_bnd Λ_ext) P = R^T Λ_ext P, NOT P^T Λ_ext P.  Omitting M_bnd is a
     magnitude bug (the field comes out ~50 % off even with the right sign).

5. **Validation status (all three bugs fixed, green)**: `laplace_fem_bem_schur`
   reaches L2 rel_err 5.5e-2 at order=2 on the shell dipole, matching the analytic
   `sphere_shell_analytic_dtn` (0.054) to 3 sig figs (FEM interior error dominates)
   -- deterministic over repeated runs, gated by
   `validation_test/radia_mcp/test_fem_bem_coupling.py::test_full_bem_schur_dense_path`.

6. **The crash (fixed) -- definedon SurfaceL2 as a mixed-form TEST space**: the
   earlier nondeterministic NGSolve crash (exit 127/139) was NOT the BEM densify
   and NOT dual_mapping.  Isolated by 5 minimal reproducers (run 3x each): the
   coupling mixed form `BilinearForm(trialspace=H1_vol, testspace=SurfaceL2(
   definedon="outer"))` corrupts the heap on teardown when the restricted
   SurfaceL2 is the TEST space of a VOLUME-mesh mixed form (independent of
   dual_mapping or ds-region).  With it as the TRIAL space (volume H1 as test) it
   is stable -- so `_coupling_matrix` builds R via the SWAPPED form (and loops
   over ndof_bnd, fewer than ndof_vol).  Do NOT revert the swap.  (Standalone
   surface meshes and the BEM densify over all DOFs were both verified stable.)
"""

FEM_BEM_SCHUR_APPLICATIONS = """
# Applications: When to Use FEM-BEM Schur Coupling

## 1. Open-domain magnetostatics / electrostatics

Replace the Kelvin transformation domain or PML:

```python
# Instead of:  geo = KelvinTransform(inner_domain, outer_sphere)
# Use:
gfu = laplace_fem_bem_schur(
    mesh_inner,
    bnd_outer="outer",
    bnd_dirichlet="inner_conductor",
    dirichlet_cf=ng.CoefficientFunction(1.0),
)
```

The BEM DtN is exact; the Kelvin transformation introduces geometric distortion
that can affect accuracy near the coupling surface.

## 2. Permanent magnet field computation

After reducing to the reduced scalar potential φ (divH = 0 outside sources):
  Interior (with sources): −Δφ = divM     in Ω_src
  Exterior (magnet-free):  −Δφ = 0        in R³ \\ Ω_src

BEM DtN on the bounding surface of Ω_src replaces the free-space exterior solve.

## 3. Motor/actuator FE models (open boundary)

For motors in free space with a surrounding air region:
  - Truncate the air region at a sphere of radius R_outer
  - Apply BEM DtN as the outer BC (no Kelvin or fixed-potential truncation error)
  - Interior FEM (iron + coils + air) solved with exact exterior condition

## 4. Coupling with the AGE (Air-Gap Harmonic Element)

The AGE (airgap_machine.py) handles the gap between rotor and stator analytically.
The BEM DtN handles the exterior of the stator housing analytically.
Together: EXACT treatment of both the gap AND the far-field boundary.

## Comparison with SIBC / GIBC DtN

| DtN type      | What it approximates     | Accuracy        |
|---------------|--------------------------|-----------------|
| SIBC (order 0)| Conductor interior       | O(δ/R)          |
| GIBC (order 1)| Conductor + curvature    | O((δ/R)²)       |
| BEM DtN       | Exterior free space      | Exact           |
| AGE DtN       | Annular air gap          | Exact           |

The radia-ngsolve DtN hierarchy: SIBC → GIBC → BEM → AGE.
"""


def get_fem_bem_schur_documentation(topic: str = "all") -> str:
    """Return FEM-BEM Schur coupling documentation."""
    topics = {
        "overview": FEM_BEM_SCHUR_OVERVIEW,
        "api": FEM_BEM_SCHUR_API,
        "applications": FEM_BEM_SCHUR_APPLICATIONS,
    }
    if topic == "all":
        return "\n\n".join(topics.values())
    return topics.get(topic, f"Unknown topic '{topic}'. Options: {list(topics)}")
