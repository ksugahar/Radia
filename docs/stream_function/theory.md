# Theory — Stream-function method + (ACA+)+TSVD

> Runnable companion: [`theory.ipynb`](theory.ipynb) executes this method live
> (the derivations/literature stay here; the notebook produces the numbers).

## The inverse problem

Given a desired magnetic field component (typically `Bz`) at a set of
**target points** `r_j ∈ R³`, design a current-carrying source surface
that produces this field as closely as possible.

For a surface current density `K` on a source surface `Γ`:

    B(r_target) = (μ₀ / 4π) ∫_Γ K(r') × (r_target − r') / |r_target − r'|³ dA(r')

If `K` is **divergence-free** (no current sources / sinks on `Γ`), then
on a 2D surface there exists a scalar **stream function** `ψ : Γ → R`
such that `K = n̂ × ∇ψ`.  The level sets of `ψ` are then the
streamlines of the surface current, which correspond directly to
**equal-current wire loops** when discretised.

## Discretisation

Two equivalent paths, both shipped:

### Path 1: basis-loop discretisation (legacy SA-25-020)

Place `N` small quad current loops on a regular `(φ, z)` or `(x, y)`
grid covering `Γ`.  Each loop carries unit current; the optimisation
finds the amplitudes `φⱼ`.  The stream function reconstruction is just
the grid sampled value `ψ(grid_pt) = φⱼ`.

Pros: straightforward, no FE machinery.
Cons: discrete grid spacing artefacts; ψ topology jumps under small
perturbations (= why Path-A iteration only oscillates with this basis).

### Path 2: FE-direct ψ (recommended; new 2026-05-30)

`ψ ∈ H¹(Γ)` is a continuous finite-element GridFunction on a Netgen
triangulation of `Γ`.  For each basis function `φᵢ` and each target
`r_j`, assemble:

    A[j, i] = (μ₀ / 4π) ∫_Γ (n̂ × ∇φᵢ(r')) × (r_j − r') / |r_j − r'|³ · ẑ dA(r')

which on a flat plane reduces to:

    A[j, i] = − (μ₀ / 4π) ∫_Γ ∇φᵢ(r') · (r_j − r')_xy / |r_j − r'|³ dA(r')

In NGSolve this is one `LinearForm` per target.  See
[`build_fem_matrix`](examples_catalog.ipynb)
for the canonical (planar) assembly.

**Γ can be ANY surface, not just a plane.** Mesh a curved former
(cylinder, **sphere**, conformal / 3D-printed) and assemble the *same*
general kernel (line above) with the surface gradient `grad(u).Trace()`
and the surface normal `specialcf.normal(3)`; surface FEM is done
robustly via a solid mesh + `H1(mesh, order=p, definedon=mesh.Boundaries(".*"))`
+ `ds` + `mesh.Curve(p)` (isoparametric).  **This is exactly the capability
the basis-loop path lacks** — it needs a structured `(φ, z)` grid and cannot
be laid on a sphere.  See
[`demo_sphere_fe_direct.py`](examples_catalog.ipynb)
(uniform Bz: cres 3e-15, single-stroke 0.24 %; Z2 shim: 0.36 % after
sheet-metal distortion).

Pros: continuous ψ, smooth contour family, Path-A iteration converges
**monotonically** (on simple topologies), arbitrary polynomial order via
`H1(mesh, order=p)`, **arbitrary former shape**.
Cons: per-entry integration is more expensive than basis-loop kernel
evaluation; without ACA+ the assembly is O(M × ndof × quad_pts).

## The (ACA+)+TSVD least-norm solver

Forming the dense `M × N` matrix `A` is `O(M N)` (basis-loop) or
`O(M × ndof × quad)` (FE-direct) — both expensive for large M or
expensive kernels (material).  Two-stage compression:

### Stage 1 — ACA+ (Adaptive Cross Approximation Plus)

Greedy pivot-based low-rank factorisation:

    A ≈ C · Dᵀ        C ∈ R^{M × k_aca},  D ∈ R^{N × k_aca}

evaluating only `~k_aca · (M + N)` entries of `A` — never the full
`M·N`.  Implemented in `HACApK` (`src/ext/HACApK/cHACApK_acaplus`), MIT
licensed, ppOpen-HPC project.

Radia delegates ACA+ entirely to HACApK via the entry-function override
`HACApK_set_entry_func` (default behaviour, the reusable material interaction matrix,
is unchanged when the override is null).  No ACA+ re-implementation in
Radia.

### Stage 2 — TSVD recompression (manuscript Method 2/3, SA-25-020)

Recompresses the rank-`k_aca` factors `C, D` to an explicit SVD of the
low-rank approximation:

    Method 3 (default, faster):
        SVD(C) → Uc, Sc, VTc;   E := diag(Sc) · VTcᵀ · Dᵀ;
        SVD(E) → UE, SE, VTE;   U := Uc · UE,  S := SE,  V := VTE.

Two SVDs of `M × k_aca` and `k_aca × N` matrices; the full `M × N` SVD
is never formed.

## Least-norm pseudo-inverse

Once we have `A ≈ U · diag(S) · Vᵀ` (truncated to `k_modes` modes):

    ψ̂ = V · diag(1/S) · Uᵀ · B_target

Truncation at `k_modes < k_aca` regularises the otherwise ill-posed
inverse (small singular values amplify noise).  Sweeping `k_modes`
traces the **L-curve** of residual `||A ψ − B||` vs solution norm
`||ψ||`.

## Path-A compensated iteration (new 2026-05-30)

Kuijpers et al. Compumag 2023 [525] observe that **the deviation of the
single-stroke chain from the iso-contour lines is collocated with the
field error in the target region**: parasitic field from the connection
segments degrades the design.  They observe and select-against the
deviation; we **fold it back into the solve**:

    ψ⁽⁰⁾ = pseudo_inverse(B_target)         # initial SF
    repeat:
      χ⁽ᵏ⁾ = single_stroke_chain(ψ⁽ᵏ⁾)     # build coil
      Bz⁽ᵏ⁾ = Bz_at_target(χ⁽ᵏ⁾, I_w)      # actual field
      r⁽ᵏ⁾  = B_target − Bz⁽ᵏ⁾             # residual
      ψ⁽ᵏ⁺¹⁾ = ψ⁽ᵏ⁾ + α · pseudo_inverse(r⁽ᵏ⁾)
    until ||r|| < tol

Each iteration RE-USES the cached TSVD factorisation (one back-
substitution + one chain rebuild) → cheap.

**Convergence depends on the ψ REPRESENTATION**:

  - Basis-loop ψ (grid-sampled, matplotlib contour): topology of the
    contour family JUMPS under small ψ perturbations → `B_c(ψ)` not
    smooth → Picard does NOT contract.  Best-effort tracking finds
    marginally better neighbourhoods.
  - FE-direct H¹ ψ: continuous GridFunction → contour family deforms
    smoothly → `B_c(ψ)` smooth → Picard CONTRACTS.  Iterations 40-47
    on the planar uniform Bz benchmark drop **monotonically**
    0.62 % → 0.49 % RMS, converging to **0.47 %**.

This is the empirical justification for the FE-direct upgrade.

## Complexity tier framework

A coil's reachable design quality (via Kuijpers chain + Path-A) is
bounded by its TOPOLOGY CLASS:

| Tier   | Topology                  | Baseline RMS  | + Path-A      | Behaviour                    |
|--------|---------------------------|---------------|---------------|------------------------------|
| EASY   | axisym. / planar uniform  | 2–3 %         | < 1 %         | Path-A useful                |
|        |   + FE-direct H¹          | 2 %           | **0.47 %**    | MONOTONE convergence         |
| MEDIUM | cylindrical Gz            | already OK    | redundant     | smooth helix natural         |
| HARD   | cylindrical Gx fingerprint| 9.3 % (field_aware) | oscillates | escape via geometry, below |
| HARDER | shielded / biplanar / 3D  | --            | --            | FE-direct on the real former |

Past EASY: FE-direct continuous ψ.  **Past HARD the practical escape is
GEOMETRY, not more current DOFs**: the single-current **sheet-metal wire
distortion** bends the manufactured wire (trade geometric DOF for current
DOF) and takes the cylinder Gx fingerprint from 8.5 % → **1.4 %** with ONE
feed (see [single_stroke.md](single_stroke.md)); a few separate-feed electric
shims then refine to ~1.0 %.  The remaining classical routes are B-spline SFD
(Kuijpers Methods 2/3) or a multivalued-potential D-path (research).  Note
the tier is set by the COIL PATTERN topology, NOT the former shape — FE-direct
solves an EASY uniform target on a *sphere* to 0.24 % single-stroke.

## Cross-reference

  - Math + literature: the paper outline (W:\02_学会資料\2025年度\2026_01_JIAM\streamfunction\)
  - Regularisation choices: [regularization.md](regularization.md)
  - Chain construction: [single_stroke.md](single_stroke.md)
  - Deformation outer loop: [deformation.md](deformation.md)
  - MCP topic: `streamfunction(topic=method)`, `streamfunction(topic=session_2026_05_30)`
