"""Non-conforming mesh coupling (mortar / Nitsche / FETI-DP / BDDC / DG)
for electromagnetic FEM — knowledge module for radia_mcp.fem.

Source folder (Sugahara Lab curated, 33+ PDFs, 6 subcategories):
    W:/03_文献・論文/00_電磁界解析/15_非接合/

The folder organizes literature for the "HFSS Mesh Fusion"-style problem:
how to independently mesh multiple regions of an EM model and couple
them weakly at non-matching interfaces, instead of forcing one global
conformal mesh.

User motivation (2026-05-25):
> HFSS Mesh Fusion を見て XFEM と混同したのを契機に、XFEM (single mesh +
> trial-space enrichment) とは別系統の multi-mesh + interface coupling
> 技術を体系的に勉強したい。

KEY DISAMBIGUATION — two senses of "nonconforming":

| Meaning | Description | Examples |
|---|---|---|
| **Nonconforming FE element** | Within a SINGLE mesh, H¹-conformity (continuity at element boundaries) is relaxed | Crouzeix-Raviart 1973, Wilson 1973, MITC plate |
| **Nonconforming / non-matching mesh** | Multiple SEPARATE meshes, coupled weakly at interfaces (meshes are not aligned at the interface) | Mortar, FETI-DP, BDDC, HFSS Mesh Fusion |

This module is about the **second** meaning (mesh-coupling).  XFEM
(trial-space enrichment within a single mesh) is a separate lineage —
see radia_mcp.fem.xfem_em_hiruma_knowledge and
W:/03_文献・論文/00_電磁界解析/Hiruma_XFEM_CLN/.
"""

# Reference: source folder organization
_FOLDER_LAYOUT = r"""
# Source folder layout

W:/03_文献・論文/00_電磁界解析/15_非接合/

  00_overview/                  4 PDFs   Concept introduction, mixed FE
  01_mortar_nitsche/           12 PDFs   Mortar + Nitsche-type coupling
  02_FETI_DP_BDDC/              9 PDFs   FETI-DP / BDDC / IETI substructuring
  03_HFSS_mesh_fusion/          2 PDFs   Sliding mesh, transfinite mortar
                                          (academic analogs of HFSS Mesh Fusion)
  04_DG_methods/                4 PDFs   Discontinuous Galerkin (mesh-coupling
                                          flavored)
  05_em_applications/          13 PDFs   EM-specific applications (eddy current,
                                          electric machines, sliding rotor-stator)

  README.md                              Hand-curated overview + reading guide

Naming convention: `[FirstAuthor]_[Year]_[Short Topic] ([Venue]).pdf`
"""


_OVERVIEW = r"""
# Overview — when do you need non-conforming mesh coupling?

Conformal mesh:  one global mesh, every element edge matches its
                 neighbor.  Easy to formulate, but:
                 - Hard to mesh heterogeneous geometry (small features
                   force tiny elements everywhere along the shared mesh)
                 - Cannot independently refine one region without
                   re-meshing the whole
                 - Sliding/rotating geometry (rotor-stator) requires
                   re-meshing each time step

Non-conforming / non-matching mesh:
                 - Mesh each region INDEPENDENTLY (best element type,
                   best density, best refinement strategy per region)
                 - Couple at the interface via a WEAK condition
                   (Lagrange multiplier, penalty term, or numerical
                   flux)
                 - Each region can be re-meshed without touching its
                   neighbors

## Decision matrix — pick a method

| Method | Couple via | Pros | Cons | Best for |
|---|---|---|---|---|
| **Mortar** | Lagrange multiplier on interface | Optimal H¹ error, well-developed theory | Saddle-point system, LBB condition, mortar space design | Rotor-stator, IGA, EM machine analysis |
| **Nitsche** | Symmetric/non-sym penalty terms | No LM space, SPD system | Penalty parameter tuning, NOT mass-conservative | CutFEM, time-dependent unfitted, multi-mesh shape opt |
| **FETI-DP / BDDC** | Iterative substructuring + coarse dofs | Massively parallel, condition-number bound | More machinery (primal/dual decomposition), library work | HPC scaling, > 1000 cores |
| **DG** | Numerical flux at element faces | Local conservation, hp adaptivity natural | More DoF per cell, jump-stabilization tuning | Hyperbolic / convection-dominated, Maxwell wave |
| **Sliding mesh / transfinite mortar** | Mortar on a rotating interface | Best for rotating geometry, conservative | Need projection / mortar update each time step | Rotating machinery (CFD, motors) |

## Position vs other lab techniques

- **XFEM** (trial enrichment, single mesh): radia_mcp.fem.xfem_em_hiruma
- **PML** (open boundary, single mesh): radia_mcp.fem.fem_gauge_open_boundary
- **Kelvin transformation** (single mesh, exterior compactified):
  radia_mcp.radia_ngsolve (kelvin_transformation)
- **FEM-BEM coupling** (one part FEM, one BEM via boundary
  integral): radia_mcp.bem.bem_fem_bem_hybrid + ngsbem
"""


_MORTAR = r"""
# Mortar element method — the standard non-conforming coupling

Original idea (Bernardi, Maday, Patera, "AsymNumMethods", 1994):
- Two subdomains Ω1, Ω2 sharing interface Γ.  Each has its own mesh.
- Trial spaces V_h^1, V_h^2 are H¹ in each subdomain.
- Add constraint: ∫_Γ (u_h^1 - u_h^2) μ_h = 0 for all μ_h ∈ M_h
  (the "mortar space" of Lagrange multipliers).
- Saddle-point problem:
      [  A    B^T ] [ u ]   [ f ]
      [  B    0   ] [ λ ] = [ 0 ]
  where A is the block-diagonal stiffness, B is the coupling.

Critical choice — the mortar space M_h:
- Must satisfy LBB (inf-sup) for stability + optimal error
- Standard choice: trace of one side's FE space, with end-point
  modifications (so the LM doesn't over-constrain corners)
- For HCurl (Maxwell): tangential traces; see Ben Belgacem-Buffa-Maday
  (2001) for 3D Maxwell mortar

## Variants

| Variant | Insight | Reference |
|---|---|---|
| **Standard mortar** | LM space = trace of one side | Bernardi-Maday-Patera 1994 |
| **Dual-mortar** | LM space biorthogonal to trace -> local LM elimination | Wohlmuth 2001 |
| **3-field mortar** | Separate trace space, allow non-matching trace approximations | Brezzi-Marini |
| **Isogeometric mortar** | NURBS / B-spline on both sides | Brivadis-Buffa-Wohlmuth-Wunderlich 2014/2015 |
| **Sliding mesh-mortar** | LM space updates at each rotor angle | Buffa-Maday-Rapetti 2001 (EM) |
| **Harmonic mortar** | Fourier basis on rotor interface | Egger et al. 2020/2021 |

## Lab usage pattern (rotor-stator, IH coil-workpiece)

1. Mesh rotor (or coil) region independently — choose density driven by
   skin depth / curvature.
2. Mesh stator (or workpiece) region independently — different density.
3. At the interface, define a mortar Lagrange multiplier space whose
   support is on the interface mesh of ONE of the two sides (typically
   the "slave" with the finer mesh).
4. Assemble the block system [A, B^T; B, 0].
5. Solve as saddle-point (e.g. block-diagonal preconditioner with A
   inverse + Schur complement on λ).

## EM mortar references (chronological)

- Rapetti-Bouillault-Santandrea-Buffa-Maday-Razek 2000 (TMag) — 2D first
- Buffa-Maday-Rapetti 2001 (M2AN, 38 pages) — 2D theory + experiments
- Ben Belgacem-Buffa-Maday 2001 (SIAM JNA) — 3D Maxwell mortar
- Bouillault-Buffa-Maday-Rapetti 2003 (SIAM JSC) — 3D magnetostatics
- Rapetti-Maday-Bouillault-Razek 2002 (TMag) — 3D moving structures
- Antunes et al. 2005 (TMag) — hierarchic + mortar electrical machines
- Bhat-DeFalco-Egger-Schöps-Wohlmuth 2019 (arXiv) — IGA mortar EM
- Egger-Harutyunyan-Loescher-Schöps-Steinbach 2020 (arXiv) — harmonic
- Egger-Harutyunyan-Merkel-Schöps 2021 (arXiv) — harmonic torque
"""


_NITSCHE = r"""
# Nitsche method — penalty-based coupling without Lagrange multipliers

Original (Nitsche 1971): a variational way to impose Dirichlet BC
without Lagrange multipliers.  Modern revival started with
Becker-Hansbo-Stenberg 2003.

Idea (for 2 subdomains Ω1, Ω2 sharing interface Γ):
    a(u, v) = sum_i a_i(u, v)
            - ∫_Γ {n.σ(u)} [v] ds      # consistency
            - ∫_Γ [u] {n.σ(v)} ds      # symmetry (gives SPD)
            + (gamma_h / h) ∫_Γ [u][v] ds  # penalty

where:
    {.}   = average across Γ
    [.]   = jump across Γ
    σ(u)  = relevant flux (gradient for Poisson, etc.)
    gamma_h = penalty parameter (must be large enough; problem-dependent)
    h     = element size on the interface

Properties:
- SPD final system (with the +symmetry term)
- No additional unknown (no LM space)
- Optimal convergence under standard regularity
- WEAK constraint -> mass not exactly conserved across interface
  (this matters for some EM problems, less for others)

## Penalty parameter tuning

The penalty gamma_h must be larger than a problem-dependent constant:
    gamma_h >= C * max_K (mu_K / h_K) * p^2
where p is the polynomial degree.  Too small -> non-coercive system;
too large -> ill-conditioned.  Practical rule: gamma_h = 10 * p^2 * mu.

## When Nitsche beats mortar

- Cut-cell / unfitted meshes (CutFEM): mortar requires the interface
  as a mesh feature; Nitsche works on cut elements directly
- Time-dependent moving interfaces (no need to track LM space at each
  step)
- Shape optimization with topology change (Nitsche tolerates degenerate
  interfaces; mortar saddle system goes singular)

## EM-relevant Nitsche references

- Becker-Hansbo-Stenberg 2003 (M2AN) — standard formulation
- Ager-Schott-Winter-Wall 2018 (arXiv) — Nitsche cut-FE flow + poroelast
- Wang-Yang-Xie 2018 (arXiv) — Nitsche-XFEM (bridges both worlds)
- Dokken-Funke-Johansson-Schmidt 2018 (arXiv) — multi-mesh Nitsche
  shape optimization
- Chernov-Hansbo-Schetzke 2026 (arXiv) — hp Nitsche FEM-BEM (most recent)
"""


_FETI_BDDC = r"""
# FETI-DP / BDDC — iterative substructuring (parallel-first)

Both methods solve the same coupled problem by IT­ER­A­TIVELY
exchanging interface data between substructure solvers.  Designed
for massive parallelism (1000s of cores).

## FETI (Farhat-Roux 1991)

Idea: enforce continuity of u across interface via a Lagrange
multiplier λ.  Eliminate λ to get a reduced system on the interface
("Schur complement").  Iterate.

For elliptic problems: condition number scales as κ = O(H/h)
(H = subdomain size, h = element size).

## FETI-DP (Farhat-Lesoinne-LeTallec-Pierson-Rixen 2001)

"Dual-Primal" variant: keep SOME interface DoFs ("primal") globally
continuous (corners, sometimes edges/faces).  The rest stay as
"dual" Lagrange multipliers.  Hugely improves the conditioning:
    κ = O((1 + log(H/h))^2)
i.e., quasi-optimal (logarithmic dependence only).

## BDDC (Dohrmann 2003)

"Balancing Domain Decomposition by Constraints" — different mechanics
but same condition number bound as FETI-DP.  Often viewed as the
PRIMAL counterpart of FETI-DP.  Unified analysis: Mandel-Sousedík 2007.

## IETI-DP / IETI for IGA

(Hofer-Langer 2015/2016) — same idea, but for isogeometric (NURBS/
B-spline) discretizations.  Each "patch" plays the role of a subdomain.

## When to use FETI-DP/BDDC instead of mortar

- > 100 subdomains, parallel solve is the bottleneck
- Iterative solver acceptable (not direct factorization)
- Substructure interior solves can be done locally + cached
- Many right-hand sides (Schur complement amortizes)

## Lab usage

- Has NOT been the lab's primary direction.  Sugahara Lab's typical
  scale (1-10M DoF) is well-served by direct factorization or
  Compact-AMS preconditioned CG.
- FETI-DP becomes attractive only at > 10M DoF AND > 64 cores.
- Reference value: read Mandel-Sousedík 2007 (unified theory) +
  Klawonn-Rheinbach 2008 (robust 3D elasticity) before considering.

## References

- Farhat-Roux 1991 (IJNME) — FETI original
- Farhat et al. 2001 — FETI-DP
- Dohrmann 2003 (SIAM JSC) — BDDC original
- Mandel-Sousedík 2007 (arXiv) — unified theory
- Šístek-Sousedík-Burda-Mandel 2011 — parallel BDDC for Stokes
- Hofer-Langer 2015/2016 — IETI-DP for IGA
- Klawonn-Rheinbach 2008 — robust FETI-DP for heterogeneous 3D
- Kim-Chung-Wang 2016 — adaptive coarse spaces
"""


_DG = r"""
# Discontinuous Galerkin (DG) — mesh-coupling viewpoint

DG was developed originally for hyperbolic problems but has been
adapted to elliptic and coupled problems.  From a mesh-coupling
perspective:
- Each element is INDEPENDENT (its own polynomial)
- Inter-element continuity is imposed via NUMERICAL FLUX
- The mesh CAN be non-matching at element interfaces

This makes DG a natural fit for the multi-mesh setting.

## Standard DG variants for elliptic (unified by Arnold-Brezzi-Cockburn-Marini 2002)

| Method | Stabilization | Symmetric? | Optimal? |
|---|---|---|---|
| Interior Penalty (IP) | gamma_h jump penalty | Yes (SIPG) / No (NIPG) | Yes |
| Local DG (LDG) | Auxiliary flux | No (typically) | Yes |
| Bassi-Rebay | Lifting operator | Yes | Yes |
| HDG (hybridizable) | Numerical trace eliminated | Yes | Yes + static condensation |

For 2nd-order elliptic, SIPG (Symmetric Interior Penalty Galerkin) is
the de facto standard.  Looks just like Nitsche's method!

## DG-IETI-DP (Hofer-Langer 2016)

Combines DG element coupling with IETI-DP substructuring.  Each
subdomain interior is DG; subdomain coupling is IETI-DP.

## Maxwell DG

Hesthaven-Warburton 2002 (Nodal hp Maxwell on unstructured grids) is
the standard reference.  Lu-Cangellaris 2006 brings it to time-domain
EM.  NGSolve has hcurl-DG built-in (rare to be needed for non-wave
problems; more often used for hyperbolic CFD).

## When to use DG for mesh coupling

- Hyperbolic / convection-dominated (transport, Maxwell wave, fluid)
- h-adaptivity with hanging nodes is natural in DG (no need to
  constrain hanging-node DoFs)
- Local conservation is required (each element conserves its own flux)

## Disadvantages

- 3-4x more DoF than continuous Galerkin (one per element-vertex
  becomes one per element-corner, replicated)
- Penalty tuning still required
- For elliptic / static problems, the DoF overhead often outweighs
  the flexibility benefit; mortar or Nitsche are preferred

## References

- Arnold-Brezzi-Cockburn-Marini 2002 (SIAM JNA) — unified analysis
  of DG for elliptic
- Hesthaven-Warburton 2002 — Nodal hp DG, Maxwell
- Cockburn-Karniadakis-Shu 2000 — DG methods review volume
- Brezzi-Cockburn-Marini-Süli 2006 — stabilization mechanisms
- Hu-Yu 2005 — Signorini FEM-BEM coupling (DG-flavored)
- Zhao-Chung 2019 (arXiv) — A posteriori error for mortar staggered DG
"""


_HFSS_MESH_FUSION = r"""
# HFSS Mesh Fusion — Ansys's product + academic analogs

Ansys HFSS Mesh Fusion (commercial, ~2018+) is the technique behind
the user's original motivation question.  Industrial workflow:

1. User selects "parts" of an HFSS model (e.g. PCB, package, antenna
   region, housing)
2. Each part is meshed INDEPENDENTLY by HFSS's adaptive mesh refiner
   (each with its own driven-frequency criterion)
3. At the interface between parts, HFSS uses a non-conforming
   coupling scheme to glue the solutions
4. Iterative refinement: HFSS can ADD elements per-part based on
   per-part error indicators

Public information on the math underneath is limited (Ansys
proprietary).  Best-guess based on the academic literature:
- Likely TFI (transfinite interpolation) mortar on the interface
- Adaptive p-enrichment per part
- Reduced-basis methods for fast frequency sweep

## Academic analogs (no public Ansys disclosure)

Zhang-Liang 2015 (arXiv 1503.02198), "High-Order Curved Sliding-Mesh
Spectral-Difference":
- High-order method on unstructured curved nonconforming sliding meshes
- Spectral difference (variant of DG) with mortar interface
- Demonstrates the geometric machinery

Zhang-Liang 2020 (arXiv 2009.04400), "Dynamic Transfinite Mortar for
Curved Sliding":
- Mortar projection between curved nonconforming meshes
- Dynamic = mortar updates as one mesh slides relative to the other
- Provably conservative, optimal-rate convergence
- Closest published academic analog of HFSS Mesh Fusion sliding

EM-specific analog: Lee-Volakis "Multimesh techniques in FE analysis"
is often cited as the academic predecessor.  Exact citation not yet
confirmed; likely a IEEE TAP or AP-S Symposium paper from ~2000-2005.

## How to build the analog in NGSolve

NGSolve has all the building blocks:
- Multi-region meshing: each region's mesh can be loaded separately
  (or combined with one mesh.AddRegion)
- Mortar: implement via FacetFESpace + Lagrange multiplier (Schöberl
  has demo notebooks)
- Nitsche: implement via BilinearForm with `dx(skeleton=True)` and
  jump/average operators
- HDG: native HCurlHDG / HDG spaces

See `ngsolve_recipe` topic for a concrete code skeleton.

For an actual lab production analog of HFSS Mesh Fusion, the right
NGSolve recipe is:
1. Use mortar (provably convergent + well-understood for EM, per
   Buffa-Maday-Rapetti lineage)
2. For sliding/rotating geometry, the transfinite-mortar recipe
   from Zhang-Liang 2020 maps to NGSolve directly (Lagrange
   multiplier space updates at each angular position)
3. For static multi-region (PCB + package + antenna), standard
   mortar suffices
"""


_NGSOLVE_RECIPE = r"""
# NGSolve recipe — mortar coupling for 2 regions

The minimum-viable NGSolve code to couple two independent meshes
via a mortar Lagrange multiplier.  Adapt for your physics by
replacing the bilinear form `a(u, v)`.

## Step 1: prepare two non-matching meshes

```python
from ngsolve import *
from netgen.occ import OCCGeometry, Box, Sphere, Glue

# Region A: a refined region (e.g. coil / rotor / antenna)
boxA = Box((0, 0, 0), (1, 1, 1))
geo_A = OCCGeometry(boxA)
mesh_A = Mesh(geo_A.GenerateMesh(maxh=0.1))

# Region B: coarser region (e.g. air / stator / housing)
boxB = Box((1, 0, 0), (2, 1, 1))
geo_B = OCCGeometry(boxB)
mesh_B = Mesh(geo_B.GenerateMesh(maxh=0.3))

# IMPORTANT: in real workflow, you'd mesh A and B with shared
# interface BOUNDARY but DIFFERENT element sizes / patterns.
# NGSolve cannot directly mortar two completely separate Mesh
# objects yet; the standard pattern is to mesh both into ONE
# Mesh object with a "soft" interface flag, then mortar the
# interface DoFs.
```

## Step 2: H¹ mortar (for Poisson / scalar Laplacian)

```python
fes_A = H1(mesh_A, order=2, dirichlet="left")
fes_B = H1(mesh_B, order=2, dirichlet="right")
# Lagrange multiplier on interface — choose ONE side's trace
# (the finer mesh's, typically).  Use a "L2 of order p-1" mortar
# space to satisfy LBB.
fes_lm = SurfaceL2(mesh_A, order=1, definedon=mesh_A.Boundaries("interface_A"))

fes = FESpace([fes_A, fes_B, fes_lm])
u_A, u_B, lam = fes.TrialFunction()
v_A, v_B, mu  = fes.TestFunction()

a = BilinearForm(fes, symmetric=True)
# Subdomain stiffness (each region has its own physics)
a += grad(u_A) * grad(v_A) * dx_A
a += grad(u_B) * grad(v_B) * dx_B
# Mortar coupling: lam * (u_A - u_B) on interface
# (in practice, project u_B's trace to the same interface mesh first)
a += lam * (u_A.Trace() - u_B.Trace()) * ds("interface")
a += mu  * (u_A.Trace() - u_B.Trace()) * ds("interface")  # symmetric

f = LinearForm(fes)
f += source_A * v_A * dx_A
f += source_B * v_B * dx_B

a.Assemble()
f.Assemble()
u_sol = GridFunction(fes)
u_sol.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
```

## Step 3: Nitsche variant (no LM space, SPD system)

```python
# Same fes_A, fes_B (no fes_lm)
fes = FESpace([fes_A, fes_B])
u_A, u_B = fes.TrialFunction()
v_A, v_B = fes.TestFunction()

h    = specialcf.mesh_size
n    = specialcf.normal(mesh_A.dim)
mu   = CoefficientFunction(1.0)         # diffusion coefficient
gamma= CoefficientFunction(10 * mu)     # Nitsche penalty (problem-tuned)

a = BilinearForm(fes, symmetric=True)
a += mu * grad(u_A) * grad(v_A) * dx_A
a += mu * grad(u_B) * grad(v_B) * dx_B
# Nitsche interface terms (symmetric IPG)
jump_u = u_A.Trace() - u_B.Trace()
jump_v = v_A.Trace() - v_B.Trace()
avg_dn_u = 0.5 * (mu * grad(u_A) * n + mu * grad(u_B) * n)
avg_dn_v = 0.5 * (mu * grad(v_A) * n + mu * grad(v_B) * n)
a += -avg_dn_u * jump_v * ds("interface")  # consistency
a += -avg_dn_v * jump_u * ds("interface")  # symmetry
a += (gamma / h) * jump_u * jump_v * ds("interface")  # penalty
```

## Step 4: HCurl mortar for Maxwell (eddy current)

```python
fes_A = HCurl(mesh_A, order=2)
fes_B = HCurl(mesh_B, order=2)
# LM space: tangential trace
fes_lm = TangentialFacetFESpace(mesh_A, order=1,
                                  definedon=mesh_A.Boundaries("interface_A"))

# ... (analog to Step 2 but with curl-curl form)
a += nu * curl(u_A) * curl(v_A) * dx_A       # 1/mu (curl curl A)
a += nu * curl(u_B) * curl(v_B) * dx_B
a += 1j * omega * sigma_A * u_A * v_A * dx_A  # eddy current
a += 1j * omega * sigma_B * u_B * v_B * dx_B
a += lam.Trace() * cross(n, (u_A.Trace() - u_B.Trace())) * ds("interface")
```

## Caveats

- NGSolve's Periodic + mortar interaction needs careful setup
- For rotating mesh (motor analysis), the LM space + assembly happens
  per time step; the AC harmonic-mortar variant (Egger 2020) is far
  more efficient
- Saddle-point preconditioning: use a block diagonal preconditioner
  with `radia.sparsesolv_ngsolve.CompactAMSPreconditioner` on the
  HCurl blocks + Schur-complement approximation for the LM block
- For verification: compare against a SINGLE-MESH conformal solution
  on a refined common mesh; mortar should agree within O(h^p) of the
  better-side discretization

## Reference implementations to study

- NGSolve official mortar demo (in i-tutorials): basic Poisson mortar
- Bhat-DeFalco-Egger-Schöps-Wohlmuth 2019 paper (isogeometric mortar
  EM) — explicit code structure for HCurl mortar
- Buffa-Maday-Rapetti 2001 — has Matlab pseudo-code in appendix for
  the 2D sliding-mesh recipe

## Worked example (VERIFIED): docs/mesh_fusion/mesh_fusion.ipynb

A self-contained, runnable 5-phase study (code + committed results_phaseN.json
+ mesh_fusion_convergence.png) in the Radia repo at
`docs/mesh_fusion/mesh_fusion.ipynb`:

- Phase 1 -- Nitsche baseline, 2D Poisson, two half-domains at different mesh
  density. VERIFIED order=2 convergence: L2 err 4.06e-04 -> 4.56e-07 over
  4 refinements (rate ~ +3.1, i.e. O(h^3)), n_dof 59 -> 3391.
- Phase 2 -- mortar / Nitsche H1 Poisson with two INDEPENDENT trial spaces.
- Phase 3 -- 2D Maxwell (A_z), Nitsche coil-air interface, complex
  j*omega*sigma eddy-current term.
- Phase 4 -- Egger (2020) harmonic-mortar (rotor-stator).
- Phase 5 -- accelerator shim-yoke: refine the shim without inflating yoke DoF,
  with mesh-independent field-at-centre recovery across the non-matching joint.

Promoted from the former `examples/mesh_fusion/` on 2026-06-26; the notebook is
the canonical runnable artifact for this topic.
"""


_EM_APPLICATIONS = r"""
# EM applications — domain-specific recipes

## 2D rotor-stator eddy current (Buffa-Maday-Rapetti 2001)

Geometry: stationary stator ring + rotating rotor core.  Interface
is the air gap at the rotor outer radius.

Formulation (2D quasi-static):
    -∇.(nu ∇A_z) + j omega sigma A_z = J_z         (TM eddy current)

Mortar coupling:
- Stator A_z and rotor A_z each on independent mesh
- Lagrange multiplier on the air-gap circle enforces continuity
- LM space updates as rotor angle changes
- Optimal H¹ error preserved across the interface

Numerical experiments in the paper:
- Up to 20 stator slots, rotor angle scan
- Mortar solution agrees with single-mesh ref to ~1e-3 in L²

## 3D eddy current in moving structures (Rapetti et al. 2002)

Extends Buffa-Maday-Rapetti 2001 to 3D:
- A-V formulation (vector potential A + scalar V)
- Tetrahedral mesh both sides; interface is a curved surface
- Hierarchic edge elements on each side; mortar in tangential trace
- 3D rotating-machine demo: induction motor with skewed rotor bars

## IGA mortar for EM (Bhat-DeFalco-Egger-Schöps-Wohlmuth 2019)

Use NURBS / B-spline trial spaces.  Advantages:
- Exact curve representation (no geometry discretization error at
  curved interface)
- Higher p comes for free
- Matches the IGA convergence theory (Pechstein et al. 2018)

Code structure: GeoPDEs (Matlab) or G+Smo (C++).  No NGSolve native
implementation yet (NGSolve has HCurl but not NURBS-HCurl out of the
box; would need custom shape functions).

## Harmonic mortar (Egger et al. 2020, 2021)

For PURELY HARMONIC rotation (constant angular velocity), the rotor's
A_z field on the interface is a Fourier series:
    A_z(r, theta, t) = sum_k a_k(r) exp(j(k theta - omega t))

Replacing the LM space with a TRUNCATED FOURIER BASIS on the interface:
- Decouples the rotor harmonics
- Each harmonic is a separate small-rank update
- Avoids per-time-step mortar reassembly
- ESPECIALLY GOOD for steady-state harmonic motor analysis

This is the modern preferred method when the rotation is uniform.

## Magnetic brake (De Gersem 2010)

Combined spectral-element / finite-element discretization for a
magnetic brake (eddy current under rotating PM):
- SE in air gap (smooth Fourier-mode treatment)
- FE in iron yoke (jagged geometry)
- Iterative coupling at the interface (similar in spirit to mortar)

## Lab tie-in scenarios

| Lab problem | Coupling at | Method |
|---|---|---|
| IH coil + workpiece (independent meshing) | Coil-workpiece air gap | Mortar (Buffa 2001 lineage) |
| BEM PEEC coil + FEM workpiece | Coil filaments don't mesh | Already FEM-BEM (radia_mcp.bem) |
| Motor rotor-stator | Air gap circle | Harmonic mortar (Egger 2020) |
| Accelerator magnet with iron + Kelvin | Iron-Kelvin interface | Single mesh today; mortar could allow refining iron only |
| Multi-conductor cable (Litz) | Strand-air interfaces | Mortar per strand (De Gersem 2002 multiconductor pattern) |
"""


_LAB_SCENARIOS = r"""
# Sugahara Lab specific scenarios

These are the cases where non-conforming mesh coupling would directly
improve the lab's current workflow.

## 1. Rotor-stator coupling for motor analysis

Status: NOT in lab pipeline yet.  Lab does PMSM cogging analysis with
single conformal mesh (calc_motor_transient.py).

What mortar would unlock:
- Independent rotor / stator mesh density (rotor often needs finer
  for surface PM)
- Easy parameter sweep over rotor angle (no remesh per angle)
- Harmonic mortar for steady-state at single rotation speed

How to build: Egger et al. 2020 harmonic mortar in NGSolve.

## 2. IH coil + workpiece independent meshing

Status: lab uses PEEC for coil (no mesh) + FEM for workpiece (mesh).
Mortar not relevant when one side is BEM/PEEC.

Where mortar HELPS: if the lab moves to volumetric FEM coil (for
strong skin/proximity), mortar lets the coil filament-bundle mesh
refine independently of the workpiece mesh.

Reference: De Gersem-Hameyer 2002 multi-conductor FE eddy current.

## 3. PCB + free-space EM analysis

Status: lab doesn't do PCB EM analysis (LAB Cubit/NGSolve workflow
is mostly motor / IH / accelerator / WPT).

If/when needed: mortar PCB-vs-air is the canonical HFSS Mesh Fusion
use case.

## 4. Moving boundary in IH (phase change / billet sliding)

Status: lab does IH on stationary workpiece.

If/when needed: sliding-mesh mortar (Zhang-Liang 2020 transfinite
mortar style) handles moving billet through induction coil.

## 5. Accelerator magnet end-pole shimming

Status: lab uses parametric remeshing for shim sensitivity (radia
Coil + Cubit hex iron).

What mortar would unlock: refine ONLY the shim region without
re-meshing the entire iron yoke.  Faster sensitivity sweeps.

## 6. BEM-CLN multi-conductor refinement

Status: lab uses single-mesh BEM-FEM coupling
(radia_mcp.bem.bem_fem_bem_hybrid).

If multi-mesh BEM coupling is needed (e.g., each transformer winding
with its own surface mesh), mortar on BEM panels is possible
(De Gersem-Weiland series).

## Priority for lab adoption

| Scenario | Effort | Payoff | Priority |
|---|---|---|---|
| Rotor-stator harmonic mortar | High (3-4 weeks NGSolve) | Big — direct motor speedup | HIGH if motor work continues |
| IH FEM coil mortar | Medium | Small — PEEC already handles this | LOW |
| Sliding-mesh IH | High (research-grade) | Niche | LOW |
| Accelerator shim mortar | Medium | Big for ATTO Tower / shimming | MEDIUM |

Recommendation: pursue rotor-stator harmonic mortar (Egger 2020)
as the first NGSolve mortar implementation when motor research
crosses the "more rotor angles needed" threshold.
"""


_BIBLIOGRAPHY = r"""
# Bibliography — 33 PDFs by category (catalog 2026-05-26)

All PDFs at W:/03_文献・論文/00_電磁界解析/15_非接合/<subfolder>/.

## 00_overview (4)
- Non-conforming Finite Elements (Compumag Newsletter 2018-11)
- Arnold — Introduction to Mixed FE Methods for Elliptic Problems
- Arnold-Falk — New Mixed Formulation of Elasticity
- Lui 2001 — Monotone Schwarz Algorithms for Nonlinear Elliptic (M2AN)

## 01_mortar_nitsche (12)
- Nitsche-Mortar for Wave Eq (JTCA 2018)
- Hoppe-Wohlmuth 1996 — Local A Posteriori Estimators (M2AN)
- Becker-Hansbo-Stenberg 2003 — FE Method for DD with Non-Matching (M2AN)
- Brivadis-Buffa-Wohlmuth-Wunderlich 2014 — Isogeometric Mortar
- Brivadis-Buffa-Wohlmuth-Wunderlich 2015 — Quadrature for IGA Mortar
- Horger-Wohlmuth-Wunderlich 2016 — RB IGA Mortar Eigenvalue
- Horger-Reali-Wohlmuth-Wunderlich 2017 — Hybrid IGA Kirchhoff Plates
- Ludescher-Gross-Reusken 2018 — Multigrid for Unfitted Interface
- Dokken-Funke-Johansson-Schmidt 2018 — Multi-Mesh Nitsche Shape Opt
- Wang-Yang-Xie 2018 — Nitsche-XFEM Optimal Control on Interfaces
- Ager-Schott-Winter-Wall 2018 — Nitsche Cut-FE Flow + Poroelasticity
- Condensed-Constrained Mortar Preconditioner (arXiv 1912.05441, 2019)
- Chernov-Hansbo-Schetzke 2026 — Nonconforming hp-FE-BE with Nitsche

## 02_FETI_DP_BDDC (9)
- Mandel-Sousedík 2007 — BDDC and FETI-DP Minimalist Assumptions
- Šístek-Sousedík-Burda-Mandel 2011 — Parallel BDDC for Stokes
- Šístek-Březina-Sousedík 2015 — BDDC Mixed-Hybrid Porous Media
- Hofer-Langer 2015 — Dual-Primal IETI Multipatch IGA
- Hofer-Langer 2016 — Dual-Primal IETI Multipatch IGA with CG
- Hofer-Langer 2016 — DG-IETI-DP
- Kim-Chung-Wang 2016 — Adaptive Coarse Spaces for BDDC/FETI-DP 3D
- Lee-Park-Park 2020 — Dual Iterative Substructuring with Small Penalty
- Unified Non-Overlapping Robin-Schwarz with Cross Points (2022)

## 03_HFSS_mesh_fusion (2)
- Zhang-Liang 2015 — High-Order Curved Sliding-Mesh Spectral-Difference
- Zhang-Liang 2020 — Dynamic Transfinite Mortar Curved Sliding

## 04_DG_methods (4)
- Arnold-Brezzi-Cockburn-Marini 2002 — Unified Analysis DG Elliptic
- Hu-Yu 2005 — Signorini FEM-BEM Coupling (M2AN)
- Entropy-Stable hp Nonconforming DG SBP (arXiv 2017)
- Zhao-Chung 2019 — A Posteriori Error for Mortar Staggered DG

## 05_em_applications (13)
arXiv:
- Buffa-Maday-Rapetti 2001 — Sliding Mesh-Mortar 2D Eddy Electric Engines (M2AN)
- Alonso Rodríguez-Meddahi 2016 — DG FEM-BEM Time-Harmonic Eddy Current
- Alonso Rodríguez-Meddahi 2017 — DG for Time-Harmonic Eddy Current
- Bhat-DeFalco-Egger-Schöps-Wohlmuth 2019 — Isogeometric Mortar EM
- Egger-Harutyunyan-Loescher-Schöps-Steinbach 2020 — Harmonic Mortar Machines
- Egger-Harutyunyan-Merkel-Schöps 2021 — Torque by Harmonic Mortar

IEEE TMag:
- Rapetti-Bouillault-Santandrea-Buffa-Maday-Razek 2000 — Non-Match Edge Moving
- De Gersem-Hameyer 2002 — Multiconductor FE Eddy Current
- Rapetti-Maday-Bouillault-Razek 2002 — Eddy Current 3D Moving
- De Gersem-Clemens-Weiland 2003 — Iterative Hybrid FE-SE
- Antunes et al. 2005 — Hierarchic Interpolation Mortar Electrical Machines
- Koch-De Gersem-Weiland 2009 — Magnetostatic Hybrid FE-SE
- De Gersem 2010 — Combined SE-FE Magnetic Brake

## Pending acquisition (priority Top 5)
1. Bernardi-Maday-Patera 1994 — Mortar original (book chapter)
2. Farhat-Roux 1991 (IJNME) — FETI original
3. Dohrmann 2003 (SIAM JSC) — BDDC original
4. Lee-Volakis — HFSS Mesh Fusion academic precursor (citation unclear)
5. Ansys HFSS Mesh Fusion product white paper (Ansys customer portal)
"""


_TOPICS = {
    "_folder": _FOLDER_LAYOUT,
    "overview": _OVERVIEW,
    "mortar": _MORTAR,
    "nitsche": _NITSCHE,
    "feti_bddc": _FETI_BDDC,
    "feti_dp": _FETI_BDDC,
    "bddc": _FETI_BDDC,
    "dg": _DG,
    "hfss_mesh_fusion": _HFSS_MESH_FUSION,
    "hfss": _HFSS_MESH_FUSION,
    "mesh_fusion": _HFSS_MESH_FUSION,
    "ngsolve_recipe": _NGSOLVE_RECIPE,
    "ngsolve": _NGSOLVE_RECIPE,
    "recipe": _NGSOLVE_RECIPE,
    "em_applications": _EM_APPLICATIONS,
    "em": _EM_APPLICATIONS,
    "applications": _EM_APPLICATIONS,
    "lab_scenarios": _LAB_SCENARIOS,
    "lab": _LAB_SCENARIOS,
    "scenarios": _LAB_SCENARIOS,
    "bibliography": _BIBLIOGRAPHY,
    "bib": _BIBLIOGRAPHY,
    "papers": _BIBLIOGRAPHY,
}


def get_nonconforming_mesh_coupling_documentation(topic: str = "overview") -> str:
    """Return knowledge for non-conforming mesh coupling methods.

    Topics:
        "overview"          - When to use; decision matrix vs XFEM/PML/Kelvin
        "mortar"            - Lagrange-multiplier coupling (Bernardi-Maday-Patera)
        "nitsche"           - Penalty-based (Nitsche 1971, Becker-Hansbo-Stenberg)
        "feti_bddc"         - Iterative substructuring (Farhat, Dohrmann)
        "dg"                - Discontinuous Galerkin as mesh-coupling
        "hfss_mesh_fusion"  - Ansys product + academic analogs (Zhang-Liang)
        "ngsolve_recipe"    - Concrete NGSolve mortar/Nitsche code
        "em_applications"   - EM-specific (rotor-stator, harmonic, magnetic brake)
        "lab_scenarios"     - Sugahara Lab adoption priority + scenarios
        "bibliography"      - 33-PDF catalog organized by subfolder
        "all"               - Everything
    """
    t = (topic or "overview").strip().lower()
    if t == "all":
        parts = ["# Non-conforming mesh coupling for EM FEM", "", _FOLDER_LAYOUT]
        seen = set()
        for k in ("overview", "mortar", "nitsche", "feti_bddc", "dg",
                   "hfss_mesh_fusion", "ngsolve_recipe", "em_applications",
                   "lab_scenarios", "bibliography"):
            if id(_TOPICS[k]) not in seen:
                seen.add(id(_TOPICS[k]))
                parts.append(_TOPICS[k])
        return "\n\n".join(parts)
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic: {topic!r}. Use 'overview', 'mortar', 'nitsche', "
            f"'feti_bddc', 'dg', 'hfss_mesh_fusion', 'ngsolve_recipe', "
            f"'em_applications', 'lab_scenarios', 'bibliography', or 'all'.")
