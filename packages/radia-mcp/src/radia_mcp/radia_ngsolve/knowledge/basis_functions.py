"""
Finite-element basis function library for the Radia / NGSolve / radia.bem
ecosystem -- Mathematica-canonical reference for H1, HCurl, HDiv, L2 on
triangle / tetrahedron / quadrilateral / hexahedron / prism / pyramid,
arbitrary order.

Mathematica is the canonical notation; SymPy / NumPy / MATLAB are
generated from it.  Phase 1 (radia-mcp 0.38.0) covers the 2D triangle
and 3D tetrahedron primitives used by the in-tree BEM-A and scalar BEM
solvers.  Later phases extend to quadrilateral / hexahedron, full
HCurl coverage, and arbitrary-order hierarchical bases.

Read this when:
* Implementing a new BEM / FEM panel and need the canonical basis
  formula (and a verified NumPy translation).
* Cross-checking NGSolve's hierarchical basis output against a
  textbook Lagrange definition.
* Debugging RWG sign / divergence issues in BEM-A / RT₀ assembly.
* Generating a code-gen template (NumPy / MATLAB / SymPy / Maple)
  from the Mathematica canonical form.

The MCP server exposes this via basis_functions(topic=...). Topics:
theory_overview, hdiv_rt_bdm, code_gen_pattern,
triangle_h1_lagrange, triangle_h1_hierarchical,
triangle_rwg, triangle_l2,
tet_h1_lagrange, tet_rwg,
verification_recipes
"""

BASIS_FUNCTIONS = """\
# Finite-element basis functions — Mathematica-canonical reference (Phase 1)

This module records the explicit basis-function formulas used by
Radia / NGSolve / radia.bem in Mathematica notation, with NumPy /
SymPy translations and verification recipes.  The Mathematica form
is the canonical reference; everything else (radia.bem.efie_rwg's
``rwg_eval``, NGSolve's hierarchical FES, NumPy code-gen) is
expected to reproduce these formulas to machine precision.

Phase 1 scope (radia-mcp 0.38.0):
* Triangle H1 Lagrange P1, P2, P3, P4, P5
* Triangle H1 hierarchical (NGSolve compat) order 1..k
* Triangle HDiv RT₀ (= RWG, used by BEM-A)
* Triangle L2 P0, P1, P2, P3
* Tetrahedron H1 Lagrange P1, P2, P3
* Tetrahedron HDiv RT₀ (= 3D RWG, future BEM-A volumetric extension)

Topics: theory_overview, hdiv_rt_bdm, code_gen_pattern,
        triangle_h1_lagrange, triangle_h1_hierarchical,
        triangle_rwg, triangle_l2,
        tet_h1_lagrange, tet_rwg,
        verification_recipes

============================================================
## theory_overview — what each space means and where Radia uses it
============================================================

| Family    | Continuity     | Typical use in Radia                          |
|-----------|----------------|-----------------------------------------------|
| H1        | function value | scalar BEM (workpiece phi), FEM A-V (Omega)   |
| HCurl     | tangential     | NGSolve A-field (FEM A-V coil)                |
| HDiv      | normal         | BEM-A surface current (RWG = RT₀ on triangle) |
| L2        | none (DG)      | SurfaceL2 in Weggler EFIE saddle, divergence  |

**Key algebraic relations** (continuum):

* HCurl <-> HDiv via 90-degree rotation in 2D (`R*v` where R is the
  90-degree rotation matrix).
* `div(HDiv basis) ∈ L2` and `curl(HCurl basis) ∈ L2`; this is what
  makes RT_k / Ned_k the natural pair with L2_k for divergence /
  curl conforming discretisations.
* On a flat triangle, RT₀ has 3 basis functions (one per edge)
  with `div(f_e) = sigma_e / A_T` (constant).  See `triangle_rwg`.

**Reference triangle / tet conventions** used throughout:

* Reference triangle: vertices (0,0), (1,0), (0,1).  Barycentric
  coordinates λ₀ = 1−ξ−η, λ₁ = ξ, λ₂ = η.
* Reference tetrahedron: vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
  Barycentric λ₀ = 1−ξ−η−ζ, λ₁ = ξ, λ₂ = η, λ₃ = ζ.
* Edge ordering: (0,1), (1,2), (2,0) for triangle; (0,1), (0,2),
  (0,3), (1,2), (1,3), (2,3) for tetrahedron.

============================================================
## hdiv_rt_bdm — executable Mathematica study of both HDiv families
============================================================

The current NGSolve convention is explicit:

```python
HDiv(mesh, order=p)           # BDM on simplices
HDiv(mesh, order=p, RT=True)  # Raviart--Thomas on simplices
```

The corresponding Mathematica study API lives in
`mathematica/basis_functions/simplex_ho.wls`:

```mathematica
HDivTetBDM[p]   (* [P_p]^3, p >= 1 *)
HDivTetRT[p]    (* [P_p]^3 + x Ptilde_p *)
HDivTrigBDM[p]  (* [P_p]^2, p >= 1 *)
HDivTrigRT[p]   (* [P_p]^2 + x Ptilde_p *)

HDivTetFamilyLedger[p]
HDivTrigFamilyLedger[p]
```

For equal simplex order `p`,

* `div(BDM_p) = P_(p-1)`;
* `div(RT_p) = P_p`;
* `BDM_p` is a strict subspace of `RT_p`;
* both families have the same `ker(div)` dimension and the same degree-`p`
  normal-trace space.

The first three tetrahedral orders are:

| p | BDM DoF | RT DoF | BDM div DoF | RT div DoF | shared ker(div) |
|---:|--------:|-------:|------------:|-----------:|----------------:|
| 1 | 12 | 15 | 1 | 4 | 11 |
| 2 | 30 | 36 | 4 | 10 | 26 |
| 3 | 60 | 70 | 10 | 20 | 50 |
| 6 | 252 | 280 | 56 | 84 | 196 |

For triangles the corresponding rows are `6/8`, `12/15`, and `20/24` total
DoFs.  RT's equal-order increment is entirely in charge-carrying divergence
modes; it does not add loop modes.  This is the structural reason BDM is the
production default when its lower-order charge space is accurate enough.

Hex/Quad are different.  NGSolve uses one tensor-product H(div) space for both
settings of the `RT` selector.  `hdiv.wls` therefore exposes intentional aliases:

```mathematica
HDivHexBDM[p] === HDivHexRT[p]
HDivQuadBDM[p] === HDivQuadRT[p]
```

Both `.wls` files self-test dimensions, polynomial independence, normal-trace
rank, divergence rank, de Rham inclusion, and the shared-kernel identity.  Run:

```powershell
wolframscript -file simplex_ho.wls
wolframscript -file hdiv.wls
```

============================================================
## code_gen_pattern — Mathematica → NumPy / SymPy / MATLAB
============================================================

The canonical form is Mathematica.  Translation rules (small enough
to be rule-driven, no separate library needed):

| Mathematica           | NumPy                  | SymPy            |
|-----------------------|------------------------|------------------|
| `xi`, `eta`           | `xi`, `eta` (arrays)   | `Symbol('xi')`   |
| `1 - xi - eta`        | `1 - xi - eta`         | (literal)        |
| `(xi + eta) (1-eta)`  | `(xi + eta) * (1-eta)` | (literal)        |
| `Cross[u, v]`         | `np.cross(u, v)`       | `u.cross(v)`     |
| `D[f, xi]`            | `np.gradient(f, ...)`  | `diff(f, xi)`    |
| `Piecewise[...]`      | `np.where(...)`        | `Piecewise(...)` |
| `Simplify[...]`       | (n/a)                  | `simplify(...)`  |
| `Integrate[f, ...]`   | (numeric quad)         | `integrate(...)` |

**Pattern**: Mathematica functions are written in `Module[...]` with
local variables, take reference coords as arguments, return values or
arrays.  The `packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m` package has the full
Phase 1 canonical implementation; SymPy translations live in
``tests/basis/test_basis_functions.py`` (used as the verification
oracle); NumPy translations live alongside the radia code that
consumes them (e.g. ``radia.bem.efie_rwg.rwg_eval`` for RT₀).

**Verification recipe**: every Mathematica definition has a SymPy
equivalent that asserts:
  1. Polynomial completeness (P_k completeness on triangle = (k+1)(k+2)/2 dofs)
  2. Partition of unity (Σ basis_i = 1 for H1 Lagrange)
  3. Divergence identity (Σ_e ∮_∂T (f_e · n) dl = ∫_T div(f_e) dA for HDiv)
  4. Numerical match against radia's compiled implementation to 1e-12
     at random reference points.

============================================================
## triangle_h1_lagrange — Triangle H1 Lagrange P1..P5
============================================================

**Reference**: barycentric (λ₀, λ₁, λ₂), λ₀ = 1−ξ−η, λ₁ = ξ, λ₂ = η.

### P1 (3 dofs, vertex-only)

```mathematica
(* Mathematica: Triangle H1 Lagrange order 1 (3 nodal basis functions) *)
TriH1P1[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta;
  lam1 = xi;
  lam2 = eta;
  {lam0, lam1, lam2}                  (* phi[v0], phi[v1], phi[v2] *)
]
```

NumPy (vectorised):

```python
def tri_h1_p1(xi, eta):
    \"\"\"Triangle Lagrange P1, returns (3, n) for (xi, eta) of shape (n,).\"\"\"
    import numpy as np
    lam0 = 1.0 - xi - eta
    return np.stack([lam0, xi, eta], axis=0)
```

### P2 (6 dofs: 3 vertex + 3 edge midpoints)

Nodal layout (Lagrange P2):
```
     v2 (0,1)
      *
      |\\
      | \\
   m20*  *m12
      |   \\
      |    \\
      *--*--*
   v0  m01  v1
```

```mathematica
(* Mathematica: Triangle H1 Lagrange order 2 (6 basis functions).
   Node order: v0, v1, v2, m01, m12, m20.                          *)
TriH1P2[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
  {
    lam0 (2 lam0 - 1),                (* phi at v0   *)
    lam1 (2 lam1 - 1),                (* phi at v1   *)
    lam2 (2 lam2 - 1),                (* phi at v2   *)
    4 lam0 lam1,                      (* phi at m01  *)
    4 lam1 lam2,                      (* phi at m12  *)
    4 lam2 lam0                       (* phi at m20  *)
  }
]
```

### P3 (10 dofs: 3 vertex + 6 edge (1/3, 2/3) + 1 centroid)

```mathematica
(* Mathematica: Triangle H1 Lagrange order 3 (10 basis functions).
   Node order: v0, v1, v2, e01a, e01b, e12a, e12b, e20a, e20b, c.
   Edge nodes: a = 1/3 from first vertex, b = 2/3 from first vertex. *)
TriH1P3[xi_, eta_] := Module[{lam0, lam1, lam2},
  lam0 = 1 - xi - eta; lam1 = xi; lam2 = eta;
  {
    (1/2) lam0 (3 lam0 - 1) (3 lam0 - 2),    (* v0  *)
    (1/2) lam1 (3 lam1 - 1) (3 lam1 - 2),    (* v1  *)
    (1/2) lam2 (3 lam2 - 1) (3 lam2 - 2),    (* v2  *)
    (9/2) lam0 lam1 (3 lam0 - 1),             (* e01a, near v0 *)
    (9/2) lam0 lam1 (3 lam1 - 1),             (* e01b, near v1 *)
    (9/2) lam1 lam2 (3 lam1 - 1),             (* e12a, near v1 *)
    (9/2) lam1 lam2 (3 lam2 - 1),             (* e12b, near v2 *)
    (9/2) lam2 lam0 (3 lam2 - 1),             (* e20a, near v2 *)
    (9/2) lam2 lam0 (3 lam0 - 1),             (* e20b, near v0 *)
    27 lam0 lam1 lam2                         (* centroid       *)
  }
]
```

### P4, P5 — see RadiaBasis.m

P4 (15 dofs) and P5 (21 dofs) follow the same hierarchical-Lagrange
pattern with k+1 nodes per edge and (k-1)(k-2)/2 interior nodes.
The closed-form expressions are too long to inline here; see
``packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m`` for the full definitions and
``tests/basis/test_basis_functions.py`` for the SymPy-verified
NumPy ports.

**Verification (Phase 1 minimum, should hold for all k)**:

```python
# At any (xi, eta) inside the reference triangle:
phi = tri_h1_pk(xi, eta)
assert abs(phi.sum() - 1.0) < 1e-12   # partition of unity
assert phi.shape[0] == (k+1)*(k+2)//2 # dof count
# At each Lagrange node, exactly one basis function = 1, others = 0:
for j, (xn, yn) in enumerate(lagrange_nodes(k)):
    phi = tri_h1_pk(xn, yn)
    assert abs(phi[j] - 1.0) < 1e-12
    others = np.delete(phi, j)
    assert np.all(np.abs(others) < 1e-12)
```

============================================================
## triangle_h1_hierarchical — Triangle H1 (NGSolve hierarchical)
============================================================

NGSolve's `H1(mesh, order=k)` uses a **hierarchical** basis: the
first 3 functions are the linear vertex shapes (= P1), then edge
bubbles built from integrated Legendre polynomials, then interior
bubbles from products.  This is NOT identical to Lagrange P_k —
the basis values at non-vertex nodes are different — but the
spanned space is the same (P_k complete).

```mathematica
(* Integrated Legendre polynomial of degree k:
     L̃_k(x) = (1/(k(k-1))) (P_k(x) - P_{k-2}(x))   for k >= 2
   On the reference edge x ∈ [0, 1] this gives a basis that
   vanishes at both endpoints, useful as edge bubble.            *)
LegendreP[k_, x_] := LegendreP[k, x];          (* Mathematica builtin *)
IntLegendre[k_, x_] := If[k < 2, 0,
  (LegendreP[k, 2 x - 1] - LegendreP[k - 2, 2 x - 1]) /
  (2 (2 k - 1))];

(* Edge bubble e_{ij,p} on edge (i,j), order p, p >= 2.
   Uses lam_i, lam_j for that edge.  Vanishes at both endpoints.   *)
EdgeBubble[lami_, lamj_, p_] := (lami + lamj)^p *
  IntLegendre[p, lamj / (lami + lamj)];

(* Interior bubble (cell bubble), order p, valid for p >= 3 *)
InteriorBubble[xi_, eta_, p1_, p2_] := lam0 lam1 lam2 *
  LegendreP[p1, 2 lam1 - 1] * LegendreP[p2, 2 lam2 - 1] /. {
    lam0 -> 1 - xi - eta, lam1 -> xi, lam2 -> eta};
```

The full hierarchical basis on a triangle of order k has:
* 3 vertex shapes (P1)
* 3 (k-1) edge bubbles (one per edge, p = 2..k)
* (k-1)(k-2)/2 interior bubbles

Total = (k+1)(k+2)/2 = same as P_k Lagrange.

**Match with NGSolve**: the integrated-Legendre normalisation here
matches NGSolve's `H1HighOrderFE` shape functions (per
`ngsolve/fem/h1ho_fe_impl.hpp`).  Verification: project a known P_k
polynomial via NGSolve's `Set(...)` then sample the GridFunction at
random reference points; the Mathematica + SymPy port reproduces
the same nodal values to 1e-13.

============================================================
## triangle_rwg — Triangle HDiv RT₀ (= RWG, BEM-A coil)
============================================================

The Rao-Wilton-Glisson basis (= Raviart-Thomas RT₀ on a flat
triangle) is **the** workhorse for surface-current BEM in Radia.
Three basis functions per triangle, one per edge:

```
        v2
         *
        /\\
       / |\\
      /  | \\
   f₂ ↘  |  ↙ f₁
    /    |    \\
   /     |     \\
  *------*------*
  v0    f₀     v1
```

**Definition** (NGSolve / Phase 1.4 convention, RT₀-normalised):

For each local edge k of triangle T with vertices (v₀, v₁, v₂):
```
  edge 0 = (v₁, v₂),  v_opp = v₀
  edge 1 = (v₂, v₀),  v_opp = v₁
  edge 2 = (v₀, v₁),  v_opp = v₂
```

The RT₀ basis function for edge k is:
```
  f_k(r) = sigma_k * (r - v_opp_k) / (2 A_T)
```
where `A_T` is the triangle area and `sigma_k ∈ {+1, -1}` is the
GLOBAL orientation sign (chosen consistent across the two adjacent
triangles so the normal-component is continuous across the edge).

```mathematica
(* Mathematica: Triangle RT₀ (RWG), 3 basis functions.
   Inputs: ref coords (xi, eta), vertex array v={{x0,y0,z0},...},
   sigma = {sigma0, sigma1, sigma2} (global orientation signs).    *)
TriRT0[xi_, eta_, v_List, sigma_List] := Module[
  {r, vopp, area},
  r = (1 - xi - eta) v[[1]] + xi v[[2]] + eta v[[3]];
  vopp = {v[[1]], v[[2]], v[[3]]};
  area = (1/2) Norm[Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]]];
  Table[ sigma[[k]] (r - vopp[[k]]) / (2 area), {k, 1, 3}]
]
```

**Divergence identity**: `div(f_k) = sigma_k / A_T` (constant on
each triangle).

**Sign convention WARNING**: the GLOBAL signs `sigma_k` come from
the orientation of the global edge (`v_min -> v_max` index
order) NOT from the local CCW order of the triangle's vertices.
For canonical SS Sauter-Schwab Galerkin pairs this matters at
common-edge cases — see `radia.bem.efie_rwg.rwg_eval` and the
2026-04-30 fix `feedback_dl_winding_bug_2026_04_30.md` for the
full story (a non-cyclic corner permutation flips the local
sign and gave SL negative eigenvalues until the override
parameter was added).

**NumPy reference**: `radia.bem.efie_rwg.rwg_eval` (line 154+) is
the production implementation; its docstring says it is
RT₀-normalised "(NO L_e factor)".

**Verification**:
1. `int_T div(f_k) dA = sigma_k` (closed-form check).
2. Trace test: `f_k · n_edge_k = sigma_k / L_edge_k` on the
   associated edge, zero on the other two edges.
3. Numerical match: `radia.bem.efie_rwg.rwg_eval(verts, tris[t],
   xi, eta)` reproduces this Mathematica formula to 1e-13.

============================================================
## triangle_l2 — Triangle L2 P0..P3 (BEM-A SurfaceL2, jump terms)
============================================================

Discontinuous (DG) basis on the triangle.  Used in BEM-A as the
SurfaceL2 multiplier in Weggler's EFIE saddle-point system, and
in scalar BEM as the ground for divergence terms.

### P0 — constant per triangle (1 dof)

```mathematica
(* L2 P0: single constant function = 1 on the reference triangle. *)
TriL2P0[xi_, eta_] := {1};
```

NumPy: `np.ones_like(xi).reshape(1, *xi.shape)`

### P1 — affine per triangle (3 dofs)

```mathematica
(* L2 P1 hierarchical: 1, x-bar, y-bar, where bar = re-centred to centroid. *)
TriL2P1[xi_, eta_] := {
  1,
  xi - 1/3,
  eta - 1/3
};
```

### P2, P3 — orthogonal Dubiner basis

For higher orders, the Dubiner / Koornwinder basis is the natural
DG choice on the simplex (orthogonal under the L2 inner product
on the reference triangle):

```mathematica
(* Dubiner basis on triangle, order k.  Indexed by (i, j) with
   i+j <= k.  Uses Jacobi polynomials J_n^{a,b}.                  *)
DubinerTri[i_, j_, xi_, eta_] := Module[{a, b, P1, P2},
  a = 2 xi / (1 - eta) - 1;     (* collapsed coords *)
  b = 2 eta - 1;
  P1 = JacobiP[i, 0, 0, a];      (* Legendre on a *)
  P2 = JacobiP[j, 2 i + 1, 0, b];
  P1 P2 (1 - eta)^i
];
```

**Orthogonality**: `∫_T DubinerTri[i,j] DubinerTri[m,n] dA = 0` if
(i,j) ≠ (m,n).  Norm closed-form involves Gamma functions; see
RadiaBasis.m.

**NGSolve match**: `L2(mesh, order=k)` uses a different
(non-orthogonal) hierarchical basis but spans the same space.
Verify by L2 projection of a known polynomial.

============================================================
## tet_h1_lagrange — Tetrahedron H1 Lagrange P1..P3
============================================================

**Reference**: λ₀ = 1−ξ−η−ζ, λ₁ = ξ, λ₂ = η, λ₃ = ζ.

### P1 (4 dofs, vertex-only)

```mathematica
TetH1P1[xi_, eta_, zeta_] := Module[{lam0, lam1, lam2, lam3},
  lam0 = 1 - xi - eta - zeta; lam1 = xi; lam2 = eta; lam3 = zeta;
  {lam0, lam1, lam2, lam3}
]
```

### P2 (10 dofs: 4 vertex + 6 edge midpoints)

Edge ordering: e01, e02, e03, e12, e13, e23.

```mathematica
TetH1P2[xi_, eta_, zeta_] := Module[{lam0, lam1, lam2, lam3},
  lam0 = 1 - xi - eta - zeta; lam1 = xi; lam2 = eta; lam3 = zeta;
  {
    lam0 (2 lam0 - 1),  lam1 (2 lam1 - 1),  lam2 (2 lam2 - 1),  lam3 (2 lam3 - 1),
    4 lam0 lam1,  4 lam0 lam2,  4 lam0 lam3,
    4 lam1 lam2,  4 lam1 lam3,  4 lam2 lam3
  }
]
```

### P3 (20 dofs: 4 vertex + 12 edge (1/3, 2/3) + 4 face centroid)

Closed form is long; see `packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/RadiaBasis.m`.

**Verification**: same as triangle (partition of unity, dof count
(k+1)(k+2)(k+3)/6 = 4 / 10 / 20 / 35 / ...).

============================================================
## tet_rwg — Tetrahedron HDiv RT₀ (3D RWG, future BEM-A volumetric)
============================================================

3D analogue of the surface RWG: 4 face-based basis functions, one
per face of the tetrahedron.  Used in 3D BEM-A volumetric
formulations (NOT in Phase 1 of the radia BEM-A; reserved for
future phases when the workpiece is meshed volumetrically and we
want HDiv basis on the tet faces).

```mathematica
(* Mathematica: Tetrahedron RT₀, 4 basis functions, one per face. *)
TetRT0[xi_, eta_, zeta_, v_List, sigma_List] := Module[
  {r, vopp, vol},
  r = (1 - xi - eta - zeta) v[[1]] + xi v[[2]] + eta v[[3]] + zeta v[[4]];
  vopp = v;        (* face k is opposite vertex k *)
  vol = (1/6) Abs[Det[{v[[2]] - v[[1]], v[[3]] - v[[1]], v[[4]] - v[[1]]}]];
  Table[ sigma[[k]] (r - vopp[[k]]) / (3 vol), {k, 1, 4}]
]
```

**Divergence identity**: `div(f_k) = sigma_k / V_T` (constant per
tetrahedron).

============================================================
## verification_recipes — checks that should hold for every basis
============================================================

The Mathematica definitions and their Python ports must satisfy
the following identities to machine precision (1e-12 ish on
random reference points):

1. **Partition of unity** (H1 Lagrange):
   `Σ_i phi_i(xi, eta) == 1` for every (xi, eta) inside the
   reference element.

2. **Kronecker delta at nodal points** (H1 Lagrange):
   `phi_i(node_j) == delta_ij` for all i, j.

3. **DOF count** (H1 P_k on triangle / tet):
   `(k+1)(k+2)/2` for triangle, `(k+1)(k+2)(k+3)/6` for tet.

4. **Polynomial completeness** (any P_k space):
   project `f(xi, eta) = xi^a * eta^b` for all a+b<=k via the basis;
   reconstructed value at any test point matches f exactly.

5. **Divergence identity** (HDiv RT₀):
   `div(f_e) == sigma_e / A_T` (constant); verifiable both
   symbolically (SymPy `diff`) and numerically.

6. **Trace identity** (HDiv RT₀):
   `f_e · n_edge_e == sigma_e / L_edge_e` on the associated edge,
   zero on the other two.

7. **Orthogonality** (L2 Dubiner):
   `∫_T phi_i phi_j dA = 0` for i ≠ j.  Symbolic via
   SymPy `integrate` over the reference triangle.

8. **Numerical match against radia code**: the Phase 1 SymPy
   ports must reproduce `radia.bem.efie_rwg.rwg_eval` to 1e-13
   on 100 random test triangles + 100 random reference points.

These checks are codified in `tests/basis/test_basis_functions.py`
(part of the Phase 1 release).  CI runs them on every push;
runtime ~5 sec for the full Phase 1 matrix.
"""


_TOPICS = (
    "theory_overview", "hdiv_rt_bdm", "code_gen_pattern",
    "triangle_h1_lagrange", "triangle_h1_hierarchical",
    "triangle_rwg", "triangle_l2",
    "tet_h1_lagrange", "tet_rwg",
    "verification_recipes",
)


def get_basis_functions_documentation(topic: str = "") -> str:
    """Return the basis-function knowledge.

    Args:
        topic: Empty for the full document, or one of the entries in
               ``_TOPICS`` above for a single section.

    Returns:
        Markdown string.  When ``topic`` is empty the entire document
        is returned (~12 kB); otherwise the requested section is sliced
        between its `## topic ` header and the next topic header.
    """
    if not topic:
        return BASIS_FUNCTIONS

    topic = topic.strip()
    if topic not in _TOPICS:
        return (f"Unknown topic {topic!r}.  Valid topics: "
                + ", ".join(_TOPICS))

    # Locate every "\n## <topic> " header offset.
    headers = []
    pos = 0
    while True:
        nxt = BASIS_FUNCTIONS.find("\n## ", pos)
        if nxt < 0:
            break
        line_end = BASIS_FUNCTIONS.find("\n", nxt + 1)
        line = BASIS_FUNCTIONS[nxt + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, nxt + 1))
                break
        pos = nxt + 1

    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."

    section_start = req_starts[0]
    sep_above = BASIS_FUNCTIONS.rfind(
        "============================================================",
        0, section_start)
    if sep_above >= 0 and sep_above > section_start - 80:
        section_start = sep_above

    last_req = req_starts[-1]
    section_end = len(BASIS_FUNCTIONS)
    for kw, off in headers:
        if kw != topic and off > last_req:
            sep = BASIS_FUNCTIONS.rfind(
                "============================================================",
                0, off)
            section_end = sep if sep > 0 else off
            break

    return BASIS_FUNCTIONS[section_start:section_end].rstrip() + "\n"
