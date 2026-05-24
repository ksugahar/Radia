# Notes: H(curl) hierarchical bases on tetrahedra — `nograds` and tree-cotree

Reference material for Mathematica-side construction & verification of
NGSolve-style hierarchical H(curl) basis functions on tetrahedra.
Used by the `mathematica_*` MCP tools when the user asks about FEM basis
functions, gauge fixing, or NGSolve element details.

Primary references:
- Zaglmayr, "High Order Finite Element Methods for Electromagnetic Field
  Computation" (PhD thesis, JKU Linz 2006), Sec. 5–6.
- Schöberl & Zaglmayr, "High order Nédélec elements with local complete
  sequence properties," COMPEL 24 (2005).
- Gross & Kotiuga, "Electromagnetic Theory and Computation: A Topological
  Approach," CUP 2004.
- Pellikka, Suuriniemi, Kettunen, Geuzaine, "Homology and cohomology
  computation in finite element modeling," SIAM J. Sci. Comp. 35 (2013).

---

## 1. The hierarchical H(curl) split (Zaglmayr)

On a tetrahedron at polynomial order p, the H(curl) shape functions split
exactly as

    Φ_HCurl^p  =  { ∇φ : φ ∈ Φ_H1^p, φ not vertex-type }   ⊕   Φ_rotational^p

The **gradient block** is in ker(curl); the **rotational block** has
non-zero curl. The basis is constructed so this split is *explicit at the
shape-function level* (you can identify each DoF as one or the other).

DoF inventory on a single tetrahedron at order p (tetrahedral H(curl)):

| Block                  | Count                 | Subset of ker(curl)? |
|------------------------|-----------------------|----------------------|
| Edge lowest-order (Nédélec/Whitney) `λᵢ∇λⱼ−λⱼ∇λᵢ` | 6           | No (lowest order)    |
| Edge gradients         `∇ L_k^s(λⱼ−λᵢ, λᵢ+λⱼ)`, k=2..p | 6·(p−1)      | **Yes**              |
| Edge "type II" rotational                              | 6·(p−1)      | No                   |
| Face gradients         `∇ (H¹ face shape)`             | 4·C(p−1,2)   | **Yes**              |
| Face type II                                            | 4·(some)     | No                   |
| Cell gradients         `∇ (H¹ cell shape)`             | C(p−1,3)     | **Yes**              |
| Cell type II                                            | (some)       | No                   |

---

## 2. `nograds=True` (NGSolve option)

```python
fes = HCurl(mesh, order=p, nograds=True)
```

**Meaning:** drop all *higher-order* gradient-block shape functions
(`Edge gradients`, `Face gradients`, `Cell gradients` rows above) from
the basis. The lowest-order Nédélec/Whitney functions and **all** type II
(rotational) functions are kept.

**Why it helps:**
1. **Condition number / solver convergence** — for `(ν curl·, curl·)`
   bilinear forms, gradients map to zero. Including them gives the
   stiffness matrix a large zero eigenspace; removing them at the basis
   level avoids this.
2. **DoF count** — at p=3 tet roughly 30% of DoFs are higher-order
   gradients; at p=5 closer to 50%.
3. **Physical interpretation** — only DoFs that change `B = curl A` are
   retained; the rest is gauge.

**What `nograds` does *not* fix:**
- The **lowest-order** gradient kernel (∇ of H¹ vertex basis on the
  Whitney complex) is still present. The basis must include these to
  preserve the discrete de Rham complex's exactness at p=1.
- Multiply connected domain (β₁ > 0) cohomology generators are unaffected.

**When NOT to use `nograds`:**
- Time-domain `∂²A/∂t² + curl curl A` if gradient components contribute
  to the mass term in a way the formulation depends on.
- Helmholtz / resonant-mode analysis where the gradient eigenspace at
  ω = 0 is part of the physical question.
- Mixed formulations using a Lagrange multiplier on the gradient space.
- NGSolve internals: `gradientrange`, `Compress`, and a few projection
  utilities expect the full basis; with `nograds=True` you may need a
  separate `HCurl_full` space alongside.

---

## 3. Tree-cotree decomposition

Classical (lowest-order, Whitney 1-forms):
- Build a spanning tree T of the edge graph of the mesh.
- **Tree edges** → set DoF = 0 (these span exactly the gradient kernel of
  the lowest-order Nédélec space on a simply connected domain).
- **Cotree edges** → keep DoF.

This is the lowest-order analogue of `nograds`. The two are
complementary:

| Strategy        | Where it acts                         |
|-----------------|---------------------------------------|
| `nograds=True`  | Higher-order gradient DoFs (basis-level removal) |
| Tree-cotree     | Lowest-order gradient DoFs (combinatorial gauge fix) |
| Cohomology cuts | Topological null space when β₁ > 0    |

Together they give a complete gauge fix on a simply connected domain at
all polynomial orders.

### Generalizations of tree-cotree

1. **High order:** in hierarchical bases, the gradient block is explicit,
   so high-order tree-cotree reduces to "zero out all higher-order
   gradient DoFs" — exactly what `nograds=True` already does.
2. **Topology:** for multiply connected domains, augment tree-cotree
   with β₁ cohomology generators (thick cuts). Algorithms: Kotiuga,
   Pellikka–Suuriniemi–Kettunen–Geuzaine 2013, Hiptmair–Ostrowski.
3. **Other complexes:** the same idea applies along the de Rham
   complex H¹ → H(curl) → H(div) → L². A **face-tree** on the
   dual graph fixes the divergence-free gauge for H(div).
4. **FEEC viewpoint (Arnold–Falk–Winther):** tree-cotree = choosing a
   combinatorial section of the exterior derivative `d`.

---

## 4. Verifying in Mathematica

A working H¹ p=2 demo lives at `C:\temp\ngsolve_tet_h1_demo.wls`
(generated from this conversation). Extension recipes for H(curl):

- **Nédélec lowest order:** `N_e = λᵢ ∇λⱼ − λⱼ ∇λᵢ` for each edge e=(i,j).
- **Edge gradient at order k:** `∇ Lob_k^scaled(λⱼ−λᵢ, λᵢ+λⱼ)` —
  these are *exactly* the DoFs `nograds=True` drops.
- **Edge type II at order k:** Zaglmayr Sec 5.4 formula
  `Lob_k^scaled(λⱼ−λᵢ, λᵢ+λⱼ) · (λᵢ ∇λⱼ − λⱼ ∇λᵢ)`.
- **Verification:** `Curl[edge_gradient_shape] == 0` should hold
  symbolically; `Curl[edge_type_II_shape]` is nonzero.

Sanity tests to run:
1. `Curl` of every gradient-type shape function is zero.
2. Mass matrix of full HCurl basis vs `nograds`-pruned basis:
   the smaller matrix should be the principal submatrix.
3. Stiffness matrix (`(curl·, curl·)`) has a much smaller kernel
   under `nograds=True` (only the lowest-order ker(curl) remains).

---

## 5. Quick decision tree

User asks about NGSolve HCurl performance / DoF count / gauge:
- If solving magnetostatics or eddy-current: recommend `nograds=True`
  + tree-cotree (or AMS preconditioner) + cohomology cuts if topology
  demands it.
- If running eigenvalue analysis or time-domain wave: think first.
- If asked to demonstrate symbolically: extend the
  `ngsolve_tet_h1_demo.wls` template with `N_e = λᵢ∇λⱼ − λⱼ∇λᵢ` and
  show `Curl` of gradient vs type II.
