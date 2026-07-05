"""NGSolve hierarchical H(curl) bases — Zaglmayr / nograds / tree-cotree.

Migrated 2026-05-23 from public-safe curated corpus
notes_fem_hcurl.md (per "all FEM knowledge lives in radia-mcp" policy).
"""

ZAGLMAYR_SPLIT = r"""
# Hierarchical H(curl) bases on tetrahedra (Zaglmayr/Schöberl, NGSolve)

Reference: Sabine Zaglmayr, "High Order Finite Element Methods for
Electromagnetic Field Computation" (PhD thesis, JKU Linz 2006), Sec. 5-6.
Schöberl-Zaglmayr, "High order Nédélec elements with local complete
sequence properties", COMPEL 24 (2005).

## The hierarchical H(curl) split

On a tetrahedron at polynomial order p, the H(curl) shape functions
split exactly as

    Φ_HCurl^p  =  { ∇φ : φ ∈ Φ_H1^p, φ not vertex-type }   ⊕   Φ_rotational^p

The **gradient block** is in ker(curl); the **rotational block** has
non-zero curl. The basis is constructed so this split is *explicit at
the shape-function level* (you can identify each DoF as one or the other).

## DoF inventory on a single tetrahedron at order p

| Block | Count | Subset of ker(curl)? |
|-------|-------|----------------------|
| Edge lowest-order (Nédélec/Whitney) `λ_i ∇λ_j − λ_j ∇λ_i` | 6 | No (lowest order) |
| Edge gradients `∇ L_k^s(λ_j-λ_i, λ_i+λ_j)`, k=2..p | 6·(p-1) | **Yes** |
| Edge "type II" rotational | 6·(p-1) | No |
| Face gradients `∇ (H¹ face shape)` | 4·C(p-1,2) | **Yes** |
| Face type II | 4·(some) | No |
| Cell gradients `∇ (H¹ cell shape)` | C(p-1,3) | **Yes** |
| Cell type II | (some) | No |

## Why the split matters

For Maxwell `(ν curl·, curl·)` bilinear forms, gradient-type DoFs map
to **zero**. The matrix has a huge kernel that hurts iterative solver
convergence. The hierarchical split lets you either:

1. **Use a preconditioner that handles the kernel** (Hiptmair-Xu AMS,
   see `radia_mcp.matrix_solvers.preconditioners.ams_hiptmair_xu`),
   OR
2. **Remove gradient-type DoFs at the basis level** (NGSolve `nograds=True`).
"""


NOGRADS_OPTION = r"""
# `nograds=True` option (NGSolve hierarchical H(curl))

```python
fes = HCurl(mesh, order=p, nograds=True)
```

## What it does

Removes the **higher-order gradient DoFs** (`Edge gradients`,
`Face gradients`, `Cell gradients` rows in the Zaglmayr inventory) from
the basis. The lowest-order Nédélec/Whitney functions and **all** type
II (rotational) functions are kept.

## Why it helps

1. **Condition number / solver convergence**: for `(ν curl·, curl·)`
   bilinear forms, gradients map to zero. Including them gives the
   stiffness matrix a large zero eigenspace; removing them at the basis
   level avoids this.
2. **DoF count**: at p=3 tet, roughly 30% of DoFs are higher-order
   gradients; at p=5 closer to 50%.
3. **Physical interpretation**: only DoFs that change `B = curl A` are
   retained; the rest is gauge.

## What `nograds` does NOT fix

- The **lowest-order gradient kernel** (∇ of H¹ vertex basis on the
  Whitney complex) is still present. The basis must include these to
  preserve the discrete de Rham complex's exactness at p=1.
- **Multiply-connected domain** (β₁ > 0) cohomology generators are
  unaffected — need topological cuts.

## When NOT to use `nograds`

- **Time-domain** `∂²A/∂t² + curl curl A` if gradient components
  contribute to the mass term in a way the formulation depends on.
- **Helmholtz / resonant-mode** analysis where the gradient eigenspace
  at ω = 0 is part of the physical question.
- **Mixed formulations** using a Lagrange multiplier on the gradient
  space.
- **NGSolve internals**: `gradientrange`, `Compress`, and a few
  projection utilities expect the full basis; with `nograds=True` you
  may need a separate `HCurl_full` space alongside.

## Lab default

The lab uses **ungauged + CompactAMS** (which handles the kernel
naturally) more often than `nograds=True`. For static magnetostatic
with `(curl·, curl·)` only, `nograds=True` is a complementary
optimization that reduces DoF count without changing physics.
"""


TREE_COTREE = r"""
# Tree-cotree decomposition (classical Whitney 1-forms gauge)

## Classical lowest-order

- Build a spanning tree T of the edge graph of the mesh.
- **Tree edges** → set DoF = 0 (these span exactly the gradient kernel
  of the lowest-order Nédélec space on a simply-connected domain).
- **Cotree edges** → keep DoF.

This is the lowest-order analogue of `nograds`. The two strategies are
**complementary**:

| Strategy | Where it acts |
|----------|---------------|
| `nograds=True` | Higher-order gradient DoFs (basis-level removal) |
| Tree-cotree | Lowest-order gradient DoFs (combinatorial gauge fix) |
| Cohomology cuts | Topological null space when β₁ > 0 |

Together they give a complete gauge fix on a simply-connected domain
at all polynomial orders.

## Generalizations of tree-cotree

1. **High order**: in hierarchical bases, the gradient block is
   explicit, so high-order tree-cotree reduces to "zero out all
   higher-order gradient DoFs" — exactly what `nograds=True` already
   does.

2. **Topology**: for multiply-connected domains, augment tree-cotree
   with β₁ cohomology generators (thick cuts). Algorithms:
   - Kotiuga 1989
   - Pellikka-Suuriniemi-Kettunen-Geuzaine 2013, SIAM J. Sci. Comp. 35
   - Hiptmair-Ostrowski

3. **Other complexes**: the same idea applies along the de Rham complex
   H¹ → H(curl) → H(div) → L². A **face-tree** on the dual graph fixes
   the divergence-free gauge for H(div).

4. **FEEC viewpoint** (Arnold-Falk-Winther): tree-cotree = choosing a
   combinatorial section of the exterior derivative `d`.

## Multiply-connected domain (β₁ > 0)

A torus, a body with a hole, or any topologically non-trivial conductor.
Tree-cotree gauge MISSES β₁ harmonic forms (loop currents). These must
be added explicitly via:
- Cohomology generators (Pellikka et al. 2013)
- "Thick cuts" through the body
- T-Ω cuts (in T-Ω formulation)

## Cross-references

- `radia_mcp.matrix_solvers.em_specific.tree_cotree` — solver-side view
- `radia_mcp.bem.surface_ie.pmchwt` — BEM analogue for multiply-connected
- `radia_mcp.differential_forms.cohomology` — exterior calculus foundation
"""


MATHEMATICA_VERIFICATION = r"""
# Mathematica verification recipes for HCurl shape functions

## H¹ p=2 demo

A working H¹ p=2 demo lives at `C:\temp\ngsolve_tet_h1_demo.wls`
(generated from earlier Mathematica session work). Use as a starting
template.

## HCurl extension recipes

- **Nédélec lowest order**:
  `N_e = λ_i ∇λ_j − λ_j ∇λ_i` for each edge e=(i,j).

- **Edge gradient at order k**:
  `∇ Lob_k^scaled(λ_j-λ_i, λ_i+λ_j)`
  — these are *exactly* the DoFs `nograds=True` drops.

- **Edge type II at order k** (Zaglmayr Sec 5.4):
  `Lob_k^scaled(λ_j-λ_i, λ_i+λ_j) · (λ_i ∇λ_j − λ_j ∇λ_i)`.

- **Verification**:
  `Curl[edge_gradient_shape] == 0` should hold symbolically;
  `Curl[edge_type_II_shape]` is nonzero.

## Sanity tests to run

1. `Curl` of every gradient-type shape function is zero.
2. Mass matrix of full HCurl basis vs `nograds`-pruned basis:
   the smaller matrix should be the principal submatrix.
3. Stiffness matrix (`(curl·, curl·)`) has a much smaller kernel under
   `nograds=True` (only the lowest-order ker(curl) remains).

## Cross-reference

- `radia_mcp.mathematica.recipes` — general Mathematica recipes
- `radia_mcp.mathematica.recipes.KELVIN_TRANSFORM` — symbolic Kelvin
"""


DECISION_TREE = r"""
# Quick decision tree: HCurl gauge / performance questions

User asks about NGSolve HCurl performance / DoF count / gauge:

```
1. Solving magnetostatic / eddy current (curl-curl + low-order mass)?
   ├── ungauged + CompactAMS (★ lab default; AMS handles kernel)
   │   → matrix_solvers_preconditioners('ams_hiptmair_xu')
   │
   ├── nograds=True (basis-level optimization)
   │   → fem_elements('ngsolve_hierarchy') ← THIS module
   │
   └── Both combined (max DoF reduction + good condition number)

2. Eigenvalue analysis or time-domain wave?
   → Do NOT use nograds=True
   → Gradient eigenspace at ω=0 is part of the problem

3. Multiply-connected domain (torus, body with hole)?
   → tree-cotree + cohomology cuts
   → fem_gauge_open_boundary('gauge')
   → bem.surface_ie.pmchwt for surface BEM equivalent

4. Need symbolic verification?
   → Use Mathematica recipes above
   → fem_elements('ngsolve_hierarchy', topic='mathematica_verification')
```

## Cross-references

- `radia_mcp.matrix_solvers.preconditioners.ams_hiptmair_xu` — AMS
  preconditioner (lab default, removes the need for nograds in most
  cases)
- `radia_mcp.fem.elements.edge` — Nédélec edge element foundations
- `radia_mcp.fem.gauge_open_boundary.gauge` — full gauge catalog
- `radia_mcp.differential_forms` — FEEC viewpoint
"""


def get_ngsolve_hierarchy_knowledge(topic: str = "decision_tree") -> str:
    """Dispatch NGSolve hierarchical H(curl) basis topics.

    Topics:
        zaglmayr_split          - Hierarchical split (gradient + rotational)
        nograds_option          - NGSolve `nograds=True` behavior
        tree_cotree             - Tree-cotree gauge + generalizations
        mathematica_verification - Symbolic verification recipes (Wolfram Lang)
        decision_tree           - Which strategy for which question (DEFAULT)
        all                     - Everything
    """
    topic = topic.lower().strip()
    if topic in ("zaglmayr_split", "zaglmayr", "split", "hierarchy"):
        return ZAGLMAYR_SPLIT
    if topic in ("nograds_option", "nograds"):
        return NOGRADS_OPTION
    if topic in ("tree_cotree", "tree-cotree", "treecotree"):
        return TREE_COTREE
    if topic in ("mathematica_verification", "mathematica", "verification",
                 "wolfram"):
        return MATHEMATICA_VERIFICATION
    if topic in ("decision_tree", "decision", "choose", "tree"):
        return DECISION_TREE
    if topic == "all":
        return "\n\n".join([ZAGLMAYR_SPLIT, NOGRADS_OPTION, TREE_COTREE,
                            MATHEMATICA_VERIFICATION, DECISION_TREE])
    return (f"Unknown topic '{topic}'. Available: zaglmayr_split, "
            "nograds_option, tree_cotree, mathematica_verification, "
            "decision_tree, all.")
