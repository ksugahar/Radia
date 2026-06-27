"""
Chain complex, boundary operator, homology, Betti numbers, and the
tree-cotree decomposition / belted-tree technique for gauging.

Distilled from:
  - Bossavit 1998 Ch.5 §5.3 (homology and chains, trees and cotrees,
    Stokes ⟨dω,c⟩ = ⟨ω,∂c⟩, belted trees for non-contractible domains)
  - Kameari 2011 §1, fig.2 (3-D tree-cotree decomposition, the
    Albanese-Rubinacci T-method, "tension" terminology)
  - Hilton-Wylie 1965, Greenberg-Harper 1981 (standard algebraic
    topology references)
"""

CHAIN_COMPLEX = """
# Chain complex: dual to the de Rham complex

## p-chains
For a simplicial mesh, a **p-chain** c is a formal linear combination

    c = ∑_{s ∈ S_p} c_s · s

where S_p is the set of p-simplices and c_s is a real (or integer)
coefficient.  The space of p-chains W_p is isomorphic to R^|S_p|.

Intuitively:
- 0-chain: weighted collection of nodes
- 1-chain: weighted directed curve (multi-edge path)
- 2-chain: oriented polyhedral surface (multi-face surface)
- 3-chain: weighted collection of cells (a region of space)

## Boundary operator ∂
Linear, ∂ : W_p → W_{p-1}, defined on individual simplices:

    ∂(node)  = 0    (or a formal (-1)-chain, ignored)
    ∂{m,n}   = n − m
    ∂{l,m,n} = {m,n} − {l,n} + {l,m}
    ∂(tet)   = signed sum of its 4 faces.

Matrix representation: G^T, R^T, D^T (transposes of the Whitney
incidence matrices) — explicit in terms of node-edge, edge-face,
face-cell incidence.

## ∂∂ = 0
"Boundary of a boundary is empty."  In coordinates:  G^T R^T = 0,
R^T D^T = 0.  This is the dual of  dd = 0  (the discrete differential
property d² = 0).

## Cycles and boundaries
- A chain c is a **cycle** if ∂c = 0 (closed chain).
- A chain c is a **boundary** if c = ∂γ for some (p+1)-chain γ.

By ∂∂=0, every boundary is a cycle.  The converse — every cycle is a
boundary — is the statement of TOPOLOGICAL TRIVIALITY (b_p = 0).

## Homology groups
    H_p := ker(∂ : W_p) / im(∂ : W_{p+1})
         = cycles / boundaries.

Their dimensions are the Betti numbers  b_p.  These count:
- b_0 = number of connected components
- b_1 = number of independent 1-loops (e.g., torus loops)
- b_2 = number of independent 2-cycles (e.g., enclosed cavities)
- b_3 = 0 for bounded domains

## de Rham theorem
The homology dimensions  b_p  COINCIDE with the dimensions of the
de Rham cohomology groups  H^p_{dR}.  This is the central insight of
de Rham (1931).

In other words: forms and chains see the SAME topology.

## Stokes as duality
The Stokes formula ∫_σ dω = ∫_∂σ ω can be rewritten as the duality

    ⟨dω, c⟩ = ⟨ω, ∂c⟩

between Whitney forms in W^p and chains in W_{p+1}.  The matrices
representing  d  (= G, R, D) and  ∂  (= G^T, R^T, D^T) are transposed
of each other.  This is **the** central duality in computational EM.
"""

TREE_COTREE = """
# Tree-cotree decomposition: discrete gauge fixing

## The problem
On a contractible domain, the kernel of  rot : W^1 → W^2  is exactly
the image of  grad : W^0 → W^1.  In matrix form,
    ker(R) = im(G).

If we wish to solve the curl-curl problem  rot(ν rot a) = j  for the
magnetic vector potential a ∈ W^1, the system matrix is singular: the
gauge  a → a + grad ϕ  leaves rot a unchanged.  We need to fix
this gauge.

## Tree of edges
A **spanning tree** T ⊂ E is a maximal subset of edges that contains
**no cycles** — visiting all N nodes with N−1 edges.  Algorithmically:
- Pick any node.
- Pick any edge from that node, add to T.  Mark its far endpoint.
- Repeat: pick any edge from a marked node to an unmarked node, add
  to T, mark the new node.
- Stop when all nodes are marked.

(Equivalent: find any submatrix of G of rank N−1.)

## Cotree
The **cotree** is E \\ T, the edges NOT in the tree.  Its size is
E − (N−1) = E − N + 1.

For a contractible domain, the Euler-Poincaré identity gives
N − E + F − T = 1, and via the exact-sequence dimension counts
the cotree size equals dim(im G^T) and provides exactly the right
number of independent constraints.

## Gauge fixing via tree
Set  a_e = 0  for every  e ∈ T.  This eliminates N−1 DoFs (the gauge
freedom).  The remaining E−(N−1) DoFs on cotree edges are independent;
the system becomes regular.

## In non-contractible domains: belted tree
If b_1 > 0 (the domain has loops), the tree alone is NOT enough.
There exist curl-free fields in W^1 that are NOT gradients (one per
loop), and these survive in the kernel even after tree-gauging.

Solution:  pick a "belt fastener" co-edge — one cotree edge for each
homology generator (loop) of D.  Add these b_1 edges to T to form
a **belted tree**.  Set a_e = 0 on the belted tree.

Total fixed DoFs: (N−1) + b_1.  The remaining DoFs are exactly
those of an irrotational + non-loop-circulation component.

## Why "tension"?
Kameari 2011 explains: the original Albanese-Rubinacci 1988 paper
used the word "tension" for the path integral of the electric vector
potential T along an edge — i.e., the edge DoF of the 1-form T.
The "tension" terminology comes from circuit theory where T plays
a role analogous to the EMF along a branch.  The modern term is
just "edge DoF" or "edge circulation".

## Practical algorithms
- BFS / DFS to build spanning tree.
- Decompose H^1 = ker(G^T) ⊕ im(G) (sum of cotree-edge basis and
  tree-edge gradient basis).
- For belt fasteners: identify representatives of H_1(D) by, e.g.,
  finding non-bounding co-edges by computing the cycle space of the
  cotree.

## Tree of faces (for div-side gauging)
By analogous duality, a "spanning tree of faces" of the dual mesh
(equivalently, a spanning tree in the cell-face graph of the primal
mesh) provides gauging for the div-side problem
    div(ε grad ϕ) = ρ
in mixed (H(div) × L^2) formulation.

## Beyond simple graph theory
Bossavit emphasizes:  these "tree" constructs are properly part of
**homology**, not graph theory.  In dim p = 1 and p = n−1 they map
to graph trees; for general p the right framework is the matrix
representation of ∂ and its rank deficiency.  In high dimensions
"trees" must be computed via algebraic methods (matrix factorization),
not graph traversal.
"""

EULER_POINCARE = """
# Euler-Poincaré formula and Betti numbers in practice

## The formula
For a finite simplicial complex (a mesh):

    χ(D) = ∑_{p=0}^n (−1)^p |S_p|  = ∑_p (−1)^p b_p

In 3-D:
    N − E + F − T = b_0 − b_1 + b_2 − b_3

For bounded D:  b_3 = 0.  For connected D:  b_0 = 1.  So:

    χ = 1 − b_1 + b_2.

## Common topologies

| D | b_0 | b_1 | b_2 | χ | Example |
|---|-----|-----|-----|---|---------|
| Ball / contractible       | 1 | 0 | 0 | 1 | Solid sphere, brick |
| Solid torus               | 1 | 1 | 0 | 0 | Donut, ring |
| Solid double torus        | 1 | 2 | 0 | −1 | Two linked rings |
| Hollow sphere shell       | 1 | 0 | 1 | 2 | Spherical shell |
| Thick torus (T² × interval)| 1 | 2 | 1 | 0 | Cooling jacket |

## How Betti numbers manifest in FEM
- b_1 > 0: scalar potential  ϕ  needs "cuts" — additional 1-cohomology
  generators that close circulation around loops.  Implementation:
  add b_1 extra DoFs per loop.
- b_2 > 0: vector potential a is not uniquely determined by b = rot a;
  one needs additional surface conditions ("solid harmonic
  representations").  Implementation: add b_2 extra DoFs per cavity.

## Cut-detection algorithms
- Make-Cuts algorithm (Kotiuga 1989): identify b_1 cuts by solving
  rank-deficient Laplacian.
- Spanning-Tree-of-Cotree algorithm: alternate spanning tree
  construction.

## radia implementation (gmsh-free) — the T-Omega cut engine
The lab realises the b_1 cuts WITHOUT gmsh.  `radia.cohomology`
(`src/radia/cohomology.py`) computes the first cohomology of a tetrahedral
NGSolve/Netgen mesh via the combinatorial Hodge Laplacian
`L1 = G G^T + Curl^T Curl`  (nullity = b_1), and returns one curl-free,
UNIT-circulation HCurl(order 0) generator `h_k` per loop — the discrete
`H^1` representative, i.e. exactly a T-Omega "cut" / one-ampere-turn loop
field.  This drops the Gmsh `computeHomology` + `.msh -> .vol` dependency
(the Pellikka 2013 path) entirely; NGSolve/Netgen ship no homology solver,
so this is the lab's standalone engine.

The total-scalar-potential (T-Omega) solver
`radia.cohomology_cut.CohomologyCutSolver` consumes it:

    H = -grad(phi) + sum_k NI_k h_k ,
    int mu grad(phi).grad(v) = sum_k NI_k int mu h_k . grad(v) ,

making the otherwise MULTI-valued magnetic scalar potential single-valued in
a current-linking (multiply-connected) region.  `oint H.dl = NI` holds
EXACTLY (by the unit-circulation cut).  This is the **primary** use of
cohomology in radia — the T-Omega cut; the Clebsch-hodograph current-linking
case is one instance of the SAME `H^1` obstruction (the scalar coordinate is
multi-valued iff a current threads a hole, period != 0).  Verified:
`validation_test/feec/test_cohomology.py` (b_1 + the curl-free unit-circulation cut
basis) and `validation_test/feec/test_tomega_cohomology_cut.py` (the straight-wire
T-Omega solve: `H_phi = I/(2 pi r)`, `oint H.dl = NI`).  Human-facing
showcase: executed/rendered `docs/cohomology/tomega_wire.ipynb`, with the
retired examples-path source preserved in
`docs/cohomology/cohomology_examples_archive_results.json`.

## Sanity check for meshes
Before running a curl-curl solver, compute χ from your mesh.  If it
doesn't match the expected topology (1 for a contractible volume),
there is a mesh defect — disconnected boundary, missing fan
tetrahedra, accidentally created handle, etc.

## Mesh-independence
Different meshes of the same continuous domain give the same Betti
numbers (provided the mesh is topologically faithful).  This is the
content of the **simplicial approximation theorem**.
"""


def get_homology_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "chain_complex"     — chains, ∂, cycles, boundaries, homology
      "tree_cotree"       — tree-cotree decomposition, belted trees
      "euler_poincare"    — N − E + F − T = χ, Betti numbers, cuts
    """
    topic = topic.lower().strip()
    if topic in ("chain_complex", "chain", "boundary", "homology"):
        return CHAIN_COMPLEX
    if topic in ("tree_cotree", "tree", "cotree", "gauge", "belted"):
        return TREE_COTREE
    if topic in ("euler_poincare", "euler", "betti"):
        return EULER_POINCARE
    if topic == "all":
        return "\n\n".join([CHAIN_COMPLEX, TREE_COTREE, EULER_POINCARE])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, chain_complex, tree_cotree, euler_poincare."
    )
