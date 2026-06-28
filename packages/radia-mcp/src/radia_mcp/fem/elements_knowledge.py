"""Element technology: edge (Nedelec/RT), high-order, XFEM, Isogeometric, DG."""

CATALOG = r"""
# Element technology for EM-FEM

| Family | Order | Continuity | Use |
|--------|-------|------------|-----|
| H1 (nodal Lagrange) | 1, 2, 3, ... p | full C^0 | scalar Omega, V, heat |
| HCurl (Nedelec / edge) | 0, 1, 2, ... | tangential | vector A, T, H |
| HDiv (RT / Raviart-Thomas) | 0, 1, 2, ... | normal | B, J (mixed formulations) |
| HDivSurface | 0, 1, 2, ... | tangential on surface | BEM surface currents (RWG) |
| L2 | discontinuous | none | DG fluxes, dual cells |

For magnetostatic / magnetodynamic EM, the natural pairing is:
- **A in HCurl, B in HDiv, ω in H1** (the de Rham complex)
- See `radia_mcp.differential_forms.complex` for the cohomology view

## Hierarchical vs nodal high-order

| Approach | Library | Pros | Cons |
|----------|---------|------|------|
| Hierarchical p (Ainsworth-Coyle) | NGSolve `order=p` | reuse low-order matrices | basis cancellation at high p |
| Modal nodal (Lagrange) | classic FEM | simpler conditioning | recompute basis for each p |
| Lobatto basis | NGSolve internal | good for p-version | requires careful implementation |

NGSolve uses **hierarchical orthogonal** basis (Schoeberl-Zaglmayr), which
is the lab production choice.

## References (this folder)

- [LOCAL] `10_FEM_定式化/12_辺要素/` (9 files) ← edge element foundations
- [LOCAL] `10_FEM_定式化/13_高次要素/` (4 files)
- [LOCAL] `10_FEM_定式化/14_xfem/` (2 files)
- [LOCAL] `10_FEM_定式化/15_Isogeometric/` (2 files)
- [LOCAL] `10_FEM_定式化/16_不連続ガレルキン法/` (7 files)
"""


EDGE_ELEMENTS = r"""
# Edge elements (Nédélec / Raviart-Thomas) — the foundation of HCurl FEM

Reference: Nédélec 1980, "Mixed finite elements in R^3", Numer. Math.
35:315-341.

## Why edge elements (vs nodal vector)

If you discretize A with **nodal** vector basis (each node has 3 DOFs Ax,
Ay, Az):
- ✗ Spurious modes appear in eigenvalue problems
- ✗ Tangential continuity is OVER-enforced (full continuity at edges/faces)
- ✗ Cannot represent gradient kernel naturally

If you use **edge** elements (tangential DOF per mesh edge):
- ✓ B = curl A is in HDiv naturally (commuting diagram)
- ✓ Gradient kernel exactly representable (∇H1 ⊂ HCurl)
- ✓ Tangential continuity preserved at material interfaces
- ✓ Natural BC: t × A = 0 (PEC) implements naturally

## Element types

Lowest-order Nédélec (NED1, k=0):
- 1 DOF per edge = ∫_edge A · t dl
- Linear basis on tetrahedra / hexahedra

Higher-order (Schoeberl-Zaglmayr in NGSolve):
- Hierarchical edge / face / cell basis
- `HCurl(mesh, order=p)` for arbitrary p

## In NGSolve

```python
from ngsolve import *
fes = HCurl(mesh, order=2)
u, v = fes.TnT()
a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
a += 1j * omega * sigma * u * v * dx("conductor")
```

## Reference papers (this folder)

[LOCAL] `10_FEM_定式化/12_辺要素/` (9 files) — Bossavit, Nédélec foundations
[LOCAL] `10_FEM_定式化/13_高次要素/Hierarchal_vector_basis_functions_*.pdf`
        — Hierarchical vector basis (key Schoeberl reference)

## Pitfalls

| Issue | Solution |
|-------|----------|
| HCurl with HDiv-conforming B problems | Use HDiv for B explicitly + Lagrange multiplier |
| HCurl with H1-conforming gauge | NGSolve `H1amg` only for H1; use CompactAMS for HCurl |
| Lowest-order NED1 with curved geometry | mesh.Curve(p) needed (Cubit netgen export order p) |
"""


HIGH_ORDER = r"""
# High-order elements (p-version)

Increase polynomial order p in same mesh → exponential convergence for
smooth problems (vs h-refinement = algebraic O(h^p)).

## When p helps

✓ Smooth solution (no singularities)
✓ Tolerable per-element cost increase (DOFs grow as p^d)
✓ Iterative solver doesn't degrade with p (Compact AMS works to ~p=10)

## When p doesn't help

✗ Solution has corner singularities (concentrate h instead → adaptive)
✗ Strong material discontinuities (interface jumps)
✗ Memory limited (p=4 is ~50x DOFs of p=1 for tet)

## hp-adaptive

Combine: h-refine near singularities, p-refine in smooth regions.
NGSolve supports this via element-wise order.

## Curved boundary geometry

For order p ≥ 2, the **mesh must be curved** to match analytic boundary,
or accuracy is capped at O(h^2) regardless of p.

```python
# NGSolve
mesh.Curve(p)   # geometric refinement matching FE order
```

For Cubit-generated mesh:
```bash
export netgen "model.vol" order 3   # exports with curved nodes
```

## References (this folder)

[LOCAL] `10_FEM_定式化/13_高次要素/` (4 files) — Schoeberl-Zaglmayr basis,
        Hierarchical vector basis, evaluation strategies

## Lab practice

- Default p=1 for prototyping
- p=2-3 for production accuracy (docs mesh-evaluation demo sweeps these)
- p=4-5 for academic validation (slower convergence due to setup cost)
"""


XFEM = r"""
# XFEM (Extended FEM) — cracks, interfaces, multiple materials

eXtended Finite Element Method enriches the standard basis with
discontinuous / singular functions to capture features that don't align
with the mesh.

## Use cases in EM

| Application | Enrichment |
|-------------|------------|
| Crack-tip singularity | r^{1/2} Heaviside |
| Material interface (μ jump) | step function across interface |
| Moving conductor / rotor | enriched local basis without remeshing |

## Lab use

- Not yet production. Reference layer only.
- Could complement `radia_mcp.ndt` (eddy current crack inspection)

## References (this folder)

[LOCAL] `10_FEM_定式化/14_xfem/` (2 files)

## Alternative: cut-FEM

Modern variant where the mesh cuts through interfaces without alignment.
NGSolve has cut-FEM support via `ngsxfem` (third-party extension).
"""


ISOGEOMETRIC = r"""
# Isogeometric Analysis (IGA) — NURBS / spline FE basis

Reference: Hughes-Cottrell-Bazilevs 2005, "Isogeometric analysis: CAD,
finite elements, NURBS, exact geometry and mesh refinement", CMAME 194.

## Idea

Use the same NURBS / B-spline basis for BOTH geometry (CAD) AND FEM
solution:
- Exact geometry preserved at all refinement levels
- Higher continuity (C^{p-1} across element boundaries)
- Better high-order approximation for smooth solutions

## EM-specific isogeometric

For HCurl / HDiv on NURBS:
- Buffa et al. constructed isogeometric de Rham complex
- Bembel library (open-source, German group) for isogeometric BEM
- Marussig "Fast Isogeometric BEM" arxiv 2014

## When IGA wins

✓ Smooth curved geometry (e.g. rotating machine rotor)
✓ Need C^1 across element boundaries (e.g. plate bending)
✓ Tight integration with CAD pipeline

## When IGA loses

✗ Sharp corners / multi-material interfaces
✗ Existing FE infrastructure assumes Lagrange
✗ Steep learning curve

## Lab use

- Not production. NGSolve does NOT have native IGA.
- Reference layer; cross-link to BEM IGA papers in `radia_mcp.bem`

## References (this folder)

[LOCAL] `10_FEM_定式化/15_Isogeometric/` (2 files)
"""


DG = r"""
# Discontinuous Galerkin (DG)

DG allows the FE function to be discontinuous across element boundaries.
Continuity enforced WEAKLY via numerical fluxes.

## Why DG for EM

✓ Local conservation (good for time-domain Maxwell)
✓ Easy h-adaptivity (no continuity constraints)
✓ Natural for shocks / discontinuous solutions

✗ Many more DOFs (no nodal sharing)
✗ Flux choice non-trivial (upwind, central, etc.)
✗ Imposes spurious modes if not careful (need consistent fluxes)

## DG-FEM for time-domain Maxwell

Nodal DG (Hesthaven-Warburton) — local time-stepping possible.
Often combined with explicit RK time integration.

## Lab use

- NOT production. Lab uses standard CG (continuous Galerkin) for MQS.
- Reference layer for full-wave time-domain (e.g. PCB simulation).

## References (this folder)

[LOCAL] `10_FEM_定式化/16_不連続ガレルキン法/` (7 files)
"""


def get_elements_knowledge(topic: str = "catalog") -> str:
    """Dispatch element technology topics.

    Topics:
        catalog        - H1/HCurl/HDiv/L2 + hierarchical (DEFAULT)
        edge           - Edge elements (Nedelec) — foundation of HCurl FEM
        high_order     - p-version (hierarchical basis)
        xfem           - eXtended FEM (cracks, interfaces)
        isogeometric   - Isogeometric Analysis (NURBS basis)
        dg             - Discontinuous Galerkin
        all            - Everything
    """
    topic = topic.lower().strip()
    if topic in ("catalog", "overview", "family"):
        return CATALOG
    if topic in ("edge", "edge_elements", "nedelec", "ned"):
        return EDGE_ELEMENTS
    if topic in ("high_order", "ho", "p_version", "hierarchical"):
        return HIGH_ORDER
    if topic == "xfem":
        return XFEM
    if topic in ("isogeometric", "iga", "nurbs"):
        return ISOGEOMETRIC
    if topic == "dg":
        return DG
    if topic == "all":
        return "\n\n".join([CATALOG, EDGE_ELEMENTS, HIGH_ORDER, XFEM,
                            ISOGEOMETRIC, DG])
    return (f"Unknown topic '{topic}'. Available: catalog, edge, high_order, "
            "xfem, isogeometric, dg, all.")
