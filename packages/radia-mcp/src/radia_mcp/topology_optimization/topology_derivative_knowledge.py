"""
Topological derivative for nonlinear magnetostatics: changing the
TOPOLOGY (adding/removing material), not just the shape boundary.
"""

TD_OVERVIEW = """
# Topological derivative — when shape derivative is not enough

Shape derivative deforms an EXISTING boundary smoothly.  But what
if the optimal design has DIFFERENT TOPOLOGY than the initial?
E.g. a magnet that needs a hole, or an iron core that should be
split into separate pieces?

**Topological derivative** measures the cost of inserting a small
inclusion (e.g. a hole) at a point x₀ ∈ Ω.

## Definition (Sokolowski-Zochowski)

For a cost functional J(Ω) and a perturbation Ω_ε = Ω \\ B(x₀, ε)
(remove a small ball of radius ε around x₀):
    dJ_top(x₀) := lim_{ε → 0+} [J(Ω_ε) - J(Ω)] / ρ(ε)
where ρ(ε) is a problem-specific scaling (typically πε² in 2D, 4πε³/3 in 3D).

Sign convention: dJ_top(x₀) < 0  ⇒  inserting a hole at x₀ DECREASES
the cost → x₀ is a good location to remove material.

## Topology optimization algorithm

Iterate:
1. Solve state equation on current Ω
2. Solve adjoint equation
3. Compute dJ_top(x) for all x ∈ Ω
4. Insert holes where dJ_top is most negative (below threshold)
5. Repeat until no further improvement

This produces designs with the RIGHT TOPOLOGY (number of holes,
disconnected components) without requiring an a priori topological
ansatz.

## Combined SD + TD methods

Most modern topology-opt codes combine:
- TD for early iterations: fast, finds the right topology
- SD for refinement: smooth boundary, fine geometry

The resulting shapes have both correct topology AND smooth contours.

## Pitfalls for nonlinear problems

For linear PDE constraints, TD formulas are well-developed
(Sokolowski-Zochowski 1999).  For nonlinear PDEs (like
magnetostatics with saturation), TD is harder:
- The asymptotic analysis as ε → 0 must account for nonlinearity
- Different B-H curves give different TD formulas
- Rigorous proofs are recent (Beck-Sturm 2018, Amstutz-Gangl 2019)

Gangl PhD thesis Part I covers the nonlinear TD theory rigorously
for magnetostatics; it's still active research.

## Connection to level-set methods

Topology-derivative methods produce a SCALAR FIELD dJ_top(x) on Ω.
Level-set methods represent the design via a level-set function
φ(x), with Ω = {φ < 0}.  The two are closely related:
- TD gives the local velocity dφ/dt = -dJ_top(x)
- Level-set evolves φ via a Hamilton-Jacobi equation
- After convergence: optimal Ω is the zero-level-set of φ

This is the de-facto standard for modern topology optimization
in mechanical and EM design (Allaire 2002, Wang-Wang-Guo 2003).

## In Radia/NGSolve practice

NGSolve doesn't have built-in TD support, but it can be implemented
on top of the standard FEM machinery:
1. Standard FEM solve for u and adjoint p
2. Hand-derived TD formula for the cost J  (Gangl thesis gives
   explicit forms for tracking-type costs in magnetostatics)
3. Evaluate dJ_top at quadrature points
4. Update topology via element-wise material switch
   (iron → air or vice versa)

For motor optimization, this typically gives a permanent-magnet
arrangement DIFFERENT from the initial — e.g. dividing one magnet
into two with optimal gap.

For Radia: not directly applicable (Radia uses fixed object lists,
not a topology-variable mesh), but can be used to suggest design
changes before encoding in Radia.
"""


def get_topology_derivative_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"  — Topological derivative overview, algorithm, nonlinear caveats
    """
    return TD_OVERVIEW
