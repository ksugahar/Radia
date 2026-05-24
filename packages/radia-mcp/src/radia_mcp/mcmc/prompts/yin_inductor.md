# Yin-Sato-Igarashi 2024: MCTS for Inductor Design

Reference: S. Yin, H. Sato, H. Igarashi, "A Comprehensive Optimal
Design of Inductors Using Monte Carlo Tree Search", IEEE Trans.
Magnetics 60(3):8400504, March 2024.
DOI: 10.1109/TMAG.2023.3308214. Hokkaido University.

Lab copy: W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/
A_Comprehensive_Optimal_Design_of_Inductors_Using_Monte_Carlo_Tree_Search.pdf

## Problem

Design inductor for switched-mode power supply:
- Variables: core material (8 options) + core length (12 options) +
  core height (6 options) + #turns horizontal (m) + #turns vertical (n)
- Per-config, optimize geometric parameters via CMA-ES
- Target: specified L_0 and I_sat; minimize R_dc and P_loss
- Considers **nonlinear core saturation** (B-H curve per material)

## Key innovations over Sato 2023

1. **Nonlinear material selection**: B-H curves of 8 candidate
   materials are part of the tree, with FEM Newton-Raphson saturation
   modeling at the simulation step.

2. **Variable-length chromosome via tree depth**: different paths
   through the tree have different numbers of geometric parameters
   (3-8). Conventional GA struggles with this; MCTS handles it
   natively.

3. **Inherited search (warm start across design targets)**: For a
   sequence of design problems with similar specifications, reuse
   the previous tree as the starting point. Empirically reaches
   the optimum in 30-50% fewer iterations.

4. **Alternative path prediction**: After search, identify the path
   along which the **mean reward at each node is maximized** (not
   just the best-ever sample). This "statistical optimum" may differ
   from the best-ever found sample → provides a robustness check.

## Objective function

    f_MCTS = -C2 * (
        w1 * |1 - L/L_0|
      + w2 * |1 - I_sat/I_0|
      + w3 * |R_dc/R_0|
      + w4 * |P_loss/P_0|
    ) -> max

L is computed from `Integrate(J_z * A_z * dV) / I^2` (energy
formulation).

## Algorithm

```
For round in 1..N_rounds:
  1. Selection: descend tree via UCT
       UCT(p_i, p) = mean(f_MCTS(p_i)) + c * sqrt(ln N(p) / N(p_i))
  2. Expansion: if not fully visited, expand new child
  3. Simulation: from leaf, random policy to terminal state → CMA-ES
     over remaining geometry → score
  4. Backpropagation: update N and mean(f_MCTS) along path from leaf
     to root  (note: ONLY along the selected path, NOT random-policy
     intermediates -- differs from Sato 2023)
```

## Inherited search recipe

```python
def search_with_inheritance(prior_tree, new_target):
    # Resume MCTS from a previous search's tree.
    tree = copy.deepcopy(prior_tree)
    # Reset visit counts? -- Yin 2024 keeps them
    # Update objective function to new_target
    for _ in range(N_rounds_new):
        path = select(tree)
        leaf_pareto = simulate(path, target=new_target)
        backpropagate(path[-1], leaf_pareto)
    return tree
```

Useful for parameter sweeps (e.g. "find best inductor at 50/100/200
kHz" — search 50 kHz first, inherit for 100, inherit again for 200).

## Lab applicability

Direct lab use case: lab works on PEEC + SIBC for IH coil design.
Yin's approach maps to:
- Tree levels: (core material, core length/height, #turns) — same
  as inductor
- Leaf optimization: CMA-ES over coil placement / winding pattern
- Objective: L + R_dc + P_wp (with workpiece loss)
- Forward: PEEC + FEM-SIBC + Kelvin (lab production pipeline)

## Implementation status in radia-mcp

- **MCTS theory** — this document
- **MCTS framework code** — TODO (could promote to lab production)
- **Outer Optuna alternative** — `radia_mcp.optuna('coil_design')`

For an MCTS-style IH coil sweep today, the user would write the
algorithm by hand using the lab PEEC + FEM panels as the inner
forward.

## Cross-links

- `radia_mcp.peec` — PEEC backend for inductor forward
- `radia_mcp.ih` — IH workpiece SIBC
- `radia_mcp.topology_optimization` — gradient-based shape opt
- `radia_mcp.optuna` — TPE / NSGA-II / CMA-ES alternatives
