# Sato-Igarashi 2023: Multi-Objective MCTS for PM Motor (★ lab focus)

Reference: H. Sato, H. Igarashi, "Multi-Objective Automatic Design of
Permanent Magnet Motor Using Monte Carlo Tree Search", IEEE Trans.
Magnetics 59(5):8201304, May 2023. DOI: 10.1109/TMAG.2023.3254510.
Hokkaido University, Graduate School of Info Sci & Tech.

Lab copy: W:/.../04_機械学習と最適化/03_最適化_モンテカルロ/
Hayaho-Multi-Objective_Automatic_Design_of_PM_Motor_Using_MCTS.pdf

## Problem

Design PMSM (permanent magnet synchronous motor) for EVs:
- Maximize average torque `T_avg`
- Minimize torque ripple `T_rip` OR iron loss `P_iron`

Design variables:
- Global configuration `s`:
  - Number of poles `Pn ∈ {4, 6, 8}`
  - Current advance angle `beta`
  - PM type `∈ {I, V, U, ∇}`
  - Number of PMs (1 or 2)
- Local shape `r`:
  - NGnet on/off shape function weights `w` (rotor topology)
  - PM placement parameters `p`

## Algorithm

```
1. Build tree: 4 levels = (Pn, beta, PM type, #PMs)
   Each node stores n_i (visits), v_i (score)
2. for t in 0..tmax:
     Selection: descend root → leaf via UCB1
       P(p, i) = v_i + C * sqrt(ln n_p / n_i),  C = 3.0
     Optimization at leaf: run NSGA-II over (w, p) to get Pareto solutions
     Backpropagation: update v_i along the path
3. Sort and return all Pareto solutions from leaf nodes
```

## Novel scoring: number of Pareto solutions per node

Standard MCTS uses single-objective score. For multi-objective:

    v_i = N_Pareto(i, parent) / sum_j N_Pareto(j, parent)

i.e. v_i is the fraction of Pareto-front solutions that came from
child i. When backpropagating, the parent's Pareto front = union of
children's fronts (re-sorted).

This biases MCTS toward branches that produce diverse, well-spread
Pareto solutions.

## Results

| Case | Objectives | Pareto solutions |
|------|------------|------------------|
| I | T_avg vs T_rip | 8-pole + various PM shapes (U for high T, I for low ripple) |
| II | T_avg vs P_iron | 8-pole + various PM types |

- Pareto front from MCTS dominates random TO at T_avg <= 260 Nm.
- Random TO occasionally found better Pareto at very high T_avg
  (stochastic property of NSGA-II).
- 30 iterations × NSGA-II inner = 3 days on Xeon Platinum 8280 (4 ×
  28 cores = 224 threads).

## Implementation pattern (Python pseudocode)

```python
import numpy as np
from collections import defaultdict

class Node:
    def __init__(self, parent=None, value=None):
        self.parent   = parent
        self.value    = value          # config value at this level
        self.children = {}
        self.visits   = 0
        self.score    = 0.0
        self.pareto   = []             # list of (Tavg, Trip) tuples

def select(node, C=3.0):
    """Return the child maximizing UCB1."""
    if not node.children:
        return node
    best_p, best_child = -np.inf, None
    for child in node.children.values():
        if child.visits == 0:
            return child   # explore unvisited first
        p = child.score + C * np.sqrt(np.log(node.visits)/child.visits)
        if p > best_p:
            best_p, best_child = p, child
    return select(best_child, C)

def optimize_at_leaf(config_path):
    """Run NSGA-II at leaf to get Pareto solutions."""
    # config_path = [(Pn, ...), (beta, ...), (PM_type, ...), (n_PM, ...)]
    Pn, beta, pm_type, n_pm = (n.value for n in config_path)
    # NSGA-II over (NGnet weights w, PM shape params p)
    pareto = nsga2(motor_fem,
                    bounds=ngnet_bounds(Pn, pm_type, n_pm),
                    objectives=["-Tavg", "Trip"],
                    n_gen=50, pop=50)
    return pareto

def backpropagate(leaf_node, new_pareto):
    """Update visits + score along path to root."""
    node = leaf_node
    node.pareto = merge_pareto(node.pareto, new_pareto)
    while node.parent is not None:
        node.parent.pareto = merge_pareto_all_children(node.parent)
        # Score = fraction of Pareto solutions contributed by this child
        N_total = sum(len(c.pareto) for c in node.parent.children.values())
        for child in node.parent.children.values():
            child.score = len(child.pareto) / N_total
        node.visits += 1
        node = node.parent
    node.visits += 1   # root
```

## How to reproduce in radia/NGSolve

Replace `motor_fem` with `calc_motor_transient.py` (already in lab
panels). Replace NGnet on/off with SIMP or level-set
(`radia_mcp.topology_optimization`).

```python
# At each MCTS leaf:
from radia.panels.calc_motor_transient import run_motor_transient

def motor_fem(Pn, beta, pm_type, n_pm, ngnet_weights, pm_params):
    cfg = build_motor_config(Pn, beta, pm_type, n_pm,
                              ngnet_weights, pm_params)
    result = run_motor_transient(cfg)
    return result["T_avg"], result["T_rip"]
```

## Cross-links

- `radia_mcp.motor` — calc_motor_transient.py + Wakao autoencoder
  alternative topology opt
- `radia_mcp.topology_optimization` — SIMP / level-set / MMA for
  the inner shape refinement
- `radia_mcp.optuna('motor_topology')` — Optuna outer + SIMP inner
  alternative (no MCTS — uses TPE conditional search space instead)

## Comparison with Optuna for the same problem

| Aspect | MCTS (Sato 2023) | Optuna conditional |
|--------|------------------|---------------------|
| Tree structure | Native | Branched via `suggest_categorical` |
| Pareto scoring | Per-branch Pareto count | NSGA-II global Pareto |
| Branch revisit memory | Strong (tree statistics) | Weak (TPE univariate per branch) |
| Wallclock for 200 trials | ~3 days (Hokkaido) | similar |
| Implementation effort | Custom MCTS code | 1-page Optuna |

Recommendation: For repeated motor design (production), invest in
MCTS code (lab lineage). For one-off design (research), Optuna's
conditional search is faster to set up.
