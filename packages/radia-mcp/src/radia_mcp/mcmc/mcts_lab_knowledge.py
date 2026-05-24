"""Monte Carlo Tree Search (MCTS) for EM design — lab lineage.

Hybrid storage: most constants are r-strings here, but the
HAYAHO_PM_MOTOR section is loaded from prompts/hayaho_pm_motor.md
as a demo of the COMSOL MCP-inspired pattern (avoids the docstring-
inside-r-string parser bug we hit twice).
"""

from ..common import load_prompt

MCTS_OVERVIEW = r"""
# Monte Carlo Tree Search (MCTS) Overview

MCTS is a tree-search algorithm with statistics. Each node represents
a partial decision sequence; each leaf represents a complete choice.
The 4 phases per iteration:

1. **Selection**: From the root, descend the tree by picking children
   that maximize an "upper confidence bound" (UCB1):

       UCB1(child) = mean_score(child) + C * sqrt( ln(parent_visits) /
                                                    child_visits )

   The first term exploits high-scoring children; the second explores
   under-visited ones. `C` (~1.4 for games, ~3.0 for engineering) sets
   the exploit/explore balance.

2. **Expansion**: At a not-fully-visited node, add a new child.

3. **Simulation** (rollout): From the new child, run a random / heuristic
   playout to a terminal state; compute the score.

4. **Backpropagation**: Update visit counts and mean scores along the
   path back to the root.

## Why MCTS for engineering design

Engineering design = sequence of discrete + continuous decisions:
- Pick motor type (PMSM / SRM / IM) [discrete]
- Pick #poles [discrete]
- Pick magnet shape (I, V, U, ∇) [discrete]
- Pick rotor geometry [continuous, via NGnet on/off or SIMP]

MCTS naturally handles the **discrete tree** of global configurations
+ **continuous leaf optimization** (CMA-ES, NSGA-II) — combining
combinatorial and shape optimization.

Compare:
- **Genetic algorithm**: works on fixed-length chromosomes; struggles
  with variable-structure designs.
- **MCTS**: each path through the tree can have different lengths
  (different #PMs → different shape params count). Native.
- **Bayesian Optimization**: same conditional space possible via
  Optuna's define-by-run, but BO doesn't separate "structural
  decision" from "shape refinement".

## MCTS in the AI canon

| Year | System | Domain |
|------|--------|--------|
| 2006 | UCT (Kocsis-Szepesvári) | Theoretical foundation |
| 2016 | AlphaGo | Go (DeepMind) — MCTS + DNN policy/value |
| 2017 | AlphaZero | Chess + Shogi + Go (one network) |
| 2020 | MuZero | Atari + Go (no rules) |
| 2021+ | EM Design | **Sato 2023, Yin 2024 (Hokkaido)** |
| 2024+ | Drug discovery, materials | Generative chemistry |

The Hokkaido lab (Sato Hayaho, Yin Shuli, advisor Igarashi Hajime)
pioneered MCTS for EM device design, transferring AlphaGo-style ideas
to motor / inductor optimization.
"""


HAYAHO_PM_MOTOR = load_prompt("mcmc", "hayaho_pm_motor")

_HAYAHO_PM_MOTOR_INLINE_DEPRECATED = r"""
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
    # Return the child maximizing UCB1.
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
    # Run NSGA-II at leaf to get Pareto solutions.
    # config_path = [(Pn, ...), (beta, ...), (PM_type, ...), (n_PM, ...)]
    Pn, beta, pm_type, n_pm = (n.value for n in config_path)
    # NSGA-II over (NGnet weights w, PM shape params p)
    pareto = nsga2(motor_fem,
                    bounds=ngnet_bounds(Pn, pm_type, n_pm),
                    objectives=["-Tavg", "Trip"],
                    n_gen=50, pop=50)
    return pareto

def backpropagate(leaf_node, new_pareto):
    # Update visits + score along path to root.
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
"""


YIN_INDUCTOR = load_prompt("mcmc", "yin_inductor")

_YIN_INDUCTOR_INLINE_DEPRECATED = r"""
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
"""


MCTS_LINEAGE = r"""
# MCTS lineage and pointers

## Pre-AlphaGo (theory)

- 1948-50: Stanislaw Ulam invents Monte Carlo (Los Alamos)
- 1973: Robbins-Lai bandit theory → UCB1
- 2006: Kocsis-Szepesvári, "Bandit based Monte-Carlo Planning"
  (the original UCT paper)
- 2008: Sylvain Gelly et al., "MoGo" — first Go program > 8-kyu
- 2012: Browne-Powley survey, "A Survey of Monte Carlo Tree Search
  Methods" (the most-cited MCTS review)

## DeepMind canon

- 2016: AlphaGo (Nature) — MCTS + DNN policy + value
- 2017: AlphaZero (Science) — same algo for Chess + Shogi + Go
- 2020: MuZero (Nature) — MCTS + learned dynamics model

## Engineering applications

- 2018-2022: Sato Hayaho early MCTS for SR motor / PM motor
  (Igarashi lab Hokkaido)
- 2023: Sato-Igarashi multi-objective MCTS for PM motor (this MCP)
- 2024: Yin-Sato-Igarashi MCTS for inductor (this MCP)
- 2023-2024: Hayaho-MCTS+TO papers (3 in `03_最適化_モンテカルロ/`)

## What to read next

1. **Browne 2012 survey** — the comprehensive MCTS reference
2. **Silver et al. 2017 AlphaZero** — modern MCTS + NN integration
3. **Sato 2023** — direct lab reference for PM motor MCTS
4. **Yin 2024** — direct lab reference for inductor MCTS + inheritance

For the underlying Bayesian / probability theory:
- Kocsis-Szepesvári 2006 — UCT proof (regret bound)
- Auer-Cesa-Bianchi-Fischer 2002 — UCB1 paper (the multi-armed
  bandit foundation)

## Lab papers also worth reading

- `Hayaho-Automatic_Design_of_PM_Motor_Using_MCTS_in_Conjunction_With_TO.pdf`
  (single-objective predecessor of Sato 2023)
"""


def get_mcts_lab_documentation(topic: str = "all") -> str:
    """Dispatch MCTS topics.

    Topics:
        overview        - MCTS algorithm + AI canon
        hayaho_pm_motor - Sato-Igarashi 2023 MO MCTS for PM motor ★
        yin_inductor    - Yin-Sato-Igarashi 2024 MCTS for inductor ★
        lineage         - MCTS history + reading recommendations
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro", "mcts_overview"):
        return MCTS_OVERVIEW
    if topic in ("hayaho_pm_motor", "hayaho", "sato_2023", "pm_motor",
                 "sato_igarashi", "multi_objective_mcts"):
        return HAYAHO_PM_MOTOR
    if topic in ("yin_inductor", "yin", "inductor_mcts", "yin_2024",
                 "inductor"):
        return YIN_INDUCTOR
    if topic in ("lineage", "history", "alphago", "reading"):
        return MCTS_LINEAGE
    if topic == "all":
        return "\n\n".join([MCTS_OVERVIEW, HAYAHO_PM_MOTOR,
                            YIN_INDUCTOR, MCTS_LINEAGE])
    return (f"Unknown topic '{topic}'. Available: overview, "
            "hayaho_pm_motor, yin_inductor, lineage, all.")
