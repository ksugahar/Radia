"""Evolutionary computation algorithms for EM optimization."""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "(see knowledge file for details)",
    "ga_de": "(see knowledge file for details)",
    "pso": "(see knowledge file for details)",
    "cma_es": "(see knowledge file for details)",
    "immune_nsga": "(see knowledge file for details)",
    "libraries": "(see knowledge file for details)",
    "all": "return \"\n\n\".join([OVERVIEW, GA_DE, PSO, CMA_ES, IMMUNE_NSGA,",
}

OVERVIEW = r"""
# Evolutionary Computation Overview

A family of stochastic, **population-based**, derivative-free
algorithms inspired by biological evolution. All share the basic
loop:

1. Maintain a population of N candidate solutions
2. Evaluate fitness (objective) of each
3. **Select** parents (better fitness more likely)
4. **Vary** (crossover + mutation) to make offspring
5. Replace some/all of the population with offspring
6. Repeat until budget exhausted

Variants differ in:
- Representation (binary, real, tree, permutation)
- Selection (tournament, roulette, rank, truncation)
- Variation (single-point crossover, SBX, UNDX, Gaussian mutation)
- Replacement (elitist, generational, steady-state, mu+lambda)

## When evolutionary computation wins

- **Discrete / mixed / structured** design spaces (graphs, trees,
  permutations)
- **Multi-modal** objective (many local optima)
- **No gradient** AND **no probabilistic structure** to exploit
- **Multi-objective** (NSGA-II, MOEA/D, SPEA2)
- **Parallel evaluation** is easy (population-based)

## When evolutionary computation loses

- Continuous + smooth + low-dim → CMA-ES / BO faster
- Differentiable → gradient descent O(d) vs EA O(d²)
- Bayesian framework available → variational inference

## EA family tree

```
Evolutionary Computation
├── Genetic Algorithm (GA) — Holland 1975
│   └── representation = binary / real, mutation + crossover
├── Genetic Programming (GP) — Koza 1992
│   └── representation = tree (programs)
├── Evolution Strategies (ES) — Rechenberg 1971
│   ├── (mu, lambda) ES, (mu + lambda) ES
│   └── CMA-ES (Hansen-Ostermeier 1996) ← gold standard
├── Differential Evolution (DE) — Storn-Price 1997
│   └── continuous, mutation = difference of population members
├── Particle Swarm Optimization (PSO) — Kennedy-Eberhart 1995
│   └── swarm-of-particles, velocity update
├── Ant Colony Optimization (ACO) — Dorigo 1992
│   └── discrete combinatorial (TSP, scheduling)
├── Artificial Immune System (AIS)
│   ├── Clonal selection, negative selection
│   └── lab papers: 免疫アルゴリズム (2 lab refs)
└── Multi-objective
    ├── NSGA-II (Deb 2002) ← industry standard
    ├── NSGA-III (Deb 2014) ← > 3 objectives
    ├── SPEA2 (Zitzler 2001)
    └── MOEA/D (Zhang 2007)
```

## Decision: which EA?

| Need | Use |
|------|-----|
| Real-valued, single-obj, < 100 dim | **CMA-ES** |
| Real-valued, multi-obj | **NSGA-II** (or NSGA-III for > 3 obj) |
| Binary / integer, single-obj | **DE** or **GA** with binary repr |
| Discrete / structured (graph, permutation) | **GA** custom operators |
| Continuous, parallel computer | **DE** (best parallel scaling) |
| Cheap evaluation + many trials | **PSO** (simplest implementation) |
| Multi-modal landscape, want all modes | **PSO with niching** or AIS |
"""


GA_DE = r"""
# Genetic Algorithm + Differential Evolution

## GA (Holland 1975) — the foundation

Population of binary strings → selection + crossover + mutation
→ next generation.

```python
import numpy as np

def genetic_algorithm(fitness_fn, pop_size=50, n_bits=20,
                      crossover_p=0.9, mutation_p=0.01,
                      n_generations=100, seed=42):
    rng = np.random.default_rng(seed)
    pop = rng.integers(0, 2, (pop_size, n_bits))

    for gen in range(n_generations):
        # Fitness
        fitness = np.array([fitness_fn(individual) for individual in pop])
        # Selection: tournament (size 3)
        parents = []
        for _ in range(pop_size):
            candidates = rng.choice(pop_size, 3, replace=False)
            winner = candidates[np.argmax(fitness[candidates])]
            parents.append(pop[winner])
        parents = np.array(parents)
        # Crossover: single point
        offspring = parents.copy()
        for i in range(0, pop_size, 2):
            if rng.uniform() < crossover_p:
                pt = rng.integers(1, n_bits)
                offspring[i,   pt:] = parents[i+1, pt:]
                offspring[i+1, pt:] = parents[i,   pt:]
        # Mutation
        mask = rng.uniform(size=(pop_size, n_bits)) < mutation_p
        offspring = np.logical_xor(offspring, mask).astype(int)
        pop = offspring

    final_fitness = np.array([fitness_fn(i) for i in pop])
    return pop[np.argmax(final_fitness)], np.max(final_fitness)
```

## Differential Evolution (Storn-Price 1997) ★ real-valued workhorse

For real-valued continuous optimization. Mutation = scaled
difference between population members:

    v_i = x_r1 + F * (x_r2 - x_r3)       (DE/rand/1)
    u_i = crossover(x_i, v_i)            (binomial)
    x_i_new = u_i if f(u_i) < f(x_i) else x_i

Knobs:
- `F` = differential weight (typical 0.5-0.8)
- `CR` = crossover probability (typical 0.7-0.9)
- Population size = 10*D (D = dim)

DE handles non-smooth, multi-modal landscapes better than CMA-ES
in some cases. Trivial to parallelize.

```python
from scipy.optimize import differential_evolution

bounds = [(-5, 5)] * 10
result = differential_evolution(my_objective, bounds,
                                 strategy='best1bin',
                                 popsize=15, mutation=(0.5, 1.0),
                                 recombination=0.7,
                                 seed=42, workers=-1)  # parallel
```

## Adaptive DE variants

- **jDE** (Brest 2006): adapt F + CR online
- **SHADE** (Tanabe-Fukunaga 2013): success-history based
- **LSHADE** (Tanabe-Fukunaga 2014): linear population size reduction
- These often top CEC competitions; available in `pyade` package.

## EM-flavored example: BH curve parameter fit via DE

```python
import numpy as np
from scipy.optimize import differential_evolution

# Frohlich-Kennelly: B(H) = mu_init * H / (1 + |H|/H_c)
def predict(H, theta):
    mu_init, H_c, B_sat = theta
    return B_sat * mu_init * H / (B_sat + mu_init * abs(H))

def loss(theta):
    H, B_meas = load_data()
    B_pred = predict(H, theta)
    return np.mean((B_pred - B_meas)**2)

bounds = [(100, 100000), (10, 5000), (0.5, 2.5)]
result = differential_evolution(loss, bounds, seed=42, workers=-1)
print(result.x, result.fun)
```
"""


PSO = r"""
# Particle Swarm Optimization (PSO)

Kennedy-Eberhart 1995. Inspired by bird flocking. Each "particle"
maintains:
- Position x_i
- Velocity v_i
- Personal best p_i (best position visited)
- Global best g (best position of swarm)

Update at each step:

    v_i ← w * v_i + c1 * r1 * (p_i - x_i) + c2 * r2 * (g - x_i)
    x_i ← x_i + v_i

where w = inertia (0.7-0.9), c1 = c2 = 2 (cognitive + social),
r1, r2 = U(0,1) random.

## Lab PSO variants (W:/.../02_最適化_進化計算/)

| Paper | Variant | Key idea |
|-------|---------|----------|
| `大域的最適化が可能なPSO` | Global PSO | Adaptive inertia → avoid premature convergence |
| `巣別れを導入したPSOアルゴリズム` | Hive-split PSO | Split swarm into sub-swarms when stagnant → niching |
| `PSO の制御設計への応用` | Control PSO | Apply PSO to PID controller tuning |

## When PSO

- Continuous, low-to-mid dim (2-50)
- Cheap-to-moderate evaluation cost
- Want simple algorithm (10 lines of code)
- Multi-modal landscape — variants like hive-split PSO find multiple
  optima

## When PSO is suboptimal

- High dim (>100): CMA-ES dominates
- Few trials (<100): BO with GP wins
- Discrete: GA or DE binary repr better
- Need theoretical guarantees: PSO has weak convergence results

## Implementation

```python
import numpy as np

def pso(fitness_fn, bounds, n_particles=30, n_iter=100,
        w=0.7, c1=1.5, c2=1.5, seed=42):
    rng = np.random.default_rng(seed)
    d = len(bounds)
    low, high = np.array(bounds).T
    x = rng.uniform(low, high, (n_particles, d))
    v = rng.uniform(-(high-low), (high-low), (n_particles, d)) * 0.1
    p_best = x.copy()
    p_best_f = np.array([fitness_fn(xi) for xi in x])
    g_best_idx = np.argmin(p_best_f)
    g_best = p_best[g_best_idx].copy()

    for it in range(n_iter):
        r1, r2 = rng.uniform(size=(2, n_particles, d))
        v = w*v + c1*r1*(p_best - x) + c2*r2*(g_best - x)
        x = np.clip(x + v, low, high)
        f = np.array([fitness_fn(xi) for xi in x])
        better = f < p_best_f
        p_best[better] = x[better]
        p_best_f[better] = f[better]
        if np.min(p_best_f) < fitness_fn(g_best):
            g_best = p_best[np.argmin(p_best_f)].copy()
    return g_best, fitness_fn(g_best)
```

## Hive-split PSO (lab variant)

When the swarm has stagnated (no improvement for K iterations),
split into 2-3 sub-swarms each with a different best particle as
attractor. Reunite if any sub-swarm finds a new global best.

This explicitly handles multi-modal landscapes that single-swarm PSO
gets stuck on. Lab paper has the Japanese implementation details.

## EM applications

- **Magnetic shield placement**: variables are (x, y) positions of N
  shields; objective is min B_leak. PSO handles non-convex.
- **Coil winding pattern**: integer #turns + continuous spacing — DE
  is better; hive-PSO competitive.
- **PID for inverter control**: lab paper PSO for control design.
"""


CMA_ES = r"""
# CMA-ES (Covariance Matrix Adaptation Evolution Strategy)

Hansen-Ostermeier 1996, refined Hansen 2003. **The gold standard for
continuous black-box optimization** in 10-300 dim.

## Algorithm

Maintain a multivariate Gaussian sampling distribution
`N(m, sigma^2 * C)`:

1. Sample lambda offspring `x_k ~ N(m, sigma^2 * C)`
2. Evaluate fitness
3. Select top mu (typically mu = lambda/2)
4. Update mean m, covariance C, step-size sigma using rank-mu and
   rank-one updates
5. Repeat

The covariance adaptation lets CMA-ES learn:
- Correlation structure of the objective
- Ill-conditioned coordinate directions (rescales axes)
- Trust region radius (sigma)

## Why CMA-ES is the standard

- **Invariant** to coordinate scaling + rotation
- **No hyperparameter tuning** — defaults work for most problems
- **Linear convergence** on convex-quadratic functions
- **Robust** to noisy / non-smooth / multi-modal landscapes
- Implementations: `pycma` (canonical), `cma` (Hansen's, on PyPI),
  built into Optuna (`optuna.samplers.CmaEsSampler`)

## Implementation via pycma

```python
import cma

es = cma.CMAEvolutionStrategy(
    x0=[0.5]*10,
    sigma0=0.1,
    inopts={"bounds": [[0]*10, [1]*10],
            "popsize": 20,
            "maxiter": 200,
            "seed": 42},
)
es.optimize(my_objective)
print(es.result.xbest, es.result.fbest)
```

## Via Optuna

```python
import optuna
study = optuna.create_study(
    sampler=optuna.samplers.CmaEsSampler(seed=42))
study.optimize(my_objective, n_trials=5000)
```

## Variants

- **BIPOP-CMA-ES**: restart with alternating small/large populations
  for multi-modal problems
- **IPOP-CMA-ES**: restart with doubled population
- **SEP-CMA-ES**: separable covariance (D variables, no off-diagonal)
  — scales to D > 1000
- **(1+1)-CMA-ES**: lightweight, mu=1, very small populations

## When CMA-ES fails

- Discrete variables (cannot diff covariance)
- Mixed types (categorical + continuous)
- Function evaluations are expensive (~$10) and budget is small
  (<100) → use BO
- Very high D (>1000) → use sep-CMA-ES or LM-CMA

## Cross-links

- official `optuna/optuna-mcp` / plain Optuna API — CmaEsSampler config
  MCTS leaves for inductor geometry optimization
"""


IMMUNE_NSGA = r"""
# Artificial Immune System + NSGA-II

## Artificial Immune System (AIS)

Inspired by biological immune system: B-cells (antibodies) bind to
antigens (objective). Clonal selection: high-affinity antibodies
clone and mutate; low-affinity die.

Lab papers (W:/.../02_最適化_進化計算/):
- `免疫アルゴリズムによる多峰性関数最適化.pdf`
- `改良型免疫アルゴリズムによる多峰性関数の大域的最大値の探索.pdf`

Both target multi-modal optimization explicitly.

## AIS algorithm (CLONALG)

```
Initialize population of N antibodies (random)
For each generation:
  Evaluate affinity (= -loss) of each
  Select top n antibodies by affinity
  Clone each selected antibody (more clones for higher affinity)
  Mutate clones (more mutation for lower affinity)
  Replace lowest-affinity antibodies with random new ones (diversity)
```

The diversity-preserving step makes AIS find MULTIPLE local optima
simultaneously — useful for design problems where you want a portfolio
of solutions, not just the best one.

## When AIS

- Multi-modal landscapes with > 5 local optima you want to find
- Diversity matters (multiple distinct designs)
- Continuous parameters

Otherwise CMA-ES + restart is usually simpler.

## NSGA-II (Multi-objective EA)

Deb-Pratap-Agarwal-Meyarivan 2002, IEEE T-EvolComp 6(2). The industry
standard for multi-objective EA.

Steps:
1. Combine parent + offspring populations (2N)
2. Non-dominated sorting → Pareto rank
3. Within same rank, **crowding distance** → prefer less crowded
4. Top N survive

Used in:
- Optuna default for `directions=[...]`
- Hayaho MCTS at the leaf NSGA-II
- Most multi-objective EM design papers

## NSGA-III for > 3 objectives

NSGA-II crowding distance breaks down at 4+ objectives (everything
becomes non-dominated). NSGA-III uses **reference points** in
objective space to maintain diversity.

```python
import optuna
sampler = optuna.samplers.NSGAIIISampler(
    reference_points=...)  # auto-generated grid in objective space
```

## Crossover operators (for real values)

| Operator | Property |
|----------|----------|
| **SBX** (Simulated Binary Crossover) | Default for Optuna NSGA-II |
| **UNDX** (Unimodal Normal Distribution) | Lab default (Hayaho MCTS) — smoother Pareto |
| **BLX-alpha** | Simple, decent |
| **PCX** (Parent-Centric Crossover) | Multi-parent, very smooth |

For Pareto fronts that need smooth coverage, UNDX > SBX > Uniform.

## Cross-links

- official `optuna/optuna-mcp` / plain Optuna API — NSGA-II / NSGA-III in Optuna
"""


LIBRARIES = r"""
# Evolutionary computation library landscape

| Library | Algorithms | Best for |
|---------|-----------|----------|
| **DEAP** | GA, DE, ES, NSGA-II, custom | Research, custom operators |
| **pymoo** | NSGA-II/III, MOEA/D, R-NSGA-II, etc. | Multi-objective, benchmarks |
| **pycma** | CMA-ES, BIPOP, IPOP | Canonical CMA-ES |
| **scipy.optimize** | DE | Quick DE without extra deps |
| **pyade** | LSHADE, jDE, SHADE | Advanced DE variants |
| **pyswarms** | PSO + variants | PSO-only, well-tested |
| **Optuna** | TPE, NSGA-II/III, CMA-ES | Production BBO interface |
| **inspyred** | GA, ES, PSO, ACO | Old but stable |
| **platypus** | NSGA-II/III, SPEA2, OMOPSO | MO research |

## DEAP skeleton (GA)

```python
import random
from deap import base, creator, tools, algorithms

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                  toolbox.attr_bool, n=20)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", lambda ind: (sum(ind),))
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

pop = toolbox.population(n=50)
algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=40)
best = tools.selBest(pop, k=1)[0]
print(best, best.fitness.values)
```

## pymoo skeleton (NSGA-II)

```python
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.problems import get_problem
from pymoo.optimize import minimize

problem = get_problem("zdt1")
algorithm = NSGA2(pop_size=100)
res = minimize(problem, algorithm, ("n_gen", 200), seed=42)
print(res.F.shape)  # (population, 2) Pareto front
```

## Optuna (recommended for lab production)

```python
import optuna
# Single-obj continuous: CMA-ES
study = optuna.create_study(
    sampler=optuna.samplers.CmaEsSampler(seed=42))
# Multi-obj: NSGA-II default
study = optuna.create_study(
    directions=["minimize", "minimize"],
    sampler=optuna.samplers.NSGAIISampler())
```

## When to use raw library vs Optuna

- **Optuna**: storage / parallelization / dashboard / standard
  workflow — usual choice
- **Raw library**: custom operators (graph repr for topology),
  custom selection (immune algorithm), research-novel variants

## Cross-references

- official `optuna/optuna-mcp` — full Optuna MCP workflow
  at leaves
"""


def get_evolutionary_documentation(topic: str = "all") -> str:
    """Dispatch evolutionary computation topics.

    Topics:
        overview      - EA family tree + when to use which
        ga_de         - Genetic Algorithm + Differential Evolution
        pso           - Particle Swarm Optimization + lab variants
        cma_es        - CMA-ES (gold standard for continuous BBO)
        immune_nsga   - Artificial Immune System + NSGA-II/III
        libraries     - DEAP/pymoo/pycma/scipy/pyade/Optuna
        all           - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("ga_de", "ga", "de", "differential_evolution",
                 "genetic_algorithm"):
        return GA_DE
    if topic in ("pso", "particle_swarm"):
        return PSO
    if topic in ("cma_es", "cma", "cmaes"):
        return CMA_ES
    if topic in ("immune_nsga", "nsga", "immune", "ais", "nsga2", "nsga3"):
        return IMMUNE_NSGA
    if topic in ("libraries", "deap", "pymoo", "library"):
        return LIBRARIES
    if topic == "all":
        return "\n\n".join([OVERVIEW, GA_DE, PSO, CMA_ES, IMMUNE_NSGA,
                            LIBRARIES])
    return (f"Unknown topic '{topic}'. Available: overview, ga_de, "
            "pso, cma_es, immune_nsga, libraries, all.")
