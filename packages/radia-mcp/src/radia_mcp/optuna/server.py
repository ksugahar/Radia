"""
Optuna MCP Server (radia_mcp.optuna)

Provides knowledge tools for Optuna black-box optimization in
the context of Radia / NGSolve EM analysis.

Knowledge distilled from:
  - Sano, Akiba, Imamura, Ota, Mizuno, Yanase,
    "Optuna によるブラックボックス最適化"
    (Ohmsha, 2023.2, ISBN 9784274230103) — the official Optuna team
    textbook. Lab copy: W:/03_文献・論文/04_機械学習と最適化/
    00_教科書/Optunaによるブラックボックス最適化/
  - Akiba et al., "Optuna: A Next-generation Hyperparameter
    Optimization Framework", KDD 2019.

Usage:
    mcp-server-optuna              # Start MCP server (stdio)
    mcp-server-optuna --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .usage_knowledge import get_usage_documentation
from .algorithm_knowledge import get_algorithm_documentation
from .lab_applications_knowledge import get_lab_applications_documentation

mcp = FastMCP("mcp-server-optuna")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def optuna_usage(topic: str = "all") -> str:
    """Optuna usage: basics, storage, visualization.

    Distilled from chapters 2 + 5 of the official Optuna team
    textbook (Sano-Akiba-Imamura-Ota-Mizuno-Yanase 2023, オーム社)
    plus the visualization patterns from chapter 3.

    Topics:
        "all"           - Everything
        "overview"      - When to use Optuna vs gradient-based;
                          position in radia-mcp universe; textbook ref
        "basic_usage"   - objective(trial) pattern, suggest API,
                          conditional search space, ML/EM examples
        "storage"       - SQLite / MySQL / PostgreSQL backends,
                          load_study / load_if_exists, heartbeat
        "visualization" - Plotly: history / importance / slice /
                          Pareto / intermediate / contour, dashboard
    """
    return get_usage_documentation(topic)


@mcp.tool()
def optuna_algorithm(topic: str = "all") -> str:
    """Optuna algorithm internals: samplers, MO, constraints, pruning.

    Distilled from chapter 3 (fluent usage) + chapter 5 (internals)
    of the textbook. Covers all production knobs the user actually
    needs to think about.

    Topics:
        "all"             - Everything
        "samplers"        - TPE / Random / CmaEs / NSGA-II / GP
                            decision table + key knobs (multivariate,
                            startup_trials, etc.)
        "multi_objective" - directions=[...], NSGA-II default,
                            UNDXCrossover for real values,
                            ZDT/Binh-Korn benchmarks
        "constraints"     - constraints_func + Deb domination rule
        "pruning"         - trial.report + should_prune, MedianPruner /
                            Hyperband / Patient, EM Karl-iteration
                            example
        "warm_start"      - enqueue_trial (params) vs add_trial
                            (params + value), transfer learning
        "parallelization" - In-process threads / multi-process via
                            RDB / distributed; heartbeat + retry
        "internals"       - DB schema, trial state machine, TPE in
                            one paragraph, NSGA-II in 4 steps
    """
    return get_algorithm_documentation(topic)


@mcp.tool()
def optuna_lab_applications(topic: str = "all") -> str:
    """Lab applications: how Optuna plugs into Radia / NGSolve work.

    Covers when to use Optuna vs gradient-based topology optimization
    (`radia_mcp.topology_optimization`), and 5 concrete recipes
    spanning IH / motor / inverse / WPT.

    Topics:
        "all"            - Everything
        "overview"       - Where Optuna fits in lab workflows;
                           gradient-vs-BBO decision table
        "coil_design"    - IH coil multi-objective (P_coil vs P_wp)
                           with PEEC + FEM-Kelvin coupling
        "motor_topology" - 2-level: Optuna outer (PMSM/SynRM/IM/SRM)
                           + SIMP inner (rotor density)
        "inverse"        - BH-curve / Play model parameter
                           identification from measured loops
        "wpt"            - WPT compensation topology (SS/SP/LCC/LCL)
                           + L/C tuning with conditional search space
        "literature"     - Finding similar BBO benchmarks in the lab
                           corpus via mcp-server-literature-index
    """
    return get_lab_applications_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_optuna_sweep(problem_kind: str) -> str:
    """Set up a new Optuna sweep for an EM design problem."""
    guidance = {
        "single_objective": (
            "Single-objective optimization:\n"
            "- sampler: TPESampler(multivariate=True, group=True)\n"
            "- start with `n_startup_trials=10` random\n"
            "- consider QMCSampler (Sobol) for the initial warmup\n"
            "- if continuous and dim>20: switch to CmaEsSampler\n"
        ),
        "multi_objective": (
            "Multi-objective optimization:\n"
            "- directions=['minimize', 'minimize'] (note: PLURAL)\n"
            "- sampler: NSGAIISampler(population_size=50, "
            "crossover=UNDXCrossover())\n"
            "- if continuous-only: UNDXCrossover\n"
            "- visualize via plot_pareto_front\n"
            "- DO NOT use a weighted sum -- defeats the purpose\n"
        ),
        "constrained": (
            "Constrained optimization:\n"
            "- For NSGA-II MO: use constraints_func + "
            "trial.set_user_attr('constraints', [c1, c2, ...])\n"
            "- For TPE single-obj: penalty in objective, OR\n"
            "  raise optuna.TrialPruned() if infeasibility detected\n"
            "- Convention: ci <= 0 means feasible\n"
        ),
        "expensive": (
            "Very expensive evaluations (hours per trial):\n"
            "- Use BoTorchSampler (GP surrogate) if dim < 20\n"
            "- Use multivariate TPE with n_startup_trials=20-50 "
            "for dim 20-100\n"
            "- ENABLE heartbeat + RetryFailedTrialCallback\n"
            "- Store in MySQL/PostgreSQL for parallel workers\n"
            "- Consider pruning if intermediate values exist\n"
        ),
        "topology": (
            "Topology selection + per-topology sizing:\n"
            "- Use conditional search space\n"
            "  (suggest_categorical('topo', [...]) + if/else branches)\n"
            "- Prefix branched params: rf_max_depth, gb_max_depth\n"
            "- group=True on TPESampler so branches don't mix\n"
            "- Consider 2-level: Optuna outer + SIMP inner\n"
        ),
    }

    pk = problem_kind.lower().strip()
    specific = guidance.get(pk, "")
    if not specific:
        types_list = ", ".join(guidance.keys())
        specific = (
            f"Unknown problem kind: '{problem_kind}'.\n"
            f"Known kinds: {types_list}\n"
            "Proceeding with general Optuna sweep guidance.\n"
        )

    return (
        f"Set up a new Optuna sweep for: {problem_kind}\n\n"
        f"{specific}\n"
        "Universal recipe:\n"
        "1. Define `objective(trial)` returning scalar (or tuple)\n"
        "2. Use `trial.suggest_float/int/categorical(...)` inside\n"
        "3. `create_study(study_name='...', storage='sqlite:///x.db',\n"
        "                  load_if_exists=True, ...)`\n"
        "4. `study.optimize(objective, n_trials=N, catch=(...))`\n"
        "5. Visualize: plot_optimization_history, plot_param_importances\n"
        "6. For long sweeps: enable heartbeat + retry\n\n"
        "Reference: optuna_usage, optuna_algorithm tools.\n"
        "Lab applications: optuna_lab_applications(topic='coil_design'/\n"
        "                  'motor_topology'/'inverse'/'wpt')\n"
    )


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-optuna',
    description='Optuna black-box optimization (Sano-Akiba-Imamura 2023 textbook)',
    subpackage='radia_mcp.optuna',
    related_servers=["bayesian-opt", "evolutionary", "mcmc"],
    optional_deps=["optuna"],
)


def main():
    if "--selftest" in sys.argv:
        print("Optuna MCP server self-test:")
        for name, fn, topics in [
            ("optuna_usage", optuna_usage,
             ["overview", "basic_usage", "storage", "visualization"]),
            ("optuna_algorithm", optuna_algorithm,
             ["samplers", "multi_objective", "constraints", "pruning",
              "warm_start", "parallelization", "internals"]),
            ("optuna_lab_applications", optuna_lab_applications,
             ["overview", "coil_design", "motor_topology",
              "inverse", "wpt", "literature"]),
        ]:
            for t in topics:
                doc = fn(t)
                if len(doc) < 200:
                    print(f"  FAIL {name}('{t}'): only {len(doc)} chars")
                    sys.exit(1)
                print(f"  {name}('{t}'): {len(doc)} chars OK")
            all_doc = fn("all")
            print(f"  {name}('all'): {len(all_doc)} chars OK")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
