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
# Advanced lab BBO recipes (2026-05-25): wire Optuna onto Stage-2
# calc_*.py scripts for PMSM cogging / WPT misalignment robustness /
# shielding placement / litz strand AC R / Karl multi-fidelity pruning.
from .recipes_advanced_knowledge import get_recipes_advanced_documentation
# Kanamori et al. (2016) continuous-optimization textbook companion
# (2026-05-26): chapter-by-chapter knowledge for gradient / 2nd-order /
# duality / KKT / SVM / sparse / matrix optimization.  The THEORETICAL
# foundation for "when gradient methods apply" (and therefore when BBO
# is NOT the right choice).
from .kanamori2016_textbook_knowledge import get_kanamori2016_documentation

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


@mcp.tool()
def optuna_recipes_advanced(topic: str = "all") -> str:
    """Advanced lab BBO recipes that wire Optuna onto a Stage-2 calc_*.py.

    Five lab-specific recipes (~4-6k chars each) with runnable code
    that drives an existing Stage-2 CLI in src/radia/panels/:

    Topics:
        "all"                - Everything (~25k chars)
        "overview"           - Recipe index + when to pick which
        "pmsm_cogging"       - Magnet alpha_p + slot b_s + skew angle:
                               cogging torque + ripple multi-objective
                               (NSGA-II), drives calc_motor_transient.py
        "wpt_misalignment"   - Compensation topology + L/C tuning;
                               objective = worst-case eta across a
                               5x3 lateral/vertical offset grid;
                               MedianPruner kills bad trials early.
                               Drives calc_inductance.py --coil-solver peec.
        "shielding_layout"   - mu-metal / Cu sheet placement (1-4 sheets)
                               with conditional dim search space;
                               Pareto: |B| at sensor zone vs shield mass.
        "litz_strand_design" - n_strands x strand_d x twist_pitch with
                               industry-standard discrete n_strands,
                               cost + DC_R pre-filter to skip expensive
                               PEEC on infeasible trials. Drives
                               calc_inductance.py --coil-solver peec.
        "karl_multifidelity" - IH sweep with Karl iter intermediate_value
                               reporting + MedianPruner; bad geometries
                               die in seconds. Drives calc_fem_kelvin.py.

    Aliases: pmsm/cogging, wpt/misalignment/robustness,
    shielding/shield, litz/strand, karl/pruning_recipe/multifidelity.

    Each recipe lists: variables + search space, objective + sampler,
    paste-runnable code, trial budget, expected outcome, lab tips,
    cross-references. Complements lab_applications_knowledge.py's
    5 pattern-level recipes; this module is the production-grade
    deep dive.
    """
    return get_recipes_advanced_documentation(topic)


@mcp.tool()
def optuna_kanamori2016_textbook(topic: str = "overview") -> str:
    """Kanamori et al. (2016) continuous-optimization textbook companion.

    Curated chapter-by-chapter summary of "機械学習のための連続最適化"
    (Kanamori-Suzuki-Takeuchi-Sato, 講談社サイエンティフィク MLP Series,
    2016, ISBN 978-4-06-152920-8).  Source PDF EasyOCR'd at 300 dpi
    (354 pages); knowledge module is hand-curated, not raw OCR.

    The textbook covers gradient / 2nd-order / duality / KKT / SVM /
    sparse / matrix optimization -- the THEORETICAL FOUNDATION for
    "when gradient methods apply" (and therefore when BBO is NOT the
    right choice).  Use this companion to decide between Optuna
    (this server's other tools) vs scipy gradient methods vs CVXPY
    convex solvers vs IPOPT.

    Topics:
        "overview"     - Chapter index + how to use the topic aliases
        "lab"          - Lab tie-ins table (textbook chapter ->
                         radia / NGSolve workflow + cross-ref MCP server)
        "all"          - Everything (~25-30 KB)
        "chapter01"    - Optimization-problem taxonomy + convergence vocab
        "chapter02"    - Math preliminaries: Taylor, convex, subdifferential
        "chapter03"    - Optimality conditions + termination criteria
        "chapter04"    - Gradient descent + Armijo/Wolfe line search
        "chapter05"    - Newton + Levenberg-Marquardt (inverse fit)
        "chapter06"    - Conjugate gradient (FR / PR / PR+)
        "chapter07"    - BFGS / L-BFGS / DFP (the lab default ~1e3-1e4 dim)
        "chapter08"    - Trust region (Cauchy / dogleg / Steihaug)
        "chapter09"    - Equality-constrained KKT + sensitivity = lambda
        "chapter10"    - Inequality KKT + constraint qualifications
        "chapter11"    - Barrier / interior-point (IPOPT, CVXPY)
        "chapter12"    - Lagrange duality + augmented-Lagrangian + ADMM
        "chapter13"    - MM / EM / IRLS (Karl iteration, hysteresis Newton)
        "chapter14"    - SVM (linear + kernel, dual problem, SMO)
        "chapter15"    - LASSO / ISTA / FISTA / ADMM (sparse learning)
        "chapter16"    - Nuclear-norm / SDP / matrix completion (MOR)

    Aliases: newton/lm/levenberg -> ch05, cg -> ch06, bfgs/lbfgs -> ch07,
    trust -> ch08, kkt/inequality -> ch10, barrier/ipm -> ch11,
    duality/lagrangian -> ch12, mm/em -> ch13, svm/smo -> ch14,
    lasso/sparse/proximal/fista -> ch15, matrix/nuclear/sdp -> ch16.
    """
    return get_kanamori2016_documentation(topic)


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
            ("optuna_recipes_advanced", optuna_recipes_advanced,
             ["overview", "pmsm_cogging", "wpt_misalignment",
              "shielding_layout", "litz_strand_design",
              "karl_multifidelity"]),
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
