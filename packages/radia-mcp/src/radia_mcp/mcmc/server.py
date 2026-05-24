"""
MCMC MCP Server (radia_mcp.mcmc)

Provides knowledge tools for Markov Chain Monte Carlo, Monte Carlo
Tree Search, and Sampled Pattern Matching in the context of Radia /
NGSolve EM analysis.

Knowledge distilled from:
  - Wada-Sugino-Mokishi, "ゼロからできるMCMC" (Asakura, 2019)
    -- statistical physics textbook
  - Suyama, "ベイズ推論による機械学習入門" + "ベイズ深層学習"
    (MLP series)
  - Shadare-Sadiku-Musa 2019 (Open J Modelling and Simulation):
    direct MCMC for Laplace's equation in axisym domain
  - ★ Sato-Igarashi 2023 (IEEE T-MAG 59(5)): multi-objective MCTS
    for PM motor (Hokkaido lab focus)
  - ★ Yin-Sato-Igarashi 2024 (IEEE T-MAG 60(3)): MCTS for inductor
    design with non-linear materials + inheritance
  - Saotome-Coulomb-Saito-Sabonnadiere 1995 (IEEE T-MAG): Sampled
    Pattern Matching method for magnetic core shape design

Usage:
    mcp-server-mcmc              # Start MCP server (stdio)
    mcp-server-mcmc --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from .algorithms_knowledge import get_algorithms_documentation
from .libraries_knowledge import get_libraries_documentation
from .em_applications_knowledge import get_em_applications_documentation
from .mcts_lab_knowledge import get_mcts_lab_documentation
from .spm_knowledge import get_spm_documentation
from ..common import register_status_tool

mcp = FastMCP("mcp-server-mcmc")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def mcmc_algorithms(topic: str = "all") -> str:
    """MCMC algorithm theory: MH / Gibbs / HMC / NUTS / PT, diagnostics.

    From the canonical MCMC textbooks (Wada-Sugino-Mokishi 2019,
    Suyama ベイズ推論). Use this for the foundational mathematics
    + algorithm selection.

    Topics:
        "all"                  - Everything
        "overview"             - MCMC family tree + when to use
        "metropolis_hastings"  - Random walk MH, proposal tuning,
                                  acceptance-rate rule (0.234 in high
                                  dim), adaptive MH
        "hamiltonian_mc"       - HMC + NUTS (default in PyMC/Stan);
                                  leapfrog integrator; when NUTS
                                  fails (discrete, discontinuity,
                                  funnel)
        "gibbs"                - Conditional sampling; conjugate
                                  cycles; block / collapsed Gibbs
        "parallel_tempering"   - Replica exchange for multi-modal
                                  posteriors; temperature ladder
                                  design
        "diagnostics"          - R-hat, ESS, trace, ACF, divergences;
                                  the "did it converge?" workflow
    """
    return get_algorithms_documentation(topic)


@mcp.tool()
def mcmc_libraries(topic: str = "decision") -> str:
    """MCMC libraries: PyMC / Stan / emcee / NumPyro / BlackJAX.

    Library selection + skeletons. Defaults to "decision" (the
    decision table) because most users ask "which library should
    I use".

    Topics:
        "all"        - Everything
        "decision"   - Decision table (DEFAULT) -- which lib for
                       which problem
        "pymc"       - PyMC (Python-first, PyTensor backend)
        "stan"       - Stan / cmdstanpy (production Bayesian)
        "emcee"      - emcee (affine-invariant; works with
                       black-box FEM forward; MPI-parallel)
        "numpyro"    - NumPyro / Pyro (JAX, GPU-friendly)
        "blackjax"   - BlackJAX (research composable, JAX scan)
    """
    return get_libraries_documentation(topic)


@mcp.tool()
def mcmc_em_applications(topic: str = "all") -> str:
    """MCMC applied to EM engineering problems.

    Concrete recipes for:
    - Direct field solution (Shadare-Sadiku 2019, MCMC for Laplace)
    - Bayesian material parameter ID (BH curve, Play model)
    - Forward uncertainty quantification (PCE, Sobol)
    - Design under uncertainty (robust optimization)

    Topics:
        "all"                       - Everything
        "direct_mcmc_laplace"       - Shadare 2019 random-walk
                                       solution of Laplace BVP;
                                       when it competes with FEM
        "bayesian_inverse"          - Material parameter ID:
                                       BH curve, Play model
                                       hysteresis parameters,
                                       permeability tomography,
                                       biomagnetic source ID
        "bayesian_uq"               - Forward UQ; PCE for cheap
                                       moments; Sobol sensitivity
                                       analysis
        "design_under_uncertainty"  - Robust optimization: outer
                                       Optuna + inner MC over
                                       uncertain inputs; chance-
                                       constrained formulations
    """
    return get_em_applications_documentation(topic)


@mcp.tool()
def mcmc_mcts_lab(topic: str = "all") -> str:
    """Monte Carlo Tree Search for EM design -- lab lineage.

    Direct reference to the Hokkaido lab Sato-Igarashi MCTS+TO
    series (★ lab focus):
    - Sato 2023 IEEE T-MAG: multi-objective MCTS for PM motor
    - Yin 2024 IEEE T-MAG: MCTS for inductor + inheritance

    Plus underlying MCTS theory (AlphaGo lineage) and pointers
    for the next thing to read.

    Topics:
        "all"             - Everything
        "overview"        - MCTS algorithm (UCB1, 4 phases);
                            AlphaGo / AlphaZero / MuZero canon;
                            why MCTS for engineering design
        "hayaho_pm_motor" - Sato-Igarashi 2023 ★: MO MCTS + NSGA-II
                            + NGnet for PM motor (T_avg vs T_rip /
                            P_iron). Novel scoring = number of
                            Pareto solutions per branch. 3-day
                            wallclock for 30 iterations on Xeon.
        "yin_inductor"    - Yin 2024 ★: MCTS + CMA-ES for inductor
                            with non-linear material selection.
                            Introduces "inherited search" warm-
                            start across similar problems.
        "lineage"         - MCTS history (Kocsis-Szepesvári 2006
                            → AlphaGo 2016 → Sato 2023 → Yin 2024);
                            reading order
    """
    return get_mcts_lab_documentation(topic)


@mcp.tool()
def mcmc_spm(topic: str = "all") -> str:
    """Sampled Pattern Matching (Saotome 1995) -- cosine-similarity
    inverse problem method for EM design.

    Topics:
        "all"               - Everything
        "overview"          - Saotome 1995 SPM algorithm; cosine
                              similarity vs least-squares; why
                              source-magnitude invariance matters
        "spm_vs_ls"         - When SPM beats LS, when LS wins;
                              hybrid SPM-then-LS workflow;
                              mathematical equivalence on unit norm
        "lab_applications"  - 3 lab papers: magnetic core (Saotome
                              1995), proton therapy dose (Mizukami
                              2010), radioactive resin shape
        "maximum_entropy"   - Sister method: MaxEnt for inverse +
                              tomographic reconstruction (lab paper
                              also in 03_最適化_モンテカルロ)
    """
    return get_spm_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_mcmc_em_problem(kind: str) -> str:
    """Set up an MCMC / MCTS / SPM analysis for an EM design problem."""
    guidance = {
        "material_id": (
            "Bayesian material parameter identification:\n"
            "- 1-10 params (mu_r, sigma, BH curve) with measured data\n"
            "- Library: PyMC with @as_op for FEM, OR emcee\n"
            "- Surrogate model (radia_mcp.pinn) if FEM > 1 min per call\n"
            "- Diagnostics: R-hat < 1.01, ESS > 400 per param\n"
            "- Report 95% credible interval, not just MAP\n"
        ),
        "design_uq": (
            "Forward uncertainty quantification:\n"
            "- < 5 inputs + smooth output: PCE (chaospy) ~ 50 FEM\n"
            "- > 10 inputs: MCMC with surrogate, 1000+ samples\n"
            "- Always start with Sobol sensitivity (~ 1000 FEM)\n"
            "  to identify which inputs matter\n"
        ),
        "robust_design": (
            "Robust design (design under uncertainty):\n"
            "- Outer Optuna (radia_mcp.optuna) over design vars\n"
            "- Inner Monte Carlo over uncertain inputs\n"
            "- Objective: E[score] - 2*Std[score] (mean - 2sigma)\n"
            "- Use surrogate for inner loop (radia_mcp.pinn)\n"
            "- 10000+ total FEM solves; budget carefully\n"
        ),
        "motor_mcts": (
            "PM motor design via MCTS + TO (★ Hokkaido Sato 2023):\n"
            "- Tree levels: (#poles, current angle, PM type, #PMs)\n"
            "- Leaf: NSGA-II over NGnet weights for rotor shape\n"
            "- UCB1 score = number of Pareto solutions per branch\n"
            "- 30 iterations ~ 3 days on 200+ cores\n"
            "- Cross-link: radia_mcp.motor for FEM backend\n"
        ),
        "inductor_mcts": (
            "Inductor design via MCTS + CMA-ES (★ Yin 2024):\n"
            "- Tree levels: material, core L/H, #turns h/v\n"
            "- Leaf: CMA-ES over geometric params (3-8 dim)\n"
            "- Objective: target L, I_sat, minimize R_dc + P_loss\n"
            "- Use 'inherited search' for parameter sweep\n"
            "  (warm-start from previous problem)\n"
        ),
        "spm_design": (
            "Sampled Pattern Matching for EM shape design:\n"
            "- Target: vector x_t at sensor locations\n"
            "- Source-magnitude invariant via cosine similarity r\n"
            "- Algorithm: greedy add-element on current contour\n"
            "- 4x better than least-squares for shape-only targets\n"
            "- Use as init for Bayesian inverse if you want UQ\n"
        ),
        "direct_laplace": (
            "Direct MCMC for Laplace's equation (Shadare 2019):\n"
            "- Only competes with FEM for SINGLE-POINT evaluation\n"
            "- Whole-field: FEM dominates (sparse solve vs O(N*N_walks))\n"
            "- Use when: open boundary + single observation point\n"
            "- Otherwise stay with NGSolve\n"
        ),
    }

    pk = kind.lower().strip()
    specific = guidance.get(pk, "")
    if not specific:
        kinds_list = ", ".join(guidance.keys())
        specific = (
            f"Unknown kind: '{kind}'.\n"
            f"Known kinds: {kinds_list}\n"
            "Falling back to general MCMC guidance.\n"
        )

    return (
        f"Set up an MCMC/MCTS/SPM analysis for: {kind}\n\n"
        f"{specific}\n"
        "Reference tools:\n"
        "- mcmc_algorithms(topic) -- MH / Gibbs / HMC / NUTS theory\n"
        "- mcmc_libraries(topic)  -- PyMC / Stan / emcee / NumPyro\n"
        "- mcmc_em_applications(topic) -- recipes for EM problems\n"
        "- mcmc_mcts_lab(topic) -- Hokkaido MCTS for motor / inductor\n"
        "- mcmc_spm(topic) -- Saotome SPM for shape design\n\n"
        "Cross-MCP links:\n"
        "- radia_mcp.optuna       -- BBO point estimates\n"
        "- radia_mcp.topology_optimization -- gradient-based shape opt\n"
        "- radia_mcp.motor / .peec / .ih   -- FEM backends\n"
        "- radia_mcp.pinn         -- surrogate models for cheap forward\n"
        "- radia_mcp.magnetic_materials    -- Play / Energy hysteresis\n"
    )


# ============================================================
# Entry point
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-mcmc",
    description=("MCMC + MCTS + SPM for EM engineering. "
                  "Hokkaido Sato 2023 MCTS PM motor + Yin 2024 inductor "
                  "+ Saotome 1995 SPM."),
    subpackage="radia_mcp.mcmc",
    related_servers=["radia-optuna", "radia-bayesian-opt",
                      "radia-evolutionary", "radia-topology-optimization"],
    optional_deps=["pymc", "emcee", "numpyro", "blackjax"],
)


def main():
    if "--selftest" in sys.argv:
        print("MCMC MCP server self-test:")
        for name, fn, topics in [
            ("mcmc_algorithms", mcmc_algorithms,
             ["overview", "metropolis_hastings", "hamiltonian_mc",
              "gibbs", "parallel_tempering", "diagnostics"]),
            ("mcmc_libraries", mcmc_libraries,
             ["pymc", "stan", "emcee", "numpyro", "blackjax",
              "decision"]),
            ("mcmc_em_applications", mcmc_em_applications,
             ["direct_mcmc_laplace", "bayesian_inverse",
              "bayesian_uq", "design_under_uncertainty"]),
            ("mcmc_mcts_lab", mcmc_mcts_lab,
             ["overview", "hayaho_pm_motor", "yin_inductor",
              "lineage"]),
            ("mcmc_spm", mcmc_spm,
             ["overview", "spm_vs_ls", "lab_applications",
              "maximum_entropy"]),
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
