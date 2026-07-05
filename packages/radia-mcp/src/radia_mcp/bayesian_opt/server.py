"""
Bayesian Optimization MCP Server (radia_mcp.bayesian_opt)

GP regression, Bayesian optimization, FMQA, physics-informed GP,
and surrogate models for CAE — distilled from 57 lab files in
public-safe curated corpus

Usage:
    mcp-server-bayesian-opt              # Start MCP server (stdio)
    mcp-server-bayesian-opt --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_bayesian_opt_documentation, TOPICS

mcp = FastMCP("mcp-server-bayesian-opt")


@mcp.tool()
def bayesian_opt(topic: str = "all") -> str:
    """Bayesian optimization, GP regression, FMQA, surrogate models.

    Topics:
        "all"                    - Everything
        "gp_fundamentals"        - GP regression theory; RBF / Matern /
                                    ARD kernels; marginal log-likelihood
                                    hyperparameter optimization
        "bayesian_optimization"  - BO algorithm; EI / UCB / KG / qNEI /
                                    qNEHVI acquisition functions;
                                    multi-fidelity, multi-objective BO
        "fmqa"                   - Factorization Machine + Quantum
                                    Annealing for binary / discrete BBO;
                                    lab metamaterial design application
        "physics_informed_gp"    - PI-GP (Raissi 2017): GP that exactly
                                    satisfies a linear PDE; Pförtner
                                    2023 generalization; wave equation
        "surrogate_models"       - Surrogate model menu (RBF / GP /
                                    RF / GBM / MLP / PINN / DeepONet /
                                    FNO); CAE懇話会 2024 workshop
        "gp_libraries"           - scikit-learn / GPy / GPyTorch /
                                    GPflow / BoTorch / Ax skeletons
        "llm_bo"                 - LLM-enhanced BO (experimental)
    """
    return get_bayesian_opt_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-bayesian-opt',
    description='BO + GP regression + FMQA + surrogate models (57 lab files; ARD kernel, PI-GP, multi-fidelity)',
    subpackage='radia_mcp.bayesian_opt',
    related_servers=["evolutionary", "pinn", "topology-optimization"],
    optional_deps=["botorch", "gpytorch", "GPy"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-bayesian-opt',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Bayesian Opt MCP server self-test:")
        topics = ["gp_fundamentals", "bayesian_optimization", "fmqa",
                  "physics_informed_gp", "surrogate_models",
                  "gp_libraries", "llm_bo", "all"]
        for t in topics:
            doc = bayesian_opt(t)
            if len(doc) < 200:
                print(f"  FAIL bayesian_opt('{t}'): only {len(doc)} chars")
                sys.exit(1)
            print(f"  bayesian_opt('{t}'): {len(doc)} chars OK")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
