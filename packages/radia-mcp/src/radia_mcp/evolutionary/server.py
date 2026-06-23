"""
Evolutionary Computation MCP Server (radia_mcp.evolutionary)

GA / DE / PSO / CMA-ES / Immune algorithm / NSGA-II/III knowledge,
distilled from W:/.../02_最適化_進化計算/.

Usage:
    mcp-server-evolutionary              # Start MCP server
    mcp-server-evolutionary --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_evolutionary_documentation, TOPICS

mcp = FastMCP("mcp-server-evolutionary")


@mcp.tool()
def evolutionary(topic: str = "all") -> str:
    """Evolutionary computation algorithms for EM optimization.

    Topics:
        "all"          - Everything
        "overview"     - EA family tree; decision table per problem type
        "ga_de"        - Genetic Algorithm + Differential Evolution
                          (DE is the real-valued workhorse; jDE/SHADE/LSHADE)
        "pso"          - Particle Swarm Optimization + lab variants
                          (global PSO, hive-split PSO, control PSO)
        "cma_es"       - CMA-ES (gold standard for 10-300 dim continuous)
        "immune_nsga"  - Artificial Immune System (lab papers) + NSGA-II/III
                          + UNDX/SBX crossover for Pareto fronts
        "libraries"    - DEAP/pymoo/pycma/scipy/pyade/Optuna skeletons
    """
    return get_evolutionary_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-evolutionary',
    description='GA / DE / PSO / CMA-ES / Immune / NSGA-II for EM',
    subpackage='radia_mcp.evolutionary',
    related_servers=["bayesian-opt", "topology-optimization"],
    optional_deps=["deap", "pymoo", "cma"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-evolutionary',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Evolutionary MCP server self-test:")
        topics = ["overview", "ga_de", "pso", "cma_es", "immune_nsga",
                  "libraries", "all"]
        for t in topics:
            doc = evolutionary(t)
            if len(doc) < 200:
                print(f"  FAIL evolutionary('{t}'): only {len(doc)} chars")
                sys.exit(1)
            print(f"  evolutionary('{t}'): {len(doc)} chars OK")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
