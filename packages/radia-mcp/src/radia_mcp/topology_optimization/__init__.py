"""radia_mcp.topology_optimization — MCP server for shape and topology optimization
in nonlinear magnetostatics.

Knowledge distilled from:
  - P. Gangl, U. Langer, A. Laurain, H. Meftahi, K. Sturm,
    "Shape optimization of an electric motor subject to nonlinear
    magnetostatics", arXiv:1501.04752 (2015)
  - Gangl PhD thesis, "Shape and Topology Optimization Subject to
    3D Nonlinear Magnetostatics", JKU Linz (~2017), Part I
    (sensitivity analysis, 31 pp) + Part II (numerics, 51 pp)
  - K. Sturm, "Minimax Lagrangian Approach to the Differentiability
    of Nonlinear PDE Constrained Shape Functions" (2015)

Tied to Radia / NGSolve practice via the IPM motor cost-tracking
formulation that appears in actual lab work.
"""
