"""radia_mcp.differential_forms — MCP server for differential forms in computational electromagnetism.

Knowledge base distilled from a curated computational-geometry bibliography,
including:
  - T. Needham, "Visual Differential Geometry and Forms" (Princeton U.P.,
    2021; Japanese translation, Maruzen, 2026)
  - A. Bossavit, "Computational Electromagnetism" (Academic Press, 1998)
  - D.N. Arnold, R.S. Falk, R. Winther, "Finite element exterior calculus,
    homological techniques, and applications", Acta Numerica 15 (2006)
  - H. Whitney, "Geometric Integration Theory" (Princeton U.P., 1957)
  - A. Kameari (亀有), "辺要素を巡って — 「張力」って？—", SA-11-013/RM-11-013 (2011)
  - 微分形式 (lecture notes, manuscript, Japanese)
  - 電磁気学とベクトル解析 (谷島賢二 編, 共立出版) — referenced as background

The server is read-only: every tool returns structured documentation strings.
No PDE solver, no symbolic computation. Use `radia_mcp.radia_ngsolve` for
concrete FEM/BEM examples on top of this theory.
"""
