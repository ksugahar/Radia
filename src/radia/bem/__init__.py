"""Radia in-tree BEM module.

HACApK-backed boundary integral solvers for Radia's SIBC / induction
heating / accelerator workflows.  Replaces the slow NGSolve.bem
Laplace path (which is dense O(N^2) at typical IH wp size, see
2026-04-29 investigation).

Phase 1 (in progress): flat-triangle Galerkin BEM with custom Laplace
kernel, dense first then HACApK ACA.  Numerical anchor is NGSolve.bem
on the same mesh.
"""
