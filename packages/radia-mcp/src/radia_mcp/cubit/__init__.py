"""Cubit MCP server, check-vol and the curved ``.vol`` export route.

The export worker (:mod:`radia_mcp.cubit.vol_export`) runs in an isolated
external CPython: it loads the fork NGSolve, Cubit's Python API and the
native ``export netgen`` command plugin built from
``packages/radia-mcp/native/cubit_vol_export``.
"""
