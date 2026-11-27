"""Cubit MCP server, the curved ``.vol`` export route and the exporter toolbar.

The export worker (:mod:`radia_mcp.cubit.vol_export`) runs in an isolated
external CPython: it loads the fork NGSolve, Cubit's Python API and the
native ``export netgen`` command plugin built from
``packages/radia-mcp/native/cubit_vol_export``.

``python -m radia_mcp.cubit.install`` deploys that plugin into Cubit and
registers the Netgen Export toolbar and menu (:mod:`radia_mcp.cubit.gui`, run
inside Cubit's private PySide6); ``python -m radia_mcp.cubit.smoke_test`` and
``python -m radia_mcp.cubit.toolbar_smoke`` are the export and GUI smoke
tests.  Kelvin helpers live in :mod:`radia_mcp.cubit.kelvin`.  This package
initializer imports nothing.
"""
