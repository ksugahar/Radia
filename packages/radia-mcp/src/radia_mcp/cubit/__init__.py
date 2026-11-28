"""Cubit MCP server, the curved ``.vol`` export route and the exporter toolbar.

The export worker (:mod:`radia_mcp.cubit.vol_export`) runs in an isolated
external CPython: it loads the fork NGSolve, Cubit's Python API and the
native ``export vol`` command plugin built from
``packages/radia-mcp/native/cubit_to_vol_export``.

The plugin is never installed into Cubit and there is no Cubit menu or
toolbar: the external export worker loads it through
``RADIA_CUBIT_PLUGIN_DIR``.  ``python -m radia_mcp.cubit.smoke_test`` is the
export smoke test.  Kelvin helpers live in :mod:`radia_mcp.cubit.kelvin`.  This package
initializer imports nothing.
"""
