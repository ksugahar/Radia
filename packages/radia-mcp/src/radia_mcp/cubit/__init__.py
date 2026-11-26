"""Cubit workflows that run in an external CPython process.

The curving worker (:mod:`radia_mcp.cubit.vol_export`) imports Cubit's
Python API and Netgen in one process and writes a curved Netgen ``.vol``.
It replaces the native Cubit plugin exporter; Cubit itself loads no Radia
or Netgen code.
"""
