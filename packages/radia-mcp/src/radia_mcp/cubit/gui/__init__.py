"""Exporter toolbar and menu assets for Cubit's private PySide6 runtime.

The modules here run inside Cubit's embedded Python, loaded by path from the
startup shim that :mod:`radia_mcp.cubit.toolbar_install` writes.  This package
initializer imports nothing: importing Qt or ``radia_mcp`` code here would
load it into the installer's interpreter or into Cubit's runtime.
"""
