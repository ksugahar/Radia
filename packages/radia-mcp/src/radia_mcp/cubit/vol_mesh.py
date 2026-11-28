"""The ``export vol`` mesh in memory: Cubit to ``netgen.meshing.Mesh`` without a file.

:func:`mesh_from_cubit` runs the same export as the ``export vol`` command
(labels, Kelvin or named periodic pairs, curving with the geometric refit) on
the current Cubit model and returns the mesh instead of writing it.  The mesh
stays attached to the Cubit geometry, so :func:`curve` can curve it again at
another order, for a p-convergence study without exporting again:

    import ngsolve                     # the fork NGSolve the module was built against
    import cubit
    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    cubit.cmd('playback "model.jou"')
    from radia_mcp.cubit import vol_mesh
    ngmesh = vol_mesh.mesh_from_cubit(order=2, plugin_dir=...)
    mesh = ngsolve.Mesh(ngmesh)
    for p in (2, 3, 4):
        vol_mesh.curve(mesh, p)
        ...

Runs only in an external CPython that imported ``cubit`` (the export worker's
environment), never inside the MCP server or Cubit's embedded Python.  The
geometry lives as long as that Cubit session: :func:`curve` refuses a mesh
whose Cubit model has changed since it was built: entity ids, vertex positions,
curve lengths, surface areas and normals, volumes, node coordinates and the
block, sideset, nodeset, surface and volume labels are compared bit for bit.
A geometry edit that keeps every one of these quantities is not detected.  Do not call NGSolve's ``mesh.Curve(p)`` on it: that also queries the
Cubit geometry but skips the model check, the geometric refit (order >= 3) and
the periodic coefficient copy.  A saved ``.vol`` keeps the order it was saved at.

The compiled module ``_cubit_to_vol_mesh`` is built by
``native/cubit_to_vol_export/build.py`` next to ``cubit_to_vol_export.dll``;
``plugin_dir`` or ``$RADIA_CUBIT_PLUGIN_DIR`` names that directory.
"""
from __future__ import annotations

import importlib.machinery
import importlib.util
import os
import sys
from pathlib import Path

MODULE = "_cubit_to_vol_mesh"


def _native(plugin_dir: str | os.PathLike[str] | None = None):
    """The compiled module from the plugin directory, loaded once."""
    if MODULE in sys.modules:
        return sys.modules[MODULE]
    if "cubit" not in sys.modules:
        raise RuntimeError("import cubit and call cubit.init(...) before building a mesh from Cubit")
    from radia_mcp.cubit.export_runner import resolve_plugin_dir
    directory = resolve_plugin_dir(plugin_dir)
    # Only the module built for this interpreter; one of another Python is ignored.
    path = directory / (MODULE + importlib.machinery.EXTENSION_SUFFIXES[0])
    if not path.is_file():
        others = sorted(p.name for p in directory.glob(MODULE + "*"))
        raise RuntimeError(f"no {path.name} for this Python in {directory} (found {others}); "
                           "build it with native/cubit_to_vol_export/build.py")
    import netgen.meshing  # noqa: F401  -- the Mesh type the module returns
    spec = importlib.util.spec_from_file_location(MODULE, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[MODULE] = module
    return module


def mesh_from_cubit(order: int = 2, *, kelvin: dict | None = None, kelvin_match: float | None = None,
                    periodic: dict | None = None,
                    plugin_dir: str | os.PathLike[str] | None = None):
    """The ``export vol`` mesh of the current Cubit model, as a ``netgen.meshing.Mesh``.

    ``kelvin`` is the dict ``add_kelvin_cubit`` returned (its ``R`` and
    ``center``), with the pairing distance ``kelvin_match`` (default 1e-6 R), as
    for :func:`radia_mcp.cubit.kelvin.export_arguments`.  ``periodic`` holds the
    nine named periodic fields of :func:`radia_mcp.cubit.vol_export.periodic_arguments`.
    """
    from radia_mcp.cubit import kelvin as _kelvin
    from radia_mcp.cubit.vol_export import periodic_arguments
    if not periodic:
        periodic = None      # as for export vol: an empty request is no periodic group
    if kelvin is not None and periodic is not None:
        raise ValueError("one named periodic group per mesh: Kelvin and periodic cannot be combined")
    if kelvin is None and kelvin_match is not None:
        raise ValueError("kelvin_match needs kelvin")
    translation = match = fields = None
    if kelvin is not None:
        translation, match = _kelvin.pair_fields(kelvin, kelvin_match)
    if periodic is not None:
        periodic_arguments(periodic)  # the same field checks as export vol
        fields = {name: str(value) for name, value in periodic.items()}
    return _native(plugin_dir).build(int(order), translation, match, fields)


def curve(mesh, order: int):
    """Curve a mesh from :func:`mesh_from_cubit` again at ``order`` (1-5); returns ``mesh``.

    ``mesh`` is the ``netgen.meshing.Mesh`` or an ``ngsolve.Mesh`` made from it.
    """
    ngmesh = getattr(mesh, "ngmesh", mesh)
    _native().curve(ngmesh, int(order))
    return mesh


def is_live(mesh) -> bool:
    """True when ``mesh`` still carries the Cubit geometry it was built on."""
    return bool(_native().is_live(getattr(mesh, "ngmesh", mesh)))
