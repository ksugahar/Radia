"""Cubit-side Python helpers for cubit-mesh-export.

These modules are intended to be invoked from inside Coreform Cubit
(via the embedded Python 3.10 interpreter) by the Cubit plugin's
APREPRO commands.  They are radia-neutral: they depend only on
``import cubit`` (provided by Cubit), the Python standard library,
and optionally ``netgen.occ`` / ``ngsolve`` for the OCC entry points.

When the cubit plugin is deployed via ``cubit-plugin-install``, these
files are copied alongside the .ccm/.pyd into Cubit's
``bin/plugins/cubit_helpers/`` directory.  The .ccm command code adds
that directory to Cubit's embedded ``sys.path`` so that
``from add_kelvin import add_kelvin_cubit`` works inside any .jou
file or panel script.

System Python (3.12) callers should use the fully qualified path::

    from cubit_mesh_export.cubit_helpers.add_kelvin import (
        add_kelvin_cubit, sym_sideset_name, parse_sym_label, ...)
"""
