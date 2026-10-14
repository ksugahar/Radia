"""Compatibility shim — add_kelvin moved to cubit-mesh-export 2026-05-05.

The canonical implementation lives in
``cubit_mesh_export.cubit_helpers.add_kelvin`` (part of the
cubit-mesh-export package).  This shim re-exports the public API so
that existing callers using ``from add_kelvin import ...`` (with
``src/radia/panels`` on ``sys.path``) keep working.

Why the move:  add_kelvin is solver-neutral infrastructure used by
*every* domain panel (IH, electromagnet, PCB, ...).  It belongs in the
shared mesh-export layer, not in radia/panels/.  Per the
target-centric architecture, ``cubit-mesh-export`` provides
mesh + Kelvin + symmetry + Dirichlet-label conventions; the per-domain
``radia-*`` tools consume those primitives.

New code should prefer the canonical path::

    from cubit_mesh_export.cubit_helpers.add_kelvin import (
        add_kelvin_cubit, add_kelvin_occ, add_kelvin_2d_axisym,
        sym_sideset_name, parse_sym_label,
        auto_add_kelvin_from_current_model)

In Cubit-embedded Python (where ``import cubit_mesh_export`` is not
available), the .ccm command code prepends
``<Cubit>/bin/plugins/cubit_helpers/`` to ``sys.path`` so that
``from add_kelvin import ...`` resolves directly to the deployed
helper.
"""

from cubit_mesh_export.cubit_helpers.add_kelvin import (  # noqa: F401
    add_kelvin_cubit,
    add_kelvin_occ,
    add_kelvin_2d_axisym,
    auto_add_kelvin_from_current_model,
    sym_sideset_name,
    parse_sym_label,
    auto_offset_direction,
    compute_offset_vector,
)
