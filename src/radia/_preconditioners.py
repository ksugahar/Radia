"""Internal NGSolve registration for Radia's HYPRE-free Compact AMG.

Importing radia registers ``compactamg`` in each Python process. This is not
NGSolve's ``h1amg`` or the external ``hypre`` backend. HCurl AMS needs a
discrete gradient, coordinates and (for complex problems) a real auxiliary
matrix; use the factories in ``radia.sparsesolv_ngsolve`` for that contract.
"""

from ngsolve.comp import RegisterPreconditioner


def _make_compact_amg(mat, freedofs, flags):
    # Load on use: registration must not depend on an earlier manual import
    # of this native module, and a broken native install must fail explicitly.
    from .sparsesolv_ngsolve import CompactAMGPreconditioner

    options = flags.ToDict()
    parameters = {
        name: options[name]
        for name in ("theta", "max_levels", "min_coarse", "num_smooth", "print_level")
        if name in options
    }
    for name in ("max_levels", "min_coarse", "num_smooth", "print_level"):
        if name in parameters:
            value = parameters[name]
            if int(value) != value:
                raise ValueError(f"compactamg {name} must be an integer")
            parameters[name] = int(value)
    return CompactAMGPreconditioner(mat, freedofs=freedofs, **parameters)


# Python imports are process-local and cached: this runs once per process,
# including fresh sessions after a wheel installation or editable relocation.
# Do not alias 'hypre': its external backend has different semantics.
RegisterPreconditioner("compactamg", _make_compact_amg)
