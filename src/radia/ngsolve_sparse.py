"""Real-double NGSolve sparse types derived from actual storage, not ABI names."""

from functools import lru_cache


@lru_cache(maxsize=2)
def real_sparse_matrix_type(*, symmetric=False):
    """Return the scalar real-double CSR type for the requested storage.

    NGSolve's exported Python class names depend on the compiler. Build a
    small mass matrix once per storage kind and use its actual type.
    This does not choose an inverse or create a TaskManager region; callers
    retain ownership of numerical execution and their existing solver route.
    """
    import ngsolve as ng
    from netgen.geom2d import unit_square

    try:
        # Use Netgen's complete public mesh builder. Inserting a raw Element2D
        # without all required region descriptors can crash in native code.
        mesh = ng.Mesh(unit_square.GenerateMesh(maxh=1))
        space = ng.L2(mesh, order=0, complex=False)
        form = ng.BilinearForm(space, symmetric=symmetric,
                              symmetric_storage=symmetric, fp32=False)
        trial, test = space.TnT()
        form += trial * test * ng.dx
        form.Assemble()
        matrix = form.mat
        cls = type(matrix)
        dtype = matrix.AsVector().FV().NumPy().dtype
        if (dtype.kind != "f" or dtype.itemsize != 8
                or not callable(getattr(cls, "CreateFromCOO", None))):
            raise TypeError(f"unsupported matrix storage: {cls!r}, dtype={dtype}")
        return cls
    except Exception as error:
        raise RuntimeError(
            "cannot resolve NGSolve real-double sparse matrix type "
            f"(symmetric={symmetric})") from error
