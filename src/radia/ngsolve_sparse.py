"""Real-double NGSolve sparse types derived from actual storage, not ABI names."""

from functools import lru_cache


@lru_cache(maxsize=2)
def real_sparse_matrix_type(*, symmetric=False):
    """Return the scalar real-double CSR type for the requested storage.

    NGSolve's exported Python class names depend on the compiler. Build a
    one-element mass matrix once per storage kind and use its actual type.
    This does not choose an inverse or create a TaskManager region; callers
    retain ownership of numerical execution and their existing solver route.
    """
    import ngsolve as ng
    from netgen.meshing import Element2D, Mesh, MeshPoint, Pnt

    try:
        mesh = Mesh(dim=2)
        vertices = [mesh.Add(MeshPoint(Pnt(*point))) for point in
                    ((0, 0, 0), (1, 0, 0), (0, 1, 0))]
        mesh.Add(Element2D(1, vertices))
        space = ng.L2(ng.Mesh(mesh), order=0, complex=False)
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
