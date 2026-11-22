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


def stable_pmis_rank(coordinates, *, fine_axis=None):
    """Return fixed PMIS ranks in the input H1 node order.

    Pass the returned list as pmis_rank to the real CPU AMS constructor.
    For an axis-aligned stretched mesh, specify its finely subdivided axis as
    0, 1, or 2. Use the default XYZ order for an unstructured mesh.
    """
    import numpy as np

    points = np.asarray(coordinates, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("coordinates must be a finite array of shape (n, 3)")
    if fine_axis is None:
        axes = (0, 1, 2)
    elif isinstance(fine_axis, (int, np.integer)) and fine_axis in (0, 1, 2):
        axes = ((fine_axis + 2) % 3, fine_axis, (fine_axis + 1) % 3)
    else:
        raise ValueError("fine_axis must be 0, 1, 2, or None")
    order = np.lexsort(tuple(points[:, i] for i in reversed(axes)))
    if np.any(np.all(points[order[1:]] == points[order[:-1]], axis=1)):
        raise ValueError("stable PMIS ranks require distinct coordinates")
    order = order[np.random.default_rng(6641).permutation(len(order))]
    rank = np.empty(len(order), dtype=np.intp)
    rank[order] = np.arange(len(order))
    return rank.tolist()
