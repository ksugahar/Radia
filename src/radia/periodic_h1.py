"""Face-consistent H1 space for one explicit 3D periodic interface."""


def periodic_h1_single_interface(base_space):
    """Return a Periodic H1 space matching faces by all their corner IDs.

    The mesh must already have one bijective periodic vertex identification.
    Coordinates and connectivity are unchanged. Base boundary flags are retained.
    Different face orders and missing face partners raise, rather than silently
    identifying unrelated face DOFs. Apply ngsolve.Compress to the result as usual.

    Rebuild this space after loading a mesh; serialization of this custom space
    is not supported. Multi-interface, adaptive, and parallel-distributed spaces
    are outside the verified scope.
    """
    from ._radia_pybind import _periodic_h1_single_interface
    return _periodic_h1_single_interface(base_space)
