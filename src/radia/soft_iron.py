"""radia.SoftIron -- the unified, intent-based soft-iron object.

"Place soft iron, solve it, read the field" -- backend-agnostic over the two demag methods
(moment-yano MSC and the FEEC HDiv-VIM).  Geometry comes from a ``.vol`` file (the canonical,
correctly-oriented netgen interchange) or an in-memory NGSolve mesh; the Radia element
representation (ObjHexahedron / ...) is an INTERNAL detail the user no longer touches.

Both demag backends read the SAME mesh -> element container -> ``rad.Fld`` (analytic open
boundary); the backend is a SELECTOR (``solve(backend=...)`` / ``radia.set_demag_backend``), not
a second API.  This is the user-facing layer of the 2-layer API (see CLAUDE.md "Reduce Proprietary
API Surface").

Example::

    import radia as rad
    iron = rad.SoftIron("yoke.vol", mu_r=1000)          # or bh_table=[[H,B],...]
    coil = rad.ObjCnt(my_coilbuilder.to_radia())        # a mesh-free Biot-Savart source
    iron.solve(source=coil, backend="auto")             # auto: mesh-backed -> HDiv-VIM
    B = iron.field("b", [[0, 0, 0.05]])                 # total (iron + source) field, exact open bdry
    iron.solve(source=coil, backend="yano")             # same object, yano-MSC instead
"""
import os

import radia as rad


class SoftIron:
    """A soft-iron body for either demag backend, built from a ``.vol`` file or an NGSolve mesh.

    Parameters
    ----------
    geometry : str | os.PathLike | ngsolve.Mesh
        A netgen ``.vol`` path (recommended -- netgen owns the orientation) or an NGSolve ``Mesh``.
    mu_r : float, optional
        Linear relative permeability.  Give exactly one of ``mu_r`` / ``bh_table``.
    bh_table : list[[H, B]], optional
        Nonlinear B-H curve.
    """

    def __init__(self, geometry, mu_r=None, bh_table=None, material_filter=None, verbose=False):
        from radia.vim import soft_iron_from_mesh, soft_iron_from_vol
        if isinstance(geometry, (str, os.PathLike)):
            self.container = soft_iron_from_vol(geometry, mu_r=mu_r, bh_table=bh_table,
                                                material_filter=material_filter, verbose=verbose)
            self._geometry = os.fspath(geometry)
        else:
            self.container = soft_iron_from_mesh(geometry, mu_r=mu_r, bh_table=bh_table,
                                                 material_filter=material_filter, verbose=verbose)
            self._geometry = "<ngsolve.Mesh>"
        self.mu_r = mu_r
        self.bh_table = bh_table
        self.result = None          # last solve() return (HDiv dict or C++ yano tuple)
        self._source = []           # last applied-field source members (for total-field queries)

    def solve(self, source=None, backend="auto", prec=1e-6, maxiter=2000, method=0, image=None):
        """Solve the demagnetization in the applied field ``source``.

        Parameters
        ----------
        source : Radia handle | sequence of handles, optional
            The applied-field object(s) -- a coil container, ``rad.ObjBckg``, a permanent magnet, ...
            Combined with the iron for the solve.  ``None`` = self-demag only (M stays 0 without a source).
        backend : {"auto", "yano", "hdiv"}
            ``"auto"`` -> a mesh-backed iron dispatches to the HDiv-VIM; ``"yano"`` / ``"hdiv"`` force
            the method.  (Per-call; the global default is restored afterwards.)
        image : str, optional
            IMA mirror symmetry (e.g. ``"+x-z"``), passed to both backends.
        """
        self._source = (list(source) if isinstance(source, (list, tuple))
                        else ([source] if source is not None else []))
        rad.set_demag_backend(backend)
        try:
            top = rad.ObjCnt([self.container] + self._source) if self._source else self.container
            kw = {} if image is None else {"image": image}
            self.result = rad.Solve(top, prec, maxiter, method, **kw)
            return self.result
        finally:
            rad.set_demag_backend("auto")

    def field(self, kind, points, include_source=True):
        """Field of the solved system (analytic open boundary).  ``kind`` = 'b'|'h'|'a'|'m'.
        By default includes the applied source from the last ``solve`` (the TOTAL field)."""
        obj = (rad.ObjCnt([self.container] + self._source)
               if include_source and self._source else self.container)
        return rad.Fld(obj, kind, points)

    def magnetization(self):
        """Per-element magnetization ``[(center, [Mx,My,Mz]), ...]`` via ``rad.ObjM``."""
        return rad.ObjM(self.container)

    def __repr__(self):
        mat = f"mu_r={self.mu_r}" if self.mu_r is not None else "bh_table=<%d pts>" % len(self.bh_table or [])
        return f"SoftIron({self._geometry!r}, {mat})"
