"""radia.SoftIron -- the unified, intent-based soft-iron object.

"Place soft iron, solve it, read the field" on the FEEC HDiv-VIM route.  Geometry comes from a
``.vol`` file (the canonical, correctly-oriented netgen interchange) or an in-memory NGSolve mesh; the
Radia element representation (ObjHexahedron / ...) is an INTERNAL detail the user no longer touches.

This is the user-facing layer of the 2-layer API (see CLAUDE.md "Reduce Proprietary API Surface").

Example::

    import radia as rad
    iron = rad.SoftIron("yoke.vol", mu_r=1000)          # or bh_table=[[H,B],...]
    coil = rad.ObjCnt(my_coilbuilder.to_radia())        # a mesh-free Biot-Savart source
    iron.solve(source=coil)                             # mesh-backed -> HDiv-VIM
    B = iron.field("b", [[0, 0, 0.05]])                 # total (iron + source) field, exact open bdry
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
    order : int, optional
        HDiv finite-element order, 1 (BDM1) or 2 (BDM2).
    """

    def __init__(self, geometry, mu_r=None, bh_table=None, material_filter=None, verbose=False,
                 order=1):
        from radia.vim import MeshSoftIron, VolSoftIron
        if isinstance(geometry, (str, os.PathLike)):
            self.container = VolSoftIron(geometry, mu_r=mu_r, bh_table=bh_table,
                                         material_filter=material_filter, verbose=verbose, order=order)
            self._geometry = os.fspath(geometry)
        else:
            self.container = MeshSoftIron(geometry, mu_r=mu_r, bh_table=bh_table,
                                          material_filter=material_filter, verbose=verbose, order=order)
            self._geometry = "<ngsolve.Mesh>"
        self.mu_r = mu_r
        self.bh_table = bh_table
        self.order = int(order)
        self.result = None          # last solve() return (HDiv dict)
        self._source = []           # last applied-field source members (for total-field queries)

    def solve(self, source=None, backend="auto", prec=1e-6, maxiter=2000, method=0, image=None):
        """Solve the demagnetization in the applied field ``source``.

        Parameters
        ----------
        source : Radia handle | sequence of handles, optional
            The applied-field object(s) -- a coil container, ``rad.ObjBckg``, a permanent magnet, ...
            Combined with the iron for the solve.  ``None`` = self-demag only (M stays 0 without a source).
        backend : {"auto", "hdiv"}
            ``"auto"`` dispatches to the HDiv-VIM for supported mesh-backed soft iron;
            ``"hdiv"`` is accepted for explicitness.  (Per-call; the global default is restored afterwards.)
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
        return f"SoftIron({self._geometry!r}, {mat}, order={self.order})"
