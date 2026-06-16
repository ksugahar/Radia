"""Bridge: dispatch ``rad.Solve`` with ``demag_backend='hdiv'`` to the FEEC HDiv-VIM.

``rad.Solve`` operates on a Radia CONTAINER of polyhedra; the HDiv-VIM
(:func:`radia.vim.hdiv_demag_solve`) operates on an NGSolve MESH.  A soft-iron
container's ``mu_r`` (applied via ``MatApl``) and its source NGSolve mesh are NOT
recoverable from a Radia container handle (Radia exposes no per-element vertex /
material getter), so this bridge keeps a registry populated at build time by
:func:`soft_iron_from_mesh`:

    {iron_container -> dict(mesh, mu_r, bh_table, handles)}

``rad.Solve(cont, demag_backend='hdiv')`` then looks up the registered iron in
``cont``, builds the applied field ``H_ext`` from the remaining (source) members
(``rad.RadiaField(.., 'h')``), runs ``hdiv_demag_solve``, and writes the per-element
magnetization back onto the iron's Radia elements via ``ObjSetM`` so that
``rad.Fld`` / ``rad.ObjM`` reflect the HDiv-VIM solution.

Element types: TET meshes use the scalable C++ charge-Gram H-matrix; HEX/WEDGE
meshes use the dense analytic polytope charge Gram (O(N^2), correct -- verified
demag_z -> 1/3 on hex/wedge cubes).  :func:`hdiv_demag_solve` (scalable=None)
auto-selects the path from the mesh element type, so this dispatch is element-agnostic.
The per-element write-back container is ObjTetrahedron / ObjHexahedron / ObjWedge
(netgen_mesh_to_radia allow_hex=allow_wedge=True), so tet/hex/wedge all round-trip.
"""
import radia as rad

_DEMAG_REGISTRY = {}   # iron container handle -> dict(mesh, mu_r, bh_table, handles)


def clear_registry():
    """Drop all mesh<->container associations (call when Radia handles are invalidated)."""
    _DEMAG_REGISTRY.clear()


def soft_iron_from_mesh(mesh, mu_r=None, bh_table=None, material_filter=None, verbose=False):
    """Build a soft-iron Radia container from an NGSolve ``mesh`` AND register it so
    ``rad.Solve`` can dispatch the FEEC HDiv-VIM backend.

    The returned container works with BOTH demag backends: ``'yano'`` (the applied
    ``MatLin`` / ``MatSatIsoTab`` MSC -- the legacy path) and ``'hdiv'`` (the FEEC
    HDiv-VIM on the registered mesh).  Exactly one of ``mu_r`` (linear) or
    ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.

    TET, HEX, and WEDGE meshes are all supported (ObjTetrahedron / ObjHexahedron /
    ObjWedge per element for the ``ObjSetM`` write-back + ``rad.Fld``).
    """
    from radia.netgen_mesh_import import netgen_mesh_to_radia
    if (mu_r is None) == (bh_table is None):
        raise ValueError("soft_iron_from_mesh: give exactly one of mu_r (linear) or bh_table (nonlinear)")
    handles = netgen_mesh_to_radia(mesh, material={'magnetization': [0.0, 0.0, 0.0]},
                                   combine=False, verbose=verbose, material_filter=material_filter,
                                   allow_hex=True, allow_wedge=True)
    # apply the soft-iron material so the legacy yano backend ALSO works on this container
    mat = rad.MatLin(float(mu_r)) if mu_r is not None else rad.MatSatIsoTab(bh_table)
    for h in handles:
        rad.MatApl(h, mat)
    cont = rad.ObjCnt(handles)
    _DEMAG_REGISTRY[cont] = dict(mesh=mesh, mu_r=mu_r, bh_table=bh_table, handles=list(handles))
    return cont


def _find_registered_iron(top):
    """Return (iron_handle, [source_handles]) for a top-level rad.Solve object: the registered
    iron container plus the other (source) members.  ``top`` may be the iron handle itself or an
    ``ObjCnt([iron, sources...])``."""
    if top in _DEMAG_REGISTRY:
        return top, []
    try:
        members = list(rad.ObjCntStuf(top))
    except Exception:
        members = []
    irons = [m for m in members if m in _DEMAG_REGISTRY]
    if len(irons) == 1:
        iron = irons[0]
        return iron, [m for m in members if m != iron]
    if not irons:
        raise NotImplementedError(
            "demag_backend='hdiv': rad.Solve received no HDiv-registered soft-iron body.  Build the "
            "iron via radia.vim.soft_iron_from_mesh(mesh, mu_r=/bh_table=) so rad.Solve can dispatch "
            "the FEEC HDiv-VIM, or call radia.vim.hdiv_demag_solve(mesh, ...) directly.")
    raise NotImplementedError(
        "demag_backend='hdiv': multiple HDiv-registered iron bodies in one rad.Solve container is not "
        "supported yet -- solve them separately or call radia.vim.hdiv_demag_solve directly.")


def dispatch(top, *solve_args, **solve_kwargs):
    """``rad.Solve(demag_backend='hdiv')`` handler.  Solves the registered iron's demag with the
    FEEC HDiv-VIM, writes per-element M back via ``ObjSetM``, and returns the hdiv_demag_solve
    result dict (note: a richer return than the legacy C++ rad.Solve tuple).  The legacy
    ``solve_args`` (prec, maxiter, method) are not used by the HDiv path."""
    import ngsolve
    from ._solve import hdiv_demag_solve

    iron, sources = _find_registered_iron(top)
    reg = _DEMAG_REGISTRY[iron]
    mesh = reg["mesh"]

    # applied field H_ext = the source members' H field (coils / ObjBckg), as an NGSolve CF
    if sources:
        H_ext = rad.RadiaField(rad.ObjCnt(list(sources)), 'h')
    else:
        H_ext = ngsolve.CoefficientFunction((0.0, 0.0, 0.0))

    res = hdiv_demag_solve(mesh, mu_r=reg["mu_r"], H_ext=H_ext, bh_table=reg["bh_table"])

    M = res["M"]
    handles = reg["handles"]
    if len(handles) != len(M):
        raise RuntimeError(
            f"demag_backend='hdiv': element/M count mismatch ({len(handles)} Radia handles vs "
            f"{len(M)} HDiv elements) -- mesh and registered handles are out of sync.")
    for h, m in zip(handles, M):
        rad.ObjSetM(h, [float(m[0]), float(m[1]), float(m[2])])
    return res
