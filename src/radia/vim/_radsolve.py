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

Element types: HDiv-VIM is TET / RT1 only (Sugahara 2026-06-29).  A TET mesh-backed iron
dispatches to the HDiv-VIM (:func:`hdiv_demag_solve`, order=1); a HEX / WEDGE mesh-backed
iron is routed by rad.Solve's 'auto' split to the collocation MMMM backend (the C++
surface-charge solve on the built ObjHexahedron / ObjWedge elements), since the HDiv-VIM
charge Gram is tet-only.  ``is_tet_registered`` makes that split; an explicit
``demag_backend='hdiv'`` on a non-tet iron fails loud (:func:`hdiv_demag_solve` raises).
The per-element write-back container is ObjTetrahedron / ObjHexahedron / ObjWedge, so tet
(HDiv-VIM) and hex/wedge (collocation MMMM) both round-trip ``ObjSetM`` / ``rad.Fld``.
"""
import radia as rad

_DEMAG_REGISTRY = {}   # iron container handle -> dict(mesh, mu_r, bh_table, handles)
_KNOWN_CONTAINER_MEMBERS = {}  # container handle -> member handles known to be safe for ObjCntStuf-free lookup


def clear_registry():
    """Drop all mesh<->container associations (call when Radia handles are invalidated)."""
    _DEMAG_REGISTRY.clear()
    _KNOWN_CONTAINER_MEMBERS.clear()


def register_container(container, members):
    """Record an ObjCnt built by the Python wrapper so HDiv dispatch can inspect it without calling
    ObjCntStuf on arbitrary handles.  ObjCntStuf segfaults on non-container handles in some Radia builds,
    so registry lookup must stay Python-side and conservative."""
    _KNOWN_CONTAINER_MEMBERS[container] = list(members)


def soft_iron_from_mesh(mesh, mu_r=None, bh_table=None, material_filter=None, verbose=False):
    """Build a soft-iron Radia container from an NGSolve ``mesh`` AND register it so
    ``rad.Solve`` can dispatch the FEEC HDiv-VIM backend.

    The returned container works with BOTH demag backends: ``'collocation_mmmm'`` (the applied
    ``MatLin`` / ``MatSatIsoTab`` canonical collocation MMMM path) and ``'hdiv'`` (the FEEC
    HDiv-VIM on the registered mesh).  Exactly one of ``mu_r`` (linear) or
    ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.

    TET, HEX, and WEDGE meshes all build a container; the demag backend then differs by
    element type: a TET iron solves with the HDiv-VIM (the default 'auto'/'hdiv', order=1),
    a HEX / WEDGE iron solves with collocation MMMM (rad.Solve's 'auto' split routes a
    non-tet iron there -- the HDiv-VIM charge Gram is tet-only).  Either way the per-element
    ObjTetrahedron / ObjHexahedron / ObjWedge round-trips ``ObjSetM`` + ``rad.Fld``.
    """
    from radia.netgen_mesh_import import netgen_mesh_to_radia
    if (mu_r is None) == (bh_table is None):
        raise ValueError("soft_iron_from_mesh: give exactly one of mu_r (linear) or bh_table (nonlinear)")
    handles = netgen_mesh_to_radia(mesh, material={'magnetization': [0.0, 0.0, 0.0]},
                                   combine=False, verbose=verbose, material_filter=material_filter,
                                   allow_hex=True, allow_wedge=True)
    # Apply the soft-iron material so the explicit collocation_mmmm backend also works on this container.
    mat = rad.MatLin(float(mu_r)) if mu_r is not None else rad.MatSatIsoTab(bh_table)
    for h in handles:
        rad.MatApl(h, mat)
    cont = rad.ObjCnt(handles)
    register_container(cont, handles)
    _DEMAG_REGISTRY[cont] = dict(mesh=mesh, mu_r=mu_r, bh_table=bh_table,
                                  handles=list(handles))
    return cont


def soft_iron_from_vol(vol_path, mu_r=None, bh_table=None, material_filter=None, verbose=False):
    """Build a soft-iron Radia container from a netgen ``.vol`` FILE -- the canonical, correctly
    oriented geometry interchange for BOTH demag backends.

    ``.vol`` is the SOLE Cubit<->NGSolve mesh interchange (Cubit ``export netgen`` / Netgen / OCC
    ``ngmesh.Save``).  Loading via NGSolve lets netgen own the mesh topology + face orientation,
    which avoids the hand-built-mesh pitfalls (e.g. inconsistent boundary-face winding that silently
    breaks the HDiv surface charge).  Both backends then read the SAME mesh: a TET iron's default/'hdiv'
    path solves on the registered mesh (FEEC HDiv-VIM, order=1); a HEX/WEDGE iron (or any iron under
    set_demag_backend('collocation_mmmm')) solves the built ObjHexahedron/Tetrahedron/Wedge elements
    (canonical collocation MMMM -- HDiv-VIM is tet-only, so the 'auto' split sends non-tet there).
    Exactly one of ``mu_r`` (linear) or ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.
    (Caller opens ``with ng.TaskManager():``.)
    """
    import ngsolve as ng
    mesh = ng.Mesh(str(vol_path))
    return soft_iron_from_mesh(mesh, mu_r=mu_r, bh_table=bh_table,
                               material_filter=material_filter, verbose=verbose)


def is_registered(top):
    """True if ``top`` (a rad.Solve object handle) IS, or CONTAINS, a soft-iron body registered via
    soft_iron_from_mesh.  Used by radia.Solve to route mesh-backed TET soft iron to the FEEC HDiv-VIM
    (RT1), route mesh-backed HEX/WEDGE soft iron to collocation MMMM via the 'auto' split, and leave
    everything else (mesh-less MSC/MMM, PM) on the C++ solve.  Read-only, never raises."""
    if top in _DEMAG_REGISTRY:
        return True
    members = _KNOWN_CONTAINER_MEMBERS.get(top, [])
    return any(m in _DEMAG_REGISTRY for m in members)


def is_tet_registered(top):
    """True if ``top``'s registered soft-iron mesh is ALL-TETRAHEDRA -- the HDiv-VIM (tet/RT1) scope.
    False for a non-tet (hex/wedge) mesh-backed iron -- rad.Solve's 'auto' split then routes it to the
    collocation MMMM backend (HDiv-VIM is tet-only) -- or if ``top`` is not registered.  Read-only, never
    raises (used by the rad.Solve wrapper to decide auto -> HDiv-VIM vs auto -> collocation MMMM)."""
    import ngsolve as ng
    iron = top if top in _DEMAG_REGISTRY else None
    if iron is None:
        irons = [m for m in _KNOWN_CONTAINER_MEMBERS.get(top, []) if m in _DEMAG_REGISTRY]
        iron = irons[0] if len(irons) == 1 else None
    if iron is None:
        return False
    try:
        mesh = _DEMAG_REGISTRY[iron]["mesh"]
        return all(len(el.vertices) == 4 for el in mesh.Elements(ng.VOL))
    except Exception:
        return False


def _find_registered_iron(top):
    """Return (iron_handle, [source_handles]) for a top-level rad.Solve object: the registered
    iron container plus the other (source) members.  ``top`` may be the iron handle itself or an
    ``ObjCnt([iron, sources...])``."""
    if top in _DEMAG_REGISTRY:
        return top, []
    members = _KNOWN_CONTAINER_MEMBERS.get(top, [])
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

    # IMA mirror symmetry is retired from HDiv-VIM.  Keep parsing the legacy image argument here only so
    # hdiv_demag_solve can fail loud with the canonical "use collocation MMMM" message.
    image = solve_kwargs.pop("image", None)
    if image is None and len(solve_args) >= 4 and solve_args[3]:
        image = solve_args[3]
    if solve_kwargs:
        raise TypeError(f"unsupported rad.Solve keyword(s) for HDiv-VIM dispatch: {sorted(solve_kwargs)}")

    iron, sources = _find_registered_iron(top)
    reg = _DEMAG_REGISTRY[iron]
    mesh = reg["mesh"]

    # applied field H_ext = the source members' H field (coils / ObjBckg), as an NGSolve CF
    if sources:
        H_ext = rad.RadiaField(rad.ObjCnt(list(sources)), 'h')
    else:
        H_ext = ngsolve.CoefficientFunction((0.0, 0.0, 0.0))

    res = hdiv_demag_solve(mesh, mu_r=reg["mu_r"], H_ext=H_ext, bh_table=reg["bh_table"], image=image)

    M = res["M"]
    handles = reg["handles"]
    if len(handles) != len(M):
        raise RuntimeError(
            f"demag_backend='hdiv': element/M count mismatch ({len(handles)} Radia handles vs "
            f"{len(M)} HDiv elements) -- mesh and registered handles are out of sync.")
    for h, m in zip(handles, M):
        rad.ObjSetM(h, [float(m[0]), float(m[1]), float(m[2])])
    return res
