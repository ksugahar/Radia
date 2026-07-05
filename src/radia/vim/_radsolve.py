"""Bridge: dispatch ``rad.Solve`` with ``demag_backend='hdiv'`` to the FEEC HDiv-VIM.

``rad.Solve`` operates on a Radia CONTAINER of polyhedra; the HDiv-VIM
(:func:`radia.vim.Solve`) operates on an NGSolve MESH.  A soft-iron
container's ``mu_r`` (applied via ``MatApl``) and its source NGSolve mesh are NOT
recoverable from a Radia container handle (Radia exposes no per-element vertex /
material getter), so this bridge keeps a registry populated at build time by
:func:`radia.vim.MeshSoftIron`:

    {iron_container -> dict(mesh, mu_r, bh_table, handles)}

``rad.Solve(cont, demag_backend='hdiv')`` then looks up the registered iron in
``cont``, builds the applied field ``H_ext`` from the remaining (source) members
(``rad.RadiaField(.., 'h')``), runs ``radia.vim.Solve``, and writes the per-element
magnetization back onto the iron's Radia elements via ``ObjSetM`` so that
``rad.Fld`` / ``rad.ObjM`` reflect the HDiv-VIM solution.

Element types: HDiv-VIM is RT1 (HDiv order 1) on a pure-TET, pure-HEX, or pure-WEDGE mesh, and is
radia's OFFICIAL soft-iron demag route (CLAUDE.md DIRECTION 2026-07-04).  rad.Solve's 'auto' split
(:func:`is_hdiv_eligible`) dispatches a mesh-backed TET / HEX / WEDGE iron to the HDiv-VIM
(:func:`radia.vim.Solve`, order=1), INCLUDING IMA image symmetry (the tet QuadDotRefl + the hex/wedge
reflected-block QuadBlockHex/Wedge(mask); validated 2026-07-04 -- a reduced 1/2,1/4,1/8 model reproduces
the full model to ~1e-4).  A MIXED / pyramid mesh-backed iron (not yet HDiv-covered), and any iron under
set_demag_backend('collocation_mmmm'), route to the collocation MMMM gated bridge.  The per-element
write-back container is ObjTetrahedron / ObjHexahedron / ObjWedge, so all three round-trip ``ObjSetM`` /
``rad.Fld``.  If ``image=`` is used, dispatch also materializes the mirror-image polyhedra after the
HDiv solve.  The solved Radia container therefore redirects ``rad.Fld`` to the full field object of that
reduced solution, while the HDiv solve itself still runs on the reduced mesh.
"""
import radia as rad

_DEMAG_REGISTRY = {}   # iron container handle -> dict(mesh, mu_r, bh_table, handles)
_KNOWN_CONTAINER_MEMBERS = {}  # container handle -> member handles known to be safe for ObjCntStuf-free lookup
_FIELD_REDIRECTS = {}  # solved reduced container/top handle -> full field container handle for rad.Fld


def clear_registry():
    """Drop all mesh<->container associations (call when Radia handles are invalidated)."""
    _DEMAG_REGISTRY.clear()
    _KNOWN_CONTAINER_MEMBERS.clear()
    _FIELD_REDIRECTS.clear()


def register_container(container, members):
    """Record an ObjCnt built by the Python wrapper so HDiv dispatch can inspect it without calling
    ObjCntStuf on arbitrary handles.  ObjCntStuf segfaults on non-container handles in some Radia builds,
    so registry lookup must stay Python-side and conservative."""
    _KNOWN_CONTAINER_MEMBERS[container] = list(members)


def _mesh_element_vertices(mesh, material_filter=None):
    """Return volume-element vertices in the same order as ``netgen_mesh_to_radia(..., combine=False)``."""
    from ngsolve import VOL

    if material_filter is None:
        allowed = None
    elif isinstance(material_filter, str):
        allowed = {material_filter}
    elif isinstance(material_filter, (list, tuple, set)):
        allowed = set(material_filter)
    else:
        raise ValueError(
            f"material_filter must be str, list, or None, got {type(material_filter)}"
        )

    vertices = []
    for el in mesh.Elements(VOL):
        if allowed is not None and el.mat not in allowed:
            continue
        vertices.append([
            [float(c) for c in mesh.vertices[v.nr].point]
            for v in el.vertices
        ])
    return vertices


def _delete_handles(handles):
    """Best-effort cleanup for explicit IMA field images from a previous solve."""
    for h in list(handles or []):
        try:
            rad.UtiDel(h)
        except Exception:
            pass


def _clear_image_field_handles(iron, reg):
    _delete_handles(reg.get("image_handles", []))
    reg["image_handles"] = []
    for key in reg.get("field_redirect_keys", []):
        _FIELD_REDIRECTS.pop(key, None)
    reg["field_redirect_keys"] = []
    reg["field_container"] = None
    reg["field_top_container"] = None
    register_container(iron, reg["handles"])


def _materialize_image_field_handles(iron, reg, image, top=None, sources=None):
    """Add explicit IMA images so ``rad.Fld(iron, ...)`` evaluates the reduced solution's full field."""
    _clear_image_field_handles(iron, reg)
    if image is None:
        return [], iron
    from radia.ima_field import add_ima_images

    images = add_ima_images(reg["handles"], reg["vertices"], image, container=None)
    field_iron = rad.ObjCnt(list(reg["handles"]) + list(images))
    reg["image_handles"] = list(images)
    reg["field_container"] = field_iron
    _FIELD_REDIRECTS[iron] = field_iron
    keys = [iron]
    field_top = field_iron
    if top is not None:
        members = [field_iron] + list(sources or [])
        field_top = rad.ObjCnt(members)
        reg["field_top_container"] = field_top
        _FIELD_REDIRECTS[top] = field_top
        keys.append(top)
    reg["field_redirect_keys"] = keys
    return images, field_iron


def field_object_for(handle):
    """Return the full-field object for a solved HDiv IMA handle, else the handle itself."""
    return _FIELD_REDIRECTS.get(handle, handle)


def soft_iron_from_mesh(mesh, mu_r=None, bh_table=None, material_filter=None, verbose=False):
    """Build a soft-iron Radia container from an NGSolve ``mesh`` AND register it so
    ``rad.Solve`` can dispatch the FEEC HDiv-VIM backend.

    The returned container works with BOTH demag backends: ``'collocation_mmmm'`` (the applied
    ``MatLin`` / ``MatSatIsoTab`` canonical collocation MMMM path) and ``'hdiv'`` (the FEEC
    HDiv-VIM on the registered mesh).  Exactly one of ``mu_r`` (linear) or
    ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.

    TET, HEX, and WEDGE meshes all build a container.  The 'auto' backend routes pure
    TET / HEX / WEDGE meshes to HDiv-VIM (order=1); mixed / pyramid meshes stay on the
    collocation MMMM bridge.  Either way the per-element
    ObjTetrahedron / ObjHexahedron / ObjWedge round-trips ``ObjSetM`` + ``rad.Fld``.
    """
    from radia.netgen_mesh_import import netgen_mesh_to_radia
    if (mu_r is None) == (bh_table is None):
        raise ValueError("vim.MeshSoftIron: give exactly one of mu_r (linear) or bh_table (nonlinear)")
    vertices = _mesh_element_vertices(mesh, material_filter=material_filter)
    handles = netgen_mesh_to_radia(mesh, material={'magnetization': [0.0, 0.0, 0.0]},
                                   combine=False, verbose=verbose, material_filter=material_filter,
                                   allow_hex=True, allow_wedge=True)
    if len(vertices) != len(handles):
        raise RuntimeError(
            f"vim.MeshSoftIron: mesh import produced {len(handles)} Radia handles but "
            f"{len(vertices)} element vertex records; material_filter/order drift?"
        )
    # Apply the soft-iron material so the explicit collocation_mmmm backend also works on this container.
    mat = rad.MatLin(float(mu_r)) if mu_r is not None else rad.MatSatIsoTab(bh_table)
    for h in handles:
        rad.MatApl(h, mat)
    cont = rad.ObjCnt(handles)
    register_container(cont, handles)
    _DEMAG_REGISTRY[cont] = dict(mesh=mesh, mu_r=mu_r, bh_table=bh_table,
                                  handles=list(handles), vertices=vertices, image_handles=[],
                                  field_redirect_keys=[], field_container=None, field_top_container=None)
    return cont


def soft_iron_from_vol(vol_path, mu_r=None, bh_table=None, material_filter=None, verbose=False):
    """Build a soft-iron Radia container from a netgen ``.vol`` FILE -- the canonical, correctly
    oriented geometry interchange for BOTH demag backends.

    ``.vol`` is the SOLE Cubit<->NGSolve mesh interchange (Cubit ``export netgen`` / Netgen / OCC
    ``ngmesh.Save``).  Loading via NGSolve lets netgen own the mesh topology + face orientation,
    which avoids the hand-built-mesh pitfalls (e.g. inconsistent boundary-face winding that silently
    breaks the HDiv surface charge).  Both backends then read the SAME mesh: pure TET / HEX / WEDGE
    irons solve on the registered mesh by default (FEEC HDiv-VIM, order=1), while
    ``demag_backend='collocation_mmmm'`` forces the built ObjHexahedron/Tetrahedron/Wedge elements through
    canonical collocation MMMM.  Mixed / pyramid mesh-backed irons remain on the collocation bridge.
    Exactly one of ``mu_r`` (linear) or ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.
    (Caller opens ``with ng.TaskManager():``.)
    """
    import ngsolve as ng
    mesh = ng.Mesh(str(vol_path))
    return soft_iron_from_mesh(mesh, mu_r=mu_r, bh_table=bh_table,
                               material_filter=material_filter, verbose=verbose)


def is_registered(top):
    """True if ``top`` (a rad.Solve object handle) IS, or CONTAINS, a soft-iron body registered via
    ``vim.MeshSoftIron``.  Used by radia.Solve to find a mesh-backed soft iron, route the default
    'auto' split (pure TET/HEX/WEDGE -> HDiv-VIM, mixed/pyramid -> collocation MMMM bridge), and leave everything else
    (mesh-less MSC/MMM, PM) on the C++ solve.  Read-only, never raises."""
    if top in _DEMAG_REGISTRY:
        return True
    members = _KNOWN_CONTAINER_MEMBERS.get(top, [])
    return any(m in _DEMAG_REGISTRY for m in members)


def is_hdiv_eligible(top):
    """True if ``top``'s registered soft-iron mesh is PURE-TET, PURE-HEX, or PURE-WEDGE -- the 'auto'
    eligibility for the FEEC HDiv-VIM, radia's OFFICIAL soft-iron demag route (CLAUDE.md DIRECTION
    2026-07-04: HDiv-VIM is the sole soft-iron method; collocation MMMM is the gated bridge -> ELF_MAGIC).
    A MIXED / pyramid mesh-backed iron is NOT yet HDiv-covered, so rad.Solve's 'auto' split routes it to
    the collocation MMMM bridge.  (Mesh-less iron is not registered, so it never reaches here and stays on
    collocation MMMM.)  Read-only, never raises."""
    import ngsolve as ng
    iron = top if top in _DEMAG_REGISTRY else None
    if iron is None:
        irons = [m for m in _KNOWN_CONTAINER_MEMBERS.get(top, []) if m in _DEMAG_REGISTRY]
        iron = irons[0] if len(irons) == 1 else None
    if iron is None:
        return False
    try:
        mesh = _DEMAG_REGISTRY[iron]["mesh"]
        vts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
        return vts in ({4}, {8}, {6})            # pure tet / hex / wedge (HDiv-VIM RT1, incl. IMA)
    except Exception:
        return False


def registered_iron_count(top):
    """Number of HDiv-registered soft-iron bodies inside ``top`` (0, 1, or >1).  rad.Solve's 'auto' split
    uses this to FAIL LOUD on a MULTI-iron container -- which is_hdiv_eligible rejects (its len!=1 guard) --
    instead of SILENTLY demoting it to collocation MMMM (a No-Fallbacks violation).  Read-only, never raises."""
    if top in _DEMAG_REGISTRY:
        return 1
    return sum(1 for m in _KNOWN_CONTAINER_MEMBERS.get(top, []) if m in _DEMAG_REGISTRY)


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
            "iron via radia.vim.MeshSoftIron(mesh, mu_r=/bh_table=) so rad.Solve can dispatch "
            "the FEEC HDiv-VIM, or call radia.vim.Solve(mesh, ...) directly.")
    raise NotImplementedError(
        "demag_backend='hdiv': multiple HDiv-registered iron bodies in one rad.Solve container is not "
        "supported yet -- solve them separately or call radia.vim.Solve directly.")


def dispatch(top, *solve_args, **solve_kwargs):
    """``rad.Solve(demag_backend='hdiv')`` handler.  Solves the registered iron's demag with the
    FEEC HDiv-VIM, writes per-element M back via ``ObjSetM``, and returns the ``radia.vim.Solve``
    result dict (note: a richer return than the legacy C++ rad.Solve tuple).  The legacy
    ``solve_args`` (prec, maxiter, method) are not used by the HDiv path."""
    import ngsolve
    from . import Solve

    # IMA mirror symmetry is WIRED in HDiv-VIM for the FLAT pure-TET / pure-HEX / pure-WEDGE paths (tet
    # QuadDotRefl + hex/wedge reflected-block).  Parse the image argument (kwarg or legacy 4th positional)
    # and pass it to Solve, which folds the mirror charges (or fails loud for the still-unwired
    # curved / mixed / pyramid cases).
    image = solve_kwargs.pop("image", None)
    if image is None and len(solve_args) >= 4 and solve_args[3]:
        image = solve_args[3]
    if solve_kwargs:
        raise TypeError(f"unsupported rad.Solve keyword(s) for HDiv-VIM dispatch: {sorted(solve_kwargs)}")

    iron, sources = _find_registered_iron(top)
    reg = _DEMAG_REGISTRY[iron]
    mesh = reg["mesh"]
    _clear_image_field_handles(iron, reg)

    # applied field H_ext = the source members' H field (coils / ObjBckg), as an NGSolve CF
    if sources:
        H_ext = rad.RadiaField(rad.ObjCnt(list(sources)), 'h')
    else:
        H_ext = ngsolve.CoefficientFunction((0.0, 0.0, 0.0))

    res = Solve(mesh, mu_r=reg["mu_r"], H_ext=H_ext, bh_table=reg["bh_table"], image=image)

    M = res["M"]
    handles = reg["handles"]
    if len(handles) != len(M):
        raise RuntimeError(
            f"demag_backend='hdiv': element/M count mismatch ({len(handles)} Radia handles vs "
            f"{len(M)} HDiv elements) -- mesh and registered handles are out of sync.")
    for h, m in zip(handles, M):
        rad.ObjSetM(h, [float(m[0]), float(m[1]), float(m[2])])
    if image is not None:
        image_handles, field_iron = _materialize_image_field_handles(iron, reg, image, top=top, sources=sources)
        res["image_field_handles"] = image_handles
        res["field_object"] = field_iron
        res["field_contract"] = "rad.Fld on the MeshSoftIron container redirects to explicit IMA mirror images of the reduced solution"
    else:
        res["field_contract"] = "rad.Fld on the MeshSoftIron container uses the solved reduced/full mesh"
    return res
