"""Bridge: dispatch ``rad.Solve`` with ``demag_backend='hdiv'`` to the FEEC HDiv-VIM.

``rad.Solve`` operates on a Radia CONTAINER of polyhedra; the HDiv-VIM
(:func:`radia.vim.Solve`) operates on an NGSolve MESH.  A soft-iron
container's ``mu_r`` (applied via ``MatApl``) and its source NGSolve mesh are NOT
recoverable from a Radia container handle (Radia exposes no per-element vertex /
material getter), so this bridge keeps a registry populated at build time by
:func:`radia.vim.MeshSoftIron`:

    {iron_container -> dict(mesh, mu_r, bh_table, handles)}

``rad.Solve(cont, demag_backend='hdiv')`` then looks up every registered iron in
``cont``, builds the applied field ``H_ext`` from the remaining source members
(``rad.RadiaField(.., 'h')``), and runs ``radia.vim.Solve`` for one body or
``radia.vim.SolveCoupled`` for several.  It writes each element magnetization
back through ``ObjSetM`` and registers the persistent BDM field evaluators so
that ``rad.Fld`` / ``rad.ObjM`` reflect the coupled HDiv-VIM solution.

Element types: HDiv-VIM is BDM1/BDM2 on a pure-TET, pure-HEX, or pure-WEDGE mesh, and is
radia's soft-iron demag route.  rad.Solve's 'auto' split
(:func:`is_hdiv_eligible`) dispatches a mesh-backed TET / HEX / WEDGE iron to the HDiv-VIM
(:func:`radia.vim.Solve`, order=1 or 2), INCLUDING IMA image symmetry (the tet QuadDotRefl + the hex/wedge
reflected-block QuadBlockHex/Wedge(mask)).  On a genuinely reflection-matched full/reduced mesh, the
material solve and reconstructed field obey the roundoff contract; ordinary non-matching mesh comparisons
also contain discretization error.  A MIXED / pyramid mesh-backed iron is not yet HDiv-covered and fails loud
instead of falling back.  The per-element write-back container is ObjTetrahedron / ObjHexahedron / ObjWedge,
so all three round-trip ``ObjSetM`` /
``rad.Fld``.  If ``image=`` is used, dispatch also materializes the mirror-image polyhedra after the
HDiv solve.  The solved Radia container therefore redirects ``rad.Fld`` to the full field object of that
reduced solution, while the HDiv solve itself still runs on the reduced mesh.
"""
import radia as rad

_DEMAG_REGISTRY = {}   # iron container handle -> dict(mesh, mu_r, bh_table, handles)
_KNOWN_CONTAINER_MEMBERS = {}  # container handle -> member handles known to be safe for ObjCntStuf-free lookup
_FIELD_SOLUTIONS = {}  # solved handle -> one/more vim results plus an optional Radia source object


def clear_registry():
    """Drop all mesh<->container associations (call when Radia handles are invalidated)."""
    _DEMAG_REGISTRY.clear()
    _KNOWN_CONTAINER_MEMBERS.clear()
    _FIELD_SOLUTIONS.clear()


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
    for key in reg.get("field_solution_keys", []):
        _FIELD_SOLUTIONS.pop(key, None)
    reg["field_solution_keys"] = []
    reg["field_container"] = None
    reg["field_top_container"] = None
    register_container(iron, reg["handles"])


def field_solution_for(handle):
    """Return the solved BDM1/BDM2 field record for ``handle``, or ``None``."""
    return _FIELD_SOLUTIONS.get(handle)


def soft_iron_from_mesh(mesh, mu_r=None, bh_table=None, material_filter=None, verbose=False,
                        order=1):
    """Build a soft-iron Radia container from an NGSolve ``mesh`` AND register it so
    ``rad.Solve`` can dispatch the FEEC HDiv-VIM backend.

    The returned container is registered for the FEEC HDiv-VIM backend.  Exactly one of ``mu_r`` (linear) or
    ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.

    TET, HEX, and WEDGE meshes all build a container.  The 'auto' backend routes pure
    TET / HEX / WEDGE meshes to HDiv-VIM (order=1 or 2); mixed / pyramid meshes fail loud until
    HDiv coverage is added.  The per-element ObjTetrahedron / ObjHexahedron / ObjWedge
    round-trips ``ObjSetM`` + ``rad.Fld``.
    """
    from radia.netgen_mesh_import import netgen_mesh_to_radia
    if (mu_r is None) == (bh_table is None):
        raise ValueError("vim.MeshSoftIron: give exactly one of mu_r (linear) or bh_table (nonlinear)")
    vertices = _mesh_element_vertices(mesh, material_filter=material_filter)
    order = int(order)
    from ._capabilities import validate_hdiv_configuration
    validate_hdiv_configuration(
        mesh.dim, {len(element) for element in vertices}, order, mesh.GetCurveOrder())
    handles = netgen_mesh_to_radia(mesh, material={'magnetization': [0.0, 0.0, 0.0]},
                                   combine=False, verbose=verbose, material_filter=material_filter,
                                   allow_hex=True, allow_wedge=True)
    if len(vertices) != len(handles):
        raise RuntimeError(
            f"vim.MeshSoftIron: mesh import produced {len(handles)} Radia handles but "
            f"{len(vertices)} element vertex records; material_filter/order drift?"
        )
    # Apply the soft-iron material to keep ObjM/Fld material metadata meaningful.
    mat = rad.MatLin(float(mu_r)) if mu_r is not None else rad.MatSatIsoTab(bh_table)
    for h in handles:
        rad.MatApl(h, mat)
    cont = rad.ObjCnt(handles)
    register_container(cont, handles)
    _DEMAG_REGISTRY[cont] = dict(mesh=mesh, mu_r=mu_r, bh_table=bh_table, order=order,
                                  handles=list(handles), vertices=vertices, image_handles=[],
                                  field_solution_keys=[], field_container=None, field_top_container=None)
    return cont


def soft_iron_from_vol(vol_path, mu_r=None, bh_table=None, material_filter=None, verbose=False,
                       order=1):
    """Build a soft-iron Radia container from a netgen ``.vol`` FILE -- the canonical, correctly
    oriented geometry interchange for the HDiv-VIM demag backend.

    ``.vol`` is the SOLE Cubit<->NGSolve mesh interchange (Cubit ``export netgen`` / Netgen / OCC
    ``ngmesh.Save``).  Loading via NGSolve lets netgen own the mesh topology + face orientation,
    which avoids the hand-built-mesh pitfalls (e.g. inconsistent boundary-face winding that silently
    breaks the HDiv surface charge).  Pure TET / HEX / WEDGE irons solve on the registered mesh by
    default (FEEC HDiv-VIM, order=1 or 2); mixed / pyramid mesh-backed irons fail loud until HDiv coverage
    is added.  Exactly one of ``mu_r`` (linear) or ``bh_table`` (nonlinear ``[[H,B],...]``) must be given.
    (Caller opens ``with ng.TaskManager():``.)
    """
    import ngsolve as ng
    mesh = ng.Mesh(str(vol_path))
    return soft_iron_from_mesh(mesh, mu_r=mu_r, bh_table=bh_table,
                               material_filter=material_filter, verbose=verbose, order=order)


def is_registered(top):
    """True if ``top`` (a rad.Solve object handle) IS, or CONTAINS, a soft-iron body registered via
    ``vim.MeshSoftIron``.  Used by radia.Solve to find a mesh-backed soft iron, route the default
    'auto' split (pure TET/HEX/WEDGE -> HDiv-VIM; unsupported meshes fail loud), and leave everything else
    (mesh-less soft iron / fixed-M PM) on the C++ solve.  Read-only, never raises."""
    if top in _DEMAG_REGISTRY:
        return True
    members = _KNOWN_CONTAINER_MEMBERS.get(top, [])
    return any(m in _DEMAG_REGISTRY for m in members)


def is_hdiv_eligible(top):
    """True if every registered soft-iron mesh is pure TET, HEX, or WEDGE.

    This is the ``auto`` eligibility for the FEEC HDiv-VIM, Radia's soft-iron
    demag route.  Multiple eligible bodies dispatch through ``vim.SolveCoupled``.
    A MIXED / pyramid mesh-backed iron is NOT yet HDiv-covered, so rad.Solve's 'auto' split rejects it.
    Mesh-less surface-charge soft iron is retired in Radia.  Read-only, never raises."""
    import ngsolve as ng
    irons = ([top] if top in _DEMAG_REGISTRY else [
        member for member in _KNOWN_CONTAINER_MEMBERS.get(top, [])
        if member in _DEMAG_REGISTRY])
    if not irons:
        return False
    try:
        return all(
            {len(el.vertices) for el in _DEMAG_REGISTRY[iron]["mesh"].Elements(ng.VOL)}
            in ({4}, {8}, {6})
            for iron in irons)
    except Exception:
        return False


def registered_iron_count(top):
    """Number of HDiv-registered soft-iron bodies inside ``top``."""
    if top in _DEMAG_REGISTRY:
        return 1
    return sum(1 for m in _KNOWN_CONTAINER_MEMBERS.get(top, []) if m in _DEMAG_REGISTRY)


def _find_registered_irons(top):
    """Return (iron handles, source handles) for a top-level solve object."""
    if top in _DEMAG_REGISTRY:
        return [top], []
    members = _KNOWN_CONTAINER_MEMBERS.get(top, [])
    irons = [m for m in members if m in _DEMAG_REGISTRY]
    if irons:
        iron_set = set(irons)
        return irons, [member for member in members if member not in iron_set]
    else:
        raise NotImplementedError(
            "demag_backend='hdiv': rad.Solve received no HDiv-registered soft-iron body.  Build the "
            "iron via radia.vim.MeshSoftIron(mesh, mu_r=/bh_table=) so rad.Solve can dispatch "
            "the FEEC HDiv-VIM, or call radia.vim.Solve(mesh, ...) directly.")


def dispatch(top, *solve_args, **solve_kwargs):
    """Solve registered HDiv iron bodies and write their magnetization back.

    One body returns ``vim.Solve``'s result; several return
    ``vim.SolveCoupled``'s result.  The legacy ``solve_args`` (prec, maxiter,
    method) are not used by the HDiv path.
    """
    import ngsolve
    from . import CoupledBody, Solve, SolveCoupled

    # IMA mirror symmetry is wired for flat/Curve(2) pure-TET / pure-HEX / pure-WEDGE paths (tet
    # QuadDotRefl + hex/wedge reflected-block).  Parse the image argument (kwarg or legacy 4th positional)
    # and pass it to Solve, which folds the mirror charges (or fails loud for the still-unwired
    # curved / mixed / pyramid cases).
    image = solve_kwargs.pop("image", None)
    if image is None and len(solve_args) >= 4 and solve_args[3]:
        image = solve_args[3]
    if solve_kwargs:
        raise TypeError(f"unsupported rad.Solve keyword(s) for HDiv-VIM dispatch: {sorted(solve_kwargs)}")

    irons, sources = _find_registered_irons(top)
    registrations = [_DEMAG_REGISTRY[iron] for iron in irons]
    for iron, reg in zip(irons, registrations):
        _clear_image_field_handles(iron, reg)

    # applied field H_ext = the source members' H field (coils / ObjBckg), as an NGSolve CF
    if sources:
        H_ext = rad.RadiaField(rad.ObjCnt(list(sources)), 'h')
    else:
        H_ext = ngsolve.CoefficientFunction((0.0, 0.0, 0.0))

    if len(irons) == 1:
        reg = registrations[0]
        result = Solve(
            reg["mesh"], mu_r=reg["mu_r"], H_ext=H_ext,
            bh_table=reg["bh_table"], image=image, order=reg["order"])
        body_results = [result]
        returned = result
    else:
        solve_options = {} if image is None else {"image": image}
        body_specs = [
            CoupledBody(
                reg["mesh"], "iron-%s" % iron,
                mu_r=reg["mu_r"], bh_table=reg["bh_table"],
                order=reg["order"], solve_options=solve_options)
            for iron, reg in zip(irons, registrations)
        ]
        returned = SolveCoupled(body_specs, H_ext=H_ext)
        body_results = list(returned["bodies"])

    for iron, reg, result in zip(irons, registrations, body_results):
        M = result["M"]
        handles = reg["handles"]
        if len(handles) != len(M):
            raise RuntimeError(
                f"demag_backend='hdiv': element/M count mismatch ({len(handles)} Radia handles vs "
                f"{len(M)} HDiv elements) -- mesh and registered handles are out of sync.")
        for handle, magnetization in zip(handles, M):
            rad.ObjSetM(handle, [float(value) for value in magnetization])
    # Register the full HDiv solution for rad.Fld.  The iron handle contributes
    # only its solved magnetization; the top handle also includes its Radia
    # source objects.  IMA is evaluated from the C++ BDM1/BDM2 field itself.
        _FIELD_SOLUTIONS[iron] = {
            "result": result, "results": (result,), "source_object": None}
        reg["field_solution_keys"] = [iron]

    if len(irons) == 1 and top == irons[0]:
        keys = [irons[0]]
    else:
        source_object = None
        if sources:
            source_object = sources[0] if len(sources) == 1 else rad.ObjCnt(list(sources))
        record = {"results": tuple(body_results), "source_object": source_object}
        if len(body_results) == 1:
            record["result"] = body_results[0]
        _FIELD_SOLUTIONS[top] = record
        keys = [top]
        for reg in registrations:
            reg["field_solution_keys"].append(top)

    returned["field_contract"] = (
        "rad.Fld evaluates the persistent C++ HDiv charge field; IMA reflects the solution without "
        "piecewise-constant image objects"
    )
    returned["hdiv_body_count"] = len(body_results)
    return returned
