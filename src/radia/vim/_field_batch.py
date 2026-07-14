"""Persistent C++ field evaluator for an HDiv RT1 solution.

The order-1 solution needs neither a piecewise-constant magnetization collapse
nor internal-face sources:
M is linear per tet (NGSolve HDiv order=1 on tets = full (P1)^3, BDM1-type), so

    * internal faces carry NO charge (HDiv conformity: M.n continuous),
    * boundary faces carry a LINEAR sigma = M.n,
    * each tet carries a CONSTANT volume charge rho = -div M,

and both pieces have closed forms in the C++ production kernel.  Solve-time
materialization stores the immutable C++ source evaluator in the result.  Its
NumPy-buffer API performs no per-call source packing, evaluates all IMA terms in
one TaskManager region, and selects a quadrupole treecode for sufficiently large
target-source work.  Tet leaves retain the analytic near kernel; hex/wedge and
curved leaves retain their NGSolve quadrature cloud.
"""
import time

import numpy as np
import radia._radia_pybind as _rp

MU0 = 4.0e-7 * np.pi
_TET_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
_FIELD_TREE_LEAF = 32
_FIELD_TREE_THETA = 0.05
_FIELD_TREE_MIN_SOURCES = 256
_FIELD_TREE_AUTO_MIN_WORK = 500_000_000
_FIELD_TREE_RELATIVE_TOLERANCE = 1.0e-5
_FIELD_TREE_PROBE_COUNT = 16


def _scalar(value):
    return float(value[0] if isinstance(value, tuple) else value)


# ------------------------------------------------------------------ solution pieces
def _linear_M_coefficients(gfM):
    """Exact per-element linear representation M(x) = a_el + G_el @ (x - c_el) of an
    HDiv order-1 (BDM1 = full P1 per tet) GridFunction, from element-wise moments.

    Uses centroid-centered moments so the per-element normal system is the
    well-conditioned block-diagonal [[V, 0], [0, S_c]].  Exact for any per-element
    linear vector field (locked against direct evaluation in the tests).
    Returns (a (n_el,3), G (n_el,3,3), c (n_el,3), V (n_el,))."""
    import ngsolve as ng
    mesh = gfM.space.mesh
    one = ng.CoefficientFunction(1.0)
    xyz = (ng.x, ng.y, ng.z)

    def ew(cf):
        return np.asarray(ng.Integrate(cf, mesh, ng.VOL, element_wise=True), float)

    V = ew(one)
    mom1 = np.stack([ew(xyz[j]) for j in range(3)], axis=1)            # INT x_j
    c = mom1 / V[:, None]
    S = np.empty((len(V), 3, 3))
    for j in range(3):
        for k in range(j, 3):
            S[:, j, k] = S[:, k, j] = ew(xyz[j] * xyz[k])              # INT x_j x_k
    Sc = S - np.einsum("e,ej,ek->ejk", V, c, c)                        # centered
    m0 = np.stack([ew(gfM[i]) for i in range(3)], axis=1)              # INT M_i
    m1 = np.empty((len(V), 3, 3))
    for i in range(3):
        for j in range(3):
            m1[:, i, j] = ew(gfM[i] * xyz[j])                          # INT M_i x_j
    a = m0 / V[:, None]                                                # M at centroid
    rhs = m1 - np.einsum("ei,ej->eij", m0, c)                          # INT M_i (x_j - c_j)
    # rhs[e,i,:] = Sc[e] @ G[e,i,:]  (Sc symmetric)  ->  batch-solve per component
    G = np.linalg.solve(Sc, rhs.transpose(0, 2, 1)).transpose(0, 2, 1)
    return a, G, c, V


def _tet_boundary_charges(mesh, a, G, c):
    """Return boundary triangles and their linear ``sigma=M.n`` coefficients."""
    import ngsolve as ng
    pts = np.array([v.point for v in mesh.vertices], float)
    faces = {}
    for e, el in enumerate(mesh.Elements(ng.VOL)):
        vid = [v.nr for v in el.vertices]
        Vv = pts[vid]
        cen = Vv.mean(axis=0)
        for tri in _TET_FACES:
            tid = tuple(sorted(vid[t] for t in tri))
            Pf = Vv[list(tri)]
            nf = np.cross(Pf[1] - Pf[0], Pf[2] - Pf[0])
            nf = nf / np.linalg.norm(nf)
            if np.dot(nf, cen - Pf[0]) > 0:
                nf = -nf                                               # outward of el
            rec = faces.get(tid)
            if rec is None:
                faces[tid] = [Pf, e, nf]
            else:
                rec[1] = None                                          # internal
    bnd_P, bnd_sigma0, bnd_s = [], [], []
    for Pf, owner, nf in faces.values():
        if owner is not None:                                          # boundary face
            # sigma(r') = M(r').n = (a + G(r'-c)).n  ->  sigma0 + s.r'
            svec = G[owner].T @ nf
            sig0 = float(np.dot(a[owner], nf) - np.dot(svec, c[owner]))
            bnd_P.append(Pf)
            bnd_sigma0.append(sig0)
            bnd_s.append(svec)
    return bnd_P, bnd_sigma0, bnd_s


def _image_arrays(res):
    image = res.get("image")
    if image is None:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=float)
    from ._image import image_group, parse_image_string
    terms = image_group(parse_image_string(image))
    masks = np.asarray([sum(1 << axis for axis in axes) for axes, _ in terms], dtype=np.int32)
    signs = np.asarray([sign for _, sign in terms], dtype=float)
    return masks, signs


def _materialize_field_evaluator(res):
    """Build and cache the immutable C++ source evaluator exactly once."""
    if not isinstance(res, dict):
        raise TypeError("vim field evaluator requires Solve's result dict")
    gfM = res.get("gfM")
    if gfM is None:
        raise ValueError(
            "vim.FieldFromSolution: res carries no 'gfM' GridFunction -- pass the dict "
            "returned by vim.Solve/rad.Solve unmodified.")
    if int(res.get("order", -1)) != 1:
        raise NotImplementedError(
            "vim.FieldFromSolution: wired for order=1 (RT1/BDM1) solutions only "
            "(got order=%r)." % (res.get("order"),))
    cached = res.get("_field_evaluator")
    if cached is not None:
        return cached
    import ngsolve as ng
    started = time.perf_counter()
    mesh = gfM.space.mesh
    vcounts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if vcounts not in ({4}, {6}, {8}):
        raise NotImplementedError(
            "vim.FieldFromSolution: pure TET/HEX/WEDGE only (element vertex counts: %s)."
            % (sorted(vcounts),))
    image_masks, image_signs = _image_arrays(res)

    if vcounts == {4} and res.get("curve_order") is None:
        a, G, c, _V = _linear_M_coefficients(gfM)
        bnd_P, bnd_sigma0, bnd_s = _tet_boundary_charges(mesh, a, G, c)
        xyz = np.asarray([v.point for v in mesh.vertices], dtype=float)
        elements = list(mesh.Elements(ng.VOL))
        volume = np.zeros((len(elements), 16), dtype=float)
        for e, el in enumerate(elements):
            volume[e, :12] = xyz[[v.nr for v in el.vertices]].reshape(-1)
            volume[e, 12] = -np.trace(G[e])
        surface = np.zeros((len(bnd_P), 22), dtype=float)
        for row, (vertices, sigma0, slope) in enumerate(zip(bnd_P, bnd_sigma0, bnd_s)):
            surface[row, :9] = np.asarray(vertices, float).reshape(-1)
            surface[row, 9] = sigma0
            surface[row, 10:13] = slope
        evaluator = _rp._HDivFieldEvaluator.from_tet(
            volume, surface, image_masks, image_signs,
            _FIELD_TREE_LEAF, _FIELD_TREE_THETA,
            _FIELD_TREE_MIN_SOURCES, _FIELD_TREE_AUTO_MIN_WORK,
            _FIELD_TREE_RELATIVE_TOLERANCE, _FIELD_TREE_PROBE_COUNT)
        source_kind = "analytic-tet"
    else:
        rho = -ng.div(gfM)
        normal = ng.specialcf.normal(mesh.dim)
        sigma = ng.InnerProduct(gfM.Trace(), normal)
        qpts, qweight = [], []
        intorder = 9
        for index in range(mesh.GetNE(ng.VOL)):
            eid = ng.ElementId(ng.VOL, index)
            trafo = mesh.GetTrafo(eid)
            for ip in ng.IntegrationRule(mesh[eid].type, intorder):
                mip = trafo(ip)
                qpts.append([float(mip.point[k]) for k in range(3)])
                qweight.append(float(ip.weight*mip.measure)*_scalar(rho(mip)))
        for index in range(mesh.GetNE(ng.BND)):
            eid = ng.ElementId(ng.BND, index)
            trafo = mesh.GetTrafo(eid)
            for ip in ng.IntegrationRule(mesh[eid].type, intorder):
                mip = trafo(ip)
                qpts.append([float(mip.point[k]) for k in range(3)])
                qweight.append(float(ip.weight*mip.measure)*_scalar(sigma(mip)))
        evaluator = _rp._HDivFieldEvaluator.from_cloud(
            np.ascontiguousarray(qpts, dtype=float), np.ascontiguousarray(qweight, dtype=float),
            image_masks, image_signs, _FIELD_TREE_LEAF, _FIELD_TREE_THETA,
            _FIELD_TREE_MIN_SOURCES, _FIELD_TREE_AUTO_MIN_WORK,
            _FIELD_TREE_RELATIVE_TOLERANCE, _FIELD_TREE_PROBE_COUNT)
        source_kind = "quadrature-cloud"

    stats = dict(evaluator.stats())
    stats["source_kind"] = source_kind
    stats["build_wall_s"] = time.perf_counter()-started
    res["_field_evaluator"] = evaluator
    res["field_evaluator_stats"] = stats
    res["field_evaluator_build_wall_s"] = stats["build_wall_s"]
    return evaluator


def field_from_solution(res, points, algorithm="auto"):
    """Demagnetizing field H_demag (A/m) of a solved HDiv-VIM magnetization at
    ``points`` (N,3), evaluated from the ORDER-1 solution directly -- no per-element
    constant-M collapse, hence none of the near-surface piecewise-constant ripple of
    ``rad.Fld`` on the write-back elements (the O(h) bumps measured at standoff ~
    element size disappear identically; see the module docstring).

    ``res`` is the dict returned by ``vim.Solve`` / ``rad.Solve`` on a MeshSoftIron
    (it must carry the ``gfM`` GridFunction; order=1 pure TET/HEX/WEDGE only).
    Iron contribution only:

        B outside the iron = MU0 * (H_ext + H_demag)
        B inside  the iron = MU0 * (H_ext + H_demag + M)

    ``algorithm="direct"`` is the exact discrete source sum.  The default
    ``"auto"`` uses direct evaluation for ordinary batches and considers the
    quadrupole tree only above the large-work threshold; representative points
    must satisfy the configured direct-reference tolerance and show a measured
    speed benefit before the full batch uses the tree."""
    pts = np.ascontiguousarray(np.asarray(points, float).reshape(-1, 3))
    evaluator = _materialize_field_evaluator(res)
    return np.asarray(evaluator.field(pts, str(algorithm)), float)/(4.0*np.pi)


def magnetization_from_solution(res, points):
    """Evaluate the RT1 magnetization inside the solved mesh, zero outside."""
    import ngsolve as ng

    gfM = res.get("gfM") if isinstance(res, dict) else None
    if gfM is None:
        raise ValueError("vim magnetization evaluation requires Solve's result dict")
    pts = np.asarray(points, float).reshape(-1, 3)
    values = np.zeros((len(pts), 3), dtype=float)
    mesh = gfM.space.mesh
    if not len(pts):
        return values
    # NGSolve's ndarray MeshPoint representation exposes nr=-1 for points
    # outside the volume, allowing one vectorized lookup and one batched GF
    # evaluation instead of N Python point-location/evaluation calls.
    with ng.TaskManager():
        mapped = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
        if isinstance(mapped, np.ndarray) and mapped.dtype.names and "nr" in mapped.dtype.names:
            valid = np.asarray(mapped["nr"] >= 0)
            if np.any(valid):
                values[valid] = np.asarray(gfM(mapped[valid]), float).reshape(-1, 3)
            return values
        for index, point in enumerate(pts):
            try:
                values[index] = np.asarray(gfM(mesh(*map(float, point))), float).reshape(3)
            except Exception:
                pass
    return values
