"""C++ batch field of an HDiv RT1 tetrahedral solution.

``rad.Fld`` on a solved ``vim.MeshSoftIron`` container evaluates the write-back
elements, i.e. the RT1/BDM1 solution COLLAPSED to per-element constant M
(``ObjSetM``).  That collapse breaks the HDiv normal continuity, so spurious
internal-face charges appear and the map at a standoff comparable to the element
size carries an O(h) piecewise-constant-M ripple (measured 2026-07-13/14 on the
parallelogram-dipole edge-focusing testbed: flat-top bumps up to +4%, fringe
integral K1g biased low by 12% -- docs/clebsch_hodograph/
edge_focusing_fem_results.json `mesh_dependence_diagnosis`).

The order-1 solution itself needs neither the collapse nor internal faces:
M is linear per tet (NGSolve HDiv order=1 on tets = full (P1)^3, BDM1-type), so

    * internal faces carry NO charge (HDiv conformity: M.n continuous),
    * boundary faces carry a LINEAR sigma = M.n,
    * each tet carries a CONSTANT volume charge rho = -div M,

and both pieces have closed forms in the C++ production kernel.
``field_from_solution`` extracts the element-local RT1 polynomial once through
NGSolve, then sends packed volume and boundary charge polynomials to the C++
analytic field kernel.  The C++ entry owns its TaskManager region and performs
the observation-point loop in parallel.  IMA images are evaluated from the
same RT1 polynomial; no piecewise-constant Radia image objects are involved.
"""
import numpy as np
import radia._radia_pybind as _rp

MU0 = 4.0e-7 * np.pi
_TET_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


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


def _tet_faces_and_charges(mesh, a, G, c):
    """Enumerate unique tet faces once and split them into the two source lists:

    * volume term: per unique face, the accumulated vector weight
      w_f = SUM_adjacent-tets rho_el * n_outward(el)  (internal faces keep the
      rho JUMP automatically; rho_el = -tr(G_el) is the constant -div M);
    * surface term: boundary faces (appearing once) with their owner element.

    Returns (face_P (list of (3,3)), face_w (list of (3,))) for the volume term and
    (bnd_P, bnd_sigma0, bnd_s) for the linear-sigma surface term."""
    import ngsolve as ng
    pts = np.array([v.point for v in mesh.vertices], float)
    rho = -np.trace(G, axis1=1, axis2=2)                               # (n_el,)
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
                faces[tid] = [Pf, rho[e] * nf, e, nf]
            else:
                rec[1] = rec[1] + rho[e] * nf
                rec[2] = None                                          # internal
    face_P, face_w, bnd_P, bnd_sigma0, bnd_s = [], [], [], [], []
    for Pf, wvec, owner, nf in faces.values():
        if float(np.dot(wvec, wvec)) > 0.0:
            face_P.append(Pf)
            face_w.append(wvec)
        if owner is not None:                                          # boundary face
            # sigma(r') = M(r').n = (a + G(r'-c)).n  ->  sigma0 + s.r'
            svec = G[owner].T @ nf
            sig0 = float(np.dot(a[owner], nf) - np.dot(svec, c[owner]))
            bnd_P.append(Pf)
            bnd_sigma0.append(sig0)
            bnd_s.append(svec)
    return face_P, face_w, bnd_P, bnd_sigma0, bnd_s


def field_from_solution(res, points):
    """Demagnetizing field H_demag (A/m) of a solved HDiv-VIM magnetization at
    ``points`` (N,3), evaluated from the ORDER-1 solution directly -- no per-element
    constant-M collapse, hence none of the near-surface piecewise-constant ripple of
    ``rad.Fld`` on the write-back elements (the O(h) bumps measured at standoff ~
    element size disappear identically; see the module docstring).

    ``res`` is the dict returned by ``vim.Solve`` / ``rad.Solve`` on a MeshSoftIron
    (it must carry the ``gfM`` GridFunction; order=1, straight TET meshes only --
    anything else raises).  Iron contribution only:

        B outside the iron = MU0 * (H_ext + H_demag)
        B inside  the iron = MU0 * (H_ext + H_demag + M)

    Exact for the discrete solution: boundary faces carry the linear sigma = M.n,
    every tet carries its constant volume charge -div M, and internal faces carry
    nothing (HDiv conformity) -- the three facts that make the order-1 evaluation
    closed-form.  Cost ~ O(N * n_faces) fully vectorized over the points."""
    gfM = res.get("gfM") if isinstance(res, dict) else None
    if gfM is None:
        raise ValueError(
            "vim.FieldFromSolution: res carries no 'gfM' GridFunction -- pass the dict "
            "returned by vim.Solve/rad.Solve (radia >= this version) unmodified.")
    if int(res.get("order", -1)) != 1:
        raise NotImplementedError(
            "vim.FieldFromSolution: wired for order=1 (RT1/BDM1) solutions only "
            "(got order=%r)." % (res.get("order"),))
    import ngsolve as ng
    mesh = gfM.space.mesh
    vcounts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if vcounts not in ({4}, {6}, {8}):
        raise NotImplementedError(
            "vim.FieldFromSolution: pure TET/HEX/WEDGE only (element vertex counts: %s)."
            % (sorted(vcounts),))
    pts = np.ascontiguousarray(np.asarray(points, float).reshape(-1, 3))
    if vcounts == {4} and res.get("curve_order") is None:
        a, G, c, _V = _linear_M_coefficients(gfM)
        _face_P, _face_w, bnd_P, bnd_sigma0, bnd_s = _tet_faces_and_charges(mesh, a, G, c)

        # C++ packing contract:
        #   volume [tet vertices 12, rho0, grad(rho) 3]
        #   surface [tri vertices 9, sigma0, grad(sigma) 3, Hessian(sigma) 9]
        # RT1 has constant rho=-div(M) and linear sigma=M.n.
        xyz = np.array([v.point for v in mesh.vertices], dtype=float)
        volume = []
        for e, el in enumerate(mesh.Elements(ng.VOL)):
            vertices = xyz[[v.nr for v in el.vertices]].reshape(-1)
            volume.extend(vertices.tolist())
            volume.extend([float(-np.trace(G[e])), 0.0, 0.0, 0.0])
        surface = []
        for vertices, sigma0, slope in zip(bnd_P, bnd_sigma0, bnd_s):
            surface.extend(np.asarray(vertices, float).reshape(-1).tolist())
            surface.append(float(sigma0))
            surface.extend(np.asarray(slope, float).tolist())
            surface.extend([0.0] * 9)

        def direct(query):
            values = _rp._hdiv_demag_field_batch(
                volume, surface, np.asarray(query, float).reshape(-1).tolist())
            return np.asarray(values, float).reshape(-1, 3) / (4.0 * np.pi)
    else:
        cloud = res.get("_rt1_field_cloud")
        if cloud is None:
            rho = -ng.div(gfM)
            normal = ng.specialcf.normal(mesh.dim)
            sigma = ng.InnerProduct(gfM.Trace(), normal)
            qpts, qweight = [], []
            # Degree 9 integrates the RT1 Q1/prism charge against the smooth
            # exterior kernel accurately while keeping the lazy field cloud
            # much smaller than the ChargeGram build.
            intorder = 9
            for i in range(mesh.GetNE(ng.VOL)):
                eid = ng.ElementId(ng.VOL, i)
                trafo = mesh.GetTrafo(eid)
                for ip in ng.IntegrationRule(mesh[eid].type, intorder):
                    mip = trafo(ip)
                    qpts.append([float(mip.point[k]) for k in range(3)])
                    qweight.append(float(ip.weight * mip.measure) * _scalar(rho(mip)))
            for i in range(mesh.GetNE(ng.BND)):
                eid = ng.ElementId(ng.BND, i)
                trafo = mesh.GetTrafo(eid)
                for ip in ng.IntegrationRule(mesh[eid].type, intorder):
                    mip = trafo(ip)
                    qpts.append([float(mip.point[k]) for k in range(3)])
                    qweight.append(float(ip.weight * mip.measure) * _scalar(sigma(mip)))
            cloud = (np.asarray(qpts, float), np.asarray(qweight, float))
            res["_rt1_field_cloud"] = cloud

        def direct(query):
            values = _rp._hdiv_charge_cloud_field(
                cloud[0].reshape(-1).tolist(), cloud[1].tolist(),
                np.asarray(query, float).reshape(-1).tolist())
            return np.asarray(values, float).reshape(-1, 3) / (4.0 * np.pi)

    field = direct(pts)
    image = res.get("image")
    if image is not None:
        from ._image import image_group, parse_image_string
        for axes, sign in image_group(parse_image_string(image)):
            reflected = pts.copy()
            reflected[:, list(axes)] *= -1.0
            contribution = direct(reflected)
            contribution[:, list(axes)] *= -1.0
            field += float(sign) * contribution
    return field


def magnetization_from_solution(res, points):
    """Evaluate the RT1 magnetization inside the solved mesh, zero outside."""
    gfM = res.get("gfM") if isinstance(res, dict) else None
    if gfM is None:
        raise ValueError("vim magnetization evaluation requires Solve's result dict")
    pts = np.asarray(points, float).reshape(-1, 3)
    values = np.zeros((len(pts), 3), dtype=float)
    mesh = gfM.space.mesh
    for index, point in enumerate(pts):
        try:
            mapped = mesh(float(point[0]), float(point[1]), float(point[2]))
            values[index] = [float(gfM[k](mapped)) for k in range(3)]
        except Exception:
            # NGSolve raises for a point outside the volume.  Its magnetization
            # contribution to B is exactly zero there.
            pass
    return values
