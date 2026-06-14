"""Correct field reconstruction for the HDiv-VIM.

Given a SOLVED magnetization gfM (an HDiv GridFunction, any order), reconstruct the magnetic field (B or H)
at arbitrary points via the EXACT analytic field of the per-element magnetization (radia.Fld -- the
documented "NGSolve magnetization -> Radia open-boundary field" pipeline).

DO NOT reconstruct the demag field as M_mass^{-1} N m (the energy operator applied + mass-inverted): the
demag operator N = B^T G B has a SOLENOIDAL NULLSPACE at order>=2 (divergence-free RT bubbles carry zero
charge, so B maps them to zero), so M_mass^{-1} N m does NOT equal the field -- it gives GARBAGE point
fields at order>=2 (it happens to work only at order 0/1).  The demag FACTOR (the energy Rayleigh quotient
m^T N m / m^T M_mass m) stays correct because it is the energy, not the field.  See
tests/feec/test_hdiv_vim_field_reconstruction.py (uniform-M sphere: B_inside = MU0*(2/3)*M at ALL orders).
"""
import numpy as np
import ngsolve as ng


def reconstruct_field(mesh, gfM, points, quantity="b", units="m"):
    """Magnetic field at `points` from a solved HDiv-VIM magnetization gfM, via the exact analytic field of
    the per-element M (radia.Fld).  CALLER wraps in TaskManager when combining with other NGSolve work.

    mesh     : the NGSolve tet mesh gfM lives on (the magnetized body, e.g. steel-only).
    gfM      : HDiv GridFunction (any order) holding the solved magnetization M(x) [A/m].
    points   : (N,3) array-like of query points (same length units as the mesh).
    quantity : 'b' (Tesla, = MU0*(H_demag + M) inside the body) or 'h' (A/m, the demag/stray H field).
    returns  : (N,3) ndarray of the field FROM THE MAGNETIZATION.  Add any SOURCE/coil field separately --
               radia.Fld here is the magnetization's contribution only.

    Per-element M is the centroid value (== the average for <= linear M; ~ for quadratic).  This is the
    field of the piecewise-constant-M body: EXACT for uniform M at ALL orders, and it carries the
    inter-element M variation the high-order solve resolved.  (Reconstructing a per-element POLYNOMIAL M
    field would need element subdivision or a polynomial-charge field kernel -- a future refinement.)
    """
    import radia as rad
    from radia.netgen_mesh_import import netgen_mesh_to_radia

    Mel = []
    for el in mesh.Elements(ng.VOL):
        c = np.array([mesh[v].point for v in el.vertices]).mean(0)   # element centroid
        mp = mesh(c[0], c[1], c[2])
        Mel.append([float(gfM[i](mp)) for i in range(3)])
    cont = netgen_mesh_to_radia(mesh, material=lambda i: {"magnetization": Mel[i]},
                                units=units, verbose=False)
    pts = np.asarray(points, float).reshape(-1, 3)
    return np.array([rad.Fld(cont, quantity, [float(p[0]), float(p[1]), float(p[2])]) for p in pts], float)


# ---------------------------------------------------------------------------------------------------
# Polynomial-charge field kernel (the order>=2 field: the centroid reconstruct_field above is only
# correct for piecewise-constant M; for a polynomial M the VOLUME charge rho = -div M is non-zero and
# its omission is a 90-230% error at div M != 0 -- measured, see tests/feec/test_hdiv_vim_poly_field.py).
# ---------------------------------------------------------------------------------------------------
def _scalar(val):
    return float(val[0] if isinstance(val, tuple) else val)


def reconstruct_field_polynomial(mesh, gfM, points, quad=4, quantity="h", include_volume=True):
    """EXTERNAL-point magnetic field of a polynomial magnetization gfM (HDiv order p), via its FULL
    polynomial charge -- volume rho = -div M (L2 order p-1) AND surface sigma = M.n (SurfaceL2 order p):

        H(r) = (1/4pi)[ INT_V (-div M)(r') (r-r')/|r-r'|^3 dV' + INT_S (M.n)(r') (r-r')/|r-r'|^3 dS' ]

    This is the order>=2 generalisation of `reconstruct_field` (which uses the per-element CENTROID M ==
    surface charge only; it DROPS the volume charge, a 90-230% error wherever div M != 0).  H = -grad phi_M
    with phi_M the magnetic scalar potential of the charges; verified on the uniform sphere (center -M/3,
    external dipole) and a linear-M body (volume-charge term essential; coarse->fine self-convergent).

    ELEMENT-TYPE AGNOSTIC: the quadrature points / weights / normals come from NGSolve's own
    ElementTransformation + IntegrationRule(el.type) + specialcf.normal, so the SAME code handles
    tetrahedral AND hexahedral (and prism) meshes, flat or curved (mip.weight carries the curved
    Jacobian, specialcf.normal the curved outward normal).

    mesh     : NGSolve mesh (tet or hex) gfM lives on.
    gfM      : HDiv GridFunction (order p) of the solved magnetization M(x) [A/m].
    points   : (N,3) query points -- EXTERNAL to the body (the integrands are 1/r^2-singular; an INTERNAL/
               near point needs the singular-aware kernel, a future step).
    quad     : controls the volume/surface integration order (intorder = 2*quad); raise for near-surface
               observation.
    quantity : 'h' (A/m, the demag/stray H) or 'b' (Tesla); 'b' = MU0*H at an EXTERNAL point (no M there).
    include_volume : drop the rho term when False (== the centroid/surface-only field) -- for diagnostics.
    returns  : (N,3) ndarray of the field FROM THE MAGNETIZATION (add any source/coil field separately).

    CALLER wraps in TaskManager when combining with other NGSolve work.  This is a reference (Python)
    implementation; a C++/H-matrix-accelerated version is the next productionisation step.
    """
    divM = ng.div(gfM)
    nsurf = ng.specialcf.normal(mesh.dim)
    sigma_cf = ng.InnerProduct(gfM.Trace(), nsurf)      # M.n on the boundary (correct outward normal)
    obs = np.asarray(points, float).reshape(-1, 3)
    H = np.zeros((len(obs), 3))
    inv4pi = 1.0 / (4.0 * np.pi)
    iord = 2 * quad

    if include_volume:                                  # volume charge rho = -div M
        for i in range(mesh.GetNE(ng.VOL)):
            ei = ng.ElementId(ng.VOL, i)
            trafo = mesh.GetTrafo(ei)
            for ip in ng.IntegrationRule(mesh[ei].type, iord):
                mip = trafo(ip)
                xq = np.array([mip.point[k] for k in range(3)])
                rho = -_scalar(divM(mip))
                d = obs - xq
                rn = np.linalg.norm(d, axis=1)
                H += inv4pi * (ip.weight * mip.measure) * rho * d / rn[:, None] ** 3

    for i in range(mesh.GetNE(ng.BND)):                 # surface charge sigma = M.n
        ei = ng.ElementId(ng.BND, i)
        trafo = mesh.GetTrafo(ei)
        for ip in ng.IntegrationRule(mesh[ei].type, iord):
            mip = trafo(ip)
            xq = np.array([mip.point[k] for k in range(3)])
            sig = _scalar(sigma_cf(mip))
            d = obs - xq
            rn = np.linalg.norm(d, axis=1)
            H += inv4pi * (ip.weight * mip.measure) * sig * d / rn[:, None] ** 3

    if quantity.lower() == "b":
        return (4e-7 * np.pi) * H                        # external point: B = MU0 * H (no M present)
    return H


# ---------------------------------------------------------------------------------------------------
# Step 2 building blocks -- the SINGULAR / near-singular field at INTERNAL points.  Each is an EXACT
# (or singularity-removed) kernel, validated standalone; assembling them into a near-surface-accurate
# internal-field call (and the curved-face case) is the next step.  See docs/hdiv_vim/POLYNOMIAL_CHARGE_FIELD.md.
# ---------------------------------------------------------------------------------------------------
def flat_triangle_charge_field(P, r):
    """EXACT field of a UNIFORM (sigma=1) charge on a flat triangle: INT_T (r-r')/|r-r'|^3 dS'
    (no 1/4pi), via the Wilton/Graglia closed form (solid-angle normal term + per-edge log tangential
    term).  Valid for ANY r (near, far, on-face principal value) -- this is the analytic surface
    near-field base.  P = (3,3) vertices; returns (3,).  Validated to machine precision vs fine Gauss
    (tests/feec/test_hdiv_vim_poly_field.py)."""
    P = np.asarray(P, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))                       # signed height of r above the plane
    r0 = r - h * n                                       # projection onto the plane
    F = np.zeros(3)
    omega = 0.0
    ah = abs(h)
    for i in range(3):
        a, b = P[i], P[(i + 1) % 3]
        u = (b - a) / np.linalg.norm(b - a)              # edge tangent
        m = np.cross(u, n)                               # in-plane edge normal
        t0 = float(np.dot(a - r0, m))
        sm = float(np.dot(a - r0, u))
        sp = float(np.dot(b - r0, u))
        R0sq = t0 * t0 + h * h
        Rm = np.sqrt(sm * sm + R0sq)
        Rp = np.sqrt(sp * sp + R0sq)
        denom_m = Rm + sm
        F += m * (-(np.log((Rp + sp) / denom_m) if denom_m > 1e-300 else 0.0))
        omega += np.arctan2(t0 * sp, R0sq + ah * Rp) - np.arctan2(t0 * sm, R0sq + ah * Rm)
    F += -n * np.sign(h) * omega
    return -F                                            # global sign fixed by the fine-Gauss validation


def reconstruct_field_internal(mesh, gfM, points, quad=4, nth=32, nph=64, quantity="h"):
    """INTERNAL/near-point magnetic field of a polynomial magnetization gfM (HDiv order p) on a TET mesh,
    via the singular-aware assembly:

        self-element VOLUME charge  -> tet_self_volume_field   (spherical ray-trace, non-singular)
        other-element VOLUME charge -> Gauss-Duffy             (non-singular: r is outside them)
        boundary SURFACE charge     -> flat_triangle_charge_field  (analytic, exact near AND far)

    This is the order>=2 field at the body's OWN quadrature points (the field Step-1
    reconstruct_field_polynomial cannot do -- its plain Gauss-Duffy is 1/r^2-singular there).  Validated
    on the uniform sphere: center -> -M/3, and near-surface the analytic surface beats Step-1 ~24x (the
    residual near-surface gap is FACETING -- the flat mesh's field genuinely differs from the smooth
    -M/3; removed by curving, a separate axis).

    TET meshes only (the self-element ray-trace uses tet face planes; boundary faces are triangles).
    Surface charge is taken CONSTANT per face (sigma = M.n at the face centroid) -- exact for uniform M,
    the constant-per-face part for polynomial M (the polynomial surface remainder is a future step).
    quantity 'h' = the demag/stray H; 'b' = MU0*(H + M(r)) (B INSIDE the body).  CALLER wraps TaskManager.
    """
    divM = ng.div(gfM)
    obs = np.asarray(points, float).reshape(-1, 3)
    inv4pi = 1.0 / (4.0 * np.pi)

    # precompute volume quadrature (point, rho, weight, element index) once -- independent of obs point
    vol_pts, vol_rho, vol_w, vol_eidx, elem_verts = [], [], [], [], []
    for i in range(mesh.GetNE(ng.VOL)):
        ei = ng.ElementId(ng.VOL, i)
        V = np.array([mesh[v].point for v in mesh[ei].vertices])
        elem_verts.append(V)
        trafo = mesh.GetTrafo(ei)
        for ip in ng.IntegrationRule(mesh[ei].type, 2 * quad):
            mip = trafo(ip)
            vol_pts.append([mip.point[k] for k in range(3)])
            vol_rho.append(-_scalar(divM(mip)))
            vol_w.append(ip.weight * mip.measure)
            vol_eidx.append(i)
    vol_pts = np.array(vol_pts); vol_rho = np.array(vol_rho); vol_w = np.array(vol_w)
    vol_eidx = np.array(vol_eidx)

    # precompute surface (triangle vertices, constant sigma = M.n at the centroid, geometric normal)
    surf = []
    for i in range(mesh.GetNE(ng.BND)):
        ei = ng.ElementId(ng.BND, i)
        P = np.array([mesh[v].point for v in mesh[ei].vertices])
        ngeom = np.cross(P[1] - P[0], P[2] - P[0]); ngeom = ngeom / np.linalg.norm(ngeom)
        c = P.mean(0)
        mc = mesh(c[0], c[1], c[2])
        sig = float(np.array([float(gfM[k](mc)) for k in range(3)]) @ ngeom)
        surf.append((P, sig))

    out = np.zeros((len(obs), 3))
    for q, r in enumerate(obs):
        # locate the self element (barycentric test)
        self_i = -1
        for i, V in enumerate(elem_verts):
            lam = np.linalg.solve(np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]).T, r - V[0])
            if lam.min() > -1e-9 and lam.sum() < 1 + 1e-9:
                self_i = i
                break
        H = np.zeros(3)
        # far volume (all quad points NOT in the self element) -- non-singular
        mask = vol_eidx != self_i
        d = r - vol_pts[mask]
        rn = np.linalg.norm(d, axis=1)
        H += inv4pi * np.sum(((vol_w[mask] * vol_rho[mask]) / rn ** 3)[:, None] * d, axis=0)
        # self volume (spherical ray-trace) -- non-singular
        if self_i >= 0:
            H += tet_self_volume_field(elem_verts[self_i], r,
                                       lambda p: -_scalar(divM(mesh(p[0], p[1], p[2]))), nth=nth, nph=nph)
        # surface (analytic triangle field, exact near and far)
        for P, sig in surf:
            H += inv4pi * sig * flat_triangle_charge_field(P, r)
        if quantity.lower() == "b":
            Mr = np.array([float(gfM[k](mesh(r[0], r[1], r[2]))) for k in range(3)])
            out[q] = (4e-7 * np.pi) * (H + Mr)
        else:
            out[q] = H
    return out


def _tet_face_planes(verts):
    """4 (n, d) half-space planes of a tet (n outward unit, interior satisfies n.x <= d)."""
    V = np.asarray(verts, float)
    c = V.mean(0)
    faces = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]   # vertex triples (opposite vertex 0,1,2,3)
    planes = []
    for (i, j, k) in faces:
        n = np.cross(V[j] - V[i], V[k] - V[i])
        n = n / np.linalg.norm(n)
        if np.dot(n, c - V[i]) > 0:                        # orient OUTWARD (interior on the n.x<=d side)
            n = -n
        planes.append((n, float(np.dot(n, V[i]))))
    return planes


def tet_self_volume_field(verts, r, rho_fn, nth=40, nph=80, ns=8):
    """SELF-element field of a tet's VOLUME charge rho at an INTERIOR point r, via the spherical
    ray-trace that removes the 1/r^2 singularity analytically:

        H = (1/4pi) INT_tet rho (r-r')/|r-r'|^3 dV' = -(1/4pi) INT_S2 s_hat [INT_0^smax rho(r+s s_hat) ds] dOmega

    (the s^2 volume Jacobian cancels 1/s^3).  smax(s_hat) = ray distance from r to the tet boundary;
    the inner s-integral is exact for polynomial rho via Gauss.  rho_fn(point)->charge density.
    Validated vs -(rho/4pi) grad(phi_tet) (the analytic Newtonian potential) to ~1e-3.

    NOTE: for a CONSTANT or LINEAR rho prefer tet_volume_field_linear -- it is a CLOSED FORM (exact to
    machine precision, no nth*nph*ns quadrature) and ~orders faster.  This spherical method is the
    fallback for genuinely higher-degree (quadratic+) rho where the closed form is not yet derived."""
    planes = _tet_face_planes(verts)
    xth, wth = np.polynomial.legendre.leggauss(nth)
    th = 0.5 * (xth + 1) * np.pi
    wth = wth * 0.5 * np.pi
    xph, wph = np.polynomial.legendre.leggauss(nph)
    ph = 0.5 * (xph + 1) * 2 * np.pi
    wph = wph * 0.5 * 2 * np.pi
    xs, ws = np.polynomial.legendre.leggauss(ns)
    r = np.asarray(r, float)
    H = np.zeros(3)
    for it in range(nth):
        sth = np.sin(th[it])
        for ip in range(nph):
            shat = np.array([sth * np.cos(ph[ip]), sth * np.sin(ph[ip]), np.cos(th[it])])
            cand = [(d - np.dot(nrm, r)) / np.dot(nrm, shat)
                    for (nrm, d) in planes if np.dot(nrm, shat) > 1e-12]
            if not cand:
                continue
            sm = min(cand)
            ss = 0.5 * (xs + 1) * sm
            wss = ws * 0.5 * sm
            integ_s = sum(wk * rho_fn(r + sk * shat) for sk, wk in zip(ss, wss))
            H += -(wth[it] * wph[ip] * sth) * shat * integ_s
    return H / (4.0 * np.pi)


# ---------------------------------------------------------------------------------------------------
# Analytic polynomial volume-charge field (EXACT, closed form) -- the order>=2 volume term.
#
# Key identity:  (r-r')/R^3 = grad'(1/R).  With the product rule + divergence theorem, the field of a
# polynomial volume charge rho reduces to LOWER-degree POTENTIAL integrals (one differential order down):
#
#     INT_V rho (r-r')/R^3 dV'  =  SUM_faces n_f INT_face rho/R dS'  -  INT_V (grad rho)/R dV'
#
# For LINEAR rho(r') = rho0 + g.r'  (grad rho = g const, == the order-2 volume charge -div M):
#
#     = SUM_faces n_f [rho0 I0_f + g . M1_f]  -  g * PhiTet        (PhiTet = INT_V 1/R dV')
#
# needing only the degree-1 TRIANGLE potential (const I0 + first moment M1) and the tet Newtonian
# potential PhiTet -- all closed form, EXACT everywhere (interior AND exterior), no quadrature.  This
# is the exact replacement for the ~1e-3 spherical tet_self_volume_field when rho is at most linear.
# See docs/hdiv_vim/POLYNOMIAL_CHARGE_FIELD.md.
# ---------------------------------------------------------------------------------------------------
def triangle_potential_const(P, r):
    """I0 = INT_T 1/R dS' on a flat triangle (Wilton 1984 closed form), pure-Python companion of
    flat_triangle_charge_field.  Bit-identical to the C++ _hdiv_tri_potential (validated)."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r0 = r - h * n
    ah = abs(h)
    I = 0.0
    omega = 0.0
    for i in range(3):
        a, b = P[i], P[(i + 1) % 3]
        u = (b - a) / np.linalg.norm(b - a)
        m = np.cross(u, n)
        t0 = float(np.dot(a - r0, m))
        sm = float(np.dot(a - r0, u))
        sp = float(np.dot(b - r0, u))
        R0sq = t0 * t0 + h * h
        Rm = np.sqrt(sm * sm + R0sq)
        Rp = np.sqrt(sp * sp + R0sq)
        dm = Rm + sm
        if dm > 1e-300:
            I += t0 * np.log((Rp + sp) / dm)
        omega += np.arctan2(t0 * sp, R0sq + ah * Rp) - np.arctan2(t0 * sm, R0sq + ah * Rm)
    return I - ah * omega


def _edge_R_integral(A, B, r):
    """INT_{edge A->B} |r - r'| dl, closed form  INT sqrt(u^2 + c^2) du."""
    t = B - A
    L = np.linalg.norm(t)
    if L < 1e-300:
        return 0.0
    that = t / L
    w = r - A
    l0 = float(np.dot(w, that))
    d2 = max(float(np.dot(w, w)) - l0 * l0, 0.0)
    u1, u2 = -l0, L - l0

    def F(u):
        if d2 < 1e-300:                                  # r on the edge line: INT |u| du
            return 0.5 * u * abs(u)
        s = np.sqrt(u * u + d2)
        return 0.5 * (u * s + d2 * np.arcsinh(u / np.sqrt(d2)))
    return F(u2) - F(u1)


def triangle_potential_moment(P, r):
    """M1 = INT_T r'/R dS' (3-vector) on a flat triangle, via the surface divergence theorem
    INT_T (r'-r_p)/R dS' = SUM_edges m_e INT_edge R dl  (m_e = in-plane outward edge normal,
    r_p = foot of the perpendicular from r).  Closed form, EXACT everywhere."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r_p = r - h * n
    cen = P.mean(axis=0)
    I0 = triangle_potential_const(P, r)
    acc = np.zeros(3)
    for i in range(3):
        A, B = P[i], P[(i + 1) % 3]
        u = (B - A) / np.linalg.norm(B - A)
        m = np.cross(u, n)
        if np.dot(m, 0.5 * (A + B) - cen) < 0:           # orient outward (away from centroid)
            m = -m
        acc += m * _edge_R_integral(A, B, r)
    return acc + r_p * I0


def tet_newtonian_potential(verts, r):
    """PhiTet = INT_V 1/R dV' over a tet, closed form  -1/2 SUM_faces h_f INT_face 1/R dS'
    (from 1/R = (1/2) lap'(R); h_f = signed height of r above outward face f).  Pure-Python
    companion of the C++ _hdiv_phi_tet (validated to machine precision); EXACT interior + exterior."""
    V = np.asarray(verts, float)
    r = np.asarray(r, float)
    c = V.mean(axis=0)
    I = 0.0
    for tri in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        P = V[list(tri)]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n = n / np.linalg.norm(n)
        if np.dot(n, c - P[0]) > 0:                      # orient outward
            n = -n
        I += -0.5 * float(np.dot(r - P[0], n)) * triangle_potential_const(P, r)
    return I


def tet_volume_field_linear(verts, r, rho0, grho):
    """EXACT field of a LINEAR volume charge rho(r') = rho0 + grho.r' on a tet:

        INT_V rho (r-r')/|r-r'|^3 dV'  (NO 1/4pi)  via the divergence-theorem recursion
        = SUM_faces n_f [rho0 I0_f + grho.M1_f]  -  grho * PhiTet

    Closed form (reuses triangle_potential_const / _moment + tet_newtonian_potential), EXACT to machine
    precision at ANY r (interior, surface, exterior) -- the exact, ~orders-faster replacement for the
    spherical tet_self_volume_field when rho is at most linear (the order-2 volume charge -div M).
    Multiply by (1/4pi) for the physical H contribution.  verts = (4,3); returns (3,)."""
    V = np.asarray(verts, float)
    r = np.asarray(r, float)
    grho = np.asarray(grho, float)
    c = V.mean(axis=0)
    acc = np.zeros(3)
    for tri in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        P = V[list(tri)]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n = n / np.linalg.norm(n)
        if np.dot(n, c - P[0]) > 0:                      # orient outward
            n = -n
        acc += n * (rho0 * triangle_potential_const(P, r) + np.dot(grho, triangle_potential_moment(P, r)))
    return acc - grho * tet_newtonian_potential(V, r)


def _edge_field_integral(A, B, r):
    """G_e = INT_{edge A->B} (r - r')/|r - r'| dl  (3-vector), closed form
    = (r-A) INT 1/R dl  -  t_hat INT l/R dl   (1D asinh / sqrt antiderivatives)."""
    t = B - A
    L = np.linalg.norm(t)
    if L < 1e-300:
        return np.zeros(3)
    that = t / L
    w = r - A
    l0 = float(np.dot(w, that))
    d2 = max(float(np.dot(w, w)) - l0 * l0, 0.0)
    d = np.sqrt(d2)
    u1, u2 = -l0, L - l0
    if d < 1e-300:                                       # r on the edge line
        asinh1 = np.sign(u2) * np.log(2 * abs(u2)) if abs(u2) > 0 else 0.0
        asinh0 = np.sign(u1) * np.log(2 * abs(u1)) if abs(u1) > 0 else 0.0
    else:
        asinh1 = np.arcsinh(u2 / d)
        asinh0 = np.arcsinh(u1 / d)
    int_1R = asinh1 - asinh0                              # INT 1/R dl
    int_lR = (np.sqrt(u2 * u2 + d2) - np.sqrt(u1 * u1 + d2)) + l0 * (asinh1 - asinh0)  # INT l/R dl
    return (r - A) * int_1R - that * int_lR


def linear_triangle_charge_field(P, r, sigma0, s):
    """EXACT field of a LINEAR surface charge sigma(r') = sigma0 + s.r' on a flat triangle:

        INT_T sigma (r-r')/|r-r'|^3 dS'  (NO 1/4pi)

    Derived as -grad_r phi_sigma with phi_sigma = INT_T sigma/R dS' = sigma0 I0 + s.M1 (the degree-1
    triangle potential), differentiated in closed form:

        = (sigma0 + s.r_p) F_const  -  SUM_edges (s.m_e) G_e  -  I0 s_par

    F_const = flat_triangle_charge_field (the constant-sigma field), G_e = _edge_field_integral,
    s_par = in-plane part of s, r_p = foot of the perpendicular.  EXACT to machine precision at ANY r
    (validated vs off-plane Gauss); s = 0 reproduces sigma0 * flat_triangle_charge_field exactly.
    This is the degree-1 SURFACE companion of tet_volume_field_linear (degree-1 VOLUME).  P = (3,3)
    vertices; returns (3,).  Multiply by 1/4pi for the physical H contribution."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    s = np.asarray(s, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r_p = r - h * n
    cen = P.mean(axis=0)
    Fc = flat_triangle_charge_field(P, r)
    I0 = triangle_potential_const(P, r)
    s_par = s - float(np.dot(s, n)) * n
    acc = (sigma0 + float(np.dot(s, r_p))) * Fc - I0 * s_par
    for i in range(3):
        A, B = P[i], P[(i + 1) % 3]
        u = (B - A) / np.linalg.norm(B - A)
        m = np.cross(u, n)
        if np.dot(m, 0.5 * (A + B) - cen) < 0:
            m = -m
        acc -= float(np.dot(s, m)) * _edge_field_integral(A, B, r)
    return acc


# ---------------------------------------------------------------------------------------------------
# Degree-2 (quadratic) charge moments + the quadratic VOLUME field.  Both new moments follow from the
# SAME identities one degree up:  V1 from 1/R = (1/2) lap'(R) weighted by r'_k;  M2 from the Hessian
# identity grad'_s grad'_s(R^3) = 3 (xi(x)xi/R + R P).  See docs/hdiv_vim/POLYNOMIAL_CHARGE_FIELD.md.
# ---------------------------------------------------------------------------------------------------
def _edge_R_xi_integral(A, B, r, r_p):
    """INT_{edge A->B} R (r' - r_p) dl  (3-vector) = (A - r_p) INT R dl + t_hat INT R l dl,
    with INT R dl and INT R l dl elementary (sqrt / asinh antiderivatives)."""
    t = B - A
    L = np.linalg.norm(t)
    if L < 1e-300:
        return np.zeros(3)
    that = t / L
    w = r - A
    l0 = float(np.dot(w, that))
    d2 = max(float(np.dot(w, w)) - l0 * l0, 0.0)
    d = np.sqrt(d2)
    u1, u2 = -l0, L - l0

    def Fsq(u):                                          # INT sqrt(u^2+d^2) du
        if d2 < 1e-300:
            return 0.5 * u * abs(u)
        return 0.5 * (u * np.sqrt(u * u + d2) + d2 * np.arcsinh(u / d))
    gR = Fsq(u2) - Fsq(u1)                               # INT R dl
    int_Rl = ((u2 * u2 + d2) ** 1.5 - (u1 * u1 + d2) ** 1.5) / 3.0 + l0 * (Fsq(u2) - Fsq(u1))  # INT R l dl
    return (A - r_p) * gR + that * int_Rl


def triangle_potential_moment2(P, r):
    """M2 = INT_T r' (x) r'/R dS' (symmetric 3x3) on a flat triangle.  Via xi (x) xi/R =
    (1/3) Hess_s(R^3) - R P  (xi = r'-r_p, P = in-plane projector), so
    INT_T xi(x)xi/R = SUM_edges (INT_edge R xi dl)(x)m_e - P INT_T R dS', with
    INT_T R dS' = (1/3)[SUM_e m_e.(INT_edge R xi dl) + h^2 I0]; then shift by r_p.  EXACT everywhere."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r_p = r - h * n
    cen = P.mean(axis=0)
    I0 = triangle_potential_const(P, r)
    Pproj = np.eye(3) - np.outer(n, n)
    Mxi1 = np.zeros(3)
    ohm_xi_m = np.zeros((3, 3))
    sum_m_dot = 0.0
    for i in range(3):
        A, B = P[i], P[(i + 1) % 3]
        u = (B - A) / np.linalg.norm(B - A)
        m = np.cross(u, n)
        if np.dot(m, 0.5 * (A + B) - cen) < 0:
            m = -m
        ARxi = _edge_R_xi_integral(A, B, r, r_p)
        Mxi1 += m * _edge_R_integral(A, B, r)
        ohm_xi_m += np.outer(ARxi, m)
        sum_m_dot += float(np.dot(m, ARxi))
    intR = (sum_m_dot + h * h * I0) / 3.0                # INT_T R dS'
    Mxi2 = ohm_xi_m - Pproj * intR
    return (np.outer(r_p, r_p) * I0 + np.outer(r_p, Mxi1) + np.outer(Mxi1, r_p) + Mxi2)


def tet_newtonian_moment(verts, r):
    """V1 = INT_V r'/R dV' over a tet (3-vector) = (1/3)[r PhiTet - SUM_faces h_f M1_f]
    (from 1/R = (1/2) lap'(R) weighted by r'_k; M1_f = surface first moment of face f, h_f = signed
    height of r above outward face f).  Closed form, EXACT interior + exterior."""
    V = np.asarray(verts, float)
    r = np.asarray(r, float)
    c = V.mean(axis=0)
    acc = r * tet_newtonian_potential(V, r)
    for tri in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        P = V[list(tri)]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n = n / np.linalg.norm(n)
        if np.dot(n, c - P[0]) > 0:
            n = -n
        acc -= float(np.dot(r - P[0], n)) * triangle_potential_moment(P, r)
    return acc / 3.0


def tet_volume_field_quadratic(verts, r, rho0, grho, Q):
    """EXACT field of a QUADRATIC volume charge rho(r') = rho0 + grho.r' + r'^T Q r' on a tet (Q
    symmetric 3x3): INT_V rho (r-r')/|r-r'|^3 dV' (NO 1/4pi), via the divergence-theorem recursion

        = SUM_faces n_f [rho0 I0_f + grho.M1_f + Q:M2_f]  -  (grho PhiTet + 2 Q.V1)

    Closed form (reuses triangle_potential_const / _moment / _moment2 + tet_newtonian_potential /
    _moment), EXACT to machine precision at ANY r.  Q = 0 reduces bit-identically to
    tet_volume_field_linear.  Multiply by 1/4pi for the physical H contribution."""
    V = np.asarray(verts, float)
    r = np.asarray(r, float)
    grho = np.asarray(grho, float)
    Q = np.asarray(Q, float)
    c = V.mean(axis=0)
    acc = np.zeros(3)
    for tri in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        P = V[list(tri)]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n = n / np.linalg.norm(n)
        if np.dot(n, c - P[0]) > 0:
            n = -n
        I0 = triangle_potential_const(P, r)
        M1 = triangle_potential_moment(P, r)
        M2 = triangle_potential_moment2(P, r)
        acc += n * (rho0 * I0 + float(np.dot(grho, M1)) + float(np.einsum("jk,jk", Q, M2)))
    return acc - (grho * tet_newtonian_potential(V, r) + 2.0 * (Q @ tet_newtonian_moment(V, r)))


def _edge_monomial_over_R(A, B, r):
    """(J0, J1, J2, t_hat, L) with Jk = INT_{edge A->B} l^k / R dl  (l = arc length from A),
    closed form (asinh / sqrt antiderivatives).  Building block for the quadratic surface field."""
    t = B - A
    L = np.linalg.norm(t)
    that = t / L
    w = r - A
    l0 = float(np.dot(w, that))
    d2 = max(float(np.dot(w, w)) - l0 * l0, 0.0)
    d = np.sqrt(d2)
    u1, u2 = -l0, L - l0
    if d < 1e-300:
        def asinh(u):
            return np.sign(u) * np.log(2 * abs(u)) if abs(u) > 0 else 0.0
    else:
        def asinh(u):
            return np.arcsinh(u / d)
    A_ = asinh(u2) - asinh(u1)                                       # INT 1/R du
    S_ = np.sqrt(u2 * u2 + d2) - np.sqrt(u1 * u1 + d2)               # INT u/R du
    U2 = (0.5 * (u2 * np.sqrt(u2 * u2 + d2) - d2 * asinh(u2))
          - 0.5 * (u1 * np.sqrt(u1 * u1 + d2) - d2 * asinh(u1)))     # INT u^2/R du
    J0 = A_
    J1 = S_ + l0 * A_                                                # l = u + l0
    J2 = U2 + 2 * l0 * S_ + l0 * l0 * A_
    return J0, J1, J2, that, L


def quadratic_triangle_charge_field(P, r, sigma0, s, S):
    """EXACT field of a QUADRATIC surface charge sigma(r') = sigma0 + s.r' + r'^T S r' on a flat
    triangle (S symmetric): INT_T sigma (r-r')/|r-r'|^3 dS' (NO 1/4pi), via the in-plane/normal split
    (r-r')/R^3 = grad'_s(1/R) + h n/R^3:

        in-plane  = SUM_e m_e INT_edge sigma/R dl  -  [P s I0 + 2 P S M1]
        normal    = h n [sigma0 J3_0 + s.J3_1 + S:J3_2]

    with the 1/R^3 moments J3_0 = INT_T 1/R^3 = (n.F_const)/h, J3_1 = INT_T r'/R^3, J3_2 = INT_T r'(x)r'/R^3
    (from xi(x)xi/R^3 = P/R - Hess_s R), all closed form.  EXACT to machine precision at ANY r (validated
    vs off-plane Gauss).  Subsumes the constant (s=S=0) and linear (S=0 -> matches
    linear_triangle_charge_field) cases.  The degree-2 SURFACE companion of tet_volume_field_quadratic.
    P = (3,3) vertices; returns (3,).  Multiply by 1/4pi for the physical H contribution."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    s = np.asarray(s, float)
    S = np.asarray(S, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r_p = r - h * n
    cen = P.mean(axis=0)
    Pproj = np.eye(3) - np.outer(n, n)
    Fc = flat_triangle_charge_field(P, r)
    I0 = triangle_potential_const(P, r)
    M1 = triangle_potential_moment(P, r)
    J3_0 = float(np.dot(n, Fc)) / h                                  # INT_T 1/R^3
    intxi1 = np.zeros(3)                                             # INT_T xi/R^3
    xixi_R3 = Pproj * I0                                             # INT_T xi(x)xi/R^3
    inplane = np.zeros(3)
    for i in range(3):
        A, B = P[i], P[(i + 1) % 3]
        J0, J1, J2, that, L = _edge_monomial_over_R(A, B, r)
        m = np.cross(that, n)
        if np.dot(m, 0.5 * (A + B) - cen) < 0:
            m = -m
        intxi1 -= m * J0
        Gxi = (A - r_p) * J0 + that * J1                            # INT_edge xi/R dl
        xixi_R3 -= np.outer(Gxi, m)
        c0 = sigma0 + float(np.dot(s, A)) + float(A @ S @ A)
        c1 = float(np.dot(s, that)) + 2.0 * float(A @ S @ that)
        c2 = float(that @ S @ that)
        inplane += m * (c0 * J0 + c1 * J1 + c2 * J2)                # m INT_edge sigma/R dl
    J3_1 = intxi1 + r_p * J3_0                                       # INT_T r'/R^3
    J3_2 = (xixi_R3 + np.outer(r_p, intxi1) + np.outer(intxi1, r_p) + np.outer(r_p, r_p) * J3_0)
    inplane -= (Pproj @ s) * I0 + 2.0 * (Pproj @ (S @ M1))          # - INT (grad_s sigma)/R
    normal = h * n * (sigma0 * J3_0 + float(np.dot(s, J3_1)) + float(np.einsum("jk,jk", S, J3_2)))
    return inplane + normal


# ---------------------------------------------------------------------------------------------------
# ARBITRARY-DEGREE polynomial charge field (the general assembler).  The degree-0/1/2 functions above
# are fast closed forms; these handle ANY polynomial degree via the general moment recursion:
#
#   master recursion  (2+p+k) INT_T xi^a R^p - p h^2 INT_T xi^a R^{p-2} = oint xi^a (xi.m) R^p dl
#   => surface 1/R moments  A_k  satisfy  A_k = (E^-1_k - h^2 B_k)/(k+1),  B_k = A_{k-2}(deriv) - edge
#      so A_k comes from A_{k-2} + edge integrals ALONE (h-safe: fold h^2 B_k = h^2[(a-1)A_{k-2} - edge]).
#   surface field  = in-plane/normal split (uses A and B moments + edge integrals).
#   volume potential moment  INT_V rp^a/R = 1/(d+2)[ -SUM_f h_f INT_face rp^a/R + SUM_i r_i a_i INT_V rp^{a-e_i}/R ]
#      (from 1/R = (1/2) lap'(R) + Euler), bottoming at PhiTet, reducing to surface potentials.
#   volume field  = SUM_f n_f INT_face rho/R - INT_V (grad rho)/R   (the divergence-theorem recursion).
#
# Validated vs Gauss to machine precision for cubic + quartic charge; reduces to the closed forms above.
# See docs/hdiv_vim/POLYNOMIAL_CHARGE_FIELD.md.  Pure-Python reference; cost ~ per (point, element).
# ---------------------------------------------------------------------------------------------------
from math import comb as _comb


def _inplane_triangle_setup(P, r):
    """In-plane orthonormal basis (e1,e2), normal n, foot r_p, signed height h, and per-edge data
    (in-plane coords of the edge endpoint and tangent/outward-normal, the u-range for l->u)."""
    P = np.asarray(P, float)
    r = np.asarray(r, float)
    n = np.cross(P[1] - P[0], P[2] - P[0])
    n = n / np.linalg.norm(n)
    h = float(np.dot(r - P[0], n))
    r_p = r - h * n
    cen = P.mean(axis=0)
    e1 = P[1] - P[0]
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    edges = []
    for i in range(3):
        A, B = P[i], P[(i + 1) % 3]
        L = np.linalg.norm(B - A)
        that = (B - A) / L
        m = np.cross(that, n)
        if np.dot(m, 0.5 * (A + B) - cen) < 0:
            m = -m
        w = r - A
        l0 = float(np.dot(w, that))
        d2 = max(float(np.dot(w, w)) - l0 * l0, 0.0)
        edges.append(dict(
            xiA=np.array([np.dot(A - r_p, e1), np.dot(A - r_p, e2)]),
            t2=np.array([np.dot(that, e1), np.dot(that, e2)]),
            m2=np.array([np.dot(m, e1), np.dot(m, e2)]),
            L=L, l0=l0, d2=d2))
    return dict(P=P, r=r, n=n, h=h, r_p=r_p, e1=e1, e2=e2, edges=edges)


def _edge_l_moments(edge, nmax):
    """Jl[n] = INT_edge l^n / R dl  (l = arc length from A), via the INT u^n/R reduction + l = u + l0."""
    d2 = edge["d2"]
    d = np.sqrt(d2)
    u1, u2 = -edge["l0"], edge["L"] - edge["l0"]

    def asinh(u):
        return np.arcsinh(u / d) if d > 1e-300 else (np.sign(u) * np.log(2 * abs(u)) if abs(u) > 0 else 0.0)
    W = np.zeros(nmax + 1)
    W[0] = asinh(u2) - asinh(u1)
    if nmax >= 1:
        W[1] = np.sqrt(u2 * u2 + d2) - np.sqrt(u1 * u1 + d2)
    for nn in range(2, nmax + 1):
        term = (u2 ** (nn - 1) * np.sqrt(u2 * u2 + d2) - u1 ** (nn - 1) * np.sqrt(u1 * u1 + d2)) / nn
        W[nn] = term - ((nn - 1) * d2 / nn) * W[nn - 2]
    l0 = edge["l0"]
    Jl = np.zeros(nmax + 1)
    for nn in range(nmax + 1):
        Jl[nn] = sum(_comb(nn, i) * l0 ** (nn - i) * W[i] for i in range(nn + 1))
    return Jl


def _edge_inplane_monomial(edge, a, b, Jl):
    """INT_edge xi1^a xi2^b / R dl  (xi1,xi2 = in-plane coords); expand (xiA + t l)^... as a poly in l."""
    p1 = np.array([edge["xiA"][0], edge["t2"][0]])
    p2 = np.array([edge["xiA"][1], edge["t2"][1]])
    poly = np.array([1.0])
    for _ in range(a):
        poly = np.convolve(poly, p1)
    for _ in range(b):
        poly = np.convolve(poly, p2)
    return float(sum(poly[nn] * Jl[nn] for nn in range(len(poly))))


def triangle_inplane_moments(P, r, degree, want_B=True):
    """The general surface moment dicts in the in-plane basis: A[(a,b)] = INT_T xi1^a xi2^b/R dS' (always,
    h-safe) and B[(a,b)] = INT_T xi1^a xi2^b/R^3 dS' (only if want_B and h != 0), for all a+b <= degree.
    A is self-contained (A_k from A_{k-2} + edge integrals).  Generalises triangle_potential_const /
    _moment / _moment2 (their r'-moments are r_p-shifts/contractions of these xi-moments)."""
    g = _inplane_triangle_setup(P, r)
    h = g["h"]
    edges = g["edges"]
    Jl = [_edge_l_moments(e, degree + 2) for e in edges]
    A = {(0, 0): triangle_potential_const(P, r)}
    B = {}

    def edge_R(e):                                       # INT_edge R dl
        u1, u2 = -e["l0"], e["L"] - e["l0"]
        d2 = e["d2"]

        def F(u):
            return (0.5 * (u * np.sqrt(u * u + d2) + d2 * np.arcsinh(u / np.sqrt(d2)))
                    if d2 > 1e-300 else 0.5 * u * abs(u))
        return F(u2) - F(u1)
    if degree >= 1:
        A1 = sum(e["m2"] * edge_R(e) for e in edges)     # INT_T xi/R = SUM m_e INT_edge R dl
        A[(1, 0)] = float(A1[0])
        A[(0, 1)] = float(A1[1])

    def Eneg1(a, b):                                     # oint xi1^a xi2^b (xi.m)/R dl
        return sum(e["m2"][0] * _edge_inplane_monomial(e, a + 1, b, Jl[i])
                   + e["m2"][1] * _edge_inplane_monomial(e, a, b + 1, Jl[i]) for i, e in enumerate(edges))

    def Eedge(j, p, q):                                  # oint xi1^p xi2^q m_j/R dl
        return sum(e["m2"][j] * _edge_inplane_monomial(e, p, q, Jl[i]) for i, e in enumerate(edges))
    for k in range(2, degree + 1):
        for a in range(k, -1, -1):
            b = k - a
            if a >= 1:                                   # h^2 B_k via the GRAD relation (finite at h=0)
                h2B = h * h * ((a - 1) * A.get((a - 2, b), 0.0) - Eedge(0, a - 1, b))
            else:
                h2B = h * h * ((b - 1) * A.get((a, b - 2), 0.0) - Eedge(1, a, b - 1))
            A[(a, b)] = (Eneg1(a, b) - h2B) / (k + 1)
    if want_B and abs(h) > 1e-12:
        Fc = flat_triangle_charge_field(P, r)
        B[(0, 0)] = float(np.dot(g["n"], Fc)) / h
        if degree >= 1:
            B[(1, 0)] = -Eedge(0, 0, 0)
            B[(0, 1)] = -Eedge(1, 0, 0)
        for k in range(2, degree + 1):
            for a in range(k, -1, -1):
                b = k - a
                if a >= 1:
                    B[(a, b)] = (a - 1) * A.get((a - 2, b), 0.0) - Eedge(0, a - 1, b)
                else:
                    B[(a, b)] = (b - 1) * A.get((a, b - 2), 0.0) - Eedge(1, a, b - 1)
    return A, B, g


def _fit_inplane_coeffs(charge_fn, g, degree):
    """Fit charge_fn(r') to in-plane monomial coeffs c[(a,b)] (xi1^a xi2^b) by sampling at barycentric
    nodes inside the triangle (exact for a degree-`degree` polynomial)."""
    P, e1, e2, r_p = g["P"], g["e1"], g["e2"], g["r_p"]
    N = max(degree, 1)
    nodes = [P[0] + (i / N) * (P[1] - P[0]) + (j / N) * (P[2] - P[0])
             for i in range(N + 1) for j in range(N + 1 - i)]
    monos = [(a, k - a) for k in range(degree + 1) for a in range(k + 1)]
    M = np.empty((len(nodes), len(monos)))
    rhs = np.empty(len(nodes))
    for ii, p in enumerate(nodes):
        x1 = float(np.dot(p - r_p, e1))
        x2 = float(np.dot(p - r_p, e2))
        for jj, (a, b) in enumerate(monos):
            M[ii, jj] = x1 ** a * x2 ** b
        rhs[ii] = charge_fn(p)
    c, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    return {monos[j]: float(c[j]) for j in range(len(monos))}


def polynomial_triangle_charge_field(P, r, sigma_fn, degree):
    """EXACT field of an ARBITRARY-degree polynomial surface charge sigma(r') on a flat triangle:
    INT_T sigma (r-r')/|r-r'|^3 dS' (NO 1/4pi), via the in-plane/normal split with the general moment
    recursion.  sigma_fn(point)->scalar (any polynomial of total degree <= `degree`).  EXACT to machine
    precision at ANY r (validated vs Gauss for cubic/quartic).  Subsumes flat_/linear_/quadratic_
    triangle_charge_field (use those closed forms when degree<=2 for speed).  Multiply by 1/4pi for H."""
    A, B, g = triangle_inplane_moments(P, r, degree + 1, want_B=True)
    e1, e2, n, h = g["e1"], g["e2"], g["n"], g["h"]
    edges = g["edges"]
    Jl = [_edge_l_moments(e, degree + 2) for e in edges]
    c = _fit_inplane_coeffs(sigma_fn, g, degree)
    normal = h * n * sum(cab * B[(a, b)] for (a, b), cab in c.items())
    inplane = np.zeros(3)
    for i, e in enumerate(edges):
        m3d = e["m2"][0] * e1 + e["m2"][1] * e2
        inplane += m3d * sum(cab * _edge_inplane_monomial(e, a, b, Jl[i]) for (a, b), cab in c.items())
    gx = sum(a * cab * A.get((a - 1, b), 0.0) for (a, b), cab in c.items() if a >= 1)
    gy = sum(b * cab * A.get((a, b - 1), 0.0) for (a, b), cab in c.items() if b >= 1)
    return inplane - (e1 * gx + e2 * gy) + normal


def _surface_potential_monomial(P, r, alpha, degree, cache):
    """INT_T r'^alpha/R dS' (global monomial alpha) -- the surface potential of the monomial on the
    triangle, via the in-plane A moments.  Cached per face."""
    key = P.tobytes()
    if key not in cache:
        A, _, g = triangle_inplane_moments(P, r, degree, want_B=False)
        cache[key] = (g, A)
    g, A = cache[key]
    c = _fit_inplane_coeffs(lambda p: p[0] ** alpha[0] * p[1] ** alpha[1] * p[2] ** alpha[2], g, sum(alpha))
    return sum(cab * A[(a, b)] for (a, b), cab in c.items())


def tet_volume_field_polynomial(verts, r, rho_fn, degree):
    """EXACT field of an ARBITRARY-degree polynomial volume charge rho(r') on a tet:
    INT_V rho (r-r')/|r-r'|^3 dV' (NO 1/4pi), via the divergence-theorem recursion

        = SUM_faces n_f INT_face rho/R dS'  -  INT_V (grad rho)/R dV'

    with the volume potential moments INT_V r'^a/R = 1/(|a|+2)[-SUM_f h_f INT_face r'^a/R
    + SUM_i r_i a_i INT_V r'^{a-e_i}/R] (bottoming at PhiTet, reducing to surface potentials).
    rho_fn(point)->scalar (any polynomial of total degree <= `degree`).  EXACT to machine precision at
    ANY r (validated vs Gauss for cubic; subsumes tet_volume_field_linear/quadratic).  Multiply by 1/4pi."""
    V = np.asarray(verts, float)
    r = np.asarray(r, float)
    c = V.mean(axis=0)
    faces = []
    for tri in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
        Pf = V[list(tri)]
        nf = np.cross(Pf[1] - Pf[0], Pf[2] - Pf[0])
        nf = nf / np.linalg.norm(nf)
        if np.dot(nf, c - Pf[0]) > 0:
            nf = -nf
        faces.append((Pf, nf))
    monos = [(ax, ay, k - ax - ay) for k in range(degree + 1)
             for ax in range(k + 1) for ay in range(k + 1 - ax)]
    N = degree + 1
    nodes = [V[0] + (i / N) * (V[1] - V[0]) + (j / N) * (V[2] - V[0]) + (k2 / N) * (V[3] - V[0])
             for i in range(N + 1) for j in range(N + 1 - i) for k2 in range(N + 1 - i - j)]
    M = np.empty((len(nodes), len(monos)))
    rhs = np.empty(len(nodes))
    for ii, p in enumerate(nodes):
        for jj, (ax, ay, az) in enumerate(monos):
            M[ii, jj] = p[0] ** ax * p[1] ** ay * p[2] ** az
        rhs[ii] = rho_fn(p)
    cg, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    coeffs = {monos[j]: float(cg[j]) for j in range(len(monos)) if abs(cg[j]) > 1e-13}
    cache = {}
    vmemo = {}

    def vol_pot_mono(alpha):
        d = sum(alpha)
        if d == 0:
            return tet_newtonian_potential(V, r)
        if alpha in vmemo:
            return vmemo[alpha]
        s = 0.0
        for (Pf, nf) in faces:
            s -= float(np.dot(r - Pf[0], nf)) * _surface_potential_monomial(Pf, r, alpha, degree, cache)
        for i in range(3):
            if alpha[i] >= 1:
                am = list(alpha)
                am[i] -= 1
                s += r[i] * alpha[i] * vol_pot_mono(tuple(am))
        vmemo[alpha] = s / (d + 2)
        return vmemo[alpha]

    F = np.zeros(3)
    for (Pf, nf) in faces:
        F += nf * sum(cc * _surface_potential_monomial(Pf, r, a, degree, cache) for a, cc in coeffs.items())
    for i in range(3):
        gi = 0.0
        for a, cc in coeffs.items():
            if a[i] >= 1:
                am = list(a)
                am[i] -= 1
                gi += cc * a[i] * vol_pot_mono(tuple(am))
        ei = np.zeros(3)
        ei[i] = 1.0
        F -= ei * gi
    return F
