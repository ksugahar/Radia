"""Loop-DOF extension of the scalar BIE + SIBC for genus-1 workpieces.

Physics
=======
On a genus-1 conductor (ring / tube) whose handle links the source flux,
the physical eddy current contains a NET circulating component -- the
shorted transformer turn.  The scalar BIE's surface current
``J_s = n x (-grad_s phi)`` with a single-valued phi carries ZERO net
current through any cut of the surface, so that component and its Lenz
screening are unrepresentable.  The plain scalar solver therefore needs
a cohomology extension whenever source flux links the surface handle.

Extension (one extra scalar DOF alpha = the net toroidal current):

    phi_multi = phi_u  (single-valued)  +  alpha * Theta

* ``Theta`` is the magnetic scalar potential of a UNIT current on a
  known ring inside the material wall (mid-wall, ray-cast from the cut).
  It is evaluated on the CUT-OPEN mesh by path integration of
  ``-H_ring . dl`` along a spanning tree; single-valuedness on the open
  mesh is guaranteed because the ring is in the same homology class as
  the cut (their surface loops have zero linking with it), and is
  VERIFIED by the uniform +-1 jump across the cut (fail-loud assert).
* The membrane (cut-disk double layer) term of the multivalued exterior
  identity cancels exactly against Theta's own identity; what survives is
  Theta's SIBC Neumann defect:

      alpha column  =  SL @ (gamma M^-1 K(Theta)  -  q_Theta),
      gamma = Z_s / (j omega mu_0),   q_Theta = -H_ring . n

  where ``K(Theta)`` is the element-local surface stiffness applied to
  the (jump-carrying) open representation of Theta -- no BEM operator is
  ever assembled on the open mesh (the duplicated cut vertices coincide
  geometrically and would poison the regular quadrature).
* Closure: Faraday's law on the cut loop (a surface loop),

      sum_edges Z_s (n x H_t) . dl  =  -j omega Phi_linked,
      Phi_linked = loop-integral of (A_inc + A_scat[J_s]) . dl,

  with the scattered vector potential from the panel currents
  (centroid-approximated single layer).

Validation
==========
* Analytic thin-wire shorted ring (torus R=30 mm / b=3 mm, copper,
  50 kHz, uniform axial field): net current alpha matches the classic
  ring circuit I = -j w Phi / (Z_s R/b + j w L_ring) to ~2 % in
  amplitude and ~0.2 deg in phase at two mesh resolutions; the frozen
  (alpha = 0) sub-system reproduces the production ScalarBIESIBCSolver
  solve to machine precision (same operators, same gauge).
* The frozen (alpha = 0) subsystem reproduces the production
  ``ScalarBIESIBCSolver`` solve to machine precision.

Scope / limitations (fail-loud, not silent)
===========================================
* genus-1 with ONE flux-linked handle (one cut, one alpha).  genus >= 2
  raises.
* The cut loop comes from ``radia.cohomology.surface_homology_loops``
  (the repo's single, gmsh-free cohomology engine: harmonic-1-cochain
  period matrix + QR-pivoted cotree selection); the Theta single-valued
  check (uniform jump) raises if the mid-wall ring construction fails
  (e.g. the ray cast finds no opposite wall, or the smoothed ring
  pierces the surface).
* P1 (order=1) intree-dense operators only -- the closed
  ``ScalarBIESIBCSolver`` must have been built with
  ``assemble_dense=True`` so ``M/K/SL/DL`` are available.

References
==========
* K. Sugahara, "Investigation of a Boundary Integral Equation
  n x H = J_s on Torus-Shaped Perfect Conductors," IEEE Trans.
  Antennas Propag., vol. 56, no. 3, pp. 722-725, Mar. 2008.  The dual
  defect of the same H^1(S) != 0 topology: on a torus the n x H = J_s
  BIE admits a spurious solution violating B . n = 0 (a harmonic
  null-space mode), closed there by ONE virtual-magnetic-current DOF +
  a one-point B . n = 0 constraint -- here the single-valued scalar
  potential LACKS the harmonic mode (representation deficit) and is
  closed by ONE loop DOF + the Faraday constraint.  Per handle, both
  formulations need exactly one extra DOF and one extra condition.
* M. Schoebinger and K. Hollaus, "An Effective Interface Approach for
  Multiply Connected Electromagnetic Shields," IEEE CEFC 2026
  conference digest, Thessaloniki, June 2026.  The FEM thin-shell
  sibling of the same topology treatment: the shielding sheet is
  reduced to a 2-D interface (1-D through-thickness analytic solution
  + nonlinear mu_eff lookup) and each hole gets ONE extra cohomology
  jump unknown (constant ``T = s_i e_z`` in the i-th hole), replacing
  the earlier non-physical auxiliary conductivity inside the holes.
* P. Dlotko, B. Kapidani, S. Pitassi, and R. Specogna, "Fake
  Conductivity or Cohomology: Which to Use When Solving Eddy Current
  Problems With h-Formulations?", IEEE Trans. Magn., vol. 55, no. 6,
  pp. 1-4, 2019.  The case AGAINST auxiliary-conductivity workarounds
  and FOR the cohomology treatment this module implements (via
  ``radia.cohomology``).

Part of the Radia project.
"""
from __future__ import annotations

import math
from collections import defaultdict, deque

import numpy as np

MU_0 = 4e-7 * math.pi


def A_from_filaments(points, filament_paths, currents):
    """Vector potential of straight filament segments (exact per segment).

    ``A_seg(x) = (mu0 I / 4 pi) t_hat * ln((L - xi + R2) / (-xi + R1))``
    with ``xi = (x - p1) . t_hat``, ``R1 = |x - p1|``, ``R2 = |x - p2|``.
    Supports complex per-filament currents.  Used as the ``A_inc_fn`` of
    ``solve_loop_extended`` for the PEEC (filament) coil source.

    Args:
        points: (n, 3) observation points [m].
        filament_paths: list of K filaments, each a list of (p1, p2)
            segment tuples (the ``coil_data["paths"]`` convention).
        currents: length-K complex per-filament currents [A].

    Returns:
        (n, 3) complex vector potential [T m].
    """
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    A = np.zeros((n, 3), dtype=complex)
    for fil, Ik in zip(filament_paths, currents):
        Ik = complex(Ik)
        if Ik == 0:
            continue
        for (p1, p2) in fil:
            p1 = np.asarray(p1, dtype=float)
            p2 = np.asarray(p2, dtype=float)
            dl = p2 - p1
            L = float(np.linalg.norm(dl))
            if L < 1e-15:
                continue
            t_hat = dl / L
            xi = (pts - p1[None, :]) @ t_hat
            R1 = np.linalg.norm(pts - p1[None, :], axis=1)
            R2 = np.linalg.norm(pts - p2[None, :], axis=1)
            num = np.maximum(L - xi + R2, 1e-30)
            den = np.maximum(-xi + R1, 1e-30)
            val = (MU_0 * Ik / (4.0 * math.pi)) * np.log(num / den)
            A += val[:, None] * t_hat[None, :]
    return A


# ----------------------------------------------------------------------
# topology
# ----------------------------------------------------------------------
def cut_open_surface(pts, tris, cut_loop):
    """Duplicate the cut-loop vertices and reattach the R-side vertex
    fans (per-vertex side resolution -- robust to zig-zag cuts).
    Requires a CONSISTENTLY oriented ``tris``.  Returns
    ``(pts_open, tris_open, dup)`` with ``dup[orig] = duplicate id``."""
    pts = np.asarray(pts, dtype=float)
    tris = np.asarray(tris, dtype=np.int64)
    nv = len(pts)
    n_c = len(cut_loop)
    ced = [(cut_loop[i], cut_loop[(i + 1) % n_c]) for i in range(n_c)]
    ces = {(min(a, b), max(a, b)) for a, b in ced}
    edge_tris = defaultdict(list)
    tov = defaultdict(list)
    for ti, t in enumerate(tris):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edge_tris[(min(a, b), max(a, b))].append(ti)
        for v in t:
            tov[v].append(ti)

    def cdir(t, a, b):
        return any(t[k] == a and t[(k + 1) % 3] == b for k in range(3))

    spv = {}
    for (a, b) in ced:
        for ti in edge_tris[(min(a, b), max(a, b))]:
            s = 'L' if cdir(tris[ti], a, b) else 'R'
            spv[(ti, a)] = s
            spv[(ti, b)] = s
    for v in cut_loop:
        ring = tov[v]
        changed = True
        while changed:
            changed = False
            for ti in ring:
                if (ti, v) in spv:
                    continue
                for tj in ring:
                    if (tj, v) not in spv:
                        continue
                    sh = set(map(int, tris[ti])) & set(map(int, tris[tj]))
                    if len(sh) == 2 and v in sh and \
                            tuple(sorted(sh)) not in ces:
                        spv[(ti, v)] = spv[(tj, v)]
                        changed = True
                        break
        missing = [ti for ti in ring if (ti, v) not in spv]
        if missing:
            raise ValueError(
                f"cut_open_surface: fan side propagation incomplete at "
                f"vertex {v} (tris {missing}) -- is the cut a simple "
                f"closed loop on a manifold surface?")

    dup = {v: nv + i for i, v in enumerate(cut_loop)}
    pts_o = np.vstack([pts, pts[cut_loop]])
    tris_o = tris.copy()
    for (ti, v), s in spv.items():
        if s == 'R':
            for k in range(3):
                if tris_o[ti, k] == v:
                    tris_o[ti, k] = dup[v]
    return pts_o, tris_o, dup


# ----------------------------------------------------------------------
# Theta carrier
# ----------------------------------------------------------------------
def _midwall_ring(pts, pts_o, tris_o, cut_loop, vnorm, n_smooth=5):
    """Ring source inside the wall: ray-cast from each cut vertex along
    the inward normal to the opposite wall, take the mid point, then
    lightly smooth.  Same homology class as the cut by construction."""
    v0 = pts_o[tris_o[:, 0]]
    e1 = pts_o[tris_o[:, 1]] - v0
    e2 = pts_o[tris_o[:, 2]] - v0

    def ray_second_hit(p, d, t_min=2e-4):
        h = np.cross(d[None, :], e2)
        a = np.einsum('ij,ij->i', e1, h)
        ok = np.abs(a) > 1e-14
        f = np.where(ok, 1.0 / np.where(ok, a, 1.0), 0.0)
        s = p[None, :] - v0
        u_b = f * np.einsum('ij,ij->i', s, h)
        q = np.cross(s, e1)
        v_b = f * np.einsum('j,ij->i', d, q)
        t = f * np.einsum('ij,ij->i', e2, q)
        hit = ok & (u_b >= -1e-9) & (v_b >= -1e-9) \
            & (u_b + v_b <= 1 + 1e-9) & (t > t_min)
        return float(np.min(t[hit])) if hit.any() else None

    ring = []
    for v in cut_loop:
        t_hit = ray_second_hit(pts[v], -vnorm[v])
        if t_hit is None:
            raise ValueError(
                f"loop extension: no opposite wall found from cut vertex "
                f"{v} along the inward normal -- cannot place the Theta "
                f"ring source inside the material.")
        ring.append(pts[v] - 0.5 * t_hit * vnorm[v])
    ring = np.array(ring)
    for _ in range(n_smooth):
        ring = 0.5 * ring + 0.25 * (np.roll(ring, 1, axis=0)
                                    + np.roll(ring, -1, axis=0))
    return ring


def _ring_H_factory(ring_pts):
    rp1 = np.asarray(ring_pts, dtype=float)
    rp2 = np.roll(rp1, -1, axis=0)

    def H_ring(robs):
        robs = np.asarray(robs, dtype=float)
        out = np.zeros((len(robs), 3))
        for s0 in range(0, len(robs), 400):
            s1 = min(len(robs), s0 + 400)
            a = rp1[None, :, :] - robs[s0:s1, None, :]
            b = rp2[None, :, :] - robs[s0:s1, None, :]
            na = np.linalg.norm(a, axis=2)
            nb = np.linalg.norm(b, axis=2)
            cr = np.cross(a, b)
            d2 = np.einsum('smd,smd->sm', cr, cr)
            ab = np.einsum('smd,smd->sm', a, b)
            coef = np.where(d2 > 1e-30,
                            (na + nb) * (1 - ab / (na * nb))
                            / np.maximum(d2, 1e-30), 0.0)
            out[s0:s1] = np.einsum('sm,smd->sd', coef, cr) / (4 * math.pi)
        return out

    return H_ring


def _theta_by_path_integration(pts_o, tris_o, dup, cut_loop, H_ring,
                               jump_tol=5e-3):
    """Theta on the open mesh: spanning-tree path integration of
    -H_ring . dl (4-pt Gauss per edge).  Verifies the uniform +-1 jump
    across the cut and raises otherwise."""
    nv_o = len(pts_o)
    adj = defaultdict(list)
    for t in tris_o:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            adj[a].append(b)
            adj[b].append(a)
    Theta = np.full(nv_o, np.nan)
    root = int(tris_o[0, 0])
    Theta[root] = 0.0
    par = {root: None}
    dq = deque([root])
    order = [root]
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in par:
                par[w] = u
                order.append(w)
                dq.append(w)
    gl_x, gl_w = np.polynomial.legendre.leggauss(4)
    gl_t = 0.5 * (gl_x + 1)
    gl_w = 0.5 * gl_w
    for w in order[1:]:
        u = par[w]
        p1, p2 = pts_o[u], pts_o[w]
        gp = p1[None, :] + gl_t[:, None] * (p2 - p1)[None, :]
        Hs = H_ring(gp)
        Theta[w] = Theta[u] - float(np.sum(gl_w * (Hs @ (p2 - p1))))
    if np.isnan(Theta).any():
        raise ValueError("loop extension: open mesh is disconnected -- "
                         "Theta path integration could not reach all "
                         "vertices.")
    jumps = np.array([Theta[dup[v]] - Theta[v] for v in cut_loop])
    if jumps.std() > jump_tol or abs(abs(jumps.mean()) - 1.0) > jump_tol:
        raise ValueError(
            f"loop extension: Theta jump across the cut is "
            f"{jumps.mean():+.4f} +- {jumps.std():.1e} (expected uniform "
            f"+-1).  The mid-wall ring source construction failed for "
            f"this geometry (it must stay inside the material wall).")
    return Theta, float(jumps.mean())


# ----------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------
def solve_loop_extended(bem_solver, phi_inc_nodal, Z_s, omega, A_inc_fn):
    """Solve the loop-extended scalar BIE + SIBC on a genus-1 workpiece.

    Args:
        bem_solver: a ``ScalarBIESIBCSolver`` built on the CLOSED
            workpiece surface mesh with ``assemble_dense=True`` and the
            intree P1 path (attributes ``M/K/SL/DL/M_inv/mesh`` used).
        phi_inc_nodal: (ndof,) complex incident scalar potential at the
            H1 nodes (e.g. the surface-Poisson reconstruction, or an
            exact expression for a uniform field).
        Z_s: complex Leontovich surface impedance (global scalar).
        omega: angular frequency [rad/s].
        A_inc_fn: callable ``A_inc_fn(points (n,3)) -> (n,3) complex`` --
            incident vector potential, used for the linked-flux term of
            the Faraday closure on the cut loop.

    Returns dict with ``alpha`` (net circulating current, complex [A]),
    ``P_total`` [W], ``H_t_rms`` [A/m], ``phi_u``, ``phi_open``,
    ``theta_jump``, ``cut_n_vertices``, plus ``P_frozen`` / ``Ht_frozen``
    (the alpha=0 sub-solve, == the plain solver, for diagnostics).
    """
    from ngsolve import BND

    mesh = bem_solver.mesh
    M, K = bem_solver.M, bem_solver.K
    SL, DL, M_inv = bem_solver.SL, bem_solver.DL, bem_solver.M_inv
    if SL is None or DL is None:
        raise ValueError("loop extension needs assemble_dense=True "
                         "(dense SL/DL) on the ScalarBIESIBCSolver.")
    gamma = Z_s / (1j * omega * MU_0) if omega > 0 else 0.0

    pts = np.array([[mesh.vertices[i].point[j] for j in range(3)]
                    for i in range(mesh.nv)])
    tris = np.array([[v.nr for v in el.vertices]
                     for el in mesh.Elements(BND)], dtype=np.int64)
    nv, nt = len(pts), len(tris)
    if bem_solver.ndof != nv:
        raise ValueError(
            f"loop extension supports the P1 nodal path only "
            f"(ndof={bem_solver.ndof} != nv={nv}).")

    # verify consistent orientation (the extractors now guarantee it; a
    # mesh from another route must be oriented BEFORE the solver is
    # built, because the BEM operators bake the winding in).
    dcount = defaultdict(int)
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            dcount[(int(a), int(b))] += 1
    conflicts = sum(1 for c in dcount.values() if c > 1)
    if conflicts:
        raise ValueError(
            f"loop extension: surface winding is inconsistent "
            f"({conflicts} directed-edge conflicts).  Orient the mesh "
            f"(surface_mesh_extract.orient_surface_triangles) BEFORE "
            f"building the BEM solver.")

    # genus / cut
    E = {(min(a, b), max(a, b)) for t in tris for a, b in
         ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
    chi = nv - len(E) + nt
    genus = (2 - chi) // 2
    if genus != 1:
        raise ValueError(
            f"loop extension supports genus-1 surfaces (got genus="
            f"{genus}, chi={chi}).  genus 0 needs no extension; "
            f"genus >= 2 needs one DOF per flux-linked handle (not "
            f"implemented).")
    # Cut selection through the repo's single cohomology engine
    # (radia.cohomology, the gmsh-free port).  The flux-linked cut must
    # be the PURE toroidal class (w_tor, w_pol) = (+-1, 0): a mixed
    # representative like (1, -1) is a valid basis element but forces
    # alpha to carry an equal poloidal component (wrong physics) and its
    # zig-zag geometry breaks the mid-wall Theta ray-cast.  A fixed
    # 2-loop generator basis is NOT guaranteed class-pure (measured on
    # Takahashi: {(1,-1), (0,1)}), so classify ALL cotree fundamental
    # cycles instead: the class map is a linear image of the harmonic
    # period vector, calibrated on one independent (QR-pivoted) pair, so
    # every cycle costs O(1) -- then take the geometrically shortest
    # pure-toroidal simple loop.  Winding numbers of a CLOSED loop are
    # exact integers: toroidal about the z axis, poloidal about the
    # wall-section centroid (rho0, z0) -- valid for any z-axis genus-1
    # workpiece (torus, tube, ring), whose section centroid lies inside
    # the material wall.
    import scipy.linalg as sla
    from radia.cohomology import surface_fundamental_cycles

    _b1, Pi, expand, _cotree = surface_fundamental_cycles(tris, nv=nv)
    rho_all = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    rho0, z0 = float(rho_all.mean()), float(pts[:, 2].mean())

    def _windings(p):
        q = pts[list(p) + [p[0]]]
        th = np.unwrap(np.arctan2(q[:, 1], q[:, 0]))
        rho = np.sqrt(q[:, 0] ** 2 + q[:, 1] ** 2)
        ph = np.unwrap(np.arctan2(q[:, 2] - z0, rho - rho0))
        return np.array([(th[-1] - th[0]) / (2 * math.pi),
                         (ph[-1] - ph[0]) / (2 * math.pi)])

    _, _, piv = sla.qr(Pi.T, pivoting=True, mode="economic")
    sel = piv[:2]
    W_sel = np.array([_windings(expand(k)) for k in sel])      # (2, 2)
    if abs(np.linalg.det(W_sel)) < 0.5:
        raise ValueError(
            f"loop extension: winding map is singular on the generator "
            f"pair (windings {np.rint(W_sel).astype(int).tolist()}) -- "
            f"the workpiece is not a z-axis genus-1 body in the assumed "
            f"sense (toroidal about z, poloidal about the wall section).")
    # windings(cycle e) = Wmap @ Pi[e]  (linear class map, calibrated once)
    Wmap = W_sel.T @ np.linalg.inv(Pi[sel].T)
    W_all = Pi @ Wmap.T                                        # (nc, 2)
    W_int = np.rint(W_all)
    cands = np.where((np.abs(W_all - W_int) < 0.2).all(axis=1)
                     & (np.abs(W_int[:, 0]) == 1)
                     & (W_int[:, 1] == 0))[0]
    if len(cands) == 0:
        raise ValueError(
            f"loop extension: no PURE toroidal fundamental cycle found "
            f"among {Pi.shape[0]} cotree candidates (generator windings "
            f"{np.rint(W_sel).astype(int).tolist()}).  Re-meshing "
            f"usually resolves this.")

    def _geo_len(p):
        q = pts[list(p) + [p[0]]]
        return float(np.sum(np.linalg.norm(np.diff(q, axis=0), axis=1)))

    cut = min((expand(int(k)) for k in cands), key=_geo_len)
    n_cut = len(cut)

    pts_o, tris_o, dup = cut_open_surface(pts, tris, cut)
    nv_o = len(pts_o)
    Tmap = np.arange(nv_o)
    for v, vd in dup.items():
        Tmap[vd] = v

    def reduce_vec(x):
        out = np.zeros(nv, dtype=complex)
        np.add.at(out, Tmap, x.astype(complex))
        return out

    # geometry tables
    areas = np.zeros(nt)
    normals = np.zeros((nt, 3))
    gvecs = np.zeros((nt, 3, 3))
    for ti, t in enumerate(tris_o):
        P = pts_o[t]
        nvec = np.cross(P[1] - P[0], P[2] - P[0])
        A2 = np.linalg.norm(nvec)
        n = nvec / A2
        areas[ti] = 0.5 * A2
        normals[ti] = n
        for k in range(3):
            opp = P[(k + 2) % 3] - P[(k + 1) % 3]
            h = np.cross(n, opp) / (2 * areas[ti])
            if h @ (P[k] - P[(k + 1) % 3]) < 0:
                h = -h
            gvecs[ti, k] = h
    cents = pts_o[tris_o].mean(axis=1)
    vnorm = np.zeros((nv, 3))
    for ti, t in enumerate(tris_o):
        for k in range(3):
            vnorm[Tmap[t[k]]] += normals[ti] * areas[ti]
    vnorm /= np.maximum(np.linalg.norm(vnorm, axis=1), 1e-30)[:, None]

    # Theta carrier
    ring = _midwall_ring(pts, pts_o, tris_o, cut, vnorm)
    H_ring = _ring_H_factory(ring)
    Theta, theta_jump = _theta_by_path_integration(
        pts_o, tris_o, dup, cut, H_ring)
    Theta = Theta.astype(complex)

    qT = -np.einsum('ij,ij->i', H_ring(pts), vnorm)
    rK_T = np.zeros(nv, dtype=complex)
    for ti, t in enumerate(tris_o):
        gTh = (gvecs[ti, 0] * Theta[t[0]] + gvecs[ti, 1] * Theta[t[1]]
               + gvecs[ti, 2] * Theta[t[2]])
        for k in range(3):
            rK_T[Tmap[t[k]]] += areas[ti] * (gvecs[ti, k] @ gTh)
    a_col = SL @ (gamma * (M_inv @ rK_T) - qT.astype(complex))

    A_sys = (0.5 * M - DL + gamma * (SL @ M_inv @ K)).astype(complex)
    RHS = (M @ np.asarray(phi_inc_nodal, dtype=complex))

    # Faraday closure on the cut loop
    edge_tris = defaultdict(list)
    for ti, t in enumerate(tris_o):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edge_tris[(min(int(a), int(b)), max(int(a), int(b)))].append(ti)
    ced = [(cut[i], cut[(i + 1) % n_cut]) for i in range(n_cut)]
    row_E = np.zeros(nv_o, dtype=complex)
    for (a, b) in ced:
        # open triangles adjacent to this cut edge (L uses originals, R
        # uses duplicates)
        tis = []
        for key in ((min(a, b), max(a, b)),
                    (min(dup[a], dup[b]), max(dup[a], dup[b]))):
            tis.extend(edge_tris.get(key, []))
        if not tis:
            raise ValueError(f"loop extension: no open triangle adjacent "
                             f"to cut edge ({a},{b}).")
        dl = pts[b] - pts[a]
        wgt = 1.0 / len(tis)
        for ti in tis:
            n = normals[ti]
            t = tris_o[ti]
            for k in range(3):
                row_E[t[k]] += wgt * Z_s * np.dot(
                    np.cross(n, -gvecs[ti, k]), dl)
    mids = np.array([0.5 * (pts[a] + pts[b]) for (a, b) in ced])
    dls = np.array([pts[b] - pts[a] for (a, b) in ced])
    row_Phi = np.zeros(nv_o, dtype=complex)
    for ti in range(nt):
        d = np.linalg.norm(mids - cents[ti][None, :], axis=1)
        w = MU_0 / (4 * math.pi) * areas[ti] / np.maximum(d, 1e-9)
        n = normals[ti]
        t = tris_o[ti]
        for k in range(3):
            Js_k = np.cross(n, -gvecs[ti, k])
            row_Phi[t[k]] += np.sum(w * (dls @ Js_k))
    A_mid = np.asarray(A_inc_fn(mids), dtype=complex)
    Phi_inc_loop = complex(np.sum(np.einsum('ij,ij->i', A_mid, dls)))

    fE = reduce_vec(row_E)
    fP = reduce_vec(row_Phi)
    f_alpha = complex((row_E + 1j * omega * row_Phi) @ Theta)

    # assemble + solve (Lagrange mean-zero gauge, production style)
    Mrow = M.sum(axis=1).astype(complex)
    N = nv + 2
    A2 = np.zeros((N, N), dtype=complex)
    b2 = np.zeros(N, dtype=complex)
    A2[:nv, :nv] = A_sys
    A2[:nv, nv] = a_col
    A2[:nv, nv + 1] = Mrow
    A2[nv, :nv] = fE + 1j * omega * fP
    A2[nv, nv] = f_alpha
    A2[nv + 1, :nv] = Mrow
    b2[:nv] = RHS
    b2[nv] = -1j * omega * Phi_inc_loop
    u = np.linalg.solve(A2, b2)
    phi_u, alpha = u[:nv], complex(u[nv])

    def _P_Ht(phi_u_v, alpha_v):
        phi_o = phi_u_v[Tmap] + alpha_v * Theta
        s2 = 0.0
        for ti, t in enumerate(tris_o):
            g = -(gvecs[ti, 0] * phi_o[t[0]] + gvecs[ti, 1] * phi_o[t[1]]
                  + gvecs[ti, 2] * phi_o[t[2]])
            s2 += areas[ti] * float(np.sum(np.abs(g) ** 2))
        return (0.5 * Z_s.real * s2, math.sqrt(s2 / float(areas.sum())))

    P_total, H_t_rms = _P_Ht(phi_u, alpha)

    # frozen sub-solve (diagnostic; == the plain production solve)
    Ng = nv + 1
    A0 = np.zeros((Ng, Ng), dtype=complex)
    b0 = np.zeros(Ng, dtype=complex)
    A0[:nv, :nv] = A_sys
    A0[:nv, nv] = Mrow
    A0[nv, :nv] = Mrow
    b0[:nv] = RHS
    phi_f = np.linalg.solve(A0, b0)[:nv]
    P_frozen, Ht_frozen = _P_Ht(phi_f, 0.0)

    return {
        "alpha": alpha,
        "P_total": float(P_total),
        "H_t_rms": float(H_t_rms),
        "P_frozen": float(P_frozen),
        "Ht_frozen": float(Ht_frozen),
        "phi_u": phi_u,
        "phi_open": phi_u[Tmap] + alpha * Theta,
        "theta_jump": float(theta_jump),
        "cut_n_vertices": int(n_cut),
        "genus": 1,
    }
