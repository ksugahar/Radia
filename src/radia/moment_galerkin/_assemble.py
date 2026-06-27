"""radia.moment_galerkin._assemble -- moment-basis assembly for the SYMMETRIC moment-Galerkin MMMM demag.

The moment-Galerkin demag operator is N = B^T G B, the SAME symmetric structure HDiv-VIM uses, but on the
MMMM moment basis instead of HDiv's RT0 flux.  G is the EXACT analytic charge-Gram H-matrix
(`radia._radia_pybind._ChargeGramHMatrix`, the same C++ kernel HDiv-VIM ships); B maps the per-element moment
amplitudes to the face surface charges sigma; M_mass is the magnetization mass.  The physical SPD demag system

    ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext

is solved by the existing C++ Krylov kernel (this file is the thin assembly layer, the radia.vim pattern).

Two moment orders (`quad=`):
  quad=False (default) -- 3 DOF/hex: a CONSTANT magnetization M_0 per hex (the dipole / standard demag).
      B col k = sigma_f = M_0 . n_f = n_f[k];  M_mass = diag(V_e).
  quad=True            -- 5 DOF/hex: 3 dipole + 2 QUAD residual-eigenmode amplitudes (the geometry-adaptive
      quadrupole patterns phi_q = the zero-monopole, zero-dipole-moment face-charge modes, the same
      Gram-Schmidt basis the C++ collocation MMMM `momentResidualEigenmodes` uses).  The quad mode's
      magnetization is the LINEAR field M_q(x) = Dm_q.(x-c), Dm_q = sum_f Ae_f n_f (x) d_f * (phi_q/Ae_f),
      so the quad MASS is the magnetization Gram W_quad[q,r] = Dm_q : M2 : Dm_r, M2 = INT (x-c)(x-c) dV.
      The 2 quad DOF buy per-element accuracy under SKEW / gradient loads (the iron-yoke case); they are
      ~null for axial / symmetric loads (correct physics, never worse).  Validated: a 3-hex bar 5-DOF field
      is strictly closer to a fine-mesh rad.Solve(yano) truth than dipole-only at every near probe under an
      oblique field, ties for axial; the cube quad amplitudes are ~0 by symmetry (the dipole 1/3 demag is
      preserved).  NOTE: on a cube / near-cube the quad mass is near-singular (the residual modes carry ~zero
      second moment) -- benign because the cube also nulls the quad amplitudes, and the diagonal-Jacobi C++
      solve (auto_prec) stays well-defined (its diagonal includes the demag self-energy N_diag).

Symmetry is by construction: G is symmetric (Coulomb reciprocity) -> N = B^T G B symmetric to machine
precision (validated 1e-16 on the loop-heavy C-yoke vs collocation MMMM's 1.56) -> loop modes field-null ->
mu_r-independent / loop-free convergence with no loop-star.

ACA defaults (leaf=40, eta=0.5): tuned so the H-matvec reproduces the exact analytic Gram to ~1e-4 on
near-heavy geometry (eta is the dominant knob; the MMMM 12-face-triangle charge group must not be split, so
leaf~40 keeps ~3 hexes intact).
"""
import numpy as np
import scipy.sparse as sp
import radia._radia_pybind as _rp

# Hexahedron face -> vertex indices.
HEX_FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [3, 2, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5]]


def _face_geom(V):
    """Per-face outward unit normal nf, area Ae, centroid fc, d = fc - element-centroid; element centroid c
    and divergence-theorem volume Ve.  Planar-faced hexes have per-face == per-triangle normals."""
    c = V.mean(0)
    fc = np.zeros((6, 3)); nf = np.zeros((6, 3)); Ae = np.zeros(6); d = np.zeros((6, 3))
    for f in range(6):
        P = V[HEX_FACES[f]]
        fcen = P.mean(0)
        n = np.cross(P[1] - P[0], P[2] - P[0]) + np.cross(P[2] - P[0], P[3] - P[0])
        nn = np.linalg.norm(n)
        if nn == 0.0:
            raise ValueError("moment_galerkin: degenerate (zero-area) hex face")
        n = n / nn
        if np.dot(n, fcen - c) < 0.0:
            n = -n
        area = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0])) + \
               0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[0]))
        fc[f] = fcen; nf[f] = n; Ae[f] = area; d[f] = fcen - c
    Ve = sum(Ae[f] * np.dot(fc[f], nf[f]) for f in range(6)) / 3.0
    return fc, nf, Ae, d, c, Ve


def _residual_eigenmodes(Ae, d, nF=6):
    """Orthonormal basis of {Ae, Ae*d_x, Ae*d_y, Ae*d_z}^perp in R^nF (the nF-4 = 2 quad eigenmodes for a hex);
    a Python mirror of the C++ momentResidualEigenmodes deterministic Gram-Schmidt."""
    basis = []

    def ortho_add(v):
        v = v.copy()
        for b in basis:
            v -= (b @ v) * b
        nrm = np.linalg.norm(v)
        if nrm > 1e-9:
            basis.append(v / nrm)
            return True
        return False

    for i in range(4):
        ortho_add(Ae.copy() if i == 0 else Ae * d[:, i - 1])
    phi = []
    for e in range(nF):
        v = np.zeros(nF); v[e] = 1.0
        if ortho_add(v):
            phi.append(basis[-1].copy())
    return np.array(phi)


def _second_moment(V, c):
    """M2 = INT (x-c)(x-c) dV over the hex by 2x2x2 Gauss on the trilinear map (exact for affine hexes)."""
    g = 1.0 / np.sqrt(3.0)
    signs = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
             (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
    nodes = [(sx * g, sy * g, sz * g) for (sx, sy, sz) in signs]
    M2 = np.zeros((3, 3))
    for (xi, et, ze) in nodes:
        N = np.array([(1 + sx * xi) * (1 + sy * et) * (1 + sz * ze) for (sx, sy, sz) in signs]) / 8.0
        dN = np.zeros((8, 3))
        for i, (sx, sy, sz) in enumerate(signs):
            dN[i, 0] = sx * (1 + sy * et) * (1 + sz * ze) / 8.0
            dN[i, 1] = sy * (1 + sx * xi) * (1 + sz * ze) / 8.0
            dN[i, 2] = sz * (1 + sx * xi) * (1 + sy * et) / 8.0
        x = N @ V
        detJ = abs(np.linalg.det(dN.T @ V))
        r = x - c
        M2 += np.outer(r, r) * detJ
    return M2


def assemble_moment_system(hexes, *, quad=False, eps=1e-9, leaf=40, eta=0.5, near_factor=1e30, build=True):
    """Build the moment-Galerkin demag pieces for a list of hexahedra.

    quad : False -> 3 DOF/hex (dipole); True -> 5 DOF/hex (3 dipole + 2 quad residual eigenmodes).
    eps, leaf, eta : ACA H-matrix parameters for the C++ charge-Gram (validated defaults).
    near_factor : analytic NEAR/FAR entry split (1e30 = all-analytic).   build : build the H-matrix now.

    Returns dict(G, B (csr n_charge x ndof_per*n_hex), M_mass (csr), vols, n_hex, n_charge, ndof_per,
                 all_tris (list of (3,3) triangle vertex arrays in B-row order, for field reconstruction)).
    """
    hexes = [np.asarray(V, float) for V in hexes]
    n = len(hexes)
    if n == 0:
        raise ValueError("moment_galerkin: empty hex list")
    ndof_per = 5 if quad else 3
    face_tris = []
    all_tris = []
    rows, cols, data = [], [], []
    tri = 0
    vols = np.zeros(n)
    Wblocks = []
    for e, V in enumerate(hexes):
        fc, nf, Ae, d, c, Ve = _face_geom(V)
        vols[e] = Ve
        col_face = np.zeros((6, ndof_per))
        for k in range(3):
            col_face[:, k] = nf[:, k]                 # dipole: sigma_f = M_0 . n_f
        if quad:
            phi = _residual_eigenmodes(Ae, d)
            col_face[:, 3] = phi[0]
            col_face[:, 4] = phi[1]
        for f in range(6):
            P = V[HEX_FACES[f]]
            for T in ((P[0], P[1], P[2]), (P[0], P[2], P[3])):
                face_tris += [T[0][0], T[0][1], T[0][2], T[1][0], T[1][1], T[1][2], T[2][0], T[2][1], T[2][2]]
                all_tris.append(np.array(T, float))
                for a in range(ndof_per):
                    if col_face[f, a] != 0.0:
                        rows.append(tri); cols.append(ndof_per * e + a); data.append(float(col_face[f, a]))
                tri += 1
        We = np.zeros((ndof_per, ndof_per))
        for k in range(3):
            We[k, k] = Ve                             # dipole mass = V_e (magnetization energy of unit M_0)
        if quad:
            M2 = _second_moment(V, c)
            Dm = []
            for q in range(2):
                Bm = phi[q] / Ae
                D = sum(Ae[f] * np.outer(nf[f], d[f]) * Bm[f] for f in range(6))
                Dm.append(D)
            for q in range(2):
                for r in range(2):
                    We[3 + q, 3 + r] = float(np.einsum("ij,ij->", Dm[q] @ M2, Dm[r]))  # Dm_q : M2 : Dm_r
        Wblocks.append(We)
    n_charge = tri
    B = sp.csr_matrix((data, (rows, cols)), shape=(n_charge, ndof_per * n))
    M_mass = sp.block_diag(Wblocks).tocsr() if quad else sp.diags(np.repeat(vols, 3)).tocsr()
    G = _rp._ChargeGramHMatrix(cell_verts=[], face_verts=face_tris, n_el=0,
                               eps=eps, leaf=int(leaf), eta=float(eta), near_factor=near_factor,
                               build=bool(build))
    return {"G": G, "B": B, "M_mass": M_mass, "vols": vols, "n_hex": n, "n_charge": n_charge,
            "ndof_per": ndof_per, "all_tris": all_tris}
