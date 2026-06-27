"""radia.moment_galerkin._assemble -- moment-basis assembly for the SYMMETRIC moment-Galerkin MMMM demag.

The moment-Galerkin demag operator is N = B^T G B, the SAME symmetric structure HDiv-VIM uses, but on the
MMMM moment basis (constant magnetization M per hexahedron) instead of HDiv's RT0 flux.  G is the EXACT
analytic charge-Gram H-matrix (`radia._radia_pybind._ChargeGramHMatrix`, the same C++ kernel HDiv-VIM ships);
B maps the per-element dipole moment to the face surface charges sigma = M . n_face; M_mass = diag(V_e) is the
magnetization mass.  The physical SPD demag system is

    ( (1/chi) M_mass + B^T G B ) m = M_mass H_ext ,      m = per-hex dipole moments (3/hex).

This file builds the SPARSE pieces (B, M_mass) + the C++ charge-Gram geometry (the hex faces as flat
triangles).  The heavy compute (G build + the CG/MINRES solve) is the existing C++ kernel; this is the thin
assembly layer (the same Python-orchestration pattern as radia.vim._solve).

Symmetry is by construction: G is symmetric (the Coulomb mutual energy is reciprocal), so N = B^T G B is
symmetric to machine precision (validated 1e-16 on the loop-heavy C-yoke vs collocation MMMM's 1.56) -> the
loop modes are field-null by construction -> mu_r-independent / loop-free convergence with no loop-star.

ACA defaults (leaf=40, eta=0.5): tuned so the H-matvec reproduces the exact analytic Gram to ~1e-4 on
near-heavy geometry.  The MMMM 6-charge/hex (= 12 face-triangles) come as a SET, so leaf must not split a
hex's charge group (leaf~40 keeps ~3 hexes intact); eta=0.5 keeps near-touching blocks dense.  See the
charge-Gram ACA tuning note (eta is the dominant knob, eps barely matters; leaf co-dominant).
"""
import numpy as np
import scipy.sparse as sp
import radia._radia_pybind as _rp

# Hexahedron face -> vertex indices (outward winding handled by per-triangle orientation below).
HEX_FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [3, 2, 6, 7], [0, 3, 7, 4], [1, 2, 6, 5]]


def _hex_volume(V):
    """Signed-robust hexahedron volume from the 8 vertices (sum of the 6 face-pyramid volumes to the
    centroid; works for affine and mildly-distorted hexes)."""
    V = np.asarray(V, float)
    c = V.mean(0)
    vol = 0.0
    for f in HEX_FACES:
        a, b, cc, d = (V[i] for i in f)
        for (p, q, r) in ((a, b, cc), (a, cc, d)):
            # tetra (centroid, p, q, r) volume, signed; abs-sum gives the hex volume
            vol += abs(np.dot(np.cross(q - p, r - p), c - p)) / 6.0
    return vol


def _tri_outward_normal(p0, p1, p2, hex_centroid):
    """Unit normal of triangle (p0,p1,p2), oriented to point AWAY from the hex centroid (outward)."""
    n = np.cross(p1 - p0, p2 - p0)
    nn = np.linalg.norm(n)
    if nn == 0.0:
        raise ValueError("moment_galerkin: degenerate (zero-area) hex face triangle")
    n = n / nn
    if np.dot(n, (0.5 * (p0 + p1 + p2)) - hex_centroid) < 0.0:
        n = -n
    return n


def assemble_moment_system(hexes, *, eps=1e-9, leaf=40, eta=0.5, near_factor=1e30, build=True):
    """Build the moment-Galerkin demag pieces for a list of hexahedra.

    Parameters
    ----------
    hexes : sequence of (8,3) arrays
        Hexahedron vertices (meters), each in the standard Radia/HEX_FACES vertex order.
    eps, leaf, eta : ACA H-matrix parameters for the C++ charge-Gram (validated defaults; eta is the
        dominant accuracy knob, leaf must not split a hex's 12-triangle charge group -> ~40).
    near_factor : analytic NEAR/FAR entry split (1e30 = all-analytic, exact entries).
    build : build the H-matrix now (True; matvec available) or defer (False; .entry() oracle only).

    Returns
    -------
    dict with
      G        : the C++ `_ChargeGramHMatrix` (the demag operator core; N v = B^T G B v)
      B        : scipy CSR (n_charge x 3*n_hex)   sigma = M . n_tri per face-triangle
      M_mass   : scipy CSR (3*n_hex)              diag(V_e) magnetization mass
      vols     : (n_hex,) hexahedron volumes
      n_hex    : number of hexes
      n_charge : number of face-triangle charges (12 per hex)
    """
    hexes = [np.asarray(V, float) for V in hexes]
    n = len(hexes)
    if n == 0:
        raise ValueError("moment_galerkin: empty hex list")
    face_tris = []
    rows, cols, data = [], [], []
    tri = 0
    for e, V in enumerate(hexes):
        c = V.mean(0)
        for f in range(6):
            a, b, cc, d = (V[i] for i in HEX_FACES[f])
            for T in ((a, b, cc), (a, cc, d)):
                p0, p1, p2 = T
                nrm = _tri_outward_normal(p0, p1, p2, c)
                face_tris += [p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]]
                for k in range(3):
                    if nrm[k] != 0.0:
                        rows.append(tri)
                        cols.append(3 * e + k)
                        data.append(float(nrm[k]))   # sigma_tri = n_tri . M  -> B[tri, 3e+k] = n_tri[k]
                tri += 1
    n_charge = tri
    B = sp.csr_matrix((data, (rows, cols)), shape=(n_charge, 3 * n))
    vols = np.array([_hex_volume(V) for V in hexes], float)
    M_mass = sp.diags(np.repeat(vols, 3)).tocsr()
    G = _rp._ChargeGramHMatrix(cell_verts=[], face_verts=face_tris, n_el=0,
                               eps=eps, leaf=int(leaf), eta=float(eta), near_factor=near_factor,
                               build=bool(build))
    return {"G": G, "B": B, "M_mass": M_mass, "vols": vols, "n_hex": n, "n_charge": n_charge}
