"""
hdiv_vim_accurate_g.py -- Phase 4 prototype: ACCURATE Coulomb Gram G (vs the centroid-monopole
placeholder), validated to give the true cube demag factor (1/3) + a positive-definite G.

The centroid-monopole G (1-point quadrature) gives uniform demag ~0.31 (7% off 1/3) and is slightly
NON-PD (spurious negative demag factors ~-0.49, shared with the NGSolve prototype).  Here every
charge-charge interaction G_ab = IntInt_{a x b} 1/(4pi r) is computed by SUB-POINT QUADRATURE:
each volume cell -> nsub^3 sub-points, each boundary face -> nsub^2 sub-points (uniform), weight =
sub-measure.  Off-diagonal: double sum over sub-points (convergent, captures the near field of
adjacent cells).  Coincident sub-points (same sub-cell) use the analytic cube/square self-energy.
This converges to the true 6D/4D Coulomb integrals as nsub grows.

Validates: uniform M_z demag factor -> 1/3, and N = B^T G B becomes POSITIVE SEMI-DEFINITE
(no spurious negative demag factors).  Compare to centroid-G.  (C++ production will use the Wilton
analytic face integral + H-matrix; this prototype confirms the accuracy target is reachable.)
"""
import json, os
import numpy as np
from math import pi
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from hdiv_vim_structured import build_structured_rt0, charge_map, assemble_mass

C_CUBE = 1.88231
C_SQ = 2.97321

def subpoints(faces, n_cell, cell_c, cell_V, nsub):
    """sub-point cloud for every charge cell (volume cells, then boundary faces).
    returns: pts (Np,3), w (Np,), charge_idx (Np,) mapping each sub-point to its charge cell."""
    # charge cell order: volume cells [0..n_cell), then boundary faces
    bnd = [f for f in range(len(faces)) if faces[f]["bnd"]]
    pts, w, cid = [], [], []
    g = (np.arange(nsub) + 0.5) / nsub
    # volume cells: nsub^3 grid in the cube (regular cell of size h around centroid)
    h = cell_V[0] ** (1/3.)
    gx = (g - 0.5) * h  # offsets from centroid in [-h/2, h/2]
    GX, GY, GZ = np.meshgrid(gx, gx, gx, indexing="ij")
    off = np.stack([GX.ravel(), GY.ravel(), GZ.ravel()], axis=1)
    wv = cell_V[0] / nsub**3
    for c in range(n_cell):
        for o in off:
            pts.append(cell_c[c] + o); w.append(wv); cid.append(c)
    # boundary faces: nsub^2 grid on the face (normal axis fixed)
    g2 = (g - 0.5)
    for bi, f in enumerate(bnd):
        fc = faces[f]; ax = fc["ax"]; A = fc["area"]; hs = A ** 0.5
        u = [0, 1, 2]; u.remove(ax)  # in-plane axes
        wa = A / nsub**2
        for a1 in g2*hs:
            for a2 in g2*hs:
                p = np.array(fc["c"], float)
                p[u[0]] += a1; p[u[1]] += a2
                pts.append(p); w.append(wa); cid.append(n_cell + bi)
    return np.array(pts), np.array(w), np.array(cid), len(bnd)

def accurate_gram(faces, n_cell, cell_c, cell_V, nsub):
    pts, w, cid, n_bnd = subpoints(faces, n_cell, cell_c, cell_V, nsub)
    n_charge = n_cell + n_bnd
    Np = len(pts)
    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    np.fill_diagonal(D, np.inf)
    K = 1.0 / (4*pi*D)                       # off (sub-point) Coulomb; diagonal handled below
    # coincident sub-point self: cube/square self-energy of the sub-cell
    self_sub = np.empty(Np)
    for p in range(Np):
        if cid[p] < n_cell: self_sub[p] = C_CUBE * w[p]**(5/3.) / (4*pi)
        else:               self_sub[p] = C_SQ   * w[p]**(1.5)  / (4*pi)
    np.fill_diagonal(K, 0.0)                  # coincident sub-points handled via self_sub below
    # G_ab = sum_{p in a, q in b}  w_p w_q K[p,q]  + coincident sub-cell self-energy
    WK = (w[:, None] * w[None, :]) * K
    np.fill_diagonal(WK, self_sub)            # coincident sub-point p==q: analytic sub-cell self
    # aggregate sub-points -> charge cells
    G = np.zeros((n_charge, n_charge))
    np.add.at(G, (cid[:, None], cid[None, :]), WK)
    return G, n_charge

def demag_check(nx, ny, nz, nsub):
    faces, nc, cc, cV = build_structured_rt0(nx, ny, nz, h=1.0)
    B, _, _ = charge_map(faces, nc, cV)
    M = assemble_mass(faces, nc, h=1.0)
    G, _ = accurate_gram(faces, nc, cc, cV, nsub)
    N = B.T @ G @ B
    df = np.sort(sla.eigh(N, M, eigvals_only=True))
    v = np.zeros(len(faces))
    for fi, f in enumerate(faces):
        if f["ax"] == 2: v[fi] = f["area"]
    dz = (v @ N @ v) / (v @ M @ v)
    return dz, df

print("=== Phase 4: accurate G via sub-point quadrature -- uniform demag -> 1/3 + PD? ===")
print(f"  true cube uniform demag factor = 1/3 = {1/3:.4f}")
res = {}
for nsub in (2, 3, 4, 6):
    dz, df = demag_check(3, 3, 3, nsub)
    nneg = int(np.sum(df < -1e-3))
    print(f"  nsub={nsub}: uniform M_z demag = {dz:.4f}, demag range [{df[0]:.3e}, {df[-1]:.3e}], "
          f"#(demag<0)={nneg}")
    res[f"nsub{nsub}"] = {"uniform_demag": float(dz), "min": float(df[0]), "max": float(df[-1]), "n_neg": nneg}
print("  (centroid-G baseline: uniform 0.310, range [-0.49, 0.89], #neg ~17)")
with open(os.path.join(HERE, "hdiv_vim_accurate_g.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "hdiv_vim_accurate_g.json"))
