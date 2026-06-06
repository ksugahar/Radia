"""
demag_jacobi_stress.py -- where does Jacobi stop working?  Stress the star-block
conditioning with ELONGATED bodies (high aspect ratio -> small longitudinal demag
factor -> small mu_min -> the Jacobi gap closes).  Radia EXACT ObjRecMag field.

The breakdown onset is mu_min < 1/chi.  For mu_r<=1e4, 1/chi=1e-4, so Jacobi holds
until the smallest demag eigenvalue drops below ~1e-4 (extreme elongation). This maps
the Jacobi-OK regime for the charge-carrying (star) demag operator.
"""
import json, os
import numpy as np
import radia as rad
from scipy.sparse.linalg import gmres, LinearOperator

HERE = os.path.dirname(os.path.abspath(__file__))

def demag_matrix(shape, a=0.01, distort=0.0, seed=1):
    """blocks at integer-grid positions in `shape` (list of (i,j,k)); exact demag N."""
    rng = np.random.default_rng(seed)
    centers = np.array([[i*a, j*a, k*a] for (i, j, k) in shape], float)
    if distort > 0:
        centers = centers + distort * a * (rng.random(centers.shape) - 0.5)
    nb = len(centers); N = np.zeros((3*nb, 3*nb)); dims = [a, a, a]
    for jb in range(nb):
        for b in range(3):
            rad.UtiDelAll()
            M = [0.0, 0.0, 0.0]; M[b] = 1.0
            src = rad.ObjRecMag(centers[jb].tolist(), dims, M)
            H = np.array(rad.Fld(src, 'h', centers.tolist())).reshape(nb, 3)
            N[:, 3*jb + b] = H.reshape(-1)
    rad.UtiDelAll()
    return N

def jac_iters(A, tol=1e-8, maxit=4000):
    d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
    Minv = LinearOperator(A.shape, matvec=lambda x: x / d)
    it = {"n": 0}
    gmres(A, np.ones(A.shape[0]), M=Minv, rtol=tol, maxiter=maxit,
          callback=lambda xk: it.__setitem__("n", it["n"]+1), callback_type="pr_norm")
    return it["n"]

cases = {
    "cube 3x3x3":   [(i, j, k) for i in range(3) for j in range(3) for k in range(3)],
    "rod 1x1x4":    [(0, 0, k) for k in range(4)],
    "rod 1x1x8":    [(0, 0, k) for k in range(8)],
    "rod 1x1x16":   [(0, 0, k) for k in range(16)],
    "slab 4x4x1":   [(i, j, 0) for i in range(4) for j in range(4)],
    "C-yoke":       ([(i, 0, 0) for i in range(5)] + [(0, j, 0) for j in range(1, 5)]
                     + [(i, 4, 0) for i in range(1, 5)]),   # a C-shaped ring (open)
}
res = {}
for name, shape in cases.items():
    for distort, dtag in [(0.0, ""), (0.3, "+dist")]:
        N = demag_matrix(shape, distort=distort)
        nd = N.shape[0]
        mu = np.sort(np.abs(np.linalg.eigvals(N).real))
        mu_min = mu[0]
        row = {"ndof": nd, "mu_min": float(mu_min)}
        line = f"{name+dtag:16s} ndof={nd:3d} |mu|_min={mu_min:.2e}"
        for mur in (1e2, 1e4):
            chi = mur - 1.0
            A = (1.0/chi)*np.eye(nd) - N
            sv = np.linalg.svd(A, compute_uv=False)
            it = jac_iters(A)
            row[f"mur_{mur:.0e}"] = {"cond": float(sv[0]/sv[-1]), "jac_iters": int(it),
                                     "gap": float(mu_min*chi)}
            line += f" | mu_r={mur:.0e}: cond={sv[0]/sv[-1]:.1e} jac={it:4d} gap={mu_min*chi:6.1f}"
        print(line)
        res[name+dtag] = row

print("\ngap = mu_min/(1/chi); gap>>1 => Jacobi bounded. Breakdown onset when gap ~ 1 (mu_min<1/chi).")
with open(os.path.join(HERE, "demag_jacobi_stress.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "demag_jacobi_stress.json"))
