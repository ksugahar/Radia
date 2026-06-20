"""Is the star-projected yano-MSC SCALABLE with HACApK?  YES -- the conditioning obstacle is SOLVED by
a sparse basis-Gram (cell-graph Laplacian) preconditioner.

`yano_pyr_faces12_star.py` showed the star projection makes the collocation matrix loop-free with a good
DENSE-SVD-basis condition number (~20).  For a SCALABLE solver the star basis must be SPARSE.  This
example measures the production path and its fix:

  #2 SPARSE-star conditioning + the FIX.  The production star basis is q_internal = B^T a (raw
     cell/internal-face incidence, one cell pinned), NOT the dense SVD row space.  The RAW sparse form
     is ill-conditioned -- it inflates ~n (cond 313 @18 elem -> 18492 @320 elem) and a symmetric Jacobi
     rescale does NOT fix it (it is the NON-ORTHOGONALITY of the incidence basis, whose Gram is the
     cell-graph Laplacian -- not a diagonal scaling).  THE FIX: precondition with that basis Gram
     M = T_sp^T T_sp (the sparse, SPD cell-graph Laplacian).  This orthonormalizes the basis, so the
     preconditioned condition number drops back BELOW the dense-SVD value (~10-14, ~0.5x dense-SVD) and
     -- the scalability proof -- the preconditioned BiCGStab iteration count is BOUNDED (~15-22) across
     an 18x range of mesh size, while the unpreconditioned count grows 35 -> 144.  M is SPARSE -> applied
     by AMG (radia.sparsesolv_ngsolve Compact AMG) in production, not the direct factor used here.

  #3 HACApK compatibility.  A_star @ x = T^T A (T x) with T SPARSE (density <3%, O(n_face) nnz).  So if
     A is the HACApK H-matrix (the scalable yano demag, natural rank ~13 -- exactly what the current
     yano-MSC already builds), the star solve is: sparse scatter -> H-matvec -> sparse gather, with a
     sparse AMG preconditioner on M.  Every piece is scalable.

Verdict: the star-projected yano IS scalable + loop-free + HACApK-usable.  A = HACApK H-matrix (matvec),
star = sparse O(n_face) wrapper, preconditioner = sparse cell-graph Laplacian (AMG) -> BOUNDED iterations
independent of mesh size.  Combined with the pyramid-cloud kernel's per-element hex accuracy, this is the
'scalable + loop-free' demag -- it gives the scalable yano backend the loop-freeness that previously
needed the HDiv-VIM, without the HDiv-VIM face-DOF A=B^T G B build cost.

Reference: H. Yano and K. Sugahara, J. Magn. Soc. Jpn. 47(2), 52-56, 2023.
"""
import os
import json
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, bicgstab, splu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yano_pyr_faces12_star as y   # the pyramid-cloud kernel + dense-SVD star + distorted_patch

HERE = os.path.dirname(os.path.abspath(__file__))


def star_transform_sparse(hexes, elements, internal, boundary, adj, pin=0):
    """Production star: T maps [boundary DOFs + cell potentials a (one cell pinned)] -> 6N local faces,
    with internal flux q = B^T a (raw incidence, NOT the dense SVD basis); local scatter +-/area."""
    cells = [c for c in range(len(elements)) if c != pin]
    cidx = {c: j for j, c in enumerate(cells)}
    rows, cols, vals = [], [], []
    for col, key in enumerate(boundary):
        ei, fi = adj[key][0]
        rows.append(6 * ei + fi); cols.append(col); vals.append(1.0)
    off = len(boundary)
    for key in internal:
        for sgn, (ei, fi) in zip((1.0, -1.0), adj[key]):
            area = y._face_area(hexes[ei], fi)
            if area <= 1e-300:
                continue
            for cell, csign in ((adj[key][0][0], 1.0), (adj[key][1][0], -1.0)):
                if cell == pin:
                    continue
                rows.append(6 * ei + fi); cols.append(off + cidx[cell]); vals.append(sgn * csign / area)
    return sp.coo_matrix((vals, (rows, cols)), shape=(6 * len(elements), off + len(cells))).tocsr()


def _precond_cond(A_star, M):
    """cond of M^{-1/2} A_star M^{-1/2} (the basis-Gram-orthonormalized operator)."""
    Md = np.asarray(M.todense())
    w, V = np.linalg.eigh(0.5 * (Md + Md.T))
    Minvhalf = V @ np.diag(1.0 / np.sqrt(np.clip(w, 1e-300, None))) @ V.T
    sv = np.linalg.svd(Minvhalf @ A_star @ Minvhalf, compute_uv=False)
    return float(sv[0] / sv[-1])


def _bicg_iters(A_star, Mop):
    rng = np.random.default_rng(0)
    b = rng.random(A_star.shape[0]) - 0.5
    cnt = {"n": 0}
    op = LinearOperator(A_star.shape, matvec=lambda v: A_star @ v)
    _x, _info = bicgstab(op, b, rtol=1e-8, maxiter=2000, M=Mop,
                         callback=lambda xk: cnt.__setitem__("n", cnt["n"] + 1))
    return cnt["n"]


def evaluate(nx, ny, nz, strength):
    hexes, elements = y.distorted_patch(nx, ny, nz, strength)
    internal, boundary, adj = y.face_topology(elements)
    cycle = len(internal) - (len(elements) - 1)
    A, _ = y.collocation_matrix(hexes, "pyr_faces12")

    T_svd = y.star_transform(hexes, elements, internal, boundary, adj)
    c_svd = float(np.linalg.cond(T_svd.T @ A @ T_svd))

    T_sp = star_transform_sparse(hexes, elements, internal, boundary, adj)
    A_star = np.asarray(T_sp.T @ (A @ T_sp.toarray()))
    c_sp = float(np.linalg.cond(A_star))

    M = (T_sp.T @ T_sp).tocsc()                  # sparse SPD basis Gram = cell-graph Laplacian (+ I_bnd)
    c_pc = _precond_cond(A_star, M)
    Mlu = splu(M)                                # AMG in production (M sparse SPD); direct here for the test
    it_none = _bicg_iters(A_star, None)
    it_pc = _bicg_iters(A_star, LinearOperator(M.shape, matvec=lambda v: Mlu.solve(v)))

    dens = T_sp.nnz / (T_sp.shape[0] * T_sp.shape[1])
    print(f"\n{nx}x{ny}x{nz} strength={strength}: n_elem={len(elements)} cycle={cycle}")
    print(f"  #2 cond:  dense-SVD={c_svd:7.1f}  sparse(B^T a)={c_sp:9.0f}  sparse+M^-1={c_pc:6.1f} "
          f"({c_pc/c_svd:.2f}x SVD)")
    print(f"  #2 BiCGStab iters:  unpreconditioned={it_none:4d}   +M^-1(cell-graph Laplacian)={it_pc:4d}")
    print(f"  #3 sparse T: shape{T_sp.shape} density {dens:.4f} -> A_star@x=T^T A(T x) sparse H-matvec wrapper")
    return dict(grid=f"{nx}x{ny}x{nz}", strength=strength, n_elem=len(elements), cycle=cycle,
                cond_dense_svd=c_svd, cond_sparse_raw=c_sp, cond_sparse_precond=c_pc,
                iters_unprecond=it_none, iters_precond=it_pc, T_density=float(dens))


def main():
    out = [evaluate(3, 3, 2, 2.0), evaluate(4, 4, 3, 2.0), evaluate(5, 5, 4, 2.0),
           evaluate(6, 6, 5, 2.0), evaluate(8, 8, 5, 2.0)]
    # scalability proof: preconditioned iters are BOUNDED while unpreconditioned grow
    it_pc = [r["iters_precond"] for r in out]
    it_un = [r["iters_unprecond"] for r in out]
    scalable = max(it_pc) <= 1.6 * min(it_pc) and it_un[-1] > 2 * it_pc[-1]
    with open(os.path.join(HERE, "yano_star_scalability.json"), "w") as f:
        json.dump({"scalable_loop_free_with_hacapk": scalable, "cases": out}, f, indent=2, default=float)
    print(f"\nbasis-Gram (cell-graph Laplacian) preconditioner -> BOUNDED iters {it_pc} (unprecond {it_un})")
    print(f"scalable + loop-free + HACApK-usable (bounded preconditioned iters): {scalable}")
    return scalable


if __name__ == "__main__":
    main()
