"""
hdiv_demag_solve.py -- assemble the ACTUAL HDiv-type demag operator (not the constant-M
proxy) and test the two halves together on the real HDiv space:
  (1) do the HDiv loops (ker of the charge map) lie in ker(N_demag)?  [exact: loops carry
      zero charge -> zero rows/cols, independent of quadrature accuracy]
  (2) is the STAR block Jacobi-solvable as mu_r grows?

N_ij = (1/4pi) Int Int charge_i(x) charge_j(x') / |x-x'| , charge = (-div M, M.n|bnd).
First cut: RT0 (HDiv order 0) so the charges are piecewise-constant (exact at element /
face centroids); the 1/r is centroid-centroid + a self-patch (r_self from cell size).
COARSE on the kernel (-> approximate conditioning numbers) but EXACT on the loop test.
"""
import json, os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import gmres, LinearOperator
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from math import sin, pi

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

def distort(x, y, z):
    return (x + 0.12*sin(pi*y)*z, y + 0.10*x*z, z + 0.08*x*sin(pi*y))

def csr(bf):
    m = bf.mat; r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width))

with ng.TaskManager():
    mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=distort)
    fes = ng.HDiv(mesh, order=0)                       # RT0: piecewise-constant charges
    ndof = fes.ndof
    nn = ng.specialcf.normal(mesh.dim)

    # ---- charge map Q (for the loop space): rows=charges, cols=HDiv dofs ----
    L2v = ng.L2(mesh, order=0); L2b = ng.SurfaceL2(mesh, order=0)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u))*L2v.TestFunction()*ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace()*nn)*L2b.TestFunction()*ng.ds; bb.Assemble()
    Bv_m = csr(bv).toarray(); Bb_m = csr(bb).toarray()
    Q = np.vstack([Bv_m, Bb_m])
    rankQ = np.linalg.matrix_rank(Q, tol=1e-9*max(1, np.linalg.norm(Q, 2)))
    n_loop = ndof - rankQ
    U, S, Vt = np.linalg.svd(Q)
    loops = Vt[rankQ:, :]                              # (n_loop x ndof) basis of ker(Q)
    print(f"distorted 3x3x3 hex, RT0: ndof={ndof}, star(rank Q)={rankQ}, loops(ker Q)={n_loop}")

    # ---- centroids (charge positions) + volumes/areas (mass diagonals) ----
    el_cent = np.array([np.array([mesh[v].point for v in el.vertices]).mean(axis=0)
                        for el in mesh.Elements(ng.VOL)])
    bf_cent = np.array([np.array([mesh[v].point for v in el.vertices]).mean(axis=0)
                        for el in mesh.Elements(ng.BND)])
    mass = ng.BilinearForm(L2v); mass += L2v.TrialFunction()*L2v.TestFunction()*ng.dx; mass.Assemble()
    massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction()*L2b.TestFunction()*ng.ds; massb.Assemble()
    el_vol = np.array(csr(mass).diagonal())            # |element| (L2-0 mass diag)
    bf_area = np.array(csr(massb).diagonal())          # |bnd face|

    # piecewise-constant charge of HDiv dof i: rho_i|elem_k = Bv[k,i]/vol_k ; sig_i|face_k = Bb[k,i]/area_k
    n_v, n_b = len(el_cent), len(bf_cent)
    C = np.zeros((n_v + n_b, ndof))
    C[:n_v, :] = Bv_m / el_vol[:, None]
    C[n_v:, :] = Bb_m / bf_area[:, None]
    pts = np.vstack([el_cent, bf_cent]); w = np.concatenate([el_vol, bf_area])

    # 1/r Gram with self-patch  K_pq = 1/(4pi |x_p-x_q|); K_pp = 1/(4pi r_self)
    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    r_self = np.empty(n_v + n_b)
    r_self[:n_v] = 0.62 * el_vol**(1/3.)               # ball-equiv radius of a cell
    r_self[n_v:] = 0.5 * np.sqrt(bf_area)
    np.fill_diagonal(D, 0.0); Dinv = np.where(D > 0, 1.0/np.maximum(D, 1e-300), 0.0)
    np.fill_diagonal(Dinv, 1.0/r_self)
    K = Dinv/(4*pi)
    WC = (w[:, None]) * C
    N = WC.T @ K @ WC                                  # ndof x ndof, symmetric demag-energy op

    # ---- (1) loops in ker(N)? ----
    Nn = np.linalg.norm(N, 2)
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / (Nn + 1e-300)
    print(f"(1) loops in ker(N): max ||N.loop||/||N|| = {loop_res:.2e}  "
          f"({'YES, field-null by construction' if loop_res < 1e-10 else 'no'})")

    # ---- (2) Jacobi on the star block: project out loops, solve A=(1/chi)I - N ----
    Pl = loops.T @ loops                               # projector onto loops (orthonormal rows from SVD)
    Ps = np.eye(ndof) - Pl                             # onto star
    def jac_iters(A):
        d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
        Minv = LinearOperator(A.shape, matvec=lambda x: x/d)
        it = {"n": 0}
        gmres(A, Ps @ np.ones(ndof), M=Minv, rtol=1e-8, maxiter=4000,
              callback=lambda xk: it.__setitem__("n", it["n"]+1), callback_type="pr_norm")
        return it["n"]
    res = {"ndof": int(ndof), "n_loop": int(n_loop), "loop_res": float(loop_res)}
    muN = np.sort(np.abs(np.linalg.eigvals(Ps @ N @ Ps).real))
    mu_star_min = muN[muN > 1e-9*muN.max()][0] if np.any(muN > 1e-9*muN.max()) else float('nan')
    print(f"    smallest star |mu| = {mu_star_min:.2e}")
    for mur in (1e2, 1e4):
        chi = mur - 1.0
        # solve only on star: A_star = Ps((1/chi)I - N)Ps + Pl (pin loops to identity to keep invertible)
        A = Ps @ ((1.0/chi)*np.eye(ndof) - N) @ Ps + Pl
        it = jac_iters(A)
        cond = np.linalg.cond(A)
        print(f"    mu_r={mur:.0e}: cond(A_star)={cond:.2e}, Jacobi-GMRES iters={it}, "
              f"gap=(star|mu|min)/(1/chi)={mu_star_min*chi:.1f}")
        res[f"mur_{mur:.0e}"] = {"cond": float(cond), "jac_iters": int(it)}
    res["mu_star_min"] = float(mu_star_min)

with open(os.path.join(HERE, "hdiv_demag_solve.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "hdiv_demag_solve.json"))
