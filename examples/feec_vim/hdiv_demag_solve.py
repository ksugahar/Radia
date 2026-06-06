"""
hdiv_demag_solve.py -- the ACTUAL HDiv-type demag operator + the two tests, refined:
  (1) do the HDiv loops (ker of the charge map) lie in ker(N_demag)?  [exact]
  (2) is the STAR block Jacobi-solvable as mu_r grows?  [now in a LOOP-FREE BASIS, with
      multi-subpoint 1/r quadrature -- removes the projection-diagonal + coarse-centroid
      artifacts of the first cut]

N_ij = (1/4pi) Int Int charge_i(x) charge_j(x')/|x-x'|, charge=(-div M, M.n|bnd).  RT0 =>
piecewise-constant charge (value rho_i|el = Bv[el,i]/vol_el, sig_i|face = Bb[face,i]/area).
Each element -> 8 interior subpoints (4 per bnd face), constant charge, sub-cell self-patch.
Star solve: change to the loop-free basis S (right-singular vecs of the charge map Q with
nonzero sigma); A_sb = (1/chi) I - S^T N S; Jacobi(diag)-GMRES.

CONVERGENCE WARNING (verify-first, 2026-06-06): the crude centroid+ball-self-patch quadrature
here is NOT convergent -- refining 8->16 subpoints/hex (NSUB_SHELLS=2) makes cond/iters WORSE
(min|mu| 3.9e-3->2.1e-3, mu_r=1e4 iters 222->1245), because the per-subcell ball self-patch is
inconsistent across refinements.  => the cond / Jacobi-iter numbers from THIS script are
quadrature ARTIFACTS, NOT the true operator's conditioning.  Only the loop-null result
(||N.loop||/||N|| ~ 4e-16) is robust (quadrature-independent).  The TRUSTWORTHY conditioning
evidence is demag_spectrum_jacobi.py (Radia EXACT field, no quadrature: Jacobi 10-84 iters,
bounded).  A reliable HDiv star-block iter count needs proper Duffy/Graglia singular quadrature
(the real remaining build for the yano-vs-HDiv benchmark).
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
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width)).toarray()

with ng.TaskManager():
    mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=distort)
    fes = ng.HDiv(mesh, order=0)
    ndof = fes.ndof
    nn = ng.specialcf.normal(mesh.dim)
    L2v = ng.L2(mesh, order=0); L2b = ng.SurfaceL2(mesh, order=0)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u))*L2v.TestFunction()*ng.dx; bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace()*nn)*L2b.TestFunction()*ng.ds; bb.Assemble()
    Bv_m, Bb_m = csr(bv), csr(bb)
    Q = np.vstack([Bv_m, Bb_m])
    Uq, Sq, Vt = np.linalg.svd(Q)
    rankQ = int(np.sum(Sq > 1e-9*Sq.max()))
    n_loop = ndof - rankQ
    loops = Vt[rankQ:, :]                 # ker(Q) basis (n_loop x ndof), orthonormal rows
    Sstar = Vt[:rankQ, :].T               # star basis (ndof x rankQ), orthonormal cols
    print(f"distorted 3x3x3 hex, RT0: ndof={ndof}, star={rankQ}, loops={n_loop}")

    mass = ng.BilinearForm(L2v); mass += L2v.TrialFunction()*L2v.TestFunction()*ng.dx; mass.Assemble()
    massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction()*L2b.TestFunction()*ng.ds; massb.Assemble()
    el_vol = np.diag(csr(mass)); bf_area = np.diag(csr(massb))
    rho = Bv_m / el_vol[:, None]          # (n_el x ndof) constant rho_i per element
    sig = Bb_m / bf_area[:, None]         # (n_bf x ndof)

    # subpoints: 8 interior per (hex) element, 4 per (quad) bnd face; constant charge each
    pts, w, rself, Crows = [], [], [], []
    NSUB = int(os.environ.get("NSUB_SHELLS", "1"))     # 1 -> 8 pts/hex; 2 -> 16; convergence knob
    shells = [0.5] if NSUB == 1 else [0.4, 0.8]
    for k, el in enumerate(mesh.Elements(ng.VOL)):
        V = np.array([mesh[v].point for v in el.vertices]); c = V.mean(0)
        sub = np.vstack([c + f*(V - c) for f in shells]); ws = el_vol[k]/len(sub)
        for s in sub:
            pts.append(s); w.append(ws); rself.append(0.62*ws**(1/3.)); Crows.append(rho[k])
    for k, el in enumerate(mesh.Elements(ng.BND)):
        V = np.array([mesh[v].point for v in el.vertices]); c = V.mean(0)
        sub = np.vstack([c + f*(V - c) for f in shells]); ws = bf_area[k]/len(sub)
        for s in sub:
            pts.append(s); w.append(ws); rself.append(0.5*ws**0.5); Crows.append(sig[k])
    pts = np.array(pts); w = np.array(w); rself = np.array(rself); C = np.array(Crows)  # (nsp x ndof)

    D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    np.fill_diagonal(D, 0.0)
    Dinv = np.where(D > 0, 1.0/np.maximum(D, 1e-300), 0.0); np.fill_diagonal(Dinv, 1.0/rself)
    K = Dinv/(4*pi)
    WC = w[:, None]*C
    N = WC.T @ K @ WC                      # ndof x ndof symmetric demag-energy operator
    Nn = np.linalg.norm(N, 2)

    # (1) loops in ker(N) (exact: loops carry zero charge)
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / (Nn + 1e-300)
    print(f"(1) loops in ker(N): max ||N.loop||/||N|| = {loop_res:.2e}  "
          f"({'YES field-null by construction' if loop_res < 1e-10 else 'NO'})")

    # (2) star block in the LOOP-FREE BASIS: A_sb = (1/chi) I - S^T N S ; Jacobi-GMRES
    Nsb = Sstar.T @ N @ Sstar             # rankQ x rankQ
    mu = np.sort(np.abs(np.linalg.eigvals(Nsb).real))
    print(f"    star demag |mu|: min={mu[0]:.3e}, max={mu[-1]:.3e}")
    res = {"ndof": int(ndof), "n_loop": int(n_loop), "star": int(rankQ),
           "loop_res": float(loop_res), "mu_min": float(mu[0]), "mu_max": float(mu[-1])}
    def jac_gmres(A):
        d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
        Minv = LinearOperator(A.shape, matvec=lambda x: x/d)
        it = {"n": 0}
        gmres(A, np.ones(A.shape[0]), M=Minv, rtol=1e-8, maxiter=3000,
              callback=lambda xk: it.__setitem__("n", it["n"]+1), callback_type="pr_norm")
        return it["n"]
    for mur in (1e2, 1e3, 1e4):
        chi = mur - 1.0
        A = (1.0/chi)*np.eye(rankQ) - Nsb
        cond = np.linalg.cond(A); it = jac_gmres(A)
        print(f"    mu_r={mur:.0e}: cond(A_star)={cond:.2e}, Jacobi-GMRES iters={it}, "
              f"gap=(mu_min)/(1/chi)={mu[0]*chi:.1f}")
        res[f"mur_{mur:.0e}"] = {"cond": float(cond), "jac_iters": int(it)}

verdict = res["loop_res"] < 1e-10 and all(res[f"mur_{m:.0e}"]["jac_iters"] < 200 for m in (1e2,1e3,1e4))
print(f"\n=> HDiv-type: loops EXACTLY field-null + star Jacobi-GMRES bounded (<200) across mu_r: "
      f"{'YES' if verdict else 'see numbers'}")
res["verdict"] = bool(verdict)
with open(os.path.join(HERE, "hdiv_demag_solve.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "hdiv_demag_solve.json"))
