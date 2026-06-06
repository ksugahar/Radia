"""
hdiv_demag_analytic_self.py -- HDiv-type demag with a CONVERGENT self term.

The previous crude self-patch (ball-equiv radius per subcell) was NON-convergent (refining
made cond/iters worse).  Fix: use the ANALYTIC self-energy of a uniform charge cell --
  cube  : Int Int_{V x V} 1/(4pi r) = c_cube  * V^{5/3} / (4pi),  c_cube  = <1/r>_unitcube
  square: Int Int_{A x A} 1/(4pi r) = c_sq    * A^{3/2} / (4pi),  c_sq    = <1/r>_unitsq
(c_cube, c_sq = the unit-cell mean reciprocal distance, computed once below).  Off-diagonal
(cell a != b) uses the centroid monopole meas_a meas_b/(4pi|c_a-c_b|).  The self term is now
FIXED (not refinement-dependent) so the operator is well-defined; accuracy vs the Radia-exact
proxy (demag_spectrum_jacobi.py) is the test, on a REGULAR mesh (cube self exact) and distorted.

N = B^T G B, B[cell,i] = charge value of HDiv dof i on cell (rho_i|el = Bv/V, sig_i|face = Bb/A);
G = charge-cell Coulomb Gram (analytic diagonal + centroid off-diagonal).  Loops from ker(charge
map).

RESULT (2026-06-07): use the NATURAL-basis Jacobi (diagonal of the full A in the HDiv-dof basis),
NOT the loop-free SVD basis (whose diagonal is meaningless -> 700 iters was that artifact).
  REGULAR mesh (cube self-energy is EXACT): natural-basis Jacobi-GMRES = 8 iters, mu_r-INDEPENDENT
    (1e2/1e3/1e4 all 8) -- matches the Radia-exact proxy. JACOBI HYPOTHESIS CONFIRMED on an accurate
    operator. Elegant: the exact de-Rham loops sit at EXACTLY eigenvalue 1/chi (perfect cluster,
    N.loop=4e-16), so GMRES absorbs them in ~1 iter -- NO explicit loop-star needed, plain
    Jacobi-GMRES on the full system works.
  DISTORTED mesh: iters explode (4200->37720->60000) -- a QUADRATURE ARTIFACT, not a breakdown:
    c_cube*V^{5/3} is exact for CUBES but WRONG for distorted hexes -> inaccurate diagonal -> bad
    Jacobi. Needs an accurate distorted-hex self-energy (Duffy), the remaining quadrature task.
"""
import json, os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import gmres, LinearOperator
from scipy import integrate
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from math import sin, pi

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- unit-cell mean reciprocal distance (computed once, reduced integrals) ----
#   <1/r>_unitsquare = 4 Int_0^1 Int_0^1 (1-u)(1-v)/sqrt(u^2+v^2) du dv
c_sq, _ = integrate.dblquad(lambda v, u: (1-u)*(1-v)/np.sqrt(u*u+v*v+1e-300),
                            0, 1, 0, 1)
c_sq *= 4.0
#   <1/r>_unitcube = 8 Int (1-u)(1-v)(1-w)/sqrt(u^2+v^2+w^2)
c_cube, _ = integrate.tplquad(lambda w, v, u: (1-u)*(1-v)*(1-w)/np.sqrt(u*u+v*v+w*w+1e-300),
                              0, 1, 0, 1, 0, 1)
c_cube *= 8.0
print(f"unit-cell mean recip dist: c_cube={c_cube:.5f} (lit ~1.88231), c_sq={c_sq:.5f}")

def csr(bf):
    m = bf.mat; r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width)).toarray()

def run(distort_amt, tag):
    def distort(x, y, z):
        d = distort_amt
        return (x + d*sin(pi*y)*z, y + 0.83*d*x*z, z + 0.67*d*x*sin(pi*y))
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=distort)
        fes = ng.HDiv(mesh, order=0); ndof = fes.ndof
        nn = ng.specialcf.normal(mesh.dim)
        L2v = ng.L2(mesh, order=0); L2b = ng.SurfaceL2(mesh, order=0)
        u = fes.TrialFunction()
        bv = ng.BilinearForm(trialspace=fes, testspace=L2v); bv += (-ng.div(u))*L2v.TestFunction()*ng.dx; bv.Assemble()
        bb = ng.BilinearForm(trialspace=fes, testspace=L2b); bb += (u.Trace()*nn)*L2b.TestFunction()*ng.ds; bb.Assemble()
        Bv_m, Bb_m = csr(bv), csr(bb)
        Q = np.vstack([Bv_m, Bb_m]); Uq, Sq, Vt = np.linalg.svd(Q)
        rankQ = int(np.sum(Sq > 1e-9*Sq.max())); n_loop = ndof - rankQ
        loops = Vt[rankQ:, :]; Sstar = Vt[:rankQ, :].T
        mass = ng.BilinearForm(L2v); mass += L2v.TrialFunction()*L2v.TestFunction()*ng.dx; mass.Assemble()
        massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction()*L2b.TestFunction()*ng.ds; massb.Assemble()
        el_vol = np.diag(csr(mass)); bf_area = np.diag(csr(massb))
        el_c = np.array([np.array([mesh[v].point for v in el.vertices]).mean(0) for el in mesh.Elements(ng.VOL)])
        bf_c = np.array([np.array([mesh[v].point for v in el.vertices]).mean(0) for el in mesh.Elements(ng.BND)])

    n_el, n_bf = len(el_c), len(bf_c)
    cent = np.vstack([el_c, bf_c]); meas = np.concatenate([el_vol, bf_area])
    B = np.vstack([Bv_m / el_vol[:, None], Bb_m / bf_area[:, None]])      # (n_cell x ndof)
    # charge-cell Coulomb Gram G
    Dd = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2)
    np.fill_diagonal(Dd, np.inf)
    G = (meas[:, None]*meas[None, :]) / (4*pi*Dd)                          # centroid monopole off-diag
    diagG = np.empty(n_el + n_bf)
    diagG[:n_el] = c_cube * el_vol**(5/3.) / (4*pi)                        # cube self (exact if regular)
    diagG[n_el:] = c_sq * bf_area**(1.5) / (4*pi)                          # square self
    np.fill_diagonal(G, diagG)
    N = B.T @ G @ B
    Nn = np.linalg.norm(N, 2)
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / (Nn + 1e-300)
    Nsb = Sstar.T @ N @ Sstar
    mu = np.sort(np.abs(np.linalg.eigvals(Nsb).real))
    def jac(A):
        d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
        it = {"n": 0}
        gmres(A, np.ones(A.shape[0]), M=LinearOperator(A.shape, matvec=lambda x: x/d),
              rtol=1e-8, maxiter=3000, callback=lambda xk: it.__setitem__("n", it["n"]+1),
              callback_type="pr_norm")
        return it["n"]
    print(f"\n[{tag}] ndof={ndof} star={rankQ} loops={n_loop}: loop_res={loop_res:.1e}, "
          f"star|mu| min={mu[0]:.3e} max={mu[-1]:.3e}")
    out = {"tag": tag, "loop_res": float(loop_res), "mu_min": float(mu[0])}
    for mur in (1e2, 1e3, 1e4):
        chi = mur - 1.0
        A_sb = (1.0/chi)*np.eye(rankQ) - Nsb            # loop-free SVD basis
        A_full = (1.0/chi)*np.eye(ndof) - N             # natural HDiv-dof basis (loops included)
        it_sb = jac(A_sb); it_nat = jac(A_full)         # Jacobi = diagonal in each basis
        cond = np.linalg.cond(A_sb)
        print(f"    mu_r={mur:.0e}: cond(star)={cond:.2e}, Jacobi iters: SVD-basis={it_sb}, "
              f"NATURAL-basis(full A)={it_nat}, gap={mu[0]*chi:.1f}")
        out[f"mur_{mur:.0e}"] = {"cond": float(cond), "jac_iters_svd": int(it_sb),
                                 "jac_iters_natural": int(it_nat)}
    return out

res = {"c_cube": float(c_cube), "c_sq": float(c_sq),
       "regular": run(0.0, "regular"), "distorted": run(0.18, "distorted")}
print("\n(regular: cube self-energy is EXACT -> conditioning should match the Radia-exact proxy"
      " ~ Jacobi 10-84 iters if the assembly is right.)")
with open(os.path.join(HERE, "hdiv_demag_analytic_self.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "hdiv_demag_analytic_self.json"))
