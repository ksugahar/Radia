"""
hdiv_demag_quad_self.py -- HDiv-type demag with an ACCURATE distorted-hex self-energy.

Closes the quadrature artifact in hdiv_demag_analytic_self.py:  there the cell self-energy
used the CUBE formula  G_aa = c_cube * V^{5/3}/(4pi)  which is exact for a cube but WRONG for
a distorted hex -> inaccurate diag(N) -> bad Jacobi (distorted exploded 4200->60000 iters).

Fix: compute G_aa = Int Int_{cell x cell} 1/(4pi|x-x'|) by a CONVERGENT subdivision of the
ACTUAL distorted geometry (trilinear hex map / bilinear quad map, NGSolve vertex order probed):
  cell -> nsub^3 (hex) / nsub^2 (quad) sub-cells (centers + |detJ| measures from the map);
  G_aa ~ Sum_{i!=j} meas_i meas_j/(4pi |c_i-c_j|)         [cross: smooth, captures the SHAPE]
       + Sum_i  c_cube * meas_i^{5/3}/(4pi)                [sub-self: cube formula on the SMALL
                                                            sub-cell -> O(nsub^-2) -> vanishes]
As nsub grows the sub-self -> 0 and the cross -> the true 6D integral, so G_aa converges to the
true self-energy of the distorted hex (the shape enters through the sub-cell centroid cloud).
Validated: regular cube -> c_cube*V^{5/3}, regular square -> c_sq*A^{3/2} (the analytic limits);
convergence in nsub shown.  Off-diagonal cell-pairs stay centroid-monopole (already gave 8 iters
on the regular mesh).  Then natural-basis Jacobi-GMRES iters on REGULAR + DISTORTED.

RESULT (2026-06-07) -- TWO artifacts found, both fixed, the true picture is HONEST and strong:
  (a) self-energy quadrature: the cube formula was only ~0.6% off even at 18% distortion (the
      subdivided self-energy converges in nsub: 6.075e-4 -> 6.052e-4).  MINOR, not the cause.
  (b) THE real artifact: scipy gmres defaults to restart=20.  A 108-dof system needing >20
      Krylov steps stagnates under GMRES(20) -> the bogus 80000-iter "explosion".  A = (1/chi)I-N
      is SYMMETRIC INDEFINITE (N=B^T G B sym PSD), so the right method is MINRES (|diag| SPD
      precond, no restart) or FULL (restart=ndof) GMRES -> the honest Krylov count.

  Corrected, honest result (Jacobi MINRES == full GMRES):
    mu_r (1e2->1e4): iters FLAT -- regular 8/8/8, distorted ~80/76/69 (NOT growing, slightly
      DECREASING).  The yano-type high-mu BREAKDOWN IS ABSENT: loops are EXACTLY field-null
      (||N.loop||/||N|| ~ 4e-16) so they sit at exactly eigenvalue 1/chi and never pollute.
      => plain Jacobi gives mu_r-INDEPENDENT convergence.  NO explicit loop-star, NO H-ILU
         needed for the high-mu regime -- the user's "Jacobi de sumu" hope HOLDS for mu_r.
    distortion (0->0.24, mu_r=1e4 fixed): iters GROW 8->57->61->68->82 (star-block cond
      71->1e5).  This is a MESH-QUALITY conditioning effect (a fixed-mu_r property), NOT a
      mu_r blowup.  Plain POINT-Jacobi is distortion-sensitive; a block-Jacobi (near off-diag)
      or AMS-type preconditioner would cut this -- but it is orthogonal to the mu_r breakdown.
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
EPS = 1e-5

# ---- analytic unit-cell mean reciprocal distance (the regular-limit reference + sub-self const) -
c_sq, _ = integrate.dblquad(lambda v, u: (1-u)*(1-v)/np.sqrt(u*u+v*v+1e-300), 0, 1, 0, 1)
c_sq *= 4.0
c_cube, _ = integrate.tplquad(lambda w, v, u: (1-u)*(1-v)*(1-w)/np.sqrt(u*u+v*v+w*w+1e-300),
                              0, 1, 0, 1, 0, 1)
c_cube *= 8.0
print(f"unit-cell mean recip dist: c_cube={c_cube:.5f} (lit ~1.88231), c_sq={c_sq:.5f}")

# ---- geometry maps in the PROBED NGSolve vertex order ----------------------------------------
#   hex local (xi=x, eta=y, zeta=z): v0(000) v1(001) v2(011) v3(010) v4(100) v5(101) v6(111) v7(110)
def hex_map(V, L):
    xi, eta, ze = L[..., 0], L[..., 1], L[..., 2]
    x0, x1, y0, y1, z0, z1 = 1-xi, xi, 1-eta, eta, 1-ze, ze
    N = np.stack([x0*y0*z0, x0*y0*z1, x0*y1*z1, x0*y1*z0,
                  x1*y0*z0, x1*y0*z1, x1*y1*z1, x1*y1*z0], axis=-1)   # (...,8)
    return N @ V                                                      # (...,3)

#   quad local (u,v): v0(00) v1(01) v2(11) v3(10) -- consecutive corners
def quad_map(V4, L):
    u, v = L[..., 0], L[..., 1]
    N = np.stack([(1-u)*(1-v), (1-u)*v, u*v, u*(1-v)], axis=-1)
    return N @ V4

def _grid(nsub, dim):
    g = (np.arange(nsub) + 0.5) / nsub
    return np.stack(np.meshgrid(*([g]*dim), indexing="ij"), -1).reshape(-1, dim)

def hex_self_energy(V, nsub):
    L = _grid(nsub, 3); C = hex_map(V, L); h = 1.0/nsub
    J = np.empty(L.shape[:-1] + (3, 3))
    for d in range(3):
        e = np.zeros(3); e[d] = EPS
        J[..., d] = (hex_map(V, L+e) - hex_map(V, L-e)) / (2*EPS)
    meas = np.abs(np.linalg.det(J)) * h**3
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2); np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(meas, meas) / (4*pi*D))
    selfsub = np.sum(c_cube * meas**(5/3.) / (4*pi))
    return cross + selfsub, meas.sum()

def quad_self_energy(V4, nsub):
    L = _grid(nsub, 2); C = quad_map(V4, L); h = 1.0/nsub
    eu = np.array([EPS, 0.0]); ev = np.array([0.0, EPS])
    xu = (quad_map(V4, L+eu) - quad_map(V4, L-eu)) / (2*EPS)
    xv = (quad_map(V4, L+ev) - quad_map(V4, L-ev)) / (2*EPS)
    meas = np.linalg.norm(np.cross(xu, xv), axis=1) * h**2
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2); np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(meas, meas) / (4*pi*D))
    selfsub = np.sum(c_sq * meas**(1.5) / (4*pi))
    return cross + selfsub, meas.sum()

def csr(bf):
    m = bf.mat; r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width)).toarray()

def jac_iters(A):
    """A = (1/chi)I - N is SYMMETRIC INDEFINITE (N = B^T G B is sym PSD).  Use MINRES with the
    SPD |diag| preconditioner (no restart -> honest Krylov iteration count).  Cross-check with
    FULL (non-restarted) GMRES so a restart-stagnation artifact cannot masquerade as a breakdown.
    """
    from scipy.sparse.linalg import minres
    n = A.shape[0]; b = np.ones(n)
    d = np.abs(np.diag(A)).copy(); d[d < 1e-30] = 1.0
    Minv = LinearOperator((n, n), matvec=lambda x: x/d)
    itm = {"n": 0}
    minres(A, b, M=Minv, rtol=1e-8, maxiter=4000,
           callback=lambda xk: itm.__setitem__("n", itm["n"]+1))
    # full GMRES: restart = n (one Krylov space, no restart) -> true count <= n
    ds = np.diag(A).copy(); ds[np.abs(ds) < 1e-30] = 1.0
    itg = {"n": 0}
    gmres(A, b, M=LinearOperator((n, n), matvec=lambda x: x/ds), rtol=1e-8,
          restart=n, maxiter=2, callback=lambda xk: itg.__setitem__("n", itg["n"]+1),
          callback_type="pr_norm")
    return itm["n"], itg["n"]

def build(distort_amt, nsub_hex, nsub_quad, n=3):
    def distort(x, y, z):
        d = distort_amt
        return (x + d*sin(pi*y)*z, y + 0.83*d*x*z, z + 0.67*d*x*sin(pi*y))
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n, mapping=distort)
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
        # sparse HDiv FEM operators (for operator-preconditioning / H-AMS feasibility tests)
        vh = fes.TestFunction()
        am = ng.BilinearForm(fes); am += u*vh*ng.dx; am.Assemble()
        ad = ng.BilinearForm(fes); ad += ng.div(u)*ng.div(vh)*ng.dx; ad.Assemble()
        mass_hdiv = csr(am); divdiv_hdiv = csr(ad)
        el_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.VOL)]
        bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
        el_dofs = [np.array([dd for dd in fes.GetDofNrs(el) if dd >= 0]) for el in mesh.Elements(ng.VOL)]
    el_c = np.array([V.mean(0) for V in el_V]); bf_c = np.array([V.mean(0) for V in bf_V])
    n_el, n_bf = len(el_c), len(bf_c)
    cent = np.vstack([el_c, bf_c]); meas = np.concatenate([el_vol, bf_area])
    B = np.vstack([Bv_m / el_vol[:, None], Bb_m / bf_area[:, None]])

    # off-diagonal: centroid monopole; diagonal: ACCURATE subdivided self-energy
    Dd = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2); np.fill_diagonal(Dd, np.inf)
    G = (meas[:, None]*meas[None, :]) / (4*pi*Dd)
    diagG = np.empty(n_el + n_bf); volerr = 0.0; areaerr = 0.0
    for k, V in enumerate(el_V):
        diagG[k], vsum = hex_self_energy(V, nsub_hex); volerr = max(volerr, abs(vsum-el_vol[k])/el_vol[k])
    for k, V in enumerate(bf_V):
        diagG[n_el+k], asum = quad_self_energy(V, nsub_quad); areaerr = max(areaerr, abs(asum-bf_area[k])/bf_area[k])
    np.fill_diagonal(G, diagG)
    N = B.T @ G @ B; Nn = np.linalg.norm(N, 2)
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / (Nn + 1e-300)
    Nsb = Sstar.T @ N @ Sstar; mu = np.sort(np.abs(np.linalg.eigvals(Nsb).real))
    return dict(ndof=ndof, rankQ=rankQ, n_loop=n_loop, N=N, Nsb=Nsb, mu=mu, el_dofs=el_dofs,
                Sstar=Sstar, loops=loops, mass_hdiv=mass_hdiv, divdiv_hdiv=divdiv_hdiv,
                loop_res=loop_res, volerr=volerr, areaerr=areaerr,
                cube_self=diagG[0], cube_self_cubeformula=c_cube*el_vol[0]**(5/3.)/(4*pi), el_vol0=el_vol[0])

def run(distort_amt, tag, nsub_hex, nsub_quad):
    d = build(distort_amt, nsub_hex, nsub_quad)
    print(f"\n[{tag}] ndof={d['ndof']} star={d['rankQ']} loops={d['n_loop']} "
          f"(nsub hex={nsub_hex} quad={nsub_quad}): subdiv vol-err={d['volerr']:.1e} area-err={d['areaerr']:.1e}")
    print(f"    cell0 self-energy: subdivided={d['cube_self']:.6e}  cube-formula={d['cube_self_cubeformula']:.6e}  "
          f"ratio={d['cube_self']/d['cube_self_cubeformula']:.4f}  (==1 iff regular)")
    print(f"    loop_res={d['loop_res']:.1e}  star|mu| min={d['mu'][0]:.3e} max={d['mu'][-1]:.3e}")
    out = {"tag": tag, "ndof": d["ndof"], "n_loop": d["n_loop"], "loop_res": float(d["loop_res"]),
           "vol_err": float(d["volerr"]), "area_err": float(d["areaerr"]), "mu_min": float(d["mu"][0]),
           "self_ratio_cell0": float(d["cube_self"]/d["cube_self_cubeformula"])}
    for mur in (1e2, 1e3, 1e4):
        chi = mur - 1.0
        A_full = (1.0/chi)*np.eye(d["ndof"]) - d["N"]
        it_minres, it_gmres = jac_iters(A_full)
        cond = np.linalg.cond((1.0/chi)*np.eye(d["rankQ"]) - d["Nsb"])
        print(f"    mu_r={mur:.0e}: Jacobi MINRES iters={it_minres}, full-GMRES iters={it_gmres}, "
              f"cond(star)={cond:.2e}, gap={d['mu'][0]*chi:.1f}")
        out[f"mur_{mur:.0e}"] = {"jac_iters_minres": int(it_minres), "jac_iters_gmres": int(it_gmres),
                                 "cond": float(cond)}
    return out

if __name__ == "__main__":
    # convergence in nsub on the DISTORTED mesh (the artifact test)
    print("\n=== self-energy convergence in nsub (distorted cell0) ===")
    def distort18(x, y, z):
        return (x + 0.18*sin(pi*y)*z, y + 0.83*0.18*x*z, z + 0.67*0.18*x*sin(pi*y))
    with ng.TaskManager():
        mz = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=distort18)
        V0 = np.array([mz[v].point for v in list(mz.Elements(ng.VOL))[0].vertices])
    prev = None
    for ns in (3, 5, 8, 12, 16):
        e, vs = hex_self_energy(V0, ns)
        chg = "" if prev is None else f"  d={abs(e-prev)/e:.1e}"
        print(f"  nsub={ns:2d}: self-energy={e:.6e}  sub-vol-sum={vs:.6e}{chg}")
        prev = e

    res = {"c_cube": float(c_cube), "c_sq": float(c_sq),
           "regular": run(0.0, "regular", 6, 8),
           "distorted": run(0.18, "distorted", 8, 10)}

    # separate the two effects: iters vs DISTORTION (fixed mu_r=1e4) vs iters vs mu_r (above, flat)
    print("\n=== Jacobi-MINRES iters vs mesh distortion (mu_r=1e4 fixed) ===")
    sweep = []
    for da in (0.0, 0.06, 0.12, 0.18, 0.24):
        d = build(da, 6, 8); chi = 1e4 - 1.0
        itm, itg = jac_iters((1.0/chi)*np.eye(d["ndof"]) - d["N"])
        cond = np.linalg.cond((1.0/chi)*np.eye(d["rankQ"]) - d["Nsb"])
        print(f"  distort={da:.2f}: MINRES iters={itm:3d}, cond(star)={cond:.2e}, "
              f"star|mu|_min={d['mu'][0]:.2e}, loop_res={d['loop_res']:.1e}")
        sweep.append({"distort": da, "iters_minres": int(itm), "cond": float(cond),
                      "mu_min": float(d["mu"][0]), "loop_res": float(d["loop_res"])})
    res["distortion_sweep_mur1e4"] = sweep
    print("  => iters FLAT in mu_r (above) but GROW with distortion -> a mesh-conditioning effect,")
    print("     NOT the yano-type high-mu breakdown (which grows with mu_r). loops stay field-null.")
    with open(os.path.join(HERE, "hdiv_demag_quad_self.json"), "w") as f:
        json.dump({k: v for k, v in res.items() if k in ("c_cube", "c_sq", "regular", "distorted")}, f, indent=2)
    print("\nsaved", os.path.join(HERE, "hdiv_demag_quad_self.json"))
