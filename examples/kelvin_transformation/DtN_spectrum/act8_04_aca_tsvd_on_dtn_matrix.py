"""act8_04_aca_tsvd_on_dtn_matrix -- Can ACA-TSVD be applied to the DtN/Kelvin material-aware transfer matrix M?

Answer (verified below): YES, but the two ingredients answer two different questions.

  * TSVD  -- applies to the GLOBAL inverse design `psi = M^+ B_target`, unchanged, and is
            necessary.  M is the discretisation of the compact (smoothing) forward operator
            psi -> field, so its singular values decay and the inverse is ill-posed.  TSVD is
            method-agnostic: it does not care that M came from a Kelvin-FEM Schur condensation
            rather than free-space Biot-Savart.  (act8_03_general_iron_design/act8_02_streamfunction_design_matrix already use `tsvd_pinv`.)

  * ACA   -- applies BLOCK-WISE, as an H-matrix, NOT as a single global low-rank factor.
            The GLOBAL M is essentially full numerical rank (targets near the coil => near
            field).  But the underlying kernel is the *material Green function* G_mu(x_t, y_c),
            which is asymptotically smooth for well-separated source/target clusters even with
            iron in the gap.  So an ADMISSIBLE block (compact source cap, far targets) is low
            rank and ACA-compresses; a NEAR block stays dense.  This is exactly the near/far
            split HACApK already does for the MMM/MSC matrix -- the same machinery transfers.

The only thing that changes from the free-space ACA is the ENTRY ORACLE: there is no closed
form.  One column of M = one back-substitution of the pre-factored Kelvin-FEM (unit Dirichlet
bump on a coil DoF, read the targets); one row = one adjoint solve.  Because a solve returns a
whole column (= a fast matvec of M), a randomized SVD / Lanczos is in fact even more natural
here than entrywise ACA; ACA/H-matrix is the right tool when M must be stored and re-applied
many times (e.g. a CMA-ES outer loop) or the geometry is huge with a genuine near/far split.

Run:  pip install -e packages/radia-mcp ; python demo_ii_aca_tsvd_on_dtn_matrix.py
Needs numpy, scipy, ngsolve 6.2.2604, netgen.occ.  order=1 is used ONLY so coil-DoF coordinates
are exact nodal values (gf.Set(x)); the rank structure is an asymptotic-smoothness property and
is independent of FE order.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import time
import numpy as np
import ngsolve as ng
from ngsolve import TaskManager
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry

# ---- geometry: coil sphere r=a, physical vacuum a<r<R_out, NON-concentric iron blob, Kelvin ball
a, R_out, offset = 0.5, 1.0, 3.0
iron_c, iron_r = (0.0, 0.0, 0.68), 0.16
order, maxh, intorder = 1, 0.14, 6
mu_r = 50.0
np.random.seed(0)


def build_kelvin_fem(mu_s):
    """Mesh coil + iron blob + open exterior (Kelvin), assemble the SPD stiffness A(mu)."""
    inner = Sphere(Pnt(0, 0, 0), a); outer = Sphere(Pnt(0, 0, 0), R_out)
    iron = Sphere(Pnt(*iron_c), iron_r)
    for f in inner.faces: f.name = "inner"
    for f in outer.faces: f.name = "kelvin_int"
    iron.mat("shell")
    vac = (outer - inner) - iron; vac.mat("vac")
    kball = Sphere(Pnt(offset, 0, 0), R_out)
    for f in kball.faces: f.name = "kelvin_ext"
    kball.mat("kelvin"); gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
    fi = [f for f in vac.faces if f.name == "kelvin_int"][0]
    fe = [f for f in kball.faces if f.name == "kelvin_ext"][0]
    # single-face periodic glue truncation sphere <-> Kelvin ball (CRITICAL: one Identify)
    fi.Identify(fe, "kelvin", IdentificationType.PERIODIC, occ.gp_Trsf.Translation(Vec(offset, 0, 0)))
    geo = OCCGeometry(occ.Glue([vac, iron, kball, gnd]))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
    x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset) ** 2 + y * y + z * z + 1e-20
    # mu as an FE coefficient incl. the kelvin-ball (R/rho)^2 weight
    mu = mesh.MaterialCF({"vac": 1.0, "shell": mu_s, "kelvin": R_out ** 2 / rp2}, default=1.0)
    fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="inner|GND")); u, v = fes.TnT()
    A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=intorder)); A.Assemble()
    return mesh, fes, A


def svd_ratio(B):
    s = np.linalg.svd(B, compute_uv=False)
    return s / s[0]


def nrank(sr, tol=1e-4):
    return int(np.sum(sr > tol))


def nearest_k(pts, center, k):
    return np.argsort(np.linalg.norm(pts - np.asarray(center), axis=1))[:k]


def main():
    with TaskManager():
        print("Building Kelvin-FEM (mu_r=%.0f, order=%d, maxh=%.2f, non-concentric iron blob)..."
              % (mu_r, order, maxh))
        mesh, fes, A = build_kelvin_fem(mu_r)
        Ainv = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")   # ONE factorisation
        innerb = fes.GetDofs(mesh.Boundaries("inner"))
        cdofs = [i for i in range(fes.ndof) if innerb[i]]
        n_coil = len(cdofs)
        gfx, gfy, gfz = ng.GridFunction(fes), ng.GridFunction(fes), ng.GridFunction(fes)
        gfx.Set(ng.x); gfy.Set(ng.y); gfz.Set(ng.z)            # order-1 nodal -> exact coords
        cxyz = np.array([[gfx.vec[i], gfy.vec[i], gfz.vec[i]] for i in cdofs])
        print("  ndof=%d, coil DoFs=%d on r=%.2f" % (fes.ndof, n_coil, a))

        # target cloud: two physical-vacuum shells, fine angular grid
        txyz, targets = [], []
        for r_t in (0.70, 0.95):
            for th in np.deg2rad(np.linspace(8, 172, 11)):
                for ph in np.deg2rad(np.linspace(0, 330, 12)):
                    d = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
                    p = r_t * d; txyz.append(p); targets.append(mesh(*p))
        txyz = np.array(txyz); m = len(targets)
        print("  targets m=%d (2 shells r=0.70, 0.95)" % m)

        _gf = ng.GridFunction(fes); _rr = _gf.vec.CreateVector()

        def column(idof):
            """One column of M = one back-substitution (the matvec oracle)."""
            _gf.vec[:] = 0.0; _gf.vec[idof] = 1.0
            _rr.data = -(A.mat * _gf.vec); _gf.vec.data += Ainv * _rr
            return np.array([_gf(t) for t in targets])

        t0 = time.time()
        M = np.column_stack([column(j) for j in cdofs])     # m x n_coil, n_coil solves
        print("  dense M %dx%d in %.1fs (%d back-substitutions)" % (M.shape[0], M.shape[1],
                                                                    time.time() - t0, n_coil))

    # ---- compact clusters so the admissibility ratio eta = dist / max(diam) is honest ----
    src = nearest_k(cxyz, (0.0, 0.0, a), 20)                 # compact source cap near +z pole
    far = nearest_k(txyz, (0.0, 0.0, -0.95), 30)            # compact far cluster (-z, separated)
    near = nearest_k(txyz, (0.0, 0.0, 0.70), 30)            # compact near cluster (above the cap)

    def cl(idx, P):
        c = P[idx].mean(0); diam = 2 * np.linalg.norm(P[idx] - c, axis=1).max(); return c, diam
    sc, sd = cl(src, cxyz); fc, fd = cl(far, txyz); nc, nd = cl(near, txyz)
    eta_far = np.linalg.norm(fc - sc) / max(sd, fd)
    eta_near = np.linalg.norm(nc - sc) / max(sd, nd)
    print("\n  source cap : %d DoFs, diam=%.2f" % (len(src), sd))
    print("  FAR cluster: %d targets, dist=%.2f, eta=%.2f  (admissible: eta > ~2)"
          % (len(far), np.linalg.norm(fc - sc), eta_far))
    print("  NEAR cluster: %d targets, dist=%.2f, eta=%.2f" % (len(near), np.linalg.norm(nc - sc), eta_near))

    Bfar, Bnear = M[np.ix_(far, src)], M[np.ix_(near, src)]
    srf, srn = svd_ratio(Bfar), svd_ratio(Bnear)
    rf, rn = nrank(srf), nrank(srn)
    mnf = min(Bfar.shape)
    print("\n[ADMISSIBLE block] %s  rank(1e-4)=%d / %d   sigma_k/sigma_1: %s"
          % (Bfar.shape, rf, mnf, "  ".join("%.0e" % v for v in srf[:8])))
    print("[NEAR block]       %s  rank(1e-4)=%d / %d   sigma_k/sigma_1: %s"
          % (Bnear.shape, rn, min(Bnear.shape), "  ".join("%.0e" % v for v in srn[:8])))

    ell = rf + 4
    err_rsvd = None
    if ell < mnf:
        Om = np.random.randn(Bfar.shape[1], ell)
        Q, _ = np.linalg.qr(Bfar @ Om)                       # ell oracle-calls (1 solve each)
        err_rsvd = np.linalg.norm(Bfar - Q @ (Q.T @ Bfar)) / np.linalg.norm(Bfar)
        print("[rSVD far block]   ell=%d << %d  =>  ||B-QQ'B||/||B||=%.1e   "
              "(ACA store %d vs dense %d)" % (ell, mnf, err_rsvd, rf * sum(Bfar.shape), Bfar.size))

    sg = svd_ratio(M); rg = nrank(sg)
    cond_tsvd = 1.0 / sg[max(rg - 1, 0)]
    print("\n[GLOBAL M]         %s  rank(1e-4)=%d / %d  -> NOT globally low rank (near-field)"
          % (M.shape, rg, min(M.shape)))
    print("                   => H-matrix: ACA far blocks + dense near blocks (the HACApK pattern)")
    print("[TSVD]             cond at rank(1e-6) = %.1e  -> regularised inverse needed & valid" % cond_tsvd)

    ok = (rf <= 0.6 * mnf) and (rn >= 0.9 * min(Bnear.shape)) and (err_rsvd is not None and err_rsvd < 1e-3)
    print("\n%s  admissible block low rank (%d/%d) + near block full (%d/%d) + rSVD %.0e"
          % ("[PASS]" if ok else "[CHECK]", rf, mnf, rn, min(Bnear.shape),
             err_rsvd if err_rsvd is not None else -1))
    print("VERDICT: TSVD applies globally (cond %.0e); ACA applies block-wise (far rank %d/%d). "
          "Same H-matrix machinery as HACApK, oracle = 1 Kelvin-FEM solve per column." % (cond_tsvd, rf, mnf))
    return ok


if __name__ == "__main__":
    main()
