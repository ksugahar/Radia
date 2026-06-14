"""demo_kk -- wire the cheap DtN material-aware M into the EXISTING stream-function solver.

This closes the Track-B loop: the DtN/Kelvin-FEM supplies the material-aware transfer matrix M
(cheap to build -- see bench_dtn_mbuild), and the *existing*, kernel-agnostic
``radia.stream_function`` (ACA+)+TSVD + ``RegularizedTSVD`` ridge solver consumes it UNCHANGED.
Nothing in the design code is DtN-specific: the free-space Biot-Savart kernel is swapped for the
DtN material Green operator just by changing the ``entry(i, j)`` callback.

What is demonstrated (all verified below):
  [1] radia.stream_function.aca_tsvd compresses the FORMED DtN M (entry = array lookup; the cheap-
      build path the matrix-construction benchmark justified) -- reproduces M to ~1e-14.
  [2] The SAME solver designs the coil: plain TSVD (min-L2) vs the RIDGE
      (RegularizedTSVD.from_stiffness, S = the surface current-density seminorm = "電流密度最小化").
      Both HIT the (iron-system) target; the ridge does so with LOWER current density ||K||.
  [3] An L-curve (alpha sweep) trading field accuracy against current density -- one cached
      factorisation, only the small k x k core re-solved per alpha.
  [4] The same solver fed the FREE-SPACE kernel designs a coil that MISSES in the iron system --
      so the material-aware DtN kernel is what makes the design correct.

The point: with M cheap to build (the benchmark), forming it and handing it to the existing
ACA+TSVD/ridge stream-function solver is the clean production path (the user's observation:
"行列がある状態で (ACA+)+TSVD は可能。構築コストが安いなら ACA+TSVD でよい").  ACA+ acts on the
formed matrix; a randomized-SVD oracle (demo_jj) is the alternative when M is not formed.

Run:  pip install -e packages/radia-mcp ; python demo_kk_streamfunction_ridge_with_dtn.py
Needs numpy, scipy, ngsolve 6.2.2604, netgen.occ, radia (radia.stream_function).
"""
import os, sys
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.sparse as sp
import ngsolve as ng
from ngsolve import TaskManager
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import demo_hh_general_iron_design as hh   # build_fem / coil_dofs / targets_on_sphere / build_transfer / fresh_forward
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, RegularizedTSVD

BLOB = ("blob", (0.0, 0.0, 0.68), 0.16)    # non-concentric iron: no closed-form Green function
MU_R, ORDER, MAXH = 50.0, 2, 0.22
ALPHAS = [0.0, 1e-6, 1e-4, 1e-2, 1e-1]


def surface_seminorm(fes, mesh, mass=1e-6):
    """S = the coil surface current-density seminorm  int_inner |grad_S psi|^2 + mass*int psi^2.
    grad_S = the TANGENTIAL trace gradient (the genuine 2D surface gradient of the coil trace;
    the normal trace component is projected out).  psi^T S psi = ||K||^2 with K = n x grad_S psi.
    The small mass term lifts the constant null mode so S is SPD (RegularizedTSVD requires it)."""
    u, v = fes.TnT()
    nrm = ng.specialcf.normal(mesh.dim)
    gu, gv = ng.grad(u).Trace(), ng.grad(v).Trace()
    gtu, gtv = gu - (gu * nrm) * nrm, gv - (gv * nrm) * nrm
    a_s = ng.BilinearForm(fes, symmetric=True, check_unused=False)
    a_s += (gtu * gtv) * ng.ds("inner")
    a_s += mass * u.Trace() * v.Trace() * ng.ds("inner")
    a_s.Assemble()
    r, c, val = a_s.mat.COO()
    return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(fes.ndof, fes.ndof))


def main():
    print("=" * 96)
    print("demo_kk -- DtN material-aware M  ->  radia.stream_function (ACA+)+TSVD + ridge design")
    print("=" * 96)
    with TaskManager():
        mesh, fes, make_A = hh.build_fem(BLOB, ORDER, MAXH)
        A_iron, A_free = make_A(MU_R), make_A(1.0)            # SAME mesh; iron vs free-space operator
        Ainv_iron = A_iron.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        Ainv_free = A_free.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        cdofs = hh.coil_dofs(mesh, fes); n_coil = len(cdofs)
        thetas = np.deg2rad(np.linspace(20, 160, 6)); phis = np.deg2rad(np.linspace(0, 270, 4))
        targets, _ = hh.targets_on_sphere(mesh, thetas, phis); n_t = len(targets)
        print("non-concentric iron blob mu_r=%.0f, order=%d; ndof=%d, coil psi DoFs=%d, targets=%d"
              % (MU_R, ORDER, fes.ndof, n_coil, n_t))
        M_iron = hh.build_transfer(mesh, fes, A_iron, Ainv_iron, cdofs, targets)   # material-aware
        M_free = hh.build_transfer(mesh, fes, A_free, Ainv_free, cdofs, targets)   # free-space kernel
        S_coil = np.array(surface_seminorm(fes, mesh)[np.ix_(cdofs, cdofs)].todense())
        # a producible target: a tilted-dipole coil trace pushed through the IRON operator
        gf = ng.GridFunction(fes); gf.vec[:] = 0.0
        gf.Set((ng.z + 0.6 * ng.x) / hh.a, ng.BND, definedon=mesh.Boundaries("inner"))
        psi_known = np.array([gf.vec[i] for i in cdofs])
        B_target = M_iron @ psi_known

        # ---- [1] ACA+ compresses the FORMED DtN M (entry = array lookup) ----
        res = aca_tsvd(n_t, n_coil, lambda i, j: float(M_iron[i, j]),
                       modes=min(n_t, n_coil), aca_eps=1e-8, method=3)
        rel_aca = np.linalg.norm(res.U @ np.diag(res.S) @ res.V.T - M_iron) / np.linalg.norm(M_iron)
        print("\n[1] aca_tsvd on the formed DtN M: k_aca=%d, ||U S V^T - M||/||M|| = %.2e"
              % (res.k_aca, rel_aca))

        # ---- [2] design through the existing solver: TSVD (min-L2) vs RIDGE (min current density) ----
        reg = RegularizedTSVD.from_stiffness(res, S_coil)
        psi_tsvd = pseudo_inverse_solve(res, B_target)
        psi_ridge = reg.solve(B_target, alpha=0.0)

        def jnorm(psi):
            return float(np.sqrt(psi @ (S_coil @ psi)))           # current-density seminorm ||K||

        def field_err(psi, Ainv, A):
            v = hh.fresh_forward(mesh, fes, A, Ainv, cdofs, psi, targets)
            return float(np.linalg.norm(v - B_target) / np.linalg.norm(B_target))

        e_tsvd, j_tsvd = field_err(psi_tsvd, Ainv_iron, A_iron), jnorm(psi_tsvd)
        e_ridge, j_ridge = field_err(psi_ridge, Ainv_iron, A_iron), jnorm(psi_ridge)
        print("\n[2] design (target hit in a fresh full Kelvin-FEM solve; ||K|| = current density):")
        print("    TSVD  (min-L2) : field_err=%.2e  ||K||=%.4e" % (e_tsvd, j_tsvd))
        print("    RIDGE (min-||K||): field_err=%.2e  ||K||=%.4e   (%.0f%% less current density)"
              % (e_ridge, j_ridge, 100.0 * (1.0 - j_ridge / j_tsvd)))

        # ---- [3] L-curve: one cached factorisation, alpha re-solves only the k x k core ----
        print("\n[3] ridge L-curve (field accuracy vs current density):")
        lcurve = []
        for alpha in ALPHAS:
            p = reg.solve(B_target, alpha=alpha)
            e, j = field_err(p, Ainv_iron, A_iron), jnorm(p)
            lcurve.append((alpha, e, j))
            print("    alpha=%.0e  field_err=%.2e  ||K||=%.4e" % (alpha, e, j))

        # ---- [4] free-space kernel through the SAME solver -> wrong in iron ----
        res_f = aca_tsvd(n_t, n_coil, lambda i, j: float(M_free[i, j]),
                         modes=min(n_t, n_coil), aca_eps=1e-8, method=3)
        psi_free = RegularizedTSVD.from_stiffness(res_f, S_coil).solve(B_target, alpha=0.0)
        e_free = field_err(psi_free, Ainv_iron, A_iron)
        print("\n[4] free-space-designed coil realised in the IRON system: field_err=%.2e"
              " (iron-aware ridge: %.2e)" % (e_free, e_ridge))

    ok = (rel_aca < 1e-8 and e_tsvd < 1e-6 and e_ridge < 1e-6
          and j_ridge < j_tsvd and e_free > 0.05)
    print("\n%s  ACA+ reproduces DtN M (%.0e); ridge hits target at LOWER ||K|| than TSVD "
          "(%.3e<%.3e); free-space design misses (%.1f%%)."
          % ("[PASS]" if ok else "[CHECK]", rel_aca, j_ridge, j_tsvd, 100.0 * e_free))
    print("VERDICT: the cheap DtN material-aware M plugs straight into radia.stream_function's "
          "(ACA+)+TSVD + ridge solver -- material-aware coil design with current-density minimisation, "
          "no DtN-specific design code.")
    return ok


if __name__ == "__main__":
    main()
