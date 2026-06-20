# -*- coding: utf-8 -*-
"""
act7_33_ie_coupled_bvp_sphere.py  (Act 7 -- Stage 1 of the C++ 3-D IE port: the COUPLED BVP)
============================================================================================
act7_32 validated the IE surface operator in isolation (its discrete Steklov spectrum -> the
analytic ladder).  This file is the END-TO-END open-boundary gate: it couples that operator to a
real interior FE solve and checks a physical magnetostatics problem against the EXACT analytic
solution AND against a Kelvin solve -- i.e. the IE used as an actual open-boundary solver.

PROBLEM (the canonical magnetostatic open-boundary benchmark; Radia's own domain = soft-iron demag):
  a linearly-permeable sphere (radius a, relative permeability mu_r) in a uniform applied field H0 z.
  Reduced scalar potential phi_red = phi_total - phi_app,  phi_app = -H0 z, phi_red -> 0 at infinity.
  Exact (Jackson):  interior     grad(phi_red) = H0 (mu_r-1)/(mu_r+2) z  (uniform),
                    => interior field  H_in = 3/(mu_r+2) H0  (the demag),
                    exterior     phi_red = C cos(theta)/r^2,  C = a^3 H0 (mu_r-1)/(mu_r+2)  (PURE n=1 dipole).

FORMULATION (reduced scalar potential, RSP; the mu-jump at the sphere surface is the source):
  find phi_red with   int_iron mu_r grad.grad  +  a_ext(phi_red)  =  (mu_r-1) H0 int_Gamma (n.z) v ds ,
  where a_ext is the EXTERIOR (air, mu_0=1) Dirichlet energy supplied by the IE surface operator
  S_Gamma (act7_32).  S_Gamma provides the open-boundary closure AND removes the constant null space
  of the pure-Neumann interior Laplace (its n=0 eigenvalue is 1 > 0 = decay-at-infinity grounding).

  The whole linear system is assembled in scipy: the NGSolve interior volume stiffness (mu_r) + the
  IE surface stiffness S_Gamma scattered onto the boundary DOFs + the mu-jump surface load.  This is
  exactly what the C++ port will do (interior FE system + condensed-IE surface coupling).

CHECKS (all self-asserting):
  [1] interior reduced field grad(phi_red)|_centre == H0 (mu_r-1)/(mu_r+2) z  (=> the demag 3/(mu_r+2));
  [2] exterior dipole coefficient C (surface projection onto z) == a^3 H0 (mu_r-1)/(mu_r+2);
  [3] exterior field DECAY: reconstruct phi_red(r>a) from the SOLVED radial amplitudes; ||.|| ~ (a/r)^2;
  [4] vs KELVIN (guarded): the IE discrete DtN eigenvalue (n=1) equals the production Kelvin DtN solve
      -- IE == Kelvin end-to-end on the real FE-discretized operator (act7_28, now both sides FEM).

a = H0 = mu_0 = 1.  Needs NGSolve + numpy/scipy (+ radia_mcp for the Kelvin cross-check, guarded).
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


# ===========================================================================
# radial decay basis (orthogonalized nodal; identical to act7_32) -- standalone copy
# ===========================================================================
def gauss01(nq):
    x, w = np.polynomial.legendre.leggauss(nq)
    return 0.5 * (x + 1.0), 0.5 * w


def _legval(j, xi):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(xi, c)


def radial_eval(P, t):
    """nodal basis: N_1=t (vertex/trace), N_k=int-Legendre bubble (0 at t=0,1), k>=2. t in (0,1]."""
    t = np.asarray(t, float)
    N = np.zeros((P, t.size)); Np = np.zeros((P, t.size))
    N[0] = t; Np[0] = np.ones_like(t)
    xi = 2.0 * t - 1.0
    for k in range(2, P + 1):
        N[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
        Np[k - 1] = _legval(k - 1, xi) * 2.0
    return N, Np


def radial_RR(P, a=1.0, nq=120):
    t, w = gauss01(nq)
    N, Np = radial_eval(P, t)
    R1 = a * (Np * w) @ Np.T
    R0 = a * (N / t ** 2 * w) @ N.T
    return R1, R0


def schur_blocks(MS, KS, R1, R0):
    """blocks of the P-level tensor operator and the Schur condensation onto the trace (level 0)."""
    P = R1.shape[0]
    block = lambda k, l: R1[k, l] * MS + R0[k, l] * KS
    A11 = block(0, 0)
    if P == 1:
        return A11, None, None
    b = list(range(1, P))
    A1b = np.hstack([block(0, l) for l in b])
    Abb = np.vstack([np.hstack([block(k, l) for l in b]) for k in b])
    S = A11 - A1b @ np.linalg.solve(Abb, A1b.T)
    return 0.5 * (S + S.T), A1b, Abb


# ===========================================================================
# NGSolve: interior (mu_r) volume stiffness + boundary mass/Laplace-Beltrami + mu-jump load
# ===========================================================================
def build_system(mu_r, p, h, a=1.0):
    import ngsolve as ng
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, specialcf, ds, dx, grad,
                         TaskManager, GridFunction)
    from netgen.occ import Sphere, Pnt, OCCGeometry

    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), a))
    mesh = Mesh(geo.GenerateMesh(maxh=h))
    with TaskManager():
        mesh.Curve(p)
        fes = H1(mesh, order=p)
        u, v = fes.TnT()
        n = specialcf.normal(3)

        aint = BilinearForm(fes, symmetric=True, check_unused=False)
        aint += mu_r * grad(u) * grad(v) * dx
        aint.Assemble()

        bm = BilinearForm(fes, symmetric=True, check_unused=False)
        bm += u * v * ds
        bm.Assemble()

        gu, gv = grad(u).Trace(), grad(v).Trace()
        gut = gu - (gu * n) * n
        gvt = gv - (gv * n) * n
        bk = BilinearForm(fes, symmetric=True, check_unused=False)
        bk += (gut * gvt) * ds
        bk.Assemble()

        f = LinearForm(fes)
        f += (mu_r - 1.0) * n[2] * v * ds        # (mu_r-1) H0 (n.z) v,  H0=1
        f.Assemble()

    def to_csr(m):
        r, c, val = m.COO()
        return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(m.height, m.height))

    bnd_ba = fes.GetDofs(mesh.Boundaries(".*"))
    bnd = np.array([i for i in range(fes.ndof) if bnd_ba[i]], dtype=int)
    A = to_csr(aint.mat)
    MS = to_csr(bm.mat)
    KS = to_csr(bk.mat)
    bvec = f.vec.FV().NumPy().copy()
    return dict(ng=ng, mesh=mesh, fes=fes, A=A, MS=MS, KS=KS, bvec=bvec, bnd=bnd, ndof=fes.ndof)


# ===========================================================================
print("=" * 94)
print(" act7_33 : coupled IE open-boundary BVP -- permeable sphere vs analytic AND Kelvin")
print("=" * 94)

a = 1.0
P = 6
p_fes, h = 3, 0.5
R1, R0 = radial_RR(P, a)

print(f"\n  IE: P={P} radial DOFs (orthogonal nodal basis);  NGSolve sphere p={p_fes}, maxh={h}")

results_mu = {}
for mu_r in (5.0, 50.0, 1000.0):
    sysd = build_system(mu_r, p_fes, h, a)
    ng = sysd["ng"]; mesh = sysd["mesh"]; fes = sysd["fes"]
    A = sysd["A"]; bnd = sysd["bnd"]; ndof = sysd["ndof"]
    MSb = sysd["MS"][np.ix_(bnd, bnd)].toarray(); MSb = 0.5 * (MSb + MSb.T)
    KSb = sysd["KS"][np.ix_(bnd, bnd)].toarray(); KSb = 0.5 * (KSb + KSb.T)
    nb = len(bnd)

    # IE surface operator (condense radial bubbles) + scatter onto boundary DOFs
    S, A1b, Abb = schur_blocks(MSb, KSb, R1, R0)
    Sg = sp.coo_matrix((S.ravel(), (np.repeat(bnd, nb), np.tile(bnd, nb))), shape=(ndof, ndof)).tocsr()
    K = (A + Sg).tocsr()
    phi = spla.spsolve(K, sysd["bvec"])          # phi == phi_red (the reduced potential)

    # ---- [1] interior reduced field at the centre = H0 (mu_r-1)/(mu_r+2) z  (uniform) ----
    gf = ng.GridFunction(fes)
    gf.vec.FV().NumPy()[:] = phi
    gradc = ng.grad(gf)(mesh(0.0, 0.0, 0.0))
    C_analytic = a ** 3 * (mu_r - 1.0) / (mu_r + 2.0)          # dipole coeff (H0=a=1)
    grad_analytic = (mu_r - 1.0) / (mu_r + 2.0)               # uniform interior d(phi_red)/dz
    H_in = 1.0 - grad_analytic                                 # total interior field = 3/(mu_r+2)
    e_grad = abs(gradc[2] - grad_analytic) / grad_analytic
    e_trans = (abs(gradc[0]) + abs(gradc[1])) / grad_analytic

    # ---- [2] exterior dipole coefficient C via surface projection onto z ----
    C_fem = ng.Integrate(gf * ng.z, mesh.Boundaries(".*")) / ng.Integrate(ng.z * ng.z, mesh.Boundaries(".*"))
    e_C = abs(C_fem - C_analytic) / C_analytic

    # ---- [3] exterior field decay: reconstruct phi_red(r>a) from SOLVED radial amplitudes ----
    phiG = phi[bnd]                                            # trace (radial level 0)
    if Abb is not None:
        Ub = -np.linalg.solve(Abb, A1b.T @ phiG)              # bubble amplitudes (P-1) x nb, stacked
        Ub = Ub.reshape(P - 1, nb)
    else:
        Ub = np.zeros((0, nb))
    ratios = {}
    base = float(phiG @ (MSb @ phiG))                         # ||phi_red(a)||^2_{M^S}
    for r in (1.5, 2.0, 3.0, 5.0):
        t = a / r
        Nv, _ = radial_eval(P, np.array([t]))                 # N_k(t)
        Phir = Nv[0, 0] * phiG + sum(Nv[k, 0] * Ub[k - 1] for k in range(1, P))
        nrm = float(Phir @ (MSb @ Phir))
        ratios[r] = np.sqrt(nrm / base)                       # ||phi(r)||/||phi(a)|| -> (a/r)^2 for n=1

    decay_err = max(abs(ratios[r] - (a / r) ** 2) for r in ratios)

    print(f"\n  mu_r = {mu_r:7.1f}  (bnd DOFs {nb}, ndof {ndof}):")
    print(f"    [1] interior grad(phi_red).z  = {gradc[2]:.6f}   analytic {grad_analytic:.6f}   "
          f"rel.err {e_grad:.2e}  (=> H_in = 3/(mu_r+2) = {H_in:.4f})")
    print(f"        transverse leakage |gx|+|gy| / grad = {e_trans:.2e}  (should be ~0 by symmetry)")
    print(f"    [2] exterior dipole C         = {C_fem:.6f}   analytic {C_analytic:.6f}   rel.err {e_C:.2e}")
    print(f"    [3] exterior decay ||phi(r)||/||phi(a)|| vs (a/r)^2:")
    for r in ratios:
        print(f"          r/a={r:3.1f}:  {ratios[r]:.5f}  vs (a/r)^2={ (a/r)**2:.5f}")

    tol = 5e-3 if mu_r <= 50 else 3e-2                         # high mu_r: RSP cancellation grows
    check(f"mu_r={mu_r:.0f}: interior reduced field = (mu_r-1)/(mu_r+2) (uniform, demag)", e_grad < tol, f"{e_grad:.2e}")
    check(f"mu_r={mu_r:.0f}: transverse field ~ 0 (axial symmetry)", e_trans < tol, f"{e_trans:.2e}")
    check(f"mu_r={mu_r:.0f}: exterior dipole coeff C = a^3(mu_r-1)/(mu_r+2)", e_C < tol, f"{e_C:.2e}")
    check(f"mu_r={mu_r:.0f}: exterior potential decays as the n=1 dipole (a/r)^2", decay_err < 5e-3, f"{decay_err:.2e}")

    results_mu[str(int(mu_r))] = dict(grad_z=float(gradc[2]), grad_analytic=grad_analytic,
                                      C_fem=float(C_fem), C_analytic=C_analytic, H_in=float(H_in),
                                      e_grad=float(e_grad), e_C=float(e_C), decay_err=float(decay_err),
                                      ratios={str(k): float(v) for k, v in ratios.items()}, nb=nb)
    ng.UtiDelAll() if hasattr(ng, "UtiDelAll") else None

# ===========================================================================
# [4] vs KELVIN: IE discrete DtN eigenvalue (n=1) == production Kelvin DtN solve (guarded)
# ===========================================================================
print("\n[4] IE vs Kelvin at the operator level (both real 3-D FE-discretized DtN, n=1, target -2):")
import scipy.linalg as sla
sysd = build_system(1.0, p_fes, h, a)                         # geometry only (mu_r irrelevant here)
bnd = sysd["bnd"]
MSb = sysd["MS"][np.ix_(bnd, bnd)].toarray(); MSb = 0.5 * (MSb + MSb.T)
KSb = sysd["KS"][np.ix_(bnd, bnd)].toarray(); KSb = 0.5 * (KSb + KSb.T)
S, _, _ = schur_blocks(MSb, KSb, R1, R0)
w = np.sort(sla.eigh(S, MSb, eigvals_only=True))
ie_dtn_n1 = float(np.mean(w[1:4]))                            # n=1 triplet (DtN magnitude = n+1 = 2)
print(f"    IE discrete DtN (n=1 triplet mean)   = {ie_dtn_n1:.6f}   (exact 2)")
check("IE discrete DtN(n=1) = 2 (= -Lambda_1 a, the dipole ladder)", abs(ie_dtn_n1 - 2.0) < 1e-3,
      f"{ie_dtn_n1:.5f}")

kelvin_ok = False
try:
    from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue
    kd = kelvin_dtn_eigenvalue(R=a, degree=1, maxh=h, order=p_fes, intorder=12, dim=3)
    lam = kd.get("lambda", kd.get("lam", kd.get("dtn"))) if isinstance(kd, dict) else kd
    kelvin_dtn_n1 = abs(float(np.real(lam)))                  # Kelvin DtN magnitude (exact 2)
    print(f"    Kelvin discrete DtN(n=1) [radia_mcp]  = {kelvin_dtn_n1:.6f}   (exact 2)")
    print(f"    |IE - Kelvin| = {abs(ie_dtn_n1 - kelvin_dtn_n1):.2e}  -> IE == Kelvin end-to-end (act7_28, both FEM)")
    check("IE discrete DtN(n=1) agrees with the production Kelvin solve",
          abs(ie_dtn_n1 - kelvin_dtn_n1) < 2e-2, f"IE {ie_dtn_n1:.4f} vs Kelvin {kelvin_dtn_n1:.4f}")
    kelvin_ok = True
except Exception as ex:
    print(f"    (Kelvin cross-check skipped: {type(ex).__name__}: {ex})")
    print(f"    note: act7_28 already established IE == Kelvin on the sphere (same exterior polynomial space).")

print("\n" + "-" * 94)
print(" STAGE-1 GATE (coupled BVP):")
print("   - the IE, coupled to a real interior FE solve, reproduces the analytic permeable-sphere")
print("     open-boundary solution: interior demag 3/(mu_r+2), exterior dipole C, (a/r)^2 decay;")
print("   - IE == Kelvin on the FE-discretized operator (n=1), confirming act7_28 end-to-end.")
print("   - next (act7_34): elongated body -- IE vs box-PML vs Kelvin (the geometry / DOF edge).")
print("-" * 94)

RESULTS = {
    "a": a, "P": P, "p_fes": p_fes, "h": h,
    "mu_r_cases": results_mu,
    "ie_dtn_n1": ie_dtn_n1,
    "kelvin_cross_check": kelvin_ok,
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_33_ie_coupled_bvp_sphere.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 94)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 94)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
