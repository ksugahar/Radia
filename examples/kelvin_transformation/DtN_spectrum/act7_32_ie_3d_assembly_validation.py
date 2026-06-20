# -*- coding: utf-8 -*-
"""
act7_32_ie_3d_assembly_validation.py  (Act 7 -- Stage 1 of the C++ 3-D infinite-element port)
=============================================================================================
acts 7_25..7_31 verified the infinite element (IE) at the MODAL / separable level (pure numpy /
sympy on analytic spherical harmonics): the radial decay basis, the de Rham complex, the
well-conditioned orthogonal basis, the two vector Steklov ladders, the curved-element requirement.
This file is the FIRST genuinely-3-D-FEM step (the algorithm gate before the C++ port): it assembles
the IE as a real surface operator on a CURVED NGSolve sphere mesh and checks the discretized Steklov
(DtN) spectrum against the exact analytic ladder.

THE ASSEMBLY (the thing the C++ port must reproduce)
----------------------------------------------------
The exterior Dirichlet energy of a function expanded as (surface FE) x (radial decay basis) is a
TENSOR PRODUCT of a surface part and a radial part:

      a_ext(u,v) = sum_{k,l}  [ R1_kl * Mtil  +  R0_kl * Ktil ]   (block (k,l), one per radial pair)

      R1_kl = int_a^inf rho_k'(r) rho_l'(r) r^2 dr      (radial STIFFNESS)
      R0_kl = int_a^inf rho_k(r)  rho_l(r)      dr      (radial MASS)
      Mtil  = int_{unit sphere} phi_i phi_j  dOmega     (surface mass,  unit sphere)
      Ktil  = int_{unit sphere} grad_S phi_i . grad_S phi_j dOmega  (Laplace-Beltrami)

For a single spherical harmonic mode n, Mtil -> 1 and Ktil -> n(n+1), so the block operator reduces
to the modal radial energy  E_kl(n) = R1_kl + n(n+1) R0_kl  of act7_25/28 -- i.e. this 3-D assembly
is EXACTLY the verified modal IE, mode by mode.  (Self-checked below.)

The radial DOFs are then statically condensed onto the surface trace level, giving the discrete DtN
surface stiffness S_Gamma; eig(S_Gamma, Mtil) is the discrete Steklov spectrum, which must approximate
the analytic ladder  +(n+1)/a  with multiplicity (2n+1):  1, 2,2,2, 3,3,3,3,3, 4(x7), ...

RADIAL BASIS -- orthogonalized, NEVER naive monomials (the act7_28 conditioning lesson)
--------------------------------------------------------------------------------------
t = a/r in (0,1] (t=1 at the surface r=a, t=0 at infinity).  N_1(t) = t  (the vertex / trace-carrying
function, value 1 at the surface, decaying); N_k(t) for k>=2 are integrated-Legendre bubbles (value 0
at BOTH t=0 and t=1).  Span{N_1..N_P} = {t, t^2, ..., t^P} = {a/r, ..., (a/r)^P} -- the SAME exterior
polynomial space as the monomials, but well-conditioned AND with a clean trace (only N_1 nonzero at
the surface, so exactly one radial DOF per surface node couples to the interior FE).

a = 1 (unit sphere) throughout, so Mtil = M^S and Ktil = K^S directly (no a-scaling); the 1/a scaling
of the ladder is recorded separately.

Needs NGSolve (surface mass + Laplace-Beltrami) + numpy/scipy.  Self-asserting; writes JSON.
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


# ===========================================================================
# 1. radial decay basis  (orthogonalized nodal: vertex t + integrated-Legendre bubbles)
# ===========================================================================
def gauss01(nq):
    x, w = np.polynomial.legendre.leggauss(nq)
    return 0.5 * (x + 1.0), 0.5 * w


def _legval(j, xi):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(xi, c)


def radial_eval(P, t, kind):
    """return N (P x len(t)) and Np (P x len(t)) for the chosen radial basis on t in (0,1]."""
    t = np.asarray(t, float)
    N = np.zeros((P, t.size))
    Np = np.zeros((P, t.size))
    if kind == "mono":
        for k in range(1, P + 1):
            N[k - 1] = t ** k
            Np[k - 1] = k * t ** (k - 1)
    elif kind == "nodal":
        # N_1 = t  (vertex / trace)
        N[0] = t
        Np[0] = np.ones_like(t)
        # N_k = integrated Legendre bubble phi_k(xi)=int_{-1}^{xi} L_{k-1}, xi=2t-1, k>=2
        xi = 2.0 * t - 1.0
        for k in range(2, P + 1):
            # phi_k(xi) = (L_k(xi) - L_{k-2}(xi)) / (2k-1)   (exact antiderivative identity)
            N[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
            Np[k - 1] = _legval(k - 1, xi) * 2.0  # d/dt int_{-1}^{2t-1} L_{k-1} = L_{k-1}(2t-1)*2
    else:
        raise ValueError(kind)
    return N, Np


def radial_RR(P, kind, a=1.0, nq=120):
    """R1_kl = a * int_0^1 N_k' N_l' dt ;  R0_kl = a * int_0^1 N_k N_l / t^2 dt."""
    t, w = gauss01(nq)
    N, Np = radial_eval(P, t, kind)
    R1 = a * (Np * w) @ Np.T
    R0 = a * (N / t ** 2 * w) @ N.T
    return R1, R0


def trace_vec(P, kind):
    """g_k = N_k(t=1) = rho_k(r=a)."""
    N, _ = radial_eval(P, np.array([1.0]), kind)
    return N[:, 0].copy()


def modal_dtn(R1, R0, g, kappa):
    """per-mode DtN(kappa) = -1 / (g^T E^{-1} g),  E = R1 + kappa R0."""
    E = R1 + kappa * R0
    return -1.0 / (g @ np.linalg.solve(E, g))


# ===========================================================================
# 2. NGSolve surface matrices on a CURVED unit sphere (M^S, K^S on boundary DOFs)
# ===========================================================================
def surface_matrices(p, h, a=1.0):
    from ngsolve import (Mesh, H1, BilinearForm, specialcf, ds, grad, TaskManager)
    from netgen.occ import Sphere, Pnt, OCCGeometry

    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), a))
    mesh = Mesh(geo.GenerateMesh(maxh=h))
    with TaskManager():
        mesh.Curve(p)
        fes = H1(mesh, order=p)
        u, v = fes.TnT()
        n = specialcf.normal(3)
        bm = BilinearForm(fes, symmetric=True, check_unused=False)
        bm += u * v * ds
        bm.Assemble()
        gu, gv = grad(u).Trace(), grad(v).Trace()
        gut = gu - (gu * n) * n
        gvt = gv - (gv * n) * n
        bk = BilinearForm(fes, symmetric=True, check_unused=False)
        bk += (gut * gvt) * ds
        bk.Assemble()

    def to_csr(m):
        r, c, val = m.COO()
        return sp.csr_matrix((np.array(val), (np.array(r), np.array(c))), shape=(m.height, m.height))

    bnd_ba = fes.GetDofs(mesh.Boundaries(".*"))
    bnd = np.array([i for i in range(fes.ndof) if bnd_ba[i]], dtype=int)
    MS = to_csr(bm.mat)[np.ix_(bnd, bnd)].toarray()
    KS = to_csr(bk.mat)[np.ix_(bnd, bnd)].toarray()
    MS = 0.5 * (MS + MS.T)
    KS = 0.5 * (KS + KS.T)
    return MS, KS, len(bnd)


# ===========================================================================
# 3. assemble the IE surface operator and condense radial DOFs  ->  S_Gamma
# ===========================================================================
def ie_surface_operator(MS, KS, R1, R0, g):
    """Build the P-level tensor operator, reorder so the TRACE level is first (g.e_k != 0),
    statically condense the bubble levels -> S_Gamma (acts on the surface trace DOFs)."""
    P = R1.shape[0]
    Ns = MS.shape[0]
    # the trace functional is g; for the nodal basis g = e_1 -> level 0 is the trace, 1..P-1 bubbles.
    # (general g: rotate so the trace is a single level -- here we rely on nodal g=e_1.)
    assert abs(g[0] - 1.0) < 1e-12 and np.allclose(g[1:], 0.0), "nodal basis expected (g=e_1)"

    def block(k, l):
        return R1[k, l] * MS + R0[k, l] * KS

    A11 = block(0, 0)
    if P == 1:
        return A11
    b = list(range(1, P))
    A1b = np.hstack([block(0, l) for l in b])               # Ns x (P-1)Ns
    Abb = np.vstack([np.hstack([block(k, l) for l in b]) for k in b])  # (P-1)Ns square
    S = A11 - A1b @ np.linalg.solve(Abb, A1b.T)
    return 0.5 * (S + S.T)


# ===========================================================================
print("=" * 94)
print(" act7_32 : 3-D IE assembly + discrete Steklov spectrum vs the analytic ladder (n+1)/a")
print("=" * 94)

P = 6
a = 1.0

# ---- [0] radial self-checks (quadrature + span + conditioning + per-mode exactness) ----
print("\n[0] radial basis self-checks (a=1):")
R1m, R0m = radial_RR(P, "mono", a)
R1n, R0n = radial_RR(P, "nodal", a)
gm, gn = trace_vec(P, "mono"), trace_vec(P, "nodal")

# quadrature vs closed forms for the monomial basis: R1=kl/((k+l)-1), R0=1/((k+l)-1)
k = np.arange(1, P + 1)
R1m_exact = np.outer(k, k) / (k[:, None] + k[None, :] - 1.0)
R0m_exact = 1.0 / (k[:, None] + k[None, :] - 1.0)
check("monomial R1 quadrature == kl/((k+l)-1)", np.max(np.abs(R1m - R1m_exact)) < 1e-9,
      f"err {np.max(np.abs(R1m - R1m_exact)):.1e}")
check("monomial R0 quadrature == 1/((k+l)-1)", np.max(np.abs(R0m - R0m_exact)) < 1e-9,
      f"err {np.max(np.abs(R0m - R0m_exact)):.1e}")
check("monomial trace g = (1,1,...,1)", np.allclose(gm, 1.0))
check("nodal trace g = (1,0,...,0) (clean trace: only N_1 at the surface)",
      abs(gn[0] - 1) < 1e-12 and np.allclose(gn[1:], 0.0))

# same space => identical per-mode DtN from monomial and nodal bases
print("\n    per-mode DtN(kappa=n(n+1)): monomial vs nodal (same space) and vs exact -(n+1):")
for n in range(0, 5):
    kap = n * (n + 1)
    dm = modal_dtn(R1m, R0m, gm, kap)
    dn = modal_dtn(R1n, R0n, gn, kap)
    print(f"      n={n}: mono {dm:8.4f}  nodal {dn:8.4f}  exact {-(n+1):4d}")
    check(f"n={n}: nodal == monomial DtN (same exterior space)", abs(dm - dn) < 1e-7)
    if n <= P - 1:
        check(f"n={n}: DtN exact = -(n+1) for P>=n+1", abs(dn - (-(n + 1))) < 1e-7, f"{dn:.5f}")

# conditioning: nodal radial energy well-conditioned vs monomial Hilbert-ill
for n in (1, 4):
    cm = np.linalg.cond(R1m + n * (n + 1) * R0m)
    cn = np.linalg.cond(R1n + n * (n + 1) * R0n)
    print(f"    cond E(n={n}, P={P}):  monomial {cm:.2e}   nodal {cn:.2e}")
    check(f"n={n}: nodal radial energy far better conditioned than monomial",
          cn < cm / 50.0, f"mono {cm:.1e} / nodal {cn:.1e}")

# ---- [1] the 3-D FEM assembly: discrete Steklov spectrum on a curved sphere ----
print("\n[1] discrete Steklov spectrum eig(S_Gamma, M^S) on a CURVED unit sphere (target = n+1):")
ladder = []                       # [1, 2,2,2, 3,3,3,3,3, ...]
for n in range(0, 4):
    ladder += [n + 1] * (2 * n + 1)
ladder = np.array(ladder, float)

p_fes, h = 3, 0.35
MS, KS, nb = surface_matrices(p_fes, h, a)
print(f"    NGSolve curved sphere: p={p_fes}, maxh={h}, boundary DOFs = {nb}")
S = ie_surface_operator(MS, KS, R1n, R0n, gn)
w = np.sort(sla.eigh(S, MS, eigvals_only=True))
nshow = len(ladder)
print("    discrete:  " + "  ".join(f"{x:6.3f}" for x in w[:nshow]))
print("    analytic:  " + "  ".join(f"{x:6.3f}" for x in ladder))
relerr = np.abs(w[:nshow] - ladder) / ladder
print("    rel.err :  " + "  ".join(f"{x:6.0e}" for x in relerr))

check("n=0 monopole eigenvalue = 1", abs(w[0] - 1.0) < 5e-3, f"{w[0]:.4f}")
check("n=1 triplet ~ 2 (mult 3)", np.max(np.abs(w[1:4] - 2.0)) < 3e-2, f"{w[1:4]}")
check("n=2 quintet ~ 3 (mult 5)", np.max(np.abs(w[4:9] - 3.0)) < 6e-2, f"max dev {np.max(np.abs(w[4:9]-3.0)):.1e}")
check("n=3 septet ~ 4 (mult 7)", np.max(np.abs(w[9:16] - 4.0)) < 2e-1, f"max dev {np.max(np.abs(w[9:16]-4.0)):.1e}")
check("spectrum is SPD (all eigenvalues > 0)", w[0] > 0)

# ---- [2] h-convergence of the discrete Steklov spectrum (angular FEM error -> 0) ----
print("\n[2] h-convergence of the n=1 triplet (lambda=2) and n=2 quintet (lambda=3), p=3:")
hs = [0.6, 0.45, 0.32]
err1, err2 = [], []
for hh in hs:
    MSx, KSx, nbx = surface_matrices(p_fes, hh, a)
    Sx = ie_surface_operator(MSx, KSx, R1n, R0n, gn)
    wx = np.sort(sla.eigh(Sx, MSx, eigvals_only=True))
    e1 = float(np.max(np.abs(wx[1:4] - 2.0)))
    e2 = float(np.max(np.abs(wx[4:9] - 3.0)))
    err1.append(e1); err2.append(e2)
    print(f"    maxh={hh:4.2f} (bnd {nbx:4d}):  |lam1-2|max={e1:.2e}   |lam2-3|max={e2:.2e}")
check("n=1 Steklov eigenvalue converges as the surface mesh refines", err1[-1] < err1[0],
      f"{err1[0]:.1e} -> {err1[-1]:.1e}")
check("n=2 Steklov eigenvalue converges as the surface mesh refines", err2[-1] < err2[0],
      f"{err2[0]:.1e} -> {err2[-1]:.1e}")

# ---- [3] radial order P: need P>=n+1 to capture mode n (radial spectral exactness in the 3-D op) ----
print("\n[3] radial order P controls mode capture (fixed fine surface): n=2 needs P>=3:")
MSf, KSf, _ = surface_matrices(p_fes, 0.32, a)
for Pr in (2, 3, 4):
    R1p, R0p = radial_RR(Pr, "nodal", a)
    gp = trace_vec(Pr, "nodal")
    Sp = ie_surface_operator(MSf, KSf, R1p, R0p, gp)
    wp = np.sort(sla.eigh(Sp, MSf, eigvals_only=True))
    e2 = float(np.max(np.abs(wp[4:9] - 3.0)))
    print(f"    P={Pr}: n=2 quintet |lam-3|max = {e2:.2e}")
    if Pr == 2:
        e2_P2 = e2
    if Pr == 3:
        e2_P3 = e2
check("radial P=3 captures n=2 far better than P=2 (P>=n+1 rule, in the full 3-D operator)",
      e2_P3 < e2_P2 / 5.0, f"P2 {e2_P2:.1e} -> P3 {e2_P3:.1e}")

print("\n" + "-" * 94)
print(" STAGE-1 GATE (assembly):")
print("   - the 3-D IE surface operator (M^S (x) R1 + K^S (x) R0, radial DOFs condensed) reproduces")
print("     the analytic Steklov ladder (n+1)/a with the correct (2n+1) multiplicities;")
print("   - it converges under surface refinement (angular FEM error) and needs radial P>=n+1 to")
print("     capture mode n -- exactly the modal IE (act7_25/28), now on a real curved 3-D FE surface.")
print("   - next (act7_33): couple to an interior FE solve (sphere BVP) vs analytic AND a Kelvin solve;")
print("     (act7_34): elongated body, IE vs box-PML vs Kelvin -- the geometry edge at production fidelity.")
print("-" * 94)

RESULTS = {
    "P": P, "a": a, "p_fes": p_fes, "h_main": h, "bnd_dofs": nb,
    "ladder": ladder.tolist(),
    "steklov_discrete": w[:nshow].tolist(),
    "steklov_relerr": relerr.tolist(),
    "h_convergence": {"hs": hs, "err_n1": err1, "err_n2": err2},
    "radial_P_capture_n2": {"P2": e2_P2, "P3": e2_P3},
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_32_ie_3d_assembly_validation.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 94)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 94)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
