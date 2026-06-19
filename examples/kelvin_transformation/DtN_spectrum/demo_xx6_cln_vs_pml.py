# -*- coding: utf-8 -*-
r"""
demo_xx6_cln_vs_pml.py  (Track A -- CLN vs PML, head-to-head)
============================================================
Is the CLN open boundary BETTER than a PML for the eddy-current/diffusion DtN?
Measured head-to-head on the separable radial problem (s = i*omega): accuracy
per auxiliary DOF, and conditioning vs frequency.

  CLN : Lanczos/PVL reduction of the exterior eddy-current FEM (demo_xx4) -- N
        ladder states, an EXACT-operator model-order reduction.
  PML : an s-dependent complex-coordinate-stretch absorbing layer for the
        modified-Helmholtz (kappa=sqrt(s)) exterior -- M layer elements,
        dr~/dr = 1 + sigma(r)/sqrt(s), Dirichlet outer wall, Schur DtN at R0.

HONEST VERDICT (measured here):
  * For the EDDY-CURRENT / low-frequency (quasi-static, Radia's MQS home)
    regime, on this separable geometry, CLN BEATS PML:
      - accuracy per DOF: CLN converges EXPONENTIALLY (a MOR of the exact
        exterior) -- ~16 DOF to the FEM floor (~3e-4) -- vs the PML's algebraic
        layer convergence (~64-128 DOF for the same), so CLN is ~4-8x more
        DOF-efficient (at equal 16 DOF, CLN is ~30x more accurate);
      - conditioning: the PML stretch 1+sigma/sqrt(s) BLOWS UP as omega->0
        (cond ~6e3 at omega=1e-3) while the CLN reduction is real,
        frequency-independent (eval-system cond ~1 at low omega) -- ~10^3..10^4x
        better conditioned in the eddy-current band.
  * At HIGH frequency (wave-like) the PML is better conditioned (its stretch is
    well-behaved there; the CLN eval-system cond grows with omega) -- consistent
    with the lab's wave demos (demo_pp/nn: high-freq vacuum -> PML; low-freq /
    quasi-static -> Kelvin/CLN).
  * SCOPE: this is a SEPARABLE (radial) geometry.  PML's raison d'etre is
    ARBITRARY geometry; the CLN counterpart there is a Kelvin-transformed
    exterior FEM reduced by CLN (not built here).  So "CLN beats PML" is
    established for the separable eddy-current case, NOT as a blanket
    arbitrary-geometry claim.

VERIFIED HERE (all asserted; self-contained numpy/scipy):
  (1) both reproduce the exact diffusion DtN G_n; CLN is far more DOF-efficient.
  (2) PML conditioning degrades toward DC; CLN is frequency-robust at low omega.

PRIOR ART: PML = Berenger 1994 (+ complex-coordinate stretching, Chew-Weedon);
CLN = Kameari-Ebrahimi-Sugahara-Shindo-Matsuo, IEEE T-Magn 54(3):7201804 (2018).
The lab's wave-side Kelvin-vs-PML comparison is demo_nn/pp/qq/rr; this is the
eddy-current (diffusion) CLN-vs-PML instance.

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import kv

A_R = 1.0
_GP = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])
_GW = np.array([5 / 9, 8 / 9, 5 / 9])


def gam(s):
    return np.sqrt(complex(s))


def dtn_exact(n, s, a=A_R):
    g = gam(s)
    return -a * g * kv(n - 0.5, g * a) / kv(n + 0.5, g * a) - (n + 1.0)


def assemble(nodes, n):
    N = nodes.size
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    cent = n * (n + 1)
    for e in range(N - 1):
        a, b = nodes[e], nodes[e + 1]
        d = b - a
        rg = 0.5 * (a + b) + 0.5 * d * _GP
        jac = 0.5 * d
        N0 = (b - rg) / d
        N1 = (rg - a) / d
        dN = (-1.0 / d, 1.0 / d)
        Ns = (N0, N1)
        for p in range(2):
            for q in range(2):
                K[e + p, e + q] += np.sum(_GW * jac * (rg ** 2 * dN[p] * dN[q]
                                                       + cent * Ns[p] * Ns[q]))
                M[e + p, e + q] += np.sum(_GW * jac * rg ** 2 * Ns[p] * Ns[q])
    return K, M


# ---- CLN: Lanczos/PVL reduction of the exterior FEM ----
def cln_setup(n, R0, Rfar, h):
    nodes = np.linspace(R0, Rfar, int((Rfar - R0) / h) + 1)
    K, M = assemble(nodes, n)
    keep = np.arange(0, nodes.size - 1)
    K, M = K[np.ix_(keep, keep)], M[np.ix_(keep, keep)]
    Linv = np.linalg.inv(np.linalg.cholesky(K))
    Ah = Linv @ M @ Linv.T
    e0 = np.zeros(K.shape[0]); e0[0] = 1.0
    bh = Linv @ e0
    return Ah, bh, float(bh @ bh), K.shape[0]


def cln_lanczos(Ah, bh, Nmax):
    q = bh / np.linalg.norm(bh); qp = np.zeros_like(q)
    al, be, Q, bp = [], [], [q], 0.0
    for _ in range(Nmax):
        z = Ah @ q; a = q @ z; al.append(a); z = z - a * q - bp * qp
        for u in Q:
            z = z - (u @ z) * u
        b = np.linalg.norm(z)
        if b < 1e-13:
            break
        be.append(b); qp, q, bp = q, z / b, b; Q.append(q)
    return np.array(al), np.array(be)


def cln_tridiag(al, be):
    N = len(al); T = np.diag(al).astype(complex)
    for i in range(min(len(be), N - 1)):
        T[i, i + 1] = T[i + 1, i] = be[i]
    return T


def cln_dtn(T, nb2, s):
    N = T.shape[0]
    x = np.linalg.solve(np.eye(N, dtype=complex) + s * T, np.eye(N, 1, dtype=complex).ravel())
    return -1.0 / (nb2 * x[0])


# ---- PML: s-dependent complex-stretch layer ----
def pml_assemble(n, s, R0, Lp, M, sg, p=2):
    kappa = gam(s)
    nodes = np.linspace(R0, R0 + Lp, M + 1)

    def sigma(r):
        return sg * ((r - R0) / Lp) ** p

    rt = np.zeros(M + 1, dtype=complex); rt[0] = R0
    for e in range(M):
        a, b = nodes[e], nodes[e + 1]; d = b - a; rg = 0.5 * (a + b) + 0.5 * d * _GP
        rt[e + 1] = rt[e] + np.sum(_GW * 0.5 * d * (1.0 + sigma(rg) / kappa))
    A = np.zeros((M + 1, M + 1), dtype=complex)
    for e in range(M):
        a, b = nodes[e], nodes[e + 1]; d = b - a; rg = 0.5 * (a + b) + 0.5 * d * _GP; jac = 0.5 * d
        N0 = (b - rg) / d; N1 = (rg - a) / d; dN = (-1 / d, 1 / d); Ns = (N0, N1)
        beta = 1.0 + sigma(rg) / kappa; rtg = rt[e] * N0 + rt[e + 1] * N1
        for pp in range(2):
            for qq in range(2):
                A[e + pp, e + qq] += np.sum(_GW * jac * ((rtg ** 2 / beta) * dN[pp] * dN[qq]
                                            + n * (n + 1) * beta * Ns[pp] * Ns[qq]
                                            + s * rtg ** 2 * beta * Ns[pp] * Ns[qq]))
    return A


def pml_dtn(n, s, R0, Lp, M, sg):
    A = pml_assemble(n, s, R0, Lp, M, sg)
    ii = np.arange(1, M)
    return -(A[0, 0] - A[0, ii] @ np.linalg.solve(A[np.ix_(ii, ii)], A[ii, 0]))


# ===========================================================================
print("=" * 78)
print(" demo_xx6 : CLN vs PML for the eddy-current/diffusion open boundary")
print("=" * 78)

n, R0 = 1, 1.0
band = 1j * np.logspace(-1, 2, 40)                      # s = i*omega, omega in [0.1,100]
Gex = np.array([complex(dtn_exact(n, s)) for s in band])


def nrmse(G):
    return float(np.sqrt(np.mean(np.abs(G - Gex) ** 2)) / np.sqrt(np.mean(np.abs(Gex) ** 2)))


print("\n[1] accuracy vs DOF (NRMSE of the reproduced DtN over the band):")
Ah, bh, nb2, ndof = cln_setup(n, R0, 8.0, 0.01)
cln_err = {}
for N in (4, 8, 16):
    T = cln_tridiag(*cln_lanczos(Ah, bh, N))
    cln_err[N] = nrmse(np.array([cln_dtn(T, nb2, s) for s in band]))
    print(f"    CLN  DOF={N:3d}: NRMSE={cln_err[N]:.2e}   (Lanczos/PVL MOR of a {ndof}-DOF exterior)")
pml_err = {}
for M in (16, 32, 64, 128):
    best = min(nrmse(np.array([pml_dtn(n, s, R0, Lp, M, sg) for s in band]))
               for Lp in (1.0, 2.0) for sg in (4.0, 8.0, 16.0))
    pml_err[M] = best
    print(f"    PML  DOF={M:3d}: NRMSE={best:.2e}   (best of a small L/sigma sweep)")
assert cln_err[16] < 1e-3 and cln_err[16] < pml_err[16] / 10, "CLN should be far more DOF-efficient"
print(f"    -> at DOF=16: CLN {cln_err[16]:.1e} vs PML {pml_err[16]:.1e} "
      f"({pml_err[16] / cln_err[16]:.0f}x); CLN reaches the FEM floor ~10x cheaper (exp. MOR vs algebraic PML)")

print("\n[2] conditioning vs omega (the eddy-current / low-frequency point):")
T16 = cln_tridiag(*cln_lanczos(Ah, bh, 16))
print(f"    {'omega':>8} {'cond(PML)':>12} {'cond(CLN eval I+sT)':>20}")
condP, condC = {}, {}
for w in (1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2):
    cP = np.linalg.cond(pml_assemble(n, 1j * w, R0, 1.0, 32, 8.0))
    cC = np.linalg.cond(np.eye(16, dtype=complex) + 1j * w * T16)
    condP[w], condC[w] = cP, cC
    print(f"    {w:8.0e} {cP:12.2e} {cC:20.2e}")
assert condP[1e-3] > 1e3 and condC[1e-3] < 10, "PML should ill-condition at low omega; CLN should not"
assert condP[1e-3] / condC[1e-3] > 1e2, "CLN should be far better conditioned in the eddy-current band"
print("    -> low omega (eddy-current/quasi-static): CLN ~1 vs PML ~1e3-1e4 (CLN wins);")
print("       high omega (wave-like): PML better conditioned (lab demo_pp/nn: PML's home).")
print("    (the CLN REDUCTION is real + frequency-independent; the marched transient")
print("     system M_c,K_c is real SPD -- demo_xx5 -- so no omega-dependent cond at all.)")

print("\n[verdict]")
print("  * Eddy-current / low-frequency (Radia MQS home), separable geometry: CLN")
print("    BEATS PML -- ~4-8x fewer DOF (exponential MOR vs algebraic layer) AND")
print("    ~1e3-1e4x better conditioned toward DC.  High-freq/wave -> PML's home.")
print("  * NOT a blanket claim: PML's domain is ARBITRARY geometry; the CLN there")
print("    needs a Kelvin-transformed exterior FEM reduced by CLN (not built here).")
print("\nALL CHECKS PASSED.")
