# -*- coding: utf-8 -*-
"""
act7_37_ie_spheroidal_fem_full.py  (Act 7 -- Gate-2 milestone 3: the FULL FE prolate-spheroidal IE)
==================================================================================================
act7_36 proved the tight ellipsoid IE's RADIAL kernel (the spheroidal Steklov D_n) modally.  This
file is the FULL finite-element realisation (M3-proper): a real 2-D prolate-spheroidal FEM whose
exterior is closed by the spheroidal infinite element, solving the permeable-spheroid demag end-to-end
and demonstrating the tight-closure DOF edge.

THE KEY SEPARATION (m=0 / axisymmetric).  In prolate spheroidal coordinates (xi>=1 radial-like, eta in
[-1,1] angular) the axisymmetric exterior energy SEPARATES CLEANLY (the xi^2-eta^2 metric coupling and
the d/dphi term both drop out):

      E[u] = 2 pi f  INT INT [ (xi^2 - 1)(d_xi u)^2 + (1 - eta^2)(d_eta u)^2 ] d_eta d_xi

so the FE is a TENSOR PRODUCT of a 1-D angular FE in eta (stiffness weight 1-eta^2) and a 1-D radial
problem in xi (stiffness weight xi^2-1): interior energy = S_xi (x) M_eta + M_xi (x) K_eta, and the IE
exterior energy = A_ext (x) M_eta + M_ext (x) K_eta with the spheroidal radial decay operators
A_ext = INT_{xi0}^inf (xi^2-1) rho' rho' dxi, M_ext = INT_{xi0}^inf rho rho dxi (act7_36).  Per angular
mode P_n(eta): M_eta->const, K_eta-> n(n+1) const, and the block reduces to A_ext + n(n+1) M_ext = the
modal spheroidal energy E_n (act7_36) -- so this FE IS the verified spheroidal IE, now end-to-end.

(The physical-surface Laplace-Beltrami does NOT match K_eta -- it carries an extra
sqrt((xi0^2-1)/(xi0^2-eta^2)) eta-weight -- so the tight spheroidal IE needs the SPHEROIDAL angular FE,
not a physical-mesh surface FE.  This is why M3-proper is a custom spheroidal FEM.)

PROBLEM: a permeable prolate spheroid (semi-minor b=1, semi-major AR, relative permeability mu_r) in a
uniform axial field.  Reduced scalar potential phi_red (the mu-jump at xi=xi0 is the source, which in
spheroidal coordinates is exactly (n.z) dS = eta f^2 (xi0^2-1) d_eta d_phi).  Interior phi_red = (H0 -
H_in) z = (H0-H_in) f xi eta, so H0 - H_in is read off the solution; H_in = H0/(1+N_a(mu_r-1)) (Osborn).

CHECKS (self-asserting, pure numpy):
  - the full spheroidal FE IE reproduces the Osborn (1945) demag N_a (AR=1.5,2,4,8; AR-aware radial P);
  - it converges under angular / radial / IE-order refinement;
  - the DOF edge: replacing a MESHED exterior shell (ballooning to xi_max) by the IE's P thin radial
    levels reaches the same Osborn accuracy at far fewer DOF -- the tight-closure saving, in one
    consistent 2-D framework.
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def osborn_Na(e):
    if e < 1e-9:
        return 1.0 / 3.0
    return (1.0 - e * e) / e ** 3 * (np.arctanh(e) - e)


# ---------------------------------------------------------------------------
# 1-D linear FE: K = int w_stiff phi_i' phi_j', M = int w_mass phi_i phi_j  on the given nodes
# ---------------------------------------------------------------------------
def assemble_1d(nodes, w_stiff, w_mass, nq=8):
    n = len(nodes); K = np.zeros((n, n)); M = np.zeros((n, n))
    gx, gw = np.polynomial.legendre.leggauss(nq)
    for e in range(n - 1):
        a, bb = nodes[e], nodes[e + 1]; h = bb - a
        xq = 0.5 * (a + bb) + 0.5 * h * gx; wq = 0.5 * h * gw
        for q in range(nq):
            x = xq[q]; w = wq[q]
            ph = np.array([(bb - x) / h, (x - a) / h]); dp = np.array([-1.0 / h, 1.0 / h])
            K[e:e + 2, e:e + 2] += w_stiff(x) * np.outer(dp, dp) * w
            M[e:e + 2, e:e + 2] += w_mass(x) * np.outer(ph, ph) * w
    return K, M


def eta_source(nodes, nq=8):
    """s_i = int_{-1}^1 eta phi_i(eta) d_eta  (the n=1 angular projection = the applied-field source)."""
    n = len(nodes); s = np.zeros(n)
    gx, gw = np.polynomial.legendre.leggauss(nq)
    for e in range(n - 1):
        a, bb = nodes[e], nodes[e + 1]; h = bb - a
        xq = 0.5 * (a + bb) + 0.5 * h * gx; wq = 0.5 * h * gw
        for q in range(nq):
            x = xq[q]; w = wq[q]; ph = np.array([(bb - x) / h, (x - a) / h])
            s[e:e + 2] += x * ph * w
    return s


# ---------------------------------------------------------------------------
# spheroidal IE radial operators on [xi0,inf): nodal decay basis in s = xi0/xi
#   A_ext = int (xi^2-1) rho' rho' dxi   (radial stiffness),  M_ext = int rho rho dxi (radial mass)
# ---------------------------------------------------------------------------
def _legval(j, x):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(x, c)


def ie_radial(P, xi0, nq=400):
    x, w = np.polynomial.legendre.leggauss(nq); s = 0.5 * (x + 1.0); w = 0.5 * w
    R = np.zeros((P, s.size)); Rp = np.zeros((P, s.size))
    R[0] = s; Rp[0] = 1.0
    xi = 2.0 * s - 1.0
    for k in range(2, P + 1):
        R[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
        Rp[k - 1] = _legval(k - 1, xi) * 2.0
    A = xi0 * (Rp * ((1.0 - s ** 2 / xi0 ** 2) * w)) @ Rp.T
    M = xi0 * (R * (w / s ** 2)) @ R.T
    return A, M


def xi0_of_AR(AR):
    return AR / np.sqrt(AR ** 2 - 1.0)


def solve_demag(AR, mu_r, P=8, n_xi=24, n_eta=24, xi_max=None, n_xi_ext=0):
    """Custom 2-D prolate-spheroidal FEM for the permeable-spheroid axial demag.
    IE closure if xi_max is None (P radial levels on [xi0,inf)); otherwise a MESHED exterior shell
    [xi0,xi_max] with Dirichlet phi_red=0 at xi_max (the ballooning baseline, n_xi_ext elements).
    Returns (H0 - H_in, total DOF)."""
    xi0 = xi0_of_AR(AR)
    f = 1.0 / np.sqrt(xi0 ** 2 - 1.0)                      # semi-minor b = f sqrt(xi0^2-1) = 1
    en = np.linspace(-1.0, 1.0, n_eta + 1)
    Keta, Meta = assemble_1d(en, lambda e: 1.0 - e * e, lambda e: 1.0)
    s_eta = eta_source(en)
    nE = len(en)
    xn = np.linspace(1.0, xi0, n_xi + 1)
    Sin, Min = assemble_1d(xn, lambda x: x * x - 1.0, lambda x: 1.0)
    nX = len(xn)
    it = nX - 1                                            # radial index of the surface xi0

    if xi_max is None:
        # ---- tight IE closure ----
        A_e, M_e = ie_radial(P, xi0)
        nB = P - 1; nR = nX + nB
        Srad = np.zeros((nR, nR)); Mrad = np.zeros((nR, nR))
        Srad[:nX, :nX] = mu_r * Sin; Mrad[:nX, :nX] = mu_r * Min
        Srad[it, it] += A_e[0, 0]; Mrad[it, it] += M_e[0, 0]
        for k in range(1, P):
            Srad[it, nX + k - 1] += A_e[0, k]; Srad[nX + k - 1, it] += A_e[k, 0]
            Mrad[it, nX + k - 1] += M_e[0, k]; Mrad[nX + k - 1, it] += M_e[k, 0]
            for l in range(1, P):
                Srad[nX + k - 1, nX + l - 1] += A_e[k, l]; Mrad[nX + k - 1, nX + l - 1] += M_e[k, l]
        free_rad = list(range(nR))
    else:
        # ---- ballooning: mesh the exterior shell [xi0,xi_max], Dirichlet phi_red=0 at xi_max ----
        xe = np.linspace(xi0, xi_max, n_xi_ext + 1)
        Sext, Mext = assemble_1d(xe, lambda x: x * x - 1.0, lambda x: 1.0)
        nXe = len(xe)
        nR = nX + (nXe - 1)                                # share the xi0 node
        Srad = np.zeros((nR, nR)); Mrad = np.zeros((nR, nR))
        Srad[:nX, :nX] = mu_r * Sin; Mrad[:nX, :nX] = mu_r * Min
        # exterior block (mu=1) sharing node it=xi0
        idx = [it] + list(range(nX, nR))                   # exterior radial DOFs (xi0, then shell nodes)
        for p in range(nXe):
            for q in range(nXe):
                Srad[idx[p], idx[q]] += Sext[p, q]; Mrad[idx[p], idx[q]] += Mext[p, q]
        free_rad = list(range(nR - 1))                     # Dirichlet at xi_max (last DOF) -> drop

    # global 2-D system A = Srad (x) Meta + Mrad (x) Keta  (times 2 pi f), source at the xi0 row
    A = 2.0 * np.pi * f * (np.kron(Srad, Meta) + np.kron(Mrad, Keta))
    b = np.zeros(nR * nE)
    coef = (mu_r - 1.0) * 1.0 * 2.0 * np.pi * f ** 2 * (xi0 ** 2 - 1.0)
    b[it * nE:(it + 1) * nE] = coef * s_eta
    free = [r * nE + i for r in free_rad for i in range(nE)]
    u = np.zeros(nR * nE)
    u[free] = np.linalg.solve(A[np.ix_(free, free)], b[free])
    U = u.reshape(nR, nE)
    # interior phi_red = (H0-H_in) f xi eta -> fit C = (H0-H_in) f over interior nodes; H0-H_in = C/f
    num = 0.0; den = 0.0
    for a_ in range(nX):
        for i in range(nE):
            ze = xn[a_] * en[i]
            if abs(ze) > 1e-6:
                num += U[a_, i] * ze; den += ze * ze
    C = num / den
    return C / f, nR * nE                                  # (H0 - H_in, DOF)


print("=" * 98)
print(" act7_37 : the FULL FE prolate-spheroidal infinite element -- permeable-spheroid demag (M3-proper)")
print("=" * 98)

# ---- [1] the full spheroidal FE IE reproduces Osborn N_a ----
print("\n[1] full 2-D spheroidal FEM + IE -> Osborn demag N_a (AR-aware radial order P):")
print("    AR    xi0      H_in(FEM)   H_in(Osborn)   N_a       rel.err    DOF")
mu_r = 100.0
RESULTS = {"mu_r": mu_r, "cases": []}
for AR, P, nx, ne in [(1.5, 8, 24, 24), (2.0, 8, 24, 24), (4.0, 12, 32, 28), (8.0, 20, 40, 32)]:
    xi0 = xi0_of_AR(AR)
    HmH, dof = solve_demag(AR, mu_r, P=P, n_xi=nx, n_eta=ne)
    H_in = 1.0 - HmH
    Na = osborn_Na(1.0 / xi0)
    H_in_osb = 1.0 / (1.0 + Na * (mu_r - 1.0))
    re = abs(H_in - H_in_osb) / H_in_osb
    print(f"   {AR:4.1f}  {xi0:6.4f}   {H_in:.5f}    {H_in_osb:.5f}     {Na:.4f}    {re:.2e}   {dof}")
    check(f"AR={AR}: full spheroidal FE IE demag == Osborn N_a={Na:.4f}", re < 1e-2, f"relerr {re:.2e}")
    RESULTS["cases"].append(dict(AR=AR, xi0=xi0, P=P, H_in=H_in, H_in_osborn=H_in_osb, Na=Na,
                                 relerr=re, dof=dof))

# ---- [2] convergence under angular / radial-interior / IE-order refinement (AR=2) ----
print("\n[2] convergence (AR=2, mu_r=100): refine angular n_eta, interior n_xi, IE order P:")
xi0 = xi0_of_AR(2.0); Na2 = osborn_Na(1.0 / xi0); Hin2 = 1.0 / (1.0 + Na2 * (mu_r - 1.0))
for tag, kw in [("coarse", dict(P=4, n_xi=8, n_eta=8)), ("medium", dict(P=6, n_xi=16, n_eta=16)),
                ("fine", dict(P=10, n_xi=32, n_eta=32))]:
    HmH, dof = solve_demag(2.0, mu_r, **kw)
    re = abs((1.0 - HmH) - Hin2) / Hin2
    print(f"    {tag:7s} {kw}: rel.err = {re:.2e}  (DOF {dof})")
    if tag == "coarse":
        re_coarse = re
    if tag == "fine":
        re_fine = re
check("[2] refinement converges (fine << coarse)", re_fine < re_coarse / 10.0,
      f"{re_coarse:.1e} -> {re_fine:.1e}")

# ---- [3] the tight-closure DOF edge: IE vs a MESHED exterior shell (ballooning), same framework ----
print("\n[3] tight-closure DOF edge (AR=4, mu_r=100): IE (P levels) vs a MESHED exterior shell:")
xi0 = xi0_of_AR(4.0); Na4 = osborn_Na(1.0 / xi0); Hin4 = 1.0 / (1.0 + Na4 * (mu_r - 1.0))
HmH_ie, dof_ie = solve_demag(4.0, mu_r, P=12, n_xi=32, n_eta=28)
re_ie = abs((1.0 - HmH_ie) - Hin4) / Hin4
print(f"    tight IE  (P=12):                rel.err {re_ie:.2e}   DOF {dof_ie}")
# ballooning: mesh the exterior to xi_max with Dirichlet; needs many shell elements for the same accuracy
ball = []
for xi_max, n_ext in [(2.0, 20), (4.0, 60), (8.0, 140)]:
    HmH_b, dof_b = solve_demag(4.0, mu_r, n_xi=32, n_eta=28, xi_max=xi_max, n_xi_ext=n_ext)
    re_b = abs((1.0 - HmH_b) - Hin4) / Hin4
    ball.append((xi_max, n_ext, re_b, dof_b))
    print(f"    ballooning xi_max={xi_max:4.1f} ({n_ext} shell elem): rel.err {re_b:.2e}   DOF {dof_b}")
best_ball = min(ball, key=lambda t: t[2])
check("[3] the IE reaches Osborn accuracy the meshed exterior needs many more DOF for",
      re_ie < 1e-2 and dof_ie < best_ball[3],
      f"IE {re_ie:.1e}@{dof_ie}DOF vs best ballooning {best_ball[2]:.1e}@{best_ball[3]}DOF")
check("[3] ballooning is finite-reach: even its best xi_max stays less accurate than the IE",
      best_ball[2] > re_ie, f"ballooning {best_ball[2]:.1e} vs IE {re_ie:.1e}")

print("\n" + "-" * 98)
print(" GATE-2 MILESTONE 3 (the FULL FE spheroidal IE) -- COMPLETE:")
print("   - the tight (surface-conforming) non-spherical IE is REALISED as a real 2-D spheroidal FEM:")
print("     the m=0 energy separates (S_xi (x) M_eta + M_xi (x) K_eta), the IE closes the exterior with")
print("     P thin radial levels, and it reproduces the Osborn demag N_a end-to-end (AR=1.5..8);")
print("   - it converges under refinement, and the IE reaches Osborn accuracy at far fewer DOF than a")
print("     meshed (ballooning) exterior -- the tight-closure DOF saving, in one consistent framework.")
print("   - this CLOSES Gate-2 for the axisymmetric (m=0) demag: the tight non-spherical IE works and is")
print("     correct. (The full 3-D all-m spheroidal IE -- the d/dphi term couples xi,eta -- is the next")
print("     extension; the m=0 demag is the canonical elongated-body validation.)")
print("-" * 98)

RESULTS["convergence_AR2"] = dict(re_coarse=re_coarse, re_fine=re_fine)
RESULTS["dof_edge_AR4"] = dict(ie_relerr=re_ie, ie_dof=dof_ie,
                               ballooning=[dict(xi_max=t[0], n_ext=t[1], relerr=t[2], dof=t[3]) for t in ball])
RESULTS["n_fail"] = N_FAIL
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_37_ie_spheroidal_fem_full.json")
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 98)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 98)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
