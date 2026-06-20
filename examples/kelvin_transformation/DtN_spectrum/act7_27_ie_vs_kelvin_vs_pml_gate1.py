# -*- coding: utf-8 -*-
"""
act7_27_ie_vs_kelvin_vs_pml_gate1.py  (Act 7 -- GATE 1: should we build a C++ infinite element?)
================================================================================================
The go / no-go benchmark BEFORE committing to a pybind11 C++ infinite element (IE).  The honest
question is NOT "IE vs Kelvin" (Kelvin is sphere-locked by Liouville, so it loses on elongated
geometry for an obvious reason) -- it is "IE vs **box-PML**", because NGSolve ALREADY ships a
non-spherical PML (`pml.Cartesian` / `BrickRadial`) that ALSO escapes the sphere-lock.  A NEW
element is only justified if the IE beats box-PML on the axes where PML is weak.

HONEST OUTCOME (this script self-asserts the MEASURED facts, not the hoped-for ones):

  (a) GEOMETRY -- exterior DOF vs aspect ratio AR = L/d (cylinder body, fixed h):  CLEAN SIGNAL.
        Kelvin  ~ AR^2  (Liouville sphere-lock: must enclose the body in a sphere, R ~ L/2)
        box-PML ~ AR    (box hugging the body)
        IE      ~ AR    (surface wrap, P radial DOFs, no volume -- leanest)
      => Kelvin is OUT for elongated/planar bodies (DOF ~ AR worse).  But IE only TIES box-PML.

  (b) STATIC/LOW-FREQ -- where the IE would have to BEAT box-PML.  TWO findings KILL the easy case:
      (b1) IE spectral accuracy: EXACT for n <= P-1 (act7_25) -- the IE's one genuine unique edge.
      (b2) BUT the naive reciprocal-power IE basis (a/r)^k is ITSELF ill-conditioned: the matrix
           A_kl = a(kl + n(n+1))/((k+l)-1) is Hilbert/Cauchy-like, so cond EXPLODES with the decay
           order P (the accuracy knob).  A production IE would need an ORTHOGONALIZED basis first.
      (b3) the cheap 1-D PML proxy did NOT reproduce the documented DC conditioning blow-up
           (act6_09/act7_24) -- it is interior-mode-dominated at this fidelity.  So the cheap
           proxy does NOT establish an IE conditioning win over box-PML.

VERDICT: the cheap Gate-1 proxies do NOT justify building a C++ IE.  Kelvin loses on geometry, but
the IE merely TIES the already-shipped box-PML there; the IE's only proven unique edge is spectral
exactness (niche), while its naive basis is ill-conditioned.  Honest call = NO clean GO: either
(i) accept box-PML (NGSolve already ships it -- "complement, do not reimplement"), or (ii) first fix
the IE basis (orthogonalize) AND settle IE-vs-box-PML on a REAL 3-D model (Gate 2) before any C++.

Pure numpy.  A measurement to ALLOCATE EFFORT honestly, not a new method.
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


# ============================================================================
# (a) GEOMETRY -- exterior DOF vs aspect ratio (fixed body diameter d, mesh size h)
# ============================================================================
def dof_kelvin(AR, d=1.0, h=0.25, nrad=4):
    """Kelvin: enclose the L x d body in a bounding SPHERE (Liouville), mesh the inverted ball.
    exterior DOF ~ sphere surface / h^2 x radial layers.  R ~ L/2 -> DOF ~ AR^2."""
    L = AR * d
    R = 0.5 * np.sqrt(L ** 2 + d ** 2)
    return (4 * np.pi * R ** 2) / h ** 2 * nrad


def dof_pml(AR, d=1.0, h=0.25, nlayer=6):
    """box-PML: square-section box hugging the body; absorbing layer = surface x nlayer.
    box surface ~ 4 d L -> DOF ~ AR."""
    L = AR * d
    return (2 * d * d + 4 * d * L) / h ** 2 * nlayer


def dof_ie(AR, d=1.0, h=0.25, P=4):
    """IE: wrap the body surface (cylinder) with P radial decay DOFs -- NO meshed volume.
    cylinder surface ~ pi d L -> DOF ~ AR (leanest)."""
    L = AR * d
    return (np.pi * d * L + 2 * np.pi * (d / 2) ** 2) / h ** 2 * P


print("=" * 92)
print(" act7_27 GATE 1 : IE vs Kelvin vs box-PML -- should we build a C++ infinite element?")
print("=" * 92)

AR = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
dK = np.array([dof_kelvin(a) for a in AR])
dP = np.array([dof_pml(a) for a in AR])
dI = np.array([dof_ie(a) for a in AR])

print("\n(a) exterior DOF vs aspect ratio AR = L/d  (cylinder body, d=1, h=0.25):")
print("   AR      Kelvin       box-PML        IE        Kelvin/IE")
for i, a in enumerate(AR):
    print(f"  {a:5.0f}  {dK[i]:11.0f}  {dP[i]:11.0f}  {dI[i]:11.0f}   {dK[i]/dI[i]:8.1f}x")

m = AR >= 5
sK = float(np.polyfit(np.log(AR[m]), np.log(dK[m]), 1)[0])
sP = float(np.polyfit(np.log(AR[m]), np.log(dP[m]), 1)[0])
sI = float(np.polyfit(np.log(AR[m]), np.log(dI[m]), 1)[0])
print(f"\n  asymptotic exponent (DOF ~ AR^p):  Kelvin {sK:.2f},  box-PML {sP:.2f},  IE {sI:.2f}")

check("(a) Kelvin DOF ~ AR^2 (Liouville sphere-lock penalty)", 1.8 <= sK <= 2.2, f"slope {sK:.2f}")
check("(a) box-PML DOF ~ AR^1 (escapes sphere-lock)", 0.8 <= sP <= 1.2, f"slope {sP:.2f}")
check("(a) IE DOF ~ AR^1 (escapes sphere-lock, leanest)", 0.8 <= sI <= 1.2, f"slope {sI:.2f}")
check("(a) Kelvin loses big at high AR (AR=50: Kelvin/IE > 20x)", dK[-1] / dI[-1] > 20, f"{dK[-1]/dI[-1]:.0f}x")
check("(a) IE only TIES box-PML on geometry (same AR^1 scaling; IE a bit leaner)",
      abs(sI - sP) < 0.2 and bool(np.all(dI < dP)), f"slopes IE {sI:.2f} vs PML {sP:.2f}")

# ============================================================================
# (b) STATIC/LOW-FREQ -- where the IE would have to BEAT box-PML
# ============================================================================
def ie_matrix(n, P, a=1.0):
    k = np.arange(1, P + 1)
    return a * (np.outer(k, k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)


def ie_dtn(n, P, a=1.0):
    A = ie_matrix(n, P, a)
    return -1.0 / (a * (np.ones(P) @ np.linalg.solve(A, np.ones(P))))


# --- (b1) IE spectral accuracy: exact for n <= P-1 (the one genuine unique edge) ---
print("\n(b1) IE spectral accuracy (static DtN, exact = -(n+1)) -- the IE's unique edge:")
P_IE = 4
dn = [abs(ie_dtn(n, P_IE) - (-(n + 1))) / (n + 1) for n in range(P_IE)]
print(f"   P={P_IE}: reldef n=0..{P_IE-1} = " + ", ".join(f"{x:.1e}" for x in dn))
check(f"(b1) IE EXACT for n <= P-1 (P={P_IE}: d_n < 1e-10)", max(dn) < 1e-10, f"max {max(dn):.1e}")

# --- (b2) BUT the naive reciprocal-power IE basis is Hilbert-like ill-conditioned (grows with P) ---
print("\n(b2) IE basis conditioning -- the naive (a/r)^k basis is Hilbert/Cauchy-like:")
Ps = [2, 4, 6, 8]
cond_ie_P = [float(np.linalg.cond(ie_matrix(1, P))) for P in Ps]
print("    P        cond(IE matrix)")
for P, c in zip(Ps, cond_ie_P):
    print(f"   {P:2d}        {c:.3e}")
print(f"  -> cond grows {cond_ie_P[-1]/cond_ie_P[0]:.1e}x from P=2 to P=8 (the decay order is the accuracy knob)")
check("(b2) naive IE basis is ILL-CONDITIONED and grows with P (P=4 cond > 1e3)", cond_ie_P[1] > 1e3,
      f"cond(P=4) {cond_ie_P[1]:.1e}")
check("(b2) IE conditioning EXPLODES with the accuracy knob P (P=8 cond > 1e6)", cond_ie_P[-1] > 1e6,
      f"cond(P=8) {cond_ie_P[-1]:.1e}  => a production IE needs an ORTHOGONALIZED basis")


def cond_pml_radial(k, n=1, a=1.0, Rin=1.7, R=2.0, N=60, sig0=8.0):
    """1-D P1 FE radial spherical Helmholtz (mode n) on [a,R] with a PML complex stretch
    gamma = 1 + i sigma/k in [Rin,R].  A cheap proxy for the PML conditioning vs frequency."""
    x = np.linspace(a, R, N)

    def sig(xx):
        return sig0 * ((xx - Rin) / (R - Rin)) ** 2 if xx >= Rin else 0.0

    cum = np.zeros(N)
    for i in range(1, N):
        cum[i] = cum[i - 1] + 0.5 * (sig(x[i]) + sig(x[i - 1])) * (x[i] - x[i - 1])
    rtil = x + 1j * cum / k
    K = np.zeros((N, N), dtype=complex)
    for e in range(N - 1):
        he = x[e + 1] - x[e]
        xm = 0.5 * (x[e] + x[e + 1])
        gm = 1 + 1j * sig(xm) / k
        rm = 0.5 * (rtil[e] + rtil[e + 1])
        Ke = (rm ** 2 / gm) / he * np.array([[1, -1], [-1, 1]], dtype=complex)
        Me = (n * (n + 1) - k ** 2 * rm ** 2) * gm * he / 6.0 * np.array([[2, 1], [1, 2]], dtype=complex)
        K[e:e + 2, e:e + 2] += Ke + Me
    return float(np.linalg.cond(K[1:, 1:]))


# --- (b3) the cheap PML proxy does NOT resolve the conditioning question (honest, inconclusive) ---
print("\n(b3) cheap 1-D PML conditioning proxy (NOT the full 3-D system) -- honest, inconclusive:")
ks = [1.0, 0.3, 0.1, 0.03, 0.01]
cond_pml = [cond_pml_radial(k) for k in ks]
for k, c in zip(ks, cond_pml):
    print(f"   k={k:5.2f}   cond(PML proxy) = {c:.3e}")
pml_growth = cond_pml[-1] / cond_pml[0]
print(f"  -> growth k=1->0.01 = {pml_growth:.2f}x : this 1-D proxy is interior-mode-dominated and")
print("     does NOT reproduce the documented full-PML DC blow-up (act6_09/act7_24).")
check("(b3) cheap PML proxy is INCONCLUSIVE (no DC blow-up at this 1-D fidelity)", pml_growth < 10,
      f"growth {pml_growth:.2f}x  => IE-vs-PML conditioning is UNRESOLVED -> needs Gate 2 (real 3-D)")

# ============================================================================
# VERDICT  (honest -- the cheap proxies do NOT greenlight a C++ IE)
# ============================================================================
kelvin_out = (1.8 <= sK <= 2.2) and (dK[-1] / dI[-1] > 20)          # geometry: Kelvin loses
ie_ties_pml_geom = abs(sI - sP) < 0.2                                # IE only ties box-PML
ie_unique_edge_spectral = max(dn) < 1e-10                            # IE's one proven edge
ie_basis_illconditioned = cond_ie_P[-1] > 1e6                        # IE weakness
pml_cond_unresolved = pml_growth < 10                               # proxy inconclusive
clean_go = kelvin_out and ie_ties_pml_geom and ie_unique_edge_spectral and \
    (not ie_basis_illconditioned) and (not pml_cond_unresolved)     # = False (honest)

print("\n" + "-" * 92)
print(" VERDICT (Gate 1) -- elongated/planar body + static/low-freq + air exterior:")
print("   Kelvin   : OUT  -- exterior DOF ~ AR^2 (Liouville sphere-lock); ~AR x worse than IE/PML.")
print("   box-PML  : escapes geometry (DOF ~ AR) AND is already shipped in NGSolve.")
print("   IE       : escapes geometry (DOF ~ AR, leanest) + spectral-exact (n<=P-1) -- BUT its naive")
print("              (a/r)^k basis is Hilbert-ill-conditioned (cond explodes with P), and the cheap")
print("              proxy does NOT prove an IE conditioning win over box-PML.")
print(f"   => {'GO' if clean_go else 'NO CLEAN GO'}: the cheap Gate-1 proxies do NOT justify a C++ IE over the shipped box-PML.")
print("      Options: (i) use box-PML (complement NGSolve, do not reimplement); OR")
print("               (ii) FIRST orthogonalize the IE basis AND settle IE-vs-box-PML on a REAL 3-D")
print("                    model (Gate 2) -- only then is a pybind11 C++ IE warranted.")
print("-" * 92)

RESULTS = {
    "geometry": {"AR": AR.tolist(), "dof_kelvin": dK.tolist(), "dof_pml": dP.tolist(),
                 "dof_ie": dI.tolist(), "exp_kelvin": sK, "exp_pml": sP, "exp_ie": sI,
                 "kelvin_over_ie_at_AR50": float(dK[-1] / dI[-1])},
    "static": {"ie_spectral_reldef": dn, "ie_basis_cond_vs_P": dict(zip([str(p) for p in Ps], cond_ie_P)),
               "pml_proxy_cond_vs_k": dict(zip([str(k) for k in ks], cond_pml)),
               "pml_proxy_growth": pml_growth},
    "findings": {
        "kelvin_out_on_geometry": bool(kelvin_out),
        "ie_only_ties_box_pml_on_geometry": bool(ie_ties_pml_geom),
        "ie_unique_edge_is_spectral_exactness": bool(ie_unique_edge_spectral),
        "naive_ie_basis_ill_conditioned_needs_orthogonalization": bool(ie_basis_illconditioned),
        "ie_vs_pml_conditioning_unresolved_by_cheap_proxy": bool(pml_cond_unresolved),
    },
    "verdict_clean_go_build_cpp_ie": bool(clean_go),
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_27_ie_vs_kelvin_vs_pml_gate1.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 92)
print(" ALL CHECKS PASSED (the asserted facts are the MEASURED ones)" if N_FAIL == 0
      else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 92)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
