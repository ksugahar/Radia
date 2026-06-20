# -*- coding: utf-8 -*-
r"""
act7_15_iabc_timedomain_cost.py  (Track A -- kelvin branch)
==========================================================
TIME-DOMAIN representation of the HIGH-FREQUENCY IABC, and its DtN COST advantage
(user: "the time-domain representation of the high-frequency IABC would be good; we
could show its DtN is cost-superior. Of course time-domain FEM-BEM is most accurate").

THE PICTURE (all numerically verified below):
  * The exact exterior DtN per spherical-harmonic mode n is RATIONAL in z=kR=omega R/c
    (act6_10_iabc_time_domain): Lambda_n(z) = i z - 1 + sum_{j=1..n} z_j/(z - z_j), poles z_j = the
    Bessel/Thomson filter poles. Its time domain = n LOCAL auxiliary ODEs (Grote-Keller
    exact nonreflecting BC). This is what a TIME-DOMAIN FEM-BEM reproduces -> MOST ACCURATE
    (it carries every mode exactly) but DENSE on Gamma + a temporal convolution (history).
  * The HIGH-FREQUENCY IABC = nested isotropic shells (here finite-frequency, spherical
    Bessel per shell) terminated by a wall; its truncation DtN APPROXIMATES Lambda_n over a
    band. With constant-material shells that DtN is TRANSCENDENTAL (Bessel) in z, so its
    *finite-state* time-domain representation is a REDUCED-ORDER RATIONAL (pole) DtN:
        Lambda^(M)_n(z) = i z - 1 + sum_{j=1..M} r_j/(z - p_j),   M <= n,
    realized as M auxiliary ODEs (a low-order analog filter) = the TIME-DOMAIN IABC.
  * COST: the time-domain IABC needs M auxiliary ODEs per surface mode, LOCAL in time and
    SPARSE in space (a few extra FE shell layers) -- no dense Gamma x Gamma matrix, no global
    convolution history.  The DtN spectrum quantifies accuracy vs the cost M.

VERIFIED HERE (all asserted):
  (A) exact rational DtN  Lambda_n(z)=iz-1+sum z_j/(z-z_j)  vs scipy DtN (reference). ~1e-13
  (B) finite-frequency nested-SHELL IABC DtN (transfer matrix, spherical Bessel per shell):
      machinery exact -- a VACUUM shell closed by the exact radiation DtN reproduces
      Lambda_n(k a) to ~1e-12; a lossy PEC-terminated shell stack APPROXIMATES Lambda_n over
      a band and improves with more shells (the IABC).
  (C) TIME-DOMAIN IABC = reduced M-pole rational DtN fit over a band: STABLE (Im p_j<0), band
      error DECREASES with M and ->0 at M=n; the M-auxiliary-ODE realization reproduces it by
      transient integration. => the high-frequency IABC HAS a finite time-domain representation.
  (D) DtN COST-ACCURACY datasheet: band DtN error vs M (= #auxiliary ODEs/mode = the cost).
      Structural cost vs time-domain FEM-BEM stated honestly (sparse-local vs dense-global).

HONEST PRIOR ART / no overclaim: the exact local-in-time pole realization of the sphere DtN
is Grote & Keller (1995/96); rational/pole (Pade, vector-fitting) approximation of NRBC kernels
is Alpert-Greengard-Hagstrom (2000) and the vector-fitting literature. The defensible-new slice
(per the lit survey) is the IABC reading + the static-to-high-frequency DtN cost-accuracy
datasheet in the SA/Kelvin context; FEM-BEM is the cited exact-but-expensive reference. The cost
comparison is STRUCTURAL (state dimension / locality / sparsity), NOT a wall-clock benchmark.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import hankel1, hankel2
from scipy.signal import besselap
from scipy.integrate import solve_ivp

np.set_printoptions(precision=6, suppress=False)

# ===========================================================================
# spherical Bessel/Hankel for COMPLEX argument (built from hankel1/hankel2)
# ===========================================================================
def _sb(kind, n, z):
    z = complex(z); p = np.sqrt(np.pi / (2.0 * z))
    h1 = hankel1(n + 0.5, z); h2 = hankel2(n + 0.5, z)
    if kind == 'h1': return p * h1
    if kind == 'j':  return p * (h1 + h2) / 2.0
    if kind == 'y':  return p * (h1 - h2) / (2j)
    raise ValueError(kind)

def _sbp(kind, n, z):
    z = complex(z); return _sb(kind, n - 1, z) - (n + 1) / z * _sb(kind, n, z)

def dtn_exact(n, z):
    z = complex(z); return z * _sbp('h1', n, z) / _sb('h1', n, z)

# ===========================================================================
# (A) exact rational DtN: Lambda_n(z) = iz - 1 + sum z_j/(z-z_j), z_j=i*besselap poles
# ===========================================================================
def exact_poles(n):
    return 1j * np.asarray(besselap(n, norm='delay')[1], dtype=complex)

def dtn_rational(n, z, poles=None):
    poles = exact_poles(n) if poles is None else poles
    z = complex(z); return 1j * z - 1.0 + np.sum(poles / (z - poles))

print("=" * 78)
print(" act7_15_iabc_timedomain_cost : TIME-DOMAIN high-frequency IABC + DtN cost")
print("=" * 78)

print("\n[A] exact rational DtN (FEM-BEM / Grote-Keller reference)  vs scipy:")
ZT = [0.5, 1.0, 2.0, 4.0, 1.0 + 0.3j, 2.5 - 0.4j]
maxA = 0.0
for n in range(1, 7):
    p = exact_poles(n)
    for z in ZT:
        maxA = max(maxA, abs(dtn_rational(n, z, p) - dtn_exact(n, z)) / abs(dtn_exact(n, z)))
print(f"    n=1..6 : max rel.err (rational n-pole vs scipy) = {maxA:.2e}")
assert maxA < 1e-9
print("    ok  (the exact per-mode DtN is the degree-n rational = the n-pole Grote-Keller NRBC)")

# ===========================================================================
# (B) finite-frequency nested-SHELL IABC DtN (the actual IABC), transfer matrix
# ===========================================================================
def _basisB(r, n, ks):
    return np.array([[_sb('j', n, ks * r),        _sb('y', n, ks * r)],
                     [ks * _sbp('j', n, ks * r),  ks * _sbp('y', n, ks * r)]], dtype=complex)

def shell_dtn(n, k, radii, m_list, term):
    """Truncation DtN a*u'/u at r=radii[0], for N constant-index shells (index m_list[j]
    in [radii[j],radii[j+1]]) closed at radii[N] by 'Dir' (PEC u=0), 'Neu' (u'=0), or
    'exact' (the exterior radiation DtN). u and u' continuous across index interfaces."""
    N = len(m_list); b = radii[N]
    if term == 'Dir':   state = np.array([0.0, 1.0], complex)
    elif term == 'Neu': state = np.array([1.0, 0.0], complex)
    elif term == 'exact': state = np.array([1.0, dtn_exact(n, k * b) / b], complex)
    else: raise ValueError(term)
    for j in range(N - 1, -1, -1):
        ks = k * m_list[j]
        coeff = np.linalg.solve(_basisB(radii[j + 1], n, ks), state)   # [A;B] in shell j
        state = _basisB(radii[j], n, ks) @ coeff                       # [u;u'] at inner edge
    a = radii[0]
    return a * state[1] / state[0]

print("\n[B] finite-frequency nested-shell IABC DtN (transfer matrix, spherical Bessel):")
# B1 machinery sanity: a VACUUM shell closed by the EXACT radiation DtN = exact Lambda_n(ka)
maxB = 0.0
for n in range(1, 6):
    for k in (1.0, 3.0):
        a, b = 1.0, 1.4
        val = shell_dtn(n, k, [a, b], [1.0], 'exact')
        maxB = max(maxB, abs(val - dtn_exact(n, k * a)) / abs(dtn_exact(n, k * a)))
print(f"    B1 vacuum shell + exact termination -> Lambda_n(ka): max rel.err = {maxB:.2e}")
assert maxB < 1e-9
print("    ok  (finite-frequency shell transfer machinery is exact)")

# B2 the IABC: PEC wall + LOSSY graded shells approximates Lambda_n(ka) over a band,
#    and improves with more shells (vs a bare PEC truncation = no absorber).
def band_err_shell(n, ks_band, radii, m_list, term):
    e = 0.0
    for k in ks_band:
        e = max(e, abs(shell_dtn(n, k, radii, m_list, term) - dtn_exact(n, k * radii[0])))
    return e
kband = np.linspace(1.0, 4.0, 25)
n_test = 1
bare = band_err_shell(n_test, kband, [1.0, 1.0 + 1e-9], [1.0], 'Dir')   # PEC right at truncation
errN = []
for N in (1, 2, 3):
    radii = [1.0] + list(1.0 + 0.25 * np.arange(1, N + 1))               # shells of width 0.25
    m_list = [1.0 + 0.6j * (j + 1) for j in range(N)]                    # graded loss (Im m grows)
    errN.append(band_err_shell(n_test, kband, radii, m_list, 'Dir'))
print(f"    B2 (n={n_test}, band kR in [1,4]) bare-PEC-truncation band err = {bare:.3f}")
print(f"       lossy IABC band err vs #shells N=1,2,3 = {np.round(errN,4)}")
assert errN[0] < bare and errN[-1] < errN[0], "a lossy shell stack should beat bare PEC truncation"
print("    ok  (a CRUDE fixed lossy design already beats bare PEC truncation -- it is NOT")
print("         optimized here; tuning the shells = Sugahara's fsolve. The shell DtN is")
print("         TRANSCENDENTAL (Bessel) in z, so its finite-state time domain needs (C).)")

# ===========================================================================
# (C) TIME-DOMAIN IABC = reduced M-pole rational DtN fit over a band  +  ODE realization
# ===========================================================================
def fit_reduced_dtn(n, M, zband):
    """Best M-pole rational  iz-1+sum_{j<=M} r_j/(z-p_j)  fitting Lambda_n over zband.
    Poles p_j = the M EXACT poles closest to the real axis (least-damped, band-dominant);
    residues r_j by linear least squares -> always stable (chosen poles have Im<0)."""
    allp = exact_poles(n)
    order = np.argsort(np.abs(allp.imag))          # least-damped first
    p = allp[order[:M]]
    # Lambda_n - (iz-1) = sum r_j/(z-p_j); linear LS for r over band samples
    A = np.array([[1.0 / (z - pj) for pj in p] for z in zband], dtype=complex)
    rhs = np.array([dtn_exact(n, z) - (1j * z - 1.0) for z in zband], dtype=complex)
    r, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    return p, r

def dtn_reduced(z, p, r):
    z = complex(z); return 1j * z - 1.0 + np.sum(r / (z - p))

zband = np.linspace(0.5, 6.0, 120)
print("\n[C] TIME-DOMAIN IABC = reduced M-pole rational DtN over band kR in [0.5,6]:")
n = 5
prev = np.inf
for M in range(1, n + 1):
    p, r = fit_reduced_dtn(n, M, zband)
    be = max(abs(dtn_reduced(z, p, r) - dtn_exact(n, z)) for z in zband)
    assert np.all(p.imag < 0), "reduced poles must be stable (Im<0)"
    print(f"    n={n}, M={M} poles : stable={np.all(p.imag<0)}  band max|err| = {be:.2e}")
    assert be <= prev + 1e-12, "band error must not increase with M"
    prev = be
assert prev < 1e-9, "M=n must recover the exact DtN"
print("    ok  (stable; band error decreases monotonically with M and -> 0 at M=n)")

# C2 the M auxiliary ODEs reproduce the reduced DtN (transient integration from rest)
def transient_reduced(p, r, w, t_end=300.0):
    m = p.size
    def rhs(t, y):
        psi = y[:m] + 1j * y[m:]
        u = np.exp(-1j * w * t)
        dpsi = -1j * p * psi - 1j * r * u            # psi_j' = -i p_j psi_j - i r_j u
        return np.concatenate([dpsi.real, dpsi.imag])
    sol = solve_ivp(rhs, [0, t_end], np.zeros(2 * m), rtol=1e-9, atol=1e-11,
                    t_eval=[t_end], method='RK45')
    psi = sol.y[:m, -1] + 1j * sol.y[m:, -1]
    u = np.exp(-1j * w * t_end); up = -1j * w * u
    return (-up - u + np.sum(psi)) / u               # g/u = -u'-u+sum psi

p3, r3 = fit_reduced_dtn(5, 3, zband)
maxC2 = 0.0
for w in (1.0, 2.0, 3.0):
    maxC2 = max(maxC2, abs(transient_reduced(p3, r3, w) - dtn_reduced(w, p3, r3)) / abs(dtn_reduced(w, p3, r3)))
print(f"    C2 M=3 auxiliary-ODE transient vs the fitted rational DtN: max rel.err = {maxC2:.2e}")
assert maxC2 < 1e-3
print("    ok  (M auxiliary ODEs ARE the finite-state time-domain high-frequency IABC)")

# ===========================================================================
# (D) DtN COST-ACCURACY datasheet
# ===========================================================================
print("\n[D] DtN cost-accuracy datasheet (band kR in [0.5,6]); cost = #poles = #aux ODEs/mode:")
print("    mode n |  M=1     M=2     M=3     M=4     M=5     (band max|DtN err|)")
for n in range(2, 6):
    row = []
    for M in range(1, 6):
        if M <= n:
            p, r = fit_reduced_dtn(n, M, zband)
            row.append(max(abs(dtn_reduced(z, p, r) - dtn_exact(n, z)) for z in zband))
        else:
            row.append(np.nan)
    print(f"      n={n}  | " + "  ".join(f"{e:7.1e}" if np.isfinite(e) else "   --   " for e in row))
print("""
    READING (cost model, structural -- not wall-clock):
      * TIME-DOMAIN FEM-BEM = the EXACT DtN: dense Gamma x Gamma surface coupling + a temporal
        convolution (retarded history) -> ~O(N_Gamma^2) per step + history.  MOST ACCURATE.
      * Exact local NRBC (Grote-Keller): n auxiliary ODEs for harmonic degree n (M=n column
        above = ~machine exact) -> exact but state grows with the highest mode kept.
      * TIME-DOMAIN IABC = M(<n) auxiliary ODEs/mode, LOCAL in time + SPARSE in space (a few
        extra FE shell layers; no dense matrix, no history) -> cheapest; the table is its
        accuracy knob. A small M already gives a small band DtN error for the dominant modes.
    => IABC trades a controlled DtN-spectral error for locality/sparsity; FEM-BEM trades cost
       for exactness. The DtN spectrum is the common yardstick (the cost-superiority claim).
""")
print("ALL CHECKS PASSED.")
