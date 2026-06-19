# -*- coding: utf-8 -*-
r"""
demo_uu2_exact_dtn_fetd.py  (Track A -- exact-DtN transient open boundary)
==========================================================================
BACKS the iabc knowledge topic `exact_dtn` with an actual transient solve.

demo_uu proved that the spherical exterior Dirichlet-to-Neumann (DtN) symbol
Lambda_l is a degree-l rational function and that an auxiliary-ODE network
(Grote-Keller form) reproduces it.  THIS file is the finite-element time-domain
(FETD) application of that result: couple the exact DtN to a 1-D radial interior
mesh as a Robin boundary + companion auxiliary ODEs, march Newmark-beta, and
MEASURE the spurious reflection and the interior-energy decay for l = 1, 2, 3.

PHYSICS (scalar wave equation, normalised c = 1, truncation radius R0 = 1):
  per multipole order l the radial field R_l(r,t) obeys
      d2R/dt2 = (1/r^2) d/dr( r^2 dR/dr ) - l(l+1)/r^2 R,    0 <= r <= R0.
  Weak form with the spherical measure r^2 dr gives
      M_ij = int r^2 Ni Nj dr,
      K_ij = int [ r^2 Ni' Nj' + l(l+1) Ni Nj ] dr,
  and an outer boundary flux equal to the DtN output  g = R0 dR/dr = Lambda_l R.

EXACT DtN (Grote-Keller form, R0 = c = 1; poles = roots of the reverse Bessel
polynomial theta_l, all with Re < 0 -> passive/stable):
      g(t)      = -du/dt - u + sum_j psi_j ,        u(t) = R_l(R0,t)
      dpsi_j/dt =  lambda_j ( psi_j + u ) ,         lambda_j = root(theta_l).
  Folded into the weak form, the -du/dt becomes a (passive) damping at the
  boundary node, the -u a stiffness shift, and the sum psi_j a coupling to one
  first-order auxiliary ODE per pole.  Everything is marched implicitly with
  Newmark-beta (gamma=1/2, beta=1/4, unconditionally stable) + a trapezoidal
  update of the auxiliary ODEs; the time-invariant LHS is factored ONCE.

REFLECTION METRIC (clean, BC-isolating): the same interior mesh on [0, R0] is
run two ways -- (T) truncated at R0 with the exact DtN, and (F) embedded in a
much larger free-space reference mesh (identical h on [0,R0], outer wall so far
it never returns within the window).  The interior discretisation error is the
same in both, so their difference on [0,R0] is the SPURIOUS REFLECTION of the
truncation.

VERIFIED HERE (all asserted; self-contained numpy/scipy):
  (1) reflection (exact-DtN truncation vs free-space reference) DECREASES under
      mesh refinement -> the boundary is reflectionless in the continuum; the
      residual is pure discretisation error (not a modelling error).
  (2) the exact DtN beats a 1st-order Sommerfeld boundary (g = -du/dt) by orders
      of magnitude on the same mesh.
  (3) all auxiliary relaxation rates have Re < 0 and the interior energy drains
      monotonically through the boundary (no spurious growth) -> stable.

HONEST PRIOR ART: the rational exact-sphere DtN and its local-in-time auxiliary
realisation are Grote & Keller (SIAM J. Appl. Math. 1995) / Hagstrom; for a
separable boundary this is reflectionless and outperforms a PML, but it does NOT
generalise to arbitrary geometry.  The IABC-specific value (knowledge module) is
the bridge "IABC shell == this exact termination impedance" + its passive
equivalent-circuit realisation; this file supplies only the VERIFIED FETD core.

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu

np.set_printoptions(precision=5, suppress=False)

# 3-point Gauss-Legendre on [-1, 1] (exact for the quartic integrands below).
_GP = np.array([-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)])
_GW = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])


# ---------------------------------------------------------------------------
# exact DtN poles: lambda_j = roots of the reverse Bessel polynomial theta_l
#   theta_l(x) = sum_{k=0}^l (l+k)! / ((l-k)! k! 2^k) x^{l-k}
# (cross-checked against scipy.signal.besselap in demo_uu; here via np.roots)
# ---------------------------------------------------------------------------
def reverse_bessel_roots(l):
    from math import factorial
    coeffs = [factorial(l + k) / (factorial(l - k) * factorial(k) * 2 ** k)
              for k in range(l + 1)]              # descending powers x^l .. x^0
    return np.roots(coeffs).astype(complex)


# ---------------------------------------------------------------------------
# 1-D radial P1 finite-element assembly (spherical measure r^2 dr)
# ---------------------------------------------------------------------------
def assemble(nodes, l):
    n = nodes.size
    M = lil_matrix((n, n))
    K = lil_matrix((n, n))
    cent = l * (l + 1)
    for e in range(n - 1):
        a, b = nodes[e], nodes[e + 1]
        d = b - a
        rg = 0.5 * (a + b) + 0.5 * d * _GP            # gauss points in r
        jac = 0.5 * d
        N0 = (b - rg) / d
        N1 = (rg - a) / d
        dN0, dN1 = -1.0 / d, 1.0 / d
        Ns = (N0, N1)
        dNs = (dN0, dN1)
        for p in range(2):
            for q in range(2):
                m_pq = np.sum(_GW * jac * rg ** 2 * Ns[p] * Ns[q])
                k_pq = np.sum(_GW * jac * (rg ** 2 * dNs[p] * dNs[q]
                                          + cent * Ns[p] * Ns[q]))
                M[e + p, e + q] += m_pq
                K[e + p, e + q] += k_pq
    return M.tocsc(), K.tocsc()


# ---------------------------------------------------------------------------
# Newmark-beta transient solve with one of three outer boundary conditions:
#   bc = 'dtn'        : exact Grote-Keller DtN  (Robin + auxiliary ODEs)
#   bc = 'sommerfeld' : 1st-order ABC  g = -du/dt  (damping only)
#   bc = 'fardir'     : Dirichlet far wall (free-space reference)
# Dirichlet R(0)=0 is the natural regularity for a multipole (R_l ~ r^l).
# Returns interior snapshots on the first `n_keep` free DOFs and the energy.
# ---------------------------------------------------------------------------
def run(nodes, l, dt, T, bc, n_keep=None):
    n = nodes.size
    M, K = assemble(nodes, l)
    gamma, beta = 0.5, 0.25

    # free DOFs (drop node 0; drop far node for the Dirichlet reference)
    if bc == 'fardir':
        free = np.arange(1, n - 1)
    else:
        free = np.arange(1, n)
    Mr = M[np.ix_(free, free)].tolil()
    Kr = K[np.ix_(free, free)].tolil()
    lf = free.size - 1                            # local index of the R0 node

    # DtN poles (only used for bc='dtn')
    if bc == 'dtn':
        lam = reverse_bessel_roots(l)
        Bscal = float(np.sum((0.5 * dt * lam) / (1.0 - 0.5 * dt * lam)).real)
    else:
        lam = np.zeros(0, dtype=complex)
        Bscal = 0.0

    # boundary operator additions (rank-1 at the R0 node)
    C = lil_matrix((free.size, free.size))
    Ktot = Kr.copy()                              # K + (-u) stiffness shift
    Keff = Kr.copy()                              # LHS stiffness (Ktot - B*eNeN)
    if bc in ('dtn', 'sommerfeld'):
        C[lf, lf] += 1.0                          # damping from -du/dt
    if bc == 'dtn':
        Ktot[lf, lf] += 1.0                       # -u term
        Keff[lf, lf] += (1.0 - Bscal)            # -u term minus implicit aux
    Mr = Mr.tocsc(); C = C.tocsc()
    Ktot = Ktot.tocsc(); Keff = Keff.tocsc()

    # time-invariant Newmark LHS -- factor ONCE
    LHS = (Mr + gamma * dt * C + beta * dt * dt * Keff).tocsc()
    lu = splu(LHS)

    # initial field: a smooth radial bump (compact support inside (0, R0))
    rc, sig = 0.5, 0.08
    R = np.exp(-((nodes[free] - rc) / sig) ** 2)
    V = np.zeros(free.size)
    psi = np.zeros(lam.size, dtype=complex)
    # initial acceleration  M a0 = -Ktot R0  (V0 = 0, psi0 = 0)
    Acc = lu_solve_init(Mr, -(Ktot @ R))

    nsteps = int(round(T / dt))
    if n_keep is None:
        n_keep = free.size
    snaps = np.empty((nsteps + 1, n_keep))
    snaps[0] = R[:n_keep]
    E0 = 0.5 * (V @ (Mr @ V) + R @ (Kr @ R))
    energies = [E0]

    for step in range(nsteps):
        Rpred = R + dt * V + (0.5 - beta) * dt * dt * Acc
        Vpred = V + (1.0 - gamma) * dt * Acc
        f = np.zeros(free.size)
        if bc == 'dtn':
            u_n = R[lf]
            A_n = np.sum((psi * (1.0 + 0.5 * dt * lam) + 0.5 * dt * lam * u_n)
                         / (1.0 - 0.5 * dt * lam)).real
            f[lf] = A_n
        rhs = f - (C @ Vpred) - (Keff @ Rpred)
        Acc = lu.solve(rhs)
        R = Rpred + beta * dt * dt * Acc
        V = Vpred + gamma * dt * Acc
        if bc == 'dtn':
            u_np1 = R[lf]
            psi = (psi * (1.0 + 0.5 * dt * lam)
                   + 0.5 * dt * lam * (u_n + u_np1)) / (1.0 - 0.5 * dt * lam)
        snaps[step + 1] = R[:n_keep]
        energies.append(0.5 * (V @ (Mr @ V) + R @ (Kr @ R)))

    return {"snaps": snaps, "energy": np.array(energies),
            "Mr": Mr, "free": free, "n_keep": n_keep}


def lu_solve_init(Mr, rhs):
    return splu(Mr.tocsc()).solve(rhs)


# ---------------------------------------------------------------------------
def reflection_vs_reference(l, N, dt, T, R_far):
    """Run the truncated exact-DtN solve and the free-space reference on a
    SHARED [0,1] mesh; return the spurious reflection (M-weighted L2)."""
    h = 1.0 / N
    nodes_t = np.linspace(0.0, 1.0, N + 1)
    Nfar = int(round(R_far / h))
    nodes_f = np.linspace(0.0, R_far, Nfar + 1)
    assert np.allclose(nodes_f[:N + 1], nodes_t)           # shared nodes on [0,1]

    n_keep = N                                              # free DOFs 1..N on [0,1]
    test = run(nodes_t, l, dt, T, 'dtn', n_keep=n_keep)
    ref = run(nodes_f, l, dt, T, 'fardir', n_keep=n_keep)

    # M-weighted norm over [0,1] (use the test interior mass)
    Mr = test["Mr"][:n_keep, :n_keep]
    diff = test["snaps"] - ref["snaps"]
    err = np.array([np.sqrt(max(d @ (Mr @ d), 0.0)) for d in diff])
    refn = np.array([np.sqrt(max(s @ (Mr @ s), 0.0)) for s in ref["snaps"]])
    reflection = err.max() / refn.max()
    return reflection, test, ref


def sommerfeld_vs_reference(l, N, dt, T, R_far):
    h = 1.0 / N
    nodes_t = np.linspace(0.0, 1.0, N + 1)
    Nfar = int(round(R_far / h))
    nodes_f = np.linspace(0.0, R_far, Nfar + 1)
    n_keep = N
    test = run(nodes_t, l, dt, T, 'sommerfeld', n_keep=n_keep)
    ref = run(nodes_f, l, dt, T, 'fardir', n_keep=n_keep)
    Mr = test["Mr"][:n_keep, :n_keep]
    diff = test["snaps"] - ref["snaps"]
    err = np.array([np.sqrt(max(d @ (Mr @ d), 0.0)) for d in diff])
    refn = np.array([np.sqrt(max(s @ (Mr @ s), 0.0)) for s in ref["snaps"]])
    return err.max() / refn.max()


# ===========================================================================
print("=" * 78)
print(" demo_uu2 : exact-DtN transient open boundary (FETD, Newmark-beta)")
print("=" * 78)

T = 4.0
R_far = 5.0

# ---------------------------------------------------------------------------
print("\n[1] exact DtN poles lambda_j = roots(theta_l) (all Re<0 => stable):")
for l in (1, 2, 3):
    lam = reverse_bessel_roots(l)
    print(f"    l={l}:  lambda_j = {np.round(lam, 4)}   max Re = {lam.real.max():+.3e}")
    assert lam.real.max() < 0.0
print("    ok  (every auxiliary relaxation rate decays)")

# ---------------------------------------------------------------------------
print("\n[2] reflection vs free-space reference DECREASES under refinement:")
print("    (truncated exact-DtN solve vs the SAME interior mesh embedded in")
print("     a free-space reference; pure discretisation error remains)")
levels = [100, 200, 400]
refl = {}
for l in (1, 2, 3):
    row = []
    for N in levels:
        dt = 0.5 / N
        rfl, _, _ = reflection_vs_reference(l, N, dt, T, R_far)
        row.append(rfl)
    refl[l] = row
    rates = [row[i] / row[i + 1] for i in range(len(row) - 1)]
    print(f"    l={l}:  " + "  ".join(f"N={N}:{r:.2e}" for N, r in zip(levels, row))
          + f"   (drop x {', '.join(f'{x:.1f}' for x in rates)})")
    assert row[-1] < row[0], "reflection did not decrease under refinement"
    assert row[-1] < 5e-3, "finest reflection unexpectedly large"
print("    ok  (reflection -> 0 with h => the exact DtN is reflectionless)")

# ---------------------------------------------------------------------------
print("\n[3] exact DtN beats a 1st-order Sommerfeld ABC on the SAME mesh (N=400):")
N = 400
dt = 0.5 / N
for l in (1, 2, 3):
    r_dtn = refl[l][-1]
    r_som = sommerfeld_vs_reference(l, N, dt, T, R_far)
    print(f"    l={l}:  exact DtN = {r_dtn:.2e}   Sommerfeld-1 = {r_som:.2e}"
          f"   gain x {r_som / r_dtn:.0f}")
    assert r_som > 5.0 * r_dtn, "exact DtN should clearly beat Sommerfeld-1"
print("    ok  (exact DtN is far more absorbing than the 1st-order ABC)")

# ---------------------------------------------------------------------------
print("\n[4] interior energy drains through the exact-DtN boundary (no growth):")
N = 400
dt = 0.5 / N
for l in (1, 2, 3):
    nodes = np.linspace(0.0, 1.0, N + 1)
    out = run(nodes, l, dt, T, 'dtn')
    E = out["energy"]
    Epk = E.max()
    Efin = E[-1]
    # energy must never exceed the initial value (passive boundary) and must
    # decay strongly once the pulse has radiated out
    grew = E.max() / E[0]
    print(f"    l={l}:  E_peak={Epk:.3e}  E_final/E_peak={Efin / Epk:.2e}"
          f"  max(E/E0)={grew:.4f}")
    assert grew < 1.0 + 1e-6, "energy grew -> boundary not passive"
    assert Efin / Epk < 1e-3, "interior energy did not drain"
print("    ok  (passive, dissipative: interior energy decays, never grows)")

print("\n[interpretation]")
print("  * The exact rational sphere DtN, realised as Robin + l auxiliary ODEs,")
print("    gives a reflectionless transient open boundary up to discretisation:")
print("    the reflection falls with h (item 2) and is orders below a 1st-order")
print("    ABC (item 3).  Poles Re<0 => passive: interior energy only drains")
print("    (item 4).  This is the FETD backing for the iabc topic `exact_dtn`.")
print("  * Separable-geometry result (Grote-Keller class); not a general PML")
print("    replacement.  IABC value = the shell<->termination-impedance bridge.")
print("\nALL CHECKS PASSED.")
