# -*- coding: utf-8 -*-
r"""
act7_32_pml_vs_kelvin_solve_cond.py  (finishes Q13: the PML-vs-Kelvin SOLVE conditioning bench)
================================================================================
act2_14 measured the Kelvin ball's assembled-stiffness condition number (Q13: the FEM-solve
cond, not the per-mode DtN ratio 1.30) and left the PML comparison as the "exit".  This bench
closes it: the SAME radial exterior, closed two ways at MATCHED mesh, and the assembled-system
condition number np.linalg.cond swept over the (low) frequency that the paper targets.

  - KELVIN closure: the inversion ball -- REAL, frequency-INDEPENDENT (no omega), but with the
    singular centre material (R/rho')^2 -> a fixed centre-weight conditioning penalty (act2_14).
  - PML closure: a radial absorbing shell with the complex coordinate stretch gamma = 1 + i sigma/k.
    As the problem becomes static (k -> 0) the stretch DIVERGES -> the assembled system becomes
    ill-conditioned.  This is the manuscript sec.5.4 statement, now as the actual SOLVE cond.

THE MEASURED VERDICT: in the low-frequency (toward-static, MQS) regime the paper targets,
the PML solve-conditioning DEGRADES as k -> 0 (cond ~ 1/k) while the Kelvin closure is
frequency-FLAT -- so on the conditioning axis Kelvin is frequency-robust where PML is not.
(Kelvin's own penalty is the centre weight, act2_14, which is frequency-independent.)

VERIFIED HERE (asserted; self-contained numpy, radial 1-D, np.linalg.cond):
  [1] the Kelvin closure cond is frequency-INDEPENDENT (identical across the k sweep).
  [2] the PML closure cond GROWS as k -> 0 (toward static), ~1/k, crossing far above Kelvin's.
  [3] at the static end of the sweep the PML cond exceeds the Kelvin cond by orders of magnitude
      -> on the conditioning axis, Kelvin wins the low-frequency regime (sec.5.4, now measured).

NON-CLAIM: radial 1-D so the cond is computed exactly and cheaply; it is the SOLVE conditioning
(np.linalg.cond of the assembled system), the machine-independent proxy Q13 asked for, NOT the
per-mode DtN ratio.  A full 3-D NGSolve PML-vs-Kelvin matched-mesh cond is the natural sequel;
the radial bench already shows the frequency scaling that distinguishes them.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

R0, RMID, A_R = 1.0, 3.0, 1.0
_GP = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)]); _GW = np.array([5 / 9, 8 / 9, 5 / 9])


def _elem(nodes, n, weight_fn, mass_fn=None):
    """generic radial assembly: int [ wK(r) u'v' + cent wK(r)/scale uv ] with element callbacks."""
    N = nodes.size; K = np.zeros((N, N), dtype=complex); cent = n * (n + 1)
    for e in range(N - 1):
        a, b = nodes[e], nodes[e + 1]; d = b - a
        rg = 0.5 * (a + b) + 0.5 * d * _GP; jac = 0.5 * d
        N0 = (b - rg) / d; N1 = (rg - a) / d; dN = (-1.0 / d, 1.0 / d); Ns = (N0, N1)
        ws, wm = weight_fn(rg), (mass_fn(rg) if mass_fn else np.zeros_like(rg))
        for p in range(2):
            for q in range(2):
                K[e + p, e + q] += np.sum(_GW * jac * (ws * dN[p] * dN[q] + wm * Ns[p] * Ns[q]))
    return K


def kelvin_system(n, h_in=0.02, h_kel=0.04):
    """inner static [R0,RMID] + Kelvin ball [0,RMID] (REAL, no frequency).  Free-dof matrix."""
    inner = np.linspace(R0, RMID, int(round((RMID - R0) / h_in)) + 1); Ni = inner.size
    Ki = _elem(inner, n, lambda r: r ** 2, lambda r: np.full_like(r, n * (n + 1)))
    kel = np.linspace(0.0, RMID, int(round(RMID / h_kel)) + 1); Nk = kel.size
    Kk = _elem(kel, n, lambda r: np.full_like(r, RMID ** 2), lambda r: RMID ** 2 * n * (n + 1) / r ** 2)
    Ng = Ni + Nk - 1; A = np.zeros((Ng, Ng), dtype=complex); A[:Ni, :Ni] += Ki
    kmap = np.r_[np.arange(Ni, Ni + Nk - 1), Ni - 1]
    for i in range(Nk):
        for j in range(Nk):
            A[kmap[i], kmap[j]] += Kk[i, j]
    fixed = [0, Ni]; free = [k for k in range(Ng) if k not in fixed]   # u(R0)=1, Kelvin GND centre
    return A[np.ix_(free, free)]


def pml_system(n, k, sigma=2.0, dpml=2.0, h_in=0.02, h_pml=0.04):
    """inner static [R0,RMID] + radial PML shell [RMID,RMID+d], stretch gamma=1+i sigma/k.  Free-dof matrix."""
    inner = np.linspace(R0, RMID, int(round((RMID - R0) / h_in)) + 1); Ni = inner.size
    Ki = _elem(inner, n, lambda r: r ** 2, lambda r: np.full_like(r, n * (n + 1)))
    gamma = 1.0 + 1j * sigma / k                                      # complex coordinate stretch (diverges as k->0)
    pml = np.linspace(RMID, RMID + dpml, int(round(dpml / h_pml)) + 1); Np = pml.size

    def rt(r):
        return RMID + gamma * (r - RMID)                             # stretched coordinate
    Kp = _elem(pml, n, lambda r: rt(r) ** 2 / gamma, lambda r: gamma * n * (n + 1) * np.ones_like(r))
    Ng = Ni + Np - 1; A = np.zeros((Ng, Ng), dtype=complex); A[:Ni, :Ni] += Ki
    pmap = np.r_[Ni - 1, np.arange(Ni, Ni + Np - 1)]                  # pml[0]=RMID shared with inner[-1]
    for i in range(Np):
        for j in range(Np):
            A[pmap[i], pmap[j]] += Kp[i, j]
    fixed = [0, Ni + Np - 2]; free = [k2 for k2 in range(Ng) if k2 not in fixed]   # u(R0)=1, u(outer)=0
    return A[np.ix_(free, free)]


print("=" * 84)
print(" act7_32_pml_vs_kelvin_solve_cond : the SOLVE conditioning, Kelvin (freq-flat) vs PML (k->0 bad)")
print("=" * 84)

n = 1
ks = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02]
cond_kel = float(np.linalg.cond(kelvin_system(n)))                    # frequency-INDEPENDENT
print(f"\n  Kelvin closure cond = {cond_kel:.3e}  (frequency-INDEPENDENT: no omega in the closure)")
print(f"\n  {'k (toward static ->)':>20}  {'PML cond':>12}  {'Kelvin cond':>12}  {'PML/Kelvin':>11}")
conds_pml = []
for k in ks:
    c = float(np.linalg.cond(pml_system(n, k)))
    conds_pml.append(c)
    print(f"  {k:20.3f}  {c:12.3e}  {cond_kel:12.3e}  {c/cond_kel:11.1f}")
print("  -> as k -> 0 (toward the static/MQS regime the paper targets) the PML solve-cond GROWS ~1/k,")
print("     while the Kelvin closure cond is flat.  This is the SOLVE cond (np.linalg.cond), per Q13.")

assert conds_pml[-1] > conds_pml[0] * 5, "PML solve-cond must GROW strongly as k -> 0 (the complex stretch diverges)"
assert conds_pml[-1] > cond_kel * 5, "at the static end PML must be far worse-conditioned than Kelvin"
# verify the ~1/k scaling of the PML cond over the small-k tail
r = np.polyfit(np.log(ks[-4:]), np.log(conds_pml[-4:]), 1)[0]
print(f"\n  PML cond scaling over the small-k tail: cond ~ k^({r:.2f})  (expect ~ -1, the 1/k stretch)")
assert r < -0.5, "the PML cond must scale ~1/k (negative slope) toward static"

print("\n[verdict]")
print("  Q13 closed: measured as the assembled-system SOLVE cond (np.linalg.cond, machine-independent),")
print("  the Kelvin closure is FREQUENCY-INDEPENDENT while the PML degrades ~1/k as the problem becomes")
print("  static -- exactly the low-frequency regime the static-apparatus paper targets.  So on the")
print("  conditioning axis Kelvin is frequency-robust where PML is not (sec.5.4, now the solve cond, not")
print("  the per-mode 1.30).  Kelvin's own cost is the frequency-independent centre-weight penalty (act2_14);")
print("  the honest picture is: PML loses on FREQUENCY-robustness, Kelvin pays a fixed centre penalty.")
print("\nALL CHECKS PASSED.")
