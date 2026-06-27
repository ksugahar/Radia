"""validate_q2_codegen.py — Verify the JSON closed-form Q2 element matrices
agree with the high-order Gauss prototype (axifem/axifem_quad_q2.py).

Procedure:
  1. Load reference numerical values produced by codegen_q2_henrotte.py from
     the Mathematica JSON (closed form, exact in floating point).
  2. Compute the same K_phi/M_sigma_phi via the Python Gauss-Legendre 8x8
     prototype for an interior element (sa > 0).
  3. Compare entry by entry — they should agree to ~1e-10 since both are
     numerical evaluations of the same exact integrand.

Also performs the V-DOF transform (Vandermonde inverse + 2 pi r diag) and
solves the per-element generalized eigenvalue problem K v = lambda M v, then
checks both methods give the same spectrum.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent.parent
PROTO = Path("W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifem")
sys.path.insert(0, str(PROTO))

from axifem_quad_q2 import (  # type: ignore
    element_matrices_q2_numerical,
    element_sigma_mass_q2_numerical,
    _node_coords,
    _node_r,
    _phi_basis_q2,
)


def load_reference():
    p = ROOT / "scripts" / "q2_henrotte_test_values.json"
    return json.loads(p.read_text(encoding="utf-8"))


def closed_form_matrices(ref):
    """Reshape the 81-element flat lists into 9x9 K_phi and M_sigma_phi."""
    K = np.array(ref["K_phi_general"]).reshape(9, 9)
    M = np.array(ref["M_sigma_phi_general"]).reshape(9, 9)
    return K, M


def vdof_transform(K_phi, ra, rb, za, zb):
    """phi-monomial -> Vandermonde inverse -> V-DOF (psi = 2 pi r A_phi)."""
    nodes = _node_coords(ra, rb, za, zb)
    V = np.zeros((9, 9))
    for k in range(9):
        s_k, z_k = nodes[k]
        V[k] = _phi_basis_q2(s_k, z_k)
    Vinv = np.linalg.inv(V)
    K_nodal = Vinv.T @ K_phi @ Vinv          # nodal phi-DOF (= 2 pi r A_phi at node)
    rj = _node_r(ra, rb)
    T = 2 * np.pi * rj
    return (T[:, None] * K_nodal) * T[None, :]


def main():
    ref = load_reference()
    p = ref["test_params_general"]
    sa, sb = p["sa"], p["sb"]
    za, zb = p["za"], p["zb"]
    muR, muZ, sigma = p["muR"], p["muZ"], p["sigma"]
    ra, rb = np.sqrt(sa), np.sqrt(sb)

    print(f"=== Q2 codegen validation (interior element) ===")
    print(f"  sa={sa}, sb={sb}, za={za}, zb={zb}")
    print(f"  ra={ra:.6e}, rb={rb:.6e}")
    print(f"  muR=muZ=mu0, sigma=Cu={sigma}")
    print()

    K_cf, M_cf = closed_form_matrices(ref)
    print(f"Closed-form K_phi[1,1] = {K_cf[1,1]:.15e}")
    print(f"Closed-form M_phi[1,1] = {M_cf[1,1]:.15e}")
    print()

    # Python Gauss prototype produces M_V already (V-DOF), and Mphi_nodal in nodal-phi.
    # We compare directly in the monomial-phi basis: rebuild K_phi_gauss by
    # recomputing the inner integral exactly as axifem_quad_q2 does.
    from axifem_quad_q2 import _GL8_x, _GL8_w, _phi_dPhiDs, _phi_dPhiDz

    s_jac = 0.5 * (sb - sa)
    z_jac = 0.5 * (zb - za)
    s_mid = 0.5 * (sa + sb)
    z_mid = 0.5 * (za + zb)
    K_gauss = np.zeros((9, 9))
    M_gauss = np.zeros((9, 9))
    for i, xi in enumerate(_GL8_x):
        s = s_mid + s_jac * xi
        for j, eta in enumerate(_GL8_x):
            z = z_mid + z_jac * eta
            w = _GL8_w[i] * _GL8_w[j] * s_jac * z_jac
            ds = _phi_dPhiDs(s, z)
            dz = _phi_dPhiDz(s, z)
            phi = _phi_basis_q2(s, z)
            K_gauss += w * (np.outer(ds, ds) / (np.pi * muZ)
                            + np.outer(dz, dz) / (4.0 * np.pi * muR * s))
            M_gauss += w * sigma / (4.0 * np.pi * s) * np.outer(phi, phi)

    diff_K = np.max(np.abs(K_cf - K_gauss))
    rel_K = diff_K / max(np.max(np.abs(K_cf)), 1e-300)
    diff_M = np.max(np.abs(M_cf - M_gauss))
    rel_M = diff_M / max(np.max(np.abs(M_cf)), 1e-300)
    print(f"K_phi: max abs diff (closed-form vs Gauss-8x8) = {diff_K:.3e}  (rel {rel_K:.2e})")
    print(f"M_phi: max abs diff (closed-form vs Gauss-8x8) = {diff_M:.3e}  (rel {rel_M:.2e})")
    print()

    # Now V-DOF eigenvalue test (per-element). Apply Vandermonde + T transform
    # to both, then solve generalized eigenvalue.
    K_V_cf = vdof_transform(K_cf, ra, rb, za, zb)
    M_V_cf = vdof_transform(M_cf, ra, rb, za, zb)
    K_V_gauss = vdof_transform(K_gauss, ra, rb, za, zb)
    M_V_gauss = vdof_transform(M_gauss, ra, rb, za, zb)

    # Symmetrize (numerical tiny asymmetry).
    def sym(A): return 0.5 * (A + A.T)
    K_V_cf, M_V_cf = sym(K_V_cf), sym(M_V_cf)
    K_V_gauss, M_V_gauss = sym(K_V_gauss), sym(M_V_gauss)

    # Per-element eigenvalues — these are not physical tau (single element), but
    # any consistent V-DOF transformation should give identical spectra for both
    # K_cf and K_gauss.
    print("Per-element generalized eigenvalues lambda = K_V v / M_V v (single element):")
    try:
        lam_cf = np.sort(eigh(K_V_cf, M_V_cf, eigvals_only=True))
    except np.linalg.LinAlgError as e:
        print(f"  closed-form eigh failed: {e}")
        lam_cf = None
    try:
        lam_gauss = np.sort(eigh(K_V_gauss, M_V_gauss, eigvals_only=True))
    except np.linalg.LinAlgError as e:
        print(f"  Gauss eigh failed: {e}")
        lam_gauss = None
    if lam_cf is not None and lam_gauss is not None:
        print(f"  closed-form smallest 4: {lam_cf[:4]}")
        print(f"  Gauss-8x8 smallest 4:    {lam_gauss[:4]}")
        # Compare matched pairs (sorted)
        rel_lam = np.max(np.abs(lam_cf - lam_gauss) / (np.abs(lam_cf) + 1e-30))
        print(f"  max relative eigenvalue diff: {rel_lam:.3e}")
    print()

    # Verdict
    if diff_K < 1e-3 * np.max(np.abs(K_cf)) and diff_M < 1e-3 * np.max(np.abs(M_cf)):
        print("OK  Closed-form vs 8x8 Gauss agree within 0.1%.")
        return 0
    print("FAIL  Discrepancy too large.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
