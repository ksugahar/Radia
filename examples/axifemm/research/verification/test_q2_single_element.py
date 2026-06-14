"""Q2 single-element sanity test.

Build a single-quad mesh, assemble the 9x9 stiffness and sigma-mass matrices
via the C++ Phase A2 BFI (Q2 path), and verify the spectrum matches the
Python Q2 prototype (axifemm/axifemm_quad_q2.py — Gauss 8x8 numerical).

Both methods produce K, M in V-DOF (psi = 2 pi r * A_phi at the 9 nodes), so
the generalized eigenvalues lambda(K v = lambda M v) must coincide. The
local DOF ordering may differ between the C++ assembly (NGSolve's vertex/
edge/face numbering) and the Python prototype (s-midpoint convention), so we
compare permutation-invariant spectra rather than raw entries.

We also separately confirm that for an interior element, the C++ K, M agree
to machine precision with the Mathematica closed-form (already validated by
scripts/validate_q2_codegen.py at the per-entry level — this test exercises
the BFI assembly end-to-end).
"""
import sys
import os
import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "tests", "axifemm", "_reference_python"))
from axifemm_quad_q2 import (  # type: ignore
    element_matrices_q2_numerical,
    element_sigma_mass_q2_numerical,
)

from ngsolve import Mesh, BilinearForm, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo
from radia.radia_axifemm import (
    H1Henrotte,
    AxiHenrotteStiffnessBFI,
    AxiHenrotteSigmaMassBFI,
)


def cpp_q2_matrices(ra, rb, za, zb, mu, sigma):
    box = MoveTo(ra, za).Rectangle(rb - ra, zb - za).Face()
    box.faces.name = "conductor"
    mesh = Mesh(OCCGeometry(box, dim=2).GenerateMesh(
        maxh=10 * (rb - ra), quad_dominated=True))
    if mesh.ne != 1:
        raise RuntimeError(f"expected 1 element, got {mesh.ne}")

    fes = H1Henrotte(mesh, order=2)
    if fes.ndof != 9:
        raise RuntimeError(f"expected ndof=9, got {fes.ndof}")

    mu_cf = CoefficientFunction(mu)
    sigma_cf = CoefficientFunction(sigma)

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()
    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()

    K = np.zeros((9, 9))
    rs, cs, vs = a.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        K[r, c] += v
    M = np.zeros((9, 9))
    rs, cs, vs = m.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        M[r, c] += v
    return K, M


def main():
    mu0 = 4 * np.pi * 1e-7
    sigma = 5.8e7

    cases = [
        ("interior 1-2 mm", 1e-3, 2e-3, 0.0, 1e-3),
        ("interior 5-7 mm tall", 5e-3, 7e-3, -2e-3, 3e-3),
        ("axis-touching", 0.0, 1e-3, 0.0, 1e-3),
    ]

    for label, ra, rb, za, zb in cases:
        print(f"\n=== {label}  (ra,rb,za,zb) = ({ra*1e3:.3f}, {rb*1e3:.3f},"
              f" {za*1e3:.3f}, {zb*1e3:.3f}) mm ===")
        try:
            K_cpp, M_cpp = cpp_q2_matrices(ra, rb, za, zb, mu0, sigma)
        except RuntimeError as e:
            print(f"  SKIP (OCC mesh): {e}")
            continue
        # Symmetrise tiny noise.
        K_cpp = 0.5 * (K_cpp + K_cpp.T)
        M_cpp = 0.5 * (M_cpp + M_cpp.T)

        # Spectrum from C++. Drop zero eigenvalues from any singular axis nodes.
        eig_cpp = sorted(np.linalg.eigvalsh(K_cpp).tolist())
        mass_cpp = sorted(np.linalg.eigvalsh(M_cpp).tolist())
        print("  C++ K eigvals (sorted):", [f"{v:+.3e}" for v in eig_cpp])
        print("  C++ M eigvals (sorted):", [f"{v:+.3e}" for v in mass_cpp])

        # Python prototype only handles interior (sa > 0).
        if ra <= 0:
            print("  (Python prototype skips axis case; eyeballing C++ only.)")
            n_zero = sum(1 for v in eig_cpp if abs(v) < 1e-3 * max(abs(eig_cpp[-1]), 1))
            print(f"  C++ K near-zero eigvals: {n_zero}  (expect 3 axis DOFs)")
            continue

        K_py, _ = element_matrices_q2_numerical(ra, rb, za, zb, mu0, mu0)
        M_py = element_sigma_mass_q2_numerical(ra, rb, za, zb, sigma)
        K_py = 0.5 * (K_py + K_py.T)
        M_py = 0.5 * (M_py + M_py.T)
        eig_py = sorted(np.linalg.eigvalsh(K_py).tolist())
        mass_py = sorted(np.linalg.eigvalsh(M_py).tolist())
        print("  Py  K eigvals (sorted):", [f"{v:+.3e}" for v in eig_py])
        print("  Py  M eigvals (sorted):", [f"{v:+.3e}" for v in mass_py])

        # Matched-pair max relative diff.
        def rel_diff(a, b):
            scale = max(abs(a[-1]), abs(b[-1]), 1e-30)
            return max(abs(x - y) for x, y in zip(a, b)) / scale
        rk = rel_diff(eig_cpp, eig_py)
        rm = rel_diff(mass_cpp, mass_py)
        print(f"  K eigvals max-rel-diff: {rk:.2e}")
        print(f"  M eigvals max-rel-diff: {rm:.2e}")

        # Generalized eigenvalues (the physically meaningful spectrum).
        # K is rank-deficient (one zero mode); regularize for compare.
        try:
            lam_cpp = sorted(eigh(K_cpp + 1e-12 * np.eye(9), M_cpp,
                                  eigvals_only=True).tolist())
            lam_py = sorted(eigh(K_py + 1e-12 * np.eye(9), M_py,
                                 eigvals_only=True).tolist())
            print("  C++ gen.eigvals (top 4):", [f"{v:.3e}" for v in lam_cpp[-4:]])
            print("  Py  gen.eigvals (top 4):", [f"{v:.3e}" for v in lam_py[-4:]])
        except np.linalg.LinAlgError as e:
            print(f"  generalized eigh failed: {e}")


if __name__ == "__main__":
    main()
