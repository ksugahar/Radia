"""cln_sibc_cuboid_3d.py -- Python port of the lab CLN-SIBC (rank-(1,1) Mixed Galerkin)
3D cuboid polarizability, for the 3D non-axisymmetric levitation route.

Ports examples/CLN/scripts/cuboid3D_schur_F.wls (Mathematica) to Python and
reproduces its documented checks, so the 3D CLN-SIBC building block is
runnable in-repo (it was Mathematica-only).  Method (Sugahara-Nagamine-Hane
2026, "Mixed Galerkin Reduction for Eddy-Current Admittance:
CLN Bulk + HOIBC Surface Envelope via Schur Composition"):

  Foster admittance (scalar 3D diffusion, odd m,n,p):
      lam   = (m*pi/a)^2 + (n*pi/b)^2 + (p*pi/c)^2
      tau   = mu*sigma/lam,   g2 = 512*sigma*a*b*c / ((m*n*p)^2 * pi^6)
      Y(s)  = sum g2 / (1 + s*tau)                       (sigma absorbed in g2)
  CLN/Cauer reduction = Krylov basis on the modal operator K(s)=nu*Lam+sigma*s*I
  Schur-F SIBC termination (the paper's key result):
      Y_R(s) = Y_CLN_N(s) + K_SIBC*sqrt(s)/(s + d),   K_SIBC = S*sqrt(sigma/mu)
  Polarizability (the per-element BEM variable):
      alpha(s) = V - Y(s)/sigma,   alpha(0)=0, alpha(inf)=V, |alpha|<=V

This is the non-axisymmetric building block: a cuboid a!=b!=c gives a
direction-split response; here we port + verify the (scalar) Schur-F model.

Run:  python cln_sibc_cuboid_3d.py
"""
import numpy as np

# Geometry: Cu cuboid 5 x 2 x 1 mm (matches the lab script)
a, b, c = 5e-3, 2e-3, 1e-3
sigma = 5.8e7
mu = 4 * np.pi * 1e-7
nu = 1 / mu
V = a * b * c
S = 2 * (a * b + b * c + c * a)
K_SIBC = S * np.sqrt(sigma / mu)
MMAX = 13   # odd m,n,p up to 13 -> 7^3 = 343 modes


def foster_spectrum():
    odd = np.arange(1, MMAX + 1, 2)
    M, N, P = np.meshgrid(odd, odd, odd, indexing="ij")
    M, N, P = M.ravel(), N.ravel(), P.ravel()
    lam = (M * np.pi / a) ** 2 + (N * np.pi / b) ** 2 + (P * np.pi / c) ** 2
    tau = mu * sigma / lam
    g2 = 512.0 * sigma * a * b * c / ((M * N * P) ** 2 * np.pi ** 6)
    return lam, tau, g2


LAM, TAU, G2 = foster_spectrum()


def Y_foster(omega):
    s = 1j * omega
    return np.sum(G2 / (1.0 + s * TAU))


def cln_basis(s0, nstage):
    """K0-orthonormal Krylov basis of K(s0)^{-1} sigma b on the modal system."""
    Lam = LAM
    b = np.sqrt(G2 / TAU)                     # modal source (slab convention)
    K0 = nu * Lam + 0.0                       # K(0) diagonal = nu*lam
    Ks0 = nu * Lam + sigma * s0               # K(s0) diagonal
    Q = []
    vk = (sigma * b) / Ks0
    vk = vk / np.sqrt(abs(np.sum(K0 * vk * vk)))
    Q.append(vk)
    for _ in range(nstage - 1):
        vk = (sigma * Q[-1]) / Ks0
        for q in Q:                            # K0 modified Gram-Schmidt
            vk = vk - (np.sum(K0 * q * vk)) * q
        nrm = np.sqrt(abs(np.sum(K0 * vk * vk)))
        if nrm < 1e-300:
            break
        Q.append(vk / nrm)
    return np.array(Q)                         # (nstage x Nmodes)


def Y_cln(Q, omega):
    Lam = LAM
    b = np.sqrt(G2 / TAU)
    Kmod = nu * Lam + sigma * (1j * omega)     # diagonal modal operator
    Kr = Q @ (Kmod[:, None] * Q.T)             # Q K Q^T  (reduced, dense)
    Qb = Q @ b
    return sigma * (Qb @ np.linalg.solve(Kr, Qb))


def Y_schurF(Q, omega, d):
    return Y_cln(Q, omega) + K_SIBC * np.sqrt(1j * omega) / (1j * omega + d)


def alpha(Yval):
    return V - Yval / sigma


def main():
    print(f"cuboid {a*1e3}x{b*1e3}x{c*1e3} mm Cu;  V={V*1e9:.3f} mm^3 "
          f"S={S*1e6:.3f} mm^2;  modes={len(LAM)}")
    print(f"K_SIBC^3D = S*sqrt(sigma/mu) = {K_SIBC:.6e}")

    # --- check 1: DC admittance + polarizability limits ---
    y0 = Y_foster(0.0).real
    print(f"\n[check] Y(0) Foster sum   = {y0:.6e}")
    print(f"[check] alpha(0) = V-Y(0)/sigma = {alpha(y0).real:.3e}  (expect ~0)")
    yhf = Y_foster(2 * np.pi * 1e12)
    print(f"[check] alpha(1THz)->        = {alpha(yhf).real:.6e} m^3  (expect V={V:.3e})")

    # --- check 2: high-freq SIBC asymptote  |Y_R| sqrt(w)/K_SIBC -> 1 ---
    d = 2 * np.pi * 1e3
    QN = cln_basis(0.0, 6)
    print(f"\n[check] SIBC asymptote ratio |Y_R| sqrt(w)/K_SIBC (Schur-F, 6-stage):")
    for f in (1e6, 1e8, 1e10, 1e12):
        w = 2 * np.pi * f
        yR = Y_schurF(QN, w, d)
        print(f"   f={f:.0e} Hz  ratio = {abs(yR)*np.sqrt(w)/K_SIBC:.4f}")

    # --- check 3: CLN converges to Foster at a mid frequency ---
    f_mid = 1e5
    w = 2 * np.pi * f_mid
    yex = Y_foster(w)
    print(f"\n[check] CLN -> Foster at {f_mid:.0e} Hz (|Y|):")
    for n in (1, 2, 3, 4, 6, 8):
        Q = cln_basis(0.0, n)
        yc = Y_cln(Q, w)
        print(f"   {n:2d} stages  |Y_CLN|={abs(yc):.5e}  rel.err={abs(yc-yex)/abs(yex)*100:7.3f}%")


if __name__ == "__main__":
    main()
