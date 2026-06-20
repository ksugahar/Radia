# -*- coding: utf-8 -*-
"""
Exterior MULTIPOLE SPECTRUM: magnetized SQUARE vs DISK (2D).

Uniform magnetization M = M * x_hat.  Magnetic surface charge sigma_m = M . n_hat.
Exterior magnetic scalar potential = 2D Coulomb potential of those surface charges:
    phi(P) = (1/2pi) * integral sigma(r') * (-ln|P - r'|) dl'

DISK radius a:   sigma = M cos(phi') on circle  -> pure 2D dipole (n=1 only).
SQUARE [-a,a]^2: sigma = +M on right edge (x=+a), -M on left edge (x=-a), 0 top/bottom.

On a bounding circle radius R, sample phi(R,theta), Fourier-decompose into cos(n theta)
amplitudes a_n.  Compare disk (pure n=1) vs square (higher odd harmonics).

KEY FINDING (corrects the naive task prediction):
  The task guessed the leading square correction is n=3 ~ (a/R)^2.  It is NOT.
  The charge config (+x-edge / -x-edge) has a D4 selection rule that KILLS n=3
  (and n=7,11,...).  The physically present harmonics are n = 1, 5, 9, 13, ...
  i.e. n == 1 (mod 4).  The LEADING correction is n=5, scaling as (a/R)^4.
  We prove n=3 is absent by showing it converges to 0 as edge resolution -> inf
  (it is pure O(h^2) trapezoid error), while n=5 is resolution-independent.

The general law |a_n/a_1| ~ (a/R)^(n-1) still holds -- but only for the harmonics
that the symmetry actually allows.
"""
import numpy as np

M = 1.0          # magnetization magnitude (ratios are scale-free)
a = 1.0          # half-size (square) / radius (disk)
N_THETA = 8192   # samples on bounding circle (power of 2 for clean FFT)
N_EDGE = 16000   # integration points per edge / arc (fine -> clean spectrum)


def _kernel(Px, Py, qx, qy):
    """2D Coulomb potential factor -ln|P - r'| (without 1/2pi or charge)."""
    dx = Px[:, None] - qx[None, :]
    dy = Py[:, None] - qy[None, :]
    return -np.log(np.sqrt(dx * dx + dy * dy))


def phi_square(R, theta, Ne=N_EDGE):
    """Exterior scalar potential of the magnetized square on circle radius R."""
    Px = R * np.cos(theta)
    Py = R * np.sin(theta)
    t = np.linspace(-a, a, Ne)
    dl = (2.0 * a) / (Ne - 1)
    w = np.full(Ne, dl); w[0] *= 0.5; w[-1] *= 0.5      # trapezoid weights
    phi = (+M) * (_kernel(Px, Py, np.full(Ne, +a), t) @ w)   # right edge +M
    phi += (-M) * (_kernel(Px, Py, np.full(Ne, -a), t) @ w)  # left edge  -M
    return phi / (2.0 * np.pi)


def phi_disk(R, theta, Ne=N_EDGE):
    """Exterior scalar potential of the magnetized disk on circle radius R.
    sigma(phi') = M cos(phi') line-charge density on circle radius a."""
    Px = R * np.cos(theta)
    Py = R * np.sin(theta)
    phip = np.linspace(0.0, 2.0 * np.pi, Ne, endpoint=False)
    dl = a * (2.0 * np.pi / Ne)                          # periodic midpoint rule
    qx = a * np.cos(phip); qy = a * np.sin(phip)
    sigma = M * np.cos(phip)
    return ((_kernel(Px, Py, qx, qy) * sigma[None, :]) @ np.full(Ne, dl)) / (2.0 * np.pi)


def cos_amps(phi, nmax=15):
    """cos(n theta) amplitudes a_n and sin amplitudes b_n from uniform samples."""
    F = np.fft.rfft(phi); Nt = len(phi)
    a = np.zeros(nmax + 1); b = np.zeros(nmax + 1)
    a[0] = F[0].real / Nt
    for n in range(1, nmax + 1):
        a[n] = 2.0 * F[n].real / Nt
        b[n] = -2.0 * F[n].imag / Nt
    return a, b


theta = np.linspace(0.0, 2.0 * np.pi, N_THETA, endpoint=False)
Rlist = [2.0 * a, 3.0 * a, 5.0 * a]

print("=" * 80)
print("DISK (M||x): exterior spectrum  -- expect PURE n=1 dipole")
print("=" * 80)
for R in Rlist:
    aN, bN = cos_amps(phi_disk(R, theta))
    a1 = aN[1]
    print(f"R={R/a:.0f}a  a1={a1: .6e}  a3/a1={abs(aN[3]/a1):.2e}  "
          f"a5/a1={abs(aN[5]/a1):.2e}  max|sin|/a1={np.max(np.abs(bN))/abs(a1):.2e}")
print("  -> disk a3/a1 ~ 1e-17 (MACHINE ZERO): pure dipole CONFIRMED.")

print()
print("=" * 80)
print("SQUARE (M||x): full odd-harmonic spectrum at R=2a")
print("=" * 80)
aN, bN = cos_amps(phi_square(2.0 * a, theta))
a1 = aN[1]
for n in range(1, 16, 2):
    tag = ""
    if n % 4 == 3:
        tag = "  <- ~noise floor (n=3 mod 4 KILLED by symmetry; ->0 as Ne->inf)"
    print(f"  n={n:2d}  a_n/a_1 = {abs(aN[n]/a1):.6e}{tag}")
print("  -> physical harmonics: n = 1, 5, 9, 13, ...  (n == 1 mod 4)")

print()
print("=" * 80)
print("PROOF n=3 is NOT physical: it is O(h^2) trapezoid error -> 0 as Ne grows")
print("=" * 80)
print("  (a5 is resolution-independent = physical; a3 falls ~4x per Ne doubling)")
for Ne in [1000, 2000, 4000, 8000, 16000]:
    aN, _ = cos_amps(phi_square(2.0 * a, theta, Ne=Ne))
    print(f"  Ne={Ne:6d}  a3/a1={abs(aN[3]/aN[1]):.4e}   a5/a1={abs(aN[5]/aN[1]):.8e}")

print()
print("=" * 80)
print("SQUARE leading-correction scaling: a5/a1 vs (a/R)^4  [law (a/R)^(n-1), n=5]")
print("=" * 80)
r51 = []
for R in Rlist:
    aN, _ = cos_amps(phi_square(R, theta))
    r = abs(aN[5] / aN[1]); r51.append(r)
    print(f"  R={R/a:.0f}a  a5/a1={r:.6e}   (4/15)*(a/R)^4={(4/15)*(a/R)**4:.6e}   "
          f"const a5/a1 / (a/R)^4 = {r/(a/R)**4:.5f}")
x = np.log(np.array([a / R for R in Rlist])); y = np.log(r51)
p, lnC = np.polyfit(x, y, 1)
print(f"  fitted power p = {p:.4f}  (analytic prediction p = 4),  "
      f"const C = {np.exp(lnC):.5f}  (= 4/15 = {4/15:.5f})")
print(f"  a5/a1 at R=2a = {r51[0]:.8f}  =  1/60 = {1/60:.8f}  (exact)")
print()
print("  ratio checks (a5/a1):")
print(f"    R=2a/R=3a = {r51[0]/r51[1]:.4f}  expect (3/2)^4 = {(3/2)**4:.4f}")
print(f"    R=2a/R=5a = {r51[0]/r51[2]:.4f}  expect (5/2)^4 = {(5/2)**4:.4f}")
print(f"    R=3a/R=5a = {r51[1]/r51[2]:.4f}  expect (5/3)^4 = {(5/3)**4:.4f}")

print()
print("=" * 80)
print("IMPLICATION FOR ELEMENT ORDER")
print("=" * 80)
print("  DISK/SPHERE  -> pure n=1 exterior field: order-1 (linear) element edge is")
print("                  exact for the boundary potential. No higher harmonics.")
print("  SQUARE/CUBOID-> exterior carries n=5,9,13,...; the boundary potential is")
print("                  NOT captured by a low order on a coarse bounding shell.")
print("                  To resolve up to harmonic n on a near boundary (small R/a)")
print("                  you need polynomial order >= ~n.  Sharp corners (square)")
print("                  inject the higher multipoles a round magnet never has.")
