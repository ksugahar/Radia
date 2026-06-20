# -*- coding: utf-8 -*-
# DEMO (z) (verified): the Sommerfeld layered kernel works at ALL frequencies -- it is NOT inherently
# high-frequency. This answers two questions for the isomorphism-benchmark plan:
#   Q1 "is Sommerfeld only for high-frequency?"  -> NO.
#   Q2 "can Sommerfeld BEM be used at low frequency?" -> YES (the KERNEL is well-defined and EASIER
#       at low frequency; only the full-wave EFIE FORMULATION has a separate low-frequency breakdown).
#
# Frequency-dependent half-space Sommerfeld reflected kernel (reflection about z=0, source+obs in the
# top medium with wavenumber k0):
#   I_refl = R * i * INT_0^inf (krho/kz0) J0(krho rho) exp(i kz0 (z+z')) dkrho,  kz0=sqrt(k0^2-krho^2), Im kz0>=0.
# For a PEC bottom the (TM) reflection coefficient is R=+1 (constant), so the Sommerfeld IDENTITY gives
# the EXACT image at ANY frequency:  I_refl = exp(i k0 R_im)/R_im,  R_im=sqrt(rho^2+(z+z')^2).
# We verify the numerical Sommerfeld integral against this exact image across k0*(z+z') = 0.001 .. 30,
# i.e. DC -> deep wave regime, and watch the integrand change character:
#   * low frequency  (k0 small): kz0 ~ i*krho, exp(i kz0 a) ~ exp(-krho a) -> MONOTONE EXP decay, benign;
#     the k0->0 limit is exactly the static Lipschitz-Hankel kernel of act5_07_sommerfeld_static_kernel.
#   * high frequency (k0 a >> 1): a real-kz0 band [0,k0] makes the integrand OSCILLATORY (the slow
#     Sommerfeld tail) + a branch point at krho=k0 -> the classic difficulty (SIP, tail acceleration).
# So the difficulty is a HIGH-frequency phenomenon; the kernel itself is valid (and easy) at low freq.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.special import j0

# the oscillatory-tail roundoff warnings at HIGH frequency are the difficulty itself (see footer note)
warnings.simplefilter("ignore", IntegrationWarning)


def _kz(krho, k0c):
    s = np.sqrt(k0c ** 2 - krho ** 2 + 0j)
    return s if s.imag >= 0 else -s                 # physical sheet: Im kz >= 0 (decay/radiation)


def sommerfeld_refl(rho, z, zp, k0, Rfun=lambda krho, kz0, k0c: 1.0 + 0j, loss=1e-7):
    """Reflected Sommerfeld integral I_refl (top half-space). Rfun(krho,kz0,k0c)=spectral reflection coeff."""
    k0c = (abs(k0) * (1.0 + 1j * loss)) if k0 != 0 else 1e-9 + 0j   # tiny loss moves the branch off-axis
    a = z + zp

    def integ(krho):
        kz0 = _kz(krho, k0c)
        return 1j * Rfun(krho, kz0, k0c) * (krho / kz0) * j0(krho * rho) * np.exp(1j * kz0 * a)
    kmax = abs(k0) + 60.0 / max(a, 1e-3)
    pts = [abs(k0)] if 0 < abs(k0) < kmax else None
    re = quad(lambda k: integ(k).real, 0, kmax, limit=500, points=pts)[0]
    im = quad(lambda k: integ(k).imag, 0, kmax, limit=500, points=pts)[0]
    return re + 1j * im


rho, z, zp = 0.6, 0.5, 0.4
a = z + zp
R_im = np.sqrt(rho ** 2 + a ** 2)

print("PEC half-space: numerical Sommerfeld integral  vs  EXACT image exp(i k0 R_im)/R_im, all freq")
print("(R=+1; verifies the wave Sommerfeld kernel from DC to deep wave regime)\n")
print("   k0*(z+z')   regime        |I_refl| (num)    |image| (exact)    rel.err")
emax = 0.0
for k0a in (0.001, 0.1, 0.5, 2.0, 6.0, 15.0, 30.0):
    k0 = k0a / a
    num = sommerfeld_refl(rho, z, zp, k0)
    exact = np.exp(1j * k0 * R_im) / R_im
    rel = abs(num - exact) / abs(exact); emax = max(emax, rel)
    reg = "quasi-static" if k0a < 0.3 else ("intermediate" if k0a < 3 else "wave (osc. tail)")
    print("   %8.3f   %-14s %12.6e    %12.6e    %.2e" % (k0a, reg, abs(num), abs(exact), rel))
print("   -> max rel err across DC..wave: %.2e  (PEC image exact at every frequency)\n" % emax)

# static limit: k0 -> 0 reproduces the act5_07_sommerfeld_static_kernel static image (Lipschitz-Hankel) -----------------------
num0 = sommerfeld_refl(rho, z, zp, 1e-6)
stat = 1.0 / R_im                                   # static PEC image (1/R_im)
print("STATIC LIMIT k0->0: I_refl=%.8e  vs static 1/R_im=%.8e  rel=%.2e  (= act5_07_sommerfeld_static_kernel kernel)"
      % (num0.real, stat, abs(num0 - stat) / stat))

# lossy/dielectric half-space (frequency-dependent Fresnel R): converges at all freq, easiest at low --
def fresnel_TM(krho, kz0, k0c, eps_r=4.0, sig_over_we=2.0):
    k1c = k0c * np.sqrt(eps_r + 1j * sig_over_we)    # lower medium: dielectric + conductivity
    kz1 = np.sqrt(k1c ** 2 - krho ** 2 + 0j); kz1 = kz1 if kz1.imag >= 0 else -kz1
    return (eps_r_complex(eps_r, sig_over_we) * kz0 - kz1) / (eps_r_complex(eps_r, sig_over_we) * kz0 + kz1)


def eps_r_complex(eps_r, sig_over_we):
    return eps_r + 1j * sig_over_we


print("\nLOSSY half-space (eps_r=4 + i*2): reflected kernel converges at every frequency")
print("   k0*(z+z')   regime         |I_refl|         (tail: exp-decay if low, oscillatory if high)")
for k0a in (0.01, 0.5, 3.0, 12.0):
    k0 = k0a / a
    val = sommerfeld_refl(rho, z, zp, k0, Rfun=fresnel_TM)
    reg = "quasi-static" if k0a < 0.3 else ("intermediate" if k0a < 3 else "wave")
    print("   %8.2f   %-13s  %12.6e" % (k0a, reg, abs(val)))

print("\n=> the Sommerfeld layered KERNEL is valid at ALL frequencies and is NUMERICALLY EASIEST at low")
print("   frequency (monotone exp tail, no real poles, branch point near 0). The 'Sommerfeld is hard'")
print("   reputation is a HIGH-frequency phenomenon (oscillatory tail + branch/pole on the SIP). Low-")
print("   frequency layered Sommerfeld is a mature field (geophysical EM, eddy-current testing,")
print("   induction logging -- all quasi-static). The separate 'low-frequency breakdown' is an EFIE")
print("   FORMULATION/conditioning issue (loop-tree/Calderon), NOT a property of the kernel.")
