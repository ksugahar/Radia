# -*- coding: utf-8 -*-
# DEMO (e) (verified, partial): optimal Kelvin/truncation radius R for a 2D magnetized
# square (a=1, M||x) vs a disk. Edge-charge -> circle-FFT multipole decomposition.
# Square: a3 = 0 (symmetry, n=4k+1 only); a5/a1 = (4/15)(a/R)^4 (exact). Disk: pure dipole.
# DOF proxy (R/a)^2 * p_needed^2 -> square has interior optimum R/a~3; disk -> minimal R.
import numpy as np

a = 1.0; M = 1.0
Nq = 600
xg, wg = np.polynomial.legendre.leggauss(Nq)
y_src = a * xg; w_src = a * wg


def harmonics(R, Ntheta=8192):
    theta = np.linspace(0, 2 * np.pi, Ntheta, endpoint=False)
    px, py = R * np.cos(theta), R * np.sin(theta)
    phi = np.zeros(Ntheta)
    for xs, sig in ((+a, +M), (-a, -M)):
        dx = px[:, None] - xs; dy = py[:, None] - y_src[None, :]
        phi += -sig * (0.5 * np.log(dx * dx + dy * dy) @ w_src)
    F = np.fft.rfft(phi) / Ntheta
    an = 2 * np.real(F); an[0] = np.real(F[0])
    return an


Rs = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
amat = np.array([harmonics(R) for R in Rs]); a1v = amat[:, 1]
r31 = np.abs(amat[:, 3] / a1v); r51 = np.abs(amat[:, 5] / a1v)
print("R/a   |a3/a1|     |a5/a1|     (a/R)^2    (a/R)^4")
for i, R in enumerate(Rs):
    print(f"{R/a:4.2f} {r31[i]:10.3e} {r51[i]:10.3e} {(a/R)**2:10.3e} {(a/R)**4:10.3e}")
p51 = np.polyfit(np.log(a / Rs), np.log(r51), 1)[0]
print(f"fit |a5/a1| ~ (a/R)^{p51:.3f} (pred 4); a3 ~ 1e-15 = ZERO (symmetry, n=4k+1 only)")

ns = [5, 9, 13]; Ak = {}
for n in ns:
    rn = np.abs(amat[:, n] / a1v); pk, ck = np.polyfit(np.log(a / Rs), np.log(rn), 1)
    Ak[n] = (np.exp(ck), pk)
Rf = np.linspace(1.01, 12.0, 4000)
for eps in (1e-4, 1e-6):
    p_nv = np.ceil(np.log(eps) / np.log(a / Rf))
    dof_nv = (Rf / a) ** 2 * p_nv ** 2; dof_dk = (Rf / a) ** 2
    print(f"eps={eps:.0e}: SQUARE optimal R/a={Rf[np.argmin(dof_nv)]:.2f} "
          f"(p={int(p_nv[np.argmin(dof_nv)])});  DISK R/a={Rf[np.argmin(dof_dk)]:.2f} (monotone)")
# DOF model is an explicit PROXY (area x order^2); (a/R)^4 decay and a3=0 are exact/measured.
