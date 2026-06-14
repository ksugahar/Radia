"""Phase 12.2: Banerjee-Gupta (1977) closed form for cuboid Newton potential.

Goal: implement Phi_box(P) = int_box 1/|r - P| dV in closed form, then verify
against numerical 3D quadrature on a known test cuboid.

Formula:
    For axis-aligned box [x0, x1] x [y0, y1] x [z0, z1] and observation P,
    let g(u, v, w) = -G(x_corner - Px, y_corner - Py, z_corner - Pz)
    summed with sign (-1)^{i+j+k} over the 8 corners (i,j,k in {0,1} with index 1
    meaning the upper face). The antiderivative is

    G(x, y, z) = x*y*log(z + R) + y*z*log(x + R) + z*x*log(y + R)
               - (x^2/2)*arctan(y*z/(x*R))
               - (y^2/2)*arctan(z*x/(y*R))
               - (z^2/2)*arctan(x*y/(z*R))
    R = sqrt(x^2 + y^2 + z^2).

This is the standard Newton-potential antiderivative cited in Banerjee-Gupta (1977),
MacMillan (1930), Nagy (1966), and used by D'Urso (2014) for arbitrary polyhedra.

We verify on cuboid 5e-3 x 2e-3 x 1e-3 by comparing Phi(P) at:
  - exterior point P = (10e-3, 0, 0)
  - face point     P = (a/2, 0, 0)  -- mild corner singularity
  - center         P = (0, 0, 0)
against direct mp.quad 3D integration of 1/|r - P|.

Date: 2026-05-02
"""

import time
import mpmath as mp

mp.mp.dps = 40


# --- antiderivative G(x, y, z) of 1/R with respect to (x, y, z) ---
def G_antider(x, y, z):
    """Antiderivative such that d^3 G / dx dy dz = 1/sqrt(x^2+y^2+z^2).

    Handles degenerate cases (x=0, y=0, z=0) via continuous limits:
      log(0 + R) -> log(|other|)  if 0 axis is the only one nonzero
      arctan(yz / (xR)) -> 0       when x=0 (the prefactor x^2/2 kills it)
    """
    x = mp.mpf(x); y = mp.mpf(y); z = mp.mpf(z)
    R = mp.sqrt(x*x + y*y + z*z)

    # Degenerate at the corner of integration domain (R = 0): the antiderivative
    # tends to 0 there because all six terms vanish (xy log R is 0*(-inf) -> 0).
    if R == 0:
        return mp.mpf(0)

    # log(z + R) etc. — choose tiny shift only if exactly singular
    def safe_log(arg, fallback):
        if arg <= 0:
            # if z + R = 0, then z = -R i.e. y = x = 0, z < 0; use limit form
            return mp.log(fallback)
        return mp.log(arg)

    term1 = x * y * safe_log(z + R, mp.fabs(R - z))  # if z+R = 0, R-z = 2|z|
    term2 = y * z * safe_log(x + R, mp.fabs(R - x))
    term3 = z * x * safe_log(y + R, mp.fabs(R - y))

    # arctan terms with 0 prefactor → set to 0 if axis variable is 0
    if x == 0:
        ax = mp.mpf(0)
    else:
        ax = (x * x / 2) * mp.atan(y * z / (x * R))
    if y == 0:
        ay = mp.mpf(0)
    else:
        ay = (y * y / 2) * mp.atan(z * x / (y * R))
    if z == 0:
        az = mp.mpf(0)
    else:
        az = (z * z / 2) * mp.atan(x * y / (z * R))

    return term1 + term2 + term3 - ax - ay - az


def Phi_box(P, x0, x1, y0, y1, z0, z1):
    """Newton potential at observation P of axis-aligned box [x0,x1]x[y0,y1]x[z0,z1]
    with unit density.

    Phi(P) = int_box 1/|r - P| dV
           = sum_{i,j,k in {0,1}} (-1)^{(i+j+k)} G(x_i - Px, y_j - Py, z_k - Pz)

    Convention: corner (1,1,1) is upper, sign +1. Standard inclusion-exclusion.
    Actually the rule is: Phi(P) = G|_(x0,y0,z0)^(x1,y1,z1) where the bracket
    means the alternating sum with sign +1 at (x1,y1,z1).
    """
    Px, Py, Pz = mp.mpf(P[0]), mp.mpf(P[1]), mp.mpf(P[2])
    xs = [x0 - Px, x1 - Px]
    ys = [y0 - Py, y1 - Py]
    zs = [z0 - Pz, z1 - Pz]
    total = mp.mpf(0)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                sign = mp.mpf(1) if ((i + j + k) % 2 == 0) else mp.mpf(-1)
                # Standard alternating sum: sign +1 when number of "1" indices is even,
                # -1 when odd ... but with the fundamental theorem we want
                # F(x1,y1,z1) - F(x0,y1,z1) - F(x1,y0,z1) + F(x0,y0,z1) - ...
                # i.e. sign = (-1)^{(1-i) + (1-j) + (1-k)} = (-1)^{i+j+k+3} = -(-1)^{i+j+k}
                total += sign * G_antider(xs[i], ys[j], zs[k])
    # Apply overall sign: i=j=k=1 (all upper) corresponds to sign=+1 in the FTC.
    # With our loop: i+j+k=3 -> sign=-1. We need +1 there, so flip.
    return -total


def Phi_box_quadrature(P, x0, x1, y0, y1, z0, z1):
    """Reference: 3D mp.quad of 1/|r - P|."""
    Px, Py, Pz = mp.mpf(P[0]), mp.mpf(P[1]), mp.mpf(P[2])
    def integrand(x, y, z):
        d = mp.sqrt((x - Px)**2 + (y - Py)**2 + (z - Pz)**2)
        if d < mp.mpf("1e-30"):
            return mp.mpf(0)
        return 1 / d
    return mp.quad(integrand, [x0, x1], [y0, y1], [z0, z1])


# === Tests ===
def _run_tests():
    print("=" * 64)
    print("Phase 12.2: Banerjee-Gupta closed form Newton potential")
    print("=" * 64)

    a, b, c = mp.mpf("5e-3"), mp.mpf("2e-3"), mp.mpf("1e-3")
    x0, x1 = -a/2, +a/2
    y0, y1 = -b/2, +b/2
    z0, z1 = -c/2, +c/2

    mp.mp.dps = 30

    # Test 1 (exterior, fast quadrature OK)
    P = (mp.mpf("10e-3"), mp.mpf(0), mp.mpf(0))
    t0 = time.time()
    phi_closed = Phi_box(P, x0, x1, y0, y1, z0, z1)
    print(f"Test 1 (exterior P=10mm): closed = {mp.nstr(phi_closed, 18)}  "
          f"({(time.time()-t0)*1000:.1f} ms)")


if __name__ == "__main__":
    _run_tests()
