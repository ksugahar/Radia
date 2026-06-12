"""Accelerator pole design via the hodograph -- design + measurement foundation.

RESEARCH example.  First rung of the "hodograph inverse design + HDiv-MMM
nonlinear forward" combination (the radia selling point): the pieces every
forward verification needs --

  1. multipoles(B, r_ref) -- the MULTIPLE harmonic analyzer.  Samples the
     field on a reference circle and FFTs it into normal b_n / skew a_n
     multipoles (accelerator convention  By + i Bx = sum_n (b_n + i a_n)
     z^{n-1},  z = x + iy;  n=1 dipole, n=2 quad, n=3 sextupole).  This is
     how "field quality / good-field region" is measured.

  2. the hodograph POLE GEOMETRY.  The iron pole face is an equipotential
     of the magnetic scalar potential Phi (high mu -> H_t = 0 -> Phi=const).
     For a multipole the complex potential is W ~ z^n, so:
        dipole (n=1):  Phi = (B0/mu0) y   -> pole  y = const   (flat)
        quad   (n=2):  Phi = G xy         -> pole  xy = r0^2/2  (hyperbola)
        sext   (n=3):  Phi = Im(z^3)*c    -> pole  Im(z^3)=const
     and the orthogonal flux net A_z = Re(W) (quad: A_z = G(x^2-y^2)/2).

The ideal infinite pole gives a PERFECT multipole (no spurious harmonics);
the design value is the FINITE pole + finite-mu + saturation, measured by
(1) on the forward solve -- the next rungs ((a)-2 HDiv-MMM, (a)-3 the
saturation loop).

Verified here (analytic): the analyzer returns ONLY b_2 for a pure quad and
detects an injected octupole b_4 to machine precision.
"""
import math
import os

import numpy as np


# ============================================================
# multipole harmonic analyzer (the field-quality measurement)
# ============================================================

def multipoles(B_func, r_ref=1.0, n_max=12, n_samples=512):
    """Normal b_n / skew a_n multipoles of a 2-D field on a circle r_ref.

    B_func(x, y) -> (Bx, By).  Convention:
        By + i Bx = sum_{n>=1} (b_n + i a_n) (z/1)^{n-1},  z = x + i y
    so b_n + i a_n is the (n-1)-th Fourier coefficient of (By + i Bx) on the
    circle of radius r_ref.  Returns {n: (b_n, a_n)} for n = 1..n_max.
    """
    th = np.linspace(0.0, 2.0 * np.pi, n_samples, endpoint=False)
    f = np.empty(n_samples, dtype=complex)
    for k, t in enumerate(th):
        bx, by = B_func(r_ref * math.cos(t), r_ref * math.sin(t))
        f[k] = by + 1j * bx                       # By + i Bx
    C = np.fft.fft(f) / n_samples                 # C[k] = b_{k+1} + i a_{k+1}
    return {n: (float(C[n - 1].real), float(C[n - 1].imag))
            for n in range(1, n_max + 1)}


def relative_harmonics(mp, main):
    """|b_n + i a_n| / |main coefficient|, the usual 'units' content table."""
    m = abs(complex(*mp[main]))
    return {n: abs(complex(*mp[n])) / m for n in mp}


# ============================================================
# hodograph pole geometry (the design output)
# ============================================================

def quad_pole_hyperbola(G=1.0, r0=1.0, n_pts=80, span=2.5):
    """The ideal quadrupole pole face = the Phi = G r0^2/2 equipotential
    (hyperbola xy = r0^2/2).  Returns the upper-right pole curve points
    (the other three follow by symmetry); r0 is the aperture (inscribed)
    radius, pole tip on the diagonal at (r0/sqrt2, r0/sqrt2)."""
    c = r0 ** 2 / 2.0
    xs = np.linspace(c / span, span, n_pts)
    ys = c / xs
    return np.column_stack([xs, ys]), c


# ============================================================
# verification + figure
# ============================================================

def solve(G=1.0, r0=1.0, r_ref=0.7, eps_oct=0.05, plot=False):
    """Verify the analyzer on analytic fields and return the dict."""
    # pure quadrupole  B = G (y, x)  ->  only b_2
    B_quad = lambda x, y: (G * y, G * x)            # noqa: E731
    mq = multipoles(B_quad, r_ref)
    main = abs(complex(*mq[2]))
    spurious = max(abs(complex(*mq[n])) / main for n in mq if n != 2)

    # quad + injected octupole  By + i Bx = G z + eps z^3  ->  detect b_4
    def B_quad_oct(x, y):
        z = complex(x, y)
        w = G * z + eps_oct * z ** 3
        return (w.imag, w.real)                     # (Bx, By)
    mqo = multipoles(B_quad_oct, r_ref)
    b4 = abs(complex(*mqo[4]))
    b4_expected = eps_oct * r_ref ** 3              # |b_4| = eps r_ref^3
    b4_rel_err = abs(b4 - b4_expected) / b4_expected

    _, c = quad_pole_hyperbola(G, r0)

    if plot:
        _plot_quad_net(G, r0)

    return {
        "quad_main_b2": float(main),
        "quad_spurious_rel": float(spurious),       # ~ machine zero
        "qo_b4_detected": float(b4),
        "qo_b4_expected": float(b4_expected),
        "qo_b4_rel_err": float(b4_rel_err),         # ~ machine zero
        "pole_hyperbola_xy_const": float(c),        # pole: xy = r0^2/2
    }


def _plot_quad_net(G, r0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = np.linspace(-2.0, 2.0, 400)
    XX, YY = np.meshgrid(g, g)
    PHI = G * XX * YY                       # pole equipotentials
    AZ = 0.5 * G * (XX ** 2 - YY ** 2)      # orthogonal flux lines
    c = r0 ** 2 / 2.0

    fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=150)
    ax.contour(XX, YY, PHI, levels=np.linspace(-3, 3, 25),
               colors="C3", linewidths=0.6)
    ax.contour(XX, YY, AZ, levels=np.linspace(-3, 3, 25),
               colors="C0", linewidths=0.6)
    # the 4 pole faces  xy = +-c  (the equipotentials bounding the aperture)
    ax.contour(XX, YY, XX * YY, levels=[-c, c], colors="k", linewidths=2.0)
    th = np.linspace(0, 2 * math.pi, 200)
    ax.plot(r0 * np.cos(th), r0 * np.sin(th), "g--", lw=1.2,
            label=f"aperture r0={r0}")
    ax.set_aspect("equal")
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Quadrupole hodograph net: pole $\\Phi=Gxy$ (red, black=poles)\n"
                 "$\\perp$ flux $A_z=G(x^2-y^2)/2$ (blue)")
    ax.legend(loc="upper right", fontsize=8)
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    r = solve(plot=True)
    print("Accelerator pole design -- multipole analyzer + hodograph geometry\n")
    print("  multipole analyzer verification:")
    print(f"    pure quad     -> main b_2 = {r['quad_main_b2']:.4f}, "
          f"max spurious |b_n|/|b_2| = {r['quad_spurious_rel']:.1e}  (->0)")
    print(f"    quad+octupole -> b_4 detected = {r['qo_b4_detected']:.5f}  "
          f"(expected {r['qo_b4_expected']:.5f}, rel err {r['qo_b4_rel_err']:.1e})")
    print(f"\n  quadrupole pole face (hodograph): hyperbola  xy = "
          f"{r['pole_hyperbola_xy_const']:.4f}  (= r0^2/2)")
    print("  -> ideal pole = perfect quad; finite pole + finite-mu + "
          "saturation harmonics come from the forward solve (next rungs).")


if __name__ == "__main__":
    main()
