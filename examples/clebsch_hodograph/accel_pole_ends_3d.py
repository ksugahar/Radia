"""(A) The magnet ENDS in 3-D: the integrated-multipole analyzer + the
equipotential-following end rule.

RESEARCH example (track A, the July main goal -- foundation rung).  Builds the
measurement tool track A needs (the INTEGRATED multipole analyzer) and proves,
analytically and to machine precision, the central design claim of
docs/clebsch_hodograph/DESIGN_METHODOLOGY.md sec 3.2:

  > The 3-D equipotential surface of the reduced-potential solution tells you
  > the optimal end-iron contour.  Shape the end iron to FOLLOW the
  > equipotential => the INTEGRATED-field harmonics are minimized.

The theorem (why it is true, and verifiable with no FEM):

1. The INTEGRATED transverse field  Bbar_perp(x,y) = INT B_perp(x,y,z) dz  is
   ALWAYS a 2-D harmonic (multipole) field.  Proof: in the current-free
   aperture div B = 0 and curl B = 0, so
       d Bbar_x/dx + d Bbar_y/dy = INT (dBx/dx + dBy/dy) dz
                                 = INT (-dBz/dz) dz = -[Bz] = 0,
       d Bbar_y/dx - d Bbar_x/dy = INT (curl B)_z dz = 0,
   i.e. Bbar_perp is 2-D div-free AND curl-free => a clean 2-D multipole.
   So the INTEGRATED multipole analyzer is well defined (it is what beam optics
   actually sees: INT B dl along the trajectory).

2. A MAXWELLIAN quad end (the 3-D-harmonic field whose gradient G(z) turns
   on/off while the cross-section keeps its m=2 quad symmetry -- i.e. the iron
   end FOLLOWS the 3-D equipotential) has TWO honest properties:
   (a) its on-axis-gradient Maxwell RADIAL corrections (the pseudo-multipoles,
       proportional to G''(z), G''''(z), ... times r^3, r^5 in the SAME m=2
       azimuthal channel) are TOTAL z-DERIVATIVES -> they integrate to ZERO
       (INT G'' dz = [G'] = 0 for a localized magnet).  So the integrated quad
       is EXACT and radially undistorted: bbar_2 = (INT G) r_ref, with no
       integrated radial pseudo-multipole.
   (b) because the symmetric end preserves the m=2 azimuthal symmetry, it
       generates NO azimuthal b_6 / b_10 at all (bbar_6 = bbar_10 = 0).

3. A NON-equipotential end (the iron deviating from the 3-D equipotential)
   BREAKS the m=2 azimuthal symmetry: it injects a GENUINE, one-signed b_6
   localized at the end whose z-integral is NONZERO => a spurious INTEGRATED
   bbar_6, growing ~linearly with the deviation.

So "follow the 3-D equipotential at the end" = keep the m=2 symmetry AND let
the gradient's radial Maxwell corrections cancel in the integral, so the
integrated field stays the pure DESIGNED multipole.  This example demonstrates
1-3 exactly: (2a) bbar_2 exact to ~1e-6, (2b)+(3) bbar_6 ~ 0 for the symmetric
end and ~ linear in the deviation for the broken end.

The 2-D cross-section (the body) is accel_pole_design.py / accel_pole_harmonics.py
(the equipotential = the hyperbola; deviation = harmonics).  THIS file is the
3-D END extension: the same lever, now on the INTEGRATED field.
"""
import math
import os

import numpy as np

ALLOWED = (2, 6, 10)            # quad-allowed normal harmonics (n = 2 mod 4)

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy 2.0 renamed trapz


# ============================================================
# the integrated multipole analyzer (the beam-optics measurement)
# ============================================================

def _field_harmonics(f_circle, n_max):
    """FFT a complex sample array (By + i Bx) on a circle into {n: C[n-1]}.
    C[n-1] = b_n + i a_n at the sampling radius (accelerator convention
    By + i Bx = sum_n (b_n + i a_n) z^{n-1})."""
    C = np.fft.fft(f_circle) / len(f_circle)
    return {n: complex(C[n - 1]) for n in range(1, n_max + 1)}


def integrated_multipoles(B_perp, r_ref=0.02, z_range=(-1.0, 1.0),
                          n_z=401, n_theta=256, n_max=12):
    """Normal/skew multipoles of the INTEGRATED transverse field
    Bbar_perp(x,y) = INT B_perp(x,y,z) dz  on a circle of radius r_ref.

    B_perp(x, y, z) -> (Bx, By).  Returns {n: (bbar_n, abar_n)} (field
    harmonics at r_ref, integrated over z) for n = 1..n_max.  Because the
    integrated field is exactly a 2-D multipole (see module docstring), this
    is well defined regardless of the per-slice fringe structure.
    """
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zs = np.linspace(z_range[0], z_range[1], n_z)
    xs = r_ref * np.cos(th)
    ys = r_ref * np.sin(th)

    fbar = np.zeros(n_theta, dtype=complex)            # INT (By + i Bx) dz
    for k in range(n_theta):
        bx = np.array([B_perp(xs[k], ys[k], z)[0] for z in zs])
        by = np.array([B_perp(xs[k], ys[k], z)[1] for z in zs])
        fbar[k] = _trapz(by, zs) + 1j * _trapz(bx, zs)

    H = _field_harmonics(fbar, n_max)
    return {n: (float(H[n].real), float(H[n].imag)) for n in H}


def per_slice_harmonic(B_perp, n, r_ref, zs, n_theta=256):
    """The n-th field harmonic |b_n + i a_n|(z) on each slice z (for the
    'pseudo-multipole oscillates but integrates to zero' picture)."""
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    xs = r_ref * np.cos(th)
    ys = r_ref * np.sin(th)
    out = np.empty(len(zs), dtype=complex)
    for j, z in enumerate(zs):
        f = np.array([B_perp(xs[k], ys[k], z)[1] + 1j * B_perp(xs[k], ys[k], z)[0]
                      for k in range(n_theta)])
        out[j] = (np.fft.fft(f) / n_theta)[n - 1]
    return out


# ============================================================
# the 3-D fields: Maxwellian (good) end vs a non-equipotential (bad) end
# ============================================================

def gaussian_gradient(G0=10.0, sigma=0.2):
    """A localized quad gradient profile G(z) = G0 exp(-(z/sigma)^2 / 2) and
    its 2nd derivative G''(z) (closed form).  INT G dz = G0 sigma sqrt(2 pi);
    INT G'' dz = 0 exactly (localized) -- the key to the theorem."""
    def G(z):
        return G0 * math.exp(-0.5 * (z / sigma) ** 2)

    def G2(z):                                          # G''(z)
        return G0 * math.exp(-0.5 * (z / sigma) ** 2) * (z * z - sigma * sigma) / sigma ** 4

    return G, G2, G0 * sigma * math.sqrt(2.0 * math.pi)


def maxwellian_quad(G, G2):
    """The 3-D HARMONIC quad whose gradient turns on/off as G(z) -- i.e. the
    iron end FOLLOWS the 3-D equipotential.  Keeps the leading n=2 term and the
    first Maxwell correction (proportional to G''(z), a total z-derivative):

        Phi = -G(z) xy + (1/12) G''(z) (x^3 y + x y^3)
        Bx  = -dPhi/dx = G y  - (1/12) G'' (3 x^2 y + y^3)
        By  = -dPhi/dy = G x  - (1/12) G'' (x^3 + 3 x y^2)

    The correction is O((r/sigma)^2) on the circle and O(G'') in z, so its
    z-INTEGRAL vanishes -> the integrated field is a pure quad.
    """
    def B_perp(x, y, z):
        g, g2 = G(z), G2(z)
        bx = g * y - (g2 / 12.0) * (3.0 * x * x * y + y ** 3)
        by = g * x - (g2 / 12.0) * (x ** 3 + 3.0 * x * y * y)
        return (bx, by)

    return B_perp


def with_bad_end(B_perp, c6=0.0, z_end=0.4, d_end=0.05):
    """Add a NON-equipotential end defect: a GENUINE, one-signed normal
    12-pole (n = 6) localized at the +z end (the iron deviating from the 3-D
    equipotential).  delta(By + i Bx) = c6 * w(z) * (x + i y)^5 with w(z) >= 0,
    so INT delta b_6 dz = c6 (INT w dz) (field) != 0 -- a spurious INTEGRATED
    harmonic that does NOT cancel (unlike the Maxwellian fringe)."""
    def B_bad(x, y, z):
        bx, by = B_perp(x, y, z)
        w = math.exp(-((z - z_end) / d_end) ** 2)       # one-signed bump at the end
        zc5 = (complex(x, y)) ** 5                       # genuine n=6 field shape
        dby = c6 * w * zc5.real
        dbx = c6 * w * zc5.imag
        return (bx + dbx, by + dby)

    return B_bad


# ============================================================
# verification + figure
# ============================================================

def solve(G0=10.0, sigma=0.2, r_ref=0.02, z_hi=1.0, n_z=801, plot=False):
    """Verify the theorem: Maxwellian end -> pure integrated quad; a
    non-equipotential end defect -> integrated b_6 growing ~linearly."""
    G, G2, intG = gaussian_gradient(G0, sigma)
    B_good = maxwellian_quad(G, G2)
    zr = (-z_hi, z_hi)

    # (1) Maxwellian (equipotential-following) end:
    #   (2a) the radial Maxwell corrections (~G'') integrate away -> bbar_2 exact;
    #   (2b) the symmetric end makes no azimuthal b_6 -> integrated bbar_6 ~ 0.
    mg = integrated_multipoles(B_good, r_ref, zr, n_z)
    main = abs(complex(*mg[2]))
    b2_expected = intG * r_ref                           # bbar_2 field at r_ref = (INT G) r_ref
    b2_rel_err = abs(main - b2_expected) / b2_expected
    spurious_good = max(abs(complex(*mg[n])) / main for n in ALLOWED if n != 2)

    # (2) non-equipotential end defect: spurious integrated b_6 grows ~linearly
    c6s = [0.0, 2.0e3, 4.0e3, 8.0e3]
    b6_rel = []
    for c6 in c6s:
        mb = integrated_multipoles(with_bad_end(B_good, c6=c6), r_ref, zr, n_z)
        b6_rel.append(abs(complex(*mb[6])) / abs(complex(*mb[2])))
    slopes = [b6_rel[i] / c6s[i] for i in range(1, len(c6s))]
    slope_spread = (max(slopes) - min(slopes)) / float(np.mean(slopes))

    if plot:
        _plot_ends(B_good, G0, sigma, r_ref, z_hi)

    return {
        "intG_expected": float(intG),
        "good_main_b2": float(main),
        "good_b2_expected": float(b2_expected),
        "good_b2_rel_err": float(b2_rel_err),            # ~0: pure integrated quad strength
        "good_spurious_b6_rel": float(spurious_good),    # ~0: Maxwellian fringe integrates away
        "bad_b6_rel_at_c6": {f"{c:.0f}": float(v) for c, v in zip(c6s, b6_rel)},
        "bad_b6_slope_spread": float(slope_spread),      # ~0: linear in the end deviation
    }


def _plot_ends(B_good, G0, sigma, r_ref, z_hi):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    G, _G2, _ = gaussian_gradient(G0, sigma)
    zs = np.linspace(-z_hi, z_hi, 400)

    # LEFT: the genuine Maxwell pseudo-multipole = the per-slice b_2 RADIAL
    # correction (b_2(z) - G(z) r_ref ~ G''(z) r_ref^3).  It oscillates in z and
    # INTEGRATES TO ZERO -> the integrated quad strength stays exact.
    b2_good = per_slice_harmonic(B_good, 2, r_ref, zs).real
    db2 = b2_good - np.array([G(z) for z in zs]) * r_ref      # radial correction ~ G''

    # RIGHT: the per-slice b_6.  Symmetric (equipotential) end -> ~0 everywhere;
    # a non-equipotential end defect -> a one-signed bump that does NOT cancel.
    b6_good = per_slice_harmonic(B_good, 6, r_ref, zs).real
    b6_bad = per_slice_harmonic(with_bad_end(B_good, c6=4.0e3), 6, r_ref, zs).real

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6), dpi=150)
    ax[0].plot(zs, db2, "C0")
    ax[0].fill_between(zs, db2, 0, alpha=0.25, color="C0")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_title("Equipotential end: radial Maxwell correction "
                    "$\\propto G''(z)$\noscillates, $\\int\\,dz=0$ "
                    "$\\Rightarrow$ integrated quad exact")
    ax[0].set_ylabel("$\\bar b_2$ radial correction at $r_{ref}$")
    ax[1].plot(zs, b6_good, "C0", label="equipotential end ($b_6\\equiv0$)")
    ax[1].plot(zs, b6_bad, "C3", label="non-equipotential defect")
    ax[1].fill_between(zs, b6_bad, 0, alpha=0.25, color="C3")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_title("$b_6(z)$: a non-equipotential end breaks symmetry\n"
                    "$\\Rightarrow$ one-signed, $\\int b_6\\,dz\\neq0$")
    ax[1].set_ylabel("$b_6(z)$ field at $r_{ref}$")
    ax[1].legend(loc="upper left", fontsize=8)
    for a in ax:
        a.set_xlabel("z [m]")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    r = solve(plot=True)
    print("(A) Magnet ENDS in 3-D -- integrated multipole analyzer + "
          "equipotential-end rule\n")
    print("  (1) MAXWELLIAN end (iron follows the 3-D equipotential):")
    print(f"      integrated quad  bbar_2 = {r['good_main_b2']:.5f}  "
          f"(expected (INT G) r_ref = {r['good_b2_expected']:.5f}, "
          f"rel err {r['good_b2_rel_err']:.1e})")
    print(f"        -> radial Maxwell corrections (~G'') integrate away: "
          f"integrated quad is exact")
    print(f"      spurious integrated |bbar_6/bbar_2| = "
          f"{r['good_spurious_b6_rel']:.2e}   (-> 0: symmetric end makes "
          f"no azimuthal b_6)\n")
    print("  (2) NON-equipotential end defect -> spurious integrated bbar_6/bbar_2:")
    for c, v in r["bad_b6_rel_at_c6"].items():
        print(f"      end-deviation c6={c:>6}:  bbar_6/bbar_2 = {v:.4e}")
    print(f"      -> grows ~linearly with the end deviation "
          f"(slope spread {r['bad_b6_slope_spread']:.1%})\n")
    print("  => 'follow the 3-D equipotential at the end' = the only end harmonic")
    print("     that does NOT integrate to zero is the desired one (track A core).")


if __name__ == "__main__":
    main()
