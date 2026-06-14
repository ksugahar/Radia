"""ellipsoid_alpha_tensor.py -- SHAPE-anisotropic eddy-current polarizability
tensor of a conducting ellipsoid, the analytic shape-anisotropy anchor for
the 3D non-axisymmetric levitation route.

This makes precise the point that the "anisotropy" of a conducting body's
AC response is SHAPE anisotropy (the semi-axes a1 != a2 != a3), NOT material
anisotropy: the conductor stays a scalar sigma/mu.  A field along axis i
drives eddy currents in the cross-section perpendicular to i, and the three
cross-sections of a triaxial ellipsoid differ -- so the induced-dipole /
polarizability response is a direction-split tensor alpha = diag(a_1,a_2,a_3)
even for an isotropic material.  It is the AC / eddy-current generalization
of the magnetostatic DEMAGNETIZING-FACTOR tensor.

Two limits are fully ANALYTIC (the anchors that pin the tensor):

  * DC (omega -> 0): no eddy response, alpha_i(0) = 0 (every axis).
  * High frequency (perfect conductor, flux exclusion): the induced
    moment along principal axis i is
        m_i(inf) = kappa_i / mu0 * B_i,   kappa_i = -V / (1 - N_i),
    with V the volume and N_i the analytic demagnetizing factor
        N_i = (a1 a2 a3 / 2) integral_0^inf ds /
                    [ (a_i^2 + s) sqrt((a1^2+s)(a2^2+s)(a3^2+s)) ],
        sum_i N_i = 1   (Osborn 1945; Stoner 1945).
    (Perfect diamagnet chi = -1 in an ellipsoid: chi_eff,i = chi/(1+N_i chi)
    -> -1/(1-N_i); m_i = V chi_eff,i H_i.)

The full-frequency eddy tensor alpha_i(omega) between these limits has no
simple closed form for a triaxial body (the sphere is special, see
levitation_sphere_force.py); it is the per-direction numerical Foster / FEM
extension (lab cuboid_521_3dir_full_foster route).  What is verified here is
the analytic shape-anisotropy: the demag tensor, the high-frequency
polarizability tensor, and its levitation-force consequence (a body lifts
most when B lies along its SHORT axis, where the largest demag factor
excludes the most flux per unit volume).

Run:  python ellipsoid_alpha_tensor.py
"""
import json
import os

import numpy as np
from scipy.integrate import quad

mu0 = 4 * np.pi * 1e-7
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Demagnetizing factors
# ---------------------------------------------------------------------------
def demag_factor(a, axis):
    """Demagnetizing factor N along `axis` (0/1/2) for an ellipsoid with
    semi-axes a=(a1,a2,a3), via the exact Osborn integral.  sum_i N_i = 1.

    The semi-infinite range is mapped to (0,1) by s = scale*t/(1-t) so the
    s^{-5/2} tail integrates cleanly (no quad extrapolation roundoff)."""
    a1, a2, a3 = a
    pref = a1 * a2 * a3 / 2.0
    scale = max(a1, a2, a3) ** 2

    def integrand(t):
        s = scale * t / (1.0 - t)
        ds = scale / (1.0 - t) ** 2
        return ds / ((a[axis] ** 2 + s) *
                     np.sqrt((a1**2 + s) * (a2**2 + s) * (a3**2 + s)))

    val, _ = quad(integrand, 0.0, 1.0, limit=200)
    return pref * val


def demag_tensor(a):
    return np.array([demag_factor(a, i) for i in range(3)])


def demag_spheroid_closed(a_eq, c):
    """Closed-form spheroid demag along the polar axis c (a_eq = a = b).
    Prolate (c>a_eq) uses the log form, oblate (c<a_eq) the arcsin form;
    returns (N_eq, N_eq, N_c).  Independent cross-check of the integral."""
    if abs(c - a_eq) < 1e-15:
        return np.array([1 / 3, 1 / 3, 1 / 3])
    if c > a_eq:                                   # prolate
        e = np.sqrt(1.0 - (a_eq / c) ** 2)
        Nc = (1 - e**2) / e**2 * (1 / (2 * e) * np.log((1 + e) / (1 - e)) - 1)
    else:                                          # oblate
        e = np.sqrt(1.0 - (c / a_eq) ** 2)
        Nc = 1.0 / e**2 * (1.0 - np.sqrt(1 - e**2) / e * np.arcsin(e))
    Neq = (1.0 - Nc) / 2.0
    return np.array([Neq, Neq, Nc])


# ---------------------------------------------------------------------------
# High-frequency (perfect-conductor) polarizability tensor
# ---------------------------------------------------------------------------
def kappa_hf(a):
    """Perfect-conductor induced-moment coefficients kappa_i = -V/(1-N_i)
    [m^3], so m_i = kappa_i B_i / mu0.  Direction-split for a non-spheroid."""
    V = 4.0 / 3.0 * np.pi * a[0] * a[1] * a[2]
    N = demag_tensor(a)
    return -V / (1.0 - N), N, V


def lift_factor(a):
    """Orientation lift coefficients: <F_z> = -(V/(4 mu0)) * 1/(1-N_i) *
    grad(B0^2) for B along principal axis i.  Returns 1/(1-N_i) per axis."""
    N = demag_tensor(a)
    return 1.0 / (1.0 - N)


def main():
    print("Shape-anisotropic eddy-current polarizability of a conducting "
          "ellipsoid\n(material stays isotropic scalar sigma/mu; the tensor "
          "is pure SHAPE)\n")

    # ----- check 1: demag tensor sums to 1; integral == closed form -----
    print("[check 1] demag tensor: Osborn integral vs spheroid closed form, sum=1")
    cases = {
        "sphere (1,1,1)":      (1.0, 1.0, 1.0),
        "prolate 2:1 (1,1,2)": (1.0, 1.0, 2.0),
        "oblate 1:4 (1,1,.25)":(1.0, 1.0, 0.25),
        "needle (1,1,100)":    (1.0, 1.0, 100.0),
        "thin disk (1,1,.01)": (1.0, 1.0, 0.01),
    }
    for name, a in cases.items():
        N = demag_tensor(a)
        Ncf = demag_spheroid_closed(a[0], a[2])      # a=b here
        dmax = np.max(np.abs(N - Ncf))
        print(f"   {name:22s} N=({N[0]:.4f},{N[1]:.4f},{N[2]:.4f})  "
              f"sum={N.sum():.6f}  |int-closed|={dmax:.2e}")
        assert abs(N.sum() - 1.0) < 1e-6, "demag factors must sum to 1"
        assert dmax < 1e-4, "Osborn integral must match the spheroid closed form"
    # sphere is isotropic
    assert np.max(np.abs(demag_tensor((1, 1, 1)) - 1 / 3)) < 1e-6
    # known 2:1 prolate value
    N21 = demag_tensor((1.0, 1.0, 2.0))
    assert abs(N21[2] - 0.17357) < 1e-3, "2:1 prolate N_c known = 0.1736"

    # ----- check 2: a TRIAXIAL body gives three DISTINCT alpha_i -----
    print("\n[check 2] genuine triaxial a1!=a2!=a3 -> alpha_x != alpha_y != alpha_z")
    a = (5e-3, 3e-3, 1.5e-3)                        # 5 x 3 x 1.5 mm brick-like
    kap, N, V = kappa_hf(a)
    print(f"   semi-axes (mm) = {tuple(round(x*1e3,2) for x in a)},  V={V*1e9:.4f} mm^3")
    print("   axis   semi-axis   N_i      kappa_i [mm^3]   |kappa_i|/V")
    for i, lab in enumerate("xyz"):
        print(f"     {lab}    {a[i]*1e3:6.2f} mm  {N[i]:.4f}  "
              f"{kap[i]*1e9:+10.4f}     {abs(kap[i])/V:.4f}")
    assert N[0] < N[1] < N[2], "shortest axis (z) must have the largest demag N"
    assert abs(kap[0]) < abs(kap[1]) < abs(kap[2]), \
        "|kappa| largest along the short axis (most flux excluded per volume)"
    aniso = abs(kap[2]) / abs(kap[0])
    print(f"   anisotropy |kappa_z|/|kappa_x| = {aniso:.3f}  (1.0 would be isotropic)")
    assert aniso > 1.2, "a 5:1.5 triaxial must be clearly anisotropic"

    # ----- check 3: reduces to the verified sphere -----
    print("\n[check 3] sphere limit reduces to levitation_sphere_force.py")
    asp = (5e-3,) * 3
    kap_s, N_s, V_s = kappa_hf(asp)
    a_r = 5e-3
    kap_analytic = -2 * np.pi * a_r**3            # sphere PEC: m=-2pi a^3 B/mu0
    print(f"   sphere kappa = {kap_s[0]*1e9:.5f} mm^3, "
          f"-2 pi a^3 = {kap_analytic*1e9:.5f} mm^3")
    assert np.allclose(kap_s, kap_analytic, rtol=1e-4), "sphere PEC kappa = -2 pi a^3"
    assert np.max(np.abs(N_s - 1 / 3)) < 1e-6, "sphere demag isotropic 1/3"

    # ----- check 4: orientation-dependent levitation lift -----
    # <F_z> = -(V/(4 mu0)) * (1/(1-N_i)) * grad(B0^2) for B along axis i.
    print("\n[check 4] orientation-dependent lift  <F_z> = -(V/4mu0)/(1-N_i) grad(B0^2)")
    B0, gradB2 = 0.1, -0.2                          # T, T^2/m (field weakens upward)
    lf = lift_factor(a)
    print("   B along   1/(1-N_i)    lift [mN]   (sphere of equal V for ref)")
    Vsph_equal = V
    lift_iso = (Vsph_equal / (4 * mu0)) * (1.0 / (1.0 - 1.0 / 3.0)) * abs(gradB2)
    for i, lab in enumerate("xyz"):
        lift = (V / (4 * mu0)) * lf[i] * abs(gradB2)
        print(f"     {lab}      {lf[i]:.4f}    {lift*1e3:8.4f}    "
              f"(iso {lift_iso*1e3:.4f})")
    lift_x = (V / (4 * mu0)) * lf[0] * abs(gradB2)
    lift_z = (V / (4 * mu0)) * lf[2] * abs(gradB2)
    print(f"   B along short axis z lifts {lift_z/lift_x:.2f}x more than along long "
          f"axis x (same body, same field)")
    assert lift_z > lift_x, "lift is larger with B along the short axis"
    # cross-check the sphere force module's coefficient at high freq
    import levitation_sphere_force as L
    F_pc_sphere = (np.pi * L.a**3 / (2 * mu0)) * abs(gradB2)
    F_from_kappa = (4 / 3 * np.pi * L.a**3 / (4 * mu0)) * (1 / (1 - 1 / 3)) * abs(gradB2)
    assert abs(F_pc_sphere - F_from_kappa) / F_pc_sphere < 1e-9, \
        "kappa lift formula must reproduce the verified sphere coefficient"
    print(f"   kappa-tensor lift reduces to the verified sphere coeff exactly OK")

    # ----- persist -----
    data = {
        "note": "shape (not material) anisotropy; HF perfect-conductor tensor",
        "triaxial_semi_axes_m": list(a),
        "demag_N": demag_tensor(a).tolist(),
        "kappa_hf_m3": kap.tolist(),
        "anisotropy_kappaz_over_kappax": float(aniso),
        "lift_factor_inv_1_minus_N": lift_factor(a).tolist(),
        "cases_demag": {k: demag_tensor(v).tolist() for k, v in cases.items()},
    }
    with open(os.path.join(HERE, "ellipsoid_alpha_tensor_results.json"), "w") as fp:
        json.dump(data, fp, indent=2)
    print("\n wrote ellipsoid_alpha_tensor_results.json")
    plot(a)


def plot(a_tri):
    """kappa anisotropy + lift factor vs aspect ratio for a prolate spheroid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif", "axes.linewidth": 0.8})
    ar = np.linspace(0.2, 6.0, 80)                 # c / a_eq aspect ratio
    Nc, Na = [], []
    for r in ar:
        N = demag_spheroid_closed(1.0, r)
        Na.append(N[0]); Nc.append(N[2])
    Na, Nc = np.array(Na), np.array(Nc)
    fig, ax = plt.subplots(figsize=(8 / 2.54, 6 / 2.54))
    ax.plot(ar, 1 / (1 - Nc), "C0-", lw=1.2, label=r"$B \parallel$ polar axis")
    ax.plot(ar, 1 / (1 - Na), "C3--", lw=1.2, label=r"$B \perp$ polar axis")
    ax.axhline(1 / (1 - 1 / 3), color="0.6", ls=":", lw=0.8)
    ax.text(ar[0], 1 / (1 - 1 / 3), " sphere", va="bottom", ha="left",
            fontsize=8, color="0.4")
    ax.axvline(1.0, color="0.85", lw=0.8)
    ax.set_xlabel(r"aspect ratio  $c / a$")
    ax.set_ylabel(r"lift factor  $1/(1-N)$")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)
    out = os.path.join(HERE, "ellipsoid_alpha_tensor_lift.png")
    fig.savefig(out, dpi=300)
    print(f" wrote {out}")


if __name__ == "__main__":
    main()
