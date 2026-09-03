"""coil_maglev_equilibrium.py -- AC magnetic levitation equilibrium of a
small conducting sphere above a coil.

This composes two already-verified pieces:

  (1) Radia (open-boundary IEM field).  A circular coil on the z=0 plane
      (axis = z) carries an AC current of PEAK amplitude I0 at angular
      frequency omega.  Radia computes the STATIC spatial field of that
      coil; because the geometry is fixed, the AC field amplitude is just
      B0(r) = (Radia field at current I0) [Tesla].  No air mesh, exact open
      boundary -- the "Radia + NGSolve maglev" policy reduced to the part
      that needs Radia (here the conductor reaction is analytic, so NGSolve
      is not invoked: the sphere's eddy response is the closed-form G).

  (2) The verified induced-dipole levitation force
      (maglev_sphere_force.py, coefficient pinned to the perfect-
      conductor energy limit to 0.09%):

          <F_z>(z) = (pi a^3 / mu0) * Re[G(omega)] * d/dz( |B0(z)|^2 )

      with B0(z) the PEAK spatial field amplitude.  Re[G] < 0 (a conducting
      sphere is an AC diamagnet) means the force points toward weaker field,
      i.e. UP, away from the coil = LIFT.

Equilibrium: a Cu sphere ABOVE the coil sits where the field weakens with
height (dB2/dz < 0); the lift balances gravity at the height z* where

          F_z(z*) = weight = rho_Cu * (4/3 pi a^3) * g.

VERTICAL STABILITY requires dF_z/dz < 0 at z* (a small upward displacement
must reduce the lift so the sphere falls back).  This selects the UPPER
branch of F_z(z): below the lift peak the force grows with height
(unstable); above it the force decays with height (stable).

Everything quantitative is backed by a hard assert (fail loud, per the
"No Fallbacks -- Fail Fast" policy).

Run:  python coil_maglev_equilibrium.py     (a few seconds)
"""
import json
import os
import sys

import numpy as np

import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from maglev_sphere_force import G_exact, delta, a, sigma  # noqa: E402
from _validation_output import validation_output  # noqa: E402

mu0 = 4 * np.pi * 1e-7
rho_Cu = 8960.0          # copper mass density [kg/m^3]
g = 9.81                 # gravity [m/s^2]

# ----- coil design (chosen so an equilibrium EXISTS with comfortable margin) -----
R_COIL = 0.03            # coil radius [m] (30 mm)
N_TURN = 100             # number of turns
I_PEAK = 100.0           # peak current per turn [A]  -> N*I = 10000 At
NI = N_TURN * I_PEAK     # total amp-turns (peak)
NSEG = 360               # filament segments approximating the circle
F_HZ = 50e3              # AC frequency [Hz]
OMEGA = 2 * np.pi * F_HZ

# sphere weight (the load to levitate)
SPHERE_VOL = 4.0 / 3.0 * np.pi * a**3
WEIGHT = rho_Cu * SPHERE_VOL * g


def Bz_loop_analytic(z, R, NI_amp):
    """On-axis field of an ideal circular current loop (radius R, NI amp-turns):
       Bz(z) = mu0 NI R^2 / (2 (R^2 + z^2)^{3/2})."""
    return mu0 * NI_amp * R**2 / (2.0 * (R**2 + z**2) ** 1.5)


def build_coil():
    """Thin circular filament loop in the z=0 plane, axis = z, carrying the
    total peak amp-turns NI.  A single multi-turn loop is the exact analogue
    of the on-axis analytic formula (turns are co-located, thin coil)."""
    phis = np.linspace(0.0, 2.0 * np.pi, NSEG + 1)
    pts = [[R_COIL * np.cos(p), R_COIL * np.sin(p), 0.0] for p in phis]
    return rad.ObjFlmCur(pts, NI)


def field_on_axis(coil, z):
    """Peak field amplitude vector B0(z) [Tesla] on the coil axis."""
    return np.asarray(rad.Fld(coil, "b", [0.0, 0.0, float(z)]), dtype=float)


def b2_on_axis(coil, zs):
    """|B0(z)|^2 [T^2] on the axis for an array of heights."""
    B = np.array([field_on_axis(coil, z) for z in zs])
    return np.sum(B * B, axis=1)


def force_on_axis(coil, zs):
    """<F_z>(z) [N] on the axis via central-difference d|B0|^2/dz.
       Positive = upward (lift)."""
    ReG = G_exact(OMEGA).real
    b2 = b2_on_axis(coil, zs)
    dz = zs[1] - zs[0]
    dB2 = np.gradient(b2, dz)                 # central differences interior
    return (np.pi * a**3 / mu0) * ReG * dB2   # ReG<0 and dB2<0 above coil -> F>0


def main():
    print("=" * 70)
    print("AC magnetic levitation equilibrium: Cu sphere above a coil")
    print("=" * 70)
    print(f"sphere:  a = {a*1e3:.1f} mm, sigma = {sigma:.2e} S/m, "
          f"rho = {rho_Cu:.0f} kg/m^3")
    print(f"weight:  {WEIGHT*1e3:.4f} mN  ({WEIGHT:.4e} N)")
    print(f"coil:    R = {R_COIL*1e3:.1f} mm, N = {N_TURN} turns, "
          f"I0 = {I_PEAK:.1f} A (peak)  ->  NI = {NI:.0f} At")
    print(f"freq:    f = {F_HZ/1e3:.1f} kHz")

    rad.UtiDelAll()
    try:
        coil = build_coil()

        # ============================================================
        # CHECK 1: Radia on-axis field vs analytic single-loop formula
        # ============================================================
        print("\n[check 1] Radia on-axis field vs analytic loop "
              "Bz = mu0 NI R^2 / (2 (R^2+z^2)^3/2)")
        print("   z[mm]    Bz_radia[T]     Bz_analytic[T]    rel.err%")
        z_check = np.array([0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.20])
        max_err = 0.0
        for z in z_check:
            Bz = field_on_axis(coil, z)[2]
            Ba = Bz_loop_analytic(z, R_COIL, NI)
            err = abs(Bz - Ba) / abs(Ba) * 100.0
            max_err = max(max_err, err)
            print(f"  {z*1e3:6.1f}   {Bz: .6e}   {Ba: .6e}   {err:.4f}")
        print(f"   max field error = {max_err:.4f} %  (require < 1 %)")
        assert max_err < 1.0, f"Radia coil field deviates {max_err:.3f}% from analytic loop"

        # ============================================================
        # CHECK 2: sphere responds at this frequency (Re[G] < 0)
        # ============================================================
        d = delta(OMEGA)
        a_over_delta = a / d
        ReG = G_exact(OMEGA).real
        print(f"\n[check 2] sphere AC response at {F_HZ/1e3:.1f} kHz")
        print(f"   skin depth delta = {d*1e6:.2f} um,  a/delta = {a_over_delta:.2f}  "
              f"(require >~ 2)")
        print(f"   Re[G(omega)] = {ReG:+.5f}  (require < 0 = diamagnetic lift)")
        assert a_over_delta >= 2.0, f"a/delta = {a_over_delta:.2f} too small (weak eddy response)"
        assert ReG < 0.0, f"Re[G] = {ReG:.4f} not negative (no lift)"

        # ============================================================
        # CHECK 3: force profile, equilibrium root, stability
        # ============================================================
        # dense axial grid from just above the coil to several radii up
        zs = np.linspace(0.004, 0.25, 2000)
        b2 = b2_on_axis(coil, zs)
        Fz = force_on_axis(coil, zs)

        i_peak = int(np.argmax(Fz))
        z_peak = zs[i_peak]
        F_peak = Fz[i_peak]
        print(f"\n[check 3] force profile  F_z(z) = (pi a^3/mu0) Re[G] d|B0|^2/dz")
        print(f"   peak lift = {F_peak*1e3:.3f} mN at z = {z_peak*1e3:.2f} mm")
        print(f"   weight    = {WEIGHT*1e3:.3f} mN   (peak/weight = {F_peak/WEIGHT:.2f})")
        assert F_peak > WEIGHT, ("peak available lift < weight: no equilibrium exists "
                                 "(raise NI or shrink R)")

        # sanity: lift -> 0 far above, and lift grows toward the coil from the top
        assert abs(Fz[-1]) < 0.05 * F_peak, "lift should decay toward zero far above coil"
        # on the STABLE (upper) branch the force decreases with height
        assert Fz[i_peak + 1] < Fz[i_peak] if i_peak + 1 < len(Fz) else True

        # equilibrium on the UPPER (stable) branch: F crosses WEIGHT going down
        # search the descending part z > z_peak for F(z*) = WEIGHT
        upper = np.arange(i_peak, len(zs))
        Fu = Fz[upper]
        # find first index where force drops below weight
        below = np.where(Fu <= WEIGHT)[0]
        assert below.size > 0, "force never falls below weight above the peak (unbounded?)"
        k = upper[below[0]]
        # linear interpolation between zs[k-1] (F>weight) and zs[k] (F<=weight)
        z0, z1 = zs[k - 1], zs[k]
        F0, F1 = Fz[k - 1], Fz[k]
        z_star = z0 + (WEIGHT - F0) * (z1 - z0) / (F1 - F0)

        # evaluate the force AT z_star directly (independent re-check of the root)
        # use a local 3-point central difference about z_star with a fine step
        h = 2e-5
        b2p = np.sum(field_on_axis(coil, z_star + h) ** 2)
        b2m = np.sum(field_on_axis(coil, z_star - h) ** 2)
        dB2_star = (b2p - b2m) / (2 * h)
        F_star = (np.pi * a**3 / mu0) * ReG * dB2_star

        print(f"\n   equilibrium height z* = {z_star*1e3:.3f} mm (stable upper branch)")
        print(f"   F_z(z*) = {F_star*1e3:.4f} mN   weight = {WEIGHT*1e3:.4f} mN   "
              f"residual = {(F_star-WEIGHT)/WEIGHT*100:+.3f} %")
        assert abs(F_star - WEIGHT) / WEIGHT < 0.02, \
            f"force at z* off by {(F_star-WEIGHT)/WEIGHT*100:.2f}% from weight (root find failed)"

        # ============================================================
        # CHECK 4: vertical stability  dF_z/dz < 0  at z*
        # ============================================================
        hh = 1e-4
        def F_at(z):
            b2p_ = np.sum(field_on_axis(coil, z + h) ** 2)
            b2m_ = np.sum(field_on_axis(coil, z - h) ** 2)
            return (np.pi * a**3 / mu0) * ReG * (b2p_ - b2m_) / (2 * h)
        dFdz = (F_at(z_star + hh) - F_at(z_star - hh)) / (2 * hh)
        print(f"\n[check 4] vertical stability at z*")
        print(f"   dF_z/dz = {dFdz:+.4e} N/m   (require < 0 = restoring)")
        # restoring stiffness magnitude and small-oscillation frequency
        mass = rho_Cu * SPHERE_VOL
        if dFdz < 0:
            f_osc = np.sqrt(-dFdz / mass) / (2 * np.pi)
            print(f"   restoring stiffness |dF/dz| = {abs(dFdz):.4e} N/m, "
                  f"mass = {mass*1e3:.3f} g  ->  f_osc ~ {f_osc:.2f} Hz")
        assert dFdz < 0.0, f"dF/dz = {dFdz:.3e} >= 0 : equilibrium is vertically UNSTABLE"

        # ----- dipole-validity caveat (reported, not asserted) -----
        # local field-gradient scale L_grad = |B0|^2 / |d|B0|^2/dz| at z*
        b2_star = np.sum(field_on_axis(coil, z_star) ** 2)
        L_grad = b2_star / abs(dB2_star)
        print(f"\n   dipole validity: sphere a = {a*1e3:.1f} mm, "
              f"field-gradient scale L = {L_grad*1e3:.1f} mm, "
              f"z*/a = {z_star/a:.2f}, a/L = {a/L_grad:.3f}")

        # ============================================================
        # CHECK 5: the dipole approximation is CONTROLLABLE.  Both lift
        # and weight scale as a^3, so the equilibrium condition reduces
        # EXACTLY to  d|B0|^2/dz = (4/3) rho g mu0 / Re[G(a/delta)]  -- the
        # a^3 cancels, so z* depends on the sphere size ONLY through
        # Re[G(a/delta)].  Shrinking the sphere drives a/L -> 0 (clean
        # point dipole); z* then shifts only modestly, and that shift is
        # purely the eddy DESATURATION of Re[G] (a/delta -> -1/2 only for
        # a/delta >> 1).  So the a=5 mm a/L~0.5 caveat above is a
        # controllable choice, not a limitation of the method.
        # ============================================================
        import maglev_sphere_force as _Lm
        dB2dz_grid = np.gradient(b2, zs[1] - zs[0])
        seg = np.arange(i_peak, len(zs))               # upper (stable) branch

        def ReG_radius(a_r):
            xr = (1 + 1j) * a_r / d
            return (-0.5 * (1 - 3 / xr**2 + 3 / xr * _Lm._cot(xr))).real

        def zstar_radius(a_r):
            target = (4.0 / 3.0) * rho_Cu * g * mu0 / ReG_radius(a_r)  # needed dB2/dz (<0)
            dd = dB2dz_grid[seg]                        # rises from very neg toward 0
            cr = np.where(dd >= target)[0]
            if cr.size == 0:
                return None
            j = seg[cr[0]]
            z0, z1 = zs[j - 1], zs[j]
            d0, d1 = dB2dz_grid[j - 1], dB2dz_grid[j]
            return z0 + (target - d0) * (z1 - z0) / (d1 - d0)

        print("\n[check 5] z* vs sphere radius (lift & weight both ~ a^3)")
        print("   a[mm]   a/delta   Re[G]     z*[mm]    a/L")
        zstars = []
        for a_r in (5e-3, 2e-3, 1e-3):
            zr = zstar_radius(a_r)
            assert zr is not None, "no equilibrium found for this radius"
            b2_r = float(np.interp(zr, zs, b2))
            dB2_r = float(np.interp(zr, zs, dB2dz_grid))
            aL = a_r * abs(dB2_r) / b2_r
            zstars.append(zr)
            print(f"   {a_r*1e3:4.1f}    {a_r/d:6.2f}   {ReG_radius(a_r):+.4f}   "
                  f"{zr*1e3:7.3f}   {aL:.3f}")
        # a/L must shrink with the sphere (the point dipole becomes clean)
        # and z* shifts only modestly (the shift is Re[G] desaturation only)
        zspread = (max(zstars) - min(zstars)) / np.mean(zstars)
        zr1 = zstar_radius(1e-3)
        aL1 = 1e-3 * abs(float(np.interp(zr1, zs, dB2dz_grid))) / float(np.interp(zr1, zs, b2))
        assert aL1 < 0.15, "a 1 mm sphere should give a clean a/L << 1"
        assert zspread < 0.20, f"z* shift over 5x radius range too large ({zspread*100:.0f}%)"
        print(f"   -> a/L is freely controllable (5 mm: 0.49 -> 1 mm: {aL1:.3f}, clean "
              f"dipole);\n      z* shifts only {zspread*100:.0f}% over a 5x radius range "
              f"(pure Re[G] desaturation)")

        print("\n" + "=" * 70)
        print("ALL CHECKS PASSED")
        print("=" * 70)

        # ----- persist durable validation evidence -----
        results = {
            "description": "AC levitation equilibrium of a Cu sphere above a coil "
                           "(Radia open-boundary field + verified induced-dipole force)",
            "force_formula": "F_z(z) = (pi a^3/mu0) Re[G(omega)] d|B0(z)|^2/dz",
            "sphere": {
                "a_m": a, "sigma_S_per_m": sigma, "rho_kg_per_m3": rho_Cu,
                "volume_m3": SPHERE_VOL, "mass_kg": mass, "weight_N": WEIGHT,
            },
            "coil": {
                "type": "filament loop (ObjFlmCur)", "R_m": R_COIL,
                "N_turn": N_TURN, "I_peak_A": I_PEAK, "NI_amp_turns": NI,
                "n_segments": NSEG, "axis": "z", "plane_z_m": 0.0,
            },
            "frequency": {
                "f_Hz": F_HZ, "omega_rad_per_s": OMEGA,
                "skin_depth_m": d, "a_over_delta": a_over_delta,
            },
            "sphere_response": {"Re_G": ReG, "G_complex": [G_exact(OMEGA).real,
                                                            G_exact(OMEGA).imag]},
            "field_check": {"max_rel_error_percent": max_err, "tol_percent": 1.0},
            "equilibrium": {
                "z_star_m": z_star, "z_star_mm": z_star * 1e3,
                "F_at_z_star_N": F_star, "weight_N": WEIGHT,
                "residual_percent": (F_star - WEIGHT) / WEIGHT * 100.0,
                "peak_lift_N": F_peak, "z_peak_lift_m": z_peak,
                "peak_over_weight": F_peak / WEIGHT,
            },
            "stability": {
                "dFdz_N_per_m": float(dFdz), "stable": bool(dFdz < 0),
                "stiffness_N_per_m": abs(float(dFdz)),
                "osc_freq_Hz": (np.sqrt(-dFdz / mass) / (2 * np.pi)
                                if dFdz < 0 else None),
            },
            "dipole_validity": {
                "grad_scale_L_m": L_grad, "a_over_L": a / L_grad,
                "z_star_over_a": z_star / a,
            },
            "profile": {
                "z_m": zs.tolist(),
                "F_z_N": Fz.tolist(),
                "B2_T2": b2.tolist(),
            },
        }
        out_json = validation_output("coil_maglev_equilibrium_results.json", HERE)
        with open(out_json, "w") as fp:
            json.dump(results, fp, indent=2, allow_nan=False)
        print(f"\nwrote {out_json}")

        plot(zs, Fz, WEIGHT, z_star, F_peak, z_peak)

    finally:
        rad.UtiDelAll()


def plot(zs, Fz, weight, z_star, F_peak, z_peak):
    """Force-vs-height figure.  NO in-figure title (lab figure policy:
    the title belongs in the LaTeX caption)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 10, "font.family": "serif", "axes.linewidth": 0.8,
        "mathtext.fontset": "dejavuserif",
    })
    fig, ax = plt.subplots(figsize=(8.0 / 2.54, 6.0 / 2.54))

    ax.plot(zs * 1e3, Fz * 1e3, "C0-", lw=1.3, label="lift $F_z(z)$")
    ax.axhline(weight * 1e3, color="C3", ls="--", lw=1.0,
               label="weight $mg$")
    ax.plot([z_star * 1e3], [weight * 1e3], "ko", ms=4.5, zorder=5)
    ax.annotate(f"$z^*$={z_star*1e3:.1f} mm",
                xy=(z_star * 1e3, weight * 1e3),
                xytext=(z_star * 1e3 + 25, weight * 1e3 + 0.18 * F_peak * 1e3),
                fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7, color="0.3"))

    ax.set_xlabel("height above coil $z$ (mm)")
    ax.set_ylabel("vertical force (mN)")
    ax.set_xlim(0, zs[-1] * 1e3)
    ax.set_ylim(0, 1.05 * F_peak * 1e3)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)

    out_png = os.path.join(HERE, "coil_maglev_equilibrium.png")
    fig.savefig(out_png, dpi=300)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
