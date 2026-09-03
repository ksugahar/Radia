r"""Combined-function magnet on a CURVED orbit: the Frenet sweep is the twist.

THE CONFLUENCE (rung 1-2 + rung 3)
----------------------------------
rung 1-2 (`ffag_sector_two_plane.py`) bent the beam with a pure DIPOLE sector;
rung 3 (`twisting_quadrupole_pole.py`) twisted a pure QUAD on a fixed station.
This file merges them: a COMBINED-FUNCTION magnet (dipole b1 + quad gradient b2
in ONE cross-section) swept along the CURVED orbit it bends.

In the Frenet frame of the orbit the cross-section is FIXED (the design spec:
b1 bends, b2 focuses).  But the Frenet frame ROTATES with the bend -- by the
bend angle theta(s) = s / rho, rho = (B rho) / b1 -- so in the LAB frame the
whole combined-function pole TWISTS by theta(s).  Rung 3's n-fold law then says:

    geometric roll theta  =>  dipole (n=1) multipole phase  psi_1 = theta,
                              quad   (n=2) multipole phase  psi_2 = 2 theta.

So BOTH the dipole and quad ORIENTATIONS track the same Frenet angle theta (the
rigid roll), while their multipole PHASES differ by the factor n (psi_2 = 2 psi_1)
-- the design-primitive surface (Frenet, fixed) reflected into a twisting lab pole.

THE COMBINED-FUNCTION CROSS-SECTION (the new piece)
---------------------------------------------------
The engineering realization is a TILTED-GAP dipole: pole faces z = +-(g/2 - t x),
so the gap g(x) = g0 - 2 t x narrows toward +x and B_z(x) ~ B0 / g(x) carries a
DIPOLE b1 plus a GRADIENT b2 = (relative) the quad.  (rung 1-2 was a flat gap =
pure dipole; this tilts it.)  A 2-D Laplace solve in the aperture (poles at +-Phi0)
recovers b1 AND b2 together.

WHAT IS VERIFIED (ngsolve only, golden-tested)
----------------------------------------------
- the combined-function cross-section: b1 (dipole) + b2 (gradient ~6%) + a small
  b3 (the 1/g curvature, a real magnet shims it);
- the Frenet sweep: rolling the magnet by theta, BOTH the dipole orientation
  (-psi_1) and the quad orientation (-psi_2/2) track theta (slope 1), and the
  quad multipole phase is exactly TWICE the dipole phase (psi_2 = 2 psi_1, the
  n-fold law) -- the lab pole twists by theta(s).

HONEST SCOPE
------------
A pure sector (rigid Frenet roll; no spiral edge, no s-varying gradient).  The
per-station 2-D design is the SLOW-BEND limit; a spiral sector (the pole twist
phi != the orbit bend theta) and an s-ramped (b1(s), b2(s)) combined function are
the extensions, and the fast-twist leaf coupling (when does the per-station 2-D
break) is the next rung.

run:  python combined_function_frenet_sweep.py          # cross-section + Frenet sweep
      python combined_function_frenet_sweep.py --fig      # + figure
"""

from _validation_output import validation_output
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))      # the analyzer

from accel_pole_design import multipoles                            # noqa: E402

G2 = 0.020                 # half-gap (flat-gap value at x=0)
HALF_W = 0.045             # pole half-width (x)
POLE_T = 0.030             # pole thickness (z)
X_AIR = 0.22               # air box half-width (x)
Z_AIR = 0.16               # air box half-height (z)
R_REF = 0.008              # multipole reference-circle radius
TILT = 0.15                # gap-tilt -> the gradient (combined-function quad)
BRHO = 1.0                 # beam rigidity B*rho (normalised) for the bend kinematics


def _pole_face_z(x, tilt):
    """Combined-function pole face: z = g/2 - tilt*x (the tilted gap -> gradient)."""
    return G2 - tilt * x


def build_combined_pole(tilt=TILT, roll_deg=0.0, n_face=40):
    """Full-aperture tilted-gap dipole (top + bottom pole islands), the whole
    magnet ROLLED by roll_deg (the Frenet frame at a station s).  Edges named on
    the unrolled shape (poleT / poleB / far), then rotated (names carry through).
    """
    from netgen.occ import (WorkPlane, OCCGeometry, Axis, Pnt, Z as ZDIR)
    xs = np.linspace(-HALF_W, HALF_W, n_face)
    zf = _pole_face_z(xs, tilt)
    wp = WorkPlane().MoveTo(float(xs[0]), float(zf[0]))
    for xi, zi in zip(xs[1:], zf[1:]):
        wp.LineTo(float(xi), float(zi))
    wp.LineTo(float(HALF_W), G2 + POLE_T)
    wp.LineTo(float(-HALF_W), G2 + POLE_T)
    wp.Close()
    top = wp.Face()
    wp2 = WorkPlane().MoveTo(float(xs[0]), float(-zf[0]))
    for xi, zi in zip(xs[1:], zf[1:]):
        wp2.LineTo(float(xi), float(-zi))
    wp2.LineTo(float(HALF_W), -(G2 + POLE_T))
    wp2.LineTo(float(-HALF_W), -(G2 + POLE_T))
    wp2.Close()
    bot = wp2.Face()
    box = WorkPlane().MoveTo(-X_AIR, -Z_AIR).Rectangle(2 * X_AIR, 2 * Z_AIR).Face()
    air = box - top - bot
    for e in air.edges:
        c = e.center
        if abs(c.x) > X_AIR - 1e-7 or abs(c.y) > Z_AIR - 1e-7:
            e.name = "far"
        elif c.y > 0:
            e.name = "poleT"
        else:
            e.name = "poleB"
    if abs(roll_deg) > 1e-12:
        air = air.Rotate(Axis(Pnt(0, 0, 0), ZDIR), float(roll_deg))
    return OCCGeometry(air, dim=2)


def solve_combined(tilt=TILT, roll_deg=0.0, order=3, maxh=0.006):
    """Laplace solve of the (rolled) combined-function pole; return the lab-frame
    multipoles + the dipole / quad phases."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, TaskManager)
    geo = build_combined_pole(tilt, roll_deg)
    with TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
        fes = H1(mesh, order=order, dirichlet="poleT|poleB|far")
        u, v = fes.TnT()
        gf = GridFunction(fes)
        gf.Set(mesh.BoundaryCF({"poleT": 1.0, "poleB": -1.0}, default=0.0),
               definedon=mesh.Boundaries("poleT|poleB|far"))
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        rhs = gf.vec.CreateVector()
        rhs.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rhs
        B = -grad(gf)

        def B_perp(x, y):
            vv = B(mesh(x, y))
            return (float(vv[0]), float(vv[1]))

        mp = multipoles(B_perp, R_REF, n_max=8, n_samples=400)
        ne = int(mesh.ne)
    c1 = complex(*mp[1])
    c2 = complex(*mp[2])
    b1 = abs(c1)
    psi1 = math.atan2(c1.imag, c1.real)            # dipole multipole phase
    psi2 = math.atan2(c2.imag, c2.real)            # quad   multipole phase
    return {
        "roll_deg": float(roll_deg), "tilt": float(tilt), "ne": ne,
        "b1": float(b1), "b2_rel": float(abs(c2) / b1),
        "b3_rel": float(abs(complex(*mp[3])) / b1),
        "psi1_deg": float(math.degrees(psi1)),
        "psi2_deg": float(math.degrees(psi2)),
        # the geometric roll recovered from each harmonic (convention sign: the
        # multipole phase turns -n*theta, so theta = -psi_n / n)
        "roll_from_dipole_deg": float(-math.degrees(psi1)),
        "roll_from_quad_deg": float(-0.5 * math.degrees(psi2)),
    }


def frenet_sweep(rolls_deg=(0.0, 10.0, 20.0, 30.0, 40.0), tilt=TILT, order=3):
    """Sweep the Frenet roll theta; confirm BOTH the dipole and quad orientations
    track theta (the rigid Frenet rotation) and the quad multipole-phase CHANGE
    is 2x the dipole's (the n-fold law).  Phases are referenced to theta=0 (the
    magnet's design orientation has its own constant offset -- the dipole points
    'down', psi ~ 180 deg) and np.unwrap'd in theta."""
    rows = [solve_combined(tilt, r, order=order) for r in rolls_deg]
    theta = np.array(rolls_deg, dtype=float)
    psi1 = np.unwrap(np.radians([r["psi1_deg"] for r in rows]))
    psi2 = np.unwrap(np.radians([r["psi2_deg"] for r in rows]))
    dpsi1 = np.degrees(psi1 - psi1[0])             # dipole phase CHANGE from 0
    dpsi2 = np.degrees(psi2 - psi2[0])             # quad   phase CHANGE from 0
    # the geometric roll recovered from each harmonic (phase turns -n*theta):
    dip = -dpsi1
    quad = -0.5 * dpsi2
    slope_dip = float(np.polyfit(theta, dip, 1)[0])
    slope_quad = float(np.polyfit(theta, quad, 1)[0])
    sp1 = float(np.polyfit(theta, dpsi1, 1)[0])
    sp2 = float(np.polyfit(theta, dpsi2, 1)[0])
    phase_ratio = sp2 / sp1 if abs(sp1) > 1e-9 else float("nan")   # n-fold: 2
    return {
        "rows": rows, "theta_deg": theta.tolist(),
        "dipole_roll_deg": dip.tolist(), "quad_roll_deg": quad.tolist(),
        "slope_dipole": slope_dip, "slope_quad": slope_quad,
        "track_err_dipole_deg": float(np.max(np.abs(dip - theta))),
        "track_err_quad_deg": float(np.max(np.abs(quad - theta))),
        "dpsi1_slope": sp1, "dpsi2_slope": sp2,
        "phase_ratio_n2_over_n1": phase_ratio,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    print("=" * 76)
    print("Combined-function magnet on a CURVED orbit -- the Frenet sweep IS the twist")
    print("=" * 76)
    base = solve_combined(roll_deg=0.0, order=args.order)
    rho = BRHO / base["b1"]
    print("Cross-section (Frenet frame) -- the combined-function pole:")
    print(f"  main dipole b1              : {base['b1']:.3f}  (ne={base['ne']})")
    print(f"  gradient (quad) b2/b1       : {base['b2_rel']:.4f}  (the tilted gap)")
    print(f"  sextupole b3/b1             : {base['b3_rel']:.3e}  (1/g curvature, shimmable)")
    print(f"  bend radius rho = Brho/b1   : {rho:.4f}  (the dipole bends; theta=s/rho)")

    print("\nThe FRENET sweep -- roll by theta, the lab pole twists:")
    sw = frenet_sweep(tilt=TILT, order=args.order)
    print(f"  {'theta':<8}{'roll<-dipole':<14}{'roll<-quad':<13}")
    for th, do, qo in zip(sw["theta_deg"], sw["dipole_roll_deg"], sw["quad_roll_deg"]):
        print(f"  {th:<8.1f}{do:<14.2f}{qo:<13.2f}")
    print(f"  => dipole orientation tracks theta: slope {sw['slope_dipole']:.3f} "
          f"(err {sw['track_err_dipole_deg']:.2f} deg)")
    print(f"  => quad   orientation tracks theta: slope {sw['slope_quad']:.3f} "
          f"(err {sw['track_err_quad_deg']:.2f} deg)")
    print(f"  => the n-fold law: quad phase change = "
          f"{sw['phase_ratio_n2_over_n1']:.3f} x dipole phase change (ideal 2.000)")
    print("  => the combined-function pole, fixed in the Frenet frame, twists by")
    print("     theta(s) in the lab -- dipole & quad roll together (the Frenet")
    print("     rotation), their multipole phases in the n:1 ratio (the n-fold law).")

    out = {
        "brho": BRHO, "rho": rho, "r_ref": R_REF,
        "base": {k: base[k] for k in base},
        "sweep": {k: sw[k] for k in sw if k != "rows"},
        "sweep_rows": [{k: r[k] for k in r} for r in sw["rows"]],
    }
    jpath = validation_output("combined_function_frenet_sweep.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(sw, base)


def _figure(sw, base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.2))
    # LEFT: the combined-function pole faces at three Frenet stations (rolled)
    xs = np.linspace(-HALF_W, HALF_W, 60)
    for th, col in ((0.0, "C0"), (20.0, "C1"), (40.0, "C2")):
        ct, st = math.cos(math.radians(th)), math.sin(math.radians(th))
        for sgn in (+1, -1):
            zf = sgn * _pole_face_z(xs, TILT)
            xr = xs * ct - zf * st
            zr = xs * st + zf * ct
            ax[0].plot(xr * 1e3, zr * 1e3, col, lw=2,
                       label=(f"theta={th:.0f} deg" if sgn > 0 else None))
    th = np.linspace(0, 2 * math.pi, 200)
    ax[0].plot(R_REF * 1e3 * np.cos(th), R_REF * 1e3 * np.sin(th), "k--", lw=0.8,
               label=f"orbit aperture r={R_REF*1e3:.0f} mm")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel("z [mm]")
    ax[0].set_title("Combined-function pole (tilted gap = dipole+quad)\n"
                    "rolling with the Frenet frame theta(s) = the lab twist")
    ax[0].legend(fontsize=8, loc="upper right")
    # RIGHT: dipole & quad orientations track theta; the phase-change ratio is 2
    th = np.array(sw["theta_deg"])
    ax[1].plot(th, sw["dipole_roll_deg"], "o-", color="C0", ms=5,
               label=f"roll <- dipole (slope {sw['slope_dipole']:.2f})")
    ax[1].plot(th, sw["quad_roll_deg"], "s-", color="C3", ms=5,
               label=f"roll <- quad (slope {sw['slope_quad']:.2f})")
    ax[1].plot(th, th, "k--", lw=1, label="ideal = theta")
    ax[1].set_xlabel("Frenet roll  theta(s) [deg]")
    ax[1].set_ylabel("recovered roll [deg]")
    ax[1].set_title(f"Both roll with the Frenet frame (slope 1);\n"
                    f"multipole phase changes in the n:1 ratio "
                    f"(n=2:n=1 = {sw['phase_ratio_n2_over_n1']:.2f})")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


if __name__ == "__main__":
    main()
