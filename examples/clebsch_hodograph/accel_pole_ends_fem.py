"""(A) The magnet ENDS, FEM rung: reduced potential + CoilBuilder, finite length.

RESEARCH example (track A, the July main project).  The next rung after the
analytic theorem in accel_pole_ends_3d.py: replace the analytic Maxwellian field
with a REAL finite-length magnet solved by the lab's forward engine --

    reduced scalar potential Omega  (H = H_s - grad Omega,  H_s = CoilBuilder
    Biot-Savart, NO coil mesh) + an iron yoke, on a self-contained netgen.occ
    mesh (no Cubit) --

and read the result through the INTEGRATED multipole analyzer (the beam-optics
quantity INT B_perp dl).  This file is the bridge "analytic -> FEM": it shows the
forward engine reproduces the integrated-field theorem on a real device, and it
exposes the magnet ENDS where the next design step (read the 3-D equipotential
surface as the end-iron contour) will act.

Convention (matches examples/c_type_electromagnet/dipole_with_coilbuilder.py):
  beam  = y   (the long / straight-section axis; INT ... dy is the beam integral)
  gap   = z   (pole faces at z = +-g/2; the dipole field is B_z in the gap)
  width = x
A C-frame dipole: top + bottom iron bars joined by a back leg (flux return),
finite length L along the beam, so the magnet has real ENDS at y = +-L/2.

Forward formulation (unified reduced potential, all regions):
    INT mu grad(Omega) . grad(v) dx = INT mu H_s . grad(v) dx
    H = H_s - grad(Omega) ,   B = mu H .

Readout: the integrated transverse multipoles of B along the beam (via
accel_pole_ends_3d.integrated_multipoles, with a coordinate adapter beam=y).
The main is the integrated dipole bbar_1 = INT B_z dy; the spurious content is
the integrated higher harmonics the ENDS contribute.

run:  python accel_pole_ends_fem.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))             # the analyzer
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "src", "radia"))         # radia_coil_builder

from accel_pole_ends_3d import integrated_multipoles                       # noqa: E402

MU0 = 4 * math.pi * 1e-7

# ---- geometry (meters); beam=y, gap=z, width=x ----
GAP = 0.040            # pole-face separation (z)
POLE_W = 0.060         # pole / bar width (x half-extent = POLE_W)
Z_OUT = 0.100          # outer z of the bars
T_LEG = 0.030          # back-leg thickness (x)
L_BEAM = 0.120         # iron length along the beam (y) -> ends at +-L/2
AIR = 0.30             # air box half-size
# ---- coil (CoilBuilder racetrack pair) ----
COIL_W = 0.015
COIL_H = 0.025
COIL_STRAIGHT = 0.220  # > iron length L_BEAM so the body sees uniform field
COIL_R = 0.020         # (coil end-turns sit OUTSIDE the iron ends)
NI = 10000.0           # ampere-turns


def build_mesh(maxh_air=0.05, maxh_iron=0.025):
    """H-frame (window) dipole (iron) + air box, netgen.occ (no Cubit).

    x-SYMMETRIC (legs on both +-x) so the dipole field is clean (B_x ~ 0 at
    centre by symmetry).  Finite length L along the beam (y) -> real ENDS.
    """
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry

    hL = L_BEAM / 2
    top = Box(Pnt(-POLE_W, -hL, GAP / 2), Pnt(POLE_W, hL, Z_OUT))
    bot = Box(Pnt(-POLE_W, -hL, -Z_OUT), Pnt(POLE_W, hL, -GAP / 2))
    leg_l = Box(Pnt(-POLE_W, -hL, -GAP / 2), Pnt(-POLE_W + T_LEG, hL, GAP / 2))
    leg_r = Box(Pnt(POLE_W - T_LEG, -hL, -GAP / 2), Pnt(POLE_W, hL, GAP / 2))
    iron = (top + bot + leg_l + leg_r)
    iron.mat("iron")
    iron.maxh = maxh_iron

    air = Box(Pnt(-AIR, -AIR, -AIR), Pnt(AIR, AIR, AIR)) - iron
    air.mat("air")
    for f in air.faces:
        c = f.center
        if max(abs(c.x), abs(c.y), abs(c.z)) > 0.9 * AIR:
            f.name = "outer"

    shape = Glue([air, iron])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh_air))
    return mesh


def build_coil():
    """Racetrack coil pair (CoilBuilder) -> Radia container.  Biot-Savart only."""
    import radia as rad
    from radia_coil_builder import CoilBuilder

    rad.UtiDelAll()
    z_coil = GAP / 2 + COIL_H / 2
    # the 180 arcs curve toward -x, so start at x = +COIL_R -> the two straights
    # land at x = +-COIL_R, i.e. the racetrack is centred on x = 0 (verified: Hs_x=0
    # and Hs_z peaks at x=0).  No spurious transverse field on the beam axis.
    upper = (CoilBuilder(current=NI)
             .set_start([COIL_R, -COIL_STRAIGHT / 2, z_coil])
             .set_cross_section(width=COIL_W, height=COIL_H)
             .add_straight(COIL_STRAIGHT)
             .add_arc(radius=COIL_R, arc_angle=180)
             .add_straight(COIL_STRAIGHT)
             .add_arc(radius=COIL_R, arc_angle=180))
    lower = upper.mirror("xy")                       # z -> -z, current reversed
    coils = rad.ObjCnt(upper.to_radia() + lower.to_radia())
    return coils


def solve(mu_r=1000.0, order=2, maxh_air=0.05, maxh_iron=0.025,
          r_ref=0.008, n_beam=121, n_theta=32, plot=False):
    """Solve reduced-Omega forward, then integrate the transverse multipoles
    along the beam.  Returns the gap field + the integrated dipole + spectrum."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    import radia as rad

    mesh = build_mesh(maxh_air, maxh_iron)
    coils = build_coil()
    Hs = rad.RadiaField(coils, "h")                  # Biot-Savart source CF
    mu = mesh.MaterialCF({"iron": mu_r * MU0}, default=MU0)

    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()

    # The RadiaField (Biot-Savart) CF is NOT thread-safe: assembling the source
    # LinearForm (and the field readout below) must run SERIALLY -- under
    # TaskManager the parallel element loop deadlocks on Radia's global state.
    # Only the pure-NGSolve stiffness assemble + solve are wrapped in TaskManager.
    f = LinearForm(fes)
    f += mu * Hs * grad(v) * dx
    f.Assemble()                                     # serial (RadiaField source)

    with TaskManager():
        a = BilinearForm(fes)
        a += mu * grad(u) * grad(v) * dx
        a.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    B = mu * (Hs - grad(gfu))                        # B = mu (Hs - grad Omega)

    # ---- field readout: SERIAL (RadiaField eval) ----
    bx_centre = float(B(mesh(0.0, 0.0, 0.0))[0])
    bz_centre = float(B(mesh(0.0, 0.0, 0.0))[2])

    # on-axis B_z(y) profile along the beam (the fringe / end falloff)
    ys = np.linspace(-1.4 * (L_BEAM / 2), 1.4 * (L_BEAM / 2), 81)
    bz_axis = np.array([float(B(mesh(0.0, y, 0.0))[2]) for y in ys])

    # INTEGRATED transverse multipoles along the beam (y).  Adapter: the
    # analyzer integrates its 3rd arg ("z") sampling a circle in (x,y);
    # map analyzer-(x,y,z) -> global (x, z=transverse, y=beam):
    #   global point = (ax, az, ay) ;  transverse field = (Bx_g, Bz_g).
    def B_perp(ax, ay, az):
        b = B(mesh(ax, az, ay))
        return (float(b[0]), float(b[2]))            # (analyzer Bx, By) = (Bx_g, Bz_g)

    y_int = 1.6 * (L_BEAM / 2)                        # integrate past the ends
    mp = integrated_multipoles(B_perp, r_ref, (-y_int, y_int),
                               n_z=n_beam, n_theta=n_theta, n_max=8)

    # ---- (B) the equipotential END contour: read Psi=Psi_pole as the end edge ----
    # In the current-free gap H = -grad(Psi_total); along z at fixed (x=0, y):
    #   Psi(y, z) = - INT_0^z H_z dz'   (Psi(y,0)=0 by the dipole's up-down antisym).
    # The iron face is an equipotential: in the body Psi_pole = Psi(0, g/2).  The
    # curve z_p(y) with Psi(y, z_p) = Psi_pole is the ideal end-iron contour --
    # in the body z_p = g/2, and past the iron end (y > L/2) it CURVES (the
    # field bows out): that curve is the end chamfer the design should follow.
    H_cf = Hs - grad(gfu)
    zc = np.linspace(0.0, 2.0 * (GAP / 2), 70)        # z: gap -> above the pole
    yc = np.linspace(0.0, 1.30 * (L_BEAM / 2), 70)    # beam: body -> past the end
    psi = np.zeros((len(yc), len(zc)))
    for i, yv in enumerate(yc):
        hz = np.array([float(H_cf(mesh(0.0, yv, zv))[2]) for zv in zc])
        psi[i, 1:] = -np.cumsum(0.5 * (hz[1:] + hz[:-1]) * np.diff(zc))
    psi_pole = float(np.interp(GAP / 2, zc, psi[0]))  # body iron-face value at z=g/2
    # Psi decreases in z (H_z > 0) -> interp on (-Psi), which increases in z.
    z_pole = np.array([float(np.interp(-psi_pole, -psi[i], zc)) for i in range(len(yc))])
    end_lift = float(np.max(z_pole) - GAP / 2)        # how far the equipotential lifts

    rad.UtiDelAll()

    main = abs(complex(*mp[1]))                       # integrated dipole bbar_1
    allowed = (1, 3, 5)                               # dipole-allowed normals
    spurious = max(abs(complex(*mp[n])) / main for n in allowed if n != 1)
    # effective magnetic length L_eff = INT B_z dy / B_z(body); > L_iron by the
    # two fringes.  B_z(body) = mean over the flat central 60% of the iron.
    body = np.abs(ys) < 0.30 * L_BEAM
    bz_body = float(np.mean(bz_axis[body]))
    l_eff = float(np.trapezoid(bz_axis, ys) / bz_body) if abs(bz_body) > 1e-9 else 0.0
    # END ENHANCEMENT: the pole-end flux concentration (peak |B_z| near the iron
    # end vs the body) -- the design-relevant end effect.
    end = np.abs(np.abs(ys) - 0.5 * L_BEAM) < 0.12 * L_BEAM
    end_overshoot = float(np.max(np.abs(bz_axis[end])) / abs(bz_body)) - 1.0

    if plot:
        _plot(ys, bz_axis, bz_body, l_eff, yc, z_pole)

    return {
        "ne": int(mesh.ne), "ndof": int(fes.ndof),
        "bz_gap_centre_T": bz_centre,
        "bz_body_T": bz_body,                         # flat-top body field
        "bx_gap_centre_T": bx_centre,                 # ~0 for the x-symmetric H-frame
        "bx_over_bz_centre": float(abs(bx_centre) / (abs(bz_centre) + 1e-30)),
        "integrated_dipole_bbar1_Tm": float(main),
        "integrated_spurious_rel": float(spurious),   # higher harmonics (ends + finite pole)
        "L_eff_m": l_eff,
        "end_overshoot": end_overshoot,               # pole-end flux concentration
        "psi_pole": psi_pole,                          # iron-face equipotential value
        "end_contour_lift_m": end_lift,               # equipotential lift at the end
        "z_pole_body_m": float(z_pole[0]),            # contour in the body (~ g/2)
    }


def _plot(ys, bz_axis, bz_body, l_eff, yc, z_pole):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.8), dpi=150)

    # LEFT: on-axis B_z(y) -- flat body + pole-end enhancement + fringe
    ax[0].plot(ys * 1e3, bz_axis, "C0")
    ax[0].fill_between(ys * 1e3, bz_axis, 0, alpha=0.2, color="C0")
    ax[0].axvspan(-L_BEAM / 2 * 1e3, L_BEAM / 2 * 1e3, color="0.85", zorder=0,
                  label="iron length L")
    ax[0].axhline(bz_body, color="C3", lw=0.8, ls="--",
                  label=f"$B_z$(body)={bz_body:.4f} T")
    ax[0].set_xlabel("beam  y [mm]")
    ax[0].set_ylabel("$B_z$ on axis [T]")
    ax[0].set_title(f"Forward solve: flat body + pole-END enhancement + fringe\n"
                    f"($L_{{eff}}={l_eff*1e3:.0f}$ mm, iron $L={L_BEAM*1e3:.0f}$ mm)")
    ax[0].legend(fontsize=8, loc="center")

    # RIGHT: the equipotential END contour z_p(y) -- the ideal end-iron edge
    ax[1].plot(yc * 1e3, z_pole * 1e3, "C0", lw=2,
               label="equipotential $\\Psi=\\Psi_{pole}$ (ideal end edge)")
    ax[1].plot([0, L_BEAM / 2 * 1e3], [GAP / 2 * 1e3, GAP / 2 * 1e3], "k--", lw=1,
               label="straight pole face $z=g/2$")
    ax[1].axvline(L_BEAM / 2 * 1e3, color="0.6", lw=0.8, ls=":")
    ax[1].annotate("iron end", (L_BEAM / 2 * 1e3, GAP / 2 * 1e3 * 1.02),
                   fontsize=8, color="0.4")
    ax[1].set_xlabel("beam  y [mm]")
    ax[1].set_ylabel("pole edge  z [mm]")
    ax[1].set_title("Read the 3-D equipotential as the end-iron contour\n"
                    "(curves past the iron end = the chamfer to follow)")
    ax[1].legend(fontsize=8, loc="upper left")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("(A) Magnet ENDS, FEM rung -- reduced potential + CoilBuilder, "
          "finite length\n")
    r = solve(plot=True)
    print(f"  mesh: ne={r['ne']}, ndof={r['ndof']}")
    print(f"  gap field  B_z(body) = {r['bz_body_T']:.5f} T  "
          f"(B_x/B_z at centre = {r['bx_over_bz_centre']:.1e} -> clean dipole)")
    print(f"  effective magnetic length L_eff = {r['L_eff_m']*1e3:.1f} mm "
          f"(iron L = {L_BEAM*1e3:.0f} mm; the 2 fringes add ~"
          f"{(r['L_eff_m']*1e3 - L_BEAM*1e3)/2:.0f} mm each)")
    print(f"  pole-END enhancement (peak/body - 1) = {r['end_overshoot']*100:.1f}% "
          f"-- the end flux concentration")
    print(f"  INTEGRATED dipole  bbar_1 = INT B_z dy = "
          f"{r['integrated_dipole_bbar1_Tm']:.5f} T.m")
    print(f"  integrated spurious |bbar_n/bbar_1| (n=3,5) = "
          f"{r['integrated_spurious_rel']:.2e}  (ends + finite pole width)")
    print(f"\n  (B) equipotential END contour (the ideal end-iron edge):")
    print(f"      body contour z_p = {r['z_pole_body_m']*1e3:.2f} mm "
          f"(= g/2 = {GAP/2*1e3:.1f} mm, recovers the flat pole face)")
    print(f"      equipotential LIFT at the end = {r['end_contour_lift_m']*1e3:.2f} mm "
          f"-- the field bows out past the iron end (the chamfer to follow)")
    print("\n  => the reduced-Omega + CoilBuilder forward engine feeds the SAME")
    print("     INTEGRATED analyzer as the analytic theorem (accel_pole_ends_3d),")
    print("     and the solved 3-D equipotential gives the end-iron contour.")
    print("     NEXT: re-shape the end iron to that contour -> re-solve -> the")
    print("     integrated spurious bbar_n drops (closes the DESIGN loop).")


if __name__ == "__main__":
    main()
