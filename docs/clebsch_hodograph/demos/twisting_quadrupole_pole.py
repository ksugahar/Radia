r"""Beam-referenced equipotential surface as the design primitive: a TWISTING
quadrupole (the curved-orbit / combined-function twist axis).

THE DESIGN PRIMITIVE (the reframing)
------------------------------------
Instead of solving a magnet and then multipole-EXPANDING the field, make the
**beam-referenced equipotential surface** the design SPEC.  In the Frenet frame
of the orbit `s` the iron pole face is a magnetic-scalar equipotential
`Omega(r, theta; s) = sum_n r^n b_n(s) sin(n theta + phi_n(s))` (high-mu =>
H_tangential = 0 => Omega = const), so the multipole `(b_n, phi_n)(s)` IS the
surface's angular Fourier mode.  Design = prescribe `(b_n, phi_n)(s)`, sweep the
equipotential surface along the orbit, place iron there.

THE TWIST (this file: the quadrupole)
-------------------------------------
The genuinely 3-D content of a CURVED-orbit / COMBINED-function magnet is that
the transverse multipole **rotates** along `s` (the Frenet frame turns with the
bend; a rotating-gradient / twisted quad turns the pole on purpose).  The key
fact is the **n-fold law**: rotating the equipotential SURFACE by `phi` rotates
the order-`n` multipole PHASE by `n phi`.  For the quadrupole (`n = 2`):

    rotate the pole by phi   <=>   (b_2, a_2) -> |b_2| (cos 2phi, sin 2phi)

so a quad twisted by 45 deg becomes a pure SKEW quad.  The recovered orientation
`alpha = (1/2) atan2(a_2, b_2)` tracks the prescribed `phi`.

WHAT IS VERIFIED (real FEM, ngsolve only)
-----------------------------------------
The quad pole face is the hyperbola `x y = +-r0^2/2` (the `Omega = const`
equipotential, `accel_pole_design.quad_pole_hyperbola`).  A 2-D Laplace solve in
the aperture with the 4 hyperbola poles at alternating `+-Omega0` (and the gaps
natural) recovers a clean quad (`a_2 ~ 0`, forbidden `n = 1, 3` at the numerical
floor, the dominant allowed spurious the finite-pole 12-pole `b_6`).  Rotating
the 4 poles by `phi` rotates the recovered orientation by exactly `phi`
(slope 1) -- the twist, measured.

HONEST SCOPE
------------
This is the per-station (Frenet cross-section) 2-D design: the SLOW-TWIST
(adiabatic) limit `d phi / ds -> 0`, where the magnet is a stack of 2-D leaves
(the rung-1 foliate-and-perturb picture, now twisting).  A fast twist / tight
bend couples adjacent leaves (a longitudinal-field correction) -- the twist rate
`d phi / ds` is a leaf-coupling perturbation parameter, the next rung; the
combined-function (dipole + quad together, a shifted+rotated hyperbola) and the
genuine curved-orbit Frenet sweep are the extensions.

run:  python twisting_quadrupole_pole.py            # design primitive + twist sweep
      python twisting_quadrupole_pole.py --fig        # + figure
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

R0 = 0.02                    # aperture (inscribed) radius -> hyperbola xy = r0^2/2
C_HYP = R0 ** 2 / 2.0        # the pole equipotential constant
BETA_DEG = 28.0             # pole angular half-width about each diagonal
R_GAP_FAC = 1.7            # gap radius = R_GAP_FAC * r0 (the inter-pole gap arc)
R_REF_FAC = 0.6            # multipole reference circle = R_REF_FAC * r0


def _sector(thd, phi_deg, beta=BETA_DEG, r_gap_fac=R_GAP_FAC):
    """At actual polar angle `thd` (deg) for a quad rotated by `phi_deg`, return
    (radius, label).  Pole sectors follow the rotated hyperbola
    r = sqrt(2c / |sin 2(thd - phi)|); the gaps are a circular arc r_gap."""
    rel = thd - phi_deg
    for d, lab in ((45, "poleP"), (135, "poleM"), (225, "poleP"), (315, "poleM")):
        if abs(((rel - d + 180) % 360) - 180) <= beta:
            s2 = abs(math.sin(2.0 * math.radians(rel)))
            r = math.sqrt(2.0 * C_HYP / max(s2, 1e-9))
            return r, lab
    return r_gap_fac * R0, "gap"


def build_quad_aperture(phi_deg=0.0, n_bnd=360, beta=BETA_DEG):
    """The 2-D quad aperture bounded by 4 hyperbola pole faces (Dirichlet
    +-Omega0) and 4 gap arcs (natural), the whole pattern rotated by `phi_deg`.
    Single-region polygon (robust); edges named poleP / poleM / gap."""
    from netgen.occ import WorkPlane, OCCGeometry
    pts = []
    for i in range(n_bnd):
        thd = 360.0 * i / n_bnd
        r, _ = _sector(thd, phi_deg, beta)
        th = math.radians(thd)
        pts.append((r * math.cos(th), r * math.sin(th)))
    wp = WorkPlane().MoveTo(*pts[0])
    for p in pts[1:]:
        wp.LineTo(*p)
    wp.Close()
    face = wp.Face()
    for e in face.edges:
        cc = e.center
        thd = math.degrees(math.atan2(cc.y, cc.x)) % 360.0
        _, lab = _sector(thd, phi_deg, beta)
        e.name = lab
    return OCCGeometry(face, dim=2)


def solve_quad(phi_deg=0.0, order=4, maxh_fac=6.0, n_bnd=360):
    """Laplace solve of the (rotated) quad pole; return the recovered multipoles
    and the quad orientation alpha = (1/2) atan2(a_2, b_2)."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, TaskManager)
    geo = build_quad_aperture(phi_deg, n_bnd)
    with TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=R0 / maxh_fac))
        fes = H1(mesh, order=order, dirichlet="poleP|poleM")
        u, v = fes.TnT()
        gf = GridFunction(fes)
        gf.Set(mesh.BoundaryCF({"poleP": 1.0, "poleM": -1.0}, default=0.0),
               definedon=mesh.Boundaries("poleP|poleM"))
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        rhs = gf.vec.CreateVector()
        rhs.data = -a.mat * gf.vec
        gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rhs
        B = -grad(gf)
        ne = int(mesh.ne)

        def B_perp(x, y):
            vv = B(mesh(x, y))
            return (float(vv[0]), float(vv[1]))

        mp = multipoles(B_perp, R_REF_FAC * R0, n_max=8, n_samples=400)
    b2, a2 = mp[2]
    main = abs(complex(b2, a2))
    # the physical pole orientation: the geometric twist phi maps to the field
    # multipole PHASE -2 phi (the By+iBx convention sign), so the pole angle is
    # alpha = -(1/2) atan2(a_2, b_2) -- it tracks the prescribed phi (slope +1).
    alpha = -0.5 * math.degrees(math.atan2(a2, b2))            # quad pole orientation
    forbidden = max(abs(complex(*mp[n])) / main for n in (1, 3, 5))
    b6_rel = abs(complex(*mp[6])) / main                       # leading allowed spurious
    return {"phi_deg": float(phi_deg), "b2": float(b2), "a2": float(a2),
            "main": float(main), "alpha_deg": float(alpha),
            "forbidden_max_rel": float(forbidden), "b6_rel": float(b6_rel),
            "ne": ne}


def twist_sweep(phis_deg=(0.0, 15.0, 30.0, 45.0, 60.0, 75.0), order=4):
    """Sweep the prescribed twist phi and confirm the recovered quad orientation
    TRACKS it (n-fold law: surface rotated by phi <=> multipole phase 2 phi, so
    alpha = phi).  Returns the rows + the tracking slope."""
    rows = [solve_quad(p, order=order) for p in phis_deg]
    # unwrap the recovered orientation (90 deg period) onto the prescribed phi
    prescribed = np.array([r["phi_deg"] for r in rows])
    rec = np.array([r["alpha_deg"] for r in rows])
    rec_unwrapped = rec.copy()
    for i in range(len(rec)):
        while rec_unwrapped[i] - prescribed[i] > 45.0:
            rec_unwrapped[i] -= 90.0
        while rec_unwrapped[i] - prescribed[i] < -45.0:
            rec_unwrapped[i] += 90.0
    slope = float(np.polyfit(prescribed, rec_unwrapped, 1)[0])
    track_err = float(np.max(np.abs(rec_unwrapped - prescribed)))
    return {"rows": rows, "prescribed_deg": prescribed.tolist(),
            "recovered_deg": rec_unwrapped.tolist(),
            "tracking_slope": slope, "tracking_err_deg": track_err}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", type=int, default=4)
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    print("=" * 76)
    print("Twisting quadrupole -- the beam-referenced equipotential surface as")
    print("the design primitive (the curved-orbit / combined-function twist)")
    print("=" * 76)
    base = solve_quad(0.0, order=args.order)
    print("Plane (transverse) -- the quad pole = the hyperbola equipotential:")
    print(f"  hyperbola pole face xy      : +-r0^2/2 = {C_HYP:.3e}")
    print(f"  recovered main quad |c_2|   : {base['main']:.3f}  (ne={base['ne']})")
    print(f"  skew a_2 / |c_2| at phi=0   : {abs(base['a2'])/base['main']:.2e}  (pure normal)")
    print(f"  forbidden n=1,3,5 (max rel) : {base['forbidden_max_rel']:.2e}  (-> floor)")
    print(f"  leading allowed spurious b_6: {base['b6_rel']:.2e}  (finite-pole 12-pole)")

    print("\nThe TWIST -- rotate the pole by phi, the multipole phase turns 2 phi:")
    sw = twist_sweep(order=args.order)
    print(f"  {'phi (deg)':<11}{'alpha_rec (deg)':<17}{'a2/b2':<11}{'b6 rel':<10}")
    for r, rec in zip(sw["rows"], sw["recovered_deg"]):
        ab = r["a2"] / r["b2"] if abs(r["b2"]) > 1e-30 else float("inf")
        print(f"  {r['phi_deg']:<11.1f}{rec:<17.2f}{ab:<+11.3f}{r['b6_rel']:<10.2e}")
    print(f"  => recovered orientation TRACKS the prescribed twist: "
          f"slope {sw['tracking_slope']:.3f} (ideal 1.000), "
          f"max err {sw['tracking_err_deg']:.2f} deg")
    print("  => the surface angular mode IS the multipole: twisting the pole")
    print("     surface by phi(s) is a controlled (b_2, a_2)(s) rotation.")

    out = {
        "hyperbola_c": C_HYP, "r0": R0,
        "base": {k: base[k] for k in base},
        "twist": {k: sw[k] for k in sw if k != "rows"},
        "twist_rows": [{k: r[k] for k in r} for r in sw["rows"]],
    }
    jpath = validation_output("twisting_quadrupole_pole.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(sw)


def _figure(sw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.2))
    # LEFT: the pole faces (hyperbola) at three twist stations -> the twist
    th = np.linspace(0, 2 * math.pi, 721)
    for phi, col in ((0.0, "C0"), (30.0, "C1"), (60.0, "C2")):
        rr = np.array([_sector(math.degrees(t), phi)[0] for t in th])
        labs = [_sector(math.degrees(t), phi)[1] for t in th]
        x = rr * np.cos(th)
        y = rr * np.sin(th)
        # plot only the pole-face arcs (mask the gaps)
        xm = np.array([xi if lb != "gap" else np.nan for xi, lb in zip(x, labs)])
        ym = np.array([yi if lb != "gap" else np.nan for yi, lb in zip(y, labs)])
        ax[0].plot(xm * 1e3, ym * 1e3, col, lw=2, label=f"phi={phi:.0f} deg")
    ax[0].plot(R0 * 1e3 * np.cos(th), R0 * 1e3 * np.sin(th), "k--", lw=0.8,
               label=f"aperture r0={R0*1e3:.0f} mm")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel("y [mm]")
    ax[0].set_title("The equipotential pole face (hyperbola) twisting:\n"
                    "rotate the surface by phi(s) along the orbit")
    ax[0].legend(fontsize=8, loc="upper right")
    # RIGHT: recovered orientation tracks the prescribed twist (slope 1)
    pre = np.array(sw["prescribed_deg"])
    rec = np.array(sw["recovered_deg"])
    ax[1].plot(pre, rec, "o", color="C3", ms=6, label="recovered alpha (FEM)")
    ax[1].plot(pre, pre, "k--", lw=1, label="ideal alpha = phi (slope 1)")
    ax[1].set_xlabel("prescribed twist  phi [deg]")
    ax[1].set_ylabel("recovered quad orientation  alpha [deg]")
    ax[1].set_title(f"The n-fold twist law (n=2), MEASURED:\n"
                    f"slope {sw['tracking_slope']:.3f}, "
                    f"max err {sw['tracking_err_deg']:.2f} deg")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


if __name__ == "__main__":
    main()
