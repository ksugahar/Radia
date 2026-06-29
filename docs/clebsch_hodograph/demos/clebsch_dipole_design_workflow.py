r"""End-to-end 3-D dipole design by a Clebsch (scalar-potential) LEVEL SET carried
from the 2-D cross-section into the 3-D magnet.

The design principle is one geometric object seen at three places: the iron pole
SURFACE is a magnetic-scalar equipotential -- a level set `Omega = const` (the
Clebsch / scalar potential whose conjugate `A_z` is the flux function of the prior
flux-line work).  An accelerator dipole is designed by making that pole surface the
RIGHT level set everywhere:

  Stage A -- 2-D cross-section (the level set in the transverse plane).
      The pole face `z_p(x)` is the equipotential.  A finite flat pole droops at its
      edges (`b_3 < 0`); the level set is corrected by (1) WIDTH and (2) a CURVATURE
      shim `z_p(x) = g/2 - delta (x/w)^2` that drives `b_3` through zero -- the
      finite-aperture analogue of the quadrupole hyperbola.  This is the cheap, exact
      instrument for the TRANSVERSE harmonics.

  Stage B -- reflect the level set into the 3-D pole SURFACE.
      The magnet BODY is the 2-D level set EXTRUDED along the beam (`y`): the body
      cross-section IS the Stage-A contour, so the body field is the 2-D field
      (`B = grad(A_z) x y_hat`, the Clebsch potential extruded).  The pole END follows
      the 3-D equipotential -- the Maxwellian end -- so the fringe pseudo-multipoles
      integrate away and the INTEGRATED field stays the designed multipole.

  Stage C -- 3-D FEM verify (reduced-Omega).
      A real finite-length H-frame solve confirms a clean flat-top dipole, a clean
      INTEGRATED dipole, and reads back the equipotential END contour `z_p(y)` (the
      level set in 3-D: `g/2` in the body, lifting past the iron end = the Maxwellian
      end the design should cut).

Honest scope (repository-first): the TRANSVERSE field quality is DESIGNED and VERIFIED
in the 2-D cross-section (Stage A) -- cheap and exact.  The 3-D solve verifies the
body-field match, the clean integrated dipole, and the longitudinal/end behaviour; the
3-D INTEGRATED transverse harmonic at a golden-feasible mesh is mesh-noise-limited
(~1e-1), well above the cross-section's intrinsic `b_3` (~1e-3..1e-5), so 2-D is the
instrument for it -- exactly how accelerator dipoles are designed (2-D cross-section
optimisation + 3-D end correction).  Reflecting a NON-negligible body shim into the
3-D mesh and resolving its (sub-mesh-noise) integrated improvement is a fine-mesh
study, deferred.  The design here chooses the WIDTH so the residual shim is negligible
(`delta_opt < 0.05 mm` at the magnet width), so the flat-body 3-D solve IS the
reflected design.

Reuses: accel_pole_dipole_body_2d (Stage A 2-D solve), accel_pole_ends_fem (Stage C
reduced-Omega FEM + integrated analyzer + equipotential end contour).

run:  python clebsch_dipole_design_workflow.py          # Stage A + B (fast)
      python clebsch_dipole_design_workflow.py --fem     # + Stage C 3-D FEM (slow)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import accel_pole_dipole_body_2d as body2d       # noqa: E402  Stage A (2-D)
import accel_pole_ends_fem as fem3d               # noqa: E402  Stage C (3-D)


# ----------------------------------------------------------------------------------
# Stage A -- the 2-D cross-section level set
# ----------------------------------------------------------------------------------
def design_cross_section(half_w=fem3d.POLE_W, order=3, maxh=0.006,
                         deltas=(0.0, 0.0005, 0.0010, 0.0020, 0.0030, 0.0040)):
    """Design the cross-section pole face (the equipotential level set) at a given
    pole half-width: sweep the curvature shim delta, find delta_opt that zeros b_3."""
    rows = []
    for d in deltas:
        r = body2d._solve_field(half_w, float(d), order, maxh)
        rows.append((float(d), r["b3_over_b1"], r["b5_over_b1"], r["spurious_rel"]))
    b3 = np.array([row[1] for row in rows]); da = np.array([row[0] for row in rows])
    flat = rows[0]
    if b3[0] * b3[-1] < 0:
        idx = np.argsort(b3)
        delta_opt = float(np.interp(0.0, b3[idx], da[idx]))
        conf = body2d._solve_field(half_w, delta_opt, order, maxh)
        spur_opt, b5_opt = conf["spurious_rel"], conf["b5_over_b1"]
    else:
        delta_opt = spur_opt = b5_opt = float("nan")
    xs = np.linspace(-half_w, half_w, 41)
    zp_flat = body2d._pole_face_z(xs, half_w, 0.0)
    zp_opt = (body2d._pole_face_z(xs, half_w, delta_opt)
              if delta_opt == delta_opt else zp_flat)
    improve = (abs(flat[3]) / spur_opt) if (spur_opt == spur_opt and spur_opt > 0) else float("nan")
    return {
        "half_w": float(half_w), "g2": float(body2d.G2), "r_ref": float(body2d.R_REF),
        "sweep": rows, "b3_flat": float(flat[1]), "spur_flat": float(flat[3]),
        "delta_opt": delta_opt, "spur_opt": spur_opt, "b5_opt": b5_opt,
        "improve_factor": float(improve),
        "contour_x": xs.tolist(),
        "contour_zp_flat": [float(v) for v in zp_flat],
        "contour_zp_opt": [float(v) for v in (zp_opt if hasattr(zp_opt, "__iter__") else zp_flat)],
    }


def width_lever(widths=(0.030, 0.040, 0.050, 0.060), order=3, maxh=0.006):
    """The width design knob: flat-pole droop b_3 and the shim delta_opt needed to
    zero it, vs pole half-width.  A wider pole flattens the field -> less shim."""
    rows = []
    for hw in widths:
        d = design_cross_section(half_w=hw, order=order, maxh=maxh)
        rows.append((float(hw), d["b3_flat"], d["delta_opt"], d["spur_opt"]))
    return rows


# ----------------------------------------------------------------------------------
# Stage B -- reflect the level set into the 3-D pole surface
# ----------------------------------------------------------------------------------
def reflect_to_3d(stage_a, length=fem3d.L_BEAM):
    """Reflect the 2-D level set into the 3-D pole SURFACE spec: BODY = the Stage-A
    contour extruded along the beam (y in [-L/2, L/2]); END = the equipotential
    (Maxwellian).  Geometry consistency: the 3-D body cross-section is exactly the
    2-D contour, and the body end-contour value is g/2 (verified in Stage C)."""
    shim_negligible = (stage_a["delta_opt"] == stage_a["delta_opt"]
                       and abs(stage_a["delta_opt"]) < 0.05e-3)
    return {
        "length_m": float(length), "half_w": stage_a["half_w"], "g2": stage_a["g2"],
        "delta_opt_m": stage_a["delta_opt"],
        "body_is_extruded_2d_contour": True,           # by construction
        "shim_negligible_at_width": bool(shim_negligible),
        "end_rule": "equipotential-following (Maxwellian): pole end follows Omega=Omega_pole",
    }


# ----------------------------------------------------------------------------------
# Stage C -- 3-D reduced-Omega FEM verify
# ----------------------------------------------------------------------------------
def verify_3d(maxh_air=0.06, maxh_iron=0.03, n_beam=61, n_theta=24, r_ref=0.008):
    """3-D reduced-Omega FEM at the designed width: clean flat-top dipole + clean
    integrated dipole + the equipotential END contour z_p(y) (the level set in 3-D)."""
    return fem3d.solve(maxh_air=maxh_air, maxh_iron=maxh_iron, n_beam=n_beam,
                       n_theta=n_theta, r_ref=r_ref, return_contour=True)


# ----------------------------------------------------------------------------------
# the whole workflow
# ----------------------------------------------------------------------------------
def run_workflow(design_half_w=fem3d.POLE_W, with_fem=False, **fem_kw):
    """Stage A (design at the magnet width) + the width lever + Stage B (reflect) +
    optionally Stage C (3-D FEM verify)."""
    a = design_cross_section(half_w=design_half_w)
    lever = width_lever()
    b = reflect_to_3d(a)
    out = {"stage_a": a, "width_lever": lever, "stage_b": b}
    if with_fem:
        out["stage_c"] = verify_3d(**fem_kw)
    return out


def _plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.6, 4.2), dpi=150)

    # Panel A: the 2-D cross-section level set (showcase a narrow width so the shim shows)
    show = design_cross_section(half_w=0.030)
    xs = np.array(show["contour_x"]) * 1e3
    axA.plot(xs, np.array(show["contour_zp_flat"]) * 1e3, "C7--", lw=1.2, label="flat pole (droops)")
    axA.plot(xs, np.array(show["contour_zp_opt"]) * 1e3, "C0-", lw=1.8,
             label=f"designed level set ($\\delta$={show['delta_opt']*1e3:.2f} mm)")
    axA.axhline(show["g2"] * 1e3, color="0.6", lw=0.6)
    axA.fill_between(xs, np.array(show["contour_zp_opt"]) * 1e3, (show["g2"] + 0.006) * 1e3,
                     color="0.85", zorder=0)
    axA.set_xlabel("x  [mm]"); axA.set_ylabel("pole face  $z_p(x)$  [mm]")
    axA.set_title("Stage A: the cross-section level set\n(pole face = scalar equipotential)")
    axA.legend(fontsize=8, loc="lower center")

    # Panel B: the width lever -- flat b_3 droop and shim delta_opt vs half-width
    lev = out["width_lever"]
    w = [row[0] * 1e3 for row in lev]
    b3 = [abs(row[1]) for row in lev]
    dopt = [row[2] * 1e3 for row in lev]
    axB.semilogy(w, b3, "C3o-", label="|$b_3/b_1$| flat (droop)")
    axB2 = axB.twinx()
    axB2.plot(w, dopt, "C0s--", label="shim $\\delta_{opt}$ [mm]")
    axB.set_xlabel("pole half-width  w  [mm]")
    axB.set_ylabel("|$b_3/b_1$| flat pole", color="C3")
    axB2.set_ylabel("$\\delta_{opt}$  [mm]  (shim to zero $b_3$)", color="C0")
    axB.set_title("the design knob: a wider pole\nneeds less shim (flatter field)")
    axB.axvline(out["stage_a"]["half_w"] * 1e3, color="0.5", ls=":", lw=1)

    # Panel C: the level set carried into 3-D (end contour) + the verified field
    sc = out.get("stage_c")
    if sc is not None and "contour_yc" in sc:
        yc = np.array(sc["contour_yc"]) * 1e3
        zp = np.array(sc["contour_zpole"]) * 1e3
        axC.plot(yc, zp, "C0-", lw=1.8, label="end contour $z_p(y)$ (level set)")
        axC.axhline(sc["z_pole_body_m"] * 1e3, color="0.6", lw=0.6, ls="--")
        axC.axvline(fem3d.L_BEAM / 2 * 1e3, color="0.5", lw=0.8, ls=":")
        axC.set_xlabel("beam  y  [mm]"); axC.set_ylabel("equipotential  $z_p(y)$  [mm]")
        axC.set_title(f"Stage C: the level set in 3-D\n(body $g/2$, lifts {sc['end_contour_lift_m']*1e3:.1f} mm "
                      f"past the end)")
        axC.legend(fontsize=8, loc="upper left")
    else:
        axC.text(0.5, 0.5, "Stage C (3-D FEM):\nrun with --fem", ha="center", va="center",
                 transform=axC.transAxes, fontsize=11)
        axC.set_xticks([]); axC.set_yticks([])
        axC.set_title("Stage C: the level set in 3-D")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    with_fem = "--fem" in sys.argv
    print("Clebsch level-set dipole design workflow: 2-D cross-section -> 3-D magnet\n")
    out = run_workflow(with_fem=with_fem)
    a, b = out["stage_a"], out["stage_b"]
    print(f"  Stage A -- 2-D cross-section level set (half-width {a['half_w']*1e3:.0f} mm, "
          f"aperture r_ref {a['r_ref']*1e3:.0f} mm):")
    print(f"    flat pole droops: b_3/b_1 = {a['b3_flat']:+.3e}  (finite-pole edge droop)")
    print(f"    shim delta_opt = {a['delta_opt']*1e3:.3f} mm zeroes b_3: spurious "
          f"{abs(a['spur_flat']):.2e} -> {a['spur_opt']:.2e}  ({a['improve_factor']:.1f}x, "
          f"residual |b_5/b_1| = {abs(a['b5_opt']):.2e})")
    print(f"  the width design knob (flat-pole droop and shim vs half-width):")
    for hw, b3f, dopt, sp in out["width_lever"]:
        print(f"    w = {hw*1e3:4.0f} mm:  flat b_3/b_1 = {b3f:+.3e}  delta_opt = {dopt*1e3:5.3f} mm")
    print(f"  Stage B -- reflect the level set into the 3-D pole surface:")
    print(f"    body = the 2-D contour EXTRUDED along the beam (body field = the 2-D field)")
    print(f"    end  = {b['end_rule']}")
    print(f"    shim negligible at the chosen width: {b['shim_negligible_at_width']} "
          f"(delta_opt = {b['delta_opt_m']*1e3:.3f} mm) -> the flat-body 3-D IS the reflected design")
    if with_fem:
        sc = out["stage_c"]
        print(f"  Stage C -- 3-D reduced-Omega FEM verify (ne={sc['ne']}, ndof={sc['ndof']}):")
        print(f"    flat-top body field B_z = {sc['bz_body_T']*1e3:.1f} mT, "
              f"B_x/B_z = {sc['bx_over_bz_centre']:.2e} (x-symmetric -> clean dipole)")
        print(f"    integrated dipole bbar_1 = {sc['integrated_dipole_bbar1_Tm']*1e3:.2f} mT*m, "
              f"L_eff = {sc['L_eff_m']*1e3:.1f} mm (> L_iron by the fringes)")
        print(f"    equipotential END contour: body z_p = {sc['z_pole_body_m']*1e3:.1f} mm (= g/2), "
              f"lifts {sc['end_contour_lift_m']*1e3:.1f} mm past the iron end (the Maxwellian end)")
        print(f"    integrated transverse spurious = {sc['integrated_spurious_rel']:.2e} "
              f"(mesh-noise-limited; the transverse quality is the 2-D instrument's job)")
    else:
        print("  Stage C -- 3-D FEM verify: run with  --fem  (slow)")
    print("\n  => the pole SURFACE is one Clebsch/scalar level set carried 2-D cross-section ->")
    print("     3-D body (extrude) -> 3-D end (Maxwellian).  2-D designs the transverse harmonics")
    print("     (cheap, exact); 3-D verifies the body field, the clean integrated dipole, and the")
    print("     equipotential end.")
    _plot(out)


if __name__ == "__main__":
    main()
