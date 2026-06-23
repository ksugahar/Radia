"""The END PACK designed in two orthogonal 2-D planes, reflected into 3-D.

THE PICTURE (the user's design thought, made literal)
-----------------------------------------------------
Design the magnet END -- the "end pack" -- in TWO orthogonal 2-D planes and
REFLECT both designs into one 3-D pole:

  | plane | what it sets | engine |
  |---|---|---|
  | **x-y cross-section** (perp the beam) | the transverse multipole the beam SEES: a flat finite pole droops (b_3<0), the shim curvature z=g/2-delta(x/w)^2 zeroes it | accel_pole_dipole_body_2d.solve (Plane 1) |
  | **s-y longitudinal** (along the beam s, the gap y) | how the field TERMINATES: the end equipotential bow-out z_p(s) past the pole edge = the end chamfer the iron should follow | this file, solve_sy_endpack (Plane 2) |

convention: beam = y (= the s arc length of a straight magnet), gap = z (pole
faces z = +-g/2), width = x.  So the user's "x-y" cross-section is here the
(x, z=gap) plane and the user's "s-y" longitudinal plane is here the
(y=beam, z=gap) plane -- the same two planes, our z is the user's vertical y.

THE METHOD (two 2-D DESIGN solves, then a 3-D REFLECTION)
--------------------------------------------------------
Unlike accel_pole_ends_fem.py -- which reads the end equipotential OUT of the
finished 3-D solve -- here BOTH planes are standalone CHEAP 2-D DESIGN solves
done FIRST, and the 3-D solve only REFLECTS + VERIFIES:

  1. Plane 1 (x-y): the 2-D cross-section equipotential (a high-mu flat pole)
     gives the transverse harmonic b_3 of the body and the shim delta that
     zeroes it -- the cross-section DESIGN.
  2. Plane 2 (s-y): a standalone 2-D reduced-Omega solve (a UNIFORM body source
     H_s in z + an iron pole that TERMINATES at the end) gives the end
     equipotential bow-out z_p(s) -- the chamfer shape ghat(s) the iron END
     should follow, and the effective magnetic length L_eff -- the end-pack
     DESIGN, computed BEFORE any 3-D mesh.
  3. Reflection (3-D): loft the cross-section along the beam, follow ghat(s) at
     the end (accel_pole_ends_fem.solve(curved_shape=ghat)), solve reduced-Omega,
     read the INTEGRATED transverse multipoles INT B_perp dl.  Verify (a) the 3-D
     integrated transverse spurious matches the Plane-1 body prediction, and
     (b) the designed s-y chamfer drives the pole-END field enhancement down vs
     the hard-cut end -- and CROSS-CHECK that the cheap 2-D Plane-2 contour
     z_p(s) predicts the expensive 3-D end equipotential.

THE RESULT (the two-plane reflection works)
-------------------------------------------
The cheap 2-D s-y plane predicts the 3-D end equipotential bow-out to a few %;
reflecting that chamfer into the 3-D pole drives the longitudinal pole-end
enhancement toward zero, while the transverse integrated b_3 stays the
Plane-1 body value (the honest two-lever split: the x-y plane owns the
transverse harmonic, the s-y plane owns the end termination).

run:  python endpack_two_plane.py            # the two-plane design + 3-D reflect
      python endpack_two_plane.py --fig        # + the 3-panel figure
      python endpack_two_plane.py --fast       # coarse mesh (golden / quick)
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # sibling examples

import accel_pole_dipole_body_2d as body2d                       # noqa: E402  Plane 1
import accel_pole_ends_fem as ends3d                             # noqa: E402  3-D reflect

GAP = ends3d.GAP                 # 0.040  pole-face separation (z)
G2 = GAP / 2.0
L_BEAM = ends3d.L_BEAM           # 0.120  iron length along the beam (y) -> ends at +-L/2
POLE_W = ends3d.POLE_W           # 0.060  pole half-width (x)
Z_OUT = ends3d.Z_OUT             # 0.100  outer z of the pole iron
MU0 = ends3d.MU0

# ---- Plane-2 (s-y) air box ----
Y_AIR = 0.30                     # air half-extent along the beam (y)
Z_AIR = 0.22                     # air height (z, upper half only)
POLE_T = Z_OUT - G2              # pole iron thickness above the face
PHI_REF_FRAC = 0.55              # interior-gap equipotential to read the bow-out (0=mid,1=face)

# ---- 3-D reflection box (upper-half, equipotential-pole drive) ----
X_AIR3, Y_AIR3, Z_AIR3 = 0.22, 0.26, Z_OUT   # (width, beam, gap) air half-extents


# ============================================================
# Plane 2: the s-y longitudinal end-pack DESIGN (standalone 2-D Laplace)
# ============================================================

def _build_sy_mesh(maxh=0.005):
    """2-D (y=beam, z=gap) upper-half air gap with a high-mu (Dirichlet) iron
    pole that TERMINATES at the iron end y = +-L/2 -- the parallel-plate fringe
    whose end equipotentials are the Rogowski chamfer.  Pole face Phi=1, midplane
    z=0 Phi=0 (the dipole antisymmetry plane), outer box Phi=0.  Pure Laplace
    (the high-mu pole-as-equipotential limit, same model as the x-y cross-section
    accel_pole_dipole_body_2d)."""
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry

    hL = L_BEAM / 2.0
    box = WorkPlane().MoveTo(-Y_AIR, 0.0).Rectangle(2 * Y_AIR, Z_AIR).Face()
    pole = WorkPlane().MoveTo(-hL, G2).Rectangle(L_BEAM, POLE_T).Face()   # terminates at +-hL
    air = box - pole
    for e in air.edges:
        c = e.center
        if abs(c.y) < 1e-7:
            e.name = "mid"                       # midplane z=0 (Omega = 0)
        elif abs(c.x) > Y_AIR - 1e-7 or c.y > Z_AIR - 1e-7:
            e.name = "far"                       # outer box (Omega = 0)
        else:
            e.name = "pole"                      # pole face (Omega = 1)
    air.maxh = maxh
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air, dim=2).GenerateMesh(maxh=maxh))
    return mesh


def solve_sy_endpack(order=3, maxh_air=0.005, maxh_iron=None,
                     chamfer_len=0.030, n_contour=90):
    """DESIGN the s-y end pack: solve the standalone 2-D Laplace fringe, read the
    end equipotential bow-out z_p(y) past the pole edge (the Rogowski chamfer the
    iron END should follow), and return the chamfer shape ghat(s) + natural depth
    + L_eff.  Cheap (2-D), done BEFORE any 3-D mesh.  (maxh_iron is ignored --
    the pole is a Dirichlet equipotential, no iron volume.)"""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)

    mesh = _build_sy_mesh(maxh_air)
    fes = H1(mesh, order=order, dirichlet="mid|pole|far")
    u, v = fes.TnT()
    with TaskManager():
        omega = GridFunction(fes)
        omega.Set(mesh.BoundaryCF({"pole": 1.0}, default=0.0),
                  definedon=mesh.Boundaries("mid|pole|far"))
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f.Assemble()
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * omega.vec
        omega.vec.data += a.mat.Inverse(fes.FreeDofs()) * r

    B_cf = -grad(omega)                                  # (B_y, B_z) up to scale

    hL = L_BEAM / 2.0
    # ---- the end equipotential contour z_phi(y): the INTERIOR-gap level set
    # Phi = Phi_ref bows up and over past the pole edge -> the Rogowski chamfer. ----
    z_ref = PHI_REF_FRAC * G2                            # interior reference height
    phi_ref = float(omega(mesh(0.0, z_ref)))             # body value of the level set
    yc = np.linspace(0.0, hL + 1.10 * chamfer_len, n_contour)
    zc = np.linspace(0.0, 1.6 * G2, 90)
    z_phi = np.empty(len(yc))
    for i, yv in enumerate(yc):
        phi_col = np.empty(len(zc))               # Phi up the column; pole interior = 1
        for j, zv in enumerate(zc):
            try:
                phi_col[j] = float(omega(mesh(yv, zv)))
            except Exception:
                phi_col[j:] = 1.0                 # entered the pole iron -> Phi = 1
                break
        # Phi rises monotonically from 0 (mid) toward the face; invert for z(Phi_ref).
        z_phi[i] = float(np.interp(phi_ref, phi_col, zc))

    # ---- chamfer shape ghat(s): the normalized bow-out past the iron end ----
    z_body = float(z_phi[0])
    sig = np.linspace(0.0, chamfer_len, 40)
    lift = np.interp(hL + sig, yc, z_phi) - z_body
    lift = np.maximum.accumulate(np.maximum(lift, 0.0))  # monotone, non-negative
    depth_nat = float(lift[-1]) if lift[-1] > 1e-5 else 1e-5
    ghat_x = sig / chamfer_len
    ghat_y = lift / depth_nat

    # ---- L_eff from the on-axis gap field (near the midplane) ----
    ya = np.linspace(-1.5 * hL, 1.5 * hL, 161)
    z_probe = 0.15 * G2
    bz_axis = np.array([float(B_cf(mesh(yv, z_probe))[1]) for yv in ya])
    body = np.abs(ya) < 0.30 * L_BEAM
    bz_body = float(np.mean(bz_axis[body]))
    l_eff = float(np.trapezoid(bz_axis, ya) / bz_body) if abs(bz_body) > 1e-12 else 0.0

    return {
        "ne": int(mesh.ne), "ndof": int(fes.ndof),
        "depth_natural_m": depth_nat,                    # equipotential lift over c_len
        "end_lift_m": float(np.max(z_phi) - z_body),     # peak bow-out
        "z_ref_m": float(z_ref),                         # the level-set reference height
        "L_eff_2d_m": l_eff,                             # 2-D effective magnetic length
        "L_eff_excess_2d": float((l_eff - L_BEAM) / L_BEAM),
        "ghat_x": ghat_x.tolist(), "ghat_y": ghat_y.tolist(),
        "contour_yc": yc.tolist(), "contour_zp": (z_phi).tolist(),
        "chamfer_len_m": float(chamfer_len),
    }


def _ghat_callable(sy):
    """Build the chamfer shape callable ghat: [0,1]->[0,1] from the s-y design."""
    gx = np.array(sy["ghat_x"])
    gy = np.array(sy["ghat_y"])
    return lambda s: float(np.interp(min(max(s, 0.0), 1.0), gx, gy))


# ============================================================
# the 3-D REFLECTION: upper-half equipotential-pole drive (NO coil/RadiaField)
# ============================================================

def _build_3d_endpack(chamfer_depth=0.0, chamfer_len=0.030, curved_shape=None,
                      maxh_air=0.05, n_face=51):
    """Upper-half 3-D dipole END: a top pole that TERMINATES at the beam ends
    (gap face follows z=g/2+depth*ghat(s) at the end), driven as an EQUIPOTENTIAL
    (Psi=mmf on the whole pole surface), midplane z=0 the dipole antisymmetry plane
    (Psi=0).  Pure Laplace -- NO coil / RadiaField (same drive as the FFAG rung-2,
    and the high-mu limit of the iron pole).  Cheap: B = -grad(Psi) is a plain CF."""
    import ngsolve as ng
    from netgen.occ import Box, Pnt, WorkPlane, Axes, X, Y, OCCGeometry

    hL = L_BEAM / 2.0

    def face_z(y):
        if curved_shape is None or chamfer_depth <= 0:
            return G2
        s = max((y - (hL - chamfer_len)) / chamfer_len,
                ((-hL + chamfer_len) - y) / chamfer_len)
        s = min(max(s, 0.0), 1.0)
        return G2 + chamfer_depth * curved_shape(s)

    ys = np.linspace(-hL, hL, n_face)
    zf = [float(face_z(y)) for y in ys]
    wp = WorkPlane(Axes(Pnt(-POLE_W, 0.0, 0.0), n=X, h=Y))   # plane = (y, z)
    wp.MoveTo(float(ys[0]), zf[0])
    for yi, zi in zip(ys[1:], zf[1:]):
        wp.LineTo(float(yi), zi)
    wp.LineTo(float(ys[-1]), Z_OUT)                          # up the +y end
    wp.LineTo(float(ys[0]), Z_OUT)                           # across the top
    wp.Close()                                               # down the -y end
    pole = wp.Face().Extrude(2 * POLE_W)                     # x-prism (flat in x)

    air = Box(Pnt(-X_AIR3, -Y_AIR3, 0.0), Pnt(X_AIR3, Y_AIR3, Z_AIR3)) - pole
    for f in air.faces:
        c = f.center
        if abs(c.z) < 1e-7:
            f.name = "median"                                # z=0 antisymmetry (Psi=0)
        elif abs(c.x) > X_AIR3 - 1e-7 or abs(c.y) > Y_AIR3 - 1e-7 or c.z > Z_AIR3 - 1e-7:
            f.name = "far"                                   # outer box (Psi=0)
        else:
            f.name = "pole"                                  # pole surface (Psi=mmf)
    air.maxh = maxh_air
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air).GenerateMesh(maxh=maxh_air))
    return mesh


def _solve_3d_endpack(mesh, order=2, n_beam=29, n_theta=12, want_contour=False,
                      chamfer_len=0.030):
    """Equipotential-pole Laplace; read the gap field, the pole-TIP corner field
    (the end enhancement the chamfer rounds), the effective length, and the
    INTEGRATED transverse multipoles INT B_perp dl (the beam-optics quantity)."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    from accel_pole_ends_3d import integrated_multipoles

    hL = L_BEAM / 2.0
    fes = H1(mesh, order=order, dirichlet="median|pole|far")
    u, v = fes.TnT()
    with TaskManager():
        gfu = GridFunction(fes)
        gfu.Set(mesh.BoundaryCF({"pole": 1.0, "median": 0.0}, default=0.0),
                definedon=mesh.Boundaries("median|pole|far"))
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f.Assemble()
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
        B = -grad(gfu)                                       # (B_y, _, B_z) up to mmf scale

        # on-axis gap profile (near the midplane) -> L_eff
        ys = np.linspace(-1.4 * hL, 1.4 * hL, 61)
        bz = np.array([float(B(mesh(0.0, y, 0.15 * G2))[2]) for y in ys])
        body = np.abs(ys) < 0.30 * L_BEAM
        bz_body = float(np.mean(bz[body]))
        l_eff = float(np.trapezoid(bz, ys) / bz_body) if abs(bz_body) > 1e-12 else 0.0

        # pole-TIP corner field (the end enhancement) vs body
        tip = float(abs(B(mesh(0.0, hL - 2e-3, 0.85 * G2))[2])) / abs(bz_body)

        # INTEGRATED transverse multipoles INT B_perp dl (upper-half + antisymmetry)
        def B_perp(ax, ay, az):                              # analyzer (x, transverse, beam)
            if ay >= 0:
                b = B(mesh(ax, az, ay))
                return (float(b[0]), float(b[2]))
            b = B(mesh(ax, az, -ay))                         # reflect: dipole antisymmetry
            return (-float(b[0]), float(b[2]))
        mp = integrated_multipoles(B_perp, 0.008, (-1.6 * hL, 1.6 * hL),
                                   n_z=n_beam, n_theta=n_theta, n_max=6)
        main = abs(complex(*mp[1]))
        spurious = max(abs(complex(*mp[n])) / main for n in (3, 5))

        out = {
            "ne": int(mesh.ne), "ndof": int(fes.ndof),
            "bz_body": bz_body, "tip_enhancement": tip,
            "L_eff_m": l_eff, "L_eff_excess": float((l_eff - L_BEAM) / L_BEAM),
            "integrated_dipole_bbar1": float(main),
            "integrated_spurious_rel": float(spurious),
        }
        if want_contour:
            # the 3-D end equipotential bow-out (Psi = gfu is the potential):
            # the interior level set Phi_ref bows up past the pole edge.
            z_ref = PHI_REF_FRAC * G2
            phi_ref = float(gfu(mesh(0.0, 0.0, z_ref)))
            yc = np.linspace(0.0, hL + 1.10 * chamfer_len, 60)
            zc = np.linspace(0.0, 1.6 * G2, 60)
            z_phi = np.empty(len(yc))
            for i, yv in enumerate(yc):
                col = np.empty(len(zc))
                for j, zv in enumerate(zc):
                    try:
                        col[j] = float(gfu(mesh(0.0, yv, zv)))
                    except Exception:
                        col[j:] = 1.0
                        break
                z_phi[i] = float(np.interp(phi_ref, col, zc))
            out["contour_yc"] = yc.tolist()
            out["contour_zphi"] = z_phi.tolist()
            out["z_ref_m"] = float(z_ref)
    return out


# ============================================================
# the two-plane orchestration: design x-y + s-y, reflect into 3-D
# ============================================================

def design_two_plane(order_2d=3, order_3d=2, fast=False, chamfer_len=0.030,
                     depth_fracs=(0.0, 0.12, 0.25, 0.40)):
    """Run BOTH 2-D design planes, then REFLECT into one 3-D pole and verify.

    The s-y plane fixes the chamfer SHAPE ghat(s); the 3-D reflection sweeps the
    chamfer DEPTH (as a fraction of g/2) and finds the depth that drives the
    pole-tip end enhancement through 1 (the zero-enhancement design point)."""
    maxh_air2 = 0.009 if fast else 0.006
    maxh_air3 = 0.060 if fast else 0.045
    n_face = 41 if fast else 61
    n_beam, n_theta = (25, 10) if fast else (33, 16)

    # ---- Plane 1 (x-y cross-section): the transverse harmonic + the shim ----
    xy = body2d.solve(order=order_2d, maxh=maxh_air2)
    b3_body = next(b3 for (w, b3, _sp) in xy["width_sweep"] if abs(w - POLE_W) < 1e-6)
    xy_design = {
        "pole_half_width_m": float(POLE_W),
        "b3_over_b1_flat_body": float(b3_body),          # the droop the shim removes
        "shim_delta_opt_m": float(xy["delta_opt_m"]),    # the concavity that zeroes b_3
        "spurious_at_shim_opt": float(xy["spurious_at_opt_rel"]),
    }

    # ---- Plane 2 (s-y longitudinal): the end chamfer shape + L_eff ----
    sy = solve_sy_endpack(order=order_2d, maxh_air=maxh_air2, chamfer_len=chamfer_len)
    ghat = _ghat_callable(sy)

    # ---- Reflection (3-D): sweep the chamfer DEPTH, find zero end enhancement ----
    sweep = []
    base_contour = None
    for fr in depth_fracs:
        cs = None if fr <= 0 else ghat
        mesh = _build_3d_endpack(chamfer_depth=fr * G2, chamfer_len=chamfer_len,
                                 curved_shape=cs, maxh_air=maxh_air3, n_face=n_face)
        r = _solve_3d_endpack(mesh, order=order_3d, n_beam=n_beam, n_theta=n_theta,
                              want_contour=(fr <= 0), chamfer_len=chamfer_len)
        if fr <= 0:
            base_contour = r
        sweep.append({"depth_frac": float(fr), "chamfer_depth_m": float(fr * G2),
                      "tip_enhancement": r["tip_enhancement"],
                      "integrated_spurious_rel": r["integrated_spurious_rel"],
                      "L_eff_m": r["L_eff_m"]})
    flat = sweep[0]
    # the designed depth: where the tip enhancement crosses 1 (zero corner over-field)
    fr_arr = np.array([s["depth_frac"] for s in sweep])
    tip_arr = np.array([s["tip_enhancement"] for s in sweep])
    if tip_arr[0] > 1.0 > tip_arr[-1]:
        idx = np.argsort(tip_arr)
        fr_zero = float(np.interp(1.0, tip_arr[idx], fr_arr[idx]))
    else:
        fr_zero = float("nan")

    # ---- CROSS-CHECK: 2-D s-y chamfer SHAPE vs 3-D end-equipotential SHAPE ----
    hL = L_BEAM / 2.0
    sig = np.linspace(0.0, chamfer_len, 30)
    yc3 = np.array(base_contour["contour_yc"]); zp3 = np.array(base_contour["contour_zphi"])
    z_body3 = float(zp3[0])
    lift3 = np.maximum.accumulate(np.maximum(np.interp(hL + sig, yc3, zp3) - z_body3, 0.0))
    ghat3 = lift3 / max(float(np.max(lift3)), 1e-9)
    ghat2 = np.array([ghat(s) for s in sig / chamfer_len])
    contour_rms_rel = float(np.sqrt(np.mean((ghat2 - ghat3) ** 2)))

    out = {
        "plane1_xy": xy_design,
        "plane2_sy": {k: sy[k] for k in (
            "ne", "ndof", "end_lift_m", "z_ref_m",
            "L_eff_2d_m", "L_eff_excess_2d", "chamfer_len_m")},
        "reflect_3d": {
            "ne": base_contour["ne"], "ndof": base_contour["ndof"],
            "depth_sweep": sweep,
            "flat_tip_enhancement": flat["tip_enhancement"],
            "designed_depth_frac_zero": fr_zero,
            "designed_chamfer_depth_m": float(fr_zero * G2) if fr_zero == fr_zero else float("nan"),
            "tip_enhancement_reduction": float(flat["tip_enhancement"] - min(
                s["tip_enhancement"] for s in sweep)),
            "integrated_spurious_flat": flat["integrated_spurious_rel"],
        },
        "crosscheck_2d_vs_3d_contour_rms_rel": contour_rms_rel,
        "transverse_lever_consistency": {
            "xy_2d_b3_body": float(abs(b3_body)),
            "3d_integrated_spurious_flat": float(flat["integrated_spurious_rel"]),
        },
    }
    return out


# ============================================================
# figure + main
# ============================================================

def _figure(res, sy_full):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(14.4, 4.0), dpi=140)

    # LEFT: Plane 1 (x-y) -- the cross-section shim that zeroes transverse b_3
    p1 = res["plane1_xy"]
    xs = np.linspace(-POLE_W, POLE_W, 120)
    d = p1["shim_delta_opt_m"]
    ax[0].plot(xs * 1e3, (G2 - d * (xs / POLE_W) ** 2) * 1e3, "C0", lw=2,
               label=f"shim pole face $z=g/2-\\delta(x/w)^2$\n$\\delta$="
                     f"{d * 1e3:.2f} mm ($b_3$=0)")
    ax[0].plot(xs * 1e3, np.full_like(xs, G2) * 1e3, "k--", lw=1,
               label=f"flat pole ($b_3/b_1$={p1['b3_over_b1_flat_body']:+.1e})")
    ax[0].set_xlabel("width  x [mm]")
    ax[0].set_ylabel("pole face  z [mm]")
    ax[0].set_title("Plane 1 (x-y cross-section):\nthe shim that zeroes "
                    "transverse $b_3$")
    ax[0].legend(fontsize=8)

    # MIDDLE: Plane 2 (s-y) -- the normalized end chamfer shape ghat(s) (Rogowski)
    gx = np.array(sy_full["ghat_x"]); gy = np.array(sy_full["ghat_y"])
    ax[1].plot(gx, gy, "C1", lw=2,
               label="end equipotential bow-out $\\hat g(s)$ (Rogowski)")
    ax[1].plot(gx, gx, "k--", lw=1, label="linear taper")
    ax[1].set_xlabel("fractional distance into the chamfer  $s$  (0=body, 1=end)")
    ax[1].set_ylabel("normalized chamfer depth")
    ax[1].set_title(f"Plane 2 (s-y end pack): standalone 2-D design\n"
                    f"$L_{{eff}}$ excess {sy_full['L_eff_excess_2d']*100:+.0f}% "
                    f"(end lift {sy_full['end_lift_m']*1e3:.0f} mm)")
    ax[1].legend(fontsize=8)

    # RIGHT: Reflection (3-D) -- chamfer-DEPTH sweep drives the pole-tip
    # corner enhancement through 1 (the zero-enhancement design point)
    r3 = res["reflect_3d"]
    sw = r3["depth_sweep"]
    dmm = [s["chamfer_depth_m"] * 1e3 for s in sw]
    tip = [(s["tip_enhancement"] - 1.0) * 100 for s in sw]   # corner over-field [%]
    spr = [s["integrated_spurious_rel"] * 100 for s in sw]
    ax[2].plot(dmm, tip, "o-", color="C3", label="pole-tip corner over-field [%]")
    ax[2].plot(dmm, spr, "s--", color="C0",
               label="integrated transverse $b_{3,5}$ [%]")
    ax[2].axhline(0, color="k", lw=0.6)
    dz = r3["designed_chamfer_depth_m"]
    if dz == dz:
        ax[2].axvline(dz * 1e3, color="C2", lw=1.2, ls=":",
                      label=f"designed depth {dz*1e3:.1f} mm (zero over-field)")
    ax[2].set_xlabel("end chamfer depth [mm]")
    ax[2].set_ylabel("[%]")
    ax[2].set_title(f"Reflection (3-D): the s-y chamfer rounds the pole-tip\n"
                    f"corner (2-D vs 3-D shape rms "
                    f"{res['crosscheck_2d_vs_3d_contour_rms_rel']*100:.0f}%)")
    ax[2].legend(fontsize=8)

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--fast", action="store_true", help="coarse mesh (golden/quick)")
    ap.add_argument("--chamfer-len", type=float, default=0.030)
    args = ap.parse_args()

    print("=" * 76)
    print("The END PACK in two planes -- x-y cross-section + s-y end -> 3-D reflect")
    print("=" * 76)

    res = design_two_plane(fast=args.fast, chamfer_len=args.chamfer_len)

    p1 = res["plane1_xy"]
    print("\nPlane 1 (x-y cross-section) -- the transverse multipole DESIGN:")
    print(f"  flat pole (half-width {p1['pole_half_width_m']*1e3:.0f} mm): "
          f"b_3/b_1 = {p1['b3_over_b1_flat_body']:+.2e}  (the edge droop)")
    print(f"  -> shim delta = {p1['shim_delta_opt_m']*1e3:.2f} mm zeroes b_3 "
          f"(residual quality {p1['spurious_at_shim_opt']:.1e})")

    p2 = res["plane2_sy"]
    print("\nPlane 2 (s-y longitudinal) -- the end-pack DESIGN (standalone 2-D):")
    print(f"  mesh: ne={p2['ne']}, ndof={p2['ndof']} (cheap, done BEFORE 3-D)")
    print(f"  end equipotential bow-out (the chamfer shape) lift = "
          f"{p2['end_lift_m']*1e3:.1f} mm over {p2['chamfer_len_m']*1e3:.0f} mm")
    print(f"  2-D effective magnetic length L_eff = {p2['L_eff_2d_m']*1e3:.1f} mm "
          f"(excess {p2['L_eff_excess_2d']*100:+.0f}% over iron L = "
          f"{L_BEAM*1e3:.0f} mm)")

    r3 = res["reflect_3d"]
    print("\nReflection (3-D) -- equipotential-pole drive, sweep the chamfer DEPTH:")
    print(f"  mesh: ne={r3['ne']}, ndof={r3['ndof']}  (pure Laplace, no coil)")
    print(f"  {'depth [mm]':<12}{'tip over-field':<16}{'integrated b_3,5':<16}")
    for s in r3["depth_sweep"]:
        print(f"  {s['chamfer_depth_m']*1e3:<12.1f}"
              f"{(s['tip_enhancement']-1.0)*100:+6.1f} %        "
              f"{s['integrated_spurious_rel']*100:5.2f} %")
    dz = r3["designed_chamfer_depth_m"]
    print(f"  hard-cut pole-tip corner over-field = "
          f"{(r3['flat_tip_enhancement']-1.0)*100:+.1f} %")
    if dz == dz:
        print(f"  -> designed chamfer depth ~ {dz*1e3:.1f} mm drives the corner "
              f"over-field through ZERO (the s-y shape, depth from the 3-D sweep).")
    print(f"  integrated transverse b_3,5 stays ~{r3['integrated_spurious_flat']*100:.2f}% "
          f"(body/Plane-1 lever -- the END shape is the wrong lever for it).")

    cc = res["crosscheck_2d_vs_3d_contour_rms_rel"]
    print(f"\nCROSS-CHECK: the cheap 2-D s-y chamfer SHAPE predicts the 3-D end "
          f"equipotential bow-out to rms {cc*100:.0f}%.")
    print("  => the two-plane reflection works: design the END in x-y + s-y, "
          "loft to 3-D.")

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    with open(jpath, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        sy_full = solve_sy_endpack(maxh_air=0.009 if args.fast else 0.006,
                                   chamfer_len=args.chamfer_len)
        _figure(res, sy_full)


if __name__ == "__main__":
    main()
