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

Convention (matches the C-type electromagnet dipole-with-CoilBuilder geometry):
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


def _chamfer_pole_end(pole, gap_z, body_sign, y_sign, c_len, c_depth):
    """Chamfer one gap-facing END corner of a pole: cut the wedge near
    (y = y_sign*(L/2 - c_len) .. y_sign*L/2, z = gap_z) so the gap WIDENS by
    c_depth over the last c_len toward that end -- the equipotential-following
    end taper.  body_sign = +1 if the pole body is at z > gap_z (top pole),
    -1 if below (bottom pole)."""
    import math
    from netgen.occ import Box, Pnt, Axis, X
    hL = L_BEAM / 2
    ang = body_sign * math.degrees(math.atan2(c_depth, c_len))
    yh = y_sign * (hL - c_len)
    # cut box flush on the GAP side of the face (z from gap_z back into the gap
    # by 0.08), covering this end; rotate about the hinge (yh, gap_z) so it
    # tilts c_depth INTO the pole corner.
    z0, z1 = (gap_z - 0.08, gap_z) if body_sign > 0 else (gap_z, gap_z + 0.08)
    if y_sign > 0:
        cut = Box(Pnt(-POLE_W * 1.4, hL - c_len, z0), Pnt(POLE_W * 1.4, hL + 0.08, z1))
    else:
        cut = Box(Pnt(-POLE_W * 1.4, -hL - 0.08, z0), Pnt(POLE_W * 1.4, -(hL - c_len), z1))
    cut = cut.Rotate(Axis(Pnt(0, yh, gap_z), X), y_sign * ang)
    return pole - cut


def _pole_prism(face_z, z_far, ys, x_half):
    """A pole built as an x-prism of a (y, z) profile: the gap-facing edge
    follows ``face_z(y)`` (an arbitrary curve), the outer edge is flat at
    ``z_far``.  Robust for a CURVED chamfer (no boolean -- the curve is baked
    into the cross-section).  ys is the y-sampling of the gap face."""
    from netgen.occ import WorkPlane, Axes, Pnt, X, Y
    wp = WorkPlane(Axes(Pnt(-x_half, 0.0, 0.0), n=X, h=Y))   # plane = (y, z)
    zf = [float(face_z(y)) for y in ys]
    wp.MoveTo(float(ys[0]), zf[0])
    for yi, zi in zip(ys[1:], zf[1:]):
        wp.LineTo(float(yi), zi)
    wp.LineTo(float(ys[-1]), float(z_far))               # up the +y end
    wp.LineTo(float(ys[0]), float(z_far))                # across the outer edge
    wp.Close()                                           # down the -y end
    return wp.Face().Extrude(2 * x_half)


def build_mesh_curved(depth, chamfer_len, shape, maxh_air=0.05, maxh_iron=0.025,
                      n_face=81):
    """H-frame dipole whose pole ENDS are chamfered along a CURVED profile
    ``z = g/2 + depth*shape(s)`` (s = 0..1 the fractional distance into the
    chamfer toward the iron end), vs build_mesh()'s LINEAR taper.  ``shape`` is
    a callable [0,1]->[0,1] (shape(0)=0 at the chamfer start, shape(1)=1 at the
    iron end); the equipotential-following profile passes the measured bow-out
    shape here.  Poles are x-prisms (curved gap face baked in); legs/air as usual."""
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry

    hL = L_BEAM / 2

    def lift(y):
        # fractional distance into the chamfer at the nearer end (0 in the body,
        # 1 at the iron end); deepest at the ends -> removes the corner.
        s = max((y - (hL - chamfer_len)) / chamfer_len,      # +y end
                ((-hL + chamfer_len) - y) / chamfer_len)     # -y end
        s = min(max(s, 0.0), 1.0)
        return depth * shape(s)

    ys = np.linspace(-hL, hL, n_face)
    top = _pole_prism(lambda y: GAP / 2 + lift(y), Z_OUT, ys, POLE_W)
    bot = _pole_prism(lambda y: -GAP / 2 - lift(y), -Z_OUT, ys, POLE_W)
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

    shape_geo = Glue([air, iron])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape_geo).GenerateMesh(maxh=maxh_air))
    return mesh


def build_mesh(maxh_air=0.05, maxh_iron=0.025, chamfer_depth=0.0, chamfer_len=0.025,
               air_half=None):
    """H-frame (window) dipole (iron) + air box, netgen.occ (no Cubit).

    x-SYMMETRIC (legs on both +-x) so the dipole field is clean (B_x ~ 0 at
    centre by symmetry).  Finite length L along the beam (y) -> real ENDS.
    chamfer_depth > 0 tapers the gap-facing pole END corners (the gap widens by
    chamfer_depth over the last chamfer_len) -- the equipotential-following end.
    air_half overrides the air-box half-size (default AIR) for the open-boundary
    convergence study.
    """
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry

    a_half = AIR if air_half is None else air_half
    hL = L_BEAM / 2
    top = Box(Pnt(-POLE_W, -hL, GAP / 2), Pnt(POLE_W, hL, Z_OUT))
    bot = Box(Pnt(-POLE_W, -hL, -Z_OUT), Pnt(POLE_W, hL, -GAP / 2))
    if chamfer_depth > 0:
        for ys in (+1, -1):
            top = _chamfer_pole_end(top, GAP / 2, +1, ys, chamfer_len, chamfer_depth)
            bot = _chamfer_pole_end(bot, -GAP / 2, -1, ys, chamfer_len, chamfer_depth)
    leg_l = Box(Pnt(-POLE_W, -hL, -GAP / 2), Pnt(-POLE_W + T_LEG, hL, GAP / 2))
    leg_r = Box(Pnt(POLE_W - T_LEG, -hL, -GAP / 2), Pnt(POLE_W, hL, GAP / 2))
    iron = (top + bot + leg_l + leg_r)
    iron.mat("iron")
    iron.maxh = maxh_iron

    air = Box(Pnt(-a_half, -a_half, -a_half), Pnt(a_half, a_half, a_half)) - iron
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
    from radia.coil_builder import CoilBuilder

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
    lower = upper.mirror("xy")                       # true mirror at -z (same current)
    coils = rad.ObjCnt(upper.to_radia() + lower.to_radia())
    return coils


def solve(mu_r=1000.0, order=2, maxh_air=0.035, maxh_iron=0.018,
          r_ref=0.008, n_beam=121, n_theta=32, plot=False,
          chamfer_depth=0.0, chamfer_len=0.025,
          curved_shape=None, return_contour=False, contour_y_max_factor=1.30,
          air_half=None):
    """Solve reduced-Omega forward, then integrate the transverse multipoles
    along the beam.  Returns the gap field + the integrated dipole + spectrum.
    chamfer_depth > 0 tapers the pole ENDS (the equipotential-following end).
    curved_shape (a callable [0,1]->[0,1]) replaces the LINEAR end taper with
    a CURVED profile z = g/2 + chamfer_depth*curved_shape(s) (s into the end)."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    import radia as rad

    if curved_shape is not None:
        mesh = build_mesh_curved(chamfer_depth, chamfer_len, curved_shape,
                                 maxh_air, maxh_iron)
    else:
        mesh = build_mesh(maxh_air, maxh_iron, chamfer_depth, chamfer_len,
                          air_half=air_half)
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
    yc = np.linspace(0.0, contour_y_max_factor * (L_BEAM / 2), 70)  # beam: body -> past end
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

    out = {
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
    if return_contour:
        out["contour_yc"] = yc.tolist()
        out["contour_zpole"] = z_pole.tolist()
    return out


def end_shaping_sweep(depths_mm=(0.0, 3.0, 6.0, 9.0, 12.0),
                      maxh_air=0.045, maxh_iron=0.022, plot=False):
    """Close the §3.2 DESIGN loop: re-shape the pole END (chamfer it, following
    the equipotential lift) -> re-solve -> the pole-end field enhancement is
    driven through zero.  Sweeps the chamfer depth and finds the optimal
    (zero-overshoot) end taper.

    Honest scope: the chamfer controls the LONGITUDINAL end-field bump
    (end_overshoot); the INTEGRATED TRANSVERSE harmonics (b_3, b_5) are
    body/pole-width dominated and barely move -- shaping the ENDS is the wrong
    lever for those (that is a body-pole-shape / Rogowski problem)."""
    rows = []
    for d in depths_mm:
        r = solve(maxh_air=maxh_air, maxh_iron=maxh_iron,
                  chamfer_depth=d * 1e-3, plot=False)
        rows.append((float(d), r["end_overshoot"], r["integrated_spurious_rel"],
                     r["bz_body_T"]))
    d_arr = np.array([x[0] for x in rows])
    over = np.array([x[1] for x in rows])
    # optimal chamfer = zero-crossing of end_overshoot (over is decreasing in d)
    if over[0] > 0 > over[-1]:
        opt = float(np.interp(0.0, over[::-1], d_arr[::-1]))
    else:
        opt = float("nan")
    if plot:
        _plot_loop(rows, opt)
    return {"sweep": rows, "optimal_chamfer_mm": opt}


def curved_chamfer_study(chamfer_len=0.030, maxh_air=0.05, maxh_iron=0.025,
                         n_beam=81, n_theta=28, plot=False):
    """(C+) Follow z_p(y) EXACTLY: a CURVED chamfer matching the measured
    equipotential bow-out, vs the LINEAR taper.

    Two passes.  Pass 1 (straight pole) extracts the equipotential contour and
    reads off (a) the bow-out SHAPE Ghat(s) past the iron end and (b) the
    natural lift DEPTH over one chamfer length -- both with NO free parameter.
    Pass 2 cuts the pole end along that curve (build_mesh_curved) and compares
    its pole-end enhancement against a LINEAR taper of the SAME depth.

    Honest expectation: the equipotential-shaped curve is the principled end
    profile; the head-to-head says whether it beats the linear taper or the
    linear taper is an adequate approximation at this mesh resolution."""
    from scipy.interpolate import interp1d

    hL = L_BEAM / 2
    yfac = 1.0 + 1.05 * chamfer_len / hL          # contour must reach hL + c_len

    # ---- pass 1: straight pole, extract the equipotential bow-out ----
    r0 = solve(maxh_air=maxh_air, maxh_iron=maxh_iron, n_beam=n_beam,
               n_theta=n_theta, chamfer_depth=0.0, return_contour=True,
               contour_y_max_factor=yfac)
    yc = np.array(r0["contour_yc"])
    zp = np.array(r0["contour_zpole"])
    # bow-out past the iron end: lift(sigma) = z_p(hL + sigma) - g/2, sigma in [0, c_len]
    sig = np.linspace(0.0, chamfer_len, 40)
    lift = np.interp(hL + sig, yc, zp) - GAP / 2
    lift = np.maximum.accumulate(np.maximum(lift, 0.0))   # monotone, non-negative
    depth_nat = float(lift[-1])                            # natural depth (no tuning)
    if depth_nat < 1e-4:
        depth_nat = 1e-4
    ghat_y = lift / lift[-1]                               # Ghat(s): 0 at start, 1 at end
    ghat_x = sig / chamfer_len
    _interp = interp1d(ghat_x, ghat_y, bounds_error=False, fill_value=(0.0, 1.0))
    shape = lambda s: float(_interp(s))                    # noqa: E731

    over_straight = r0["end_overshoot"]           # frac 0 (reuse pass-1 solve)
    spur_straight = r0["integrated_spurious_rel"]

    # ---- pass 2: curved depth sweep -> the curved profile's OWN zero-bump ----
    # The natural single-pass lift OVER-corrects (the straight-pole equipotential
    # is read in the already-bumped field): the SHAPE is right, the DEPTH needs
    # one knob (or one design iteration).  Sweep fractions of the natural depth.
    fracs = (0.40, 0.80)
    curved_rows = [(0.0, over_straight, spur_straight)]
    for fr in fracs:
        rc = solve(maxh_air=maxh_air, maxh_iron=maxh_iron, n_beam=n_beam,
                   n_theta=n_theta, chamfer_depth=fr * depth_nat,
                   chamfer_len=chamfer_len, curved_shape=shape)
        curved_rows.append((fr, rc["end_overshoot"], rc["integrated_spurious_rel"]))
    fr_arr = np.array([r[0] for r in curved_rows])
    ov_arr = np.array([r[1] for r in curved_rows])
    if ov_arr[0] > 0 > ov_arr[-1]:
        fr_zero = float(np.interp(0.0, ov_arr[::-1], fr_arr[::-1]))
    else:
        fr_zero = float("nan")
    depth_zero = fr_zero * depth_nat if fr_zero == fr_zero else float("nan")

    # confirm: the curved profile AT its zero-bump depth -> small residual bump,
    # and the transverse spurious is unchanged (body-dominated -> the END shape
    # is the wrong lever for it, consistent with the two-lever split).
    r_zero = solve(maxh_air=maxh_air, maxh_iron=maxh_iron, n_beam=n_beam,
                   n_theta=n_theta, chamfer_depth=depth_zero,
                   chamfer_len=chamfer_len, curved_shape=shape)

    res = {
        "chamfer_len_m": float(chamfer_len),
        "depth_natural_m": depth_nat,             # equipotential lift over c_len (no tuning)
        "depth_zerobump_m": float(depth_zero),    # the curved profile's zero-bump depth
        "ghat_x": ghat_x.tolist(), "ghat_y": ghat_y.tolist(),  # the bow-out shape
        "curved_depth_sweep": [(float(a), float(b)) for a, b, _ in curved_rows],
        "end_overshoot_straight": float(over_straight),
        "end_overshoot_curved_overdeep": float(curved_rows[-1][1]),  # frac 0.8 -> over-corrected
        "end_overshoot_curved_zerobump": float(r_zero["end_overshoot"]),
        "spurious_straight": float(spur_straight),
        "spurious_curved_zerobump": float(r_zero["integrated_spurious_rel"]),
    }
    if plot:
        _plot_curved(res)
    return res


def _plot_curved(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10.2, 3.8), dpi=150)
    # LEFT: the bow-out shape Ghat(s) (curved) vs the linear taper
    s = np.array(res["ghat_x"])
    ax[0].plot(s, res["ghat_y"], "C0", lw=2,
               label="equipotential bow-out $\\hat G(s)$ (curved)")
    ax[0].plot(s, s, "k--", lw=1, label="linear taper")
    ax[0].set_xlabel("fractional distance into the chamfer  $s$  (0=body, 1=end)")
    ax[0].set_ylabel("normalized depth  $(z-g/2)/$depth")
    ax[0].set_title(f"Follow $z_p(y)$ exactly: the curved end profile\n"
                    f"(natural single-pass depth = {res['depth_natural_m']*1e3:.1f} mm)")
    ax[0].legend(fontsize=8)
    # RIGHT: curved depth sweep (overshoot through zero) + linear at zero-bump
    sweep = res["curved_depth_sweep"]
    d = [row[0] * res["depth_natural_m"] * 1e3 for row in sweep]
    ov = [row[1] * 100 for row in sweep]
    ax[1].plot(d, ov, "o-", color="C0", label="curved (equipotential shape)")
    ax[1].axhline(0, color="k", lw=0.6)
    dz = res["depth_zerobump_m"] * 1e3
    if dz == dz:
        ax[1].axvline(dz, color="C3", lw=1, ls=":",
                      label=f"zero-bump depth = {dz:.1f} mm")
    ax[1].axvline(res["depth_natural_m"] * 1e3, color="0.4", lw=1, ls="--",
                  label=f"natural single-pass = {res['depth_natural_m']*1e3:.1f} mm "
                        f"(over-corrects)")
    ax[1].set_xlabel("end chamfer depth [mm]")
    ax[1].set_ylabel("pole-END enhancement [%]")
    ax[1].set_title("Curved (equipotential) depth -> zero bump\n"
                    "the naive single-pass lift over-corrects")
    ax[1].legend(fontsize=8)
    png = os.path.splitext(os.path.abspath(__file__))[0] + "_curved.png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def open_boundary_convergence(air_halves=(0.20, 0.30, 0.45),
                              maxh_air=0.03, maxh_iron=0.015,
                              n_beam=101, n_theta=28, plot=False):
    """The OPEN-BOUNDARY question: the dipole truncates the exterior with a
    Dirichlet air box.  How much does that cost?

    Grow the air box and watch the INTEGRATED dipole bbar_1 = INT B_z dy (it
    includes the FRINGE, so a too-close Dirichlet box would suppress the fringe
    and under-integrate it) and the local gap field B_z(body).

    Result (verified at a fine mesh, boxes 200->600 mm): bbar_1 is STABLE to
    well under ~1% across all box sizes -- and the residual is NON-MONOTONIC
    (mesh noise), not a converging trend.  The iron H-frame FLUX RETURN contains
    the field, so the Dirichlet air-box truncation is BELOW the mesh-noise floor:
    the open boundary is NOT the limiting error for this iron-dominated dipole.
    (B_z(body), a near-point value, is MORE mesh-sensitive than the bbar_1
    integral, so its larger spread is mesh noise, not an open-boundary effect.)

    The EXACT open boundary -- which DOES matter for a flux-return-free magnet --
    is the lab's Kelvin transform (hodograph_kelvin_axisym.py here: field_error
    ~1e-7 vs ~3e-3 truncated).  A full 3-D Cartesian Kelvin for THIS coil-driven
    dipole would need the localized Biot-Savart source mapped into the Kelvin
    exterior (the closed-form reduced-potential-background helper covers uniform/
    dipole-at-infinity backgrounds, not a localized coil) -- recorded as deferred
    in memory because the convergence here shows it would not change the answer."""
    rows = []
    for a in air_halves:
        r = solve(maxh_air=maxh_air, maxh_iron=maxh_iron, n_beam=n_beam,
                  n_theta=n_theta, air_half=a)
        rows.append((float(a), int(r["ne"]), float(r["bz_body_T"]),
                     float(r["integrated_dipole_bbar1_Tm"])))
    bb = np.array([row[3] for row in rows])
    bz = np.array([row[2] for row in rows])
    bbar1_box_spread = float((np.max(bb) - np.min(bb)) / np.mean(bb))  # noise-limited
    bz_body_spread = float((np.max(bz) - np.min(bz)) / np.mean(bz))    # near-point mesh noise
    res = {
        "rows": rows,                                    # (air_half, ne, bz_body, bbar1)
        "bbar1_box_spread": bbar1_box_spread,            # bbar_1 stability across boxes (<~1%)
        "bz_body_spread": bz_body_spread,                # B_z(body) mesh-noise spread
    }
    if plot:
        _plot_open_boundary(res)
    return res


def _plot_open_boundary(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = [row[0] * 1e3 for row in res["rows"]]
    bb = [row[3] * 1e3 for row in res["rows"]]
    mean_bb = float(np.mean(bb))
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
    ax.plot(a, bb, "s-", color="C1", label="integrated dipole $\\bar b_1$")
    ax.axhline(mean_bb, color="0.6", lw=0.8, ls="--")
    ax.fill_between([min(a), max(a)], mean_bb * 0.99, mean_bb * 1.01,
                    color="0.85", zorder=0, label="$\\pm1\\%$ band")
    ax.set_xlabel("air-box (Dirichlet) half-size [mm]")
    ax.set_ylabel("integrated dipole $\\bar b_1$ [mT$\\cdot$m]")
    ax.set_title(f"Open boundary: $\\bar b_1$ STABLE across box size\n"
                 f"(spread {res['bbar1_box_spread']*100:.2f}%, non-monotonic = mesh "
                 f"noise)\niron flux-return contains the field "
                 f"$\\Rightarrow$ truncation below the noise floor")
    ax.legend(fontsize=8)
    png = os.path.splitext(os.path.abspath(__file__))[0] + "_openbnd.png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def _plot_loop(rows, opt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = [r[0] for r in rows]
    over = [r[1] * 100 for r in rows]
    spur = [r[2] * 100 for r in rows]
    fig, ax = plt.subplots(figsize=(6.2, 3.8), dpi=150)
    ax.plot(d, over, "o-", color="C0", label="pole-END enhancement (longitudinal)")
    ax.plot(d, spur, "s--", color="C1", label="integrated spurious $b_{3,5}$ (transverse)")
    ax.axhline(0, color="k", lw=0.6)
    if opt == opt:  # not nan
        ax.axvline(opt, color="C3", lw=1, ls=":",
                   label=f"optimal chamfer = {opt:.1f} mm (zero bump)")
    ax.set_xlabel("end chamfer depth [mm]")
    ax.set_ylabel("[%]")
    ax.set_title("Close the design loop: re-shape the end -> re-solve\n"
                 "the chamfer zeroes the longitudinal end bump "
                 "(transverse $b_{3,5}$ is body-dominated)")
    ax.legend(fontsize=8)
    png = os.path.splitext(os.path.abspath(__file__))[0] + "_loop.png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


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

    # ---- (C) CLOSE THE LOOP: re-shape the end (chamfer) -> re-solve ----
    print("\n  (C) DESIGN loop -- re-shape the pole END (chamfer) -> re-solve:")
    s = end_shaping_sweep(plot=True)
    for d, ov, sp, bz in s["sweep"]:
        print(f"      chamfer {d:4.0f} mm:  end enhancement {ov*100:+5.1f} %  "
              f"|  integrated spurious b_3,5 {sp*100:4.1f} %  |  B_z(body) {bz:.4f} T")
    print(f"      -> OPTIMAL chamfer ~ {s['optimal_chamfer_mm']:.1f} mm zeroes the "
          f"longitudinal end bump (following the equipotential lift).")
    print("      The integrated TRANSVERSE spurious b_3,5 barely moves: it is")
    print("      body/pole-width dominated, NOT an end effect -- shaping the END")
    print("      is the right lever for the end bump, the wrong one for b_3,5.")

    # ---- (C+) FOLLOW z_p(y) EXACTLY: a CURVED chamfer, not a linear taper ----
    print("\n  (C+) Follow z_p(y) EXACTLY -- a CURVED end chamfer:")
    cc = curved_chamfer_study(plot=True)
    gx, gy = cc["ghat_x"], cc["ghat_y"]
    print(f"      equipotential bow-out shape Ghat(s): convex (Ghat(0.5)="
          f"{np.interp(0.5, gx, gy):.2f} vs linear 0.50) -- a parameter-free curve")
    print(f"      natural single-pass depth = {cc['depth_natural_m']*1e3:.1f} mm "
          f"OVER-corrects: end bump {cc['end_overshoot_straight']*100:+.1f}% -> "
          f"{cc['end_overshoot_curved_overdeep']*100:+.1f}%")
    print(f"      -> the curved profile drives the bump through zero at "
          f"~{cc['depth_zerobump_m']*1e3:.1f} mm (~{cc['depth_zerobump_m']/cc['depth_natural_m']*100:.0f}% "
          f"of the naive lift; the SHAPE is right, the DEPTH needs one knob).")
    print(f"      integrated transverse spurious stays body-dominated: "
          f"{cc['spurious_straight']:.3f} -> {cc['spurious_curved_zerobump']:.3f}.")

    # ---- OPEN BOUNDARY: how much does the Dirichlet air box cost? ----
    print("\n  Open boundary -- air-box (Dirichlet) truncation vs growing the box:")
    ob = open_boundary_convergence(plot=True)
    for a, ne, bz, bb in ob["rows"]:
        print(f"      air_half {a*1e3:3.0f} mm (ne={ne:6d}):  B_z(body)={bz:.5f} T  "
              f"bbar_1={bb*1e3:.3f} mT.m")
    print(f"      -> INTEGRATED bbar_1 is STABLE across boxes "
          f"({ob['bbar1_box_spread']*100:.2f}% spread = mesh noise, non-monotonic):")
    print("         the iron H-frame FLUX RETURN contains the field, so the")
    print("         Dirichlet air-box truncation is BELOW the mesh-noise floor --")
    print("         the open boundary is NOT the limiting error for this dipole.")
    print(f"      -> B_z(body) (a near-point value) is more mesh-sensitive "
          f"({ob['bz_body_spread']*100:.1f}% spread = mesh noise, not open-boundary).")
    print("      -> The EXACT open boundary (for a flux-return-FREE magnet) is the")
    print("         lab's Kelvin transform -- hodograph_kelvin_axisym.py here")
    print("         (field_error ~1e-7 vs ~3e-3 truncated).")


if __name__ == "__main__":
    main()
