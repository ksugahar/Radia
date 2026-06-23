r"""FFAG scaling SECTOR cell: the TWO-PLANE -> 3D design method (iron poles).

THE METHOD (what this file is the first rung of)
------------------------------------------------
An accelerator (bending / focusing) magnet is designed in TWO orthogonal 2-D
planes and the two designs are REFLECTED into one 3-D pole:

  * TRANSVERSE plane -- here the (r, z) radial-vertical cross-section: sets the
    field the beam SEES.  For a scaling FFAG the median-plane radial profile is
    B_z(r) ~ r^k (field index k = const => geometrically-similar orbits =>
    momentum-independent tune).  The pole gap g(r) = g0 (r/r0)^{-k} encodes it.
    Solved by scaling_ffag_pole_2d.py (reused here as Plane A).

  * AZIMUTHAL plane -- here the (s, z) plane along the orbit at a fixed radius:
    sets how the field TURNS ON / OFF along the path -- the sector ENDS, the
    fringe, and the effective magnetic length L_eff = INT B_z ds / B_z(body).
    A finite-length iron pole of gap g(r0) over the sector arc s in [-L/2, L/2].
    Solved here (Plane B, a new ngsolve scalar-potential solve).

The 3-D pole is the (r, z) gap profile SWEPT around the sector arc and
TRUNCATED at the azimuthal ends shaped by Plane B.  The two 2-D designs are
EXACT in the body and couple only at the ends; the prior leaf-coupling study
(leaf_coupling_perturbation_3d.py) measured that coupling ~ gap / L, so the
two-plane reflection is controlled by the aspect ratio L_sector / g(r).  This
file reports that aspect and the fringe excess -- the lever that says where the
2-plane stack is exact and where the genuine 3-D end is needed.

SCOPE (rung 1, this file -- ngsolve only, no radia)
---------------------------------------------------
  Plane A : reuse scaling_ffag_pole_2d -> k(r), g(r) over the momentum aperture.
  Plane B : NEW (s, z) finite-length pole -> B_z(s), L_eff, fringe at r0.
  Compose : tie Plane A's g(r0) to Plane B's sector ends; report the validity
            aspect L_sector / g and the fringe excess (the leaf-coupling lever).
  The 3-D curved-sector forward solve (CoilBuilder arc + reduced-Omega, radia)
  is rung 2 (needs the radia binaries; runs in the main tree).

run:  python ffag_sector_two_plane.py            # Plane A + Plane B + compose
      python ffag_sector_two_plane.py --fig       # + figure
      python ffag_sector_two_plane.py --dtheta 0.5  # wider sector (aspect sweep)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))      # scaling_ffag

from scaling_ffag_pole_2d import (                                   # noqa: E402
    K_INDEX, G0, R0, MU0,
    aperture_radii, scaling_gap, field_index, solve_complementary,
)


# --------------------------------------------------------------------------- #
# Plane B: the AZIMUTHAL (s, z) plane -- a finite-length iron pole at radius r0.
# s = arc length along the orbit (s = r0 * theta); z = the gap height.  The pole
# spans the sector arc s in [-L/2, L/2]; past the ends the gap field fringes.
# This is the SAME end physics as accel_pole_ends_fem.py, in 2-D, ngsolve only
# (scalar magnetic potential, no coil -- the iron back is driven at Psi = mmf).
# --------------------------------------------------------------------------- #
def _azimuthal_geometry(s_half, z_top, sector_half, g_ref, pole_t):
    """Upper-half (s, z) air box [-s_half, s_half] x [0, z_top] with an iron pole
    block s in [-sector_half, sector_half], z in [g_ref/2, g_ref/2 + pole_t].
    Boundaries: 'median' (z=0), 'irontop' (the driven pole back), 's_wall'
    (|s| = s_half), 'top' (z = z_top); the gap/iron interface is interior."""
    from netgen.occ import WorkPlane, OCCGeometry, Glue
    air = WorkPlane().MoveTo(-s_half, 0.0).Rectangle(2.0 * s_half, z_top).Face()
    iron = (WorkPlane().MoveTo(-sector_half, g_ref / 2.0)
            .Rectangle(2.0 * sector_half, pole_t).Face())
    gap_air = air - iron
    gap_air.faces.name = "air"
    iron.faces.name = "iron"
    shape = Glue([gap_air, iron])
    z_back = g_ref / 2.0 + pole_t
    for e in shape.edges:
        c = e.center
        if abs(c.y) < 1e-9:
            e.name = "median"
        elif abs(c.y - z_top) < 1e-9:
            e.name = "top"
        elif abs(abs(c.x) - s_half) < 1e-9:
            e.name = "s_wall"
        elif abs(c.y - z_back) < 1e-9 and abs(c.x) < sector_half + 1e-9:
            e.name = "irontop"
        else:
            e.name = "poleface"                                # interior interface
    return OCCGeometry(shape, dim=2)


def solve_azimuthal_end(r_ref, sector_arclen, k=K_INDEX, g0=G0, r0=R0,
                        B_design=1.0, mu_r=2000.0, pole_t=0.03, order=3,
                        win_g=6.0, z_buffer=3.0):
    """Plane B: solve the (s, z) finite-length pole at radius r_ref with gap
    g_ref = g(r_ref) over the sector arc `sector_arclen` (= r0 * dtheta).
    Returns the on-axis B_z(s) profile, the effective magnetic length L_eff,
    and the fringe excess (L_eff - L_sector) / L_sector -- the azimuthal-end
    design quantities.  Linear iron (mu_r); the body field is set to ~ B_design
    by the mmf so the numbers read in Tesla.

    The magnetic-length integral uses a FIXED beyond-end window W = win_g * g_ref
    (independent of the sector length) so the fringe is captured WITHOUT the
    window growing with L -- otherwise a length-proportional window pollutes the
    fringe-vs-aspect scaling (the end effect is local, ~g, not ~L)."""
    from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, TaskManager,
                         Mesh)
    g_ref = float(scaling_gap(r_ref, k, g0, r0))
    sector_half = 0.5 * sector_arclen
    win = win_g * g_ref                                        # fixed end window
    s_half = sector_half + win + z_buffer * g_ref              # domain past window
    z_top = 0.5 * g_ref + pole_t + z_buffer * g_ref
    maxh = g_ref / 3.0
    mmf = B_design * g_ref / (2.0 * MU0)                        # B_body ~ B_design
    n_eval = max(161, int(2.0 * (sector_half + win) / (g_ref / 20.0)))
    ss = np.linspace(-(sector_half + win), sector_half + win, n_eval)
    y_probe = 0.02 * g_ref
    geo = _azimuthal_geometry(s_half, z_top, sector_half, g_ref, pole_t)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        mu = mesh.MaterialCF({"iron": mu_r}, default=1.0)
        a = BilinearForm(fes)
        a += mu * grad(u) * grad(v) * dx
        a.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0),
                definedon=mesh.Boundaries("median|irontop"))
        rhs = gfu.vec.CreateVector()
        rhs.data = -a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * rhs
        g = grad(gfu)
        # gap field on the median: air mu_r=1 => B_z = -mu0 dPsi/dz
        Bz = np.array([-MU0 * g(mesh(float(s), y_probe))[1] for s in ss])
        ndof = int(fes.ndof)
        ne = int(mesh.ne)
    body = np.abs(ss) < 0.3 * sector_half
    Bz_body = float(np.mean(Bz[body]))
    L_eff = (float(np.trapezoid(Bz, ss) / Bz_body)
             if abs(Bz_body) > 1e-30 else 0.0)
    fringe_excess = (L_eff - sector_arclen) / sector_arclen
    # pole-end enhancement: peak |B_z| near the sector end vs the body
    end = np.abs(np.abs(ss) - sector_half) < 0.12 * sector_arclen
    end_overshoot = (float(np.max(np.abs(Bz[end])) / abs(Bz_body)) - 1.0
                     if end.any() and abs(Bz_body) > 1e-30 else 0.0)
    return {
        "r_ref": float(r_ref), "g_ref": g_ref,
        "sector_arclen": float(sector_arclen),
        "aspect_L_over_g": float(sector_arclen / g_ref),
        "ss": ss, "Bz": Bz, "Bz_body": Bz_body,
        "L_eff": L_eff, "fringe_excess": float(fringe_excess),
        "end_overshoot": float(end_overshoot),
        "ne": ne, "ndof": ndof,
    }


# --------------------------------------------------------------------------- #
# Compose the two planes into the two-plane spec.
# --------------------------------------------------------------------------- #
def design_two_plane(k=K_INDEX, g0=G0, r0=R0, dtheta=0.30, B_design=1.0,
                     order_A=4, maxh_A=0.012, order_B=3):
    """Two-plane FFAG-sector design.

    Plane A : scaling_ffag k(r) over the momentum aperture (transverse, (r,z)).
    Plane B : azimuthal end at the reference radius r0, sector arc r0 * dtheta.
    Compose : the validity aspect L_sector / g(r0) and the fringe excess.
    """
    r_min, r_max = aperture_radii(k=k, r0=r0)
    # ---- Plane A (transverse): field index k(r) with the A/phi bracket ----
    A = solve_complementary(r_min, r_max, k=k, g0=g0, r0=r0,
                            order=order_A, maxh=maxh_A)
    m = slice(2, -2)                                            # interior window
    k_mid = 0.5 * (A["k_phi"] + A["k_A"])
    k_dev = float(np.max(np.abs(k_mid[m] - k)))                 # flatness vs design
    bracket = float(np.max(np.abs(A["k_phi"] - A["k_A"])))
    # ---- Plane B (azimuthal): sector ends at r0 ----
    sector_arclen = r0 * dtheta
    B = solve_azimuthal_end(r0, sector_arclen, k=k, g0=g0, r0=r0,
                            B_design=B_design, order=order_B)
    return {
        "k_design": float(k), "aperture": (float(r_min), float(r_max)),
        "dtheta_rad": float(dtheta), "r0": float(r0),
        "plane_A": {
            "k_interior_mean": float(np.mean(k_mid[m])),
            "k_interior_dev_from_design": k_dev,
            "bracket_gap_max": bracket,
            "ndof": A["ndof"],
            "_sol": A,
        },
        "plane_B": B,
        "compose": {
            "g_r0": B["g_ref"],
            "L_sector": B["sector_arclen"],
            "aspect_L_over_g": B["aspect_L_over_g"],
            "L_eff": B["L_eff"],
            "fringe_excess": B["fringe_excess"],
            "end_overshoot": B["end_overshoot"],
        },
    }


def aspect_sweep(dthetas=(0.20, 0.35, 0.60, 1.00), k=K_INDEX, g0=G0, r0=R0,
                 B_design=1.0, order_B=3):
    """Plane B across sector angles: confirm the fringe excess falls as the
    aspect L_sector / g grows -- the same leaf-coupling ~ gap/L law, now on the
    azimuthal (sector) plane.  Reports the log-log slope of fringe vs aspect."""
    rows = []
    for dt in dthetas:
        B = solve_azimuthal_end(r0, r0 * dt, k=k, g0=g0, r0=r0,
                                B_design=B_design, order=order_B)
        rows.append({"dtheta": float(dt),
                     "aspect_L_over_g": B["aspect_L_over_g"],
                     "fringe_excess": B["fringe_excess"],
                     "end_overshoot": B["end_overshoot"],
                     "L_eff": B["L_eff"]})
    asp = np.array([r["aspect_L_over_g"] for r in rows])
    fr = np.array([abs(r["fringe_excess"]) for r in rows])
    ok = (asp > 0) & (fr > 0)
    slope = (float(np.polyfit(np.log(asp[ok]), np.log(fr[ok]), 1)[0])
             if ok.sum() >= 2 else float("nan"))
    return {"rows": rows, "fringe_vs_aspect_slope": slope,
            "monotone_decay": bool(np.all(np.diff(fr) < 0))}


# --------------------------------------------------------------------------- #
# Rung 2: the 3-D REFLECTION.  Sweep the (r,z) scaling gap profile around the
# sector arc into a 3-D iron pole and forward-solve, verifying the orbit sees
# B_z(r) ~ r^k.  Iron-pole (equipotential) drive: an upper-half model (median
# z=0 = the up-down antisymmetry plane, Psi=0) with the high-mu pole back driven
# at Psi=mmf -- the 3-D form of Plane A/B's scalar-potential pole.  The field
# index is set by the pole GEOMETRY (B ~ 1/g(r) under a high-mu equipotential),
# so it is drive-independent; the CoilBuilder + reduced-Omega coil drive sets
# the AMPLITUDE (a further step), not the index.  ngsolve only.
# --------------------------------------------------------------------------- #
def build_sector_pole_3d(dtheta, k=K_INDEX, g0=G0, r0=R0, r_lo=None, r_hi=None,
                         z_pole=0.0, n_face=33, maxh_iron=None, maxh_air=None):
    """Revolve the (r,z) scaling gap cross-section (gap face z=g(r)/2, top flat
    at z_pole) around the Z axis by `dtheta` (rad) -> the 3-D TOP sector pole
    (upper half-model).  Returns (mesh, sector_center_azimuth, (r_lo, r_hi))."""
    import ngsolve as ng
    from netgen.occ import (WorkPlane, Axes, Pnt, X, Y, Z, Axis, Box, Glue,
                            OCCGeometry)
    # the pole spans a radial buffer BEYOND the measured momentum aperture so the
    # radial-wall fringing stays out of the field-index window (cf. Plane A).
    if r_lo is None:
        r_lo = aperture_radii(k=k, r0=r0)[0] * 0.86
    if r_hi is None:
        r_hi = aperture_radii(k=k, r0=r0)[1] * 1.14
    g_lo = float(scaling_gap(r_lo, k, g0, r0))                # widest gap (inner r)
    g_hi = float(scaling_gap(r_hi, k, g0, r0))                # narrowest gap (outer r, highest B)
    if maxh_iron is None:
        maxh_iron = g_hi / 1.3                                # resolve the smallest gap
    if maxh_air is None:
        maxh_air = g_hi * 2.2
    if z_pole <= 0.0:
        z_pole = g_lo / 2.0 + 0.05                            # inner gap face + iron
    rs = np.linspace(r_lo, r_hi, n_face)
    zf = scaling_gap(rs, k, g0, r0) / 2.0
    wp = WorkPlane(Axes(Pnt(0, 0, 0), n=-Y, h=X))             # MoveTo(r,z)->(r,0,z)
    wp.MoveTo(float(rs[0]), float(zf[0]))
    for r, z in zip(rs[1:], zf[1:]):
        wp.LineTo(float(r), float(z))                          # gap face z=g(r)/2
    wp.LineTo(float(r_hi), z_pole)
    wp.LineTo(float(r_lo), z_pole)
    wp.Close()
    pole = wp.Face().Revolve(Axis(Pnt(0, 0, 0), Z), float(np.degrees(dtheta)))
    pole.faces.Max(Z).name = "irontop"                         # the driven pole back
    pole.mat("iron")
    pole.maxh = maxh_iron
    mx = 0.5 * g_lo + 0.06
    az_pad = 0.20                                              # azimuthal fringe room
    x_lo, x_hi = r_lo * np.cos(dtheta) - mx, r_hi + mx
    y_lo, y_hi = -az_pad, r_hi * np.sin(dtheta) + az_pad
    box = Box(Pnt(x_lo, y_lo, 0.0), Pnt(x_hi, y_hi, z_pole))
    air = box - pole
    air.mat("air")
    air.maxh = maxh_air
    air.faces.Min(Z).name = "median"                           # z=0 antisymmetry plane
    shape = Glue([air, pole])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh_air))
    return mesh, 0.5 * dtheta, (float(r_lo), float(r_hi))


def solve_sector_field_index(dtheta=0.30, k=K_INDEX, g0=G0, r0=R0, B_design=1.0,
                             mu_r=2000.0, order=2, n_r=31, n_az=41, **geo):
    """3-D reflection forward-solve: B_z(r) along the median orbit at the sector
    center azimuth recovers the field index k(r); B_z(theta) along the arc at r0
    cross-checks Plane B's azimuthal end (effective magnetic length)."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, TaskManager)
    mesh, theta_c, (r_lo, r_hi) = build_sector_pole_3d(dtheta, k, g0, r0, **geo)
    mmf = B_design * float(scaling_gap(r0, k, g0, r0)) / (2.0 * MU0)
    zc = 0.02 * float(scaling_gap(r0, k, g0, r0))
    with TaskManager():
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        mu = mesh.MaterialCF({"iron": mu_r}, default=1.0)
        a = BilinearForm(fes)
        a += mu * grad(u) * grad(v) * dx
        a.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0),
                definedon=mesh.Boundaries("median|irontop"))
        rhs = gfu.vec.CreateVector()
        rhs.data = -a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * rhs
        g = grad(gfu)
        # (i) field index: B_z(r) on the median at the sector center azimuth,
        # measured over the momentum APERTURE (inside the wider pole, so the
        # radial-wall fringing is out of the window).
        r_meas_lo, r_meas_hi = aperture_radii(k=k, r0=r0)
        rr = np.linspace(r_meas_lo, r_meas_hi, n_r)
        Bz_r = np.array([-MU0 * g(mesh(float(rv * np.cos(theta_c)),
                                       float(rv * np.sin(theta_c)), zc))[2]
                         for rv in rr])
        # (ii) azimuthal end: B_z(theta) along the arc at r0, sampled PAST the
        # sector ends (fixed angular window ~1.6 g/r0 each side) so the fringe is
        # captured for the effective magnetic length.
        half = 0.5 * dtheta
        w_ang = 1.6 * float(scaling_gap(r0, k, g0, r0)) / r0
        ths = np.linspace(theta_c - half - w_ang, theta_c + half + w_ang, n_az)
        Bz_az = np.array([-MU0 * g(mesh(float(r0 * np.cos(t)),
                                        float(r0 * np.sin(t)), zc))[2]
                          for t in ths])
        ne, ndof = int(mesh.ne), int(fes.ndof)
    rmid, k_r = field_index(rr, Bz_r)
    m = slice(1, -1) if len(k_r) > 2 else slice(None)
    # azimuthal effective length along the arc (s = r0*theta), body-normalised
    s = r0 * ths
    body = np.abs(ths - theta_c) < 0.25 * dtheta
    Bz_body = float(np.mean(Bz_az[body]))
    L_eff_arc = (float(np.trapezoid(Bz_az, s) / Bz_body)
                 if abs(Bz_body) > 1e-30 else 0.0)
    return {
        "k_design": float(k), "dtheta": float(dtheta), "theta_c": float(theta_c),
        "k_mean": float(np.mean(k_r[m])), "k_min": float(np.min(k_r)),
        "k_max": float(np.max(k_r)), "k_dev_from_design": float(np.max(np.abs(k_r[m] - k))),
        "rr": rr, "Bz_r": Bz_r, "rmid": rmid, "k_r": k_r,
        "ths": ths, "Bz_az": Bz_az, "Bz_body_T": Bz_body,
        "L_eff_arc": L_eff_arc, "L_sector": float(r0 * dtheta),
        "ne": ne, "ndof": ndof,
    }


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dtheta", type=float, default=0.30,
                    help="sector half-cell angle (rad); L_sector = r0 * dtheta")
    ap.add_argument("--B-design", type=float, default=1.0, dest="B_design")
    ap.add_argument("--sweep", action="store_true",
                    help="sweep sector angle -> fringe vs aspect (leaf-coupling)")
    ap.add_argument("--rung2", action="store_true",
                    help="rung 2: the 3-D sector reflection (revolve g(r), recover "
                         "B(r) ~ r^k along the orbit); ~1-2 min")
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    print("=" * 76)
    print("FFAG scaling SECTOR cell -- the TWO-PLANE -> 3D design method (rung 1)")
    print("=" * 76)
    d = design_two_plane(dtheta=args.dtheta, B_design=args.B_design)
    A, B, C = d["plane_A"], d["plane_B"], d["compose"]
    rmin, rmax = d["aperture"]
    print("Plane A  (r,z) transverse -- the scaling field index:")
    print(f"  field index k_design        : {d['k_design']:.3f}")
    print(f"  radial aperture [r_min,r_max]: [{rmin:.4f}, {rmax:.4f}]")
    print(f"  k interior mean             : {A['k_interior_mean']:.3f}")
    print(f"  k interior dev from design  : {A['k_interior_dev_from_design']:.2e}")
    print(f"  A/phi bracket gap (max)     : {A['bracket_gap_max']:.2e}")
    print("Plane B  (s,z) azimuthal -- the sector ends at r0:")
    print(f"  gap at r0  g(r0)            : {C['g_r0']:.4f}")
    print(f"  sector arc L_sector=r0*dth  : {C['L_sector']:.4f}  (dtheta={d['dtheta_rad']:.3f})")
    print(f"  B_z(body)                   : {B['Bz_body']:.4f} T")
    print(f"  effective magnetic length   : {C['L_eff']:.4f}  (excess vs L_sector)")
    print(f"  pole-end enhancement        : {C['end_overshoot']*100:+.1f} %")
    print("Compose -- the two-plane validity lever (leaf-coupling ~ gap/L):")
    print(f"  aspect  L_sector / g(r0)    : {C['aspect_L_over_g']:.2f}")
    print(f"  fringe excess (L_eff/L - 1) : {C['fringe_excess']*100:+.1f} %"
          f"  ({'PERTURBATIVE' if abs(C['fringe_excess']) < 0.10 else 'NON-perturbative (3D end needed)'})")
    print("  => the (r,z) and (s,z) 2-D designs are EXACT in the body; the")
    print("     sector ENDS need the genuine 3-D end (rung 2: --rung2).")

    out = {
        "k_design": d["k_design"], "aperture": d["aperture"],
        "dtheta_rad": d["dtheta_rad"], "r0": d["r0"],
        "plane_A": {kk: A[kk] for kk in A if kk != "_sol"},
        "plane_B": {kk: (B[kk] if kk not in ("ss", "Bz") else None) for kk in B},
        "compose": C,
    }

    rung2 = None
    if args.rung2:
        print("\n" + "-" * 76)
        print("Rung 2 -- the 3-D REFLECTION: revolve g(r) around the arc, recover B(r)~r^k")
        print("-" * 76)
        rung2 = solve_sector_field_index(dtheta=args.dtheta, B_design=args.B_design)
        ex = 100.0 * (rung2["L_eff_arc"] - rung2["L_sector"]) / rung2["L_sector"]
        print(f"  3-D mesh ne / ndof          : {rung2['ne']} / {rung2['ndof']}")
        print(f"  field index k(r) recovered  : mean {rung2['k_mean']:.3f}  "
              f"[{rung2['k_min']:.2f}, {rung2['k_max']:.2f}]  (design {rung2['k_design']:.1f})")
        print(f"  => the swept g(r)~r^(-k) pole gives B_z(r)~r^k along the 3-D orbit")
        print(f"  body field B_z(body)        : {rung2['Bz_body_T']:.4f} T")
        print(f"  azimuthal L_eff (arc)       : {rung2['L_eff_arc']:.4f} "
              f"(excess {ex:+.1f}% vs L_sector; cross-checks Plane B)")
        out["rung2"] = {kk: rung2[kk] for kk in rung2
                        if kk not in ("rr", "Bz_r", "rmid", "k_r", "ths", "Bz_az")}

    if args.sweep:
        print("\n" + "-" * 76)
        print("Aspect sweep -- fringe excess vs L_sector/g (the azimuthal leaf-coupling)")
        print("-" * 76)
        sw = aspect_sweep(B_design=args.B_design)
        print(f"{'dtheta':<9}{'aspect L/g':<12}{'fringe %':<11}{'end %':<9}{'L_eff':<8}")
        for r in sw["rows"]:
            print(f"{r['dtheta']:<9.2f}{r['aspect_L_over_g']:<12.2f}"
                  f"{r['fringe_excess']*100:<+11.1f}{r['end_overshoot']*100:<+9.1f}"
                  f"{r['L_eff']:<8.4f}")
        print(f"  fringe ~ aspect^{sw['fringe_vs_aspect_slope']:.2f}"
              f"  (monotone decay: {sw['monotone_decay']})")
        out["aspect_sweep"] = sw

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(d, args)
        if rung2 is not None:
            _figure_rung2(rung2)


def _figure(d, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A, B, C = d["plane_A"]["_sol"], d["plane_B"], d["compose"]
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.0))
    # LEFT: Plane A -- field index k(r) with the A/phi bracket
    ri = A["r_index"]
    ax[0].plot(ri, A["k_phi"], "o-", ms=3, label="k_phi")
    ax[0].plot(ri, A["k_A"], "s-", ms=3, label="k_A")
    ax[0].fill_between(ri, np.minimum(A["k_phi"], A["k_A"]),
                       np.maximum(A["k_phi"], A["k_A"]), alpha=0.2,
                       label="A/phi bracket")
    ax[0].axhline(d["k_design"], color="k", ls="--", lw=1,
                  label=f"k_design={d['k_design']:.0f}")
    ax[0].set_xlabel("r"); ax[0].set_ylabel("field index k(r)")
    ax[0].set_title("Plane A  (r,z) transverse: the scaling field index")
    ax[0].legend(fontsize=8)
    # RIGHT: Plane B -- B_z(s)/B_body along the orbit (body flat + sector fringe)
    ss, Bz = B["ss"], B["Bz"] / B["Bz_body"]
    sh = 0.5 * C["L_sector"]
    ax[1].plot(ss, Bz, "C0")
    ax[1].fill_between(ss, Bz, 0, alpha=0.2, color="C0")
    ax[1].axvspan(-sh, sh, color="0.85", zorder=0, label="sector iron L")
    ax[1].axhline(1.0, color="C3", lw=0.8, ls="--", label="body")
    ax[1].set_xlabel("orbit  s = r0*theta")
    ax[1].set_ylabel("B_z(s) / B_z(body)")
    ax[1].set_title(f"Plane B  (s,z) azimuthal: sector ends\n"
                    f"L_eff/L_sector-1={C['fringe_excess']*100:+.1f}%, "
                    f"aspect L/g={C['aspect_L_over_g']:.1f}")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


def _figure_rung2(r2):
    """Rung 2 (3-D reflection): the recovered field index along the orbit + the
    azimuthal end profile from the full 3-D sector solve."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.0))
    # LEFT: the field index k(r) recovered from the 3-D solve vs design
    ax[0].plot(r2["rmid"], r2["k_r"], "o-", ms=3, color="C2",
               label=f"3-D solve  (mean {r2['k_mean']:.2f})")
    ax[0].axhline(r2["k_design"], color="k", ls="--", lw=1,
                  label=f"k_design={r2['k_design']:.0f}")
    ax[0].set_xlabel("r (orbit, sector center)")
    ax[0].set_ylabel("field index k(r)")
    ax[0].set_title("Rung 2: B_z(r)~r^k recovered from the 3-D swept-gap pole")
    ax[0].legend(fontsize=8)
    # RIGHT: the azimuthal end profile B_z(theta)/body along the arc
    th, bz = r2["ths"], r2["Bz_az"] / r2["Bz_body_T"]
    th_c, half = r2["theta_c"], 0.5 * r2["dtheta"]
    ax[1].plot(th, bz, "C2")
    ax[1].fill_between(th, bz, 0, alpha=0.2, color="C2")
    ax[1].axvspan(th_c - half, th_c + half, color="0.85", zorder=0,
                  label="sector iron")
    ax[1].axhline(1.0, color="C3", lw=0.8, ls="--", label="body")
    ax[1].set_xlabel("azimuth  theta (orbit arc)")
    ax[1].set_ylabel("B_z(theta) / B_z(body)")
    ax[1].set_title("Rung 2: azimuthal sector ends (3-D)\n"
                    f"L_eff/L_sector-1 = "
                    f"{100*(r2['L_eff_arc']-r2['L_sector'])/r2['L_sector']:+.1f}%")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(os.path.abspath(__file__))[0] + "_rung2.png"
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


if __name__ == "__main__":
    main()
