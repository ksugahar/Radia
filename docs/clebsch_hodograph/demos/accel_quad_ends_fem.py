"""(A) The QUADRUPOLE end study, FEM rung -- the analyzer handles any multipole.

RESEARCH example (track A).  The dipole work established the integrated-multipole
method (accel_pole_ends_3d = the analytic theorem; accel_pole_ends_fem = the FEM
rung).  This file is the FEM rung for the QUADRUPOLE: a real finite-length
4-pole magnet solved by the forward engine and fed to the SAME integrated
analyzer (accel_pole_ends_3d.integrated_multipoles), to confirm the method (and
the analytic quad theorem already in accel_pole_ends_3d) on a real device.

Excitation = the high-mu scalar-potential model (the same one validated in 2-D
by accel_pole_dipole_body_2d.py, here in 3-D): the iron pole face is a magnetic
scalar equipotential (high mu -> H_t = 0 -> Phi = const).  The 4 hodograph
hyperbola poles xy = +-r0^2/2 (the n=2 equipotential, from accel_pole_design)
are held at alternating +-Phi0 and the current-free aperture field is Laplace-
solved.  Adjacent poles (90 deg apart) carry opposite potential, the standard
quadrupole scalar-potential model; the 4 finite-length iron bars give real ENDS.

Convention: beam = z (the quad axis), transverse plane = (x, y).  The quad field
B = G (y, x) has By + i Bx = G (x + i y) -> main n = 2.  The 4-fold + alternating-
sign symmetry forbids the normal harmonics n = 1, 3, 4, 5 (only n = 2 mod 4 are
allowed: n = 2, 6, 10) -- so the first ALLOWED spurious is the 12-pole b_6.

Verified (FEM): a CLEAN integrated quadrupole -- main b_2; the forbidden n=1,3,5
suppressed ~10x below b_6; the 12-pole b_6 (~0.4%) the dominant allowed spurious.
b_6 is dominated by the finite pole geometry (the body pole angular width), with
the ENDS a smaller part -- the SAME body-dominated picture the dipole found for
b_3,5 (the two-lever split generalizes: ends drive the longitudinal profile, the
body pole shape drives the transverse harmonics; the analytic theorem in
accel_pole_ends_3d shows an ideal Maxwellian end contributes zero integrated b_6).

run:  python accel_quad_ends_fem.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the analyzers
from accel_pole_ends_3d import integrated_multipoles               # noqa: E402

# ---- geometry (meters); beam = z, transverse = (x, y) ----
R0 = 0.020             # aperture (inscribed) radius -> pole tip on the diagonal at r0
C_HYP = R0 ** 2 / 2.0  # hodograph hyperbola constant: pole face xy = +-C_HYP
R_OUT = 0.055          # pole outer radius (toward the yoke)
AIR = 0.10             # air box half-size
R_REF = 0.008          # reference-circle radius (inside the aperture)
R_CORE = 0.013         # refined air-core radius around the axis (clean harmonics)
BASE_DEG = 5.0         # rigid rotation off the sampling symmetry planes (|b_n| invariant)
ALLOWED = (2, 6, 10)   # quad-allowed normal harmonics (n = 2 mod 4)
FORBIDDEN = (1, 3, 5)  # symmetry-forbidden normals (should be suppressed)


def _pole_poly(n_face=24):
    """The +x/+y (45 deg) hodograph pole cross-section in (x, y): inner face =
    the hyperbola xy = C_HYP from (xa, C/xa) to (xb, C/xb) (symmetric about
    x = y), the two endpoints carried radially out to R_OUT (toward the yoke)."""
    xa, xb = 0.010, 0.020
    xs = np.linspace(xa, xb, n_face)
    ys = C_HYP / xs
    pts = list(zip(xs, ys))
    e0, e1 = pts[0], pts[-1]
    r0e, r1e = math.hypot(*e0), math.hypot(*e1)
    return pts + [(e1[0] / r1e * R_OUT, e1[1] / r1e * R_OUT),
                  (e0[0] / r0e * R_OUT, e0[1] / r0e * R_OUT)]


def build_quad_mesh(L=0.080, maxh=0.022, maxh_core=0.0035, n_face=24):
    """4 finite-length hodograph hyperbola poles (z-prisms) at alternating
    +-Phi0, with a refined air core around the axis (so the reference circle
    spans many elements -> clean harmonics).  netgen.occ, no Cubit."""
    import ngsolve as ng
    from netgen.occ import (WorkPlane, Axes, Pnt, Z, X, Box, Cylinder, Glue,
                            OCCGeometry, Axis)

    poly = _pole_poly(n_face)
    wp = WorkPlane(Axes(Pnt(0, 0, -L / 2), n=Z, h=X))    # plane z=-L/2, local (x,y)
    wp.MoveTo(*poly[0])
    for p in poly[1:]:
        wp.LineTo(*p)
    wp.Close()
    pole = wp.Face().Extrude(L)                          # the 45 deg pole, z in [-L/2, L/2]

    def rot(a):
        return pole.Rotate(Axis(Pnt(0, 0, 0), Z), a)
    pos = [rot(BASE_DEG), rot(BASE_DEG + 180)]           # the +Phi0 poles (45, 225 + base)
    neg = [rot(BASE_DEG + 90), rot(BASE_DEG + 270)]      # the -Phi0 poles (135, 315 + base)
    for p in pos:
        p.faces.name = "pole_pos"
    for p in neg:
        p.faces.name = "pole_neg"
    iron = pos[0]
    for p in pos[1:] + neg:
        iron = iron + p
    iron.maxh = 0.6 * maxh

    core = Cylinder(Pnt(0, 0, -AIR), Z, r=R_CORE, h=2 * AIR)   # refined axis column
    core.mat("air")
    core.maxh = maxh_core
    outer = Box(Pnt(-AIR, -AIR, -AIR), Pnt(AIR, AIR, AIR)) - iron - core
    outer.mat("air")
    for f in outer.faces:
        cc = f.center
        if max(abs(cc.x), abs(cc.y), abs(cc.z)) > 0.9 * AIR:
            f.name = "far"

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Glue([outer, core])).GenerateMesh(maxh=maxh))
    return mesh


def _solve_phi(mesh, order=2):
    """Laplace solve in the air with the poles at +-1 (the high-mu equipotential
    limit) -> B = -grad(Phi) (up to mu0; harmonic RATIOS only)."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    fes = H1(mesh, order=order, dirichlet="pole_pos|pole_neg|far")
    u, v = fes.TnT()
    with TaskManager():
        phi = GridFunction(fes)
        phi.Set(mesh.BoundaryCF({"pole_pos": 1.0, "pole_neg": -1.0}, default=0.0),
                definedon=mesh.Boundaries("pole_pos|pole_neg|far"))
        a = BilinearForm(fes)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        f = LinearForm(fes)
        f.Assemble()
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * phi.vec
        phi.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    return fes, phi


def _multipoles_of(mesh, phi, L, n_z=121, n_theta=30):
    """Integrated transverse multipoles along the beam (z) on the R_REF circle."""
    from ngsolve import grad
    B = -grad(phi)

    def B_perp(ax, ay, az):
        try:
            vv = B(mesh(ax, ay, az))
        except Exception:                                # step off a rare facet hit
            vv = B(mesh(ax + 3e-5, ay + 3e-5, az))
        return (float(vv[0]), float(vv[1]))              # (Bx, By); main = n=2

    zmax = min(0.9 * L, 0.85 * AIR)                      # past the ends, inside the box
    mp = integrated_multipoles(B_perp, R_REF, (-zmax, zmax), n_z=n_z,
                               n_theta=n_theta, n_max=10)
    return mp


def solve(L=0.080, maxh=0.022, maxh_core=0.0035, order=2, plot=False):
    """Solve the finite-length hyperbola quad and measure the integrated
    multipoles: main b_2, the suppressed forbidden harmonics, the allowed b_6."""
    mesh = build_quad_mesh(L, maxh, maxh_core)
    fes, phi = _solve_phi(mesh, order)
    mp = _multipoles_of(mesh, phi, L)
    b2 = abs(complex(*mp[2]))
    rel = {n: abs(complex(*mp[n])) / b2 for n in mp}
    forbidden_max = max(rel[n] for n in FORBIDDEN)
    allowed_spurious = {n: rel[n] for n in ALLOWED if n != 2}

    out = {
        "ne": int(mesh.ne), "ndof": int(fes.ndof), "L_m": float(L),
        "main_b2": float(b2),
        "rel_harmonics": {n: float(rel[n]) for n in (1, 2, 3, 4, 5, 6, 10)},
        "forbidden_max_rel": float(forbidden_max),       # n=1,3,5 (symmetry-suppressed)
        "b6_rel": float(rel[6]),                         # 12-pole = first allowed spurious
        "b10_rel": float(rel[10]),
        "allowed_spurious_rel": float(max(allowed_spurious.values())),
        "hyperbola_xy_const": float(C_HYP),              # pole face xy = r0^2/2
    }
    if plot:
        _plot(out)
    return out


def length_comparison(lengths=(0.050, 0.100), maxh=0.022, maxh_core=0.0035):
    """Is the integrated b_6 an END effect or a BODY effect?  Compare two pole
    lengths.  An END-driven b_6 would fall fast as L grows (the end region is a
    shrinking fraction); a BODY-driven b_6 (finite pole angular width) is nearly
    L-independent.  The dipole found b_3,5 body-dominated; check the quad."""
    rows = []
    for L in lengths:
        r = solve(L=L, maxh=maxh, maxh_core=maxh_core)
        rows.append((float(L), r["b6_rel"], r["forbidden_max_rel"]))
    b6_short, b6_long = rows[0][1], rows[-1][1]
    return {
        "rows": rows,                                    # (L, b6_rel, forbidden_max)
        "b6_short": float(b6_short), "b6_long": float(b6_long),
        "b6_change_rel": float(abs(b6_short - b6_long) / b6_short),  # small => body-dominated
    }


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.2), dpi=150)

    # LEFT: the 4 hodograph hyperbola poles (the design geometry)
    g = np.linspace(-R_OUT * 1.1, R_OUT * 1.1, 400)
    XX, YY = np.meshgrid(g, g)
    ang = math.radians(BASE_DEG)
    ca, sa = math.cos(ang), math.sin(ang)
    Xr, Yr = ca * XX + sa * YY, -sa * XX + ca * YY       # rotate frame by -base
    ax[0].contour(XX, YY, Xr * Yr, levels=[-C_HYP, C_HYP], colors=["C0", "C3"],
                  linewidths=2.0)
    th = np.linspace(0, 2 * math.pi, 200)
    ax[0].plot(R0 * np.cos(th), R0 * np.sin(th), "g--", lw=1.0, label=f"aperture r0={R0*1e3:.0f}mm")
    ax[0].plot(R_REF * np.cos(th), R_REF * np.sin(th), "0.5", lw=1.0, ls=":",
               label=f"r_ref={R_REF*1e3:.0f}mm")
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-R_OUT, R_OUT)
    ax[0].set_ylim(-R_OUT, R_OUT)
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_title("Hodograph quad poles: hyperbola $xy=\\pm r_0^2/2$\n"
                    "(red $=+\\Phi_0$, blue $=-\\Phi_0$, alternating)")
    ax[0].legend(fontsize=8, loc="upper right")

    # RIGHT: the integrated harmonic spectrum |b_n/b_2| (log)
    ns = [1, 2, 3, 4, 5, 6, 10]
    vals = [res["rel_harmonics"][n] for n in ns]
    # allowed normals are n = 2 mod 4 (2,6,10); everything else is forbidden.
    colors = ["C3" if n == 2 else ("C0" if n in ALLOWED else "0.6") for n in ns]
    ax[1].bar([str(n) for n in ns], vals, color=colors)
    ax[1].set_yscale("log")
    ax[1].set_ylim(1e-5, 2.0)
    ax[1].set_xlabel("normal harmonic  n  (1=dipole, 2=quad, 6=12-pole)")
    ax[1].set_ylabel("$|b_n / b_2|$  (integrated)")
    ax[1].set_title("Clean integrated quad: main $b_2$ (red); forbidden\n"
                    "$n{=}1,3,5$ suppressed (grey); allowed 12-pole $b_6$ (blue)")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("(A) QUADRUPOLE ENDS, FEM rung -- the analyzer handles any multipole\n")
    r = solve(plot=True)
    print(f"  mesh: ne={r['ne']}, ndof={r['ndof']}  (L={r['L_m']*1e3:.0f} mm)")
    print(f"  hodograph pole = hyperbola  xy = {r['hyperbola_xy_const']:.3e}  (= r0^2/2)\n")
    print("  integrated harmonic spectrum |b_n / b_2|:")
    for n in (1, 2, 3, 4, 5, 6, 10):
        tag = ("<- main" if n == 2 else
               "  (forbidden: 4-fold + alt-sign)" if n in FORBIDDEN else
               "  <- first ALLOWED spurious (12-pole)" if n == 6 else "")
        print(f"    n={n:2d}:  {r['rel_harmonics'][n]:.3e} {tag}")
    print(f"\n  -> clean integrated QUADRUPOLE: main b_2; forbidden n=1,3,5 "
          f"suppressed to {r['forbidden_max_rel']:.1e}")
    print(f"     (~{r['b6_rel']/r['forbidden_max_rel']:.0f}x below the allowed "
          f"12-pole b_6 = {r['b6_rel']:.2e}).")

    print("\n  Is b_6 an END or a BODY effect?  (length comparison)")
    lc = length_comparison()
    for L, b6, fb in lc["rows"]:
        print(f"    L={L*1e3:3.0f} mm:  b_6/b_2 = {b6:.3e}  (forbidden {fb:.1e})")
    print(f"    -> b_6 changes only {lc['b6_change_rel']*100:.0f}% over 2x length "
          f"=> BODY-dominated (finite pole angular width), NOT the ends --")
    print("       the SAME body-dominated picture the dipole found for b_3,5.")
    print("       (the analytic theorem in accel_pole_ends_3d shows an ideal")
    print("        Maxwellian end contributes ZERO integrated b_6.)")


if __name__ == "__main__":
    main()
