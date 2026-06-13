"""(A) The dipole BODY lever: transverse b_3,5 = pole-shape, not ends.

RESEARCH example (track A).  The 3-D end study (accel_pole_ends_fem.py) found
the honest TWO-LEVER split: the pole-end chamfer drives the LONGITUDINAL end
bump, but the INTEGRATED TRANSVERSE spurious harmonics b_3, b_5 are
body/pole-width dominated and the END chamfer does not move them.  This script
isolates and exercises the OTHER lever -- the BODY pole-face shape -- in the
2-D cross-section (the body is translationally invariant along the beam, so
its cross-section has NO ends; this is the clean body model).

Physics = the same design principle as accel_pole_harmonics.py but for the
DIPOLE and in the *fix* direction.  The iron pole face is a magnetic-scalar
equipotential (high mu -> H_t = 0 -> Phi = const).  An infinite flat pole
gives a perfect dipole; a FINITE flat pole droops at its edges, so the
midplane field expands as

    B_z(x) = b_1 + b_3 x^2 + b_5 x^4 + ...        (even in x by symmetry)

with b_3 < 0 (the field falls off away from centre).  Two body levers fix it:

  (1) WIDTH   -- a wider pole flattens the field over a fixed aperture
                 -> |b_3|, |b_5| fall.
  (2) CURVATURE (the dipole "shim"/Rogowski analog) -- narrow the gap slightly
     toward the pole edges, z_face(x) = g/2 - delta (x/w)^2, to BOOST the field
     where the flat pole drooped.  Tuning delta drives b_3 THROUGH ZERO -- the
     finite-aperture analog of the quad hyperbola: the ideal finite-width
     dipole pole is slightly concave, not flat.

The field is a real FEM Laplace solve in the air gap with the pole face held
at a fixed potential (the high-mu equipotential limit) and the midplane at 0
(the dipole up-down antisymmetry).  The transverse multipoles are measured on
a reference circle with the SAME analyzer as the 3-D study
(accel_pole_design.multipoles); the lower half of the circle is reconstructed
from the upper half by the dipole antisymmetry B_z(x,-z)=B_z(x,z),
B_x(x,-z)=-B_x(x,z).

run:  python accel_pole_dipole_body_2d.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the analyzer
from accel_pole_design import multipoles                          # noqa: E402

# ---- geometry (meters); cross-section is (x = width, z = gap) ----
GAP = 0.040            # pole-face separation -> half-gap g/2 = 0.020
G2 = GAP / 2.0
POLE_T = 0.040         # pole thickness (z extent of the iron above the face)
X_AIR = 0.24           # air box half-width (x)
Z_AIR = 0.18           # air box height   (z, upper half only)
R_REF = 0.008          # reference-circle radius (inside the gap)

ALLOWED = (1, 3, 5)    # dipole-allowed normal harmonics


def _pole_face_z(x, half_w, delta):
    """Pole-face profile z_face(x): flat g/2, narrowed by delta (x/half_w)^2
    toward the edges (delta > 0 = the compensating concave 'shim')."""
    return G2 - delta * (x / half_w) ** 2


def build_mesh(half_w=0.040, delta=0.0, maxh=0.006, n_face=40):
    """Air gap (upper half) with a floating iron pole island whose bottom face
    is the (possibly curved) equipotential.  netgen.occ 2-D, no Cubit."""
    import ngsolve as ng
    from netgen.occ import WorkPlane, OCCGeometry

    # outer air box: x in [-X_AIR, X_AIR], z in [0, Z_AIR]
    box = WorkPlane().MoveTo(-X_AIR, 0.0).Rectangle(2 * X_AIR, Z_AIR).Face()

    # pole island: shaped bottom face (left->right), up the right, across the
    # top, down the left (Close).
    xs = np.linspace(-half_w, half_w, n_face)
    zf = _pole_face_z(xs, half_w, delta)
    z_top = G2 + POLE_T
    wp = WorkPlane().MoveTo(float(xs[0]), float(zf[0]))
    for xi, zi in zip(xs[1:], zf[1:]):
        wp.LineTo(float(xi), float(zi))
    wp.LineTo(float(half_w), z_top)
    wp.LineTo(float(-half_w), z_top)
    wp.Close()
    pole = wp.Face()

    air = box - pole
    for e in air.edges:
        c = e.center
        if abs(c.y) < 1e-7:
            e.name = "mid"                       # midplane z = 0 (Omega = 0)
        elif abs(c.x) > X_AIR - 1e-7 or c.y > Z_AIR - 1e-7:
            e.name = "far"                       # outer box (Omega = 0)
        else:
            e.name = "pole"                      # pole face (Omega = 1)
    air.maxh = maxh

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air, dim=2).GenerateMesh(maxh=maxh))
    return mesh


def _solve_field(half_w, delta, order=3, maxh=0.006):
    """Laplace solve; return the multipole dict measured on the R_REF circle."""
    import ngsolve as ng
    from ngsolve import H1, BilinearForm, LinearForm, GridFunction, grad, dx, TaskManager

    mesh = build_mesh(half_w, delta, maxh)
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

    Bcf = -grad(omega)                           # (Bx, Bz) up to mu0 (ratios only)

    def B_perp(ax, ay):                          # analyzer (x, y) = phys (x, z)
        if ay >= 0:
            vv = Bcf(mesh(ax, ay))
            return (float(vv[0]), float(vv[1]))  # (Bx, Bz=main)
        vv = Bcf(mesh(ax, -ay))                  # reflect: dipole antisymmetry
        return (-float(vv[0]), float(vv[1]))

    mp = multipoles(B_perp, R_REF, n_max=8, n_samples=512)
    main = abs(complex(*mp[1]))
    b3 = complex(*mp[3]).real / (complex(*mp[1]).real)   # signed b_3/b_1
    b5 = complex(*mp[5]).real / (complex(*mp[1]).real)
    spurious = max(abs(complex(*mp[n])) / main for n in ALLOWED if n != 1)
    return {"ne": int(mesh.ne), "main_b1": float(main),
            "b3_over_b1": float(b3), "b5_over_b1": float(b5),
            "spurious_rel": float(spurious)}


def solve(order=3, maxh=0.006):
    """Exercise the two body levers; return the verification dict."""
    # ---- (1) WIDTH lever: wider pole -> flatter field -> smaller |b_3| ----
    widths = (0.030, 0.045, 0.060)
    width_rows = []
    for w in widths:
        r = _solve_field(w, 0.0, order, maxh)
        width_rows.append((w, r["b3_over_b1"], r["spurious_rel"]))
    narrow_spur = width_rows[0][2]
    wide_spur = width_rows[-1][2]

    # ---- (2) CURVATURE lever at fixed width: tune delta -> b_3 through zero ----
    half_w = 0.040
    deltas = (0.0, 0.0006, 0.0012, 0.0018, 0.0024)
    curv_rows = []
    for d in deltas:
        r = _solve_field(half_w, d, order, maxh)
        curv_rows.append((d, r["b3_over_b1"], r["spurious_rel"]))
    b3_flat = curv_rows[0][1]
    b3_arr = np.array([row[1] for row in curv_rows])
    d_arr = np.array([row[0] for row in curv_rows])
    # b3 is monotone decreasing through zero as delta grows: find the crossing.
    if b3_arr[0] > 0 > b3_arr[-1] or b3_arr[0] < 0 < b3_arr[-1]:
        order_idx = np.argsort(b3_arr)
        delta_opt = float(np.interp(0.0, b3_arr[order_idx], d_arr[order_idx]))
    else:
        delta_opt = float("nan")

    # confirm AT the optimum: with b_3 ~ 0 the residual spurious is the
    # next harmonic b_5 -- the genuine field-quality the curvature lever buys.
    if delta_opt == delta_opt:
        conf = _solve_field(half_w, delta_opt, order, maxh)
        spur_at_opt = conf["spurious_rel"]
        b5_at_opt = conf["b5_over_b1"]
    else:
        spur_at_opt = float("nan")
        b5_at_opt = float("nan")

    return {
        "width_sweep": width_rows,           # (w, b3/b1, spurious) -- |b3| falls with w
        "narrow_spurious_rel": float(narrow_spur),
        "wide_spurious_rel": float(wide_spur),
        "curvature_sweep": curv_rows,        # (delta, b3/b1, spurious)
        "b3_flat": float(b3_flat),           # flat finite pole: b3 < 0 (droop)
        "delta_opt_m": delta_opt,            # the concavity that zeroes b3
        "spurious_at_opt_rel": float(spur_at_opt),   # residual quality at b3=0 (~|b5|)
        "b5_at_opt": float(b5_at_opt),
        "g2_m": float(G2), "r_ref_m": float(R_REF),
    }


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(10.2, 3.8), dpi=150)

    w = [row[0] * 1e3 for row in res["width_sweep"]]
    b3w = [row[1] * 1e4 for row in res["width_sweep"]]
    ax[0].plot(w, b3w, "o-", color="C0")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xlabel("pole half-width  w [mm]")
    ax[0].set_ylabel("$b_3/b_1$  [units, $\\times10^{-4}$]")
    ax[0].set_title("WIDTH lever: a wider flat pole flattens the field\n"
                    "(|$b_3$| falls; the finite pole droops at its edges)")

    d = [row[0] * 1e3 for row in res["curvature_sweep"]]
    b3d = [row[1] * 1e4 for row in res["curvature_sweep"]]
    ax[1].plot(d, b3d, "s-", color="C1")
    ax[1].axhline(0, color="k", lw=0.6)
    if res["delta_opt_m"] == res["delta_opt_m"]:
        ax[1].axvline(res["delta_opt_m"] * 1e3, color="C3", lw=1, ls=":",
                      label=f"$\\delta$={res['delta_opt_m']*1e3:.2f} mm ($b_3$=0)")
        ax[1].legend(fontsize=8)
    ax[1].set_xlabel("edge gap narrowing  $\\delta$ [mm]")
    ax[1].set_ylabel("$b_3/b_1$  [units, $\\times10^{-4}$]")
    ax[1].set_title("CURVATURE lever (the dipole shim): $z_{face}=g/2-\\delta(x/w)^2$\n"
                    "tuning the concavity drives $b_3$ through zero")

    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("(A) Dipole BODY lever -- transverse b_3,5 is pole-shape, not ends\n")
    r = solve()
    print("  (1) WIDTH lever (flat pole, fixed aperture r_ref="
          f"{r['r_ref_m']*1e3:.0f} mm):")
    for w, b3, sp in r["width_sweep"]:
        print(f"      half-width {w*1e3:4.0f} mm:  b_3/b_1 = {b3:+.3e}  "
              f"spurious |b_n/b_1| = {sp:.3e}")
    print(f"      -> a wider pole flattens the field: spurious "
          f"{r['narrow_spurious_rel']:.2e} -> {r['wide_spurious_rel']:.2e}\n")
    print("  (2) CURVATURE lever (the dipole shim, fixed width 40 mm):")
    for d, b3, sp in r["curvature_sweep"]:
        print(f"      delta {d*1e3:5.2f} mm:  b_3/b_1 = {b3:+.3e}  "
              f"spurious = {sp:.3e}")
    print(f"      flat pole b_3/b_1 = {r['b3_flat']:+.3e} (< 0: edge droop)")
    print(f"      -> concavity delta = {r['delta_opt_m']*1e3:.2f} mm zeroes b_3: "
          f"spurious {abs(r['b3_flat']):.2e} -> {r['spurious_at_opt_rel']:.2e} "
          f"(residual = |b_5/b_1| = {abs(r['b5_at_opt']):.2e})\n")
    print("  => the transverse b_3,5 are a BODY pole-shape lever (width +")
    print("     curvature), the lever the END chamfer cannot move -- the honest")
    print("     other half of the two-lever split in accel_pole_ends_fem.py.")
    _plot(r)


if __name__ == "__main__":
    main()
