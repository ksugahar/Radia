r"""3-D Clebsch pole-face shape OPTIMIZATION: null b3 AND b5 together.

The accelerator dipole pole face is a magnetic-scalar EQUIPOTENTIAL (high mu => the
tangential field vanishes => Phi = const on the iron) -- i.e. a CLEBSCH LEVEL SET.  The
3-D body pole is the translationally-invariant extrusion of its 2-D cross-section
contour, so OPTIMIZING THE 2-D CONTOUR IS optimizing the 3-D pole's Clebsch level set
(`accel_pole_dipole_body_2d.py` established that the integrated transverse spurious
harmonics b3, b5 are a BODY/pole-shape lever, not an end effect).

A finite flat pole droops at its edges, so the midplane field expands as

    B_z(x) = b1 + b3 x^2 + b5 x^4 + ...        (even in x by dipole symmetry).

A SINGLE quadratic shim (`accel_pole_dipole_body_2d`) has ONE knob, so it can null b3
but it leaves b5.  Here a TWO-parameter Clebsch contour

    z_face(x) = g/2 - d2 (x/w)^2 - d4 (x/w)^4

gives TWO knobs, optimized by a 2-D Newton on the residual (b3(d2,d4), b5(d2,d4)) = 0
(the harmonics respond near-linearly to the shim coefficients, so the Jacobian solve
lands close and one refinement tightens it).  Both leading spurious harmonics are
nulled SIMULTANEOUSLY -- a genuine multi-parameter shape optimization, each evaluation
a 2-D Laplace solve in the high-mu equipotential limit, verified by the harmonic
content of the optimized pole.

run:  python clebsch_pole_shape_optimization_2d.py   [--fig]
"""

from _validation_output import validation_output

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from accel_pole_design import multipoles                                # noqa: E402

# ---- geometry (meters); cross-section is (x = width, z = gap) ----
GAP = 0.040            # pole-face separation -> half-gap g/2
G2 = GAP / 2.0
POLE_T = 0.040         # pole thickness (z extent of iron above the face)
X_AIR = 0.24           # air box half-width (x)
Z_AIR = 0.18           # air box height (z, upper half only)
R_REF = 0.008          # reference-circle radius (inside the gap)
HALF_W = 0.040         # pole half-width
ALLOWED = (1, 3, 5)    # dipole-allowed normal harmonics


def face_z(x, d2, d4, half_w=HALF_W):
    """Two-parameter Clebsch pole-face contour: flat g/2, narrowed by a quadratic
    (d2) and a quartic (d4) shim toward the edges."""
    xr = np.asarray(x) / half_w
    return G2 - d2 * xr ** 2 - d4 * xr ** 4


def field_harmonics(d2, d4, half_w=HALF_W, order=3, maxh=0.007, n_face=46):
    """2-D Laplace solve (high-mu equipotential pole), return the transverse multipoles
    measured on the reference circle: signed b3/b1, b5/b1, and the worst spurious."""
    import ngsolve as ng
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx,
                         TaskManager)
    from netgen.occ import WorkPlane, OCCGeometry

    xs = np.linspace(-half_w, half_w, n_face)
    zf = face_z(xs, d2, d4, half_w)
    box = WorkPlane().MoveTo(-X_AIR, 0.0).Rectangle(2 * X_AIR, Z_AIR).Face()
    wp = WorkPlane().MoveTo(float(xs[0]), float(zf[0]))
    for xi, zi in zip(xs[1:], zf[1:]):
        wp.LineTo(float(xi), float(zi))
    wp.LineTo(float(half_w), G2 + POLE_T)
    wp.LineTo(float(-half_w), G2 + POLE_T)
    wp.Close()
    air = box - wp.Face()
    for e in air.edges:
        c = e.center
        if abs(c.y) < 1e-7:
            e.name = "mid"                       # midplane z = 0 (Omega = 0)
        elif abs(c.x) > X_AIR - 1e-7 or c.y > Z_AIR - 1e-7:
            e.name = "far"                       # outer box (Omega = 0)
        else:
            e.name = "pole"                      # pole face (Omega = 1)
    air.maxh = maxh

    with TaskManager():
        mesh = ng.Mesh(OCCGeometry(air, dim=2).GenerateMesh(maxh=maxh))
        fes = H1(mesh, order=order, dirichlet="mid|pole|far")
        u, v = fes.TnT()
        omega = GridFunction(fes)
        omega.Set(mesh.BoundaryCF({"pole": 1.0}, default=0.0),
                  definedon=mesh.Boundaries("mid|pole|far"))
        a = BilinearForm(grad(u) * grad(v) * dx)
        a.Assemble()
        f = LinearForm(fes)
        f.Assemble()
        r = f.vec.CreateVector()
        r.data = f.vec - a.mat * omega.vec
        omega.vec.data += a.mat.Inverse(fes.FreeDofs()) * r

    Bcf = -grad(omega)

    def B_perp(ax, ay):
        if ay >= 0:
            vv = Bcf(mesh(ax, ay))
            return (float(vv[0]), float(vv[1]))
        vv = Bcf(mesh(ax, -ay))                  # dipole antisymmetry
        return (-float(vv[0]), float(vv[1]))

    mp = multipoles(B_perp, R_REF, n_max=8, n_samples=512)
    b1 = complex(*mp[1]).real
    main = abs(complex(*mp[1]))
    return {
        "b3": complex(*mp[3]).real / b1,
        "b5": complex(*mp[5]).real / b1,
        "spurious": float(max(abs(complex(*mp[n])) / main for n in ALLOWED if n != 1)),
        "ne": int(mesh.ne),
    }


def optimize_shape(half_w=HALF_W, order=3, maxh=0.007, h_fd=2e-4, n_refine=1):
    """Optimize (d2, d4) to null (b3, b5) by a 2-D Newton (Jacobian by finite diff at
    the origin, then n_refine refinements).  Returns the flat / 1-param / 2-param
    harmonics and the optimum shim."""
    def resid(d2, d4):
        h = field_harmonics(d2, d4, half_w, order, maxh)
        return np.array([h["b3"], h["b5"]]), h

    r0, flat = resid(0.0, 0.0)
    r2, _ = resid(h_fd, 0.0)
    r4, _ = resid(0.0, h_fd)
    J = np.column_stack([(r2 - r0) / h_fd, (r4 - r0) / h_fd])   # d(b3,b5)/d(d2,d4)
    step = np.linalg.solve(J, -r0)
    d2, d4 = float(step[0]), float(step[1])
    r, opt = resid(d2, d4)
    for _ in range(n_refine):                                  # Newton refinement (reuse J)
        step = np.linalg.solve(J, -r)
        d2 += float(step[0]); d4 += float(step[1])
        r, opt = resid(d2, d4)

    # the single-shim baseline: d2 alone nulls b3 (linearised), d4 = 0.
    d2_1 = float(-r0[0] / J[0, 0])
    _, one = resid(d2_1, 0.0)
    return {
        "flat": flat, "one_param": one, "two_param": opt,
        "d2_opt": d2, "d4_opt": d4, "d2_one": d2_1,
        "jacobian": J.ravel().tolist(),
    }


def run(order=3, maxh=0.007):
    out = optimize_shape(order=order, maxh=maxh)
    return out


def _make_figure(out, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    xs = np.linspace(-HALF_W, HALF_W, 200)
    axL.plot(xs * 1e3, face_z(xs, 0, 0) * 1e3, "-", label="flat")
    axL.plot(xs * 1e3, face_z(xs, out["d2_one"], 0) * 1e3, "--", label="1 shim (b3)")
    axL.plot(xs * 1e3, face_z(xs, out["d2_opt"], out["d4_opt"]) * 1e3, "-",
             label="2 shim (b3,b5)")
    axL.set_xlabel("x [mm]"); axL.set_ylabel("pole face z [mm]")
    axL.legend(frameon=False, fontsize=9)
    cases = [("flat", out["flat"]), ("1 shim", out["one_param"]), ("2 shim", out["two_param"])]
    x = np.arange(3); w = 0.36
    axR.bar(x - w / 2, [abs(c[1]["b3"]) for c in cases], w, label="|b3/b1|", color="#0072B2")
    axR.bar(x + w / 2, [abs(c[1]["b5"]) for c in cases], w, label="|b5/b1|", color="#D55E00")
    axR.set_yscale("log"); axR.set_xticks(x); axR.set_xticklabels([c[0] for c in cases])
    axR.set_ylabel("spurious harmonic")
    axR.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return True


def main():
    import json
    print("=" * 70)
    print("3-D Clebsch pole-face shape optimization: null b3 AND b5")
    print("=" * 70)
    out = run()
    for name in ("flat", "one_param", "two_param"):
        h = out[name]
        print("  %-10s b3/b1 = %+.2e  b5/b1 = %+.2e  spurious = %.2e"
              % (name, h["b3"], h["b5"], h["spurious"]))
    print("  optimum shim: d2 = %.3f mm, d4 = %.3f mm  (1-shim d2 = %.3f mm)"
          % (out["d2_opt"] * 1e3, out["d4_opt"] * 1e3, out["d2_one"] * 1e3))

    here = os.path.dirname(os.path.abspath(__file__))
    res = {k: (v if not isinstance(v, dict) else v) for k, v in out.items()}
    with validation_output("clebsch_pole_shape_optimization_2d.json").open("w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nsaved clebsch_pole_shape_optimization_2d.json")
    if "--fig" in sys.argv:
        if _make_figure(out, os.path.join(here, "clebsch_pole_shape_optimization_2d.png")):
            print("saved clebsch_pole_shape_optimization_2d.png")


if __name__ == "__main__":
    main()
