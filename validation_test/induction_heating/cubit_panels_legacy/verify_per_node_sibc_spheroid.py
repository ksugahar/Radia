"""Per-node Z_s SIBC validation on a PROLATE SPHEROID.

The sphere validation in ``verify_per_node_sibc_sphere.py`` only proves
the per-node path WORKS — every panel sees the same R, so any uniform
ndarray would have done. The real test is a workpiece with **spatially
varying curvature**, where the per-node R must change panel by panel.

A prolate spheroid (cigar) ``x^2/b^2 + y^2/b^2 + z^2/a^2 = 1`` with
``a > b`` has principal radii of curvature at parametric angle
``t`` (measured from the +z polar axis) given by::

    rho = sqrt(a^2 sin^2 t + b^2 cos^2 t)
    R_meridional  = rho^3 / (a*b)
    R_azimuthal   = rho      * b / a   (   b * rho / a   )

derived from the standard Gauss / second-fundamental-form result for
a surface of revolution. At the equator (t = pi/2): R_m = a^2/b,
R_a = b. At the poles (t = 0): R_m = R_a = b^2/a (the "tip" radius).

The percentile-10 estimator in ``_compute_panel_local_radii`` picks
the **smaller** principal radius (i.e. the maximum local curvature
direction). For a prolate (a > b) cigar:

    smaller-of-two:    pole = b^2/a   equator = b
    extractor expects: pole = b^2/a   equator = b

These should agree to within mesh discretization error.

Three SIBC paths are compared:

  1. **scalar(R = mean)** — the legacy path, where the user supplies
     a single global half_thickness equal to the mean principal R
     (the most-favorable scalar choice).
  2. **per-node(mesh R)** — the production path; per-panel R from the
     mesh extractor, projected to H1 nodes by vertex averaging.
  3. **per-node(analytical R)** — the GROUND TRUTH; per-panel R from
     the analytical spheroid principal-curvature formula evaluated at
     each panel centroid. No mesh discretization error.

If (2) and (3) agree to a few percent, and both differ from (1) in
the regime where Z_s depends on R, the per-node implementation is
verified.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src", "radia"))
sys.path.insert(0, os.path.join(REPO_ROOT, "src", "radia", "panels"))

from ngsolve import Mesh, BND
from ngsolve import TaskManager
from netgen.occ import (Sphere, Pnt, OCCGeometry, Glue, Axes, Dir,
                         Vec)

from radia.bem_sibc_solver import ScalarBIESIBCSolver
from calc_inductance import (_extract_panels_from_mesh,
                              _compute_panel_local_radii,
                              _build_per_node_Zs, _dowell_Zs)

MU_0 = 4e-7 * math.pi


# ---------------------------------------------------------------------
# Analytical spheroid principal radii of curvature
# ---------------------------------------------------------------------

def spheroid_principal_radii(x, y, z, a, b):
    """Analytical principal radii of curvature on a prolate spheroid.

    Surface: x^2/b^2 + y^2/b^2 + z^2/a^2 = 1, a > b > 0.

    Returns ``(R_meridional, R_azimuthal)`` at the surface point
    ``(x, y, z)``. Sign convention: both positive (convex outward).

    Derivation: parameterise as x = b sin(t) cos(phi), y = b sin(t)
    sin(phi), z = a cos(t). The standard surface-of-revolution
    formulas give::

        rho^2(t)      = a^2 sin^2(t) + b^2 cos^2(t)
        R_meridional  = rho^3 / (a * b)
        R_azimuthal   = b * rho / a
    """
    r_xy = math.sqrt(x * x + y * y)
    cos_t = z / a
    cos_t = max(-1.0, min(1.0, cos_t))
    sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
    rho2 = a * a * sin_t * sin_t + b * b * cos_t * cos_t
    rho = math.sqrt(rho2)
    R_m = rho ** 3 / (a * b)
    R_a = b * rho / a
    return R_m, R_a


def build_prolate_spheroid_mesh(a, b, maxh):
    """Surface mesh of x^2/b^2 + y^2/b^2 + z^2/a^2 = 1.

    Approach: build a fine **polyline approximation** of the half
    ellipse (z from -a to +a, r = b * sqrt(1 - z^2/a^2)) in the
    (x, z) plane via WorkPlane.Spline, then revolve 360 around z.
    Revolving a face whose boundary touches the rotation axis only
    at single points works robustly with this construction.
    """
    from netgen.occ import WorkPlane, Axis

    # Half-ellipse polyline: from (0, +a) along the right side of
    # the ellipse to (0, -a). Use 64 sample points for smoothness.
    n = 64
    pts = []
    for i in range(n + 1):
        t = math.pi * i / n  # 0 .. pi
        x = b * math.sin(t)
        z = a * math.cos(t)
        pts.append((x, z))
    # Workplane in the (x, z) plane: origin at (0, 0, 0), normal +y,
    # local x-axis along world +x, local y-axis along world +z.
    from netgen.occ import gp_Pnt2d

    wp = WorkPlane(Axes(p=Pnt(0, 0, 0), n=Dir(0, 1, 0), h=Dir(1, 0, 0)))
    wp = wp.MoveTo(pts[0][0], pts[0][1])
    # Spline through the half-ellipse points (2D in the workplane).
    wp = wp.Spline([gp_Pnt2d(p[0], p[1]) for p in pts[1:]])
    wp = wp.LineTo(pts[0][0], pts[0][1])  # close: along the z-axis
    half_face = wp.Face()
    spheroid = half_face.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    for f in spheroid.faces:
        f.name = "sibc"
    geo = OCCGeometry(Glue(spheroid.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    with TaskManager():
        mesh.Curve(2)
        return mesh


def uniform_field_phi_inc(mesh, H0):
    """phi_inc = -H0 * z (uniform H = H0 z_hat)."""
    nodes = np.array([mesh.vertices[i].point for i in range(mesh.nv)],
                     dtype=float)
    return (-H0 * nodes[:, 2]).astype(complex)


def analytical_R_per_panel(panels, a, b):
    """Analytical SMALLER principal R at each panel centroid."""
    R = np.empty(len(panels))
    for i, p in enumerate(panels):
        x, y, z = p["center"]
        R_m, R_a = spheroid_principal_radii(x, y, z, a, b)
        R[i] = min(R_m, R_a)
    return R


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------

def run_case(a, b, sigma, mu_r, frequency, B0, maxh_factor=3):
    """Compare scalar vs per-node-mesh vs per-node-analytical."""
    omega = 2 * math.pi * frequency
    mu_eff = MU_0 * mu_r
    rho = 1.0 / sigma
    delta = math.sqrt(2 * rho / (omega * mu_eff))

    print(f"--- prolate a={a*1e3:.1f}mm b={b*1e3:.1f}mm  "
          f"sigma={sigma:.2e}  mu_r={mu_r}  f={frequency} Hz ---")
    print(f"  delta = {delta*1e3:.4f} mm")
    print(f"  pole R (b^2/a) = {b*b/a*1e3:.3f} mm  "
          f"equator b = {b*1e3:.3f} mm")
    print(f"  ratio        = {b/(b*b/a):.2f}x  "
          f"(this is the spread of R across the surface)")
    print()

    mesh = build_prolate_spheroid_mesh(a, b, b / maxh_factor)
    nse = mesh.GetNE(BND)
    nv = mesh.nv
    print(f"  mesh: nv={nv}, nse={nse}")

    panels = _extract_panels_from_mesh(mesh, "sibc")
    n_panels = len(panels)
    if n_panels != nse:
        print(f"  WARN: extracted {n_panels} panels but nse={nse}")

    # ---- Step 1: extractor vs analytical principal R ----
    R_mesh = _compute_panel_local_radii(mesh, panels, default_R=10 * a)
    R_analytic = analytical_R_per_panel(panels, a, b)
    rel_err = (R_mesh - R_analytic) / R_analytic
    print(f"  R_mesh:     min={R_mesh.min()*1e3:.3f}  "
          f"mean={R_mesh.mean()*1e3:.3f}  "
          f"max={R_mesh.max()*1e3:.3f} mm")
    print(f"  R_analytic: min={R_analytic.min()*1e3:.3f}  "
          f"mean={R_analytic.mean()*1e3:.3f}  "
          f"max={R_analytic.max()*1e3:.3f} mm")
    print(f"  R rel err:  mean={rel_err.mean()*100:+.1f}%  "
          f"std={rel_err.std()*100:.1f}%  "
          f"max|.|={np.max(np.abs(rel_err))*100:.1f}%")
    print()

    # ---- Step 2: SIBC solver paths ----
    solver = ScalarBIESIBCSolver(mesh, order=1)
    H0 = B0 / MU_0
    phi_inc = uniform_field_phi_inc(mesh, H0)

    # Path A: scalar Z_s using mean analytical R as the "best scalar"
    R_scalar = float(R_analytic.mean())
    Z_s_scalar = _dowell_Zs(R_scalar, rho, omega, mu_eff)
    sol_scalar = solver.solve(phi_inc, Z_s=Z_s_scalar, omega=omega)

    # Path B: per-node Z_s from mesh extractor
    Zs_node_mesh, _ = _build_per_node_Zs(
        mesh, panels, R_mesh, rho, omega, mu_eff, solver.ndof)
    sol_node_mesh = solver.solve(phi_inc, Z_s=Zs_node_mesh, omega=omega)

    # Path C: per-node Z_s from ANALYTICAL R per panel (ground truth)
    Zs_node_ana, _ = _build_per_node_Zs(
        mesh, panels, R_analytic, rho, omega, mu_eff, solver.ndof)
    sol_node_ana = solver.solve(phi_inc, Z_s=Zs_node_ana, omega=omega)

    # H_t_rms comparison
    Hrms_s = sol_scalar["H_t_rms"]
    Hrms_m = sol_node_mesh["H_t_rms"]
    Hrms_a = sol_node_ana["H_t_rms"]
    print(f"  {'method':<28} {'H_t_rms':>14} {'vs analytical':>16}")
    print(f"  {'-'*28} {'-'*14} {'-'*16}")
    print(f"  {'scalar(R = mean R)':<28} {Hrms_s:>14.4f} "
          f"{(Hrms_s - Hrms_a) / Hrms_a * 100:>+15.2f}%")
    print(f"  {'per-node(mesh extractor)':<28} {Hrms_m:>14.4f} "
          f"{(Hrms_m - Hrms_a) / Hrms_a * 100:>+15.2f}%")
    print(f"  {'per-node(analytical R)':<28} {Hrms_a:>14.4f} "
          f"{'(ground truth)':>16}")
    print()

    return {
        "a": a, "b": b, "frequency": frequency,
        "Hrms_scalar": Hrms_s,
        "Hrms_per_node_mesh": Hrms_m,
        "Hrms_per_node_ana": Hrms_a,
        "R_extractor_mean_err_pct": float(rel_err.mean() * 100),
        "R_extractor_max_err_pct": float(np.max(np.abs(rel_err)) * 100),
    }


def main():
    print("=" * 72)
    print("Prolate spheroid SIBC validation: per-node R extractor + solver")
    print("=" * 72)
    print()

    cases = [
        # (a, b, sigma, mu_r, frequency, B0)
        # Cu cigar (a/b = 2), 10 Hz: xi << 1, R-sensitive
        (0.020, 0.010, 5.8e7, 1, 10, 1.0),
        # Cu cigar (a/b = 2), 1 MHz: xi >> 1, asymptotic
        (0.020, 0.010, 5.8e7, 1, 1e6, 1.0),
        # More elongated cigar (a/b = 4) at low freq — bigger spread
        # of principal R, stronger per-panel test
        (0.040, 0.010, 5.8e7, 1, 100, 1.0),
    ]

    results = []
    for c in cases:
        results.append(run_case(*c))

    print("=" * 72)
    print("Summary")
    print("=" * 72)
    print(f"  {'case':<46} {'mesh err':>11} {'mesh-vs-ana':>13}")
    for r in results:
        label = f"a={r['a']*1e3:.0f}mm b={r['b']*1e3:.0f}mm f={r['frequency']:.0f}Hz"
        ratio = (r['Hrms_per_node_mesh'] - r['Hrms_per_node_ana']) \
            / r['Hrms_per_node_ana'] * 100
        print(f"  {label:<46} {r['R_extractor_mean_err_pct']:>+10.1f}% "
              f"{ratio:>+12.2f}%")

    print()
    print("Pass criteria:")
    print("  1. extractor mean error < 15% (chord vs arc, mesh refinement)")
    print("  2. per-node(mesh) within 5% of per-node(analytical)")
    print("  3. per-node beats scalar in low-freq / xi~1 regime")


if __name__ == "__main__":
    main()
