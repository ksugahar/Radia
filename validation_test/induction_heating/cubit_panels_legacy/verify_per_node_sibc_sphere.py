"""Validate per-node Z_s SIBC against the analytical sphere SIBC.

The canonical analytical reference is the conducting sphere (radius R,
conductivity sigma, optionally permeability mu_r) immersed in a uniform
external magnetic field H_0 = H0 * z_hat. With surface impedance Z_s
applied on r=R, the surface tangential current is::

    J_s(theta) = (3/2) * j*omega*B0*R / (j*omega*mu_0*R + Z_s) * sin(theta)
    H_t(theta) = J_s(theta)
    H_t_rms    = J_s_max * sqrt(2/3)

This script feeds **the same uniform-field phi_inc** to two SIBC paths:

  1. **Scalar Z_s**: legacy path, uses the global Dowell tanh value
     evaluated at ``R = half_thickness`` (the user input).
  2. **Per-node Z_s**: per-node Dowell evaluated with the per-panel
     local R extracted from the mesh (this is the path Phase 5 wires
     into the coupled BEM solver).

Both should reproduce the analytical answer. The per-node path tests
that:

  - The mesh-driven local R extraction recovers the sphere radius.
  - ``ScalarBIESIBCSolver.solve(phi_inc, Z_s=ndarray, omega)`` is
    numerically correct.
  - The per-node and scalar paths agree to numerical precision when
    R_local is uniform (regression).
  - When the user supplies a WRONG global half_thickness, the per-node
    path still recovers the right answer (the headline result).

Reference: Smythe, "Static and Dynamic Electricity" 3rd ed, Section 11
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
from netgen.occ import Sphere, Pnt, OCCGeometry, Glue

from radia.bem_sibc_solver import ScalarBIESIBCSolver
from calc_inductance import (_extract_panels_from_mesh,
                              _compute_panel_local_radii,
                              _build_per_node_Zs, _dowell_Zs)

MU_0 = 4e-7 * math.pi


def analytical_sphere_sibc(R, Z_s, B0, omega):
    """H_t_rms for a conducting sphere in a uniform field with SIBC.

    From Smythe Section 11. Returns the RMS tangential field over the
    sphere surface (accounting for the sin^2 weight)::

        J_s(theta) = (3/2) * j*omega*B0*R / (j*omega*mu_0*R + Z_s) * sin(theta)
        H_t_rms    = max(|J_s|) * sqrt(2/3)
    """
    beta = 1j * omega * MU_0 * R
    J_s_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    H_t_rms = J_s_max * math.sqrt(2.0 / 3.0)
    return {"J_s_max": J_s_max, "H_t_rms": H_t_rms}


def build_sphere_mesh(R, maxh):
    """Surface-only mesh of a sphere with the BND label "sibc"."""
    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sibc"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    with TaskManager():
        mesh.Curve(2)
        return mesh


def uniform_field_phi_inc(mesh, H0):
    """phi_inc nodal vector for a uniform H0 z-direction field.

    For an external H = H0 z_hat, the scalar magnetic potential is
    phi = -H0 * z (so that H = -grad phi). The H1 nodal values are
    just the z-coordinate at each vertex multiplied by -H0.
    """
    nodes = np.array([mesh.vertices[i].point for i in range(mesh.nv)],
                     dtype=float)
    return (-H0 * nodes[:, 2]).astype(complex)


def run_case(R, sigma, mu_r, frequency, B0, half_thickness_guess,
             maxh_factor=4):
    """Run analytical / scalar / per-node and report."""
    omega = 2 * math.pi * frequency
    mu_eff = MU_0 * mu_r
    rho = 1.0 / sigma
    delta = math.sqrt(2 * rho / (omega * mu_eff))

    # Analytical
    Z_s_true = _dowell_Zs(R, rho, omega, mu_eff)
    ana = analytical_sphere_sibc(R, Z_s_true, B0, omega)

    # Build sphere mesh
    mesh = build_sphere_mesh(R, R / maxh_factor)
    nse = mesh.GetNE(BND)
    nv = mesh.nv
    H0 = B0 / MU_0

    # ScalarBIESIBCSolver
    solver = ScalarBIESIBCSolver(mesh, order=1)
    phi_inc = uniform_field_phi_inc(mesh, H0)

    # Scalar path: Dowell at user-supplied half_thickness_guess
    Z_s_scalar = _dowell_Zs(half_thickness_guess, rho, omega, mu_eff)
    sol_scalar = solver.solve(phi_inc, Z_s=Z_s_scalar, omega=omega)

    # Per-node path: extract per-panel R from the sphere mesh, then
    # build the H1-node Z_s array
    panels = _extract_panels_from_mesh(mesh, "sibc")
    R_panel = _compute_panel_local_radii(mesh, panels, default_R=10 * R)
    Zs_node, Zs_panel = _build_per_node_Zs(
        mesh, panels, R_panel, rho, omega, mu_eff, solver.ndof)
    sol_node = solver.solve(phi_inc, Z_s=Zs_node, omega=omega)

    # Per-node "RIGHT" baseline: per-node Z_s with R = R (uniform global)
    Zs_uniform = np.full(solver.ndof, Z_s_true, dtype=complex)
    sol_uniform = solver.solve(phi_inc, Z_s=Zs_uniform, omega=omega)

    print(f"--- R={R*1e3:.1f} mm  sigma={sigma:.2e}  mu_r={mu_r}  "
          f"f={frequency} Hz  B0={B0} T ---")
    print(f"  delta = {delta*1e3:.4f} mm  xi=R/delta={R/delta:.3f}")
    print(f"  Z_s(true)  = {Z_s_true.real:.4e} + {Z_s_true.imag:.4e}j")
    print(f"  half_thickness_guess (user input) = "
          f"{half_thickness_guess*1e3:.1f} mm")
    print(f"  Z_s(guess) = {Z_s_scalar.real:.4e} + {Z_s_scalar.imag:.4e}j")
    print()
    print(f"  R_local (mesh): mean={R_panel.mean()*1e3:.3f}  "
          f"min={R_panel.min()*1e3:.3f}  "
          f"max={R_panel.max()*1e3:.3f} mm  (target {R*1e3:.1f})")
    print()
    print(f"  {'method':<22} {'H_t_rms':>10} {'err':>10}")
    print(f"  {'-'*22} {'-'*10} {'-'*10}")
    print(f"  {'analytical':<22} {ana['H_t_rms']:>10.4f}      ---")
    err_s = (sol_scalar['H_t_rms'] - ana['H_t_rms']) / ana['H_t_rms'] * 100
    err_n = (sol_node['H_t_rms'] - ana['H_t_rms']) / ana['H_t_rms'] * 100
    err_u = (sol_uniform['H_t_rms'] - ana['H_t_rms']) / ana['H_t_rms'] * 100
    print(f"  {'BEM scalar(guess R)':<22} {sol_scalar['H_t_rms']:>10.4f} "
          f"{err_s:>+9.2f}%")
    print(f"  {'BEM per-node(mesh R)':<22} {sol_node['H_t_rms']:>10.4f} "
          f"{err_n:>+9.2f}%")
    print(f"  {'BEM per-node(true R)':<22} {sol_uniform['H_t_rms']:>10.4f} "
          f"{err_u:>+9.2f}%")
    print(f"  (mesh: nv={nv}, nse={nse}, ndof={solver.ndof})")
    print()

    return {
        "R": R, "frequency": frequency,
        "ana": ana['H_t_rms'],
        "scalar_guess": sol_scalar['H_t_rms'],
        "per_node_mesh": sol_node['H_t_rms'],
        "per_node_true": sol_uniform['H_t_rms'],
        "err_scalar_guess_pct": err_s,
        "err_per_node_mesh_pct": err_n,
        "err_per_node_true_pct": err_u,
    }


def main():
    print("=" * 70)
    print("Per-node Z_s SIBC validation: conducting sphere in uniform H field")
    print("=" * 70)
    print()

    # Three regimes:
    #   1. xi >> 1 (skin << R): Z_s asymptotic, R-independent.
    #      scalar(wrong R) and per-node should both match.
    #   2. xi ~ 1 (skin ~ R): Z_s sensitive to R.
    #      scalar(wrong R) WRONG, per-node CORRECT.
    #   3. xi << 1 (skin >> R, "transparent"): Z_s ~ jw mu R / 2.
    #      Strongly R-dependent.
    cases = [
        # (R, sigma, mu_r, frequency, B0, half_thickness_guess)
        # Cu, asymptotic
        (0.010, 5.8e7, 1, 1e6, 1.0, 0.005),
        # Cu, intermediate (skin depth ~ 2 mm at 1 kHz)
        (0.010, 5.8e7, 1, 1000, 1.0, 0.005),
        # Cu, transparent (skin depth ~ 21 mm at 10 Hz, R = 10 mm)
        (0.010, 5.8e7, 1, 10, 1.0, 0.005),
        # Steel mu_r=100, intermediate
        (0.010, 2e6, 100, 1000, 1.0, 0.005),
    ]

    results = []
    for c in cases:
        results.append(run_case(*c))

    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  {'case':<48} {'scalar':>10} {'per-node':>10}")
    for r in results:
        label = (f"R={r['R']*1e3:.0f}mm f={r['frequency']:.0f}Hz")
        print(f"  {label:<48} {r['err_scalar_guess_pct']:>+9.2f}% "
              f"{r['err_per_node_mesh_pct']:>+9.2f}%")

    print()
    print("Pass criteria:")
    print("  - per-node always within 5% of analytical")
    print("  - per-node beats scalar in the xi~1 and xi<<1 regimes")
    print("  - scalar = per-node (regression) in the xi>>1 regime")


if __name__ == "__main__":
    main()
