"""Kameari-CLN v22e: axisymmetric Cu disk WITHOUT Kelvin (simple outer Dirichlet).

Diagnostic: if v22b's τ_dot vs τ_B2 mismatch (15-43% from n=1) is caused by
Kelvin pullback, then this version (no Kelvin) should give <1% match.
If still mismatched, Kelvin is innocent and J leakage / iteration is suspect.

Setup identical to v22b except:
- No outer Kelvin half-circle, no GND vertex
- Single half-circle of radius A_outer = 5*R_disk = 50mm
- Outer arc has Dirichlet u=0 (truncation BC)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v22e (no Kelvin, outer Dirichlet)...", flush=True)

from netgen.occ import (WorkPlane, Pnt, Axes, X, Z, MoveTo,
                         Glue, OCCGeometry)
from ngsolve import (Mesh, H1, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, grad, dx,
                     x, y, IfPos, sqrt, Integrate, TaskManager, ngsglobals)
from math import pi
import json, time
from pathlib import Path

# === Parameters (matching v22b for direct comparison) ===
R_disk = 10e-3
half_thick = 1e-3
sigma_Cu = 5.8e7
mu0 = 4 * pi * 1e-7
nu0 = 1.0 / mu0

A_outer = 50e-3  # outer truncation radius (= a_kelvin in v22b)
maxh = 2e-3
maxh_disk = 0.1e-3
ORDER = 3
N_STAGES = 6
B0 = 1.0

TAU_INF_CYL_US = (mu0 * sigma_Cu * R_disk**2) / 3.83171**2 * 1e6


def build_geo():
    """Single half-circle: disk + air, with outer Dirichlet truncation."""
    # Outer half-circle (replaces the Kelvin two-sphere setup)
    wp_outer = WorkPlane()
    outer_full = wp_outer.Circle(A_outer).Face()
    outer_full.name = "air"

    # Disk
    wp_disk = WorkPlane(Axes((0, -half_thick, 0), n=Z, h=X))
    disk_full = wp_disk.Rectangle(R_disk, 2 * half_thick).Face()
    disk_full.name = "conductor"
    disk_full.maxh = maxh_disk

    # Cut to half-circle (r >= 0)
    cutter = MoveTo(-A_outer - 0.1, -A_outer - 0.1).Rectangle(
        A_outer + 0.1, 2 * A_outer + 0.2).Face()
    outer_half = outer_full - cutter
    disk_half = disk_full - cutter

    # Air region
    air_only = outer_half - disk_half
    air_only.name = "air"

    # Edge naming
    for edge in air_only.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + v0.p.y ** 2)
            d1 = sqrt(v1.p.x ** 2 + v1.p.y ** 2)
            is_outer_arc = (abs(d0 - A_outer) < 0.01 * A_outer
                             and abs(d1 - A_outer) < 0.01 * A_outer
                             and cx > 0.001 * A_outer)
        except Exception:
            is_outer_arc = False

        if cx < 0.001 * A_outer:
            edge.name = "axis"
        elif is_outer_arc:
            edge.name = "outer_arc"
        else:
            edge.name = "interface"

    for edge in disk_half.edges:
        cx = edge.center.x
        if cx < 0.001 * A_outer:
            edge.name = "axis"
        else:
            edge.name = "disk_bnd"

    shape = Glue([air_only, disk_half])
    return OCCGeometry(shape, dim=2)


def main():
    ngsglobals.msg_level = 0
    print(f"=== v22e: axisymmetric Cu disk WITHOUT Kelvin ===", flush=True)
    print(f"  Disk: R={R_disk*1000:.1f}mm, t={2*half_thick*1000:.1f}mm",
          flush=True)
    print(f"  Outer truncation radius: {A_outer*1000:.0f}mm",
          flush=True)
    print(f"  N_STAGES = {N_STAGES}\n", flush=True)

    geo = build_geo()
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh, grading=0.5)
    mesh = Mesh(ngmesh)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}", flush=True)
    print(f"  boundaries = {mesh.GetBoundaries()}\n", flush=True)

    r_weight = IfPos(x - 1e-12, x, 1e-12)
    nu_cf = CoefficientFunction(nu0)  # uniform vacuum (no Kelvin)
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu}, default=0.0)
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0 / sigma_Cu}, default=0.0)

    # FE spaces — no Periodic (no Kelvin); Dirichlet on axis + outer_arc
    fes_u = H1(mesh, order=ORDER, dirichlet="axis|outer_arc")
    fes_J = H1(mesh, order=ORDER, definedon=mesh.Materials("conductor"),
                dirichlet="axis")
    print(f"  fes_u (H1, no Periodic) ndof = {fes_u.ndof}", flush=True)
    print(f"  fes_J (H1 on conductor) ndof = {fes_J.ndof}", flush=True)

    u, v = fes_u.TnT()

    a_form = BilinearForm(fes_u)
    a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_form.Assemble()
        try:
            inv_u = a_form.mat.Inverse(fes_u.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_u = a_form.mat.Inverse(fes_u.FreeDofs(),
                                        inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    A_phi_imposed_cf = B0 * x / 2
    gfJ = GridFunction(fes_J)
    gfJ.Set(sigma_Cu * A_phi_imposed_cf,
             definedon=mesh.Materials("conductor"))

    proj_err = float(Integrate(
        (gfJ - sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    proj_norm = float(Integrate(
        (sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    print(f"  ||gfJ - σA||/||σA|| = {(proj_err/proj_norm)**0.5:.4e}",
          flush=True)

    twopi = 2 * pi
    R_inv_check = twopi * float(Integrate(
        gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
    R0_test = 1.0 / R_inv_check
    R_inv_anal = pi * sigma_Cu * B0**2 * half_thick * R_disk**4 / 4
    print(f"  R_0: {R0_test:.6e} vs analytical {1/R_inv_anal:.6e}, "
          f"ratio {R0_test * R_inv_anal:.6f}\n", flush=True)

    raise RuntimeError("Historical CLN field-stage reduction is retired.")


if __name__ == "__main__":
    main()
