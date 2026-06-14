"""Kameari + Kelvin v22: Axisymmetric Cu disk Kameari-CLN proof of concept (2026-05-05).

Sugahara guidance: 3D HCurl Kameari-CLN (v15-v21) struggled to reach
canonical answers due to J·n leakage at conductor surfaces. In axisymmetric
the J = J_φ ê_φ has only azimuthal component, so J·n = 0 is satisfied
by geometry — no leakage possible. v22 demonstrates Kameari-CLN works
cleanly in axisymmetric.

Setup:
- Cu disk (R=10mm, full thickness 2mm) at origin
- Inner air half-circle (radius a_kelvin = 50mm) on r >= 0 side of (r,z) plane
- Outer Kelvin half-circle (same radius) at z = z_offset (Sugahara two-sphere)
- Periodic identification on Kelvin arcs
- GND vertex at exterior center (image of infinity)

Discretization:
- u = r × A_φ in H1 (Periodic, full domain, dirichlet="axis|axis_ext|GND")
- J_φ in H1 with definedon="conductor" — DISCONTINUOUS at conductor surface
  (J = 0 outside, naturally enforced — no current leak into air)

Weak form (axisymmetric A-formulation):
- a(u, v) = ∫_2D (ν/r) ∇u · ∇v dr dz
- f(v) = ∫_2D J_φ × v dr dz (in conductor)

Inner products (3D = 2π × 2D-with-r-weight):
- R_inv = 2π × ∫_2D (J_φ)²/σ × r dr dz (in conductor)
- L_n via dot: τ = 2π × ∫_2D J_φ × u_pot dr dz (no r — A_φ=u/r cancels with r)
- L_n via B²: 2π × ∫_2D (ν/r) × |∇u|² dr dz (per stage)

The historical CLN stage update has been retired.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting v22 (axisymmetric Cu disk Kameari-CLN)...", flush=True)

from netgen.occ import (WorkPlane, Vertex, Pnt, Axes, X, Y, Z, MoveTo,
                         Glue, OCCGeometry, IdentificationType, Dir)
from ngsolve import (Mesh, H1, Periodic, BilinearForm, LinearForm,
                     GridFunction, CoefficientFunction, grad, dx,
                     x, y, IfPos, sqrt, Integrate, TaskManager, ngsglobals)
from kelvin_source import kelvin_nu_factor_axisym_cf, build_material_cf
from math import pi
import json, time
from pathlib import Path
import numpy as np

# === Parameters ===
R_disk = 10e-3       # disk radius [m]
half_thick = 1e-3    # disk half-thickness [m] (full 2*half_thick = 2mm)
sigma_Cu = 5.8e7
mu0 = 4 * pi * 1e-7
nu0 = 1.0 / mu0

a_kelvin = 50e-3     # Kelvin boundary radius
z_offset = 5 * a_kelvin
maxh = 5e-3
maxh_disk = 0.5e-3
ORDER = 2
N_STAGES = 6
B0 = 1.0

# Analytical reference (infinite cylinder eddy mode):
# τ_n = μ_0 σ R² / j_{1,n}²,  j_{1,1} ≈ 3.83171
TAU_INF_CYL_US = (mu0 * sigma_Cu * R_disk**2) / 3.83171**2 * 1e6


def build_geo():
    """Two-half-circle axisymmetric geometry per sphere_with_Kelvin pattern."""
    # Inner half-circle
    wp1 = WorkPlane()
    inner_full = wp1.Circle(a_kelvin).Face()
    inner_full.name = "air"

    # Disk (rectangle in r-z plane, r in [0, R_disk], z in [-half_thick, half_thick])
    wp_disk = WorkPlane(Axes((0, -half_thick, 0), n=Z, h=X))
    disk_full = wp_disk.Rectangle(R_disk, 2 * half_thick).Face()
    disk_full.name = "conductor"
    disk_full.maxh = maxh_disk

    # Outer Kelvin half-circle
    wp3 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
    outer_full = wp3.Circle(a_kelvin).Face()
    outer_full.name = "kelvin"

    # GND vertex
    gnd_point = Vertex(Pnt(0, z_offset, 0))
    gnd_point.name = "GND"

    # Cut to half (keep r >= 0)
    cutter = MoveTo(-a_kelvin - 0.1, -a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    cutter_ext = MoveTo(-a_kelvin - 0.1, z_offset - a_kelvin - 0.1).Rectangle(
        a_kelvin + 0.1, 2 * a_kelvin + 0.2).Face()
    inner_half = inner_full - cutter
    disk_half = disk_full - cutter
    outer_half = outer_full - cutter_ext

    # Inner air = inner_half - disk_half
    inner_air = inner_half - disk_half
    inner_air.name = "air"

    # Edge naming (axis, kelvin_int, kelvin_ext, GND adjacent edges)
    kelvin_inner_edges = []
    for edge in inner_air.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + v0.p.y ** 2)
            d1 = sqrt(v1.p.x ** 2 + v1.p.y ** 2)
            is_kelvin_arc = (abs(d0 - a_kelvin) < 0.01 * a_kelvin
                              and abs(d1 - a_kelvin) < 0.01 * a_kelvin
                              and cx > 0.001 * a_kelvin)
        except Exception:
            is_kelvin_arc = False

        if cx < 0.001 * a_kelvin:
            edge.name = "axis"
        elif is_kelvin_arc:
            edge.name = "kelvin_int"
            kelvin_inner_edges.append(edge)
        else:
            edge.name = "interface"

    for edge in disk_half.edges:
        cx = edge.center.x
        if cx < 0.001 * a_kelvin:
            edge.name = "axis"
        else:
            edge.name = "disk_bnd"

    kelvin_outer_edges = []
    for edge in outer_half.edges:
        cx = edge.center.x
        try:
            v0, v1 = edge.vertices
            d0 = sqrt(v0.p.x ** 2 + (v0.p.y - z_offset) ** 2)
            d1 = sqrt(v1.p.x ** 2 + (v1.p.y - z_offset) ** 2)
            is_kelvin_arc = (abs(d0 - a_kelvin) < 0.01 * a_kelvin
                              and abs(d1 - a_kelvin) < 0.01 * a_kelvin
                              and cx > 0.001 * a_kelvin)
        except Exception:
            is_kelvin_arc = False

        if cx < 0.001 * a_kelvin:
            edge.name = "axis_ext"
        elif is_kelvin_arc:
            edge.name = "kelvin_ext"
            kelvin_outer_edges.append(edge)
        else:
            edge.name = "default"

    # Periodic identification of Kelvin arcs (match by y sign)
    matched = 0
    for ie in kelvin_inner_edges:
        iy = ie.center.y
        for oe in kelvin_outer_edges:
            oy = oe.center.y - z_offset
            if (iy > 0 and oy > 0) or (iy < 0 and oy < 0):
                ie.Identify(oe, "kelvin", IdentificationType.PERIODIC)
                matched += 1
                break
    print(f"  Periodic edge pairs matched: {matched}", flush=True)

    shape = Glue([inner_air, disk_half, outer_half, gnd_point])
    return OCCGeometry(shape, dim=2)


def main():
    ngsglobals.msg_level = 0
    print(f"=== v22: axisymmetric Cu disk Kameari-CLN ===", flush=True)
    print(f"  Disk: R={R_disk*1000:.1f}mm, thickness={2*half_thick*1000:.1f}mm",
          flush=True)
    print(f"  Kelvin: a={a_kelvin*1000:.0f}mm, z_offset={z_offset*1000:.0f}mm",
          flush=True)
    print(f"  Analytical (∞-cylinder leading) τ_1 = {TAU_INF_CYL_US:.2f} μs",
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

    # === Material with Kelvin pullback ===
    r_weight = IfPos(x - 1e-12, x, 1e-12)
    nu_kelvin_factor = kelvin_nu_factor_axisym_cf(z_offset=z_offset, R=a_kelvin)
    nu_cf = build_material_cf(mesh, nu0, nu_kelvin_factor)
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu}, default=0.0)
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0 / sigma_Cu}, default=0.0)

    # === FE spaces ===
    # u = r * A_phi (continuous, full domain)
    fes_u = Periodic(H1(mesh, order=ORDER, dirichlet="axis|axis_ext|GND"))
    # J_phi (discontinuous at conductor surface: P2 on conductor, 0 outside)
    fes_J = H1(mesh, order=ORDER, definedon=mesh.Materials("conductor"),
                dirichlet="axis")
    print(f"  fes_u (Periodic H1) ndof = {fes_u.ndof}", flush=True)
    print(f"  fes_J (H1 on conductor) ndof = {fes_J.ndof}", flush=True)

    u, v = fes_u.TnT()

    # === Bilinear form (axisymmetric A-formulation) ===
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

    # === Initial J_φ = σ × A_φ_imposed = σ × B_0 × r/2 (in conductor) ===
    A_phi_imposed_cf = B0 * x / 2  # = B_0 × r / 2
    gfJ = GridFunction(fes_J, name="J_imp")
    gfJ.Set(sigma_Cu * A_phi_imposed_cf,
             definedon=mesh.Materials("conductor"))

    # Diagnostic: how well does the H1 P2 projection match analytical?
    proj_err_int = float(Integrate(
        (gfJ - sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    proj_norm_int = float(Integrate(
        (sigma_Cu * A_phi_imposed_cf) ** 2 * dx("conductor"), mesh))
    rel_proj_err = (proj_err_int / proj_norm_int) ** 0.5 if proj_norm_int else 0
    print(f"  ||gfJ - σA_phys||/||σA_phys|| in conductor = "
          f"{rel_proj_err:.4e}", flush=True)

    # R_0 verification (Sugahara CLN convention with 3D = 2π × 2D-with-r):
    twopi = 2 * pi
    R_inv_check = twopi * float(Integrate(
        gfJ * gfJ * sigma_inv_cf * r_weight * dx("conductor"), mesh))
    R0_test = 1.0 / R_inv_check if R_inv_check > 1e-30 else float('inf')
    # Analytical R_0 for cylinder (sigma A_phys with A_phys = B_0 r/2):
    # R_inv_anal = 2π × σ × (B_0/2)² × ∫ r³ × 2 half_thick dr
    #            = 2π × σ × B_0²/4 × (R^4/4) × (2 half_thick)
    #            = π × σ × B_0² × half_thick × R^4 / 4
    R_inv_anal = pi * sigma_Cu * B0**2 * half_thick * R_disk**4 / 4
    R0_anal = 1.0 / R_inv_anal
    print(f"  R_0 verification: {R0_test:.6e} vs analytical {R0_anal:.6e}, "
          f"ratio {R0_test/R0_anal:.6f}\n", flush=True)

    raise RuntimeError("The historical CLN field-reduction implementation was retired.")


if __name__ == "__main__":
    main()
