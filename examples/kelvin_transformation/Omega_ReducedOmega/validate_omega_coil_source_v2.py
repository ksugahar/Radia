"""validate_omega_coil_source_v2.py — EMPY total+reduced split for
Omega + coil source + Kelvin (Plan A v2, no compromise).

This implementation mirrors the EMPY Omega_ReducedOmega.py method
exactly (Static/Omega_ReducedOmega.py and the working uniform-field
reference 3D_cylinder_with_Kelvin.py). The fully-reduced
ScalarPotentialSolver.solve_single_potential is NOT used because it
implicitly assumes H_s defined globally with proper Kelvin pullback,
which we cannot satisfy from a physical-bbox voxel.

Method (EMPY):
    omega = Omega_t in iron, Omega_r in air_inner+air_outer
    H_s = grad(Omega_s)  (Hs = coil Biot-Savart, Omega_s = scalar
        potential of coil computed by line-integral)
    Bs = mu0 * Hs in air, 0 in iron and Kelvin

    a = Mu * grad(omega) * grad(psi) * dx(all)         + Kelvin metric
    gfOmega.Set(Omega_s, BND, iron-air boundary)        (Dirichlet lift)
    f = Mu * grad(gfOmega) * grad(psi) * dx(air_inner)  (lift contrib)
        + (n . Bs) * psi * ds(iron-air boundary)        (Neumann jump)

    Bt(iron) = Mu * grad(omega_in_iron)
    Br(air)  = mu0 * (grad(omega_in_air) - grad(Omega_s_lifted))
    BField(air)   = Br + Bs
    BField(iron)  = Bt
    BField(Kelvin) is purely the Kelvin pullback of air contribution

Acceptance: < 5 % max |B| diff on sample points.

BUGS HUNTED AND FIXED during this implementation:
  1. `rad.Fld(coil, 'phi')` is NOT the magnetic scalar potential for
     arc coils (rad_arc_current.cpp uses an incorrect 1/r^3 integrand).
     Fix: compute Omega_s by line-integrating rad.Fld('h') from a
     below-iron reference point along a path that doesn't cross the
     coil cut.
  2. Mesh inversion from overlapping Kelvin spheres. inner@origin,
     R=0.18 and outer@(0.30,0,0), R=0.18 overlap (2R > offset). OCC
     Glue assigned iron tets to the air_outer material, giving wrong
     Mu inside iron. Fix: kelvin_center=(0.50,0,0) → offset 2.78×R.
  3. Bs voxel bbox too small — clamped values outside iron bbox gave
     constant `|B|_FEM` for all probe points in air_inner beyond
     z>0.05. Fix: BS_BBOX = full physical-region bbox ±phys_R.
  4. strict Dirichlet on iron_surf freezes omega to harmonic-extension
     of Omega_s in iron → grad(omega) = Hs → Bt = mu_iron*Hs (the
     long-cylinder limit, no demagnetization). Fix: use
     `dirichlet="GND"` only; iron_surf BC enters WEAKLY via the air
     RHS lift integral. This matches EMPY/3D_cylinder_with_Kelvin.py
     and allows demagnetization to develop.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

import radia as rad

from ngsolve import (Mesh, TaskManager, CF, GridFunction, BilinearForm,
                     LinearForm, H1, Periodic, grad, dx, ds, x as ngx,
                     y as ngy, z as ngz, sqrt as ngsqrt, specialcf,
                     VoxelCoefficient, Integrate, BND, VOL)
from netgen.occ import (Cylinder, Sphere, Pnt, Z, Vertex, Glue, OCCGeometry,
                         IdentificationType)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src', 'radia'))
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)
from netgen_mesh_import import netgen_mesh_to_radia
from kelvin_source import kelvin_mu_factor_3d_cf

MU_0 = 4 * math.pi * 1e-7

# --- Geometry / physics ---------------------------------------------
cylinder_radius = 0.03
cylinder_height = 0.06
mu_r = 100

R_coil = 0.08
a_coil = 0.008
gap_deg = 5.0
arc_deg = 360.0 - gap_deg
z_coil = 0.10
I_total = 1000.0

phys_R = 0.18
# kelvin_center offset MUST be > 2 * phys_R to avoid sphere overlap
# (working ref 3D_cylinder_with_Kelvin.py uses offset_z=4.0 with R=1.5 → 2.67×).
# Earlier offset 0.30 (1.67×) caused mesh inversion: iron at origin was
# being meshed AS air_outer because of OCC Glue overlap behaviour.
kelvin_center = (0.50, 0.0, 0.0)   # offset/R = 2.78 — safely non-overlapping

maxh_iron = 0.008
maxh_air = 0.025
maxh_kelvin = 0.03
fe_order = 2
radia_maxh = 0.004   # tighter than FEM mesh for convergence check

# Reference point for Omega_s line integral: well below the cylinder,
# along axis. The straight-line path from this reference to any
# evaluation point in the iron region (z in [-0.03, 0.03], r<0.03) does
# NOT cross the coil at z=0.10 (coil cut surface).
OMEGA_REF = np.array([0.0, 0.0, -1.0])

# Voxel grid for Omega_s — covers the iron + small margin.
# Omega_s is only needed on the iron-air boundary (Dirichlet lift).
OS_BBOX = [[-0.05, 0.05], [-0.05, 0.05], [-0.05, 0.05]]
OS_RES = 41     # 41^3 = 68921 voxel pts (~2.4mm spacing)
OS_LINE_N = 100  # samples per line integral

# Voxel grid for B_s — needs to cover all of air_inner (r < phys_R = 0.18)
# plus the iron surface integration. Bs is NOT needed in the Kelvin region
# (masked by per-material CF in build_BField).
BS_BBOX = [[-0.18, 0.18], [-0.18, 0.18], [-0.18, 0.18]]
BS_RES = 51


# ---------------------------------------------------------------
# Build coil
# ---------------------------------------------------------------
def build_radia_coil():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    coil = rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, 200, 'man', 'z', J0)
    return coil


# ---------------------------------------------------------------
# Reference: Radia BEM with iron
# ---------------------------------------------------------------
def build_and_solve_radia_reference():
    rad.UtiDelAll()
    J0 = I_total / (2 * a_coil) ** 2
    coil = rad.ObjArcCur(
        [0, 0, z_coil],
        [R_coil - a_coil, R_coil + a_coil],
        [0, math.radians(arc_deg)],
        2 * a_coil, 200, 'man', 'z', J0)
    half_h = cylinder_height / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("magnetic")
    mag_cyl.maxh = radia_maxh
    geo = OCCGeometry(mag_cyl)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=radia_maxh)
    mesh_rad = Mesh(ngmesh)
    cyl_obj = netgen_mesh_to_radia(mesh_rad,
                                    material={'magnetization': [0, 0, 0]},
                                    units='m', material_filter='magnetic',
                                    verbose=False)
    rad.MatApl(cyl_obj, rad.MatLin(mu_r))
    grp = rad.ObjCnt([coil, cyl_obj])
    rad.SolverConfig(bicgstab_tol=1e-8, relax_param=0.0)
    print(f"  Radia tetra: {mesh_rad.ne} elements")
    t0 = time.perf_counter()
    # Tighter tol + more iterations for convergence reference
    result = rad.Solve(grp, 1e-7, 2000, 0)
    print(f"  Radia Solve: {time.perf_counter()-t0:.1f}s, "
          f"residual={result[0]:.3e}")
    return coil, cyl_obj, grp


# ---------------------------------------------------------------
# Mesh build (EMPY layout: total + reduced + Kelvin, periodic BC)
# ---------------------------------------------------------------
def build_kelvin_mesh():
    half_h = cylinder_height / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z,
                        r=cylinder_radius, h=cylinder_height)
    mag_cyl.mat("iron")
    mag_cyl.maxh = maxh_iron
    for f in mag_cyl.faces:
        f.name = "iron_surf"
    inner = Sphere(Pnt(0, 0, 0), phys_R)
    inner.maxh = maxh_air
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - mag_cyl
    inner_air.mat("air_inner")
    outer = Sphere(Pnt(*kelvin_center), phys_R)
    outer.maxh = maxh_kelvin
    outer.mat("air_outer")
    for f in outer.faces:
        f.name = "kelvin_ext"
    gnd = Vertex(Pnt(*kelvin_center))
    gnd.name = "GND"
    geo = Glue([inner_air, mag_cyl, outer, gnd])
    geo.solids[0].name = "air_inner"
    geo.solids[1].name = "iron"
    geo.solids[2].name = "air_outer"
    k_int = k_ext = None
    for solid in geo.solids:
        for f in solid.faces:
            if f.name == "kelvin_int" and k_int is None:
                k_int = f
            elif f.name == "kelvin_ext" and k_ext is None:
                k_ext = f
    if k_int is None or k_ext is None:
        raise RuntimeError("Could not find Kelvin periodic faces")
    k_int.Identify(k_ext, "periodic", IdentificationType.PERIODIC)
    with TaskManager():
        ngmesh = OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(fe_order)
    return mesh


# ---------------------------------------------------------------
# Build Omega_s as a voxel CF via line integral of H_s
# ---------------------------------------------------------------
def build_Omega_s_voxel(coil, bbox, resolution, n_samples):
    """Compute Omega_s(p) = integral_ref^p H_s . dl on a voxel grid.

    H_s is the coil's vacuum Biot-Savart H field. The straight path from
    OMEGA_REF to each grid point must not cross the coil's cut surface
    (here the disk at z=z_coil). Reference is at (0,0,-1.0); all grid
    points have z in [-0.05, 0.05] < z_coil = 0.1 -- safe.

    Uses grad(Omega_s) = H_s convention (matches EMPY).
    """
    nx = ny = nz = resolution
    xs = np.linspace(bbox[0][0], bbox[0][1], nx)
    ys = np.linspace(bbox[1][0], bbox[1][1], ny)
    zs = np.linspace(bbox[2][0], bbox[2][1], nz)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
    grid_pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    n_pts = len(grid_pts)
    print(f"  Omega_s voxel grid: {nx}x{ny}x{nz} = {n_pts} pts")
    print(f"  line samples per pt: {n_samples}")

    # Trapezoid rule along straight line ref -> p
    # path(t) = ref + t * (p - ref),  t in [0, 1]
    # dl = (p - ref) dt
    # Omega_s(p) = sum_i H(path(t_i)) . (p-ref) * w_i  (trapezoid weights)
    ts = np.linspace(0.0, 1.0, n_samples)   # shape (M,)
    # Trapezoid weights
    dt = ts[1] - ts[0]
    w = np.full(n_samples, dt)
    w[0] *= 0.5
    w[-1] *= 0.5

    # Build all sample points: shape (n_pts, M, 3)
    delta = grid_pts - OMEGA_REF[None, :]      # (n_pts, 3)
    sample_pts = (OMEGA_REF[None, None, :] +
                   ts[None, :, None] * delta[:, None, :])  # (n_pts, M, 3)
    sample_flat = sample_pts.reshape(-1, 3)             # (n_pts*M, 3)

    print(f"  evaluating H on {len(sample_flat)} batch points...")
    t0 = time.perf_counter()
    H_flat = np.asarray(rad.Fld(coil, 'h', sample_flat))   # (n_pts*M, 3)
    print(f"    rad.Fld batch H: {time.perf_counter()-t0:.1f}s")
    H = H_flat.reshape(n_pts, n_samples, 3)

    # H . delta  (delta repeated along axis 1)
    integrand = np.einsum('pmi,pi->pm', H, delta)         # (n_pts, M)
    Omega_s_arr = np.einsum('pm,m->p', integrand, w)      # (n_pts,)
    Omega_s_grid = Omega_s_arr.reshape(nx, ny, nz)

    # Build VoxelCoefficient (note: NGSolve VoxelCoefficient expects
    # data layout with last axis being x for whatever reason — see
    # set_source_from_radia in scalar_potential_solver.py for the
    # transpose convention)
    data = np.ascontiguousarray(Omega_s_grid.transpose(2, 1, 0))
    start = (bbox[0][0], bbox[1][0], bbox[2][0])
    end = (bbox[0][1], bbox[1][1], bbox[2][1])
    return VoxelCoefficient(start, end, data, linear=True)


# ---------------------------------------------------------------
# Build Bs as voxel CF (vector)
# ---------------------------------------------------------------
def build_Bs_voxel(coil, bbox, resolution):
    nx = ny = nz = resolution
    xs = np.linspace(bbox[0][0], bbox[0][1], nx)
    ys = np.linspace(bbox[1][0], bbox[1][1], ny)
    zs = np.linspace(bbox[2][0], bbox[2][1], nz)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
    pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    print(f"  Bs voxel: {nx}x{ny}x{nz} = {len(pts)} pts")
    t0 = time.perf_counter()
    B_flat = np.asarray(rad.Fld(coil, 'b', pts))   # (N, 3) in Tesla
    print(f"    rad.Fld batch B: {time.perf_counter()-t0:.1f}s")
    start = (bbox[0][0], bbox[1][0], bbox[2][0])
    end = (bbox[0][1], bbox[1][1], bbox[2][1])
    cfs = []
    for c in range(3):
        d = B_flat[:, c].reshape(nx, ny, nz)
        d = np.ascontiguousarray(d.transpose(2, 1, 0))
        cfs.append(VoxelCoefficient(start, end, d, linear=True))
    return CF(tuple(cfs))


# ---------------------------------------------------------------
# Sanity check: numerical grad(Omega_s) vs Hs at probe points
# ---------------------------------------------------------------
def check_Omega_s_consistency(coil, Omega_s_cf, mesh, eps=1e-4):
    print("\n[Sanity] grad(Omega_s) finite-diff vs rad.Fld(coil,'h')")
    for pt in [(0.0, 0.0, 0.0), (0.02, 0.0, 0.0), (0.0, 0.0, 0.02)]:
        try:
            Op = Omega_s_cf(mesh(pt[0] + eps, pt[1], pt[2]))
            Om = Omega_s_cf(mesh(pt[0] - eps, pt[1], pt[2]))
            dOdx = (Op - Om) / (2 * eps)
            Op = Omega_s_cf(mesh(pt[0], pt[1] + eps, pt[2]))
            Om = Omega_s_cf(mesh(pt[0], pt[1] - eps, pt[2]))
            dOdy = (Op - Om) / (2 * eps)
            Op = Omega_s_cf(mesh(pt[0], pt[1], pt[2] + eps))
            Om = Omega_s_cf(mesh(pt[0], pt[1], pt[2] - eps))
            dOdz = (Op - Om) / (2 * eps)
            H_rad = rad.Fld(coil, 'h', list(pt))
            print(f"  pt={pt}  grad(Omega_s)=({dOdx:.3e},{dOdy:.3e},{dOdz:.3e})  "
                  f"H_rad={tuple(f'{v:.3e}' for v in H_rad)}")
        except Exception as e:
            print(f"  pt={pt}  ERR {e}")


# ---------------------------------------------------------------
# EMPY-style total+reduced solve
# ---------------------------------------------------------------
def solve_empy_pattern(mesh, Omega_s_cf, Bs_cf):
    print("\n[Solve] EMPY total+reduced pattern with Kelvin")
    # Material per region: iron has mu_r * mu0; air_inner has mu0;
    # air_outer (Kelvin) has mu0 * (R/r')^2 (Nagamine canonical).
    mu_kelvin_factor = kelvin_mu_factor_3d_cf(center=kelvin_center, R=phys_R)
    mat_dict = {
        "iron": mu_r * MU_0,
        "air_inner": MU_0,
        "air_outer": MU_0 * mu_kelvin_factor,
    }
    Mu = CF([mat_dict[m] for m in mesh.GetMaterials()])

    # H1 + Periodic. CRITICAL: iron_surf must NOT be Dirichlet — strict
    # Dirichlet pins omega_iron = harmonic_extension(Omega_s) which gives
    # the long-cylinder limit (no demagnetization). The EMPY/working-ref
    # pattern uses dirichlet="GND" only; iron_surf BC enters WEAKLY via the
    # RHS lift integral, allowing demagnetization to develop in iron.
    fes_pre = H1(mesh, order=fe_order, dirichlet="GND")
    fes = Periodic(fes_pre)
    print(f"  ndof = {fes.ndof}  free = {sum(1 for d in fes.FreeDofs() if d)}")
    omega = fes.TrialFunction()
    psi = fes.TestFunction()

    # Bilinear form (Mu * grad . grad over all regions)
    a = BilinearForm(fes)
    a += Mu * grad(omega) * grad(psi) * dx
    with TaskManager():
        a.Assemble()

    # Lift: gfOmega = Omega_s on iron-air boundary, 0 elsewhere
    gfOmega = GridFunction(fes)
    gfOmega.Set(Omega_s_cf, BND, mesh.Boundaries("iron_surf"))
    print(f"  lift max |gfOmega| = {max(abs(np.asarray(gfOmega.vec.FV()))):.3e}")

    # RHS: Mu * grad(gfOmega) * grad(psi) * dx(air_inner)
    #    + (n . Bs) * psi * ds(iron_surf)
    f = LinearForm(fes)
    f += Mu * grad(gfOmega) * grad(psi) * dx("air_inner")
    # Normal points OUTWARD from physical (cf 3D_cylinder_with_Kelvin sign note)
    normal = -specialcf.normal(mesh.dim)
    f += (normal * Bs_cf) * psi * ds("iron_surf")
    with TaskManager():
        f.Assemble()

    # Direct solve (working-ref pattern). iron_surf is NOT Dirichlet, so
    # all relevant DoFs are FreeDofs and no inhomogeneous correction needed.
    print("  inverting (sparsecholesky)...")
    t0 = time.perf_counter()
    sol = GridFunction(fes)
    sol.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                   inverse="sparsecholesky") * f.vec
    print(f"    inverse: {time.perf_counter()-t0:.1f}s")
    return Mu, sol, gfOmega, mu_kelvin_factor


# ---------------------------------------------------------------
# Build B-field for evaluation
# ---------------------------------------------------------------
def build_BField(mesh, Mu, sol, Omega_s_cf, mu_kelvin_factor, Bs_cf):
    # Simplified post-processing: use `sol` directly (no Ot/Orr extraction
    # — those L2 projections produced near-zero Ot on iron for unclear
    # NGSolve reasons). Compute B per material:
    #   iron:       Bt  = mu_iron * grad(sol)           (total)
    #   air_inner:  Br+Bs = mu0 * (grad(sol) - grad(Oxr)) + Bs
    #               = mu0 * grad(sol) - mu0 * grad(Oxr) + mu0 * Hs
    #               Oxr is the BND lift of Omega_s; in air_inner the lift
    #               is a thin boundary skin, grad(Oxr) is non-zero only in
    #               elements adjacent to iron surface.
    #   air_outer:  mu_kelvin * (grad(sol) - grad(Oxr))  (no Bs in Kelvin)
    fes_air = H1(mesh, order=fe_order, definedon="air_inner|air_outer")
    Oxr = GridFunction(fes_air)
    Oxr.Set(Omega_s_cf, BND, mesh.Boundaries("iron_surf"))

    zero_vec = CF((0.0, 0.0, 0.0))
    # Per-material B assembly
    B_per_mat = []
    for m in mesh.GetMaterials():
        if m == "iron":
            B_per_mat.append(Mu * grad(sol))
        elif m == "air_inner":
            B_per_mat.append(Mu * (grad(sol) - grad(Oxr)) + Bs_cf)
        elif m == "air_outer":
            B_per_mat.append(Mu * (grad(sol) - grad(Oxr)))
        else:
            B_per_mat.append(zero_vec)
    BField = CF(B_per_mat, dims=(3,))
    return BField


def compare_B(mesh, BField, grp_radia, sample_points, label):
    print(f"\n{label}")
    print(f"  {'point':>22s}  {'|B|_FEM':>12s}  {'|B|_Radia':>12s}  "
          f"{'rel':>9s}")
    rels = []
    for pt in sample_points:
        try:
            mip = mesh(pt[0], pt[1], pt[2])
            B_fem = np.array([BField(mip)[i] for i in range(3)])
        except Exception:
            B_fem = np.full(3, np.nan)
        B_rad = np.array(rad.Fld(grp_radia, 'b', list(pt)))
        m_fem = float(np.linalg.norm(B_fem))
        m_rad = float(np.linalg.norm(B_rad))
        if np.isnan(m_fem):
            continue
        rel = (m_fem - m_rad) / max(m_rad, 1e-12)
        rels.append(rel)
        ps = f"({pt[0]:+.3f},{pt[1]:+.3f},{pt[2]:+.3f})"
        print(f"  {ps:>22s}  {m_fem:11.4e}  {m_rad:11.4e}  {rel*100:+7.2f}%")
    return rels


def main():
    print("=" * 72)
    print("validate_omega_coil_source v2 — EMPY total+reduced + Kelvin")
    print("=" * 72)

    print("\n[1] Build Radia coil (for source field sampling)")
    coil = build_radia_coil()

    print("\n[2] Build Omega_s voxel CF via line-integral of H_s")
    # Need a (temporary) mesh just to evaluate the CF — use the one
    # we'll solve on later, save by building up front
    print("\n[3] Build mesh (iron, air_inner, air_outer with periodic Kelvin BC)")
    mesh = build_kelvin_mesh()
    print(f"  ne={mesh.ne}  nv={mesh.nv}  mats={mesh.GetMaterials()}")

    Omega_s_cf = build_Omega_s_voxel(coil, OS_BBOX, OS_RES, OS_LINE_N)
    Bs_cf = build_Bs_voxel(coil, BS_BBOX, BS_RES)

    check_Omega_s_consistency(coil, Omega_s_cf, mesh)

    print("\n[4] Build Radia BEM reference (coil + iron, Solve)")
    coil_ref, cyl_ref, grp_ref = build_and_solve_radia_reference()

    # Re-build coil for the FEM source extraction since UtiDelAll() in
    # the reference build wiped it. We use the ref's coil handle for
    # everything afterwards. But Omega_s_cf and Bs_cf are now numpy
    # voxels — they're not tied to Radia handles. Good.

    Mu, sol, gfOmega, mu_kelvin_factor = solve_empy_pattern(
        mesh, Omega_s_cf, Bs_cf)

    # Post-solve diagnostics: probe sol at iron interior
    print("\n[Post-solve diagnostic] sol & Omega_s at iron points")
    for pt in [(0.0, 0.0, 0.0), (0.01, 0.0, 0.0), (0.0, 0.0, 0.02),
                (0.03, 0.0, 0.0)]:
        try:
            mip = mesh(*pt)
            print(f"  pt={pt}  sol={float(sol(mip)):+.4e}  "
                  f"Omega_s={float(Omega_s_cf(mip)):+.4e}")
        except Exception as e:
            print(f"  pt={pt}  ERR {e}")

    BField = build_BField(mesh, Mu, sol, Omega_s_cf, mu_kelvin_factor, Bs_cf)

    # Raw grad(sol) + Mu*grad(sol) probe, incl material
    print("\n[Raw grad(sol) probe inside iron]")
    gs = grad(sol)
    Mu_gs = Mu * grad(sol)
    mats = mesh.GetMaterials()
    for pt in [(0.0, 0.0, 0.0), (0.0, 0.0, -0.01), (0.0, 0.0, 0.02),
                (0.0, 0.0, -0.02), (0.1, 0.0, 0.0), (0.3, 0.0, 0.0)]:
        try:
            mip = mesh(*pt)
            el_idx = mip.nr if hasattr(mip, 'nr') else None
            el_mat = None
            try:
                el_mat = mesh[mip.mip.GetElement(VOL)].mat if hasattr(mip, 'mip') else None
            except Exception:
                pass
            g_v = np.array([gs(mip)[i] for i in range(3)])
            Mu_val = float(Mu(mip))
            sol_v = float(sol(mip))
            # try to get material name via another route
            print(f"  pt={pt}  sol={sol_v:+.2e}  Mu={Mu_val:.3e}  "
                  f"grad_z={g_v[2]:+.2e}  el_nr={el_idx} mat={el_mat}")
        except Exception as e:
            print(f"  pt={pt}  ERR {e}")

    print(f"\n  mesh.GetMaterials() indexed:")
    for i, m in enumerate(mats):
        print(f"    idx {i}: {m}")

    # Probe Bt (iron total) vs Br (air reduced) vs Bs individually
    print("\n[Bt/Br/Bs component probes]")
    fes_iron_p = H1(mesh, order=fe_order, definedon="iron")
    fes_air_p = H1(mesh, order=fe_order, definedon="air_inner|air_outer")
    Ot_p = GridFunction(fes_iron_p)
    Orr_p = GridFunction(fes_air_p)
    Oxr_p = GridFunction(fes_air_p)
    Ot_p.Set(sol)   # simplest form: fes_iron_p is already restricted to "iron"
    Orr_p.Set(sol)
    print(f"  Ot_p.vec nnz={np.count_nonzero(np.abs(np.asarray(Ot_p.vec.FV()))>1e-10)}"
          f"/{len(np.asarray(Ot_p.vec.FV()))}  "
          f"max|val|={np.max(np.abs(np.asarray(Ot_p.vec.FV()))):.3e}")
    print(f"  Orr_p.vec nnz={np.count_nonzero(np.abs(np.asarray(Orr_p.vec.FV()))>1e-10)}"
          f"/{len(np.asarray(Orr_p.vec.FV()))}  "
          f"max|val|={np.max(np.abs(np.asarray(Orr_p.vec.FV()))):.3e}")
    Oxr_p.Set(Omega_s_cf, BND, mesh.Boundaries("iron_surf"))
    Bt_cf = grad(Ot_p) * Mu
    Br_cf = (grad(Orr_p) - grad(Oxr_p)) * Mu
    for pt in [(0.0, 0.0, 0.0), (0.0, 0.0, -0.01), (0.0, 0.0, 0.05),
                (0.0, 0.0, 0.10)]:
        try:
            mip = mesh(*pt)
            Bt_v = np.array([Bt_cf(mip)[i] for i in range(3)])
            Br_v = np.array([Br_cf(mip)[i] for i in range(3)])
            try:
                Bs_v = np.array([Bs_cf(mip)[i] for i in range(3)])
            except Exception:
                Bs_v = np.zeros(3)
            Ot_v = float(Ot_p(mip))
            print(f"  pt={pt}  Ot={Ot_v:+.3e}  "
                  f"Bt_z={Bt_v[2]:+.3e}  Br_z={Br_v[2]:+.3e}  "
                  f"Bs_z={Bs_v[2]:+.3e}  sum={Bt_v[2]+Br_v[2]+Bs_v[2]:+.3e}")
        except Exception as e:
            print(f"  pt={pt}  ERR {e}")

    half_h = cylinder_height / 2
    z_inside = [(0.0, 0.0, zv) for zv in
                 np.linspace(-half_h + 0.01, half_h - 0.01, 5)]
    x_inside = [(xv, 0.0, 0.0) for xv in
                 np.linspace(0.005, cylinder_radius - 0.005, 4)]
    z_outside = [(0.0, 0.0, zv) for zv in
                  np.linspace(half_h + 0.02, phys_R - 0.04, 4)]

    rels_zi = compare_B(mesh, BField, grp_ref, z_inside,
                          "Bz along z-axis (inside cylinder)")
    rels_xi = compare_B(mesh, BField, grp_ref, x_inside,
                          "Bz along x-axis (inside cylinder, z=0)")
    rels_zo = compare_B(mesh, BField, grp_ref, z_outside,
                          "Bz along z-axis (outside cylinder, towards coil)")

    all_rels = np.array(rels_zi + rels_xi + rels_zo)
    print()
    print("=" * 72)
    print("SUMMARY (target |max rel| < 5 %)")
    print("=" * 72)
    print(f"  samples  = {len(all_rels)}")
    print(f"  max|rel| = {np.max(np.abs(all_rels))*100:.3f} %")
    print(f"  RMS rel  = {np.sqrt(np.mean(all_rels**2))*100:.3f} %")
    ok = np.max(np.abs(all_rels)) < 0.05
    print(f"  acceptance: {'OK' if ok else 'FAIL'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
