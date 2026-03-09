#!/usr/bin/env python
"""
C-Type Electromagnet 1/4 Model -- Simkin vs Cohomology Comparison

Quarter model (x>0, z>0) with symmetry boundary conditions:
  - x=0: Neumann (dPhi/dn=0, natural BC) -- B_n=Bx=0 on this plane
  - z=0: Dirichlet (phi=0) -- B_n=Bz!=0 on this plane
  - outer: Dirichlet (phi=0) -- Kelvin transform

Method A: Simkin reduced potential
  H = H_s - grad(phi), H_s from Radia coil (full coil, Biot-Savart)

Method B: Cohomology total scalar potential
  H = -grad(phi) + NI * h_k, h_k = cohomology basis (relative to Dirichlet bdry)

Both compared against Radia direct calculation.
"""

import sys
import os
import time
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(repo_root, 'src'))

MU_0 = 4 * np.pi * 1e-7

# =========================================================================
# Geometry constants (meters) -- full model
# =========================================================================
IRON_PARTS = [
    ('main_leg',     0.0, 0.13125, 0.0,    0.050, 0.0625, 0.210),
    ('yoke_back_lo', 0.0, 0.060,  -0.080,  0.050, 0.080,  0.050),
    ('yoke_back_hi', 0.0, 0.060,   0.080,  0.050, 0.080,  0.050),
    ('pole_bottom',  0.0, -0.000, -0.055,   0.050, 0.040,  0.100),
    ('pole_top',     0.0, -0.000,  0.055,   0.050, 0.040,  0.100),
]

TUBE_OUTER_CENTER = (0.0, 0.1325, 0.0)
TUBE_OUTER_SIZE = (0.110, 0.115, 0.008)
TUBE_INNER_CENTER = (0.0, 0.1325, 0.0)
TUBE_INNER_SIZE = (0.100, 0.105, 0.008)

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300
KELVIN_R = 0.600

NI = 2000.0
MU_R_IRON = 1000.0


# =========================================================================
# 1/4 Geometry builder
# =========================================================================

def build_quarter_geometry(gb, include_tube=False,
                           maxh_gap=0.003, maxh_iron=0.008, maxh_air=0.060):
    """Build 1/4 C-type geometry (x>0, z>0) using fragment_tracked.

    Parameters
    ----------
    gb : GmshBuilder
        Active builder context.
    include_tube : bool
        If True, include coil hole tube (for cohomology method).

    Returns
    -------
    dict with keys: iron, air, kelvin, (tube if include_tube)
    """
    import gmsh

    # --- Full iron (5 boxes fused) ---
    iron_ids = []
    for name, cx, cy, cz, sx, sy, sz in IRON_PARTS:
        iron_ids.append(gb.add_box([cx, cy, cz], [sx, sy, sz]))
    iron = gb.fuse(iron_ids)

    # --- Coil tube (cohomology only) ---
    tube = None
    if include_tube:
        outer = gb.add_box(list(TUBE_OUTER_CENTER), list(TUBE_OUTER_SIZE))
        inner = gb.add_box(list(TUBE_INNER_CENTER), list(TUBE_INNER_SIZE))
        tube = gb.cut(outer, inner)

    # --- Full spheres ---
    air = gb.add_sphere(SPHERE_CENTER, AIR_R)
    kelvin = gb.add_sphere(SPHERE_CENTER, KELVIN_R)

    # --- Quarter cutting box (x>=0, z>=0) ---
    L = KELVIN_R * 2.5
    quarter = gb.add_box([L/2, 0, L/2], [L, 4*L, L])

    # --- Fragment all ---
    all_vols = [iron, air, kelvin, quarter]
    if tube is not None:
        all_vols.insert(1, tube)

    frag_map = gb.fragment_tracked(all_vols)

    iron_set = set(frag_map[iron])
    air_set = set(frag_map[air])
    kelvin_set = set(frag_map[kelvin])
    quarter_set = set(frag_map[quarter])
    tube_set = set(frag_map[tube]) if tube is not None else set()

    # --- Classify quarter volumes ---
    iron_q = sorted(iron_set & quarter_set)
    tube_q = sorted((tube_set & quarter_set) - iron_set) if tube else []
    air_q = sorted((air_set & quarter_set) - iron_set - tube_set)
    kelvin_q = sorted((kelvin_set & quarter_set) - air_set)

    # --- Remove outside volumes ---
    keep = set(iron_q + air_q + kelvin_q + tube_q)
    all_frags = set()
    for v in frag_map.values():
        all_frags.update(v)
    remove = all_frags - keep
    for vol in remove:
        try:
            gmsh.model.occ.remove([(3, vol)], recursive=False)
        except Exception:
            pass
    gmsh.model.occ.synchronize()

    print(f"  Iron: {len(iron_q)}, Air: {len(air_q)}, "
          f"Kelvin: {len(kelvin_q)}", end="")
    if tube:
        print(f", Tube: {len(tube_q)}", end="")
    print()

    # --- Physical groups ---
    gb.add_block_by_name(iron_q, 'iron')
    gb.add_block_by_name(air_q, 'air')
    gb.add_block_by_name(kelvin_q, 'kelvin')
    if tube_q:
        gb.add_block_by_name(tube_q, 'coil_hole')

    domain_ids = iron_q + air_q + kelvin_q
    gb.add_block_by_name(domain_ids, 'domain')

    # --- Identify boundary surfaces ---
    all_surfs = set()
    for vid in keep:
        all_surfs.update(gb.get_surfaces(vid))

    sym_z_surfs = []
    sym_x_surfs = []
    kelvin_outer_surfs = set()
    for vid in kelvin_q:
        kelvin_outer_surfs.update(gb.get_surfaces(vid))
    air_surfs_set = set()
    for vid in air_q + iron_q + tube_q:
        air_surfs_set.update(gb.get_surfaces(vid))

    # Outer kelvin = kelvin surfaces not shared with interior
    outer_surfs = sorted(kelvin_outer_surfs - air_surfs_set)

    # Classify flat faces by bounding box
    tol = 1e-6
    for sid in all_surfs:
        try:
            bbox = gmsh.model.occ.getBoundingBox(2, sid)
        except Exception:
            continue
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        if abs(zmin) < tol and abs(zmax) < tol:
            sym_z_surfs.append(sid)
        elif abs(xmin) < tol and abs(xmax) < tol:
            sym_x_surfs.append(sid)

    # Remove sym_z and sym_x from outer
    outer_surfs = sorted(set(outer_surfs) - set(sym_z_surfs) - set(sym_x_surfs))

    if sym_z_surfs:
        gb.add_sideset_by_name(sym_z_surfs, 'sym_z')
    if sym_x_surfs:
        gb.add_sideset_by_name(sym_x_surfs, 'sym_x')
    if outer_surfs:
        gb.add_sideset_by_name(outer_surfs, 'outer')

    # NOTE: do NOT create a separate 'dirichlet' physical group because
    # surfaces in multiple groups may lose their label during export.
    # Instead, use 'sym_z|outer' in H1 FES dirichlet specification.

    print(f"  sym_z: {len(sym_z_surfs)}, sym_x: {len(sym_x_surfs)}, "
          f"outer: {len(outer_surfs)} surfaces")

    # --- Mesh control ---
    f_gap = gb.add_field_box(
        0.0, 0.030, -0.025, 0.025, 0.0, 0.010,
        size_in=maxh_gap, size_out=maxh_air)

    f_iron = gb.add_field_box(
        0.0, 0.030, -0.025, 0.170, 0.0, 0.110,
        size_in=maxh_iron, size_out=maxh_air)

    f_min = gb.add_field_min([f_gap, f_iron])
    gb.set_background_field(f_min)
    gb.set_min_max_size(maxh_gap / 3, maxh_air)
    gb.set_algorithm_3d(1)

    result = {
        'iron': iron_q, 'air': air_q, 'kelvin': kelvin_q,
    }
    if tube:
        result['tube'] = tube_q
    return result


# =========================================================================
# Radia coil (FULL coil for H_s)
# =========================================================================

def build_radia_coil():
    """Build full Radia racetrack coil (not 1/4)."""
    import radia as rad

    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'mu=1000')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil

    return create_racetrack_coil(NI)


# (Simkin solver now uses ScalarPotentialSolver from scalar_potential_solver.py)


# =========================================================================
# Radia reference (full model)
# =========================================================================

def build_radia_reference():
    """Build Radia full model for comparison."""
    import radia as rad

    rad.UtiDelAll()
    blocks = []
    for name, cx, cy, cz, dx, dy, dz in IRON_PARTS:
        x0, y0, z0 = cx - dx/2, cy - dy/2, cz - dz/2
        nx, ny, nz = 2, 2, max(2, int(round(dz / dx * 2)))
        sx, sy, sz = dx / nx, dy / ny, dz / nz
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    bx = x0 + (ix + 0.5) * sx
                    by = y0 + (iy + 0.5) * sy
                    bz = z0 + (iz + 0.5) * sz
                    block = rad.ObjRecMag([bx, by, bz],
                                         [sx, sy, sz], [0, 0, 0])
                    mat = rad.MatLin(MU_R_IRON)
                    rad.MatApl(block, mat)
                    blocks.append(block)

    iron = rad.ObjCnt(blocks)
    coil = build_radia_coil()
    model = rad.ObjCnt([iron, coil])

    result = rad.Solve(model, 0.0001, 1000, 0)
    print(f"   Radia: {len(blocks)} blocks, max|dM|={result[1]:.2e}")
    return model


# =========================================================================
# Comparison table
# =========================================================================

def compare_fields(mesh, B_cf, radia_model, label=""):
    """Print B comparison at test points."""
    import radia as rad

    test_points = [
        ([0, 0, 0],      "gap center"),
        ([0, 0, 0.003],  "z=3mm"),
        ([0.010, 0, 0],  "x=10mm"),
        ([0, 0.010, 0],  "y=10mm"),
    ]

    print(f"   {'Point':16s} {'FEM (mT)':>10s} {'Radia (mT)':>12s} {'Diff':>8s}")
    print("   " + "-" * 52)

    for pt, desc in test_points:
        B_radia = np.array(rad.Fld(radia_model, 'b', pt))
        B_r_mag = np.linalg.norm(B_radia)
        try:
            mip = mesh(*pt)
            B_fem_pt = np.array(B_cf(mip))
            B_f_mag = np.linalg.norm(B_fem_pt)
            pct = (B_f_mag - B_r_mag) / B_r_mag * 100 if B_r_mag > 1e-6 else 0
            print(f"   {desc:16s} {B_f_mag*1e3:8.2f}   {B_r_mag*1e3:10.2f}   {pct:+6.1f}%")
        except Exception:
            print(f"   {desc:16s} {'(outside)':>8s}   {B_r_mag*1e3:10.2f}")


# =========================================================================
# Main
# =========================================================================

def main():
    print("=" * 70)
    print("  C-Type Electromagnet: 1/4 Model Comparison")
    print("  Simkin (reduced potential) vs Cohomology (total potential)")
    print("=" * 70)
    print()
    print(f"NI = {NI:.0f} AT, mu_r = {MU_R_IRON:.0f}")
    print(f"Symmetry: x=0 (Neumann), z=0 (Dirichlet phi=0)")
    print()

    import radia as rad
    from radia.gmsh_builder import GmshBuilder
    from radia.scalar_potential_solver import ScalarPotentialSolver

    fem_order = 3

    # =================================================================
    # Method A: Simkin 1/4 (using ScalarPotentialSolver)
    # =================================================================
    print("=" * 70)
    print("  Method A: Simkin Reduced Potential (1/4)")
    print("=" * 70)
    print()

    print("A1. Building 1/4 geometry (no coil hole)...")
    with GmshBuilder(model_name='quarter_simkin', verbose=False) as gb:
        vols = build_quarter_geometry(gb, include_tube=False,
                                      maxh_gap=0.003,
                                      maxh_iron=0.008,
                                      maxh_air=0.060)
        gb.generate(element_type='tet')
        mesh_s = gb.to_ngsolve_volume()

    mesh_s.Curve(fem_order)
    print(f"   Mesh: {mesh_s.nv} vertices, {mesh_s.ne} elements")
    print(f"   Materials: {sorted(set(mesh_s.GetMaterials()))}")
    print(f"   Boundaries: {sorted(set(mesh_s.GetBoundaries()))}")
    print()

    print("A2. Computing H_s from full Radia coil...")
    rad.UtiDelAll()
    coil = build_radia_coil()

    simkin = ScalarPotentialSolver(
        mesh_s, mu_r_dict={'iron': MU_R_IRON}, order=fem_order,
        kelvin_region='kelvin', kelvin_radius=AIR_R,
        kelvin_center=SPHERE_CENTER)
    simkin.set_source_from_radia(coil, resolution=61)
    print()

    print("A3. Solving Simkin + Kelvin (1/4)...")
    t0 = time.time()
    phi_s = simkin.solve(dirichlet='sym_z|outer')
    t_simkin = time.time() - t0
    B_s = simkin.get_B()
    print(f"   Solve: {t_simkin:.1f}s")
    print(f"   phi range: [{min(phi_s.vec):.2f}, {max(phi_s.vec):.2f}]")
    print(f"   ndof: {phi_s.space.ndof}, free: {sum(phi_s.space.FreeDofs())}")

    try:
        mip = mesh_s(0.001, 0, 0.001)
        B_val = B_s(mip)
        print(f"   B at (1,0,1)mm: Bz={B_val[2]*1e3:.2f} mT, |B|={np.sqrt(sum(b**2 for b in B_val))*1e3:.2f} mT")
    except Exception as e:
        print(f"   Cannot evaluate at origin: {e}")
    print()

    # =================================================================
    # Method B: Cohomology 1/4
    # =================================================================
    print("=" * 70)
    print("  Method B: Cohomology Total Potential (1/4)")
    print("=" * 70)
    print()

    print("B1. Building 1/4 geometry (with coil hole)...")
    from radia.cohomology_cut import CohomologyCutSolver

    with GmshBuilder(model_name='quarter_cohomology', verbose=False) as gb:
        vols = build_quarter_geometry(gb, include_tube=True,
                                      maxh_gap=0.003,
                                      maxh_iron=0.008,
                                      maxh_air=0.060)

        print()
        print("B2. Computing cohomology (relative to Dirichlet bdry)...")
        solver = CohomologyCutSolver()
        try:
            n_coils = solver.setup_from_gmsh(
                domain_physical_name='domain',
                boundary_physical_name='outer')
        except RuntimeError as e:
            print(f"   Cohomology failed: {e}")
            n_coils = 0

    if n_coils >= 1:
        mesh_c = solver.get_mesh()
        mesh_c.Curve(fem_order)
        print(f"   Cohomology generators: {n_coils}")
        print(f"   Mesh: {mesh_c.nv} vertices, {mesh_c.ne} elements")
    else:
        mesh_c = None

    if n_coils >= 1:
        print()
        print("B3. Solving total potential + Kelvin (1/4)...")
        t0 = time.time()
        phi_c = solver.solve(
            [NI],
            mu_r_dict={'iron': MU_R_IRON},
            order=fem_order,
            dirichlet='sym_z|outer',
            kelvin_region='kelvin',
            kelvin_radius=AIR_R,
            kelvin_center=SPHERE_CENTER)
        t_cohom = time.time() - t0
        print(f"   Solve: {t_cohom:.1f}s")
        B_c = solver.get_B()

        try:
            mip = mesh_c(0, 0, 0)
            B_val = B_c(mip)
            print(f"   B at origin: Bz={B_val[2]*1e3:.2f} mT, "
                  f"|B|={np.sqrt(sum(b**2 for b in B_val))*1e3:.2f} mT")
        except Exception as e:
            print(f"   Cannot evaluate at origin: {e}")
    else:
        print(f"   Cohomology: {n_coils} generators (1/4 tube is not a closed loop)")
        print("   1/4 cohomology requires relative cohomology (future work)")
        B_c = None
        t_cohom = 0
    print()

    # =================================================================
    # Radia reference
    # =================================================================
    print("=" * 70)
    print("  Radia Reference (full model)")
    print("=" * 70)
    print()

    radia_model = build_radia_reference()
    print()

    print("--- Simkin 1/4 vs Radia ---")
    compare_fields(mesh_s, B_s, radia_model, "Simkin")
    print()

    if B_c is not None:
        print("--- Cohomology 1/4 vs Radia ---")
        compare_fields(mesh_c, B_c, radia_model, "Cohomology")
        print()

    rad.UtiDelAll()

    # =================================================================
    # Summary
    # =================================================================
    print("=" * 70)
    print("  Summary: 1/4 Model Comparison")
    print("=" * 70)
    print(f"  Simkin:     {mesh_s.ne:,} elements, {t_simkin:.1f}s solve")
    if B_c is not None:
        print(f"  Cohomology: {mesh_c.ne:,} elements, {t_cohom:.1f}s solve")
    print()
    print("  Symmetry BCs (phi formulation):")
    print("    x=0: Neumann (natural) -- reversed from A formulation")
    print("    z=0: Dirichlet (phi=0) -- reversed from A formulation")
    print("    outer: Dirichlet (phi=0) -- Kelvin transform")
    print("=" * 70)


if __name__ == '__main__':
    main()
