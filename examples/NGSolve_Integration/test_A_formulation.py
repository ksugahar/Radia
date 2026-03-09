#!/usr/bin/env python
"""
Test: NGSolve A-formulation (vector potential) for C-type electromagnet.

Solves: curl(nu(|curl A|) * curl A) = J
where nu = 1/mu is reluctivity, derived from B-H curve.

This avoids the Simkin reduced scalar potential decomposition entirely.
No H_s computation, no cancellation issues.

Reference: BEM (Tosca-verified) Bz = -995 mT at origin, NI=20000 AT.
"""

import sys
import os
import time
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
sys.path.insert(0, os.path.join(repo_root, 'src'))

from ngsolve import *
from ngsolve.solvers import Newton
import ngsolve

MU_0 = 4 * np.pi * 1e-7

# =========================================================================
# Geometry: Same ELF polygon yoke as test_elf_yoke_simkin.py
# =========================================================================
YOKE_POLYGON_YZ = [
    (-0.020, -0.105),
    ( 0.1625, -0.105),
    ( 0.1625,  0.105),
    (-0.020,  0.105),
    (-0.020,  0.005),
    ( 0.020,  0.005),
    ( 0.020,  0.055),
    ( 0.100,  0.055),
    ( 0.100, -0.055),
    ( 0.020, -0.055),
    ( 0.020, -0.005),
    (-0.020, -0.005),
]

YOKE_X_MIN = -0.025
YOKE_X_MAX =  0.025
YOKE_X_EXTENT = YOKE_X_MAX - YOKE_X_MIN

SPHERE_CENTER = [0.0, 0.070, 0.0]
AIR_R = 0.300

# Coil parameters (racetrack coil, same as ELF model)
# Coil cross-section: 35mm x 105mm
# Coil center Y=131.25mm (between y=20mm arm inner and y=100mm arm outer,
#   coil spans y=76.25 to y=186.25 mm... no.
# From coil_model.py: center (0, 0.13125, 0), cross-section 35x105 mm
# So coil goes: Y = 131.25-52.5 to 131.25+52.5 = 78.75 to 183.75 mm
#               Z = 0-17.5 to 0+17.5 = -17.5 to 17.5 mm
# X extent same as yoke: -25 to 25 mm
# For FEM, we only need the straight portions (Y-direction current)
# Upper leg: X=[-25,25], Y=[78.75,183.75], Z=[17.5] (current in -Y... no)
# Actually racetrack has current flowing in Y direction along two sides.
# The two sides at Z=+17.5 and Z=-17.5:
#   Side 1 (Z=+17.5): current in +Y direction
#   Side 2 (Z=-17.5): current in -Y direction
# But this is racetrack - current goes around, so:
#   Front (closest to gap): Z=[-17.5, +17.5], current in +X or -X?
# Let me re-read the coil model.

COIL_Y_MIN = 0.07875    # 78.75 mm
COIL_Y_MAX = 0.18375    # 183.75 mm
COIL_Z_MIN = -0.0175    # -17.5 mm
COIL_Z_MAX =  0.0175    # +17.5 mm
COIL_X_MIN = YOKE_X_MIN
COIL_X_MAX = YOKE_X_MAX


def load_bh_curve():
    """Load B-H curve from BH.txt. Returns (H_array, B_array)."""
    bh_file = os.path.join(repo_root, 'examples',
                           'c_type_electromagnet', 'nonlinear', 'BH.txt')
    data = np.loadtxt(bh_file)
    return data[:, 0], data[:, 1]


def build_geometry_and_mesh(maxh_gap=0.004, maxh_iron=0.010, maxh_air=0.050):
    """Build geometry with GmshBuilder, including coil region."""
    from radia.gmsh_builder import GmshBuilder

    with GmshBuilder(model_name='a_form_test', verbose=False) as gb:
        # Iron yoke (polygon extrusion)
        n = len(YOKE_POLYGON_YZ)
        point_tags = []
        for y, z in YOKE_POLYGON_YZ:
            pt = gb.add_point(YOKE_X_MIN, y, z)
            point_tags.append(pt)

        line_tags = []
        for i in range(n):
            j = (i + 1) % n
            line = gb.add_line(point_tags[i], point_tags[j])
            line_tags.append(line)

        wire = gb.add_wire(line_tags)
        surf = gb.add_plane_surface([wire])
        extruded = gb.extrude_surface(surf, YOKE_X_EXTENT, 0, 0)
        iron_vol = extruded[0]

        # Coil regions (two rectangular cross-section volumes in coil window)
        # Upper coil (current in +X): Z=[COIL_Z_MAX-small, COIL_Z_MAX] ... no
        # Actually for A-formulation, we model the coil as a volume with J.
        # The racetrack coil in ELF/Radia has current in Y direction
        # on the two Z-legs, and current in X direction on the two end caps.
        # For a 3D model, let's define the full rectangular coil cross-section.
        # Coil cross-section area: 35mm(Z) x 50mm(X) = 1750 mm^2
        # But coil_model.py says 35x105mm cross-section...
        # Let me just use the coil as J source directly (no coil mesh region).
        # Simpler: define J as a CoefficientFunction.

        # Air sphere
        air = gb.add_sphere(SPHERE_CENTER, AIR_R)
        frag_map = gb.fragment_tracked([iron_vol, air])

        iron_ids = set(frag_map[iron_vol])
        air_ids = set(frag_map[air])

        iron_final = sorted(iron_ids)
        air_final = sorted(air_ids - iron_ids)

        gb.add_block_by_name(iron_final, 'iron')
        gb.add_block_by_name(air_final, 'air')

        air_surfs = set()
        for vid in air_final:
            air_surfs.update(gb.get_surfaces(vid))
        iron_surfs = set()
        for vid in iron_final:
            iron_surfs.update(gb.get_surfaces(vid))
        outer_surfs = sorted(air_surfs - iron_surfs)
        if outer_surfs:
            gb.add_sideset_by_name(outer_surfs, 'outer')

        print(f"  Iron: {len(iron_final)}, Air: {len(air_final)}")
        print(f"  Outer boundary: {len(outer_surfs)} surfaces")

        # Mesh size fields
        f_gap = gb.add_field_box(
            -0.030, 0.030, -0.025, 0.025, -0.010, 0.010,
            size_in=maxh_gap, size_out=maxh_air)
        f_iron = gb.add_field_box(
            -0.030, 0.030, -0.025, 0.170, -0.110, 0.110,
            size_in=maxh_iron, size_out=maxh_air)
        f_min = gb.add_field_min([f_gap, f_iron])
        gb.set_background_field(f_min)
        gb.set_min_max_size(maxh_gap / 3, maxh_air)
        gb.set_algorithm_3d(1)

        gb.generate(element_type='tet')
        mesh = gb.to_ngsolve_volume()

    return mesh


def create_coil_source_cf(ni):
    """Create current density J as CoefficientFunction for the racetrack coil.

    The racetrack coil has rectangular cross-section and carries current
    around the coil window. For 3D FEM, we model J as a volume current.

    Racetrack coil (from coil_model.py):
      Center: (0, 0.13125, 0)
      Cross-section: 35mm (Z) x 105mm (Y/along-coil)
      But actually it's a RACETRACK (ring), so current flows around.

    For the A-formulation, we use the Biot-Savart H_s from Radia
    coil and compute curl(H_s) = J_free. But that's circular.

    Better approach: define J directly in the coil volume.
    The coil is a rectangular ring in the Y-Z plane:
      Outer: Y=[78.75, 183.75], Z=[-17.5, 17.5]
      But it's a racetrack, so current flows in closed loop.

    Actually, let's use VoxelCoefficient for H_s and compute
    J = curl(H_s). No - that would couple back to Simkin.

    Simplest correct approach: use Radia's ObjRaceTrk coil to compute
    H_s on a grid, then use VoxelCoefficient as RHS for A-formulation:
        curl(nu * curl A) = curl(H_s) in air
    No - that's wrong. The RHS should be J_free (the free current density).

    For A-formulation with iron:
        curl(nu * curl A) = J_coil

    where J_coil is the current density in the coil region.
    The coil is NOT in the mesh (it's in air), so J_coil is a body force.

    Let me define J_coil as rectangular current loops.
    Racetrack coil has current I flowing through cross-section A_coil.
    J = NI / A_coil (magnitude).

    From coil_model.py: create_racetrack_coil creates a Radia ObjRaceTrk.
    The coil has inner/outer dimensions and current.
    We need to know: which direction does J flow?

    Looking at Radia ObjRaceTrk: current flows in a loop in the coil.
    For a racetrack in the Y-Z plane:
    - Two straight legs parallel to Y (at Z=+17.5mm and Z=-17.5mm)
    - Two semicircular ends

    Wait - the coil_model.py says:
    coil = rad.ObjRaceTrk(center, radii, straights, height, nseg, cur_dens, 'man', 'y')
    The 'y' means the symmetry axis is Y. So current flows in loops
    around Y axis. The straight sections are parallel to Y.

    For this coil geometry with 'y' axis:
    - Current loops are in the X-Z plane (perpendicular to Y)
    - Wait, that doesn't make sense for a magnet...

    Let me check the actual coil_model.py parameters.
    """
    # We'll use VoxelCoefficient from Radia coil for the source.
    # For A-formulation: curl(nu * curl A) = J
    # We cannot easily define J_coil analytically for a racetrack.
    #
    # Alternative approach: Use a "two-potential" method instead.
    # Or: compute H_s from Radia, then solve for A using:
    #   curl(nu * curl A) = curl(mu_0 * H_s)  ... NO, that's wrong.
    #
    # The correct A-formulation with source field:
    #   B_total = B_s + curl(A_r)
    #   where B_s = mu_0 * H_s (Biot-Savart field of coil alone in vacuum)
    #   and A_r is the "reduced" vector potential due to iron.
    #
    #   Weak form: (nu * (B_s + curl A_r), curl v) = 0  for all v
    #   => (nu * curl A_r, curl v) = -(nu * B_s, curl v)
    #
    # This IS the A-formulation equivalent of Simkin!
    # H_total = nu * B_total = nu * (B_s + curl A_r)
    # In air:  nu = 1/mu_0, so H = H_s + curl(A_r)/mu_0
    # In iron: nu = nu(|B|), nonlinear
    #
    # The advantage: we work with B directly, no cancellation in H.
    pass


def solve_A_formulation(mesh, ni, bh_data_H, bh_data_B, order=2,
                        maxiter=50, tol=1e-4):
    """Solve nonlinear magnetostatic problem with A-formulation.

    Uses reduced vector potential: B = B_s + curl(A_r)
    where B_s = mu_0 * H_s from Radia coil.

    Weak form (energy minimization):
      W(A_r) = integral nu(|B_s + curl A_r|) * |B_s + curl A_r|^2 / 2 dx
      or equivalently via co-energy:
      W(A_r) = integral w*(|B_s + curl A_r|) dx

    where w*(B) = integral_0^B H(B') dB'  (co-energy density in terms of B).

    For iron: H(B) from inverted B-H curve.
    For air:  H(B) = B / mu_0.
    """
    import radia as rad

    print("\n--- A-formulation: B = B_s + curl(A_r) ---")

    # Compute B_s from Radia coil
    print("Computing B_s from Radia coil...")
    rad.UtiDelAll()
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'nonlinear')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    coil = create_racetrack_coil(ni)

    # Use VoxelCoefficient for B_s = mu_0 * H_s
    from ngsolve import VoxelCoefficient, CF as CF_
    resolution = 51

    # Bounding box from mesh
    pmin, pmax = mesh.ngmesh.bounding_box
    pmin_l = [pmin[i] for i in range(3)]
    pmax_l = [pmax[i] for i in range(3)]
    margin = 0.05 * max(pmax_l[i] - pmin_l[i] for i in range(3))
    bbox = [[pmin_l[i] - margin, pmax_l[i] + margin] for i in range(3)]

    nx = ny = nz = resolution
    xs = np.linspace(bbox[0][0], bbox[0][1], nx)
    ys = np.linspace(bbox[1][0], bbox[1][1], ny)
    zs = np.linspace(bbox[2][0], bbox[2][1], nz)

    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    print(f"  Computing H_s on {nx}x{ny}x{nz} grid...")
    Hs_field = np.asarray(rad.Fld(coil, 'h', points))
    Bs_field = MU_0 * Hs_field

    start = (bbox[0][0], bbox[1][0], bbox[2][0])
    end = (bbox[0][1], bbox[1][1], bbox[2][1])

    cfs_n = []
    for comp in range(3):
        data = Bs_field[:, comp].reshape(nx, ny, nz)
        data = np.ascontiguousarray(data.transpose(2, 1, 0))
        cfs_n.append(VoxelCoefficient(start, end, data, linear=True))

    Bs_cf = CF_(tuple(cfs_n))

    # Check B_s at origin
    mip_origin = mesh(0, 0, 0)
    Bs_origin = Bs_cf(mip_origin)
    print(f"  B_s at origin: ({Bs_origin[0]*1e3:.2f}, {Bs_origin[1]*1e3:.2f}, "
          f"{Bs_origin[2]*1e3:.2f}) mT")

    rad.UtiDelAll()

    # Build nu(B) from B-H curve
    # nu = H / B (reluctivity)
    # For co-energy: w*(B) = integral_0^B H(B') dB'
    print("Building nu(B) from B-H curve...")

    # Invert B-H curve to get H(B)
    # bh_data_H, bh_data_B are arrays of H and B values
    # We need H as function of B (nu = H/B)
    # Sort by B
    sort_idx = np.argsort(bh_data_B)
    B_sorted = bh_data_B[sort_idx]
    H_sorted = bh_data_H[sort_idx]

    # Remove duplicate B values
    mask = np.diff(B_sorted, prepend=-1) > 0
    B_unique = B_sorted[mask]
    H_unique = H_sorted[mask]

    # nu(B) = H(B) / B
    # For B=0: nu = 1/(mu_0 * mu_r_initial)
    # Use mu_r from initial slope: mu_r = B[1]/(mu_0*H[1])
    if H_unique[0] == 0:
        mu_r_init = B_unique[1] / (MU_0 * H_unique[1]) if H_unique[1] > 0 else 1000.0
    else:
        mu_r_init = B_unique[0] / (MU_0 * H_unique[0])
    nu_0_iron = 1.0 / (MU_0 * mu_r_init)
    print(f"  Initial mu_r = {mu_r_init:.1f}, nu_0_iron = {nu_0_iron:.2e}")

    # Co-energy density w*(B) = integral_0^B H(B') dB' (for iron)
    # Numerically integrate
    n_bh = len(B_unique)
    w_star = np.zeros(n_bh)
    for i in range(1, n_bh):
        # Trapezoidal rule
        dB = B_unique[i] - B_unique[i-1]
        w_star[i] = w_star[i-1] + 0.5 * (H_unique[i] + H_unique[i-1]) * dB

    # Extend to high B values (linear extrapolation with nu_0 = 1/mu_0)
    B_max_data = B_unique[-1]
    nu_max = H_unique[-1] / B_unique[-1] if B_unique[-1] > 0 else 1.0/MU_0
    print(f"  B-H data range: B = [0, {B_max_data:.3f}] T")
    print(f"  nu at B_max = {nu_max:.2e} (cf. 1/mu_0 = {1/MU_0:.2e})")

    # Create BSpline for w*(B) in iron
    # NGSolve BSpline: BSpline(order, knots, values)
    # We need w* as function of B^2 for easier chain rule? No, use B directly.

    # For NGSolve SymbolicEnergy, we need w*(|B|) as a function of |B|.
    # BSpline maps scalar -> scalar.

    # Extend data if needed
    B_ext = list(B_unique)
    w_ext = list(w_star)
    if B_max_data < 5.0:
        # Extend to 5T with linear H(B) slope
        dH_dB_end = (H_unique[-1] - H_unique[-2]) / (B_unique[-1] - B_unique[-2])
        for B_val in np.linspace(B_max_data + 0.1, 5.0, 30):
            H_val = H_unique[-1] + dH_dB_end * (B_val - B_max_data)
            dB = B_val - B_ext[-1]
            H_prev = H_unique[-1] + dH_dB_end * (B_ext[-1] - B_max_data)
            w_val = w_ext[-1] + 0.5 * (H_val + H_prev) * dB
            B_ext.append(B_val)
            w_ext.append(w_val)

    B_ext = np.array(B_ext)
    w_ext = np.array(w_ext)

    # Create NGSolve BSpline for w*(B)
    # BSpline needs: order, knots (with multiplicity), control points
    # Simpler: use BSpline.Integrate() ... no, we already have w*(B).
    # Use BSpline(order=3, knot_vector, values) where values are at knot points.
    from ngsolve import BSpline

    # BSpline(order, [t0, t1, ...], [v0, v1, ...])
    w_star_bspline = BSpline(3, list(B_ext), list(w_ext))

    # For air: w*(B) = B^2 / (2*mu_0)
    # For iron: w*(B) = w_star_bspline(B)

    # Setup FE space: HCurl (Nedelec) for vector potential A
    fes = HCurl(mesh, order=order, dirichlet='outer')
    print(f"  HCurl DOFs: {fes.ndof}")

    A = fes.TrialFunction()
    v = fes.TestFunction()

    # Domain indicators
    iron_indicator = CoefficientFunction([1 if mat == 'iron' else 0
                                           for mat in mesh.GetMaterials()])
    air_indicator = 1.0 - iron_indicator

    # B_total = B_s + curl(A_r)
    B_total = Bs_cf + curl(A)

    # |B_total|
    B_mag = Norm(B_total)
    B_mag_safe = IfPos(B_mag - 1e-15, B_mag, 1e-15)

    # Energy functional:
    # W = integral_iron w*(|B_total|) dx + integral_air |B_total|^2/(2*mu_0) dx

    # Iron co-energy
    W_iron = w_star_bspline(B_mag_safe) * iron_indicator * dx

    # Air co-energy: B^2 / (2*mu_0)
    W_air = (1.0 / (2.0 * MU_0)) * InnerProduct(B_total, B_total) * air_indicator * dx

    # Total energy
    W_total = W_iron + W_air

    # GridFunction for solution
    A_gf = GridFunction(fes)
    A_gf.vec[:] = 0

    # Newton solve using SymbolicEnergy
    print("\nNewton solve...")
    t0 = time.time()

    try:
        Newton(W_total, A_gf, maxerr=tol, maxit=maxiter, printing=True,
               dampfactor=0.25, inverse='pardiso')
        converged = True
    except Exception as e:
        print(f"  Newton exception: {e}")
        converged = False

    t_solve = time.time() - t0
    print(f"  Solve time: {t_solve:.1f}s")

    # Evaluate B at origin
    B_total_cf = Bs_cf + curl(A_gf)
    B_origin = B_total_cf(mip_origin)
    Bz = B_origin[2]
    B_mag_origin = np.sqrt(sum(b**2 for b in B_origin))

    print(f"\n  B at origin: Bx={B_origin[0]*1e3:.2f}, "
          f"By={B_origin[1]*1e3:.2f}, Bz={B_origin[2]*1e3:.2f} mT")
    print(f"  |B| = {B_mag_origin*1e3:.2f} mT")
    print(f"  Expected (BEM/Tosca): Bz = -995 mT")
    if abs(Bz) > 1e-6:
        print(f"  Error: {abs(Bz*1e3 - (-995)) / 995 * 100:.1f}%")

    return Bz, A_gf, B_total_cf


def solve_A_formulation_picard(mesh, ni, bh_data_H, bh_data_B, order=2,
                                maxiter=80, tol=1e-4):
    """Picard iteration for A-formulation (more robust than Newton).

    Solves: (nu_k * curl A_{k+1}, curl v) = -(nu_k * B_s, curl v)
    where nu_k = H(|B_k|) / |B_k| is updated from previous iteration.
    """
    import radia as rad
    from scipy.interpolate import interp1d

    print("\n--- A-formulation Picard: B = B_s + curl(A_r) ---")

    # Compute B_s from Radia coil
    print("Computing B_s from Radia coil...")
    rad.UtiDelAll()
    coil_dir = os.path.join(repo_root, 'examples',
                            'c_type_electromagnet', 'nonlinear')
    if coil_dir not in sys.path:
        sys.path.insert(0, coil_dir)
    from coil_model import create_racetrack_coil
    coil = create_racetrack_coil(ni)

    from ngsolve import VoxelCoefficient, CF
    resolution = 51

    # Bounding box from mesh
    pmin, pmax = mesh.ngmesh.bounding_box
    pmin = [pmin[i] for i in range(3)]
    pmax = [pmax[i] for i in range(3)]
    margin = 0.05 * max(pmax[i] - pmin[i] for i in range(3))
    bbox = [[pmin[i] - margin, pmax[i] + margin] for i in range(3)]

    nx = ny = nz = resolution
    xs = np.linspace(bbox[0][0], bbox[0][1], nx)
    ys = np.linspace(bbox[1][0], bbox[1][1], ny)
    zs = np.linspace(bbox[2][0], bbox[2][1], nz)

    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    print(f"  Computing H_s on {nx}x{ny}x{nz} grid...")
    Hs_field = np.asarray(rad.Fld(coil, 'h', points))

    # B_s = mu_0 * H_s
    Bs_field = MU_0 * Hs_field

    start = (bbox[0][0], bbox[1][0], bbox[2][0])
    end = (bbox[0][1], bbox[1][1], bbox[2][1])

    cfs = []
    for comp in range(3):
        data = Bs_field[:, comp].reshape(nx, ny, nz)
        data = np.ascontiguousarray(data.transpose(2, 1, 0))
        cfs.append(VoxelCoefficient(start, end, data, linear=True))

    Bs_cf = CF(tuple(cfs))

    mip_origin = mesh(0, 0, 0)
    Bs_origin = Bs_cf(mip_origin)
    print(f"  B_s at origin: ({Bs_origin[0]*1e3:.2f}, {Bs_origin[1]*1e3:.2f}, "
          f"{Bs_origin[2]*1e3:.2f}) mT")

    rad.UtiDelAll()

    # Build nu(B) interpolation from B-H curve
    # nu = H/B (reluctivity)
    sort_idx = np.argsort(bh_data_B)
    B_sorted = bh_data_B[sort_idx]
    H_sorted = bh_data_H[sort_idx]

    mask = np.diff(B_sorted, prepend=-1) > 0
    B_unique = B_sorted[mask]
    H_unique = H_sorted[mask]

    # nu(B) = H(B)/B, handle B=0
    nu_arr = np.zeros_like(B_unique)
    if B_unique[0] == 0 and H_unique[0] == 0:
        # Use initial slope
        nu_arr[0] = H_unique[1] / B_unique[1] if B_unique[1] > 0 else 1.0/(MU_0*1000)
    else:
        nu_arr[0] = H_unique[0] / B_unique[0] if B_unique[0] > 0 else 1.0/(MU_0*1000)

    for i in range(1, len(B_unique)):
        nu_arr[i] = H_unique[i] / B_unique[i]

    # Extend to high B
    B_ext = list(B_unique)
    nu_ext = list(nu_arr)
    if B_unique[-1] < 5.0:
        dH_dB_end = (H_unique[-1] - H_unique[-2]) / (B_unique[-1] - B_unique[-2])
        for B_val in np.linspace(B_unique[-1] + 0.1, 10.0, 50):
            H_val = H_unique[-1] + dH_dB_end * (B_val - B_unique[-1])
            nu_ext.append(H_val / B_val)
            B_ext.append(B_val)

    nu_interp = interp1d(B_ext, nu_ext, kind='linear', fill_value='extrapolate')

    mu_r_init = 1.0 / (MU_0 * nu_ext[0])
    print(f"  Initial mu_r = {mu_r_init:.1f}")
    print(f"  nu(B=0) = {nu_ext[0]:.2e}, nu(B_max) = {nu_ext[-1]:.2e}")

    # B_sat for convergence check
    B_sat = B_unique[-1]
    print(f"  B_sat = {B_sat:.3f} T")

    # FE space: HCurl for A_r
    fes = HCurl(mesh, order=order, dirichlet='outer')
    print(f"  HCurl DOFs: {fes.ndof}")

    A = fes.TrialFunction()
    v = fes.TestFunction()

    # nu GridFunction (L2, order=0, per element)
    fes_nu = L2(mesh, order=0)
    nu_gf = GridFunction(fes_nu)

    # Initialize nu: iron gets nu_iron, air gets 1/mu_0
    nu_air = 1.0 / MU_0
    nu_iron_init = nu_ext[0]
    print("  Initializing nu per element...")
    for el in mesh.Elements(VOL):
        if str(el.mat) == 'iron':
            nu_gf.vec[el.nr] = nu_iron_init
        else:
            nu_gf.vec[el.nr] = nu_air

    A_gf = GridFunction(fes)
    A_gf.vec[:] = 0

    B_old_arr = None

    use_iterative = fes.ndof > 500000  # Use pardiso up to 500k DOF

    print(f"\nPicard iteration (solver: {'BDDC+CG' if use_iterative else 'pardiso'}):")
    for it in range(maxiter):
        t_iter = time.time()
        # Correct weak form for reduced A:
        #   (nu * curl A_r, curl v) = ((1/mu_0 - nu) * B_s, curl v)
        # In air: 1/mu_0 - nu = 0 -> no source (correct)
        # In iron: 1/mu_0 - nu != 0 -> source from iron permeability contrast
        #
        # Also add small regularization for gauge uniqueness:
        #   + eps * (A, v)
        eps = 1e-10
        a = BilinearForm(fes)
        a += nu_gf * InnerProduct(curl(A), curl(v)) * dx
        a += eps * InnerProduct(A, v) * dx

        f = LinearForm(fes)
        nu_contrast = 1.0/MU_0 - nu_gf
        f += nu_contrast * InnerProduct(Bs_cf, curl(v)) * dx

        if use_iterative:
            from ngsolve.krylovspace import CGSolver
            c = Preconditioner(a, 'bddc')
            a.Assemble()
            f.Assemble()
            inv = CGSolver(a.mat, c.mat, maxiter=2000, tol=1e-8,
                           printrates=False)
            A_gf.vec.data = inv * f.vec
        else:
            a.Assemble()
            f.Assemble()
            A_gf.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                           inverse='pardiso') * f.vec

        t_solve_iter = time.time() - t_iter

        # Compute B_total = B_s + curl(A_r)
        B_total_cf = Bs_cf + curl(A_gf)

        # Quick check B at origin
        B_orig = B_total_cf(mip_origin)
        Bz_now = B_orig[2]

        # Update nu in iron elements using L2 projection approach
        # Instead of element loop, use GridFunction Set
        B_mag_cf = Norm(B_total_cf)

        # For iron: nu = H(B)/B from interp, for air: nu = 1/mu_0
        # We need per-element evaluation, so loop is necessary
        n_iron = 0
        B_new_list = []

        for el in mesh.Elements(VOL):
            if str(el.mat) == 'iron':
                # Element centroid from vertices
                verts = el.vertices
                nv = len(verts)
                cx, cy, cz = 0.0, 0.0, 0.0
                for v_idx in verts:
                    p = mesh[v_idx].point
                    cx += p[0]; cy += p[1]; cz += p[2]
                cx /= nv; cy /= nv; cz /= nv

                try:
                    mip_el = mesh(cx, cy, cz)
                    B_val = B_total_cf(mip_el)
                    B_mag = float(np.sqrt(B_val[0]**2 + B_val[1]**2 + B_val[2]**2))
                except:
                    B_mag = 0.0

                if B_mag > 1e-10:
                    nu_new = float(nu_interp(B_mag))
                else:
                    nu_new = nu_iron_init

                # Under-relaxation: nu = alpha*nu_new + (1-alpha)*nu_old
                alpha = 0.3
                nu_old = nu_gf.vec[el.nr]
                nu_gf.vec[el.nr] = alpha * nu_new + (1.0 - alpha) * nu_old
                B_new_list.append(B_mag)
                n_iron += 1

        B_new_arr = np.array(B_new_list)

        # Convergence check
        if B_old_arr is not None and len(B_old_arr) == len(B_new_arr):
            max_change = np.max(np.abs(B_new_arr - B_old_arr)) / B_sat
            print(f"   iter {it}: max |dB|/B_sat = {max_change:.2e}, "
                  f"Bz={Bz_now*1e3:.1f} mT, {t_solve_iter:.1f}s")
            if max_change < tol:
                print(f"   Converged at iteration {it}")
                break
        else:
            print(f"   iter {it}: initial solve, "
                  f"Bz={Bz_now*1e3:.1f} mT, {t_solve_iter:.1f}s")

        B_old_arr = B_new_arr.copy()

    # Final result
    B_total_cf = Bs_cf + curl(A_gf)
    B_origin = B_total_cf(mip_origin)
    Bz = B_origin[2]
    B_mag_origin = np.sqrt(sum(b**2 for b in B_origin))

    print(f"\n  B at origin: Bx={B_origin[0]*1e3:.2f}, "
          f"By={B_origin[1]*1e3:.2f}, Bz={B_origin[2]*1e3:.2f} mT")
    print(f"  |B| = {B_mag_origin*1e3:.2f} mT")
    print(f"  Expected (BEM/Tosca): Bz = -995 mT")
    if abs(Bz) > 1e-6:
        print(f"  Error: {abs(Bz*1e3 - (-995)) / 995 * 100:.1f}%")

    return Bz


def main():
    print("=" * 70)
    print("  A-formulation test: curl(nu*curl A) nonlinear magnetostatic")
    print("  Reference: BEM (Tosca-verified) Bz = -995 mT at NI=20000")
    print("=" * 70)

    # Load B-H curve
    H_arr, B_arr = load_bh_curve()
    print(f"B-H curve: {len(H_arr)} points, H=[0, {H_arr[-1]:.0f}] A/m, "
          f"B=[0, {B_arr[-1]:.3f}] T")

    # Build geometry + mesh
    print("\n1. Building geometry (coarse mesh for speed)...")
    mesh = build_geometry_and_mesh(maxh_gap=0.006, maxh_iron=0.020,
                                   maxh_air=0.060)
    mesh.Curve(2)
    print(f"   Mesh: {mesh.nv} vertices, {mesh.ne} elements")
    print(f"   Materials: {sorted(set(mesh.GetMaterials()))}")

    # Run Picard A-formulation (order=1 for faster iteration)
    print("\n2. Picard A-formulation solve...")
    Bz = solve_A_formulation_picard(mesh, ni=20000.0,
                                     bh_data_H=H_arr, bh_data_B=B_arr,
                                     order=1, maxiter=80, tol=1e-4)

    print("\n" + "=" * 70)
    print(f"  RESULT: A-formulation Bz = {Bz*1e3:.2f} mT")
    print(f"  Expected (BEM/Tosca):       -995 mT")
    print(f"  Simkin result:              -834 mT")
    print("=" * 70)


if __name__ == '__main__':
    main()
