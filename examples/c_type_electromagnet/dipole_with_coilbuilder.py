"""
C-Type Dipole Magnet with CoilBuilder

Demonstrates the complete accelerator magnet workflow:
1. Define coil with CoilBuilder (racetrack, upper + lower via mirror)
2. Create yoke geometry in Cubit
3. Mesh yoke + air with Cubit
4. Compute coil field with Radia Biot-Savart
5. Solve iron magnetization with NGSolve FEM
6. Visualize with GmshPostExport + coil STEP

Coil is NOT meshed - it's a field source only.
Yoke + air are meshed and solved with FEM.

Run with system Python (Cubit + NGSolve access):
    python dipole_with_coilbuilder.py
"""

import sys
import os
import math
import numpy as np

# Paths
# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
if cubit_path not in sys.path:
    sys.path.insert(0, cubit_path)

radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..', '..', 'src', 'radia')
sys.path.insert(0, os.path.abspath(radia_src))

work_dir = os.path.dirname(os.path.abspath(__file__))

# Import NGSolve BEFORE Cubit
from ngsolve import Mesh, Integrate, CF, BND, VOL, dx, GridFunction, H1, BilinearForm, LinearForm
from ngsolve import grad
from ngsolve import TaskManager

print("=" * 60)
print("C-TYPE DIPOLE MAGNET WITH COILBUILDER")
print("=" * 60)

# ============================================================
# Parameters (all in meters)
# ============================================================
mm = 1e-3

# Yoke dimensions
YOKE_WIDTH = 100 * mm       # X direction
YOKE_HEIGHT = 200 * mm      # Y direction (total, including gap)
YOKE_DEPTH = 100 * mm       # Z direction (beam direction)
POLE_WIDTH = 60 * mm        # Pole face width
POLE_HEIGHT = 30 * mm       # Pole height
GAP = 40 * mm               # Air gap between poles

# Coil dimensions
COIL_WIDTH = 15 * mm         # Radial width
COIL_HEIGHT = 25 * mm        # Axial height
COIL_STRAIGHT = 80 * mm      # Straight section length (along beam)
COIL_BEND_R = 20 * mm        # Bend radius
NI = 10000                   # Ampere-turns

# Derived
COIL_CURRENT = NI  # Total current (single-turn model)
COIL_Z_CENTER = GAP / 2 + COIL_HEIGHT / 2  # Coil center z-position

print(f"Gap:      {GAP*1e3:.0f} mm")
print(f"NI:       {NI} A-turns")
print(f"Coil z:   +/-{COIL_Z_CENTER*1e3:.1f} mm")
print()

# ============================================================
# Step 1: Define coil with CoilBuilder
# ============================================================
print("Step 1: Define coil with CoilBuilder")

from coil_builder import CoilBuilder

upper_coil = (CoilBuilder(current=COIL_CURRENT)
    .set_start([0, -COIL_STRAIGHT / 2, COIL_Z_CENTER])
    .set_cross_section(width=COIL_WIDTH, height=COIL_HEIGHT)
    .add_straight(COIL_STRAIGHT)
    .add_arc(radius=COIL_BEND_R, arc_angle=180)
    .add_straight(COIL_STRAIGHT)
    .add_arc(radius=COIL_BEND_R, arc_angle=180))

print(f"  Upper coil: {len(upper_coil.segments)} segments, closed={upper_coil.is_closed}")

lower_coil = upper_coil.mirror('xy')
print(f"  Lower coil: I={lower_coil.current} A (reversed)")

# Export coil STEP for visualization
coil_step = os.path.join(work_dir, "dipole_coils.step")
combined_shape = upper_coil.combined_occ([lower_coil])
combined_shape.WriteStep(coil_step)
print(f"  Coil STEP: {coil_step}")

# ============================================================
# Step 2: Create yoke + air in Cubit
# ============================================================
print("\nStep 2: Create yoke + air geometry in Cubit")

import contextlib, io
with contextlib.redirect_stderr(io.StringIO()):
    import cubit
    cubit.init(['cubit', '-nojournal', '-nographics'])

cubit.cmd("reset")

# Yoke (C-shape): outer box minus gap region
cubit.cmd(f"brick x {YOKE_WIDTH} y {YOKE_DEPTH} z {YOKE_HEIGHT}")
cubit.cmd(f"brick x {POLE_WIDTH} y {YOKE_DEPTH*1.1} z {GAP}")

# Air region (box around everything)
air_size = max(YOKE_WIDTH, YOKE_HEIGHT, YOKE_DEPTH) * 3
cubit.cmd(f"brick x {air_size} y {air_size} z {air_size}")

# Imprint and merge to create shared surfaces
cubit.cmd("imprint all")
cubit.cmd("merge all")

# Subtract gap from yoke (volume 1 - volume 2)
cubit.cmd("subtract volume 2 from volume 1")

# Now: yoke = what's left of vol 1, air = vol 3 minus yoke
# Get volume IDs
all_vols = list(cubit.get_entities("volume"))
print(f"  Volumes after boolean: {all_vols}")

# Identify yoke and air by bounding box size
def vol_bbox_size(vid):
    bb = cubit.get_bounding_box("volume", vid)
    # bb = [xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, diag]
    return bb[2] * bb[5] * bb[8]  # dx * dy * dz

vol_sizes = [(v, vol_bbox_size(v)) for v in all_vols]
vol_sizes.sort(key=lambda x: x[1])
yoke_id = vol_sizes[0][0]  # smaller bbox = yoke
air_id = vol_sizes[-1][0]  # larger bbox = air

print(f"  Yoke = {yoke_id}, Air = {air_id}")

# ============================================================
# Step 3: Mesh
# ============================================================
print("\nStep 3: Mesh")

mesh_size = GAP / 2  # coarser for demo speed
cubit.cmd(f"volume all size {mesh_size}")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

cubit.cmd(f'block 1 add volume {yoke_id}')
cubit.cmd('block 1 name "yoke"')
if air_id:
    cubit.cmd(f'block 2 add volume {air_id}')
    cubit.cmd('block 2 name "air"')
cubit.cmd('block 3 add tri all')
cubit.cmd('block 3 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 4: Export to NGSolve
# ============================================================
print("\nStep 4: Export to NGSolve")

from cubit_mesh_export import extract_curved_mesh
try:
    mesh = Mesh(extract_curved_mesh(cubit, order=2))
except (ImportError, AttributeError):
    # Fallback: order=1 without CallbackGeometry (pip netgen)
    print("  Note: CallbackGeometry not available, using order=1")
    # Use Netgen .vol export as fallback
    vol_tmp = os.path.join(work_dir, "_temp_dipole.vol")
    cubit.cmd(f'export netgen "{vol_tmp}" order 1 overwrite')
    mesh = Mesh(vol_tmp)
    os.remove(vol_tmp)

print(f"  NGSolve mesh: ne={mesh.ne}, nv={mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")
print(f"  Boundaries: {set(mesh.GetBoundaries())}")

# ============================================================
# Step 5: Solve magnetostatic problem
# ============================================================
print("\nStep 5: Solve magnetostatic (scalar potential)")

mu0 = 4 * math.pi * 1e-7
mu_r_iron = 1000

# Material properties
mu = mesh.MaterialCF({"yoke": mu_r_iron * mu0}, default=mu0)

# Coil source field via Radia Biot-Savart
print("  Computing coil source field with Radia...")
import radia as rad
rad.UtiDelAll()

# Create Radia coil objects
radia_upper = upper_coil.to_radia()
radia_lower = lower_coil.to_radia()
all_coils = rad.ObjCnt(radia_upper + radia_lower)

# Source field as CoefficientFunction
Hs = rad.RadiaField(all_coils, 'h')  # H-field from coils

# Scalar potential formulation (reduced region = everywhere, no iron source)
fes = H1(mesh, order=2, dirichlet="boundary")
u, v = fes.TnT()

a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx
with TaskManager():
    a.Assemble()

    f = LinearForm(fes)
    f += mu * Hs * grad(v) * dx
    f.Assemble()

    phi = GridFunction(fes)
    phi.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    print(f"  Solved: {fes.ndof} DOFs")

    # Total H = Hs - grad(phi)
    H_total = Hs - grad(phi)
    B_total = mu * H_total

    # ============================================================
    # Step 6: Evaluate field in gap
    # ============================================================
    print("\nStep 6: Field in gap center")

    # Field at gap center (0, 0, 0)
    B_gap = B_total(mesh(0, 0, 0))
    print(f"  B_gap = ({B_gap[0]:.6f}, {B_gap[1]:.6f}, {B_gap[2]:.6f}) T")
    print(f"  |B_gap| = {math.sqrt(sum(b**2 for b in B_gap)):.6f} T")

    # Expected: B ~ mu0 * NI / gap (rough estimate for C-magnet)
    B_est = mu0 * NI / GAP
    print(f"  Estimate (mu0*NI/gap): {B_est:.6f} T")

    # ============================================================
    # Step 7: Export to GMSH for visualization
    # ============================================================
    print("\nStep 7: GMSH visualization")

    from gmsh_post_export import GmshPostExport

    post = GmshPostExport(mesh)
    post.add_field("|B|", CF(B_total.Norm()), ncomp=1)
    msh_file = os.path.join(work_dir, "dipole_field.msh")
    post.write(msh_file)

    print(f"\nTo visualize:")
    print(f"  gmsh {msh_file} {coil_step}")
    print(f"  (Load both files in GMSH to see field + coil geometry)")

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
