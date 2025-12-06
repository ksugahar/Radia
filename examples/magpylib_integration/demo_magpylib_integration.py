"""
Magpylib + Radia Integration Demo

Simple demonstration of using magpylib permanent magnet as background field
for Radia soft magnetic material computation.

Requirements:
- pip install magpylib>=5.0
- pip install radia (or use local build)
"""

import sys
import os
import numpy as np

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad

try:
    import magpylib as magpy
    print(f"magpylib version: {magpy.__version__}")
except ImportError:
    print("ERROR: magpylib not installed. Run: pip install magpylib")
    sys.exit(1)


def main():
    """Main demo function."""
    print("=" * 70)
    print("Magpylib + Radia Integration Demo")
    print("=" * 70)

    # Use meters for SI consistency
    rad.FldUnits('m')
    rad.UtiDelAll()

    # =========================================================================
    # Step 1: Create magpylib permanent magnet
    # =========================================================================
    print("\n1. Creating magpylib permanent magnet...")

    # Large NdFeB magnet very close to iron cube
    # Magnetization = Br / mu_0 (A/m)
    Br = 1.4  # Tesla (strong NdFeB)
    M_mag = Br / (4 * np.pi * 1e-7)  # Convert to A/m

    magnet = magpy.magnet.Cuboid(
        magnetization=(0, 0, M_mag),  # Magnetized in +z
        dimension=(0.1, 0.1, 0.03),   # 100x100x30mm magnet
        position=(0, 0, 0.035)        # Just 5mm gap above iron cube top
    )

    # Check field at center of iron cube location
    B_at_center = magnet.getB([0, 0, 0])
    print(f"   Magpylib field at origin: [{B_at_center[0]:.4f}, {B_at_center[1]:.4f}, {B_at_center[2]:.4f}] T")
    print(f"   |B| = {np.linalg.norm(B_at_center):.4f} T")

    # =========================================================================
    # Step 2: Create Radia soft iron cube
    # =========================================================================
    print("\n2. Creating Radia soft iron cube...")

    # Iron cube: 40mm side, centered at origin
    cube_size = 0.04  # 40mm
    cube = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])

    # Subdivide for better accuracy
    n_div = 3
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    print(f"   Cube size: {cube_size*1000:.0f}mm, Elements: {n_div**3}")

    # Apply linear material (mu_r = 1000)
    mu_r = 1000
    mat = rad.MatLin(mu_r - 1)  # chi = mu_r - 1
    rad.MatApl(cube, mat)
    print(f"   Material: mu_r = {mu_r}")

    # =========================================================================
    # Step 3: Create background field from magpylib
    # =========================================================================
    print("\n3. Creating Radia background field from magpylib...")

    # Define callback function
    # CRITICAL: ObjBckgCF passes positions in mm (Radia internal units),
    # regardless of FldUnits() setting. magpylib expects meters.
    def magpylib_field(pos):
        """Callback that evaluates magpylib field.

        IMPORTANT: pos is in mm (Radia internal units), must convert to m for magpylib.
        """
        pos_m = [p * 0.001 for p in pos]  # mm -> m
        B = magnet.getB(pos_m)
        return [float(B[0]), float(B[1]), float(B[2])]

    # Test callback (pass mm coordinates as Radia would)
    test_B = magpylib_field([0, 0, 0])  # Origin in mm
    print(f"   Callback test at origin: {test_B}")

    # Create background field object
    bg_field = rad.ObjBckgCF(magpylib_field)

    # =========================================================================
    # Step 4: Assemble and solve
    # =========================================================================
    print("\n4. Solving...")

    # Create system container
    system = rad.ObjCnt([cube, bg_field])

    # Solve the magnetization problem
    result = rad.Solve(system, 0.0001, 1000)

    print(f"   Max |M| change: {result[0]:.6f}")
    print(f"   Iterations: {int(result[3])}")

    # =========================================================================
    # Step 5: Get results
    # =========================================================================
    print("\n5. Results:")

    # Get magnetization of each element
    M_data = rad.ObjM(cube)

    # Compute average magnetization
    M_sum = np.array([0.0, 0.0, 0.0])
    for elem in M_data:
        M_sum += np.array(elem[1])
    M_avg = M_sum / len(M_data)

    print(f"   Average M: [{M_avg[0]:.0f}, {M_avg[1]:.0f}, {M_avg[2]:.0f}] A/m")
    print(f"   |M_avg| = {np.linalg.norm(M_avg):.0f} A/m")

    # =========================================================================
    # Step 6: Compare fields at external points
    # =========================================================================
    print("\n6. Field at external observation points:")

    obs_points = [
        [0, 0, 0.05],     # 50mm above center (above cube)
        [0.05, 0, 0],     # 50mm to the right
        [0, 0, -0.05],    # 50mm below center
    ]

    for pt in obs_points:
        # Magpylib field only
        B_magpy = magnet.getB(pt)

        # Radia iron contribution
        B_iron = rad.Fld(cube, 'b', pt)

        # Total field
        B_total = np.array(B_magpy) + np.array(B_iron)

        print(f"   At {pt}:")
        print(f"      Magnet only: [{B_magpy[0]:.4f}, {B_magpy[1]:.4f}, {B_magpy[2]:.4f}] T")
        print(f"      Iron effect: [{B_iron[0]:.4f}, {B_iron[1]:.4f}, {B_iron[2]:.4f}] T")
        print(f"      Total:       [{B_total[0]:.4f}, {B_total[1]:.4f}, {B_total[2]:.4f}] T")

    # =========================================================================
    # Comparison test with ObjBckg
    # =========================================================================
    print("\n" + "=" * 70)
    print("Comparison: Using ObjBckg (uniform field) instead of ObjBckgCF")
    print("=" * 70)

    rad.UtiDelAll()

    # Recreate cube
    cube2 = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])
    rad.ObjDivMag(cube2, [n_div, n_div, n_div])
    mat2 = rad.MatLin(mu_r - 1)
    rad.MatApl(cube2, mat2)

    # Use ObjBckg with the same field strength as magpylib gives at origin
    B_uniform = B_at_center[2]  # Use z-component
    print(f"\nUsing uniform field: B = {B_uniform:.4f} T in z-direction")

    bg_uniform = rad.ObjBckg([0, 0, B_uniform])
    system2 = rad.ObjCnt([cube2, bg_uniform])

    result2 = rad.Solve(system2, 0.0001, 1000)

    print(f"Max |M| change: {result2[0]:.6f}")
    print(f"Iterations: {int(result2[3])}")

    M_data2 = rad.ObjM(cube2)
    M_sum2 = np.array([0.0, 0.0, 0.0])
    for elem in M_data2:
        M_sum2 += np.array(elem[1])
    M_avg2 = M_sum2 / len(M_data2)

    print(f"Average M: [{M_avg2[0]:.0f}, {M_avg2[1]:.0f}, {M_avg2[2]:.0f}] A/m")
    print(f"|M_avg| = {np.linalg.norm(M_avg2):.0f} A/m")

    print("\n" + "=" * 70)
    print("Demo completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
