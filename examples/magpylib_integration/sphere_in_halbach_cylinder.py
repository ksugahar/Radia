"""
Soft Iron Sphere in External Magnetic Field (magpylib + Radia integration)

This example demonstrates using magpylib to create various magnetic field
configurations as background fields for Radia MMM computation.

Three scenarios are demonstrated:
1. Single magnet - Simple permanent magnet field
2. Helmholtz coil - Uniform axial field
3. Quadrupole - Field gradient configuration

Physical Setup:
- Soft iron cube: 40mm side, centered at origin
- Material: Linear isotropic (mu_r = 1000)

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

# Check magpylib availability
try:
    import magpylib as magpy
    print(f"magpylib version: {magpy.__version__}")
except ImportError:
    print("ERROR: magpylib not installed. Run: pip install magpylib")
    sys.exit(1)


def create_single_magnet(position=(0, 0, 0.1), magnetization_dir=(0, 0, 1),
                         remanence=1.2, size=(0.05, 0.05, 0.02)):
    """
    Create a single permanent magnet using magpylib.

    Parameters:
    -----------
    position : tuple
        Magnet center position [m]
    magnetization_dir : tuple
        Magnetization direction (normalized)
    remanence : float
        Remanent magnetization [T]
    size : tuple
        (width, depth, height) of magnet [m]

    Returns:
    --------
    magpy.magnet.Cuboid
        Single magnet object
    """
    # Convert remanence to magnetization (A/m): M = Br / mu_0
    M = remanence / (4 * np.pi * 1e-7)

    # Normalize direction
    mag_dir = np.array(magnetization_dir)
    mag_dir = mag_dir / np.linalg.norm(mag_dir)

    magnetization = tuple(M * mag_dir)

    magnet = magpy.magnet.Cuboid(
        magnetization=magnetization,
        dimension=size,
        position=position
    )

    return magnet


def create_helmholtz_pair(radius=0.1, separation=0.1, n_turns=100, current=10.0):
    """
    Create a Helmholtz coil pair using magpylib current loops.

    The Helmholtz condition (separation = radius) produces a uniform
    field in the center region.

    Parameters:
    -----------
    radius : float
        Coil radius [m]
    separation : float
        Distance between coils (typically = radius for Helmholtz) [m]
    n_turns : int
        Number of turns per coil
    current : float
        Current per turn [A]

    Returns:
    --------
    magpy.Collection
        Helmholtz coil pair
    """
    coils = []

    # Upper coil
    upper = magpy.current.Circle(
        current=current * n_turns,  # Total ampere-turns
        diameter=2 * radius,
        position=(0, 0, separation / 2)
    )
    coils.append(upper)

    # Lower coil (same current direction for Helmholtz)
    lower = magpy.current.Circle(
        current=current * n_turns,
        diameter=2 * radius,
        position=(0, 0, -separation / 2)
    )
    coils.append(lower)

    return magpy.Collection(coils)


def create_quadrupole_magnets(bore_radius=0.08, magnet_size=(0.04, 0.04, 0.1),
                              remanence=1.2):
    """
    Create a quadrupole magnet configuration using 4 permanent magnets.

    The magnets are arranged to produce a field gradient (dBy/dx = dBx/dy).

    Parameters:
    -----------
    bore_radius : float
        Distance from center to magnet center [m]
    magnet_size : tuple
        (width, depth, height) of each magnet [m]
    remanence : float
        Remanent magnetization [T]

    Returns:
    --------
    magpy.Collection
        Quadrupole magnet array
    """
    M = remanence / (4 * np.pi * 1e-7)  # A/m

    magnets = []

    # Quadrupole configuration:
    # Top: magnetized toward center
    # Bottom: magnetized toward center
    # Right: magnetized toward center
    # Left: magnetized toward center

    # Top magnet (+y position, magnetized toward center = -y)
    magnets.append(magpy.magnet.Cuboid(
        magnetization=(0, -M, 0),
        dimension=magnet_size,
        position=(0, bore_radius, 0)
    ))

    # Bottom magnet (-y position, magnetized toward center = +y)
    magnets.append(magpy.magnet.Cuboid(
        magnetization=(0, M, 0),
        dimension=magnet_size,
        position=(0, -bore_radius, 0)
    ))

    # Right magnet (+x position, magnetized toward center = -x)
    magnets.append(magpy.magnet.Cuboid(
        magnetization=(-M, 0, 0),
        dimension=magnet_size,
        position=(bore_radius, 0, 0)
    ))

    # Left magnet (-x position, magnetized toward center = +x)
    magnets.append(magpy.magnet.Cuboid(
        magnetization=(M, 0, 0),
        dimension=magnet_size,
        position=(-bore_radius, 0, 0)
    ))

    return magpy.Collection(magnets)


def magpylib_to_radia_callback(magpy_source):
    """
    Create a Radia ObjBckgCF callback function from a magpylib source.

    This adapter function converts magpylib field evaluation to the
    format expected by Radia's ObjBckgCF.

    CRITICAL: ObjBckgCF passes positions in mm (Radia internal units),
    regardless of the FldUnits() setting. magpylib expects meters.
    This function handles the mm -> m conversion automatically.

    Parameters:
    -----------
    magpy_source : magpylib source object
        Any magpylib source (Cuboid, Collection, Circle, etc.)

    Returns:
    --------
    callable
        Function: pos [x, y, z] in mm -> B [Bx, By, Bz] in Tesla
    """
    def field_callback(pos):
        """Compute B field at position using magpylib.

        Args:
            pos: Position in mm (Radia internal units)

        Returns:
            B field in Tesla [Bx, By, Bz]
        """
        # Convert mm (Radia internal) to m (magpylib)
        pos_m = np.array(pos) * 0.001
        B = magpy_source.getB(pos_m)
        return [float(B[0]), float(B[1]), float(B[2])]

    return field_callback


def run_example(name, magpy_source, description):
    """
    Run a single example with given magpylib source.

    Parameters:
    -----------
    name : str
        Example name
    magpy_source : magpylib source
        The magnetic field source
    description : str
        Description of the example
    """
    print("\n" + "=" * 70)
    print(f"Example: {name}")
    print(f"Description: {description}")
    print("=" * 70)

    # Use meters for NGSolve compatibility
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Observation points OUTSIDE the iron cube (cube is 40mm = 0.04m centered at origin)
    # Cube extends from -0.02m to +0.02m in each direction
    external_points = [
        [0, 0, 0.05],      # 50mm above center (outside cube)
        [0.05, 0, 0],      # 50mm to the right (outside cube)
        [0, 0.05, 0],      # 50mm forward (outside cube)
        [0.03, 0.03, 0],   # 30mm diagonal (outside cube)
    ]

    # Check field at external points
    print("\n1. Magpylib field at external observation points (no iron):")
    for pt in external_points:
        B = magpy_source.getB(pt)
        B_mag = np.linalg.norm(B)
        print(f"   B at {pt}: [{B[0]:.4f}, {B[1]:.4f}, {B[2]:.4f}] T (|B|={B_mag:.4f} T)")

    # Create Radia background field
    print("\n2. Creating Radia system...")
    field_callback = magpylib_to_radia_callback(magpy_source)
    bg_field = rad.ObjBckgCF(field_callback)

    # Create soft iron cube
    cube_size = 0.04  # 40mm
    cube = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])

    # Subdivide
    n_div = 4
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    print(f"   Iron cube: {cube_size*1000:.0f}mm side, {n_div**3} elements")

    # Apply material
    mu_r = 1000
    mat = rad.MatLin(mu_r)  # relative permeability
    rad.MatApl(cube, mat)
    print(f"   Material: mu_r = {mu_r}")

    # Create system and solve
    system = rad.ObjCnt([cube, bg_field])

    print("\n3. Solving...")
    result = rad.Solve(system, 0.0001, 1000)
    print(f"   Max |M| change: {result[0]:.6f}")
    print(f"   Iterations: {int(result[3])}")

    # Get magnetization
    M_data = rad.ObjM(cube)
    M_sum = np.array([0.0, 0.0, 0.0])
    for elem in M_data:
        M_sum += np.array(elem[1])
    M_avg = M_sum / len(M_data)

    print(f"\n4. Results:")
    print(f"   Average M: [{M_avg[0]:.0f}, {M_avg[1]:.0f}, {M_avg[2]:.0f}] A/m")
    print(f"   |M_avg| = {np.linalg.norm(M_avg):.0f} A/m")

    # Field comparison at EXTERNAL points (Radia field is accurate outside magnets)
    print("\n5. Field comparison at external observation points:")
    print("   (Radia field is most accurate OUTSIDE magnetic materials)")
    for pt in external_points:
        B_with = rad.Fld(cube, 'b', pt)  # Only iron contribution
        B_without = magpy_source.getB(pt)
        B_total = np.array(B_with) + np.array(B_without)
        print(f"   At {pt}:")
        print(f"      Magpylib only:     [{B_without[0]:.4f}, {B_without[1]:.4f}, {B_without[2]:.4f}] T")
        print(f"      Iron contribution: [{B_with[0]:.4f}, {B_with[1]:.4f}, {B_with[2]:.4f}] T")
        print(f"      Total field:       [{B_total[0]:.4f}, {B_total[1]:.4f}, {B_total[2]:.4f}] T")

    return system


def main():
    """Main function running all examples."""
    print("=" * 70)
    print("Magpylib + Radia Integration Examples")
    print("Soft Iron Cube in Various External Fields")
    print("=" * 70)

    # Example 1: Single magnet above iron cube
    magnet = create_single_magnet(
        position=(0, 0, 0.06),  # 60mm above center
        magnetization_dir=(0, 0, -1),  # Magnetized downward
        remanence=1.2,
        size=(0.04, 0.04, 0.02)  # 40x40x20mm
    )
    run_example(
        "Single Permanent Magnet",
        magnet,
        "NdFeB magnet (1.2T) 60mm above soft iron cube"
    )

    # Example 2: Helmholtz coil pair
    helmholtz = create_helmholtz_pair(
        radius=0.1,      # 100mm radius
        separation=0.1,  # Helmholtz condition
        n_turns=100,
        current=5.0      # 5A -> 500 A-turns
    )
    run_example(
        "Helmholtz Coil Pair",
        helmholtz,
        "Uniform axial field from Helmholtz coils (100mm radius, 500 A-turns)"
    )

    # Example 3: Quadrupole magnets
    quadrupole = create_quadrupole_magnets(
        bore_radius=0.08,  # 80mm from center
        magnet_size=(0.03, 0.03, 0.08),  # 30x30x80mm
        remanence=1.0
    )
    run_example(
        "Quadrupole Magnets",
        quadrupole,
        "Field gradient configuration (zero field at center, gradient away)"
    )

    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
