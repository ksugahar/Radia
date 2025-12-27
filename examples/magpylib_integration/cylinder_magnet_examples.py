"""
Cylindrical Magnet Examples (magpylib + Radia integration)

This example demonstrates using magpylib's cylindrical magnet models
as background fields for Radia MMM computation.

Examples included:
1. Axially magnetized cylinder - Simple axial field
2. Diametrically magnetized cylinder - Transverse field
3. Halbach cylinder - Strong uniform field in bore
4. Ring magnet pair - Axial gradient field

Physical Setup:
- Soft iron sphere/cube centered at origin
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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
from netgen_mesh_import import HEX_FACES

# Check magpylib availability
try:
    import magpylib as magpy
    print(f"magpylib version: {magpy.__version__}")
except ImportError:
    print("ERROR: magpylib not installed. Run: pip install magpylib")
    sys.exit(1)


def magpylib_to_radia_callback(magpy_source):
    """
    Create a Radia ObjBckgCF callback function from a magpylib source.

    CRITICAL: ObjBckgCF passes positions in mm (Radia internal units),
    regardless of the FldUnits() setting. magpylib expects meters.
    This function handles the mm -> m conversion automatically.

    Parameters:
    -----------
    magpy_source : magpylib source object
        Any magpylib source (Cylinder, CylinderSegment, Collection, etc.)

    Returns:
    --------
    callable
        Function: pos [x, y, z] in mm -> B [Bx, By, Bz] in Tesla
    """
    def field_callback(pos):
        # Convert mm (Radia internal) to m (magpylib)
        pos_m = np.array(pos) * 0.001
        B = magpy_source.getB(pos_m)
        return [float(B[0]), float(B[1]), float(B[2])]

    return field_callback


def create_axial_cylinder(diameter=0.05, height=0.03, remanence=1.2, position=(0, 0, 0.05)):
    """
    Create an axially magnetized cylinder (magnetization along cylinder axis).

    This is the most common permanent magnet configuration.

    Parameters:
    -----------
    diameter : float
        Cylinder diameter [m]
    height : float
        Cylinder height [m]
    remanence : float
        Remanent magnetization Br [T]
    position : tuple
        Center position [m]

    Returns:
    --------
    magpy.magnet.Cylinder
        Axially magnetized cylinder
    """
    # Convert Br to magnetization M = Br / mu_0
    M = remanence / (4 * np.pi * 1e-7)

    cylinder = magpy.magnet.Cylinder(
        magnetization=(0, 0, M),  # Axial magnetization (+z)
        dimension=(diameter, height),
        position=position
    )

    return cylinder


def create_diametric_cylinder(diameter=0.05, height=0.03, remanence=1.2,
                               position=(0, 0, 0), mag_direction=(1, 0, 0)):
    """
    Create a diametrically magnetized cylinder (magnetization perpendicular to axis).

    Used for transverse field generation and magnetic couplings.

    Parameters:
    -----------
    diameter : float
        Cylinder diameter [m]
    height : float
        Cylinder height [m]
    remanence : float
        Remanent magnetization Br [T]
    position : tuple
        Center position [m]
    mag_direction : tuple
        Magnetization direction in xy-plane (will be normalized)

    Returns:
    --------
    magpy.magnet.Cylinder
        Diametrically magnetized cylinder
    """
    M = remanence / (4 * np.pi * 1e-7)

    # Normalize direction
    mag_dir = np.array(mag_direction)
    mag_dir = mag_dir / np.linalg.norm(mag_dir)

    magnetization = tuple(M * mag_dir)

    cylinder = magpy.magnet.Cylinder(
        magnetization=magnetization,
        dimension=(diameter, height),
        position=position
    )

    return cylinder


def create_halbach_cylinder(inner_radius=0.03, outer_radius=0.06, height=0.08,
                            n_segments=8, remanence=1.2, k=2):
    """
    Create a Halbach cylinder using magpylib CylinderSegment objects.

    Halbach cylinders produce a strong, uniform field inside the bore
    with minimal external field (k=2 for internal dipole field).

    Parameters:
    -----------
    inner_radius : float
        Inner radius (bore radius) [m]
    outer_radius : float
        Outer radius [m]
    height : float
        Cylinder height [m]
    n_segments : int
        Number of segments (more = more uniform field, typically 8-16)
    remanence : float
        Remanent magnetization Br [T]
    k : int
        Halbach parameter (k=2 for dipole, k=4 for quadrupole)

    Returns:
    --------
    magpy.Collection
        Collection of CylinderSegment objects forming Halbach array
    """
    M = remanence / (4 * np.pi * 1e-7)

    segments = []
    angle_step = 2 * np.pi / n_segments

    for i in range(n_segments):
        # Angular position of segment center
        theta = i * angle_step

        # Halbach magnetization direction
        # For k=2: magnetization rotates twice as fast as position angle
        mag_angle = k * theta

        # Magnetization vector (in global coordinates)
        mx = M * np.cos(mag_angle)
        my = M * np.sin(mag_angle)

        # Create segment
        # CylinderSegment: dimension = (r1, r2, h, phi1, phi2)
        phi1 = np.degrees(theta - angle_step/2)
        phi2 = np.degrees(theta + angle_step/2)

        segment = magpy.magnet.CylinderSegment(
            magnetization=(mx, my, 0),
            dimension=(inner_radius, outer_radius, height, phi1, phi2),
            position=(0, 0, 0)
        )
        segments.append(segment)

    return magpy.Collection(segments)


def create_ring_magnet_pair(inner_radius=0.02, outer_radius=0.05, height=0.02,
                             separation=0.06, remanence=1.2, mode='attract'):
    """
    Create a pair of ring magnets with adjustable separation.

    Can be configured for attractive (opposing magnetization) or
    repulsive (same magnetization) configuration.

    Parameters:
    -----------
    inner_radius : float
        Ring inner radius [m]
    outer_radius : float
        Ring outer radius [m]
    height : float
        Ring height [m]
    separation : float
        Center-to-center separation [m]
    remanence : float
        Remanent magnetization Br [T]
    mode : str
        'attract' - opposing magnetization (strong gradient)
        'repel' - same magnetization (weaker gradient)

    Returns:
    --------
    magpy.Collection
        Pair of ring magnets
    """
    M = remanence / (4 * np.pi * 1e-7)

    # Upper ring (magnetized downward for attraction mode)
    if mode == 'attract':
        upper_mag = (0, 0, -M)
        lower_mag = (0, 0, M)
    else:  # repel
        upper_mag = (0, 0, M)
        lower_mag = (0, 0, M)

    upper_ring = magpy.magnet.CylinderSegment(
        magnetization=upper_mag,
        dimension=(inner_radius, outer_radius, height, 0, 360),
        position=(0, 0, separation/2)
    )

    lower_ring = magpy.magnet.CylinderSegment(
        magnetization=lower_mag,
        dimension=(inner_radius, outer_radius, height, 0, 360),
        position=(0, 0, -separation/2)
    )

    return magpy.Collection([upper_ring, lower_ring])


def run_example(name, magpy_source, description, iron_size=0.03):
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
    iron_size : float
        Size of soft iron cube [m]
    """
    print("\n" + "=" * 70)
    print(f"Example: {name}")
    print(f"Description: {description}")
    print("=" * 70)

    # Use meters for consistency
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Observation points outside the iron cube
    half_size = iron_size / 2
    external_points = [
        [0, 0, iron_size],           # Above
        [0, 0, -iron_size],          # Below
        [iron_size, 0, 0],           # Side
        [iron_size/2, iron_size/2, iron_size/2],  # Diagonal
    ]

    # Check field at observation points (without iron)
    print("\n1. Magpylib field at observation points (no iron):")
    for pt in external_points:
        B = magpy_source.getB(pt)
        B_mag = np.linalg.norm(B)
        print(f"   B at [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}]: "
              f"[{B[0]:.4f}, {B[1]:.4f}, {B[2]:.4f}] T (|B|={B_mag:.4f} T)")

    # Field at center (inside bore for Halbach)
    B_center = magpy_source.getB([0, 0, 0])
    print(f"\n   B at center: [{B_center[0]:.4f}, {B_center[1]:.4f}, {B_center[2]:.4f}] T")
    print(f"   |B_center| = {np.linalg.norm(B_center):.4f} T")

    # Create Radia background field
    print("\n2. Creating Radia system...")
    field_callback = magpylib_to_radia_callback(magpy_source)
    bg_field = rad.ObjBckgCF(field_callback)

    # Create soft iron cube: iron_size side, centered at origin
    half = iron_size / 2
    # Hexahedron vertices for cube centered at [0, 0, 0] with dimensions [iron_size, iron_size, iron_size]
    vertices = [
        [-half, -half, -half], [half, -half, -half], [half, half, -half], [-half, half, -half],
        [-half, -half, half], [half, -half, half], [half, half, half], [-half, half, half]
    ]
    cube = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])

    # Subdivide
    n_div = 4
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    print(f"   Iron cube: {iron_size*1000:.0f}mm side, {n_div**3} elements")

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
    M_vectors = np.array([m[1] for m in M_data])
    M_avg = np.mean(M_vectors, axis=0)

    print(f"\n4. Results:")
    print(f"   Average M: [{M_avg[0]:.0f}, {M_avg[1]:.0f}, {M_avg[2]:.0f}] A/m")
    print(f"   |M_avg| = {np.linalg.norm(M_avg):.0f} A/m")

    # Field comparison at external points
    print("\n5. Field with iron vs without:")
    for pt in external_points[:2]:  # Show first two points
        B_iron = rad.Fld(cube, 'b', pt)
        B_magpy = magpy_source.getB(pt)
        B_total = np.array(B_iron) + np.array(B_magpy)

        print(f"   At [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}]:")
        print(f"      Without iron: [{B_magpy[0]:.4f}, {B_magpy[1]:.4f}, {B_magpy[2]:.4f}] T")
        print(f"      With iron:    [{B_total[0]:.4f}, {B_total[1]:.4f}, {B_total[2]:.4f}] T")

    return system


def main():
    """Main function running all cylindrical magnet examples."""
    print("=" * 70)
    print("Cylindrical Magnet Examples")
    print("magpylib + Radia Integration")
    print("=" * 70)

    # Example 1: Axially magnetized cylinder
    axial_cyl = create_axial_cylinder(
        diameter=0.04,     # 40mm diameter
        height=0.03,       # 30mm height
        remanence=1.2,     # NdFeB
        position=(0, 0, 0.05)  # 50mm above center
    )
    run_example(
        "Axially Magnetized Cylinder",
        axial_cyl,
        "40mm dia x 30mm NdFeB cylinder, axially magnetized, 50mm above iron",
        iron_size=0.03
    )

    # Example 2: Diametrically magnetized cylinder
    diametric_cyl = create_diametric_cylinder(
        diameter=0.05,
        height=0.04,
        remanence=1.0,
        position=(0.06, 0, 0),  # 60mm to the side
        mag_direction=(1, 0, 0)  # Magnetized in x-direction
    )
    run_example(
        "Diametrically Magnetized Cylinder",
        diametric_cyl,
        "50mm dia x 40mm cylinder, diametrically magnetized, beside iron",
        iron_size=0.03
    )

    # Example 3: Halbach cylinder
    halbach = create_halbach_cylinder(
        inner_radius=0.025,  # 25mm bore
        outer_radius=0.05,   # 50mm outer
        height=0.06,         # 60mm height
        n_segments=12,       # 12 segments for good uniformity
        remanence=1.2,
        k=2                  # Dipole configuration
    )
    run_example(
        "Halbach Cylinder (k=2 Dipole)",
        halbach,
        "Halbach array: 25mm bore, 50mm OD, 60mm height, 12 segments",
        iron_size=0.02  # Smaller iron inside bore
    )

    # Example 4: Ring magnet pair (magnetic gradient)
    ring_pair = create_ring_magnet_pair(
        inner_radius=0.015,
        outer_radius=0.04,
        height=0.015,
        separation=0.05,
        remanence=1.0,
        mode='attract'
    )
    run_example(
        "Ring Magnet Pair (Gradient Field)",
        ring_pair,
        "Two opposing ring magnets creating axial field gradient",
        iron_size=0.02
    )

    print("\n" + "=" * 70)
    print("All cylindrical magnet examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
