"""
C-Type Electromagnet Nonlinear Simulation
Comparison with ELF_MAGIC: nonlinear_20000AT (20000 AT)

ELF_MAGIC reference:
  S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1

Model:
  - C-type yoke with nonlinear B-H curve
  - SQRING coil: 10000 AT (with symmetry -> 20000 AT effective)
  - Symmetry: Image symmetry (+x, -z) using new API
  - Quarter model solved, 13 yoke elements

Updated 2026-01-31: Using new Image symmetry API (TrfMlt removed)
"""

import sys
import os
import numpy as np

# Add Radia path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))
import radia as rad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from meg_to_vol import parse_meg_file


def load_bh_curve(filepath):
    """Load B-H curve from text file."""
    data = np.loadtxt(filepath, comments='#')
    return data.tolist()  # [[H1, B1], [H2, B2], ...]


def get_element_centroid(nodes, node_ids):
    """Calculate centroid of hexahedron."""
    coords = np.array([nodes[nid] for nid in node_ids])
    return coords.mean(axis=0)


def elf_hex_to_radia_vertices(nodes, elf_node_ids):
    """
    Convert ELF hex node list to Radia vertex list.

    ELF node ordering is directly compatible with Radia ObjHexahedron.
    No reordering needed.

    Returns list of 8 vertex coordinates.
    """
    # ELF node ordering is compatible with Radia - use directly
    vertices = [list(nodes[nid]) for nid in elf_node_ids]
    return vertices


def main():
    print("=" * 70)
    print("C-Type Electromagnet Nonlinear Simulation")
    print("Comparison with ELF_MAGIC nonlinear_20000AT")
    print("=" * 70)

    # Initialize Radia
    rad.UtiDelAll()
    rad.FldUnits('m')  # Use meters (SI units)

    # File paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/ELF_magic.meg"
    bh_file = os.path.join(script_dir, "BH.txt")

    # Load B-H curve
    print("\nLoading B-H curve from:", bh_file)
    bh_data = load_bh_curve(bh_file)
    print(f"  {len(bh_data)} data points")
    print(f"  H range: {bh_data[0][0]:.1f} - {bh_data[-1][0]:.1f} A/m")
    print(f"  B range: {bh_data[0][1]:.3f} - {bh_data[-1][1]:.3f} T")

    # Parse MEG file
    print("\nParsing MEG file:", meg_file)
    nodes, elements, scale = parse_meg_file(meg_file)
    print(f"  Nodes: {len(nodes)}")
    print(f"  Elements: {len(elements)}")
    print(f"  Scale: {scale} (mm -> m)")

    # Note: After parsing, coordinates are in meters (MEG file scale applied)
    # No conversion needed since we use FldUnits('m')

    # Separate yoke (MMB8T) and coil (MCL8T) elements
    yoke_elements = [e for e in elements if e['type'] == 'MMB8T']
    coil_elements = [e for e in elements if e['type'] == 'MCL8T']

    print(f"\n  Yoke elements (MMB8T): {len(yoke_elements)}")
    print(f"  Coil elements (MCL8T): {len(coil_elements)}")

    # Create nonlinear material for yoke
    print("\nCreating nonlinear material...")
    mat = rad.MatSatIsoTab(bh_data)

    # Create yoke hexahedral elements
    print("\nCreating yoke geometry...")
    yoke_objs = []
    elem_info = []

    for elem in yoke_elements:
        elem_id = elem['id']
        node_ids = elem['nodes']

        # Get vertices in correct order
        vertices = elf_hex_to_radia_vertices(nodes, node_ids)
        centroid = get_element_centroid(nodes, node_ids)

        print(f"  Creating element {elem_id}...")
        print(f"    Centroid: [{centroid[0]*1000:.1f}, {centroid[1]*1000:.1f}, {centroid[2]*1000:.1f}] mm")

        try:
            hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
            rad.MatApl(hex_obj, mat)
            yoke_objs.append(hex_obj)
            elem_info.append({
                'id': elem_id,
                'obj': hex_obj,
                'centroid': centroid,
                'vertices': vertices
            })
            print(f"    OK (Radia object: {hex_obj})")
        except Exception as e:
            print(f"    FAILED: {e}")

    print(f"\nCreated {len(yoke_objs)} hexahedral yoke elements")

    if len(yoke_objs) == 0:
        print("ERROR: No yoke elements created!")
        return None, None, None

    # Create container for yoke (quarter model only)
    print("\nCreating yoke container...")
    yoke = rad.ObjCnt(yoke_objs)
    print(f"  Yoke container: {yoke}")
    print(f"  Quarter model: {len(yoke_objs)} elements")
    print("  Image symmetry will be applied during Solve: '+x-z'")

    # Create coil using ObjRaceTrk
    # Note: Coil is analytical (not discretized) and provides background field
    # The coil wraps around the yoke with axis along Y (perpendicular to yoke)
    print("\nCreating coil...")
    # Convert mm to m for all dimensions
    # ELF: AA O 0 131.25 26.25, AA SQRING 60 72.5 52.5 35 5 3
    coil_center = [0, 0.13125, 0.02625]  # m (was 131.25 mm, 26.25 mm)
    r_inner = 0.060   # m (was 60 mm)
    r_outer = 0.0725  # m (was 72.5 mm)
    straight_half = 0.0525  # m (was 52.5 mm)
    height = 0.070    # m (was 70 mm = 2 * 35 mm)
    current = 10000.0  # A (10000 AT for quarter model -> 20000 AT effective with symmetry)
    nseg = 8  # number of segments

    print(f"  Coil center: {coil_center} m")
    print(f"  R_inner: {r_inner*1000:.1f} mm, R_outer: {r_outer*1000:.1f} mm")
    print(f"  Straight section: {straight_half*2*1000:.1f} mm")
    print(f"  Height: {height*1000:.1f} mm")
    print(f"  Current: {current} A")
    print(f"  Axis: Z")

    try:
        # ObjRaceTrk(center, radii, lengths, h, nseg, man_auto, axis, j)
        # For Z-axis coil: lengths[0] is in X, lengths[1] is in Y
        coil = rad.ObjRaceTrk(
            coil_center,
            [r_inner, r_outer],
            [straight_half * 2, height],  # [X-extent, Y-extent]
            height,  # h: coil thickness in Z direction
            nseg,    # number of segments
            'man',   # manual mode
            'z',     # axis along Z
            current / (height * (r_outer - r_inner))  # j: current density
        )
        print(f"  Coil created: {coil}")
    except Exception as e:
        print(f"  ObjRaceTrk failed: {e}")
        print("  Trying ObjArcCur...")

        # Use ObjArcCur with 4 radial segments
        r_avg = (r_inner + r_outer) / 2
        try:
            # ObjArcCur(center, radii, phi, h, nseg, man_auto, axis, j)
            coil = rad.ObjArcCur(
                coil_center,
                [r_avg, r_avg],  # radii (inner, outer)
                [0, 360],        # phi angles
                height,          # h: height
                nseg,            # number of phi segments
                'man',           # manual mode
                'z',             # axis along Z
                current / (height * (r_outer - r_inner))  # j: current density
            )
            print(f"  ObjArcCur coil created: {coil}")
        except Exception as e2:
            print(f"  ObjArcCur also failed: {e2}")
            coil = None

    # Create model containing both yoke and coil
    if coil is not None:
        model = rad.ObjCnt([yoke, coil])
        print(f"\nModel container: {model} (yoke + coil)")
    else:
        model = yoke
        print(f"\nModel container: {model} (yoke only, no coil)")

    # Solve with Image symmetry
    # The model must contain the coil for it to contribute field to the yoke
    print("\nSolving nonlinear problem with Image symmetry...")
    print("  Method: BiCGSTAB (1)")
    print("  Precision: 0.001")
    print("  Max iterations: 500")
    print("  Image: '+x-z' (MIMA X symmetric, MIMA -Z antisymmetric)")

    try:
        # Solve the model (not just yoke) so coil field is included
        # Use BiCGSTAB (method 1) with relaxed precision for nonlinear problems
        result = rad.Solve(model, 0.001, 500, 1, image='+x-z')
        print(f"\nSolution result:")
        print(f"  Max |M|: {result[0]:.2f} A/m")
        print(f"  Max |dM|/|M|: {result[1]:.6f}")
    except Exception as e:
        print(f"  Solve failed: {e}")
        print("  Continuing with current magnetization state...")
        result = [0, 0]

    # ELF_MAGIC reference
    elf_results = [
        (1, 1.8836e-02, 2.3029e-02, 2.2495e+00, 2.2497e+00, 19.305),
        (2, 6.3480e-02, -4.5267e-02, -1.4775e+00, 1.4795e+00, 2597.522),
        (3, 4.6444e-02, -1.0082e+00, 1.1721e+00, 1.5468e+00, 1586.327),
        (4, 4.1590e-02, -1.9292e+00, 1.0014e-02, 1.9297e+00, 138.400),
        (5, 4.7680e-02, -9.0421e-01, -7.5272e-01, 1.1775e+00, 9627.043),
        (6, -2.3166e-01, 3.1331e-02, -8.9152e-01, 9.2166e-01, 11391.015),
        (7, -1.9244e-01, 1.7487e-01, -8.1971e-01, 8.5997e-01, 11619.973),
        (8, -5.3199e-01, 1.2006e-01, -1.3566e+00, 1.4621e+00, 2977.398),
        (9, -4.7339e-01, -3.3145e-02, -1.3056e+00, 1.3892e+00, 5225.880),
        (10, 2.7923e-01, 9.6955e-02, -1.0071e+00, 1.0496e+00, 10711.152),
        (11, 2.1037e-01, -8.4530e-02, -1.0593e+00, 1.0833e+00, 10462.834),
        (12, 4.1241e-01, 5.9766e-02, -1.3042e+00, 1.3692e+00, 5897.984),
        (13, 3.5043e-01, -2.7636e-01, -1.3379e+00, 1.4103e+00, 4502.731),
    ]

    # Compare results
    print("\n" + "=" * 70)
    print("Comparison with ELF_MAGIC results")
    print("=" * 70)

    print("\nRadia field at element centroids:")
    print("-" * 90)
    print(f"{'Elem':>4} {'Bx(T)':>10} {'By(T)':>10} {'Bz(T)':>10} {'|B|(T)':>10} {'ELF|B|':>10} {'Diff%':>8}")
    print("-" * 90)

    total_error = 0
    for i, info in enumerate(elem_info):
        elem_id = info['id']
        centroid = info['centroid']

        try:
            B = rad.Fld(model, 'b', list(centroid))
            B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)

            elf_B_mag = elf_results[i][4] if i < len(elf_results) else 0
            diff_pct = 100 * abs(B_mag - elf_B_mag) / elf_B_mag if elf_B_mag > 0 else 0
            total_error += diff_pct

            print(f"{elem_id:>4} {B[0]:>10.4f} {B[1]:>10.4f} {B[2]:>10.4f} {B_mag:>10.4f} {elf_B_mag:>10.4f} {diff_pct:>7.1f}%")
        except Exception as e:
            print(f"{elem_id:>4} Error: {e}")

    print("-" * 90)
    if len(elem_info) > 0:
        print(f"Average |B| difference: {total_error / len(elem_info):.1f}%")

    print("\nELF_MAGIC reference:")
    print("-" * 90)
    print(f"{'Elem':>4} {'Bx(T)':>10} {'By(T)':>10} {'Bz(T)':>10} {'|B|(T)':>10} {'perm':>10}")
    print("-" * 90)
    for elem_id, Bx, By, Bz, B_mag, perm in elf_results:
        print(f"{elem_id:>4} {Bx:>10.4f} {By:>10.4f} {Bz:>10.4f} {B_mag:>10.4f} {perm:>10.1f}")
    print("-" * 90)

    print("\n" + "=" * 70)
    print("Simulation complete")
    print("=" * 70)

    return model, elem_info, elf_results


if __name__ == "__main__":
    main()
